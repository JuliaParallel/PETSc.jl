# EXCLUDE FROM TESTING
#=
  ex62.jl — Isoviscous Stokes Problem (2D/3D) with finite elements.
  Julia port of PETSc snes/tutorials/ex62.c.

  Solves  −∇·(2μ ε(u)) + ∇p = f   (momentum)
           ∇·u                = 0   (incompressibility)
  on Ω = (0,1)ⁿ (n=2 or 3) with constant viscosity μ and Dirichlet BCs
  on all walls using manufactured exact solutions.

  Solution types (-sol):
    quadratic  — polynomial MMS, exact for P2/Q2 velocity  (default)
    trig       — trigonometric MMS, tests convergence rate

  Quadratic MMS:
    2D: u = [x²+y², 2x²−2xy],  p = x+y−1,       f = (4μ−1) [1,1]
    3D: u = [2x²+y²+z², 2x²−2xy, 2x²−2xz],  p = x+y+z−3/2,
        f = [8μ−1, 4μ−1, 4μ−1]

  Trigonometric MMS:
    2D: u = [sin(πx)+sin(πy), −π cos(πx)y],  p = sin(2πx)+sin(2πy)
    3D: u = [2sin(πx)+sin(πy)+sin(πz), −π cos(πx)y, −π cos(πx)z],
        p = sin(2πx)+sin(2πy)+sin(2πz)

  Parameters:
    -mu val    dynamic shear viscosity (default 1.0)
    -sol name  solution type: quadratic (default) or trig

  ── Run examples ─────────────────────────────────────────────────────────────

  Common solver options (Schur-complement fieldsplit):
    SOLVER="-pc_use_amat -pc_type fieldsplit -pc_fieldsplit_type schur \
            -pc_fieldsplit_schur_factorization_type full \
            -pc_fieldsplit_schur_precondition a11 \
            -fieldsplit_velocity_pc_type lu \
            -fieldsplit_pressure_ksp_rtol 1e-9 -fieldsplit_pressure_pc_type lu \
            -ksp_rtol 1e-9"

  Note: ex62 uses the single "marker" label (id=1) covering the whole
  boundary — do NOT pass -dm_plex_separate_marker.

  2D P2/P1 simplex, quadratic MMS (L² ≈ machine eps — exact for P2):
    julia --project examples/ex62.jl -sol quadratic \
          -vel_petscspace_degree 2 -pres_petscspace_degree 1 \
          $SOLVER

  2D Q2/Q1 quad, trig MMS (L² ≈ 0.65 on default 2×2 mesh):
    julia --project examples/ex62.jl -sol trig \
          -dm_plex_simplex 0 \
          -vel_petscspace_degree 2 -pres_petscspace_degree 1 \
          $SOLVER

  3D Q2/Q1 hex, quadratic MMS (L² ≈ machine eps):
    julia --project examples/ex62.jl -sol quadratic \
          -dm_plex_dim 3 -dm_plex_simplex 0 -dm_plex_box_faces 3,3,3 \
          -vel_petscspace_degree 2 -pres_petscspace_degree 1 \
          $SOLVER

  Refinement convergence study — trig MMS, Q2/Q1, dm_refine N:
    julia --project examples/ex62.jl -sol trig \
          -dm_plex_simplex 0 -dm_refine 2 \
          -vel_petscspace_degree 2 -pres_petscspace_degree 1 \
          $SOLVER

  ── C-code test equivalents (from PETSc src/snes/tutorials/ex62.c) ──────────

  2d_q2_q1_conv  (L₂ convergence study, 2 uniform refinements, Q2/Q1 quad):
    julia --project examples/ex62.jl -sol trig \
          -dm_plex_simplex 0 -dm_plex_separate_marker \
          -vel_petscspace_degree 2 -pres_petscspace_degree 1 \
          -snes_convergence_estimate -convest_num_refine 2 \
          -ksp_atol 1e-10 -pc_use_amat \
          -pc_type fieldsplit -pc_fieldsplit_type schur \
          -pc_fieldsplit_schur_factorization_type full \
          -pc_fieldsplit_schur_precondition a11 \
          -fieldsplit_velocity_pc_type lu \
          -fieldsplit_pressure_ksp_rtol 1e-10 -fieldsplit_pressure_pc_type lu

  2d_q2_q1_check  (Jacobian/residual check with quadratic MMS):
    julia --project examples/ex62.jl -sol quadratic \
          -dm_plex_simplex 0 -dm_plex_separate_marker \
          -vel_petscspace_degree 2 -pres_petscspace_degree 1 \
          -dmsnes_check 0.0001 \
          -pc_use_amat -pc_type fieldsplit -pc_fieldsplit_type schur \
          -pc_fieldsplit_schur_factorization_type full \
          -pc_fieldsplit_schur_precondition a11 \
          -fieldsplit_velocity_pc_type lu \
          -fieldsplit_pressure_ksp_rtol 1e-9 -fieldsplit_pressure_pc_type lu \
          -ksp_rtol 1e-9

  References:
    - PETSc 3.23  src/snes/tutorials/ex62.c
=#

using MPI
using PETSc
using PETSc: LibPETSc

# ── PETSc / MPI init ─────────────────────────────────────────────────────────
_args    = isinteractive() ? String[] : collect(String, ARGS)
opts     = PETSc.parse_options(_args)
petsclib = PETSc.getlib(; PetscScalar = Float64)
MPI.Initialized() || MPI.Init()
PETSc.initialize(petsclib)

const PetscInt    = petsclib.PetscInt
const PetscScalar = Float64
const PetscReal   = Float64

comm = MPI.COMM_WORLD

# ── Option parsing ────────────────────────────────────────────────────────────
_o       = NamedTuple(pairs(opts))
sol_type = string(get(_o, :sol,  "quadratic"))
mu_val   = parse(Float64, string(get(_o, :mu, "1.0")))
_simp    = get(_o, Symbol("dm_plex_simplex"), nothing)
simplex  = _simp === nothing ? true : (parse(Int, string(_simp)) != 0)
vel_degree  = parse(Int, string(get(_o, Symbol("vel_petscspace_degree"),  "1")))
pres_degree = parse(Int, string(get(_o, Symbol("pres_petscspace_degree"), "1")))

@assert sol_type in ("quadratic", "trig") "-sol must be quadratic or trig"

# ── Pointwise functions: shared momentum flux f1_u ────────────────────────────
#
# f1[c*dim+d] = μ(∂u_c/∂x_d + ∂u_d/∂x_c) − p δ_{cd}
# constants[1] = μ
#
function f1_u(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
              aOff, aOff_x, a, a_t, a_x, t, x, nC, cst, f1)
    mu = real(cst[1])
    for c in 0:dim_-1
        for d in 0:dim_-1
            f1[c*dim_+d+1] = mu * (u_x[c*dim_+d+1] + u_x[d*dim_+c+1])
        end
        f1[c*dim_+c+1] -= u[dim_+1]   # subtract pressure (field 1)
    end
end
const f1_u_ptr = PETSc.@petsc_residual_fn(f1_u, dim_*dim_)

# Continuity: f0 = −∇·u
function f0_p(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
              aOff, aOff_x, a, a_t, a_x, t, x, nC, cst, f0)
    f0[1] = 0.0
    for d in 0:dim_-1
        f0[1] -= u_x[d*dim_+d+1]
    end
end
const f0_p_ptr = PETSc.@petsc_residual_fn(f0_p, 1)

# ── Body-force functions (one per MMS solution type) ──────────────────────────
#
# Quadratic MMS body force:
#   2D: f = [1 − 4μ, 1 − 4μ]   (negated → added as rhs, sign convention)
#   3D: f = [1 − 8μ, 1 − 4μ, 1 − 4μ]
# Note: PETSc weak form sign: < v, f > on the RHS means f0 = −f_body
#   so f0[1] = (dim−1)*4μ − 1, f0[d≥2] = 4μ − 1
#
function f0_quadratic_u(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                         aOff, aOff_x, a, a_t, a_x, t, x, nC, cst, f0)
    mu = real(cst[1])
    f0[1] = (dim_ - 1) * 4mu - 1.0
    for d in 2:dim_
        f0[d] = 4mu - 1.0
    end
end
const f0_quadratic_u_ptr = PETSc.@petsc_residual_fn(f0_quadratic_u, dim_)

# Trigonometric MMS body force:
#   f_x = 2π cos(2πx) + μ(dim−1)π² sin(πx) + μπ² Σ_{d>0} sin(πx_d)
#   f_d = 2π cos(2πx_d) − μπ³ cos(πx) x_d   for d > 0
#
function f0_trig_u(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                   aOff, aOff_x, a, a_t, a_x, t, x, nC, cst, f0)
    mu = real(cst[1])
    f0[1] = -2π*cos(2π*x[1]) - (dim_-1)*mu*π^2*sin(π*x[1])
    for d in 2:dim_
        f0[1] -= mu*π^2*sin(π*x[d])
        f0[d]  = -2π*cos(2π*x[d]) + mu*π^3*cos(π*x[1])*x[d]
    end
end
const f0_trig_u_ptr = PETSc.@petsc_residual_fn(f0_trig_u, dim_)

# ── Jacobians ─────────────────────────────────────────────────────────────────

# J_uu (g3): μ(∇v : ∇u + ∇v : ∇uᵀ)
function g3_uu(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
               aOff, aOff_x, a, a_t, a_x, t, utShift, x, nC, cst, g3)
    mu = real(cst[1])
    for c in 0:dim_-1, d in 0:dim_-1
        g3[((c*dim_+c)*dim_+d)*dim_+d+1] += mu
        g3[((c*dim_+d)*dim_+d)*dim_+c+1] += mu
    end
end
const g3_uu_ptr = PETSc.@petsc_jacobian_fn(g3_uu, dim_*dim_*dim_*dim_)

# J_up (g2): −⟨∇·v, p⟩
function g2_up(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
               aOff, aOff_x, a, a_t, a_x, t, utShift, x, nC, cst, g2)
    for d in 0:dim_-1; g2[d*dim_+d+1] = -1.0; end
end
const g2_up_ptr = PETSc.@petsc_jacobian_fn(g2_up, dim_*dim_)

# J_pu (g1): ⟨q, −∇·u⟩
function g1_pu(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
               aOff, aOff_x, a, a_t, a_x, t, utShift, x, nC, cst, g1)
    for d in 0:dim_-1; g1[d*dim_+d+1] = -1.0; end
end
const g1_pu_ptr = PETSc.@petsc_jacobian_fn(g1_pu, dim_*dim_)

# J_pp preconditioner (g0): (1/μ)⟨q, p⟩
function g0_pp(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
               aOff, aOff_x, a, a_t, a_x, t, utShift, x, nC, cst, g0)
    g0[1] = 1.0 / real(cst[1])
end
const g0_pp_ptr = PETSc.@petsc_jacobian_fn(g0_pp, 1)

# ── Exact solutions ───────────────────────────────────────────────────────────
#
# Quadratic MMS:
#   2D: u=[x²+y², 2x²−2xy],  p=x+y−1
#   3D: u=[2x²+y²+z², 2x²−2xy, 2x²−2xz],  p=x+y+z−3/2
#
function quadratic_vel(t, x, u, ctx)
    dim = length(x)
    u[1] = (dim - 1) * x[1]^2
    for c in 2:dim
        u[1] += x[c]^2
        u[c]  = 2*x[1]^2 - 2*x[1]*x[c]
    end
end
const quadratic_vel_ptr = PETSc.@petsc_simple_fn(quadratic_vel)

function quadratic_pres(t, x, u, ctx)
    u[1] = sum(x) - 0.5*length(x)
end
const quadratic_pres_ptr = PETSc.@petsc_simple_fn(quadratic_pres)

#
# Trigonometric MMS:
#   2D: u=[sin(πx)+sin(πy), −π cos(πx)y],  p=sin(2πx)+sin(2πy)
#   3D: u=[2sin(πx)+sin(πy)+sin(πz), −π cos(πx)y, −π cos(πx)z],
#       p=sin(2πx)+sin(2πy)+sin(2πz)
#
function trig_vel(t, x, u, ctx)
    dim = length(x)
    u[1] = (dim - 1) * sin(π*x[1])
    for c in 2:dim
        u[1] += sin(π*x[c])
        u[c]  = -π*cos(π*x[1])*x[c]
    end
end
const trig_vel_ptr = PETSc.@petsc_simple_fn(trig_vel)

function trig_pres(t, x, u, ctx)
    u[1] = sum(sin(2π*xi) for xi in x)
end
const trig_pres_ptr = PETSc.@petsc_simple_fn(trig_pres)

# Helper functions for the pressure null space constructor
function zero_vel(t, x, u, ctx)
    for c in 1:length(u); u[c] = 0.0; end
end
const zero_vel_ptr = PETSc.@petsc_simple_fn(zero_vel)

function one_pres(t, x, u, ctx)
    u[1] = 1.0
end
const one_pres_ptr = PETSc.@petsc_simple_fn(one_pres)

# ── Select exact solution functions based on sol_type ─────────────────────────
exact_vel_ptr  = sol_type == "quadratic" ? quadratic_vel_ptr  : trig_vel_ptr
exact_pres_ptr = sol_type == "quadratic" ? quadratic_pres_ptr : trig_pres_ptr

# ── Pressure null space constructor ───────────────────────────────────────────
function pressure_nsp_constructor(
    dm_ptr::Ptr{Cvoid}, ::PetscInt, ::PetscInt,
    nsp_pp::Ptr{LibPETSc.MatNullSpace},
)::PetscInt
    PL   = typeof(petsclib)
    dm_w = LibPETSc.PetscDM{PL}(dm_ptr)
    nvec = PETSc.DMGlobalVec(dm_w)
    PETSc.dm_project_function!(petsclib, dm_w, 0.0,
        [zero_vel_ptr, one_pres_ptr], nothing, LibPETSc.INSERT_ALL_VALUES, nvec)
    LibPETSc.VecNormalize(petsclib, nvec)
    GC.@preserve nvec begin
        nsp = PETSc.mat_null_space_create(petsclib, MPI.COMM_WORLD, (nvec,))
        unsafe_store!(nsp_pp, nsp)
    end
    PETSc.destroy(nvec)
    return PetscInt(0)
end
const pressure_nsp_ptr = Base.@cfunction(pressure_nsp_constructor, PetscInt,
    (Ptr{Cvoid}, PetscInt, PetscInt, Ptr{LibPETSc.MatNullSpace}))

# ── Mesh ──────────────────────────────────────────────────────────────────────
dm = PETSc.DMPlex(petsclib, comm; opts...)

# ── Discretisation ────────────────────────────────────────────────────────────
dim     = LibPETSc.DMGetDimension(petsclib, dm)
simplex = PETSc.isplexsimplex(dm)

fe_vel  = PETSc.fe_create_default(petsclib, MPI.COMM_SELF, dim, dim, simplex;
                                   degree = vel_degree, prefix = "vel_")
LibPETSc.PetscObjectSetName(petsclib, convert(Ptr{Cvoid}, fe_vel), "velocity")

fe_pres = PETSc.fe_create_default(petsclib, MPI.COMM_SELF, dim, 1, simplex;
                                   degree = pres_degree, prefix = "pres_")
LibPETSc.PetscObjectSetName(petsclib, convert(Ptr{Cvoid}, fe_pres), "pressure")

# Copy quadrature from velocity to pressure for consistent integration
PETSc.fe_copy_quadrature!(petsclib, fe_vel, fe_pres)

PETSc.setfield!(dm, 0, fe_vel)
PETSc.setfield!(dm, 1, fe_pres)
PETSc.createds!(dm)

# ── PetscDS setup ─────────────────────────────────────────────────────────────
ds = PETSc.getds(dm)

PETSc.set_constants!(ds, [mu_val])

if sol_type == "quadratic"
    PETSc.set_residual!(ds, 0, f0_quadratic_u_ptr, f1_u_ptr)
else
    PETSc.set_residual!(ds, 0, f0_trig_u_ptr, f1_u_ptr)
end
PETSc.set_residual!(ds, 1, f0_p_ptr, C_NULL)

PETSc.set_jacobian!(ds, 0, 0, C_NULL, C_NULL, C_NULL, g3_uu_ptr)
PETSc.set_jacobian!(ds, 0, 1, C_NULL, C_NULL, g2_up_ptr, C_NULL)
PETSc.set_jacobian!(ds, 1, 0, C_NULL, g1_pu_ptr, C_NULL, C_NULL)
PETSc.set_jacobian_preconditioner!(ds, 0, 0, C_NULL, C_NULL, C_NULL, g3_uu_ptr)
PETSc.set_jacobian_preconditioner!(ds, 1, 1, g0_pp_ptr, C_NULL, C_NULL, C_NULL)

PETSc.set_exact_solution!(ds, 0, exact_vel_ptr)
PETSc.set_exact_solution!(ds, 1, exact_pres_ptr)

# ── Boundary condition: Dirichlet on all walls (marker label, id=1) ───────────
# All velocity components are constrained; pressure is free.
label = PETSc.getlabel(dm, "marker")
PETSc.add_boundary!(petsclib, dm, LibPETSc.DM_BC_ESSENTIAL, "wall", label,
                    PetscInt[1], 0, PetscInt[], exact_vel_ptr)

# ── Propagate discretisation and null space constructor to coarser levels ─────
let cdm = dm
    while convert(Ptr{Cvoid}, cdm) != C_NULL
        PETSc.dm_copy_disc!(dm, cdm)
        LibPETSc.DMSetNullSpaceConstructor(petsclib, cdm, PetscInt(1), pressure_nsp_ptr)
        cdm = PETSc.dm_get_coarse(cdm)
    end
end

# Attach constant null space to the pressure FE for the fieldsplit preconditioner
PETSc.fe_compose_constant_null_space!(petsclib, comm, fe_pres)

# ── Pressure null space (normalized constant-pressure mode) ───────────────────
null_vec = PETSc.DMGlobalVec(dm)
PETSc.dm_project_function!(petsclib, dm, 0.0,
    [zero_vel_ptr, one_pres_ptr], nothing, LibPETSc.INSERT_ALL_VALUES, null_vec)
LibPETSc.VecNormalize(petsclib, null_vec)
nullspace = GC.@preserve null_vec PETSc.mat_null_space_create(petsclib, comm, (null_vec,))

# ── SNES + solve ─────────────────────────────────────────────────────────────
snes = PETSc.SNES(petsclib, comm; opts...)
PETSc.setDM!(snes, dm)
u = PETSc.DMGlobalVec(dm)
PETSc.plex_set_snes_local_fem!(petsclib, dm)

push!(snes.opts)
try
    LibPETSc.SNESSetFromOptions(petsclib, snes)
    LibPETSc.SNESSetUp(petsclib, snes)
    PETSc.snes_set_jacobian_null_space!(snes, nullspace)
    LibPETSc.SNESSolve(petsclib, snes, C_NULL, u)
finally
    pop!(snes.opts)
end

if MPI.Comm_rank(comm) == 0
    its = LibPETSc.SNESGetIterationNumber(petsclib, snes)
    println("SNES converged in $its Newton iteration(s).")
end

# ── L² errors ────────────────────────────────────────────────────────────────
l2err = PETSc.dm_compute_l2diff(petsclib, dm, 0.0,
    [exact_vel_ptr, exact_pres_ptr], nothing, u)

if MPI.Comm_rank(comm) == 0
    println("L2 error: $l2err")
end

# ── Cleanup ───────────────────────────────────────────────────────────────────
GC.gc(true)
MPI.Barrier(comm)
PETSc.mat_null_space_destroy!(petsclib, nullspace)
PETSc.destroy(snes)
PETSc.destroy(u)
PETSc.destroy(null_vec)
PETSc.destroy(dm)

if !isinteractive()
    flush(stdout); flush(stderr)
    PETSc.finalize(petsclib)
    MPI.Barrier(comm)
    MPI.Finalize()
    ccall(:quick_exit, Cvoid, (Cint,), 0)
end
