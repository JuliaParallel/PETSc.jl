# EXCLUDE FROM TESTING
#=
  ex69.jl — Variable-viscosity Stokes Problem (2D) with finite elements.
  Julia port of PETSc snes/tutorials/ex69.c.

  Solves  −∇·(2μ ε(u)) + ∇p = f   (momentum)
           ∇·u                  = 0   (incompressibility)
  on Ω = (0,1)², with free-slip BCs on all walls and zero-mean pressure.
  The saddle-point system is solved with a Schur-complement fieldsplit
  preconditioner; the constant-pressure null space is handled automatically
  via DMSetNullSpaceConstructor.

  Solution types (-sol_type):
    solkx   — exponentially varying viscosity μ = exp(2Bx)   (default)
    solcx   — piecewise-constant μ, jump at x = xc

  Parameters:
    -n N      x-wavenumber for forcing/exact solution   (default 1)
    -m N      z-wavenumber (may be non-integer)         (default 1.0)
    -B val    viscosity exponent for SolKx              (default 1.0)
    -etaA v   low-side viscosity for SolCx              (default 1.0)
    -etaB v   high-side viscosity for SolCx             (default 1.0)
    -xc  v    viscosity jump location for SolCx         (default 0.5)

  Convergence rates (L² combined velocity+pressure, SolKx default params):
    Q2/Q1 quad:  ~4× per refinement  (O(h²) dominated by pressure)
    Q1/P0 quad:  ~2× per refinement  (O(h)  dominated by velocity)

  ── Run examples ─────────────────────────────────────────────────────────────

  Common solver options (Schur-complement fieldsplit, quad mesh):
    SOLVER="-pc_use_amat -pc_type fieldsplit -pc_fieldsplit_type schur \
            -pc_fieldsplit_schur_factorization_type full \
            -pc_fieldsplit_schur_precondition a11 \
            -fieldsplit_velocity_pc_type lu \
            -fieldsplit_pressure_ksp_rtol 1e-9 -fieldsplit_pressure_pc_type lu \
            -ksp_rtol 1e-9"

  Q2/Q1 quad, SolKx (L² ≈ 0.014 on default 2×2 mesh):
    julia --project examples/ex69.jl \
          -dm_plex_simplex 0 -dm_plex_separate_marker \
          -vel_petscspace_degree 2 -pres_petscspace_degree 1 \
          $SOLVER

  Q1/P0 quad, SolKx (L² ≈ 0.062 on default 2×2 mesh):
    julia --project examples/ex69.jl \
          -dm_plex_simplex 0 -dm_plex_separate_marker \
          -vel_petscspace_degree 1 -pres_petscspace_degree 0 \
          $SOLVER

  SolCx with viscosity jump etaB=1000 (L² ≈ 0.039 on 2×2 mesh):
    julia --project examples/ex69.jl \
          -dm_plex_simplex 0 -dm_plex_separate_marker \
          -sol_type solcx -etaB 1e3 \
          -vel_petscspace_degree 2 -pres_petscspace_degree 1 \
          $SOLVER

  Refinement convergence study — Q2/Q1, -dm_refine N gives 2^N × 2^N cells:
    julia --project examples/ex69.jl \
          -dm_plex_simplex 0 -dm_plex_separate_marker -dm_refine 2 \
          -vel_petscspace_degree 2 -pres_petscspace_degree 1 \
          $SOLVER

  MPI (2 ranks, Q2/Q1, L² ≈ 0.0029 after 1 uniform refinement):
    mpiexec -n 2 julia --project examples/ex69.jl \
            -dm_plex_simplex 0 -dm_plex_separate_marker -dm_refine 1 \
            -vel_petscspace_degree 2 -pres_petscspace_degree 1 \
            $SOLVER

  MPI (4 ranks, Q2/Q1, SolCx etaB=1000):
    mpiexec -n 4 julia --project examples/ex69.jl \
            -dm_plex_simplex 0 -dm_plex_separate_marker -dm_refine 1 \
            -sol_type solcx -etaB 1e3 \
            -vel_petscspace_degree 2 -pres_petscspace_degree 1 \
            $SOLVER

  P2/P1 simplex, direct fieldsplit (requires triangle mesher in PETSc_jll):
    julia --project examples/ex69.jl -dm_plex_separate_marker \
          -vel_petscspace_degree 2 -pres_petscspace_degree 1 \
          -pc_use_amat -pc_type fieldsplit -pc_fieldsplit_type schur \
          -pc_fieldsplit_schur_factorization_type full \
          -pc_fieldsplit_schur_precondition a11 \
          -fieldsplit_velocity_pc_type lu \
          -fieldsplit_pressure_ksp_rtol 1e-9 -fieldsplit_pressure_pc_type lu \
          -ksp_rtol 1e-9

  ── C-code test equivalents (from PETSc src/snes/tutorials/ex69.c) ──────────

  q2q1_conv  (L² convergence study, 2 uniform refinements, Q2/Q1 quad):
    julia --project examples/ex69.jl \
          -dm_plex_simplex 0 -dm_plex_separate_marker \
          -vel_petscspace_degree 2 -pres_petscspace_degree 1 \
          -snes_convergence_estimate -convest_num_refine 2 \
          -ksp_rtol 1e-9 -pc_use_amat \
          -pc_type fieldsplit -pc_fieldsplit_type schur \
          -pc_fieldsplit_schur_factorization_type full \
          -pc_fieldsplit_schur_precondition a11 \
          -fieldsplit_velocity_pc_type lu \
          -fieldsplit_pressure_ksp_rtol 1e-9 -fieldsplit_pressure_pc_type lu

  q1p0_conv  (L² convergence study, 2 uniform refinements, Q1/P0 quad):
    julia --project examples/ex69.jl \
          -dm_plex_simplex 0 -dm_plex_separate_marker \
          -vel_petscspace_degree 1 -pres_petscspace_degree 0 \
          -snes_convergence_estimate -convest_num_refine 2 \
          -ksp_rtol 1e-9 -pc_use_amat \
          -pc_type fieldsplit -pc_fieldsplit_type schur \
          -pc_fieldsplit_schur_factorization_type full \
          -pc_fieldsplit_schur_precondition a11 \
          -fieldsplit_velocity_pc_type lu \
          -fieldsplit_pressure_ksp_rtol 1e-9 -fieldsplit_pressure_pc_type lu

  q2p1_conv  (Q2 velocity, discontinuous P1 pressure, 2 uniform refinements):
    julia --project examples/ex69.jl \
          -dm_plex_simplex 0 -dm_plex_separate_marker \
          -vel_petscspace_degree 2 -pres_petscspace_degree 1 \
          -pres_petscspace_poly_tensor 0 \
          -pres_petscdualspace_lagrange_continuity 0 \
          -snes_convergence_estimate -convest_num_refine 2 \
          -ksp_rtol 1e-9 -pc_use_amat \
          -pc_type fieldsplit -pc_fieldsplit_type schur \
          -pc_fieldsplit_schur_factorization_type full \
          -pc_fieldsplit_schur_precondition a11 \
          -fieldsplit_velocity_pc_type lu \
          -fieldsplit_pressure_ksp_rtol 1e-9 -fieldsplit_pressure_pc_type lu

  References:
    - PETSc 3.23  src/snes/tutorials/ex69.c
    - Velic et al. (2006), "A Fast and Robust Algorithm for Solving the
      Stokes Equations with Strong Discontinuous Viscosity"
=#

include("ex69_exact.jl")

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
libname = petsclib.petsc_library   # library path for direct ccall

# ── Option parsing ────────────────────────────────────────────────────────────
_o       = NamedTuple(pairs(opts))
sol_type = string(get(_o, :sol_type, "solkx"))
n_val    = parse(Int,     string(get(_o, :n,    "1")))
m_val    = parse(Float64, string(get(_o, :m,    "1.0")))
B_val    = parse(Float64, string(get(_o, :B,    "1.0")))
etaA_val = parse(Float64, string(get(_o, :etaA, "1.0")))
etaB_val = parse(Float64, string(get(_o, :etaB, "1.0")))
xc_val   = parse(Float64, string(get(_o, :xc,   "0.5")))
_simp       = get(_o, Symbol("dm_plex_simplex"), nothing)
simplex     = _simp === nothing ? true : (parse(Int, string(_simp)) != 0)
vel_degree  = parse(Int, string(get(_o, Symbol("vel_petscspace_degree"),  "1")))
pres_degree = parse(Int, string(get(_o, Symbol("pres_petscspace_degree"), "1")))

@assert sol_type in ("solkx", "solcx") "-sol_type must be solkx or solcx"

# ── Pointwise functions: residuals ────────────────────────────────────────────
#
# Field layout in u[]:  u[1..dim] = velocity (field 0), u[dim+1] = pressure (field 1).
# u_x layout: u_x[c*dim+d+1] = ∂u_c/∂x_d  (all fields concatenated)
#
# f0_u: body force (rhs of momentum).  cst[1]=m, cst[2]=n.
function f0_stokes_u(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                     aOff, aOff_x, a, a_t, a_x, t, x, nC, cst, f0)
    f0[1] = 0.0
    f0[2] = -sin(real(cst[2]) * pi * x[2]) * cos(real(cst[1]) * pi * x[1])
end
const f0_stokes_u_ptr = PETSc.@petsc_residual_fn(f0_stokes_u, dim_)

# Momentum flux σ = 2μ ε(u) − pI  for SolKx (μ = exp(2Bx), cst[3]=B)
function stokes_momentum_kx(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                             aOff, aOff_x, a, a_t, a_x, t, x, nC, cst, f1)
    mu = exp(2.0 * real(cst[3]) * x[1])
    for c in 0:dim_-1
        for d in 0:dim_-1
            f1[c*dim_+d+1] = mu * (u_x[c*dim_+d+1] + u_x[d*dim_+c+1])
        end
        f1[c*dim_+c+1] -= u[dim_+1]      # subtract pressure
    end
end
const stokes_momentum_kx_ptr = PETSc.@petsc_residual_fn(stokes_momentum_kx, dim_*dim_)

# Momentum flux for SolCx (μ = etaA or etaB, cst[3]=etaA, cst[4]=etaB, cst[5]=xc)
function stokes_momentum_cx(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                             aOff, aOff_x, a, a_t, a_x, t, x, nC, cst, f1)
    mu = x[1] < real(cst[5]) ? real(cst[3]) : real(cst[4])
    for c in 0:dim_-1
        for d in 0:dim_-1
            f1[c*dim_+d+1] = mu * (u_x[c*dim_+d+1] + u_x[d*dim_+c+1])
        end
        f1[c*dim_+c+1] -= u[dim_+1]
    end
end
const stokes_momentum_cx_ptr = PETSc.@petsc_residual_fn(stokes_momentum_cx, dim_*dim_)

# Continuity: f0 = −∇·u
function stokes_mass(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                     aOff, aOff_x, a, a_t, a_x, t, x, nC, cst, f0)
    f0[1] = 0.0
    for d in 0:dim_-1
        f0[1] -= u_x[d*dim_+d+1]
    end
end
const stokes_mass_ptr = PETSc.@petsc_residual_fn(stokes_mass, 1)

# Zero f1 for pressure field
function f1_zero_p(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                   aOff, aOff_x, a, a_t, a_x, t, x, nC, cst, f1)
    for d in 1:dim_; f1[d] = 0.0; end
end
const f1_zero_p_ptr = PETSc.@petsc_residual_fn(f1_zero_p, dim_)

# ── Pointwise functions: Jacobians ────────────────────────────────────────────

# J_uu (g3): 2μ (∇v : ∇u + ∇v : ∇uᵀ)
function stokes_vel_J_kx(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                          aOff, aOff_x, a, a_t, a_x, t, utShift, x, nC, cst, g3)
    mu = exp(2.0 * real(cst[3]) * x[1])
    for cI in 0:dim_-1, d in 0:dim_-1
        g3[((cI*dim_+cI)*dim_+d)*dim_+d+1] += mu
        g3[((cI*dim_+d)*dim_+d)*dim_+cI+1] += mu
    end
end
const stokes_vel_J_kx_ptr = PETSc.@petsc_jacobian_fn(stokes_vel_J_kx, dim_*dim_*dim_*dim_)

function stokes_vel_J_cx(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                          aOff, aOff_x, a, a_t, a_x, t, utShift, x, nC, cst, g3)
    mu = x[1] < real(cst[5]) ? real(cst[3]) : real(cst[4])
    for cI in 0:dim_-1, d in 0:dim_-1
        g3[((cI*dim_+cI)*dim_+d)*dim_+d+1] += mu
        g3[((cI*dim_+d)*dim_+d)*dim_+cI+1] += mu
    end
end
const stokes_vel_J_cx_ptr = PETSc.@petsc_jacobian_fn(stokes_vel_J_cx, dim_*dim_*dim_*dim_)

# J_up (g2): −⟨∇·v, p⟩  →  g2[d*dim+d] = -1
function stokes_pres_J(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                        aOff, aOff_x, a, a_t, a_x, t, utShift, x, nC, cst, g2)
    for d in 0:dim_-1; g2[d*dim_+d+1] = -1.0; end
end
const stokes_pres_J_ptr = PETSc.@petsc_jacobian_fn(stokes_pres_J, dim_*dim_)

# J_pu (g1): ⟨q, ∇·u⟩  →  g1[d*dim+d] = -1
function stokes_mass_J(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                        aOff, aOff_x, a, a_t, a_x, t, utShift, x, nC, cst, g1)
    for d in 0:dim_-1; g1[d*dim_+d+1] = -1.0; end
end
const stokes_mass_J_ptr = PETSc.@petsc_jacobian_fn(stokes_mass_J, dim_*dim_)

# J_pp preconditioner (g0): (1/μ)⟨q, p⟩
function stokes_id_J_kx(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                          aOff, aOff_x, a, a_t, a_x, t, utShift, x, nC, cst, g0)
    g0[1] = 1.0 / exp(2.0 * real(cst[3]) * x[1])
end
const stokes_id_J_kx_ptr = PETSc.@petsc_jacobian_fn(stokes_id_J_kx, 1)

function stokes_id_J_cx(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                          aOff, aOff_x, a, a_t, a_x, t, utShift, x, nC, cst, g0)
    mu = x[1] < real(cst[5]) ? real(cst[3]) : real(cst[4])
    g0[1] = 1.0 / mu
end
const stokes_id_J_cx_ptr = PETSc.@petsc_jacobian_fn(stokes_id_J_cx, 1)

# ── Exact solution / BC callbacks ─────────────────────────────────────────────
# PetscSimplePointFn: f(t, x, u, ctx) → nothing

function zero_vel(t, x, u, ctx)
    for c in 1:length(u); u[c] = 0.0; end
end
const zero_vel_ptr = PETSc.@petsc_simple_fn(zero_vel)

function one_pres(t, x, u, ctx)
    u[1] = 1.0
end
const one_pres_ptr = PETSc.@petsc_simple_fn(one_pres)

# SolKx exact solution (captures parameters via module-level globals)
function exact_vel_kx(t, x, u, ctx)
    vx, vz, _ = SolKxSolution(x, m_val, n_val, B_val)
    u[1] = vx; u[2] = vz
end
const exact_vel_kx_ptr = PETSc.@petsc_simple_fn(exact_vel_kx)

function exact_pres_kx(t, x, u, ctx)
    _, _, p = SolKxSolution(x, m_val, n_val, B_val)
    u[1] = p
end
const exact_pres_kx_ptr = PETSc.@petsc_simple_fn(exact_pres_kx)

# SolCx exact solution
function exact_vel_cx(t, x, u, ctx)
    vx, vz, _ = SolCxSolution(x, m_val, n_val, xc_val, etaA_val, etaB_val)
    u[1] = vx; u[2] = vz
end
const exact_vel_cx_ptr = PETSc.@petsc_simple_fn(exact_vel_cx)

function exact_pres_cx(t, x, u, ctx)
    _, _, p = SolCxSolution(x, m_val, n_val, xc_val, etaA_val, etaB_val)
    u[1] = p
end
const exact_pres_cx_ptr = PETSc.@petsc_simple_fn(exact_pres_cx)

exact_vel_ptr  = sol_type == "solkx" ? exact_vel_kx_ptr  : exact_vel_cx_ptr
exact_pres_ptr = sol_type == "solkx" ? exact_pres_kx_ptr : exact_pres_cx_ptr

# ── Pressure null space constructor ───────────────────────────────────────────
# Registered via DMSetNullSpaceConstructor so PETSc automatically attaches the
# constant-pressure null space to the pressure block during matrix setup.
# Signature: (DM, orig_field, field, MatNullSpace*) → PetscErrorCode
function pressure_nsp_constructor(
    dm_ptr::Ptr{Cvoid}, ::PetscInt, ::PetscInt,
    nsp_pp::Ptr{LibPETSc.MatNullSpace},
)::PetscInt
    PL   = typeof(petsclib)
    dm_w = LibPETSc.PetscDM{PL}(dm_ptr)
    # Project constant pressure = 1 onto the DM's global vector
    nvec = PETSc.DMGlobalVec(dm_w)
    PETSc.dm_project_function!(petsclib, dm_w, 0.0,
        [zero_vel_ptr, one_pres_ptr], nothing, LibPETSc.INSERT_ALL_VALUES, nvec)
    LibPETSc.VecNormalize(petsclib, nvec)
    # Build MatNullSpace from the normalized vector
    vec_ptrs = [nvec.ptr]
    GC.@preserve nvec begin
        nsp_ref = Ref{LibPETSc.MatNullSpace}(C_NULL)
        LibPETSc.@chk ccall((:MatNullSpaceCreate, libname),
            LibPETSc.PetscErrorCode,
            (LibPETSc.MPI_Comm, LibPETSc.PetscBool, PetscInt,
             Ptr{LibPETSc.CVec}, Ptr{LibPETSc.MatNullSpace}),
            MPI.COMM_WORLD, LibPETSc.PETSC_FALSE, PetscInt(1), vec_ptrs, nsp_ref)
        unsafe_store!(nsp_pp, nsp_ref[])
    end
    PETSc.destroy(nvec)
    return PetscInt(0)
end
const pressure_nsp_ptr = Base.@cfunction(pressure_nsp_constructor, PetscInt,
    (Ptr{Cvoid}, PetscInt, PetscInt, Ptr{LibPETSc.MatNullSpace}))

# ── create_split_labels! ──────────────────────────────────────────────────────
# Splits the single "marker" label into one per wall so corners can belong to
# multiple labels.  Face marker IDs: bottom=1, right=2, top=3, left=4.
function create_split_labels!(dm_ptr::Ptr{Cvoid})
    names = ["markerBottom", "markerRight", "markerTop", "markerLeft"]
    ids   = PetscInt[1, 2, 3, 4]
    for (name, id) in zip(names, ids)
        LibPETSc.@chk ccall((:DMCreateLabel, libname),
            LibPETSc.PetscErrorCode, (Ptr{Cvoid}, Cstring), dm_ptr, name)

        is_ref = Ref{Ptr{Cvoid}}(C_NULL)
        LibPETSc.@chk ccall((:DMGetStratumIS, libname),
            LibPETSc.PetscErrorCode,
            (Ptr{Cvoid}, Cstring, PetscInt, Ptr{Ptr{Cvoid}}),
            dm_ptr, "marker", id, is_ref)
        is_ptr = is_ref[]
        is_ptr == C_NULL && continue

        label_ref = Ref{Ptr{Cvoid}}(C_NULL)
        LibPETSc.@chk ccall((:DMGetLabel, libname),
            LibPETSc.PetscErrorCode,
            (Ptr{Cvoid}, Cstring, Ptr{Ptr{Cvoid}}),
            dm_ptr, name, label_ref)

        LibPETSc.@chk ccall((:DMLabelInsertIS, libname),
            LibPETSc.PetscErrorCode,
            (Ptr{Cvoid}, Ptr{Cvoid}, PetscInt),
            label_ref[], is_ptr, PetscInt(1))

        LibPETSc.@chk ccall((:ISDestroy, libname),
            LibPETSc.PetscErrorCode, (Ptr{Ptr{Cvoid}},), is_ref)
    end
end

# ── Mesh ──────────────────────────────────────────────────────────────────────
dm = PETSc.DMPlex(petsclib, comm; opts...)

# Create split per-wall labels on the main DM and all coarser DMs
let cdm = dm
    while convert(Ptr{Cvoid}, cdm) != C_NULL
        create_split_labels!(convert(Ptr{Cvoid}, cdm))
        cdm = PETSc.dm_get_coarse(cdm)
    end
end

# ── Discretisation ────────────────────────────────────────────────────────────
dim = LibPETSc.DMGetDimension(petsclib, dm)
@assert dim == 2 "ex69.jl only supports 2D"

fe_vel  = PETSc.fe_create_default(petsclib, MPI.COMM_SELF, dim, dim, simplex;
                                   degree = vel_degree, prefix = "vel_")
LibPETSc.PetscObjectSetName(petsclib, convert(Ptr{Cvoid}, fe_vel), "velocity")

fe_pres = PETSc.fe_create_default(petsclib, MPI.COMM_SELF, dim, 1, simplex;
                                   degree = pres_degree, prefix = "pres_")
LibPETSc.PetscObjectSetName(petsclib, convert(Ptr{Cvoid}, fe_pres), "pressure")

# Copy quadrature from velocity to pressure for consistent integration
LibPETSc.@chk ccall((:PetscFECopyQuadrature, libname),
    LibPETSc.PetscErrorCode,
    (Ptr{Cvoid}, Ptr{Cvoid}),
    convert(Ptr{Cvoid}, fe_vel), convert(Ptr{Cvoid}, fe_pres))

PETSc.setfield!(dm, 0, fe_vel)
PETSc.setfield!(dm, 1, fe_pres)
PETSc.createds!(dm)

# ── PetscDS setup ─────────────────────────────────────────────────────────────
ds = PETSc.getds(dm)

if sol_type == "solkx"
    PETSc.set_constants!(ds, [m_val, Float64(n_val), B_val])
    PETSc.set_residual!(ds, 0, f0_stokes_u_ptr,   stokes_momentum_kx_ptr)
    PETSc.set_residual!(ds, 1, stokes_mass_ptr,    f1_zero_p_ptr)
    PETSc.set_jacobian!(ds, 0, 0, C_NULL, C_NULL, C_NULL, stokes_vel_J_kx_ptr)
    PETSc.set_jacobian!(ds, 0, 1, C_NULL, C_NULL, stokes_pres_J_ptr, C_NULL)
    PETSc.set_jacobian!(ds, 1, 0, C_NULL, stokes_mass_J_ptr, C_NULL, C_NULL)
    PETSc.set_jacobian_preconditioner!(ds, 0, 0, C_NULL, C_NULL, C_NULL, stokes_vel_J_kx_ptr)
    PETSc.set_jacobian_preconditioner!(ds, 1, 1, stokes_id_J_kx_ptr, C_NULL, C_NULL, C_NULL)
    PETSc.set_exact_solution!(ds, 0, exact_vel_kx_ptr)
    PETSc.set_exact_solution!(ds, 1, exact_pres_kx_ptr)
else
    PETSc.set_constants!(ds, [m_val, Float64(n_val), etaA_val, etaB_val, xc_val])
    PETSc.set_residual!(ds, 0, f0_stokes_u_ptr,   stokes_momentum_cx_ptr)
    PETSc.set_residual!(ds, 1, stokes_mass_ptr,    f1_zero_p_ptr)
    PETSc.set_jacobian!(ds, 0, 0, C_NULL, C_NULL, C_NULL, stokes_vel_J_cx_ptr)
    PETSc.set_jacobian!(ds, 0, 1, C_NULL, C_NULL, stokes_pres_J_ptr, C_NULL)
    PETSc.set_jacobian!(ds, 1, 0, C_NULL, stokes_mass_J_ptr, C_NULL, C_NULL)
    PETSc.set_jacobian_preconditioner!(ds, 0, 0, C_NULL, C_NULL, C_NULL, stokes_vel_J_cx_ptr)
    PETSc.set_jacobian_preconditioner!(ds, 1, 1, stokes_id_J_cx_ptr, C_NULL, C_NULL, C_NULL)
    PETSc.set_exact_solution!(ds, 0, exact_vel_cx_ptr)
    PETSc.set_exact_solution!(ds, 1, exact_pres_cx_ptr)
end

# ── Boundary conditions (free-slip) ───────────────────────────────────────────
# Free-slip = zero normal velocity on each wall.  The exact solution has exactly
# zero normal velocity at every wall, so the exact solution function is used as BC.
# Component index (0-based): 0 = v_x,  1 = v_z.
for (wall, comp) in (("markerBottom", 1), ("markerRight", 0),
                     ("markerTop",    1), ("markerLeft",  0))
    label = PETSc.getlabel(dm, wall)
    PETSc.add_boundary!(petsclib, dm, LibPETSc.DM_BC_ESSENTIAL, wall, label,
                        PetscInt[1], 0, PetscInt[comp], exact_vel_ptr)
end

# Propagate discretisation and pressure null space constructor to coarser levels
let cdm = dm
    while convert(Ptr{Cvoid}, cdm) != C_NULL
        PETSc.dm_copy_disc!(dm, cdm)
        LibPETSc.DMSetNullSpaceConstructor(petsclib, cdm, PetscInt(1), pressure_nsp_ptr)
        cdm = PETSc.dm_get_coarse(cdm)
    end
end

# Attach a constant (trivial) null space to the pressure FE object.
# This informs the fieldsplit preconditioner that the pressure block has a
# constant null space, enabling proper Schur complement approximations.
let fe_pres_ptr = convert(Ptr{Cvoid}, fe_pres)
    nsp_ref = Ref{LibPETSc.MatNullSpace}(C_NULL)
    LibPETSc.@chk ccall((:MatNullSpaceCreate, libname),
        LibPETSc.PetscErrorCode,
        (LibPETSc.MPI_Comm, LibPETSc.PetscBool, PetscInt,
         Ptr{LibPETSc.CVec}, Ptr{LibPETSc.MatNullSpace}),
        comm, LibPETSc.PETSC_TRUE, PetscInt(0), C_NULL, nsp_ref)
    LibPETSc.@chk ccall((:PetscObjectCompose, libname),
        LibPETSc.PetscErrorCode,
        (Ptr{Cvoid}, Cstring, Ptr{Cvoid}),
        fe_pres_ptr, "nullspace", nsp_ref[])
    LibPETSc.@chk ccall((:MatNullSpaceDestroy, libname),
        LibPETSc.PetscErrorCode,
        (Ptr{LibPETSc.MatNullSpace},), nsp_ref)
end

# ── Pressure null space (normalized constant-pressure mode) ──────────────────
null_vec = PETSc.DMGlobalVec(dm)
PETSc.dm_project_function!(petsclib, dm, 0.0,
    [zero_vel_ptr, one_pres_ptr], nothing, LibPETSc.INSERT_ALL_VALUES, null_vec)
LibPETSc.VecNormalize(petsclib, null_vec)

# Create null space from the normalized vector via direct ccall
nullspace_ref = Ref{LibPETSc.MatNullSpace}(C_NULL)
GC.@preserve null_vec begin
    vec_ptrs = [null_vec.ptr]   # Vector{CVec} = Vector{Ptr{Cvoid}}
    LibPETSc.@chk ccall((:MatNullSpaceCreate, libname),
        LibPETSc.PetscErrorCode,
        (LibPETSc.MPI_Comm, LibPETSc.PetscBool, PetscInt,
         Ptr{LibPETSc.CVec}, Ptr{LibPETSc.MatNullSpace}),
        comm, LibPETSc.PETSC_FALSE, PetscInt(1), vec_ptrs, nullspace_ref)
end
nullspace = nullspace_ref[]

# ── SNES + linear algebra ────────────────────────────────────────────────────
# Do NOT call SNESSetJacobian with J,J — that forces a single matrix for both
# Amat and Pmat, preventing PETSc from separately assembling the preconditioner
# Jacobian (PetscDSSetJacobianPreconditioner).  Instead, let DMPlexSetSNESLocalFEM
# + SNESSetUp create separate Amat and Pmat internally from the DM.
snes = PETSc.SNES(petsclib, comm; opts...)
PETSc.setDM!(snes, dm)
u = PETSc.DMGlobalVec(dm)
PETSc.plex_set_snes_local_fem!(petsclib, dm)

# ── Initial guess: zero ───────────────────────────────────────────────────────
PETSc.dm_project_function!(petsclib, dm, 0.0,
    [zero_vel_ptr, zero_vel_ptr], nothing, LibPETSc.INSERT_VALUES, u)

# ── Solve ─────────────────────────────────────────────────────────────────────
# The Stokes system has a constant-pressure null space.  We must attach it to
# the Jacobian AFTER SNESSetUp assembles the matrix but BEFORE SNESSolve.
# PETSc.solve! does SetFromOptions + SNESSolve in one shot, so we replicate the
# relevant part of solve! manually to insert the MatSetNullSpace in between.
push!(snes.opts)
try
    LibPETSc.SNESSetFromOptions(petsclib, snes)
    LibPETSc.SNESSetUp(petsclib, snes)
    # Get the actual Jacobian PETSc will use after setup and attach the null space.
    # CMat = Ptr{Cvoid}; SNESGetJacobian takes Mat* (i.e. Ptr{CMat} = Ptr{Ptr{Cvoid}}).
    J_ref = Ref{LibPETSc.CMat}(C_NULL)
    LibPETSc.@chk ccall((:SNESGetJacobian, libname),
        LibPETSc.PetscErrorCode,
        (LibPETSc.CSNES, Ptr{LibPETSc.CMat}, Ptr{LibPETSc.CMat},
         Ptr{Cvoid}, Ptr{Ptr{Cvoid}}),
        snes, J_ref, C_NULL, C_NULL, C_NULL)
    LibPETSc.@chk ccall((:MatSetNullSpace, libname),
        LibPETSc.PetscErrorCode,
        (LibPETSc.CMat, LibPETSc.MatNullSpace),
        J_ref[], nullspace)
    LibPETSc.SNESSolve(petsclib, snes, C_NULL, u)
finally
    pop!(snes.opts)
end

if MPI.Comm_rank(comm) == 0
    its = LibPETSc.SNESGetIterationNumber(petsclib, snes)
    println("SNES converged in $its Newton iteration(s).")
end

# ── L² errors ────────────────────────────────────────────────────────────────
# Combined L²-norm of the error across both velocity and pressure fields.
l2err = PETSc.dm_compute_l2diff(petsclib, dm, 0.0,
    [exact_vel_ptr, exact_pres_ptr], nothing, u)

if MPI.Comm_rank(comm) == 0
    println("L2 error: $l2err")
end

# ── Cleanup ───────────────────────────────────────────────────────────────────
GC.gc(true)
MPI.Barrier(comm)
LibPETSc.@chk ccall((:MatNullSpaceDestroy, libname),
    LibPETSc.PetscErrorCode, (Ptr{LibPETSc.MatNullSpace},), nullspace_ref)
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
