# EXCLUDE FROM TESTING
#=
  ex12.jl — Julia port of PETSc snes/tutorials/ex12.c

  Solves:  -∇·(κ(·)∇u) = f   in  Ω = (0,1)^dim
                      u = g   on  ∂Ω  (Dirichlet)  or
          κ(·)∇u · n = h   on  ∂Ω  (Neumann)

  Exact solution (any dim):  u(x) = Σᵢ xᵢ²

  Options:
    -dim N       spatial dimension, 2 or 3            (default: 3)
    -coeff MODE  coefficient mode:
                   field      κ(x) = x₂+1   linear PDE
                   nonlinear  κ(u) = u+1    nonlinear PDE  (default: field)
    -bc MODE     boundary-condition mode:
                   dirichlet  Dirichlet on all faces  (default)
                   neumann    pure Neumann; constant null space removed

  Body-force derivation (f = ∇·(κ∇u)):
    field:      f = (2·dim+2)·x₂ + 2·dim
    nonlinear:  f = (2·dim+4)·u + 2·dim

  ── Basic usage ──────────────────────────────────────────────────────────────

  Serial, 2-D, linear coefficient (triangles):
    julia --project examples/ex12.jl -dim 2 -coeff field

  Serial, 3-D, default 8³ hex mesh:
    julia --project examples/ex12.jl

  Finer mesh with refinement:
    julia --project examples/ex12.jl -dim 2 -dm_plex_box_faces 16,16 -dm_refine 1

  ── MPI ──────────────────────────────────────────────────────────────────────

  4 ranks, 3-D, field coefficient:
    mpiexec -n 4 julia --project examples/ex12.jl -dim 3 -dm_plex_box_faces 8,8,8

  4 ranks, nonlinear coefficient:
    mpiexec -n 4 julia --project examples/ex12.jl -dim 2 -coeff nonlinear \
            -dm_plex_box_faces 32,32

  ── Boundary conditions ───────────────────────────────────────────────────────

  Dirichlet (default):
    julia --project examples/ex12.jl -dim 2 -bc dirichlet

  Pure Neumann (singular system, solution unique up to a constant):
    julia --project examples/ex12.jl -dim 2 -coeff field -bc neumann

  ── Solver variants (all pass through as PETSc options) ─────────────────────

  GAMG algebraic multigrid (good general choice):
    julia --project examples/ex12.jl -dim 2 -dm_plex_box_faces 32,32 \
          -pc_type gamg -ksp_type cg

  Matrix-free Jacobian (no matrix assembled):
    julia --project examples/ex12.jl -dim 2 -coeff nonlinear -snes_mf

  Matrix-free preconditioner (matrix used only for PC, not for matvec):
    julia --project examples/ex12.jl -dim 2 -coeff nonlinear -snes_mf_operator

  Periodic mesh (no Dirichlet, singular — add null space or use neumann mode):
    julia --project examples/ex12.jl -dim 2 -dm_plex_box_bd periodic,periodic \
          -bc neumann

  BDDC domain decomposition (4 ranks):
    mpiexec -n 4 julia --project examples/ex12.jl -dim 2 \
            -dm_plex_box_faces 6,6 -dm_refine 2 \
            -dm_mat_type is -pc_type bddc -ksp_type gmres

  ASM with ILU subsolvers:
    julia --project examples/ex12.jl -dim 2 -dm_plex_box_faces 16,16 \
          -pc_type asm -pc_asm_blocks 4 -sub_pc_type ilu -ksp_type gmres

  ── What requires an external (non-jll) PETSc build ─────────────────────────

  The following ex12.c tests need external libraries not compiled into PETSc_jll:

    p4est   — adaptive mesh refinement (AMR), p4est-based partitioning,
              BDDC/FAS on forest meshes.
              All tests with -dm_plex_convert_type p4est or -dm_plex_adapt*.

    FAS nonlinear multigrid (-snes_type fas -snes_fas_levels N) requires that
    the FE field and DS be propagated to each coarsened DM via a DM coarsen
    callback (DMCoarsenHookAdd).  This is not implemented in this example;
    see PETSc ex12.c SetupProblem for the full C implementation.

    ctetgen/triangle — 3-D simplex (tetrahedral) mesh generation.
              ctetgen/TetGen is present in the jll but crashes; use
              -dim 2 (triangles via DMPlex) or hex meshes (-dim 3) instead.

    parmetis — parallel graph partitioning (-petscpartitioner_type parmetis).

    hpddm / slepc — multilevel HPDDM preconditioner tests.

    pragmatic / bamg — anisotropic mesh adaptation.

    cgns / hdf5 — mesh I/O tests (-vec_view cgns:... or hdf5:...).
              hdf5 may be available depending on the Julia artifact version.

  To use these, build PETSc from source with the relevant options and load it
  via the PETSC_OPT / PETSC_DEB environment variables (see PETSc.jl docs).

  References:
    - PETSc 3.23  src/snes/tutorials/ex12.c
    - PetscFE/PetscDS weak-form conventions
=#

using MPI
using PETSc
using PETSc: LibPETSc

# ── PETSc setup ───────────────────────────────────────────────────────────────
_args    = isinteractive() ? String[] : collect(String, ARGS)
opts     = PETSc.parse_options(["-snes_monitor"; _args])
petsclib = PETSc.getlib(; PetscScalar = Float64)
MPI.Initialized() || MPI.Init()
PETSc.initialize(petsclib)

const PetscInt    = petsclib.PetscInt
const PetscScalar = Float64
const PetscReal   = Float64

comm = MPI.COMM_WORLD

let _d = get(NamedTuple(pairs(opts)), :dim,   nothing); global dim   = _d === nothing ? 3 : parse(Int, string(_d)); end
let _c = get(NamedTuple(pairs(opts)), :coeff, nothing); global coeff = _c === nothing ? "field" : string(_c);       end
let _b = get(NamedTuple(pairs(opts)), :bc,    nothing); global bc    = _b === nothing ? "dirichlet" : string(_b);   end
@assert dim in (2, 3)                             "-dim must be 2 or 3"
@assert coeff in ("field", "nonlinear")           "-coeff must be field or nonlinear"
@assert bc    in ("dirichlet", "neumann")         "-bc must be dirichlet or neumann"

# ── Common simple-point functions ─────────────────────────────────────────────

# Exact solution  u = Σ xᵢ²  — Dirichlet seed and L² reference.
function quadratic_u(t, x, u, ctx)
    u[1] = sum(xi^2 for xi in x)
end
const quadratic_u_ptr = PETSc.@petsc_simple_fn(quadratic_u)

# Zero initial guess for the interior DOFs.
function zero_u(t, x, u, ctx)
    fill!(u, 0)
end
const zero_u_ptr = PETSc.@petsc_simple_fn(zero_u)

# ── COEFF_FIELD:  κ(x) = x₂ + 1 ─────────────────────────────────────────────
#
# f0 = (2·dim+2)·x₂ + 2·dim   (body force)
# f1 = κ(x)·∇u                 (diffusion flux)
# g3 = κ(x)·I                  (Jacobian of f1 w.r.t. ∇u)

function f0_field(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                  aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, f0)
    f0[1] = (2*dim_+2)*x[2] + 2*dim_
end
const f0_field_ptr = PETSc.@petsc_residual_fn(f0_field, Nf)

function f1_field(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                  aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, f1)
    kap = x[2] + 1.0
    @inbounds for d in 1:dim_; f1[d] = kap * u_x[d]; end
end
const f1_field_ptr = PETSc.@petsc_residual_fn(f1_field, dim_)

function g3_field(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                  aOff, aOff_x, a, a_t, a_x, t, u_tShift, x, numConstants, constants, g3)
    kap = x[2] + 1.0
    fill!(g3, 0)
    @inbounds for d in 1:dim_; g3[(d-1)*dim_ + d] = kap; end
end
const g3_field_ptr = PETSc.@petsc_jacobian_fn(g3_field, dim_*dim_)

# Neumann flux for COEFF_FIELD:  h = κ(x)·∇u·n = (x₂+1)·2(x·n)
# Convention: f0_bd = −h  (PETSc boundary residual sign)
function f0_bd_field(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                     aOff, aOff_x, a, a_t, a_x, t, x, n, numConstants, constants, f0)
    kap  = x[2] + 1.0
    xdotn = zero(eltype(x))
    @inbounds for d in 1:dim_; xdotn += x[d]*n[d]; end
    f0[1] = -kap * 2.0 * xdotn
end
const f0_bd_field_ptr = PETSc.@petsc_bd_fn(f0_bd_field, Nf)

# ── COEFF_NONLINEAR:  κ(u) = u + 1 ──────────────────────────────────────────
#
# PDE is now nonlinear: −∇·((u+1)∇u) = f
#
# Exact solution u = Σxᵢ²  →  κ = Σxᵢ² + 1
# Body force  f = ∇·(κ∇u) = (2·dim+4)·u + 2·dim
#
# Weak-form residual and Jacobian terms:
#   f0[1]    = (2·dim+4)·u + 2·dim        ← body force (depends on u)
#   f1[d]    = (u+1)·∂u/∂xd               ← diffusion flux
#
#   g0[1]    = ∂f0/∂u         = 2·dim+4
#   g2[d]    = ∂f1[d]/∂u      = ∂u/∂xd   (κ′=1)
#   g3[d,d]  = ∂f1[d]/∂(∂u/∂xd) = u+1

function f0_nl(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
               aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, f0)
    u_val = u[Int(uOff[1])+1]
    f0[1] = (2*dim_+4)*u_val + 2*dim_
end
const f0_nl_ptr = PETSc.@petsc_residual_fn(f0_nl, Nf)

function f1_nl(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
               aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, f1)
    kap = u[Int(uOff[1])+1] + 1.0
    @inbounds for d in 1:dim_; f1[d] = kap * u_x[d]; end
end
const f1_nl_ptr = PETSc.@petsc_residual_fn(f1_nl, dim_)

# g0: ∂f0/∂u = 2·dim+4
function g0_nl(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
               aOff, aOff_x, a, a_t, a_x, t, u_tShift, x, numConstants, constants, g0)
    g0[1] = 2*dim_ + 4
end
const g0_nl_ptr = PETSc.@petsc_jacobian_fn(g0_nl, Nf*Nf)

# g2: ∂f1[d]/∂u = u_x[d]   (κ′(u)=1, so ∂((u+1)·u_x[d])/∂u = u_x[d])
function g2_nl(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
               aOff, aOff_x, a, a_t, a_x, t, u_tShift, x, numConstants, constants, g2)
    @inbounds for d in 1:dim_; g2[d] = u_x[d]; end
end
const g2_nl_ptr = PETSc.@petsc_jacobian_fn(g2_nl, Nf*Nf*dim_)

# g3: ∂f1[d]/∂(∂u/∂xd′) = (u+1)·δ_{dd′}
function g3_nl(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
               aOff, aOff_x, a, a_t, a_x, t, u_tShift, x, numConstants, constants, g3)
    kap = u[Int(uOff[1])+1] + 1.0
    fill!(g3, 0)
    @inbounds for d in 1:dim_; g3[(d-1)*dim_ + d] = kap; end
end
const g3_nl_ptr = PETSc.@petsc_jacobian_fn(g3_nl, dim_*dim_)

# Neumann flux for COEFF_NONLINEAR:  h = (u+1)·2(x·n)
# f0_bd = −h
function f0_bd_nl(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                  aOff, aOff_x, a, a_t, a_x, t, x, n, numConstants, constants, f0)
    kap   = u[Int(uOff[1])+1] + 1.0
    xdotn = zero(eltype(x))
    @inbounds for d in 1:dim_; xdotn += x[d]*n[d]; end
    f0[1] = -kap * 2.0 * xdotn
end
const f0_bd_nl_ptr = PETSc.@petsc_bd_fn(f0_bd_nl, Nf)

# ── Build the mesh ────────────────────────────────────────────────────────────
# 3-D simplex (TetGen) is broken in the current PETSc_jll; use hexahedra for 3D.
simplex = dim < 3
let _s = get(NamedTuple(pairs(opts)), :dm_plex_box_faces, nothing)
    global faces = _s === nothing ? fill(8, dim) : parse.(Int, split(string(_s), ","))
end
dm = PETSc.DMPlex(petsclib, comm, dim, simplex, faces; opts...)

# ── Discretization ────────────────────────────────────────────────────────────
fe = PETSc.fe_create_default(petsclib, MPI.COMM_SELF, dim, 1, simplex; prefix = "")
LibPETSc.PetscObjectSetName(petsclib, convert(Ptr{Cvoid}, fe), "potential")
PETSc.setfield!(dm, 0, fe)
PETSc.createds!(dm)

# ── DS: residual / Jacobian / exact solution ──────────────────────────────────
ds = PETSc.getds(dm)

if coeff == "field"
    PETSc.set_residual!(ds, 0, f0_field_ptr, f1_field_ptr)
    PETSc.set_jacobian!(ds, 0, 0, C_NULL, C_NULL, C_NULL, g3_field_ptr)
else  # nonlinear
    PETSc.set_residual!(ds, 0, f0_nl_ptr, f1_nl_ptr)
    PETSc.set_jacobian!(ds, 0, 0, g0_nl_ptr, C_NULL, g2_nl_ptr, g3_nl_ptr)
end
PETSc.set_exact_solution!(ds, 0, quadratic_u_ptr)

# ── Boundary conditions ───────────────────────────────────────────────────────
label = PETSc.getlabel(dm, "marker")

if bc == "dirichlet"
    PETSc.add_boundary!(petsclib, dm, LibPETSc.DM_BC_ESSENTIAL, "wall", label,
                        PetscInt[1], 0, PetscInt[], quadratic_u_ptr)
else  # neumann
    bd_ptr = coeff == "field" ? f0_bd_field_ptr : f0_bd_nl_ptr
    PETSc.add_boundary!(petsclib, dm, LibPETSc.DM_BC_NATURAL, "flux", label,
                        PetscInt[1], 0, PetscInt[], bd_ptr)
end

# ── SNES + matrix + vectors ───────────────────────────────────────────────────
snes = PETSc.SNES(petsclib, comm; opts...)
PETSc.setDM!(snes, dm)
u = PETSc.DMGlobalVec(dm)
J = PETSc.MatAIJ(dm)
PETSc.plex_set_snes_local_fem!(petsclib, dm)
LibPETSc.SNESSetJacobian(petsclib, snes, J, J, C_NULL, C_NULL)

# For pure Neumann: the assembled matrix is singular (null space = constants).
# MatSetNullSpace        → projects null space from the solution each KSP solve.
# MatSetTransposeNullSpace → projects null space from the RHS (b) each KSP solve.
# Both are required for the residual to converge to zero.
if bc == "neumann"
    nsp = LibPETSc.MatNullSpaceCreate(petsclib, comm,
              LibPETSc.PetscBool(true), PetscInt(0), LibPETSc.PetscVec[])
    LibPETSc.MatSetNullSpace(petsclib, J, nsp)
    LibPETSc.MatSetTransposeNullSpace(petsclib, J, nsp)
end

# ── Initial guess ─────────────────────────────────────────────────────────────
# Dirichlet: seed BCs with exact solution, then reset free DOFs to zero so SNES
# has a non-trivial problem to solve while BC nodes are already correct.
# Neumann:   no constrained DOFs exist, so start directly from the L²-projected
# exact solution.  The solution is unique only up to a constant; starting near
# the answer avoids the singular-system pitfall of a pure-zero initial guess.
PETSc.dm_project_function!(petsclib, dm, 0.0, [quadratic_u_ptr], nothing,
                            LibPETSc.INSERT_ALL_VALUES, u)
if bc == "dirichlet"
    PETSc.dm_project_function!(petsclib, dm, 0.0, [zero_u_ptr], nothing,
                                LibPETSc.INSERT_VALUES, u)
end

# ── Solve ─────────────────────────────────────────────────────────────────────
PETSc.solve!(u, snes)

if MPI.Comm_rank(comm) == 0
    its = LibPETSc.SNESGetIterationNumber(petsclib, snes)
    println("SNES converged in $its Newton iterations.")
end

# ── L² error ──────────────────────────────────────────────────────────────────
l2err = PETSc.dm_compute_l2diff(petsclib, dm, 0.0, [quadratic_u_ptr], nothing, u)
if MPI.Comm_rank(comm) == 0
    println("L2 error: $l2err")
end

# ── Cleanup ───────────────────────────────────────────────────────────────────
GC.gc(true)
MPI.Barrier(comm)
PETSc.destroy(snes)
PETSc.destroy(J)
PETSc.destroy(u)
PETSc.destroy(dm)
if !isinteractive()
    PETSc.finalize(petsclib)
    MPI.Barrier(comm)
    MPI.Finalize()
    ccall(:quick_exit, Cvoid, (Cint,), 0)
end
