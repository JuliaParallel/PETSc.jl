# EXCLUDE FROM TESTING
#=
  ex12.jl — Julia port of PETSc snes/tutorials/ex12.c (simplest variant)

  2D Poisson with Dirichlet BCs on the unit square using DMPLEX + P1 simplex
  finite elements:

      -Δu + f = 0   in Ω = (0,1)²,   i.e.  Δu = f = 4
             u = g  on ∂Ω

  with exact solution  u(x,y) = x² + y²,  f = 4,  g = u|_∂Ω.

  This is the minimal slice of ex12.c covering bcType = DIRICHLET,
  variableCoefficient = COEFF_NONE, runType = RUN_FULL.  Variable coefficients,
  Neumann/Periodic BCs, FAS, mat-free Jacobian etc. are deliberately omitted.

  Usage:
    julia --project ex12.jl
    julia --project ex12.jl -dm_plex_box_faces 16,16 -snes_monitor \
                            -ksp_monitor -dm_refine 1

  References:
    - PETSc 3.24  src/snes/tutorials/ex12.c
    - PetscFE/PetscDS weak-form conventions
=#

using MPI
using PETSc
using PETSc: LibPETSc

# ── PETSc setup ──────────────────────────────────────────────────────────────
_args    = isinteractive() ? String[] : collect(String, ARGS)
opts     = PETSc.parse_options(["-snes_monitor"; _args])
petsclib = PETSc.getlib(; PetscScalar = Float64)
MPI.Initialized() || MPI.Init()
PETSc.initialize(petsclib)

const PetscInt    = petsclib.PetscInt
const PetscScalar = Float64
const PetscReal   = Float64

comm = MPI.COMM_WORLD
dim  = 2

# ── Pointwise functions ──────────────────────────────────────────────────────
#
# PETSc's PetscSimplePointFn signature (used by DMProjectFunction / DMAddBoundary):
#   PetscErrorCode f(PetscInt dim, PetscReal t, const PetscReal x[],
#                    PetscInt Nc, PetscScalar u[], void *ctx)

# Exact solution  u(x,y) = x² + y²  — used as Dirichlet BC, initial seed,
# and reference for the L² error.
function quadratic_u_2d(t, x, u, ctx)
    u[1] = x[1]^2 + x[2]^2
end
const quadratic_u_2d_ptr = PETSc.@petsc_simple_fn(quadratic_u_2d)

# Zero function — initial guess for the free (interior) DOFs.
function zero_u(t, x, u, ctx)
    fill!(u, 0)
end
const zero_u_ptr = PETSc.@petsc_simple_fn(zero_u)

# ── Weak-form point functions (Poisson, f = 4) ───────────────────────────────
#
# PETSc residual convention: F_i = ∫ φ_i f0 dx + ∫ ∇φ_i · f1 dx = 0
# For -Δu + f = 0  (i.e. ∫∇v·∇u + ∫v·f = 0 after IBP):
#   f0 = +4  (the body force term, same sign as in ex12.c)
#   f1 = ∇u  (the diffusion flux)
#
# Jacobian:  g3_{ij} = ∂f1_i/∂(∂u/∂x_j) = δ_{ij}  (identity)
#
# Pointer signature: PetscPointFn (no u_tShift) for residual,
#                    PetscPointJacFn (has u_tShift) for Jacobian.

# Residual f0 = 4 (body force), f1 = ∇u (diffusion flux)
function f0_u(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
              aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, f0)
    f0 .= 4
end
const f0_u_ptr = PETSc.@petsc_residual_fn(f0_u, Nf)

function f1_u(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
              aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, f1)
    f1 .= u_x[1:dim_]
end
const f1_u_ptr = PETSc.@petsc_residual_fn(f1_u, dim_)

function g3_uu(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
               aOff, aOff_x, a, a_t, a_x, t, u_tShift, x, numConstants, constants, g3)
    # g3 is a (dim × dim) tensor stored row-major; identity matrix.
    fill!(g3, 0)
    @inbounds for d in 1:dim_
        g3[(d-1)*dim_ + d] = 1
    end
end
const g3_uu_ptr = PETSc.@petsc_jacobian_fn(g3_uu, dim_*dim_)

# ── Build the mesh ───────────────────────────────────────────────────────────
# Allow -dm_plex_box_faces n,n from the command line to override the resolution.
# C equivalent: DMCreate + DMSetType(DMPLEX) + DMSetFromOptions.
let _s = get(NamedTuple(pairs(opts)), :dm_plex_box_faces, nothing)
    global faces = _s === nothing ? [8, 8] : parse.(Int, split(string(_s), ","))
end
dm = PETSc.DMPlex(petsclib, comm, dim, true, faces; opts...)

# ── Discretization: P1 Lagrange on simplices, single field ──────────────────
# NOTE: PetscFECreateDefault must be given PETSC_COMM_SELF (= MPI.COMM_SELF),
# not the world communicator.  The FE basis is a purely local object — using
# COMM_WORLD produces wrong DOF layouts in some PETSc configurations.
fe = PETSc.fe_create_default(petsclib, MPI.COMM_SELF, dim, 1, true; prefix = "")
LibPETSc.PetscObjectSetName(petsclib, convert(Ptr{Cvoid}, fe), "potential")

PETSc.setfield!(dm, 0, fe)
PETSc.createds!(dm)

# ── Set residual / Jacobian point functions on the DS ────────────────────────
ds = PETSc.getds(dm)
PETSc.set_residual!(ds, 0, f0_u_ptr, f1_u_ptr)
PETSc.set_jacobian!(ds, 0, 0, C_NULL, C_NULL, C_NULL, g3_uu_ptr)
# Register exact solution so DMComputeL2Diff / DMComputeExactSolution work.
PETSc.set_exact_solution!(ds, 0, quadratic_u_2d_ptr)

# ── Dirichlet BC on the "marker" boundary label ─────────────────────────────
# DMPlexCreateBoxMesh tags the entire domain boundary with label "marker", id 1.
let
    label = PETSc.getlabel(dm, "marker")
    PETSc.add_boundary!(
        petsclib, dm,
        LibPETSc.DM_BC_ESSENTIAL,
        "wall",
        label,
        PetscInt[1],   # boundary value id(s)
        0,             # field index (0-based)
        PetscInt[],    # comps — empty = all components
        quadratic_u_2d_ptr,
    )
end

# ── SNES + solution vector ───────────────────────────────────────────────────
snes = PETSc.SNES(petsclib, comm; opts...)
PETSc.setDM!(snes, dm)

u = PETSc.DMGlobalVec(dm)

J = PETSc.MatAIJ(dm)

# Wire DMPlex's built-in FEM residual/Jacobian assembly into SNES, then set
# the Jacobian matrices.  (NULL function in SNESSetJacobian preserves the
# callback already installed by plex_set_snes_local_fem!.)
PETSc.plex_set_snes_local_fem!(petsclib, dm)
LibPETSc.SNESSetJacobian(petsclib, snes, J, J, C_NULL, C_NULL)

# ── Initial guess: BC values pre-seeded, interior = 0 ────────────────────────
# Step 1: project the exact solution into u with INSERT_ALL_VALUES — this fills
#         both the free DOFs and the constrained (Dirichlet) DOFs.  PETSc then
#         knows the correct BC values even before the first Newton iteration.
PETSc.dm_project_function!(petsclib, dm, 0.0, [quadratic_u_2d_ptr], nothing,
                           LibPETSc.INSERT_ALL_VALUES, u)
# Step 2: reset only the free (interior) DOFs to zero so that SNES actually has
#         something to solve — INSERT_VALUES leaves constrained DOFs untouched.
PETSc.dm_project_function!(petsclib, dm, 0.0, [zero_u_ptr], nothing,
                           LibPETSc.INSERT_VALUES, u)

# ── Solve ────────────────────────────────────────────────────────────────────
# SNESSetFromOptions + SNESSolve; b = nothing → PETSc assembles RHS internally.
PETSc.solve!(u, snes)

if MPI.Comm_rank(comm) == 0
    its = LibPETSc.SNESGetIterationNumber(petsclib, snes)
    println("SNES converged in $its Newton iterations.")
end

# ── L² error ─────────────────────────────────────────────────────────────────
l2err = PETSc.dm_compute_l2diff(petsclib, dm, 0.0, [quadratic_u_2d_ptr], nothing, u)
if MPI.Comm_rank(comm) == 0
    println("L2 error: $l2err")
end

# ── Cleanup ──────────────────────────────────────────────────────────────────
GC.gc(true)
MPI.Barrier(comm)

PETSc.destroy(snes)
PETSc.destroy(J)
PETSc.destroy(u)
# `fe` was handed off to dm via DMSetField — DMDestroy releases it.
PETSc.destroy(dm)
if !isinteractive()
    PETSc.finalize(petsclib)
    MPI.Barrier(comm)
    MPI.Finalize()
    ccall(:quick_exit, Cvoid, (Cint,), 0)
end
