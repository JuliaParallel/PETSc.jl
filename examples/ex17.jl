# EXCLUDE FROM TESTING
#=
  ex17.jl — Julia port of PETSc snes/tutorials/ex17.c
  Linear elasticity in 2D and 3D with finite elements (DMPLEX + PetscFE).

  Solves: −∇·σ(u) = f   in Ω = (0,1)^dim
  where σ = λ tr(ε) I + 2μ ε  (linear isotropic elasticity)
        ε = ½(∇u + ∇uᵀ)       (small-strain tensor)

  Material parameters (command-line overridable):
    -mu     N   shear modulus λ              (default 1.0)
    -lambda N   Lamé's first parameter μ    (default 1.0)
    -N      N   tension on right wall (axial_disp only) (default −1.0)

  Solution types (-sol_type):
    vlap_quad      — vector Laplacian, quadratic exact solution
    elas_quad      — linear elasticity, quadratic exact solution (default)
    vlap_trig      — vector Laplacian, trigonometric exact solution
    elas_trig      — linear elasticity, trig exact solution
    elas_axial_disp  — axial displacement; mixed Dirichlet/Neumann BCs
                       requires -dm_plex_separate_marker
    elas_uniform_strain — uniform strain; Dirichlet BCs with known exact solution
    elas_ge        — geological engineering; Dirichlet shift on right wall,
                     no analytical exact solution

  Near-null space for AMG:
    -near_nullspace false   disable rigid-body-mode near-null space (default: true)

  ── Basic usage ────────────────────────────────────────────────────────────────

  2D Q1 elasticity, quadratic exact solution:
    julia --project examples/ex17.jl \
          -dm_plex_simplex 0 -dm_plex_box_faces 4,4 \
          -sol_type elas_quad -displacement_petscspace_degree 1

  Check convergence (compares FEM solution to exact solution):
    julia --project examples/ex17.jl \
          -dm_plex_simplex 0 -sol_type elas_quad \
          -displacement_petscspace_degree 2 -dmsnes_check .0001

  2D P1 triangle mesh (requires the default triangle mesher in PETSc_jll):
    julia --project examples/ex17.jl \
          -sol_type vlap_quad -displacement_petscspace_degree 1 \
          -dm_refine 2

  3D Q1 hex mesh:
    julia --project examples/ex17.jl \
          -dm_plex_dim 3 -dm_plex_simplex 0 -dm_plex_box_faces 3,3,3 \
          -sol_type elas_quad -displacement_petscspace_degree 1

  ── Solver variants ────────────────────────────────────────────────────────────

  Direct solver (exact):
    julia --project examples/ex17.jl \
          -dm_plex_simplex 0 -dm_plex_box_faces 4,4 \
          -sol_type elas_quad -pc_type lu

  GAMG (algebraic multigrid) with rigid-body near-null space:
    julia --project examples/ex17.jl \
          -dm_plex_simplex 0 -dm_plex_box_faces 8,8 \
          -sol_type elas_quad -pc_type gamg -ksp_type cg

  Geometric multigrid (3 refinement levels):
    julia --project examples/ex17.jl \
          -dm_plex_dim 3 -dm_plex_simplex 0 -dm_plex_box_faces 2,2,2 \
          -dm_refine_hierarchy 3 -sol_type elas_quad \
          -pc_type mg -pc_mg_type full \
          -mg_levels_ksp_type chebyshev -mg_levels_pc_type jacobi \
          -ksp_type cg

  ── Specific solution types ─────────────────────────────────────────────────────

  Axial displacement (tension on right wall):
    julia --project examples/ex17.jl \
          -dm_plex_simplex 0 -dm_plex_box_faces 4,4 \
          -sol_type elas_axial_disp -dm_plex_separate_marker \
          -displacement_petscspace_degree 1 -pc_type lu

  Uniform strain:
    julia --project examples/ex17.jl \
          -dm_plex_simplex 0 -dm_plex_box_faces 4,4 \
          -sol_type elas_uniform_strain -pc_type lu

  Geological engineering (slab with imposed right-wall displacement):
    julia --project examples/ex17.jl \
          -dm_plex_simplex 0 -dm_plex_box_faces 4,4 \
          -sol_type elas_ge -dm_plex_separate_marker -pc_type gamg -ksp_type cg

  ── MPI ────────────────────────────────────────────────────────────────────────

  4 ranks, 3D elasticity with GAMG:
    mpiexec -n 4 julia --project examples/ex17.jl \
            -dm_plex_dim 3 -dm_plex_simplex 0 -dm_plex_box_faces 4,4,4 \
            -sol_type elas_quad -pc_type gamg -ksp_type cg

  ── Notes ──────────────────────────────────────────────────────────────────────

  - 3D simplex (tetrahedral) meshes require ctetgen/TetGen; PETSc_jll does not
    include it.  Use -dm_plex_simplex 0 (hex meshes) for 3D runs.
  - elas_axial_disp and elas_ge require -dm_plex_separate_marker so that each
    face of the box gets its own marker ID.  Face marker IDs:
      2D: bottom=1, right=2, top=3, left=4
      3D: bottom=1, right=2, top=3, left=4, front=5, back=6
  - elas_ge has no analytical exact solution; the printed L² error is meaningless.
  - VTK output (-vtk_output file.vtu) writes displacement as a 3-component vector
    and the Cauchy stress as a 9-component 3×3 tensor (row-major, zero-padded in 2D).
    PETSc.vtk_merge_tensor! post-processes the file so ParaView sees a proper tensor.

  References:
    - PETSc 3.23  src/snes/tutorials/ex17.c
    - https://en.wikipedia.org/wiki/Linear_elasticity
=#

using MPI
using PETSc
using PETSc: LibPETSc

# ── PETSc setup ─────────────────────────────────────────────────────────────────
_args    = isinteractive() ? String[] : collect(String, ARGS)
opts     = PETSc.parse_options(["-snes_monitor"; _args])
petsclib = PETSc.getlib(; PetscScalar = Float64)
MPI.Initialized() || MPI.Init()
PETSc.initialize(petsclib)

const PetscInt    = petsclib.PetscInt
const PetscScalar = Float64
const PetscReal   = Float64

comm = MPI.COMM_WORLD

# ── Option parsing ─────────────────────────────────────────────────────────────
_o = NamedTuple(pairs(opts))
dim        = parse(Int,    string(get(_o, :dm_plex_dim,    "2")))
sol_type   = string(get(_o, :sol_type, "elas_quad"))
_nns       = string(get(_o, :near_nullspace, "true"))
use_near_nullspace = _nns != "false" && _nns != "0"
mu_val     = parse(Float64, string(get(_o, :mu,     "1.0")))
lambda_val = parse(Float64, string(get(_o, :lambda, "1.0")))
N_val      = parse(Float64, string(get(_o, :N,      "-1.0")))

# Simplex flag: default triangles for 2D, hexahedra for 3D (TetGen not in jll).
_simp = get(_o, Symbol("dm_plex_simplex"), nothing)
simplex = _simp === nothing ? (dim < 3) : (parse(Int, string(_simp)) != 0)

VALID_SOL = ("vlap_quad", "elas_quad", "vlap_trig", "elas_trig",
             "elas_axial_disp", "elas_uniform_strain", "elas_ge")
@assert dim in (2, 3)        "-dm_plex_dim must be 2 or 3"
@assert sol_type in VALID_SOL  "-sol_type must be one of: $(join(VALID_SOL, ", "))"

# ── Exact solution / BC functions ──────────────────────────────────────────────
# Each wraps a PetscSimplePointFn: f(t, x::Vec[dim], u::Vec[Nc], ctx) → nothing
# Nc = dim (one component per spatial dimension for the displacement field).

function zero_u(t, x, u, ctx)
    fill!(u, 0)
end
const zero_u_ptr = PETSc.@petsc_simple_fn(zero_u)

function ge_shift_u(t, x, u, ctx)
    u[1] = 0.1
    for d in 2:length(u); u[d] = 0.0; end
end
const ge_shift_u_ptr = PETSc.@petsc_simple_fn(ge_shift_u)

function quadratic_2d_u(t, x, u, ctx)
    u[1] = x[1]^2
    u[2] = x[2]^2 - 2*x[1]*x[2]
end
const quadratic_2d_u_ptr = PETSc.@petsc_simple_fn(quadratic_2d_u)

function quadratic_3d_u(t, x, u, ctx)
    u[1] = x[1]^2
    u[2] = x[2]^2 - 2*x[1]*x[2]
    u[3] = x[3]^2 - 2*x[2]*x[3]
end
const quadratic_3d_u_ptr = PETSc.@petsc_simple_fn(quadratic_3d_u)

function trig_2d_u(t, x, u, ctx)
    u[1] = sin(2π * x[1])
    u[2] = sin(2π * x[2]) - 2*x[1]*x[2]
end
const trig_2d_u_ptr = PETSc.@petsc_simple_fn(trig_2d_u)

function trig_3d_u(t, x, u, ctx)
    u[1] = sin(2π * x[1])
    u[2] = sin(2π * x[2]) - 2*x[1]*x[2]
    u[3] = sin(2π * x[3]) - 2*x[2]*x[3]
end
const trig_3d_u_ptr = PETSc.@petsc_simple_fn(trig_3d_u)

function axial_disp_u(t, x, u, ctx)
    mu = mu_val; lambda = lambda_val; N = N_val
    u[1] = (3*lambda^2 + 8*lambda*mu + 4*mu^2) /
           (4*mu*(3*lambda^2 + 5*lambda*mu + 2*mu^2)) * N * x[1]
    u[2] = 0.25*lambda / (mu*(lambda+mu)) * N * x[2]
    for d in 3:length(u); u[d] = 0.0; end
end
const axial_disp_u_ptr = PETSc.@petsc_simple_fn(axial_disp_u)

function uniform_strain_u(t, x, u, ctx)
    eps_xx = 0.1; eps_xy = 0.3; eps_yy = 0.25
    u[1] = eps_xx * x[1] + eps_xy * x[2]
    u[2] = eps_xy * x[1] + eps_yy * x[2]
    for d in 3:length(u); u[d] = 0.0; end
end
const uniform_strain_u_ptr = PETSc.@petsc_simple_fn(uniform_strain_u)

# ── Residual functions f0 (zeroth-order) ───────────────────────────────────────
# f0 output size = Nc = dim (one entry per displacement component).
# constants[1]=mu, constants[2]=lambda, constants[3]=N  (1-indexed in Julia).

# Vector Laplacian, quadratic solution: body force = -Δu = −2 in each component.
function f0_vlap_quadratic(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                           aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, f0)
    for d in 1:dim_; f0[d] += 2.0; end
end
const f0_vlap_quadratic_ptr = PETSc.@petsc_residual_fn(f0_vlap_quadratic, dim_)

# Linear elasticity, quadratic solution: body force from −∇·σ.
function f0_elas_quadratic(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                           aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, f0)
    mu = real(constants[1]); lambda = real(constants[2])
    for d in 1:dim_-1; f0[d] += 2.0*mu; end
    f0[dim_] += 2.0*lambda + 4.0*mu
end
const f0_elas_quadratic_ptr = PETSc.@petsc_residual_fn(f0_elas_quadratic, dim_)

# Vector Laplacian, trig solution: body force = 4π² sin(2πx_d) per component.
function f0_vlap_trig(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                      aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, f0)
    for d in 1:dim_; f0[d] += -4.0*π^2 * sin(2π * x[d]); end
end
const f0_vlap_trig_ptr = PETSc.@petsc_residual_fn(f0_vlap_trig, dim_)

# Linear elasticity, trig solution: body force from −∇·σ.
function f0_elas_trig(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                      aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, f0)
    mu = real(constants[1]); lambda = real(constants[2])
    fact = 4.0*π^2
    for d in 1:dim_
        f0[d] += -(2.0*mu + lambda) * fact * sin(2π * x[d]) -
                 (d < dim_ ? 2.0*(mu + lambda) : 0.0)
    end
end
const f0_elas_trig_ptr = PETSc.@petsc_residual_fn(f0_elas_trig, dim_)

# ── Residual functions f1 (first-order, diffusion flux) ────────────────────────
# f1 output size = Nc * dim = dim * dim (gradient of Nc-component field).
# Indexing (0-based c,d): f1[c*dim+d+1] = flux in direction d of component c.

function f1_vlap(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                 aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, f1)
    for c in 0:dim_-1, d in 0:dim_-1
        f1[c*dim_ + d + 1] += u_x[c*dim_ + d + 1]
    end
end
const f1_vlap_ptr = PETSc.@petsc_residual_fn(f1_vlap, dim_*dim_)

function f1_elas(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                 aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, f1)
    mu = real(constants[1]); lambda = real(constants[2])
    for c in 0:dim_-1, d in 0:dim_-1
        f1[c*dim_ + d + 1] += mu * (u_x[c*dim_ + d + 1] + u_x[d*dim_ + c + 1])
        f1[c*dim_ + c + 1] += lambda * u_x[d*dim_ + d + 1]
    end
end
const f1_elas_ptr = PETSc.@petsc_residual_fn(f1_elas, dim_*dim_)

# ── Jacobian functions ─────────────────────────────────────────────────────────
# g3 output size = Nc*Nc*dim*dim = dim^4 (full stiffness tensor).
# Indexing (0-based c,gc,df,dg): g3[((c*Nc+gc)*dim+df)*dim+dg + 1].

function g3_vlap(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                 aOff, aOff_x, a, a_t, a_x, t, u_tShift, x, numConstants, constants, g3)
    fill!(g3, 0)
    for c in 0:dim_-1, d in 0:dim_-1
        g3[((c*dim_ + c)*dim_ + d)*dim_ + d + 1] = 1.0
    end
end
const g3_vlap_ptr = PETSc.@petsc_jacobian_fn(g3_vlap, dim_*dim_*dim_*dim_)

function g3_elas(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                 aOff, aOff_x, a, a_t, a_x, t, u_tShift, x, numConstants, constants, g3)
    mu = real(constants[1]); lambda = real(constants[2])
    fill!(g3, 0)
    for c in 0:dim_-1, d in 0:dim_-1
        g3[((c*dim_ + c)*dim_ + d)*dim_ + d + 1] += mu
        g3[((c*dim_ + d)*dim_ + d)*dim_ + c + 1] += mu
        g3[((c*dim_ + d)*dim_ + c)*dim_ + d + 1] += lambda
    end
end
const g3_elas_ptr = PETSc.@petsc_jacobian_fn(g3_elas, dim_*dim_*dim_*dim_)

# ── VTK post-processing: displacement (3-vector) and stress (3×3 tensor) ──────
#
# PETSc's VTK writer requires all fields to share the same DM.  We therefore
# create a combined output DM with two fields:
#   field 0 — displacement  (3 components, zero-padded in 2D)
#   field 1 — Cauchy stress (9 components, 3×3 row-major, zero-padded in 2D)
#
# Both are filled by DMProjectField from the solved displacement vector.
# Using 3 / 9 components gives ParaView proper Vector / Tensor data arrays.

# Copy displacement into a 3-component output (pad z=0 in 2D).
function copy_displacement_3(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                              aOff, aOff_x, a, a_t, a_x, t, x, nC, cst, out)
    for c in 1:dim_; out[c] = u[c]; end
    for c in dim_+1:3; out[c] = 0.0; end
end
const copy_displacement_3_ptr = PETSc.@petsc_residual_fn(copy_displacement_3, 3)

# Compute Cauchy stress σ_{cd} = μ(∂u_c/∂x_d + ∂u_d/∂x_c) + λ tr(ε) δ_{cd}
# and place it in a 3×3 output (9 components, row-major).
# Rows/columns beyond dim are zeroed (2D → z-row/column = 0).
function compute_stress_3x3(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                             aOff, aOff_x, a, a_t, a_x, t, x, nC, cst, out)
    mu = mu_val; lambda = lambda_val
    for i in 1:9; out[i] = 0.0; end   # zero-pad (covers z-row/col in 2D)
    div_u = zero(eltype(u_x))
    for d in 0:dim_-1; div_u += u_x[d*dim_+d+1]; end
    for c in 0:dim_-1, d in 0:dim_-1
        out[c*3+d+1] = mu*(u_x[c*dim_+d+1] + u_x[d*dim_+c+1])
    end
    for c in 0:dim_-1; out[c*3+c+1] += lambda * div_u; end
end
const compute_stress_3x3_ptr = PETSc.@petsc_residual_fn(compute_stress_3x3, 9)

# ── Neumann BC for elas_axial_disp: traction N on right wall ──────────────────
# f0_bd output size = Nc = dim.
function f0_elas_axial_disp_bd(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
                                aOff, aOff_x, a, a_t, a_x, t, x, n,
                                numConstants, constants, f0)
    f0[1] = real(constants[3])   # N — tension force, constants[3] = N (1-indexed)
end
const f0_elas_axial_disp_bd_ptr = PETSc.@petsc_bd_fn(f0_elas_axial_disp_bd, dim_)

# ── Near-null space constructor for AMG ────────────────────────────────────────
# Called by PETSc for each level when setting up the preconditioner.
# Registers rigid body modes (translations + rotations) as the near-null space.
function elasticity_nsp_constructor(
    dm_ptr::Ptr{Cvoid}, ::PetscInt, ::PetscInt,
    nsp_pp::Ptr{LibPETSc.MatNullSpace},
)::PetscInt
    PL   = typeof(petsclib)
    dm_w = LibPETSc.PetscDM{PL}(dm_ptr)
    nsp  = LibPETSc.DMPlexCreateRigidBody(LibPETSc.getlib(PL), dm_w, PetscInt(0))
    unsafe_store!(nsp_pp, nsp)
    return PetscInt(0)
end
const elasticity_nsp_ptr = Base.@cfunction(
    elasticity_nsp_constructor, PetscInt,
    (Ptr{Cvoid}, PetscInt, PetscInt, Ptr{LibPETSc.MatNullSpace}))

# ── Build mesh ─────────────────────────────────────────────────────────────────
# DMCreate + DMSetType(DMPLEX) + DMSetFromOptions — all mesh parameters
# (-dm_plex_dim, -dm_plex_simplex, -dm_plex_box_faces, -dm_refine_hierarchy, …)
# are passed through the options database.
dm = PETSc.DMPlex(petsclib, comm; opts...)

# ── Discretisation: displacement field with dim components ──────────────────────
# The FE prefix "displacement" means -displacement_petscspace_degree N controls
# the polynomial degree (default 1 = Q1/P1).
fe = PETSc.fe_create_default(petsclib, MPI.COMM_SELF, dim, dim, simplex;
                              prefix = "displacement")
PETSc.petsc_setname!(petsclib, fe, "displacement")
PETSc.setfield!(dm, 0, fe)
PETSc.createds!(dm)

# ── PetscDS: residual, Jacobian, exact solution, material constants ─────────────
ds = PETSc.getds(dm)

# Constants: [mu, lambda, N] — accessible in all pointwise functions.
PETSc.set_constants!(ds, [mu_val, lambda_val, N_val])

# Choose residual/Jacobian/exact-solution based on sol_type.
exact_ptr = if sol_type == "vlap_quad"
    PETSc.set_residual!(ds, 0, f0_vlap_quadratic_ptr, f1_vlap_ptr)
    PETSc.set_jacobian!(ds, 0, 0, C_NULL, C_NULL, C_NULL, g3_vlap_ptr)
    dim == 2 ? quadratic_2d_u_ptr : quadratic_3d_u_ptr

elseif sol_type == "elas_quad"
    PETSc.set_residual!(ds, 0, f0_elas_quadratic_ptr, f1_elas_ptr)
    PETSc.set_jacobian!(ds, 0, 0, C_NULL, C_NULL, C_NULL, g3_elas_ptr)
    dim == 2 ? quadratic_2d_u_ptr : quadratic_3d_u_ptr

elseif sol_type == "vlap_trig"
    PETSc.set_residual!(ds, 0, f0_vlap_trig_ptr, f1_vlap_ptr)
    PETSc.set_jacobian!(ds, 0, 0, C_NULL, C_NULL, C_NULL, g3_vlap_ptr)
    dim == 2 ? trig_2d_u_ptr : trig_3d_u_ptr

elseif sol_type == "elas_trig"
    PETSc.set_residual!(ds, 0, f0_elas_trig_ptr, f1_elas_ptr)
    PETSc.set_jacobian!(ds, 0, 0, C_NULL, C_NULL, C_NULL, g3_elas_ptr)
    dim == 2 ? trig_2d_u_ptr : trig_3d_u_ptr

elseif sol_type == "elas_axial_disp"
    # Zero body force; traction applied on right wall via Neumann BC below.
    PETSc.set_residual!(ds, 0, C_NULL, f1_elas_ptr)
    PETSc.set_jacobian!(ds, 0, 0, C_NULL, C_NULL, C_NULL, g3_elas_ptr)
    axial_disp_u_ptr

elseif sol_type == "elas_uniform_strain"
    PETSc.set_residual!(ds, 0, C_NULL, f1_elas_ptr)
    PETSc.set_jacobian!(ds, 0, 0, C_NULL, C_NULL, C_NULL, g3_elas_ptr)
    uniform_strain_u_ptr

else  # elas_ge
    PETSc.set_residual!(ds, 0, C_NULL, f1_elas_ptr)
    PETSc.set_jacobian!(ds, 0, 0, C_NULL, C_NULL, C_NULL, g3_elas_ptr)
    zero_u_ptr  # no exact solution
end

PETSc.set_exact_solution!(ds, 0, exact_ptr)

# ── Boundary conditions ─────────────────────────────────────────────────────────
label = PETSc.getlabel(dm, "marker")

# Face marker IDs with -dm_plex_separate_marker:
#   2D: bottom=1, right=2, top=3, left=4
#   3D: bottom=1, right=2, top=3, left=4, front=5, back=6
right_id = dim == 3 ? PetscInt(5) : PetscInt(2)
left_id  = dim == 3 ? PetscInt(6) : PetscInt(4)

if sol_type in ("vlap_quad", "elas_quad", "vlap_trig", "elas_trig")
    # Dirichlet on entire boundary (marker 1 covers all faces without separate markers).
    PETSc.add_boundary!(petsclib, dm, LibPETSc.DM_BC_ESSENTIAL, "wall", label,
                        PetscInt[1], 0, PetscInt[], exact_ptr)

elseif sol_type == "elas_axial_disp"
    # Neumann traction on right wall; zero Dirichlet on left/bottom (some components).
    # Requires -dm_plex_separate_marker.
    PETSc.add_natural_boundary!(petsclib, dm, ds, "right", label, right_id, 0,
                                f0_elas_axial_disp_bd_ptr)
    PETSc.add_boundary!(petsclib, dm, LibPETSc.DM_BC_ESSENTIAL, "left", label,
                        PetscInt[left_id], 0, PetscInt[0], zero_u_ptr)
    bottom_cmp = dim == 3 ? PetscInt[2] : PetscInt[1]
    PETSc.add_boundary!(petsclib, dm, LibPETSc.DM_BC_ESSENTIAL, "bottom", label,
                        PetscInt[1], 0, bottom_cmp, zero_u_ptr)
    if dim == 3
        PETSc.add_boundary!(petsclib, dm, LibPETSc.DM_BC_ESSENTIAL, "front", label,
                            PetscInt[3], 0, PetscInt[1], zero_u_ptr)
    end

elseif sol_type == "elas_uniform_strain"
    PETSc.add_boundary!(petsclib, dm, LibPETSc.DM_BC_ESSENTIAL, "wall", label,
                        PetscInt[1], 0, PetscInt[], exact_ptr)

else  # elas_ge — requires -dm_plex_separate_marker
    PETSc.add_boundary!(petsclib, dm, LibPETSc.DM_BC_ESSENTIAL, "left", label,
                        PetscInt[left_id], 0, PetscInt[], zero_u_ptr)
    PETSc.add_boundary!(petsclib, dm, LibPETSc.DM_BC_ESSENTIAL, "right", label,
                        PetscInt[right_id], 0, PetscInt[0], ge_shift_u_ptr)
end

# ── Propagate disc + near-null space to coarser DMs (for GMG / GAMG hierarchy) ──
let cdm = PETSc.dm_get_coarse(dm)
    while convert(Ptr{Cvoid}, cdm) != C_NULL
        PETSc.dm_copy_disc!(dm, cdm)
        if use_near_nullspace
            LibPETSc.DMSetNearNullSpaceConstructor(
                petsclib, cdm, PetscInt(0), elasticity_nsp_ptr)
        end
        cdm = PETSc.dm_get_coarse(cdm)
    end
end
if use_near_nullspace
    LibPETSc.DMSetNearNullSpaceConstructor(
        petsclib, dm, PetscInt(0), elasticity_nsp_ptr)
end

# ── SNES + matrix ───────────────────────────────────────────────────────────────
snes = PETSc.SNES(petsclib, comm; opts...)
PETSc.setDM!(snes, dm)
u = PETSc.DMGlobalVec(dm)
J = PETSc.MatAIJ(dm)
PETSc.plex_set_snes_local_fem!(petsclib, dm)
LibPETSc.SNESSetJacobian(petsclib, snes, J, J, C_NULL, C_NULL)

# ── Initial guess ───────────────────────────────────────────────────────────────
# Seed constrained DOFs with the exact solution (or zero for elas_ge).
PETSc.dm_project_function!(petsclib, dm, 0.0, [exact_ptr], nothing,
                            LibPETSc.INSERT_ALL_VALUES, u)
# Reset free DOFs to zero so SNES has a well-posed starting point.
PETSc.dm_project_function!(petsclib, dm, 0.0, [zero_u_ptr], nothing,
                            LibPETSc.INSERT_VALUES, u)

# ── Solve ───────────────────────────────────────────────────────────────────────
PETSc.solve!(u, snes)

if MPI.Comm_rank(comm) == 0
    its = LibPETSc.SNESGetIterationNumber(petsclib, snes)
    println("SNES converged in $its Newton iteration(s).")
end

# ── L² error ────────────────────────────────────────────────────────────────────
l2err = PETSc.dm_compute_l2diff(petsclib, dm, 0.0, [exact_ptr], nothing, u)
if MPI.Comm_rank(comm) == 0
    if sol_type == "elas_ge"
        println("(elas_ge has no exact solution; L² comparison is against zero.)")
    end
    println("L2 error: $l2err")
end

# ── VTK output ──────────────────────────────────────────────────────────────────
let _vtk = get(NamedTuple(pairs(opts)), :vtk_output, nothing)
    if _vtk !== nothing
        fname = string(_vtk)

        # ── Combined output DM: displacement (field 0) + stress (field 1) ──────
        # PETSc's VTK viewer requires all fields in one DM.  We clone the mesh and
        # attach both fields, then write a single combined vector.
        #
        # VTK arrays produced (row-major 3×3, zero-padded in 2D):
        #   "displacement"  NumberOfComponents=3  → Vector in ParaView
        #   "stress"        NumberOfComponents=9  → 3×3 Tensor in ParaView
        dm_out = PETSc.dmclone(dm)

        fe_disp = PETSc.fe_create_default(petsclib, MPI.COMM_SELF,
                      dim, dim, simplex; prefix = "displacement")
        PETSc.petsc_setname!(petsclib, fe_disp, "displacement")

        fe_stress = PETSc.fe_create_default(petsclib, MPI.COMM_SELF,
                        dim, 9, simplex; prefix = "stress")
        PETSc.petsc_setname!(petsclib, fe_stress, "stress")

        PETSc.setfield!(dm_out, 0, fe_disp)
        PETSc.setfield!(dm_out, 1, fe_stress)
        PETSc.createds!(dm_out)

        out_vec = PETSc.dm_create_global_vec(dm_out)
        # Empty vector name: VTK array names come purely from the FE names above.
        PETSc.petsc_setname!(petsclib, out_vec, "")

        # DMProjectField projects functions of the displacement solution (u, 1 field)
        # into both output fields simultaneously.
        PETSc.dm_project_field!(petsclib, dm_out, 0.0, u,
                                [copy_displacement_3_ptr, compute_stress_3x3_ptr],
                                LibPETSc.INSERT_ALL_VALUES, out_vec)

        PETSc.vtk_save!(petsclib, comm, fname, out_vec)

        MPI.Comm_rank(comm) == 0 && PETSc.vtk_merge_tensor!(fname, "stress")

        MPI.Comm_rank(comm) == 0 &&
            println("Displacement (vector) and stress (tensor) written to $fname")
    end
end

# ── Cleanup ─────────────────────────────────────────────────────────────────────
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
