# EXCLUDE FROM TESTING
#=
  ex19.jl — 2D Driven Cavity: velocity–vorticity–temperature formulation
  Port of PETSc snes/tutorials/ex19.c

  Solves on the unit square [0,1]²:
    −∇²u  − ∂ω/∂y               = 0
    −∇²v  + ∂ω/∂x               = 0
    −∇²ω  + ∇·(uω, vω) − Gr·∂T/∂x = 0
    −∇²T  + Pr·∇·(uT, vT)       = 0

  4 DOFs per node: (u, v, ω, T).  Convective terms are upwinded.

  Boundary conditions:
    All walls:   u = v = 0  (no-slip),  except top lid: u = lidvelocity
    Left:        T = 0  (cold, Dirichlet)
    Right:       T = 1  (hot, Dirichlet, only when Gr > 0)
    Top/bottom:  ∂T/∂n = 0  (insulated, Neumann)
    ω:           derived from the no-slip condition at each wall

  Usage (from the examples/ directory):
    # Basic run (4×4 default grid)
    julia --project ex19.jl

    # Larger grid with SNES convergence output
    julia --project ex19.jl -snes_monitor -snes_converged_reason -da_grid_x 129 -da_grid_y 129

    # With PETSc performance log
    julia --project ex19.jl -da_grid_x 129 -da_grid_y 129 -log_view

    # Multigrid preconditioner (3 levels, Chebyshev+Jacobi smoothers)
    julia --project ex19.jl -da_grid_x 125 -da_grid_y 125 \\
        -pc_type mg -pc_mg_levels 3 \\
        -mg_levels_ksp_type chebyshev -mg_levels_pc_type jacobi \\
        -snes_monitor -ksp_monitor

    # MPI parallel (4 ranks)
    mpiexec -n 4 julia --project ex19.jl \\
        -da_grid_x 256 -da_grid_y 256 \\
        -pc_type mg -pc_mg_levels 4 \\
        -mg_levels_ksp_type chebyshev -mg_levels_pc_type jacobi

    # GPU (CUDA) — set  useCUDA = true  at the top of the file, then:
    julia --project ex19.jl -da_grid_x 256 -da_grid_y 256 \\
        -pc_type mg -pc_mg_levels 4 \\
        -mg_levels_ksp_type chebyshev -mg_levels_pc_type jacobi

  Jacobian strategy:
    Fine-grid level: manual FD coloring via PETSc.dmda_star_fd_coloring +
      MatSetPreallocationCOOLocal + MatSetValuesCOO.  On GPU the entire
      perturb → F(x+h) → accumulate loop runs on-device with no host copies.
    Coarser MG levels: fall back to SNESComputeJacobianDefaultColor
      (correct for each level's DM, CPU-only).

  Requires: LocalPreferences.toml in examples/ with PetscInt = "Int32" matching
    the PETSc build (check with: grep sizeof_PetscInt petscconf.h).

  GPU usage: set  useCUDA = true  at the top of this file.
    Requires PETSc built with --with-cuda, and CUDA.jl in the environment.
=#

# ── GPU switch ────────────────────────────────────────────────────────────────
const useCUDA = true

using MPI
using PETSc
using KernelAbstractions

if useCUDA
    using CUDA
    import CUDA: CuArray, CuPtr, unsafe_wrap
    const backend = CUDABackend()
else
    const backend = KernelAbstractions.CPU()
end


# ── Physical parameters ──────────────────────────────────────────────────────
Base.@kwdef mutable struct AppCtx{T}
    lidvelocity :: T = one(T)
    prandtl     :: T = one(T)
    grashof     :: T = one(T)
end

# ── Residual kernel ───────────────────────────────────────────────────────────
#
# Iterates over the locally-owned grid: (li, lj) ∈ 1..nx_own × 1..ny_own.
#
# Array layout (plain 1-based, no OffsetArray — GPU-compatible):
#   x_par[dof, xi, xj]   ghost input,  xi = li + ox, xj = lj + oy
#   f_par[dof, li, lj]   owned output
#
# ox = xs - xsg, oy = ys - ysg  (ghost offsets; 0 at a domain boundary, 1 otherwise)
#
# Corner priority matches the C code: left/right walls processed last, so they
# take precedence over bottom/top at corners.  Achieved here by checking i==1
# and i==mx before j==1 and j==my.
#
@kernel function cavity_residual_kernel!(
    f_par,
    x_par,
    dhx, dhy,           # mx-1, my-1
    hx,  hy,            # 1/dhx, 1/dhy
    hydhx, hxdhy,       # metric factors for Laplacian scaling
    grashof, prandtl, lid,
    mx :: Int, my :: Int,
    xs :: Int, ys :: Int,
    ox :: Int, oy :: Int,
)
    li, lj = @index(Global, NTuple)

    i  = xs + li - 1   # global 1-based x-coordinate
    j  = ys + lj - 1   # global 1-based y-coordinate
    xi = li + ox        # ghost-array x-index for this point
    xj = lj + oy        # ghost-array y-index for this point

    # ── Boundary conditions ───────────────────────────────────────────────────

    if i == 1                           # left wall — cold (T = 0)
        f_par[1, li, lj] = x_par[1, xi,   xj  ]
        f_par[2, li, lj] = x_par[2, xi,   xj  ]
        f_par[3, li, lj] = x_par[3, xi,   xj  ] -
            (x_par[2, xi+1, xj] - x_par[2, xi, xj]) * dhx
        f_par[4, li, lj] = x_par[4, xi,   xj  ]

    elseif i == mx                      # right wall — hot (T = 1 when Gr > 0)
        f_par[1, li, lj] = x_par[1, xi,   xj  ]
        f_par[2, li, lj] = x_par[2, xi,   xj  ]
        f_par[3, li, lj] = x_par[3, xi,   xj  ] -
            (x_par[2, xi, xj] - x_par[2, xi-1, xj]) * dhx
        f_par[4, li, lj] = x_par[4, xi,   xj  ] -
            (grashof > zero(grashof) ? one(grashof) : zero(grashof))

    elseif j == 1                       # bottom wall — no-slip, insulated
        f_par[1, li, lj] = x_par[1, xi,   xj  ]
        f_par[2, li, lj] = x_par[2, xi,   xj  ]
        f_par[3, li, lj] = x_par[3, xi,   xj  ] +
            (x_par[1, xi, xj+1] - x_par[1, xi, xj]) * dhy
        f_par[4, li, lj] = x_par[4, xi,   xj  ] - x_par[4, xi, xj+1]

    elseif j == my                      # top wall — moving lid, insulated
        f_par[1, li, lj] = x_par[1, xi,   xj  ] - lid
        f_par[2, li, lj] = x_par[2, xi,   xj  ]
        f_par[3, li, lj] = x_par[3, xi,   xj  ] +
            (x_par[1, xi, xj] - x_par[1, xi, xj-1]) * dhy
        f_par[4, li, lj] = x_par[4, xi,   xj  ] - x_par[4, xi, xj-1]

    else                                # ── interior point ───────────────────

        # Upwind split of advecting velocities
        vx  = x_par[1, xi, xj];  avx = abs(vx)
        vxp = oftype(vx, 0.5) * (vx + avx)   # max(vx, 0)
        vxm = oftype(vx, 0.5) * (vx - avx)   # min(vx, 0)

        vy  = x_par[2, xi, xj];  avy = abs(vy)
        vyp = oftype(vy, 0.5) * (vy + avy)
        vym = oftype(vy, 0.5) * (vy - avy)

        # u-equation:  −∇²u − ∂ω/∂y = 0
        cu  = x_par[1, xi, xj]
        uxx = (2cu - x_par[1, xi-1, xj] - x_par[1, xi+1, xj]) * hydhx
        uyy = (2cu - x_par[1, xi, xj-1] - x_par[1, xi, xj+1]) * hxdhy
        f_par[1, li, lj] = uxx + uyy -
            oftype(cu, 0.5) * (x_par[3, xi, xj+1] - x_par[3, xi, xj-1]) * hx

        # v-equation:  −∇²v + ∂ω/∂x = 0
        cv  = x_par[2, xi, xj]
        vxx = (2cv - x_par[2, xi-1, xj] - x_par[2, xi+1, xj]) * hydhx
        vyy = (2cv - x_par[2, xi, xj-1] - x_par[2, xi, xj+1]) * hxdhy
        f_par[2, li, lj] = vxx + vyy +
            oftype(cv, 0.5) * (x_par[3, xi+1, xj] - x_par[3, xi-1, xj]) * hy

        # ω-equation:  −∇²ω + ∇·(uω, vω) − Gr·∂T/∂x = 0
        cω  = x_par[3, xi, xj]
        wxx = (2cω - x_par[3, xi-1, xj] - x_par[3, xi+1, xj]) * hydhx
        wyy = (2cω - x_par[3, xi, xj-1] - x_par[3, xi, xj+1]) * hxdhy
        f_par[3, li, lj] = wxx + wyy +
            (vxp * (cω - x_par[3, xi-1, xj]) + vxm * (x_par[3, xi+1, xj] - cω)) * hy +
            (vyp * (cω - x_par[3, xi, xj-1]) + vym * (x_par[3, xi, xj+1] - cω)) * hx -
            oftype(cω, 0.5) * grashof * (x_par[4, xi+1, xj] - x_par[4, xi-1, xj]) * hy

        # T-equation:  −∇²T + Pr·∇·(uT, vT) = 0
        cT  = x_par[4, xi, xj]
        txx = (2cT - x_par[4, xi-1, xj] - x_par[4, xi+1, xj]) * hydhx
        tyy = (2cT - x_par[4, xi, xj-1] - x_par[4, xi, xj+1]) * hxdhy
        f_par[4, li, lj] = txx + tyy + prandtl * (
            (vxp * (cT - x_par[4, xi-1, xj]) + vxm * (x_par[4, xi+1, xj] - cT)) * hy +
            (vyp * (cT - x_par[4, xi, xj-1]) + vym * (x_par[4, xi, xj+1] - cT)) * hx)
    end
end

# ── FD-coloring GPU helper kernels ────────────────────────────────────────────
#
# scatter_perturb_kernel!: add h to selected entries of a vector.
#   cols[k] is the 1-based Julia index to perturb.
#
@kernel function scatter_perturb_kernel!(x, cols, h)
    k = @index(Global)
    @inbounds x[cols[k]] += h
end

# fd_accumulate_kernel!: write (f1 − f0)/h into val at COO indices.
#   coo_idxs[k] and row_idxs[k] are 1-based Julia indices.
#
@kernel function fd_accumulate_kernel!(val, f0, f1, coo_idxs, row_idxs, inv_h)
    k = @index(Global)
    @inbounds val[coo_idxs[k]] = (f1[row_idxs[k]] - f0[row_idxs[k]]) * inv_h
end

# ── FD-coloring Jacobian fill ─────────────────────────────────────────────────
#
# Fills val_dev[k] with the forward-difference Jacobian value for COO slot k.
# Loops over colors: scatter +h onto the owned columns of that color, evaluate
# F(x + h·eₖ), then accumulate (F1 − F0)/h into val_dev at the matching slots.
# Does NOT call MatSetValuesCOO or MatAssembly; those remain with the caller.
#
# Captures from module scope: useCUDA, backend, CuArray, CuPtr,
#   scatter_perturb_kernel!, fd_accumulate_kernel!, KernelAbstractions.
#
function maybe_wrap_device(arr, mtype, n, ::Type{T}, useCUDA) where T
    if useCUDA && mtype == LibPETSc.PETSC_MEMTYPE_DEVICE
        return unsafe_wrap(CuArray, CuPtr{T}(UInt64(pointer(arr))), n)
    else
        return arr
    end
end

function fd_coloring_jac!(
    petsclib,
    snes,
    g_x,
    f0_vec, f1_vec, x_pert_vec,
    val_dev  :: AbstractVector{T},
    n_colors    :: Int,
    n_local_dofs :: Int,
    perturb_cols_dev,
    coo_idxs_dev,
    local_rows_dev,
    h_eps :: T,
    inv_h :: T,
) where T
    LibPETSc.SNESComputeFunction(petsclib, snes, g_x, f0_vec)
    f0_arr, f0_mtype = LibPETSc.VecGetArrayReadAndMemType(petsclib, f0_vec)
    f0_dev = maybe_wrap_device(f0_arr, f0_mtype, n_local_dofs, T, useCUDA)

    for c in 1:n_colors
        isempty(perturb_cols_dev[c]) && continue

        LibPETSc.VecCopy(petsclib, g_x, x_pert_vec)
        xp_arr, xp_mtype = LibPETSc.VecGetArrayAndMemType(petsclib, x_pert_vec)
        xp_dev = maybe_wrap_device(xp_arr, xp_mtype, n_local_dofs, T, useCUDA)
        scatter_perturb_kernel!(backend, 64)(
            xp_dev, perturb_cols_dev[c], h_eps;
            ndrange = length(perturb_cols_dev[c]))
        KernelAbstractions.synchronize(backend)
        LibPETSc.VecRestoreArrayAndMemType(petsclib, x_pert_vec, xp_arr)

        LibPETSc.SNESComputeFunction(petsclib, snes, x_pert_vec, f1_vec)

        f1_arr, f1_mtype = LibPETSc.VecGetArrayReadAndMemType(petsclib, f1_vec)
        f1_dev = maybe_wrap_device(f1_arr, f1_mtype, n_local_dofs, T, useCUDA)
        fd_accumulate_kernel!(backend, 64)(
            val_dev, f0_dev, f1_dev,
            coo_idxs_dev[c], local_rows_dev[c], inv_h;
            ndrange = length(coo_idxs_dev[c]))
        KernelAbstractions.synchronize(backend)
        LibPETSc.VecRestoreArrayReadAndMemType(petsclib, f1_vec, f1_arr)
    end

    LibPETSc.VecRestoreArrayReadAndMemType(petsclib, f0_vec, f0_arr)
    return nothing
end

# ── Setup ─────────────────────────────────────────────────────────────────────
opts     = isinteractive() ? NamedTuple() : PETSc.parse_options(filter(a -> a != "-log_view", ARGS))
log_view = "-log_view" in ARGS

petsclib = PETSc.getlib(; PetscScalar = Float64, PetscInt = Int32)
PETSc.initialize(petsclib; log_view)

_T        = Float64
PetscInt = petsclib.PetscInt
comm     = MPI.COMM_WORLD

# DMDA: 4×4 default (matches ex19.c); override via -da_grid_x / -da_grid_y
da = PETSc.DMDA(
    petsclib, comm,
    (PETSc.DM_BOUNDARY_NONE, PETSc.DM_BOUNDARY_NONE),
    (4, 4),
    4,   # DOFs per node: (u, v, ω, T)
    1,   # stencil width
    PETSc.DMDA_STENCIL_STAR;
    opts...,
)

# GPU vecs and GPU matrix enable a fully GPU-resident FD coloring path via
# COO preallocation.  No host↔device bouncing in residual or Jacobian.
if useCUDA
    LibPETSc.DMSetVecType(petsclib, da, "cuda")
    LibPETSc.DMSetMatType(petsclib, da, "aijcusparse")
end

snes = PETSc.SNES(petsclib, comm; opts...)
PETSc.setDM!(snes, da)

# Actual grid size after setfromoptions (may differ from the 4×4 default)
info = PETSc.getinfo(da)
mx   = Int(info.global_size[1])
my   = Int(info.global_size[2])

user = AppCtx{_T}(
    lidvelocity = _T(1) / (mx - 1),
    prandtl     = _T(1),
    grashof     = _T(1),
)

# Precomputed grid metrics
dhx   = _T(mx - 1);   dhy   = _T(my - 1)
hx    = one(_T) / dhx; hy    = one(_T) / dhy
hydhx = hy * dhx;    hxdhy = hx * dhy

# ── Initial condition: u = v = ω = 0, T linear in x ─────────────────────────
x = PETSc.DMGlobalVec(da)

PETSc.withlocalarray!(x; read = false) do x_arr
    corners = PETSc.getcorners(da)
    xs = corners.lower[1];  ys = corners.lower[2]
    xe = corners.upper[1];  ye = corners.upper[2]
    nx_own = xe - xs + 1;   ny_own = ye - ys + 1
    dx = one(_T) / (mx - 1)
    x_par = reshape(x_arr, 4, nx_own, ny_own)
    for lj in 1:ny_own, li in 1:nx_own
        ig = xs + li - 1
        x_par[1, li, lj] = zero(_T)
        x_par[2, li, lj] = zero(_T)
        x_par[3, li, lj] = zero(_T)
        x_par[4, li, lj] = user.grashof > 0 ? _T(ig - 1) * dx : zero(_T)
    end
end

# ── Residual callback ─────────────────────────────────────────────────────────
r = similar(x)

PETSc.setfunction!(snes, r) do g_fx, snes, g_x
    da = PETSc.getDM(snes)

    l_x = PETSc.DMLocalVec(da)
    PETSc.dm_global_to_local!(g_x, l_x, da, PETSc.INSERT_VALUES)

    # Get arrays for the output (g_fx) and ghost-padded input (l_x) Vecs.
    # On GPU, returns CuArray wrappers (zero-copy when both Vecs are device-
    # resident) together with the raw PETSc handles needed for the restore call.
    # On CPU, returns plain Array views backed by VecGetArray.
    # fx_bounce is a GPU scratch buffer used when g_fx is host-resident; it is
    # copied back D2H by restore_petsc_arrays after the kernel completes.
    fx, lx, fx_arr, lx_arr, fx_bounce = PETSc.get_petsc_arrays(petsclib, g_fx, l_x)

    corners       = PETSc.getcorners(da)
    ghost_corners = PETSc.getghostcorners(da)

    xs  = corners.lower[1];        ys  = corners.lower[2]
    xe  = corners.upper[1];        ye  = corners.upper[2]
    xsg = ghost_corners.lower[1];  ysg = ghost_corners.lower[2]
    xeg = ghost_corners.upper[1];  yeg = ghost_corners.upper[2]

    nx_own = xe  - xs  + 1;  ny_own = ye  - ys  + 1
    nx_g   = xeg - xsg + 1;  ny_g   = yeg - ysg + 1

    # Plain [dof, x, y] arrays — no OffsetArray, safe for KA on GPU
    x_par = reshape(lx, 4, nx_g,   ny_g)
    f_par = reshape(fx, 4, nx_own, ny_own)

    # Ghost offset: ghost-array index for owned start = 1 + ox (0 at domain wall)
    ox = xs - xsg
    oy = ys - ysg

    # Recompute grid metrics from the DM so this callback is correct on every
    # MG level (coarsen/refine changes mx/my; capturing outer-scope values
    # would give wrong stencil weights and lid velocity on coarse grids).
    info_  = PETSc.getinfo(da)
    mx_    = Int(info_.global_size[1])
    my_    = Int(info_.global_size[2])
    dhx_   = _T(mx_ - 1);   dhy_   = _T(my_ - 1)
    hx_    = one(_T) / dhx_; hy_   = one(_T) / dhy_
    hydhx_ = hy_ * dhx_;    hxdhy_ = hx_ * dhy_
    lid_   = _T(1) / dhx_    # lidvelocity = 1/(mx-1)

    cavity_residual_kernel!(backend, 64)(
        f_par, x_par,
        dhx_, dhy_, hx_, hy_, hydhx_, hxdhy_,
        user.grashof, user.prandtl, lid_,
        mx_, my_, xs, ys, ox, oy;
        ndrange = (nx_own, ny_own),
    )
    KernelAbstractions.synchronize(backend)

    PETSc.restore_petsc_arrays(petsclib, g_fx, l_x, fx, lx, fx_arr, lx_arr, fx_bounce)
    PETSc.destroy(l_x)
    return PetscInt(0)
end

# ── Jacobian: manual FD coloring with GPU-efficient COO matrix assembly ───────
#
# Setup (one-time):
#   dmda_star_fd_coloring  — builds IS_COLORING_LOCAL coloring and analytically
#     derives the STAR-stencil COO (row, col) triplets; returns per-color index
#     arrays for the owned columns and COO slots to fill.
#   MatSetPreallocationCOOLocal  — registers the COO pattern (enables on-device
#     scatter on GPU; avoids per-entry hash-table lookups on CPU).
#
# Each Newton step (fd_coloring_jac!):
#   For each color: copy x → x_pert, scatter +h to owned cols of that color
#   (GPU kernel), evaluate F(x_pert), accumulate (F1−F0)/h into val[] (GPU
#   kernel).  Then MatSetValuesCOO assembles J from val[] in one call.
#
# NOTE: For MG, coarser levels fall back to SNESComputeJacobianDefaultColor
#       (correct, but CPU-only FD coloring for those levels).

# ── Coloring + COO index setup ────────────────────────────────────────────────
# Builds the IS_COLORING_LOCAL coloring for da's STAR stencil, ghost-local COO
# (row, col) pairs, and per-color owned-column / COO-entry index arrays.
# 2-D DMDA STAR stencil only; see PETSc.dmda_star_fd_coloring for 3-D notes.
coloring      = PETSc.dmda_star_fd_coloring(petsclib, da)
n_colors      = coloring.n_colors
n_local_dofs  = coloring.n_local_dofs
nnz_coo       = coloring.nnz_coo
row_coo_local = coloring.row_coo_local
col_coo_local = coloring.col_coo_local

# ── Create J ─────────────────────────────────────────────────────────────────
J = LibPETSc.DMCreateMatrix(petsclib, da)
# Register the COO pattern on both CPU and GPU.  This allows MatSetValuesCOO
# to be used for assembly in both cases, avoiding per-entry hash-table lookups
# that MatSetValuesLocal incurs.  On GPU it also enables device-side scatter.
LibPETSc.MatSetPreallocationCOOLocal(petsclib, J, LibPETSc.PetscCount(nnz_coo), row_coo_local, col_coo_local)

# ── Per-color index arrays for the FD loop ───────────────────────────────────
perturb_cols_1b = coloring.perturb_cols
coo_idxs_1b     = coloring.coo_idxs
local_rows_1b   = coloring.local_rows

if useCUDA
    perturb_cols_dev = [CuArray(v) for v in perturb_cols_1b]
    coo_idxs_dev     = [CuArray(v) for v in coo_idxs_1b]
    local_rows_dev   = [CuArray(v) for v in local_rows_1b]
    val_dev          = CUDA.zeros(_T, nnz_coo)
else
    perturb_cols_dev = perturb_cols_1b
    coo_idxs_dev     = coo_idxs_1b
    local_rows_dev   = local_rows_1b
    val_dev          = zeros(_T, nnz_coo)
end

# ── Scratch vectors for the FD loop ────────────────────────────────────────
x_pert_vec = LibPETSc.VecDuplicate(petsclib, x)
f0_vec     = LibPETSc.VecDuplicate(petsclib, x)
f1_vec     = LibPETSc.VecDuplicate(petsclib, x)
h_eps      = _T(sqrt(eps(_T)))
inv_h      = _T(1) / h_eps

# ── Jacobian callback ────────────────────────────────────────────────────────
PETSc.setjacobian!(snes, J) do Jmat, actual_snes, g_x
    # For MG: if this is a coarser level (different DM than the fine-grid da),
    # fall back to PETSc's built-in FD coloring (correct for that level's DM).
    da_level = PETSc.getDM(actual_snes)
    if da_level.ptr != da.ptr
        LibPETSc.SNESComputeJacobianDefaultColor(petsclib, actual_snes, g_x, Jmat, Jmat, C_NULL)
        return PetscInt(0)
    end

    # ── FD coloring: fill val_dev with Jacobian entries ────────────────────────
    fd_coloring_jac!(
        petsclib, actual_snes, g_x,
        f0_vec, f1_vec, x_pert_vec, val_dev,
        n_colors, n_local_dofs,
        perturb_cols_dev, coo_idxs_dev, local_rows_dev,
        h_eps, inv_h,
    )

    # ── Assemble J via COO ─────────────────────────────────────────────────────
    # On GPU: pass a raw device pointer so cuSPARSE scatters directly on device.
    # On CPU: pass the Vector{T} directly (uses the vector overload).
    # Both paths use the COO pattern registered by MatSetPreallocationCOOLocal.
    if useCUDA
        LibPETSc.MatSetValuesCOO(petsclib, Jmat, Ptr{_T}(UInt64(pointer(val_dev))), LibPETSc.INSERT_VALUES)
    else
        LibPETSc.MatSetValuesCOO(petsclib, Jmat, val_dev, LibPETSc.INSERT_VALUES)
    end
    LibPETSc.MatAssemblyBegin(petsclib, Jmat, LibPETSc.MAT_FINAL_ASSEMBLY)
    LibPETSc.MatAssemblyEnd(petsclib, Jmat, LibPETSc.MAT_FINAL_ASSEMBLY)
    if useCUDA
        # Force GPU→CPU sync so both copies are valid.
        # MatBindToCPU(PETSC_TRUE) triggers MatSeqAIJCUSPARSECopyFromGPU when
        # offloadmask==PETSC_OFFLOAD_GPU, making the CPU CSR correct.
        # MatBindToCPU(PETSC_FALSE) then releases the CPU-only restriction,
        # leaving offloadmask==PETSC_OFFLOAD_BOTH so that:
        #   MatGetDiagonal (Jacobi smoother in MG)  → reads CPU copy ✓
        #   MatPtAP (Galerkin coarse-op formation)  → uses GPU copy  ✓
        LibPETSc.MatBindToCPU(petsclib, Jmat, LibPETSc.PETSC_TRUE)
        LibPETSc.MatBindToCPU(petsclib, Jmat, LibPETSc.PETSC_FALSE)
    end
    return PetscInt(0)
end

# ── Solve ─────────────────────────────────────────────────────────────────────
@show Threads.nthreads()
PETSc.solve!(x, snes)

if MPI.Comm_rank(comm) == 0
    its = LibPETSc.SNESGetIterationNumber(petsclib, snes)
    println("SNES converged in $its Newton iterations.")
end

# ── Cleanup ───────────────────────────────────────────────────────────────────
# Run a full GC now so any lingering VecRestoreArray finalizers from
# withlocalarray! run while PETSc is still valid, then barrier all ranks.
GC.gc(true)
MPI.Barrier(comm)

# SNES holds internal PETSc references to J, da, and r — destroy it first so
# those reference counts are decremented before we explicitly free the objects.
PETSc.destroy(snes)
PETSc.destroy(J)
PETSc.destroy(x_pert_vec)
PETSc.destroy(f0_vec)
PETSc.destroy(f1_vec)
PETSc.destroy(x)
PETSc.destroy(r)
PETSc.destroy(da)
PETSc.finalize(petsclib)

# On macOS ARM64 with MPICH ch4:ofi, MPICH's C atexit handler crashes during
# process teardown (SIGSEGV in libfabric/OFI cleanup).  Using quick_exit(0)
# after explicitly finalizing PETSc and MPI bypasses all C atexit() handlers
# (while still running at_quick_exit() handlers) and avoids the crash.
# All MPI communication is already complete at this point.
MPI.Barrier(comm)
MPI.Finalize()
ccall(:quick_exit, Cvoid, (Cint,), 0)