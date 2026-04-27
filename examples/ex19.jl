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
    julia --project ex19.jl
    julia --project ex19.jl -snes_monitor -da_grid_x 129 -da_grid_y 129
    julia --project ex19.jl -snes_monitor -da_grid_x 129 -da_grid_y 129 -log_view
    mpiexec -n 4 julia --project ex19.jl -snes_monitor -pc_type mg -da_grid_x 64 -da_grid_y 64

  Requires: LocalPreferences.toml in examples/ with PetscInt = "Int32" matching the
  PETSc build (check with: grep sizeof_PetscInt petscconf.h).

  GPU usage: set  useCUDA = true  then run as above.
    Requires PETSc built with --with-cuda, and CUDA.jl in the environment.
=#

# ── GPU switch ────────────────────────────────────────────────────────────────
const useCUDA = false

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

# ── Setup ─────────────────────────────────────────────────────────────────────
opts     = isinteractive() ? NamedTuple() : PETSc.parse_options(filter(a -> a != "-log_view", ARGS))
log_view = "-log_view" in ARGS

petsclib = PETSc.getlib(; PetscScalar = Float64, PetscInt = Int32)
PETSc.initialize(petsclib; log_view)

T        = Float64
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

# Stage 2: GPU vecs and GPU matrix enable a fully GPU-resident FD coloring
# path via COO preallocation.  No host↔device bouncing in residual or Jacobian.
# NOTE: MG with GPU vecs needs per-level COO rebuild (not yet implemented);
#       for MG tests use -pc_type lu or default ILU.
if useCUDA
    GC.@preserve begin
        vt = "cuda"
        mt = "aijcusparse"
        LibPETSc.DMSetVecType(petsclib, da, pointer(vt))
        LibPETSc.DMSetMatType(petsclib, da, pointer(mt))
    end
end

snes = PETSc.SNES(petsclib, comm; opts...)
PETSc.setDM!(snes, da)

# Actual grid size after setfromoptions (may differ from the 4×4 default)
info = PETSc.getinfo(da)
mx   = Int(info.global_size[1])
my   = Int(info.global_size[2])

user = AppCtx{T}(
    lidvelocity = T(1) / (mx - 1),
    prandtl     = T(1),
    grashof     = T(1),
)

# Precomputed grid metrics
dhx   = T(mx - 1);   dhy   = T(my - 1)
hx    = one(T) / dhx; hy    = one(T) / dhy
hydhx = hy * dhx;    hxdhy = hx * dhy

# ── Initial condition: u = v = ω = 0, T linear in x ─────────────────────────
x = PETSc.DMGlobalVec(da)

PETSc.withlocalarray!(x; read = false) do x_arr
    corners = PETSc.getcorners(da)
    xs = corners.lower[1];  ys = corners.lower[2]
    xe = corners.upper[1];  ye = corners.upper[2]
    nx_own = xe - xs + 1;   ny_own = ye - ys + 1
    dx = one(T) / (mx - 1)
    x_par = reshape(x_arr, 4, nx_own, ny_own)
    for lj in 1:ny_own, li in 1:nx_own
        ig = xs + li - 1
        x_par[1, li, lj] = zero(T)
        x_par[2, li, lj] = zero(T)
        x_par[3, li, lj] = zero(T)
        x_par[4, li, lj] = user.grashof > 0 ? T(ig - 1) * dx : zero(T)
    end
end

# ── Helpers for PETSc array access with optional GPU residual kernel ──────────
#
# Three cases:
#
#   useCUDA=false  →  plain CPU arrays via unsafe_localarray.
#
#   useCUDA=true, vecs on GPU (PETSC_MEMTYPE_CUDA)
#              →  zero-copy CuArray wraps; no host↔device transfer.
#
#   useCUDA=true, vecs on CPU (PETSC_MEMTYPE_HOST, Stage-1 coloring path)
#              →  allocate GPU scratch buffers, copy lx H2D before kernel,
#                 copy fx D2H after kernel, then let PETSc see the result in
#                 the HOST fx_arr.  FD coloring arithmetic stays on CPU (BLAS)
#                 which avoids the VecPlaceArray+cuBLAS bug.
#
# Returns (fx, lx, fx_arr, lx_arr, fx_bounce)
#   fx / lx       — arrays passed to the kernel (CuArray or plain Array)
#   fx_arr/lx_arr — raw PETSc handles for Restore calls (nothing on CPU path)
#   fx_bounce     — CuArray whose contents must be copied back to fx_arr after
#                   the kernel; nothing when no copy is needed.
#
function get_petsc_arrays(petsclib, g_fx, l_x)
    if useCUDA
        fx_arr, fx_mtype = LibPETSc.VecGetArrayAndMemType(petsclib, g_fx)
        lx_arr, lx_mtype = LibPETSc.VecGetArrayReadAndMemType(petsclib, l_x)

        if lx_mtype == LibPETSc.PETSC_MEMTYPE_DEVICE
            # Native GPU vecs: zero-copy wrap, no bounce needed.
            fx = unsafe_wrap(CuArray, CuPtr{T}(UInt64(pointer(fx_arr))), length(fx_arr))
            lx = unsafe_wrap(CuArray, CuPtr{T}(UInt64(pointer(lx_arr))), length(lx_arr))
            return fx, lx, fx_arr, lx_arr, nothing
        else
            # CPU vecs (FD coloring path): bounce residual through GPU.
            lx_gpu = CuArray{T}(undef, length(lx_arr))
            fx_gpu = CuArray{T}(undef, length(fx_arr))
            copyto!(lx_gpu, lx_arr)          # H2D: send ghost input to GPU
            return fx_gpu, lx_gpu, fx_arr, lx_arr, fx_gpu
        end
    else
        fx = PETSc.unsafe_localarray(g_fx; read = true,  write = true)
        lx = PETSc.unsafe_localarray(l_x;  read = true,  write = false)
        return fx, lx, nothing, nothing, nothing
    end
end

function restore_petsc_arrays(petsclib, g_fx, l_x, fx, lx, fx_arr, lx_arr, fx_bounce)
    if useCUDA
        if fx_bounce !== nothing
            # D2H: copy GPU residual result back to the HOST PETSc array.
            CUDA.synchronize()
            copyto!(fx_arr, fx_bounce)
        end
        LibPETSc.VecRestoreArrayAndMemType(petsclib, g_fx, fx_arr)
        LibPETSc.VecRestoreArrayReadAndMemType(petsclib, l_x, lx_arr)
    else
        Base.finalize(fx)
        Base.finalize(lx)
    end
end

# ── Residual callback ─────────────────────────────────────────────────────────
r = similar(x)

PETSc.setfunction!(snes, r) do g_fx, snes, g_x
    da = PETSc.getDM(snes)

    l_x = PETSc.DMLocalVec(da)
    PETSc.dm_global_to_local!(g_x, l_x, da, PETSc.INSERT_VALUES)

    fx, lx, fx_arr, lx_arr, fx_bounce = get_petsc_arrays(petsclib, g_fx, l_x)

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
    dhx_   = T(mx_ - 1);    dhy_   = T(my_ - 1)
    hx_    = one(T) / dhx_; hy_    = one(T) / dhy_
    hydhx_ = hy_ * dhx_;    hxdhy_ = hx_ * dhy_
    lid_   = T(1) / dhx_    # lidvelocity = 1/(mx-1)

    cavity_residual_kernel!(backend, 64)(
        f_par, x_par,
        dhx_, dhy_, hx_, hy_, hydhx_, hxdhy_,
        user.grashof, user.prandtl, lid_,
        mx_, my_, xs, ys, ox, oy;
        ndrange = (nx_own, ny_own),
    )
    KernelAbstractions.synchronize(backend)

    restore_petsc_arrays(petsclib, g_fx, l_x, fx, lx, fx_arr, lx_arr, fx_bounce)
    PETSc.destroy(l_x)
    return PetscInt(0)
end

# ── Jacobian: manual FD coloring with GPU-efficient COO matrix assembly ───────
#
# 1) Obtain the DM's ISColoring (PETSc determines which columns can be
#    perturbed simultaneously without touching the same nonzero row twice).
# 2) Build the sparse (row, col) COO triplets analytically from the DMDA STAR
#    stencil; store per-COO the 0-based color of its column.
# 3) Pre-allocate J via MatSetPreallocationCOO (one-time GPU setup).
# 4) Each Newton step: loop over colors, scatter +h to owned cols of that
#    color (GPU kernel), evaluate F(x_pert), accumulate (F1−F0)/h into
#    val[] (GPU kernel), assemble via MatSetValuesCOO(J, val_dev).
#
# NOTE: Assumes serial (single MPI rank).  For parallel runs ghost-column
#       colors must be communicated; see ISColoringGetColors documentation.
# NOTE: For MG, coarser levels fall back to SNESComputeJacobianDefaultColor
#       (correct, but CPU-only FD coloring for those levels).

# ── 1. ISColoring ─────────────────────────────────────────────────────────────
iscoloring = LibPETSc.DMCreateColoring(petsclib, da, LibPETSc.IS_COLORING_GLOBAL)

# ── 2. Per-column color via raw ISColoringGetColors call ──────────────────────
#   C API:  ISColoringGetColors(iscoloring, PetscInt *n, PetscInt *nc,
#                               const ISColoringValue **colors)
#   ISColoringValue = unsigned short (UInt16) per petscconf.h PETSC_IS_COLORING_VALUE_TYPE=short
n_cols_ref     = Ref{PetscInt}(0)
nc_ref         = Ref{PetscInt}(0)
colors_ptr_ref = Ref{Ptr{UInt16}}(C_NULL)
LibPETSc.@chk ccall(
    (:ISColoringGetColors, petsclib.petsc_library), PetscInt,
    (LibPETSc.ISColoring, Ptr{PetscInt}, Ptr{PetscInt}, Ptr{Ptr{UInt16}}),
    iscoloring, n_cols_ref, nc_ref, colors_ptr_ref)
n_cols   = Int(n_cols_ref[])
n_colors = Int(nc_ref[])
# Copy colors to an owned Julia array before we destroy the ISColoring.
col_colors_host = copy(unsafe_wrap(Vector{UInt16}, colors_ptr_ref[], n_cols; own = false))

# ── 3. Ownership range (0-based PETSc indices) ────────────────────────────────
row_start, row_end = LibPETSc.VecGetOwnershipRange(petsclib, x)
col_start = row_start
@assert n_cols == Int(row_end - row_start) "serial assumption: n_cols ($n_cols) != n_owned_dofs ($(Int(row_end-row_start)))"
n_local_dofs = Int(row_end - row_start)

# ── 4. Build COO from DMDA STAR stencil ───────────────────────────────────────
#  For each owned node (ii, jj) and each stencil neighbor, emit dof×dof
#  (row, col) pairs.  Every entry also records the 0-based color of its column
#  and the 0-based local row index (row_global − row_start).
dof_per_node = 4
coo_corners = PETSc.getcorners(da)
xs_da = coo_corners.lower[1];  ys_da = coo_corners.lower[2]
xe_da = coo_corners.upper[1];  ye_da = coo_corners.upper[2]

CPetscInt = petsclib.PetscInt           # matches the actual C sizeof(PetscInt)
row_coo_host      = CPetscInt[]
col_coo_host      = CPetscInt[]
local_row_per_coo = CPetscInt[]   # 0-based local row  (row_global − row_start)
color_per_coo     = CPetscInt[]   # 0-based color of each COO entry's column

for jj in ys_da:ye_da, ii in xs_da:xe_da
    ig = ii - 1;  jg = jj - 1           # 0-based global node coords
    row_base = (jg * mx + ig) * dof_per_node

    neighbors = Tuple{Int,Int}[(ii, jj)]
    ii > 1  && push!(neighbors, (ii-1, jj))
    ii < mx && push!(neighbors, (ii+1, jj))
    jj > 1  && push!(neighbors, (ii, jj-1))
    jj < my && push!(neighbors, (ii, jj+1))

    for (ni, nj) in neighbors
        nig = ni - 1;  njg = nj - 1
        col_base = (njg * mx + nig) * dof_per_node
        for d_row in 0:dof_per_node-1, d_col in 0:dof_per_node-1
            r_g = row_base + d_row
            c_g = col_base + d_col
            push!(row_coo_host, CPetscInt(r_g))
            push!(col_coo_host, CPetscInt(c_g))
            # Serial: local col index = c_g  (col_start == 0)
            push!(color_per_coo,     CPetscInt(col_colors_host[c_g - Int(col_start) + 1]))
            push!(local_row_per_coo, CPetscInt(r_g - Int(row_start)))
        end
    end
end
nnz_coo = length(row_coo_host)

# ── 5. Create J with COO preallocation ────────────────────────────────────────
J = LibPETSc.DMCreateMatrix(petsclib, da)
# Direct ccall: use CPetscInt (= petsclib.PetscInt) so this works for both
# 32-bit and 64-bit PETSc builds.  PetscCount = ptrdiff_t = Int64 always.
LibPETSc.@chk ccall(
    (:MatSetPreallocationCOO, petsclib.petsc_library), Cint,
    (LibPETSc.CMat, Int64, Ptr{CPetscInt}, Ptr{CPetscInt}),
    J, Int64(nnz_coo), row_coo_host, col_coo_host)

# ── 6. Per-color index arrays for the FD loop ─────────────────────────────────
# perturb_cols_1b[c]: 1-based local x-indices of owned columns with color c-1.
# coo_idxs_1b[c]:    1-based COO entry indices whose column color == c-1.
# local_rows_1b[c]:  1-based local residual-row indices for those COO entries.
perturb_cols_1b = [Int32[] for _ in 1:n_colors]
for k_local in 1:n_cols
    c = Int(col_colors_host[k_local]) + 1   # 1-based color
    push!(perturb_cols_1b[c], Int32(k_local))
end

coo_idxs_1b   = [Int32[] for _ in 1:n_colors]
local_rows_1b = [Int32[] for _ in 1:n_colors]
for k in 1:nnz_coo
    c = Int(color_per_coo[k]) + 1            # 1-based color
    push!(coo_idxs_1b[c],   Int32(k))
    push!(local_rows_1b[c], Int32(local_row_per_coo[k] + 1))  # 0→1-based
end

if useCUDA
    perturb_cols_dev = [CuArray(v) for v in perturb_cols_1b]
    coo_idxs_dev     = [CuArray(v) for v in coo_idxs_1b]
    local_rows_dev   = [CuArray(v) for v in local_rows_1b]
    val_dev          = CUDA.zeros(T, nnz_coo)
else
    perturb_cols_dev = perturb_cols_1b
    coo_idxs_dev     = coo_idxs_1b
    local_rows_dev   = local_rows_1b
    val_dev          = zeros(T, nnz_coo)
end

LibPETSc.ISColoringDestroy(petsclib, iscoloring)

# ── 7. Scratch vectors for the FD loop ────────────────────────────────────────
x_pert_vec = LibPETSc.VecDuplicate(petsclib, x)
f0_vec     = LibPETSc.VecDuplicate(petsclib, x)
f1_vec     = LibPETSc.VecDuplicate(petsclib, x)
h_eps      = T(sqrt(eps(T)))
inv_h      = T(1) / h_eps

# ── 8. Custom Jacobian callback ───────────────────────────────────────────────
PETSc.setjacobian!(snes, J) do Jmat, actual_snes, g_x
    # For MG: if this is a coarser level (grid size differs from fine grid),
    # fall back to PETSc's built-in FD coloring (correct for that level's DM).
    da_level   = PETSc.getDM(actual_snes)
    info_level = PETSc.getinfo(da_level)
    if info_level.global_size[1] != mx || info_level.global_size[2] != my
        LibPETSc.@chk ccall(
            (:SNESComputeJacobianDefaultColor, petsclib.petsc_library), PetscInt,
            (LibPETSc.CSNES, LibPETSc.CVec, LibPETSc.CMat, LibPETSc.CMat, Ptr{Cvoid}),
            actual_snes.ptr, g_x.ptr, Jmat.ptr, Jmat.ptr, C_NULL)
        return PetscInt(0)
    end

    # ── Evaluate F(x) → f0 ────────────────────────────────────────────────────
    LibPETSc.SNESComputeFunction(petsclib, actual_snes, g_x, f0_vec)
    f0_arr, f0_mtype = LibPETSc.VecGetArrayReadAndMemType(petsclib, f0_vec)
    f0_dev = (useCUDA && f0_mtype == LibPETSc.PETSC_MEMTYPE_DEVICE) ?
        unsafe_wrap(CuArray, CuPtr{T}(UInt64(pointer(f0_arr))), n_local_dofs) :
        f0_arr

    # ── FD loop over colors ────────────────────────────────────────────────────
    for c in 1:n_colors
        isempty(perturb_cols_dev[c]) && continue

        # Copy x → x_pert, then scatter +h to owned cols of color c.
        LibPETSc.VecCopy(petsclib, g_x, x_pert_vec)
        xp_arr, xp_mtype = LibPETSc.VecGetArrayAndMemType(petsclib, x_pert_vec)
        xp_dev = (useCUDA && xp_mtype == LibPETSc.PETSC_MEMTYPE_DEVICE) ?
            unsafe_wrap(CuArray, CuPtr{T}(UInt64(pointer(xp_arr))), n_local_dofs) :
            xp_arr
        scatter_perturb_kernel!(backend, 64)(
            xp_dev, perturb_cols_dev[c], h_eps;
            ndrange = length(perturb_cols_dev[c]))
        KernelAbstractions.synchronize(backend)
        LibPETSc.VecRestoreArrayAndMemType(petsclib, x_pert_vec, xp_arr)

        # Evaluate F(x_pert) → f1.
        LibPETSc.SNESComputeFunction(petsclib, actual_snes, x_pert_vec, f1_vec)

        # Accumulate (f1 − f0)/h into val[] at the COO indices of color c.
        f1_arr, f1_mtype = LibPETSc.VecGetArrayReadAndMemType(petsclib, f1_vec)
        f1_dev = (useCUDA && f1_mtype == LibPETSc.PETSC_MEMTYPE_DEVICE) ?
            unsafe_wrap(CuArray, CuPtr{T}(UInt64(pointer(f1_arr))), n_local_dofs) :
            f1_arr
        fd_accumulate_kernel!(backend, 64)(
            val_dev, f0_dev, f1_dev,
            coo_idxs_dev[c], local_rows_dev[c], inv_h;
            ndrange = length(coo_idxs_dev[c]))
        KernelAbstractions.synchronize(backend)
        LibPETSc.VecRestoreArrayReadAndMemType(petsclib, f1_vec, f1_arr)
    end

    LibPETSc.VecRestoreArrayReadAndMemType(petsclib, f0_vec, f0_arr)

    # ── Assemble J via COO ─────────────────────────────────────────────────────
    # For GPU matrix (aijcusparse): pass device pointer so PETSc's CUDA kernel
    # scatters val[] into CSR storage entirely on device (no D2H transfer).
    if useCUDA
        LibPETSc.@chk ccall(
            (:MatSetValuesCOO, petsclib.petsc_library), PetscInt,
            (LibPETSc.CMat, Ptr{T}, LibPETSc.InsertMode),
            Jmat.ptr, Ptr{T}(UInt64(pointer(val_dev))), LibPETSc.INSERT_VALUES)
    else
        LibPETSc.MatSetValuesCOO(petsclib, Jmat, val_dev, LibPETSc.INSERT_VALUES)
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
# Explicitly destroy the PetscOptions stored on the SNES before finalization.
# Its GC finalizer calls PetscOptionsDestroy, which can use MPI internally.
# If GC runs it after MPI is alive but in a different collective-sync state
# across ranks, it triggers intermittent crashes.  Destroying it explicitly
# here (while all ranks are synchronized and PETSc/MPI are still fully active)
# is safe and prevents any later GC-driven call.
if !isnothing(snes.opts)
    PETSc.destroy(snes.opts)
    snes.opts = nothing
end

# Run a full GC now so any lingering VecRestoreArray finalizers from
# withlocalarray! run while PETSc is still valid, then barrier all ranks.
GC.gc(true)
MPI.Barrier(comm)

# SNES holds internal PETSc references to J, da, and r — destroy it first so
# those reference counts are decremented before we explicitly free the objects.
PETSc.destroy(snes)
PETSc.destroy(J)
LibPETSc.VecDestroy(petsclib, x_pert_vec)
LibPETSc.VecDestroy(petsclib, f0_vec)
LibPETSc.VecDestroy(petsclib, f1_vec)
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