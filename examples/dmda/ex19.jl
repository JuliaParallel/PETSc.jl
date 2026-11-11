#=
  ex19.jl: 2-D driven cavity in the velocity–vorticity–temperature formulation
  Port of PETSc snes/tutorials/ex19.c

  On the unit square, with the unknowns (u, v, ω, T) at every grid point:
    −∇²u − ∂ω/∂y                 = 0
    −∇²v + ∂ω/∂x                 = 0
    −∇²ω + ∇·(uω, vω) − Gr ∂T/∂x = 0
    −∇²T + Pr ∇·(uT, vT)         = 0
  with upwinded advection. All walls are no-slip, and the top lid moves with
  u = 1/(mx − 1). The left wall is cold (T = 0), the right one hot (T = 1 when
  Gr > 0), top and bottom are insulated, and ω follows from no-slip at each wall.

  Usage (from the examples/ directory):
    # Multigrid with Chebyshev/Jacobi smoothers on 4 ranks, with PETSc's log
    mpiexec -n 4 julia --project dmda/ex19.jl -da_grid_x 257 -da_grid_y 257 \
        -pc_type mg -pc_mg_levels 4 -mg_levels_ksp_type chebyshev \
        -mg_levels_pc_type jacobi -snes_monitor -log_view

    # One GPU per rank: PETSc built with CUDA (e.g. PETSc.set_library!(:cuda)) and
    # CUDA.jl, or `-gpu hip` with HIP and AMDGPU.jl
    julia --project dmda/ex19.jl -gpu cuda -da_grid_x 257 -da_grid_y 257 \
        -pc_type mg -pc_mg_levels 4 -mg_levels_ksp_type chebyshev -mg_levels_pc_type jacobi

  The residual is a KernelAbstractions kernel, so it runs on the CPU or the GPU.
  The finest grid's Jacobian comes from finite differences with a coloring of
  the stencil, on the device; coarser levels use PETSc's own coloring.
=#
using PETSc, MPI, KernelAbstractions
using PETSc: LibPETSc

# ============================================================================
#   Physics
# ============================================================================
const U, V, Ω, T = 1, 2, 3, 4   # the four unknowns at every grid point
const params = (; prandtl = 1.0, grashof = 1.0)

# The residual (F_u, F_v, F_ω, F_T) at global grid point (i, j), with every
# equation multiplied by hx hy, as in ex19.c. The point is at (a, b) in the
# ghosted array x_ghost, which also holds its neighbors.
# The helpers are anonymous functions (`f = q -> ...`): Julia would create a named
# one (`f(q) = ...`) before its variables are set and box them, which GPUs reject.
@inline function cavity_residual(x_ghost, a, b, i, j, grid, params)
    hx, hy = grid.hx, grid.hy
    # X(q) is unknown q at this point, X(q, di, dj) at the neighbor (i + di, j + dj)
    X = (q, di = 0, dj = 0) -> x_ghost[q, a + di, b + dj]

    # Walls: u = v = 0 (u = lid speed on top), ω from no-slip, T fixed left and
    # right, insulated bottom and top. Corners count as left or right wall.
    T_right = params.grashof > 0 ? 1.0 : 0.0
    if i == 1               # left
        return (X(U), X(V), X(Ω) - (X(V, 1, 0) - X(V)) / hx, X(T))
    elseif i == grid.mx     # right
        return (X(U), X(V), X(Ω) - (X(V) - X(V, -1, 0)) / hx, X(T) - T_right)
    elseif j == 1           # bottom
        return (X(U), X(V), X(Ω) + (X(U, 0, 1) - X(U)) / hy, X(T) - X(T, 0, 1))
    elseif j == grid.my     # top
        return (X(U) - hx, X(V), X(Ω) + (X(U) - X(U, 0, -1)) / hy, X(T) - X(T, 0, -1))  # lid speed 1/(mx − 1)
    end

    # Interior. −∇²q, ∂q/∂x and ∂q/∂y by centered differences, times hx hy
    laplacian = q -> (2X(q) - X(q, -1, 0) - X(q, 1, 0)) * hy / hx +
                     (2X(q) - X(q, 0, -1) - X(q, 0, 1)) * hx / hy
    d_dx = q -> (X(q, 1, 0) - X(q, -1, 0)) / 2 * hy
    d_dy = q -> (X(q, 0, 1) - X(q, 0, -1)) / 2 * hx

    # ∇·(uq, vq) times hx hy, upwinded: a backward difference where the
    # velocity is positive, a forward one where it is negative
    u, v = X(U), X(V)
    u⁺, u⁻ = max(u, 0), min(u, 0)
    v⁺, v⁻ = max(v, 0), min(v, 0)
    advection = q -> (u⁺ * (X(q) - X(q, -1, 0)) + u⁻ * (X(q, 1, 0) - X(q))) * hy +
                     (v⁺ * (X(q) - X(q, 0, -1)) + v⁻ * (X(q, 0, 1) - X(q))) * hx

    return (laplacian(U) - d_dy(Ω),
            laplacian(V) + d_dx(Ω),
            laplacian(Ω) + advection(Ω) - params.grashof * d_dx(T),
            laplacian(T) + params.prandtl * advection(T))
end

# One thread per grid point (li, lj) of this rank: writes its residual to F[:, li, lj]
@kernel function cavity_residual_kernel!(F, x_ghost, grid, params)
    li, lj = @index(Global, NTuple)
    i = grid.i_first + li - 1               # global grid index
    j = grid.j_first + lj - 1
    r = cavity_residual(x_ghost, li + grid.i_shift, lj + grid.j_shift, i, j, grid, params)
    for q in 1:4
        F[q, li, lj] = r[q]
    end
end

# The grid of one multigrid level, as this rank sees it
function grid_metrics(dm)
    mx, my = PETSc.info(dm).global_size
    own = PETSc.corners(dm)                 # the points this rank owns
    ghost = PETSc.ghost_corners(dm)         # the same, plus one layer of neighbors
    hx, hy = 1 / (mx - 1), 1 / (my - 1)
    return (;
        mx, my,                             # points in the whole grid
        hx, hy,                             # grid spacing
        i_first = own.lower[1], j_first = own.lower[2],   # global index of this rank's first point
        nx = own.size[1], ny = own.size[2],               # number of points this rank owns
        # The ghosted array starts with the neighbor layer, except at a wall:
        # this rank's first point sits at (1 + i_shift, 1 + j_shift) in it.
        i_shift = own.lower[1] - ghost.lower[1], j_shift = own.lower[2] - ghost.lower[2],
        nx_ghost = ghost.size[1], ny_ghost = ghost.size[2],
    )
end

# ============================================================================
#   Finite-difference Jacobian with a stencil coloring
# ============================================================================
# Columns of one color share no row, so one residual evaluation per color gives
# all their Jacobian entries: (F(x + h eₖ) − F(x)) / h for the columns k.

@kernel function perturb!(x, cols, h)
    k = @index(Global)
    x[cols[k]] += h
end

@kernel function difference!(vals, f0, f1, slots, rows, inv_h)
    k = @index(Global)
    vals[slots[k]] = (f1[rows[k]] - f0[rows[k]]) * inv_h
end

# The Jacobian matrix, the coloring of the grid of da and the work Vecs
function fd_jacobian_setup(da, x)
    coloring = PETSc.star_fd_coloring(da)
    # The Jacobian with the sparsity pattern of the stencil, filled in COO order
    J = LibPETSc.DMCreateMatrix(petsclib, da)
    LibPETSc.MatSetPreallocationCOOLocal(petsclib, J, LibPETSc.PetscCount(coloring.nnz_coo),
                                         coloring.row_coo_local_0b, coloring.col_coo_local_0b)

    # The Jacobian values, on the device for a device matrix. MatSetValuesCOO
    # takes an Array: on the GPU, one over the device memory, which PETSc reads
    # as device memory for a device matrix.
    vals = to_device(zeros(coloring.nnz_coo))
    coo_vals = gpu == "none" ? vals :
        unsafe_wrap(Array, Ptr{Float64}(UInt(pointer(vals))), length(vals))

    return (; J, vals, coo_vals,
              # for each color: the columns to perturb, the COO entries they
              # fill, and the residual rows those entries come from
              cols  = map(to_device, coloring.perturb_cols_1b),
              slots = map(to_device, coloring.coo_idxs_1b),
              rows  = map(to_device, coloring.local_rows_1b),
              # x + h eₖ, F(x) and F(x + h eₖ)
              x1 = LibPETSc.VecDuplicate(petsclib, x),
              f0 = LibPETSc.VecDuplicate(petsclib, x),
              f1 = LibPETSc.VecDuplicate(petsclib, x),
              h  = sqrt(eps()))
end

# Fill fd.vals with the finite-difference Jacobian at x
function fd_jacobian!(fd, snes, x)
    LibPETSc.SNESComputeFunction(petsclib, snes, x, fd.f0)
    PETSc.with_local_array!(fd.f0; write = false) do f0
        for c in eachindex(fd.cols)
            if isempty(fd.cols[c])          # a color with no columns on this rank
                continue
            end
            LibPETSc.VecCopy(petsclib, x, fd.x1)
            PETSc.with_local_array!(fd.x1) do x1
                backend = get_backend(x1)
                perturb!(backend, 64)(x1, fd.cols[c], fd.h; ndrange = length(fd.cols[c]))
                KernelAbstractions.synchronize(backend)
            end
            LibPETSc.SNESComputeFunction(petsclib, snes, fd.x1, fd.f1)
            PETSc.with_local_array!(fd.f1; write = false) do f1
                backend = get_backend(f1)
                difference!(backend, 64)(fd.vals, f0, f1, fd.slots[c], fd.rows[c], 1 / fd.h;
                                         ndrange = length(fd.slots[c]))
                KernelAbstractions.synchronize(backend)
            end
        end
    end
    return nothing
end

# ============================================================================
#   GPU selection (-gpu cuda or -gpu hip)
# ============================================================================
# Rank number within its node, to give each rank of a node its own GPU
function node_rank()
    MPI.Initialized() || MPI.Init()
    node = MPI.Comm_split_type(MPI.COMM_WORLD, MPI.COMM_TYPE_SHARED, 0)
    r = MPI.Comm_rank(node)
    MPI.free(node)
    return r
end

const gpu = get(PETSc.parse_options(ARGS), :gpu, "none")
if gpu == "cuda"
    using CUDA
elseif gpu == "hip"
    using AMDGPU
elseif gpu != "none"
    error("-gpu takes cuda or hip, got $gpu")
end

# Each rank picks a GPU, and the DM gets GPU Vec and Mat types, which every
# multigrid level inherits. PETSc uses the GPU selected here, so this comes
# before PETSc.initialize. (A separate block from the `using` above, so that
# the functions of CUDA.jl or AMDGPU.jl are available.)
if gpu == "cuda"
    CUDA.device!(node_rank() % length(CUDA.devices()))
    const gpu_options = (; dm_vec_type = "cuda", dm_mat_type = "aijcusparse")
elseif gpu == "hip"
    AMDGPU.device!(AMDGPU.devices()[node_rank() % length(AMDGPU.devices()) + 1])
    const gpu_options = (; dm_vec_type = "hip", dm_mat_type = "aijhipsparse")
else
    const gpu_options = (;)
end

# A copy of the array `a` on the GPU, or `a` itself on the CPU
function to_device(a)
    if gpu == "cuda"
        return CuArray(a)
    elseif gpu == "hip"
        return ROCArray(a)
    else
        return a
    end
end

petsclib = PETSc.getlib(; PetscScalar = Float64)
# The command line goes into the options database, which the DMDA and the SNES
# read: -da_grid_x, -snes_monitor, -log_view and the like
PETSc.initialize(petsclib; options = isinteractive() ? String[] : collect(String, ARGS))

# ============================================================================
#   Solve
# ============================================================================
comm = MPI.COMM_WORLD
# 4 × 4 by default, as in ex19.c
da = PETSc.DMDA(petsclib, comm, (PETSc.DM_BOUNDARY_NONE, PETSc.DM_BOUNDARY_NONE), (4, 4), 4;
                stencil_width = 1, stencil_type = PETSc.DMDA_STENCIL_STAR, gpu_options...)
snes = PETSc.SNES(petsclib, comm)
PETSc.set_dm!(snes, da)

# Initial guess: u = v = ω = 0, and T linear in x when Gr > 0
x = PETSc.global_vec(da)
PETSc.with_local_array!(x; read = false) do x_array
    grid = grid_metrics(da)
    X = reshape(x_array, 4, grid.nx, grid.ny)    # X[q, i, j]: unknown q at point (i, j)
    fill!(X, 0)
    if params.grashof > 0
        x_coordinates = (grid.i_first - 1 .+ (0:grid.nx - 1)) .* grid.hx
        X[T, :, :] .= x_coordinates             # the same in every row j
    end
end

# The residual F(x), on any multigrid level. The kernel needs the neighbors of
# this rank's points too, so x is first copied into a ghosted local Vec.
F = similar(x)
PETSc.set_function!(snes, F) do F, snes, x
    dm = PETSc.dm(snes)
    grid = grid_metrics(dm)
    x_local = PETSc.local_vec(dm)
    PETSc.global_to_local!(x_local, dm, x; mode = PETSc.INSERT_VALUES)
    # Both arrays come out on the same side: on the GPU when either Vec is there
    # (PETSc's coloring on coarser levels can pass a host Vec)
    PETSc.with_local_arrays!(F, x_local) do F_array, x_array
        backend = get_backend(F_array)
        cavity_residual_kernel!(backend, 64)(
            reshape(F_array, 4, grid.nx, grid.ny),
            reshape(x_array, 4, grid.nx_ghost, grid.ny_ghost),
            grid, params;
            ndrange = (grid.nx, grid.ny),       # one thread per owned point
        )
        KernelAbstractions.synchronize(backend)
    end
    PETSc.destroy!(x_local)
    return 0
end

# The Jacobian. Multigrid calls this for every level, each with its own SNES and
# DM: on the finest grid (da) our finite differences fill it, on coarser grids
# PETSc's own finite-difference coloring does.
fd = fd_jacobian_setup(da, x)
PETSc.set_snes_jacobian!(snes, fd.J) do J, snes, x
    on_finest_grid = PETSc.dm(snes).ptr == da.ptr
    if on_finest_grid
        fd_jacobian!(fd, snes, x)
        LibPETSc.MatSetValuesCOO(petsclib, J, fd.coo_vals, LibPETSc.INSERT_VALUES)
        PETSc.assemble!(J)
    else
        LibPETSc.SNESComputeJacobianDefaultColor(petsclib, snes, x, J, J, C_NULL)
    end
    return 0
end

PETSc.solve!(x, snes)
if MPI.Comm_rank(comm) == 0
    println("SNES converged in $(LibPETSc.SNESGetIterationNumber(petsclib, snes)) Newton iterations.")
end

foreach(PETSc.destroy!, (snes, fd.J, fd.x1, fd.f0, fd.f1, x, F, da))
PETSc.finalize(petsclib)
