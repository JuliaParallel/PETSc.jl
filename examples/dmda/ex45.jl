# INCLUDE IN MPI TEST
#=
  ex45.jl: −∇²u = f on the unit cube with CG and geometric multigrid
  Port of PETSc ksp/tutorials/ex45.c

  Usage (from the examples/ directory):
    mpiexec -n 4 julia --project dmda/ex45.jl -N 33 -ksp_type cg -pc_type mg \
        -pc_mg_levels 3 -log_view

    # One GPU per rank: PETSc built with CUDA (e.g. PETSc.set_library!(:cuda)) and
    # CUDA.jl, or `-gpu hip` with HIP and AMDGPU.jl. Smoothers must run on the GPU.
    julia --project dmda/ex45.jl -gpu cuda -N 129 -ksp_type cg -pc_type mg \
        -pc_mg_levels 5 -mg_levels_ksp_type chebyshev -mg_levels_pc_type jacobi

  `Setup time` covers assembling every multigrid level and the multigrid setup,
  `Solve time` only the CG iterations.
=#
using PETSc, MPI, Printf, OffsetArrays
using LinearAlgebra: norm, axpy!
using PETSc: LibPETSc

# ============================================================================
#   Physics
# ============================================================================
# Exact solution, and the forcing that gives it: −∇²u = 3(2π)² u = 12π² u
u_exact(x, y, z) = sin(2π * x) * sin(2π * y) * sin(2π * z)
forcing(x, y, z) = 12π^2 * u_exact(x, y, z)

# 7-point finite differences with spacings hx, hy, hz, multiplied by the cell
# volume V = hx hy hz so that the matrix entries are O(h) on every multigrid
# level instead of O(1/h²). With wx = V/hx² (and wy, wz alike), row (i, j, k) is
#   wx (2u[i,j,k] − u[i-1,j,k] − u[i+1,j,k]) + (same along j and k) = V f[i,j,k]
# Boundary rows are the identity, with u_exact on the right-hand side.
rhs(x, y, z, on_boundary, V) = on_boundary ? u_exact(x, y, z) : V * forcing(x, y, z)

# ============================================================================
#   Discretization on the DMDA
# ============================================================================
# Grid index i (1-based, global) of an N-point grid on the unit cube
on_boundary(i, N) = any(i .== 1) || any(i .== N)
coordinates(i, N) = (i .- 1) ./ (N .- 1)

# Global matrix row (0-based) of every point in this rank's ghosted box,
# indexed by grid index, so row[i, j, k] also finds neighbors on other ranks
function global_rows(dm)
    ghost = PETSc.ghost_corners(dm)
    ltog = LibPETSc.DMGetLocalToGlobalMapping(petsclib, dm)
    indices = LibPETSc.ISLocalToGlobalMappingGetIndices(petsclib, ltog)
    rows = OffsetArray(reshape(copy(indices), ghost.size), Tuple(ghost.lower) .- 1)
    LibPETSc.ISLocalToGlobalMappingRestoreIndices(petsclib, ltog, indices)
    return rows
end

# The Laplacian rows this rank owns, as CSR arrays (0-based row offsets and
# column indices)
function laplacian_csr(dm)
    own = PETSc.corners(dm)
    N = PETSc.info(dm).global_size
    hx, hy, hz = 1 ./ (N .- 1)
    V = hx * hy * hz
    wx, wy, wz = V / hx^2, V / hy^2, V / hz^2
    row = global_rows(dm)
    PetscInt = petsclib.PetscInt
    rowptr, cols, vals = PetscInt[0], PetscInt[], Float64[]
    for I in own.lower:own.upper
        i, j, k = Tuple(I)
        if on_boundary((i, j, k), N)
            entries = [row[i, j, k] => 1.0]
        else
            entries = [
                row[i, j, k] => 2 * (wx + wy + wz),
                row[i-1, j, k] => -wx, row[i+1, j, k] => -wx,
                row[i, j-1, k] => -wy, row[i, j+1, k] => -wy,
                row[i, j, k-1] => -wz, row[i, j, k+1] => -wz,
            ]
        end
        # PETSc's CSR input wants the columns of a row in increasing order
        for (col, val) in sort!(entries)
            push!(cols, col); push!(vals, val)
        end
        push!(rowptr, length(cols))
    end
    return rowptr, cols, vals
end

# Fill the global vector v of the DMDA dm with f(x, y, z, on_boundary) at
# this rank's grid points
function fill_grid!(v, dm, f)
    own = PETSc.corners(dm)
    N = PETSc.info(dm).global_size
    function value_at(I)                    # I: global grid index of a point
        i = Tuple(I)
        x, y, z = coordinates(i, N)
        return f(x, y, z, on_boundary(i, N))
    end
    PETSc.with_local_array!(v; read = false) do v_array
        # A broadcast over this rank's points, so for a device Vec it runs on the GPU
        reshape(v_array, own.size) .= value_at.(own.lower:own.upper)
    end
    return v
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

petsclib = PETSc.getlib(; PetscScalar = Float64)
# The command line goes into the options database, so -log_view and the like apply
PETSc.initialize(petsclib; options = isinteractive() ? String[] : collect(String, ARGS))

# ============================================================================
#   Solve
# ============================================================================
toint(n) = n isa Integer ? Int(n) : parse(Int, n)

"""
    solve_ex45(N = 7; Nx = N, Ny = N, Nz = N, opts...)

Solve the Poisson problem on an `Nx × Ny × Nz` grid; `opts` are PETSc options
for the DMDA and the KSP. Sizes may be given as strings, as on the command line.
Returns `(; norm, final_grid, niter, setup_time, solve_time, L2, global_max)`.
"""
function solve_ex45(N = 7; Nx = N, Ny = N, Nz = N, opts...)
    comm = MPI.COMM_WORLD
    opts = (; opts..., gpu_options...)
    da = PETSc.DMDA(
        petsclib, comm, ntuple(_ -> PETSc.DM_BOUNDARY_NONE, 3), toint.((Nx, Ny, Nz)), 1;
        stencil_width = 1, stencil_type = PETSc.DMDA_STENCIL_STAR, opts...,
    )
    ksp = PETSc.KSP(da; opts...)

    # Called once per multigrid level, with that level's DM. Of the two CSR
    # calls, only the one matching the matrix (one rank or several) acts.
    PETSc.set_compute_operators!(ksp) do _, jac, ksp
        rowptr, cols, vals = laplacian_csr(PETSc.dm(ksp))
        LibPETSc.MatSeqAIJSetPreallocationCSR(petsclib, jac, rowptr, cols, vals)
        LibPETSc.MatMPIAIJSetPreallocationCSR(petsclib, jac, rowptr, cols, vals)
        return 0
    end

    final_grid = PETSc.info(da).global_size
    V = prod(1 ./ (final_grid .- 1))
    b = fill_grid!(PETSc.global_vec(da), da, (x, y, z, bnd) -> rhs(x, y, z, bnd, V))
    u = PETSc.global_vec(da)

    setup_time = @elapsed PETSc.setup!(ksp)
    solve_time = @elapsed PETSc.solve!(u, ksp, b)
    niter = LibPETSc.KSPGetIterationNumber(petsclib, ksp)
    res_norm = LibPETSc.KSPGetResidualNorm(petsclib, ksp)

    # err = u_exact − u; ‖err‖₂ √V approximates the L² norm on the cube
    err = fill_grid!(PETSc.global_vec(da), da, (x, y, z, _) -> u_exact(x, y, z))
    axpy!(-1, u, err)
    L2, global_max = sqrt(V) * norm(err, 2), norm(err, Inf)

    # Read by examples/scalability_tests/parse_scaling.jl; flushed so it comes
    # before PETSc's -log_view summary
    if MPI.Comm_rank(comm) == 0
        @printf("""
            Fine grid resolution: %d × %d × %d
            Outer KSP iterations: %d
            Final residual norm %g
            Setup time: %.6f seconds
            Solve time: %.6f seconds
            L2 error: %.6e
            Max error: %.6e
            """, final_grid..., niter, res_norm, setup_time, solve_time, L2, global_max)
        flush(stdout)
    end

    foreach(PETSc.destroy!, (err, u, b, ksp, da))
    return (; norm = res_norm, final_grid, niter, setup_time, solve_time, L2, global_max)
end

if !isinteractive() && abspath(PROGRAM_FILE) == @__FILE__
    opts = PETSc.parse_options(ARGS)
    solve_ex45(get(opts, :N, 7); opts...)
    PETSc.finalize(petsclib)
end
