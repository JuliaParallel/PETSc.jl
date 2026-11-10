# INCLUDE IN MPI TEST
#=
  ex45.jl: −∇²u = f on the unit cube with CG and geometric multigrid
  Port of PETSc ksp/tutorials/ex45.c

  The exact solution is u = sin(2πx) sin(2πy) sin(2πz), so f = 12π² u, with
  u = 0 on the boundary. The 7-point Laplacian is assembled on every multigrid
  level, and the error against u is printed at the end.

  Usage (from the examples/ directory):
    julia --project dmda/ex45.jl -N 33 -ksp_type cg -pc_type mg -pc_mg_levels 3

    # Non-cubic grid on 4 ranks, with PETSc's performance log
    mpiexec -n 4 julia --project dmda/ex45.jl -Nx 65 -Ny 65 -Nz 129 \
        -ksp_type cg -pc_type mg -pc_mg_levels 4 -log_view

    # One GPU per rank. Needs PETSc built with CUDA (e.g. PETSc.set_library!(:cuda))
    # and CUDA.jl, or `-gpu hip` with HIP and AMDGPU.jl. The smoothers must run
    # on the GPU too, so Chebyshev/Jacobi instead of SOR.
    julia --project dmda/ex45.jl -gpu cuda -N 129 \
        -ksp_type cg -pc_type mg -pc_mg_levels 5 \
        -mg_levels_ksp_type chebyshev -mg_levels_pc_type jacobi

  The right-hand side and the exact solution are filled by a broadcast, so the
  same code fills host and device vectors. `Setup time` covers the operator assembly on all levels and the multigrid
  setup, `Solve time` only the CG iterations.
=#
using PETSc, MPI, Printf
using LinearAlgebra: norm, axpy!
using PETSc: LibPETSc

const gpu = Symbol(get(PETSc.parse_options(ARGS), :gpu, "none"))
gpu in (:none, :cuda, :hip) || throw(ArgumentError("-gpu takes cuda or hip, got $gpu"))

if gpu === :cuda
    using CUDA
elseif gpu === :hip
    using AMDGPU
end

# Vec and Mat types for the DM and, through it, every multigrid level
const gpu_options = (;
    none = (;),
    cuda = (; dm_vec_type = "cuda", dm_mat_type = "aijcusparse" ),
    hip  = (; dm_vec_type = "hip" , dm_mat_type = "aijhipsparse"),
)[gpu]

exact(x, y, z) = sin(2π * x) * sin(2π * y) * sin(2π * z)
forcing(x, y, z) = 12π^2 * exact(x, y, z)

petsclib = PETSc.getlib(; PetscScalar = Float64)

# Each rank uses the GPU (node-local rank) mod (GPUs per node). PETSc picks up
# the device selected here, so this runs before PETSc.initialize.
function select_device()
    MPI.Initialized() || MPI.Init()
    node = MPI.Comm_split_type(MPI.COMM_WORLD, MPI.COMM_TYPE_SHARED, 0)
    r = MPI.Comm_rank(node)
    MPI.free(node)
    if gpu === :cuda
        CUDA.device!(r % length(CUDA.devices()))
    else
        devs = AMDGPU.devices()
        AMDGPU.device!(devs[r % length(devs) + 1])
    end
    return nothing
end
gpu === :none || select_device()

# The command line goes into the options database, so -log_view and the like apply
PETSc.initialize(petsclib; options = isinteractive() ? String[] : collect(String, ARGS))

on_boundary(g, N) = any(g .== 1) || any(g .== N)

# The 7-point Laplacian of the rows this rank owns, as COO triplets in the DM's
# ghosted local numbering (0-based), with identity rows on the boundary.
# Scaled by the cell volume HxHyHz so that the entries are O(h) on every
# multigrid level instead of O(1/h²).
function laplacian_coo(dm)
    own = PETSc.corners(dm)
    ghost = PETSc.ghost_corners(dm)
    N = PETSc.info(dm).global_size
    H = 1 ./ (N .- 1)
    c = prod(H) ./ H .^ 2    # (HyHz/Hx, HxHz/Hy, HxHy/Hz)
    neighbors = (
        (CartesianIndex(-1, 0, 0), -c[1]), (CartesianIndex(1, 0, 0), -c[1]),
        (CartesianIndex(0, -1, 0), -c[2]), (CartesianIndex(0, 1, 0), -c[2]),
        (CartesianIndex(0, 0, -1), -c[3]), (CartesianIndex(0, 0, 1), -c[3]),
    )
    # Global grid index → 0-based position in the ghosted local array
    ghost_linear = LinearIndices(ghost.size)
    local_index(I) = ghost_linear[I - ghost.lower + oneunit(I)] - 1

    PetscInt = petsclib.PetscInt
    rows, cols, vals = PetscInt[], PetscInt[], Float64[]
    function add!(i, j, v)
        push!(rows, local_index(i))
        push!(cols, local_index(j))
        push!(vals, v)
    end
    for I in own.lower:own.upper
        if on_boundary(Tuple(I), N)
            add!(I, I, 1.0)
        else
            for (d, v) in neighbors
                add!(I, I + d, v)
            end
            add!(I, I, 2 * sum(c))
        end
    end
    return rows, cols, vals
end

# f(x, y, z, on_boundary) at the local grid point I; `lo` is the global grid
# index of local point (1, 1, 1) and H the grid spacing
function grid_value(f, I, lo, N, H)
    g = Tuple(I) .+ lo .- 1
    x, y, z = (g .- 1) .* H
    return f(x, y, z, on_boundary(g, N))
end

# Fill the global vector `v` of the DMDA `dm` with f(x, y, z, on_boundary).
# A broadcast, so it runs on the device for a device Vec.
function grid_eval!(v, dm, f)
    own = PETSc.corners(dm)
    N = PETSc.info(dm).global_size
    H = 1 ./ (N .- 1)
    # A device array for a device Vec once CUDA.jl or AMDGPU.jl is loaded
    PETSc.with_local_array!(v; read = false) do a
        a3 = reshape(a, own.size)
        a3 .= grid_value.(Ref(f), CartesianIndices(a3), Ref(Tuple(own.lower)), Ref(N), Ref(H))
    end
    return v
end

toint(n::Integer) = Int(n)
toint(n::AbstractString) = parse(Int, n)

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

    # Called once per multigrid level, with that level's DM
    PETSc.set_compute_operators!(ksp) do _, jac, ksp
        rows, cols, vals = laplacian_coo(PETSc.dm(ksp))
        LibPETSc.MatSetPreallocationCOOLocal(petsclib, jac, LibPETSc.PetscCount(length(rows)), rows, cols)
        LibPETSc.MatSetValuesCOO(petsclib, jac, vals, LibPETSc.INSERT_VALUES)
        PETSc.assemble!(jac)
        return 0
    end

    final_grid = PETSc.info(da).global_size
    vol = prod(1 ./ (final_grid .- 1))
    # Boundary rows are identity rows, so b holds u there; inside it is f scaled like the matrix
    b = grid_eval!(PETSc.global_vec(da), da, (x, y, z, bnd) -> bnd ? exact(x, y, z) : vol * forcing(x, y, z))
    x = PETSc.global_vec(da)

    setup_time = @elapsed PETSc.setup!(ksp)
    solve_time = @elapsed PETSc.solve!(x, ksp, b)
    niter = LibPETSc.KSPGetIterationNumber(petsclib, ksp)
    res_norm = LibPETSc.KSPGetResidualNorm(petsclib, ksp)

    # err = u − x; ‖err‖₂ √vol approximates the L² norm on the cube
    err = grid_eval!(PETSc.global_vec(da), da, (x, y, z, _) -> exact(x, y, z))
    axpy!(-1, x, err)
    L2 = sqrt(vol) * norm(err, 2)
    global_max = norm(err, Inf)

    # One value per line: examples/scalability_tests/parse_scaling.jl reads these
    if MPI.Comm_rank(comm) == 0
        @printf("Fine grid resolution: %d × %d × %d\n", final_grid...)
        @printf("Outer KSP iterations: %d\n", niter)
        @printf("Final residual norm %g\n", res_norm)
        @printf("Setup time: %.6f seconds\n", setup_time)
        @printf("Solve time: %.6f seconds\n", solve_time)
        @printf("L2 error: %.6e\n", L2)
        @printf("Max error: %.6e\n", global_max)
        flush(stdout)    # before PETSc's -log_view summary, which the parser expects after these
    end

    PETSc.destroy!.((err, x, b))
    PETSc.destroy!(ksp)
    PETSc.destroy!(da)
    return (; norm = res_norm, final_grid, niter, setup_time, solve_time, L2, global_max)
end

if !isinteractive() && abspath(PROGRAM_FILE) == @__FILE__
    opts = PETSc.parse_options(ARGS)
    solve_ex45(get(opts, :N, 7); opts...)
    PETSc.finalize(petsclib)
end
