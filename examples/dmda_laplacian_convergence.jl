# EXCLUDE FROM TESTING
# # Finite Difference Laplacian - Convergence Analysis
#
# This example performs convergence analysis for a 2-D Laplacian solver with zero 
# Dirichlet boundary conditions using the PETSc DMDA interface.
#
# The approach taken is that the 5-point Laplacian is applied in the interior
# and boundary points are included in the solve, but the forcing for the
# boundary data is the Dirichlet boundary condition and the matrix entry for
# these points is set to the identity.
#
# By default, runs convergence tests on grids 17×17, 33×33, 65×65, and 129×129.
# Use --single_solve flag for a single solve on a 65×65 grid instead.
#
# ## Command Line Options:
#
# Grid configuration:
#   --grids <values>         Comma-separated list of grid sizes (e.g., 17,65,257,1025)
#   --N_min <value>          Minimum grid size for auto-generation (default: 17)
#   --N_max <value>          Maximum grid size for auto-generation (default: 129)
#                            Auto-generates: 17, 33, 65, 129, 257, 513, 1025, 2049, ...
#
# Multigrid configuration:
#   --mg_levels_list <vals>  Comma-separated MG levels per grid (e.g., 3,5,7,9)
#   --mg_levels <value>      Fixed number of MG levels for all grids
#                            If neither specified: auto-scales with grid (N=17→3, N=33→4, ...)
#
# Other options:
#   -single_solve           Run single solve instead of convergence test
#
# PETSc options (use single dash):
#   -ksp_type <type>         KSP solver type (preonly, gmres, richardson, etc.)
#   -pc_type <type>          Preconditioner type (lu, mg, sor, ilu, etc.)
#   -ksp_rtol <value>        Relative convergence tolerance
#   -mg_coarse_pc_type <type> Coarse grid preconditioner for multigrid
#   -ksp_monitor             Monitor convergence
#
# ## Usage Examples:
#
# 1. Default convergence analysis (N=17,33,65,129 with LU solver):
#
#     julia --project=.. dmda_laplacian_convergence.jl
#
# 2. Custom grid sizes with geometric multigrid:
#
#     julia --project=.. dmda_laplacian_convergence.jl --grids 17,65,257,1025 -pc_type mg -mg_coarse_pc_type lu
#
# 3. Custom grids with specific MG levels for each:
#
#     julia --project=.. dmda_laplacian_convergence.jl --grids 33,129,513 --mg_levels_list 4,6,8 -pc_type mg -mg_coarse_pc_type lu
#
# 4. Auto-generate grids with fixed MG levels:
#
#     julia --project=.. dmda_laplacian_convergence.jl --N_max 513 --mg_levels 7 -pc_type mg -mg_coarse_pc_type lu
#
# 5. Direct solver with preonly KSP:
#
#     julia --project=.. dmda_laplacian_convergence.jl --N_max 257 -ksp_type preonly -pc_type lu
#
# 6. Iterative SOR solver:
#
#     julia --project=.. dmda_laplacian_convergence.jl --grids 17,65,257,1025 -pc_type sor -ksp_rtol 1e-10
#
# 7. Single solve with multigrid:
#
#     julia --project=.. dmda_laplacian_convergence.jl --single_solve -pc_type mg -ksp_rtol 1e-12 -mg_coarse_pc_type lu
#
# 8. Run in parallel with MPI (4 processes):
#
#     julia> using MPI
#     julia> mpiexec(cmd -> run(`$cmd -n 4 julia --project dmda_laplacian_convergence.jl -pc_type mg`))
#
# Notes:
# - Geometric multigrid (MG) uses the DMDA grid hierarchy for automatic coarsening
# - By default, MG levels automatically scale with grid size (N=17→3, N=33→4, N=65→5, N=129→6)
# - Use --mg_levels to fix the number of levels for all grids
# - PETSc options use single dash (-), Julia script options use double dash (--)

using MPI
using PETSc
using UnicodePlots: heatmap
using ForwardDiff: derivative
using Printf
using Statistics: mean

# Set up the problem by using an exact solution and then using this to set the
# boundary data and forcing  
exact(x, y) = sin(2π * x) * sin(2π * y)
# The discrete operator is -∇², so forcing = -∇²(exact)
# For exact(x,y) = sin(2πx)sin(2πy), we have ∇²u = -2(2π)²sin(2πx)sin(2πy)
# So forcing = -∇²u = 8π²sin(2πx)sin(2πy)
forcing(x, y) = 8 * π^2 * sin(2π * x) * sin(2π * y)

"""
    solve_laplacian(petsclib, comm, N::Int, opts)

Solve the 2D Laplacian problem on an N×N grid and return the L2 and max errors
(computed on interior points only).
"""
function solve_laplacian(petsclib, comm, N::Int, opts; mg_levels=nothing)
    PetscScalar = petsclib.PetscScalar
    Nq = (N, N)
    
    # Add mg_levels to opts if using multigrid and mg_levels is specified
    if mg_levels !== nothing && get(opts, :pc_type, nothing) == "mg"
        opts = merge(opts, (var"-pc_mg_levels" = mg_levels,))
    end
    
    # Add direct coarse grid solver for multigrid if not already specified
    if get(opts, :pc_type, nothing) == "mg" && !haskey(opts, :mg_coarse_pc_type)
        opts = merge(opts, (var"-mg_coarse_pc_type" = "lu",))
    end
    
    # Create the DMDA
    da = PETSc.DMDA(
        petsclib,
        comm,
        (PETSc.DM_BOUNDARY_NONE, PETSc.DM_BOUNDARY_NONE),
        Nq,
        1,
        1,
        PETSc.DMDA_STENCIL_STAR;
        opts...,
    )
    
    # Setup the Krylov solver
    ksp = PETSc.KSP(da; opts...)
    
    # Define the operator assembly function
    function assemble_operator!(A, _, ksp)
        da = PETSc.getDMDA(ksp)
        corners = PETSc.getcorners(da)
        Nq = PETSc.getinfo(da).global_size[1:2]
        
        Δx = PetscScalar(1 / (Nq[1] - 1))
        Δy = PetscScalar(1 / (Nq[2] - 1))
        
        interior = (CartesianIndex(2, 2, 1)):(CartesianIndex(Nq[1] - 1, Nq[2] - 1, 1))
        
        sten = (
            CartesianIndex(-1, 0, 0),
            CartesianIndex(1, 0, 0),
            CartesianIndex(0, -1, 0),
            CartesianIndex(0, 1, 0),
            CartesianIndex(0, 0, 0),
        )
        vals = (-1 / Δx^2, -1 / Δx^2, -1 / Δy^2, -1 / Δy^2, 2 / Δx^2 + 2 / Δy^2)
        
        for i in ((corners.lower):(corners.upper))
            if i ∈ interior
                for (s, v) in zip(sten, vals)
                    A[i, i + s] = v
                end
            else
                A[i, i] = 1
            end
        end
        
        PETSc.assemble!(A)
        return 0
    end
    
    PETSc.setcomputeoperators!(ksp, assemble_operator!)
    
    # Set the right-hand side
    PETSc.setcomputerhs!(ksp) do petsc_b, ksp
        da = PETSc.getDMDA(ksp)
        corners = PETSc.getcorners(da)
        Nq = PETSc.getinfo(da).global_size[1:2]
        
        g_x = range(PetscScalar(0), length = Nq[1], stop = 1)
        g_y = range(PetscScalar(0), length = Nq[2], stop = 1)
        
        l_x = g_x[(corners.lower[1]):(corners.upper[1])]
        l_y = g_y[(corners.lower[2]):(corners.upper[2])]
        
        PETSc.withlocalarray!(petsc_b; read = false) do b
            b = reshape(b, Int64(corners.size[1]), Int64(corners.size[2]))
            b .= forcing.(l_x, l_y')
            
            if corners.lower[1] == 1
                b[1, :] .= exact.(l_x[1], l_y)
            end
            if corners.lower[2] == 1
                b[:, 1] .= exact.(l_x, l_y[1])
            end
            if corners.upper[1] == Nq[1]
                b[end, :] .= exact.(l_x[end], l_y)
            end
            if corners.upper[2] == Nq[2]
                b[:, end] .= exact.(l_x, l_y[end])
            end
        end
        
        return 0
    end
    
    # Solve the problem and time it
    t_start = time()
    PETSc.solve!(ksp)
    t_solve = time() - t_start
    
    # Check convergence and get iteration count
    reason = PETSc.LibPETSc.KSPGetConvergedReason(petsclib, ksp)
    niter = PETSc.LibPETSc.KSPGetIterationNumber(petsclib, ksp)
    if Integer(reason) < 0 && MPI.Comm_rank(comm) == 0
        @warn "KSP diverged for N=$N with reason: $reason"
    end
    
    # Get the solution and compute error
    sol = PETSc.get_solution(ksp)
    corners = PETSc.getcorners(da)
    Nq = PETSc.getinfo(da).global_size[1:2]
    
    g_x = range(PetscScalar(0), length = Nq[1], stop = 1)
    g_y = range(PetscScalar(0), length = Nq[2], stop = 1)
    
    l_x = g_x[(corners.lower[1]):(corners.upper[1])]
    l_y = g_y[(corners.lower[2]):(corners.upper[2])]
    
    Δx = PetscScalar(l_x.step)
    
    x = sol[:]
    u = reshape(x, corners.size[1], corners.size[2])
    
    # Compute error only in interior (excluding boundaries)
    loc_err2_sum = 0.0
    loc_max_err = 0.0
    loc_npts = 0
    
    for i in 1:length(l_x)
        for j in 1:length(l_y)
            gi = corners.lower[1] + i - 1
            gj = corners.lower[2] + j - 1
            
            if gi > 1 && gi < Nq[1] && gj > 1 && gj < Nq[2]
                err = abs(u[i,j] - exact(l_x[i], l_y[j]))
                loc_err2_sum += err^2
                loc_max_err = max(loc_max_err, err)
                loc_npts += 1
            end
        end
    end
    
    root = 0
    err2_sum = MPI.Reduce(loc_err2_sum, +, root, comm)
    npts_total = MPI.Reduce(loc_npts, +, root, comm)
    max_err = MPI.Reduce(loc_max_err, max, root, comm)
    
    L2_error = nothing
    if root == MPI.Comm_rank(comm)
        L2_error = sqrt(err2_sum / npts_total)
    end
    
    # Return u for visualization (only needed for single solve mode)
    return L2_error, max_err, Δx, reason, niter, t_solve, u, l_x, l_y, ksp, da
end

"""
    run_convergence_analysis(petsclib, comm, grid_sizes, opts; mg_levels_list=nothing, mg_levels_override=nothing)

Run a convergence analysis for the specified grid sizes and return errors and mesh spacings.
For multigrid solvers:
- If mg_levels_list is provided, uses the specified level for each grid
- If mg_levels_override is provided, uses that fixed level for all grids
- Otherwise, automatically scales levels with grid size
"""
function run_convergence_analysis(petsclib, comm, grid_sizes, opts; mg_levels_list=nothing, mg_levels_override=nothing)
    L2_errors = Float64[]
    max_errors = Float64[]
    h_values = Float64[]
    iterations = Int[]
    solve_times = Float64[]
    mg_levels_used = Union{Int,Nothing}[]
    
    # Check if using multigrid
    using_mg = get(opts, :pc_type, nothing) == "mg"
    
    for (i, N) in enumerate(grid_sizes)
        # Determine MG levels for this grid
        mg_levels = if mg_levels_list !== nothing
            mg_levels_list[i]
        elseif mg_levels_override !== nothing
            mg_levels_override
        elseif using_mg
            i + 2  # Auto-scale: N=17→3, N=33→4, N=65→5, etc.
        else
            nothing
        end
        
        if using_mg && MPI.Comm_rank(comm) == 0
            println("  Solving N=$N with $mg_levels MG levels...")
        end
        
        # For the first grid, run twice to exclude compilation time
        if i == 1
            # Warmup run (discard results)
            L2_err_warmup, max_err_warmup, h_warmup, reason_warmup, niter_warmup, t_warmup, u_warmup, l_x_warmup, l_y_warmup, ksp_warmup, da_warmup = solve_laplacian(petsclib, comm, N, opts; mg_levels=mg_levels)
            PETSc.destroy(ksp_warmup)
            PETSc.destroy(da_warmup)
        end
        
        # Actual run (record results)
        L2_err, max_err, h, reason, niter, t_solve, u, l_x, l_y, ksp, da = solve_laplacian(petsclib, comm, N, opts; mg_levels=mg_levels)
        
        if MPI.Comm_rank(comm) == 0
            push!(L2_errors, L2_err)
            push!(max_errors, max_err)
            push!(h_values, h)
            push!(iterations, niter)
            push!(solve_times, t_solve)
            push!(mg_levels_used, mg_levels)
        end
        
        # Clean up
        PETSc.destroy(ksp)
        PETSc.destroy(da)
    end
    
    return L2_errors, max_errors, h_values, iterations, solve_times, mg_levels_used
end

"""
    report_convergence_results(grid_sizes, L2_errors, max_errors, h_values, iterations, solve_times, mg_levels_used)

Print a formatted table of convergence results and compute convergence orders.
"""
function report_convergence_results(grid_sizes, L2_errors, max_errors, h_values, iterations, solve_times, mg_levels_used)
    println("\n" * "="^115)
    println("Convergence Analysis Results")
    println("="^115)
    @printf("%-8s %-8s %-12s %-12s %-10s %-12s %-6s %-8s %-12s\n", 
            "N", "MG Lvls", "h", "L2 Error", "L2 Order", "Max Error", "Iters", "Time(s)", "Time Order")
    println("-"^115)
    
    for i in 1:length(grid_sizes)
        N = grid_sizes[i]
        h = h_values[i]
        L2_err = L2_errors[i]
        max_err = max_errors[i]
        niter = iterations[i]
        t_solve = solve_times[i]
        mg_lvl = mg_levels_used[i]
        mg_str = isnothing(mg_lvl) ? "---" : string(mg_lvl)
        
        if i == 1
            @printf("%-8d %-8s %-12.6e %-12.6e %-10s %-12.6e %-6d %-8.4f %-12s\n",
                N, mg_str, h, L2_err, "---", max_err, niter, t_solve, "---")
        else
            L2_order = log(L2_errors[i-1] / L2_errors[i]) / log(h_values[i-1] / h_values[i])
            # Compute time scaling relative to expected complexity for 2D problems
            # Total unknowns = N² (2D grid)
            # Optimal multigrid: O(N²) work, so when N doubles, expect time × 4
            # Sparse direct solver (2D): O(N³) for optimal ordering, so when N doubles, expect time × 8
            # time_order = actual_ratio / expected_O(N²)_ratio
            # For MG: Time Order ≈ 1.0 means optimal O(N²) scaling
            # For sparse direct: Time Order ≈ 2.0 means optimal O(N³) scaling
            N_ratio = grid_sizes[i] / grid_sizes[i-1]
            expected_time_ratio = N_ratio^2  # Always compare to O(N²)
            actual_time_ratio = solve_times[i] / solve_times[i-1]
            time_order = actual_time_ratio / expected_time_ratio
            @printf("%-8d %-8s %-12.6e %-12.6e %-10.6f %-12.6e %-6d %-8.4f %-12.6f\n",
                N, mg_str, h, L2_err, L2_order, max_err, niter, t_solve, time_order)
        end
    end
    
    println("-"^115)
    
    if length(L2_errors) > 1
        # Compute average convergence order
        avg_L2_order = mean([log(L2_errors[i-1] / L2_errors[i]) / log(h_values[i-1] / h_values[i]) 
                             for i in 2:length(L2_errors)])
        
        println("\nSummary:")
        println(@sprintf("  Average L2 convergence order: %.4f", avg_L2_order))
        println("  Expected order for 2nd-order FD: 2.0000")
        
        if avg_L2_order > 1.95 && avg_L2_order < 2.05
            println("  ✓ Convergence test PASSED (within 2.5% of expected)")
        else
            println("  ✗ Convergence test FAILED (not within 2.5% of expected)")
        end
        
        println("\nNote:")
        println("  Time Order shows actual time ratio / expected O(N²) ratio for 2D problem (N×N grid)")
        println("  Optimal multigrid: Time Order ≈ 1.0 (O(N²) = O(total unknowns))")
        println("  Optimal sparse direct (2D): Time Order ≈ 2.0 (O(N³) with nested dissection)")
    end
    println("="^90)
end

# Parse options
opts = if !isinteractive()
    PETSc.parse_options(ARGS)
else
    (ksp_monitor = false, ksp_type = "preonly", pc_type = "lu", ksp_rtol = 1e-12)
end

# Check if single solve is requested (default is convergence test)
run_single_solve = haskey(opts, :single_solve) || (!isinteractive() && "single_solve" in ARGS)

# Parse grid sizes
# Option 1: -grids 17,33,65,129 (comma-separated list)
# Option 2: -N_min 17 -N_max 129 (auto-generate sequence)
grid_sizes = nothing
if haskey(opts, :grids)
    grids_str = opts[:grids]
    grid_sizes = parse.(Int, split(grids_str, ','))
else
    # Default: N=17 to N=129
    N_min = parse(Int, get(opts, :N_min, "17"))
    N_max = parse(Int, get(opts, :N_max, "129"))
end

# Parse MG levels
# Option 1: -mg_levels_list 3,4,5,6 (one per grid, comma-separated)
# Option 2: -mg_levels 5 (fixed for all grids)
# Option 3: auto-scale (default if using MG)
mg_levels_list = nothing
mg_levels_override = nothing
if haskey(opts, :mg_levels_list)
    levels_str = opts[:mg_levels_list]
    mg_levels_list = parse.(Int, split(levels_str, ','))
elseif haskey(opts, :mg_levels)
    mg_levels_override = parse(Int, opts[:mg_levels])
end

# Set our MPI communicator
comm = MPI.COMM_WORLD

# Set our PETSc Scalar Type
PetscScalar = Float64

# get the PETSc lib with our chosen `PetscScalar` type
petsclib = PETSc.petsclibs[1]
#petsclib = PETSc.getlib(; PetscScalar = PetscScalar)

# Initialize PETSc
PETSc.initialize(petsclib)

if !run_single_solve
    # Run convergence analysis
    if MPI.Comm_rank(comm) == 0
        println("\nRunning convergence analysis...")
        solver_type = get(opts, :pc_type, "lu")
        println("Using solver: $solver_type")
        if solver_type == "mg"
            if mg_levels_list !== nothing
                println("  MG levels: $mg_levels_list (per-grid)")
            elseif mg_levels_override !== nothing
                println("  MG levels: $mg_levels_override (fixed)")
            else
                println("  MG levels will scale with grid refinement")
            end
        end
        if haskey(opts, :pc_factor_mat_solver_type)
            println("  Matrix solver: $(opts[:pc_factor_mat_solver_type])")
        end
    end
    
    # Generate or use provided grid sizes
    if grid_sizes === nothing
        # Auto-generate: 17, 33, 65, 129, 257, 513, ...
        # Pattern: 2^n + 1 for n = 4, 5, 6, ...
        grid_sizes = Int[]
        let N = 17
            while N <= N_max
                if N >= N_min
                    push!(grid_sizes, N)
                end
                N = 2*N - 1  # Next value in sequence: 17->33, 33->65, etc.
            end
        end
    end
    
    # Validate mg_levels_list length if provided
    if mg_levels_list !== nothing && length(mg_levels_list) != length(grid_sizes)
        error("Number of MG levels ($(length(mg_levels_list))) must match number of grids ($(length(grid_sizes)))")
    end
    
    if MPI.Comm_rank(comm) == 0
        println("Grid sizes: $grid_sizes")
    end
    
    # Pass mg_levels info to the convergence analysis
    L2_errors, max_errors, h_values, iterations, solve_times, mg_levels_used = run_convergence_analysis(
        petsclib, comm, grid_sizes, opts; mg_levels_list=mg_levels_list, mg_levels_override=mg_levels_override)
    
    if MPI.Comm_rank(comm) == 0
        report_convergence_results(grid_sizes, L2_errors, max_errors, h_values, iterations, solve_times, mg_levels_used)
    end
else
    # Run single solve with visualization
    N = 65
    L2_err, max_err, h, reason, niter, t_solve, u, l_x, l_y, ksp, da = solve_laplacian(petsclib, comm, N, opts)
    
    # Report convergence status
    if MPI.Comm_rank(comm) == 0
        if Integer(reason) > 0
            println("KSP converged with reason: $reason")
        elseif Integer(reason) < 0
            println("WARNING: KSP diverged with reason: $reason")
        else
            println("KSP still iterating: $reason")
        end
        
        println(@sprintf("L2-error (interior): %.6e", L2_err))
        println(@sprintf("Max-error (interior): %.6e", max_err))
        println(@sprintf("Mesh spacing h: %.6e", h))
    end
    
    # If we only have 1 MPI rank we plot the solution and error
    if isinteractive() && MPI.Comm_size(comm) == 1
        display(heatmap(u, zlabel = "sol"))
        display(heatmap(u - exact.(l_x, l_y'), zlabel = "error"))
    end
    
    # Clean up
    PETSc.destroy(ksp)
    PETSc.destroy(da)
end

PETSc.finalize(petsclib)

nothing

