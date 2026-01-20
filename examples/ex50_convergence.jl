# This example demonstrates solving a 2D Poisson equation with Neumann boundary
# conditions using a DMDA and KSP solver, and performs convergence analysis.
# 
# Based on PETSc's ex50.c:
# https://petsc.org/main/src/ksp/ksp/tutorials/ex50.c.html
#
# Poisson equation: div(grad p) = f on [0,1]^2
# with Neumann boundary conditions: dp/dx = 0 at x=0,1 and dp/dy = 0 at y=0,1
# Forcing function: f = -cos(π*x)*cos(π*y)
#
# ## Usage Examples:
#
# 1. Run single solve with default settings (100×100 grid, multigrid):
#    julia --project=. examples/ex50_convergence.jl
#
# 2. Set grid size with -N option (creates N×N grid):
#    julia --project=. examples/ex50_convergence.jl -N 50 -ksp_monitor
#
# 3. Run convergence analysis with custom resolutions:
#    julia --project=. examples/ex50_convergence.jl -convergence_analysis -resolutions 9,17,33
#
# ## Convergence Test Examples:
#
# 1. Direct solver (LU factorization):
#    julia --project=. examples/ex50_convergence.jl -convergence_analysis -pc_type lu -ksp_type preonly
#
# 2. AMG solver (algebraic multigrid):
#    julia --project=. examples/ex50_convergence.jl -convergence_analysis -pc_type gamg -ksp_type cg
#
# 3. Geometric multigrid:
#    julia --project=. examples/ex50_convergence.jl -convergence_analysis -pc_type mg -ksp_type cg
#
# 4. Geometric Galerkin multigrid:
#    julia --project=. examples/ex50_convergence.jl -convergence_analysis -pc_type mg -pc_mg_galerkin -ksp_type cg
#
# 3. Monitor convergence with custom grid:
#    julia --project=. examples/ex50.jl -da_grid_x 50 -da_grid_y 50 -ksp_monitor
#
# 4. Use direct LU solver:
#    julia --project=. examples/ex50.jl -pc_type lu -ksp_monitor -ksp_converged_reason
#
# 5. Multigrid with refinement:
#    julia --project=. examples/ex50.jl -da_grid_x 3 -da_grid_y 3 -da_refine 5 -pc_type mg -ksp_monitor
#
# 6. Programmatic usage with custom solver options:
#    julia -e "include(\"examples/ex50.jl\"); solve_poisson(50, 0; pc_mg_levels=3, ksp_monitor=true)"
#
# 7. Run convergence analysis for different resolutions:
#    julia --project=. examples/ex50.jl -convergence_analysis
#
# 8. Run convergence analysis with custom MG levels cap:
#    julia --project=. examples/ex50_convergence.jl -convergence_analysis -resolutions 33,65,129,257,513,1025,2049 -pc_type mg -ksp_type cg -mg_levels_cap 9

using PETSc, MPI, Printf

# Get the PETSc lib with our chosen `PetscScalar` type
petsclib = PETSc.getlib(; PetscScalar = Float64)

# Initialize PETSc 
PETSc.initialize(petsclib)


"""
    analytical_solution(x, y)

Compute the analytical solution for the Poisson equation with forcing function
f = -cos(π*x)*cos(π*y).

The solution is p(x,y) = (1)/(2*π²) * cos(π*x)*cos(π*y) + C,
where C is an arbitrary constant due to Neumann boundary conditions.
Here, C=0.
"""
function analytical_solution(x, y)
    K = 1.0 / (2 * π^2)
    return -K * cos(π * (x - 0.5)) * cos(π * (y - 0.5))
end

"""
    solve_poisson(N=100, da_refine=0; solver_opts...)

Solve the 2D Poisson equation in serial or parallel with Neumann boundary conditions using PETSc.

# Arguments
- `N::Int`: Grid size (N×N grid). Default is 100.
- `da_refine::Int`: Refinement level for the DMDA. Default is 0.
- `solver_opts`: Additional keyword options passed to the PETSc solver (e.g., `pc_mg_levels=3`, `ksp_monitor=true`).

# Returns
- `solve_time::Float64`: Time taken to solve the system.
- `niter::Int`: Number of iterations taken by the solver.
- `final_grid_size::Tuple{Int, Int}`: The final grid size after refinement.
- `sol2D::Array`: The 2D solution array.
- `l2_error::Float64`: The L2 norm of the error compared to analytical solution.

# Examples
```julia
# Solve with default settings
time, iters, grid = solve_poisson()

# Solve with custom grid size and multigrid levels
time, iters, grid = solve_poisson(50, 0; pc_mg_levels=3, ksp_monitor=true)

# Command-line usage
# julia examples/ex50.jl -N 50 -pc_mg_levels 3
```
"""
function solve_poisson(N=100, da_refine=0; solver_opts...)

    # Set our PETSc Scalar Type
    PetscScalar = Float64

    # Get MPI communicator
    comm = MPI.COMM_WORLD

    # Parse command-line options
    opts = PETSc.parse_options(ARGS)

    # Use parameter values as defaults, but allow options to override
    N = parse(Int, get(opts, Symbol("N"), string(N)))
    da_refine = parse(Int, get(opts, Symbol("da_refine"), string(da_refine)))
    da_grid_x = parse(Int, get(opts, Symbol("da_grid_x"), string(N)))
    da_grid_y = parse(Int, get(opts, Symbol("da_grid_y"), string(N)))

    # Use CG solver with AMG preconditioner as defaults
    ksp_type = get(opts, Symbol("ksp_type"), "cg")
    pc_type = get(opts, Symbol("pc_type"), "gamg")
    ksp_rtol = parse(Float64, get(opts, Symbol("ksp_rtol"), "1e-12"))

    # Set the grid options
    opts = merge(opts, (da_grid_x = da_grid_x, da_grid_y = da_grid_y, da_refine = da_refine))

    # Set solver options
    opts = merge(opts, (ksp_type = ksp_type, pc_type = pc_type, ksp_rtol = ksp_rtol))

    # Merge in any additional solver options passed as keyword arguments
    opts = merge(opts, Dict(solver_opts))

    boundary_type = (PETSc.DM_BOUNDARY_NONE, PETSc.DM_BOUNDARY_NONE)
    stencil_type = PETSc.DMDA_STENCIL_STAR
    global_size = (11, 11)
    procs = (PETSc.PETSC_DECIDE, PETSc.PETSC_DECIDE)
    dof_per_node = 1
    stencil_width = 1

    da = PETSc.DMDA(
        petsclib,
        comm,
        boundary_type,
        global_size,
        dof_per_node,
        stencil_width,
        stencil_type;
        opts...,
    )

    ksp = PETSc.KSP(da; opts...)

    # Print the final grid size after any refinement (only on rank 0)
    final_grid_size = PETSc.getinfo(da).global_size[1:2]
    if MPI.Comm_rank(comm) == 0
        @printf("Solving on %d × %d grid\n", final_grid_size[1], final_grid_size[2])
    end

    # Set the Jacobian/operator
    PETSc.setcomputeoperators!(ksp) do J, jac, ksp
        dm = PETSc.getDMDA(ksp)
        corners = PETSc.getcorners(dm)
        global_size = PETSc.getinfo(dm).global_size[1:2]

        # Grid spacing in each direction
        h = PetscScalar(1) ./ global_size
        HxdHy = h[1] / h[2]
        HydHx = h[2] / h[1]
        
        # Assemble the matrix - Neumann BCs with variable stencil at boundaries
        for j in corners.lower[2]:corners.upper[2]
            for i in corners.lower[1]:corners.upper[1]
                idx = CartesianIndex(i, j, 1)
                
                # Check if we're on boundary
                is_boundary = (i == 1 || j == 1 || i == global_size[1] || j == global_size[2])
                
                if is_boundary
                    # NOTE: we use a 1th order accurate stencil at boundaries for Neumann BCs
                    # This mimics the C-example but limits the overall accuracy of the scheme to first order.

                    # Boundary point - variable stencil (only interior neighbors)
                    diag_val = PetscScalar(0)
                    numi = 0  # count of i-direction neighbors
                    numj = 0  # count of j-direction neighbors
                    
                    if j > 1  # not on bottom boundary
                        jac[idx, idx + CartesianIndex(0, -1, 0)] = -HxdHy
                        numj += 1
                    end
                    if i > 1  # not on left boundary
                        jac[idx, idx + CartesianIndex(-1, 0, 0)] = -HydHx
                        numi += 1
                    end
                    if i < global_size[1]  # not on right boundary
                        jac[idx, idx + CartesianIndex(1, 0, 0)] = -HydHx
                        numi += 1
                    end
                    if j < global_size[2]  # not on top boundary
                        jac[idx, idx + CartesianIndex(0, 1, 0)] = -HxdHy
                        numj += 1
                    end
                    
                    # Diagonal: sum of neighbor coefficients
                    diag_val = numj * HxdHy + numi * HydHx
                    jac[idx, idx] = diag_val
                else
                    # Interior point - full 5-point stencil
                    jac[idx, idx + CartesianIndex(0, -1, 0)] = -HxdHy
                    jac[idx, idx + CartesianIndex(-1, 0, 0)] = -HydHx
                    jac[idx, idx] = 2.0 * (HxdHy + HydHx)
                    jac[idx, idx + CartesianIndex(1, 0, 0)] = -HydHx
                    jac[idx, idx + CartesianIndex(0, 1, 0)] = -HxdHy
                end
            end
        end

        PETSc.assemble!(jac)
        
        # Set a constant nullspace on the matrix for Neumann BCs
        nullspace = LibPETSc.MatNullSpaceCreate(petsclib, MPI.COMM_WORLD, LibPETSc.PETSC_TRUE, 0, LibPETSc.PetscVec[])
        LibPETSc.MatSetNullSpace(petsclib, jac, nullspace)
        LibPETSc.MatNullSpaceDestroy(petsclib, nullspace)
        # Don't destroy nullspace here - let the matrix manage it
        
        return 0
    end

    # Set the right-hand side
    PETSc.setcomputerhs!(ksp) do b_vec, ksp
        dm = PETSc.getDMDA(ksp)
        corners = PETSc.getcorners(dm)
        global_size = PETSc.getinfo(dm).global_size[1:2]
        
        # Grid spacing in each direction
        h = PetscScalar(1) ./ global_size

        # Build the RHS vector with forcing function
        PETSc.withlocalarray!(b_vec; read = false) do b
            b = reshape(b, Int64(corners.size[1]), Int64(corners.size[2]))
            
            for (iy, y) in enumerate(corners.lower[2]:corners.upper[2])
                # Cell-centered y-coordinate
                y_coord = (y - 0.5) * h[2]
                
                for (ix, x) in enumerate(corners.lower[1]:corners.upper[1])
                    # Cell-centered x-coordinate
                    x_coord = (x - 0.5) * h[1]
                    # Forcing function: -cos(π*x)*cos(π*y) * Hx * Hy
                    b[ix, iy] = -cos(π * x_coord) * cos(π * y_coord) * h[1] * h[2]
                end
            end
        end
        
        # Assemble the vector before modifying it
        PETSc.assemble!(b_vec)
        
        # Remove the nullspace from the RHS to make it consistent for the singular matrix
        LibPETSc.VecAssemblyBegin(petsclib, b_vec)
        LibPETSc.VecAssemblyEnd(petsclib, b_vec)
        nullspace = LibPETSc.MatNullSpaceCreate(petsclib, MPI.COMM_WORLD, LibPETSc.PETSC_TRUE, 0, LibPETSc.PetscVec[])
        LibPETSc.MatNullSpaceRemove(petsclib, nullspace, b_vec)
        LibPETSc.MatNullSpaceDestroy(petsclib, nullspace)

        return 0
    end

    # Solve the problem
    solve_time = @elapsed PETSc.solve!(ksp)

    # Get iteration count and print results (only on rank 0)
    niter = LibPETSc.KSPGetIterationNumber(petsclib, ksp)
    if MPI.Comm_rank(comm) == 0
        @printf("KSP converged in %d iterations in %.4f seconds\n", niter, solve_time)
    end

    # Compute L2 norm of the error and extrema
    dm = PETSc.getDMDA(ksp)
    corners = PETSc.getcorners(dm)
    global_size = PETSc.getinfo(dm).global_size[1:2]
    h = PetscScalar(1) ./ global_size

    sol = PETSc.get_solution(ksp)
    
    # Get local solution array for analysis
    sol2D = nothing
    PETSc.withlocalarray!(sol; read=true) do s
        sol2D = copy(s)
        sol2D = reshape(sol2D, Int64(corners.size[1]), Int64(corners.size[2]))
    end
    
    l2_error = 0.0
    PETSc.withlocalarray!(sol; read=true) do s
        nx, ny = corners.size[1:2]
        s2D = reshape(s, Int64(corners.size[1]), Int64(corners.size[2]))
        
        x_coords = range(-0.5, 0.5, length=global_size[1])
        y_coords = range(-0.5, 0.5, length=global_size[2])
        x_coord = x_coords[corners.lower[1]:corners.upper[1]]
        y_coord = y_coords[corners.lower[2]:corners.upper[2]]
        
        l2_local = 0.0
        sum_diff_local = 0.0
        count_local = 0
        for iy=1:ny, ix=1:nx
            p_ana = analytical_solution(x_coord[ix], y_coord[iy])
            p_num = s2D[ix, iy]
            diff = p_num - p_ana
            sum_diff_local += diff
            count_local += 1
        end
        # Compute global average difference to remove constant shift
        sum_diff_global = MPI.Allreduce(sum_diff_local, MPI.SUM, comm)
        count_global = MPI.Allreduce(count_local, MPI.SUM, comm)
        avg_diff = sum_diff_global / count_global
        # Now compute L2 error with constant removed
        l2_local = 0.0
        for iy=1:ny, ix=1:nx
            p_ana = analytical_solution(x_coord[ix], y_coord[iy])
            p_num = s2D[ix, iy]
            diff = (p_num - p_ana) - avg_diff
            l2_local += diff^2 * h[1] * h[2]
        end
        l2_global = MPI.Allreduce(l2_local, MPI.SUM, comm)
        l2_error = sqrt(l2_global)
        
        if MPI.Comm_rank(comm) == 0
            @printf("L2 norm of error: %.6e\n", l2_error)
        end
    end

    # Clean up
    # Note: When using KSP(da), the KSP takes ownership of the DM reference
    # so we only need to destroy the KSP, not the DA
    PETSc.destroy(ksp)

    #PETSc.finalize(petsclib)

    return solve_time, niter, final_grid_size, sol2D, l2_error
end

"""
    convergence_test(; solver_opts...)

Run convergence test for different grid resolutions and compute order of convergence.

# Arguments
- `solver_opts`: Additional keyword options passed to the PETSc solver.

# Returns
- `results::Vector{Tuple{Int, Float64, Float64, Float64}}`: Vector of (N, h, l2_error, order) for each resolution.
"""
function convergence_test(; solver_opts...)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    
    # Grid sizes for convergence test
    Ns = [10, 20, 40, 80, 160]
    results = []
    
    for (i, N) in enumerate(Ns)
        if rank == 0
            @printf("Running convergence test for N = %d\n", N)
        end
        
        solve_time, niter, grid_size, sol2D, l2_error = solve_poisson(N, 0; solver_opts...)
        
        h = 1.0 / N
        push!(results, (N, h, l2_error, 0.0))  # order will be computed later
    end
    
    # Compute orders
    for i in 1:length(results)-1
        N1, h1, e1, _ = results[i]
        N2, h2, e2, _ = results[i+1]
        if e1 > 0 && e2 > 0
            order = log(e1 / e2) / log(h1 / h2)
            results[i] = (N1, h1, e1, order)
        end
    end
    results[end] = (results[end][1], results[end][2], results[end][3], NaN)  # last one has no order
    
    # Print table
    if rank == 0
        @printf("\nConvergence Test Results:\n")
        @printf("%-5s %-10s %-12s %-8s\n", "N", "h", "L2 Error", "Order")
        @printf("%-5s %-10s %-12s %-8s\n", "-"^5, "-"^10, "-"^12, "-"^8)
        for (N, h, err, order) in results
            if isnan(order)
                @printf("%-5d %-10.6f %-12.6e %-8s\n", N, h, err, "N/A")
            else
                @printf("%-5d %-10.6f %-12.6e %-8.2f\n", N, h, err, order)
            end
        end
    end
    
    return results
end

"""
    convergence_analysis(resolutions=[10, 20, 40, 80, 160], base_mg_levels=1, mg_levels_cap=7; solver_opts...)

Run a convergence analysis by solving the Poisson equation on different grid resolutions
and computing the order of convergence.

# Arguments
- `resolutions::Vector{Int}`: List of grid sizes (N×N) to test. Default is [10, 20, 40, 80, 160].
- `base_mg_levels::Int`: Number of MG levels at the coarsest grid. Default is 1.
- `mg_levels_cap::Int`: Maximum number of MG levels to use. Default is 7.
- `solver_opts`: Additional keyword options passed to the PETSc solver.

# Returns
- `results::Vector{Tuple{Int, Float64, Float64, Float64, Float64, String, Float64, Int}}`: Vector of (N, h, l2_error, order, solve_time, levels_str, time_order, niter) for each resolution.
"""
function convergence_analysis(resolutions=[10, 20, 40, 80, 160], base_mg_levels=1, mg_levels_cap=7; solver_opts...)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    
    min_N = minimum(resolutions)
    pc_type = get(solver_opts, :pc_type, "gamg")
    ksp_type = get(solver_opts, :ksp_type, "cg")
    results = []
    
    for (i, N) in enumerate(resolutions)
        # Compute MG levels: start with base_mg_levels for first resolution, increase by 1 for each subsequent
        levels = base_mg_levels + (i - 1)
        levels_str = if pc_type in ["mg", "gamg"] string(levels) else "--" end
        opts_dict = Dict(solver_opts)
        if pc_type == "mg"
            # Use proper MG configuration for optimal multigrid behavior
            opts_dict[:pc_mg_levels] = string(levels)
            opts_dict[:mg_coarse_ksp_type] = "preonly"
            opts_dict[:mg_coarse_pc_type] = "lu"
            opts_dict[:mg_levels_ksp_type] = "chebyshev"
            opts_dict[:mg_levels_pc_type] = "jacobi"
            opts_dict[:pc_mg_smoothup] = "2"
            opts_dict[:pc_mg_smoothdown] = "2"
        end
        
        if rank == 0
            @printf("Running convergence analysis for N = %d (MG levels = %s)\n", N, levels_str)
        end
        
        # Run twice for the first grid to exclude compilation time
        if i == 1
            solve_poisson(N, 0; opts_dict...)  # discard first run
        end
        solve_time, niter, grid_size, sol2D, l2_error = solve_poisson(N, 0; opts_dict...)
        
        h = 1.0 / N
        push!(results, (N, levels_str, h, niter, l2_error, 0.0, solve_time, 0.0))  # order and time_order will be computed later
    end
    
    # Compute orders
    for i in 1:length(results)-1
        N1, levels_str1, h1, niter1, e1, _, t1, _ = results[i]
        N2, levels_str2, h2, niter2, e2, _, t2, _ = results[i+1]
        if e1 > 0 && e2 > 0 && h1 > 0 && h2 > 0
            order = log(e1 / e2) / log(h1 / h2)
        else
            order = NaN
        end
        if t1 > 0 && t2 > 0 && N1 > 0 && N2 > 0
            time_order = log(t2 / t1) / (2 * log(N2 / N1))  # Divide by 2 for 2D scaling
        else
            time_order = NaN
        end
        results[i] = (N1, levels_str1, h1, niter1, e1, order, t1, time_order)
    end
    if length(results) > 0
        results[end] = (results[end][1], results[end][2], results[end][3], results[end][4], results[end][5], NaN, results[end][7], NaN)  # last ones have no order
    end
    
    # Print table
    if rank == 0
        @printf("\nConvergence Analysis Results (KSP: %s, PC: %s):\n", ksp_type, pc_type)
        @printf("%-5s %-8s %-10s %-10s %-12s %-8s %-10s %-8s\n", "N", "MG Lvl", "h", "KSP Iters", "L2 Error", "𝒪(N)", "Time (s)", "𝒪(Time)")
        @printf("%-5s %-8s %-10s %-10s %-12s %-8s %-10s %-8s\n", "-"^5, "-"^8, "-"^10, "-"^10, "-"^12, "-"^8, "-"^10, "-"^8)
        for (N, levels_str, h, niter, err, order, time, time_order) in results
            if isnan(order) && isnan(time_order)
                @printf("%-5d %-8s %-10.6f %-10d %-12.6e %-8s %-10.4f %-8s\n", N, levels_str, h, niter, err, "N/A", time, "N/A")
            elseif isnan(order)
                @printf("%-5d %-8s %-10.6f %-10d %-12.6e %-8s %-10.4f %-8.2f\n", N, levels_str, h, niter, err, "N/A", time, time_order)
            elseif isnan(time_order)
                @printf("%-5d %-8s %-10.6f %-10d %-12.6e %-8.2f %-10.4f %-8s\n", N, levels_str, h, niter, err, order, time, "N/A")
            else
                @printf("%-5d %-8s %-10.6f %-10d %-12.6e %-8.2f %-10.4f %-8.2f\n", N, levels_str, h, niter, err, order, time, time_order)
            end
        end
    end
    
    return results
end

# If run as script, call the function with defaults
if !isinteractive() && abspath(PROGRAM_FILE) == @__FILE__
    if "-convergence_analysis" in ARGS
        # Remove the flag from ARGS
        args = filter(x -> x != "-convergence_analysis", ARGS)
        opts = PETSc.parse_options(args)
        
        # Parse resolutions from -resolutions option
        resolutions_str = get(opts, Symbol("resolutions"), "10,20,40,80,160")
        resolutions = parse.(Int, split(resolutions_str, ","))
        
        # Parse base_mg_levels from -base_mg_levels option
        base_mg_levels_str = get(opts, Symbol("base_mg_levels"), "1")
        base_mg_levels = parse(Int, base_mg_levels_str)
        
        # Parse mg_levels_cap from -mg_levels_cap option
        mg_levels_cap_str = get(opts, Symbol("mg_levels_cap"), "7")
        mg_levels_cap = parse(Int, mg_levels_cap_str)
        
        convergence_analysis(resolutions, base_mg_levels, mg_levels_cap; opts...)
    elseif "-convergence_test" in ARGS
        # Remove the flag from ARGS
        args = filter(x -> x != "-convergence_test", ARGS)
        opts = PETSc.parse_options(args)
        convergence_test(; opts...)
    else
        solve_poisson()
    end
else
    solve_time, niter, final_grid_size, sol2D, l2_error = solve_poisson();
end

