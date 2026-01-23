# INCLUDE IN MPI TEST
#
# This example demonstrates solving a 2D Poisson equation with Neumann boundary
# conditions using a DMDA and KSP solver.
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
# 1. Run with default settings (100×100 grid, LU direct solver):
#    julia --project=. examples/ex50.jl
#
# 2. Set grid size with -N option (creates N×N grid):
#    julia --project=. examples/ex50.jl -N 50 -ksp_monitor
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

    # Use CG solver with AMG preconditioner
    ksp_type = get(opts, Symbol("ksp_type"), "cg")
    pc_type = get(opts, Symbol("pc_type"), "gamg")
    ksp_rtol = parse(Float64, get(opts, Symbol("ksp_rtol"), "1e-12"))

    # Boundary condition type: "neumann" (default) or "dirichlet"
    bc_type = get(opts, Symbol("bc_type"), "neumann")

    # Set the grid options
    opts = merge(opts, (da_grid_x = da_grid_x, da_grid_y = da_grid_y, da_refine = da_refine))

    # Set solver options
    opts = merge(opts, (ksp_type = ksp_type, pc_type = pc_type, ksp_rtol = ksp_rtol))

    # Merge in any additional solver options passed as keyword arguments
    opts = merge(opts, Dict(solver_opts))

    boundary_type = (PETSc.DM_BOUNDARY_NONE, PETSc.DM_BOUNDARY_NONE)
    stencil_type = PETSc.DMDA_STENCIL_STAR
    global_size = (da_grid_x, da_grid_y)
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
        
        for j in corners.lower[2]:corners.upper[2]
            for i in corners.lower[1]:corners.upper[1]
                idx = CartesianIndex(i, j, 1)
                
                # Check if we're on boundary
                is_boundary = (i == 1 || j == 1 || i == global_size[1] || j == global_size[2])
                
                if bc_type == "dirichlet" && is_boundary
                    # Dirichlet boundary point: enforce u = analytic value -> set row to identity
                    jac[idx, idx] = 1.0
                elseif is_boundary
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
        comm = PETSc.getcomm(ksp)
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
                    is_boundary = (x == 1 || y == 1 || x == global_size[1] || y == global_size[2])

                    if bc_type == "dirichlet" && is_boundary
                        # Dirichlet boundary value evaluated at boundary coordinate
                        xb = x == 1 ? 0.0 : (x == global_size[1] ? 1.0 : x_coord)
                        yb = y == 1 ? 0.0 : (y == global_size[2] ? 1.0 : y_coord)
                        b[ix, iy] = analytical_solution(xb, yb)
                    else
                        # Forcing function: -cos(π*x)*cos(π*y) * Hx * Hy
                        b[ix, iy] = -cos(π * x_coord) * cos(π * y_coord) * h[1] * h[2]
                    end
                end
            end
        end
        
        # Assemble the vector before modifying it
        PETSc.assemble!(b_vec)
        
        # Remove the nullspace from the RHS only for Neumann BCs (singular matrix)
        if bc_type == "neumann"
            LibPETSc.VecAssemblyBegin(petsclib, b_vec)
            LibPETSc.VecAssemblyEnd(petsclib, b_vec)
            nullspace = LibPETSc.MatNullSpaceCreate(petsclib, MPI.COMM_WORLD, LibPETSc.PETSC_TRUE, 0, LibPETSc.PetscVec[])
            LibPETSc.MatNullSpaceRemove(petsclib, nullspace, b_vec)
            LibPETSc.MatNullSpaceDestroy(petsclib, nullspace)
        end

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

        if bc_type == "neumann"
            # Compute global average difference to remove constant shift
            sum_diff_global = MPI.Allreduce(sum_diff_local, MPI.SUM, comm)
            count_global = MPI.Allreduce(count_local, MPI.SUM, comm)
            avg_diff = sum_diff_global / count_global
        else
            avg_diff = 0.0
        end

        # Now compute L2 error (remove avg only for Neumann)

    solve_poisson(N, da_refine; pc_type=pc_type, ksp_type=ksp_type, ksp_rtol=ksp_rtol)
end

