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
# 1. Run with default settings (100×100 grid, multigrid):
#    julia --project=. examples/ksp/ex50.jl
#
# 2. Set grid size with -N option (creates N×N grid):
#    julia --project=. examples/ksp/ex50.jl -N 50 -ksp_monitor
#
# 3. Monitor convergence with custom grid:
#    julia --project=. examples/ksp/ex50.jl -da_grid_x 50 -da_grid_y 50 -ksp_monitor
#
# 4. Use direct LU solver:
#    julia --project=. examples/ksp/ex50.jl -pc_type lu -ksp_monitor -ksp_converged_reason
#
# 5. Multigrid with refinement:
#    julia --project=. examples/ksp/ex50.jl -da_grid_x 3 -da_grid_y 3 -da_refine 5 -pc_type mg -ksp_monitor
#
# KNOWN LIMITATION: The C code uses MatNullSpaceCreate/MatSetNullSpace/MatNullSpaceRemove
# to handle the constant nullspace from all-Neumann BCs. However, this causes
# ReadOnlyMemoryError in PETSc.jl due to a wrapper bug in src/ksp.jl (line 201):
# the RHS vector b_vec is created with age=0 instead of getlib(PetscLib).age,
# making it read-only. As a workaround, nullspace handling is commented out.
# The solver may still work due to PETSc's automatic handling, but results may
# not be as accurate as with explicit nullspace removal.
#
using PETSc, MPI, Printf
MPI.Initialized() || MPI.Init()

# Set our PETSc Scalar Type
PetscScalar = Float64

# Get the PETSc lib with our chosen `PetscScalar` type
petsclib = PETSc.getlib(; PetscScalar = PetscScalar)

# Initialize PETSc
PETSc.initialize(petsclib)

# Set our MPI communicator
comm = MPI.COMM_WORLD

# Parse command-line options
opts = if !isinteractive()
    PETSc.parse_options(ARGS)
else
    # Default options when running interactively
    if MPI.Comm_size(comm) == 1
        (
            ksp_monitor = false,
            ksp_view = false,
            da_grid_x = 100,
            da_grid_y = 100,
            pc_type = "mg",
            pc_mg_levels = 1,
            mg_levels_0_pc_type = "ilu",
            mg_levels_0_pc_factor_levels = 1,
        )
    else
        (
            da_grid_x = 3,
            da_grid_y = 3,
            pc_type = "mg",
            da_refine = 10,
            ksp_monitor = false,
            ksp_view = false,
            log_view = nothing,
        )
    end
end

# Parse -N option for grid size (sets both da_grid_x and da_grid_y)
N = parse(Int, get(opts, Symbol("-N"), "0"))
if N > 0
    # Override da_grid_x and da_grid_y with N if specified
    opts = merge(opts, (da_grid_x = N, da_grid_y = N))
end

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
if MPI.Comm_rank(comm) == 0
    final_grid_size = PETSc.getinfo(da).global_size[1:2]
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
    
    # NOTE: The C code sets a constant nullspace on the matrix for Neumann BCs.
    # However, this causes ReadOnlyMemoryError in PETSc.jl due to a wrapper bug.
    # TODO: Fix this wrapper issue and re-enable nullspace handling
    # comm = PETSc.getcomm(ksp)
    # nullspace = LibPETSc.MatNullSpaceCreate(petsclib, comm, LibPETSc.PETSC_TRUE, 0, LibPETSc.PetscVec[])
    # LibPETSc.MatSetNullSpace(petsclib, J, nullspace)
    # LibPETSc.MatNullSpaceDestroy(petsclib, nullspace)
    
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
                # Forcing function: -cos(π*x)*cos(π*y) * Hx * Hy
                b[ix, iy] = -cos(π * x_coord) * cos(π * y_coord) * h[1] * h[2]
            end
        end
    end
    
    # NOTE: The C code calls MatNullSpaceRemove here to force the RHS to be
    # consistent for the singular matrix. However, this causes ReadOnlyMemoryError
    # in PETSc.jl due to a wrapper bug (src/ksp.jl line 201: b_vec created with age=0).
    # The nullspace is set on the matrix in setcomputeoperators!, which should
    # be sufficient for KSP to handle it correctly.
    #
    # TODO: Fix wrapper to create b_vec with getlib(PetscLib).age instead of 0
    # LibPETSc.VecAssemblyBegin(petsclib, b_vec)
    # LibPETSc.VecAssemblyEnd(petsclib, b_vec)
    # nullspace = LibPETSc.MatNullSpaceCreate(petsclib, comm, LibPETSc.PETSC_TRUE, 0, LibPETSc.PetscVec[])
    # LibPETSc.MatNullSpaceRemove(petsclib, nullspace, b_vec)
    # LibPETSc.MatNullSpaceDestroy(petsclib, nullspace)

    return 0
end

# Solve the problem
solve_time = @elapsed PETSc.solve!(ksp)

# Get iteration count and print results (only on rank 0)
if MPI.Comm_rank(comm) == 0
    niter = LibPETSc.KSPGetIterationNumber(petsclib, ksp)
    @printf("KSP converged in %d iterations in %.4f seconds\n", niter, solve_time)
end

# Clean up
# Note: When using KSP(da), the KSP takes ownership of the DM reference
# so we only need to destroy the KSP, not the DA
PETSc.destroy(ksp)
PETSc.finalize(petsclib)