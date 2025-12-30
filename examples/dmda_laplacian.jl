# INCLUDE IN MPI TEST
# # Finite Difference Laplacian
#
# In this example a simple finite difference, 2-D laplacian with zero Dirichlet
# boundary conditions is solved using the PETSc DMDA interface
#
# The approach taken is that the 5-point Laplacian is applied in the interior
# and boundary points are included in the solve, but the forcing for the
# boundary data is the Dirichlet boundary condition and the matrix entry for
# these points is set to the identity.
#
# When run from the commandline this script will parse the commandline arguements:
#
#     julia --project dmda_laplacian.jl -da_refine 4 -ksp_monitor
#
# The script can also be run in parallel:
#
#     julia> using MPI
#     julia> mpiexec(cmd -> run(`$cmd -n 4 julia --project dmda_laplacian.jl -da_refine 4 -ksp_monitor -pc_type asm`))

using MPI
using PETSc
using SparseArrays: spdiagm
using UnicodePlots: heatmap
using ForwardDiff: derivative

opts = if !isinteractive()
    PETSc.parse_options(ARGS)
else
    (ksp_monitor = true, ksp_view = true, pc_type = "mg", pc_mg_levels = 2)
end

# Set our MPI communicator
comm = MPI.COMM_WORLD

# Set our PETSc Scalar Type
PetscScalar = Float64

# get the PETSc lib with our chosen `PetscScalar` type
petsclib = PETSc.getlib(; PetscScalar = PetscScalar)

# Initialize PETSc
PETSc.initialize(petsclib)

# Set the total number of grid points in each direction
Nq = (11, 13)

# Set up the problem by using an exact solution and then using this to set the
# boundary data and forcing
exact(x, y) = cos(π * x) * sin(π * y)
forcing(x, y) =
    derivative(x -> exact(x, y), x) * derivative(x -> exact(x, y), x) +
    derivative(x -> exact(x, y), y) * derivative(x -> exact(x, y), y)

# Create the PETSC dmda object
da = PETSc.DMDA(
    petsclib,
    comm,
    (PETSc.DM_BOUNDARY_NONE, PETSc.DM_BOUNDARY_NONE), # Use no ghost nodes
    Nq,                                               # Global grid size
    1,                                                # Number of DOF per node
    1,                                                # Stencil width
    PETSc.DMDA_STENCIL_STAR;                          # Stencil type
    opts...,
)

# Setup the Krylov solver for the distributed array
ksp = PETSc.KSP(da; opts...)

#
# Set the right-hand side of the ksp
#
PETSc.setcomputerhs!(ksp) do petsc_b, ksp
    # Get the distributed array from the ksp.
    # The reason we need to query the ksp for the `da` is that this ksp may be
    # different than the one we started with (e.g., multi-grid)
    da = PETSc.getDMDA(ksp)

    # Get the corners of the box that this processor is responsible for
    corners = PETSc.getcorners(da)

    # get the global grid dimension
    Nq = PETSc.getinfo(da).global_size[1:2]

    # Determine the global grid in each direction.
    g_x = range(PetscScalar(0), length = Nq[1], stop = 1)
    g_y = range(PetscScalar(0), length = Nq[2], stop = 1)

    # Get the local grid in each direction
    l_x = g_x[(corners.lower[1]):(corners.upper[1])]
    l_y = g_y[(corners.lower[2]):(corners.upper[2])]
   
    # Build the RHS vector for the finite difference grid. Since `petsc_b` is a
    # PETSc vectors, we covert it to a julia vector before using it
    PETSc.withlocalarray!(petsc_b; read = false) do b

        # Reshape the vector into a 2-D vector of the size of the data we own
        b = reshape(b, Int64(corners.size[1]), Int64(corners.size[2]))

        # Set the RHS forcing
        b .= forcing.(l_x, l_y')

        # If on the boundary overwrite the forcing the with the Dirichlet data
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
    
    # 0 is the PETSc success value
    return 0
end

#
# Set the left-hand side operators of the ksp
#
PETSc.setcomputeoperators!(ksp) do A, _, ksp
    # Get the distributed array from the ksp.
    da = PETSc.getDMDA(ksp)

    # Get the corners of the box that this processor is responsible for
    corners = PETSc.getcorners(da)

    # get the global grid dimension
    Nq = PETSc.getinfo(da).global_size[1:2]

    # Get the grid spacing
    Δx = PetscScalar(range(PetscScalar(0), length = Nq[1] + 2, stop = 1).step)
    Δy = PetscScalar(range(PetscScalar(0), length = Nq[2] + 2, stop = 1).step)

    # Interior Points
    interior =
        (CartesianIndex(2, 2, 1)):(CartesianIndex(Nq[1] - 1, Nq[2] - 1, 1))

    # Computational stencil for the interior points
    sten = (
        CartesianIndex(-1, 0, 0),
        CartesianIndex(1, 0, 0),
        CartesianIndex(0, -1, 0),
        CartesianIndex(0, 1, 0),
        CartesianIndex(0, 0, 0),
    )
    vals = (1 / Δx^2, 1 / Δx^2, 1 / Δy^2, 1 / Δy^2, -2 / Δx^2 - 2 / Δy^2)

    # Loop over all the points on the processory
    for i in ((corners.lower):(corners.upper))
        # for the interior points set the stencil otherwise just set to identity
        if i ∈ interior
            for (s, v) in zip(sten, vals)
                A[i, i + s] = v
            end
        else
            A[i, i] = 1
        end
    end

    # Assemble the matrix after inserting values
    PETSc.assemble!(A)

    # 0 is the PETSc success value
    return 0
end

#
# Solve the problem
#
PETSc.solve!(ksp)

#
# Analyze the error in the solution
#

# Get the solution
sol = PETSc.get_solution(ksp)

# Get the corners we own
corners = PETSc.getcorners(da)

# Get the true global size
Nq = PETSc.getinfo(da).global_size[1:2]

# Determine the global grid in each direction.
g_x = range(PetscScalar(0), length = Nq[1], stop = 1)
g_y = range(PetscScalar(0), length = Nq[2], stop = 1)

# Get the local grid in each direction
l_x = g_x[(corners.lower[1]):(corners.upper[1])]
l_y = g_y[(corners.lower[2]):(corners.upper[2])]

# Determine the grid spacing
Δx = PetscScalar(l_x.step)
Δy = PetscScalar(l_y.step)

# Calculate the error
PETSc.withlocalarray!(sol; read = false) do x
    # Change to a 2-D array
    u = reshape(x, corners.size[1], corners.size[2])

    # Compute the local error squared on the mpi rank
    loc_err2 = mapreduce(+, x, l_x, l_y') do u, x, y
        (u - exact(x, y))^2 * Δx * Δy
    end

    # reduce the error to the root rank and display
    root = 0
    ϵ = MPI.Reduce(loc_err2, +, root, comm)
    if root == MPI.Comm_rank(comm)
        println("L2-error is $ϵ")
    end

    # If we only have 1 MPI rank we plot the solution and error
    if isinteractive() && MPI.Comm_size(comm) == 1
        display(heatmap(u, zlabel = "sol"))
        display(heatmap(u - exact.(l_x, l_y'), zlabel = "error"))
    end
end

PETSc.destroy(da)
PETSc.destroy(ksp)
PETSc.finalize(petsclib)

nothing

