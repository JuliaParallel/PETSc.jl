# INCLUDE IN MPI TEST
#=
In this example we solve the [Liouville–Bratu–Gelfand
equation](https://en.wikipedia.org/wiki/Liouville%E2%80%93Bratu%E2%80%93Gelfand_equation):
```math
∇² ψ + λ exp(ψ)
```
with zero Dirichlet boundary conditions.

To solve the problem we use the standard central, finite difference
approximation of the Laplacian in `d`-dimensions

This example is motivated by the following PETSc examples:
- [`snes/ex5`](https://petsc.org/release/src/snes/tutorials/ex5.c.html)
- [`snes/ex14`](https://petsc.org/release/src/snes/tutorials/ex14.c.html)
=#

using MPI
using PETSc
using OffsetArrays: OffsetArray
using LinearAlgebra: norm

opts = if !isinteractive()
    PETSc.parse_options(ARGS)
else
    (ksp_monitor = true, ksp_view = true)
end

# Set our MPI communicator
comm = MPI.COMM_WORLD

# Set our PETSc Scalar Type
PetscScalar = Float64

# get the PETSc lib with our chosen `PetscScalar` type
petsclib = PETSc.getlib(; PetscScalar = PetscScalar)

# Initialize PETSc
PETSc.initialize(petsclib)

# dimensionality of the problem
dim = haskey(opts, :dim) ? opts.dim : 3

# Set the total number of grid points in each direction
Nq = ntuple(_ -> 10, dim)

# Set the boundary conditions on each side
bcs = ntuple(_ -> PETSc.DM_BOUNDARY_NONE, dim)

# Set parameter
λ = PetscScalar(6)

# Create the PETSC dmda object
da = PETSc.DMDA(
    petsclib,
    comm,
    bcs,                     # boundary conditions
    Nq,                      # Global grid size
    1,                       # Number of DOF per node
    1,                       # Stencil width
    PETSc.DMDA_STENCIL_STAR; # Stencil type
    opts...,
)

# Create the PETSC snes object
snes = PETSc.SNES(petsclib, comm; opts...)

# add the da to the snes
PETSc.setDM!(snes, da)

# Set up the initial guess
x = PETSc.DMGlobalVec(da)
xl = PETSc.DMLocalVec(da)
PETSc.withlocalarray!(xl; read = false) do l_x
    corners = PETSc.getcorners(da)

    # Get the global grid dimensions
    Nq = PETSc.getinfo(da).global_size

    # Figure out the interior points 
    int_min = min(CartesianIndex(corners.size), CartesianIndex(2, 2, 2))
    int_max = max(CartesianIndex(corners.size .- 1), CartesianIndex(1, 1, 1))
    interior = (int_min):(int_max)

    # Allows us to adress the local array with global indexing
    ox = @view PETSc.reshapelocalarray(l_x, da)[1, :, :, :]

    # Set up the global coordinates in each direction
    # -1 to 1 when Nq > 1 and 0 otherwise
    coords = map(
        Nq ->
            Nq == 1 ? range(PetscScalar(0), stop = 0, length = Nq) :
            range(-PetscScalar(1), stop = 1, length = Nq),
        Nq,
    )

    scaling = λ / (λ + 1)

    # Loop over all the points on the processor and set the initial condition to
    # be a hat function
    for i in ((corners.lower):(corners.upper))
        ox[i] =
            scaling * sqrt(
                min(
                    (1 - abs(coords[1][i[1]])),
                    (1 - abs(coords[2][i[2]])),
                    (1 - abs(coords[3][i[3]])),
                ),
            )
    end
end
PETSc.update!(x, xl, PETSc.INSERT_VALUES)    # update local -> global

# Set up the nonlinear function
r = similar(x)
PETSc.setfunction!(snes, r) do g_fx, snes, g_x
    # Get the DMDA associated with the snes
    da = PETSc.getDMDA(snes)

    # Get a local vector and transfer the data from the global vector into it
    l_x = PETSc.DMLocalVec(da)
    PETSc.update!(l_x, g_x, PETSc.INSERT_VALUES)

    ghostcorners = PETSc.getghostcorners(da)
    corners = PETSc.getcorners(da)

    # Global grid size
    Nq = PETSc.getinfo(da).global_size

    # grid spacing in each dimension
    Δx, Δy, Δz = PetscScalar(1) ./ Nq

    # Get local arrays
    PETSc.withlocalarray!(
        (g_fx, l_x);
        read = (false, true),
        write = (true, false),
    ) do fx, x

        # reshape the array and allow for global indexing
        x = @view PETSc.reshapelocalarray(x, da)[1, :, :, :]
        fx = @view PETSc.reshapelocalarray(fx, da)[1, :, :, :]

        # Store a tuple of stencils in each direction
        stencils = (
            (
                CartesianIndex(-1, 0, 0),
                CartesianIndex(0, 0, 0),
                CartesianIndex(1, 0, 0),
            ),
            (
                CartesianIndex(0, -1, 0),
                CartesianIndex(0, 0, 0),
                CartesianIndex(0, 1, 0),
            ),
            (
                CartesianIndex(0, 0, -1),
                CartesianIndex(0, 0, 0),
                CartesianIndex(0, 0, 1),
            ),
        )
        # Weights for each direction
        weights = (Δy * Δz / Δx, Δx * Δz / Δy, Δx * Δy / Δz)

        # loop over indices and set the function value
        for ind in ((corners.lower):(corners.upper))
            # If on the boundary just set equal to the incoming data
            # otherwise apply the finite difference operator
            if any(ntuple(j -> ind[j] == 1 || ind[j] == Nq[j], dim))
                fx[ind] = x[ind]
            else
                # Apply the source
                u = -Δx * Δy * Δz * λ * exp(x[ind])

                # Apply the finite diffference stencil
                for (s, w) in zip(stencils[1:dim], weights[1:dim])
                    u +=
                        w * (-x[ind + s[1]] + 2 * x[ind + s[2]] - x[ind + s[3]])
                end
                fx[ind] = u
            end
        end
    end

    # Clean up the local vector
    PETSc.destroy(l_x)
    return 0
end

J = PETSc.MatAIJ(da)
PETSc.setjacobian!(snes, J) do J, snes, g_x
    # Get the DMDA associated with the snes
    da = PETSc.getDMDA(snes)

    # Get the corners of the points we own
    corners = PETSc.getcorners(da)

    # Global grid size
    Nq = PETSc.getinfo(da).global_size

    # grid spacing in each dimension
    Δx, Δy, Δz = PetscScalar(1) ./ Nq

    # Store a tuple of stencils in each direction
    stencils = (
        (
            CartesianIndex(-1, 0, 0),
            CartesianIndex(0, 0, 0),
            CartesianIndex(1, 0, 0),
        ),
        (
            CartesianIndex(0, -1, 0),
            CartesianIndex(0, 0, 0),
            CartesianIndex(0, 1, 0),
        ),
        (
            CartesianIndex(0, 0, -1),
            CartesianIndex(0, 0, 0),
            CartesianIndex(0, 0, 1),
        ),
    )
    # Weights for each direction
    weights = (Δy * Δz / Δx, Δx * Δz / Δy, Δx * Δy / Δz)

    # Get a local array of the solution vector
    PETSc.withlocalarray!(g_x; write = false) do l_x
        # reshape so we can use multi-D indexing
        x = @view PETSc.reshapelocalarray(l_x, da)[1, :, :, :]

        # loop over indices and set the function value
        for ind in ((corners.lower):(corners.upper))
            # If on the boundary just set equal to the incoming data
            # otherwise apply the finite difference operator
            if any(ntuple(j -> ind[j] == 1 || ind[j] == Nq[j], dim))
                J[ind, ind] = 1
            else
                # We accumulate the diagonal and add it at the end
                Jii = -Δx * Δy * Δz * λ * exp(x[ind]) # Apply the source

                # Apply the finite diffference stencil
                for (s, w) in zip(stencils[1:dim], weights[1:dim])
                    Jii += w * 2
                    J[ind, ind + s[1]] = -w
                    J[ind, ind + s[3]] = -w
                end
                J[ind, ind] = Jii
            end
        end
    end

    # Assemble the Jacobian matrix
    PETSc.assemble!(J)
    return 0
end

if MPI.Comm_rank(comm) == 0
    println(@elapsed(PETSc.solve!(x, snes)))
else
    PETSc.solve!(x, snes)
end
g = similar(x)
snes.f!(g, snes, x)
nm = norm(g)
if MPI.Comm_rank(comm) == 0
    @show nm
end

# Do some clean up
PETSc.destroy(J)
PETSc.destroy(x)
PETSc.destroy(g)
PETSc.destroy(r)
PETSc.destroy(da)
PETSc.destroy(snes)

PETSc.finalize(petsclib)
