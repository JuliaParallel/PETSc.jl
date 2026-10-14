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

# `corners`, `info` and `reshape_local_array` are dimension-correct
# (docs/src/man/naming.md §12): on a `dim`-dimensional DMDA they answer with
# `dim`-tuples and `CartesianIndex{dim}`, and the reshaped local array has
# `1 + dim` axes. These three helpers are what the example needs as a result.

# The single-dof slice of a reshaped local array.
dof_slice(A) = view(A, 1, ntuple(_ -> Colon(), dim)...)

# The unit offset along each axis, for the finite-difference stencil.
const units = ntuple(
    j -> CartesianIndex(ntuple(k -> k == j ? 1 : 0, dim)),
    dim,
)

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
PETSc.set_dm!(snes, da)

# Set up the initial guess
x = PETSc.global_vec(da)
xl = PETSc.local_vec(da)
PETSc.with_local_array!(xl; read = false) do l_x
    corners = PETSc.corners(da)

    # Get the global grid dimensions
    Nq = PETSc.info(da).global_size

    # Allows us to adress the local array with global indexing
    ox = dof_slice(PETSc.reshape_local_array(l_x, da))

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
            scaling *
            sqrt(minimum(ntuple(j -> 1 - abs(coords[j][i[j]]), dim)))
    end
end
PETSc.local_to_global!(x, da, xl, PETSc.INSERT_VALUES)

# Set up the nonlinear function
r = similar(x)
PETSc.set_function!(snes, r) do g_fx, snes, g_x
    # Get the DMDA associated with the snes
    da = PETSc.dm(snes)

    # Get a local vector and transfer the data from the global vector into it
    l_x = PETSc.local_vec(da)
    PETSc.global_to_local!(l_x, da, g_x, PETSc.INSERT_VALUES)

    ghostcorners = PETSc.ghost_corners(da)
    corners = PETSc.corners(da)

    # Global grid size
    Nq = PETSc.info(da).global_size

    # grid spacing in each dimension
    Δ = PetscScalar(1) ./ Nq

    # Get local arrays
    PETSc.with_local_array!(
        (g_fx, l_x);
        read = (false, true),
        write = (true, false),
    ) do fx, x

        # reshape the array and allow for global indexing
        x = dof_slice(PETSc.reshape_local_array(x, da))
        fx = dof_slice(PETSc.reshape_local_array(fx, da))

        # Weights for each direction
        weights = ntuple(j -> prod(Δ) / Δ[j]^2, dim)

        # loop over indices and set the function value
        for ind in ((corners.lower):(corners.upper))
            # If on the boundary just set equal to the incoming data
            # otherwise apply the finite difference operator
            if any(ntuple(j -> ind[j] == 1 || ind[j] == Nq[j], dim))
                fx[ind] = x[ind]
            else
                # Apply the source
                u = -prod(Δ) * λ * exp(x[ind])

                # Apply the finite diffference stencil
                for (e, w) in zip(units, weights)
                    u += w * (-x[ind - e] + 2 * x[ind] - x[ind + e])
                end
                fx[ind] = u
            end
        end
    end

    # Clean up the local vector
    PETSc.destroy!(l_x)
    return 0
end

J = LibPETSc.DMCreateMatrix(petsclib, da)
PETSc.set_snes_jacobian!(snes, J) do J, snes, g_x
    # Get the DMDA associated with the snes
    da = PETSc.dm(snes)

    # Get the corners of the points we own
    corners = PETSc.corners(da)

    # Global grid size
    Nq = PETSc.info(da).global_size

    # grid spacing in each dimension
    Δ = PetscScalar(1) ./ Nq

    # Weights for each direction
    weights = ntuple(j -> prod(Δ) / Δ[j]^2, dim)

    # Get a local array of the solution vector
    PETSc.with_local_array!(g_x; write = false) do l_x
        # reshape so we can use multi-D indexing
        x = dof_slice(PETSc.reshape_local_array(l_x, da))

        # loop over indices and set the function value
        for ind in ((corners.lower):(corners.upper))
            # If on the boundary just set equal to the incoming data
            # otherwise apply the finite difference operator
            if any(ntuple(j -> ind[j] == 1 || ind[j] == Nq[j], dim))
                J[ind, ind] = 1
            else
                # We accumulate the diagonal and add it at the end
                Jii = -prod(Δ) * λ * exp(x[ind]) # Apply the source

                # Apply the finite diffference stencil
                for (e, w) in zip(units, weights)
                    Jii += w * 2
                    J[ind, ind - e] = -w
                    J[ind, ind + e] = -w
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
PETSc.destroy!(J)
PETSc.destroy!(x)
PETSc.destroy!(g)
PETSc.destroy!(r)
PETSc.destroy!(da)
PETSc.destroy!(snes)

PETSc.finalize(petsclib)

