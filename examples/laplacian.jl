# # Finite Difference Laplacian
#
# In this example a simple finite difference, 1-D laplacian with Dirichlet
# boundary conditions is solved first constructing the operator as Julia sparse
# matrix then solving using PETsc

using PETSc
using SparseArrays: spdiagm
using UnicodePlots: lineplot, lineplot!

# Initialize PETSc
PETSc.initialize()

# Set our PETSc Scalar Type
PetscScalar = Float64

# Set the total number of grid points
Nq = 100

# number of interior points
n = Nq - 2

# create the interior grid and get the grid spacing
x = range(PetscScalar(0), length = Nq, stop = 1)[2:(end - 1)]
Δx = PetscScalar(x.step)

# Create the finite difference operator with zero Dirichlet boundary conditions
A =
    spdiagm(
        -1 => -ones(PetscScalar, n - 1),
        0 => 2ones(PetscScalar, n),
        1 => -ones(PetscScalar, n - 1),
    ) / Δx^2

# Setup the exact solution and the forcing function
κ = 2PetscScalar(π)
u(x) = sin(κ * x)
f(x) = κ^2 * sin(κ * x)

# Set up the PETSc solver
ksp = PETSc.KSP(A; ksp_monitor = true);

# Solve the problem
v = ksp \ f.(x);

# Compute the L2-error
ϵ = sqrt(sum((v - u.(x)) .^ 2 * Δx))
println("L2-error is $ϵ")

# Plot the solution and error
display(lineplot(x, v, xlabel = "x", ylabel = "solution"))
display(lineplot(x, v - u.(x), xlabel = "x", ylabel = "error"))
