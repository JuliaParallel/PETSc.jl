# Newton's method for u'' + u² = f on [0, 1], u(0) = 0, u(1) = 1, after PETSc's
# src/snes/tutorials/ex2.c. Serial.
#
# The Jacobian is coded by hand, or with -ad_jacobian computed by ForwardDiff from the
# residual. The ForwardDiff Jacobian is dense, which is fine for 21 points; see
# examples/dmstag/dmstag_stokes.jl for one that uses the stencil.
#
#   julia --project=examples examples/basics/snes_ex2.jl
#   julia --project=examples examples/basics/snes_ex2.jl -ad_jacobian

using PETSc, MPI, LinearAlgebra, SparseArrays, ForwardDiff, UnicodePlots

MPI.Initialized() || MPI.Init()
petsclib = first(PETSc.petsclibs)
PETSc.isinitialized(petsclib) || PETSc.initialize(petsclib)
comm = MPI.COMM_WORLD

const n = 21
const dx = 1 / (n - 1)
const source = [6x + (x + 1e-12)^6 for x in range(0, 1; length = n)]

# The residual: the boundary values, and u'' + u² − f inside
function residual!(f, u)
    f[1] = u[1] - 0
    for i in 2:(n - 1)
        f[i] = (u[i - 1] - 2u[i] + u[i + 1]) / dx^2 + u[i]^2 - source[i]
    end
    f[n] = u[n] - 1
    return f
end

function hand_jacobian!(J, u)
    J[1, 1] = 1
    for i in 2:(n - 1)
        J[i, i - 1] = 1 / dx^2
        J[i, i] = -2 / dx^2 + 2u[i]
        J[i, i + 1] = 1 / dx^2
    end
    J[n, n] = 1
    return J
end

ad_jacobian!(J, u) = copyto!(J, sparse(ForwardDiff.jacobian(u -> residual!(similar(u), u), u)))

jacobian! = haskey(PETSc.parse_options(ARGS), :ad_jacobian) ? ad_jacobian! : hand_jacobian!

# The matrix gets the tridiagonal pattern of the hand-coded Jacobian
pattern = sparse(hand_jacobian!(zeros(n, n), zeros(n)) .!= 0)
J = PETSc.PetscMat(petsclib, comm, Float64.(pattern))

u = PETSc.PetscVec(petsclib, fill(0.5, n))      # the initial guess, then the solution
r = PETSc.PetscVec(petsclib, zeros(n))

snes = PETSc.SNES(petsclib, comm; snes_rtol = 1e-12, snes_monitor = true)
PETSc.set_function!((r, _, u) -> (residual!(r, u); nothing), snes, r)
PETSc.set_snes_jacobian!(snes, J) do J, _, u
    jacobian!(J, u[:])          # u[:] copies the iterate into a Vector
    PETSc.assemble!(J)
end
PETSc.solve!(u, snes)

solution = u[:]
@show norm(residual!(zeros(n), solution))
foreach(PETSc.destroy!, (snes, J, r, u))

lineplot(range(0, 1; length = n), solution; xlabel = "x", ylabel = "u")
