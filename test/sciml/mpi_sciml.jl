# Standalone script: run under `mpiexec -n N julia --project=... mpi_sciml.jl`.
# Exercises distributed (MPI) explicit time integration through the SciML
# extension. The RHS is purely *local* (decoupled componentwise decay), which
# is the milestone-1 contract: each rank owns its block of the global state and
# no cross-rank ghost exchange is required.
using Test
using MPI
MPI.Initialized() || MPI.Init()
using PETSc
using SciMLBase

const comm = MPI.COMM_WORLD
const mpisize = MPI.Comm_size(comm)
const mpirank = MPI.Comm_rank(comm)

@testset "MPI explicit integration (rank $mpirank / $mpisize)" begin
    # Each rank owns `nloc` components of the global state. Decay rates and
    # initial values differ per rank so a wrong data layout would be caught.
    nloc = 3
    k = fill(1.0 + mpirank, nloc)                 # local decay rates
    u0 = Float64[mpirank + 1 + i for i in 1:nloc] # local initial block

    # Purely local RHS: du[i] = -k[i] * u[i]. No neighbour coupling.
    decay!(du, u, p, t) = (@. du = -k * u; nothing)

    tspan = (0.0, 1.0)
    prob = ODEProblem(decay!, u0, tspan)

    # `init` first so we can inspect the distributed solution vector before it
    # is torn down by `solve!`.
    integ = init(prob, TSRK("5dp"); dt = 0.05, comm = comm)

    # The solution vector must be genuinely distributed: its global length is
    # the sum of the local blocks, and this rank owns exactly `nloc` of them.
    @test length(integ.u_petsc) == nloc * mpisize
    rng = PETSc.ownershiprange(integ.u_petsc, false)
    @test length(rng) == nloc

    sol = solve!(integ)
    @test sol.retcode == ReturnCode.Success

    # Each rank holds only its local block of the trajectory.
    @test length(sol.u[end]) == nloc
    uref = u0 .* exp.(-k .* tspan[2])
    @test sol.u[end] ≈ uref atol = 1e-4

    # `saveat` exercises the `TSInterpolate` path, which must duplicate the
    # *distributed* solution vector (a serial `VecSeq` of the local length would
    # mismatch the TS). Save times that fall strictly between steps so the
    # interpolation branch — not just the step-endpoint save — is hit.
    saveat = [0.1, 0.3, 0.7]
    sol_sa = solve(prob, TSRK("5dp"); dt = 0.05, comm = comm, saveat = saveat)
    @test sol_sa.retcode == ReturnCode.Success
    for (i, ts) in enumerate(saveat)
        idx = findfirst(t -> isapprox(t, ts; atol = 1e-10), sol_sa.t)
        @test idx !== nothing
        # Each saved sample is this rank's local block, interpolated correctly.
        # `saveat` values come from `TSInterpolate` (dense output), which is a
        # touch less accurate than the step-endpoint solution, so allow a looser
        # tolerance than the `sol.u[end]` check above.
        @test length(sol_sa.u[idx]) == nloc
        @test sol_sa.u[idx] ≈ u0 .* exp.(-k .* ts) atol = 1e-3
    end
end

# A failing top-level @testset throws at exit, so `success(cmd)` in the parent
# reflects the result. Make a clean run exit 0 explicitly for clarity.
MPI.Barrier(comm)
nothing
