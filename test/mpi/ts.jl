# test/mpi/ts.jl
# A TS on a distributed vector: forward Euler on du/dt = -u, checked against the
# discrete solution, and a right-hand side that throws on every rank inside
# `with_destroy`. Runs on any number of ranks.

using Test
using MPI
MPI.Initialized() || MPI.Init()
using PETSc

isdefined(Main, :PETScTestUtils) || include(joinpath(@__DIR__, "..", "testutils.jl"))
using .PETScTestUtils: without_petsc_traceback

@testset "TS on $(MPI.Comm_size(MPI.COMM_WORLD)) ranks" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    comm = MPI.COMM_WORLD
    PetscScalar = petsclib.PetscScalar
    N = 4 * MPI.Comm_size(comm)
    dt, steps = 0.1, 10
    da() = PETSc.DMDA(petsclib, comm, (PETSc.DM_BOUNDARY_NONE,), (N,), 1; stencil_width = 1)

    # u(0) = i at global index i, so every rank holds different values
    function initial!(u)
        rng = PETSc.ownership_range(u)
        PETSc.with_local_array!(a -> (a .= rng), u; read = false)
        return u
    end
    function decay!(F, x)
        PETSc.with_local_array!((x, F); read = (true, false), write = (false, true)) do xa, Fa
            Fa .= .-xa
        end
    end
    function explicit_ts(dm, rhs!)
        ts = PETSc.TS(petsclib, comm)
        PETSc.set_type!(ts, :euler)
        PETSc.set_adapt_type!(ts, :none)
        PETSc.set_dm!(ts, dm)
        PETSc.set_rhs_function!(rhs!, ts)
        PETSc.set_time!(ts, 0.0)
        PETSc.set_timestep!(ts, dt)
        PETSc.set_max_time!(ts, steps * dt)
        return ts
    end

    @testset "forward Euler matches the discrete solution" begin
        PETSc.with_destroy() do list
            dm = PETSc.track!(list, da())
            u = PETSc.track!(list, initial!(PETSc.global_vec(dm)))
            ts = PETSc.track!(list, explicit_ts(dm, (F, _ts, _t, x) -> decay!(F, x)))
            PETSc.solve!(u, ts)
            @test PETSc.step_number(ts) == steps
            # u_n = (1 - dt)^n u_0 for forward Euler, entry by entry
            rng = PETSc.ownership_range(u)
            got = PETSc.with_local_array!(copy, u; write = false)
            @test got ≈ (1 - dt)^steps .* rng
        end
    end

    @testset "a right-hand side that throws on every rank" begin
        objects = Any[]
        @test_throws "rhs failed" without_petsc_traceback(petsclib) do
            PETSc.with_destroy() do list
                dm = PETSc.track!(list, da())
                u = PETSc.track!(list, initial!(PETSc.global_vec(dm)))
                ts = PETSc.track!(list, explicit_ts(dm, (F, _ts, _t, x) -> error("rhs failed")))
                append!(objects, (dm, u, ts))
                PETSc.solve!(u, ts)
            end
        end
        @test length(objects) == 3
        @test all(o -> o.ptr == C_NULL, objects)
        # the ranks are still in step: a collective constructor and a solve go through
        PETSc.with_destroy() do list
            dm = PETSc.track!(list, da())
            u = PETSc.track!(list, initial!(PETSc.global_vec(dm)))
            ts = PETSc.track!(list, explicit_ts(dm, (F, _ts, _t, x) -> decay!(F, x)))
            PETSc.solve!(u, ts)
            @test PETSc.step_number(ts) == steps
        end
        MPI.Barrier(comm)
    end

    PETSc.finalize(petsclib)
end
