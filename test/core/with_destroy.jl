# test/core/with_destroy.jl
# with_destroy, track! and untrack!: tracked objects are destroyed when the block
# ends, in reverse order, also after an exception. Runs on any number of ranks.

using Test
using MPI
MPI.Initialized() || MPI.Init()
using PETSc

# Records when `destroy!` reaches it, to check the order
struct Probe
    name::Symbol
    log::Vector{Symbol}
end
PETSc.destroy!(p::Probe) = (push!(p.log, p.name); nothing)

@testset "with_destroy" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    comm = MPI.COMM_WORLD
    da() = PETSc.DMDA(petsclib, comm, (PETSc.DM_BOUNDARY_NONE,), (16,), 1, 1)

    @testset "destroys what it tracks, last tracked first" begin
        log = Symbol[]
        result = PETSc.with_destroy() do list
            for name in (:a, :b, :c)
                @test PETSc.track!(list, Probe(name, log)).name === name
            end
            :result
        end
        @test result === :result
        @test log == [:c, :b, :a]

        objects = PETSc.with_destroy() do list
            dm = PETSc.track!(list, da())
            x = PETSc.track!(list, PETSc.global_vec(dm))
            ksp = PETSc.track!(list, PETSc.KSP(dm))
            (dm, x, ksp)
        end
        @test all(o -> o.ptr == C_NULL, objects)
    end

    @testset "after an exception on every rank" begin
        log = Symbol[]
        dm = nothing
        @test_throws "inside the block" PETSc.with_destroy() do list
            dm = PETSc.track!(list, da())
            PETSc.track!(list, Probe(:probe, log))
            error("inside the block")
        end
        @test dm.ptr == C_NULL
        @test log == [:probe]
        # the ranks are still in step: a collective constructor goes through
        PETSc.destroy!(da())
        MPI.Barrier(comm)
    end

    @testset "borrowed and untracked objects survive" begin
        dm = da()
        ksp = PETSc.KSP(dm)
        x = PETSc.with_destroy() do list
            p = PETSc.track!(list, PETSc.pc(ksp))   # borrowed: not tracked
            @test PETSc.untrack!(list, p) === p
            x = PETSc.track!(list, PETSc.global_vec(dm))
            PETSc.untrack!(list, x)
        end
        @test x.ptr != C_NULL
        @test PETSc.pc(ksp).ptr != C_NULL
        @test PETSc.type_name(PETSc.pc(ksp)) isa Union{Nothing, Symbol}
        foreach(PETSc.destroy!, (x, ksp, dm))
    end

    @testset "misuse is an error" begin
        PETSc.with_destroy() do list
            dm = PETSc.track!(list, da())
            @test_throws "already tracked" PETSc.track!(list, dm)
            @test_throws "has no method for a Int64" PETSc.track!(list, 1)
            @test_throws "not tracked" PETSc.untrack!(list, Probe(:x, Symbol[]))
        end
    end

    PETSc.finalize(petsclib)
end
