# test/core/finalizer_queue.jl
# A handle the garbage collector drops is destroyed by the next PETSc.jl call that
# creates or destroys an object, not by the finalizer itself: the finalizer runs on
# whichever thread triggered the collection, and the shipped PETSc is not thread safe.

using Test
using PETSc
using MPI

MPI.Initialized() || MPI.Init()

# Reference count of the object behind `ptr` after the Julia handle was dropped and
# the collector ran. A second reference is held, so the object is still alive and
# the count says how many references remain.
@noinline function dropped_vec_pointer(petsclib)
    v = PETSc.PetscVec(petsclib, zeros(petsclib.PetscScalar, 4))
    PETSc.LibPETSc.PetscObjectReference(petsclib, v)
    return v.ptr
end
references(petsclib, ptr) =
    PETSc.LibPETSc.PetscObjectGetReference(petsclib, PETSc.PetscVec(petsclib, ptr; own = false))

@testset "the finalizer queues, a PETSc.jl call destroys" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    PETSc.drain_destroy_queue!()

    @testset "nothing is destroyed inside the finalizer" begin
        ptr = dropped_vec_pointer(petsclib)
        GC.gc(); GC.gc()
        # the finalizer ran and queued the handle; the object has not been destroyed
        @test references(petsclib, ptr) == 2
        PETSc.drain_destroy_queue!()
        @test references(petsclib, ptr) == 1
        PETSc.destroy!(PETSc.PetscVec(petsclib, ptr; own = true))
    end

    @testset "a constructor destroys what is queued" begin
        ptr = dropped_vec_pointer(petsclib)
        GC.gc(); GC.gc()
        @test references(petsclib, ptr) == 2
        PETSc.destroy!(PETSc.PetscVec(petsclib, zeros(petsclib.PetscScalar, 2)))
        @test references(petsclib, ptr) == 1
        PETSc.destroy!(PETSc.PetscVec(petsclib, ptr; own = true))
    end

    @testset "finalize destroys what is queued" begin
        ptr = dropped_vec_pointer(petsclib)
        GC.gc(); GC.gc()
        @test PETSc.destroy_queue_length[] >= 1
        PETSc.finalize(petsclib)
        @test PETSc.destroy_queue_length[] == 0
        @test isempty(PETSc.destroy_queue)
    end
end

# The threaded case, in a process of its own: it needs several threads and PETSc's
# debug malloc, which reports a destroy that overlaps another PETSc call
@testset "collections on another thread during solves" begin
    if Sys.CPU_THREADS >= 4
        script = joinpath(@__DIR__, "finalizer_stress.jl")
        julia = `$(Base.julia_cmd()) --project=$(Base.active_project()) -t 4 $script 15`
        cmd = addenv(julia, "PETSC_OPTIONS" => "-malloc_debug -malloc_test")
        log = tempname()
        process = run(pipeline(ignorestatus(cmd); stdout = log, stderr = log); wait = false)
        # a destroy that overlaps another PETSc call can also hang the process
        finished = timedwait(() -> !process_running(process), 180) === :ok
        finished || kill(process)
        output = read(log, String)
        rm(log)
        ok = finished && success(process) && !occursin("PETSC ERROR", output)
        @test ok
        ok || print(output)
    else
        @info "skipped: fewer than 4 CPU threads"
    end
end
