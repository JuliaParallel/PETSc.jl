# test/test_destroy.jl
# Destroying PETSc objects must stay safe in the awkward cases: a second
# explicit destroy, and an object left over from an earlier
# initialize/finalize cycle.
#
# The stale-cycle case is the one that used to crash. `initialize` and
# `finalize` each bump `petsclib.age`, and `PetscFinalize` frees the inner
# communicator. An object created before that still holds a non-null pointer,
# so without an age check `destroy` called `XXXDestroy` on a dead communicator
# and PETSc aborted the process. Because leaked objects are destroyed by GC
# finalizers, the abort landed at an arbitrary later point, typically inside
# whichever test file happened to trigger a collection.

using Test
using PETSc
using MPI

if !Sys.iswindows()
    MPI.Initialized() || MPI.Init()
end

@testset "destroy" begin

for petsclib in PETSc.petsclibs
    PetscScalar = PETSc.scalartype(petsclib)
    PetscInt    = PETSc.inttype(petsclib)

    @testset "$(PetscScalar)/$(PetscInt)" begin

        # ── objects that outlive their initialize/finalize cycle ─────────────
        @testset "stale cycle" begin
            PETSc.initialize(petsclib)
            v = PETSc.VecSeq(petsclib, PetscScalar[1, 2, 3, 4])
            m = PETSc.MatSeqAIJ(petsclib, 4, 4, 1)
            age_created = v.age
            PETSc.finalize(petsclib)

            # A new cycle: the library is live again, so a `finalized` check
            # alone would wrongly conclude these are safe to destroy.
            PETSc.initialize(petsclib)
            @test PETSc.LibPETSc.getlib(typeof(petsclib)).age > age_created
            @test !PETSc.isdestroyable(v, typeof(petsclib))
            @test !PETSc.isdestroyable(m, typeof(petsclib))

            # Must be a no-op rather than a call into the dead communicator.
            @test PETSc.destroy(v) === nothing
            @test PETSc.destroy(m) === nothing
            @test v.ptr == C_NULL
            @test m.ptr == C_NULL

            PETSc.finalize(petsclib)
        end

        # ── ordinary lifetime, and destroying twice ──────────────────────────
        @testset "double destroy" begin
            PETSc.initialize(petsclib)

            v = PETSc.VecSeq(petsclib, PetscScalar[1, 2, 3, 4])
            @test PETSc.isdestroyable(v, typeof(petsclib))
            PETSc.destroy(v)
            @test v.ptr == C_NULL
            @test !PETSc.isdestroyable(v, typeof(petsclib))
            # The finalizer will reach this object again after the explicit
            # destroy, so a repeat call has to stay harmless.
            @test PETSc.destroy(v) === nothing

            m = PETSc.MatSeqAIJ(petsclib, 4, 4, 1)
            PETSc.destroy(m)
            @test m.ptr == C_NULL
            @test PETSc.destroy(m) === nothing

            PETSc.finalize(petsclib)
        end

        # ── after the library is finalized ───────────────────────────────────
        @testset "after finalize" begin
            PETSc.initialize(petsclib)
            v = PETSc.VecSeq(petsclib, PetscScalar[1, 2, 3, 4])
            PETSc.finalize(petsclib)

            @test !PETSc.isdestroyable(v, typeof(petsclib))
            @test PETSc.destroy(v) === nothing
        end

    end # @testset "$(PetscScalar)/$(PetscInt)"
end # for petsclib

end # @testset "destroy"
