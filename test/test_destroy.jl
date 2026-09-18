# test/test_destroy.jl
# Destroying PETSc objects must stay safe in the awkward cases: a second
# explicit destroy!, and an object left over from an earlier
# initialize/finalize cycle.
#
# The stale-cycle case is the one that used to crash. `initialize` and
# `finalize` each bump `petsclib.age`, and `PetscFinalize` frees the inner
# communicator. An object created before that still holds a non-null pointer,
# so without an age check `destroy!` called `XXXDestroy` on a dead communicator
# and PETSc aborted the process. Because leaked objects are destroyed by GC
# finalizers, the abort landed at an arbitrary later point, typically inside
# whichever test file happened to trigger a collection.

using Test
using PETSc
using MPI

MPI.Initialized() || MPI.Init()

@testset "destroy!" begin

for petsclib in PETSc.petsclibs
    PetscScalar = PETSc.scalartype(petsclib)
    PetscInt    = PETSc.inttype(petsclib)

    @testset "$(PetscScalar)/$(PetscInt)" begin

        # ── objects that outlive their initialize/finalize cycle ─────────────
        @testset "stale cycle" begin
            PETSc.initialize(petsclib)
            v = PETSc.PetscVec(petsclib, PetscScalar[1, 2, 3, 4])
            m = PETSc.PetscMat(petsclib, 4, 4, 1)
            o = PETSc.PetscOptions(petsclib; ksp_monitor = "")
            age_created = v.age
            PETSc.finalize(petsclib)

            # A new cycle: the library is live again, so a `finalized` check
            # alone would wrongly conclude these are safe to destroy!.
            PETSc.initialize(petsclib)
            @test PETSc.LibPETSc.getlib(typeof(petsclib)).age > age_created
            @test !PETSc.isdestroyable(v, typeof(petsclib))
            @test !PETSc.isdestroyable(m, typeof(petsclib))
            @test !PETSc.isdestroyable(o, typeof(petsclib))

            # Must be a no-op rather than a call into the dead communicator.
            @test PETSc.destroy!(v) === nothing
            @test PETSc.destroy!(m) === nothing
            @test PETSc.destroy!(o) === nothing
            @test v.ptr == C_NULL
            @test m.ptr == C_NULL
            @test o.ptr == C_NULL

            PETSc.finalize(petsclib)
        end

        # ── ordinary lifetime, and destroying twice ──────────────────────────
        @testset "double destroy" begin
            PETSc.initialize(petsclib)

            v = PETSc.PetscVec(petsclib, PetscScalar[1, 2, 3, 4])
            @test PETSc.isdestroyable(v, typeof(petsclib))
            PETSc.destroy!(v)
            @test v.ptr == C_NULL
            @test !PETSc.isdestroyable(v, typeof(petsclib))
            # The finalizer will reach this object again after the explicit
            # destroy!, so a repeat call has to stay harmless.
            @test PETSc.destroy!(v) === nothing

            m = PETSc.PetscMat(petsclib, 4, 4, 1)
            PETSc.destroy!(m)
            @test m.ptr == C_NULL
            @test PETSc.destroy!(m) === nothing

            o = PETSc.PetscOptions(petsclib; ksp_monitor = "")
            @test PETSc.isdestroyable(o, typeof(petsclib))
            PETSc.destroy!(o)
            @test o.ptr == C_NULL
            @test PETSc.destroy!(o) === nothing

            PETSc.finalize(petsclib)
        end

        # ── objects built empty and filled in through an out-parameter ───────
        # DMClone and friends take a pre-allocated object and write the pointer
        # into it, so the object is built by the empty constructor rather than
        # from a pointer. That path has to stamp the age too, or destroy
        # silently skips the object and it leaks.
        @testset "empty constructor stamps age" begin
            PETSc.initialize(petsclib)
            libage = PETSc.LibPETSc.getlib(typeof(petsclib)).age

            for obj in (PETSc.LibPETSc.PetscDM(petsclib),
                        PETSc.LibPETSc.PetscVec(petsclib),
                        PETSc.LibPETSc.PetscMat(petsclib))
                @test obj.age == libage
            end

            PETSc.finalize(petsclib)
        end

        # ── wrappers that borrow a handle instead of owning it ───────────────
        # VecPtr and MatPtr wrap a pointer PETSc still owns, which is how the
        # TS and SNES callbacks hand their arguments to Julia. Destroying one of
        # those frees an object the solver is still using, so `destroy!` consults
        # `own` and leaves the wrapper untouched.
        @testset "borrowed handles" begin
            PETSc.initialize(petsclib)

            v = PETSc.PetscVec(petsclib, PetscScalar[1, 2, 3, 4])
            borrowed_v = PETSc.VecPtr(petsclib, v.ptr, false)
            @test PETSc.owns(v)
            @test !PETSc.owns(borrowed_v)
            @test PETSc.destroy!(borrowed_v) === nothing
            @test borrowed_v.ptr == v.ptr
            @test PETSc.LibPETSc.VecGetSize(petsclib, v) == 4

            m = PETSc.PetscMat(petsclib, 4, 4, 1)
            borrowed_m = PETSc.MatPtr(petsclib, m.ptr, false)
            @test !PETSc.owns(borrowed_m)
            @test PETSc.destroy!(borrowed_m) === nothing
            @test borrowed_m.ptr == m.ptr

            PETSc.destroy!(v)
            PETSc.destroy!(m)
            PETSc.finalize(petsclib)
        end

        # ── borrowed DM handles ──────────────────────────────────────────────
        # `narrow` and every reader returning a DM hand back a second handle
        # onto one PETSc object (docs/src/man/naming.md §3.3, §5.4). Destroying
        # it would invalidate the owner's copy, so `destroy!` consults `own`.
        # A constructor result owns its handle and really is destroyed.
        @testset "borrowed DM handles" begin
            PETSc.initialize(petsclib)
            comm = MPI.COMM_SELF

            da = PETSc.DMDA(
                petsclib, comm, (PETSc.DM_BOUNDARY_NONE,), (8,), 1, 1,
            )
            @test da isa PETSc.DMDA{typeof(petsclib), 1}
            @test PETSc.owns(da)

            borrowed = PETSc.narrow(da)
            @test borrowed isa PETSc.DMDA{typeof(petsclib), 1}
            @test !PETSc.owns(borrowed)
            @test PETSc.destroy!(borrowed) === nothing
            # The no-op leaves both the wrapper and the owner's object usable.
            @test borrowed.ptr == da.ptr
            @test PETSc.ndims(da) == 1

            # A reader hands back a borrowed handle too.
            ksp = PETSc.KSP(da)
            d = PETSc.dm(ksp)
            @test d isa PETSc.DMDA{typeof(petsclib), 1}
            @test !PETSc.owns(d)
            @test PETSc.destroy!(d) === nothing
            @test PETSc.ndims(da) == 1

            # `clone` is the other side: a new object the caller owns.
            c = PETSc.clone(da)
            @test PETSc.owns(c)
            @test c.ptr != da.ptr
            PETSc.destroy!(c)
            @test c.ptr == C_NULL

            PETSc.destroy!(ksp)
            PETSc.destroy!(da)
            @test da.ptr == C_NULL
            PETSc.finalize(petsclib)
        end

        # ── an owning VecPtr still destroys ──────────────────────────────────
        # It records an age like every other wrapper, so `isdestroyable` reads
        # the field instead of throwing on a type that never had one.
        @testset "owned raw pointer" begin
            PETSc.initialize(petsclib)
            libage = PETSc.LibPETSc.getlib(typeof(petsclib)).age

            raw = PETSc.LibPETSc.VecCreateSeq(
                petsclib, PETSc.LibPETSc.PETSC_COMM_SELF, PetscInt(4),
            )
            owned = PETSc.VecPtr(petsclib, raw.ptr, true)
            @test owned.age == libage
            @test PETSc.owns(owned)
            @test PETSc.isdestroyable(owned, typeof(petsclib))

            PETSc.destroy!(owned)
            @test owned.ptr == C_NULL
            @test PETSc.destroy!(owned) === nothing

            raw.ptr = C_NULL  # the VecPtr freed it; keep the finalizer off it
            PETSc.finalize(petsclib)
        end

        # ── after the library is finalized ───────────────────────────────────
        @testset "after finalize" begin
            PETSc.initialize(petsclib)
            v = PETSc.PetscVec(petsclib, PetscScalar[1, 2, 3, 4])
            PETSc.finalize(petsclib)

            @test !PETSc.isdestroyable(v, typeof(petsclib))
            @test PETSc.destroy!(v) === nothing
        end

    end # @testset "$(PetscScalar)/$(PetscInt)"
end # for petsclib

end # @testset "destroy!"
