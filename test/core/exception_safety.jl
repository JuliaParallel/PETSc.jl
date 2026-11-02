# test/core/exception_safety.jl
# A throw half way through a checkout or a constructor hands back what was
# already taken: no Vec stays checked out, and a solver that failed to set up
# drops its reference on the matrix or DM it was given.

using Test
using PETSc
using MPI

MPI.Initialized() || MPI.Init()

# A backend that records the Vecs it hands back
struct Recorder <: PETSc.AbstractPetscMemBackend end
const released = Any[]
PETSc.release_local_array(cpu_arr, ::Recorder, vec::PETSc.LibPETSc.AbstractPetscVec; read, write) =
    push!(released, vec)

@testset "exception safety" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    comm = MPI.COMM_SELF
    LibPETSc = PETSc.LibPETSc
    PetscScalar = petsclib.PetscScalar
    references(obj) = LibPETSc.PetscObjectGetReference(petsclib, obj)

    @testset "a failing second checkout hands back the first" begin
        # This PETSc build tracks no checkouts, so a getter that fails on the
        # second Vec and a backend that records each hand-back stand in for it
        acquire(v; read, write) = v === w ? error("second checkout fails") : (nothing, nothing, Recorder())
        v = PETSc.PetscVec(petsclib, PetscScalar[1, 2])
        w = PETSc.PetscVec(petsclib, PetscScalar[3, 4])
        @test_throws "second checkout fails" PETSc.checkout_arrays(
            (_, _) -> nothing, nothing, (v, w), (true, true), (false, false), acquire)
        @test released == [v]
        foreach(PETSc.destroy!, (v, w))
    end

    @testset "a solver that fails to set up is destroyed" begin
        A = PETSc.PetscMat(petsclib, 4, 4, 1)
        before = references(A)
        @test_throws "options_at must be" PETSc.KSP(A; options_at = :never, ksp_type = "cg")
        @test references(A) == before

        dm = PETSc.DMDA(petsclib, comm, (PETSc.DM_BOUNDARY_NONE,), (8,), 1; stencil_width = 1)
        before = references(dm)
        @test_throws "options_at must be" PETSc.KSP(dm; options_at = :never, ksp_type = "cg")
        @test references(dm) == before

        @test_throws "options_at must be" PETSc.SNES(petsclib, comm; options_at = :never, snes_type = "newtonls")
        @test_throws "options_at must be" PETSc.TS(petsclib, comm; options_at = :never, ts_type = "rk")
        foreach(PETSc.destroy!, (dm, A))
    end

    PETSc.finalize(petsclib)
end
