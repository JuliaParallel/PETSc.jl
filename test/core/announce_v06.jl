using Test
using MPI
MPI.Initialized() || MPI.Init()
using PETSc
using PETSc: LibPETSc

# Calls whose behavior changes in v0.6 warn once and still return what they did.

@testset "v0.6 announcements" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    PetscScalar = petsclib.PetscScalar
    comm = LibPETSc.PETSC_COMM_SELF

    v = PETSc.PetscVec(petsclib, PetscScalar[1, 2, 3])
    m = PETSc.PetscMat(petsclib, 3, 3, 1)

    @testset "VecPtr and MatPtr" begin
        vp = @test_logs (:warn, r"VecPtr goes in v0.6.*PetscVec\(petsclib, ptr; own\)") PETSc.VecPtr(petsclib, v.ptr)
        @test vp isa PETSc.VecPtr
        @test vp.ptr == v.ptr && !PETSc.owns(vp)
        mp = @test_logs (:warn, r"MatPtr goes in v0.6.*PetscMat\(petsclib, ptr; own\)") PETSc.MatPtr(petsclib, m.ptr)
        @test mp isa PETSc.MatPtr
        @test mp.ptr == m.ptr && !PETSc.owns(mp)

        # The replacements give the same wrappers and do not warn
        wv = @test_logs PETSc.PetscVec(petsclib, v.ptr; own = false)
        @test wv isa PETSc.VecPtr && wv.ptr == v.ptr && !PETSc.owns(wv)
        wm = @test_logs PETSc.PetscMat(petsclib, m.ptr; own = false)
        @test wm isa PETSc.MatPtr && wm.ptr == m.ptr && !PETSc.owns(wm)
        @test PETSc.destroy!(wm) === nothing
        @test LibPETSc.MatGetSize(petsclib, m) == (3, 3)
    end

    foreach(PETSc.destroy!, (v, m))
    PETSc.finalize(petsclib)
end
