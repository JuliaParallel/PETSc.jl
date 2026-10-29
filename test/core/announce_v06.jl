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

    @testset "isapprox on two PetscVecs" begin
        w = PETSc.PetscVec(petsclib, PetscScalar[1, 2, 3 + 1e-12])
        @test @test_logs (:warn, r"isapprox on two PetscVecs.*v0.6") isapprox(v, w)
        # Entry by entry, 1e-9 is not ≈ 0, though the norm of the difference is
        # within rtol of the norms, so Julia arrays compare as approximately equal
        x = PETSc.PetscVec(petsclib, PetscScalar[1e-9, 1])
        y = PETSc.PetscVec(petsclib, PetscScalar[0, 1])
        @test !(@test_logs (:warn, r"isapprox on two PetscVecs") isapprox(x, y))
        @test isapprox(x[:], y[:])
        foreach(PETSc.destroy!, (w, x, y))
    end

    foreach(PETSc.destroy!, (v, m))
    PETSc.finalize(petsclib)
end
