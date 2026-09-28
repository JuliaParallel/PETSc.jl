using Test
using PETSc

@testset "lib" begin
    for petsclib in PETSc.petsclibs
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt
        @test petsclib ==
              PETSc.getlib(; PetscScalar = PetscScalar, PetscInt = PetscInt)
        if PetscInt == Int64
            @test petsclib == PETSc.getlib(; PetscScalar = PetscScalar)
        end
        if PetscScalar == Float64
            @test petsclib == PETSc.getlib(; PetscInt = PetscInt)
        end
        if PetscInt == Int64 && PetscScalar == Float64
            @test petsclib == PETSc.getlib()
        end
    end
end
