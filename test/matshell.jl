using Test
using PETSc
using MPI

@testset "MatShell" begin
    for petsclib in PETSc.petsclibs
        petsclib = PETSc.petsclibs[1]
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt

        local_rows = 10
        local_cols = 5
        function f!(p_x, p_y)
            PETSc.withlocalarray!((p_x, p_y), write=(true,false)) do x, y
                x .= [2y; 3y]
            end
        end

        matshell =  PETSc.MatShell(petsclib, f!, MPI.COMM_SELF, local_rows, local_cols)
        x = PetscScalar.(collect(1:5))
        petsc_x = LibPETSc.VecCreateSeqWithArray(petsclib,MPI.COMM_SELF, PetscInt(1), PetscInt(local_cols), x)
        petsc_y = LibPETSc.VecCreateSeqWithArray(petsclib,MPI.COMM_SELF, PetscInt(1), PetscInt(local_rows), zeros(PetscScalar, local_rows))
        
        LibPETSc.MatMult(petsclib, matshell, petsc_x, petsc_y)
        @test petsc_y[:] == [2x; 3x]
        @test matshell * x == [2x; 3x]

        PETSc.destroy(matshell)
        PETSc.finalize(petsclib)
    end
end
