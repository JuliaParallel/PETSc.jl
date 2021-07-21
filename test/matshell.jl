using Test
using PETSc
using MPI

@testset "MatShell" begin
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar

        local_rows = 10
        local_cols = 5
        function f!(p_x, p_y)
            PETSc.withlocalarray!((p_x, p_y)) do x, y
                x .= [2y; 3y]
            end
        end
        x_jl = collect

        matshell =
            PETSc.MatShell(petsclib, f!, MPI.COMM_SELF, local_rows, local_cols)
        x = PetscScalar.(collect(1:5))
        @test matshell * x == [2x; 3x]

        PETSc.destroy(matshell)
        PETSc.finalize(petsclib)
    end
end
