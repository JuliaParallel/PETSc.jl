using Test
using PETSc
using LinearAlgebra: norm

@testset "VecSeq" begin
    N = 10
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        x = rand(PetscScalar, N)
        petsc_x = PETSc.VecSeq(petsclib, x)
        @test length(petsc_x) == N
        @test norm(petsc_x) ≈ norm(x)

        # make sure the viewer works
        _stdout = stdout
        (rd, wr) = redirect_stdout()
        @show petsc_x
        @test readline(rd) == "petsc_x = Vec Object: 1 MPI processes"
        @test readline(rd) == "  type: seq"
        redirect_stdout(_stdout)

        _stdout = stdout
        (rd, wr) = redirect_stdout()
        show(stdout, "text/plain", petsc_x)
        @test readline(rd) == "Vec Object: 1 MPI processes"
        @test readline(rd) == "  type: seq"
        redirect_stdout(_stdout)

        @test PETSc.ownershiprange(petsc_x) == 1:N
        @test PETSc.ownershiprange(petsc_x, false) == 0:(N - 1)

        PETSc.with_unsafe_localarray!(petsc_x) do x2
            @test x2 == x
        end

        PETSc.destroy(petsc_x)
        PETSc.finalize(petsclib)
    end
end
