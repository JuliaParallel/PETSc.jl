using Test
using PETSc
using LinearAlgebra: norm

@testset "VecSeqWithArray" begin
    N = 10
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        x = rand(PetscScalar, N)
        petsc_x = PETSc.VecSeqWithArray(petsclib, x)
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

        PETSc.withlocalarray!(petsc_x) do x2
            @test x2 == x
        end

        PETSc.destroy(petsc_x)
        PETSc.finalize(petsclib)
    end
end

@testset "VecSeq" begin
    N = 10
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        petsc_x = PETSc.VecSeq(petsclib, N)
        @test length(petsc_x) == N

        @test PETSc.ownershiprange(petsc_x) == 1:N
        @test PETSc.ownershiprange(petsc_x, false) == 0:(N - 1)

        x = rand(PetscScalar, N)
        PETSc.withlocalarray!(petsc_x) do x2
            x2 .= x
        end
        @test norm(petsc_x) ≈ norm(x)

        PETSc.destroy(petsc_x)
        PETSc.finalize(petsclib)
    end
end
