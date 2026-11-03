using Test
using MPI
MPI.Initialized() || MPI.Init()
using PETSc
using LinearAlgebra: mul!, I
using SparseArrays: sparse

# Index ranges and vector lengths that do not match throw before PETSc reads past
# a buffer, and range assignment converts to the scalar type like Julia's setindex!
@testset "Index and length checks" begin
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        comm = MPI.COMM_SELF
        n = 5

        @testset "Vec range assignment ($PetscScalar)" begin
            v = PETSc.PetscVec(petsclib, comm, n)
            fill!(v, 0)
            v[2:3] = 7
            v[[1, 5]] = [1, 2]
            PETSc.assemble!(v)
            @test v[:] == PetscScalar[1, 7, 7, 0, 2]
            @test_throws DimensionMismatch v[1:3] = [1, 2]
            @test_throws DimensionMismatch v[1:2] = [1, 2, 3]
            PETSc.destroy!(v)
        end

        @testset "Mat row and column assignment ($PetscScalar)" begin
            A = PETSc.PetscMat(petsclib, comm, sparse(PetscScalar(1) * I, n, n))
            @test_throws DimensionMismatch A[1, 1:3] = [1.0, 2.0]
            @test_throws DimensionMismatch A[1:3, 1] = [1.0, 2.0, 3.0, 4.0]
            PETSc.destroy!(A)
        end

        @testset "MatShell product lengths ($PetscScalar)" begin
            shell = PETSc.MatShell(petsclib, comm, n, n) do y, x
                y .= x
            end
            x = PETSc.PetscVec(petsclib, comm, n)
            short = PETSc.PetscVec(petsclib, comm, n - 1)
            fill!(x, 1)
            fill!(short, 1)
            @test_throws DimensionMismatch shell * short
            @test_throws DimensionMismatch mul!(short, shell, x)
            @test_throws DimensionMismatch shell' * short
            y = shell * x
            @test y[:] == fill(PetscScalar(1), n)
            foreach(PETSc.destroy!, (y, x, short, shell))
        end

        @testset "Null space from Vecs ($PetscScalar)" begin
            v1 = PETSc.PetscVec(petsclib, comm, n)
            v2 = PETSc.PetscVec(petsclib, comm, n)
            fill!(v1, 1)
            fill!(v2, 2)
            nsp = PETSc.mat_nullspace_create(petsclib, comm, [v1, v2])
            @test nsp isa PETSc.LibPETSc.AbstractMatNullSpace
            PETSc.destroy!(nsp)
            PETSc.destroy!(v2)
            @test_throws ArgumentError PETSc.mat_nullspace_create(petsclib, comm, [v1, v2])
            PETSc.destroy!(v1)
        end
    end
end
