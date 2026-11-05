using Test
using MPI
MPI.Initialized() || MPI.Init()
using PETSc
using LinearAlgebra: mul!, I
using SparseArrays: sparse
using OffsetArrays: OffsetVector

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

        # PETSc_jll has no bounds checks of its own: an index outside the array read
        # or wrote arbitrary memory, and a Bool mask was read as the integers 1 and 0
        @testset "Vec bounds, masks and offset indices ($PetscScalar)" begin
            v = PETSc.PetscVec(petsclib, PetscScalar[1, 2, 3])
            @test_throws BoundsError v[0]
            @test_throws BoundsError v[4]
            @test_throws BoundsError v[[1, 4]]
            @test_throws BoundsError v[0] = 1
            @test_throws BoundsError v[2:4] = 1
            @test_throws BoundsError v[[true, false]]
            @test v[[true, false, true]] == PetscScalar[1, 3]
            v[[false, true, false]] = 9
            PETSc.assemble!(v)
            @test v[:] == PetscScalar[1, 9, 3]
            # the result takes the axes of the index vector, as for any AbstractVector
            r = OffsetVector([3, 1], 0:1)
            @test v[r] == OffsetVector(PetscScalar[3, 1], 0:1)
            v[r] = [7, 8]
            PETSc.assemble!(v)
            @test v[:] == PetscScalar[8, 9, 7]
            PETSc.destroy!(v)
        end

        @testset "Mat bounds, masks and offset indices ($PetscScalar)" begin
            A = PETSc.PetscMat(petsclib, comm, sparse(PetscScalar[1 2; 3 4]))
            @test_throws BoundsError A[0, 1]
            @test_throws BoundsError A[1, 3]
            @test_throws BoundsError A[0, 1] = 5
            @test_throws BoundsError A[1:3, 1]
            @test_throws BoundsError A[1, [true, false, true]]
            @test A[[false, true], [true, true]] == PetscScalar[3 4]
            @test A[1, [false, true]] == PetscScalar[2]
            @test A[OffsetVector([2, 1], 0:1), 1] == OffsetVector(PetscScalar[3, 1], 0:1)
            PETSc.destroy!(A)
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
