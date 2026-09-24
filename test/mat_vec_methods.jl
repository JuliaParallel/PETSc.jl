# test/mat_vec_methods.jl
# Base and Mat methods: copyto! between vectors, fill! on a matrix, row zeroing
# for Dirichlet conditions, matrix options, the diagonal and isassembled.

using Test
using PETSc
using MPI
using SparseArrays: sparse

MPI.Initialized() || MPI.Init()

@testset "Mat and Vec methods" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    comm = MPI.COMM_SELF
    PetscScalar = petsclib.PetscScalar
    LibPETSc = PETSc.LibPETSc

    # the 1D Laplacian on three points
    laplacian() = PETSc.PetscMat(petsclib, sparse(PetscScalar[2 -1 0; -1 2 -1; 0 -1 2]))

    @testset "copyto! between vectors" begin
        v = PETSc.PetscVec(petsclib, PetscScalar[1, 2, 3])
        w = PETSc.PetscVec(petsclib, PetscScalar[0, 0, 0])
        @test copyto!(w, v) === w
        @test w[:] == PetscScalar[1, 2, 3]
        u = PETSc.PetscVec(petsclib, PetscScalar[0, 0])
        @test_throws DimensionMismatch copyto!(u, v)
        foreach(PETSc.destroy!, (u, w, v))
    end

    @testset "isassembled, fill!, set_option!" begin
        A = PETSc.PetscMat(petsclib, 3, 3, 1)
        @test !PETSc.isassembled(A)
        for i in 1:3
            A[i, i] = 2
        end
        PETSc.assemble!(A)
        @test PETSc.isassembled(A)

        @test fill!(A, 0) === A
        @test A[:, :] == zeros(PetscScalar, 3, 3)
        @test_throws ArgumentError fill!(A, 1)

        # one nonzero per row was preallocated; a second needs the option off
        @test PETSc.set_option!(A, LibPETSc.MAT_NEW_NONZERO_ALLOCATION_ERR, false) === A
        A[1, 2] = 5
        PETSc.assemble!(A)
        @test A[1, 2] == 5
        PETSc.destroy!(A)
    end

    @testset "diagonal!" begin
        A = laplacian()
        d = PETSc.PetscVec(petsclib, PetscScalar[0, 0, 0])
        @test PETSc.diagonal!(d, A) === d
        @test d[:] == PetscScalar[2, 2, 2]
        PETSc.destroy!(d)
        PETSc.destroy!(A)
    end

    @testset "zero_rows!" begin
        A = laplacian()
        @test PETSc.zero_rows!(A, [0], 3) === A
        @test A[1, :] == PetscScalar[3, 0, 0]
        @test A[2, :] == PetscScalar[-1, 2, -1]

        # with x and b, b takes diag * x on the zeroed rows
        x = PETSc.PetscVec(petsclib, PetscScalar[4, 0, 7])
        b = PETSc.PetscVec(petsclib, PetscScalar[1, 1, 1])
        PETSc.zero_rows!(A, [0, 2]; x, b)
        @test A[3, :] == PetscScalar[0, 0, 1]
        @test b[:] == PetscScalar[4, 1, 7]
        @test_throws ArgumentError PETSc.zero_rows!(A, [0]; x)
        foreach(PETSc.destroy!, (b, x, A))
    end

    @testset "zero_rows_local!" begin
        # a DMDA matrix carries the local-to-global mapping zero_rows_local! needs
        da = PETSc.DMDA(petsclib, comm, (PETSc.DM_BOUNDARY_NONE,), (4,), 1, 1)
        A = LibPETSc.DMCreateMatrix(petsclib, da)
        for i in 1:4
            A[i, i] = 2
        end
        PETSc.assemble!(A)
        @test PETSc.zero_rows_local!(A, [1], 5) === A
        @test A[2, 2] == 5
        @test A[3, 3] == 2
        PETSc.destroy!(A)
        PETSc.destroy!(da)
    end

    PETSc.finalize(petsclib)
end
