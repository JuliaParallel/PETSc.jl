# test/vecmat/mat_vec_methods.jl
# Base and Mat methods: copyto! between vectors, fill! on a matrix, row zeroing
# for Dirichlet conditions, matrix options, the diagonal and isassembled; the
# LinearAlgebra and array-interface contracts (transposed products, norms,
# isapprox, size and axes).

using Test
using PETSc
using MPI
using LinearAlgebra: norm, opnorm, mul!, issymmetric, ishermitian
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

    # a nonsymmetric matrix, so a product with A and one with its transpose differ
    nonsymmetric() = PETSc.PetscMat(petsclib, sparse(PetscScalar[1 2; 3 4]))

    @testset "transposed products" begin
        A = nonsymmetric()
        x = PETSc.PetscVec(petsclib, PetscScalar[1, 1])
        y = PETSc.PetscVec(petsclib, PetscScalar[0, 0])
        @test mul!(y, A', x) === y
        @test y[:] == PetscScalar[4, 6]
        @test mul!(y, transpose(A), x) === y
        @test y[:] == PetscScalar[4, 6]
        @test (A' * x)[:] == PetscScalar[4, 6]
        @test (transpose(A) * x)[:] == PetscScalar[4, 6]
        @test (A * x)[:] == PetscScalar[3, 7]
        foreach(PETSc.destroy!, (y, x, A))
    end

    @testset "issymmetric, ishermitian" begin
        L, A = laplacian(), nonsymmetric()
        @test issymmetric(L) && ishermitian(L)
        @test !issymmetric(A) && !ishermitian(A)
        foreach(PETSc.destroy!, (A, L))
    end

    @testset "norms" begin
        v = PETSc.PetscVec(petsclib, PetscScalar[3, -4])
        @test norm(v) == 5
        @test norm(v, 1) == 7
        @test norm(v, 2) == 5
        @test norm(v, Inf) == 4
        @test norm(v, LibPETSc.NORM_1) == 7
        @test_throws "only for p = 1, 2 and Inf, got p = 3" norm(v, 3)

        A = nonsymmetric()
        @test norm(A) ≈ sqrt(30)
        @test norm(A, 2) ≈ sqrt(30)
        @test opnorm(A, 1) == 6       # largest column sum
        @test opnorm(A, Inf) == 7     # largest row sum
        @test_throws "use opnorm(A, 1) or opnorm(A, Inf)" norm(A, 1)
        @test_throws "only for p = 1 and p = Inf, got p = 2" opnorm(A)
        foreach(PETSc.destroy!, (A, v))
    end

    @testset "isapprox keywords" begin
        v = PETSc.PetscVec(petsclib, PetscScalar[1, 2])
        w = PETSc.PetscVec(petsclib, PetscScalar[1, 2 + 1e-6])
        @test isapprox(v, v)
        @test !isapprox(v, w)
        @test isapprox(v, w; rtol = 1e-5)
        @test isapprox(v, w; atol = 1e-5)
        foreach(PETSc.destroy!, (w, v))
    end

    @testset "size, axes, ndims" begin
        v = PETSc.PetscVec(petsclib, PetscScalar[1, 2, 3])
        @test size(v, 1) === 3
        @test size(v, 2) === 1
        @test_throws ArgumentError size(v, 0)
        @test ndims(v) == 1

        A = PETSc.PetscMat(petsclib, sparse(PetscScalar[1 2 3; 4 5 6]))
        @test size(A, 1) === 2
        @test size(A, 2) === 3
        @test size(A, 3) === 1
        @test_throws ArgumentError size(A, 0)
        @test axes(A) == (Base.OneTo(2), Base.OneTo(3))
        @test axes(A, 2) == Base.OneTo(3)
        @test ndims(A) == 2
        @test ndims(typeof(A)) == 2
        foreach(PETSc.destroy!, (A, v))
    end

    @testset "copyto! from a SparseMatrixCSC" begin
        S = sparse(PetscScalar[1 0 2; 0 3 0])
        M = PETSc.PetscMat(petsclib, 2, 3, 3)
        @test copyto!(M, S) === M
        PETSc.assemble!(M)
        @test M[:, :] == Matrix(S)
        @test_throws "does not fit" copyto!(M, sparse(PetscScalar[1 2 3 4; 5 6 7 8]))
        PETSc.destroy!(M)
    end

    PETSc.finalize(petsclib)
end
