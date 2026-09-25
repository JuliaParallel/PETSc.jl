# test/handle_arrays.jl
# C arrays of handles: the vector a wrapper returns goes back to the matching
# release function, and a vector built in Julia is refused.

using Test
using PETSc
using MPI
using SparseArrays

MPI.Initialized() || MPI.Init()

@testset "arrays of handles" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    comm = MPI.COMM_SELF
    PetscInt = petsclib.PetscInt
    PetscScalar = petsclib.PetscScalar
    LibPETSc = PETSc.LibPETSc
    stride(n, first) = LibPETSc.ISCreateStride(petsclib, comm, PetscInt(n), PetscInt(first), PetscInt(1))

    @testset "MatCreateSubMatrices and MatDestroySubMatrices" begin
        A = PETSc.PetscMat(petsclib, sparse(PetscScalar[1 2 0; 3 4 5; 0 6 7]))
        rows, cols = [stride(2, 0)], [stride(2, 1)]

        sub = LibPETSc.MatCreateSubMatrices(petsclib, A, PetscInt(1), rows, cols, LibPETSc.MAT_INITIAL_MATRIX)
        @test length(sub) == 1
        @test all(PETSc.owns, sub)
        @test sub[1][:, :] == PetscScalar[2 0; 4 5]

        # MAT_REUSE_MATRIX refills the same vector
        A[2, 3] = 50
        PETSc.assemble!(A)
        @test LibPETSc.MatCreateSubMatrices(petsclib, A, PetscInt(1), rows, cols, LibPETSc.MAT_REUSE_MATRIX, sub) === sub
        @test sub[1][:, :] == PetscScalar[2 0; 4 50]
        @test_throws ArgumentError LibPETSc.MatCreateSubMatrices(petsclib, A, PetscInt(1), rows, cols, LibPETSc.MAT_REUSE_MATRIX)

        # only the returned vector is accepted, and only once
        @test_throws ArgumentError LibPETSc.MatDestroySubMatrices(petsclib, PetscInt(1), copy(sub))
        @test_throws DimensionMismatch LibPETSc.MatDestroySubMatrices(petsclib, PetscInt(2), sub)
        LibPETSc.MatDestroySubMatrices(petsclib, PetscInt(1), sub)
        @test sub[1].ptr == C_NULL
        @test_throws ArgumentError LibPETSc.MatDestroySubMatrices(petsclib, PetscInt(1), sub)

        foreach(PETSc.destroy!, (rows..., cols..., A))
    end

    @testset "VecNestGetSubVecsRead and VecNestRestoreSubVecsRead" begin
        blocks = [PETSc.PetscVec(petsclib, PetscScalar[1, 2]), PETSc.PetscVec(petsclib, PetscScalar[3])]
        index_sets = [stride(2, 0), stride(1, 2)]
        X = LibPETSc.VecCreateNest(petsclib, comm, PetscInt(2), index_sets, blocks)

        N, sx = LibPETSc.VecNestGetSubVecsRead(petsclib, X)
        @test N == 2
        @test !any(PETSc.owns, sx)
        @test [v.ptr for v in sx] == [v.ptr for v in blocks]
        @test_throws ArgumentError LibPETSc.VecNestRestoreSubVecsRead(petsclib, X, N, copy(sx))
        LibPETSc.VecNestRestoreSubVecsRead(petsclib, X, N, sx)
        @test all(v -> v.ptr == C_NULL, sx)

        # the read lock is released: X can be written again
        LibPETSc.VecSet(petsclib, X, PetscScalar(0))
        @test blocks[1][:] == PetscScalar[0, 0]

        foreach(PETSc.destroy!, (X, index_sets..., blocks...))
    end

    @testset "PCASMCreateSubdomains2D and PCASMDestroySubdomains" begin
        nsub, is, is_local = LibPETSc.PCASMCreateSubdomains2D(
            petsclib, PetscInt(4), PetscInt(4), PetscInt(2), PetscInt(2), PetscInt(1), PetscInt(1),
        )
        @test nsub == length(is) == length(is_local) == 4
        @test all(PETSc.owns, is)
        @test all(x -> LibPETSc.PetscObjectGetReference(petsclib, x) == 1, is)
        @test_throws DimensionMismatch LibPETSc.PCASMDestroySubdomains(petsclib, PetscInt(3), is, is_local)
        @test_throws ArgumentError LibPETSc.PCASMDestroySubdomains(petsclib, nsub, is, C_NULL)
        LibPETSc.PCASMDestroySubdomains(petsclib, nsub, is, is_local)
        @test all(x -> x.ptr == C_NULL, is)
        @test all(x -> x.ptr == C_NULL, is_local)
    end

    PETSc.finalize(petsclib)
end
