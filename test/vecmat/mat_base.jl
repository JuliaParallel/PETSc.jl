using Test
using MPI
MPI.Initialized() || MPI.Init()
using PETSc
using LinearAlgebra: mul!
using SparseArrays: sparse

# A PetscMat as an AbstractMatrix, on COMM_WORLD: equality, conversions and reductions
# give the global answer on every rank.

# A distributed AIJ matrix from a list of (i, j, value), 1-based; each rank writes the
# rows it owns. `local_rows` sets this rank's row count instead of PETSc's choice.
function petsc_mat(petsclib, comm, n, entries; local_rows = LibPETSc.PETSC_DECIDE)
    PetscInt = petsclib.PetscInt
    A = LibPETSc.MatCreateAIJ(petsclib, comm, PetscInt(local_rows),
                              PetscInt(LibPETSc.PETSC_DECIDE), PetscInt(n), PetscInt(n),
                              PetscInt(4), C_NULL, PetscInt(4), C_NULL)
    rows = PETSc.ownership_range(A)
    for (i, j, a) in entries
        i in rows && (A[i, j] = a)
    end
    return PETSc.assemble!(A)
end

dense(n, entries, T) = (D = zeros(T, n, n); foreach(((i, j, a),) -> D[i, j] = a, entries); D)

@testset "PetscMat as AbstractMatrix" begin
    comm = MPI.COMM_WORLD
    nranks = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        n = 4 * nranks + 3
        z = PetscScalar <: Complex ? im : 0
        # tridiagonal, with a zero stored explicitly on the last diagonal entry
        entries = [(i, j, PetscScalar(i == j == n ? 0 : 10i + j + z * (i - j)))
                   for i in 1:n for j in max(1, i - 1):min(n, i + 1)]
        ref = dense(n, entries, PetscScalar)
        A = petsc_mat(petsclib, comm, n, entries)
        # the same matrix without the stored zero
        B = petsc_mat(petsclib, comm, n, filter(e -> !(e[1] == e[2] == n), entries))

        @testset "type ($PetscScalar)" begin
            @test A isa AbstractMatrix{PetscScalar}
            @test eltype(A) === PetscScalar
            @test size(A) === (n, n)
            @test axes(A) === (Base.OneTo(n), Base.OneTo(n))
            @test size(A, 3) == 1
            @test startswith(sprint(show, MIME("text/plain"), A), "PETSc")
            r = PETSc.ownership_range(A)
            @test A[first(r), first(r)] == ref[first(r), first(r)]
            @test A[r, 1:n] == ref[r, :]
            if nranks > 1
                foreign = first(r) == 1 ? n : 1
                @test_throws LibPETSc.PetscError A[foreign, foreign]
            end
        end

        @testset "== and equal ($PetscScalar)" begin
            @test A == A
            @test A == B                 # a stored zero makes no difference
            @test isequal(A, B)
            @test !PETSc.equal(A, B)     # MatEqual compares what is stored
            @test PETSc.equal(A, A)
            C = copy(A)
            r = PETSc.ownership_range(C)
            rank == nranks - 1 && (C[last(r), last(r)] = 1)
            PETSc.assemble!(C)
            @test A != C
            @test C != A                 # the same answer on every rank
            short = petsc_mat(petsclib, comm, n - 1, [(1, 1, PetscScalar(1))])
            @test A != short
            N = petsc_mat(petsclib, comm, n, [(i, i, PetscScalar(NaN)) for i in 1:n])
            @test N != N
            @test isequal(N, N)
            @test PETSc.equal(N, N)
            @test_throws "has no hash" hash(A)
            @test_throws ArgumentError Dict(A => 1)
            @test IdDict(A => 1)[A] == 1
            foreach(PETSc.destroy!, (C, short, N))
        end

        @testset "isapprox ($PetscScalar)" begin
            @test A ≈ B
            C = copy(A)
            rank == 0 && (C[1, 1] = C[1, 1] + 1)
            PETSc.assemble!(C)
            @test !(A ≈ C)
            @test isapprox(A, C; atol = 2)
            PETSc.destroy!(C)
        end

        @testset "conversions and sum ($PetscScalar)" begin
            @test collect(A) == ref
            @test Matrix(A) == ref
            D = zeros(PetscScalar, n, n)
            @test copyto!(D, A) === D
            @test D == ref
            @test_throws DimensionMismatch copyto!(zeros(PetscScalar, n, n + 1), A)
            @test sum(A) ≈ sum(ref)
            W = copy(B)
            @test copyto!(W, A) === W
            @test W == A
            S = petsc_mat(petsclib, comm, n - 1, [])
            @test_throws DimensionMismatch copyto!(S, A)
            foreach(PETSc.destroy!, (W, S))
        end

        @testset "sparse reductions ($PetscScalar)" begin
            # row 2 empty, and the explicitly stored zero at (n, n) kept, as `sparse` keeps it
            es = filter(e -> e[1] != 2, entries)
            S = sparse(first.(es), getindex.(es, 2), last.(es), n, n)
            P = petsc_mat(petsclib, comm, n, es)
            @test sum(abs2, P) ≈ sum(abs2, S)
            @test sum(x -> real(x) > 0, P) == sum(x -> real(x) > 0, S)
            @test prod(x -> 1 + real(x) / 1000, P) ≈ prod(x -> 1 + real(x) / 1000, S)
            @test maximum(real, P) == maximum(real, S)
            @test minimum(abs, P) == minimum(abs, S)
            @test extrema(real, P) == extrema(real, S)
            @test reduce(max, real.(collect(P))) == maximum(real, S)
            @test mapreduce(abs, +, P; init = 5) ≈ 5 + mapreduce(abs, +, S)
            @test count(x -> real(x) > 0, P) == count(x -> real(x) > 0, S)
            @test count(iszero, P) == count(iszero, S)
            @test any(x -> real(x) > 10n, P) == any(x -> real(x) > 10n, S)
            @test any(iszero, P)
            @test !all(x -> real(x) > 0, P)
            @test all(x -> abs(x) < 100n, P)
            if PetscScalar <: Real
                @test maximum(P) == maximum(S)
                @test minimum(P) == minimum(S)
            end
            # nothing stored: every entry is an implicit zero
            Z = petsc_mat(petsclib, comm, n, [])
            @test maximum(abs, Z) == 0
            @test sum(x -> x + 1, Z) == n^2
            @test count(iszero, Z) == n^2
            foreach(PETSc.destroy!, (P, Z))
        end

        @testset "products ($PetscScalar)" begin
            x = PETSc.PetscVec(petsclib, comm, n)
            fill!(x, 1)
            y = A * x
            @test y isa PETSc.PetscVec
            @test collect(y) ≈ ref * ones(PetscScalar, n)
            mul!(y, A, x)
            @test collect(y) ≈ ref * ones(PetscScalar, n)
            foreach(PETSc.destroy!, (x, y))
        end

        if nranks > 1
            @testset "different layouts ($PetscScalar)" begin
                # Same size, one more row on rank 0: every rank throws
                nloc = length(PETSc.ownership_range(A))
                shift = rank == 0 ? 1 : rank == 1 ? -1 : 0
                L = petsc_mat(petsclib, comm, n, entries; local_rows = nloc + shift)
                @test_throws "different row layouts" A == L
                PETSc.destroy!(L)
            end
        end

        foreach(PETSc.destroy!, (A, B))
    end
end
