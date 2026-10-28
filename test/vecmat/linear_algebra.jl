using Test
using MPI
MPI.Initialized() || MPI.Init()
using PETSc
using LinearAlgebra: dot, axpy!, axpby!, rmul!, lmul!

# Vec and Mat arithmetic on COMM_WORLD, so on several ranks it checks that the
# operations are collective and give the global result. Each rank compares the
# entries it owns against the same computation on Julia arrays.

# The entries of `v` this rank owns, and the matching slice of a Julia vector
owned(v::PETSc.AbstractPetscVec) = v[PETSc.ownership_range(v)]
owned(v::PETSc.AbstractPetscVec, ref::AbstractVector) = ref[PETSc.ownership_range(v)]

function petsc_vec(petsclib, comm, vals)
    v = PETSc.PetscVec(petsclib, comm, length(vals))
    r = PETSc.ownership_range(v)
    v[r] = vals[r]
    return PETSc.assemble!(v)
end

# Rows this rank owns of a distributed AIJ matrix filled from `entries`, a list of
# (i, j, value) with 1-based indices
function petsc_mat(petsclib, comm, n, entries)
    PetscInt = petsclib.PetscInt
    A = LibPETSc.MatCreateAIJ(petsclib, comm, PetscInt(LibPETSc.PETSC_DECIDE),
                              PetscInt(LibPETSc.PETSC_DECIDE), PetscInt(n), PetscInt(n),
                              PetscInt(3), C_NULL, PetscInt(3), C_NULL)
    rows = PETSc.ownership_range(A)
    for (i, j, a) in entries
        i in rows && (A[i, j] = a)
    end
    return PETSc.assemble!(A)
end

@testset "Vec and Mat linear algebra" begin
    comm = MPI.COMM_WORLD
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        n = 4 * MPI.Comm_size(comm) + 3
        # Complex entries in a complex build, so `dot` has something to conjugate
        z = PetscScalar <: Complex ? im : 0
        xs = PetscScalar[i + z * (n - i) for i in 1:n]
        ys = PetscScalar[2i - 1 - z * i for i in 1:n]
        α = PetscScalar(3 + 2z)
        β = PetscScalar(-0.5)

        @testset "Vec ($PetscScalar)" begin
            x = petsc_vec(petsclib, comm, xs)
            y = petsc_vec(petsclib, comm, ys)

            w = copy(x)
            @test w isa PETSc.PetscVec
            @test w.ptr != x.ptr
            @test owned(w) == owned(x, xs)
            fill!(w, 0)
            @test owned(x) == owned(x, xs)

            @test dot(x, y) ≈ dot(xs, ys)
            @test dot(x, y) ≈ conj(dot(y, x))
            @test dot(x, x) ≈ sum(abs2, xs)

            results = [
                x + y => xs + ys,
                x - y => xs - ys,
                -x => -xs,
                α * x => α * xs,
                x * α => xs * α,
                x / α => xs / α,
                2 * x => 2 * xs,
            ]
            for (v, ref) in results
                @test v isa PETSc.PetscVec
                @test owned(v) ≈ owned(v, ref)
            end
            # Division is entry by entry, as for a Julia vector, not a scale by 1/α
            @test owned(results[6].first) == owned(x, xs ./ α)

            u = copy(y)
            @test axpy!(α, x, u) === u
            @test owned(u) ≈ owned(u, α * xs + ys)
            @test axpby!(α, x, β, u) === u
            @test owned(u) ≈ owned(u, α * xs + β * (α * xs + ys))
            @test rmul!(u, β) === u
            @test lmul!(α, u) === u
            @test owned(u) ≈ owned(u, α * β * (α * xs + β * (α * xs + ys)))

            # The same vector on both sides, as for Julia arrays
            u = copy(x)
            axpy!(2, u, u)
            @test owned(u) ≈ owned(u, 3 * xs)
            axpby!(1, u, 1, u)
            @test owned(u) ≈ owned(u, 6 * xs)

            short = PETSc.PetscVec(petsclib, comm, n - 1)
            @test_throws DimensionMismatch x + short
            @test_throws DimensionMismatch dot(x, short)
            @test_throws DimensionMismatch axpy!(1, x, short)
            if PetscScalar <: Real
                @test_throws InexactError 2im * x
            end

            foreach(PETSc.destroy!, (x, y, w, u, short, first.(results)...))
        end

        @testset "Mat ($PetscScalar)" begin
            # Upper and lower bidiagonal, so A ± B has the union of the patterns
            a = [(i, j, PetscScalar(i + 2j + z)) for i in 1:n for j in (i, i + 1) if j <= n]
            b = [(i, j, PetscScalar(3i - j)) for i in 1:n for j in (i - 1, i) if j >= 1]
            dense(entries) = (M = zeros(PetscScalar, n, n); foreach(((i, j, v),) -> M[i, j] = v, entries); M)
            As, Bs = dense(a), dense(b)
            A = petsc_mat(petsclib, comm, n, a)
            B = petsc_mat(petsclib, comm, n, b)
            rows = PETSc.ownership_range(A)

            C = copy(A)
            @test C isa PETSc.PetscMat
            @test C.ptr != A.ptr
            @test C == A
            @test C[rows, 1:n] == As[rows, :]
            fill!(C, 0)
            @test A[rows, 1:n] == As[rows, :]

            results = [
                A + B => As + Bs,
                A - B => As - Bs,
                -A => -As,
                α * A => α * As,
                A * α => As * α,
            ]
            for (M, ref) in results
                @test M isa PETSc.PetscMat
                @test M[rows, 1:n] ≈ ref[rows, :]
            end
            @test A[rows, 1:n] == As[rows, :]
            @test B[rows, 1:n] == Bs[rows, :]

            short = petsc_mat(petsclib, comm, n - 1, [(1, 1, PetscScalar(1))])
            @test_throws DimensionMismatch A + short

            foreach(PETSc.destroy!, (A, B, C, short, first.(results)...))
        end
    end
end
