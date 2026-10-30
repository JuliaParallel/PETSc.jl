using Test
using MPI
MPI.Initialized() || MPI.Init()
using PETSc

isdefined(Main, :PETScTestUtils) || include(joinpath(@__DIR__, "..", "testutils.jl"))
using .PETScTestUtils: without_petsc_traceback
using LinearAlgebra: dot, I
using SparseArrays: sparse

# A PetscVec as an AbstractVector, on COMM_WORLD: generic Julia methods either reach a
# PETSc override and give the global answer on every rank, or throw at the first
# scalar access instead of making one C call per entry.

function petsc_vec(petsclib, comm, vals)
    v = PETSc.PetscVec(petsclib, comm, length(vals))
    r = PETSc.ownership_range(v)
    v[r] = vals[r]
    return PETSc.assemble!(v)
end

@testset "PetscVec as AbstractVector" begin
    comm = MPI.COMM_WORLD
    nranks = MPI.Comm_size(comm)
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        n = 4 * nranks + 3
        xs = PetscScalar[(-1)^i * i for i in 1:n]
        x = petsc_vec(petsclib, comm, xs)
        y = copy(x)

        @testset "type ($PetscScalar)" begin
            @test x isa AbstractVector{PetscScalar}
            @test eltype(x) === PetscScalar
            @test size(x) === (n,)
            @test axes(x) === (Base.OneTo(n),)
            @test lastindex(x) == n
            @test size(x, 2) == 1
            @test similar(x, PetscScalar) isa PETSc.PetscVec
            @test similar(x, Int) isa Vector{Int}
            # The REPL shows the handle, not the entries
            @test startswith(sprint(show, MIME("text/plain"), x), "PETSc")
        end

        @testset "scalar access ($PetscScalar)" begin
            r = PETSc.ownership_range(x)
            i = first(r)
            @test x[i] == xs[i]
            @test x[i] isa PetscScalar
            @test first(r) == i && x[last(r)] == xs[last(r)]
            x[i] = 7
            PETSc.assemble!(x)
            @test x[i] == 7
            x[i] = xs[i]
            PETSc.assemble!(x)
            @test x[[i, last(r)]] == xs[[i, last(r)]]
            if nranks > 1
                foreign = first(PETSc.ownership_range(x)) == 1 ? n : 1
                @test_throws LibPETSc.PetscError without_petsc_traceback(() -> x[foreign], petsclib)
            end
        end

        @testset "== and hash ($PetscScalar)" begin
            @test x == y
            @test isequal(x, y)
            short = petsc_vec(petsclib, comm, xs[1:(n - 1)])
            @test x != short
            z = copy(x)
            r = PETSc.ownership_range(z)
            MPI.Comm_rank(comm) == nranks - 1 && (z[last(r):last(r)] = [-xs[end]])
            z = PETSc.assemble!(z)
            @test x != z
            @test y == x    # the same answer on every rank
            nanv = petsc_vec(petsclib, comm, fill(PetscScalar(NaN), n))
            @test nanv != nanv
            @test isequal(nanv, nanv)
            zs = petsc_vec(petsclib, comm, zeros(PetscScalar, n))
            ms = petsc_vec(petsclib, comm, fill(PetscScalar(-0.0), n))
            @test zs == ms
            @test !isequal(zs, ms)
            @test_throws "has no hash" hash(x)
            @test_throws ArgumentError Dict(x => 1)
            @test_throws ArgumentError Set([x])
            @test IdDict(x => 1)[x] == 1
            foreach(PETSc.destroy!, (short, z, nanv, zs, ms))
        end

        @testset "reductions ($PetscScalar)" begin
            @test sum(x) ≈ sum(xs)
            @test any(v -> real(v) > n - 1, x) == any(v -> real(v) > n - 1, xs)
            @test any(v -> real(v) > n, x) == false
            @test all(v -> abs(v) >= 1, x)
            @test !all(v -> real(v) > 0, x)
            @test count(v -> real(v) < 0, x) == count(v -> real(v) < 0, xs)
            @test count(v -> real(v) < 0, x; init = 10) == 10 + count(v -> real(v) < 0, xs)
            @test sum(abs2, x) ≈ sum(abs2, xs)
            @test sum(v -> real(v) > 0, x) == sum(v -> real(v) > 0, xs)   # Bool results sum to an Int
            @test sum(abs, PETSc.PetscVec(petsclib, comm, 0)) == 0               # Base's neutral element
            @test mapreduce(abs, +, x) ≈ mapreduce(abs, +, xs)
            @test mapreduce(abs, +, x; init = 5) ≈ 5 + mapreduce(abs, +, xs)
            @test reduce(+, x) ≈ sum(xs)
            @test extrema(real, x) == extrema(real, xs)
            @test prod(v -> v / 2, x) ≈ prod(v -> v / 2, xs)
            # the partial results combine in rank order, also for an op that is not commutative
            @test mapreduce(v -> (real(v), real(v)), (a, b) -> (a[1], b[2]), x) == (real(xs[1]), real(xs[end]))
            if PetscScalar <: Real
                @test maximum(x) == maximum(xs)
                @test minimum(x) == minimum(xs)
                e = PETSc.PetscVec(petsclib, comm, 0)
                @test_throws ArgumentError maximum(e)
                @test mapreduce(abs, +, e) == mapreduce(abs, +, PetscScalar[])   # as Base: zero
                @test mapreduce(abs, +, e; init = 3) == 3
                @test sum(e) == 0
                PETSc.destroy!(e)
            end
        end

        @testset "conversions ($PetscScalar)" begin
            @test collect(x) == xs
            @test Vector(x) == xs
            @test Array(x) == xs
            dest = zeros(PetscScalar, n)
            @test copyto!(dest, x) === dest
            @test dest == xs
            @test_throws DimensionMismatch copyto!(zeros(PetscScalar, n + 1), x)
        end

        @testset "generic algorithms ($PetscScalar)" begin
            @test dot(x, y) ≈ dot(xs, xs)
            @test isapprox(x, y)
            w = similar(x)
            @test copyto!(w, x) === w
            @test w == x
            if nranks == 1
                x .= 1:n
                @test collect(x) == PetscScalar.(1:n)
                @test sum(x) ≈ n * (n + 1) / 2
                x .= xs
            else
                # a Julia array in a broadcast must have the local length
                @test_throws DimensionMismatch x .= 1:n
            end
            PETSc.destroy!(w)
        end

        if nranks > 1
            @testset "different layouts ($PetscScalar)" begin
                # Same global length, one more entry on rank 0: every rank throws
                nloc = length(PETSc.ownership_range(x))
                rank = MPI.Comm_rank(comm)
                shift = rank == 0 ? 1 : rank == 1 ? -1 : 0
                v = PETSc.PetscVec(petsclib, comm, LibPETSc.PETSC_DETERMINE; local_size = nloc + shift)
                fill!(v, 0)
                @test_throws "different parallel layouts" v == x
                PETSc.destroy!(v)
            end
        end

        PETSc.destroy!(x)
        PETSc.destroy!(y)
    end
end

# On one rank a vector that `similar`, `copy`, a broadcast or `ksp \ b` returns has a finalizer,
# so the garbage collector frees it. The check holds a second PETSc reference on the object,
# drops the Julia handle and collects: the count is back to 1 only if the finalizer ran.
@noinline function references_after_collection(make, petsclib, settle = () -> nothing)
    LibPETSc = PETSc.LibPETSc
    ptr = let w = make()
        LibPETSc.PetscObjectReference(petsclib, w)
        w.ptr
    end
    GC.gc()
    GC.gc()
    settle()
    held = PETSc.PetscVec(petsclib, ptr; own = true)
    refs = LibPETSc.PetscObjectGetReference(petsclib, held)
    PETSc.destroy!(held)
    return refs
end

@testset "vectors from similar, copy, a broadcast and ksp \\ b are freed by the GC" begin
    if MPI.Comm_size(MPI.COMM_WORLD) == 1
        petsclib = PETSc.petsclibs[1]
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        n = 10
        x = PETSc.PetscVec(petsclib, PetscScalar[1:n;])
        A = PETSc.PetscMat(petsclib, sparse(PetscScalar(2) * I, n, n))
        ksp = PETSc.KSP(A; ksp_type = "preonly", pc_type = "jacobi")

        @test references_after_collection(() -> similar(x), petsclib) == 1
        @test references_after_collection(() -> similar(x, PetscScalar), petsclib) == 1
        @test references_after_collection(() -> copy(x), petsclib) == 1
        @test references_after_collection(() -> 2 .* x, petsclib) == 1
        @test references_after_collection(() -> x / 2, petsclib) == 1
        # a KSP keeps a reference to its last solution until the next solve
        scratch = similar(x)
        @test references_after_collection(() -> ksp \ x, petsclib, () -> PETSc.solve!(scratch, ksp, x)) == 1
        PETSc.destroy!(scratch)
        foreach(PETSc.destroy!, (ksp, A, x))
    end
end
