using Test
using PETSc

const LEAKY = joinpath(@__DIR__, "fixtures", "leaky.jl")
const CONSTRUCTS = joinpath(@__DIR__, "fixtures", "constructs.jl")

@testset "audit_petsc_file" begin
    report = PETSc.audit_petsc_file(LEAKY; verbose = false)

    created = Set(v for (_, v, _) in report.created if v !== nothing)
    destroyed = Set(v for (_, v) in report.destroyed)

    @testset "creations" begin
        # high-level constructor, low-level creator, and tuple assignment
        for v in (:ksp, :v1, :v2, :x, :b, :dm, :mat)
            @test v in created
        end
    end

    @testset "releases" begin
        # both PETSc.destroy! and LibPETSc.VecDestroy count
        for v in (:ksp, :v1, :v2, :x)
            @test v in destroyed
        end
    end

    @testset "leaks" begin
        @test report.leaked == [:b, :dm, :mat]
    end

    @testset "two statements on one line" begin
        # v2 is created and destroyed on the same source line; the previous
        # text-matching implementation stopped after the first match
        line = only(l for (l, v, _) in report.created if v === :v2)
        @test any(report.destroyed) do (l, v)
            v === :v2 && l == line
        end
    end

    @testset "comments and strings are not code" begin
        # the fixture mentions destroy!(dm) in a comment and MatDestroy(mat)
        # inside a string; neither may register as a release
        @test !(:dm in destroyed)
        @test !(:mat in destroyed)
    end

    @testset "finalizer counts as a release" begin
        mktemp() do path, io
            write(
                io,
                """
                v = LibPETSc.VecCreateSeq(petsclib, comm, 10)
                finalizer(destroy, v)

                mat = PETSc.MatSeqAIJ(petsclib, 10, 10, 3)
                finalizer(m -> (destroy(m); data), mat)

                dm = PETSc.DMStag(petsclib, comm, bt, sz, 1, 1)
                """,
            )
            close(io)
            # v and mat are handed to the garbage collector; only dm leaks
            @test PETSc.audit_petsc_file(path; verbose = false).leaked == [:dm]
        end
    end

    @testset "borrowed pointers are not creations" begin
        mktemp() do path, io
            # VecPtr/MatPtr wrap a handle the caller may not own, so they must
            # not be reported as objects needing a destroy
            write(io, "f = PETSc.VecPtr(petsclib, f_ptr, false)\n")
            close(io)
            @test isempty(PETSc.audit_petsc_file(path; verbose = false).leaked)
        end
    end

    @testset "clean file has no leaks" begin
        mktemp() do path, io
            write(
                io,
                """
                v = LibPETSc.VecCreateSeq(petsclib, comm, 10)
                LibPETSc.VecDestroy(petsclib, v)
                """,
            )
            close(io)
            @test isempty(PETSc.audit_petsc_file(path; verbose = false).leaked)
        end
    end

    @testset "language constructs" begin
        # Every construct in the fixture is released except the two named here.
        # A failure names the construct that broke.
        report = PETSc.audit_petsc_file(CONSTRUCTS; verbose = false)
        @test report.leaked == [:never_freed, :quoted_v]
    end

    @testset "does not parse" begin
        mktemp() do path, io
            write(io, "v = LibPETSc.VecCreateSeq(petsclib, comm, 10\nfunction broken(\n")
            close(io)
            # a syntax error must not silently read as a clean bill of health
            @test_logs (:warn, r"does not parse") PETSc.audit_petsc_file(
                path;
                verbose = false,
            )
        end
    end

    @testset "helpers" begin
        @test PETSc.audit_creator(:VecCreateSeq) == "Vec"
        @test PETSc.audit_creator(:MatSeqAIJWithArrays) == "Mat"
        @test PETSc.audit_creator(:DMStag) == "DM"
        @test PETSc.audit_creator(:solve!) === nothing

        @test PETSc.audit_destroyer(:destroy)
        @test PETSc.audit_destroyer(:destroy!)
        @test PETSc.audit_destroyer(:VecDestroy)
        @test PETSc.audit_destroyer(:finalizer)
        @test PETSc.audit_creator(:VecPtr) === nothing
        @test PETSc.audit_creator(:MatPtr) === nothing
        @test !PETSc.audit_destroyer(:Destroy)
        @test !PETSc.audit_destroyer(:assemble!)

        @test PETSc.audit_callee(:(PETSc.LibPETSc.f(x))) === :f
        @test PETSc.audit_callee(:(PETSc.f.(xs))) === :f
        @test !PETSc.audit_isbroadcast(:(a.b))
        @test PETSc.audit_isbroadcast(:(f.(x)))
        @test PETSc.audit_argnames(:(destroy!(v))) == [:v]
        @test PETSc.audit_argnames(:(destroy!(v...))) == [:v]
        @test PETSc.audit_argnames(:(destroy!.([v, w]))) == [:v, :w]
        @test PETSc.audit_argnames(:(VecDestroy(lib, v))) == [:lib, :v]
        @test PETSc.audit_callee(:(f(x))) === :f
        @test PETSc.audit_callee(:(x + 1)) === :+
    end
end
