using Test

include(joinpath(@__DIR__, "..", "..", "scripts", "audit.jl"))

const LEAKY = joinpath(dirname(@__DIR__), "fixtures", "leaky.jl")
const CONSTRUCTS = joinpath(dirname(@__DIR__), "fixtures", "constructs.jl")

@testset "audit_file" begin
    report = Audit.audit_file(LEAKY; verbose = false)

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
                finalizer(destroy!, v)

                mat = PETSc.PetscMat(petsclib, 10, 10, 3)
                finalizer(m -> (destroy!(m); data), mat)

                dm = PETSc.DMStag(petsclib, comm, bt, sz, 1; stencil_width = 1)
                """,
            )
            close(io)
            # v and mat are handed to the garbage collector; only dm leaks
            @test Audit.audit_file(path; verbose = false).leaked == [:dm]
        end
    end

    @testset "borrowed pointers are not creations" begin
        mktemp() do path, io
            # PetscVec/PetscMat(petsclib, ptr; own) wrap a handle the caller may not own, so they must
            # not be reported as objects needing a destroy
            write(io, "g = PETSc.PetscVec(petsclib, g_ptr; own = false)\n")
            write(io, "h = PETSc.PetscMat(petsclib, h_ptr; own = false)\n")
            close(io)
            @test isempty(Audit.audit_file(path; verbose = false).leaked)
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
            @test isempty(Audit.audit_file(path; verbose = false).leaked)
        end
    end

    @testset "language constructs" begin
        # Every construct in the fixture is released except the two named here.
        # A failure names the construct that broke.
        report = Audit.audit_file(CONSTRUCTS; verbose = false)
        @test report.leaked == [:never_freed, :quoted_v]
    end

    @testset "does not parse" begin
        mktemp() do path, io
            write(io, "v = LibPETSc.VecCreateSeq(petsclib, comm, 10\nfunction broken(\n")
            close(io)
            # a syntax error must not silently read as a clean bill of health
            @test_logs (:warn, r"does not parse") Audit.audit_file(
                path;
                verbose = false,
            )
        end
    end

    @testset "helpers" begin
        @test Audit.audit_creator(:VecCreateSeq) == "Vec"
        @test Audit.audit_creator(:MatSeqAIJWithArrays) == "Mat"
        @test Audit.audit_creator(:DMStag) == "DM"
        @test Audit.audit_creator(:solve!) === nothing

        @test Audit.audit_destroyer(:destroy!)
        @test Audit.audit_destroyer(:destroy!)
        @test Audit.audit_destroyer(:VecDestroy)
        @test Audit.audit_destroyer(:finalizer)
        @test !Audit.audit_destroyer(:Destroy)
        @test !Audit.audit_destroyer(:assemble!)

        @test Audit.audit_callee(:(PETSc.LibPETSc.f(x))) === :f
        @test Audit.audit_callee(:(PETSc.f.(xs))) === :f
        @test !Audit.audit_isbroadcast(:(a.b))
        @test Audit.audit_isbroadcast(:(f.(x)))
        @test Audit.audit_argnames(:(destroy!(v))) == [:v]
        @test Audit.audit_argnames(:(destroy!(v...))) == [:v]
        @test Audit.audit_argnames(:(destroy!.([v, w]))) == [:v, :w]
        @test Audit.audit_argnames(:(VecDestroy(lib, v))) == [:lib, :v]
        @test Audit.audit_callee(:(f(x))) === :f
        @test Audit.audit_callee(:(x + 1)) === :+
    end
end

@testset "audit script command line" begin
    script = joinpath(@__DIR__, "..", "..", "scripts", "audit.jl")
    out = read(`$(Base.julia_cmd()) --startup-file=no $script $LEAKY`, String)
    @test occursin("POSSIBLY UNDESTROYED", out)
    @test occursin("  dm\n", out)
end
