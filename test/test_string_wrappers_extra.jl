using Test
using PETSc

@testset "String wrappers for TS subtype setters" begin
    petsclib = PETSc.getlib(PetscScalar = Float64)
    PETSc.initialize(petsclib)

    @testset "TSRKSetType(::String) round-trips" begin
        for sub in ["3bs", "5dp", "4"]
            ts = PETSc.LibPETSc.TSCreate(petsclib, PETSc.LibPETSc.PETSC_COMM_SELF)
            PETSc.LibPETSc.TSSetType(petsclib, ts, "rk")
            @test_nowarn PETSc.LibPETSc.TSRKSetType(petsclib, ts, sub)
            @test PETSc.LibPETSc.TSRKGetType(petsclib, ts) == sub
            PETSc.LibPETSc.TSDestroy(petsclib, ts)
        end
    end

    @testset "TSRKSetType(::String) rejects unknown subtype" begin
        ts = PETSc.LibPETSc.TSCreate(petsclib, PETSc.LibPETSc.PETSC_COMM_SELF)
        PETSc.LibPETSc.TSSetType(petsclib, ts, "rk")
        @test_throws Exception PETSc.LibPETSc.TSRKSetType(petsclib, ts, "this-subtype-does-not-exist")
        PETSc.LibPETSc.TSDestroy(petsclib, ts)
    end

    @testset "TSRosWSetType(::String) round-trips" begin
        for sub in ["ra34pw2", "rodas3", "2m"]
            ts = PETSc.LibPETSc.TSCreate(petsclib, PETSc.LibPETSc.PETSC_COMM_SELF)
            PETSc.LibPETSc.TSSetType(petsclib, ts, "rosw")
            @test_nowarn PETSc.LibPETSc.TSRosWSetType(petsclib, ts, sub)
            @test PETSc.LibPETSc.TSRosWGetType(petsclib, ts) == sub
            PETSc.LibPETSc.TSDestroy(petsclib, ts)
        end
    end

    @testset "TSRosWSetType(::String) rejects unknown subtype" begin
        ts = PETSc.LibPETSc.TSCreate(petsclib, PETSc.LibPETSc.PETSC_COMM_SELF)
        PETSc.LibPETSc.TSSetType(petsclib, ts, "rosw")
        @test_throws Exception PETSc.LibPETSc.TSRosWSetType(petsclib, ts, "this-subtype-does-not-exist")
        PETSc.LibPETSc.TSDestroy(petsclib, ts)
    end

    @testset "TSARKIMEXSetType(::String) round-trips" begin
        for sub in ["2e", "3", "4"]
            ts = PETSc.LibPETSc.TSCreate(petsclib, PETSc.LibPETSc.PETSC_COMM_SELF)
            PETSc.LibPETSc.TSSetType(petsclib, ts, "arkimex")
            @test_nowarn PETSc.LibPETSc.TSARKIMEXSetType(petsclib, ts, sub)
            @test PETSc.LibPETSc.TSARKIMEXGetType(petsclib, ts) == sub
            PETSc.LibPETSc.TSDestroy(petsclib, ts)
        end
    end

    @testset "TSARKIMEXSetType(::String) rejects unknown subtype" begin
        ts = PETSc.LibPETSc.TSCreate(petsclib, PETSc.LibPETSc.PETSC_COMM_SELF)
        PETSc.LibPETSc.TSSetType(petsclib, ts, "arkimex")
        @test_throws Exception PETSc.LibPETSc.TSARKIMEXSetType(petsclib, ts, "this-subtype-does-not-exist")
        PETSc.LibPETSc.TSDestroy(petsclib, ts)
    end

    PETSc.finalize(petsclib)
end
