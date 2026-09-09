using Test
using PETSc, MPI

if !Sys.iswindows()
    MPI.Initialized() || MPI.Init()
end

@testset "argument validation" begin
    @testset "PetscNotInitialized" begin
        e = PETSc.PetscNotInitialized(PETSc.petsclibs[1])
        msg = sprint(showerror, e)
        # names the library by its element types, not by dumping the handle
        @test occursin("PetscScalar", msg)
        @test occursin("PetscInt", msg)
        @test occursin("initialize", msg)
        # a handle that answers neither accessor still produces a message
        @test occursin("not initialized", sprint(showerror, PETSc.PetscNotInitialized("x")))
    end

    @testset "parse_options" begin
        @test PETSc.parse_options(["-ksp_monitor", "-pc_type", "mg"]) ==
              (ksp_monitor = nothing, pc_type = "mg")

        # each of these used to trip an @assert, and the empty string raised a
        # BoundsError before the length was checked first
        for bad in (["notanoption"], ["-"], [""])
            @test_throws ArgumentError PETSc.parse_options(bad)
        end
    end

    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt

        @testset "element type mismatch ($PetscScalar)" begin
            # a Julia array whose element type differs from the library's
            wrong = PetscScalar === Float64 ? Float32 : Float64
            @test_throws ArgumentError PETSc.VecSeq(petsclib, wrong[1, 2, 3])
            @test_throws ArgumentError PETSc.MatSeqDense(petsclib, wrong[1 2; 3 4])
        end

        @testset "size mismatch ($PetscScalar)" begin
            # nonzeros shorter than the number of rows
            @test_throws DimensionMismatch PETSc.MatSeqAIJ(
                petsclib,
                PetscInt(4),
                PetscInt(4),
                PetscInt[1, 1],
            )
        end

        PETSc.finalize(petsclib)
    end
end
