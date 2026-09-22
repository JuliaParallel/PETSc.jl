using Test
using PETSc, MPI
using SparseArrays: sparse

MPI.Initialized() || MPI.Init()

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

    # every high-level constructor that creates a PETSc object checks the library first,
    # instead of failing inside PETSc with a bare PetscError
    @testset "constructors on an uninitialized library" begin
        petsclib = PETSc.petsclibs[1]
        PetscScalar = petsclib.PetscScalar
        comm = MPI.COMM_SELF
        B = PETSc.DM_BOUNDARY_NONE

        # a DMStag to derive from, created before the library goes down
        PETSc.initialize(petsclib)
        stag = PETSc.DMStag(petsclib, comm, (B,), (5,), (1, 1), 1)
        PETSc.finalize(petsclib)

        S = sparse(PetscScalar[1 0; 0 1])
        @test_throws PETSc.PetscNotInitialized PETSc.PetscVec(petsclib, 3)
        @test_throws PETSc.PetscNotInitialized PETSc.PetscVec(petsclib, comm, PetscScalar[1, 2])
        @test_throws PETSc.PetscNotInitialized PETSc.PetscMat(petsclib, S)
        @test_throws PETSc.PetscNotInitialized PETSc.PetscMat(petsclib, comm, S; with_arrays = true)
        @test_throws PETSc.PetscNotInitialized PETSc.PetscMat(petsclib, [0, 1, 2], [0, 1], PetscScalar[1, 1])
        @test_throws PETSc.PetscNotInitialized PETSc.KSP(petsclib, comm, S)
        @test_throws PETSc.PetscNotInitialized PETSc.SNES(petsclib, comm)
        @test_throws PETSc.PetscNotInitialized PETSc.TS(petsclib, comm)
        @test_throws PETSc.PetscNotInitialized PETSc.DMDA(petsclib, comm, (B,), (5,), 1, 1)
        @test_throws PETSc.PetscNotInitialized PETSc.DMStag(petsclib, comm, (B,), (5,), (1, 1), 1)
        @test_throws PETSc.PetscNotInitialized PETSc.DMStag(stag, (1, 1))
        @test_throws PETSc.PetscNotInitialized PETSc.DMPlex(petsclib, comm)
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
            @test_throws ArgumentError PETSc.PetscVec(petsclib, wrong[1, 2, 3])
            @test_throws ArgumentError PETSc.PetscMat(petsclib, wrong[1 2; 3 4])
        end

        @testset "size mismatch ($PetscScalar)" begin
            # nonzeros shorter than the number of rows
            @test_throws DimensionMismatch PETSc.PetscMat(
                petsclib,
                PetscInt(4),
                PetscInt(4),
                PetscInt[1, 1],
            )
        end

        PETSc.finalize(petsclib)
    end
end
