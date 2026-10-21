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
        @test_throws PETSc.PetscNotInitialized PETSc.PetscOptions(petsclib; ksp_monitor = "")
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
            @test_throws "array has element type" PETSc.PetscVec(petsclib, MPI.COMM_SELF, wrong[1, 2, 3])
            @test_throws ArgumentError PETSc.PetscVec(petsclib, MPI.COMM_SELF, wrong[1, 2, 3])
        end

        @testset "DMs have 1 to 3 dimensions ($PetscScalar)" begin
            none = ntuple(_ -> PETSc.DM_BOUNDARY_NONE, 4)
            @test_throws "a DMDA has 1, 2 or 3 dimensions, got 4" PETSc.DMDA(
                petsclib, MPI.COMM_SELF, none, (3, 3, 3, 3), 1, 1,
            )
            @test_throws "a DMStag has 1, 2 or 3 dimensions, got 4" PETSc.DMStag(
                petsclib, MPI.COMM_SELF, none, (3, 3, 3, 3), (1, 0, 0, 0, 0), 1,
            )
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

@testset "PetscError names the C function and PETSc's message" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    PetscInt = petsclib.PetscInt
    LibPETSc = PETSc.LibPETSc

    v = PETSc.PetscVec(petsclib, 3)
    err = try
        LibPETSc.VecSetSizes(petsclib, v, PetscInt(4), PetscInt(4))
    catch e
        e
    end
    @test err isa LibPETSc.PetscError
    @test err.code == 56               # PETSC_ERR_SUP
    @test err.func === :VecSetSizes
    @test occursin("Cannot change/reset vector sizes", err.message)
    @test_throws "PetscError(56) from VecSetSizes: No support for this operation" throw(err)

    # An error raised deep inside PETSc is named after the call Julia made
    err = try
        PETSc.DMDA(petsclib, MPI.COMM_SELF, (PETSc.DM_BOUNDARY_NONE,), (10,), 1, 1; da_refine = "abc")
    catch e
        e
    end
    @test err.func === :DMSetFromOptions
    @test occursin("Input string abc has no integer value", err.message)

    PETSc.destroy!(v)
    PETSc.finalize(petsclib)

    # The message is captured when the error is raised, so it still prints
    # once the library is finalized
    @test occursin("Input string abc", sprint(showerror, err))
    # A code alone prints as it always has
    @test sprint(showerror, LibPETSc.PetscError(63)) == "PetscError(63)"
    @test LibPETSc.PetscError(63).code == 63
end
