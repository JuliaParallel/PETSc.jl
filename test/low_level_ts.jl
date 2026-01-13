using Test
using PETSc
using MPI

@testset "Low-level TS (Time Stepping) functions" begin
    petsclib = PETSc.getlib(PetscScalar=Float64)
    PETSc.initialize(petsclib)
    
    @testset "TS object creation and destruction" begin
        ts = PETSc.LibPETSc.TSCreate(petsclib, MPI.COMM_SELF)
        @test ts isa PETSc.LibPETSc.TS
        @test ts.ptr != C_NULL
        
        PETSc.LibPETSc.TSDestroy(petsclib, ts)
        @test ts.ptr == C_NULL
    end
    
    @testset "TS problem type and solver type" begin
        ts = PETSc.LibPETSc.TSCreate(petsclib, MPI.COMM_SELF)
        
        # Set problem type
        @test_nowarn PETSc.LibPETSc.TSSetProblemType(petsclib, ts, PETSc.LibPETSc.TS_LINEAR)
        
        # Set TS type
        @test_nowarn PETSc.LibPETSc.TSSetType(petsclib, ts, Base.unsafe_convert(Ptr{Int8}, "bdf"))
        
        # Get TS type back
        tstype = PETSc.LibPETSc.TSGetType(petsclib, ts)
        @test tstype == "bdf"
        
        PETSc.LibPETSc.TSDestroy(petsclib, ts)
    end
    
    @testset "TS time parameters" begin
        ts = PETSc.LibPETSc.TSCreate(petsclib, MPI.COMM_SELF)
        PETSc.LibPETSc.TSSetType(petsclib, ts, Base.unsafe_convert(Ptr{Int8}, "bdf"))
        
        # Set time parameters
        @test_nowarn PETSc.LibPETSc.TSSetTime(petsclib, ts, 0.0)
        @test_nowarn PETSc.LibPETSc.TSSetMaxTime(petsclib, ts, 1.0)
        @test_nowarn PETSc.LibPETSc.TSSetTimeStep(petsclib, ts, 0.1)
        
        # Get time parameters back
        current_time = PETSc.LibPETSc.TSGetTime(petsclib, ts)
        @test current_time ≈ 0.0
        
        time_step = PETSc.LibPETSc.TSGetTimeStep(petsclib, ts)
        @test time_step ≈ 0.1
        
        PETSc.LibPETSc.TSDestroy(petsclib, ts)
    end
    
    @testset "TS with different solver types" begin
        for tstype in ["euler", "bdf", "rk"]
            ts = PETSc.LibPETSc.TSCreate(petsclib, MPI.COMM_SELF)
            @test_nowarn PETSc.LibPETSc.TSSetType(petsclib, ts, Base.unsafe_convert(Ptr{Int8}, tstype))
            retrieved_type = PETSc.LibPETSc.TSGetType(petsclib, ts)
            @test retrieved_type == tstype
            PETSc.LibPETSc.TSDestroy(petsclib, ts)
        end
    end
    
    PETSc.finalize(petsclib)
end
