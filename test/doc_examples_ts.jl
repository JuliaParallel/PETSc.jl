using Test
using PETSc
using MPI

@testset "Documentation examples for TS" begin
    petsclib = PETSc.getlib(PetscScalar=Float64)
    PETSc.initialize(petsclib)
        test_comm = Sys.iswindows() ? LibPETSc.PETSC_COMM_SELF : MPI.COMM_SELF
    
    @testset "Basic Usage" begin
        # Create a TS object
        ts = PETSc.LibPETSc.TSCreate(petsclib, test_comm)
        @test ts isa PETSc.LibPETSc.TS
        @test ts.ptr != C_NULL
        
        # Set the problem type (ODE or DAE)
        PETSc.LibPETSc.TSSetProblemType(petsclib, ts, PETSc.LibPETSc.TS_NONLINEAR)
        
        # Set the time stepping method (e.g., BDF, RK, Theta)
        PETSc.LibPETSc.TSSetType(petsclib, ts, "bdf")  # String convenience wrapper
        
        # Set time span
        PETSc.LibPETSc.TSSetTime(petsclib, ts, 0.0)  # Initial time
        PETSc.LibPETSc.TSSetMaxTime(petsclib, ts, 1.0)  # Final time
        PETSc.LibPETSc.TSSetExactFinalTime(petsclib, ts, PETSc.LibPETSc.TS_EXACTFINALTIME_STEPOVER)
        
        # Set initial time step
        PETSc.LibPETSc.TSSetTimeStep(petsclib, ts, 0.01)
        
        # Set options from command line/options database
        PETSc.LibPETSc.TSSetFromOptions(petsclib, ts)
        
        # Get solution time (returns value directly, not through Ref)
        final_time = PETSc.LibPETSc.TSGetSolveTime(petsclib, ts)
        @test final_time isa Float64
        
        # Cleanup
        PETSc.LibPETSc.TSDestroy(petsclib, ts)
        @test ts.ptr == C_NULL
    end
    
    @testset "Different Time Integration Schemes" begin
        # Test different TS types mentioned in docs
        ts_types = ["euler", "rk", "bdf", "theta", "arkimex", "ssp"]
        
        for ts_type in ts_types
            ts = PETSc.LibPETSc.TSCreate(petsclib, test_comm)
            
            # Set type using string
            PETSc.LibPETSc.TSSetType(petsclib, ts, ts_type)
            
            # Verify type was set (TSGetType already returns a String)
            type_str = PETSc.LibPETSc.TSGetType(petsclib, ts)
            @test type_str == ts_type
            
            PETSc.LibPETSc.TSDestroy(petsclib, ts)
        end
    end
    
    @testset "Time Parameter Setters/Getters" begin
        ts = PETSc.LibPETSc.TSCreate(petsclib, test_comm)
        
        # Set and get time
        PETSc.LibPETSc.TSSetTime(petsclib, ts, 0.5)
        time = PETSc.LibPETSc.TSGetTime(petsclib, ts)
        @test time ≈ 0.5
        
        # Set and get time step
        PETSc.LibPETSc.TSSetTimeStep(petsclib, ts, 0.01)
        dt = PETSc.LibPETSc.TSGetTimeStep(petsclib, ts)
        @test dt ≈ 0.01
        
        # Set max time
        PETSc.LibPETSc.TSSetMaxTime(petsclib, ts, 10.0)
        max_time = PETSc.LibPETSc.TSGetMaxTime(petsclib, ts)
        @test max_time ≈ 10.0
        
        # Set max steps
        PETSc.LibPETSc.TSSetMaxSteps(petsclib, ts, 1000)
        max_steps = PETSc.LibPETSc.TSGetMaxSteps(petsclib, ts)
        @test max_steps == 1000
        
        PETSc.LibPETSc.TSDestroy(petsclib, ts)
    end
    
    PETSc.finalize(petsclib)
end
