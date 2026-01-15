using Test
using PETSc
using MPI

@testset "Documentation examples for Tao" begin
    petsclib = PETSc.getlib(PetscScalar=Float64)
    PETSc.initialize(petsclib)
        test_comm = Sys.iswindows() ? LibPETSc.PETSC_COMM_SELF : MPI.COMM_SELF
    PetscInt = PETSc.LibPETSc.PetscInt
    
    @testset "Basic Usage" begin
        # Create a Tao object
        tao = PETSc.LibPETSc.TaoCreate(petsclib, test_comm)
        @test tao isa PETSc.LibPETSc.Tao
        @test tao.ptr != C_NULL
        
        # Set the optimization algorithm (e.g., LMVM, BLMVM, NLS)
        PETSc.LibPETSc.TaoSetType(petsclib, tao, Base.unsafe_convert(Ptr{Int8}, "lmvm"))
        
        # Set convergence tolerances
        PETSc.LibPETSc.TaoSetTolerances(petsclib, tao, 1e-8, 1e-8, 1e-8)
        
        # Set maximum iterations
        PETSc.LibPETSc.TaoSetMaximumIterations(petsclib, tao, 1000)
        
        # Set options from command line/options database
        PETSc.LibPETSc.TaoSetFromOptions(petsclib, tao)
        
        # Get iteration information (returns value directly)
        iter = PETSc.LibPETSc.TaoGetIterationNumber(petsclib, tao)
        @test iter == 0  # Should be 0 before solving
        
        # Cleanup
        PETSc.LibPETSc.TaoDestroy(petsclib, tao)
        @test tao.ptr == C_NULL
    end
    
    @testset "Different Optimization Algorithms" begin
        # Test unconstrained algorithms
        unconstrained_types = ["lmvm", "cg", "nm", "nls", "ntr"]
        
        for tao_type in unconstrained_types
            tao = PETSc.LibPETSc.TaoCreate(petsclib, test_comm)
            
            # Set type using string
            PETSc.LibPETSc.TaoSetType(petsclib, tao, Base.unsafe_convert(Ptr{Int8}, tao_type))
            
            # Verify type was set (TaoGetType already returns a String)
            type_str = PETSc.LibPETSc.TaoGetType(petsclib, tao)
            @test type_str == tao_type
            
            PETSc.LibPETSc.TaoDestroy(petsclib, tao)
        end
        
        # Test bound-constrained algorithms
        bound_types = ["blmvm", "bncg", "bqnls", "bntl", "tron"]
        
        for tao_type in bound_types
            tao = PETSc.LibPETSc.TaoCreate(petsclib, test_comm)
            
            # Set type using string
            PETSc.LibPETSc.TaoSetType(petsclib, tao, Base.unsafe_convert(Ptr{Int8}, tao_type))
            
            # Verify type was set (TaoGetType already returns a String)
            type_str = PETSc.LibPETSc.TaoGetType(petsclib, tao)
            @test type_str == tao_type
            
            PETSc.LibPETSc.TaoDestroy(petsclib, tao)
        end
    end
    
    @testset "Tolerance and Iteration Settings" begin
        tao = PETSc.LibPETSc.TaoCreate(petsclib, test_comm)
        PETSc.LibPETSc.TaoSetType(petsclib, tao, Base.unsafe_convert(Ptr{Int8}, "lmvm"))
        
        # Set convergence tolerances (gatol, grtol, gttol)
        PETSc.LibPETSc.TaoSetTolerances(petsclib, tao, 1e-6, 1e-6, 1e-6)
        
        # Get tolerances back
        gatol, grtol, gttol = PETSc.LibPETSc.TaoGetTolerances(petsclib, tao)
        @test gatol ≈ 1e-6
        @test grtol ≈ 1e-6
        @test gttol ≈ 1e-6
        
        # Set and get maximum iterations
        PETSc.LibPETSc.TaoSetMaximumIterations(petsclib, tao, 500)
        max_iter = PETSc.LibPETSc.TaoGetMaximumIterations(petsclib, tao)
        @test max_iter == 500
        
        # Get current iteration number
        iter = PETSc.LibPETSc.TaoGetIterationNumber(petsclib, tao)
        @test iter == 0  # Should be 0 before solving
        
        # Set and get maximum function evaluations
        PETSc.LibPETSc.TaoSetMaximumFunctionEvaluations(petsclib, tao, 10000)
        max_funcs = PETSc.LibPETSc.TaoGetMaximumFunctionEvaluations(petsclib, tao)
        @test max_funcs == 10000
        
        PETSc.LibPETSc.TaoDestroy(petsclib, tao)
    end
    
    @testset "Solution Vector Setup" begin
        tao = PETSc.LibPETSc.TaoCreate(petsclib, test_comm)
        PETSc.LibPETSc.TaoSetType(petsclib, tao, Base.unsafe_convert(Ptr{Int8}, "lmvm"))
        
        # Create a solution vector
        x = PETSc.LibPETSc.VecCreate(petsclib, test_comm)
        PETSc.LibPETSc.VecSetSizes(petsclib, x, 10, 10)
        PETSc.LibPETSc.VecSetFromOptions(petsclib, x)
        PETSc.LibPETSc.VecSet(petsclib, x, 0.0)
        
        # Set initial solution
        PETSc.LibPETSc.TaoSetSolution(petsclib, tao, x)
        
        # Note: TaoGetSolution modifies the Vec parameter in-place
        # Skip this test as the API pattern is complex
        
        PETSc.LibPETSc.VecDestroy(petsclib, x)
        PETSc.LibPETSc.TaoDestroy(petsclib, tao)
    end
    
    PETSc.finalize(petsclib)
end
