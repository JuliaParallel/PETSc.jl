using Test
using PETSc
using MPI

@testset "Low-level Tao (optimization) functions" begin
    petsclib = PETSc.getlib(PetscScalar=Float64)
    PETSc.initialize(petsclib)
        test_comm = Sys.iswindows() ? LibPETSc.PETSC_COMM_SELF : MPI.COMM_SELF
    
    @testset "Tao object creation and destruction" begin
        tao = PETSc.LibPETSc.TaoCreate(petsclib, test_comm)
        @test tao isa PETSc.LibPETSc.Tao
        @test tao.ptr != C_NULL
        
        PETSc.LibPETSc.TaoDestroy(petsclib, tao)
        @test tao.ptr == C_NULL
    end
    
    @testset "Tao solver type" begin
        tao = PETSc.LibPETSc.TaoCreate(petsclib, test_comm)
        
        # Set Tao type (using NLS - Newton Line Search, always available)
        @test_nowarn PETSc.LibPETSc.TaoSetType(petsclib, tao, Base.unsafe_convert(Ptr{Int8}, "nls"))
        
        # Get Tao type back
        taotype = PETSc.LibPETSc.TaoGetType(petsclib, tao)
        @test taotype == "nls"
        
        PETSc.LibPETSc.TaoDestroy(petsclib, tao)
    end
    
    @testset "Tao with different solver types" begin
        # Test various Tao solver types (using types that don't require special builds)
        for taotype in ["nls", "ntr", "cg", "nm"]
            tao = PETSc.LibPETSc.TaoCreate(petsclib, test_comm)
            @test_nowarn PETSc.LibPETSc.TaoSetType(petsclib, tao, Base.unsafe_convert(Ptr{Int8}, taotype))
            retrieved_type = PETSc.LibPETSc.TaoGetType(petsclib, tao)
            @test retrieved_type == taotype
            PETSc.LibPETSc.TaoDestroy(petsclib, tao)
        end
    end
    
    @testset "Tao tolerances" begin
        tao = PETSc.LibPETSc.TaoCreate(petsclib, test_comm)
        PETSc.LibPETSc.TaoSetType(petsclib, tao, Base.unsafe_convert(Ptr{Int8}, "nls"))
        
        # Set tolerances (gatol, grtol, gttol)
        @test_nowarn PETSc.LibPETSc.TaoSetTolerances(petsclib, tao, 1e-8, 1e-8, 1e-8)
        
        # Get tolerances back
        gatol, grtol, gttol = PETSc.LibPETSc.TaoGetTolerances(petsclib, tao)
        @test gatol ≈ 1e-8
        @test grtol ≈ 1e-8
        @test gttol ≈ 1e-8
        
        PETSc.LibPETSc.TaoDestroy(petsclib, tao)
    end
    
    @testset "Tao max iterations and function evaluations" begin
        tao = PETSc.LibPETSc.TaoCreate(petsclib, test_comm)
        
        # Set max iterations
        @test_nowarn PETSc.LibPETSc.TaoSetMaximumIterations(petsclib, tao, 100)
        
        # Get max iterations back
        maxits = PETSc.LibPETSc.TaoGetMaximumIterations(petsclib, tao)
        @test maxits == 100
        
        # Set max function evaluations
        @test_nowarn PETSc.LibPETSc.TaoSetMaximumFunctionEvaluations(petsclib, tao, 1000)
        
        # Get max function evaluations back
        maxfuncs = PETSc.LibPETSc.TaoGetMaximumFunctionEvaluations(petsclib, tao)
        @test maxfuncs == 1000
        
        PETSc.LibPETSc.TaoDestroy(petsclib, tao)
    end
    
    @testset "Tao iteration count" begin
        tao = PETSc.LibPETSc.TaoCreate(petsclib, test_comm)
        PETSc.LibPETSc.TaoSetType(petsclib, tao, Base.unsafe_convert(Ptr{Int8}, "nls"))
        
        # Get iteration count (should be 0 before solving)
        iter = PETSc.LibPETSc.TaoGetIterationNumber(petsclib, tao)
        @test iter == 0
        
        PETSc.LibPETSc.TaoDestroy(petsclib, tao)
    end
    
    PETSc.finalize(petsclib)
end
