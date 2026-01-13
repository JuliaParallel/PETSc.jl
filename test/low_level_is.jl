using Test
using PETSc
using MPI

@testset "Low-level IS (Index Set) functions" begin
    petsclib = PETSc.getlib(PetscScalar=Float64)
    PETSc.initialize(petsclib)
    PetscInt = PETSc.LibPETSc.PetscInt
    
    @testset "ISCreateGeneral" begin
        indices = PetscInt[0, 2, 4, 6, 8]
        is = PETSc.LibPETSc.ISCreateGeneral(petsclib, MPI.COMM_SELF, 5, indices, PETSc.LibPETSc.PETSC_COPY_VALUES)
        @test is isa PETSc.LibPETSc.IS
        @test is.ptr != C_NULL
        
        # Get size
        size = PETSc.LibPETSc.ISGetSize(petsclib, is)
        @test size == 5
        
        local_size = PETSc.LibPETSc.ISGetLocalSize(petsclib, is)
        @test local_size == 5
        
        PETSc.LibPETSc.ISDestroy(petsclib, is)
        @test is.ptr == C_NULL
    end
    
    @testset "ISCreateStride" begin
        # Create stride: 0, 2, 4, 6, 8 (first=0, step=2, n=5)
        is = PETSc.LibPETSc.ISCreateStride(petsclib, MPI.COMM_SELF, 5, 0, 2)
        @test is isa PETSc.LibPETSc.IS
        
        size = PETSc.LibPETSc.ISGetSize(petsclib, is)
        @test size == 5
        
        PETSc.LibPETSc.ISDestroy(petsclib, is)
    end
    
    @testset "ISCreateBlock" begin
        # Create block IS
        indices = PetscInt[0, 2, 4]
        is = PETSc.LibPETSc.ISCreateBlock(petsclib, MPI.COMM_SELF, 2, 3, indices, PETSc.LibPETSc.PETSC_COPY_VALUES)
        @test is isa PETSc.LibPETSc.IS
        
        # Get block size
        bs = PETSc.LibPETSc.ISGetBlockSize(petsclib, is)
        @test bs == 2
        
        PETSc.LibPETSc.ISDestroy(petsclib, is)
    end
    
    @testset "IS properties" begin
        indices = PetscInt[0, 1, 2, 3, 4]
        is = PETSc.LibPETSc.ISCreateGeneral(petsclib, MPI.COMM_SELF, 5, indices, PETSc.LibPETSc.PETSC_COPY_VALUES)
        
        # Check if sorted
        sorted = PETSc.LibPETSc.ISSorted(petsclib, is)
        @test sorted == PETSc.LibPETSc.PETSC_TRUE
        
        # Check if identity
        is_identity = PETSc.LibPETSc.ISIdentity(petsclib, is)
        @test is_identity == PETSc.LibPETSc.PETSC_TRUE
        
        # Check if permutation
        is_perm = PETSc.LibPETSc.ISPermutation(petsclib, is)
        @test is_perm == PETSc.LibPETSc.PETSC_TRUE
        
        PETSc.LibPETSc.ISDestroy(petsclib, is)
    end
    
    @testset "IS operations" begin
        indices = PetscInt[4, 2, 8, 1, 5]
        is = PETSc.LibPETSc.ISCreateGeneral(petsclib, MPI.COMM_SELF, 5, indices, PETSc.LibPETSc.PETSC_COPY_VALUES)
        
        # Sort the IS
        @test_nowarn PETSc.LibPETSc.ISSort(petsclib, is)
        
        # Check it's sorted now
        sorted = PETSc.LibPETSc.ISSorted(petsclib, is)
        @test sorted == PETSc.LibPETSc.PETSC_TRUE
        
        # Get min/max
        min_val, max_val = PETSc.LibPETSc.ISGetMinMax(petsclib, is)
        @test min_val == 1
        @test max_val == 8
        
        # Duplicate
        is_copy = PETSc.LibPETSc.ISDuplicate(petsclib, is)
        @test is_copy isa PETSc.LibPETSc.IS
        @test is_copy.ptr != is.ptr
        
        PETSc.LibPETSc.ISDestroy(petsclib, is)
        PETSc.LibPETSc.ISDestroy(petsclib, is_copy)
    end
    
    @testset "IS with viewer convenience functions" begin
        indices = PetscInt[0, 2, 4, 6, 8]
        is = PETSc.LibPETSc.ISCreateGeneral(petsclib, MPI.COMM_SELF, 5, indices, PETSc.LibPETSc.PETSC_COPY_VALUES)
        
        # Redirect output to suppress IS viewer output during testing
        redirect_stdout(devnull) do
            redirect_stderr(devnull) do
                # Test PETSC_VIEWER_STDOUT_SELF
                viewer_stdout_self = PETSc.LibPETSc.PETSC_VIEWER_STDOUT_SELF(petsclib)
                @test viewer_stdout_self isa PETSc.LibPETSc.PetscViewer
                @test viewer_stdout_self != C_NULL
                @test PETSc.LibPETSc.ISView(petsclib, is, viewer_stdout_self) === nothing
                
                # Test PETSC_VIEWER_STDOUT_WORLD
                viewer_stdout_world = PETSc.LibPETSc.PETSC_VIEWER_STDOUT_WORLD(petsclib)
                @test viewer_stdout_world isa PETSc.LibPETSc.PetscViewer
                @test viewer_stdout_world != C_NULL
                @test PETSc.LibPETSc.ISView(petsclib, is, viewer_stdout_world) === nothing
                
                # Test PETSC_VIEWER_STDERR_SELF
                viewer_stderr_self = PETSc.LibPETSc.PETSC_VIEWER_STDERR_SELF(petsclib)
                @test viewer_stderr_self isa PETSc.LibPETSc.PetscViewer
                @test viewer_stderr_self != C_NULL
                @test PETSc.LibPETSc.ISView(petsclib, is, viewer_stderr_self) === nothing
                
                # Test PETSC_VIEWER_STDERR_WORLD
                viewer_stderr_world = PETSc.LibPETSc.PETSC_VIEWER_STDERR_WORLD(petsclib)
                @test viewer_stderr_world isa PETSc.LibPETSc.PetscViewer
                @test viewer_stderr_world != C_NULL
                @test PETSc.LibPETSc.ISView(petsclib, is, viewer_stderr_world) === nothing
            end
        end
        
        PETSc.LibPETSc.ISDestroy(petsclib, is)
    end
    
    PETSc.finalize(petsclib)
end
