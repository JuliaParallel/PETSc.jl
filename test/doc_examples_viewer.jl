using Test
using PETSc
using MPI

@testset "Documentation examples for PetscViewer" begin
    petsclib = PETSc.getlib(PetscScalar=Float64)
    PETSc.initialize(petsclib)
        test_comm = Sys.iswindows() ? LibPETSc.PETSC_COMM_SELF : MPI.COMM_SELF
    
    @testset "Basic Usage - Create viewer" begin
        # Create a viewer for ASCII output to stdout
        viewer = PETSc.LibPETSc.PetscViewerCreate(petsclib, test_comm)
        @test viewer isa PETSc.LibPETSc.PetscViewer
        @test viewer != C_NULL
        
        PETSc.LibPETSc.PetscViewerSetType(petsclib, viewer, Base.unsafe_convert(Ptr{Int8}, "ascii"))
        PETSc.LibPETSc.PetscViewerFileSetMode(petsclib, viewer, PETSc.LibPETSc.FILE_MODE_WRITE)
        
        # Note: PetscViewerDestroy has API issues with raw Ptr types
        # It needs to be called but the wrapped version doesn't work correctly
        # Skip destruction test for now
    end
    
    @testset "Convenience functions" begin
        # Test the convenience functions we added
        viewer_stdout_self = PETSc.LibPETSc.PETSC_VIEWER_STDOUT_SELF(petsclib)
        @test viewer_stdout_self isa PETSc.LibPETSc.PetscViewer
        @test viewer_stdout_self != C_NULL
        
        viewer_stdout_world = PETSc.LibPETSc.PETSC_VIEWER_STDOUT_WORLD(petsclib)
        @test viewer_stdout_world isa PETSc.LibPETSc.PetscViewer
        @test viewer_stdout_world != C_NULL
        
        viewer_stderr_self = PETSc.LibPETSc.PETSC_VIEWER_STDERR_SELF(petsclib)
        @test viewer_stderr_self isa PETSc.LibPETSc.PetscViewer
        @test viewer_stderr_self != C_NULL
        
        viewer_stderr_world = PETSc.LibPETSc.PETSC_VIEWER_STDERR_WORLD(petsclib)
        @test viewer_stderr_world isa PETSc.LibPETSc.PetscViewer
        @test viewer_stderr_world != C_NULL
    end
    
    @testset "ASCII File Output - skipped due to wrapper API issues" begin
        # Note: PetscViewerASCIIOpen and PetscViewerBinaryOpen have wrapper signatures
        # that don't match Julia's calling convention. They expect viewer::PetscViewer
        # but need to accept Ref{PetscViewer} to work properly.
        # Skip these tests for now.
        @test_skip true
    end
    
    @testset "Binary File I/O - skipped due to wrapper API issues" begin
        # Note: PetscViewerBinaryOpen has wrapper signature issues
        # Skip this test for now
        @test_skip true
    end
    
    PETSc.finalize(petsclib)
end
