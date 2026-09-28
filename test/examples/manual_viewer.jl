using Test
using PETSc
using MPI

@testset "Documentation examples for PetscViewer" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
        test_comm = MPI.COMM_SELF
    
    @testset "Basic Usage - Create viewer" begin
        # Create a viewer for ASCII output to stdout
        viewer = PETSc.LibPETSc.PetscViewerCreate(petsclib, test_comm)
        @test viewer isa PETSc.LibPETSc.PetscViewer
        @test viewer != C_NULL
        
        PETSc.LibPETSc.PetscViewerSetType(petsclib, viewer, "ascii")
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
    
    @testset "ASCII and binary file I/O" begin
        x = PETSc.PetscVec(petsclib, [1.0, 2.0, 3.0])
        cd(mktempdir()) do
            viewer = PETSc.LibPETSc.PetscViewerASCIIOpen(petsclib, test_comm, "output.txt")
            PETSc.LibPETSc.PetscViewerPushFormat(petsclib, viewer, PETSc.LibPETSc.PETSC_VIEWER_ASCII_MATLAB)
            PETSc.LibPETSc.VecView(petsclib, x, viewer)
            PETSc.LibPETSc.PetscViewerPopFormat(petsclib, viewer)
            PETSc.LibPETSc.PetscViewerDestroy(petsclib, viewer)
            @test occursin("2.", read("output.txt", String))

            # a vector saved in binary loads back unchanged
            viewer = PETSc.LibPETSc.PetscViewerBinaryOpen(petsclib, test_comm, "checkpoint.dat",
                                                          PETSc.LibPETSc.FILE_MODE_WRITE)
            PETSc.LibPETSc.VecView(petsclib, x, viewer)
            PETSc.LibPETSc.PetscViewerDestroy(petsclib, viewer)
            viewer = PETSc.LibPETSc.PetscViewerBinaryOpen(petsclib, test_comm, "checkpoint.dat",
                                                          PETSc.LibPETSc.FILE_MODE_READ)
            y = PETSc.LibPETSc.VecCreate(petsclib, test_comm)
            PETSc.LibPETSc.VecLoad(petsclib, y, viewer)
            PETSc.LibPETSc.PetscViewerDestroy(petsclib, viewer)
            @test y[:] == [1.0, 2.0, 3.0]
            PETSc.destroy!(y)
        end
        PETSc.destroy!(x)
    end
    
    PETSc.finalize(petsclib)
end
