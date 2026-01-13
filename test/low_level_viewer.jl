using Test
using PETSc
using MPI

@testset "Low-level PetscViewer convenience functions" begin
    petsclib = PETSc.getlib(PetscScalar=Float64)
    PETSc.initialize(petsclib)
    
    @testset "PETSC_VIEWER_STDOUT_SELF" begin
        viewer = PETSc.LibPETSc.PETSC_VIEWER_STDOUT_SELF(petsclib)
        @test viewer isa PETSc.LibPETSc.PetscViewer
        @test viewer != C_NULL
    end
    
    @testset "PETSC_VIEWER_STDOUT_WORLD" begin
        viewer = PETSc.LibPETSc.PETSC_VIEWER_STDOUT_WORLD(petsclib)
        @test viewer isa PETSc.LibPETSc.PetscViewer
        @test viewer != C_NULL
    end
    
    @testset "PETSC_VIEWER_STDERR_SELF" begin
        viewer = PETSc.LibPETSc.PETSC_VIEWER_STDERR_SELF(petsclib)
        @test viewer isa PETSc.LibPETSc.PetscViewer
        @test viewer != C_NULL
    end
    
    @testset "PETSC_VIEWER_STDERR_WORLD" begin
        viewer = PETSc.LibPETSc.PETSC_VIEWER_STDERR_WORLD(petsclib)
        @test viewer isa PETSc.LibPETSc.PetscViewer
        @test viewer != C_NULL
    end
    
    @testset "Using viewers with ISView" begin
        PetscInt = PETSc.LibPETSc.PetscInt
        is = PETSc.LibPETSc.ISCreateGeneral(petsclib, MPI.COMM_SELF, 5, PetscInt[0,2,4,6,8], PETSc.LibPETSc.PETSC_COPY_VALUES)
        
        # Redirect output to suppress IS viewer output during testing
        redirect_stdout(devnull) do
            redirect_stderr(devnull) do
                # Test STDOUT viewer - it will print to stdout
                viewer_out = PETSc.LibPETSc.PETSC_VIEWER_STDOUT_SELF(petsclib)
                @test PETSc.LibPETSc.ISView(petsclib, is, viewer_out) === nothing
                
                # Test STDERR viewer - it will print to stderr
                viewer_err = PETSc.LibPETSc.PETSC_VIEWER_STDERR_SELF(petsclib)
                @test PETSc.LibPETSc.ISView(petsclib, is, viewer_err) === nothing
            end
        end
        
        PETSc.LibPETSc.ISDestroy(petsclib, is)
    end
    
    PETSc.finalize(petsclib)
end
