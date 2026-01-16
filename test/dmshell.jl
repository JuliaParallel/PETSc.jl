using Test
using PETSc, MPI

# Initialize MPI & PETSc only if needed
 mpi_started = false
 if !MPI.Initialized()
     MPI.Init()
     mpi_started = true
 end
 try
     local petsclib = PETSc.getlib()
     petsc_started = false
     if !PETSc.initialized(petsclib)
         PETSc.initialize(petsclib)
         petsc_started = true
     end

    @testset "DMShell basic" begin
        shell = LibPETSc.DMShellCreate(petsclib, MPI.COMM_WORLD)
        @test shell.ptr != C_NULL

        LibPETSc.DMSetUp(petsclib, shell)

        LibPETSc.DMDestroy(petsclib, shell)
    end

    if petsc_started
        PETSc.finalize(petsclib)
    end
catch e
    try if petsc_started; PETSc.finalize(petsclib); end catch end
    rethrow(e)
finally
    # Do not call MPI.Finalize() from tests; test harness manages MPI lifecycle
    nothing
end