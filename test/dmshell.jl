using Test
using PETSc, MPI

MPI.Init()
try
    petsclib = PETSc.getlib()
    PETSc.initialize(petsclib)

    @testset "DMShell basic" begin
        shell = LibPETSc.DMShellCreate(petsclib, MPI.COMM_WORLD)
        @test shell.ptr != C_NULL

        LibPETSc.DMSetUp(petsclib, shell)

        LibPETSc.DMDestroy(petsclib, shell)
    end

    PETSc.finalize(petsclib)
catch e
    try PETSc.finalize(petsclib) catch end
    rethrow(e)
finally
    try MPI.Finalize() catch end
end