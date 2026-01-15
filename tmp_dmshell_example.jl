using PETSc, MPI

MPI.Init()
try
    petsclib = PETSc.getlib()
    PETSc.initialize(petsclib)

    shell = LibPETSc.DMShellCreate(petsclib, MPI.COMM_WORLD)
    println("Created shell: ", shell)

    LibPETSc.DMSetUp(petsclib, shell)
    println("DMShell set up successfully")

    LibPETSc.DMDestroy(petsclib, shell)

    PETSc.finalize(petsclib)
    println("DMShell example completed")
catch e
    @error "DMShell example failed" exception=(e, catch_backtrace())
    try PETSc.finalize(petsclib) catch end
    rethrow(e)
finally
    try MPI.Finalize() catch end
end