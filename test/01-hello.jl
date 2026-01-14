using MPI
using PETSc
MPI.Init()

# Windows PETSc binaries are built without MPI support, use PETSC_COMM_SELF instead
comm = Sys.iswindows() ? LibPETSc.PETSC_COMM_SELF : MPI.COMM_WORLD
println("Hello world, I am $(MPI.Comm_rank(comm)) of $(MPI.Comm_size(comm))")
MPI.Barrier(comm)