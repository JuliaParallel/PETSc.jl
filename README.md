# PETSc
This package provides thin wrappers for PETSc, as well as a few convience functions that take advantage of multiple dispatch.

This package requires the MPI.jl package be installed.  Once it is installed you should be able to run both Julia and Petsc in parallel using MPI for all communcation.  Initial some initial experiments indicate this works, although it is not thoroughly tested yet.

To run in parallel with 4 processors, do:

mpirun -np 4 julia ./name_of_file

To run in serial, do:

julia ./name_of_file

The only currently working example is  test/exKSP3.jl, which solves a simple system with a Krylov Subspace method and compares the result with a direct solve using Julia's backslash operator



To do:
  * Use Petsc typealiases for all arguments
  * Generate file at installation that defines the correct typealiases, taking into account the options Petsc was built with (the Petsc function PetscDataTypeGetSize()  should help)
  * Handle Petsc error codes properly
  * Make the script for building Petsc more flexible, eg. allowing more configuration options like building blas or lapack, while ensure it remains completely autonomous (needed for Travis testig)
  * Wrap more PetscVec functions
  * Wrap more PetscMat functions
  * Wrap more KSP function
  * Determine priorities for wrapping additional functions



  
