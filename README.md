# PETSc
This package provides a high level interface for PETSc, enabling the use of PETSc as an `AbstractArray.  A low level interface is also available in the submodule PETSc.C.

This package requires the MPI.jl package be installed.  Once it is installed you should be able to run both Julia and Petsc in parallel using MPI for all communication.  The testing verifies that PETSc can be used both serially and in parallel.

To use the package, simply put `using PETSc` at the top of your Julia source file.  The module exports the names of all the functions, as well as the Petsc data type aliases and constants such as `PETSC_DECIDE`.

In general, it is possible to run PETSC in parallel. To do so with 4 processors, do:

```
mpirun -np 4 julia ./name_of_file
```

Note that this launches 4 independent Julia processes.  They are not aware of each other using Julia's built-in parallelism, and MPI is used for all communications.  

To run in serial, do:
```
julia ./name_of_file
```

Even when running serially, the MPI.jl package must be installed.


An example of using a Krylov Sub-Space method to solve a linear system is in  `test/test_ksp.jl`, which solves a simple system with a Krylov Subspace method and compares the result with a direct solve using Julia's backslash operator.  This works in serial and in parallel.  It requires some variables declared at the top of `runtests.jl` to work.



## To do:
  * Make the script for building Petsc more flexible, eg. allowing more configuration options like building BLAS or LAPCK, while ensure it remains completely autonomous (needed for Travis testing)
  * Wrap more KSP function

## Status
### Vector
  The AbstractArray for PetscVec is implemented.  Some additional PETSc 
  BLAS functions are wrapped as well.
### Matrix
 The AbstractArray interface for PetscMat is implemented.  Preallocation 
 is supported through optional keyword arguments to the matrix constructor or
 the `setpreallocation` function.  It possible to set multiple values in the 
  matrix without intermediate assembly using the `assemble` function or by 
 setting the `Mat` object field `assembling` to `false` and calling `setindex`
 repeatedly.

### KSP
 Just enough KSP functions are implimented to do a GMRES solve.  After the vector and matrix functions I will focus on KSP.

## Directory Structure
  `/src` : source files.  PETSc.jl is the main file containing initialization, with the functions for each type of Petsc object in its own file.  All constants are declared in `petsc_constants.jl`.

  `/src/generated`: auto generated wrappers from Clang.jl.  Not directly useful, but easy to modify to make useful

  `/test` : contains `runtest.jl`, which does some setup and runs all tests on all three version of Petsc currently supported.  Tests for each type of Petsc object (mirroring the files in `/src`) are contained in separate files.

  `/deps` : builds Petsc if needed.  See description below


## Building Petsc
Building the package will build build the 3 versions of PETSc in the `/deps` 
 directory, and writes the file `lib_locations.jl` to the `/src/generated` 
 directory to tell the package the location of the libraries.  Note that 
this builds the debug versions of PETSc, which are recommended to use for all 
development.  If you wish to do high performance computations, you should 
build the optimized versions of the library.  See the PETSc website for 
details.

## Installing MPI.jl
This package requires MPI.jl, although it is not listed in the REQUIRE file because that would download the release version of MPI.jl, which does not work.  Instead, you must use the master branch.  After you have an MPI implementation installed, Pkg.build("Petsc") will install it and then Petsc, according to the description above.  If you wish to install it manually, do:
```
  Pkg.clone("MPI")
  Pkg.build("MPI")
```

Currently, only MPI implimentations where the Fortran communicator is the same as the C communictor are supported.  This is due to the current implimentation of the MPI.jl package.  MPICH is one MPI implimentation that satisfies this requirement.  Note that OPENMPI does not.


## Auto Generation Notes
PETSc uses preprocessor variables to decide what code to include when compiling 
the library.  Clang does not know what preprocessor variables were defined at 
compile time, so it does not correctly detect the typealiases PetscScalar, PetscReal, etc.  To correctly autogenerate wrappers, the proper variables must be passed to Clang with the -D switch.  Note that users will not need to generate their own wrappers because they have already been generated and commit to the repo.

[![Build Status](https://travis-ci.org/JaredCrean2/PETSc.jl.svg?branch=master)](https://travis-ci.org/JaredCrean2/PETSc.jl)
