# PETSc

[![Build Status](https://travis-ci.org/JuliaParallel/PETSc.jl.svg?branch=master)](https://travis-ci.org/JuliaParallel/PETSc.jl)
[![codecov.io](http://codecov.io/github/JuliaParallel/PETSc.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaParallel/PETSc.jl?branch=master)
[![Coverage Status](https://coveralls.io/repos/JuliaParallel/PETSc.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/JuliaParallel/PETSc.jl?branch=master)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://JuliaParallel.github.io/PETSc.jl/latest)

This package provides a high level interface for PETSc, enabling the use of PETSc as an `AbstractArray`.  
A low level interface is also available in the submodule `PETSc.C`.
The package supports 64-bit integers the `PetscInt` type described in 
the PETSc documentation, and `Float64`, `Float32`, and `Complex128` for the 
`PetscScalar` type.  In a default build of the package, all types can be used
simultaneously, using multiple dispatch to determine which version of PETSc
to use.

This package requires the [MPI.jl package](https://github.com/JuliaParallel/MPI.jl) be installed.  Once it is installed you should be able to run both Julia and Petsc in parallel using MPI for all communication.  The testing verifies that PETSc can be used both serially and in parallel.

To use the package, simply put `using PETSc` at the top of your Julia source file.  The module exports the names of all the functions, as well as the PETSc data type aliases and constants such as `PETSC_DECIDE`.

In general, it is possible to run PETSc in parallel. To do so with 4 processors, do:

```
mpirun -np 4 julia ./name_of_file
```

Note that this launches 4 independent Julia processes.  They are not aware of each other using Julia's built-in parallelism, and MPI is used for all communications.  

To run in serial, do:
```
julia ./name_of_file
```

Even when running serially, the [MPI.jl package](https://github.com/JuliaParallel/MPI.jl) must be installed.


An example of using a Krylov subspace method to solve a linear system is in  `test/test_ksp.jl`, which solves a simple system with a Krylov subspace method and compares the result with a direct solve using Julia's backslash operator.  This works in serial and in parallel.  It requires some variables declared at the top of `runtests.jl` to work.



## To do:
  * Make the script for building PETSc more flexible, e.g. allowing more configuration options like building BLAS or LAPCK, while ensure it remains completely autonomous (needed for Travis testing)
  * Wrap more KSP functions

## Status
### Vector
  The `AbstractArray` for `PetscVec` is implemented.  Some additional PETSc 
  BLAS functions are wrapped as well.
### Matrix
 The AbstractArray interface for `PetscMat` is implemented.  Preallocation 
 is supported through optional keyword arguments to the matrix constructor or
 the `setpreallocation` function.  It possible to set multiple values in the 
  matrix without intermediate assembly using the `assemble` function or by 
 setting the `Mat` object field `assembling` to `false` and calling `setindex`
 repeatedly.

### KSP
 Just enough KSP functions are implimented to do a GMRES solve.  Adding more 
functionality is the current priority.

## Directory Structure
  `/src` : source files.  PETSc.jl is the main file containing initialization, with the functions for each type of Petsc object in its own file.  All constants are declared in `petsc_constants.jl`.

  `/src/generated`: auto generated wrappers from Clang.jl.  Not directly useful, but easy to modify to make useful

  `/test` : contains `runtest.jl`, which does some setup and runs all tests on all three version of Petsc currently supported.  Tests for each type of Petsc object (mirroring the files in `/src`) are contained in separate files.

  `/deps` : builds Petsc if needed.  See description below


## Building PETSc
By default, building the package will build 3 versions of PETSc in the `/deps` 
 directory, and writes the file `lib_locations.jl` to the `/deps` 
 directory to tell the package the location of the libraries.  Note that 
this builds the debug versions of PETSc, which are recommended to use for all 
development.  If you wish to do high performance computations, you should 
build the optimized versions of the library.  See the PETSc website for 
details.

If you wish to build fewer than 3 version of PETSc or to use your own build 
of PETSc rather than having the package build it for you, there a several 
environmental variables that control what the build system will do.
For all the variables listed below, `name` is one of `RealDouble`, `RealSingle`,
or `ComplexDouble`, and specifies which version of the library the variable
describes.

### What to build
If the varibles `JULIA_PETSC_name_DIR` and `JULIA_PETSC_name_ARCH` are set to 
the `PETSC_DIR` and `PETSC_ARCH` of an existing PETSc installation, the build 
system will use that PETSc installation for the version of PETSc specified by
`name`.

If the variable `JULIA_PETSC_name_NOBUILD` exists (the value does not matter),
then the package will not build a version the `name`d version of PETSc.

### How to build it
If the variable `JULIA_PETSC_OPT` exists (the value does not matter), then 
a set of default optimization flags are passed to the PETSc `configure` 
script.

If the variable `JULIA_PETSC_FLAGS` exists and `JULIA_PETSC_OPT` does not, 
its value is used passed to the 
PETSc configure script (for all builds).  The user should *never* specify `--with-64-bit-indices`, `--with-scalar-type` or `--with-precision`, because this 
would break the build process for the different version of PETSc.

If neither of the above variables exist, a standard build is performed.


## Auto Generation Notes
PETSc uses preprocessor variables to decide what code to include when compiling 
the library.  Clang does not know what preprocessor variables were defined at 
compile time, so it does not correctly detect the typealiases `PetscScalar`, `PetscReal`, etc.  To correctly autogenerate wrappers, the proper variables must be passed to Clang with the -D switch.  Note that users will not need to generate their own wrappers because they have already been generated and commit to the repo.
