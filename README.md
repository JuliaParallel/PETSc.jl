# PETSc

[![Build Status](https://github.com/JuliaParallel/PETSc.jl/workflows/CI/badge.svg)](https://github.com/JuliaParallel/PETSc.jl/actions/workflows/ci.yml)
[![doc dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaparallel.github.io/PETSc.jl/dev/)

This package provides a low level interface for [PETSc](https://www.mcs.anl.gov/petsc/) and allows combining julia features (such as automatic differentiation) with the PETSc infrastructure and nonlinear solvers.

## Installation

This package can be added with the julia command:
```julia
]add PETSc
```
The installation can be tested with
```julia
]test PETSc
```

## PETSc binaries

By default, the package uses a pre-built binary of (see [PETSc_jll](https://github.com/JuliaBinaryWrappers/PETSc_jll.jl)) along with a default installation of `MPI.jl`, so you don't have to install it on your machine.

If you want to use the package with custom builds of the PETSc library, this can be done by specifying the environment variable `JULIA_PETSC_LIBRARY`. This is a colon separated list of paths to custom builds of PETSc; the reason for using multiple builds is to enable single, double, and complex numbers in the same julia session. These should be built against the same version of MPI as used with `MPI.jl`

After setting the variable you should
```julia
]build PETSc
```
and the library will be persistently set until the next time the build command is issued.

To see the currently set library use
```julia
using PETSc
PETSc.libs
```
## Windows users 
The package currently does not work on windows, mainly because `MicrosoftMPI_jll` does not function when used along with the precompiled version used in `PETSc_jll`. Windows users are therefore advised to install the Windows Subshell for Linux (WSL) and run PETSc through there. 

## Getting started
The documentation is currently rather minimalistic; yet, if you want to see what is possible have a look at the [examples](./examples/) directory or at the tests in the [test](./test) directory. We do keep the tests up to date, so that is a good starting point.
