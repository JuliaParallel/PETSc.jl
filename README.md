# PETSc

[![Build Status](https://github.com/JuliaParallel/PETSc.jl/workflows/CI/badge.svg)](https://github.com/JuliaParallel/PETSc.jl/actions/workflows/ci.yml)
[![doc stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaparallel.github.io/PETSc.jl/stable/)
[![doc dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaparallel.github.io/PETSc.jl/dev/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18274810.svg)](https://zenodo.org/badge/latestdoi/18274810/JuliaParallel/PETSc.jl)


`PETSc.jl` provides an interface to the Portable, Extensible Toolkit for Scientific Computation ([PETSc](https://petsc.org)) library, allowing the combination of Julia features (such as automatic differentiation) with the PETSc's infrastructure, including linear, nonlinear, and optimization solvers, timesteppers, domain management (DM), and more, in a distributed-memory (MPI) environment. 

This package comprises two main components:

1. An automatically generated, low-level interface for large parts of the PETSc API (see `PETSc.LibPETSc`).
2. A curated, high-level, more Julianic interface for selected functionality.

The low-level interface covers nearly the entire PETSc API, but may be awkward to work with and likely requires previous experience with PETSc to use effectively. The high level interface is designed to be more familiar and convenient for Julia users, and allows, for example, to set matrix entries with `A[1,2] = 3.0`, rather than having to call `LibPETSc.MatSetValue`. It, however, exposes only a small portion of the functionality of the underlying library. 

## Installation
This package can be added with the julia command:
```julia
julia>]add PETSc
```
The installation can be tested with
```julia
julia>]test PETSc
```

## PETSc binaries

By default, the package uses a pre-built binary of PETSc (see [PETSc_jll](https://github.com/JuliaBinaryWrappers/PETSc_jll.jl)) along with a default installation of `MPI.jl`, so you don't have to install it on your machine.

If you want to use the package with custom builds of the PETSc library, this can be done by using the function `set_petsclib` which requires you to point to the correct dynamic library (which should be compatible with the MPI version used by `MPI.jl`)

```julia
using PETSc

# Create custom library instance
petsclib = set_petsclib("/path/to/custom/libpetsc.so"; 
                       PetscScalar=Float64, PetscInt=Int64)
# Use it like any precompiled library
PETSc.initialize(petsclib, log_view=true)
# ... your code ...
PETSc.finalize(petsclib)
```
To get an overview of available precompiled libraries:
```julia
julia>using PETSc
julia>[PETSc.petsclibs...]
```

## Windows users 
The package currently does not work on windows, mainly because `MicrosoftMPI_jll` does not function when used along with the precompiled version used in `PETSc_jll`. Windows users are therefore advised to install the [Windows Subsystem for Linux](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux) (WSL) and run PETSc through there. 

## Getting started
The documentation is currently rather minimalistic; yet, if you want to see what is possible have a look at the [examples](./examples/) directory or at the tests in the [test](./test) directory. We do keep the tests up to date, so that is a good starting point. 

Note, that we do not have tests in place for the whole library at this stage. The best supported parts are `DMDA`,`DMStag`, `KSP`,`SNES`,`Vec` and `Mat` interfaces, while other parts such as `DMPlex` do not have a high-level interface or tests yet. Users will thus have to rely on the low-level interface.