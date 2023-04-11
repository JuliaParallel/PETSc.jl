# PETSc

[![Build Status](https://github.com/JuliaParallel/PETSc.jl/workflows/CI/badge.svg)](https://github.com/JuliaParallel/PETSc.jl/actions/workflows/ci.yml)
[![doc dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaparallel.github.io/PETSc.jl/dev/)

PETSc.jl provides an interface to the [PETSc](https://www.mcs.anl.gov/petsc/) library,
allowing the combination of Julia features (such as automatic differentiation) with the PETSc's infrastructure, including
linear, nonlinear, and optimization solvers, timesteppers, domain management (DM), and more, in a distributed-memory (MPI) environment.

This package comprises two main components:

1. An automatically generated, low-level interface for large parts of the PETSc API.
2. A curated, high-level, more Julianic interface for selected functionality.

The low-level interface covers the entire PETSc API, but may be awkward to work with and likely requires
previous experience with PETSc to use effectively. The high level interface is designed to be
more familiar and convenient for Julia users but exposes only a small portion of the functionality
of the underlying library. This high-level interface should be considered unstable, as its
implementation involves design decisions which may be revisited - note the low version number
of this package if relying on these high-level features.


## Installation

This package can be added with the julia command:
```julia
]add PETSc
```
The installation can be tested with
```julia
]test PETSc
```

## BinaryBuilder Version

By default, the package uses a pre-built binary of
[`PETSc`](https://github.com/JuliaBinaryWrappers/PETSc_jll.jl) along with a
default installation of `MPI.jl`. Note that the distributed version of PETSc is using real,
`Float64` numbers; build details can be found
[here](https://github.com/JuliaPackaging/Yggdrasil/blob/master/P/PETSc/build_tarballs.jl)

## System Builds

If you want to use the package with custom builds of the PETSc library, this can
be done by specifying the environment variable `JULIA_PETSC_LIBRARY`. This is a
colon separated list of paths to custom builds of PETSc; the reason for using
multiple builds is to enable single, double, and complex numbers in the same
julia session. These should be built against the same version of MPI as used
with `MPI.jl`

After setting the variable you should
```julia
]build PETSc
```
and the library will be persistently set until the next time the build command
is issued.

To see the currently set library use
```julia
using PETSc
PETSc.libs
```
