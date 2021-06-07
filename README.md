# PETSc

[![Build Status](https://github.com/jkozdon/PETSc.jl/workflows/CI/badge.svg)](https://github.com/jkozdon/PETSc.jl/actions/workflows/ci.yml)
[![codecov.io](http://codecov.io/github/JuliaParallel/PETSc.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaParallel/PETSc.jl?branch=master)

This package provides a low level interface for
[PETSc](https://www.mcs.anl.gov/petsc/)


## Installation

This package can be added with the julia command:
```julia
]add https://github.com/JuliaParallel/PETSc.jl
```

## Caveats

The package requires the uses a pre-build binary of
[`PETSc`](https://github.com/JuliaBinaryWrappers/PETSc_jll.jl) along with a
default installation of `MPI.jl`; use of system install MPI and PETSc is not
currently supported.
