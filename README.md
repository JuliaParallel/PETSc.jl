# PETSc

[![Build Status](https://github.com/JuliaParallel/PETSc.jl/workflows/CI/badge.svg)](https://github.com/JuliaParallel/PETSc.jl/actions/workflows/ci.yml)

This package provides a low level interface for [PETSc](https://www.mcs.anl.gov/petsc/) and allows combining julia features (such as automatic differentiation) with the PETSc infrastructure and nonlinear solvers.

## Installation

This package can be added with the julia command:
```julia
]add https://github.com/JuliaParallel/PETSc.jl
```
The installation can be tested with
```julia
]test PETSc
```

## Caveats

The package uses a pre-build binary of [`PETSc`](https://github.com/JuliaBinaryWrappers/PETSc_jll.jl) that is automatically installed, along with a default installation of `MPI.jl`; use of system install MPI and PETSc is not currently supported. Also note that `PETSc_jll.jl` is currently not working on Windows for julia > 1.3. 
