# DMDA

The `DMDA` (Distributed Array) module provides functionality for creating and managing structured grids in 1D, 2D, or 3D.

## Overview

`DMDA` is ideal for problems on regular structured grids where:
- The grid is logically rectangular
- Each grid point has the same number of degrees of freedom
- Stencil operations follow a regular pattern (star or box stencils)

## Creating a DMDA

```julia
# 2D grid example
da = DMDA(
    petsclib,
    MPI.COMM_WORLD,
    (DM_BOUNDARY_NONE, DM_BOUNDARY_NONE),  # boundary types
    (nx, ny),                               # global dimensions
    1,                                      # degrees of freedom per node
    1,                                      # stencil width
    DMDA_STENCIL_STAR                      # stencil type
)
```

## Functions

```@autodocs
Modules = [PETSc]
Pages   = ["dmda.jl"]
```
