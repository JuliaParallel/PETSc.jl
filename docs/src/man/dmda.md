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
    (PETSc.DM_BOUNDARY_NONE, PETSc.DM_BOUNDARY_NONE),  # boundary types
    (nx, ny),                                          # global dimensions
    1,                                                 # degrees of freedom per node
    1,                                                 # stencil width
    PETSc.DMDA_STENCIL_STAR,                           # stencil type
)

da isa PETSc.DMDA{typeof(petsclib), 2}   # true: the flavour and the dimension
```

`DMDA` is a Julia type rather than a factory function, so the dimension is in
the type and the returns are dimension-correct: `PETSc.corners(da).size` is a
2-tuple for the grid above, and `PETSc.info(da)` reports `ndofs` and `procs`
([naming conventions](naming.md), §5.3 and §12).

## Reading the local extent

```julia
c  = PETSc.corners(da)         # (; lower, upper, size), 1-based and inclusive
gc = PETSc.ghost_corners(da)   # the same, with the ghost points
i  = PETSc.info(da)            # (; dim, global_size, procs, ndofs, …)
ndims(da)                      # the dimension, as Base spells it
```

## Colouring for a finite-difference Jacobian

```julia
# Takes no `petsclib`, and every index field says its base (naming.md §12.1)
c = PETSc.star_fd_coloring(da)
c.row_coo_local_0b     # straight to MatSetPreallocationCOOLocal
c.perturb_cols_1b      # indexes Julia arrays
```

## Functions

```@autodocs
Modules = [PETSc]
Pages   = ["dmda.jl"]
```
