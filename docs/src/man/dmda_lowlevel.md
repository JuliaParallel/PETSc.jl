# DMDA (Structured Grids)

DMDA manages structured grids with regular topology in 1D, 2D, and 3D. It's designed for finite difference and finite volume methods on Cartesian grids.

## Overview

DMDA provides:
- Regular grid topology in 1D, 2D, or 3D
- Efficient stencil-based communication
- Natural ordering and ghost point management
- Support for multiple degrees of freedom per grid point
- Boundary type specification (periodic, ghosted, none)

## Basic Usage Pattern

```julia
using PETSc, MPI

petsclib = PETSc.getlib()

# Create a 2D structured grid: 10x10 global size, 1 dof per point
dm = Ref{LibPETSc.CDM}()
LibPETSc.DMDACreate2d(
    petsclib,
    MPI.COMM_WORLD,
    LibPETSc.DM_BOUNDARY_NONE,  # x boundary
    LibPETSc.DM_BOUNDARY_NONE,  # y boundary
    LibPETSc.DMDA_STENCIL_STAR, # stencil type
    10, 10,                      # global dimensions
    LibPETSc.PETSC_DECIDE,       # processes in x
    LibPETSc.PETSC_DECIDE,       # processes in y
    1,                           # dof per node
    1,                           # stencil width
    C_NULL, C_NULL,              # nodes per process
    dm
)

# Set up and use
LibPETSc.DMSetFromOptions(petsclib, dm[])
LibPETSc.DMSetUp(petsclib, dm[])

# Get local vector with ghost points
local_vec = Ref{LibPETSc.CVec}()
LibPETSc.DMCreateLocalVector(petsclib, dm[], local_vec)

# Get global vector
global_vec = Ref{LibPETSc.CVec}()
LibPETSc.DMCreateGlobalVector(petsclib, dm[], global_vec)

# Cleanup
LibPETSc.VecDestroy(petsclib, local_vec)
LibPETSc.VecDestroy(petsclib, global_vec)
LibPETSc.DMDestroy(petsclib, dm)
```

## DMDA Functions

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["DM_wrappers.jl"]
Order   = [:function]
Filter = t -> startswith(string(nameof(t)), "DMDA")
```
