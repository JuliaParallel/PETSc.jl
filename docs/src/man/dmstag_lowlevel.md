# DMStag (Staggered Grids)

DMStag manages staggered grid discretizations, commonly used in computational fluid dynamics and geophysics. It handles variables defined at different locations (vertices, edges, faces, element centers).

## Overview

DMStag provides:
- Staggered grid support for MAC/marker-and-cell schemes
- Multiple degrees of freedom at vertices, edges, faces, and element centers
- Efficient communication for staggered variables
- 1D, 2D, and 3D support
- Boundary type specification

## Basic Usage Pattern

```julia
using PETSc, MPI

petsclib = PETSc.getlib()

# Create a 2D staggered grid
# Velocity components on edges, pressure at cell centers
dm = Ref{LibPETSc.CDM}()
LibPETSc.DMStagCreate2d(
    petsclib,
    MPI.COMM_WORLD,
    LibPETSc.DM_BOUNDARY_NONE,
    LibPETSc.DM_BOUNDARY_NONE,
    10, 10,                      # global dimensions
    LibPETSc.PETSC_DECIDE,
    LibPETSc.PETSC_DECIDE,
    0,                           # dof per vertex
    1,                           # dof per edge (velocity)
    1,                           # dof per element (pressure)
    LibPETSc.DMSTAG_STENCIL_BOX,
    1,                           # stencil width
    C_NULL, C_NULL,
    dm
)

# Set up
LibPETSc.DMSetFromOptions(petsclib, dm[])
LibPETSc.DMSetUp(petsclib, dm[])

# Create vectors
x = Ref{LibPETSc.CVec}()
LibPETSc.DMCreateGlobalVector(petsclib, dm[], x)

# Access staggered components
# Use DMStagGetLocationSlot to get indices for each component

# Cleanup
LibPETSc.VecDestroy(petsclib, x)
LibPETSc.DMDestroy(petsclib, dm)
```

## DMStag Functions

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["DM_wrappers.jl"]
Order   = [:function]
Filter = t -> startswith(string(nameof(t)), "DMStag")
```
