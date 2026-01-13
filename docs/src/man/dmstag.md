# DMStag

The DMStag (Staggered Grid DM) module provides data management for staggered grids, commonly used in finite difference/volume methods for fluid dynamics and similar applications.

## Overview

DMStag is designed for problems where:
- Variables live at different grid locations (vertices, edges, faces, cell centers)
- Staggered grids provide better stability for incompressible flow
- Multiple degrees of freedom per grid location are needed

### Staggered Grid Layout

In a staggered grid, different physical quantities are stored at different locations:

**1D**: 
- Vertices: Scalar quantities (pressure, temperature)
- Elements: Flux quantities

**2D**:
- Vertices: Corner values
- Edges: Face-normal velocities (u on vertical edges, v on horizontal edges)  
- Elements: Cell-centered values (pressure)

**3D**:
- Vertices: Corner values
- Edges: Edge-centered values
- Faces: Face-normal quantities
- Elements: Cell-centered values

## Creating a DMStag

```julia
# 2D staggered grid
dm = DMStag(
    petsclib,
    MPI.COMM_WORLD,
    (DM_BOUNDARY_NONE, DM_BOUNDARY_NONE),  # boundary types
    (nx, ny),                               # global dimensions
    (dof_vertex, dof_edge, dof_element),   # DOF at each location
    1,                                      # stencil width
    DMSTAG_STENCIL_BOX                     # stencil type
)

# 3D staggered grid
dm = DMStag(
    petsclib,
    MPI.COMM_WORLD,
    (DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE),
    (nx, ny, nz),
    (dof_vertex, dof_edge, dof_face, dof_element),
    1,
    DMSTAG_STENCIL_BOX
)
```

## Accessing Data

### Grid Corners and Sizes

```julia
# Get local grid extent (without ghost points)
corners = getcorners_dmstag(dm)
# Returns (lower=CartesianIndex, upper=CartesianIndex, size=Tuple)

# Get local grid extent (with ghost points)  
ghost_corners = getghostcorners_dmstag(dm)
```

### Working with Vectors

```julia
# Create global and local vectors
global_vec = DMGlobalVec(dm)
local_vec = DMLocalVec(dm)

# Transfer data between global and local
dm_global_to_local!(dm, global_vec, INSERT_VALUES, local_vec)
dm_local_to_global!(dm, local_vec, ADD_VALUES, global_vec)
```

### Getting Location Indices

```julia
# Get indices for accessing specific DOF locations
indices = DMStagGetIndices(dm)
# Use indices to access vertex, edge, face, or element DOFs
```

## Setting Coordinates

```julia
# Set uniform coordinates
setuniformcoordinates_stag!(dm, xmin, xmax)           # 1D
setuniformcoordinates_stag!(dm, xmin, xmax, ymin, ymax)  # 2D
setuniformcoordinates_stag!(dm, xmin, xmax, ymin, ymax, zmin, zmax)  # 3D

# Get local coordinate array
coords = getlocalcoordinatearray(dm)
```

## Stencil Types

- `DMSTAG_STENCIL_BOX` - Full box stencil (includes diagonals)
- `DMSTAG_STENCIL_STAR` - Star stencil (axis-aligned neighbors only)

## Example: 2D Stokes Flow Setup

```julia
# Create staggered grid for Stokes: velocity on edges, pressure in cells
dm = DMStag(
    petsclib,
    MPI.COMM_WORLD,
    (DM_BOUNDARY_NONE, DM_BOUNDARY_NONE),
    (64, 64),      # 64x64 grid
    (0, 1, 1),     # 0 DOF at vertices, 1 at edges (velocity), 1 in elements (pressure)
    1,
    DMSTAG_STENCIL_BOX
)

setuniformcoordinates_stag!(dm, 0.0, 1.0, 0.0, 1.0)

# Create vectors and matrix
x = DMGlobalVec(dm)
b = DMGlobalVec(dm)
A = DMCreateMatrix(dm)
```

## Functions

```@autodocs
Modules = [PETSc]
Pages   = ["dmstag.jl"]
```
