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
c = PETSc.corners(dm)
# Returns (lower=CartesianIndex, upper=CartesianIndex, size=Tuple, nextra=Tuple)

# Get local grid extent (with ghost points)
gc = PETSc.ghost_corners(dm)
```

### Working with Vectors

```julia
# Create global and local vectors
gvec = PETSc.global_vec(dm)
lvec = PETSc.local_vec(dm)

# Transfer data between global and local
PETSc.global_to_local!(gvec, lvec, dm, PETSc.INSERT_VALUES)
PETSc.local_to_global!(lvec, gvec, dm, PETSc.ADD_VALUES)
```

### Getting Location Indices

```julia
# Get indices (ghost-aware) for accessing specific DOF locations in a local array
indices = local_indices(dm)
# Use indices to access vertex, edge, face, or element DOFs

# Get indices (no ghosts) for accessing specific DOF locations in a global array
indices = global_indices(dm)
```

## Setting Coordinates

```julia
# Set uniform coordinates
set_uniform_coordinates!(dm, xmin, xmax)           # 1D
set_uniform_coordinates!(dm, xmin, xmax, ymin, ymax)  # 2D
set_uniform_coordinates!(dm, xmin, xmax, ymin, ymax, zmin, zmax)  # 3D

# Get local coordinate array
coords = local_coordinate_array(dm)
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

set_uniform_coordinates!(dm, 0.0, 1.0, 0.0, 1.0)

# Create vectors and matrix
x = global_vec(dm)
b = global_vec(dm)
A = DMCreateMatrix(dm)
```

## Functions

```@autodocs
Modules = [PETSc]
Pages   = ["dmstag.jl"]
```
