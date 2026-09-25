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
    (PETSc.DM_BOUNDARY_NONE, PETSc.DM_BOUNDARY_NONE),  # boundary types
    (nx, ny),                                          # global dimensions
    (dof_vertex, dof_edge, dof_element),               # DOF at each location
    1,                                                 # stencil width
    PETSc.DMSTAG_STENCIL_BOX,                          # stencil type
)

dm isa PETSc.DMStag{typeof(petsclib), 2}   # true: the flavour and the dimension

# 3D staggered grid
dm = DMStag(
    petsclib,
    MPI.COMM_WORLD,
    (PETSc.DM_BOUNDARY_NONE, PETSc.DM_BOUNDARY_NONE, PETSc.DM_BOUNDARY_NONE),
    (nx, ny, nz),
    (dof_vertex, dof_edge, dof_face, dof_element),
    1,
    PETSc.DMSTAG_STENCIL_BOX,
)
```

## Accessing Data

### Grid Corners and Sizes

```julia
# Get local grid extent (without ghost points)
c = PETSc.corners(dm)
# (; lower, upper, size, nextra). `lower` and `upper` are CartesianIndex{N},
# `size` and `nextra` are NTuple{N,Int}: the returns are dimension-correct, so a
# 2D DM answers with 2-tuples and `c.size[3]` is a BoundsError (naming.md §12)

# Get local grid extent (with ghost points): (; lower, upper, size), no nextra
gc = PETSc.ghost_corners(dm)
```

### Working with Vectors

```julia
# Create global and local vectors
gvec = PETSc.global_vec(dm)
lvec = PETSc.local_vec(dm)

# Transfer data between global and local
# The written vector comes first, the DM follows it (naming.md §8)
PETSc.global_to_local!(lvec, dm, gvec, PETSc.INSERT_VALUES)
PETSc.local_to_global!(gvec, dm, lvec, PETSc.ADD_VALUES)

# Refresh the ghost points of a local vector from its neighbours, in place
PETSc.local_to_local!(lvec, dm)
```

### Getting Location Indices

```julia
# Get indices (ghost-aware) for accessing specific DOF locations in a local array
indices = PETSc.local_indices(dm)
# `indices.center` and `indices.vertex` are NamedTuples keyed by axis, and a 2D
# DM yields (x = …, y = …) with no `z` (naming.md §12)

# Get indices (no ghosts) for accessing specific DOF locations in a global array
indices = PETSc.global_indices(dm)
```

### Locations and stencils

DMStag names a point by where it sits on an element: `DMSTAG_LEFT`, `DMSTAG_DOWN`, `DMSTAG_BACK_DOWN_LEFT`, and so on. Those names shift meaning between dimensions (`DOWN` is the second axis, whatever it's called in your model), so the location functions count axes instead:

```julia
PETSc.vertex_location(dm)        # the element's lower corner
PETSc.face_location(dm, axis)    # the face normal to `axis`, on the lower side
PETSc.edge_location(dm, a, b)    # the edge touching the lower faces of `a` and `b` (3D; the vertex in 2D)
PETSc.element_location(dm)       # the element interior
```

`face_location(dm, ndims(dm))` is the last axis in any dimension. An axis outside `1:ndims(dm)` throws an `ArgumentError`.

A stencil addresses one unknown: a location, the element index and a component. [`stencil`](@ref) takes the 1-based element index that `corners` and `ghost_corners` use, and PETSc's 0-based component, as [`dof_slot`](@ref) does. Indices are not bounds-checked, so ghost elements work, and the call allocates nothing:

```julia
I = corners(dm).lower
row  = PETSc.stencil(dm, PETSc.face_location(dm, 1), I)
cols = [PETSc.stencil(dm, PETSc.element_location(dm), I),
        PETSc.stencil(dm, PETSc.element_location(dm), I + CartesianIndex(1, 0))]
```

### Assembling with stencils

```julia
J = PETSc.PetscMat(dm)
PETSc.set_values!(J, dm, [row], cols, [-1.0, 1.0], PETSc.ADD_VALUES)   # row-major block
PETSc.assemble!(J)

PETSc.zero_rows_local!(J, dm, wall_rows, 1.0)   # Dirichlet rows, given as stencils
PETSc.set_values!(b, dm, wall_rows, wall_values)
```

Under `ADD_VALUES`, repeated (row, column) entries in one call are summed. `zero_rows_local!` is collective: a rank that owns no wall passes an empty vector.

An index set of whole fields, for a field split, takes (location, component) pairs:

```julia
flow = LibPETSc.IS(dm, PETSc.face_location(dm, 1) => 0,
                       PETSc.face_location(dm, 2) => 0,
                       PETSc.element_location(dm) => 0)
PETSc.set_fieldsplit_is!(PETSc.pc(ksp), "flow", flow)
```

## Setting Coordinates

```julia
# Set uniform coordinates
PETSc.set_uniform_coordinates!(dm, xmin, xmax)                           # 1D
PETSc.set_uniform_coordinates!(dm, xmin, xmax, ymin, ymax)               # 2D
PETSc.set_uniform_coordinates!(dm, xmin, xmax, ymin, ymax, zmin, zmax)   # 3D

# Get local coordinate array
coords = PETSc.local_coordinate_array(dm)

# Read the per-axis coordinates: x[i, 1] is the lower face of element i, x[i, 2] its centre
PETSc.with_product_coordinates(dm) do x, y
    x[i, 2], y[j, 2]
end
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
    (PETSc.DM_BOUNDARY_NONE, PETSc.DM_BOUNDARY_NONE),
    (64, 64),      # 64x64 grid
    (0, 1, 1),     # 0 DOF at vertices, 1 at edges (velocity), 1 in elements (pressure)
    1,
    PETSc.DMSTAG_STENCIL_BOX,
)

PETSc.set_uniform_coordinates!(dm, 0.0, 1.0, 0.0, 1.0)

# Create vectors and matrix
x = PETSc.global_vec(dm)
b = PETSc.global_vec(dm)
A = PETSc.PetscMat(dm)
```

`PetscMat(dm)` stores every coupling the stencil allows, as explicit zeros. To let the first assembly define the pattern instead, call `PETSc.set_matrix_preallocate_only!(dm, true)` before creating the matrix.

## Functions

```@autodocs
Modules = [PETSc]
Pages   = ["dmstag.jl"]
```
