# PetscSF (Star Forest) - Low-level Interface

The PetscSF (Star Forest) component provides efficient parallel communication patterns for distributed data structures. A star forest is a specialized graph structure optimized for scatter/gather operations in parallel computing.

## Overview

PetscSF enables:
- **Point-to-point communication**: Efficient MPI communication patterns
- **Scatter/gather operations**: Move data between processors
- **Halo exchange**: Update ghost/boundary values
- **Reduction operations**: Parallel sums, max, min across shared data
- **Irregular communication**: Handle non-uniform data distributions

A star forest consists of:
- **Roots**: Data owned locally
- **Leaves**: Data needed from remote processes (or local)
- **Communication pattern**: Which leaves come from which roots

PetscSF is the underlying communication layer for DM ghost point updates and other parallel operations.

## Basic Usage

```julia
using PETSc, MPI

# Initialize MPI and PETSc
MPI.Init()
petsclib = PETSc.getlib()
PETSc.initialize(petsclib)
PetscInt = petsclib.PetscInt

# Create a star forest
sf = LibPETSc.PetscSFCreate(petsclib, MPI.COMM_WORLD)

# Define communication pattern
# nleaves: number of leaves (data items we need)
# ilocal: local indices for leaves (can be C_NULL if identity)
# iremote: (rank, index) pairs specifying which process/index to get from

nleaves = 5
# number of roots owned locally (for this simple example set equal to nleaves)
nroots = 5
ilocal = [0, 1, 2, 3, 4]  # Local indices where data will be stored
iremote = [
    LibPETSc.PetscSFNode(0, 0),
    LibPETSc.PetscSFNode(0, 1),
    LibPETSc.PetscSFNode(0, 2),
    LibPETSc.PetscSFNode(0, 3),
    LibPETSc.PetscSFNode(0, 4),
]

LibPETSc.PetscSFSetGraph(petsclib, sf, nroots, nleaves, ilocal, LibPETSc.PETSC_COPY_VALUES,
                         iremote, LibPETSc.PETSC_COPY_VALUES)

# Setup
LibPETSc.PetscSFSetUp(petsclib, sf)

# Cleanup
LibPETSc.PetscSFDestroy(petsclib, sf)

# Finalize PETSc and MPI
PETSc.finalize(petsclib)
MPI.Finalize()
```

## Communication Operations

### Broadcast (Scatter)

Send data from roots to leaves:

```julia
# Root data: data we own
root_data = Float64[1.0, 2.0, 3.0, 4.0, 5.0]

# Leaf data: buffer to receive data
leaf_data = zeros(Float64, nleaves)

# Broadcast: send root data to leaves
LibPETSc.PetscSFBcastBegin(petsclib, sf, LibPETSc.MPI_DOUBLE, root_data, leaf_data,
                           LibPETSc.MPI_REPLACE)
LibPETSc.PetscSFBcastEnd(petsclib, sf, LibPETSc.MPI_DOUBLE, root_data, leaf_data,
                         LibPETSc.MPI_REPLACE)
```

### Reduce

Accumulate data from leaves back to roots:

```julia
# Leaf contributions
leaf_data = Float64[0.1, 0.2, 0.3, 0.4, 0.5]

# Root accumulator
root_data = zeros(Float64, nroots)

# Reduce: accumulate leaf data to roots
LibPETSc.PetscSFReduceBegin(petsclib, sf, LibPETSc.MPI_DOUBLE, leaf_data, root_data,
                            LibPETSc.MPI_SUM)
LibPETSc.PetscSFReduceEnd(petsclib, sf, LibPETSc.MPI_DOUBLE, leaf_data, root_data,
                          LibPETSc.MPI_SUM)
```

### Fetch and Operations

Atomic operations for concurrent updates:

```julia
# Fetch data and apply operation
LibPETSc.PetscSFFetchAndOpBegin(petsclib, sf, LibPETSc.MPI_DOUBLE, root_data,
                                leaf_data, leaf_updates, LibPETSc.MPI_SUM)
LibPETSc.PetscSFFetchAndOpEnd(petsclib, sf, LibPETSc.MPI_DOUBLE, root_data,
                              leaf_data, leaf_updates, LibPETSc.MPI_SUM)
```

## MPI Operations

Supported MPI operations for reduce:
- `MPI_SUM`: Sum values
- `MPI_MAX`: Maximum value
- `MPI_MIN`: Minimum value
- `MPI_REPLACE`: Replace (last write wins)
- `MPI_PROD`: Product

## Star Forest Types

Available through `PetscSFSetType`:
- **PETSCSFBASIC**: Basic implementation
- **PETSCSFNEIGHBOR**: MPI neighborhood collectives (efficient for structured patterns)
- **PETSCSFALLGATHERV**: All-gather based
- **PETSCSFALLGATHER**: All-gather for small data
- **PETSCSFGATHERV**: Gather-based
- **PETSCSFGATHER**: Simple gather
- **PETSCSFALLTOALL**: All-to-all based

## Graph Queries

```julia
# Get number of roots (locally owned data)
nroots = Ref{PetscInt}()
LibPETSc.PetscSFGetGraph(petsclib, sf, nroots, C_NULL, C_NULL, C_NULL)

# Get number of leaves
nleaves = Ref{PetscInt}()
LibPETSc.PetscSFGetGraph(petsclib, sf, C_NULL, nleaves, C_NULL, C_NULL)

# Get full graph
ilocal_ptr = Ref{Ptr{PetscInt}}()
iremote_ptr = Ref{Ptr{LibPETSc.PetscSFNode}}()
LibPETSc.PetscSFGetGraph(petsclib, sf, nroots, nleaves, ilocal_ptr, iremote_ptr)
```

## Multi-Root Support

Handle communication with multiple root data per point:

```julia
# Create multi-SF for multiple DOFs per point
nroots_mult = nroots * num_components
multi_sf = LibPETSc.PetscSFCreateEmbeddedRootSF(petsclib, sf, nroots_mult, iroot_indices)
```

## Common Use Cases

### 1. Ghost Point Updates (Halo Exchange)

```julia
# After modifying owned data, update ghost points
# 1. Pack local data
# 2. Broadcast to leaves (ghost points)
LibPETSc.PetscSFBcastBegin(petsclib, sf, datatype, local_data, ghost_data, op)
LibPETSc.PetscSFBcastEnd(petsclib, sf, datatype, local_data, ghost_data, op)
```

### 2. Parallel Assembly

```julia
# After local assembly, accumulate contributions from other processes
# 1. Each process computes local contributions
# 2. Reduce to accumulate at owners
LibPETSc.PetscSFReduceBegin(petsclib, sf, datatype, local_contrib, global_data, MPI_SUM)
LibPETSc.PetscSFReduceEnd(petsclib, sf, datatype, local_contrib, global_data, MPI_SUM)
```

### 3. DM Point Communication

```julia
# Get natural SF for a DM (describes point distribution)
dm_sf = Ref{LibPETSc.PetscSF}()
# LibPETSc.DMGetPointSF(petsclib, dm, dm_sf)

# Use to communicate point-based data
```

## Performance Considerations

- **Choose appropriate type**: `PETSCSFNEIGHBOR` is often best for structured grids
- **Reuse SF objects**: Creating the communication pattern is expensive
- **Batch communications**: Combine multiple small messages when possible
- **Alignment**: Use properly aligned data types for better performance

## Function Reference

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscSF_wrappers.jl"]
Order   = [:function]
```
