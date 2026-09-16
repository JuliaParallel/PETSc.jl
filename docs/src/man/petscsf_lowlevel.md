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

# Broadcast: send root data to leaves (the datatype and op are MPI.jl objects)
LibPETSc.PetscSFBcastBegin(petsclib, sf, MPI.Datatype(Float64), root_data, leaf_data, MPI.REPLACE)
LibPETSc.PetscSFBcastEnd(petsclib, sf, MPI.Datatype(Float64), root_data, leaf_data, MPI.REPLACE)
```

### Reduce

Accumulate data from leaves back to roots:

```julia
# Leaf contributions
leaf_data = Float64[0.1, 0.2, 0.3, 0.4, 0.5]

# Root accumulator
root_data = zeros(Float64, nroots)

# Reduce: accumulate leaf data to roots
LibPETSc.PetscSFReduceBegin(petsclib, sf, MPI.Datatype(Float64), leaf_data, root_data, MPI.SUM)
LibPETSc.PetscSFReduceEnd(petsclib, sf, MPI.Datatype(Float64), leaf_data, root_data, MPI.SUM)
```

### Fetch and Operations

Atomic operations for concurrent updates:

```julia
# Fetch the old root value into leaf_updates, then apply op(root, leaf) at the root
leaf_updates = zeros(Float64, nleaves)
LibPETSc.PetscSFFetchAndOpBegin(petsclib, sf, MPI.Datatype(Float64), root_data, leaf_data, leaf_updates, MPI.SUM)
LibPETSc.PetscSFFetchAndOpEnd(petsclib, sf, MPI.Datatype(Float64), root_data, leaf_data, leaf_updates, MPI.SUM)
```

## MPI Operations

Supported MPI operations for reduce (MPI.jl objects):
- `MPI.SUM`: Sum values
- `MPI.MAX`: Maximum value
- `MPI.MIN`: Minimum value
- `MPI.REPLACE`: Replace (last write wins)
- `MPI.PROD`: Product

The communication routines (`PetscSFBcastBegin/End`, `PetscSFReduceBegin/End`,
`PetscSFFetchAndOpBegin/End`) are hand-written wrappers: PETSc's API extractor skips
functions taking an `MPI_Datatype`. They accept `Array`s or raw pointers and keep the
arrays alive for the duration of the call; the arrays must not be freed between `Begin`
and `End`.

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
# Get the graph: root count, leaf count, leaf locations and remote (rank, index) pairs.
# `ilocal` is `nothing` when the leaves are contiguous [0, nleaves); the arrays are
# owned by the SF and valid until it changes.
nroots, nleaves, ilocal, iremote = LibPETSc.PetscSFGetGraph(petsclib, sf)
for leaf in 1:nleaves
    node = iremote[leaf]          # node.rank, node.index
end
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
# After modifying owned data, update ghost points: broadcast roots to leaves
LibPETSc.PetscSFBcastBegin(petsclib, sf, MPI.Datatype(Float64), local_data, ghost_data, MPI.REPLACE)
LibPETSc.PetscSFBcastEnd(petsclib, sf, MPI.Datatype(Float64), local_data, ghost_data, MPI.REPLACE)
```

### 2. Parallel Assembly

```julia
# After local assembly, accumulate contributions from other processes at the owners
LibPETSc.PetscSFReduceBegin(petsclib, sf, MPI.Datatype(Float64), local_contrib, global_data, MPI.SUM)
LibPETSc.PetscSFReduceEnd(petsclib, sf, MPI.Datatype(Float64), local_contrib, global_data, MPI.SUM)
```

### 3. DM Point Communication

```julia
# Get the point SF of a DM (describes the point distribution); the DM owns it
dm_sf = LibPETSc.DMGetPointSF(petsclib, dm)

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

### Communication (hand-written wrappers)

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/extra_wrappers.jl"]
Order   = [:function]
Filter  = t -> startswith(string(nameof(t)), "PetscSF")
```
