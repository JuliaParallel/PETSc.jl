# IS (Index Sets) - Low-level Interface

The IS (Index Set) component provides data structures and operations for managing sets of indices, which are fundamental for parallel data distribution, matrix/vector assembly, and mapping between different numbering schemes.

## Overview

Index sets are used throughout PETSc for:
- **Parallel data distribution**: Specifying which indices are owned by each processor
- **Submatrix/subvector extraction**: Selecting specific rows/columns
- **Scatter/gather operations**: Defining communication patterns
- **Field splitting**: Identifying DOFs belonging to different fields
- **Reordering**: Specifying permutations for bandwidth reduction

PETSc provides several IS types:
- **General IS**: Arbitrary list of indices
- **Stride IS**: Regularly spaced indices (start, step, length)
- **Block IS**: Block-structured indices
- **Complement IS**: All indices except specified ones

## Basic Usage

```julia
using PETSc

# Initialize PETSc
petsclib = PETSc.getlib()

# Create an index set from an array of indices (0-based)
# Note: indices should be PetscInt type (typically Int64, not Int32)
indices = PetscInt[0, 2, 4, 6, 8]
is = LibPETSc.ISCreateGeneral(petsclib, MPI.COMM_SELF, length(indices), indices, LibPETSc.PETSC_COPY_VALUES)

# Create a stride index set: indices = first:step:(first + step*(n-1))
is_stride = LibPETSc.ISCreateStride(petsclib, MPI.COMM_SELF, 10, 0, 2)  # 0, 2, 4, ..., 18

# Get the size of an index set (returns value directly)
n = LibPETSc.ISGetSize(petsclib, is)

# Get local size (returns value directly)
local_n = LibPETSc.ISGetLocalSize(petsclib, is)

# Get indices as an array
indices_ptr = Ref{Ptr{PetscInt}}()
LibPETSc.ISGetIndices(petsclib, is, indices_ptr)
# ... use indices ...
LibPETSc.ISRestoreIndices(petsclib, is, indices_ptr)

# Cleanup
LibPETSc.ISDestroy(petsclib, is)
LibPETSc.ISDestroy(petsclib, is_stride)
```

## Common Operations

### Creating Index Sets

```julia
# General index set from array
LibPETSc.ISCreateGeneral(petsclib, comm, n, indices, copymode, is)

# Stride index set: first, first+step, first+2*step, ...
LibPETSc.ISCreateStride(petsclib, comm, n, first, step, is)

# Block index set: block-structured indices
LibPETSc.ISCreateBlock(petsclib, comm, blocksize, n, indices, copymode, is)
```

### Set Operations

```julia
# Union of two index sets
LibPETSc.ISSum(petsclib, is1, is2, is_union)

# Difference: is1 - is2
LibPETSc.ISDifference(petsclib, is1, is2, is_diff)

# Intersection
LibPETSc.ISIntersect(petsclib, is1, is2, is_intersect)

# Complement: all indices in [0, n) not in is
LibPETSc.ISComplement(petsclib, is, nmin, nmax, is_complement)
```

### Querying Properties

```julia
# Check if index set is sorted
is_sorted = Ref{PetscBool}()
LibPETSc.ISSorted(petsclib, is, is_sorted)

# Check if identity permutation
is_identity = Ref{PetscBool}()
LibPETSc.ISIdentity(petsclib, is, is_identity)

# Check if a permutation
is_perm = Ref{PetscBool}()
LibPETSc.ISPermutation(petsclib, is, is_perm)
```

## Index Set Types

Available through `ISSetType`:
- **ISGENERAL**: General list of indices
- **ISSTRIDE**: Arithmetic sequence
- **ISBLOCK**: Block-structured indices

## Scatter Context

Index sets are used to create scatter contexts for moving data between vectors:

```julia
# Create scatter context
scatter = Ref{LibPETSc.VecScatter}()
LibPETSc.VecScatterCreate(petsclib, vec_from, is_from, vec_to, is_to, scatter)

# Perform scatter operation
LibPETSc.VecScatterBegin(petsclib, scatter[], vec_from, vec_to, INSERT_VALUES, SCATTER_FORWARD)
LibPETSc.VecScatterEnd(petsclib, scatter[], vec_from, vec_to, INSERT_VALUES, SCATTER_FORWARD)
```

## Parallel Considerations

- Index sets use global indexing by default
- Each processor owns a portion of the index set
- Use `ISGetLocalSize` to get the local portion
- Scatter operations handle parallel communication automatically

## Function Reference

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/IS_wrappers.jl"]
Order   = [:function]
```

## IS Add-ons

Additional IS utilities and helper functions:

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/ISaddons_wrappers.jl"]
Order   = [:function]
```
