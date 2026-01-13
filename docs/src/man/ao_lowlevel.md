# AO (Application Ordering) - Low-level Interface

The AO (Application Ordering) component provides mappings between the natural "application" ordering of unknowns and the PETSc parallel ordering, which is optimized for parallel sparse matrix operations.

## Overview

Application orderings are used when:
- **Natural ordering differs from parallel ordering**: Your application uses a different numbering scheme
- **I/O operations**: Reading/writing data in application order
- **User interaction**: Displaying results in familiar ordering
- **Legacy code integration**: Interfacing with existing applications

The AO object maintains bidirectional mappings:
- **Application → PETSc**: Convert from your numbering to PETSc's
- **PETSc → Application**: Convert from PETSc's numbering to yours

## Basic Usage

```julia
using PETSc

# Initialize PETSc
petsclib = PETSc.getlib()

# Define the mapping
# application[i] is the application index for PETSc index i
n = 10
application = Int32[9, 8, 7, 6, 5, 4, 3, 2, 1, 0]  # Reverse ordering
petsc = Int32[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]        # Natural PETSc ordering

# Create AO object
ao = Ref{LibPETSc.AO}()
LibPETSc.AOCreateBasic(petsclib, MPI.COMM_SELF, n, application, petsc, ao)

# Convert application indices to PETSc indices
app_indices = Int32[0, 5, 9]
LibPETSc.AOApplicationToPetsc(petsclib, ao[], length(app_indices), app_indices)
# app_indices now contains corresponding PETSc indices

# Convert PETSc indices to application indices  
petsc_indices = Int32[0, 1, 2]
LibPETSc.AOPetscToApplication(petsclib, ao[], length(petsc_indices), petsc_indices)
# petsc_indices now contains corresponding application indices

# Cleanup
LibPETSc.AODestroy(petsclib, ao)
```

## Creating AO Objects

### Basic AO

Most common: explicit mapping arrays:

```julia
# Create from application and PETSc index arrays
LibPETSc.AOCreateBasic(petsclib, comm, n, app_indices, petsc_indices, ao)
```

### Identity Mapping

When application and PETSc orderings are the same:

```julia
# Create identity mapping (no conversion needed)
LibPETSc.AOCreateIdentity(petsclib, comm, n, ao)
```

### Memory Scalable

For large problems where storing full mapping is expensive:

```julia
# Create memory-scalable AO (uses hash table)
LibPETSc.AOCreateMemoryScalable(petsclib, comm, n, app_indices, petsc_indices, ao)
```

## Conversion Operations

### Forward Conversion (Application → PETSc)

```julia
# Convert array of application indices to PETSc indices
indices = Int32[10, 20, 30, 40]
LibPETSc.AOApplicationToPetsc(petsclib, ao[], length(indices), indices)
# indices are now in PETSc ordering (modified in-place)
```

### Reverse Conversion (PETSc → Application)

```julia
# Convert array of PETSc indices to application indices
indices = Int32[0, 5, 10, 15]
LibPETSc.AOPetscToApplication(petsclib, ao[], length(indices), indices)
# indices are now in application ordering (modified in-place)
```

### Index Set Conversion

```julia
# Convert IS (index set) from application to PETSc ordering
is_app = Ref{LibPETSc.IS}()
# ... create IS with application indices ...

LibPETSc.AOApplicationToPetscIS(petsclib, ao[], is_app[])
# is_app now contains PETSc indices
```

## Parallel Considerations

- Each processor has its own local application ordering
- AO objects handle parallel communication automatically
- Indices not owned locally are mapped via MPI communication

```julia
# Create parallel AO
LibPETSc.AOCreateBasic(petsclib, MPI.COMM_WORLD, local_n, 
                       local_app_indices, local_petsc_indices, ao)
```

## Common Use Cases

### 1. Reading Data in Application Order

```julia
# Read matrix entries in application ordering
# Then convert to PETSc ordering for assembly

# Application-ordered rows/cols
app_rows = Int32[...]
app_cols = Int32[...]

# Convert to PETSc ordering
LibPETSc.AOApplicationToPetsc(petsclib, ao[], length(app_rows), app_rows)
LibPETSc.AOApplicationToPetsc(petsclib, ao[], length(app_cols), app_cols)

# Now use app_rows, app_cols (which are now in PETSc ordering) for MatSetValues
```

### 2. Displaying Results

```julia
# After solve, convert solution indices for output
solution_indices = Int32[0, 1, 2, 3, 4]  # PETSc ordering

# Convert to application ordering for display
LibPETSc.AOPetscToApplication(petsclib, ao[], length(solution_indices), solution_indices)

# Display in application order
# for i in solution_indices
#     println("Application DOF $i: value = ", solution[i])
# end
```

### 3. Integration with DM

```julia
# Get AO from DM
dm_ao = Ref{LibPETSc.AO}()
# LibPETSc.DMGetAO(petsclib, dm, dm_ao)

# Use for converting between natural and distributed orderings
```

## Viewing and Debugging

```julia
# View the mapping (for debugging)
viewer = LibPETSc.PETSC_VIEWER_STDOUT_SELF(petsclib)
LibPETSc.AOView(petsclib, ao[], viewer)
```

## AO Types

Available types (set with `AOSetType`):
- **AOBASIC**: Hash table based, good for general use
- **AOMEMSCALABLE**: Memory efficient for large problems
- **AOMAPPING1TO1**: Optimized for one-to-one mappings

## Performance Tips

- **Reuse AO objects**: Creation can be expensive for large problems
- **Batch conversions**: Convert arrays rather than individual indices
- **Use identity when possible**: Skip AO entirely if orderings match
- **Choose appropriate type**: `AOMEMSCALABLE` for very large problems

## Function Reference

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/AO_wrappers.jl"]
Order   = [:function]
```
