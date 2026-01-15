# Mat - Low-level Interface

This page documents the low-level, automatically wrapped PETSc Mat (matrix) functions available through `LibPETSc`. These functions provide direct access to the PETSc C API.

For the high-level Julia interface, see [Mat](@ref).

## Overview

The Mat interface includes:
- **Core Mat operations**: Creation, manipulation, and mathematical operations on matrices (~632 functions in `Mat_wrappers.jl`)
- **Mat utilities**: Additional matrix operations and conversions (~124 functions in `Mataddons_wrappers.jl`)

## Usage

All low-level Mat functions require a `petsclib` parameter as the first argument:

```julia
using PETSc

# Get the library instance
petsclib = PETSc.petsclibs[1]
PETSc.initialize(petsclib)

# Get PETSc types
PetscInt = petsclib.PetscInt
PetscScalar = petsclib.PetscScalar

# Create a sparse matrix
mat = LibPETSc.MatCreate(petsclib, LibPETSc.PETSC_COMM_SELF)
LibPETSc.MatSetSizes(petsclib, mat, PetscInt(10), PetscInt(10), PetscInt(10), PetscInt(10))
LibPETSc.MatSetType(petsclib, mat, "seqaij")
LibPETSc.MatSetUp(petsclib, mat)

# Set values (0-based indexing!)
row = PetscInt[0]
cols = PetscInt[0, 1]
vals = PetscScalar[2.0, -1.0]
LibPETSc.MatSetValues(petsclib, mat, PetscInt(1), row, PetscInt(2), cols, vals, LibPETSc.INSERT_VALUES)

# Assemble
LibPETSc.MatAssemblyBegin(petsclib, mat, LibPETSc.MAT_FINAL_ASSEMBLY)
LibPETSc.MatAssemblyEnd(petsclib, mat, LibPETSc.MAT_FINAL_ASSEMBLY)

# Clean up
LibPETSc.MatDestroy(petsclib, mat)
PETSc.finalize(petsclib)
```

## Function Reference

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/Mat_wrappers.jl", "autowrapped/Mataddons_wrappers.jl"]
Order   = [:function]
```
