# Vec - Low-level Interface

This page documents the low-level, automatically wrapped PETSc Vec functions available through `LibPETSc`. These functions provide direct access to the PETSc C API.

For the high-level Julia interface, see [Vec](@ref).

## Overview

The Vec interface includes:
- **Core Vec operations**: Creation, manipulation, and mathematical operations on vectors (~186 functions in `Vec_wrappers.jl`)
- **Vec utilities**: Scatter operations, FFTW integration, and advanced vector operations (~76 functions in `Vecs_wrappers.jl`)
- **VecTagger**: Tools for tagging/marking vector elements based on criteria (~31 functions in `VecTagger_wrappers.jl`)

## Usage

All low-level Vec functions require a `petsclib` parameter as the first argument:

```julia
using PETSc

# Get the library instance
petsclib = PETSc.petsclibs[1]
PETSc.initialize(petsclib)

# Get PETSc types
PetscInt = petsclib.PetscInt

# Create a vector
vec = LibPETSc.VecCreate(petsclib, LibPETSc.PETSC_COMM_SELF)
LibPETSc.VecSetSizes(petsclib, vec, PetscInt(10), PetscInt(10))
LibPETSc.VecSetType(petsclib, vec, "seq")  # Set vector type
LibPETSc.VecSetFromOptions(petsclib, vec)

# Set values
LibPETSc.VecSet(petsclib, vec, 1.0)

# Query the vector
size = LibPETSc.VecGetSize(petsclib, vec)
println("Vector size: $size")

# Get values (0-based indexing!)
idx = PetscInt(0)
val = LibPETSc.VecGetValues(petsclib, vec, PetscInt(1), [idx])
println("Value at index 0: $(val[1])")

# Clean up
LibPETSc.VecDestroy(petsclib, vec)
PETSc.finalize(petsclib)
```

## Function Reference

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/Vec_wrappers.jl", "autowrapped/Vecs_wrappers.jl", "autowrapped/VecTagger_wrappers.jl"]
Order   = [:function]
```
