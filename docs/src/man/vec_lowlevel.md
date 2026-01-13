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

# Use low-level functions
vec = LibPETSc.VecCreate(petsclib, MPI.COMM_SELF)
LibPETSc.VecSetSizes(petsclib, vec, 10, 10)
LibPETSc.VecSetFromOptions(petsclib, vec)

# ... work with vec ...

LibPETSc.VecDestroy(petsclib, vec)
PETSc.finalize(petsclib)
```

## Function Reference

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/Vec_wrappers.jl", "autowrapped/Vecs_wrappers.jl", "autowrapped/VecTagger_wrappers.jl"]
Order   = [:function]
```
