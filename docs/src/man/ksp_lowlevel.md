# KSP - Low-level Interface

This page documents the low-level, automatically wrapped PETSc KSP (Krylov Subspace Methods) functions available through `LibPETSc`. These functions provide direct access to the PETSc C API.

For the high-level Julia interface, see [KSP](@ref).

## Overview

The KSP interface includes:
- **Core KSP operations**: Linear solver creation, configuration, and solution (~244 functions in `KSP_wrappers.jl`)
- **KSPGuess**: Initial guess generation for iterative solvers (~12 functions in `KSPGuess_wrappers.jl`)

## Usage

All low-level KSP functions require a `petsclib` parameter as the first argument:

```julia
using PETSc

# Get the library instance
petsclib = PETSc.petsclibs[1]
PETSc.initialize(petsclib)

# Assume mat, b, x are already created

# Create KSP solver
ksp = LibPETSc.KSPCreate(petsclib, MPI.COMM_SELF)
LibPETSc.KSPSetOperators(petsclib, ksp, mat, mat)
LibPETSc.KSPSetFromOptions(petsclib, ksp)

# Solve Ax = b
LibPETSc.KSPSolve(petsclib, ksp, b, x)

# Get convergence info
reason = Ref{PetscInt}()
LibPETSc.KSPGetConvergedReason(petsclib, ksp, reason)

iterations = Ref{PetscInt}()
LibPETSc.KSPGetIterationNumber(petsclib, ksp, iterations)

# Clean up
LibPETSc.KSPDestroy(petsclib, ksp)
PETSc.finalize(petsclib)
```

## Function Reference

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/KSP_wrappers.jl", "autowrapped/KSPGuess_wrappers.jl"]
Order   = [:function]
```
