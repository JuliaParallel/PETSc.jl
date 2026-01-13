# SNES - Low-level Interface

This page documents the low-level, automatically wrapped PETSc SNES (Scalable Nonlinear Equations Solvers) functions available through `LibPETSc`. These functions provide direct access to the PETSc C API.

For the high-level Julia interface, see [SNES](@ref).

## Overview

The SNES interface includes:
- **Core SNES operations**: Nonlinear solver creation, configuration, and solution (~283 functions in `SNES_wrappers.jl`)
- **SNESLineSearch**: Line search methods for nonlinear solvers (~50 functions in `SNESLineSearch_wrappers.jl`)

## Usage

All low-level SNES functions require a `petsclib` parameter as the first argument:

```julia
using PETSc

# Get the library instance
petsclib = PETSc.petsclibs[1]
PETSc.initialize(petsclib)

# Create SNES solver
snes = LibPETSc.SNESCreate(petsclib, MPI.COMM_SELF)

# Set function and Jacobian (callbacks would need to be defined)
# LibPETSc.SNESSetFunction(petsclib, snes, r, compute_function, ctx)
# LibPETSc.SNESSetJacobian(petsclib, snes, J, J, compute_jacobian, ctx)

# Configure from command line
LibPETSc.SNESSetFromOptions(petsclib, snes)

# Solve
# LibPETSc.SNESSolve(petsclib, snes, C_NULL, x)

# Get convergence info
reason = Ref{PetscInt}()
LibPETSc.SNESGetConvergedReason(petsclib, snes, reason)

iterations = Ref{PetscInt}()
LibPETSc.SNESGetIterationNumber(petsclib, snes, iterations)

# Clean up
LibPETSc.SNESDestroy(petsclib, snes)
PETSc.finalize(petsclib)
```

## Function Reference

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/SNES_wrappers.jl", "autowrapped/SNESLineSearch_wrappers.jl"]
Order   = [:function]
```
