# PetscOptions (Options Database) - Low-level Interface

The options database holds the `-ksp_type gmres`-style settings PETSc objects read in their
`SetFromOptions` calls. The high-level `PETSc.PetscOptions` type (see [Utilities](utilities.md))
covers the common cases; the functions here give full control.

A `PetscOptions` handle with a NULL pointer denotes the global options database.

## Basic Usage

```julia
using PETSc, MPI
MPI.Init()
petsclib = PETSc.getlib()
PETSc.initialize(petsclib)

opts = LibPETSc.PetscOptions{typeof(petsclib)}()    # NULL handle: the global database

LibPETSc.PetscOptionsSetValue(petsclib, opts, "-ksp_type", "gmres")
LibPETSc.PetscOptionsSetValue(petsclib, opts, "-ksp_max_it", "50")

maxit, set = LibPETSc.PetscOptionsGetInt(petsclib, opts, "", "-ksp_max_it")   # (50, PETSC_TRUE)
has = LibPETSc.PetscOptionsHasName(petsclib, opts, "", "-ksp_type")           # PETSC_TRUE

LibPETSc.PetscOptionsClearValue(petsclib, opts, "-ksp_max_it")
PETSc.finalize(petsclib)
```

`pre` is an optional prefix (`"mg_levels_"`), `""` for none. The `PetscOptionsGet*` functions
return the value and a `PetscBool` telling whether the option was set.

## Function Reference

### PetscOptions

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscOptions_wrappers.jl"]
Order   = [:function]
```

### Options helpers (PetscOptionsGetViewer, ...)

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscOptions_addons_wrappers.jl"]
Order   = [:function]
```
