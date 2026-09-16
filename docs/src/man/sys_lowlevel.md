# Sys (Runtime, Memory, Strings, Errors) - Low-level Interface

`Sys_wrappers.jl` collects the PETSc functions that do not belong to a class: initialization
and version queries, memory and timing utilities, error handling, string and sorting helpers,
MPI helpers (`PetscSplitOwnership`, `PetscGlobalMinMax*`, ...), binary I/O and the
`PetscInfo` diagnostics.

PETSc.jl calls `PetscInitialize`/`PetscFinalize` for you through `PETSc.initialize(petsclib)`
and `PETSc.finalize(petsclib)`; do not call them directly.

## Basic Usage

```julia
using PETSc, MPI
MPI.Init()
petsclib = PETSc.getlib()
PETSc.initialize(petsclib)

major, minor, subminor, release = LibPETSc.PetscGetVersionNumber(petsclib)
mem = LibPETSc.PetscMemoryGetCurrentUsage(petsclib)       # bytes used by this process
t0 = LibPETSc.PetscGetCPUTime(petsclib)                    # PetscTime is a C inline; not wrapped
LibPETSc.PetscSleep(petsclib, 0.01)
elapsed = LibPETSc.PetscGetCPUTime(petsclib) - t0

# Split N unknowns over the processes of a communicator (n = PETSC_DECIDE)
n, N = LibPETSc.PetscSplitOwnership(petsclib, MPI.COMM_WORLD, LibPETSc.PETSC_DECIDE, 100)

PETSc.finalize(petsclib)
```

## Function Reference

### Sys

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/Sys_wrappers.jl"]
Order   = [:function]
```
