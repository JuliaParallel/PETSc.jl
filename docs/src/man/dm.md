# DM

The `DM` module provides the base functionality for managing distributed data structures in PETSc. It serves as a foundation for various grid managers.

## Overview

A `DM` object encapsulates the topology and data layout of a computational grid, enabling:
- Parallel data distribution across MPI processes
- Ghost point management for communication
- Creation of vectors and matrices with appropriate parallel layouts
- Multigrid hierarchy management

## DM Types in PETSc

PETSc provides several DM implementations for different mesh types:

### High-Level Interface Available in PETSc.jl

| DM Type | Description | Status |
|---------|-------------|--------|
| **DMDA** | Distributed arrays for structured grids (1D/2D/3D) | ✅ Full support |
| **DMStag** | Staggered grids for finite volume/difference methods | ✅ Full support |
| **DMPlex** | Unstructured meshes + full FEM workflow (Gmsh, FE spaces, callbacks, VTK) | ✅ Full support |

### Flavour is a type

`DMDA`, `DMStag` and `DMPlex` are Julia types, and constructing one returns that
type:

```julia
da = PETSc.DMDA(petsclib, comm, (PETSc.DM_BOUNDARY_NONE,), (10,), 1, 1)
da isa PETSc.DMDA{typeof(petsclib), 1}     # true: flavour and dimension
```

All three are subtypes of `LibPETSc.AbstractPetscDM`, which is what a function
working on any DM takes. `corners`, `ghost_corners`, `local_indices`,
`set_uniform_coordinates!` and `Base.size` differ by flavour, so they are
ordinary methods on `DMDA` and `DMStag` rather than one function comparing the
string `DMGetType` returns.

`LibPETSc.PetscDM` remains the low-level handle: `LibPETSc.DMCreate` and friends
have to return something before the flavour is known. Turn one into a typed
handle with [`PETSc.narrow`](@ref):

```julia
d = PETSc.narrow(dm)          # DMDA{L,N}, DMStag{L,N}, DMPlex{L}, or dm unchanged
```

`narrow` queries PETSc, so its return type is a wide `Union` and the call is a
dynamic dispatch. One dispatch is cheap; propagating an abstractly-typed DM
through a hot loop is not, so narrow once behind a function barrier. The
accessors that hand back a DM PETSc owns — `dm(ksp)`, `dm(snes)`, `dm(ts)`,
`coarse_dm` — narrow for you, and what they return is a **borrowed** handle: it
belongs to the object it was asked of, and `destroy!` on it is a no-op.

### Low-Level Interface Only (via LibPETSc)

The following DM types are available through the low-level `LibPETSc` wrapper but do not yet have a convenient high-level Julia interface:

| DM Type | Description | Use Case |
|---------|-------------|----------|
| **DMForest** | Adaptive mesh refinement (AMR) via p4est/p8est | Octree-based adaptivity |
| **DMNetwork** | Graph/network structures | Power grids, pipe networks |
| **DMSwarm** | Particle data management | PIC methods, Lagrangian particles |
| **DMProduct** | Tensor product of DMs | Semi-structured problems |
| **DMSliced** | Sliced representation | Legacy, specialized uses |
| **DMShell** | User-defined DM | Custom implementations |
| **DMComposite** | Composition of multiple DMs | Multi-physics coupling |
| **DMRedundant** | Redundant storage on all ranks | Small coupled systems |

### Using Low-Level DM Types

For DM types without high-level wrappers, you can use the `LibPETSc` module directly:

```julia
using PETSc
using PETSc.LibPETSc

# Example: Create a DMPlex (low-level)
petsclib = PETSc.petsclibs[1]
dm = LibPETSc.DMPlexCreate(petsclib, MPI.COMM_WORLD)
# ... configure using LibPETSc functions ...
LibPETSc.DMDestroy(petsclib, dm)
```

!!! note "Contributing"
    Contributions to add high-level interfaces for additional DM types are welcome! 
    See the [Contributing](@ref) page for guidelines.

## Functions

```@autodocs
Modules = [PETSc]
Pages   = ["dm.jl"]
```
