# DMForest (Adaptive Mesh Refinement)

DMForest manages adaptive mesh refinement (AMR) using forest-of-octrees (2D: quadtrees, 3D: octrees). It provides hierarchical mesh adaptation with efficient parallel implementation.

## Overview

DMForest provides:
- Adaptive mesh refinement and coarsening
- Forest-of-octrees topology
- Integration with p4est library
- Multi-level mesh hierarchies for multigrid
- Metric-based adaptation
- Label-based adaptation

## Basic Usage Pattern

```julia
using PETSc, MPI

petsclib = PETSc.getlib()

# Create a DMForest
forest = Ref{LibPETSc.CDM}()
LibPETSc.DMCreate(petsclib, MPI.COMM_WORLD, forest)
LibPETSc.DMSetType(petsclib, forest[], LibPETSc.DMFOREST)

# Set forest topology (e.g., unit square/cube)
LibPETSc.DMForestSetTopology(
    petsclib, forest[],
    LibPETSc.DM_FOREST_TOPOLOGY_BRICK
)

# Set base DM (initial coarse mesh)
base_dm = Ref{LibPETSc.CDM}()
# ... create base DMPlex ...
LibPETSc.DMForestSetBaseDM(petsclib, forest[], base_dm[])

# Set initial refinement level
LibPETSc.DMForestSetInitialRefinement(petsclib, forest[], 2)

# Set up
LibPETSc.DMSetFromOptions(petsclib, forest[])
LibPETSc.DMSetUp(petsclib, forest[])

# Adapt based on some criterion
# Create adapted forest
adapted = Ref{LibPETSc.CDM}()
# ... set adaptation criterion ...
LibPETSc.DMForestTemplate(petsclib, forest[], MPI.COMM_WORLD, adapted)

# Cleanup
LibPETSc.DMDestroy(petsclib, adapted)
LibPETSc.DMDestroy(petsclib, forest)
LibPETSc.DMDestroy(petsclib, base_dm)
```

## DMForest Functions

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["DM_wrappers.jl"]
Order   = [:function]
Filter = t -> startswith(string(nameof(t)), "DMForest")
```
