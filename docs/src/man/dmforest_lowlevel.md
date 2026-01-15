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

# Initialize MPI and PETSc
MPI.Init()
petsclib = PETSc.getlib()
PETSc.initialize(petsclib)

# Create a DMForest
forest = LibPETSc.DMCreate(petsclib, MPI.COMM_WORLD)
LibPETSc.DMSetType(petsclib, forest, "forest")  # String convenience wrapper
# Set^ the geometric/topological dimension (e.g., 2 for a surface)
LibPETSc.DMSetDimension(petsclib, forest, 2)

# Set forest topology (e.g., unit square/cube)
# Set topology by name (e.g., "brick")
LibPETSc.DMForestSetTopology(petsclib, forest, "brick")

# Set base DM (initial coarse mesh) using a simple DMPlex box mesh
base_dm = LibPETSc.DMPlexCreateBoxMesh(
    petsclib, MPI.COMM_WORLD, 2, LibPETSc.PETSC_FALSE,
    [2, 2], [0.0, 0.0], [1.0, 1.0], [LibPETSc.DM_BOUNDARY_NONE, LibPETSc.DM_BOUNDARY_NONE], LibPETSc.PETSC_TRUE, 0, LibPETSc.PETSC_FALSE
)
LibPETSc.DMForestSetBaseDM(petsclib, forest, base_dm)

# Set initial refinement level
LibPETSc.DMForestSetInitialRefinement(petsclib, forest, 2)

# Ensure adjacency dimension and partition overlap are non-negative (some builds may leave them unset)
LibPETSc.DMForestSetAdjacencyDimension(petsclib, forest, 0)
LibPETSc.DMForestSetPartitionOverlap(petsclib, forest, 0)

# Set up
# (Skip DMSetFromOptions in this simple example to avoid parsing unexpected runtime options)
LibPETSc.DMSetUp(petsclib, forest)

# Adapt based on some criterion
# Create adapted forest (returns via out-parameter)
tdm = LibPETSc.PetscDM(C_NULL, petsclib)
LibPETSc.DMForestTemplate(petsclib, forest, MPI.COMM_WORLD, tdm)
adapted = tdm
# `adapted` now points to the adapted DM (if any) and should be checked before use.

# Cleanup
LibPETSc.DMDestroy(petsclib, adapted)
LibPETSc.DMDestroy(petsclib, forest)
LibPETSc.DMDestroy(petsclib, base_dm)

# Finalize PETSc and MPI
PETSc.finalize(petsclib)
MPI.Finalize()
```

## DMForest Functions

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["DM_wrappers.jl"]
Order   = [:function]
Filter = t -> startswith(string(nameof(t)), "DMForest")
```
