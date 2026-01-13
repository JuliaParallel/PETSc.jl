# DM (Domain Management)

The DM (Domain Management) object encapsulates the relationship between a mesh data structure and the algebraic objects (vectors, matrices). It provides a common interface for structured and unstructured meshes, particles, networks, and other geometric/topological structures.

## Overview

The DM abstraction in PETSc manages:
- Topology and geometry of computational domains
- Mapping between geometric entities and degrees of freedom
- Distribution of data across MPI processes
- Multi-level representations for multigrid
- Field definitions and discretizations

## Basic Usage Pattern

```julia
using PETSc, MPI

# Initialize
petsclib = PETSc.getlib()

# Create a DM (specific type created via DMDACreate, DMPlexCreate, etc.)
# See subtype-specific pages for creation functions

# Setup the DM
LibPETSc.DMSetFromOptions(petsclib, dm)
LibPETSc.DMSetUp(petsclib, dm)

# Create global vectors
x = Ref{LibPETSc.CVec}()
LibPETSc.DMCreateGlobalVector(petsclib, dm, x)

# Create matrices
A = Ref{LibPETSc.CMat}()
LibPETSc.DMCreateMatrix(petsclib, dm, A)

# Use the DM...

# Cleanup
LibPETSc.VecDestroy(petsclib, x)
LibPETSc.MatDestroy(petsclib, A)
LibPETSc.DMDestroy(petsclib, dm)
```

## DM Subtypes

PETSc provides several DM implementations for different problem types:

- **[DMDA](dmda_lowlevel.md)**: Structured grids (1D, 2D, 3D) with regular topology
- **[DMPlex](dmplex_lowlevel.md)**: Unstructured meshes using Sieve/DMPlex topology
- **[DMStag](dmstag_lowlevel.md)**: Staggered grids for staggered finite differences
- **[DMSwarm](dmswarm_lowlevel.md)**: Particle methods and particle-in-cell
- **[DMForest](dmforest_lowlevel.md)**: Adaptive mesh refinement with forest-of-octrees
- **[DMNetwork](dmnetwork_lowlevel.md)**: Network/graph-based problems
- **[DMShell and others](dmshell_lowlevel.md)**: User-defined and composite DM types

## General DM Functions

The following functions work across all DM types:

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["DM_wrappers.jl", "DMaddons_wrappers.jl"]
Order   = [:function]
Filter = t -> begin
    name = string(nameof(t))
    # Include only general DM functions, exclude subtype-specific ones
    startswith(name, "DM") && 
    !startswith(name, "DMDA") &&
    !startswith(name, "DMPlex") &&
    !startswith(name, "DMStag") &&
    !startswith(name, "DMSwarm") &&
    !startswith(name, "DMForest") &&
    !startswith(name, "DMNetwork") &&
    !startswith(name, "DMShell") &&
    !startswith(name, "DMProduct") &&
    !startswith(name, "DMRedundant") &&
    !startswith(name, "DMComposite")
end
```
