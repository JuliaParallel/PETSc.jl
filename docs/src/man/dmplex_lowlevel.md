# DMPlex (Unstructured Meshes)

DMPlex manages unstructured meshes using a flexible topology representation. It supports finite element methods, finite volume methods, and general mesh-based computations.

## Overview

DMPlex provides:
- General unstructured mesh topology (simplices, tensor products, polyhedra)
- Point-based topological representation (vertices, edges, faces, cells)
- Mesh partitioning and distribution
- Overlap and ghost cell management
- Section-based field layout
- Support for FEM assembly via PetscDS and PetscFE
- Mesh refinement and adaptation

## Basic Usage Pattern

```julia
using PETSc, MPI

MPI.Init()
petsclib = PETSc.getlib()
PETSc.initialize(petsclib)

# Create a simple 2D quadrilateral mesh
dm = LibPETSc.DMPlexCreateBoxMesh(
    petsclib,
    MPI.COMM_WORLD,
    2,                           # dimension
    LibPETSc.PETSC_FALSE,        # simplex (false = tensor/quad)
    [5, 5],                      # faces per dimension
    [0.0, 0.0],                  # lower bounds
    [1.0, 1.0],                  # upper bounds
    [LibPETSc.DM_BOUNDARY_NONE, LibPETSc.DM_BOUNDARY_NONE],  # boundary types
    LibPETSc.PETSC_TRUE,         # interpolate
    0,                           # localizationHeight
    LibPETSc.PETSC_FALSE         # sparseLocalize
)

# Distribute mesh across processes (for serial, this is a no-op)
dmParallel = LibPETSc.PetscDM(C_NULL, petsclib)
LibPETSc.DMPlexDistribute(petsclib, dm, 0, C_NULL, dmParallel)
if dmParallel.ptr != C_NULL
    LibPETSc.DMDestroy(petsclib, dm)
    dm = dmParallel
end

# Set up
LibPETSc.DMSetFromOptions(petsclib, dm)
LibPETSc.DMSetUp(petsclib, dm)

# Create section to define field layout
section = Ref{LibPETSc.PetscSection}()
LibPETSc.DMGetLocalSection(petsclib, dm, section)

# Create vectors and matrices
x = LibPETSc.DMCreateGlobalVector(petsclib, dm)

# Cleanup
LibPETSc.VecDestroy(petsclib, x)
LibPETSc.DMDestroy(petsclib, dm)
PETSc.finalize(petsclib)
MPI.Finalize()
```

## DMPlex Functions

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["DM_wrappers.jl"]
Order   = [:function]
Filter = t -> startswith(string(nameof(t)), "DMPlex")
```
