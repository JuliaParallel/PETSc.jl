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

petsclib = PETSc.getlib()

# Create a simple 2D triangular mesh
dm = Ref{LibPETSc.CDM}()
LibPETSc.DMPlexCreateBoxMesh(
    petsclib,
    MPI.COMM_WORLD,
    2,                           # dimension
    LibPETSc.PETSC_FALSE,        # simplex (false = tensor/quad)
    [5, 5],                      # faces per dimension
    C_NULL, C_NULL,              # lower/upper bounds (NULL for unit box)
    C_NULL,                      # boundary types
    LibPETSc.PETSC_TRUE,         # interpolate
    dm
)

# Distribute mesh across processes
dmDist = Ref{LibPETSc.CDM}()
LibPETSc.DMPlexDistribute(petsclib, dm[], 0, C_NULL, dmDist)
if dmDist[] != C_NULL
    LibPETSc.DMDestroy(petsclib, dm)
    dm[] = dmDist[]
end

# Set up
LibPETSc.DMSetFromOptions(petsclib, dm[])
LibPETSc.DMSetUp(petsclib, dm[])

# Create section to define field layout
section = Ref{LibPETSc.PetscSection}()
LibPETSc.DMGetLocalSection(petsclib, dm[], section)

# Create vectors and matrices
x = Ref{LibPETSc.CVec}()
LibPETSc.DMCreateGlobalVector(petsclib, dm[], x)

# Cleanup
LibPETSc.VecDestroy(petsclib, x)
LibPETSc.DMDestroy(petsclib, dm)
```

## DMPlex Functions

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["DM_wrappers.jl"]
Order   = [:function]
Filter = t -> startswith(string(nameof(t)), "DMPlex")
```
