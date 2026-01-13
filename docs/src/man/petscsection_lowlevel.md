# PetscSection - Low-level Interface

The PetscSection component provides a flexible mechanism for describing the layout of degrees of freedom (DOFs) in a discretization, particularly for finite element methods and other multi-field problems where different points have different numbers of DOFs.

## Overview

PetscSection is used to:
- **Define DOF layouts**: Specify how many DOFs exist at each point (vertex, edge, face, cell)
- **Multi-field problems**: Handle systems with multiple variables per point
- **Local-to-global mapping**: Convert between local and global numbering
- **FEM discretizations**: Manage DOF distributions for finite elements
- **DM integration**: Work with DM objects for mesh-based computations

A section maps from points (mesh entities like vertices, cells) to DOF ranges, answering:
- How many DOFs are at point `p`?
- What is the offset of the DOFs for point `p`?
- Which field do specific DOFs belong to?

## Basic Usage

```julia
using PETSc

# Initialize PETSc
petsclib = PETSc.getlib()

# Create a section
section = Ref{LibPETSc.PetscSection}()
LibPETSc.PetscSectionCreate(petsclib, MPI.COMM_SELF, section)

# Set chart: range of valid point indices [pStart, pEnd)
LibPETSc.PetscSectionSetChart(petsclib, section[], 0, 10)

# Set DOF count for each point
for p in 0:9
    num_dofs = (p < 4) ? 1 : 2  # Different DOFs per point
    LibPETSc.PetscSectionSetDof(petsclib, section[], p, num_dofs)
end

# Setup: compute offsets
LibPETSc.PetscSectionSetUp(petsclib, section[])

# Query the section
dof = Ref{PetscInt}()
LibPETSc.PetscSectionGetDof(petsclib, section[], 5, dof)  # Get DOF count for point 5

offset = Ref{PetscInt}()
LibPETSc.PetscSectionGetOffset(petsclib, section[], 5, offset)  # Get offset for point 5

# Get total storage size
storage_size = Ref{PetscInt}()
LibPETSc.PetscSectionGetStorageSize(petsclib, section[], storage_size)

# Cleanup
LibPETSc.PetscSectionDestroy(petsclib, section)
```

## Multi-Field Sections

For problems with multiple fields (e.g., velocity + pressure):

```julia
# Create section with 2 fields
section = Ref{LibPETSc.PetscSection}()
LibPETSc.PetscSectionCreate(petsclib, MPI.COMM_SELF, section)
LibPETSc.PetscSectionSetNumFields(petsclib, section[], 2)

# Set field names
LibPETSc.PetscSectionSetFieldName(petsclib, section[], 0, "velocity")
LibPETSc.PetscSectionSetFieldName(petsclib, section[], 1, "pressure")

# Set chart
LibPETSc.PetscSectionSetChart(petsclib, section[], 0, 10)

# Set field components: velocity has 3 components (vx, vy, vz), pressure has 1
LibPETSc.PetscSectionSetFieldComponents(petsclib, section[], 0, 3)
LibPETSc.PetscSectionSetFieldComponents(petsclib, section[], 1, 1)

# Set DOFs per field per point
for p in 0:9
    LibPETSc.PetscSectionSetFieldDof(petsclib, section[], p, 0, 3)  # 3 velocity DOFs
    LibPETSc.PetscSectionSetFieldDof(petsclib, section[], p, 1, 1)  # 1 pressure DOF
    LibPETSc.PetscSectionSetDof(petsclib, section[], p, 4)          # Total: 4 DOFs
end

LibPETSc.PetscSectionSetUp(petsclib, section[])
```

## Constrained DOFs

Mark certain DOFs as constrained (e.g., for boundary conditions):

```julia
# Set chart and DOFs...

# Mark 2 DOFs at point 0 as constrained
LibPETSc.PetscSectionSetConstraintDof(petsclib, section[], 0, 2)

# Specify which DOFs are constrained (e.g., DOFs 0 and 2)
constrained = Int32[0, 2]
LibPETSc.PetscSectionSetConstraintIndices(petsclib, section[], 0, constrained)

LibPETSc.PetscSectionSetUp(petsclib, section[])

# Query unconstrained DOF count
unconstrained_dof = Ref{PetscInt}()
LibPETSc.PetscSectionGetConstrainedStorageSize(petsclib, section[], unconstrained_dof)
```

## Integration with DM

Sections are commonly used with DM objects:

```julia
# Get section from DM
dm_section = Ref{LibPETSc.PetscSection}()
# LibPETSc.DMGetSection(petsclib, dm, dm_section)

# Set section on DM
# LibPETSc.DMSetSection(petsclib, dm, section[])

# Get local section (ghosted)
# local_section = Ref{LibPETSc.PetscSection}()
# LibPETSc.DMGetLocalSection(petsclib, dm, local_section)
```

## Common Workflows

### 1. FEM Discretization

```julia
# Create section for Q1 finite elements on structured grid
# - Vertices: 1 DOF each
# - Cells: 0 DOFs

section = Ref{LibPETSc.PetscSection}()
LibPETSc.PetscSectionCreate(petsclib, comm, section)
LibPETSc.PetscSectionSetChart(petsclib, section[], vStart, cEnd)

# Set DOFs (vStart to vEnd are vertices, vEnd to cEnd are cells)
for v in vStart:vEnd-1
    LibPETSc.PetscSectionSetDof(petsclib, section[], v, 1)
end

LibPETSc.PetscSectionSetUp(petsclib, section[])
```

### 2. Point Closure

Get all DOFs in the closure of a point (point + its boundary):

```julia
# Get closure DOFs for a cell
closure_size = Ref{PetscInt}()
closure = Ref{Ptr{PetscInt}}()
# LibPETSc.DMPlexGetTransitiveClosure(petsclib, dm, cell, PETSC_TRUE, closure_size, closure)

# Map closure points to DOF offsets using section
# ... use section to get DOF offsets for each point in closure ...

# LibPETSc.DMPlexRestoreTransitiveClosure(petsclib, dm, cell, PETSC_TRUE, closure_size, closure)
```

## Querying Section Properties

```julia
# Get total number of fields
num_fields = Ref{PetscInt}()
LibPETSc.PetscSectionGetNumFields(petsclib, section[], num_fields)

# Get field components
components = Ref{PetscInt}()
LibPETSc.PetscSectionGetFieldComponents(petsclib, section[], field, components)

# Get field name
name = Ref{Ptr{Cchar}}()
LibPETSc.PetscSectionGetFieldName(petsclib, section[], field, name)

# Get maximum DOF count across all points
max_dof = Ref{PetscInt}()
LibPETSc.PetscSectionGetMaxDof(petsclib, section[], max_dof)
```

## Cloning and Permutation

```julia
# Clone a section
new_section = Ref{LibPETSc.PetscSection}()
LibPETSc.PetscSectionClone(petsclib, section[], new_section)

# Create permuted section
perm = Int32[3, 1, 2, 0]  # Permutation
perm_section = Ref{LibPETSc.PetscSection}()
LibPETSc.PetscSectionPermute(petsclib, section[], perm, perm_section)
```

## Function Reference

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/PetscSection_wrappers.jl"]
Order   = [:function]
```
