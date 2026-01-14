# DMSwarm (Particle Methods)

DMSwarm manages particle-based methods including particle-in-cell (PIC), smoothed particle hydrodynamics (SPH), and general particle simulations.

## Overview

DMSwarm provides:
- Dynamic particle creation and migration
- Field registration for particle data
- Particle-mesh coupling
- Particle migration across MPI processes
- Support for PIC and general particle methods
- Cell-based particle management

## Basic Usage Pattern

```julia
using PETSc, MPI

petsclib = PETSc.getlib()

# Create a DMSwarm in PIC mode
swarm = Ref{LibPETSc.CDM}()
LibPETSc.DMCreate(petsclib, MPI.COMM_WORLD, swarm)
LibPETSc.DMSetType(petsclib, swarm[], "swarm")  # String convenience wrapper

# Set swarm type to PIC
LibPETSc.DMSwarmSetType(petsclib, swarm[], LibPETSc.DMSWARM_PIC)

# Create background mesh (e.g., DMDA or DMPlex)
mesh = Ref{LibPETSc.CDM}()
# ... create mesh ...

# Set the cell DM for PIC
LibPETSc.DMSwarmSetCellDM(petsclib, swarm[], mesh[])

# Register particle fields
LibPETSc.DMSwarmRegisterPetscDatatypeField(
    petsclib, swarm[],
    "velocity", 3, LibPETSc.PETSC_REAL
)

# Finalize field registration
LibPETSc.DMSwarmFinalizeFieldRegister(petsclib, swarm[])

# Set number of particles
nparticles = 1000
LibPETSc.DMSwarmSetLocalSizes(petsclib, swarm[], nparticles, 0)

# Access and set particle data
# Use DMSwarmGetField and DMSwarmRestoreField

# Cleanup
LibPETSc.DMDestroy(petsclib, swarm)
LibPETSc.DMDestroy(petsclib, mesh)
```

## DMSwarm Functions

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["DM_wrappers.jl"]
Order   = [:function]
Filter = t -> startswith(string(nameof(t)), "DMSwarm")
```
