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

# Initialize MPI and PETSc
MPI.Init()
petsclib = PETSc.getlib()
PETSc.initialize(petsclib)

# Create a DMSwarm in PIC mode
swarm = LibPETSc.DMCreate(petsclib, MPI.COMM_WORLD)
LibPETSc.DMSetType(petsclib, swarm, LibPETSc.DMSWARM)
# Set the geometric/topological dimension for the swarm (required)
LibPETSc.DMSetDimension(petsclib, swarm, 1)

# Set swarm type to PIC
LibPETSc.DMSwarmSetType(petsclib, swarm, LibPETSc.DMSWARM_PIC)

# Create background mesh (here we use a 1D DMDA for simplicity)
da = LibPETSc.DMDACreate1d(petsclib, MPI.COMM_WORLD, LibPETSc.DM_BOUNDARY_NONE, 10, 1, 1, C_NULL)

# Set the cell DM for PIC
LibPETSc.DMSwarmSetCellDM(petsclib, swarm, da)

# Register a particle field `velocity` with blocksize 3 (x,y,z components)
LibPETSc.DMSwarmRegisterPetscDatatypeField(
    petsclib, swarm,
    "velocity", 3, LibPETSc.PETSC_DOUBLE
)

# Finalize field registration
LibPETSc.DMSwarmFinalizeFieldRegister(petsclib, swarm)

# Set number of particles
nparticles = 100
LibPETSc.DMSwarmSetLocalSizes(petsclib, swarm, nparticles, 0)

# Access and set particle data using DMSwarmGetField / DMSwarmRestoreField.
# `DMSwarmGetField` returns the block size, the data type and a raw pointer to the
# field storage (blocksize * nparticles entries); wrap it with `unsafe_wrap`.
blocksize, dtype, ptr = LibPETSc.DMSwarmGetField(petsclib, swarm, "velocity")
@assert dtype == LibPETSc.PETSC_DOUBLE
data = unsafe_wrap(Array, Ptr{Float64}(ptr), blocksize * nparticles)
# Initialize velocity values
for i in eachindex(data)
    data[i] = 0.1 * i
end
# Restore the field when done (unlocks internal storage)
LibPETSc.DMSwarmRestoreField(petsclib, swarm, "velocity", blocksize, dtype)

# Cleanup
LibPETSc.DMDestroy(petsclib, swarm)
LibPETSc.DMDestroy(petsclib, da)

# Finalize PETSc and MPI
PETSc.finalize(petsclib)
MPI.Finalize()
```

## DMSwarm Functions

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["DM_wrappers.jl"]
Order   = [:function]
Filter = t -> startswith(string(nameof(t)), "DMSwarm")
```
