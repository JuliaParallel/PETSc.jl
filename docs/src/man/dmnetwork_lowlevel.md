# DMNetwork (Network/Graph Problems)

DMNetwork manages network and graph-based problems, such as power grids, transportation networks, and general graph computations with edge and vertex data.

## Overview

DMNetwork provides:
- Network topology with edges and vertices
- Component registration for edges and vertices
- Coupling of network problems with other DM types
- Distribution across MPI processes
- Support for algebraic constraints at vertices

## Basic Usage Pattern

```julia
using PETSc, MPI

petsclib = PETSc.getlib()

# Create a DMNetwork
network = Ref{LibPETSc.CDM}()
LibPETSc.DMNetworkCreate(petsclib, MPI.COMM_WORLD, network)

# Set network sizes
nedges = 10
nvertices = 8
LibPETSc.DMNetworkSetNumSubNetworks(petsclib, network[], 1, 1)
LibPETSc.DMNetworkAddSubnetwork(
    petsclib, network[],
    C_NULL, nedges, C_NULL, nvertices
)

# Register components
# e.g., for power grid: buses, generators, transmission lines
compkey = Ref{LibPETSc.PetscInt}()
LibPETSc.DMNetworkRegisterComponent(
    petsclib, network[],
    "bus", sizeof(BusData), compkey
)

# Add components to vertices/edges
# Use DMNetworkAddComponent

# Finalize network
LibPETSc.DMNetworkLayoutSetUp(petsclib, network[])

# Set up
LibPETSc.DMSetUp(petsclib, network[])

# Create vectors
x = Ref{LibPETSc.CVec}()
LibPETSc.DMCreateGlobalVector(petsclib, network[], x)

# Cleanup
LibPETSc.VecDestroy(petsclib, x)
LibPETSc.DMDestroy(petsclib, network)
```

## DMNetwork Functions

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["DM_wrappers.jl"]
Order   = [:function]
Filter = t -> startswith(string(nameof(t)), "DMNetwork")
```
