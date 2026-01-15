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

# Initialize MPI and PETSc
MPI.Init()
petsclib = PETSc.getlib()
PETSc.initialize(petsclib)

# Create a DMNetwork (returns a PetscDM)
network = LibPETSc.DMNetworkCreate(petsclib, MPI.COMM_WORLD)

# Set network sizes
# For simple examples, provide a name and an explicit edge list vector
nedges = 2
nvertices = 3
LibPETSc.DMNetworkSetNumSubNetworks(petsclib, network, 1, 1)
edgelist = [1, 2, 2, 3]  # pairs of vertex global indices
LibPETSc.DMNetworkAddSubnetwork(petsclib, network, "main", nedges, edgelist)

# Register components (size in bytes; use appropriate struct size in real examples)
# Register a component; the function returns a component key (PetscInt)
compkey = LibPETSc.DMNetworkRegisterComponent(petsclib, network, "bus", Csize_t(0))

# Add components to vertices/edges using DMNetworkAddComponent (omitted here)

# Finalize network and set up
LibPETSc.DMNetworkLayoutSetUp(petsclib, network)
LibPETSc.DMSetUp(petsclib, network)

# Create vectors
x = LibPETSc.DMCreateGlobalVector(petsclib, network)

# Cleanup
LibPETSc.VecDestroy(petsclib, x)
LibPETSc.DMDestroy(petsclib, network)

# Finalize PETSc and MPI
PETSc.finalize(petsclib)
MPI.Finalize()
```

## DMNetwork Functions

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["DM_wrappers.jl"]
Order   = [:function]
Filter = t -> startswith(string(nameof(t)), "DMNetwork")
```
