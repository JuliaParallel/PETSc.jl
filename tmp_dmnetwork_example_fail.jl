using PETSc, MPI

# Running the DMNetwork minimal example from docs (as-is) to reproduce crash

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
compkey = Ref{LibPETSc.PetscInt}()
LibPETSc.DMNetworkRegisterComponent(
    petsclib, network[],
    "bus", 0, compkey
)

# Finalize network
LibPETSc.DMNetworkLayoutSetUp(petsclib, network[])

# Set up
LibPETSc.DMSetUp(petsclib, network[])

# Create vectors
x = Ref{LibPETSc.CVec}()
LibPETSc.DMCreateGlobalVector(petsclib, network[], x)

println("Finished DMNetwork example (no init).")