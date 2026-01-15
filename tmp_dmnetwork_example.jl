using PETSc, MPI

MPI.Init()
try
    petsclib = PETSc.getlib()
    PETSc.initialize(petsclib)

    # Create a DMNetwork (returns a PetscDM)
    network = LibPETSc.DMNetworkCreate(petsclib, MPI.COMM_WORLD)

    # Set network sizes
    nedges = 2
    nvertices = 3
    LibPETSc.DMNetworkSetNumSubNetworks(petsclib, network, 1, 1)
    edgelist = [1, 2, 2, 3]  # pairs of vertex global indices
    LibPETSc.DMNetworkAddSubnetwork(petsclib, network, "main", nedges, edgelist)

    # Register components (size in bytes; use appropriate struct size in real examples)
    compkey = LibPETSc.DMNetworkRegisterComponent(petsclib, network, "bus", Csize_t(0))
    println("Registered component key: ", compkey)

    LibPETSc.DMNetworkLayoutSetUp(petsclib, network)
    LibPETSc.DMSetUp(petsclib, network)

    # Create vectors
    x = LibPETSc.DMCreateGlobalVector(petsclib, network)

    println("DMNetwork example completed successfully: ", x)

    # Cleanup
    LibPETSc.VecDestroy(petsclib, x)
    LibPETSc.DMDestroy(petsclib, network)

    PETSc.finalize(petsclib)
catch e
    @error "DMNetwork example failed" exception=(e, catch_backtrace())
    try PETSc.finalize(petsclib) catch end
    rethrow(e)
finally
    try MPI.Finalize() catch end
end