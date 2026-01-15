using Test
using PETSc, MPI

# Initialize MPI & PETSc only if needed
mpi_started = false
if !MPI.Initialized()
    MPI.Init()
    mpi_started = true
end
try
    petsclib = PETSc.getlib()
    petsc_started = false
    if !PETSc.initialized(petsclib)
        PETSc.initialize(petsclib)
        petsc_started = true
    end

    @testset "DMNetwork basic" begin
        network = LibPETSc.DMNetworkCreate(petsclib, MPI.COMM_WORLD)
        @test network.ptr != C_NULL

        nedges = 2
        nvertices = 3
        LibPETSc.DMNetworkSetNumSubNetworks(petsclib, network, 1, 1)
        edgelist = [1,2,2,3]
        netnum = LibPETSc.DMNetworkAddSubnetwork(petsclib, network, "main", nedges, edgelist)
        @test netnum >= 0

        compkey = LibPETSc.DMNetworkRegisterComponent(petsclib, network, "bus", Csize_t(0))
        @test compkey >= 0

        LibPETSc.DMNetworkLayoutSetUp(petsclib, network)
        LibPETSc.DMSetUp(petsclib, network)

        x = LibPETSc.DMCreateGlobalVector(petsclib, network)
        @test x.ptr != C_NULL

        LibPETSc.VecDestroy(petsclib, x)
        LibPETSc.DMDestroy(petsclib, network)
    end

    if petsc_started
        PETSc.finalize(petsclib)
    end
catch e
    try if petsc_started; PETSc.finalize(petsclib); end catch end
    rethrow(e)
finally
    # Do not call MPI.Finalize() from tests; test harness manages MPI lifecycle
    nothing
end