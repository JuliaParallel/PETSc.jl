using Test
using PETSc, MPI

MPI.Init()
try
    petsclib = PETSc.getlib()
    PETSc.initialize(petsclib)

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

    PETSc.finalize(petsclib)
catch e
    try PETSc.finalize(petsclib) catch end
    rethrow(e)
finally
    try MPI.Finalize() catch end
end