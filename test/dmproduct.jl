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

    @testset "DMProduct basic" begin
        dm1 = LibPETSc.DMPlexCreateBoxMesh(
            petsclib, MPI.COMM_WORLD, 2, LibPETSc.PETSC_FALSE,
            [1,1], [0.0,0.0], [1.0,1.0], [LibPETSc.DM_BOUNDARY_NONE, LibPETSc.DM_BOUNDARY_NONE], LibPETSc.PETSC_TRUE, 0, LibPETSc.PETSC_FALSE
        )
        dm2 = LibPETSc.DMPlexCreateBoxMesh(
            petsclib, MPI.COMM_WORLD, 2, LibPETSc.PETSC_FALSE,
            [1,1], [0.0,0.0], [1.0,1.0], [LibPETSc.DM_BOUNDARY_NONE, LibPETSc.DM_BOUNDARY_NONE], LibPETSc.PETSC_TRUE, 0, LibPETSc.PETSC_FALSE
        )

        product = LibPETSc.DMCreate(petsclib, MPI.COMM_WORLD)
        LibPETSc.DMSetType(petsclib, product, "product")
        # set dimension to number of component DMs and map dimension indices
        LibPETSc.DMSetDimension(petsclib, product, 2)
        LibPETSc.DMProductSetDimensionIndex(petsclib, product, 0, 0)
        LibPETSc.DMProductSetDimensionIndex(petsclib, product, 1, 0)

        @test product.ptr != C_NULL

        LibPETSc.DMProductSetDM(petsclib, product, 0, dm1)
        LibPETSc.DMProductSetDM(petsclib, product, 1, dm2)

        LibPETSc.DMSetUp(petsclib, product)

        LibPETSc.DMDestroy(petsclib, product)
        LibPETSc.DMDestroy(petsclib, dm1)
        LibPETSc.DMDestroy(petsclib, dm2)
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