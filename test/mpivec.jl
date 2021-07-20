using Test
using MPI
MPI.Initialized() || MPI.Init()
using PETSc
using LinearAlgebra: norm

@testset "VecMPI" begin
    comm = MPI.COMM_WORLD
    mpisize = MPI.Comm_size(comm)
    mpirank = MPI.Comm_rank(comm)

    # local first and last value
    n0 = sum(10 .+ (0:mpirank)) - 10
    n1 = sum(10 .+ (0:(mpirank + 1))) - 11

    # global last value
    ne = sum(10 .+ (0:mpisize)) - 11
    exact_length = ne + 1
    exact_norm = norm(0:ne)

    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt

        for version in 1:2
            if version == 1
                # Create using local size
                petsc_x = PETSc.VecMPI(petsclib, comm, n1 - n0 + 1)

                # check the data ownership
                rng = n0:n1
                @test rng == PETSc.ownershiprange(petsc_x, false)
            else
                # Create using local size
                petsc_x = PETSc.VecMPI(
                    petsclib,
                    comm,
                    PETSc.PETSC_DECIDE;
                    global_length = exact_length,
                )

                # get the data ownership
                rng = PETSc.ownershiprange(petsc_x, false)
            end

            # insert some values
            julia_x = PetscScalar.(rng)
            inds = PetscInt.(rng)
            # 0-based
            PETSc.setvalues!(petsc_x, inds, julia_x, PETSc.INSERT_VALUES)

            PETSc.assemblybegin!(petsc_x)
            PETSc.assemblyend!(petsc_x)

            @test length(petsc_x) == exact_length
            vec_norm = norm(petsc_x)
            @test exact_norm ≈ vec_norm
            # 1-based
            @test petsc_x[rng .+ 1] == julia_x

            PETSc.with_unsafe_localarray!(petsc_x) do x
                @test x == julia_x
            end

            PETSc.destroy(petsc_x)

        end
        PETSc.finalize(petsclib)
    end
end

@testset "VecGhost" begin
    comm = MPI.COMM_WORLD
    mpisize = MPI.Comm_size(comm)
    mpirank = MPI.Comm_rank(comm)

    # local first and last value
    n0 = sum(10 .+ (0:mpirank)) - 10
    n1 = sum(10 .+ (0:(mpirank + 1))) - 11
    local_length = n1 - n0 + 1

    # global last value
    ne = sum(10 .+ (0:mpisize)) - 11
    exact_length = ne + 1

    # exact norm
    exact_norm = norm(0:ne)

    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt

        # two ghost to left and right
        ghost = PetscInt.([])
        if mpirank != 0
            ghost = PetscInt.([ghost..., n0 - 2, n0 - 1])
        end
        if mpirank != mpisize - 1
            ghost = PetscInt.([ghost..., n1 + 1, n1 + 2])
        end

        for version in 1:2
            if version == 1
                # Create using local size
                petsc_x = PETSc.VecGhost(petsclib, comm, local_length, ghost)

                rng = PetscInt.(PETSc.ownershiprange(petsc_x))
                julia_x = PetscScalar.(rng)
                PETSc.setvalues!(petsc_x, rng, julia_x, PETSc.INSERT_VALUES)
            else
                # Create using local size
                continue
            end

            PETSc.assemblybegin!(petsc_x)
            PETSc.assemblyend!(petsc_x)

            @test length(petsc_x) == exact_length

            vec_norm = norm(petsc_x)
            @test exact_norm ≈ vec_norm

            PETSc.withlocalform(petsc_x) do l_x
                @test length(l_x) == local_length + length(ghost)

                # Check the ghost has propogated
                if length(ghost) > 0
                    vals = zeros(PetscScalar, length(ghost))
                    inds = PetscInt.(local_length - 1 .+ (1:length(ghost)))
                    # Initially we have pushed the numbers so shouldn't match
                    PETSc.getvalues!(vals, l_x, inds)
                    @test !(vals == ghost)

                    # propagate the ghost
                    PETSc.ghostupdatebegin!(petsc_x)
                    PETSc.ghostupdateend!(petsc_x)

                    # Recheck the numbers
                    PETSc.getvalues!(vals, l_x, inds)
                    @test vals == ghost
                end
            end

            PETSc.destroy(petsc_x)
        end
        PETSc.finalize(petsclib)
    end
end
