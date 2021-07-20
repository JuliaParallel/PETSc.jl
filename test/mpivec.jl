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
