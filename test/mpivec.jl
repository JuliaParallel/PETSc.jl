using Test
using MPI
if !Sys.iswindows()
    MPI.Initialized() || MPI.Init()
end
using PETSc
using LinearAlgebra: norm


@testset "VecMPI" begin
    comm = Sys.iswindows() ? LibPETSc.PETSC_COMM_SELF : MPI.COMM_WORLD
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
        #petsclib = PETSc.petsclibs[1]
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt

        for version in 1:2
            if version == 1
                # Create using local size
                petsc_x = LibPETSc.VecCreateMPI(petsclib, comm, PetscInt(n1 - n0 + 1), PetscInt(LibPETSc.PETSC_DECIDE))

                # check the data ownership
                rng = n0:n1
                @test rng == PETSc.ownershiprange(petsc_x, false)
            else
                # Create using global size
                petsc_x = LibPETSc.VecCreateMPI(
                    petsclib,
                    comm,
                    PetscInt(LibPETSc.PETSC_DECIDE),
                    PetscInt(exact_length),
                )

                # get the data ownership
                rng = PETSc.ownershiprange(petsc_x, false)
            end

            # insert some values
            julia_x = PetscScalar.(rng)
            inds = PetscInt.(rng)

            # NOTE: PETSc is 0-based
            LibPETSc.VecSetValues(petsclib,petsc_x, PetscInt(length(inds)), inds, julia_x, PETSc.INSERT_VALUES)
            PETSc.assemble!(petsc_x)

            @test length(petsc_x) == exact_length
            
            # use this instead of norm(petsc_x), as it can deal with parallel vectors
            vec_norm1 = LibPETSc.VecNorm(petsclib,petsc_x, LibPETSc.NORM_2)
            @test exact_norm ≈ vec_norm1
            
            vec_norm = norm(petsc_x)    # multiple dispatch version of the line above
            @test exact_norm ≈ vec_norm
            
            # 1-based
            @test petsc_x[rng .+ 1] == julia_x

            PETSc.withlocalarray!(petsc_x) do x
                @test x == julia_x
            end

            @test "mpi" == PETSc.type(petsc_x)
            PETSc.destroy(petsc_x)

        end
        PETSc.finalize(petsclib)
    end
end


@testset "VecGhost" begin
    comm = Sys.iswindows() ? LibPETSc.PETSC_COMM_SELF : MPI.COMM_WORLD
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
    exact_norm = norm(1:ne+1)

    for petsclib in PETSc.petsclibs
        #petsclib = PETSc.petsclibs[1]
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt

        # two ghost to left and right
        ghost = PetscInt.([])
        if mpirank != 0
            ghost = PetscInt.([ghost..., n0 - 2, n0 - 1]  )
        end
        if mpirank != mpisize - 1
            ghost = PetscInt.([ghost..., n1 + 1, n1 + 2]  )
        end
        
        for version in 1:2
            if version == 1
                # Create using local size
                petsc_x = LibPETSc.VecCreateGhost(petsclib, comm, PetscInt(local_length), PetscInt(LibPETSc.PETSC_DECIDE), PetscInt(length(ghost)), ghost)



                rng = PetscInt.(PETSc.ownershiprange(petsc_x))
                julia_x = PetscScalar.(rng)
                LibPETSc.VecSetValues(petsclib,petsc_x, PetscInt(length(rng)), PetscInt.(rng .- 1), julia_x, PETSc.INSERT_VALUES)

            else
                # Create using local size
                continue
            end

            PETSc.assemble!(petsc_x)
            @test length(petsc_x) == exact_length

            vec_norm = norm(petsc_x)
            @test exact_norm ≈ vec_norm

            petsc_x_local = LibPETSc.VecGhostGetLocalForm(petsclib,petsc_x) 
            
          #  LibPETSc.VecGetLocalVector(petsclib,petsc_x, petsc_x_local)

            PETSc.withlocalarray!(petsc_x_local) do l_x
                @test length(l_x) == local_length + length(ghost)

                # Check the ghost has propagated
                if length(ghost) > 0
                    vals = zeros(PetscScalar, length(ghost))
                    inds = PetscInt.(local_length - 1 .+ (1:length(ghost)))
                   
                    # Initially we have pushed the numbers so shouldn't match
                    vals = l_x[inds .+ 1]  # +1 for 1-based indexing
                    @test !(vals == ghost)

                    # propagate the ghost
                    PETSc.ghostupdatebegin!(petsc_x)
                    PETSc.ghostupdateend!(petsc_x)

                    # Recheck the numbers
                    vals = l_x[inds .+ 1]  
                    @test vals == (ghost .+1)
                end
            end
            

            PETSc.destroy(petsc_x)
        end
        PETSc.finalize(petsclib)
    end
end
