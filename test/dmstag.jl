using Test
using PETSc, MPI
MPI.Initialized() || MPI.Init()

@testset "DMStagCreate1d" begin
    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt

        dof_per_node = (3, 4)
        stencil_width = 5

        for boundary_type in instances(PETSc.DMBoundaryType)

            # DMStag cannot be used with these boundary types
            boundary_type == PETSc.DM_BOUNDARY_MIRROR && continue
            boundary_type == PETSc.DM_BOUNDARY_TWIST && continue

            @testset "$boundary_type" begin

                # We test both setting and not setting the point distribution
                points_per_proc = [PetscInt(10 + i) for i in 1:mpisize]
                proc_global_offsets =
                    [PetscInt(0), accumulate(+, points_per_proc)...]

                global_size = proc_global_offsets[end]

                dm = PETSc.DMStag(
                    petsclib,
                    comm,
                    (boundary_type,),
                    (global_size,),
                    dof_per_node,
                    stencil_width;
                    points_per_proc = (points_per_proc,),
                )

                @test PETSc.gettype(dm) == "stag"
                @test PETSc.getdimension(dm) == 1
                @test PETSc.globalsize(dm) ===
                      (global_size, PetscInt(0), PetscInt(0))
                @test size(dm) === (global_size, PetscInt(0), PetscInt(0))
                @test PETSc.localsize(dm) ===
                      (points_per_proc[mpirank + 1], PetscInt(0), PetscInt(0))

                corners = PETSc.getcorners(dm)
                @test corners.lower ==
                      CartesianIndex(proc_global_offsets[mpirank + 1] + 1, 1, 1)
                @test corners.upper ==
                      CartesianIndex(proc_global_offsets[mpirank + 2], 1, 1)
                @test corners.size == (points_per_proc[mpirank + 1], 1, 1)

                # Check the extra and first / last rank
                map(PETSc.islastrank(dm), corners.nextra) do b, n
                    @test b && boundary_type != PETSc.DM_BOUNDARY_PERIODIC ?
                          n == 1 : n == 0
                end

                gl =
                    boundary_type == PETSc.DM_BOUNDARY_NONE && mpirank == 0 ?
                    0 : stencil_width
                gr =
                    boundary_type == PETSc.DM_BOUNDARY_NONE &&
                    mpirank == mpisize - 1 ? 1 : stencil_width
                ghost_corners = PETSc.getghostcorners(dm)

                @test ghost_corners.lower == CartesianIndex(
                    proc_global_offsets[mpirank + 1] + 1 - gl,
                    1,
                    1,
                )
                @test ghost_corners.upper == CartesianIndex(
                    proc_global_offsets[mpirank + 2] + gr,
                    1,
                    1,
                )
                @test ghost_corners.size ==
                      (points_per_proc[mpirank + 1] + gl + gr, 1, 1)

                @test PETSc.boundarytypes(dm) === (
                    boundary_type,
                    PETSc.DM_BOUNDARY_NONE,
                    PETSc.DM_BOUNDARY_NONE,
                )
                PETSc.destroy(dm)
            end
        end
    end
end
