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
        dof_per_nodec= (4, 3)
        stencil_width = 3

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

                dmnew = PETSc.DMStag(
                    petsclib,
                    dm,
                    dof_per_nodec,
                )

                @test PETSc.gettype(dm) == "stag"
                @test PETSc.gettype(dmnew) == "stag"
                @test PETSc.getdimension(dm) == 1
                @test PETSc.getdof(dm) == (3, 4)
                @test PETSc.getdof(dmnew) == (4, 3)
                @test PETSc.globalsize(dm) ===
                      (global_size, PetscInt(1), PetscInt(1))
                @test size(dm) === (global_size, PetscInt(1), PetscInt(1))
                @test PETSc.localsize(dm) ===
                      (points_per_proc[mpirank + 1], PetscInt(1), PetscInt(1))

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

@testset "DMStagCreate2d" begin
    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)
    global_size_x = 100
    global_size_y = 45
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt

        dof_per_node = (3, 4, 5)
        dof_per_nodec= (4,3)
        stencil_width = 5

        for boundary_type_y in instances(PETSc.DMBoundaryType),
            boundary_type_x in instances(PETSc.DMBoundaryType)

            # DMStag cannot be used with these boundary types
            boundary_type_x ∈
            (PETSc.DM_BOUNDARY_MIRROR, PETSc.DM_BOUNDARY_TWIST) && continue
            boundary_type_y ∈
            (PETSc.DM_BOUNDARY_MIRROR, PETSc.DM_BOUNDARY_TWIST) && continue

            @testset "$boundary_type_x, $boundary_type_y" begin
                dm = PETSc.DMStag(
                    petsclib,
                    comm,
                    (boundary_type_x, boundary_type_y),
                    (global_size_x, global_size_y),
                    dof_per_node,
                    stencil_width,
                )

                dmnew = PETSc.DMStag(
                    petsclib,
                    dm,
                    dof_per_nodec,
                )

                @test PETSc.gettype(dm) == "stag"
                @test PETSc.gettype(dmnew) == "stag"
                @test PETSc.getdimension(dm) == 2
                @test PETSc.getdof(dm) == (3, 4, 5)
                @test PETSc.getdof(dmnew) == (4, 3, 0)
                @test PETSc.globalsize(dm) === (
                    PetscInt(global_size_x),
                    PetscInt(global_size_y),
                    PetscInt(1),
                )
                @test size(dm) === (
                    PetscInt(global_size_x),
                    PetscInt(global_size_y),
                    PetscInt(1),
                )

                corners = PETSc.getcorners(dm)
                ghost_corners = PETSc.getghostcorners(dm)
                isfirst = PETSc.isfirstrank(dm)
                islast = PETSc.islastrank(dm)

                bt = (boundary_type_x, boundary_type_y, PETSc.DM_BOUNDARY_NONE)
                for d in 1:3
                    # Check left side ghost_corners and corner
                    if d == 3
                        @test corners.lower[d] == ghost_corners.lower[d]
                    elseif !isfirst[d] || bt[d] != PETSc.DM_BOUNDARY_NONE
                        @test corners.lower[d] ==
                              ghost_corners.lower[d] + stencil_width
                    else
                        @test corners.lower[d] == ghost_corners.lower[d]
                    end

                    # Check right side ghost_corners and corner
                    if d == 3
                        @test corners.upper[d] == ghost_corners.upper[d]
                    elseif bt[d] != PETSc.DM_BOUNDARY_NONE
                        @test corners.upper[d] ==
                              ghost_corners.upper[d] - stencil_width
                    elseif !islast[d]
                        @test corners.upper[d] ==
                              ghost_corners.upper[d] - stencil_width
                    else
                        @test corners.upper[d] == ghost_corners.upper[d] - 1
                    end

                    # Check the nextra
                    if islast[d] && bt[d] != PETSc.DM_BOUNDARY_PERIODIC
                        @test corners.nextra[d] == 1
                    else
                        @test corners.nextra[d] == 0
                    end
                end

                @test PETSc.boundarytypes(dm) ===
                      (boundary_type_x, boundary_type_y, PETSc.DM_BOUNDARY_NONE)

                PETSc.destroy(dm)
            end
        end
    end
end

@testset "DMStagCreate3d" begin
    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)
    global_size_x = 20
    global_size_y = 25
    global_size_z = 30
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt

        dof_per_node = (2, 3, 4, 5)
        dof_per_nodec= (4,3)
        stencil_width = 5

        for boundary_type_z in instances(PETSc.DMBoundaryType),
            boundary_type_y in instances(PETSc.DMBoundaryType),
            boundary_type_x in instances(PETSc.DMBoundaryType)

            # DMStag cannot be used with these boundary types
            boundary_type_x ∈
            (PETSc.DM_BOUNDARY_MIRROR, PETSc.DM_BOUNDARY_TWIST) && continue
            boundary_type_y ∈
            (PETSc.DM_BOUNDARY_MIRROR, PETSc.DM_BOUNDARY_TWIST) && continue
            boundary_type_z ∈
            (PETSc.DM_BOUNDARY_MIRROR, PETSc.DM_BOUNDARY_TWIST) && continue

            @testset "$boundary_type_x, $boundary_type_y, $boundary_type_z" begin
                dm = PETSc.DMStag(
                    petsclib,
                    comm,
                    (boundary_type_x, boundary_type_y, boundary_type_z),
                    (global_size_x, global_size_y, global_size_z),
                    dof_per_node,
                    stencil_width,
                )

                dmnew = PETSc.DMStag(
                    petsclib,
                    dm,
                    dof_per_nodec,
                )

                @test PETSc.gettype(dm) == "stag"
                @test PETSc.gettype(dmnew) == "stag"
                @test PETSc.getdimension(dm) == 3
                @test PETSc.getdof(dm) == (2, 3, 4, 5)
                @test PETSc.getdof(dmnew) == (4, 3, 0, 0)
                @test PETSc.globalsize(dm) === (
                    PetscInt(global_size_x),
                    PetscInt(global_size_y),
                    PetscInt(global_size_z),
                )
                @test size(dm) === (
                    PetscInt(global_size_x),
                    PetscInt(global_size_y),
                    PetscInt(global_size_z),
                )

                corners = PETSc.getcorners(dm)
                ghost_corners = PETSc.getghostcorners(dm)
                isfirst = PETSc.isfirstrank(dm)
                islast = PETSc.islastrank(dm)

                bt = (boundary_type_x, boundary_type_y, boundary_type_z)
                for d in 1:3
                    # Check left side ghost_corners and corner
                    if !isfirst[d] || bt[d] != PETSc.DM_BOUNDARY_NONE
                        @test corners.lower[d] ==
                              ghost_corners.lower[d] + stencil_width
                    else
                        @test corners.lower[d] == ghost_corners.lower[d]
                    end

                    # Check right side ghost_corners and corner
                    if bt[d] != PETSc.DM_BOUNDARY_NONE
                        @test corners.upper[d] ==
                              ghost_corners.upper[d] - stencil_width
                    elseif !islast[d]
                        @test corners.upper[d] ==
                              ghost_corners.upper[d] - stencil_width
                    else
                        @test corners.upper[d] == ghost_corners.upper[d] - 1
                    end

                    # Check the nextra
                    if islast[d] && bt[d] != PETSc.DM_BOUNDARY_PERIODIC
                        @test corners.nextra[d] == 1
                    else
                        @test corners.nextra[d] == 0
                    end
                end

                @test PETSc.boundarytypes(dm) ===
                      (boundary_type_x, boundary_type_y, boundary_type_z)

                PETSc.destroy(dm)
            end
        end
    end
end
