using Test
using PETSc, MPI
MPI.Initialized() || MPI.Init()

@testset "DMDACreate1D" begin
    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt

        # Loop over all boundary types and try to use them
        for boundary_type in instances(PETSc.DMBoundaryType)
            @testset "$boundary_type" begin
                dof_per_node = 4
                stencil_width = 5

                # We test both setting and not setting the point distribution
                points_per_proc = [PetscInt(10 + i) for i in 1:mpisize]
                proc_global_offsets =
                    [PetscInt(0), accumulate(+, points_per_proc)...]

                global_size = proc_global_offsets[end]

                # left and right ghost points
                gl =
                    boundary_type == PETSc.DM_BOUNDARY_NONE && mpirank == 0 ?
                    0 : stencil_width
                gr =
                    boundary_type == PETSc.DM_BOUNDARY_NONE &&
                    mpirank == mpisize - 1 ? 0 : stencil_width

                # Set the points
                da = PETSc.DMDA(
                    petsclib,
                    comm,
                    (boundary_type,),
                    (global_size,),
                    dof_per_node,
                    stencil_width;
                    points_per_proc = (points_per_proc,),
                )

                @test PETSc.gettype(da) == "da"
                @test PETSc.getdimension(da) == 1

                da_info = PETSc.getinfo(da)

                @test da_info.dim == 1
                @test da_info.global_size == (global_size, 1, 1)
                @test da_info.procs_per_dim == (mpisize, 1, 1)
                @test da_info.boundary_type == (
                    boundary_type,
                    PETSc.DM_BOUNDARY_NONE,
                    PETSc.DM_BOUNDARY_NONE,
                )
                @test da_info.stencil_type == PETSc.DMDA_STENCIL_BOX
                @test da_info.stencil_width == stencil_width

                corners = PETSc.getcorners(da)
                @test corners.lower ==
                      CartesianIndex(proc_global_offsets[mpirank + 1] + 1, 1, 1)
                @test corners.upper ==
                      CartesianIndex(proc_global_offsets[mpirank + 2], 1, 1)
                @test corners.size == (points_per_proc[mpirank + 1], 1, 1)

                ghost_corners = PETSc.getghostcorners(da)
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

                PETSc.destroy(da)

                # Do not set the points and test option parsing
                da_refine = 2
                da = PETSc.DMDA(
                    petsclib,
                    comm,
                    (boundary_type,),
                    (global_size,),
                    dof_per_node,
                    stencil_width;
                    da_refine = da_refine,
                )
                @test PETSc.gettype(da) == "da"
                @test PETSc.getdimension(da) == 1

                da_info = PETSc.getinfo(da)

                @test da_info.dim == 1
                if boundary_type == PETSc.DM_BOUNDARY_PERIODIC
                    @test da_info.global_size ==
                          (2^da_refine * global_size, 1, 1)
                else
                    @test da_info.global_size ==
                          (2^da_refine * (global_size - 1) + 1, 1, 1)
                end
                @test da_info.procs_per_dim == (mpisize, 1, 1)
                @test da_info.boundary_type == (
                    boundary_type,
                    PETSc.DM_BOUNDARY_NONE,
                    PETSc.DM_BOUNDARY_NONE,
                )
                @test da_info.stencil_type == PETSc.DMDA_STENCIL_BOX
                @test da_info.stencil_width == stencil_width
                # In this case we cannot check the numbers locally
                PETSc.destroy(da)
                #=
                # TODO: Need a better test?
                ksp = PETSc.KSP(da)
                @test PETSc.gettype(ksp) == "gmres"
                =#
            end
        end
        PETSc.finalize(petsclib)
    end
end

@testset "DMDACreate2D" begin
    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)
    global_size_x = 100
    global_size_y = 45
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt

        # Loop over all boundary types and stencil types
        for stencil_type in instances(PETSc.DMDAStencilType),
            boundary_type_y in instances(PETSc.DMBoundaryType),
            boundary_type_x in instances(PETSc.DMBoundaryType)

            # skip unsupported stencils
            stencil_type == PETSc.DMDA_STENCIL_BOX &&
                (
                    boundary_type_x == PETSc.DM_BOUNDARY_MIRROR ||
                    boundary_type_y == PETSc.DM_BOUNDARY_MIRROR
                ) &&
                continue

            @testset "$boundary_type_x, $boundary_type_y, $stencil_type" begin
                dof_per_node = 4
                stencil_width = 5

                # Set the points
                da = PETSc.DMDA(
                    petsclib,
                    comm,
                    (boundary_type_x, boundary_type_y),
                    (global_size_x, global_size_y),
                    dof_per_node,
                    stencil_width,
                    stencil_type,
                )
                @test PETSc.gettype(da) == "da"
                @test PETSc.getdimension(da) == 2

                da_info = PETSc.getinfo(da)

                @test da_info.global_size == (global_size_x, global_size_y, 1)
                @test da_info.dim == 2
                @test prod(da_info.procs_per_dim) == mpisize
                @test da_info.boundary_type ==
                      (boundary_type_x, boundary_type_y, PETSc.DM_BOUNDARY_NONE)
                @test da_info.stencil_type == stencil_type
                @test da_info.stencil_width == stencil_width

                # test refinement
                da_refine = 2
                da = PETSc.DMDA(
                    petsclib,
                    comm,
                    (boundary_type_x, boundary_type_y),
                    (global_size_x, global_size_y),
                    dof_per_node,
                    stencil_width,
                    stencil_type;
                    da_refine = da_refine,
                )
                @test PETSc.gettype(da) == "da"
                @test PETSc.getdimension(da) == 2

                da_info = PETSc.getinfo(da)

                # Compute refined global size
                ref_global_size_x =
                    boundary_type_x == PETSc.DM_BOUNDARY_PERIODIC ?
                    2^da_refine * global_size_x :
                    2^da_refine * (global_size_x - 1) + 1
                ref_global_size_y =
                    boundary_type_y == PETSc.DM_BOUNDARY_PERIODIC ?
                    2^da_refine * global_size_y :
                    2^da_refine * (global_size_y - 1) + 1

                @test da_info.global_size ==
                      (ref_global_size_x, ref_global_size_y, 1)
                @test prod(da_info.procs_per_dim) == mpisize
                @test da_info.boundary_type ==
                      (boundary_type_x, boundary_type_y, PETSc.DM_BOUNDARY_NONE)
                @test da_info.stencil_type == stencil_type
                @test da_info.stencil_width == stencil_width

                # TODO: Test with specific distribution of processors and sizes

                # TODO: Need a better test?
                #=
                ksp = PETSc.KSP(da)
                @test PETSc.gettype(ksp) == "gmres"
                =#

                PETSc.destroy(da)
            end
        end
        PETSc.finalize(petsclib)
    end
end

@testset "DMDACreate3D" begin
    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)
    global_size_x = 12
    global_size_y = 13
    global_size_z = 14
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt
        # Loop over all boundary types and stencil types
        for stencil_type in instances(PETSc.DMDAStencilType),
            boundary_type_z in instances(PETSc.DMBoundaryType),
            boundary_type_y in instances(PETSc.DMBoundaryType),
            boundary_type_x in instances(PETSc.DMBoundaryType)

            stencil_type == PETSc.DMDA_STENCIL_BOX &&
                (
                    boundary_type_x == PETSc.DM_BOUNDARY_MIRROR ||
                    boundary_type_y == PETSc.DM_BOUNDARY_MIRROR ||
                    boundary_type_z == PETSc.DM_BOUNDARY_MIRROR
                ) &&
                continue

            @testset "$boundary_type_x, $boundary_type_y, $boundary_type_z, $stencil_type" begin
                dof_per_node = 4
                stencil_width = 2

                # Set the points
                da = PETSc.DMDA(
                    petsclib,
                    comm,
                    (boundary_type_x, boundary_type_y, boundary_type_z),
                    (global_size_x, global_size_y, global_size_z),
                    dof_per_node,
                    stencil_width,
                    stencil_type,
                )
                @test PETSc.gettype(da) == "da"
                @test PETSc.getdimension(da) == 3

                da_info = PETSc.getinfo(da)

                @test da_info.global_size ==
                      (global_size_x, global_size_y, global_size_z)
                @test da_info.dim == 3
                @test prod(da_info.procs_per_dim) == mpisize
                @test da_info.boundary_type ==
                      (boundary_type_x, boundary_type_y, boundary_type_z)
                @test da_info.stencil_type == stencil_type
                @test da_info.stencil_width == stencil_width

                # test refinement
                da_refine = 2
                da = PETSc.DMDA(
                    petsclib,
                    comm,
                    (boundary_type_x, boundary_type_y, boundary_type_z),
                    (global_size_x, global_size_y, global_size_z),
                    dof_per_node,
                    stencil_width,
                    stencil_type;
                    da_refine = da_refine,
                )
                @test PETSc.gettype(da) == "da"
                @test PETSc.getdimension(da) == 3

                da_info = PETSc.getinfo(da)

                # Compute refined global size
                ref_global_size_x =
                    boundary_type_x == PETSc.DM_BOUNDARY_PERIODIC ?
                    2^da_refine * global_size_x :
                    2^da_refine * (global_size_x - 1) + 1
                ref_global_size_y =
                    boundary_type_y == PETSc.DM_BOUNDARY_PERIODIC ?
                    2^da_refine * global_size_y :
                    2^da_refine * (global_size_y - 1) + 1
                ref_global_size_z =
                    boundary_type_z == PETSc.DM_BOUNDARY_PERIODIC ?
                    2^da_refine * global_size_z :
                    2^da_refine * (global_size_z - 1) + 1

                @test da_info.global_size ==
                      (ref_global_size_x, ref_global_size_y, ref_global_size_z)
                @test prod(da_info.procs_per_dim) == mpisize
                @test da_info.boundary_type ==
                      (boundary_type_x, boundary_type_y, boundary_type_z)
                @test da_info.stencil_type == stencil_type
                @test da_info.stencil_width == stencil_width

                # TODO: Test with specific distribution of processors and sizes

                # TODO: Need a better test?
                #=
                ksp = PETSc.KSP(da)
                @test PETSc.gettype(ksp) == "gmres"
                =#
            end
        end
        PETSc.finalize(petsclib)
    end
end

@testset "DM MatAIJ" begin
    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)
    # Just check a couple libraries
    for petsclib in PETSc.petsclibs[1:2]
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt
        for dim in 1:3
            boundary_type = ntuple(_ -> PETSc.DM_BOUNDARY_NONE, dim)
            dof_per_node = 1
            stencil_width = 1
            number_points = 10
            global_size = ntuple(i -> (10 + i) * mpisize, 3)
            # Set the points
            da = PETSc.DMDA(
                petsclib,
                comm,
                boundary_type,
                global_size[1:dim],
                dof_per_node,
                stencil_width,
                PETSc.DMDA_STENCIL_STAR;
                # Force processor distribution so that we know how to go from
                # Stencil to Linear easily
                processors = (ntuple(i -> 1, dim - 1)..., mpisize),
            )
            mat = PETSc.MatAIJ(da)

            # Build the dim-dimensional Laplacian FD matrix
            corners = PETSc.getcorners(da)

            for i in (corners.lower):(corners.upper)
                for d in 1:dim
                    if i[d] - 1 > 1
                        j = i - CartesianIndex(d == 1, d == 2, d == 3)
                        mat[i, j] = 1
                    end
                    if i[d] + 1 < global_size[d]
                        j = i + CartesianIndex(d == 1, d == 2, d == 3)
                        mat[i, j] = 1
                    end
                end
                mat[i, i] = -2dim
            end

            PETSc.assemble!(mat)

            ind = LinearIndices(ntuple(i -> 1:global_size[i], 3))
            for ci in (corners.lower):(corners.upper)
                i = ind[ci]
                @test mat[i, i] == -2dim
                for d in 1:dim
                    if ci[d] - 1 > 1
                        j = ind[ci - CartesianIndex(d == 1, d == 2, d == 3)]
                        @test mat[i, j] == 1
                    end
                    if ci[d] + 1 < global_size[d]
                        j = ind[ci + CartesianIndex(d == 1, d == 2, d == 3)]
                        @test mat[i, j] == 1
                    end
                end
            end

            PETSc.destroy(mat)
            PETSc.destroy(da)
        end
        PETSc.finalize(petsclib)
    end
end

@testset "DM Vectors and Coordinates" begin
    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscInt = petsclib.PetscInt

        boundary_type = PETSc.DM_BOUNDARY_NONE
        dof_per_node = 1
        stencil_width = 1
        number_points = 10
        points_per_proc = [PetscInt(10) for i in 1:mpisize]
        global_size = sum(points_per_proc)

        # Set the points
        da = PETSc.DMDA(
            petsclib,
            comm,
            (boundary_type,),
            (global_size,),
            dof_per_node,
            stencil_width;
            points_per_proc = (points_per_proc,),
        )

        corners = PETSc.getcorners(da)

        # Create the local and global vectors
        local_vec = PETSc.DMLocalVec(da)
        global_vec = PETSc.DMGlobalVec(da)
        bot_val = 0
        top_val = 0

        # Fill everything with some data
        fill!(local_vec, mpirank)
        fill!(global_vec, mpisize)

        # add the local values to the global values
        PETSc.update!(global_vec, local_vec, PETSc.ADD_VALUES)

        # end points added with neighbor due to ghost of size 1
        bot_val = mpisize + mpirank + (mpirank == 0 ? 0 : mpirank - 1)
        top_val = mpisize + mpirank + (mpirank == mpisize - 1 ? 0 : mpirank + 1)
        @test global_vec[corners.lower[1]] == bot_val
        @test global_vec[corners.upper[1]] == top_val

        # Center is just self plus the global
        for i in (corners.lower[1] + 1):(corners.upper[1] - 1)
            @test global_vec[i] == mpisize + mpirank
        end

        # reset the local values with the global values
        PETSc.update!(local_vec, global_vec, PETSc.INSERT_VALUES)

        # My first value and my ghost should be the bot/top values
        @test local_vec[1] == bot_val
        @test local_vec[2] == bot_val
        @test local_vec[end - 1] == top_val
        @test local_vec[end] == top_val

        # interior is just self plus the global
        for i in 3:(length(local_vec) - 2)
            @test local_vec[i] == mpisize + mpirank
        end

        #=
        # Test DM Coordinates
        coord_da = PETSc.getcoordinateDM(da)
        # Crank it up to 11!
        xmin, xmax = 0, 11
        PETSc.setuniformcoordinates!(coord_da, (xmin,), (xmax,))
        coord_vec = PETSc.getcoordinateslocal(coord_da)
        Δx = (xmax - xmin) / (global_size - 1)

        # Figure out the values we should have in the coordinate vector
        ghost_lower = corners.lower[1] - (mpirank == 0 ? 0 : 1)
        ghost_upper = corners.upper[1] + (mpirank == mpisize - 1 ? 0 : 1)
        for (loc, glo) in enumerate(ghost_lower:ghost_upper)
            @test coord_vec[loc] ≈ (glo - 1) * Δx
        end
        =#
        PETSc.finalize(petsclib)
    end
end

nothing
