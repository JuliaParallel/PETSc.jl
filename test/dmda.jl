using Test
using PETSc, MPI
MPI.Initialized() || MPI.Init()
PETSc.initialize()

@testset "DMDACreate1D" begin
    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)
    for petsclib in PETSc.petsclibs
        PetscScalar = PETSc.scalartype(petsclib)
        PetscInt = PETSc.inttype(petsclib)
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
                da = PETSc.DMDACreate1d(
                    petsclib,
                    comm,
                    boundary_type,
                    global_size,
                    dof_per_node,
                    stencil_width,
                    points_per_proc,
                )
                PETSc.DMSetUp!(da)

                @test PETSc.gettype(da) == "da"
                @test PETSc.DMGetDimension(da) == 1

                da_info = PETSc.DMDAGetInfo(da)
                corners = PETSc.DMDAGetCorners(da)
                ghost_corners = PETSc.DMDAGetGhostCorners(da)

                @test da_info.dim == 1
                @test da_info.global_size == [global_size, 1, 1]
                @test da_info.procs_per_dim == [mpisize, 1, 1]
                @test da_info.boundary_type == [
                    boundary_type,
                    PETSc.DM_BOUNDARY_NONE,
                    PETSc.DM_BOUNDARY_NONE,
                ]
                @test da_info.stencil_type == PETSc.DMDA_STENCIL_BOX
                @test da_info.stencil_width == stencil_width
                @test corners.lower ==
                      [proc_global_offsets[mpirank + 1] + 1, 1, 1]
                @test corners.upper == [proc_global_offsets[mpirank + 2], 1, 1]
                @test corners.size == [points_per_proc[mpirank + 1], 1, 1]
                @test ghost_corners.lower ==
                      [proc_global_offsets[mpirank + 1] + 1 - gl, 1, 1]
                @test ghost_corners.upper ==
                      [proc_global_offsets[mpirank + 2] + gr, 1, 1]
                @test ghost_corners.size ==
                      [points_per_proc[mpirank + 1] + gl + gr, 1, 1]

                # Do not set the points and test option parsing
                da_refine = 2
                da = PETSc.DMDACreate1d(
                    petsclib,
                    comm,
                    boundary_type,
                    global_size,
                    dof_per_node,
                    stencil_width,
                    nothing;
                    da_refine = da_refine,
                )
                PETSc.DMSetUp!(da)
                @test PETSc.gettype(da) == "da"
                @test PETSc.DMGetDimension(da) == 1

                da_info = PETSc.DMDAGetInfo(da)

                @test da_info.dim == 1
                if boundary_type == PETSc.DM_BOUNDARY_PERIODIC
                    @test da_info.global_size ==
                          [2^da_refine * global_size, 1, 1]
                else
                    @test da_info.global_size ==
                          [2^da_refine * (global_size - 1) + 1, 1, 1]
                end
                @test da_info.procs_per_dim == [mpisize, 1, 1]
                @test da_info.boundary_type == [
                    boundary_type,
                    PETSc.DM_BOUNDARY_NONE,
                    PETSc.DM_BOUNDARY_NONE,
                ]
                @test da_info.stencil_type == PETSc.DMDA_STENCIL_BOX
                @test da_info.stencil_width == stencil_width
                # In this case we cannot check the numbers locally

                # TODO: Need a better test?
                ksp = PETSc.KSP(da)
                @test PETSc.gettype(ksp) == "gmres"
            end
        end
    end
end
@testset "DMDACreate2D" begin
    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)
    global_size_x = 100
    global_size_y = 45
    for petsclib in PETSc.petsclibs
        PetscScalar = PETSc.scalartype(petsclib)
        PetscInt = PETSc.inttype(petsclib)
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
                da = PETSc.DMDACreate2d(
                    petsclib,
                    comm,
                    boundary_type_x,
                    boundary_type_y,
                    stencil_type,
                    global_size_x,
                    global_size_y,
                    PETSc.PETSC_DECIDE,
                    PETSc.PETSC_DECIDE,
                    dof_per_node,
                    stencil_width,
                    nothing,
                    nothing,
                )
                PETSc.DMSetUp!(da)
                @test PETSc.gettype(da) == "da"
                @test PETSc.DMGetDimension(da) == 2

                da_info = PETSc.DMDAGetInfo(da)

                @test da_info.global_size == [global_size_x, global_size_y, 1]
                @test da_info.dim == 2
                @test prod(da_info.procs_per_dim) == mpisize
                @test da_info.boundary_type ==
                      [boundary_type_x, boundary_type_y, PETSc.DM_BOUNDARY_NONE]
                @test da_info.stencil_type == stencil_type
                @test da_info.stencil_width == stencil_width

                # test refinement
                da_refine = 2
                da = PETSc.DMDACreate2d(
                    petsclib,
                    comm,
                    boundary_type_x,
                    boundary_type_y,
                    stencil_type,
                    global_size_x,
                    global_size_y,
                    PETSc.PETSC_DECIDE,
                    PETSc.PETSC_DECIDE,
                    dof_per_node,
                    stencil_width,
                    nothing,
                    nothing;
                    da_refine = da_refine,
                )
                PETSc.DMSetUp!(da)
                @test PETSc.gettype(da) == "da"
                @test PETSc.DMGetDimension(da) == 2

                da_info = PETSc.DMDAGetInfo(da)

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
                      [ref_global_size_x, ref_global_size_y, 1]
                @test prod(da_info.procs_per_dim) == mpisize
                @test da_info.boundary_type ==
                      [boundary_type_x, boundary_type_y, PETSc.DM_BOUNDARY_NONE]
                @test da_info.stencil_type == stencil_type
                @test da_info.stencil_width == stencil_width

                # TODO: Test with specific distribution of processors and sizes

                # TODO: Need a better test?
                ksp = PETSc.KSP(da)
                @test PETSc.gettype(ksp) == "gmres"
            end
        end
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
        PetscScalar = PETSc.scalartype(petsclib)
        PetscInt = PETSc.inttype(petsclib)
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
                da = PETSc.DMDACreate3d(
                    petsclib,
                    comm,
                    boundary_type_x,
                    boundary_type_y,
                    boundary_type_z,
                    stencil_type,
                    global_size_x,
                    global_size_y,
                    global_size_z,
                    PETSc.PETSC_DECIDE,
                    PETSc.PETSC_DECIDE,
                    PETSc.PETSC_DECIDE,
                    dof_per_node,
                    stencil_width,
                    nothing,
                    nothing,
                    nothing,
                )
                PETSc.DMSetUp!(da)
                @test PETSc.gettype(da) == "da"
                @test PETSc.DMGetDimension(da) == 3

                da_info = PETSc.DMDAGetInfo(da)

                @test da_info.global_size ==
                      [global_size_x, global_size_y, global_size_z]
                @test da_info.dim == 3
                @test prod(da_info.procs_per_dim) == mpisize
                @test da_info.boundary_type ==
                      [boundary_type_x, boundary_type_y, boundary_type_z]
                @test da_info.stencil_type == stencil_type
                @test da_info.stencil_width == stencil_width

                # test refinement
                da_refine = 2
                da = PETSc.DMDACreate3d(
                    petsclib,
                    comm,
                    boundary_type_x,
                    boundary_type_y,
                    boundary_type_z,
                    stencil_type,
                    global_size_x,
                    global_size_y,
                    global_size_z,
                    PETSc.PETSC_DECIDE,
                    PETSc.PETSC_DECIDE,
                    PETSc.PETSC_DECIDE,
                    dof_per_node,
                    stencil_width,
                    nothing,
                    nothing,
                    nothing;
                    da_refine = da_refine,
                )
                PETSc.DMSetUp!(da)
                @test PETSc.gettype(da) == "da"
                @test PETSc.DMGetDimension(da) == 3

                da_info = PETSc.DMDAGetInfo(da)

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
                      [ref_global_size_x, ref_global_size_y, ref_global_size_z]
                @test prod(da_info.procs_per_dim) == mpisize
                @test da_info.boundary_type ==
                      [boundary_type_x, boundary_type_y, boundary_type_z]
                @test da_info.stencil_type == stencil_type
                @test da_info.stencil_width == stencil_width

                # TODO: Test with specific distribution of processors and sizes

                # TODO: Need a better test?
                ksp = PETSc.KSP(da)
                @test PETSc.gettype(ksp) == "gmres"
            end
        end
    end
end

@testset "DMCreateMatrix" begin
    for petsclib in PETSc.petsclibs
        comm = MPI.Comm_dup(MPI.COMM_SELF)
        PetscScalar = PETSc.scalartype(petsclib)
        PetscInt = PETSc.inttype(petsclib)
        boundary_type = PETSc.DM_BOUNDARY_NONE
        dof_per_node = 1
        stencil_width = 1
        number_points = 10
        points_per_proc = [PetscInt(number_points)]
        global_size = points_per_proc[end]
        gl = 0
        gr = 0
        # Set the points
        da = PETSc.DMDACreate1d(
            petsclib,
            comm,
            boundary_type,
            global_size,
            dof_per_node,
            stencil_width,
            points_per_proc,
        )
        PETSc.DMSetUp!(da)
        mat = PETSc.DMCreateMatrix(da)

        # Build the 1-D Laplacian FD matrix
        Sten = PETSc.MatStencil{PetscInt}
        col = Vector{Sten}(undef, 3)
        row = Vector{Sten}(undef, 1)
        val = Vector{PetscScalar}(undef, 3)
        corners = PETSc.DMDAGetCorners(da)

        for i in corners.lower[1]:corners.upper[1]
            row[1] = Sten(i = i)
            num = 1
            fill!(val, 0)
            if i > 1
                val[num] = -1
                col[num] = Sten(i = i - 1)
                num += 1
            end
            if i < global_size[1]
                val[num] = -1
                col[num] = Sten(i = i + 1)
                num += 1
            end
            val[num] = -sum(val)
            col[num] = Sten(i = i)
            PETSc.MatSetValuesStencil!(
                mat,
                row,
                col,
                val,
                PETSc.INSERT_VALUES;
                num_cols = num,
            )
        end

        PETSc.assemblybegin(mat)
        PETSc.assemblyend(mat)

        # Check the the entries in the matrix are correct
        @test mat[1, 1:2] == [1, -1]
        for r in 2:(number_points - 1)
            @test mat[r, (r - 1):(r + 1)] == [-1, 2, -1]
        end
        @test mat[end, (end - 1):end] == [-1, 1]
    end
end

nothing
