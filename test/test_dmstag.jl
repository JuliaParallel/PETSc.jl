using Test
using PETSc, MPI
using ForwardDiff, SparseArrays
MPI.Initialized() || MPI.Init()


@testset "DMStagCreate1d" begin

    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = PETSc.scalartype(petsclib)
        PetscInt    = PETSc.inttype(petsclib)

        # Create 1D DMStag
        dm = PETSc.DMStagCreate1d(
                petsclib,
                comm,
                PETSc.DM_BOUNDARY_PERIODIC,
                20,
                2,
                2,
                PETSc.DMSTAG_STENCIL_BOX,
                2)
        @test PETSc.DMStagGetBoundaryTypes(dm)==PETSc.DM_BOUNDARY_PERIODIC
        PETSc.destroy(dm)

        # Create 1D DMStag with array of local @ of points
        dm = PETSc.DMStagCreate1d(petsclib,comm,PETSc.DM_BOUNDARY_NONE,20,2,2,PETSc.DMSTAG_STENCIL_BOX,2,[20])

        # Test get size
        @test PETSc.DMStagGetGlobalSizes(dm) == 20
        @test PETSc.DMStagGetLocalSizes(dm) == 20

        # Test
        @test PETSc.gettype(dm) == "stag"
        @test PETSc.getdimension(dm) == 1

        # Info about ranks
        @test PETSc.DMStagGetIsFirstRank(dm) == (true,false,false)
        @test PETSc.DMStagGetIsLastRank(dm) == (true,false,false)

        # Boundary
        @test PETSc.DMStagGetBoundaryTypes(dm)==PETSc.DM_BOUNDARY_NONE

        # Corners
        corners         = PETSc.getcorners(dm)
        ghost_corners   = PETSc.getghostcorners(dm)

        @test corners.lower[1] == 1
        @test corners.upper[1] == 20
        @test corners.size[1]  == 20
        @test corners.extra[1] == 1

        @test ghost_corners.lower[1] == 1
        @test ghost_corners.upper[1] == 21
        @test ghost_corners.size[1]  == 21
        @test ghost_corners.extra[1]  == 0

        # DOF
        @test PETSc.DMStagGetDOF(dm) == (2,2)
        PETSc.destroy(dm)

        # Create new struct and pass keyword arguments
        dm_1D = PETSc.DMStagCreate1d(petsclib,comm,PETSc.DM_BOUNDARY_NONE,200,2,2; stag_grid_x=10);
        @test PETSc.DMStagGetGlobalSizes(dm_1D) == 10
        @test PETSc.DMStagGetEntriesPerElement(dm_1D)==4

        # Stencil width & type
        @test  PETSc.DMStagGetStencilWidth(dm_1D)==2
        @test  PETSc.DMStagGetBoundaryTypes(dm_1D) == PETSc.DM_BOUNDARY_NONE

        PETSc.destroy(dm_1D)

        # test ghosted array set using keywords
        dm_ghosted = PETSc.DMStagCreate1d(petsclib,comm,PETSc.DM_BOUNDARY_GHOSTED,200,2,2; stag_grid_x=10);
        @test  PETSc.DMStagGetStencilWidth(dm_ghosted)==2
        corners         = PETSc.getcorners(dm_ghosted)

        @test corners.size[1]==10       # keyword overrides the specified value
        @test  PETSc.DMStagGetBoundaryTypes(dm_ghosted) == PETSc.DM_BOUNDARY_GHOSTED

        ind = PETSc.DMStagGetIndices(dm_ghosted);
        @test ind.center.x[3] == 5

        # simple test to retrieve the KSP object
        # NOTE: need to implement a similar SNES routine
        ksp = PETSc.KSP(dm_ghosted)
        @test PETSc.gettype(ksp)=="gmres"

        PETSc.destroy(dm_ghosted)

        PETSc.finalize(petsclib)
    end
end


@testset "DMStagCreate2d" begin

    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = PETSc.scalartype(petsclib)
        PetscInt    = PETSc.inttype(petsclib)

        # Create 2D DMStag
        dm_2D = PETSc.DMStagCreate2d(petsclib, comm,
                PETSc.DM_BOUNDARY_NONE,
                PETSc.DM_BOUNDARY_NONE,
                20,
                21,
                1,
                1,
                1,
                1,
                1,
                PETSc.DMSTAG_STENCIL_BOX,
                2)

        @test PETSc.DMStagGetGlobalSizes(dm_2D) == (20,21)
        corners = PETSc.getcorners(dm_2D)
        @test corners.size  == [20,21,0]
        @test corners.extra == [1, 1, 0]
        @test corners.lower == [1, 1, 1]
        @test corners.upper == [20,21,0]


        PETSc.finalize(petsclib)
    end
end


@testset "DMStagCreate3d" begin

    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = PETSc.scalartype(petsclib)
        PetscInt    = PETSc.inttype(petsclib)

        # Create 3D DMStag
        dm_3D = PETSc.DMStagCreate3d(petsclib,comm,
                PETSc.DM_BOUNDARY_NONE,
                PETSc.DM_BOUNDARY_NONE,
                PETSc.DM_BOUNDARY_NONE,
                20,
                21,
                22,
                1,
                1,
                1,
                2,
                2,
                2,
                2,
                PETSc.DMSTAG_STENCIL_BOX,
                1,
                [],
                [],
                [])
        @test PETSc.DMStagGetGlobalSizes(dm_3D) == (20, 21, 22)

        # copy struct
        dmnew = PETSc.DMStagCreateCompatibleDMStag(dm_3D,1,1,2,2)
        @test PETSc.DMStagGetGlobalSizes(dmnew) == (20, 21, 22)


        @test PETSc.DMStagGetGlobalSizes(dm_3D) == (20,21,22)
        corners = PETSc.getcorners(dm_3D)
        @test corners.size  == [20,21,22]
        @test corners.extra == [1, 1, 1]
        @test corners.lower == [1, 1, 1]
        @test corners.upper == [20,21,22]


        PETSc.finalize(petsclib)
    end
end



@testset "DMStag Vectors and Coordinates" begin
    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = PETSc.scalartype(petsclib)
        PetscInt    = PETSc.inttype(petsclib)

        # Create 1D DMStag
        dm_1D = PETSc.DMStagCreate1d(petsclib,comm,PETSc.DM_BOUNDARY_NONE,200,2,2; stag_grid_x=10);
        @test PETSc.DMStagGetGlobalSizes(dm_1D) == 10

        # Set coordinates
        PETSc.setuniformcoordinates!(dm_1D, (0,), (10,))

        DMcoord = PETSc.getcoordinateDM(dm_1D)
        @test PETSc.gettype(DMcoord)=="stag"

        # Retrieve array with staggered coordinates
        coord_vec   = PETSc.getcoordinateslocal(dm_1D)
        X_coord     = PETSc.DMStagVecGetArray(DMcoord, coord_vec);
        @test  X_coord[1,2] == 0.5

        # Set coordinates using product (1D) arrays
        dm_1D = PETSc.DMStagCreate1d(petsclib,comm,PETSc.DM_BOUNDARY_NONE,200,2,2; stag_grid_x=10);
        PETSc.setuniformcoordinatesproduct!(dm_1D, (0,), (10,))

        # retrieve DM with coordinates
        DMcoord = PETSc.getcoordinateDM(dm_1D)
        @test PETSc.gettype(DMcoord)=="product"
        coord_vec   = PETSc.getcoordinateslocal(dm_1D)

        # Note: retrieving 1D coordinate vectors using the "product" type appears broken
        #  This is something to be looked at later
        #x_coord,y_coord,z_coord = PETSc.DMStagGetProductCoordinateArraysRead(dm_1D); # BROKEN

        PETSc.DMStagGetLocationSlot(dm_1D, PETSc.DMSTAG_RIGHT, 0) ==4
        @test PETSc.DMStagGetProductCoordinateLocationSlot(dm_1D, PETSc.DMSTAG_RIGHT) == 2

        global_vec      = PETSc.createglobalvector(dm_1D)
        local_vec       = PETSc.createlocalvector(dm_1D)

        # Fill everything with some data
        fill!(local_vec, mpisize)
        fill!(global_vec, mpisize)
        @test global_vec[3] == 1.0

        # Add the local values to the global values
        PETSc.update!(global_vec, local_vec, PETSc.ADD_VALUES)
        @test global_vec[3] == 2.0



        # PETSc.destroy(DMcoord);


        # Do 2D tests
        dm_2D = PETSc.DMStagCreate2d(petsclib, comm,
            PETSc.DM_BOUNDARY_NONE,
            PETSc.DM_BOUNDARY_NONE,
            3,
            4,
            1,
            1,
            1,
            1,
            1,
            PETSc.DMSTAG_STENCIL_BOX,
            2)

        PETSc.setuniformcoordinates!(dm_2D, (1,3), (10,11))
        coord_vec   = PETSc.getcoordinateslocal(dm_2D)

        # Retrieve array with staggered coordinates
        DMcoord_2D = PETSc.getcoordinateDM(dm_2D)
        X_coord_2D = PETSc.DMStagVecGetArray(DMcoord_2D, coord_vec);

        @test X_coord_2D[3,3,1] ≈ 7.0
        @test X_coord_2D[3,4,2] ≈ 9.0
        @test X_coord_2D[3,3,3] ≈ 8.5
        @test X_coord_2D[3,3,4] ≈ 7.0

        vec_test_2D     = PETSc.createlocalvector(dm_2D)
        X               = PETSc.DMStagVecGetArray(dm_2D,vec_test_2D);
        X[end,end,end] = 111;                   # modify 3D array @ some point and DOF
        @test vec_test_2D[end]==111.0           # verify that this modified the vector as well
        Base.finalize(X)                        # release from memory

        #test stencil locations
        pos1 = PETSc.DMStagStencil{PetscInt}(PETSc.DMSTAG_LEFT,1,0,0,1)
        @test pos1.c == 1
        pos2 = PETSc.DMStagStencil{PetscInt}(PETSc.DMSTAG_RIGHT,4,0,0,0)
        pos  = [pos1, pos2]
        @test pos2.loc == PETSc.DMSTAG_RIGHT
        @test pos2.i == 4

        # Retrieve value from stencil
        vec_test       = PETSc.createlocalvector(dm_1D)
        vec_test      .= 1:length(vec_test)                 # point wise copy of data to PetscVec
        val             = PETSc.DMStagVecGetValuesStencil(dm_1D, vec_test, pos1) # this gets a single value
        @test val ==6
        vals            = PETSc.DMStagVecGetValuesStencil(dm_1D, vec_test, 2, pos)     # this gets an array of values
        @test vals[1] == 6

        X_1D            = PETSc.DMStagVecGetArray(dm_1D,vec_test);
        @test X_1D[2,3] == 7.0

        # Set values using stencils
        vec_test_global = PETSc.createglobalvector(dm_1D)
        val1 = PetscScalar.([2222.2, 3.2]);
        PETSc.DMStagVecSetValuesStencil(dm_1D, vec_test_global,    pos1, val1[1],   PETSc.INSERT_VALUES)
        @test vec_test_global[6] ≈ 2222.2
        PETSc.DMStagVecSetValuesStencil(dm_1D, vec_test_global, 2, pos,  val1,      PETSc.INSERT_VALUES)
        @test vec_test_global[21] ≈ 3.2

        pos3 = PETSc.DMStagStencil{PetscInt}(PETSc.DMSTAG_LEFT,1,0,0,1)
        val = PETSc.DMStagVecGetValuesStencil(dm_1D, vec_test, 2, [pos3; pos3])
        @test val[2] == 6.0
        PETSc.destroy(dm_1D);

        PETSc.finalize(petsclib)
    end
end

@testset "DMStag create matrixes" begin
    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = PETSc.scalartype(petsclib)
        PetscInt    = PETSc.inttype(petsclib)

        dm_1D = PETSc.DMStagCreate1d(petsclib,comm,PETSc.DM_BOUNDARY_NONE,200,2,2; stag_grid_x=10);
        PETSc.setuniformcoordinatesproduct!(dm_1D, (0,), (10,))

        A = PETSc.creatematrix(dm_1D);  #
        PETSc.MatSetOption(A, PETSc.MAT_NEW_NONZERO_ALLOCATION_ERR, false)
        @test size(A) == (42,42)
      

        # set some values using normal indices:
        A[1,1]  = 1.0
        A[1,10] = 1.0

        pos1 = PETSc.DMStagStencil{PetscInt}(PETSc.DMSTAG_LEFT,1,0,0,1)
        pos2 = PETSc.DMStagStencil{PetscInt}(PETSc.DMSTAG_RIGHT,4,0,0,0)
        pos  = [pos1, pos2]
        val1 = PetscScalar.([2222.2, 3.2]);
        PETSc.DMStagMatSetValuesStencil(dm_1D, A, pos1, pos1, 11.1, PETSc.INSERT_VALUES)
        PETSc.DMStagMatSetValuesStencil(dm_1D, A, 1, [pos2], 2, pos, val1, PETSc.INSERT_VALUES)

        @test PETSc.assembled(A) == false
        PETSc.assemble(A)
        @test PETSc.assembled(A) == true
        @test A[1,10] == 1.0

        # Reads a value from the matrix, using the stencil structure
        @test PETSc.DMStagMatGetValuesStencil(dm_1D, A, pos1, pos1)== PetscScalar(11.1)
        @test PETSc.DMStagMatGetValuesStencil(dm_1D, A, 1, [pos2], 2, pos)==val1
            
        #PETSc.destroy(A);
        PETSc.destroy(dm_1D);
            
        dofCenter       =   1;
        dofEdge         =   1;
        dofVertex       =   0
        stencilWidth    =   1;
        dm_2D = PETSc.DMStagCreate2d(petsclib,comm,
                                        PETSc.DM_BOUNDARY_GHOSTED,
                                        PETSc.DM_BOUNDARY_GHOSTED,
                                        10,11,
                                        PETSc.PETSC_DECIDE,PETSc.PETSC_DECIDE,
                                        dofVertex,dofEdge,dofCenter,
                                        PETSc.DMSTAG_STENCIL_BOX,stencilWidth)

        vec_test_2D_global      =   PETSc.createglobalvector(dm_2D)
        vec_test_2D_local       =   PETSc.createlocalvector(dm_2D)

        corners                 =   PETSc.getcorners(dm_2D)
        ghost_corners           =   PETSc.getghostcorners(dm_2D)


        for ix=corners.lower[1]:corners.upper[1]
            for iy=corners.lower[2]:corners.upper[2]
                local dof
                # DOF at the center point
                dof     = 0;
                posA    = PETSc.DMStagStencil{PetscInt}(PETSc.DMSTAG_DOWN,ix,iy,0,dof)
                value   = PetscScalar(ix+10);
                PETSc.DMStagVecSetValuesStencil(dm_2D, vec_test_2D_global, posA, value, PETSc.INSERT_VALUES)
                dof     = 0;
                posB    = PETSc.DMStagStencil{PetscInt}(PETSc.DMSTAG_LEFT,ix,iy,0,dof)
                value   = PetscScalar(33);
                PETSc.DMStagVecSetValuesStencil(dm_2D, vec_test_2D_global, posB, value, PETSc.INSERT_VALUES)
                dof     = 0;
                posC    = PETSc.DMStagStencil{PetscInt}(PETSc.DMSTAG_ELEMENT,ix,iy,0,dof)
                value   = PetscScalar(44);
                PETSc.DMStagVecSetValuesStencil(dm_2D, vec_test_2D_global, posC, value, PETSc.INSERT_VALUES)
            end
        end
        PETSc.assemble(vec_test_2D_global) # assemble global vector

        # Add the global values to the local values
        PETSc.update!(vec_test_2D_local, vec_test_2D_global,PETSc.INSERT_VALUES)

        # retrieve value back from the local array and check that it agrees with global one
        dof     = 0;
        pos     = PETSc.DMStagStencil{PetscInt}(PETSc.DMSTAG_DOWN,2,2,0,dof)
        @test PETSc.DMStagVecGetValuesStencil(dm_2D, vec_test_2D_local, pos) == 12.0

        # Extract an array that holds all DOF's
        X2D_dofs  = PETSc.DMStagVecGetArray(dm_2D,vec_test_2D_local)           # extract arrays with all DOF (mostly for visualizing)
        @test X2D_dofs[4,4,1] ≈ PetscScalar(12.0)
        @test X2D_dofs[4,4,2] ≈ PetscScalar(33.0)
        @test X2D_dofs[4,4,3] ≈ PetscScalar(44.0)

        # Extract an array of a specific DOF (here a face velocity @ the left)
        Xarray = PETSc.DMStagGetGhostArrayLocationSlot(dm_2D,vec_test_2D_local, PETSc.DMSTAG_LEFT, 0)
        @test sum(X2D_dofs[:,:,2]-Xarray)==0                # check if the local array is identical to the full array

        Xarray .= 111.                                      # Set a value @ a specific location
        @test vec_test_2D_local[2] ≈ PetscScalar(111)       # verify that this is changed in the PETSc Vec

          
        # cleanup
        #PETSc.destroy(vec_test_2D_global);
        #PETSc.destroy(vec_test_2D_local);
        PETSc.destroy(dm_2D);

           
        PETSc.finalize(petsclib)
    end
end


# -----------------
# Example of DMStag & SNES with AD jacobian
@testset "DMStag: 1D SNES AD" begin

    # Tell AD that it can handle Complex as scalars
    ForwardDiff.can_dual(::Type{ComplexF64}) = true
    ForwardDiff.can_dual(::Type{ComplexF32}) = true

    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)

    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = PETSc.scalartype(petsclib)
        PetscInt    = PETSc.inttype(petsclib)
        if PetscScalar == Float64 || PetscScalar == Float32
            # Define a struct that holds data we need in the local SNES routines below
            mutable struct Data_1{PetscScalar,PetscInt}
                dm
                x_l
                f_l
            end

            user_ctx = Data_1{PetscScalar,PetscInt}(nothing, nothing, nothing);  # holds data we need in the local

            function FormRes!(ptr_fx_g, ptr_x_g, user_ctx)

                # Note that in PETSc, ptr_x_g and ptr_fx_g are pointers to global vectors.
                # Copy global to local vectors that are stored in user_ctx
                PETSc.update!(user_ctx.x_l, ptr_x_g,   PETSc.INSERT_VALUES)
                PETSc.update!(user_ctx.f_l, ptr_fx_g,  PETSc.INSERT_VALUES)

                # Retrieve arrays from the local vectors
                ArrayLocal_x     =   PETSc.DMStagVecGetArrayRead(user_ctx.dm, user_ctx.x_l);  # array with all local x-data
                ArrayLocal_f     =   PETSc.DMStagVecGetArray(user_ctx.dm, user_ctx.f_l);      # array with all local residual

                # Compute local residual
                ComputeLocalResidual(user_ctx.dm, ArrayLocal_x, ArrayLocal_f, user_ctx)

                # Finalize local arrays
                Base.finalize(ArrayLocal_x)
                Base.finalize(ArrayLocal_f)

                # Copy local into global residual vector
                PETSc.update!(ptr_fx_g, user_ctx.f_l,   PETSc.INSERT_VALUES)

            end

            function ComputeLocalResidual(dm, ArrayLocal_x, ArrayLocal_f, user_ctx)
                # Compute the local residual. The vectors include ghost points

                T              =   PETSc.DMStagGetGhostArrayLocationSlot(dm,ArrayLocal_x, PETSc.DMSTAG_LEFT,    0);
                fT             =   PETSc.DMStagGetGhostArrayLocationSlot(dm,ArrayLocal_f, PETSc.DMSTAG_LEFT,    0);

                P              =   PETSc.DMStagGetGhostArrayLocationSlot(dm,ArrayLocal_x, PETSc.DMSTAG_ELEMENT, 0);
                fP             =   PETSc.DMStagGetGhostArrayLocationSlot(dm,ArrayLocal_f, PETSc.DMSTAG_ELEMENT, 0);

                # compute the FD stencil
                indices         =     PETSc.DMStagGetIndices(dm);      # indices of (center/element) points, not including ghost values.
                gc              =     PETSc.getghostcorners(user_ctx.dm);   # start and end of loop including ghost points
                c               =     PETSc.getcorners(user_ctx.dm);        # start and end of loop including ghost points

                nT             =     length(T);                             # array length
                dx             =     1.0/(c.size[1]-1);
                xp             =     (gc.lower[1]:gc.upper[1]).*dx;         # coordinates including ghost points (to define source term)
                F              =     6.0.*xp .+ (xp .+1.e-12).^6.0;         # define source term function

                # Nonlinear equation @ nodal points
                ind            =     indices.vertex.x;                         #  There is one more "vertex" point
                i              =     ind[2:end-1]
                fT[ind[1]]     =     T[ind[1]  ]-0.5;                       # left BC
                fT[ind[end]]   =     T[ind[end]]-2.0;                       # right BC
                fT[i]          =     (T[i .+ 1] - 2*T[i] + T[i .- 1])/dx^2  + T[i].*T[i] - F[i] # NL diffusion with source term

                # second, non-coupled, equation @ center points
                ind            =     indices.center.x;                             #  There is one more "vertex" point
                i              =     ind[2:end-1];
                fP[ind[1]]     =     P[ind[1]]-30.;                             # left BC
                fP[ind[end]]   =     P[ind[end]]-20.;                           # right BC
                fP[i]          =     (P[i .+ 1] - 2*P[i] + P[i .- 1])/dx^2      # steady state diffusion

            end

            function  ForwardDiff_res(x, user_ctx)

                f   = zero(x)               # vector of zeros, of same type as x (local vector)

                ArrayLocal_x     =   PETSc.DMStagVecGetArray(user_ctx.dm, x);        # array with all local x-data
                ArrayLocal_f     =   PETSc.DMStagVecGetArray(user_ctx.dm, f);        # array with all local residual

                ComputeLocalResidual(user_ctx.dm, ArrayLocal_x, ArrayLocal_f, user_ctx);
                # As the residual vector f is linked with ArrayLocal_f, we don't need to pass ArrayLocal_f back

                return f;
            end

            function FormJacobian!(ptr_x_g, J, P, user_ctx)

                # This requires several steps:
                #
                #   1) Extract local vector from global solution (x) vector
                #   2) Compute local jacobian from the residual routine (note that
                #       this routine requires julia vectors as input)

                # Extract the local vector
                PETSc.update!(user_ctx.x_l, ptr_x_g,  PETSc.INSERT_VALUES)
                x               =   PETSc.unsafe_localarray(PetscScalar, user_ctx.x_l.ptr;  write=false, read=true)
                f_Residual      =   (x -> ForwardDiff_res(x, user_ctx));        # pass additional arguments into the routine
                J_julia         =   ForwardDiff.jacobian(f_Residual,x);

                # Note: since x is the LOCAL vector, J_julia also ends up having the same size.
                ind             =   PETSc.LocalInGlobalIndices(user_ctx.dm);
                if PETSc.assembled(J) == false
                    J           =   PETSc.MatSeqAIJ(sparse(J_julia[ind,ind]));
                else
                    J           .=   sparse(J_julia[ind,ind]);
                end

                Base.finalize(ptr_x_g)
                return sparse(J_julia[ind,ind]), ind
            end

            # Main part

            # Construct a 1D test case for a coupled P-T diffusion solver, with 1 DOF @ the center & 1 DOF @ faces
            nx              =   21;
            user_ctx.dm     =   PETSc.DMStagCreate1d(petsclib,comm,
                                    PETSc.DM_BOUNDARY_GHOSTED,
                                    nx,
                                    1,                              # DOF @ vertex
                                    1,                              # DOF @ center
                                    PETSc.DMSTAG_STENCIL_BOX,
                                    1);                             # Stencil width


            x_g             =   PETSc.createglobalvector(user_ctx.dm)
            f_g             =   PETSc.createglobalvector(user_ctx.dm)
            user_ctx.x_l    =   PETSc.createlocalvector(user_ctx.dm)
            user_ctx.f_l    =   PETSc.createlocalvector(user_ctx.dm)

            PJ           =      PETSc.creatematrix(user_ctx.dm);                  # extract (global) matrix from DMStag
            PETSc.MatSetOption(PJ, PETSc.MAT_NEW_NONZERO_ALLOCATION_ERR, false)

            J_julia, ind =      FormJacobian!(x_g, PJ, PJ, user_ctx)
            PJ           =      PETSc.MatSeqAIJ(J_julia)                # assemble non-zero structure
            PETSc.MatSetOption(PJ, PETSc.MAT_NEW_NONZERO_ALLOCATION_ERR, false)

            S = PETSc.SNES{PetscScalar}(petsclib, comm;
                    snes_rtol=1e-12,
                    snes_view=false,
                    snes_monitor=false,
                    ksp_view=false,     # set this to true if you want to see output
                    snes_monitor_true_residual=false,
                    snes_converged_reason=false);
            S.user_ctx  =       user_ctx;

            PETSc.setfunction!(S, FormRes!, f_g)
            PETSc.setjacobian!(S, FormJacobian!, PJ, PJ)

            # Solv
            PETSc.solve!(x_g, S);

            # check
            @test x_g[4] ≈ 29.5
            @test x_g[11] ≈ 0.63797 rtol=1e-4

            # cleanup
            PETSc.destroy(PJ);
            PETSc.destroy(user_ctx.dm);

        end
        PETSc.finalize(petsclib)
    end
end

@testset "DMStag: 2D SNES AD"  begin

    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)

    # Tell AD that it can handle Complex as scalars
    ForwardDiff.can_dual(::Type{ComplexF64}) = true
    ForwardDiff.can_dual(::Type{ComplexF32}) = true

    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = PETSc.scalartype(petsclib)
        PetscInt    = PETSc.inttype(petsclib)
        if PetscScalar == Float64 || PetscScalar == Float32

            mutable struct Data_2D{PetscScalar,PetscInt}
                dm
                x_l
                f_l
            end
            user_ctx = Data_2D{PetscScalar,PetscInt}(nothing, nothing, nothing);  # holds data we need in the local

            function FormRes!(ptr_fx_g, ptr_x_g, user_ctx)
                # Note that in PETSc, ptr_x_g and ptr_fx_g are pointers to global vectors.

                # Copy global to local vectors
                PETSc.update!(user_ctx.x_l, ptr_x_g,  PETSc.INSERT_VALUES)
                PETSc.update!(user_ctx.f_l, ptr_fx_g, PETSc.INSERT_VALUES)

                # Retrieve arrays from the local vectors
                ArrayLocal_x     =   PETSc.DMStagVecGetArrayRead(user_ctx.dm, user_ctx.x_l);  # array with all local x-data
                ArrayLocal_f     =   PETSc.DMStagVecGetArray(user_ctx.dm, user_ctx.f_l);      # array with all local residual

                # Compute local residual
                ComputeLocalResidual(user_ctx.dm, ArrayLocal_x, ArrayLocal_f, user_ctx)

                # Finalize local arrays
                Base.finalize(ArrayLocal_x)
                Base.finalize(ArrayLocal_f)

                # Copy local into global residual vector
                PETSc.update!(ptr_fx_g, user_ctx.f_l, PETSc.INSERT_VALUES)

            end

            function  ForwardDiff_res(x, user_ctx)
                f   = zero(x)               # vector of zeros, of same type as x (local vector)

                ArrayLocal_x     =   PETSc.DMStagVecGetArray(user_ctx.dm, x);        # array with all ocal x-data
                ArrayLocal_f     =   PETSc.DMStagVecGetArray(user_ctx.dm, f);        # array with all ocal residual

                ComputeLocalResidual(user_ctx.dm, ArrayLocal_x, ArrayLocal_f, user_ctx);

                # As the residual vector f is linked with ArrayLocal_f, we don't need to pass ArrayLocal_f back

                return f;
            end

            function ComputeLocalResidual(dm, ArrayLocal_x, ArrayLocal_f, user_ctx)
                # Compute the local residual. The vectors include ghost points

                # Important! Make sure you retrieve the values from the correct locations. In this example we have a
                T              =   PETSc.DMStagGetGhostArrayLocationSlot(dm,ArrayLocal_x, PETSc.DMSTAG_ELEMENT, 0);
                fT             =   PETSc.DMStagGetGhostArrayLocationSlot(dm,ArrayLocal_f, PETSc.DMSTAG_ELEMENT, 0);

                # compute the FD stencil
                indices         =     PETSc.DMStagGetIndices(dm);      # indices of (center/element) points, not including ghost values.

                sz             =     size(user_ctx.dm);                                 # array length
                dx             =     1.0/(sz[1]-1);
                dz             =     1.0/(sz[2]-1);

                # set ghost points for BC'S
                bnd            =    PETSc.DMStagGetBoundaryTypes(user_ctx.dm)
                if bnd[1] == PETSc.DM_BOUNDARY_GHOSTED
                    T[1,:]     =    T[2,:];        # zero flux; dT/dx=0
                    T[end,:]   =    T[end-1,:];    # zero flux
                end

                # Diffusion @ center points
                indx        = indices.center.x;
                indz        = indices.center.y;

                ix          =     indx[1:end]                             # use ghost points in x  (required GHOSTED x-boundary)
                iz          =     indz[2:end-1]                           # center points

                # upper and lower BC (including corners)
                fT[indx,indz[1]]   =    T[indx,indz[1]] .- 0.5;                             # bottom BC
                fT[indx,indz[end]] =    T[indx,indz[end]] .- 2.0;                             # top BC


                fT[ix,iz]       =    (T[ix .+ 1,iz] - 2*T[ix,iz] + T[ix .- 1,iz])/dx^2   +
                                    (T[ix,iz .+ 1] - 2*T[ix,iz] + T[ix,iz .- 1])/dz^2

            end

            function FormJacobian!(ptr_x_g, J, P, user_ctx)
                # This requires several steps:
                #
                #   1) Extract local vector from global solution (x) vector
                #   2) Compute local jacobian from the residual routine (note that
                #       this routine requires julia vectors as input)

                # Extract the local vector
                PETSc.update!(user_ctx.x_l, ptr_x_g,  PETSc.INSERT_VALUES)
                x               =   PETSc.unsafe_localarray(PetscScalar, user_ctx.x_l.ptr;  write=false, read=true)

                f_Residual      =   (x -> ForwardDiff_res(x, user_ctx));        # pass additional rguments into the routine
                J_julia         =   ForwardDiff.jacobian(f_Residual,x);

                # Note: since x is the LOCAL vector, J_julia also ends up having the same size.
                ind             =   PETSc.LocalInGlobalIndices(user_ctx.dm);
                if PETSc.assembled(J) == false
                    J           =   PETSc.MatSeqAIJ(sparse(J_julia[ind,ind]));
                else
                    J           .=   sparse(J_julia[ind,ind]);
                end

                return sparse(J_julia[ind,ind]), ind, sparse(J_julia)
            end

            # Main routine starts here ----

            dofVertex   =   0
            dofEdge     =   0
            dofCenter   =   1
            nx,nz       =   6,25
            user_ctx.dm =   PETSc.DMStagCreate2d(petsclib,comm,
                        PETSc.DM_BOUNDARY_GHOSTED,
                        PETSc.DM_BOUNDARY_NONE,
                        nx,
                        nz,
                        1,
                        1,
                        dofVertex,
                        dofEdge,
                        dofCenter,
                        PETSc.DMSTAG_STENCIL_BOX,
                        1)
            PJ           =      PETSc.creatematrix(user_ctx.dm);                  # extract global) matrix from DMStag
            PETSc.MatSetOption(PJ, PETSc.MAT_NEW_NONZERO_ALLOCATION_ERR, false)

            x_g             =   PETSc.createglobalvector(user_ctx.dm)
            f_g             =   PETSc.createglobalvector(user_ctx.dm)
            user_ctx.x_l    =   PETSc.createlocalvector(user_ctx.dm)
            user_ctx.f_l    =   PETSc.createlocalvector(user_ctx.dm)

            S = PETSc.SNES{PetscScalar}(petsclib, comm;
                    snes_rtol=1e-12,
                    snes_monitor=false, # set to true if you convergence information
                    pc_type="none",
                    snes_monitor_true_residual=true,
                    snes_converged_reason=false);
            S.user_ctx  =       user_ctx;

            J_julia, ind, J_full =      FormJacobian!(x_g, PJ, PJ, user_ctx)
            PJ           =      PETSc.MatSeqAIJ(J_julia)                # assemble non-zero structure
            PETSc.MatSetOption(PJ, PETSc.MAT_NEW_NONZERO_ALLOCATION_ERR, false)

            PETSc.setfunction!(S, FormRes!, f_g)
            PETSc.setjacobian!(S, FormJacobian!, PJ, PJ)

            # Solve 2D system
            sol = PETSc.solve!(x_g, S);

            PETSc.update!(user_ctx.x_l,sol, PETSc.INSERT_VALUES);   # copy global solution -> local vector
            T2d =   PETSc.DMStagGetGhostArrayLocationSlot(user_ctx.dm,user_ctx.x_l, PETSc.DMSTAG_ELEMENT,    0);

            @test T2d[5,5] ≈ 0.75 rtol=1e-3
            #
            # -----------------

            # cleanup
            PETSc.destroy(PJ);
            PETSc.destroy(user_ctx.dm);

        end
        PETSc.finalize(petsclib)

    end
end
