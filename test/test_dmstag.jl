using Test
using PETSc, MPI
#using ForwardDiff
using SparseArrays
MPI.Initialized() || MPI.Init()


#@testset "DMStag All" begin

    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)
    #for petsclib in PETSc.petsclibs
        petsclib = PETSc.petsclibs[1]
        PETSc.initialize(petsclib)
        PetscScalar = PETSc.scalartype(petsclib)
        PetscInt    = PETSc.inttype(petsclib)

        # Create 1D DMStag with new unified constructor
        dm_1D = PETSc.DMStag(petsclib,
                                comm,
                                (PETSc.DM_BOUNDARY_NONE,),
                                (20,),
                                (1,1),
                                2,
                                PETSc.DMSTAG_STENCIL_BOX;
                                points_per_proc=([20],))

        # Create 2D DMStag with new unified constructor
        dm_2D = PETSc.DMStag(petsclib,
                                comm,
                                (PETSc.DM_BOUNDARY_NONE, PETSc.DM_BOUNDARY_NONE),
                                (20, 21),
                                (1,1,1),
                                2,
                                PETSc.DMSTAG_STENCIL_BOX;
                                processors=(1,1),
                                points_per_proc=([20],[21]))

      
        # Create 3D DMStag
        dm_3D = PETSc.DMStag(
            petsclib,
            comm,
            (PETSc.DM_BOUNDARY_NONE,PETSc.DM_BOUNDARY_NONE, PETSc.DM_BOUNDARY_NONE),
            (20,21,22),    # global size
            (1,1,1,2),     # dof_per_node (4 in 3D)
            1,             # stencil_width
            PETSc.DMSTAG_STENCIL_BOX    # stencil type
        )        

        dmTo = PETSc.DMStag(
            petsclib,
            comm,
            (PETSc.DM_BOUNDARY_NONE,PETSc.DM_BOUNDARY_NONE, PETSc.DM_BOUNDARY_NONE),
            (20,21,22),    # global size
            (1,1,2,2),     # dof_per_node (4 in 3D)
            1,             # stencil_width
            PETSc.DMSTAG_STENCIL_BOX    # stencil type
        )   
        @test size(dm_3D) == (20, 21, 22)
        
        # copy struct - using new interface
        dmnew = LibPETSc.DMStagCreateCompatibleDMStag(petsclib,dm_3D,1,1,2,2)
        @test size(dmnew) == (20, 21, 22)


        corners = PETSc.getcorners(dm_3D)
        @test corners.size   == (20,21,22)
        @test corners.nextra == (1, 1, 1)
        @test corners.lower  == CartesianIndex(1, 1, 1)
        @test corners.upper  == CartesianIndex(20,21,22)


        bound = LibPETSc.DMStagGetBoundaryTypes(petsclib, dm_3D) 
        @test bound == (PETSc.LibPETSc.DM_BOUNDARY_NONE, PETSc.LibPETSc.DM_BOUNDARY_NONE, PETSc.LibPETSc.DM_BOUNDARY_NONE)

        corners = LibPETSc.DMStagGetCorners(petsclib, dm_3D) 
        @test corners == (0, 0, 0, 20, 21, 22, 1, 1, 1)

        dof = LibPETSc.DMStagGetDOF(petsclib, dm_3D)
        @test dof == (1, 1, 1, 2)

        entries = LibPETSc.DMStagGetEntries(petsclib, dm_3D)
        @test entries == 88575

        entries = LibPETSc.DMStagGetEntriesLocal(petsclib, dm_3D)
        @test entries == 95634

        out = LibPETSc.DMStagGetEntriesPerElement(petsclib, dm_3D)
        @test out == 9

        dm = dm_3D
        out = LibPETSc.DMStagGetGhostCorners(petsclib, dm)
        @test out == (0, 0, 0, 21, 22, 23)

        out = LibPETSc.DMStagGetGlobalSizes(petsclib, dm)
        @test out == (20, 21, 22)

        #out = LibPETSc.DMStagGetIsFirstRank(petsclib, dm)
        #@test out == (false, false, true)

        #out = LibPETSc.DMStagGetIsLastRank(petsclib, dm)
        #@test out == (false, false, true)

        out = LibPETSc.DMStagGetLocalSizes(petsclib, dm)
        @test out == (20, 21, 22)

        loc = PETSc.LibPETSc.DMSTAG_DOWN
        out = LibPETSc.DMStagGetLocationDOF(petsclib, dm, loc)
        @test out == 1

        loc = PETSc.LibPETSc.DMSTAG_DOWN
        c = 1
        out = LibPETSc.DMStagGetLocationSlot(petsclib, dm, loc, 0)
        @test out == 5

        out = LibPETSc.DMStagGetNumRanks(petsclib, dm)
        @test out == (1, 1, 1)

        out = LibPETSc.DMStagGetStencilType(petsclib, dm)
        @test out == LibPETSc.DMSTAG_STENCIL_BOX

        out = LibPETSc.DMStagGetStencilWidth(petsclib, dm)
        @test out == 1

        out = LibPETSc.DMStagGetRefinementFactor(petsclib, dm)
        @test out == (2, 2, 2)

        out = LibPETSc.DMStagPopulateLocalToGlobalInjective(petsclib, dm)
        @test out == nothing

        v = PETSc.DMLocalVec(dm);
        v[10] = 10.1;
        PETSc.assemble!(v);
        
        #X,X_ptr = LibPETSc.DMStagVecGetArray(petsclib, dm,v); # doesn't crash but gives wrong result
        #@test X[10] == PetscScalar(10.1)        

        dm_new = LibPETSc.DMStagCreateCompatibleDMStag(petsclib, dm, 1,2,3,4)
        @test LibPETSc.DMStagGetDOF(petsclib, dm_new) == (1,2,3,4)

        # this needs fixing!
        #out = LibPETSc.DMStagGetOwnershipRanges(petsclib, dm)
        #@test out == ([20], [21], [22])
        
        # Not sure we should test that here as it should be called before DMSetup
        #b0 = PETSc.LibPETSc.DM_BOUNDARY_NONE
        #b1 = PETSc.LibPETSc.DM_BOUNDARY_GHOSTED
        #b2 = PETSc.LibPETSc.DM_BOUNDARY_MIRROR
        #out = LibPETSc.DMStagSetBoundaryTypes(petsclib, dm, b0, b1, b2)
        #@test out == 
        
        dwn = PETSc.LibPETSc.DMSTAG_DOWN
        loc  = PETSc.LibPETSc.DMStagStencil(dwn,1,1,1,0)
        loc2 = PETSc.LibPETSc.DMStagStencil(dwn,1,1,2,0)
        
        vg = PETSc.DMGlobalVec(dm);
        LibPETSc.DMStagVecSetValuesStencil(petsclib, dm,vg,1,[loc],[1.0],PETSc.LibPETSc.INSERT_VALUES) 
        @test vg[4145] == 1.0

        LibPETSc.DMStagVecSetValuesStencil(petsclib, dm,vg,2, [loc, loc2],[1.0,2.0],PETSc.LibPETSc.INSERT_VALUES) 
        @test extrema(vg) == (0.0,2.0)

        out = LibPETSc.DMStagVecGetValuesStencil(petsclib, dm,v, 1, [loc])
        @test out == PetscScalar[0.0]
        
        # results in a segfault for some libs
        if isa(PetscScalar, AbstractFloat)
            out = LibPETSc.DMStagVecGetValuesStencil(petsclib, dm,v, 2,[loc, loc2])
            @test out == [0.0, 0.0]

            # results in segfault for some libs
            out = LibPETSc.DMStagStencilToIndexLocal(petsclib, dm,3,2,[loc, loc2])
            ##@test out == PetscInt[0, 4361] #wrong??
            @test out ==   PetscInt[4361, 8519]   
        end

        # NOT YET WORKING
        #pda = LibPETSc.DMDAFromDMStag(petsclib,dm)
        #pda, pdavec = LibPETSc.DMStagVecSplitToDMDA(petsclib, dm,v,dwn,0)
        #@test length(pdavec) == 9680


        vecTo = PETSc.DMGlobalVec(dmTo);
        out = LibPETSc.DMStagMigrateVec(petsclib, dm,vg,dmTo,vecTo)
        @test isnothing(out)

        v1D = PETSc.DMLocalVec(dm_1D);

        # Set values using the 2D array interface
        # DMStagVecGetArray does this internally
        x,_,_,m,_,_ = LibPETSc.DMStagGetGhostCorners(petsclib, dm_1D)
        entriesEl = LibPETSc.DMStagGetEntriesPerElement(petsclib,dm_1D)

        array2D = LibPETSc.VecGetArray2d(petsclib, v1D, m, entriesEl, 0, 0 )
        array2D .= 1.0
        array2D[2,2] = 2.2
        LibPETSc.VecRestoreArray2d(petsclib, v1D, m, entriesEl, 0, 0, array2D)
        @test extrema(v1D) == (1.0, 2.2)

        # Now lets use DMStagVecGetArray. Note that we had to do manual modifcations to that routine
        array2D_1 = LibPETSc.DMStagVecGetArray(petsclib, dm_1D, v1D) 
        array2D_1 .= 2.0
        array2D_1[2,2] = 4.2
        LibPETSc.DMStagVecRestoreArray(petsclib, dm_1D, v1D, array2D_1) 
        @test extrema(v1D) == (2.0, 4.2)


#=

        # Test product coordinates and friends 
        # Converted from old DMStagCreate3d interface to new unified DMStag constructor
        dm_3D_pd = PETSc.DMStag(
            petsclib,
            comm,
            (PETSc.DM_BOUNDARY_NONE, PETSc.DM_BOUNDARY_NONE, PETSc.DM_BOUNDARY_NONE),
            (20, 21, 22),               # global sizes
            (1, 1, 1, 1),               # dof_per_node (vertex, edge, face, element)
            2,                          # stencil width
            PETSc.DMSTAG_STENCIL_BOX;   # stencil type
            processors = (1, 1, 1),
            points_per_proc = ([20], [21], [22])
        )

        @test LibPETSc.DMStagGetGlobalSizes(petsclib, dm_3D_pd) == (20, 21, 22)
         
        LibPETSc.DMStagSetUniformCoordinatesProduct(petsclib, dm_3D_pd, 0.0, 1.0, 0.0, 2.0, 0.0, 3.0)

        # Retrieve 1D coordinate arrays
        x,y,z,px,py,pz = LibPETSc.DMStagGetProductCoordinateArrays(petsclib, dm_3D_pd)
        dims = LibPETSc.DMStagGetGhostCorners(petsclib, dm_3D_pd)[4:6] 
        @test length(x) == dims[1]
        @test length(z) == dims[3]
        @test x[10] ≈ PetscScalar(0.225)
        @test y[10] == PetscScalar(0.42857142857142855)
        @test z[10] == PetscScalar(0.6136363636363635)

        x[10] = 0.230
        
        LibPETSc.DMStagRestoreProductCoordinateArrays(petsclib, dm_3D_pd, px,py,pz)
        

        x,y,z,px,py,pz = LibPETSc.DMStagGetProductCoordinateArraysRead(petsclib, dm_3D_pd)
        LibPETSc.DMStagRestoreProductCoordinateArraysRead(petsclib, dm_3D_pd, px,py,pz)



        slot = LibPETSc.DMStagGetProductCoordinateLocationSlot(petsclib, dm_3D_pd, PETSc.LibPETSc.DMSTAG_ELEMENT)
        @test slot==1
        slot = LibPETSc.DMStagGetProductCoordinateLocationSlot(petsclib, dm_3D_pd, PETSc.LibPETSc.DMSTAG_LEFT)
        @test slot==0

        # Converted from old DMStagCreate2d interface (fine grid) to new DMStag constructor
        dmf = PETSc.DMStag(
            petsclib,
            comm,
            (PETSc.DM_BOUNDARY_NONE, PETSc.DM_BOUNDARY_NONE),
            (32, 32),            # global sizes
            (1, 0, 0),           # dof_per_node (vertex, edge, element)
            2,                   # stencil width
            PETSc.DMSTAG_STENCIL_BOX;
            processors = (1, 1),
            points_per_proc = ([32], [32])
        )

        # Converted from old DMStagCreate2d interface (coarse grid) to new DMStag constructor
        dmc = PETSc.DMStag(
            petsclib,
            comm,
            (PETSc.DM_BOUNDARY_NONE, PETSc.DM_BOUNDARY_NONE),
            (16, 16),            # global sizes
            (1, 0, 0),           # dof_per_node (vertex, edge, element)
            2,                   # stencil width
            PETSc.DMSTAG_STENCIL_BOX;
            processors = (1, 1),
            points_per_proc = ([16], [16])
        )    

        xf = PETSc.DMLocalVec(dmf)      
        xf .= 1.0
        
        
        xc = LibPETSc.DMStagRestrictSimple(petsclib, dmf,xf,dmc)
        @test xc[5] == PetscScalar(1.0)


        # To be added as test: get matrix first
        #out = LibPETSc.DMStagMatGetValuesStencil(petsclib, dm,mat,nRow,posRow,nCol,posCol)
        #@test out == 

        # To be added as test: get matrix first
        #out = LibPETSc.DMStagMatSetValuesStencil(petsclib, dm,mat,nRow,posRow,nCol,posCol,val,insertMode)
        #@test out == 

        #out = LibPETSc.DMStagSetCoordinateDMType(petsclib, dm)
        #@test out == 

        
        #out = LibPETSc.DMStagSetDOF(petsclib, dm,1,2,1,1) # needs to be called before DMSetup
        #@test out == 

        #out = LibPETSc.DMStagSetGlobalSizes(petsclib, dm,10,11,12) # needs to be called before DMSetup
        #@test out == 

        #out = LibPETSc.DMStagSetNumRanks(petsclib, dm,1,1,1) # needs to be called before DMSetup
        #@test out == 

        #out = LibPETSc.DMStagSetOwnershipRanges(petsclib, dm,lx,Int64,Int64)
        #@test out == 

        #out = LibPETSc.DMStagSetStencilType(petsclib, dm,stencilType)
        #@test out == 
#
        #out = LibPETSc.DMStagSetStencilWidth(petsclib, dm,stencilWidth)
        #@test out == 
#
        #out = LibPETSc.DMStagSetRefinementFactor(petsclib, dm,refine_x,refine_y,refine_z)
        #@test out == 
#
        #out = LibPETSc.DMStagSetUniformCoordinates(petsclib, dm,0.0,1.0,1.1,2.1,3.1,4.1)
        #@test out == 

        out = LibPETSc.DMStagSetUniformCoordinatesExplicit(petsclib, dm,0.1,1.1,1.3,3.1,2.1,2.8)
        @test isnothing(out)

    
        xmin,xmax,ymin,ymax,zmin,zmax = 0.1,1.1,1.3,3.1,2.1,2.8
        out = LibPETSc.DMStagSetUniformCoordinatesProduct(petsclib, dmc,xmin,xmax,ymin,ymax,zmin,zmax)
        @test isnothing(out)

        PETSc.finalize(petsclib)
    #end
#end
=#

#=
@testset "DMStagCreate1d" begin

    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = PETSc.scalartype(petsclib)
        PetscInt    = PETSc.inttype(petsclib)

        # Create 1D DMStag
        dm = PETSc.DMStag(
                petsclib,
                comm,
                (PETSc.DM_BOUNDARY_PERIODIC,),
                (20,),
                (2,0),
                2,
                PETSc.DMSTAG_STENCIL_BOX)

        @test LibPETSc.DMStagGetBoundaryTypes(petsclib, dm)== (PETSc.LibPETSc.DM_BOUNDARY_PERIODIC, PETSc.LibPETSc.DM_BOUNDARY_NONE, PETSc.LibPETSc.DM_BOUNDARY_NONE)
        PETSc.destroy(dm)

        # Create 1D DMStag with array of local @ of points
        dm = PETSc.DMStag(
            petsclib,
            comm,
            (PETSc.DM_BOUNDARY_NONE,),
            (20,),
            (2,2),
            2;
            points_per_proc = ([20],),
        )

        # Test get size
        @test LibPETSc.DMStagGetGlobalSizes(petsclib, dm) == (20,0,0)
        @test LibPETSc.DMStagGetLocalSizes(petsclib, dm) == (20,0,0)

        # Test
       # @test PETSc.gettype(dm) == "stag"
        @test PETSc.getdimension(dm) == 1

        # Info about ranks  
        @test LibPETSc.DMStagGetIsFirstRank(petsclib, dm) == (true,false,false)
        @test LibPETSc.DMStagGetIsFirstRank(petsclib, dm) == (true,false,false)

        # Boundary
        @test LibPETSc.DMStagGetBoundaryTypes(petsclib, dm)[1]==PETSc.LibPETSc.DM_BOUNDARY_NONE

        # Corners
        corners         = PETSc.getcorners(dm)
        ghost_corners   = PETSc.getghostcorners(dm)

        @test corners.lower[1] == 1
        @test corners.upper[1] == 20
        @test corners.size[1]  == 20
        @test corners.nextra[1] == 1

        @test ghost_corners.lower[1] == 1
        @test ghost_corners.upper[1] == 21
        @test ghost_corners.size[1]  == 21

        # DOF
        @test LibPETSc.DMStagGetDOF(petsclib, dm) == (2,2,0,0)
        PETSc.destroy(dm)

        # Create new struct and pass keyword arguments
        dm_1D = PETSc.DMStag(petsclib,comm,(PETSc.DM_BOUNDARY_NONE,),(200,),(2,2),2; stag_grid_x=10);
        @test PETSc.globalsize(dm_1D)[1] == 10
        @test PETSc.getentriesperelement(dm_1D)==4

        # Stencil width & type
        @test  LibPETSc.DMStagGetStencilWidth(petsclib, dm_1D)==2
        @test  LibPETSc.DMStagGetBoundaryTypes(petsclib, dm_1D)[1] == PETSc.DM_BOUNDARY_NONE

        PETSc.destroy(dm_1D)

        # test ghosted array set using keywords
        dm_ghosted = PETSc.DMStag(
            petsclib,
            comm,
            (PETSc.DM_BOUNDARY_GHOSTED,),
            (200,),
            (2,2),
            2;
            stag_grid_x = 10,
        )


        @test  LibPETSc.DMStagGetStencilWidth(petsclib, dm_ghosted)==2
        corners         = PETSc.getcorners(dm_ghosted)

        @test corners.size[1]==10       # keyword overrides the specified value
        @test  LibPETSc.DMStagGetBoundaryTypes(petsclib, dm_ghosted)[1] == PETSc.DM_BOUNDARY_GHOSTED

        ind = LibPETSc.DMStagGetIndices(petsclib, dm_ghosted);
        @test ind.center.x[3] == 5

        # simple test to retrieve the KSP object
        # NOTE: need to implement a similar SNES routine
        ksp = PETSc.KSP(dm_ghosted)
        @test PETSc.KSPGetType(ksp)=="gmres"

        PETSc.destroy(dm_ghosted)

        PETSc.finalize(petsclib)
    end
end


@testset "DMStagCreate2d" begin

    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)
    for petsclib in PETSc.petsclibs
        petsclib =  PETSc.petsclibs[1]
        
        PETSc.initialize(petsclib)
        PetscScalar = PETSc.scalartype(petsclib)
        PetscInt    = PETSc.inttype(petsclib)

        # Create 2D DMStag
        dm_2D = PETSc.DMStag(
            petsclib,
            comm,
            (PETSc.DM_BOUNDARY_NONE,PETSc.DM_BOUNDARY_NONE),
            (20,21),    # global size
            (1,1,1),    # dof_per_node (3 in 2D)
            1,          # stencil_width
            PETSc.DMSTAG_STENCIL_BOX    # stencil type
        )

        @test LibPETSc.DMStagGetGlobalSizes(petsclib, dm_2D) == (20,21,0)
        corners = PETSc.getcorners(dm_2D)
        @test corners.size   == (20,21,0)
        @test corners.nextra == (1, 1, 0)
        @test corners.lower  == CartesianIndex(1, 1, 1)
        @test corners.upper  == CartesianIndex(20,21,0)

        PETSc.finalize(petsclib)
    end
end


#=

#@testset "DMStag Vectors and Coordinates" begin
    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)
    #for petsclib in PETSc.petsclibs
        petsclib = PETSc.petsclibs[1]
    
        PETSc.initialize(petsclib)
        PetscScalar = PETSc.scalartype(petsclib)
        PetscInt    = PETSc.inttype(petsclib)

        # Create 1D DMStag
        dm_1D = PETSc.DMStag(
            petsclib,
            comm,
            (PETSc.DM_BOUNDARY_NONE,),
            (200,),
            (2,2),
            1,
            PETSc.DMSTAG_STENCIL_BOX,
            stag_grid_x=10)

        @test LibPETSc.DMStagGetGlobalSizes(petsclib, dm_1D) == (10,0,0)

        # Set coordinates using product (1D) arrays
        PETSc.setuniformcoordinates!(dm_1D, (0,), (10,))

        #DMcoord = PETSc.getcoordinateDM(dm_1D)
        @test PETSc.gettype(DMcoord)=="product"

        # Retrieve array with staggered coordinates
        X_coord = PETSc.getcoordinatearray(dm_1D)[1]
        @test  X_coord[2,1] == 0.5

        LibPETSc.DMStagGetLocationSlot(petsclib, dm_1D, LibPETSc.DMSTAG_RIGHT, 0) ==4
     
        global_vec      = PETSc.DMLocalVec(dm_1D)
        local_vec       = PETSc.DMGlobalVec(dm_1D)

        # Fill everything with some data
        fill!(local_vec, mpisize)
        fill!(global_vec, mpisize)
        @test global_vec[3] == 1.0

        # Add the local values to the global values
        PETSc.update!(global_vec, local_vec, PETSc.ADD_VALUES)
        @test global_vec[3] == 2.0

        # Do 2D tests
        dm_2D = PETSc.DMStag(
            petsclib,
            comm,
            (PETSc.DM_BOUNDARY_NONE,PETSc.DM_BOUNDARY_NONE),
            (3,4),      # global size
            (1,1,1),    # dof_per_node (3 in 2D)
            1,          # stencil_width
            PETSc.DMSTAG_STENCIL_BOX    # stencil type
        )

        PETSc.setuniformcoordinates!(dm_2D, (1,3), (10,11))
       
        
        # Retrieve array with staggered coordinates
        X_coord_2D = PETSc.getcoordinatearray(dm_2D)
        @test X_coord_2D[1][1,3] ≈ 7.0
        @test X_coord_2D[2][1,4] ≈ 9.0
        @test X_coord_2D[1][2,3] ≈ 8.5
        @test X_coord_2D[1][1,3] ≈ 7.0

        vec_test_2D     = PETSc.DMLocalVec(dm_2D)
        X               = LibPETSc.DMStagVecGetArray(petsclib, dm_2D,vec_test_2D);
        X[end,end,end] = 111;                   # modify 3D array @ some point and DOF
        @test vec_test_2D[end]==111.0           # verify that this modified the vector as well
        Base.finalize(X)                        # release from memory

#### finished until here 
#=
        #test stencil locations
        pos1 = PETSc.DMStagStencil{PetscInt}(PETSc.DMSTAG_LEFT,1,0,0,1)
        @test pos1.c == 1
        pos2 = PETSc.DMStagStencil{PetscInt}(PETSc.DMSTAG_RIGHT,4,0,0,0)
        pos  = [pos1, pos2]
        @test pos2.loc == PETSc.DMSTAG_RIGHT
        @test pos2.i == 4

        # Retrieve value from stencil
        vec_test       = PETSc.DMLocalVec(dm_1D)
        vec_test      .= 1:length(vec_test)                 # point wise copy of data to PetscVec
        val             = LibPETSc.DMStagVecGetValuesStencil(petsclib, dm_1D, vec_test, pos1) # this gets a single value
        @test val ==6
        vals            = LibPETSc.DMStagVecGetValuesStencil(petsclib, dm_1D, vec_test, 2, pos)     # this gets an array of values
        @test vals[1] == 6

        X_1D            = LibPETSc.DMStagVecGetArray(petsclib, dm_1D,vec_test);
        @test X_1D[2,3] == 7.0

        # Set values using stencils
        vec_test_global = PETSc.DMGlobalVec(dm_1D)
        val1 = PetscScalar.([2222.2, 3.2]);
        LibPETSc.DMStagVecSetValuesStencil(petsclib, dm_1D, vec_test_global,    pos1, val1[1],   PETSc.INSERT_VALUES)
        @test vec_test_global[6] ≈ 2222.2
        LibPETSc.DMStagVecSetValuesStencil(petsclib, dm_1D, vec_test_global, 2, pos,  val1,      PETSc.INSERT_VALUES)
        @test vec_test_global[21] ≈ 3.2

        pos3 = PETSc.DMStagStencil{PetscInt}(PETSc.DMSTAG_LEFT,1,0,0,1)
        val = LibPETSc.DMStagVecGetValuesStencil(petsclib, dm_1D, vec_test, 2, [pos3; pos3])
        @test val[2] == 6.0
        #PETSc.destroy(dm_1D);
=#
        #PETSc.finalize(petsclib)
    #end
#end


# FIXME: part below that is commented segfaults on linux
@testset "DMStag create matrixes" begin
    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)
    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = PETSc.scalartype(petsclib)
        PetscInt    = PETSc.inttype(petsclib)

    # Converted from old DMStagCreate1d to new DMStag constructor
    dm_1D = PETSc.DMStag(petsclib, comm, (PETSc.DM_BOUNDARY_NONE,), (200,), (2, 2), 2;
                  stag_grid_x = 10)
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
        LibPETSc.DMStagMatSetValuesStencil(petsclib, dm_1D, A, pos1, pos1, 11.1, PETSc.INSERT_VALUES)
        LibPETSc.DMStagMatSetValuesStencil(petsclib, dm_1D, A, 1, [pos2], 2, pos, val1, PETSc.INSERT_VALUES)

        @test PETSc.assembled(A) == false
        PETSc.assemble(A)
        @test PETSc.assembled(A) == true
        @test A[1,10] == 1.0

        # Reads a value from the matrix, using the stencil structure
        @test LibPETSc.DMStagMatGetValuesStencil(petsclib, dm_1D, A, pos1, pos1)== PetscScalar(11.1)
        @test LibPETSc.DMStagMatGetValuesStencil(petsclib, dm_1D, A, 1, [pos2], 2, pos)==val1
            
        #PETSc.destroy(A);
        PETSc.destroy(dm_1D);
            
        dofCenter       =   1;
        dofEdge         =   1;
        dofVertex       =   0
        stencilWidth    =   1;
        # Converted from old DMStagCreate2d (matrix test) to new DMStag constructor
        dm_2D = PETSc.DMStag(
            petsclib,
            comm,
            (PETSc.DM_BOUNDARY_GHOSTED, PETSc.DM_BOUNDARY_GHOSTED),
            (10, 11),                        # global sizes
            (dofVertex, dofEdge, dofCenter), # dof_per_node
            stencilWidth,
            PETSc.DMSTAG_STENCIL_BOX;
            processors = (PETSc.PETSC_DECIDE, PETSc.PETSC_DECIDE)
        )

        vec_test_2D_global      =   PETSc.createglobalvector(dm_2D)
        vec_test_2D_local       =   PETSc.createlocalvector(dm_2D)

        corners                 =   PETSc.getcorners(dm_2D)
        ghost_corners           =   PETSc.getghostcorners(dm_2D)

        # ----
        # FIXME: 
        # the commented lines below result in a segfault on linux
        # To be checked whether this is still the case for the auto-wrapped library
        for ix=corners.lower[1]:corners.upper[1]
            for iy=corners.lower[2]:corners.upper[2]
                local dof
                # DOF at the center point
                dof     = 0;
                posA    = PETSc.DMStagStencil{PetscInt}(PETSc.DMSTAG_DOWN,ix,iy,0,dof)
                value   = PetscScalar(ix+10);
                #PETSc.DMStagVecSetValuesStencil(dm_2D, vec_test_2D_global, 1, [posA], [value], PETSc.INSERT_VALUES)
                dof     = 0;
                posB    = PETSc.DMStagStencil{PetscInt}(PETSc.DMSTAG_LEFT,ix,iy,0,dof)
                value   = PetscScalar(33);
                #PETSc.DMStagVecSetValuesStencil(dm_2D, vec_test_2D_global, posB, value, PETSc.INSERT_VALUES)
                dof     = 0;
                posC    = PETSc.DMStagStencil{PetscInt}(PETSc.DMSTAG_ELEMENT,ix,iy,0,dof)
                value   = PetscScalar(44);
                #PETSc.DMStagVecSetValuesStencil(dm_2D, vec_test_2D_global, posC, value, PETSc.INSERT_VALUES)
            end
        end
        
        PETSc.assemble(vec_test_2D_global) # assemble global vector

        # Add the global values to the local values
        PETSc.update!(vec_test_2D_local, vec_test_2D_global,PETSc.INSERT_VALUES)

        
        # retrieve value back from the local array and check that it agrees with global one
        dof     = 0;
        pos     = PETSc.DMStagStencil{PetscInt}(PETSc.DMSTAG_DOWN,2,2,0,dof)
#        @test LibPETSc.DMStagVecGetValuesStencil(petsclib, dm_2D, vec_test_2D_local, pos) == 12.0

#        # Extract an array that holds all DOF's
        X2D_dofs  = LibPETSc.DMStagVecGetArray(petsclib, dm_2D,vec_test_2D_local)           # extract arrays with all DOF (mostly for visualizing)
#        @test X2D_dofs[4,4,1] ≈ PetscScalar(12.0)
#        @test X2D_dofs[4,4,2] ≈ PetscScalar(33.0)
#        @test X2D_dofs[4,4,3] ≈ PetscScalar(44.0)

        # ----

        # Extract an array of a specific DOF (here a face velocity @ the left)
        Xarray = LibPETSc.DMStagGetGhostArrayLocationSlot(petsclib, dm_2D,vec_test_2D_local, PETSc.DMSTAG_LEFT, 0)
        @test sum(X2D_dofs[:,:,2]-Xarray)==0                # check if the local array is identical to the full array

        Xarray .= 111.                                      # Set a value @ a specific location
        @test vec_test_2D_local[2] ≈ PetscScalar(111)       # verify that this is changed in the PETSc Vec

          
        # cleanup
        PETSc.destroy(dm_2D);
        
           
        PETSc.finalize(petsclib)
    end
end
=#

#=
# -----------------
# Example of DMStag & SNES with AD jacobian
#@testset "DMStag: 1D SNES AD" begin

    # Tell AD that it can handle Complex as scalars
    ForwardDiff.can_dual(::Type{ComplexF64}) = true
    ForwardDiff.can_dual(::Type{ComplexF32}) = true

    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)

    mutable struct Data_1{PetscScalar,PetscInt}
        dm
        x_l
        f_l
    end

    #for petsclib in PETSc.petsclibs
    petsclib = PETSc.petsclibs[1]
        PETSc.initialize(petsclib)
        PetscScalar = PETSc.scalartype(petsclib)
        PetscInt    = PETSc.inttype(petsclib)
        ##if PetscScalar == Float64 || PetscScalar == Float32
            # Define a struct that holds data we need in the local SNES routines below
         
            user_ctx = Data_1{PetscScalar,PetscInt}(nothing, nothing, nothing);  # holds data we need in the local

            function FormRes!(ptr_fx_g, ptr_x_g, user_ctx)

                # Note that in PETSc, ptr_x_g and ptr_fx_g are pointers to global vectors.
                # Copy global to local vectors that are stored in user_ctx
                PETSc.update!(user_ctx.x_l, ptr_x_g,   PETSc.INSERT_VALUES)
                PETSc.update!(user_ctx.f_l, ptr_fx_g,  PETSc.INSERT_VALUES)

                # Retrieve arrays from the local vectors
                ArrayLocal_x     =   LibPETSc.DMStagVecGetArrayRead(petsclib, user_ctx.dm, user_ctx.x_l);  # array with all local x-data
                ArrayLocal_f     =   LibPETSc.DMStagVecGetArray(petsclib, user_ctx.dm, user_ctx.f_l);      # array with all local residual

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

                T              =   LibPETSc.DMStagGetGhostArrayLocationSlot(petsclib, dm,ArrayLocal_x, PETSc.DMSTAG_LEFT,    0);
                fT             =   LibPETSc.DMStagGetGhostArrayLocationSlot(petsclib, dm,ArrayLocal_f, PETSc.DMSTAG_LEFT,    0);

                P              =   LibPETSc.DMStagGetGhostArrayLocationSlot(petsclib, dm,ArrayLocal_x, PETSc.DMSTAG_ELEMENT, 0);
                fP             =   LibPETSc.DMStagGetGhostArrayLocationSlot(petsclib, dm,ArrayLocal_f, PETSc.DMSTAG_ELEMENT, 0);

                # compute the FD stencil
                indices         =     LibPETSc.DMStagGetIndices(petsclib, dm);      # indices of (center/element) points, not including ghost values.
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
            user_ctx.dm     =   PETSc.DMStag(petsclib,comm,
                                    (PETSc.DM_BOUNDARY_GHOSTED,),
                                    (nx,),
                                    (1,1),                              # (DOF @ vertex, DOF @ center)
                                    PETSc.DMSTAG_STENCIL_BOX,
                                    1; stag_grid_x=[nx]);                             # Stencil width


            x_g             =   PETSc.DMGlobalVec(user_ctx.dm)
            f_g             =   PETSc.DMGlobalVec(user_ctx.dm)
            user_ctx.x_l    =   PETSc.DMLocalVec(user_ctx.dm)
            user_ctx.f_l    =   PETSc.DMLocalVec(user_ctx.dm)

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

        ##end
        #PETSc.finalize(petsclib)
    #end
#end

@testset "DMStag: 2D SNES AD"  begin

    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)

    # Tell AD that it can handle Complex as scalars
#    ForwardDiff.can_dual(::Type{ComplexF64}) = true
#    ForwardDiff.can_dual(::Type{ComplexF32}) = true
    mutable struct Data_2D{PetscScalar,PetscInt}
        dm
        x_l
        f_l
    end

    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = PETSc.scalartype(petsclib)
        PetscInt    = PETSc.inttype(petsclib)
        if PetscScalar == Float64 || PetscScalar == Float32

          
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
            user_ctx.dm =   PETSc.DMStag(petsclib,comm,
                        (PETSc.DM_BOUNDARY_GHOSTED, PETSc.DM_BOUNDARY_NONE),
                        (nx, nz),
                        (1, 1),
                        (dofVertex, dofEdge, dofCenter),
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

=#


=#