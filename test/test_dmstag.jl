using Test
using PETSc, MPI
#using SparseArrays
MPI.Initialized() || MPI.Init()

@testset "DMStag All" begin

    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)
    for petsclib in PETSc.petsclibs[1:4]
        #petsclib = PETSc.petsclibs[1]
        PETSc.initialize(petsclib)
        PetscScalar = PETSc.scalartype(petsclib)
        PetscInt    = PETSc.inttype(petsclib)
        PetscReal   = real(PetscScalar)
        # Create 1D DMStag with new unified constructor
        dm_1D = PETSc.DMStag(petsclib,
                                comm,
                                (PETSc.DM_BOUNDARY_NONE,),
                                (20,),
                                (1,1),
                                2,
                                PETSc.DMSTAG_STENCIL_BOX;
                                points_per_proc=([PetscInt(20)],))

        # Create 2D DMStag with new unified constructor
        dm_2D = PETSc.DMStag(petsclib,
                                comm,
                                (PETSc.DM_BOUNDARY_NONE, PETSc.DM_BOUNDARY_NONE),
                                (20, 21),
                                (1,1,1),
                                2,
                                PETSc.DMSTAG_STENCIL_BOX;
                                processors=(1,1),
                                points_per_proc=([PetscInt(20)],[PetscInt(21)]))

      
        # Create 3D DMStag
        dm_3D = PETSc.DMStag(
            petsclib,
            comm,
            (PETSc.DM_BOUNDARY_NONE,PETSc.DM_BOUNDARY_NONE, PETSc.DM_BOUNDARY_NONE),
            (PetscInt(20),PetscInt(21),PetscInt(22)),    # global size
            (PetscInt(1),PetscInt(1),PetscInt(1),PetscInt(2)),     # dof_per_node (4 in 3D)
            PetscInt(1),             # stencil_width
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
        dmnew = LibPETSc.DMStagCreateCompatibleDMStag(petsclib,dm_3D,PetscInt(1),PetscInt(1),PetscInt(2),PetscInt(2))
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
        out = LibPETSc.DMStagGetLocationSlot(petsclib, dm, loc, PetscInt(0))
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

        dm_new = LibPETSc.DMStagCreateCompatibleDMStag(petsclib, dm, PetscInt(1),PetscInt(2),PetscInt(3),PetscInt(4))
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
        loc  = PETSc.LibPETSc.DMStagStencil(dwn,PetscInt(1),PetscInt(1),PetscInt(1),PetscInt(0))
        loc2 = PETSc.LibPETSc.DMStagStencil(dwn,PetscInt(1),PetscInt(1),PetscInt(2),PetscInt(0))
        
        #loc  = PETSc.LibPETSc.DMStagStencil(dwn,1,1,1,0)
        #loc2 = PETSc.LibPETSc.DMStagStencil(dwn,1,1,2,0)
        

        vg = PETSc.DMGlobalVec(dm);
        LibPETSc.DMStagVecSetValuesStencil(petsclib, dm,vg,PetscInt(1),[loc],[PetscScalar(1.0)],PETSc.LibPETSc.INSERT_VALUES) 
        @test vg[4145] == PetscScalar(1.0)

        LibPETSc.DMStagVecSetValuesStencil(petsclib, dm,vg,PetscInt(2), [loc, loc2],[PetscScalar(1.0),PetscScalar(2.0)],PETSc.LibPETSc.INSERT_VALUES) 
        @test extrema(PetscReal.(vg)) == (0.0,2.0)

        out = LibPETSc.DMStagVecGetValuesStencil(petsclib, dm,v, PetscInt(1), [loc])
        @test out == PetscScalar[0.0]
        
        # results in a segfault for some libs
        if isa(PetscScalar, AbstractFloat)
            out = LibPETSc.DMStagVecGetValuesStencil(petsclib, dm,v, PetscInt(2),[loc, loc2])
            @test out == [0.0, 0.0]

            # results in a segfault for some libs
            out = LibPETSc.DMStagStencilToIndexLocal(petsclib, dm,PetscInt(3),PetscInt(2),[loc, loc2])
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
        @test extrema(PetscReal.(v1D)) == (PetscReal(1.0), PetscReal(2.2))

        # Now lets use DMStagVecGetArray. Note that we had to do manual modifcations to that routine
        array2D_1 = LibPETSc.DMStagVecGetArray(petsclib, dm_1D, v1D) 
        array2D_1 .= 2.0
        array2D_1[2,2] = 4.2
        LibPETSc.DMStagVecRestoreArray(petsclib, dm_1D, v1D, array2D_1) 
        @test extrema(PetscReal.(v1D)) == (PetscScalar(2.0), PetscScalar(4.2))
        
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
            points_per_proc = ([PetscInt(20)], [PetscInt(21)], [PetscInt(22)])
        )

        @test LibPETSc.DMStagGetGlobalSizes(petsclib, dm_3D_pd) == (20, 21, 22)
         
        LibPETSc.DMStagSetUniformCoordinatesProduct(petsclib, dm_3D_pd, PetscReal(0.0), PetscReal(1.0), PetscReal(0.0), PetscReal(2.0), PetscReal(0.0), PetscReal(3.0))
#        LibPETSc.DMStagSetUniformCoordinatesProduct(petsclib, dm_3D_pd, 0.0, 1.0, 0.0, 2.0,0.0, 3.0)

        # Retrieve 1D coordinate arrays
        x,y,z = LibPETSc.DMStagGetProductCoordinateArrays(petsclib, dm_3D_pd)
        dims = LibPETSc.DMStagGetGhostCorners(petsclib, dm_3D_pd)[4:6] 
        @test size(x,1) == dims[1]
        @test size(z,1) == dims[3]
#        @test x[10] ≈ PetscScalar(0.225)
#        @test y[10] == PetscScalar(0.42857142857142855)
#        @test z[10] == PetscScalar(0.6136363636363635)

        x[10] = 0.230
     
        LibPETSc.DMStagRestoreProductCoordinateArrays(petsclib, dm_3D_pd, x,y,z)
        
   
        x,y,z = LibPETSc.DMStagGetProductCoordinateArraysRead(petsclib, dm_3D_pd)
        LibPETSc.DMStagRestoreProductCoordinateArraysRead(petsclib, dm_3D_pd, x,y,z)


        slot = LibPETSc.DMStagGetLocationSlot(petsclib, dm_3D_pd, PETSc.LibPETSc.DMSTAG_ELEMENT,0)
        @test slot==7
        slot = LibPETSc.DMStagGetLocationSlot(petsclib, dm_3D_pd, PETSc.LibPETSc.DMSTAG_LEFT,0)
        @test slot==6

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
            points_per_proc = ([PetscInt(32)], [PetscInt(32)])
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
            points_per_proc = ([PetscInt(16)], [PetscInt(16)])
        )    

        xf = PETSc.DMLocalVec(dmf)      
        xf .= 1.0
        
        xc = PETSc.DMLocalVec(dmc)      
        LibPETSc.DMStagRestrictSimple(petsclib, dmf,xf,dmc, xc)
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
    end
end

@testset "DMStagCreate1d" begin

    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)
    for petsclib in PETSc.petsclibs
        #petsclib = PETSc.petsclibs[8]
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

        @test LibPETSc.DMStagGetBoundaryTypes(petsclib, dm)[1]== PETSc.LibPETSc.DM_BOUNDARY_PERIODIC
        PETSc.destroy(dm)

        # Create 1D DMStag with array of local @ of points
        dm = PETSc.DMStag(
            petsclib,
            comm,
            (PETSc.DM_BOUNDARY_NONE,),
            (20,),
            (2,2),
            2;
            points_per_proc = ([PetscInt(20)],),
        )

        # Test get size
        @test LibPETSc.DMStagGetGlobalSizes(petsclib, dm) == (20,0,0)
        @test LibPETSc.DMStagGetLocalSizes(petsclib, dm) == (20,0,0)

        # Test
        @test PETSc.gettype(dm) == "stag"
        @test PETSc.getdimension(dm) == 1

        # Info about ranks  
        @test LibPETSc.DMStagGetIsFirstRank(petsclib, dm) == (false,false,false)
        @test LibPETSc.DMStagGetIsFirstRank(petsclib, dm) == (false,false,false)

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

        ##
        # Create new struct and pass keyword arguments
        dm_1D = PETSc.DMStag(petsclib,comm,(PETSc.DM_BOUNDARY_NONE,),(200,),(2,2),PetscInt(2); stag_grid_x=PetscInt(10));
       
        @test  LibPETSc.DMStagGetGlobalSizes(petsclib,dm_1D)[1] == 10
        @test LibPETSc.DMStagGetEntriesPerElement(petsclib, dm_1D)==4

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
            stag_grid_x = PetscInt(10),
        )


        @test  LibPETSc.DMStagGetStencilWidth(petsclib, dm_ghosted)==2
        corners         = PETSc.getcorners(dm_ghosted)

        @test corners.size[1]==10       # keyword overrides the specified value
        @test  LibPETSc.DMStagGetBoundaryTypes(petsclib, dm_ghosted)[1] == PETSc.DM_BOUNDARY_GHOSTED

        ind = PETSc.DMStagGetIndices(dm_ghosted);
        @test ind.center.x[3] == 5

    
        # simple test to retrieve the KSP object
        # NOTE: need to implement a similar SNES routine
        ksp = PETSc.KSP(dm_ghosted, ksp_type="gmres")
        @test LibPETSc.KSPGetType(petsclib, ksp)=="gmres"

        PETSc.destroy(dm_ghosted)

        PETSc.finalize(petsclib)
    end
end

@testset "DMStagCreate2d" begin

    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)
    for petsclib in PETSc.petsclibs
        #petsclib =  PETSc.petsclibs[1]
        
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


@testset "DMStag Vectors and Coordinates" begin
    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)
    for petsclib in PETSc.petsclibs[1:4]
        #petsclib = PETSc.petsclibs[5]
    
        PETSc.initialize(petsclib)
        PetscScalar = PETSc.scalartype(petsclib)
        PetscInt    = PETSc.inttype(petsclib)
        PetscReal   = real(PetscScalar)
        # Create 1D DMStag
        dm_1D = PETSc.DMStag(
            petsclib,
            comm,
            (PETSc.DM_BOUNDARY_NONE,),
            (200,),
            (2,2),
            1,
            PETSc.DMSTAG_STENCIL_BOX,
            stag_grid_x=PetscInt(10))

        @test LibPETSc.DMStagGetGlobalSizes(petsclib, dm_1D) == (10,0,0)

        # Set coordinates using product (1D) arrays
        #PETSc.setuniformcoordinates_dmstag!(dm_1D, (0,), (10,))
        
        LibPETSc.DMStagSetUniformCoordinatesProduct(petsclib, dm_1D,PetscReal(0.),PetscReal(10.),PetscReal(0.),PetscReal(1.),PetscReal(0.),PetscReal(1.) )


        DMcoord = LibPETSc.DMGetCoordinateDM(petsclib,dm_1D)
        @test PETSc.gettype(DMcoord)=="product"

        # Retrieve array with staggered coordinates
        X_coord,_,_ = LibPETSc.DMStagGetProductCoordinateArrays(petsclib, dm_1D)
        @test  X_coord[1,2] == 0.5

        LibPETSc.DMStagGetLocationSlot(petsclib, dm_1D, LibPETSc.DMSTAG_RIGHT, 0) ==4
     
        global_vec      = PETSc.DMLocalVec(dm_1D)
        local_vec       = PETSc.DMGlobalVec(dm_1D)

        # Fill everything with some data
        fill!(local_vec, mpisize)
        fill!(global_vec, mpisize)
        @test global_vec[3] == 1.0

        # Add the local values to the global values
        LibPETSc.DMLocalToGlobalBegin(petsclib, dm_1D, local_vec, PETSc.ADD_VALUES, global_vec)
        LibPETSc.DMLocalToGlobalEnd(petsclib, dm_1D, local_vec, PETSc.ADD_VALUES, global_vec)

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

        LibPETSc.DMStagSetUniformCoordinatesProduct(petsclib, dm_2D,PetscReal(1.),PetscReal(3.),PetscReal(10.),PetscReal(11.),PetscReal(0.),PetscReal(0.))

        #PETSc.setuniformcoordinates!(dm_2D, (1,3), (10,11))
       
        
        # Retrieve array with staggered coordinates
        X_coord,Y_coord,_ = LibPETSc.DMStagGetProductCoordinateArrays(petsclib, dm_2D)


        vec_test_2D     = PETSc.DMLocalVec(dm_2D)
        X               = LibPETSc.DMStagVecGetArray(petsclib, dm_2D,vec_test_2D);
        X[end,end,end] = 111;                   # modify 3D array @ some point and DOF

        LibPETSc.DMStagVecRestoreArray(petsclib, dm_2D,vec_test_2D,X);
        @test vec_test_2D[end]==111.0           # verify that this modified the vector as well
        Base.finalize(X)                        # release from memory

        #test stencil locations
        pos1 = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_LEFT,1,0,0,1)
        @test pos1.c == 1
        pos2 = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_RIGHT,4,0,0,0)
        pos  = [pos1, pos2]
        @test pos2.loc == LibPETSc.DMSTAG_RIGHT
        @test pos2.i == 4

        # Retrieve value from stencil
        vec_test       = PETSc.DMLocalVec(dm_1D)
        vec_test      .= 1:length(vec_test)                 # point wise copy of data to PetscVec
        val             = LibPETSc.DMStagVecGetValuesStencil(petsclib, dm_1D, vec_test, PetscInt(1), [pos1]) # this gets a single value
        @test val[1] ==6
        vals            = LibPETSc.DMStagVecGetValuesStencil(petsclib, dm_1D, vec_test, PetscInt(2), pos)     # this gets an array of values
        @test vals[1] == 6

        X_1D            = LibPETSc.DMStagVecGetArray(petsclib, dm_1D,vec_test);
#        @test X_1D[2,3] == 14.0


        # Set values using stencils
        vec_test_global = PETSc.DMGlobalVec(dm_1D)
        val1 = PetscScalar.([2222.2, 3.2]);
        LibPETSc.DMStagVecSetValuesStencil(petsclib, dm_1D, vec_test_global, PetscInt(1),  [pos1], val1,   PETSc.INSERT_VALUES)
        @test vec_test_global[6] ≈ 2222.2
        LibPETSc.DMStagVecSetValuesStencil(petsclib, dm_1D, vec_test_global, PetscInt(2), pos,  val1,      PETSc.INSERT_VALUES)
        @test vec_test_global[21] ≈ 3.2

        pos3 = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_LEFT,1,0,0,1)
        val = LibPETSc.DMStagVecGetValuesStencil(petsclib, dm_1D, vec_test, PetscInt(2), [pos3; pos3])
        @test val[2] == 6.0
        PETSc.destroy(dm_1D);

        PETSc.finalize(petsclib)
    end
end



@testset "DMStag create matrixes" begin
    comm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(comm)
    mpisize = MPI.Comm_size(comm)
    for petsclib in PETSc.petsclibs[1:4]
        #petsclib = PETSc.petsclibs[4]
        PETSc.initialize(petsclib)
        PetscScalar = PETSc.scalartype(petsclib)
        PetscInt    = PETSc.inttype(petsclib)
        PetscReal   = real(PetscScalar)

        # Converted from old DMStagCreate1d to new DMStag constructor
        dm_1D = PETSc.DMStag(petsclib, comm, (PETSc.DM_BOUNDARY_NONE,), (200,), (2, 2), 2;
                  stag_grid_x = 10)
        #PETSc.setuniformcoordinates!(dm_1D, (0,), (10,))
        LibPETSc.DMStagSetUniformCoordinatesProduct(petsclib, dm_1D,PetscReal(0.),PetscReal(10.),PetscReal(0.),PetscReal(1.),PetscReal(0.),PetscReal(1.) )

        A = LibPETSc.DMCreateMatrix(petsclib,dm_1D)

        #PETSc.MatSetOption(A, PETSc.MAT_NEW_NONZERO_ALLOCATION_ERR, false)
        LibPETSc.MatSetOption(petsclib,A, LibPETSc.MAT_NEW_NONZERO_ALLOCATION_ERR, false)
        @test size(A) == (42,42)
      

        # set some values using normal indices:
        A[1,1]  = 1.0
        A[1,10] = 1.0

        pos1 = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_LEFT,1,0,0,1)
        pos2 = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_RIGHT,4,0,0,0)
        pos  = [pos1, pos2]
        val1 = PetscScalar.([2222.2 3.2]);
        LibPETSc.DMStagMatSetValuesStencil(petsclib, dm_1D, A, PetscInt(1), [pos1], PetscInt(1), [pos1], [PetscScalar(11.1)], PETSc.INSERT_VALUES)
        LibPETSc.DMStagMatSetValuesStencil(petsclib, dm_1D, A, PetscInt(1), [pos2], PetscInt(2), pos, val1[:], PETSc.INSERT_VALUES)

        @test LibPETSc.MatAssembled(petsclib, A) == false
        PETSc.assemble!(A);
        @test LibPETSc.MatAssembled(petsclib, A) == true
        @test A[1,10] == 1.0


        # Reads a value from the matrix, using the stencil structure
        @test LibPETSc.DMStagMatGetValuesStencil(petsclib, dm_1D, A, PetscInt(1), [pos1], PetscInt(1), [pos1])[1]== PetscScalar(11.1)
   
     
        # result is a 1x2 matrix
        @test LibPETSc.DMStagMatGetValuesStencil(petsclib, dm_1D, A, PetscInt(1), [pos2], PetscInt(2), pos)==val1

            
        PETSc.destroy(A);
        PETSc.destroy(dm_1D);
            
        dofCenter       =   1;
        dofEdge         =   1;
        dofVertex       =   1
        stencilWidth    =   1;
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

        vec_test_2D_global      =   PETSc.DMGlobalVec(dm_2D)
        vec_test_2D_local       =   PETSc.DMLocalVec(dm_2D)
        fill!(vec_test_2D_global, 0.0)
        fill!(vec_test_2D_local, 0.0)
        corners                 =   PETSc.getcorners(dm_2D)
        ghost_corners           =   PETSc.getghostcorners(dm_2D)


        for ix=corners.lower[1]:corners.upper[1]
            for iy=corners.lower[2]:corners.upper[2]
                local dof
                # DOF at the center point
                dof     = 0;
                posA    = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_ELEMENT,ix,iy,0,dof)
                value   = PetscScalar(ix);
                LibPETSc.DMStagVecSetValuesStencil(petsclib, dm_2D, vec_test_2D_global, PetscInt(1), [posA], [value], PETSc.INSERT_VALUES)
                
                dof     = 0;
                posB    = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_LEFT,ix,iy,0,dof)
                value   = PetscScalar(33);
                LibPETSc.DMStagVecSetValuesStencil(petsclib, dm_2D, vec_test_2D_global, PetscInt(1), [posB], [value], PETSc.INSERT_VALUES)
                
                dof     = 0;
                posC    = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_DOWN,ix,iy,0,dof)
                value   = PetscScalar(44);
                LibPETSc.DMStagVecSetValuesStencil(petsclib, dm_2D, vec_test_2D_global, PetscInt(1), [posC], [value], PETSc.INSERT_VALUES)

                dof     = 0;
                posC    = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_BACK_UP_RIGHT,ix,iy,0,dof)
                value   = PetscScalar(55);
                LibPETSc.DMStagVecSetValuesStencil(petsclib, dm_2D, vec_test_2D_global, PetscInt(1), [posC], [value], PETSc.INSERT_VALUES)
                
            end
        end
        
        PETSc.assemble!(vec_test_2D_global) # assemble global vector

        
        # Add the global values to the local values
        LibPETSc.DMGlobalToLocalBegin(petsclib, dm_2D, vec_test_2D_global, LibPETSc.INSERT_VALUES, vec_test_2D_local)
        LibPETSc.DMGlobalToLocalEnd(petsclib, dm_2D, vec_test_2D_global, LibPETSc.INSERT_VALUES, vec_test_2D_local)
        
        # retrieve value back from the local array and check that it agrees with global one
        dof     = 0;
        pos     = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_ELEMENT,3,2,0,dof)
        @test LibPETSc.DMStagVecGetValuesStencil(petsclib, dm_2D, vec_test_2D_local, PetscInt(1), [pos])[1] == 3.0

        # Extract an array that holds all DOF's
        X2D_dofs  = LibPETSc.DMStagVecGetArray(petsclib, dm_2D,vec_test_2D_local)           # extract arrays with all DOF (mostly for visualizing)
        
        @test X2D_dofs[4,4,1] ≈ PetscScalar(55.0)
        @test X2D_dofs[4,4,2] ≈ PetscScalar(44.0)
        @test X2D_dofs[4,4,3] ≈ PetscScalar(33.0)
        @test X2D_dofs[4,4,4] ≈ PetscScalar(3.0)


          
        # cleanup
        PETSc.destroy(dm_2D);
        
           
        PETSc.finalize(petsclib)
    end
end

