# Attempt at include dmstag functions
const CDMStag = Ptr{Cvoid}
const CDMStagType = Cstring

mutable struct DMStag{T} <: Factorization{T}
    ptr::CDMStag
    comm::MPI.Comm
    dim::Int64
    opts::Options{T}
end

# allows us to pass XXMat objects directly into CMat ccall signatures
Base.cconvert(::Type{CDMStag}, obj::DMStag) = obj.ptr
# allows us to pass XXMat objects directly into Ptr{CMat} ccall signatures
Base.unsafe_convert(::Type{Ptr{CDMStag}}, obj::DMStag) =
    convert(Ptr{CDMStag}, pointer_from_objref(obj))

Base.eltype(::DMStag{T}) where {T} = T

@for_libpetsc begin

    """
        Creates a 1D DMStag object
        
        Usage:

            dm = DMStagCreate1d(comm::MPI.Comm, bndx::DMBoundaryType, M, dofVertex, dofCenter, stencilType::DMStagStencilType=DMSTAG_STENCIL_BOX, stencilWidth=1, lx::Vector=[]; kwargs...)

                comm            -   MPI communicator
                bndx            -   boundary type: DM_BOUNDARY_NONE, DM_BOUNDARY_PERIODIC, or DM_BOUNDARY_GHOSTED. 
                M               -   global number of grid points
                dofVertex       -   [=1] number of degrees of freedom per vertex/point/node/0-cell
                dofCenter       -   [=1] number of degrees of freedom per element/edge/1-cell
                stencilType     -   ghost/halo region type: DMSTAG_STENCIL_BOX or DMSTAG_STENCIL_NONE
                stencilWidth    -   width, in elements, of halo/ghost region
                lx              -   [Optional] Vector of local sizes, of length equal to the comm size, summing to M
                kwargs...       -   [Optional] keyword arguments (see PETSc webpage), specifiable as stag_grid_x=100, etc. 
    """
    function DMStagCreate1d(comm::MPI.Comm, bndx::DMBoundaryType, M, dofVertex=1,dofCenter=1,stencilType::DMStagStencilType=DMSTAG_STENCIL_BOX,stencilWidth=2, lx=C_NULL; kwargs...)

        if isempty(lx); lx = C_NULL; end
        opts = Options{$PetscScalar}(kwargs...)

        dm  = DMStag{$PetscScalar}(C_NULL, comm, 1, opts)   # retrieve options
        
        @chk ccall((:DMStagCreate1d, $libpetsc), PetscErrorCode,
                (MPI.MPI_Comm, DMBoundaryType, $PetscInt, $PetscInt, $PetscInt, DMStagStencilType, $PetscInt,  Ptr{$PetscInt}, Ptr{CDMStag}),
                comm, bndx, M,dofVertex,dofCenter,stencilType,stencilWidth,lx, dm )

        with(dm.opts) do
            setfromoptions!(dm)
        end

        DMSetUp(dm);

        if comm == MPI.COMM_SELF
            finalizer(destroy, dm)
        end
        
        return dm
    end

    """
        Creates a 2D DMStag object
        
        Usage:

            dm = DMStagCreate2d(comm::MPI.Comm, bndx::DMBoundaryType, bndy::DMBoundaryType, M, N, m, n, dofVertex, dofEdge, dofElement, stencilType::DMStagStencilType=DMSTAG_STENCIL_BOX, stencilWidth, lx, ly; kwargs...)

                comm            -   MPI communicator
                bndx,bndy       -   boundary type: DM_BOUNDARY_NONE, DM_BOUNDARY_PERIODIC, or DM_BOUNDARY_GHOSTED. 
                M,N             -   global number of grid points
                m,n             -   number of ranks in the x,y directions (may be PETSC_DECIDE TO do) 
                dofVertex       -   [=1] number of degrees of freedom per vertex/point/node/0-cell
                dofEdge         -   [=1] number of degrees of freedom per edge/1-cell 
                dofElement      -   [=1] number of degrees of freedom per element/2-cell 
                stencilType     -   ghost/halo region type: DMSTAG_STENCIL_BOX or DMSTAG_STENCIL_NONE
                stencilWidth    -   width, in elements, of halo/ghost region
                lx,ly           -   [Optional] arrays of local x,y element counts, of length equal to m,n, summing to M,N 
                kwargs...       -   [Optional] keyword arguments (see PETSc webpage), specifiable as stag_grid_x=100, etc. 
    """
    function DMStagCreate2d(comm::MPI.Comm, bndx::DMBoundaryType, bndy::DMBoundaryType, M, N, m=C_NULL, n=C_NULL, dofVertex=1, dofEdge=1, dofElement=1, stencilType::DMStagStencilType=DMSTAG_STENCIL_BOX, stencilWidth=2, lx=C_NULL, ly=C_NULL; kwargs...)
        
        if isempty(lx); lx = C_NULL; end
        if isempty(ly); ly = C_NULL; end
        opts = Options{$PetscScalar}(kwargs...)
        
        dm = DMStag{$PetscScalar}(C_NULL, comm, 2, opts)

        @chk ccall((:DMStagCreate2d, $libpetsc), PetscErrorCode,
                (MPI.MPI_Comm, DMBoundaryType, DMBoundaryType, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, DMStagStencilType, $PetscInt, Ptr{$PetscInt},  Ptr{$PetscInt}, Ptr{CDMStag}),
                comm, bndx, bndy, M, N, m, n ,dofVertex ,dofEdge ,dofElement ,stencilType ,stencilWidth ,lx ,ly ,dm )
        
        with(dm.opts) do
            setfromoptions!(dm)
        end

        DMSetUp(dm);

        if comm == MPI.COMM_SELF
            finalizer(destroy, dm)
        end
        
        return dm
    end

    """
        Creates a 3D DMStag object
        
        Usage:

            dm = DMStagCreate3d(comm::MPI.Comm, bndx::DMBoundaryType, bndy::DMBoundaryType, bndz::DMBoundaryType, M, N, P, m, n, p, dofVertex, dofEdge, dofFace, dofElement, stencilType::DMStagStencilType=DMSTAG_STENCIL_BOX, stencilWidth, lx, ly, lz; kwargs...)

                comm            -   MPI communicator
                bndx,bndy,bndz  -   boundary type: DM_BOUNDARY_NONE, DM_BOUNDARY_PERIODIC, or DM_BOUNDARY_GHOSTED. 
                M,N,P           -   global number of grid points
                m,n,p           -   number of ranks in the x,y directions (may be PETSC_DECIDE TO do) 
                dofVertex       -   [=1] number of degrees of freedom per vertex/point/node/0-cell
                dofEdge         -   [=1] number of degrees of freedom per edge/1-cell 
                dofFace         -   [=1] number of degrees of freedom per face/2-cell 
                dofElement      -   [=1] number of degrees of freedom per element/3-cell 
                stencilType     -   ghost/halo region type: DMSTAG_STENCIL_BOX or DMSTAG_STENCIL_NONE
                stencilWidth    -   width, in elements, of halo/ghost region
                lx,ly,lz        -   [Optional] arrays of local x,y element counts, of length equal to m,n, summing to M,N 
                kwargs...       -   [Optional] keyword arguments (see PETSc webpage), specifiable as stag_grid_x=100, etc. 
    """
    function DMStagCreate3d(comm::MPI.Comm, bndx::DMBoundaryType, bndy::DMBoundaryType, bndz::DMBoundaryType, M, N, P, m=C_NULL, n=C_NULL, p=C_NULL, dofVertex=1, dofEdge=1, dofFace=1, dofElement=1, stencilType::DMStagStencilType=DMSTAG_STENCIL_BOX, stencilWidth=2, lx=C_NULL, ly=C_NULL, lz=C_NULL; kwargs...)
        
        if isempty(lx); lx = C_NULL; end
        if isempty(ly); ly = C_NULL; end
        if isempty(lz); lz = C_NULL; end
        opts = Options{$PetscScalar}(kwargs...)
        
        dm = DMStag{$PetscScalar}(C_NULL, comm, 3, opts)

        @chk ccall((:DMStagCreate3d, $libpetsc), PetscErrorCode,
                (MPI.MPI_Comm, DMBoundaryType, DMBoundaryType, DMBoundaryType, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, DMStagStencilType, $PetscInt, Ptr{$PetscInt},  Ptr{$PetscInt},  Ptr{$PetscInt}, Ptr{CDMStag}),
                comm, bndx, bndy, bndz, M, N, P, m, n ,p ,dofVertex ,dofEdge ,dofFace ,dofElement ,stencilType ,stencilWidth ,lx ,ly ,lz ,dm )
        
        with(dm.opts) do
            setfromoptions!(dm)
        end

        DMSetUp(dm);

        if comm == MPI.COMM_SELF
            finalizer(destroy, dm)
        end
        
        return dm
    end

    """
        sets up the data structures inside a DM object 
        
        Usage:

            DMSetUp(dm::DMStag)

                dm              -   the DMStag object 
    """
    function DMSetUp(dm::DMStag{$PetscScalar})

        @chk ccall((:DMSetUp, $libpetsc), PetscErrorCode, (CDMStag, ), dm )

        return nothing
    end

    """
        sets parameters in a DM from the options database 
        
        Usage:

            setfromoptions!(dm::DMStag)

                dm              -   the DMStag object 
    """
    function setfromoptions!(dm::DMStag{$PetscScalar})

        @chk ccall((:DMSetFromOptions, $libpetsc), PetscErrorCode, (CDMStag, ), dm )

        return nothing
    end


    
    """
        Creates a compatible DMStag with different dof/stratum 
        
        Usage:

            dm = DMStagCreateCompatibleDMStag(dm::DMStag, dofVertex, dofEdge, dofFace, dofElement; kwargs...)

                dm              -   the DMStag object 
                dofVertex       -   [=0] number of degrees of freedom per vertex/point/node/0-cell
                dofEdge         -   [=0] number of degrees of freedom per edge/1-cell 
                dofFace         -   [=0] number of degrees of freedom per face/2-cell 
                dofElement      -   [=0] number of degrees of freedom per element/3-cell 
                kwargs...       -   [Optional] keyword arguments (see PETSc webpage), specifiable as stag_grid_x=100, etc. 
    """
    function DMStagCreateCompatibleDMStag(dm::DMStag{$PetscScalar}, dofVertex=0, dofEdge=0, dofFace=0, dofElement=0; kwargs...)

        comm  = MPI.COMM_SELF

        dim   = DMGetDimension(dm)

        opts  = Options{$PetscScalar}(kwargs...)

        dmnew = DMStag{$PetscScalar}(C_NULL, comm, dim, opts)

        @chk ccall((:DMStagCreateCompatibleDMStag, $libpetsc), PetscErrorCode, 
        (CDMStag, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{CDMStag}), 
        dm, dofVertex, dofEdge, dofFace, dofElement, dmnew)

        with(dm.opts) do
            setfromoptions!(dmnew)
        end

        DMSetUp(dmnew);

        if comm == MPI.COMM_SELF
            finalizer(destroy, dmnew)
        end
        
        return dmnew

    end

    """
        Get number of DOF associated with each stratum of the grid 
        
        Usage:

            dofVertex, dofEdge, dofFace, dofElement = DMStagGetDOF(dm::DMStag, dofVertex, dofEdge, dofFace, dofElement; kwargs...)

                dm              -   the DMStag object 
                dof0 	- the number of points per 0-cell (vertex/node)
	            dof1 	- the number of points per 1-cell (element in 1D, edge in 2D and 3D)
	            dof2 	- the number of points per 2-cell (element in 2D, face in 3D)
	            dof3 	- the number of points per 3-cell (element in 3D) 

    """
    function DMStagGetDOF(dm::DMStag{$PetscScalar})

        dof0 = Ref{$PetscInt}()
        dof1 = Ref{$PetscInt}()
        dof2 = Ref{$PetscInt}()
        dof3 = Ref{$PetscInt}()

        @chk ccall((:DMStagGetDOF, $libpetsc), PetscErrorCode, 
        (CDMStag, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}), 
        dm, dof0, dof1, dof2, dof3)

        if dm.dim==1
            return  dof0[],dof1[]   
        elseif dm.dim==2
            return dof0[],dof1[],dof2[]
        elseif dm.dim==3
            return dof0[],dof1[],dof2[],dof3[]
        end

    end


    """
        Gets the global size of the DMStag object
            M,N,P = DMStagGetGlobalSizes(dm::DMStag)
    """
    function DMStagGetGlobalSizes(dm::DMStag{$PetscScalar})

        M = Ref{$PetscInt}()
        N = Ref{$PetscInt}()
        P = Ref{$PetscInt}()

        @chk ccall((:DMStagGetGlobalSizes, $libpetsc), PetscErrorCode,
            (CDMStag, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}), 
            dm, M, N, P )
        
        if dm.dim==1    
            return M[]
        elseif dm.dim==2
            return M[], N[] 
        elseif dm.dim==3
            return M[], N[], P[]    
        end
    end


    """
        Sets coordinates for a DMStag object using the Product method to specify coordinates (1D arrays)
    """
    function DMStagSetUniformCoordinatesProduct(dm::DMStag, xmin, xmax, ymin, ymax, zmin, zmax)
        
        @chk ccall((:DMStagSetUniformCoordinatesProduct, $libpetsc), PetscErrorCode,
                    ( CDMStag,   $PetscScalar, $PetscScalar, $PetscScalar, 
                                $PetscScalar, $PetscScalar, $PetscScalar), 
                            dm, xmin, xmax, ymin, ymax, zmin, zmax)

        return nothing
    end

    """
        Sets uniform coordinates for a 1D DMStag object 
            DMStagSetUniformCoordinates(dm::DMStag, xmin, xmax)
    """
    function DMStagSetUniformCoordinates(dm::DMStag, xmin, xmax)
        DMStagSetUniformCoordinatesProduct(dm::DMStag, xmin, xmax, 0.0, 0.0, 0.0, 0.0);
        return nothing
    end

    """
        Sets uniform coordinates for a 2D DMStag object 
            DMStagSetUniformCoordinates(dm::DMStag, xmin, xmax, ymin, ymax)
    """
    function DMStagSetUniformCoordinates(dm::DMStag, xmin, xmax, ymin, ymax)
        DMStagSetUniformCoordinatesProduct(dm::DMStag, xmin, xmax, ymin, ymax, 0.0, 0.0);
        return nothing
    end

    """
        Sets uniform coordinates for a 3D DMStag object 
            DMStagSetUniformCoordinates(dm::DMStag, xmin, xmax, ymin, ymax, zmin, zmax)
    """
    function DMStagSetUniformCoordinates(dm::DMStag, xmin, xmax, ymin, ymax, zmin, zmax)
        DMStagSetUniformCoordinatesProduct(dm::DMStag, xmin, xmax, ymin, ymax, zmin, zmax);
        return nothing
    end



    # NOT WORKING YET
    function DMStagGetProductCoordinateArrays(dm::DMStag)
        
        arrX = Ref{$PetscScalar}()
        arrY = Ref{$PetscScalar}()
        arrZ = Ref{$PetscScalar}()

        @chk ccall((:DMStagGetProductCoordinateArrays, $libpetsc), PetscErrorCode,
            ( CDMStag,   Ptr{$PetscScalar}, Ptr{$PetscScalar}, Ptr{$PetscScalar}), 
                dm, arrX, arrY, arrZ)


        
        return arrX, arrY, arrZ        

    end

    """
        This extracts a global vector from the DMStag object
            NOTE: for now this is initialized sequentially; MPI should be added
    """
    function DMCreateGlobalVector(dm::DMStag; write_val=true, read_val=true)

        v = VecSeq(C_NULL, dm.comm, [0.0])  # empty vector
        
        ccall((:DMCreateGlobalVector, $libpetsc), PetscErrorCode, (CDMStag, Ptr{CVec}), dm, v)

        # Link a julia array to the values from the new vector
        # If we modify values here, it will automatically be changed in the PetcVec as well
        v.array = unsafe_localarray($PetscScalar, v.ptr; write=write_val, read=read_val)
        
        return v
    end

    """
    This extracts a local vector from the DMStag object
            NOTE: for now this is initialized sequentially; MPI should be added
    """
    function DMCreateLocalVector(dm::DMStag; write_val=true, read_val=true)

        v = VecSeq(C_NULL, dm.comm, [0.0])  # empty vector
        
        ccall((:DMCreateLocalVector, $libpetsc), PetscErrorCode, (CDMStag, Ptr{CVec}), dm, v)

        # Link a julia array to the values from the new vector
        # If we modify values here, it will automatically be changed in the PetcVec as well
        v.array = unsafe_localarray($PetscScalar, v.ptr; write=write_val, read=read_val)
        
        return v
    end

    function DMStagVecGetArray(dm::DMStag, cv::CVec)
        ## NOT WORKING YET!

        # the following works, so lets reproduce it first
        #X = PETSc.unsafe_localarray(Float64, cv; read=true, write=false)
        if 1==0
            # copy of what we do before, but spelled out (not for array)
            r_pv = Ref{Ptr{$PetscScalar}}()
     
            @chk ccall((:DMStagVecGetArrayRead, $libpetsc), PetscErrorCode,
                (CDMStag, CVec, Ptr{$PetscScalar}), dm, cv, r_pv)
            
            len_array = localsize(cv);    

            X = unsafe_wrap(Array, r_pv[], len_array*2; own = false)


        end

        if 1==1
        
            #  r_pv = Ref{Ptr{$PetscScalar}}()
            
            len_array = localsize(cv);    

            #VecGetArray2d(Vec x,PetscInt m,PetscInt n,PetscInt mstart,PetscInt nstart,PetscScalar **a[])

            
            m = len_array;
            n = 2;
            mstart = 0
            nstart = 0


            entriesPerElement = PETSc.DMStagGetEntriesPerElement(dm)
            nGhost = PETSc.DMStagGetGhostCorners(dm)
            @show nGhost, entriesPerElement
            X = zeros(nGhost[2],entriesPerElement);

            @chk ccall((:VecGetArray2d, $libpetsc), PetscErrorCode,
                (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscScalar}), 
                    cv, nGhost[2],entriesPerElement, nGhost[1], 0, X)
            
            
            #X = unsafe_wrap(Array, r_pv[], len_array; own = false)


        end
      

        return X

    end

    """
        Retrieves a coordinate slot from a DMStag object, if the coordinates are set as ProductCoordinate 
            slot = DMStagGetProductCoordinateLocationSlot(dm::DMStag,loc::DMStagStencilLocation)
    """
    function DMStagGetProductCoordinateLocationSlot(dm::DMStag,loc::DMStagStencilLocation)
        slot = Ref{$PetscInt}()
        @chk ccall((:DMStagGetProductCoordinateLocationSlot, $libpetsc), PetscErrorCode,
                    ( CDMStag,   DMStagStencilLocation, Ptr{$PetscInt}), dm, loc, slot)

        return slot[]
    end


    function DMStagGetEntriesPerElement(dm::DMStag)
        entriesPerElement = Ref{$PetscInt}()
        @chk ccall((:DMStagGetEntriesPerElement, $libpetsc), PetscErrorCode,
                    ( CDMStag,  Ptr{$PetscInt}), dm,  entriesPerElement)

        return entriesPerElement[]
    end

    """
    Retrieves a coordinate slot from a DMStag object, if the coordinates are set as ProductCoordinate 

        slot = DMStagGetLocationSlot(dm::DMStag,loc::DMStagStencilLocation, c)
        
        Input Parameters
            dm	    - the DMStag object
            loc	    - location relative to an element
            c	    - component
        
        Output Parameter

            slot	- index to use

    """
    function DMStagGetLocationSlot(dm::DMStag,loc::DMStagStencilLocation, c)
        
        slot = Ref{$PetscInt}()
        @chk ccall((:DMStagGetLocationSlot, $libpetsc), PetscErrorCode,
                    ( CDMStag,   DMStagStencilLocation, $PetscInt, Ptr{$PetscInt}), dm, loc, c, slot)

        return slot[]
    end

    """
        Destroys the DMStag object
    """
    function destroy(dm::DMStag{$PetscScalar})
        finalized($PetscScalar) ||
            @chk ccall((:DMDestroy, $libpetsc), PetscErrorCode, (Ptr{CDMStag},), dm)
        return nothing
    end

    """
        Retrieves the Type of the DMStag object
    """
    function gettype(dm::DMStag{$PetscScalar})
        t_r = Ref{CDMStagType}()
        @chk ccall((:DMGetType, $libpetsc), PetscErrorCode, (CDMStag, Ptr{CDMStagType}), dm, t_r)
        return unsafe_string(t_r[])
    end

    function view(dm::DMStag{$PetscScalar}, viewer::Viewer{$PetscScalar}=ViewerStdout{$PetscScalar}(dm.comm))
        @chk ccall((:DMView, $libpetsc), PetscErrorCode, 
                    (CDMStag, CPetscViewer),
                dm, viewer);
        return nothing
    end

    """ 
        Gets the corners of the DMStag grid
            x,m,nExtrax = DMStagGetCorners(dm:DMStag)   in 1D
    """
    function  DMStagGetCorners(dm::DMStag)

        x = Ref{$PetscInt}()
        y = Ref{$PetscInt}()
        z = Ref{$PetscInt}()
        m = Ref{$PetscInt}()
        n = Ref{$PetscInt}()
        p = Ref{$PetscInt}()
        nExtrax = Ref{$PetscInt}()
        nExtray = Ref{$PetscInt}()
        nExtraz = Ref{$PetscInt}()
        
        @chk ccall((:DMStagGetCorners, $libpetsc), PetscErrorCode,
            (CDMStag,   Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt},
                        Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt},
                        Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}), 
            dm, x,y,z, m,n,p, nExtrax,nExtray,nExtraz )

            if dm.dim==1
                return x[], m[], nExtrax[]    
            elseif dm.dim==2
                return x[], y[], m[],n[], nExtrax[],nExtray[]    
            elseif dm.dim==3
                return x[], y[], z[], m[],n[],p[], nExtrax[],nExtray[],nExtraz[]    
            end
    end
    
    """ 
    Gets the corners of the DMStag grid including the ghost nodes
        x,m = DMStagGetGhostCorners(dm:DMStag)   in 1D
    """
    function  DMStagGetGhostCorners(dm::DMStag)

        x = Ref{$PetscInt}()
        y = Ref{$PetscInt}()
        z = Ref{$PetscInt}()
        m = Ref{$PetscInt}()
        n = Ref{$PetscInt}()
        p = Ref{$PetscInt}()
        
        @chk ccall((:DMStagGetGhostCorners, $libpetsc), PetscErrorCode,
            (CDMStag,   Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt},
                        Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}), 
            dm, x,y,z, m,n,p)

            if dm.dim==1
                return x[], m[]  
            elseif dm.dim==2
                return x[], y[], m[],n[]    
            elseif dm.dim==3
                return x[], y[], z[], m[],n[],p[]   
            end
    end

    """
        returns the types of the boundary of the DMStag object in x/y/z direction 
            Bx = DMStagGetBoundaryTypes(dm::DMStag) in 1D
    """
    function  DMStagGetBoundaryTypes(dm::DMStag)

        Bx = Ref{$DMBoundaryType}()
        By = Ref{$DMBoundaryType}()
        Bz = Ref{$DMBoundaryType}()
      
        @chk ccall((:DMStagGetBoundaryTypes, $libpetsc), PetscErrorCode,
            (CDMStag,   Ptr{$DMBoundaryType}, Ptr{$DMBoundaryType}, Ptr{$DMBoundaryType}), dm, Bx,By,Bz)

            if dm.dim==1
                return Bx[]    
            elseif dm.dim==2
                return Bx[], By[]
            elseif dm.dim==3
                return Bx[], By[], Bz[]
            end
    end

    function  DMStagSetStencilWidth(dm::DMStag, stencilWidth::Int64)

        @chk ccall((:DMStagSetStencilWidth, $libpetsc), PetscErrorCode,
             (CDMStag,  $PetscInt), dm, stencilWidth)

        return nothing
    end


 

    """
        returns the # of dimensions of the DMStag object
    """
    function DMGetDimension(dm::DMStag)
        dim = Ref{$PetscInt}()

        @chk ccall((:DMGetDimension, $libpetsc), PetscErrorCode, (CDMStag,Ptr{$PetscInt}), dm, dim )

        return dim[]
    end

end

Base.show(io::IO, dm::DMStag) = _show(io, dm)
