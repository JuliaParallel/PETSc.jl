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

    function DMSetUp(dm::DMStag{$PetscScalar})

        @chk ccall((:DMSetUp, $libpetsc), PetscErrorCode, (CDMStag, ), dm )

        return nothing
    end
   
    function setfromoptions!(dm::DMStag{$PetscScalar})

        @chk ccall((:DMSetFromOptions, $libpetsc), PetscErrorCode, (CDMStag, ), dm )

        return nothing
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
    function DMCreateGlobalVector(dm::DMStag)


        v = VecSeq(C_NULL, dm.comm, [0.0])  # empty vector
        
        ccall((:DMCreateGlobalVector, $libpetsc), PetscErrorCode, (CDMStag, Ptr{CVec}), dm, v)

        # extract array of values from new vector
        v.array = unsafe_localarray($PetscScalar, v.ptr; write=true)
        
        return v
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
