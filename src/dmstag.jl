# Attempt at include dmstag functions
const CDMStag = Ptr{Cvoid}
const CDMStagType = Cstring

mutable struct DMStag{T} <: Factorization{T}
    ptr::CDMStag
    comm::MPI.Comm
    dim::Int64
    opts::Options{T}
end

mutable struct DMSTAGSTENCIL
    loc::DMStagStencilLocation
    i::Int64
    j::Int64
    k::Int64
    c::Int64
end

mutable struct DMSTAGSTENCIL_C
    loc::DMStagStencilLocation
    i::Cint
    j::Cint
    k::Cint
    c::Cint
end

const DMStagStencil     = DMSTAGSTENCIL
const DMStagStencil_c   = DMSTAGSTENCIL_C


# allows us to pass XXMat objects directly into CMat ccall signatures
Base.cconvert(::Type{CDMStag}, obj::DMStag) = obj.ptr
# allows us to pass XXMat objects directly into Ptr{CMat} ccall signatures
Base.unsafe_convert(::Type{Ptr{CDMStag}}, obj::DMStag) =
    convert(Ptr{CDMStag}, pointer_from_objref(obj))


Base.eltype(::DMStag{T}) where {T} = T

# allows us to pass XXMat objects directly into CMat ccall signatures
#Base.cconvert(::Type{DMStagStencil_c}, obj::Ref{DMStagStencil}) = obj
#Base.cconvert(::Type{DMStagStencil_c}, v::DMStagStencil) = DMStagStencil_c(v.loc, v.i, v.j,v.k, v.c)

Base.convert(::Type{DMStagStencil_c}, v::DMStagStencil) = DMStagStencil_c(v.loc, v.i, v.j,v.k, v.c)
#Base.unsafe_convert(::Type{DMStagStencil_c}, v::Tuple) = DMStagStencil_c(v[1], v[2], v[3], v[4], v[5]);


@for_libpetsc begin

    """
    Creates a 1D DMStag object

        dm = DMStagCreate1d(comm::MPI.Comm, bndx::DMBoundaryType, M, dofVertex, dofCenter, stencilType::DMStagStencilType=DMSTAG_STENCIL_BOX, stencilWidth=2, lx=C_NULL; kwargs...)

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

        DMSetUp(dm::DMStag)

            dm    -   the DMStag object 
    """
    function DMSetUp(dm::DMStag{$PetscScalar})

        @chk ccall((:DMSetUp, $libpetsc), PetscErrorCode, (CDMStag, ), dm )

        return nothing
    end

    """
    sets parameters in a DM from the options database 

        setfromoptions!(dm::DMStag)

            dm              -   the DMStag object 
    """
    function setfromoptions!(dm::DMStag{$PetscScalar})

        @chk ccall((:DMSetFromOptions, $libpetsc), PetscErrorCode, (CDMStag, ), dm )

        return nothing
    end


    
    """
    Creates a compatible DMStag with different dof/stratum 

        dm = DMStagCreateCompatibleDMStag(dm::DMStag, dofVertex, dofEdge, dofFace, dofElement; kwargs...)

            dm              -   the DMStag object 
            dofVertex       -   [=0] number of degrees of freedom per vertex/point/node/0-cell
            dofEdge         -   [=0] number of degrees of freedom per edge/1-cell 
            dofFace         -   [=0] number of degrees of freedom per face/2-cell 
            dofElement      -   [=0] number of degrees of freedom per element/3-cell 
            kwargs...       -   [Optional] keyword arguments (see PETSc webpage), specifiable as stag_grid_x=100, etc. 
    """
    function DMStagCreateCompatibleDMStag(dm::DMStag{$PetscScalar}, dofVertex=0, dofEdge=0, dofFace=0, dofElement=0; kwargs...)

        comm  = dm.comm

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

        dof0, dof1, dof2, dof3 = DMStagGetDOF(dm::DMStag)

            dm      - the DMStag object 
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

            dm      - the DMStag object 
            M,N,P   - size in x,y,z
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


    function Base.size(dm::DMStag{$PetscScalar})
        size = DMStagGetGlobalSizes(dm)
        return size
    end


    """
    Gets the local size of the DMStag object

        M,N,P = DMStagGetLocalSizes(dm::DMStag)

            dm      - the DMStag object 
            M,N,P   - size in x,y,z
    """
    function DMStagGetLocalSizes(dm::DMStag{$PetscScalar})

        M = Ref{$PetscInt}()
        N = Ref{$PetscInt}()
        P = Ref{$PetscInt}()

        @chk ccall((:DMStagGetLocalSizes, $libpetsc), PetscErrorCode,
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
    get number of native entries in the global representation 

        entries = DMStagGetEntries(dm::DMStag)

            dm      - the DMStag object 
            entries - number of rank-native entries in the global representation 
    """

    function DMStagGetEntries(dm::DMStag)
        # doesn't work, as it cannot find the name of the routine. 
        # Changed in more recent PETSc version?

        entries = Ref{$PetscInt}()
        @chk ccall((:DMStagGetEntries, $libpetsc), PetscErrorCode,
                    ( CDMStag,  Ptr{$PetscInt}), dm,  entries)

        return entries[]
    end


    """
    Set the coordinate DM to be a DMProduct of 1D DMStag objects, each of which have a coordinate DM (also a 1d DMStag) holding uniform coordinates. 

        DMStagSetUniformCoordinatesProduct(dm::DMStag, xmin, xmax, ymin=0, ymax=0, zmin=0, zmax=0)

            dm 	                           - the DMStag object 
            xmin,xmax,ymin,ymax,zmin,zmax  - maximum and minimum global coordinate values
    """
    function DMStagSetUniformCoordinatesProduct(dm::DMStag, xmin, xmax, ymin=0, ymax=0, zmin=0, zmax=0)
        
        @chk ccall((:DMStagSetUniformCoordinatesProduct, $libpetsc), PetscErrorCode,
                    ( CDMStag,   $PetscScalar, $PetscScalar, $PetscScalar, 
                                $PetscScalar, $PetscScalar, $PetscScalar), 
                            dm, xmin, xmax, ymin, ymax, zmin, zmax)

        return nothing
    end

    """
    set DMStag coordinates to be a uniform grid, storing all values

        DMStagSetUniformCoordinatesExplicit(dm::DMStag, xmin, xmax, ymin=0, ymax=0, zmin=0, zmax=0)

            dm 	                          - the DMStag object
	        xmin,xmax,ymin,ymax,zmin,zmax - maximum and minimum global coordinate values 
    """
    function DMStagSetUniformCoordinatesExplicit(dm::DMStag, xmin, xmax, ymin=0, ymax=0, zmin=0, zmax=0)
        
        @chk ccall((:DMStagSetUniformCoordinatesExplicit, $libpetsc), PetscErrorCode,
                    ( CDMStag,   $PetscScalar, $PetscScalar, $PetscScalar, 
                                $PetscScalar, $PetscScalar, $PetscScalar), 
                            dm, xmin, xmax, ymin, ymax, zmin, zmax)

        return nothing
    end

    """
    Creates a global vector from a DM object 
    NOTE: for now this is initialized sequentially; MPI should be added

        vec = DMCreateGlobalVector(dm::DMStag; write_val=true, read_val=true)

            dm 	- the DM object
            vec - the global vector
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
    Creates a local vector from a DM object 
    NOTE: for now this is initialized sequentially; MPI should be added

        vec = DMCreateLocalVector(dm::DMStag; write_val=true, read_val=true)

            dm 	- the DM object
            vec - the local vector
    """
    function DMCreateLocalVector(dm::DMStag; write_val=true, read_val=true)

        v = VecSeq(C_NULL, dm.comm, [0.0])  # empty vector
        
        ccall((:DMCreateLocalVector, $libpetsc), PetscErrorCode, (CDMStag, Ptr{CVec}), dm, v)

        # Link a julia array to the values from the new vector
        # If we modify values here, it will automatically be changed in the PetcVec as well
        v.array = unsafe_localarray($PetscScalar, v.ptr; write=write_val, read=read_val)
        
        return v
    end 


    """
    get access to local array (including ghost points) of the DMStag
            
        Array =  DMStagVecGetArray(dm::DMStag, v::AbstractVec)

            dm 	  - the DMStag object
            vec   - the Vec object 
            Array - the array

    Once you are done with work on the array, you MUST release the memory with
                
        Base.finalize(Array)

    Otherwise the values are not returned correctly to v    

    """
    function DMStagVecGetArray(dm::DMStag, v::AbstractVec)
        # Note: there is actually no need to call PETSc again, as Julia has the possibility 
        # to wrap an existing array into another one. Our vec already has the array wrapper, 
        # so we reshape that 

        # Extract array from vector. Note: we need to release this by calling 
        # Base.finalize on X1!
        v.array     =   unsafe_localarray($PetscScalar, v.ptr;  write=true, read=true)

        X1          =   DMStagVecGetArray(dm, v.array) 
        
        return X1
    end

    """
    get read-only access to a local array (including ghost points) of the DMStag

        Array =  DMStagVecGetArrayRead(dm::DMStag, v::AbstractVec)

            dm 	  - the DMStag object
            vec   - the Vec object 
            Array - the read-only array
       
    """
    function DMStagVecGetArrayRead(dm::DMStag, v::AbstractVec)
        # Note: there is actually no need to call PETSc again, as Julia has the possibility 
        # to wrap an existing array into another one. Our vec already has the array wrapper, 
        # so we reshape that 

        # Extract array from vector. Note: we need to release this by calling 
        # finalize on X1!
        v.array     =   unsafe_localarray($PetscScalar, v.ptr;  write=false, read=true)

        X1          =   DMStagVecGetArray(dm, v.array) 
        
        return X1
    end

    function DMStagVecGetArray(dm::DMStag, v::Vector)

        entriesPerElement   =   DMStagGetEntriesPerElement(dm)
        nGhost              =   DMStagGetGhostCorners(dm)
        dim                 =   DMGetDimension(dm);         
      
        # Dimensions of new array (see the PETSc DMStagVecGetArrayRead routine)
        dim_vec             =   [entriesPerElement; collect(nGhost[2])];  

        # Wrap julia vector to new vector.
        X                   =    Base.view(v,:);
        
        # reshape to correct format
        X                   =   reshape(v, Tuple(dim_vec))
        X1                  =   PermutedDimsArray(X, Tuple([2:dim+1;1]));   # permute to take care of different array ordering in C/Julia
       
        return X1
    end


    

    """
        Julia routine that extracts an array related to a certain DOF. 
        Modifying values in the array will change them in the local PetscVec. Use LocalToGlobal to update global vector values

        Usage:

            Array = DMStagGetGhostArrayLocationSlot(dm::DMStag, v::AbstractVec              , loc::DMStagStencilLocation, dof::Int)
            Array = DMStagGetGhostArrayLocationSlot(dm::DMStag, ArrayFull::PermutedDimsArray, loc::DMStagStencilLocation, dof::Int)

        Input:
            dm           -   the DMStag object 
            v,ArrayFull  -   the local vector as obtained with DMCreateLocalVector, can also be a local array
            loc  	     -   a DMStagStencilLocation
            dof 	     -   the degree of freedom on loc, which you want to extracts
        
        Output:

            Array 	-   local array that includes the ghost points, that is linked to the vector v. 
                        Modifying values in Array will update v

    """
    function DMStagGetGhostArrayLocationSlot(dm::DMStag, v::AbstractVec{$PetscScalar}, loc::DMStagStencilLocation, dof::Int)
        entriesPerElement   =   DMStagGetEntriesPerElement(dm)
        dim                 =   DMGetDimension(dm);  
        slot                =   DMStagGetLocationSlot(dm, loc, dof); 
        slot_start          =   mod(slot,entriesPerElement);          # figure out which component we are interested in

        ArrayFull           =   DMStagVecGetArray(dm, v);             # obtain access to full array

        # now extract only the dimension belonging to the current point
        Array               =   selectdim(ArrayFull,dim+1, slot_start+1);

        return Array
    end

    function DMStagGetGhostArrayLocationSlot(dm::DMStag, ArrayFull::PermutedDimsArray, loc::DMStagStencilLocation, dof::Int)
        entriesPerElement   =   DMStagGetEntriesPerElement(dm)
        dim                 =   DMGetDimension(dm);  
        slot                =   DMStagGetLocationSlot(dm, loc, dof); 
        slot_start          =   mod(slot,entriesPerElement);          # figure out which component we are interested in

        # now extract only the dimension belonging to the current point
        Array               =   selectdim(ArrayFull,dim+1, slot_start+1);

        return Array
    end

    """
    get slot for use with local product coordinate arrays 

        slot = DMStagGetProductCoordinateLocationSlot(dm::DMStag,loc::DMStagStencilLocation)

            dm 	    - the DMStag object
            loc 	- the grid location 
            slot 	- the index to use in local arrays
    """
    function DMStagGetProductCoordinateLocationSlot(dm::DMStag,loc::DMStagStencilLocation)
        slot = Ref{$PetscInt}()
        @chk ccall((:DMStagGetProductCoordinateLocationSlot, $libpetsc), PetscErrorCode,
                    ( CDMStag,   DMStagStencilLocation, Ptr{$PetscInt}), dm, loc, slot)

        return slot[]
    end

    """
    get number of entries per element in the local representation 

        entriesPerElement = DMStagGetEntriesPerElement(dm::DMStag)

            dm 	              - the DMStag objects
            entriesPerElement - number of entries associated with each element in the local representation
    """
    function DMStagGetEntriesPerElement(dm::DMStag)
        entriesPerElement = Ref{$PetscInt}()
        @chk ccall((:DMStagGetEntriesPerElement, $libpetsc), PetscErrorCode,
                    ( CDMStag,  Ptr{$PetscInt}), dm,  entriesPerElement)

        return entriesPerElement[]
    end

    """
    get elementwise stencil width 

        stencilWidth = DMStagGetStencilWidth(dm::DMStag)

            dm 	         - the DMStag objects
            stencilWidth - stencil/halo/ghost width in elements
    """
    function DMStagGetStencilWidth(dm::DMStag)
        stencilWidth = Ref{$PetscInt}()
        @chk ccall((:DMStagGetStencilWidth, $libpetsc), PetscErrorCode,
                    ( CDMStag,  Ptr{$PetscInt}), dm,  stencilWidth)

        return stencilWidth[]
    end

    """
    get index to use in accessing raw local arrays 

        slot = DMStagGetLocationSlot(dm::DMStag,loc::DMStagStencilLocation, c)
        
            dm	    - the DMStag object
            loc	    - location relative to an element
            c	    - component ( the degree of freedom)
            slot	- index to use

    """
    function DMStagGetLocationSlot(dm::DMStag,loc::DMStagStencilLocation, c)
        
        slot = Ref{$PetscInt}()
        @chk ccall((:DMStagGetLocationSlot, $libpetsc), PetscErrorCode,
                    ( CDMStag,   DMStagStencilLocation, $PetscInt, Ptr{$PetscInt}), dm, loc, c, slot)

        return slot[]
    end

    """
    Destroys a DM

        destroy(dm::DMStag)

            dm 	- the DM object to destroy
    """
    function destroy(dm::DMStag{$PetscScalar})
        finalized($PetscScalar) ||
            @chk ccall((:DMDestroy, $libpetsc), PetscErrorCode, (Ptr{CDMStag},), dm)
        return nothing
    end

    """
    Gets the DM type name (as a string) from the DM

        type = gettype(dm::DMStag)

            dm 	 - The DM
            type - The DM type name 
    """
    function gettype(dm::DMStag{$PetscScalar})
        t_r = Ref{CDMStagType}()
        @chk ccall((:DMGetType, $libpetsc), PetscErrorCode, (CDMStag, Ptr{CDMStagType}), dm, t_r)
        return unsafe_string(t_r[])
    end

    """
    Views a DM

    view(dm::DMStag, viewer::Viewer{$PetscScalar}=ViewerStdout{$PetscScalar}(dm.comm))

            dm 	   - the DM object to view 
            viewer - the viewer 
    """
    function view(dm::DMStag{$PetscScalar}, viewer::Viewer{$PetscScalar}=ViewerStdout{$PetscScalar}(dm.comm))
        @chk ccall((:DMView, $libpetsc), PetscErrorCode, 
                    (CDMStag, CPetscViewer),
                dm, viewer);
        return nothing
    end

    """ 
    return global element indices of the local region (excluding ghost points) 
        
        x,m,nExtrax = DMStagGetCorners(dm:DMStag)   in 1D
        x[],m[],nExtrax[] = DMStagGetCorners(dm:DMStag)   in 2D or 3D

            dm 	    - the DMStag object
            x,y,z 	- starting element indices in each direction
            m,n,p 	- element widths in each direction
            nExtrax,nExtray,nExtraz 	- number of extra partial elements in each direction. 
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
                X = (x[],)
                M = (m[],)
                NEXTRA = (nExtrax[],)
                return X[1], M[1], NEXTRA[1]    
            elseif dm.dim==2
                return (x[], y[]), (m[],n[]), (nExtrax[],nExtray[])    
            elseif dm.dim==3
                return (x[], y[], z[]), (m[],n[],p[]), (nExtrax[],nExtray[],nExtraz[])    
            end
    end
    
    """ 
    return global element indices of the local region (including ghost points) 
        
        x,m,nExtrax = DMStagGetCorners(dm:DMStag)   in 1D
        x[],m[],nExtrax[] = DMStagGetCorners(dm:DMStag)   in 2D or 3D

            dm 	    - the DMStag object
            x,y,z 	- starting element indices in each direction
            m,n,p 	- element widths in each direction
            nExtrax,nExtray,nExtraz 	- number of extra partial elements in each direction. 
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
                X = (x[],)
                M = (m[],)
                NEXTRA = (nExtrax[],)
                return X[1], M[1], NEXTRA[1]
            elseif dm.dim==2
                return (x[], y[]), (m[],n[])    
            elseif dm.dim==3
                return (x[], y[], z[]), (m[],n[],p[])   
            end
    end

    """
    Return indices of start and end of the central nodes of a local array built from the input dm (excluding ghost nodes)    
    
        Cen_start, Cen_end = DMStagGetCentralNodes(dm::DMStag)

            dm 	                - the DMStag object
            Cen_start, Cen_end 	- indices of start and finish of central nodes
    """
    function DMStagGetCentralNodes(dm::DMStag)
        # in Julia, indices in arrays start @ 1, whereas they can go negative in C
        # This routine  

        g_start, g_N    =   DMStagGetGhostCorners(dm);
        g_width         =   DMStagGetStencilWidth(dm);
        start,N, nExtra =   DMStagGetCorners(dm);
        
        Cen_start       =   zeros(Int64,dm.dim)
        for i=1:length(g_start)
            Cen_start[i] = -g_start[i] + 1;
        end

        Cen_end         =   Cen_start  .+ N .- 1;
        return Cen_start, Cen_end
    end


    """
    get boundary types
        
        Bx = DMStagGetBoundaryTypes(dm::DMStag) in 1D
        Bx,By,Bz = DMStagGetBoundaryTypes(dm::DMStag) in 3D

            dm 	     - the DMStag object 
            Bx,By,Bz - boundary types
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

    """
    set elementwise stencil width in dm
        
        DMStagSetStencilWidth(dm::DMStag, stencilWidth::Int64)
        Bx,By,Bz = DMStagGetBoundaryTypes(dm::DMStag) in 3D

            dm 	         - the DMStag object 
            stencilWidth - stencil/halo/ghost width in elements 
    """
    function  DMStagSetStencilWidth(dm::DMStag, stencilWidth::Int64)

        @chk ccall((:DMStagSetStencilWidth, $libpetsc), PetscErrorCode,
             (CDMStag,  $PetscInt), dm, stencilWidth)

        return nothing
    end

    """
    get number of ranks in each direction in the global grid decomposition

        nRanks0 = DMStagGetNumRanks(dm::DMStag) in 1D
        nRanks0,nRanks1,nRanks2 = DMStagGetNumRanks(dm::DMStag) in 3D

            dm 	                     - the DMStag object 
            nRanks0,nRanks1,nRanks2  - number of ranks in each direction in the grid decomposition
    """
    function  DMStagGetNumRanks(dm::DMStag)

        nRanks0 = Ref{$PetscInt}()
        nRanks1 = Ref{$PetscInt}()
        nRanks2 = Ref{$PetscInt}()
        
        @chk ccall((:DMStagGetNumRanks, $libpetsc), PetscErrorCode,
            (CDMStag, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}), dm, nRanks0,nRanks1,nRanks2)
            
        if dm.dim==1
            return nRanks0[]
        elseif dm.dim==2
            return nRanks0[], nRanks1[]
        elseif dm.dim==3
            return nRanks0[], nRanks1[], nRanks2[]
        end
    end

    """
    set ranks in each direction in the global rank grid. Arguments corresponding to higher dimensions are ignored for 1D and 2D grids.  

        DMStagSetNumRanks(dm::DMStag, nRanks0::Int32, nRanks1=1::Int32, nRanks2=1::Int32)

            	dm 	                        - the DMStag object
	            nRanks0,nRanks1,nRanks2 	- number of ranks in each direction 
    """
    function  DMStagSetNumRanks(dm::DMStag, nRanks0::Int32, nRanks1=1::Int32, nRanks2=1::Int32)

        
        @chk ccall((:DMStagSetNumRanks, $libpetsc), PetscErrorCode,
            (CDMStag, $PetscInt, $PetscInt, $PetscInt), dm, nRanks0,nRanks1,nRanks2)
            
        return nothing
    end

    # NOT WORKING!
    function  DMStagGetOwnershipRanges(dm::DMStag, lx::Array{$PetscInt},  ly::Array{$PetscInt}, lz::Array{$PetscInt})

        nx = DMStagGetNumRanks(dm)

        
        @chk ccall((:DMStagGetOwnershipRanges, $libpetsc), PetscErrorCode,
            (CDMStag, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), dm, lx,ly,lz)


       # lx_vec     = unsafe_wrap(Vector{$PetscInt},lx[],  nx)

    #    return lx
     
    end


    """
        This puts a single value inside a global vector using DMStagStencil
        
            DMStagVecSetValueStencil(dm::DMStag,vec::Abstractvec, pos::DMStagStencil, val::Float64, insertMode::InsertMode)


                dm	- the DMStag object
                vec	- the Vec
                pos	- the location of the set values, given by a DMStagStencil struct
                val	- the value to be set
                insertMode	- INSERT_VALUES or ADD_VALUES
    """
    function  DMStagVecSetValueStencil(dm::DMStag, vec::AbstractVec{$PetscScalar}, pos::DMStagStencil, val, insertMode::InsertMode)

        n=1;
        @chk ccall((:DMStagVecSetValuesStencil, $libpetsc), PetscErrorCode,
                    (CDMStag, CVec, $PetscInt, Ptr{DMStagStencil_c}, Ptr{$PetscScalar}, InsertMode), 
                        dm, vec.ptr, n, Ref{DMStagStencil_c}(pos), Ref{$PetscScalar}(val), insertMode)

        return nothing
    end


    """
    get vector values using grid indexing (currently only valid for a single value)

        val = DMStagVecGetValueStencil(dm::DMStag, vec::AbstractVec{$PetscScalar}, pos::DMStagStencil)

            dm 	- the DMStag object
            vec - the vector object
            n 	- the number of values to obtain
            pos - locations to obtain values from (as an array of DMStagStencil values) 
            val - value at the point 
    """
    function  DMStagVecGetValueStencil(dm::DMStag, vec::AbstractVec{$PetscScalar}, pos::DMStagStencil)

        n=1;
        val = Ref{$PetscScalar}()
        @chk ccall((:DMStagVecGetValuesStencil, $libpetsc), PetscErrorCode,
                    (CDMStag, CVec, $PetscInt, Ptr{DMStagStencil_c}, Ptr{$PetscScalar}), 
                        dm, vec.ptr, n, Ref{DMStagStencil_c}(pos), val)
    
        return val[]
    end

    """
    This reads a single value from a matrix DMStagStencil
        
        val =  DMStagMatGetValueStencil(dm::DMStag,mat::AbstractMat, posRow::DMStagStencil,  posCol::DMStagStencil)

            dm	    - the DMStag object
            mat	    - the Mat
            posRow	- the location of the row of the set value, given by a DMStagStencil struct
            posCol	- the location of the row of the set value, given by a DMStagStencil struct
            val	    - the value
    """
    function  DMStagMatGetValueStencil(dm::DMStag, mat::AbstractMat{$PetscScalar},  posRow::DMStagStencil, posCol::DMStagStencil)

        nRow= 1;
        nCol= 1;
        val = Ref{$PetscScalar}()
        @chk  ccall((:DMStagMatGetValuesStencil, $libpetsc), PetscErrorCode,
                    (CDMStag, CMat, $PetscInt, Ptr{DMStagStencil_c}, $PetscInt, Ptr{DMStagStencil_c}, Ptr{$PetscScalar}), 
                        dm, mat.ptr, nRow, Ref{DMStagStencil_c}(posRow), nCol, Ref{DMStagStencil_c}(posCol), val)
    
        return val[]
    end

    """ 
        Give the non-ghosted indices in the local vector that contribute to the global vector

            indices = LocalInGlobalIndices(dm::DMStag)

                dm	    - the DMStag object
                indices - local indices
    """
    function LocalInGlobalIndices(dm::DMStag)
        # note: this can likely be done more efficiently and will have to be modified in parallel
        ind_g   =   DMCreateGlobalVector(dm)
        v_ind_l =   DMCreateLocalVector(dm)

        ind_l   = unsafe_localarray(Float64, v_ind_l.ptr);
        for i=1:length(ind_l)
            ind_l[i] = i
        end
        
        DMLocalToGlobal(dm,v_ind_l, INSERT_VALUES, ind_g);

        return Int64.(ind_g.array)

    end
    
    """
    This puts a single value inside a matrix using DMStagStencil position
    
        DMStagMatSetValueStencil(dm::DMStag,mat::AbstractMat, posRow::DMStagStencil,  posCol::DMStagStencil, val::Float64, insertMode::InsertMode)

            dm	        - the DMStag object
            mat	        - the Mat
            posRow	    - the location of the row of the set value, given by a DMStagStencil struct
            posCol	    - the location of the row of the set value, given by a DMStagStencil struct
            val	        - the value to be set
            insertMode	- INSERT_VALUES or ADD_VALUES
    """
    function  DMStagMatSetValueStencil(dm::DMStag, mat::AbstractMat{$PetscScalar},  posRow::DMStagStencil, posCol::DMStagStencil, val, insertMode::InsertMode)

        nRow= 1;
        nCol= 1;
        @chk ccall((:DMStagMatSetValuesStencil, $libpetsc), PetscErrorCode,
                    (CDMStag, CMat, $PetscInt, Ptr{DMStagStencil_c},  $PetscInt, Ptr{DMStagStencil_c}, Ptr{$PetscScalar}, InsertMode), 
                        dm, mat.ptr, nRow, Ref{DMStagStencil_c}(posRow), nCol, Ref{$DMStagStencil_c}(posCol), Ref{$PetscScalar}(val), insertMode)

        return nothing
    end



  
    """
    Return the topological dimension of the DM

        dim = DMGetDimension(dm::DMStag)

            dm 	- The DM
            dim - dimensions
    """
    function DMGetDimension(dm::DMStag)
        dim = Ref{$PetscInt}()

        @chk ccall((:DMGetDimension, $libpetsc), PetscErrorCode, (CDMStag,Ptr{$PetscInt}), dm, dim )

        return dim[]
    end


    """
    updates global vectors from local vectors 

        DMLocalToGlobal(dm::DMStag,l::AbstractVec{$PetscScalar}, mode::InsertMode,g::AbstractVec{$PetscScalar})
        DMLocalToGlobal(dm::DMStag,l::AbstractVec{$PetscScalar}, mode::InsertMode,g::CVec)

            dm 	 - the DM object
	        l 	 - the local vector
	        mode - if INSERT_VALUES then no parallel communication is used, if ADD_VALUES then all ghost points from the same base point accumulate into that base point.
	        g 	 - the global vector 
    """
    function DMLocalToGlobal(dm::DMStag,l::AbstractVec{$PetscScalar}, mode::InsertMode,g::AbstractVec{$PetscScalar})

        DMLocalToGlobal(dm,l.ptr, mode,g.ptr)

        return nothing
    end
    function DMLocalToGlobal(dm::DMStag,l::AbstractVec{$PetscScalar}, mode::InsertMode,g::CVec)

        DMLocalToGlobal(dm,l.ptr, mode,g)

        return nothing
    end
    function DMLocalToGlobal(dm::DMStag,l::CVec, mode::InsertMode,g::CVec)

        @chk ccall((:DMLocalToGlobal, $libpetsc), PetscErrorCode,
        (CDMStag, CVec, InsertMode, CVec), 
            dm, l, mode, g)

        return nothing
    end

    
    """
    update local vectors from global vector 

        DMGlobalToLocal(dm::DMStag,g::AbstractVec{$PetscScalar}, mode::InsertMode,l::AbstractVec{$PetscScalar})
        DMGlobalToLocal(dm::DMStag,g::CVec, mode::InsertMode,l::AbstractVec{$PetscScalar})

            dm 	 - the DM object
	        g 	 - the global vector
	        mode - INSERT_VALUES or ADD_VALUES
	        l 	 - the local vector 
    """
    function DMGlobalToLocal(dm::DMStag,g::AbstractVec{$PetscScalar}, mode::InsertMode,l::AbstractVec{$PetscScalar})

        DMGlobalToLocal(dm,g.ptr, mode::InsertMode,l.ptr)

        return nothing
    end

    function DMGlobalToLocal(dm::DMStag,g::CVec, mode::InsertMode,l::AbstractVec{$PetscScalar})

        DMGlobalToLocal(dm,g, mode::InsertMode,l.ptr)

        return nothing
    end

    function DMGlobalToLocal(dm::DMStag,g::CVec, mode::InsertMode,l::CVec)

        @chk ccall((:DMGlobalToLocal, $libpetsc), PetscErrorCode,
        (CDMStag, CVec, InsertMode, CVec), 
            dm, g, mode, l)

        return nothing
    end


    """
    Generates a matrix from a DMStag object. The type is a MatSeqAIJ is we are on 1 core 

        mat = DMCreateMatrix(dm::DMStag)

            dm	- the DMStag object
            mat	- the matrix of type MatSeqAIJ (on 1 core) 
    """
    function DMCreateMatrix(dm::DMStag)
        # Note: the matrix cannot be viewed yet, as it remains unassembled
        #  ideally, we should modify the viewer to take care of this case

        if dm.comm==MPI.COMM_SELF
            mat = MatSeqAIJ{$PetscScalar}(C_NULL, dm.comm)
        elseif dm.comm==MPI.COMM_WORLD
            error("MatMPIAIJ still to be implemented")
        end

        @chk ccall((:DMCreateMatrix, $libpetsc), PetscErrorCode,
                    (CDMStag, Ptr{CMat}), dm, mat)
        
        return mat
    end

    """
    get elementwise ghost/halo stencil type

        stencilType = DMStagGetStencilType(dm::DMStag)

            dm 	        - the DMStag object
            stencilType - the elementwise ghost stencil type: DMSTAG_STENCIL_BOX, DMSTAG_STENCIL_STAR, or DMSTAG_STENCIL_NONE
    """
    function DMStagGetStencilType(dm::DMStag)
        stencilType =  Ref{DMStagStencilType}()

        @chk ccall((:DMStagGetStencilType, $libpetsc), PetscErrorCode, (CDMStag, Ptr{DMStagStencilType}), dm, stencilType)

        return stencilType[]
    end

    """
    get boolean value for whether this rank is first in each direction in the rank grid. Arguments corresponding to higher dimensions are ignored for 1D and 2D grids. These arguments may be set to NULL in this case.

        fr_X,fr_Y,fr_Z = DMStagGetIsFirstRank(dm::DMStag)

            dm 	           - the DMStag object
            fr_X,fr_Y,fr_Z - whether this rank is first in each direction
    """
    function DMStagGetIsFirstRank(dm::DMStag)
        fr_X = Ref{PetscBool}()
        fr_Y = Ref{PetscBool}()
        fr_Z = Ref{PetscBool}()
        
        @chk ccall((:DMStagGetIsFirstRank, $libpetsc), PetscErrorCode, (CDMStag, Ptr{PetscBool}, Ptr{PetscBool}, Ptr{PetscBool}), dm, fr_X, fr_Y, fr_Z)
        
        return fr_X[]== PETSC_TRUE, fr_Y[]== PETSC_TRUE, fr_Z[]== PETSC_TRUE
    end



    """
    get boolean value for whether this rank is last in each direction in the rank grid 

        fr_X,fr_Y,fr_Z = DMStagGetIsLastRank(dm::DMStag)

            dm 	           - the DMStag object
            fr_X,fr_Y,fr_Z - whether this rank is last in each direction
    """
    function DMStagGetIsLastRank(dm::DMStag)
        fr_X = Ref{PetscBool}()
        fr_Y = Ref{PetscBool}()
        fr_Z = Ref{PetscBool}()
        
        @chk ccall((:DMStagGetIsLastRank, $libpetsc), PetscErrorCode, (CDMStag, Ptr{PetscBool}, Ptr{PetscBool}, Ptr{PetscBool}), dm, fr_X, fr_Y, fr_Z)
        
        return fr_X[]== PETSC_TRUE, fr_Y[]== PETSC_TRUE, fr_Z[]== PETSC_TRUE
    end

    # NOT WORKING YET!
    function DMStagGetProductCoordinateArrays(dm::DMStag)

        #  arrX = Ref{Ptr{$PetscScalar}}()
          arrY = Ref{Ptr{$PetscScalar}}()
          arrZ = Ref{Ptr{$PetscScalar}}()
  
          #arrX = zeros(20,2)
  
          array_ref = Ref{Ptr{$PetscScalar}}()
  
          @chk ccall((:DMStagGetProductCoordinateArrays, $libpetsc), PetscErrorCode,
              ( CDMStag,   Ref{Ptr{$PetscScalar}}, Ptr{Ptr{$PetscScalar}}, Ptr{Ptr{$PetscScalar}}), 
                  dm, array_ref, arrY, arrZ)
       
          arrX = unsafe_wrap(Array, array_ref[], 20)
        #same  r_sz = Ref{$PetscInt}()
          
          #@chk ccall((:VecGetLocalSize, $libpetsc), PetscErrorCode,
          #            (CVec, Ptr{$PetscInt}), arrX, r_sz)
        #  x_local = PETSc.unsafe_localarray(Float64, arrX)
  
          #v = unsafe_wrap(Array, arrX[], r_sz[]; own = false)
  
          #x = unsafe_localarray($PetscScalar, arrX; write=false)
          
          return arrX
  
      end

    """
    Gets the DM that prescribes coordinate layout and scatters between global and local coordinates

        dmnew = DMGetCoordinateDM(dm::DMStag; kwargs...)

            dm 	           - the DMStag object
            kwargs...      -   [Optional] keyword arguments (see PETSc webpage), specifiable as stag_grid_x=100, etc. 
            dmnew          - Coordinate dm
    """
      function DMGetCoordinateDM(dm::DMStag; kwargs...)

        comm  = dm.comm

        dim   = DMGetDimension(dm)

        opts  = Options{$PetscScalar}(kwargs...)

        dmnew = DMStag{$PetscScalar}(C_NULL, comm, dim, opts)

        @chk ccall((:DMGetCoordinateDM, $libpetsc), PetscErrorCode,
        (CDMStag, Ptr{CDMStag}), 
         dm, dmnew)

        return dmnew
    end

    """
    Gets a local vector with the coordinates associated with the DM.

        v = DMGetCoordinatesLocal(dm::DMStag; write_val=true, read_val=true)

            dm 	- the DMStag object
            v   - coordinate local vector
    """
    function DMGetCoordinatesLocal(dm::DMStag; write_val=true, read_val=true)

        v = VecSeq(C_NULL, dm.comm, [0.0])  # empty vector
        
        ccall((:DMGetCoordinatesLocal, $libpetsc), PetscErrorCode, (CDMStag, Ptr{CVec}), dm, v)

        # Link a julia array to the values from the new vector
        # If we modify values here, it will automatically be changed in the PetcVec as well
        v.array = unsafe_localarray($PetscScalar, v.ptr; write=write_val, read=read_val)
        
        return v
    end
 

end

Base.show(io::IO, dm::DMStag) = _show(io, dm)
