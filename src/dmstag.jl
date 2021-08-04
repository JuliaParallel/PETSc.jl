# Attempt at include dmstag functions
const CDMStag = Ptr{Cvoid}
const CDMStagType = Cstring

mutable struct DMStag{PetscLib} <: AbstractDM{PetscLib}
    ptr::CDMStag
    opts::Options{PetscLib}
    
    DMStag{PetscLib}(ptr, opts = Options(PetscLib)) where {PetscLib} =
    new{PetscLib}(ptr, opts)
end

"""
    empty(dm::DMStag)

return an uninitialized `DMStag` struct.
"""
Base.empty(::DMStag{PetscLib}) where {PetscLib} = DMStag{PetscLib}(C_NULL)


mutable struct DMSTAGSTENCIL{PetscInt}
    loc::DMStagStencilLocation
    i::PetscInt
    j::PetscInt
    k::PetscInt
    c::PetscInt
end

const DMStagStencil     = DMSTAGSTENCIL

# allows us to pass XXMat objects directly into CMat ccall signatures
Base.cconvert(::Type{CDMStag}, obj::DMStag) = obj.ptr
# allows us to pass XXMat objects directly into Ptr{CMat} ccall signatures
Base.unsafe_convert(::Type{Ptr{CDMStag}}, obj::DMStag) =
    convert(Ptr{CDMStag}, pointer_from_objref(obj))

"""
    dm = DMStagCreate1d(::PetscLib,
        comm::MPI.Comm, 
        bndx::DMBoundaryType, 
        M, 
        dofVertex, 
        dofCenter, 
        stencilType::DMStagStencilType=DMSTAG_STENCIL_BOX, 
        stencilWidth=2, 
        lx=C_NULL; 
        dmsetfromoptions=true,
        dmsetup=true,
        options...
        )

Creates a 1D DMStag object.
        ::PetscLib      -   PETSc library,
        comm            -   MPI communicator
        bndx            -   boundary type: DM_BOUNDARY_NONE, DM_BOUNDARY_PERIODIC, or DM_BOUNDARY_GHOSTED. 
        M               -   global number of grid points
        dofVertex       -   [=1] number of degrees of freedom per vertex/point/node/0-cell
        dofCenter       -   [=1] number of degrees of freedom per element/edge/1-cell
        stencilType     -   ghost/halo region type: DMSTAG_STENCIL_BOX or DMSTAG_STENCIL_NONE
        stencilWidth    -   width, in elements, of halo/ghost region
        lx              -   [Optional] Vector of local sizes, of length equal to the comm size, summing to M
        options...      -   [Optional] keyword arguments (see PETSc webpage), specifiable as stag_grid_x=100, etc. 

Creates a 1-D distributed staggered array with the options specified using keyword
arguments.

If keyword argument `dmsetfromoptions == true` then `setfromoptions!` called.
If keyword argument `dmsetup == true` then `setup!` is called.

# External Links
$(_doc_external("DMSTAG/DMStagCreate1d"))

"""
function DMStagCreate1d end

@for_petsc function DMStagCreate1d(
    ::$UnionPetscLib,
    comm::MPI.Comm, 
    bndx::DMBoundaryType, 
    M, 
    dofVertex=1,
    dofCenter=1,
    stencilType::DMStagStencilType=DMSTAG_STENCIL_BOX,
    stencilWidth=2, 
    lx=C_NULL;
    dmsetfromoptions = true,
    dmsetup = true, 
    options...,
)

    if isempty(lx); lx = C_NULL; end
    opts = Options($petsclib, options...)
    dm   = DMStag{$PetscLib}(C_NULL, opts)   # retrieve options
    with(dm.opts) do

        @chk ccall(
                (:DMStagCreate1d, $petsc_library), 
                PetscErrorCode,
                (
                    MPI.MPI_Comm, 
                    DMBoundaryType, 
                    $PetscInt, 
                    $PetscInt, 
                    $PetscInt, 
                    DMStagStencilType, 
                    $PetscInt,  
                    Ptr{$PetscInt}, 
                    Ptr{CDMStag}
                ),
                comm, 
                bndx, 
                M,
                dofVertex,
                dofCenter,
                stencilType,
                stencilWidth,
                lx, 
                dm 
                )
    end
    dmsetfromoptions && setfromoptions!(dm)
    dmsetup && setup!(dm)

    if comm == MPI.COMM_SELF
        finalizer(destroy, dm)
    end
        
    return dm
end

"""
    dm = DMStagCreate2d(
        ::PetscLib,
        comm::MPI.Comm, 
        bndx::DMBoundaryType, 
        bndy::DMBoundaryType, 
        M, N, 
        m=C_NULL, n=C_NULL, 
        dofVertex=1, 
        dofEdge=1, 
        dofElement=1, 
        stencilType::DMStagStencilType=DMSTAG_STENCIL_BOX, 
        stencilWidth=2, 
        lx=C_NULL, ly=C_NULL; 
        dmsetfromoptions=true,
        dmsetup=true,
        options...
        )

Creates a 2D DMStag object.

If keyword argument `dmsetfromoptions == true` then `setfromoptions!` called.
If keyword argument `dmsetup == true` then `setup!` is called.

# External Links
$(_doc_external("DMSTAG/DMStagCreate2d"))

"""
function DMStagCreate2d end

@for_petsc function DMStagCreate2d(
    ::$UnionPetscLib,
    comm::MPI.Comm, 
    bndx::DMBoundaryType, 
    bndy::DMBoundaryType, 
    M, N, 
    m=C_NULL, n=C_NULL, 
    dofVertex=1, 
    dofEdge=1, 
    dofElement=1, 
    stencilType::DMStagStencilType=DMSTAG_STENCIL_BOX, 
    stencilWidth=2, 
    lx=C_NULL, ly=C_NULL; 
    dmsetfromoptions=true,
    dmsetup=true,
    options...,
)
        
    if isempty(lx); lx = C_NULL; end
    if isempty(ly); ly = C_NULL; end
    opts = Options($petsclib, options...)

    dm = DMStag{$PetscLib}(C_NULL, opts)

    with(dm.opts) do
        @chk ccall(
            (:DMStagCreate2d, $petsc_library), 
            PetscErrorCode,
                (
                    MPI.MPI_Comm, 
                    DMBoundaryType, DMBoundaryType, 
                    $PetscInt, $PetscInt, 
                    $PetscInt, $PetscInt, 
                    $PetscInt, 
                    $PetscInt, 
                    $PetscInt, 
                    DMStagStencilType, 
                    $PetscInt, 
                    Ptr{$PetscInt},  Ptr{$PetscInt}, 
                    Ptr{CDMStag}
                ),
                comm, 
                bndx, bndy, 
                M, N, 
                m, n ,
                dofVertex ,
                dofEdge ,
                dofElement ,
                stencilType ,
                stencilWidth ,
                lx ,ly ,
                dm 
            )
    end
    dmsetfromoptions && setfromoptions!(dm)
    dmsetup && setup!(dm)

    if comm == MPI.COMM_SELF
        finalizer(destroy, dm)
    end
        
    return dm
end

"""

    dm = DMStagCreate3d(
        ::PetscLib,
        comm::MPI.Comm, 
        bndx::DMBoundaryType, bndy::DMBoundaryType, bndz::DMBoundaryType, 
        M, N, P, 
        m, n, p, 
        dofVertex, 
        dofEdge, 
        dofFace, 
        dofElement, 
        stencilType::DMStagStencilType=DMSTAG_STENCIL_BOX, 
        stencilWidth, 
        lx, ly, lz; 
        dmsetfromoptions=true,
        dmsetup=true,
        options...
        )

Creates a 3D DMStag object.

If keyword argument `dmsetfromoptions == true` then `setfromoptions!` called.
If keyword argument `dmsetup == true` then `setup!` is called.

# External Links
$(_doc_external("DMSTAG/DMStagCreate3d"))
"""
function DMStagCreate3d end

@for_petsc function DMStagCreate3d(
    ::$UnionPetscLib,
    comm::MPI.Comm, 
    bndx::DMBoundaryType, bndy::DMBoundaryType, bndz::DMBoundaryType, 
    M, N, P, 
    m=C_NULL, n=C_NULL, p=C_NULL, 
    dofVertex=1, 
    dofEdge=1, 
    dofFace=1, 
    dofElement=1, 
    stencilType::DMStagStencilType=DMSTAG_STENCIL_BOX, 
    stencilWidth=2, 
    lx=C_NULL, ly=C_NULL, lz=C_NULL; 
    dmsetfromoptions=true,
    dmsetup=true,
    options...,
    )
        
    if isempty(lx); lx = C_NULL; end
    if isempty(ly); ly = C_NULL; end
    if isempty(lz); lz = C_NULL; end
    opts = Options($petsclib, options...)
        
    dm = DMStag{$PetscLib}(C_NULL, opts)

    with(dm.opts) do
        @chk ccall((:DMStagCreate3d, $petsc_library), PetscErrorCode,
                (
                    MPI.MPI_Comm, 
                    DMBoundaryType, DMBoundaryType, DMBoundaryType, 
                    $PetscInt, $PetscInt, $PetscInt, 
                    $PetscInt, $PetscInt, $PetscInt, 
                    $PetscInt, 
                    $PetscInt, 
                    $PetscInt, 
                    $PetscInt, 
                    DMStagStencilType, 
                    $PetscInt, 
                    Ptr{$PetscInt},  Ptr{$PetscInt},  Ptr{$PetscInt}, 
                    Ptr{CDMStag}
                ),
                comm, 
                bndx, bndy, bndz, 
                M, N, P, 
                m, n ,p ,
                dofVertex ,
                dofEdge ,
                dofFace ,
                dofElement ,
                stencilType ,
                stencilWidth ,
                lx ,ly ,lz ,
                dm 
                )
    end
    dmsetfromoptions && setfromoptions!(dm)
    dmsetup && setup!(dm)

    # We can only let the garbage collect finalize when we do not need to
    # worry about MPI (since garbage collection is asyncronous)
    if comm == MPI.COMM_SELF
        finalizer(destroy, dm)
    end
        
    return dm
end

    
"""
    dm = DMStagCreateCompatibleDMStag(
        dm::DMStag, 
        dofVertex, 
        dofEdge, 
        dofFace, 
        dofElement; 
        kwargs...)

Creates a compatible DMStag with different dof/stratum 

        dm              -   the DMStag object 
        dofVertex       -   [=0] number of degrees of freedom per vertex/point/node/0-cell
        dofEdge         -   [=0] number of degrees of freedom per edge/1-cell 
        dofFace         -   [=0] number of degrees of freedom per face/2-cell 
        dofElement      -   [=0] number of degrees of freedom per element/3-cell 
        kwargs...       -   [Optional] keyword arguments (see PETSc webpage), specifiable as stag_grid_x=100, etc. 

# External Links
$(_doc_external("DMSTAG/DMStagCreateCompatibleDMStag"))
"""
function DMStagCreateCompatibleDMStag end

@for_petsc function DMStagCreateCompatibleDMStag(
    dm::DMStag{$PetscLib}, 
    dofVertex=0, 
    dofEdge=0, 
    dofFace=0, 
    dofElement=0; 
    dmsetfromoptions=true,
    dmsetup=true,
    options...
    )

    opts = Options($petsclib, options...)

    dmnew = DMStag{$PetscLib}(C_NULL, opts)
    comm  = getcomm(dm);

    with(dm.opts) do
        @chk ccall((:DMStagCreateCompatibleDMStag, $petsc_library), PetscErrorCode, 
        (
            CDMStag, 
            $PetscInt, 
            $PetscInt, 
            $PetscInt, 
            $PetscInt, 
            Ptr{CDMStag}
        ), 
        dm, 
        dofVertex, 
        dofEdge, 
        dofFace, 
        dofElement, 
        dmnew
        )
    end

    dmsetfromoptions && setfromoptions!(dm)
    dmsetup && setup!(dm)

    if comm == MPI.COMM_SELF
        finalizer(destroy, dmnew)
    end
        
    return dmnew

end

"""
    dof0,dof1,dof2,dof3 = DMStagGetDOF(dm::DMStag)
        
Get number of DOF associated with each stratum of the grid. 
        
    dm      - the DMStag object 
    dof0 	- the number of points per 0-cell (vertex/node)
    dof1 	- the number of points per 1-cell (element in 1D, edge in 2D and 3D)
    dof2 	- the number of points per 2-cell (element in 2D, face in 3D)
    dof3 	- the number of points per 3-cell (element in 3D)     

# External Links
$(_doc_external("DMSTAG/DMStagGetDOF"))
"""
function DMStagGetDOF end


@for_petsc function DMStagGetDOF(dm::DMStag{$PetscLib})

    dof0 = Ref{$PetscInt}()
    dof1 = Ref{$PetscInt}()
    dof2 = Ref{$PetscInt}()
    dof3 = Ref{$PetscInt}()

    @chk ccall((:DMStagGetDOF, $petsc_library), PetscErrorCode, 
    (
        CDMStag, 
        Ptr{$PetscInt}, 
        Ptr{$PetscInt}, 
        Ptr{$PetscInt}, 
        Ptr{$PetscInt}), 
    dm, 
    dof0, 
    dof1, 
    dof2, 
    dof3
    )

    dim   = getdimension(dm)

    if dim==1
        return  dof0[],dof1[]   
    elseif dim==2
        return dof0[],dof1[],dof2[]
    elseif dim==3
        return dof0[],dof1[],dof2[],dof3[]
    end

end


"""
    M,N,P = DMStagGetGlobalSizes(dm::DMStag)

Gets the global size of the DMStag object

    dm      - the DMStag object 
    M,N,P   - size in x,y,z

# External Links
$(_doc_external("DMSTAG/DMStagGetGlobalSizes"))
"""
function DMStagGetGlobalSizes end

@for_petsc function DMStagGetGlobalSizes(dm::DMStag{$PetscLib})

    M = Ref{$PetscInt}()
    N = Ref{$PetscInt}()
    P = Ref{$PetscInt}()

    @chk ccall((:DMStagGetGlobalSizes, $petsc_library), PetscErrorCode,
        (
            CDMStag, 
            Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}
        ), 
        dm,
        M, N, P 
        )

        dim   = getdimension(dm)
        
    if dim==1    
        return M[]
    elseif dim==2
        return M[], N[] 
    elseif dim==3
        return M[], N[], P[]    
    end
end


@for_petsc function Base.size(dm::DMStag{$PetscLib})
    size = DMStagGetGlobalSizes(dm)
    return size
end


"""
    M,N,P = DMStagGetLocalSizes(dm::DMStag)

Gets the local size of the DMStag object

    dm      - the DMStag object 
    M,N,P   - size in x,y,z

# External Links
$(_doc_external("DMSTAG/DMStagGetLocalSizes"))
"""
function DMStagGetLocalSizes end

@for_petsc function DMStagGetLocalSizes(dm::DMStag{$PetscLib})

    M = Ref{$PetscInt}()
    N = Ref{$PetscInt}()
    P = Ref{$PetscInt}()

    @chk ccall((:DMStagGetLocalSizes, $petsc_library), PetscErrorCode,
        (
            CDMStag, 
            Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}
        ), 
        dm, 
        M, N, P 
        )

        dim   = getdimension(dm)
        
    if dim==1    
        return M[]
    elseif dim==2
        return M[], N[] 
    elseif dim==3
        return M[], N[], P[]    
    end
end


"""
    setuniformcoordinatesproduct!(
        dm::DMStag,
        xyzmin::NTuple{N, Real},
        xyzmax::NTuple{N, Real},
        )

Set uniform coordinates for the `dmstag` using the lower and upper corners defined
    by the `NTuple`s `xyzmin` and `xyzmax`. If `N` is less than the dimension of the
    `dmstag` then the value of the trailing coordinates is set to `0`.
    
# External Links
$(_doc_external("DMSTAG/DMStagSetUniformCoordinatesProduct"))
    
"""
function setuniformcoordinatesproduct! end

@for_petsc function setuniformcoordinatesproduct!(
    dm::DMStag{$PetscLib},
    xyzmin::NTuple{N, Real},
    xyzmax::NTuple{N, Real},
    ) where {N}

    xmin = $PetscReal(xyzmin[1])
    xmax = $PetscReal(xyzmax[1])

    ymin = (N > 1) ? $PetscReal(xyzmin[2]) : $PetscReal(0)
    ymax = (N > 1) ? $PetscReal(xyzmax[2]) : $PetscReal(0)

    zmin = (N > 2) ? $PetscReal(xyzmin[3]) : $PetscReal(0)
    zmax = (N > 2) ? $PetscReal(xyzmax[3]) : $PetscReal(0)

    @chk ccall((:DMStagSetUniformCoordinatesProduct, $petsc_library), PetscErrorCode,
                ( 
                    CDMStag,   
                    $PetscReal,
                    $PetscReal,
                    $PetscReal,
                    $PetscReal,
                    $PetscReal,
                    $PetscReal,
                ), 
                dm, 
                xmin, xmax, 
                ymin, ymax, 
                zmin, zmax)

    return nothing
end

# NEED TO BE REPAIRED
@for_petsc function DMStagGetProductCoordinateArraysRead(dm::DMStag{$PetscLib})

    Arrx = Ref{$PetscScalar}()
    Arry = Ref{$PetscScalar}()
    Arrz = Ref{$PetscScalar}()

    @chk ccall((:DMStagGetProductCoordinateArraysRead, $petsc_library), PetscErrorCode, 
        (
            CDMStag, 
            Ptr{$PetscScalar}, Ptr{$PetscScalar}, Ptr{$PetscScalar}
        ), 
        dm, 
        Arrx, Arry, Arrz
        )
        
    return Arrx[],Arry[],Arrz[]
end


# TO BE REMOVED, AS setuniformcoordinates! DOES THE SAME?
"""
    DMStagSetUniformCoordinatesExplicit(
        dm::DMStag, 
        xmin, xmax, 
        ymin=0, ymax=0, 
        zmin=0, zmax=0
        )

Set DMStag coordinates to be a uniform grid, storing all values.

    dm                            - the DMStag object
    xmin,xmax,ymin,ymax,zmin,zmax - maximum and minimum global coordinate values 

# External Links
$(_doc_external("DMSTAG/DMStagSetUniformCoordinatesExplicit"))
"""
function DMStagSetUniformCoordinatesExplicit end

@for_petsc function DMStagSetUniformCoordinatesExplicit(
    dm::DMStag{$PetscLib}, 
    xmin, xmax, 
    ymin=0, ymax=0, 
    zmin=0, zmax=0
    )
        
    @chk ccall((:DMStagSetUniformCoordinatesExplicit, $petsc_library), PetscErrorCode,
                ( 
                    CDMStag,   
                    $PetscScalar, $PetscScalar, 
                    $PetscScalar, $PetscScalar, 
                    $PetscScalar, $PetscScalar
                ), 
                dm, 
                xmin, xmax, 
                ymin, ymax, 
                zmin, zmax
                )

    return nothing
end


"""
    Array =  DMStagVecGetArray(dm::DMStag, v::AbstractVec)
    
Get access to local array (including ghost points) of the DMStag.

    dm 	  - the DMStag object
    vec   - the Vec object 
    Array - the array
    
Once you are done with work on the array, you MUST release the memory with
    
    Base.finalize(Array)
    
Otherwise the values are not returned correctly to v    
"""
function DMStagVecGetArray end

@for_petsc function DMStagVecGetArray(dm::DMStag{$PetscLib}, v::AbstractVec)
    # Note: there is actually no need to call PETSc again, as Julia has the possibility 
    # to wrap an existing array into another one. Our vec already has the array wrapper, 
    # so we reshape that 

    # Extract array from vector. Note: we need to release this by calling 
    # Base.finalize on X1!
    array       =   unsafe_localarray($PetscScalar, v.ptr;  write=true, read=true)

    X1          =   DMStagVecGetArray(dm, array) 
        
    return X1
end

"""
    Array =  DMStagVecGetArrayRead(dm::DMStag, v::AbstractVec)
    
Get read-only access to a local array (including ghost points) of the DMStag

    dm    - the DMStag object
    vec   - the Vec object 
    Array - the read-only array
"""
function DMStagVecGetArrayRead end

@for_petsc function DMStagVecGetArrayRead(dm::DMStag{$PetscLib}, v::AbstractVec)
    # Note: there is actually no need to call PETSc again, as Julia has the possibility 
    # to wrap an existing array into another one. Our vec already has the array wrapper, 
    # so we reshape that 

    # Extract array from vector. Note: we need to release this by calling 
    # finalize on X1!
    array       =   unsafe_localarray($PetscScalar, v.ptr;  write=false, read=true)

    X1          =   DMStagVecGetArray(dm, array) 
        
    return X1
end

"""
    X1 = DMStagVecGetArray(dm::DMStag, v::Vector)

Returns a julia array from a vector `v`, in the same shape as the DMSTAG, which can be used to set values.
"""
function DMStagVecGetArray end

@for_petsc function DMStagVecGetArray(dm::DMStag{$PetscLib}, v::Vector)

    entriesPerElement   =   DMStagGetEntriesPerElement(dm)
    ghost_corners       =   getghostcorners(dm);
    dim                 =   getdimension(dm);         

    # Dimensions of new array (see the PETSc DMStagVecGetArrayRead routine)
    dim_vec             =   [entriesPerElement; collect(ghost_corners.size[1:dim])];  

    # Wrap julia vector to new vector.
    X                   =    Base.view(v,:);
        
    # reshape to correct format
    X                   =   reshape(v, Tuple(dim_vec))
    X1                  =   PermutedDimsArray(X, Tuple([2:dim+1;1]));   # permute to take care of different array ordering in C/Julia
       
    return X1
end  

"""
    Array = DMStagGetGhostArrayLocationSlot(
        dm::DMStag, 
        v::AbstractVec, 
        loc::DMStagStencilLocation, 
        dof::Int
        )
    
    
Julia routine that extracts an array related to a certain DOF. Modifying values in the array will change them in the local PetscVec. Use LocalToGlobal to update global vector values.

Usage:
            
    Input:
        dm           -   the DMStag object 
        v,ArrayFull  -   the local vector as obtained with DMCreateLocalVector, can also be a local array
        loc          -   a DMStagStencilLocation
        dof          -   the degree of freedom on loc, which you want to extracts
        
    Output:
        Array       -   local array that includes the ghost points, that is linked to the vector `v`. 
                            Modifying values in Array will update `v`

"""
function DMStagGetGhostArrayLocationSlot end

@for_petsc function DMStagGetGhostArrayLocationSlot(
    dm::DMStag{$PetscLib}, 
    v::AbstractVec{$PetscScalar}, 
    loc::DMStagStencilLocation, 
    dof::Int
    )

    entriesPerElement   =   DMStagGetEntriesPerElement(dm)
    dim                 =   getdimension(dm);  
    slot                =   DMStagGetLocationSlot(dm, loc, dof); 
    slot_start          =   mod(slot,entriesPerElement);          # figure out which component we are interested in

    ArrayFull           =   DMStagVecGetArray(dm, v);             # obtain access to full array

    # now extract only the dimension belonging to the current point
    Array               =   selectdim(ArrayFull,dim+1, slot_start+1);

    return Array
end

@for_petsc function DMStagGetGhostArrayLocationSlot(
    dm::DMStag{$PetscLib}, 
    ArrayFull::PermutedDimsArray, 
    loc::DMStagStencilLocation, 
    dof::Int
    )

    entriesPerElement   =   DMStagGetEntriesPerElement(dm)
    dim                 =   getdimension(dm);  
    slot                =   DMStagGetLocationSlot(dm, loc, dof); 
    slot_start          =   mod(slot,entriesPerElement);          # figure out which component we are interested in

    # now extract only the dimension belonging to the current point
    Array               =   selectdim(ArrayFull,dim+1, slot_start+1);

    return Array
end

"""
    slot = DMStagGetProductCoordinateLocationSlot(
        dm::DMStag,
        loc::DMStagStencilLocation
        )
    
Get slot for use with local product coordinate arrays.

    dm      - the DMStag object
    loc     - the grid location 
    slot    - the index to use in local arrays

# External Links
$(_doc_external("DMSTAG/DMStagGetProductCoordinateLocationSlot"))
"""
function DMStagGetProductCoordinateLocationSlot end

#REPAIR THAT AS WELL (OR GET RID OF...)
@for_petsc function DMStagGetProductCoordinateLocationSlot(
    dm::DMStag{$PetscLib},
    loc::DMStagStencilLocation
    )

    slot = Ref{$PetscInt}()
    @chk ccall((:DMStagGetProductCoordinateLocationSlot, $petsc_library), PetscErrorCode,
                ( 
                    CDMStag,   
                    DMStagStencilLocation, 
                    Ptr{$PetscInt}
                ), 
                dm, 
                loc, 
                slot
                )

    return slot[]
end

"""
    entriesPerElement = DMStagGetEntriesPerElement(dm::DMStag)
    
Get number of entries per element in the local representation. 

    dm                - the DMStag objects
    entriesPerElement - number of entries associated with each element in the local representation

# External Links
$(_doc_external("DMSTAG/DMStagGetEntriesPerElement"))
"""
function DMStagGetEntriesPerElement end

@for_petsc function DMStagGetEntriesPerElement(dm::DMStag)
    entriesPerElement = Ref{$PetscInt}()
    @chk ccall((:DMStagGetEntriesPerElement, $petsc_library), PetscErrorCode,
                ( 
                    CDMStag,  
                    Ptr{$PetscInt}
                ), 
                dm,  
                entriesPerElement
                )

    return entriesPerElement[]
end

"""
    stencilWidth = DMStagGetStencilWidth(dm::DMStag)
    
Get elementwise stencil width. 

    dm           - the DMStag objects
    stencilWidth - stencil/halo/ghost width in elements

# External Links
$(_doc_external("DMSTAG/DMStagGetStencilWidth"))
"""
function DMStagGetStencilWidth end

@for_petsc function DMStagGetStencilWidth(dm::DMStag)
    stencilWidth = Ref{$PetscInt}()
    @chk ccall((:DMStagGetStencilWidth, $petsc_library), PetscErrorCode,
                ( CDMStag,  Ptr{$PetscInt}), dm,  stencilWidth)

    return stencilWidth[]
end

"""

    slot = DMStagGetLocationSlot(
        dm::DMStag,
        loc::DMStagStencilLocation, 
        c
        )
            
Get index to use in accessing raw local arrays.

    dm      - the DMStag object
    loc     - location relative to an element
    c       - component ( the degree of freedom)
    slot    - index to use

# External Links
$(_doc_external("DMSTAG/DMStagGetLocationSlot"))
"""
function DMStagGetLocationSlot end

@for_petsc function DMStagGetLocationSlot(dm::DMStag{$PetscLib},loc::DMStagStencilLocation, c)
        
    slot = Ref{$PetscInt}()
    @chk ccall((:DMStagGetLocationSlot, $petsc_library), PetscErrorCode,
                ( 
                    CDMStag,   
                    DMStagStencilLocation, 
                    $PetscInt, 
                    Ptr{$PetscInt}
                ), 
                dm, 
                loc, 
                c, 
                slot
                )

    return slot[]
end

"""
    destroy(dm::DMStag)

Destroys a DMSTAG object and releases the memory

    dm 	- the DM object to destroy

# External Links
$(_doc_external("DM/DMDestroy"))
"""
function destroy(dm::DMStag) end

@for_petsc function destroy(dm::DMStag{$PetscLib})
    finalized($petsclib) ||
        @chk ccall((:DMDestroy, $petsc_library), PetscErrorCode, (Ptr{CDMStag},), dm)
    return nothing
end

"""
    Indices = DMStagGetCentralNodes(dm::DMStag)
    
Return indices of start and end of the central nodes of a local array built from the input `dm` (excluding ghost nodes).   

    dm 	        - the DMStag object
    Indices 	- indices of lower and upper range of center and vertex nodes
"""
function DMStagGetCentralNodes end

@for_petsc function DMStagGetCentralNodes(dm::DMStag)
    # In Julia, indices in arrays start @ 1, whereas they can go negative in C
    
    #g_start, g_N    =   getghostcorners(dm);  # MODIFY
    gc              =   getghostcorners(dm);  
  
    #start,N, nExtra =   getcorners(dm);       # MODIFY
    c               =   getcorners(dm); 

    # Note that we add the +1 for julia/petsc consistency
    center = (  x= (c.lower[1]+1):(c.upper[1]+1),
                y= (c.lower[2]+1):(c.upper[2]+1),  
                z= (c.lower[3]+1):(c.upper[3]+1) )

    vertex = (  x=c.lower[1]+1:c.upper[1]+2 ,
                y=c.lower[2]+1:c.upper[2]+2 ,  
                z=c.lower[3]+1:c.upper[3]+2 )

    return (center=center, vertex=vertex)
            
end


"""
    Bx = DMStagGetBoundaryTypes(dm::DMStag) in 1D
    Bx,By,Bz = DMStagGetBoundaryTypes(dm::DMStag) in 3D

Get boundary types.

        dm 	     - the DMStag object 
        Bx,By,Bz - boundary types

# External Links
$(_doc_external("DMSTAG/DMStagGetBoundaryTypes"))
"""
function DMStagGetBoundaryTypes end

@for_petsc function DMStagGetBoundaryTypes(dm::DMStag)

    Bx = Ref{$DMBoundaryType}()
    By = Ref{$DMBoundaryType}()
    Bz = Ref{$DMBoundaryType}()
      
    @chk ccall((:DMStagGetBoundaryTypes, $petsc_library), PetscErrorCode,
        (
            CDMStag,   
            Ptr{$DMBoundaryType}, Ptr{$DMBoundaryType}, Ptr{$DMBoundaryType}
        ), 
        dm, 
        Bx,By,Bz
        )

        dim = getdimension(dm); 

        if dim==1
            return Bx[]    
        elseif dim==2
            return Bx[], By[]
        elseif dim==3
            return Bx[], By[], Bz[]
        end
end


"""
    nRanks0 = DMStagGetNumRanks(dm::DMStag) in 1D
    nRanks0,nRanks1,nRanks2 = DMStagGetNumRanks(dm::DMStag) in 3D

Get number of ranks in each direction in the global grid decomposition.

        dm 	                     - the DMStag object 
        nRanks0,nRanks1,nRanks2  - number of ranks in each direction in the grid decomposition

# External Links
$(_doc_external("DMSTAG/DMStagGetNumRanks"))
"""
function DMStagGetNumRanks end

@for_petsc function  DMStagGetNumRanks(dm::DMStag)

    nRanks0 = Ref{$PetscInt}()
    nRanks1 = Ref{$PetscInt}()
    nRanks2 = Ref{$PetscInt}()
        
    @chk ccall((:DMStagGetNumRanks, $petsc_library), PetscErrorCode,
        (
            CDMStag, 
            Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}
        ), 
        dm, 
        nRanks0,nRanks1,nRanks2
        )

        dim = getdimension(dm); 
            
    if dim==1
        return nRanks0[]
    elseif dim==2
        return nRanks0[], nRanks1[]
    elseif dim==3
        return nRanks0[], nRanks1[], nRanks2[]
    end
end


"""
    DMStagVecSetValuesStencil(
        dm::DMStag,
        vec::Abstractvec, 
        pos::DMStagStencil, 
        val::Float64, 
        insertMode::InsertMode
        )
    
This puts a value inside a global vector using DMStagStencil.
    
    dm  - the DMStag object
    vec - the Vec
    n   - the number of values (do not fill if only 1)
    pos - the location of the set values, given by a DMStagStencil struct
    val - the value to be set
    insertMode  - INSERT_VALUES or ADD_VALUES

# External Links
$(_doc_external("DMSTAG/DMStagVecSetValuesStencil"))
"""
function DMStagVecSetValuesStencil end

@for_petsc function  DMStagVecSetValuesStencil(
        dm::DMStag{$PetscLib}, 
        vec::AbstractVec{$PetscScalar}, 
        pos::DMStagStencil{$PetscInt}, 
        val, 
        insertMode::InsertMode)

    n=1;
    @chk ccall((:DMStagVecSetValuesStencil, $petsc_library), PetscErrorCode,
                (
                    CDMStag, 
                    CVec, 
                    $PetscInt, 
                    Ptr{DMStagStencil{$PetscInt}}, 
                    Ptr{$PetscScalar}, 
                    InsertMode
                ), 
                dm, 
                vec.ptr, 
                n, 
                Ref{DMStagStencil{$PetscInt}}(pos), 
                Ref{$PetscScalar}(val), 
                insertMode
                )

    return nothing
end

"""

    DMStagVecSetValuesStencil(
        dm::DMStag, 
        vec::AbstractVec{PetscScalar}, 
        n, 
        pos::Vector{DMStagStencil}, 
        values::Vector{PetscScalar}, 
        insertMode::InsertMode
        )    

This puts values inside a global vector using DMStagStencil
    
    dm  - the DMStag object
    vec - the Vec
    n   - the number of values (do not fill if only 1)
    pos - the location of the set values, given by a DMStagStencil struct
    val - the value to be set
    insertMode  - INSERT_VALUES or ADD_VALUES

# External Links
$(_doc_external("DMSTAG/DMStagVecGetValuesStencil"))
"""
function DMStagVecSetValuesStencil end

@for_petsc function  DMStagVecSetValuesStencil(
    dm::DMStag{$PetscLib}, 
    vec::AbstractVec{$PetscScalar}, 
    n::Integer,
    pos::Vector{$DMStagStencil{$PetscInt}}, 
    values::Vector{$PetscScalar}, 
    insertMode::InsertMode
    )

    i = 1;
    while i <= n
        pos0 = pos[i];
        val  = values[i];
        m=1;
        @chk ccall((:DMStagVecSetValuesStencil, $petsc_library), PetscErrorCode,
                    (
                        CDMStag, 
                        CVec, 
                        $PetscInt, 
                        Ptr{DMStagStencil{$PetscInt}}, 
                        Ptr{$PetscScalar}, 
                        InsertMode
                    ), 
                    dm, 
                    vec.ptr, 
                    m, 
                    Ref{DMStagStencil{$PetscInt}}(pos0), 
                    Ref{$PetscScalar}(val), 
                    insertMode
                    )
    i += 1;
    end
    return nothing
end


"""
    val = DMStagVecGetValuesStencil(
        dm::DMStag, 
        vec::AbstractVec, 
        pos::DMStagStencil
        )

Get vector values using grid indexing

    dm  - the DMStag object
    vec - the vector object
    n   - the number of values to obtain (do not fill if only one)
    pos - locations to obtain values from (as an array of DMStagStencil values) 
    val - value at the point 

# External Links
$(_doc_external("DMSTAG/DMStagVecGetValuesStencil"))
"""
function DMStagVecGetValuesStencil end

@for_petsc function  DMStagVecGetValuesStencil(
    dm::DMStag{$PetscLib}, 
    vec::AbstractVec{$PetscScalar}, 
    pos::DMStagStencil
    )

    n=1;
    val = Ref{$PetscScalar}()
    @chk ccall((:DMStagVecGetValuesStencil, $petsc_library), PetscErrorCode,
                (
                    CDMStag, 
                    CVec, 
                    $PetscInt, 
                    Ptr{DMStagStencil{$PetscInt}}, 
                    Ptr{$PetscScalar}
                ), 
                dm, 
                vec.ptr, 
                n, 
                Ref{DMStagStencil{$PetscInt}}(pos), 
                val
                )
    
    return val[]
end


"""
    val = DMStagVecGetValuesStencil(
        dm::DMStag, 
        vec::AbstractVec, 
        n, 
        pos::Vector{DMStagStencil}
        )

Get vector values using grid indexing.

    dm  - the DMStag object
    vec - the vector object
    n   - the number of values to obtain (do not fill if only one)
    pos - locations to obtain values from (as an array of DMStagStencil values) 
    val - value at the point 

# External Links
$(_doc_external("DMSTAG/DMStagVecGetValuesStencil"))
"""
function DMStagVecGetValuesStencil end

@for_petsc function  DMStagVecGetValuesStencil(
    dm::DMStag{$PetscLib}, 
    vec::AbstractVec{$PetscScalar}, 
    n,
    pos::Vector{$DMStagStencil{$PetscInt}}
    )

    i = 1;
    values = zeros(n);
    while i <= n
        pos0 = pos[i];
        m=1;
        val = Ref{$PetscScalar}()
        @chk ccall((:DMStagVecGetValuesStencil, $petsc_library), PetscErrorCode,
                    (
                        CDMStag, 
                        CVec, 
                        $PetscInt, 
                        Ptr{DMStagStencil{$PetscInt}}, 
                        Ptr{$PetscScalar}
                    ), 
                    dm, 
                    vec.ptr, 
                    m, 
                    Ref{DMStagStencil{$PetscInt}}(pos0), 
                    val
                    )
            
        values[i] = val[];
        i += 1;
    end
    return values
end

"""
    val =  DMStagMatGetValuesStencil(
        dm::DMStag,
        mat::AbstractMat, 
        posRow::DMStagStencil,  
        posCol::DMStagStencil
        )

This reads a single value from a matrix DMStagStencil
      
    dm      - the DMStag object
    mat     - the Mat
    posRow  - the location of the row of the set value, given by a DMStagStencil struct (as a vector)
    posCol  - the location of the row of the set value, given by a DMStagStencil struct (as a vector)
    val     - the value

# External Links
$(_doc_external("DMSTAG/DMStagMatGetValuesStencil"))
"""
function DMStagMatGetValuesStencil end

@for_petsc function  DMStagMatGetValuesStencil(
    dm::DMStag{$PetscLib}, 
    mat::AbstractMat{$PetscScalar},  
    posRow::DMStagStencil, 
    posCol::DMStagStencil
    )

    nRow= 1;
    nCol= 1;
    val = Ref{$PetscScalar}()
    @chk  ccall((:DMStagMatGetValuesStencil, $petsc_library), PetscErrorCode,
                (
                    CDMStag, 
                    CMat, 
                    $PetscInt, 
                    Ptr{DMStagStencil{$PetscInt}}, 
                    $PetscInt, 
                    Ptr{DMStagStencil{$PetscInt}}, 
                    Ptr{$PetscScalar}
                ), 
                dm, 
                mat.ptr, 
                nRow, 
                Ref{DMStagStencil{$PetscInt}}(posRow), 
                nCol, 
                Ref{DMStagStencil{$PetscInt}}(posCol), 
                val
                )

    return val[]
end

"""
    val =  DMStagMatGetValuesStencil(
        dm::DMStag, 
        mat::AbstractMat{PetscScalar}, 
        nRow,  
        posRow::Vector{DMStagStencil}, 
        nCol, 
        posCol::Vector{DMStagStencil}
        )

This reads a single value from a matrix DMStagStencil.

    dm      - the DMStag object
    mat     - the Mat
    posRow  - the location of the row of the set value, given by a DMStagStencil struct (as a vector)
    posCol  - the location of the row of the set value, given by a DMStagStencil struct (as a vector)
    val     - the value


# External Links
$(_doc_external("DMSTAG/DMStagMatGetValuesStencil"))
"""
function DMStagMatGetValuesStencil end

@for_petsc function  DMStagMatGetValuesStencil(
    dm::DMStag{$PetscLib},
    mat::AbstractMat{$PetscScalar}, 
    nRow,  
    posRow::Vector{$DMStagStencil{$PetscInt}}, 
    nCol, 
    posCol::Vector{$DMStagStencil{$PetscInt}}
    )

    i = 1;
    j = 1;
    values = zeros(nRow*nCol);
    while i <= nRow
        while j <= nCol
            posr = posRow[i];
            posc = posCol[j];
            n_Row= 1;
            n_Col= 1;
            val = Ref{$PetscScalar}()
            @chk  ccall((:DMStagMatGetValuesStencil, $petsc_library), PetscErrorCode,
                        (
                            CDMStag, 
                            CMat, 
                            $PetscInt, 
                            Ptr{DMStagStencil{$PetscInt}}, 
                            $PetscInt, 
                            Ptr{DMStagStencil{$PetscInt}}, 
                            Ptr{$PetscScalar}
                        ), 
                        dm, 
                        mat.ptr, 
                        n_Row, 
                        Ref{DMStagStencil{$PetscInt}}(posr), 
                        n_Col, 
                        Ref{DMStagStencil{$PetscInt}}(posc), 
                        val
                        )
            values[i*j] = val[]
            j += 1;
        end
        i += 1;
    end

    return values
end

""" 
    indices = LocalInGlobalIndices(dm::DMStag)

Give the non-ghosted indices in the local vector that contribute to the global vector.

    dm      - the DMStag object
    indices - local indices
"""
function LocalInGlobalIndices end

function LocalInGlobalIndices(dm::DMStag)
    # note: this can likely be done more efficiently and will have to be modified in parallel
    ind_g   =   createglobalvector(dm)
    v_ind_l =   createlocalvector(dm)

    # Set indices in local vector
    for i=1:length(v_ind_l)
        v_ind_l[i] = i
    end
    
    update!(ind_g, v_ind_l, INSERT_VALUES); # update global vector

    return Int64.(ind_g)
end
    
"""

    DMStagMatSetValuesStencil(
        dm::DMStag, 
        mat::AbstractMat, 
        nRow,  
        posRow::Vector{DMStagStencil}, 
        nCol, 
        posCol::Vector{DMStagStencil}, 
        values::Vector{PetscScalar}, 
        insertMode::InsertMode
        )

This puts values inside a matrix using DMStagStencil position.
    
    dm	        - the DMStag object
    mat	        - the Mat
    posRow	    - the location of the row of the set value, given by a DMStagStencil struct (as a vector)
    posCol	    - the location of the row of the set value, given by a DMStagStencil struct (as a vector)
    val	        - the value to be set
    insertMode	- INSERT_VALUES or ADD_VALUES

# External Links
$(_doc_external("DMSTAG/DMStagMatSetValuesStencil"))
"""
function DMStagMatSetValuesStencil end

@for_petsc function  DMStagMatSetValuesStencil(
    dm::DMStag{$PetscLib}, 
    mat::AbstractMat{$PetscScalar},  
    posRow::DMStagStencil, 
    posCol::DMStagStencil, 
    val, 
    insertMode::InsertMode
    )

    nRow= 1;
    nCol= 1;
    @chk ccall((:DMStagMatSetValuesStencil, $petsc_library), PetscErrorCode,
                (
                    CDMStag, 
                    CMat, 
                    $PetscInt, 
                    Ptr{DMStagStencil{$PetscInt}},  
                    $PetscInt, 
                    Ptr{DMStagStencil{$PetscInt}}, 
                    Ptr{$PetscScalar}, 
                    InsertMode
                ), 
                dm, 
                mat.ptr, 
                nRow, 
                Ref{DMStagStencil{$PetscInt}}(posRow), 
                nCol, 
                Ref{$DMStagStencil{$PetscInt}}(posCol), 
                Ref{$PetscScalar}(val), 
                insertMode
                )

    return nothing
end

"""
    DMStagMatSetValuesStencil(
        dm::DMStag, 
        mat::AbstractMat, 
        nRow,  
        posRow::Vector{DMStagStencil}, 
        nCol, 
        posCol::Vector{DMStagStencil}, 
        values::Vector{PetscScalar}, 
        insertMode::InsertMode
        )

This puts values inside a matrix using DMStagStencil position

    dm	        - the DMStag object
    mat	        - the Mat
    posRow	    - the location of the row of the set value, given by a DMStagStencil struct (as a vector)
    posCol	    - the location of the row of the set value, given by a DMStagStencil struct (as a vector)
    val	        - the value to be set
    insertMode	- INSERT_VALUES or ADD_VALUES

# External Links
$(_doc_external("DMSTAG/DMStagMatSetValuesStencil"))
"""
function DMStagMatSetValuesStencil end

@for_petsc function  DMStagMatSetValuesStencil(
    dm::DMStag{$PetscLib}, 
    mat::AbstractMat{$PetscScalar}, 
    nRow,  
    posRow::Vector{$DMStagStencil{$PetscInt}}, 
    nCol, 
    posCol::Vector{$DMStagStencil{$PetscInt}}, 
    values::Vector{$PetscScalar}, 
    insertMode::InsertMode
    )


    i = 1;
    j = 1;
    while i <= nRow
        while j <= nCol
            posr = posRow[i];
            posc = posCol[j];
            val  = values[i*j];
            n_Row= 1;
            n_Col= 1;
            @chk ccall((:DMStagMatSetValuesStencil, $petsc_library), PetscErrorCode,
                        (
                            CDMStag, 
                            CMat, 
                            $PetscInt, 
                            Ptr{DMStagStencil{$PetscInt}},  
                            $PetscInt, 
                            Ptr{DMStagStencil{$PetscInt}}, 
                            Ptr{$PetscScalar}, 
                            InsertMode
                        ), 
                        dm, 
                        mat.ptr, 
                        n_Row, 
                        Ref{DMStagStencil{$PetscInt}}(posr), 
                        n_Col, 
                        Ref{$DMStagStencil{$PetscInt}}(posc), 
                        Ref{$PetscScalar}(val), 
                        insertMode
                        )
        j += 1;
        end
    i += 1;
    end
    return nothing
end

# NOTE: We should likely rewrite this to make it consistent with update!, which has a slighlt different calling sequence
"""
    DMLocalToGlobal(
        dm::DMStag,
        l::AbstractVec, 
        mode::InsertMode,
        g::AbstractVec
        )

Updates global vectors from local vectors. 

    dm 	 - the DM object
    l 	 - the local vector
    mode - if INSERT_VALUES then no parallel communication is used, if ADD_VALUES then all ghost points from the same base point accumulate into that base point.
    g 	 - the global vector    
    
# External Links
$(_doc_external("DM/DMLocalToGlobal"))

"""
function DMLocalToGlobal end

@for_petsc function DMLocalToGlobal(
    l::AbstractVec{$PetscScalar}, 
    mode::InsertMode,
    g::AbstractVec{$PetscScalar}
    )

    update!(l,g,mode);

    return nothing
end

@for_petsc function DMLocalToGlobal(
    l::AbstractVec{$PetscScalar}, 
    mode::InsertMode,
    g::CVec
    )

    update!(l,g,mode);

    return nothing
end

@for_petsc function DMGlobalToLocal(
    g::AbstractVec{$PetscScalar}, 
    mode::InsertMode,
    l::AbstractVec{$PetscScalar}
    )

    update!(g,l,mode)

    return nothing
end


@for_petsc function DMGlobalToLocal(
    g::CVec, 
    mode::InsertMode,
    l::AbstractVec{$PetscScalar}
    )
    update!(g,l,mode)

    return nothing
end


"""
    stencilType = DMStagGetStencilType(dm::DMStag)

Get elementwise ghost/halo stencil type.

    dm          - the DMStag object
    stencilType - the elementwise ghost stencil type: DMSTAG_STENCIL_BOX, DMSTAG_STENCIL_STAR, or DMSTAG_STENCIL_NONE

# External Links
$(_doc_external("DMSTAG/DMStagGetStencilType"))
"""
function DMStagGetStencilType end

@for_petsc function DMStagGetStencilType(dm::DMStag)
    stencilType =  Ref{DMStagStencilType}()

    @chk ccall((:DMStagGetStencilType, $petsc_library), PetscErrorCode, 
    (
        CDMStag, 
        Ptr{DMStagStencilType}
    ), 
    dm, 
    stencilType
    )

    return stencilType[]
end

"""
    fr_X,fr_Y,fr_Z = DMStagGetIsFirstRank(dm::DMStag)
    
Returns boolean value to indicate whether this rank is first in each direction in the rank grid. Arguments corresponding to higher dimensions are ignored for 1D and 2D grids. These arguments may be set to NULL in this case.

    dm             - the DMStag object
    fr_X,fr_Y,fr_Z - whether this rank is first in each direction

# External Links
$(_doc_external("DMSTAG/DMStagGetIsFirstRank"))
"""
function DMStagGetIsFirstRank end

@for_petsc function DMStagGetIsFirstRank(dm::DMStag)
    fr_X = Ref{PetscBool}()
    fr_Y = Ref{PetscBool}()
    fr_Z = Ref{PetscBool}()
        
    @chk ccall((:DMStagGetIsFirstRank, $petsc_library), PetscErrorCode, 
    (
        CDMStag, 
        Ptr{PetscBool}, Ptr{PetscBool}, Ptr{PetscBool}
    ), 
    dm, 
    fr_X, fr_Y, fr_Z
    )
        
    return fr_X[]== PETSC_TRUE, fr_Y[]== PETSC_TRUE, fr_Z[]== PETSC_TRUE
end



"""
    fr_X,fr_Y,fr_Z = DMStagGetIsLastRank(dm::DMStag)

Returns boolean value to indicate whether this rank is last in each direction in the rank grid.

    dm             - the DMStag object
    fr_X,fr_Y,fr_Z - whether this rank is last in each direction

# External Links
$(_doc_external("DMSTAG/DMStagGetIsLastRank"))
"""
function DMStagGetIsLastRank end

@for_petsc function DMStagGetIsLastRank(dm::DMStag)
    fr_X = Ref{PetscBool}()
    fr_Y = Ref{PetscBool}()
    fr_Z = Ref{PetscBool}()
        
    @chk ccall((:DMStagGetIsLastRank, $petsc_library), PetscErrorCode, 
    (
        CDMStag, 
        Ptr{PetscBool}, Ptr{PetscBool}, Ptr{PetscBool}
    ), 
    dm, 
    fr_X, fr_Y, fr_Z
    )
        
    return fr_X[]== PETSC_TRUE, fr_Y[]== PETSC_TRUE, fr_Z[]== PETSC_TRUE
end

"""
    setuniformcoordinates!(
        dm::DMStag,
        xyzmin::NTuple{N, Real},
        xyzmax::NTuple{N, Real},
    ) where {N}

Set uniform coordinates for the `dmstag` using the lower and upper corners defined
by the `NTuples` `xyzmin` and `xyzmax`. If `N` is less than the dimension of the
`dm` then the value of the trailing coordinates is set to `0`.

# External Links
$(_doc_external("DMSTAG/DMStagSetUniformCoordinatesExplicit"))
"""
function setuniformcoordinates!(dm::DMStag,xyzmin,xyzmax) end

@for_petsc function setuniformcoordinates!(
    dm::DMStag{$PetscLib},
    xyzmin::NTuple{N, Real},
    xyzmax::NTuple{N, Real},
) where {N}
    xmin = $PetscReal(xyzmin[1])
    xmax = $PetscReal(xyzmax[1])

    ymin = (N > 1) ? $PetscReal(xyzmin[2]) : $PetscReal(0)
    ymax = (N > 1) ? $PetscReal(xyzmax[2]) : $PetscReal(0)

    zmin = (N > 2) ? $PetscReal(xyzmin[3]) : $PetscReal(0)
    zmax = (N > 2) ? $PetscReal(xyzmax[3]) : $PetscReal(0)

    @chk ccall(
        (:DMStagSetUniformCoordinatesExplicit, $petsc_library),
        PetscErrorCode,
        (
            CDMStag,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscReal,
        ),
        dm,
        xmin,
        xmax,
        ymin,
        ymax,
        zmin,
        zmax,
    )
    return nothing
end

"""
    getcorners(dm::DMSTAG)

Returns a `NamedTuple` with the global indices (excluding ghost points) of the
`lower` and `upper` corners as well as the `size`.

# External Links
$(_doc_external("DMSTAG/DMStagGetCorners"))
"""
function getcorners(dm::DMStag) end

@for_petsc function getcorners(dm::DMStag{$PetscLib})
    corners     = [$PetscInt(0), $PetscInt(0), $PetscInt(0)]
    local_size  = [$PetscInt(0), $PetscInt(0), $PetscInt(0)]
    nExtra      = [$PetscInt(0), $PetscInt(0), $PetscInt(0)]   
    @chk ccall(
        (:DMStagGetCorners, $petsc_library),
        PetscErrorCode,
        (
            CDMStag,
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
        ),
        dm,
        Ref(corners, 1),
        Ref(corners, 2),
        Ref(corners, 3),
        Ref(local_size, 1),
        Ref(local_size, 2),
        Ref(local_size, 3),
        Ref(nExtra,     1),
        Ref(nExtra,     2),
        Ref(nExtra,     3),
    )
    corners .+= 1
    return (
        lower = corners,
        upper = corners .+ local_size .- $PetscInt(1),
        size  = local_size,
        extra = nExtra
    )
end


"""
    getghostcorners(dm::DMStag)

Returns a `NamedTuple` with the global indices (including ghost points) of the
`lower` and `upper` corners as well as the `size`.

# External Links
$(_doc_external("DMSTAG/DMStagGetGhostCorners"))
"""
function getghostcorners(dm::DMStag) end

@for_petsc function getghostcorners(dm::DMStag)
    corners     = [$PetscInt(0), $PetscInt(0), $PetscInt(0)]
    local_size  = [$PetscInt(0), $PetscInt(0), $PetscInt(0)]
    nExtra      = [$PetscInt(0), $PetscInt(0), $PetscInt(0)]  
    @chk ccall(
        (:DMStagGetGhostCorners, $petsc_library),
        PetscErrorCode,
        (
            CDMStag,
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
            Ref{$PetscInt},
        ),
        dm,
        Ref(corners, 1),
        Ref(corners, 2),
        Ref(corners, 3),
        Ref(local_size, 1),
        Ref(local_size, 2),
        Ref(local_size, 3),
        Ref(nExtra,     1),
        Ref(nExtra,     2),
        Ref(nExtra,     3),
    )
    corners .+= 1
    return (
        lower = corners,
        upper = corners .+ local_size .- 1,
        size  = local_size,
        extra = nExtra
    )
end


Base.show(io::IO, dm::DMStag) = _show(io, dm)
