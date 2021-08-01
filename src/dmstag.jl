# Attempt at include dmstag functions
const CDMStag = Ptr{Cvoid}
const CDMStagType = Cstring

mutable struct DMStag{PetscLib} <: AbstractDM{PetscLib}
    ptr::CDMStag
    comm::MPI.Comm
    #dim::PetscInt
    opts::Options{PetscLib}
end

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


#Base.eltype(::DMStag{T,PetscInt}) where {T,PetscInt} = T,PetscInt


"""
    dm = DMStagCreate1d(::PetscLib,
        comm::MPI.Comm, 
        bndx::DMBoundaryType, 
        M, dofVertex, 
        dofCenter, 
        stencilType::DMStagStencilType=DMSTAG_STENCIL_BOX, 
        stencilWidth=2, 
        lx=C_NULL; 
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

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMSTAG/DMStagCreate1d.html)
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
    options...,
)

    if isempty(lx); lx = C_NULL; end
    #opts = Options{$PetscScalar}(kwargs...)
    opts = Options($petsclib, options...)

    #dm  = DMStag{$PetscScalar,$PetscInt}(C_NULL, comm, 1, opts)   # retrieve options
    dm  = DMStag{$PetscLib}(C_NULL, comm, opts)   # retrieve options
        
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
        options...
        )

Creates a 2D DMStag object.
        ::PetscLib      -   PETSc library
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

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMSTAG/DMStagCreate2d.html)
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
    options...,
)
        
    if isempty(lx); lx = C_NULL; end
    if isempty(ly); ly = C_NULL; end
    opts = Options($petsclib, options...)

    dm = DMStag{$PetscLib}(C_NULL, comm, opts)

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
        options...
        )

Creates a 3D DMStag object.

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
        options...      -   [Optional] keyword arguments (see PETSc webpage), specifiable as stag_grid_x=100, etc. 

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMSTAG/DMStagCreate3d.html)
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
    options...,
    )
        
    if isempty(lx); lx = C_NULL; end
    if isempty(ly); ly = C_NULL; end
    if isempty(lz); lz = C_NULL; end
    opts = Options($petsclib, options...)
        
    dm = DMStag{$PetscLib}(C_NULL, comm, opts)

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
    DMSetUp(dm::DMStag)

Sets up the data structures inside a DM object (automatically called in the DMStagCreate routines). 

        dm    -   the DMStag object 

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DM/DMSetUp.html)
"""
function DMSetUp end

@for_petsc function DMSetUp(dm::DMStag{$PetscLib})

    @chk ccall((:DMSetUp, $petsc_library), PetscErrorCode, (CDMStag, ), dm )

    return nothing
end

"""
    setfromoptions!(dm::DMStag)

Sets parameters in a DM from the options database (automatically called in the DMStagCreate routines).

    dm              -   the DMStag object 

More info on [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DM/DMSetFromOptions.html)
"""
function setfromoptions! end

@for_petsc function setfromoptions!(dm::DMStag{$PetscLib})

    @chk ccall((:DMSetFromOptions, $petsc_library), PetscErrorCode, (CDMStag, ), dm )

    return nothing
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

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMSTAG/DMStagCreateCompatibleDMStag.html)
"""
function DMStagCreateCompatibleDMStag end

@for_petsc function DMStagCreateCompatibleDMStag(
    dm::DMStag{$PetscLib}, 
    dofVertex=0, 
    dofEdge=0, 
    dofFace=0, 
    dofElement=0; 
    kwargs...
    )

    comm  = dm.comm

    opts = Options($petsclib, kwargs...)

    dmnew = DMStag{$PetscLib}(C_NULL, comm, opts)

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
    dof0,dof1,dof2,dof3 = DMStagGetDOF(dm::DMStag)
        
Get number of DOF associated with each stratum of the grid. 
        
    dm      - the DMStag object 
    dof0 	- the number of points per 0-cell (vertex/node)
    dof1 	- the number of points per 1-cell (element in 1D, edge in 2D and 3D)
    dof2 	- the number of points per 2-cell (element in 2D, face in 3D)
    dof3 	- the number of points per 3-cell (element in 3D)     

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMSTAG/DMStagGetDOF.html)
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

    dim   = DMGetDimension(dm)

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

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMSTAG/DMStagGetGlobalSizes.html)
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

        dim   = DMGetDimension(dm)
        
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

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMSTAG/DMStagGetLocalSizes.html)
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

        dim   = DMGetDimension(dm)
        
    if dim==1    
        return M[]
    elseif dim==2
        return M[], N[] 
    elseif dim==3
        return M[], N[], P[]    
    end
end


"""
    DMStagSetUniformCoordinatesProduct(
        dm::DMStag, 
        xmin, xmax,
        ymin=0, ymax=0, 
        zmin=0, zmax=0)

Set the coordinate DM to be a DMProduct of 1D DMStag objects, each of which have a coordinate DM (also a 1d DMStag) holding uniform coordinates. 

    dm                             - the DMStag object 
    xmin,xmax,ymin,ymax,zmin,zmax  - maximum and minimum global coordinate values

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMSTAG/DMStagSetUniformCoordinatesProduct.html)
"""
function DMStagSetUniformCoordinatesProduct end

@for_petsc function DMStagSetUniformCoordinatesProduct(
    dm::DMStag{$PetscLib}, 
    xmin, xmax, 
    ymin=0, ymax=0, 
    zmin=0, zmax=0
    )
        
    @chk ccall((:DMStagSetUniformCoordinatesProduct, $petsc_library), PetscErrorCode,
                ( 
                    CDMStag,   
                    $PetscScalar, $PetscScalar, 
                    $PetscScalar, $PetscScalar, 
                    $PetscScalar, $PetscScalar
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

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMSTAG/DMStagSetUniformCoordinatesExplicit.html)
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
    vec = DMCreateGlobalVector(
        dm::DMStag; 
        write_val=true, read_val=true
        )
    
Creates a global vector from a DM object. 

NOTE: for now this is initialized sequentially; MPI should be added

        dm 	- the DM object
        vec - the global vector

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DM/DMCreateGlobalVector.html) 
"""
function DMCreateGlobalVector end

@for_petsc function DMCreateGlobalVector(
    dm::DMStag{$PetscLib}; 
    write_val=true, read_val=true
    )

    v = VecSeq(C_NULL, [0.0])  # empty vector
        
    @chk ccall((:DMCreateGlobalVector, $petsc_library), PetscErrorCode, (CDMStag, Ptr{CVec}), dm, v)

    # Link a julia array to the values from the new vector
    # If we modify values here, it will automatically be changed in the PetcVec as well
    v.array = unsafe_localarray($PetscScalar, v.ptr; write=write_val, read=read_val)
        
    return v
end

"""
    vec = DMCreateLocalVector(dm::DMStag; 
    write_val=true, read_val=true
    )
    
Creates a local vector from a DM object.

NOTE: for now this is initialized sequentially; MPI should be added

        dm 	- the DM object
        vec - the local vector

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DM/DMCreateLocalVector.html)
"""
function DMCreateLocalVector end

@for_petsc function DMCreateLocalVector(
    dm::DMStag{$PetscLib}; 
    write_val=true, read_val=true
    )

    v = VecSeq(C_NULL, [0.0])  # empty vector
        
    @chk ccall((:DMCreateLocalVector, $petsc_library), PetscErrorCode, (CDMStag, Ptr{CVec}), dm, v)

    # Link a julia array to the values from the new vector
    # If we modify values here, it will automatically be changed in the PetcVec as well
    v.array = unsafe_localarray($PetscScalar, v.ptr; write=write_val, read=read_val)
        
    return v
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
    v.array     =   unsafe_localarray($PetscScalar, v.ptr;  write=true, read=true)

    X1          =   DMStagVecGetArray(dm, v.array) 
        
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
    v.array     =   unsafe_localarray($PetscScalar, v.ptr;  write=false, read=true)

    X1          =   DMStagVecGetArray(dm, v.array) 
        
    return X1
end

"""
    X1 = DMStagVecGetArray(dm::DMStag, v::Vector)

Returns a julia array from a vector `v`, in the same shape as the DMSTAG, which can be used to set values.
"""
function DMStagVecGetArray end

@for_petsc function DMStagVecGetArray(dm::DMStag{$PetscLib}, v::Vector)

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
    dim                 =   DMGetDimension(dm);  
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
    dim                 =   DMGetDimension(dm);  
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

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMSTAG/DMStagGetProductCoordinateLocationSlot.html)
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

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-master/docs/manualpages/DMSTAG/DMStagGetEntriesPerElement.html)
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

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMSTAG/DMStagGetStencilWidth.html)
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

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMSTAG/DMStagGetLocationSlot.html)
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

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DM/DMDestroy.html)
"""
function destroy end

@for_petsc function destroy(dm::DMStag{$PetscLib})
    finalized($petsclib) ||
        @chk ccall((:DMDestroy, $petsc_library), PetscErrorCode, (Ptr{CDMStag},), dm)
    return nothing
end

"""
    type = gettype(dm::DMStag)

Gets the DM type name (as a string) from the DM.

    dm  - The DM
    type- The DM type name
        
From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DM/DMGetType.html)
"""
function gettype end

@for_petsc function gettype(dm::DMStag{$PetscLib})
    t_r = Ref{CDMStagType}()
    @chk ccall((:DMGetType, $petsc_library), PetscErrorCode, (CDMStag, Ptr{CDMStagType}), dm, t_r)
    return unsafe_string(t_r[])
end

"""
    view(dm::DMStag, viewer::Viewer=ViewerStdout(dm.comm))
    
Views a DMSTAG object.

    dm      - the DM object to view 
    viewer  - the viewer 

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DM/DMView.html)
"""
function view end

@for_petsc function view(dm::DMStag{$PetscLib}, viewer::AbstractViewer{$PetscLib}=ViewerStdout($petsclib, dm.comm))
    @chk ccall((:DMView, $petsc_library), PetscErrorCode, 
                (CDMStag, CPetscViewer),
            dm, viewer);
    return nothing
end

""" 
        
    x,m,nExtrax = DMStagGetCorners(dm:DMStag)   in 1D
    x[],m[],nExtrax[] = DMStagGetCorners(dm:DMStag)   in 2D or 3D

Returns the global element indices of the local region (excluding ghost points). 

    dm 	    - the DMStag object
    x,y,z 	- starting element indices in each direction
    m,n,p 	- element widths in each direction
    nExtrax,nExtray,nExtraz 	- number of extra partial elements in each direction. 

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMSTAG/DMStagGetCorners.html)
"""
function DMStagGetCorners end

@for_petsc function  DMStagGetCorners(dm::DMStag)

    x = Ref{$PetscInt}()
    y = Ref{$PetscInt}()
    z = Ref{$PetscInt}()
    m = Ref{$PetscInt}()
    n = Ref{$PetscInt}()
    p = Ref{$PetscInt}()
    nExtrax = Ref{$PetscInt}()
    nExtray = Ref{$PetscInt}()
    nExtraz = Ref{$PetscInt}()
        
    @chk ccall((:DMStagGetCorners, $petsc_library), PetscErrorCode,
        (
            CDMStag,   
            Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt},
            Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt},
            Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}
        ), 
        dm, 
        x,y,z, 
        m,n,p, 
        nExtrax,nExtray,nExtraz 
        )

        dim = DMGetDimension(dm);

        if dim==1
            X = (x[],)
            M = (m[],)
            NEXTRA = (nExtrax[],)
            return X[1], M[1], NEXTRA[1]    
        elseif dim==2
            return (x[], y[]), (m[],n[]), (nExtrax[],nExtray[])    
        elseif dim==3
            return (x[], y[], z[]), (m[],n[],p[]), (nExtrax[],nExtray[],nExtraz[])    
        end
end
    
""" 
    x,m = DMStagGetGhostCorners(dm:DMStag)   in 1D
    x[],m[] = DMStagGetGhostCorners(dm:DMStag)   in 2D or 3D
    
Return global element indices of the local region (including ghost points). 
       
    dm 	    - the DMStag object
    x[] 	- starting element indices in each direction
    m[] 	- element widths in each direction

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMSTAG/DMStagGetGhostCorners.html)
"""
function DMStagGetGhostCorners end

@for_petsc function  DMStagGetGhostCorners(dm::DMStag)

    x = Ref{$PetscInt}()
    y = Ref{$PetscInt}()
    z = Ref{$PetscInt}()
    m = Ref{$PetscInt}()
    n = Ref{$PetscInt}()
    p = Ref{$PetscInt}()
        
    @chk ccall((:DMStagGetGhostCorners, $petsc_library), PetscErrorCode,
        (
            CDMStag,   
            Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt},
            Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}
        ), 
        dm, 
        x,y,z, 
        m,n,p
        )

        dim =   DMGetDimension(dm);

        if dim==1
            X = (x[],)
            M = (m[],)
            return X[1], M[1]
        elseif dim==2
            return (x[], y[]), (m[],n[])    
        elseif dim==3
            return (x[], y[], z[]), (m[],n[],p[])   
        end
end

"""
    Cen_start, Cen_end = DMStagGetCentralNodes(dm::DMStag)
    
Return indices of start and end of the central nodes of a local array built from the input `dm` (excluding ghost nodes).   

    dm 	                - the DMStag object
    Cen_start, Cen_end 	- indices of start and finish of central nodes
"""
function DMStagGetCentralNodes end

@for_petsc function DMStagGetCentralNodes(dm::DMStag)
    # in Julia, indices in arrays start @ 1, whereas they can go negative in C
    # This routine  

    g_start, g_N    =   DMStagGetGhostCorners(dm);
    g_width         =   DMStagGetStencilWidth(dm);
    start,N, nExtra =   DMStagGetCorners(dm);
    dim             =   DMGetDimension(dm); 
        
    Cen_start       =   zeros(Int64,dim)
    for i=1:length(g_start)
        Cen_start[i] = -g_start[i] + 1;
    end

    Cen_end         =   Cen_start  .+ N .- 1;
    return Cen_start, Cen_end
end


"""
    Bx = DMStagGetBoundaryTypes(dm::DMStag) in 1D
    Bx,By,Bz = DMStagGetBoundaryTypes(dm::DMStag) in 3D

Get boundary types.

        dm 	     - the DMStag object 
        Bx,By,Bz - boundary types

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMSTAG/DMStagGetBoundaryTypes.html)
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

        dim = DMGetDimension(dm); 

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

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMSTAG/DMStagGetNumRanks.html)
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

        dim = DMGetDimension(dm); 
            
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

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMSTAG/DMStagVecSetValuesStencil.html)
"""
function DMStagVecSetValuesStencil end

@for_petsc function  DMStagVecSetValuesStencil(dm::DMStag{$PetscLib}, vec::AbstractVec{$PetscScalar}, pos::DMStagStencil{$PetscInt}, val, insertMode::InsertMode)

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

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMSTAG/DMStagVecSetValuesStencil.html)
"""
function DMStagVecSetValuesStencil end

@for_petsc function  DMStagVecSetValuesStencil(
    dm::DMStag{$PetscLib}, 
    vec::AbstractVec{$PetscScalar}, 
    n,
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

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMSTAG/DMStagVecGetValuesStencil.html)
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

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMSTAG/DMStagVecGetValuesStencil.html)
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

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMSTAG/DMStagMatGetValuesStencil.html)
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

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMSTAG/DMStagMatGetValuesStencil.html)
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

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMSTAG/DMStagMatSetValuesStencil.html)
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

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMSTAG/DMStagMatSetValuesStencil.html)
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

  
"""
    dim = DMGetDimension(dm::DMStag)

Return the topological dimension of the DM.

    dm 	- The DM
    dim - dimensions

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DM/DMGetDimension.html)
"""
function DMGetDimension end

@for_petsc function DMGetDimension(dm::DMStag)
    dim = Ref{$PetscInt}()

    @chk ccall((:DMGetDimension, $petsc_library), PetscErrorCode, (CDMStag,Ptr{$PetscInt}), dm, dim )

    return dim[]
end


"""
    DMLocalToGlobal(
        dm::DMStag,
        l::AbstractVec, 
        mode::InsertMode,
        g::AbstractVec
        )
    DMLocalToGlobal(
        dm::DMStag,
        l::AbstractVec, 
        mode::InsertMode,
        g::CVec
        )
    DMLocalToGlobal(
        dm::DMStag,
        l::CVec, 
        mode::InsertMode,
        g::CVec
        )

Updates global vectors from local vectors. 

    dm 	 - the DM object
    l 	 - the local vector
    mode - if INSERT_VALUES then no parallel communication is used, if ADD_VALUES then all ghost points from the same base point accumulate into that base point.
    g 	 - the global vector    

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DM/DMLocalToGlobal.html)
"""
function DMLocalToGlobal end

@for_petsc function DMLocalToGlobal(
    dm::DMStag{$PetscLib},
    l::AbstractVec{$PetscScalar}, 
    mode::InsertMode,
    g::AbstractVec{$PetscScalar}
    )

    DMLocalToGlobal(dm,l.ptr, mode,g.ptr)

    return nothing
end

@for_petsc function DMLocalToGlobal(
    dm::DMStag{$PetscLib},
    l::AbstractVec{$PetscScalar}, 
    mode::InsertMode,
    g::CVec
    )

    DMLocalToGlobal(dm,l.ptr, mode,g)

    return nothing
end

@for_petsc function DMLocalToGlobal(
    dm::DMStag{$PetscLib},
    l::CVec, 
    mode::InsertMode,
    g::CVec
    )

    @chk ccall((:DMLocalToGlobal, $petsc_library), PetscErrorCode,
    (
        CDMStag, 
        CVec, 
        InsertMode, 
        CVec
    ), 
    dm, 
    l, 
    mode, 
    g)


    return nothing
end


@for_petsc function DMGlobalToLocal(
    dm::DMStag{$PetscLib},
    g::AbstractVec{$PetscScalar}, 
    mode::InsertMode,
    l::AbstractVec{$PetscScalar}
    )

    DMGlobalToLocal(dm,g.ptr, mode::InsertMode,l.ptr)

    return nothing
end


@for_petsc function DMGlobalToLocal(
    dm::DMStag{$PetscLib},
    g::CVec, 
    mode::InsertMode,
    l::AbstractVec{$PetscScalar}
    )

    DMGlobalToLocal(dm,g, mode::InsertMode,l.ptr)

    return nothing
end

"""
    DMGlobalToLocal(
        dm::DMStag,
        g::CVec, 
        mode::InsertMode,
        l::CVec
        )

Update local vectors from global vector.

    dm   - the DM object
    g    - the global vector
    mode - INSERT_VALUES or ADD_VALUES
    l    - the local vector 

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DM/DMGlobalToLocal.html)
"""
function DMGlobalToLocal end

@for_petsc function DMGlobalToLocal(
    dm::DMStag{$PetscLib},
    g::CVec, 
    mode::InsertMode,
    l::CVec
    )

    @chk ccall((:DMGlobalToLocal, $petsc_library), PetscErrorCode,
    (
        CDMStag, 
        CVec, 
        InsertMode, 
        CVec
    ), 
    dm, 
    g, 
    mode, 
    l
    )

    return nothing
end


"""
    mat = DMCreateMatrix(dm::DMStag)

Generates a matrix from a DMStag object. The type is a MatSeqAIJ is we are on 1 core.

    dm  - the DMStag object
    mat - the matrix of type MatSeqAIJ (on 1 core) 

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DM/DMCreateMatrix.html)
"""
function DMCreateMatrix end

@for_petsc function DMCreateMatrix(dm::DMStag)
    # Note: the matrix cannot be viewed yet, as it remains unassembled
    #  ideally, we should modify the viewer to take care of this case

    if dm.comm==MPI.COMM_SELF
        mat = MatSeqAIJ{$PetscScalar}(C_NULL)
    elseif dm.comm==MPI.COMM_WORLD
        error("MatMPIAIJ still to be implemented")
    end

    @chk ccall((:DMCreateMatrix, $petsc_library), PetscErrorCode,
                (CDMStag, Ptr{CMat}), dm, mat)
        
    return mat
end

"""
    stencilType = DMStagGetStencilType(dm::DMStag)

Get elementwise ghost/halo stencil type.

    dm          - the DMStag object
    stencilType - the elementwise ghost stencil type: DMSTAG_STENCIL_BOX, DMSTAG_STENCIL_STAR, or DMSTAG_STENCIL_NONE

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMSTAG/DMStagGetStencilType.html)
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

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMSTAG/DMStagGetIsFirstRank.html)
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

From [PETSc Manual])(https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMSTAG/DMStagGetIsLastRank.html)
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
    dmnew = DMGetCoordinateDM(dm::DMStag; kwargs...)

Gets the DM that prescribes coordinate layout and scatters between global and local coordinates.

    dm          -    the DMStag object
    kwargs...   -   [Optional] keyword arguments (see PETSc webpage), specifiable as stag_grid_x=100, etc. 
    dmnew       -    Coordinate DM

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DM/DMGetCoordinateDM.html)
"""
function DMGetCoordinateDM end

@for_petsc function DMGetCoordinateDM(
    ::$UnionPetscLib,
    dm::DMStag; kwargs...)

    comm  = dm.comm

    opts = Options($petsclib, kwargs...)

    dmnew = DMStag{$PetscLib}(C_NULL, comm, opts)   # retrieve options

    @chk ccall((:DMGetCoordinateDM, $petsc_library), PetscErrorCode,
    (CDMStag, Ptr{CDMStag}), 
     dm, dmnew)

    return dmnew
end

"""
    v = DMGetCoordinatesLocal(dm::DMStag; write_val=true, read_val=true)

Gets a local vector with the coordinates associated with the DM.

    dm 	- the DMStag object
    v   - coordinate local vector

From [PETSc Manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DM/DMGetCoordinatesLocal.html)
"""
function DMGetCoordinatesLocal end

@for_petsc function DMGetCoordinatesLocal(
    ::$UnionPetscLib,
    dm::DMStag; write_val=true, read_val=true)

    # Should probably be a parallel vector
    v = VecSeq(C_NULL, [0.0])  # empty vector
        
    @chk ccall((:DMGetCoordinatesLocal, $petsc_library), PetscErrorCode, (CDMStag, Ptr{CVec}), dm, v)

    # Link a julia array to the values from the new vector
    # If we modify values here, it will automatically be changed in the PetcVec as well
    v.array = unsafe_localarray($PetscScalar, v.ptr; write=write_val, read=read_val)
        
    return v
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
function setuniformcoordinates!(da::DMStag,xyzmin,xyzmax) end

@for_petsc function setuniformcoordinates!(
    da::DMStag{$PetscLib},
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
            CDM,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscReal,
            $PetscReal,
        ),
        da,
        xmin,
        xmax,
        ymin,
        ymax,
        zmin,
        zmax,
    )
    return nothing
end

Base.show(io::IO, dm::DMStag) = _show(io, dm)
