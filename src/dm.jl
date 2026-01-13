import .LibPETSc: AbstractPetscDM, PetscDM, CDM

# Custom display for REPL
function Base.show(io::IO, v::AbstractPetscDM{PetscLib}) where {PetscLib}
    if v.ptr == C_NULL
        print(io, "PETSc DM (null pointer)")
        return
    else
        ty = LibPETSc.DMGetType(PetscLib,v)
        di = LibPETSc.DMGetDimension(PetscLib,v)
        print(io, "PETSc DM $ty object in $di dimensions")
    end
    return nothing
end


"""
    destroy(dm::AbstractPetscDM)

Destroy a DM object and release associated resources.

This function is typically called automatically via finalizers when the object
is garbage collected, but can be called explicitly to free resources immediately.

# External Links
$(_doc_external("DM/DMDestroy"))
"""
function destroy(dm::AbstractPetscDM{PetscLib}) where {PetscLib}
    if !(finalized(PetscLib)) && dm.ptr != C_NULL 
        LibPETSc.DMDestroy(PetscLib, dm)
    end
    dm.ptr = C_NULL
    return nothing
end



"""
    getinfo(dm::AbstractPetscDM)

Get information about a DMDA.

# Returns

A `NamedTuple` with the following fields:
- `dim`: Dimension of the `DMDA` (1, 2, or 3)
- `global_size`: Tuple with global dimensions in each direction of the array
- `mpi_proc_size`: Tuple with number of MPI processes in each direction
- `dof`: Degrees of freedom per node
- `stencil_width`: Width of the stencil
- `boundary_type`: Tuple with boundary types in each direction
- `stencil_type`: Stencil type, either `DMDA_STENCIL_STAR` or `DMDA_STENCIL_BOX`

# External Links
$(_doc_external("DMDA/DMDAGetInfo"))
"""
function getinfo(dm::AbstractPetscDM{PetscLib}) where {PetscLib}

    dim, M, N, P, m, n, p, dof, s, bx, by, bz, st = LibPETSc.DMDAGetInfo(PetscLib, dm)
    global_size   = (M,N,P)
    mpi_proc_size = (m,n,p)
    boundary_type = (bx,by,bz)
    stencil_width = s
    stencil_type  = st
               

	return (;dim,global_size,mpi_proc_size,dof,s,boundary_type,stencil_width,stencil_type)
end

"""
    lower, upper, size = getcorners_dmda(da::AbstractDMDA)

Returns a `NamedTuple` with the global indices (excluding ghost points) of the
`lower` and `upper` corners as well as the `size`.


Calls `LibPETSc.DMDAGetCorners`.
"""
function getcorners_dmda(dm::AbstractPetscDM{PetscLib}) where {PetscLib}
    PetscInt = inttype(PetscLib)
    xs, ys, zs, xm, ym, zm = LibPETSc.DMDAGetCorners(PetscLib, dm)
    corners = [PetscInt(xs), PetscInt(ys), PetscInt(zs)]
    local_size = [PetscInt(xm), PetscInt(ym), PetscInt(zm)]
    
    corners .+= 1
    upper = corners .+ local_size .- PetscInt(1)

    return (
        lower = CartesianIndex(corners...),
        upper = CartesianIndex(upper...),
        size = (local_size...,),
    )
end


"""
    lower, upper, size = getcorners(da::AbstractDMDA)

Returns a `NamedTuple` with the global indices (excluding ghost points) of the
`lower` and `upper` corners as well as the `size`. 
Works for both a DMDA and DMStag object
"""
function getcorners(dm::AbstractPetscDM{PetscLib}) where {PetscLib}
    type = gettype(dm)
    if type == "da"
        return getcorners_dmda(dm)
    elseif type == "stag"
        return getcorners_dmstag(dm)
    else
        error("getcorners only works for DMDA and DMStag objects")
    end
end

"""
    lower, upper, size = getghostcorners(da::AbstractDMDA)

Returns a `NamedTuple` with the global indices (including ghost points) of the
`lower` and `upper` corners as well as the `size`. 
Works for both a `DMDA` and `DMStag` object
"""
function getghostcorners(dm::AbstractPetscDM{PetscLib}) where {PetscLib}
    type = gettype(dm)
    if type == "da"
        return getghostcorners_dmda(dm)
    elseif type == "stag"
        return getghostcorners_dmstag(dm)
    else
        error("getghostcorners only works for DMDA and DMStag objects")
    end
end

"""
    lower, upper, size = getghostcorners_dmda(da::AbstractDMDA)

Returns a `NamedTuple` with the global indices (including ghost points) of the
`lower` and `upper` corners as well as the `size` of the local part of the domain.

Calls `LibPETSc.DMDAGetCorners`.
"""
function getghostcorners_dmda(dm::AbstractPetscDM{PetscLib}) where {PetscLib}
    PetscInt = inttype(PetscLib)
    xs, ys, zs, xm, ym, zm = LibPETSc.DMDAGetGhostCorners(PetscLib, dm)
    corners = [PetscInt(xs), PetscInt(ys), PetscInt(zs)]
    local_size = [PetscInt(xm), PetscInt(ym), PetscInt(zm)]
    
    corners .+= 1
    upper = corners .+ local_size .- PetscInt(1)

    return (
        lower = CartesianIndex(corners...),
        upper = CartesianIndex(upper...),
        size = (local_size...,),
    )
end


"""
    setup!(dm::DM)

# External Links
$(_doc_external("DM/DMSetUp"))
"""
setup!(dm::AbstractPetscDM{PetscLib}) where {PetscLib} = LibPETSc.DMSetUp(PetscLib, dm)


"""
    setfromoptions!(dm::AbstractPetscDM)

Sets the global options to the `dm`    
# External Links
$(_doc_external("DM/DMSetFromOptions"))
"""
setfromoptions!(dm::AbstractPetscDM{PetscLib}) where {PetscLib} = LibPETSc.DMSetFromOptions(PetscLib, dm)



"""
    empty(da::AbstractPetscDM)

return an uninitialized `DMDA` struct.
"""
Base.empty(da::AbstractPetscDM{PetscLib}) where {PetscLib} = PetscDM{PetscLib}(C_NULL, da.age)


"""
    v::PetscVec = DMLocalVec(dm::AbstractPetscDM{PetscLib}) where {PetscLib}

Returns a local vector `v` from the `dm` object.
"""
DMLocalVec(dm::AbstractPetscDM{PetscLib}) where {PetscLib} = LibPETSc.DMCreateLocalVector(getlib(PetscLib), dm)

"""
    v::PetscVec = DMGlobalVec(dm::AbstractPetscDM{PetscLib}) where {PetscLib}

Returns a global vector `v` from the `dm` object.
"""
DMGlobalVec(dm::AbstractPetscDM{PetscLib}) where {PetscLib} = LibPETSc.DMCreateGlobalVector(getlib(PetscLib), dm)

"""
    dm_local_to_global!(local_vec, global_vec, dm, mode = INSERT_VALUES)

Transfer values from the `local_vec` to the `global_vec` associated with the `dm` object.

# Arguments
- `local_vec::AbstractPetscVec`: Local vector (source)
- `global_vec::AbstractPetscVec`: Global vector (destination)
- `dm::AbstractPetscDM`: DM object
- `mode::InsertMode`: Insert mode, either `INSERT_VALUES` or `ADD_VALUES`

# External Links
$(_doc_external("DM/DMLocalToGlobal"))

"""
function dm_local_to_global!(   local_vec::AbstractPetscVec{PetscLib},
                                global_vec::AbstractPetscVec{PetscLib},
                                   dm::AbstractPetscDM{PetscLib},
                                   mode::InsertMode = INSERT_VALUES) where {PetscLib}

    LibPETSc.DMLocalToGlobalBegin(PetscLib, dm, local_vec, mode, global_vec)
    LibPETSc.DMLocalToGlobalEnd(PetscLib, dm, local_vec,  mode, global_vec)
    return nothing
end


"""
    dm_global_to_local!(global_vec, local_vec, dm, mode = INSERT_VALUES)

Transfer values from the `global_vec` to the `local_vec` associated with the `dm` object,
including ghost point values from neighboring processes.

# Arguments
- `global_vec::AbstractPetscVec`: Global vector (source)
- `local_vec::AbstractPetscVec`: Local vector (destination)
- `dm::AbstractPetscDM`: DM object
- `mode::InsertMode`: Insert mode, either `INSERT_VALUES` or `ADD_VALUES`

# External Links
$(_doc_external("DM/DMLocalToGlobal"))
"""
function dm_global_to_local!(global_vec::AbstractPetscVec{PetscLib},
                              local_vec::AbstractPetscVec{PetscLib},
                                     dm::AbstractPetscDM{PetscLib},
                                   mode::InsertMode = INSERT_VALUES) where {PetscLib}

    LibPETSc.DMGlobalToLocalBegin(getlib(PetscLib), dm, global_vec, mode, local_vec)
    LibPETSc.DMGlobalToLocalEnd(getlib(PetscLib), dm, global_vec,  mode, local_vec)
    return nothing
end


"""
    setuniformcoordinates!(
        da::DMDA
        xyzmin::NTuple{N, Real},
        xyzmax::NTuple{N, Real},
    ) where {N}

Set uniform coordinates for the `da` using the lower and upper corners defined
by the `NTuple`s `xyzmin` and `xyzmax`. If `N` is less than the dimension of the
`da` then the value of the trailing coordinates is set to `0`.

# External Links
$(_doc_external("DMDA/DMDASetUniformCoordinates"))
"""
function setuniformcoordinates_dmda!(
    da::PetscDM{PetscLib},
    xyzmin::NTuple{N, Real},
    xyzmax::NTuple{N, Real},
) where {N, PetscLib}
    @assert gettype(da) == "da" "setuniformcoordinates_dmda! only works for DMDA objects"
    PetscReal = PetscLib.PetscReal
    xmin = PetscReal(xyzmin[1])
    xmax = PetscReal(xyzmax[1])

    ymin = (N > 1) ? PetscReal(xyzmin[2]) : PetscReal(0)
    ymax = (N > 1) ? PetscReal(xyzmax[2]) : PetscReal(0)

    zmin = (N > 2) ? PetscReal(xyzmin[3]) : PetscReal(0)
    zmax = (N > 2) ? PetscReal(xyzmax[3]) : PetscReal(0)


    LibPETSc.DMDASetUniformCoordinates(
        PetscLib,
        da,
        xmin,
        xmax,
        ymin,
        ymax,
        zmin,
        zmax,
    )
    return da
end

"""
    coordinatesDMLocalVec(dm::AbstractDM)

Gets a local vector with the coordinates associated with `dm`.

Note that the returned vector is borrowed from the `dm` and is not a new vector.

# External Links
$(_doc_external("DM/DMGetCoordinatesLocal"))
"""
function coordinatesDMLocalVec(dm::AbstractPetscDM{PetscLib}) where {PetscLib}
    petsclib = getlib(PetscLib)
    coord_vec = DMLocalVec(dm)
    LibPETSc.DMGetCoordinatesLocal(PetscLib, dm, coord_vec)

    return coord_vec
end

"""
    getlocalcoordinatearray(da::AbstractPetscDM)

Return coordinate arrays for the local portion of the domain.

The returned arrays are `OffsetArray`s that can be addressed using global indices,
accounting for ghost points.

# External Links
$(_doc_external("DM/DMGetCoordinatesLocal"))
"""
function getlocalcoordinatearray(da::AbstractPetscDM{PetscLib}) where {PetscLib}
    # retrieve local coordinates
    coord_vec = coordinatesDMLocalVec(da)
    # array
    array1D = unsafe_localarray(coord_vec; read = true, write = false)
    dim = [PetscLib.PetscInt(0)]
    dim = LibPETSc.DMGetCoordinateDim(PetscLib, da)
    dim = dim[1]
    corners = getghostcorners(da)

    return reshapelocalarray(array1D, da, dim)
end


gettype(dm::PetscDM{PetscLib}) where {PetscLib} = LibPETSc.DMGetType(PetscLib,dm)

"""
    getdimension(dm::AbstractPetscDM)

Return the topological dimension of the `dm`

# External Links
$(_doc_external("DM/DMGetDimension"))
"""
getdimension(dm::AbstractPetscDM{PetscLib}) where PetscLib = LibPETSc.DMGetDimension(PetscLib,dm)


"""
    size(dm::AbstractPetscDM)

Return the global size of a DM object as a tuple.

For DMDA and DMStag, returns `(M, N, P)` where unused dimensions are 1.
"""
function Base.size(dm::AbstractPetscDM{PetscLib}) where PetscLib
    if gettype(dm) == "stag"
        size = LibPETSc.DMStagGetGlobalSizes(PetscLib,dm)
    elseif gettype(dm) == "da"
        dim, M,N,P,_ = LibPETSc.DMDAGetInfo(PetscLib, dm)
        size = (M,N,P)
    else
        error("Size not defined for DMStag objects. Use getinfo(dm).global_size instead.")
    end
    return size
end


"""
    dm_local_to_global(dm, x_L, x_G, mode = INSERT_VALUES)

Transfer values from the local vector `x_L` to the global vector `x_G`.

# Arguments
- `dm`: The DM object
- `x_L`: Local vector (source)
- `x_G`: Global vector (destination)
- `mode`: `INSERT_VALUES` (default) or `ADD_VALUES`

# External Links
$(_doc_external("DM/DMLocalToGlobal"))
"""
function dm_local_to_global(dm::PetscDM{PetscLib},
                             x_L::AbstractPetscVec{PetscLib},
                             x_G::AbstractPetscVec{PetscLib}, 
                             mode=LibPETSc.INSERT_VALUES) where {PetscLib}
    
    petsclib = getlib(PetscLib)
    LibPETSc.DMLocalToGlobalBegin(petsclib, dm, x_L, mode, x_G)
    LibPETSc.DMLocalToGlobalEnd(petsclib, dm, x_L, mode, x_G)
    
    return nothing
end

"""
    dm_global_to_local(dm, x_G, x_L, mode = INSERT_VALUES)

Transfer values from the global vector `x_G` to the local vector `x_L`,
including ghost point values from neighboring processes.

# Arguments
- `dm`: The DM object
- `x_G`: Global vector (source)
- `x_L`: Local vector (destination)
- `mode`: `INSERT_VALUES` (default) or `ADD_VALUES`

# External Links
$(_doc_external("DM/DMGlobalToLocal"))
"""
function dm_global_to_local(dm::PetscDM{PetscLib},
                             x_G::AbstractPetscVec{PetscLib},
                             x_L::AbstractPetscVec{PetscLib}, 
                             mode=LibPETSc.INSERT_VALUES) where {PetscLib}
    
    petsclib = getlib(PetscLib)
    LibPETSc.DMGlobalToLocalBegin(petsclib, dm, x_G, mode, x_L)
    LibPETSc.DMGlobalToLocalEnd(petsclib, dm, x_G, mode, x_L)

    return nothing
end

"""
    MatAIJ(da::AbstractPetscDM)

Create a sparse matrix (AIJ format) with sparsity pattern determined by the DM.

# Returns

A `PetscMat` object compatible with vectors from the DM.

# External Links
$(_doc_external("DM/DMCreateMatrix"))
"""
function MatAIJ(da::AbstractPetscDM{PetscLib}) where {PetscLib}
    J = LibPETSc.DMCreateMatrix(getlib(PetscLib), da)
    return J
end
