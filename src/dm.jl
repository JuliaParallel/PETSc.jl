import .LibPETSc: AbstractPetscDM, PetscDM, CDM

# Custom display for REPL
function Base.show(io::IO, v::AbstractPetscDM{PetscLib}) where {PetscLib}
    if v.ptr == C_NULL
        print(io, "PETSc DM (null pointer)")
        return
    end
    # DMGetType internally calls DMInitializePackage which queries the PETSc
    # options database.  Calling it before PETSc is initialised causes a C-level
    # SIGSEGV that cannot be caught with try/catch.
    if !isinitialized(PetscLib)
        print(io, "PETSc DM (PETSc not initialized)")
        return
    end
    try
        ty = LibPETSc.DMGetType(PetscLib, v)
        di = LibPETSc.DMGetDimension(PetscLib, v)
        print(io, "PETSc DM $ty object in $di dimensions")
    catch
        print(io, "PETSc DM (type not set)")
    end
    return nothing
end


"""
    destroy!(dm::AbstractPetscDM)

Destroy a DM object and release associated resources.

This function is typically called automatically via finalizers when the object
is garbage collected, but can be called explicitly to free resources immediately.

# External Links
$(doc_external("DM/DMDestroy"))
"""
function destroy!(dm::AbstractPetscDM{PetscLib}) where {PetscLib}
    if isdestroyable(dm, PetscLib)
        LibPETSc.DMDestroy(PetscLib, dm)
    end
    dm.ptr = C_NULL
    return nothing
end



"""
    info(dm::AbstractPetscDM)

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
$(doc_external("DMDA/DMDAGetInfo"))
"""
function info(dm::AbstractPetscDM{PetscLib}) where {PetscLib}

    dim, M, N, P, m, n, p, dof, s, bx, by, bz, st = LibPETSc.DMDAGetInfo(PetscLib, dm)
    global_size   = (M,N,P)
    mpi_proc_size = (m,n,p)
    boundary_type = (bx,by,bz)
    stencil_width = s
    stencil_type  = st
               

	return (;dim,global_size,mpi_proc_size,dof,s,boundary_type,stencil_width,stencil_type)
end

"""
    lower, upper, size = corners_dmda(da::AbstractDMDA)

Returns a `NamedTuple` with the global indices (excluding ghost points) of the
`lower` and `upper` corners as well as the `size`.


Calls `LibPETSc.DMDAGetCorners`.
"""
function corners_dmda(dm::AbstractPetscDM{PetscLib}) where {PetscLib}
    PetscInt = inttype(PetscLib)
    xs, ys, zs, xm, ym, zm = LibPETSc.DMDAGetCorners(PetscLib, dm)
    lo = [PetscInt(xs), PetscInt(ys), PetscInt(zs)]
    local_size = [PetscInt(xm), PetscInt(ym), PetscInt(zm)]

    lo .+= 1
    upper = lo .+ local_size .- PetscInt(1)

    return (
        lower = CartesianIndex(lo...),
        upper = CartesianIndex(upper...),
        size = (local_size...,),
    )
end


"""
    lower, upper, size = corners(da::AbstractDMDA)

Returns a `NamedTuple` with the global indices (excluding ghost points) of the
`lower` and `upper` corners as well as the `size`. 
Works for both a DMDA and DMStag object
"""
function corners(dm::AbstractPetscDM{PetscLib}) where {PetscLib}
    tname = type_name(dm)
    if tname == "da"
        return corners_dmda(dm)
    elseif tname == "stag"
        return corners_dmstag(dm)
    else
        error("corners only works for DMDA and DMStag objects")
    end
end

"""
    lower, upper, size = ghost_corners(da::AbstractDMDA)

Returns a `NamedTuple` with the global indices (including ghost points) of the
`lower` and `upper` corners as well as the `size`. 
Works for both a `DMDA` and `DMStag` object
"""
function ghost_corners(dm::AbstractPetscDM{PetscLib}) where {PetscLib}
    tname = type_name(dm)
    if tname == "da"
        return ghost_corners_dmda(dm)
    elseif tname == "stag"
        return ghost_corners_dmstag(dm)
    else
        error("ghost_corners only works for DMDA and DMStag objects")
    end
end

"""
    lower, upper, size = ghost_corners_dmda(da::AbstractDMDA)

Returns a `NamedTuple` with the global indices (including ghost points) of the
`lower` and `upper` corners as well as the `size` of the local part of the domain.

Calls `LibPETSc.DMDAGetCorners`.
"""
function ghost_corners_dmda(dm::AbstractPetscDM{PetscLib}) where {PetscLib}
    PetscInt = inttype(PetscLib)
    xs, ys, zs, xm, ym, zm = LibPETSc.DMDAGetGhostCorners(PetscLib, dm)
    lo = [PetscInt(xs), PetscInt(ys), PetscInt(zs)]
    local_size = [PetscInt(xm), PetscInt(ym), PetscInt(zm)]

    lo .+= 1
    upper = lo .+ local_size .- PetscInt(1)

    return (
        lower = CartesianIndex(lo...),
        upper = CartesianIndex(upper...),
        size = (local_size...,),
    )
end


"""
    setup!(dm::DM)

# External Links
$(doc_external("DM/DMSetUp"))
"""
setup!(dm::AbstractPetscDM{PetscLib}) where {PetscLib} = LibPETSc.DMSetUp(PetscLib, dm)


"""
    set_from_options!(dm::AbstractPetscDM)

Sets the global options to the `dm`    
# External Links
$(doc_external("DM/DMSetFromOptions"))
"""
set_from_options!(dm::AbstractPetscDM{PetscLib}) where {PetscLib} = LibPETSc.DMSetFromOptions(PetscLib, dm)



"""
    v::PetscVec = local_vec(dm::AbstractPetscDM{PetscLib}) where {PetscLib}

Returns a local vector `v` from the `dm` object.
"""
local_vec(dm::AbstractPetscDM{PetscLib}) where {PetscLib} = LibPETSc.DMCreateLocalVector(getlib(PetscLib), dm)

"""
    v::PetscVec = global_vec(dm::AbstractPetscDM{PetscLib}) where {PetscLib}

Returns a global vector `v` from the `dm` object.
"""
global_vec(dm::AbstractPetscDM{PetscLib}) where {PetscLib} = LibPETSc.DMCreateGlobalVector(getlib(PetscLib), dm)

"""
    local_to_global!(local_vec, global_vec, dm, mode = INSERT_VALUES)

Transfer values from the `local_vec` to the `global_vec` associated with the `dm` object.

# Arguments
- `local_vec::AbstractPetscVec`: Local vector (source)
- `global_vec::AbstractPetscVec`: Global vector (destination)
- `dm::AbstractPetscDM`: DM object
- `mode::InsertMode`: Insert mode, either `INSERT_VALUES` or `ADD_VALUES`

# External Links
$(doc_external("DM/DMLocalToGlobal"))

"""
function local_to_global!(   local_vec::AbstractPetscVec{PetscLib},
                                global_vec::AbstractPetscVec{PetscLib},
                                   dm::AbstractPetscDM{PetscLib},
                                   mode::InsertMode = INSERT_VALUES) where {PetscLib}

    LibPETSc.DMLocalToGlobalBegin(PetscLib, dm, local_vec, mode, global_vec)
    LibPETSc.DMLocalToGlobalEnd(PetscLib, dm, local_vec,  mode, global_vec)
    return nothing
end


"""
    global_to_local!(global_vec, local_vec, dm, mode = INSERT_VALUES)

Transfer values from the `global_vec` to the `local_vec` associated with the `dm` object,
including ghost point values from neighboring processes.

# Arguments
- `global_vec::AbstractPetscVec`: Global vector (source)
- `local_vec::AbstractPetscVec`: Local vector (destination)
- `dm::AbstractPetscDM`: DM object
- `mode::InsertMode`: Insert mode, either `INSERT_VALUES` or `ADD_VALUES`

# External Links
$(doc_external("DM/DMLocalToGlobal"))
"""
function global_to_local!(global_vec::AbstractPetscVec{PetscLib},
                              local_vec::AbstractPetscVec{PetscLib},
                                     dm::AbstractPetscDM{PetscLib},
                                   mode::InsertMode = INSERT_VALUES) where {PetscLib}

    LibPETSc.DMGlobalToLocalBegin(getlib(PetscLib), dm, global_vec, mode, local_vec)
    LibPETSc.DMGlobalToLocalEnd(getlib(PetscLib), dm, global_vec,  mode, local_vec)
    return nothing
end


"""
    set_uniform_coordinates!(
        dm::AbstractPetscDM,
        xyzmin::NTuple{N, Real},
        xyzmax::NTuple{N, Real},
    ) where {N}

Set uniform coordinates on `dm` using the lower and upper corners defined by the
`NTuple`s `xyzmin` and `xyzmax`. If `N` is less than the dimension of the `dm`
then the value of the trailing coordinates is set to `0`.

Works for both a `DMDA` and a `DMStag`. The flavour is still resolved at runtime
from `type_name(dm)`, because both are the same concrete type until the typed DM
hierarchy lands.

# External Links
$(doc_external("DMDA/DMDASetUniformCoordinates"))
$(doc_external("DMSTAG/DMStagSetUniformCoordinatesProduct"))
"""
function set_uniform_coordinates!(
    dm::AbstractPetscDM{PetscLib},
    xyzmin::NTuple,
    xyzmax::NTuple,
) where {PetscLib}
    tname = type_name(dm)
    if tname == "da"
        return set_uniform_coordinates_dmda!(dm, xyzmin, xyzmax)
    elseif tname == "stag"
        return set_uniform_coordinates_stag!(dm, xyzmin, xyzmax)
    else
        throw(
            ArgumentError(
                "set_uniform_coordinates! only works for DMDA and DMStag objects, " *
                "got a DM of type \"$tname\"",
            ),
        )
    end
end

"""
    set_uniform_coordinates_dmda!(da, xyzmin, xyzmax)

The `DMDA` method behind [`set_uniform_coordinates!`](@ref).

# External Links
$(doc_external("DMDA/DMDASetUniformCoordinates"))
"""
function set_uniform_coordinates_dmda!(
    da::AbstractPetscDM{PetscLib},
    xyzmin::NTuple{N, Real},
    xyzmax::NTuple{N, Real},
) where {N, PetscLib}
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
    local_coordinates(dm::AbstractDM)

Gets a local vector with the coordinates associated with `dm`.

Note that the returned vector is borrowed from the `dm` and is not a new vector.

# External Links
$(doc_external("DM/DMGetCoordinatesLocal"))
"""
function local_coordinates(dm::AbstractPetscDM{PetscLib}) where {PetscLib}
    petsclib = getlib(PetscLib)
    coord_vec = LibPETSc.DMGetCoordinatesLocal(petsclib, dm)
    # borrowed from the DM: `destroy!` on the returned handle is a no-op
    return VecPtr(petsclib, coord_vec.ptr, false)
end

"""
    local_coordinate_array(da::AbstractPetscDM)

Return coordinate arrays for the local portion of the domain.

The returned arrays are `OffsetArray`s that can be addressed using global indices,
accounting for ghost points.

# External Links
$(doc_external("DM/DMGetCoordinatesLocal"))
"""
function local_coordinate_array(da::AbstractPetscDM{PetscLib}) where {PetscLib}
    # retrieve local coordinates
    coord_vec = local_coordinates(da)
    # array
    array1D = unsafe_local_array(coord_vec; read = true, write = false)
    dim = [PetscLib.PetscInt(0)]
    dim = LibPETSc.DMGetCoordinateDim(PetscLib, da)
    dim = dim[1]

    return reshape_local_array(array1D, da, dim)
end


"""
    type_name(dm::PetscDM)

The name PETSc knows this DM's flavour by, as a `String` (`"da"`, `"stag"`, `"plex"`, …).

# External Links
$(doc_external("DM/DMGetType"))
"""
type_name(dm::PetscDM{PetscLib}) where {PetscLib} = LibPETSc.DMGetType(PetscLib, dm)

"""
    ndims(dm::AbstractPetscDM)

Return the topological dimension of the `dm`

# External Links
$(doc_external("DM/DMGetDimension"))
"""
Base.ndims(dm::AbstractPetscDM{PetscLib}) where {PetscLib} =
    LibPETSc.DMGetDimension(PetscLib, dm)


"""
    size(dm::AbstractPetscDM)

Return the global size of a DM object as a tuple.

For DMDA and DMStag, returns `(M, N, P)` where unused dimensions are 1.
"""
function Base.size(dm::AbstractPetscDM{PetscLib}) where PetscLib
    tname = type_name(dm)
    if tname == "stag"
        sz = LibPETSc.DMStagGetGlobalSizes(PetscLib, dm)
    elseif tname == "da"
        dim, M, N, P, _ = LibPETSc.DMDAGetInfo(PetscLib, dm)
        sz = (M, N, P)
    else
        throw(
            ArgumentError(
                "size is only defined for DMDA and DMStag objects, got a DM of type \"$tname\"",
            ),
        )
    end
    return sz
end

#=
"""
    dm_local_to_global(dm, x_L, x_G, mode = INSERT_VALUES)

Transfer values from the local vector `x_L` to the global vector `x_G`.

# Arguments
- `dm`: The DM object
- `x_L`: Local vector (source)
- `x_G`: Global vector (destination)
- `mode`: `INSERT_VALUES` (default) or `ADD_VALUES`

# External Links
$(doc_external("DM/DMLocalToGlobal"))
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
=#
#=
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
$(doc_external("DM/DMGlobalToLocal"))
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
=#

"""
    MatAIJ(da::AbstractPetscDM)

Create a sparse matrix (AIJ format) with sparsity pattern determined by the DM.

# Returns

A `PetscMat` object compatible with vectors from the DM.

# External Links
$(doc_external("DM/DMCreateMatrix"))
"""
function MatAIJ(da::AbstractPetscDM{PetscLib}) where {PetscLib}
    J = LibPETSc.DMCreateMatrix(getlib(PetscLib), da)
    return J
end

