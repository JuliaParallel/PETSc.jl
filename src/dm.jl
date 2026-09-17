import .LibPETSc: AbstractPetscDM, PetscDM, CDM

# ─────────────────────────────────────────────────────────────────────────────
# The typed DM hierarchy (docs/src/man/naming.md §5.3)
#
# `PetscDM{PetscLib}` stays the low-level handle, because `LibPETSc.DMCreate`
# and friends must return something before the flavour is known (§5.4). The
# high-level API uses the three types below, which carry the flavour — and, for
# DMDA and DMStag, the dimension — in the type, so that `corners`,
# `local_indices` and the coordinate handling dispatch instead of comparing the
# string `DMGetType` returns. Giving flavour a type is what deletes the nine
# `@assert type_name(dm) == "stag"` checks (§14).
#
# DMPlex has no dimension parameter: nothing in dmplex.jl dispatches on it, and
# neither constructor could supply one honestly, since `DMPlex(petsclib, comm)`
# leaves the dimension unset until setup (§5.3).
#
# `own` is the ownership flag of §3.3: `true` for a constructor result, `false`
# for a borrowed handle produced by `narrow`.
# ─────────────────────────────────────────────────────────────────────────────

"""
    DMDA{PetscLib, N}

A `DMDA`: a `N`-dimensional PETSc distributed array.

Build one with [`DMDA(petsclib, comm, boundary_type, global_dim, dof_per_node,
stencil_width, stencil_type)`](@ref). A handle onto a DMDA that PETSc owns comes
from [`narrow`](@ref).

# External Links
$(doc_external("DMDA/DMDA"))
"""
mutable struct DMDA{PetscLib, N} <: AbstractPetscDM{PetscLib}
    ptr::CDM
    age::Int
    own::Bool
end

"""
    DMStag{PetscLib, N}

A `DMSTAG`: a `N`-dimensional PETSc staggered-grid distributed array.

Build one with [`DMStag(petsclib, comm, boundary_type, global_dim, dof_per_node,
stencil_width, stencil_type)`](@ref). A handle onto a DMStag that PETSc owns
comes from [`narrow`](@ref).

# External Links
$(doc_external("DMSTAG/DMSTAG"))
"""
mutable struct DMStag{PetscLib, N} <: AbstractPetscDM{PetscLib}
    ptr::CDM
    age::Int
    own::Bool
end

"""
    DMPlex{PetscLib}

A `DMPLEX`: an unstructured PETSc mesh.

Dimension is a runtime property here rather than a type parameter, because
nothing dispatches on it and `DMPlex(petsclib, comm)` leaves it unset until
setup. Ask for it with `ndims(dm)`.

# External Links
$(doc_external("DMPLEX/DMPLEX"))
"""
mutable struct DMPlex{PetscLib} <: AbstractPetscDM{PetscLib}
    ptr::CDM
    age::Int
    own::Bool
end

const TypedPetscDM{PetscLib} =
    Union{DMDA{PetscLib}, DMStag{PetscLib}, DMPlex{PetscLib}}

owns(dm::TypedPetscDM) = dm.own

# Wrap a handle a LibPETSc creator returned, taking ownership of it. The
# finalizer is only attached on a serial communicator: `DMDestroy` is
# collective and a GC finalizer runs at an arbitrary point.
function own_dm!(dm::TypedPetscDM{PetscLib}, comm) where {PetscLib}
    if MPI.Comm_size(comm) == 1
        finalizer(destroy!, dm)
    end
    return dm
end

"""
    narrow(dm::AbstractPetscDM; own = false)

A handle onto the same PETSc object whose Julia type carries the DM's flavour
and, for a `DMDA` or a `DMStag`, its dimension.

`LibPETSc.DMGetType` and `LibPETSc.DMGetDimension` are queried, so the return
type is a wide `Union` and the call is a dynamic dispatch. One dispatch is
cheap; propagating an abstractly-typed DM through a hot loop is not, so narrow
once behind a function barrier. A flavour with no type of its own comes back
unchanged.

$(doc_borrowed())

The exception is `own = true`, which the callers that really do hand over a new
object use: [`clone`](@ref) and [`distribute!`](@ref).

# External Links
$(doc_external("DM/DMGetType"))
$(doc_external("DM/DMGetDimension"))
"""
function narrow(dm::AbstractPetscDM{PetscLib}; own::Bool = false) where {PetscLib}
    dm.ptr == C_NULL && return dm
    tname = LibPETSc.DMGetType(PetscLib, dm)
    if tname == "da"
        N = Int(LibPETSc.DMGetDimension(PetscLib, dm))
        out = DMDA{PetscLib, N}(dm.ptr, dm.age, own)
    elseif tname == "stag"
        N = Int(LibPETSc.DMGetDimension(PetscLib, dm))
        out = DMStag{PetscLib, N}(dm.ptr, dm.age, own)
    elseif tname == "plex"
        out = DMPlex{PetscLib}(dm.ptr, dm.age, own)
    else
        return dm
    end
    own && own_dm!(out, LibPETSc.PetscObjectGetComm(getlib(PetscLib), out))
    return out
end


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
Does nothing on a DM that only borrows its handle: see [`owns`](@ref).

# External Links
$(doc_external("DM/DMDestroy"))
"""
function destroy!(dm::AbstractPetscDM{PetscLib}) where {PetscLib}
    owns(dm) || return nothing
    if isdestroyable(dm, PetscLib)
        LibPETSc.DMDestroy(PetscLib, dm)
    end
    dm.ptr = C_NULL
    return nothing
end



"""
    info(dm::DMDA)

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
function info(dm::DMDA{PetscLib}) where {PetscLib}

    dim, M, N, P, m, n, p, dof, s, bx, by, bz, st = LibPETSc.DMDAGetInfo(PetscLib, dm)
    global_size   = (M,N,P)
    mpi_proc_size = (m,n,p)
    boundary_type = (bx,by,bz)
    stencil_width = s
    stencil_type  = st
               

	return (;dim,global_size,mpi_proc_size,dof,s,boundary_type,stencil_width,stencil_type)
end

"""
    lower, upper, size = corners(da::DMDA)

Returns a `NamedTuple` with the global indices (excluding ghost points) of the
`lower` and `upper` corners as well as the `size`.

# External Links
$(doc_external("DMDA/DMDAGetCorners"))
"""
function corners(dm::DMDA{PetscLib}) where {PetscLib}
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
    lower, upper, size = ghost_corners(da::DMDA)

Returns a `NamedTuple` with the global indices (including ghost points) of the
`lower` and `upper` corners as well as the `size` of the local part of the domain.

# External Links
$(doc_external("DMDA/DMDAGetGhostCorners"))
"""
function ghost_corners(dm::DMDA{PetscLib}) where {PetscLib}
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
    corners(dm)

Returns a `NamedTuple` with the global indices (excluding ghost points) of the
`lower` and `upper` corners as well as the `size`. A `DMStag` also reports
`nextra`, the number of extra partial elements in each direction.

Defined for [`DMDA`](@ref) and [`DMStag`](@ref); the flavour is a type
parameter, so any other DM is a `MethodError` rather than a runtime check.
"""
function corners end

"""
    ghost_corners(dm)

Returns a `NamedTuple` with the global indices (including ghost points) of the
`lower` and `upper` corners as well as the `size`.

Defined for [`DMDA`](@ref) and [`DMStag`](@ref).
"""
function ghost_corners end

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

Defined for [`DMDA`](@ref) and [`DMStag`](@ref), which reach different PETSc
calls; the flavour is a type parameter, so the two are ordinary methods.

# External Links
$(doc_external("DMDA/DMDASetUniformCoordinates"))
$(doc_external("DMSTAG/DMStagSetUniformCoordinatesProduct"))
"""
function set_uniform_coordinates! end

"""
    set_uniform_coordinates!(da::DMDA, xyzmin, xyzmax)

The [`DMDA`](@ref) method of [`set_uniform_coordinates!`](@ref).

# External Links
$(doc_external("DMDA/DMDASetUniformCoordinates"))
"""
function set_uniform_coordinates!(
    da::DMDA{PetscLib},
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

$(doc_borrowed())

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
    local_coordinate_array(da::Union{DMDA, DMStag})

Return coordinate arrays for the local portion of the domain.

The returned arrays are `OffsetArray`s that can be addressed using global indices,
accounting for ghost points.

# External Links
$(doc_external("DM/DMGetCoordinatesLocal"))
"""
function local_coordinate_array(
    da::Union{DMDA{PetscLib}, DMStag{PetscLib}},
) where {PetscLib}
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
    type_name(dm::AbstractPetscDM)

The name PETSc knows this DM's flavour by, as a `String` (`"da"`, `"stag"`, `"plex"`, …).

# External Links
$(doc_external("DM/DMGetType"))
"""
type_name(dm::AbstractPetscDM{PetscLib}) where {PetscLib} =
    LibPETSc.DMGetType(PetscLib, dm)

"""
    ndims(dm::AbstractPetscDM)

Return the topological dimension of the `dm`

# External Links
$(doc_external("DM/DMGetDimension"))
"""
Base.ndims(dm::AbstractPetscDM{PetscLib}) where {PetscLib} =
    LibPETSc.DMGetDimension(PetscLib, dm)


"""
    size(dm::DMDA)
    size(dm::DMStag)

Return the global size of the DM as a tuple.

Returns `(M, N, P)`, with unused dimensions set to `1`.

# External Links
$(doc_external("DMDA/DMDAGetInfo"))
$(doc_external("DMSTAG/DMStagGetGlobalSizes"))
"""
function Base.size(dm::DMDA{PetscLib}) where {PetscLib}
    _, M, N, P, _ = LibPETSc.DMDAGetInfo(PetscLib, dm)
    return (M, N, P)
end

Base.size(dm::DMStag{PetscLib}) where {PetscLib} =
    LibPETSc.DMStagGetGlobalSizes(PetscLib, dm)

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
    PetscMat(da::AbstractPetscDM)

Create a sparse matrix (AIJ format) with sparsity pattern determined by the DM.

Replaces v0.4's `MatAIJ`: construction goes through the type
(docs/src/man/naming.md §5.1).

# Returns

A `PetscMat` object compatible with vectors from the DM.

# External Links
$(doc_external("DM/DMCreateMatrix"))
"""
function LibPETSc.PetscMat(da::AbstractPetscDM{PetscLib}) where {PetscLib}
    J = LibPETSc.DMCreateMatrix(getlib(PetscLib), da)
    return J
end

