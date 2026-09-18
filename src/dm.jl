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

# ─────────────────────────────────────────────────────────────────────────────
# Dimension-correct returns (docs/src/man/naming.md §12)
#
# Every index-shaped return is built from the DM's type parameter `N` with
# `ntuple(..., Val(N))`, never by splatting a `Vector`: splatting hides the
# length from the compiler, so `lower` and `upper` used to infer as `Any` and
# cost three allocations per call.
#
# `axis_names` is three methods rather than `(:x, :y, :z)[1:N]` so the tuple is
# a compile-time constant and `local_indices(dm2d).center` infers concretely as
# `@NamedTuple{x::UnitRange{Int}, y::UnitRange{Int}}`.
# ─────────────────────────────────────────────────────────────────────────────

axis_names(::Val{1}) = (:x,)
axis_names(::Val{2}) = (:x, :y)
axis_names(::Val{3}) = (:x, :y, :z)

@inline axis_ranges(f::F, ::Val{N}) where {F, N} =
    NamedTuple{axis_names(Val(N))}(ntuple(f, Val(N)))

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

A `NamedTuple` with the following fields, all of them `N`-dimensional for an
`N`-dimensional `DMDA` (docs/src/man/naming.md §12):
- `dim`: Dimension of the `DMDA` (1, 2, or 3)
- `global_size`: `NTuple{N,Int}` with the global dimensions in each direction
- `procs`: `NTuple{N,Int}` with the number of MPI processes in each direction
- `ndofs`: Degrees of freedom per node
- `stencil_width`: Width of the stencil
- `boundary_type`: `NTuple{N,DMBoundaryType}` with the boundary type per direction
- `stencil_type`: Stencil type, either `DMDA_STENCIL_STAR` or `DMDA_STENCIL_BOX`

v0.4 reported the stencil width twice, as `s` and as `stencil_width`, spelled
the dof count `dof` and the process grid `mpi_proc_size`, and padded every
tuple to three entries. All four are gone; this is a break with no shim (§16).

# External Links
$(doc_external("DMDA/DMDAGetInfo"))
"""
function info(dm::DMDA{PetscLib, N}) where {PetscLib, N}
    dim, M, Ny, P, m, n, p, dof, s, bx, by, bz, st =
        LibPETSc.DMDAGetInfo(PetscLib, dm)
    gsize = (Int(M), Int(Ny), Int(P))
    nproc = (Int(m), Int(n), Int(p))
    btype = (bx, by, bz)

    return (;
        dim = Int(dim),
        global_size = ntuple(i -> gsize[i], Val(N)),
        procs = ntuple(i -> nproc[i], Val(N)),
        ndofs = Int(dof),
        stencil_width = Int(s),
        boundary_type = ntuple(i -> btype[i], Val(N)),
        stencil_type = st,
    )
end

"""
    lower, upper, size = corners(da::DMDA{PetscLib, N})

Returns a `NamedTuple` with the global indices (excluding ghost points) of the
`lower` and `upper` corners as well as the `size`.

The result is dimension-correct (§12): `lower` and `upper` are
`CartesianIndex{N}` and `size` is an `NTuple{N,Int}`, with no padding to three
entries. This is a break with no shim (§16).

# External Links
$(doc_external("DMDA/DMDAGetCorners"))
"""
function corners(dm::DMDA{PetscLib, N}) where {PetscLib, N}
    xs, ys, zs, xm, ym, zm = LibPETSc.DMDAGetCorners(PetscLib, dm)
    lo = (Int(xs), Int(ys), Int(zs))
    sz = (Int(xm), Int(ym), Int(zm))

    return (
        lower = CartesianIndex(ntuple(i -> lo[i] + 1, Val(N))),
        upper = CartesianIndex(ntuple(i -> lo[i] + sz[i], Val(N))),
        size = ntuple(i -> sz[i], Val(N)),
    )
end

"""
    lower, upper, size = ghost_corners(da::DMDA{PetscLib, N})

Returns a `NamedTuple` with the global indices (including ghost points) of the
`lower` and `upper` corners as well as the `size` of the local part of the domain.

Dimension-correct like [`corners`](@ref): `CartesianIndex{N}` and `NTuple{N,Int}`.

# External Links
$(doc_external("DMDA/DMDAGetGhostCorners"))
"""
function ghost_corners(dm::DMDA{PetscLib, N}) where {PetscLib, N}
    xs, ys, zs, xm, ym, zm = LibPETSc.DMDAGetGhostCorners(PetscLib, dm)
    lo = (Int(xs), Int(ys), Int(zs))
    sz = (Int(xm), Int(ym), Int(zm))

    return (
        lower = CartesianIndex(ntuple(i -> lo[i] + 1, Val(N))),
        upper = CartesianIndex(ntuple(i -> lo[i] + sz[i], Val(N))),
        size = ntuple(i -> sz[i], Val(N)),
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
    local_to_global!(gvec, dm, lvec, mode = INSERT_VALUES)

Transfer values from the local vector `lvec` to the global vector `gvec`
associated with the `dm` object.

The written vector comes first and the `dm` follows it (§8); v0.4 took the DM
last. A shim forwards the v0.4 name `dm_local_to_global!` from the old order
(§17.2); the new name only accepts the new one (§16).

# Arguments
- `gvec::AbstractPetscVec`: Global vector (destination, written)
- `dm::AbstractPetscDM`: DM object
- `lvec::AbstractPetscVec`: Local vector (source)
- `mode::InsertMode`: Insert mode, either `INSERT_VALUES` or `ADD_VALUES`

# External Links
$(doc_external("DM/DMLocalToGlobal"))
"""
function local_to_global!(
    gvec::AbstractPetscVec{PetscLib},
    dm::AbstractPetscDM{PetscLib},
    lvec::AbstractPetscVec{PetscLib},
    mode::InsertMode = INSERT_VALUES,
) where {PetscLib}
    LibPETSc.DMLocalToGlobalBegin(PetscLib, dm, lvec, mode, gvec)
    LibPETSc.DMLocalToGlobalEnd(PetscLib, dm, lvec, mode, gvec)
    return nothing
end


"""
    global_to_local!(lvec, dm, gvec, mode = INSERT_VALUES)

Transfer values from the global vector `gvec` to the local vector `lvec`
associated with the `dm` object, including ghost point values from neighboring
processes.

The written vector comes first and the `dm` follows it (§8); v0.4 had two
spellings of this call, one taking the DM last and one taking it first, and
they collapse onto this one. A shim forwards the v0.4 name
`dm_global_to_local!` from the old order (§17.2).

# Arguments
- `lvec::AbstractPetscVec`: Local vector (destination, written)
- `dm::AbstractPetscDM`: DM object
- `gvec::AbstractPetscVec`: Global vector (source)
- `mode::InsertMode`: Insert mode, either `INSERT_VALUES` or `ADD_VALUES`

# External Links
$(doc_external("DM/DMGlobalToLocal"))
"""
function global_to_local!(
    lvec::AbstractPetscVec{PetscLib},
    dm::AbstractPetscDM{PetscLib},
    gvec::AbstractPetscVec{PetscLib},
    mode::InsertMode = INSERT_VALUES,
) where {PetscLib}
    LibPETSc.DMGlobalToLocalBegin(getlib(PetscLib), dm, gvec, mode, lvec)
    LibPETSc.DMGlobalToLocalEnd(getlib(PetscLib), dm, gvec, mode, lvec)
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

The name PETSc knows this DM's flavour by, as a `Symbol` (`:da`, `:stag`,
`:plex`, …), or `nothing` when no type has been set yet
(docs/src/man/naming.md §3.1). v0.4 answered with a `String`; that is a break
with no shim (§16).

The flavour is a type parameter in v0.5 ([`DMDA`](@ref), [`DMStag`](@ref),
[`DMPlex`](@ref)), so dispatch, not this reader, is the way to branch on it.

# External Links
$(doc_external("DM/DMGetType"))
"""
type_name(dm::AbstractPetscDM{PetscLib}) where {PetscLib} =
    type_name_symbol(LibPETSc.DMGetType(PetscLib, dm))

"""
    set_type!(dm::AbstractPetscDM, type::Symbol)

Set the DM flavour, for example `:da`, `:stag` or `:plex`.

This configures a raw `PetscDM` handle; the high-level constructors set the
flavour themselves, and [`narrow`](@ref) is what puts it into the Julia type.

# External Links
$(doc_external("DM/DMSetType"))
"""
function set_type!(dm::AbstractPetscDM{PetscLib}, type::Symbol) where {PetscLib}
    LibPETSc.DMSetType(getlib(PetscLib), dm, String(type))
    return nothing
end

"""
    ndims(dm::AbstractPetscDM)

Return the topological dimension of the `dm`

# External Links
$(doc_external("DM/DMGetDimension"))
"""
Base.ndims(dm::AbstractPetscDM{PetscLib}) where {PetscLib} =
    LibPETSc.DMGetDimension(PetscLib, dm)


"""
    size(dm::DMDA{PetscLib, N})
    size(dm::DMStag{PetscLib, N})

Return the global size of the DM as an `NTuple{N,Int}`.

v0.4 returned `(M, N, P)` with the unused dimensions padded to `1`; the result
is now dimension-correct (§12), which is a break with no shim (§16).

# External Links
$(doc_external("DMDA/DMDAGetInfo"))
$(doc_external("DMSTAG/DMStagGetGlobalSizes"))
"""
function Base.size(dm::DMDA{PetscLib, N}) where {PetscLib, N}
    _, M, Ny, P, _ = LibPETSc.DMDAGetInfo(PetscLib, dm)
    gsize = (Int(M), Int(Ny), Int(P))
    return ntuple(i -> gsize[i], Val(N))
end

function Base.size(dm::DMStag{PetscLib, N}) where {PetscLib, N}
    M, Ny, P = LibPETSc.DMStagGetGlobalSizes(PetscLib, dm)
    gsize = (Int(M), Int(Ny), Int(P))
    return ntuple(i -> gsize[i], Val(N))
end

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

