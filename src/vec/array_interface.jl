# =============================================================================
# Multiple dispatch to make PetscVec behave like Julia Vector
# =============================================================================

# Treat PETSc vectors as 1-D array-like objects for broadcasting/indexing
Base.ndims(::Type{<:AbstractPetscVec}) = 1
Base.ndims(v::AbstractPetscVec) = ndims(typeof(v))
Base.IndexStyle(::Type{<:AbstractPetscVec}) = IndexLinear()
Base.axes(v::AbstractPetscVec) = (Base.OneTo(length(v)),)
# A Vec takes part in a broadcast as itself; the default would `collect` it
# entry by entry, since it iterates but is not an AbstractArray
Base.Broadcast.broadcastable(v::AbstractPetscVec) = v

"""
    PetscVecStyle

The broadcast style of a `PetscVec`: an expression with one in it runs on the
local arrays rather than entry by entry. It wins over an ordinary `Vector`, so
`x .+ v` uses it too.
"""
struct PetscVecStyle <: Broadcast.AbstractArrayStyle{1} end
Base.BroadcastStyle(::Type{<:AbstractPetscVec}) = PetscVecStyle()
PetscVecStyle(::Val{1}) = PetscVecStyle()
PetscVecStyle(::Val{N}) where {N} = Broadcast.DefaultArrayStyle{N}()

# The library a wrapper carries in its type parameter, as a value. Used where a
# call has to recover `petsclib` from an object rather than take it as an
# argument (§8), for instance `save_vtk(vecs, filename)`.
petsclib_of(::AbstractPetscVec{PetscLib}) where {PetscLib} = getlib(PetscLib)

# Array interface - size and length
Base.size(v::AbstractPetscVec{PetscLib}) where {PetscLib} = LibPETSc.VecGetSize(PetscLib,v)
Base.length(v::AbstractPetscVec{PetscLib}) where {PetscLib} = prod(size(v))
Base.lastindex(v::AbstractPetscVec{PetscLib}) where {PetscLib} = length(v)

# As for any vector: the global length along dimension 1, and 1 beyond it
function Base.size(v::AbstractPetscVec, d::Integer)
    d >= 1 || throw(ArgumentError("dimension must be ≥ 1, got $d"))
    return d == 1 ? Int(length(v)) : 1
end

"""
    similar(v::AbstractPetscVec)

A new vector with the same layout as `v` (size, parallel distribution, type), its
entries not set. There is no `similar(v, T)` or `similar(v, dims)`: every vector of a
library holds that library's scalar type, and a different size has no layout to copy.

# External Links
$(doc_external("Vec/VecDuplicate"))
"""
Base.similar(v::AbstractPetscVec{PetscLib}) where {PetscLib} =  LibPETSc.VecDuplicate(getlib(PetscLib), v)

function Base.getindex(v::AbstractPetscVec{PetscLib}, i::Integer) where {PetscLib}
    PetscInt = inttype(PetscLib)
    val = LibPETSc.VecGetValues(PetscLib, v, PetscInt(1), PetscInt[i - 1])
    return val[1]
end

# Range indexing
function Base.getindex(v::AbstractPetscVec{PetscLib}, r::AbstractRange) where {PetscLib}
    PetscInt = inttype(PetscLib)
    val = LibPETSc.VecGetValues(PetscLib, v, PetscInt(length(r)), PetscInt[i - 1 for i in r])
    return val
end

# Get all values
function Base.getindex(v::AbstractPetscVec{PetscLib}, ::Colon) where {PetscLib}
    n  = length(v)
    return getindex(v, 1:n)
end

"""
    isapprox(v::AbstractPetscVec, w::AbstractPetscVec; atol, rtol, nans)

Whether every entry of `v` is approximately equal to the matching entry of `w`, with
the keywords of `isapprox` for numbers applied entry by entry. Unlike `isapprox` on two
Julia arrays, this does not compare the norm of the difference. It warns once: from v0.6 it
compares norms, ‖v − w‖ ≤ max(atol, rtol·max(‖v‖, ‖w‖)), as for Julia arrays.
"""
function Base.isapprox(v::AbstractPetscVec{PetscLib}, w::AbstractPetscVec{PetscLib}; kwargs...) where {PetscLib}
    @warn "isapprox on two PetscVecs compares entry by entry; from v0.6 it compares norms as for Julia arrays, ‖v - w‖ ≤ max(atol, rtol * max(‖v‖, ‖w‖))" maxlog = 1
    return all(isapprox.(v[:], w[:]; kwargs...))
end

function Base.setindex!(v::AbstractPetscVec{PetscLib}, val, i::Integer) where {PetscLib}
    PetscInt = inttype(PetscLib)
    PetscScalar = scalartype(PetscLib)
    LibPETSc.VecSetValues(PetscLib, v, PetscInt(1), PetscInt[i - 1],
                          PetscScalar[val], PETSc.INSERT_VALUES)
    return v
end

function Base.setindex!(v::AbstractPetscVec{PetscLib}, vals, r::AbstractRange) where {PetscLib}
    PetscInt = inttype(PetscLib)
    LibPETSc.VecSetValues(PetscLib, v, PetscInt(length(r)), PetscInt[i - 1 for i in r],
                          vals, PETSc.INSERT_VALUES)
    return v
end

function Base.fill!(v::AbstractPetscVec{PetscLib}, val) where {PetscLib}
    LibPETSc.VecSet(PetscLib, v, PetscLib.PetscScalar(val))
    return v
end

"""
    copyto!(dst::AbstractPetscVec, src::AbstractPetscVec)

Copy the entries of `src` into `dst` and return `dst`. The two vectors must
have the same global length, which throws a `DimensionMismatch` otherwise, and
the same parallel layout, which PETSc checks.

# External Links
$(doc_external("Vec/VecCopy"))
"""
function Base.copyto!(
    dst::AbstractPetscVec{PetscLib},
    src::AbstractPetscVec{PetscLib},
) where {PetscLib}
    length(dst) == length(src) || throw(
        DimensionMismatch("destination has length $(length(dst)), source has length $(length(src))"),
    )
    LibPETSc.VecCopy(PetscLib, src, dst)
    return dst
end

# In-place broadcasting, `dest .= f.(args...)`, runs on the local arrays:
# `dest` and every `PetscVec` among the arguments are checked out, the broadcast runs
# over the arrays, and all are handed back. So it is elementwise over the entries
# this rank owns, which on one rank is every entry.
# Other array arguments must have the local length.
function Base.copyto!(dest::AbstractPetscVec, bc::Base.Broadcast.Broadcasted)
    flat = Base.Broadcast.flatten(bc)
    arr, raw, backend = acquire_local_array(dest; read = true, write = true)
    try
        with_broadcast_arrays(dest, arr, flat.args) do args
            copyto!(arr, Base.Broadcast.instantiate(Base.Broadcast.Broadcasted(flat.f, args)))
        end
    finally
        release_local_array(raw, backend, dest; read = true, write = true)
    end
    return dest
end

# Out-of-place broadcasting, `w = f.(args...)`, runs on the local arrays too and
# gives a `Vector`, so `2 .* x` is what it was before, without the one C call per
# entry. A `Vector` of the whole vector only means something where one rank holds
# all of it, so a distributed vector throws instead.
function Base.copy(bc::Broadcast.Broadcasted{PetscVecStyle})
    flat = Broadcast.flatten(bc)
    foreach(check_undistributed, flat.args)
    return with_broadcast_arrays(nothing, nothing, flat.args) do args
        copy(Broadcast.instantiate(Broadcast.Broadcasted(flat.f, args)))
    end
end

check_undistributed(_) = nothing
function check_undistributed(v::AbstractPetscVec{PetscLib}) where {PetscLib}
    lib = getlib(PetscLib)
    LibPETSc.VecGetLocalSize(lib, v) == LibPETSc.VecGetSize(lib, v) || throw(ArgumentError(
        "a broadcast of a distributed Vec cannot give a Vector of the whole vector. " *
        "Broadcast into a Vec of the same layout, or use with_local_array! for the " *
        "entries this rank owns."
    ))
    return nothing
end

# Calls `f` with `args`, each `PetscVec` replaced by its local array: `dest_arr`
# for `dest` itself, a read-only checkout for any other. `dest` is `nothing` when
# the broadcast is out of place and nothing can alias.
with_broadcast_arrays(f, dest, dest_arr, ::Tuple{}) = f(())
function with_broadcast_arrays(f, dest, dest_arr, args::Tuple)
    a = first(args)
    rest = Base.tail(args)
    a isa AbstractPetscVec || return with_broadcast_arrays(r -> f((a, r...)), dest, dest_arr, rest)
    aliases_dest(a, dest) && return with_broadcast_arrays(r -> f((dest_arr, r...)), dest, dest_arr, rest)
    arr, raw, backend = acquire_local_array(a; read = true, write = false)
    try
        return with_broadcast_arrays(r -> f((arr, r...)), dest, dest_arr, rest)
    finally
        release_local_array(raw, backend, a; read = true, write = false)
    end
end

aliases_dest(::AbstractPetscVec, ::Nothing) = false
aliases_dest(a::AbstractPetscVec, dest::AbstractPetscVec) = a.ptr == dest.ptr

Base.materialize!(dest::AbstractPetscVec, bc::Base.Broadcast.Broadcasted) = (Base.copyto!(dest, bc); dest)

# Iterator interface
Base.iterate(v::AbstractPetscVec{PetscLib})  where {PetscLib} = iterate(v, 1)
function Base.iterate(v::AbstractPetscVec{PetscLib}, state)  where {PetscLib}
    if state > length(v)
        return nothing
    end
    return (v[state], state + 1)
end
