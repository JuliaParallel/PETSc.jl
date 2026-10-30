# =============================================================================
# Multiple dispatch to make PetscVec behave like Julia Vector
# =============================================================================

# A PETSc vector is an `AbstractVector` of its library's scalar type. `v[i]` is one C call
# and reaches only the entries this rank owns; PETSc raises an error for any other.
# The reductions and conversions below cover what generic code would otherwise do
# entry by entry.
Base.IndexStyle(::Type{<:AbstractPetscVec}) = IndexLinear()

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
# Global length, as `Int` and as a tuple, like any `AbstractVector`
Base.size(v::AbstractPetscVec{PetscLib}) where {PetscLib} = (Int(LibPETSc.VecGetSize(PetscLib, v)),)
Base.length(v::AbstractPetscVec) = only(size(v))

"""
    similar(v::AbstractPetscVec)

A new vector with the same layout as `v` (size, parallel distribution, type), its
entries not set. `similar(v, T)` with `T = eltype(v)` gives the same. With another
element type or a size it gives a Julia `Array`, as for any `AbstractArray`: every
vector of a library holds that library's scalar type, and a different size has no
layout to copy.

As from a constructor, the new vector has a finalizer only on a one-rank
communicator; on more ranks call [`destroy!`](@ref) on it on every rank.

# External Links
$(doc_external("Vec/VecDuplicate"))
"""
Base.similar(v::AbstractPetscVec{PetscLib}) where {PetscLib} =
    finalize_if_serial!(LibPETSc.VecDuplicate(getlib(PetscLib), v))
Base.similar(v::AbstractPetscVec{PetscLib, T}, ::Type{T}) where {PetscLib, T} = similar(v)

function Base.getindex(v::AbstractPetscVec{PetscLib}, i::Integer) where {PetscLib}
    PetscInt = inttype(PetscLib)
    val = LibPETSc.VecGetValues(PetscLib, v, PetscInt(1), PetscInt[i - 1])
    return val[1]
end

# Range or index-vector indexing
function Base.getindex(v::AbstractPetscVec{PetscLib}, r::AbstractVector{<:Integer}) where {PetscLib}
    PetscInt = inttype(PetscLib)
    val = LibPETSc.VecGetValues(PetscLib, v, PetscInt(length(r)), PetscInt[i - 1 for i in r])
    return val
end

# Get all values
function Base.getindex(v::AbstractPetscVec{PetscLib}, ::Colon) where {PetscLib}
    n  = length(v)
    return getindex(v, 1:n)
end

function Base.setindex!(v::AbstractPetscVec{PetscLib}, val, i::Integer) where {PetscLib}
    PetscInt = inttype(PetscLib)
    PetscScalar = scalartype(PetscLib)
    LibPETSc.VecSetValues(PetscLib, v, PetscInt(1), PetscInt[i - 1],
                          PetscScalar[val], PETSc.INSERT_VALUES)
    return v
end

function Base.setindex!(v::AbstractPetscVec{PetscLib}, vals, r::AbstractVector{<:Integer}) where {PetscLib}
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
# Other array arguments must have the local length. Base's `materialize!` gives the
# broadcast `PetscVecStyle` from `dest`, unless an argument has more dimensions.
function Base.copyto!(dest::AbstractPetscVec, bc::Broadcast.Broadcasted{PetscVecStyle})
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

# Out-of-place broadcasting, `w = f.(args...)`, gives a new `PetscVec` with the layout
# of the first `PetscVec` among the arguments, filled by the in-place broadcast above.
# So it works on a distributed vector, and the caller owns the result: on more than
# one rank it has no finalizer (as `copy(v)`), and `destroy!` must be called on every rank.
function Base.copy(bc::Broadcast.Broadcasted{PetscVecStyle})
    flat = Broadcast.flatten(bc)
    check_scalar_result(Broadcast.combine_eltypes(flat.f, flat.args))
    w = similar(first_petscvec(flat.args))
    try
        copyto!(w, bc)
    catch
        destroy!(w)
        rethrow()
    end
    return w
end

first_petscvec(::Tuple{}) = throw(ArgumentError("a PetscVec broadcast needs a PetscVec argument"))
first_petscvec(args::Tuple) = first_petscvec(first(args), Base.tail(args))
first_petscvec(a::AbstractPetscVec, ::Tuple) = a
first_petscvec(::Any, rest::Tuple) = first_petscvec(rest)

# A PetscVec holds only its library's scalar type, so a broadcast whose entries are
# not numbers (`x .> 0` gives `Bool`) has no PetscVec to land in. Converting a mask
# to 0 and 1 silently would be a different result from Julia's, so it throws.
function check_scalar_result(::Type{T}) where {T}
    (T === Any || (T <: Number && !(T <: Bool))) && return nothing
    throw(ArgumentError(
        "the broadcast gives entries of type $T, which a PetscVec cannot hold. " *
        "Broadcast over x[:] (or with_local_array!) for the entries as a Julia array."
    ))
end

# Calls `f` with `args`, each `PetscVec` replaced by its local array: `dest_arr`
# for `dest` itself, a read-only checkout for any other.
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

aliases_dest(a::AbstractPetscVec, dest::AbstractPetscVec) = a.ptr == dest.ptr

# `dest .= f.(args...)` takes its style from the arguments alone, so `x .= 1:n` would be
# DefaultArrayStyle and run entry by entry. Run it as a PetscVecStyle broadcast.
Base.materialize!(dest::AbstractPetscVec, bc::Broadcast.Broadcasted) =
    (copyto!(dest, Broadcast.Broadcasted{PetscVecStyle}(bc.f, bc.args, bc.axes)); dest)
