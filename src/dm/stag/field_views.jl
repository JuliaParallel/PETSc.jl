# ============================================================================
#   Views of a local vector by field
# ============================================================================

"""
    with_field_views!(f, dm::DMStag, vecs::AbstractPetscVec...; fields = nothing, read = true, write = true, memtype = :native)

Check out the local vectors `vecs` of `dm`, call `f` with one argument per vector,
hand the vectors back when `f` returns or throws, and return what `f` returns.

`fields` is a tuple of `location => dof` pairs, as [`IS`](@ref LibPETSc.IS) takes
them, with the 0-based component of [`stencil`](@ref). For each vector `f` gets a
tuple of views, one per field and in that order. Without `fields`, `f` gets the
whole array per vector, indexed `[I..., slot]` with `slot` from [`dof_slot`](@ref).
Either way the element index `I` is 1-based and covers the ghost elements, as in
[`ghost_corners`](@ref), so a view and [`stencil`](@ref) address the same point.

`read` and `write` work as in [`with_local_array!`](@ref): a `Bool`, or one per
vector. `write = false` checks out read-only; the views have the same types either
way. The vectors must be local vectors of `dm` ([`local_vec`](@ref)); anything else
throws an `ArgumentError` before any vector is checked out.

Vectors on a GPU, with its extension loaded, give views of device arrays, for use in
a kernel. Broadcast over `parent(view)`: an `OffsetArray` broadcasts element by
element on the host whatever its parent. `memtype = :host` gives views of host
arrays instead, as `memtype` in [`with_local_array!`](@ref), so host code runs on
device vectors; that is also what they give while no GPU extension is loaded.

```julia
flow = (face_location(dm, 1) => 0, face_location(dm, 2) => 0, element_location(dm) => 0)
c = corners(dm)
with_field_views!(dm, x, r; fields = flow, write = (false, true)) do (Vx, Vy, P), (Rx, Ry, Rp)
    for I in c.lower:c.upper      # the owned elements
        Rp[I] = Vx[I + CartesianIndex(1, 0)] - Vx[I] + Vy[I + CartesianIndex(0, 1)] - Vy[I]
    end
end
```

Vectors of two DMs are checked out by nesting, and a vector is handed back early by
closing its block. Here the stress `τ`, on a second DM `dm_τ` with the same elements,
is computed, handed back for [`local_to_local!`](@ref) to fill its ghosts, then read
by the momentum equation, all inside one checkout of the flow:

```julia
xx = (element_location(dm_τ) => 0,)
with_field_views!(dm, x, r; fields = flow, write = (false, true)) do (Vx, Vy, P), (Rx, Ry, Rp)
    with_field_views!(dm_τ, τ; fields = xx) do (Txx,)
        for I in c.lower:c.upper
            Txx[I] = 2η * (Vx[I + CartesianIndex(1, 0)] - Vx[I]) / h
        end
    end
    local_to_local!(τ, dm_τ)
    with_field_views!(dm_τ, τ; fields = xx, write = false) do (Txx,)
        for I in c.lower:c.upper
            Rx[I] = (Txx[I] - Txx[I - CartesianIndex(1, 0)] - P[I] + P[I - CartesianIndex(1, 0)]) / h
        end
    end
end
```

# External Links
$(doc_external("Vec/VecGetArrayAndMemType"))
$(doc_external("DMStag/DMStagGetLocationSlot"))
"""
function with_field_views!(
    f,
    dm::DMStag{PetscLib, N},
    vecs::Vararg{AbstractPetscVec{PetscLib}, M};
    fields = nothing,
    read::Union{Bool, NTuple{M, Bool}} = true,
    write::Union{Bool, NTuple{M, Bool}} = true,
    memtype::Symbol = :native,
) where {PetscLib, N, M}
    return with_field_views!(
        f, dm, vecs, fields, flags(read, Val(M)), flags(write, Val(M)), array_getter(memtype),
    )
end

function with_field_views!(
    f,
    dm::DMStag{PetscLib, N},
    vecs::NTuple{M, AbstractPetscVec{PetscLib}},
    fields,
    read::NTuple{M, Bool},
    write::NTuple{M, Bool},
    acquire::G,
) where {PetscLib, N, M, G}
    lib = getlib(PetscLib)
    gc = ghost_corners(dm)
    q = Int(LibPETSc.DMStagGetEntriesPerElement(lib, dm))
    dims = (q, gc.size...)
    offsets = ntuple(d -> gc.lower[d] - 1, Val(N))
    slots = field_slots(dm, fields)
    for v in vecs
        LibPETSc.VecGetLocalSize(lib, v) == prod(dims) || throw(
            ArgumentError("with_field_views! takes local vectors of the DMStag (local_vec(dm))"),
        )
    end
    acquired = map(vecs, read, write) do v, r, w
        acquire(v; read = r, write = w)
    end
    try
        return call_with_views(f, map(first, acquired), dims, offsets, slots)
    finally
        foreach(vecs, acquired, read, write) do v, (_, cpu_arr, backend), r, w
            release_local_array(cpu_arr, backend, v; read = r, write = w)
        end
    end
end

# The 1-based slot of each field, or `nothing` for the whole array
field_slots(dm, ::Nothing) = nothing
field_slots(dm::DMStag{PetscLib}, fields::Tuple) where {PetscLib} = map(fields) do (loc, dof)
    dof_slot(dm, loc, Int(dof))
end

# A function barrier: the array types come from a runtime memory type, and
# everything below, `f` included, compiles for the concrete ones
function call_with_views(f, arrays, dims, offsets, slots)
    return f(map(a -> field_views(reshape(a, dims), offsets, slots), arrays)...)
end

# The array is stored with the component first; the views put it last, or drop it
function field_views(A::AbstractArray{T, D}, offsets::NTuple{N, Int}, ::Nothing) where {T, D, N}
    whole = PermutedDimsArray(A, (ntuple(d -> d + 1, Val(N))..., 1))
    return OffsetArray(whole, offsets..., 0)
end

field_views(A::AbstractArray{T, D}, offsets::NTuple{N, Int}, slots::Tuple) where {T, D, N} =
    map(slot -> OffsetArray(view(A, slot, ntuple(_ -> Colon(), Val(N))...), offsets...), slots)
