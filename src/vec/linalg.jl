# =============================================================================
# Linear algebra on Vecs
# =============================================================================

"""
    norm(v::AbstractPetscVec, p::Real = 2)
    norm(v::AbstractPetscVec, normtype::NormType)

The `p`-norm of `v`, for `p` = 1, 2 or `Inf`; other values throw an `ArgumentError`.
The second form passes a PETSc `NormType` straight to `VecNorm`.

# External Links
$(doc_external("Vec/VecNorm"))
"""
function LinearAlgebra.norm(
    v::AbstractPetscVec{PetscLib},
    normtype::LibPETSc.NormType = LibPETSc.NORM_2,
) where {PetscLib}
    return LibPETSc.VecNorm(PetscLib, v, normtype)
end

function LinearAlgebra.norm(v::AbstractPetscVec, p::Real)
    p == 1 && return norm(v, LibPETSc.NORM_1)
    p == 2 && return norm(v, LibPETSc.NORM_2)
    p == Inf && return norm(v, LibPETSc.NORM_INFINITY)
    throw(ArgumentError("PETSc computes the vector norm only for p = 1, 2 and Inf, got p = $p"))
end

# A new Vec or Mat gets a finalizer only on one rank, as from its constructor:
# its Destroy is collective, and a GC finalizer runs at an arbitrary point
function finalize_if_serial!(
    obj::Union{AbstractPetscVec{PetscLib}, LibPETSc.AbstractPetscMat{PetscLib}},
) where {PetscLib}
    if MPI.Comm_size(LibPETSc.PetscObjectGetComm(getlib(PetscLib), obj)) == 1
        finalizer(destroy!, obj)
    end
    return obj
end

function check_same_length(x::AbstractPetscVec, y::AbstractPetscVec)
    length(x) == length(y) ||
        throw(DimensionMismatch("vectors have lengths $(length(x)) and $(length(y))"))
    return nothing
end

"""
    copy(v::AbstractPetscVec)

A new `PetscVec` with the layout and entries of `v`. Collective. As from a
constructor, the copy has a finalizer only on a one-rank communicator; on more
ranks call [`destroy!`](@ref) on it on every rank.

# External Links
$(doc_external("Vec/VecDuplicate"))
$(doc_external("Vec/VecCopy"))
"""
function Base.copy(v::AbstractPetscVec{PetscLib}) where {PetscLib}
    w = LibPETSc.VecDuplicate(getlib(PetscLib), v)
    LibPETSc.VecCopy(PetscLib, v, w)
    return finalize_if_serial!(w)
end

"""
    dot(x::AbstractPetscVec, y::AbstractPetscVec)

`xᴴy`, conjugating `x` as `dot` does for Julia arrays. Collective: every rank
gets the global result.

# External Links
$(doc_external("Vec/VecDot"))
"""
function LinearAlgebra.dot(
    x::AbstractPetscVec{PetscLib},
    y::AbstractPetscVec{PetscLib},
) where {PetscLib}
    check_same_length(x, y)
    # VecDot(x, y) is yᴴx
    return LibPETSc.VecDot(PetscLib, y, x)
end

"""
    axpy!(α, x::AbstractPetscVec, y::AbstractPetscVec)

`y = αx + y`, returning `y`. Collective.

# External Links
$(doc_external("Vec/VecAXPY"))
"""
function LinearAlgebra.axpy!(
    α::Number,
    x::AbstractPetscVec{PetscLib},
    y::AbstractPetscVec{PetscLib},
) where {PetscLib}
    check_same_length(x, y)
    LibPETSc.VecAXPY(PetscLib, y, scalartype(PetscLib)(α), x)
    return y
end

"""
    axpby!(α, x::AbstractPetscVec, β, y::AbstractPetscVec)

`y = αx + βy`, returning `y`. Collective.

# External Links
$(doc_external("Vec/VecAXPBY"))
"""
function LinearAlgebra.axpby!(
    α::Number,
    x::AbstractPetscVec{PetscLib},
    β::Number,
    y::AbstractPetscVec{PetscLib},
) where {PetscLib}
    check_same_length(x, y)
    PetscScalar = scalartype(PetscLib)
    LibPETSc.VecAXPBY(PetscLib, y, PetscScalar(α), PetscScalar(β), x)
    return y
end

"""
    rmul!(v::AbstractPetscVec, α)
    lmul!(α, v::AbstractPetscVec)

Scale `v` by `α` in place and return it. Collective.

# External Links
$(doc_external("Vec/VecScale"))
"""
function LinearAlgebra.rmul!(v::AbstractPetscVec{PetscLib}, α::Number) where {PetscLib}
    LibPETSc.VecScale(PetscLib, v, scalartype(PetscLib)(α))
    return v
end
LinearAlgebra.lmul!(α::Number, v::AbstractPetscVec) = rmul!(v, α)

"""
    x + y
    x - y
    -x
    α * x, x * α
    x / α

Vector arithmetic on `PetscVec`s, each giving a new `PetscVec` with the layout of
`x`. Collective. The result has a finalizer only on a one-rank communicator, as
[`copy`](@ref) does. For an update in place use [`axpy!`](@ref), [`axpby!`](@ref)
and [`rmul!`](@ref), or broadcast into an existing vector, `y .= α .* x .+ y`.

# External Links
$(doc_external("Vec/VecWAXPY"))
$(doc_external("Vec/VecScale"))
"""
function Base.:+(x::AbstractPetscVec{PetscLib}, y::AbstractPetscVec{PetscLib}) where {PetscLib}
    return waxpy(one(scalartype(PetscLib)), x, y)
end

function Base.:-(x::AbstractPetscVec{PetscLib}, y::AbstractPetscVec{PetscLib}) where {PetscLib}
    return waxpy(-one(scalartype(PetscLib)), y, x)
end

Base.:-(x::AbstractPetscVec) = rmul!(copy(x), -1)
Base.:*(α::Number, x::AbstractPetscVec) = rmul!(copy(x), α)
Base.:*(x::AbstractPetscVec, α::Number) = rmul!(copy(x), α)
# Through the broadcast rather than VecScale with 1/α, so each entry is divided
# exactly as `x ./ α` would
function Base.:/(x::AbstractPetscVec{PetscLib}, α::Number) where {PetscLib}
    w = similar(x)
    w .= x ./ scalartype(PetscLib)(α)
    return w
end

# A new vector αx + y with the layout of x
function waxpy(
    α,
    x::AbstractPetscVec{PetscLib},
    y::AbstractPetscVec{PetscLib},
) where {PetscLib}
    check_same_length(x, y)
    w = similar(x)
    LibPETSc.VecWAXPY(PetscLib, w, α, x, y)
    return w
end

"""
    isapprox(v::AbstractPetscVec, w::AbstractPetscVec; atol = 0, rtol, nans = false, norm = norm)

Whether `v ≈ w` as for Julia arrays: ‖v − w‖ ≤ max(atol, rtol·max(‖v‖, ‖w‖)), where
`rtol` defaults to `√eps` of the real type when `atol` is zero. Collective: every rank
gets the same answer. If ‖v − w‖ is not finite, because of an `Inf` or `NaN` entry, the
entries are compared one by one, as Base does, and `nans = true` counts two `NaN`s as
equal.
"""
function Base.isapprox(
    v::AbstractPetscVec{PetscLib},
    w::AbstractPetscVec{PetscLib};
    atol::Real = 0,
    rtol::Real = default_rtol(real(scalartype(PetscLib)), atol),
    nans::Bool = false,
    norm::Function = LinearAlgebra.norm,
) where {PetscLib}
    diff = v - w
    d = try
        norm(diff)
    finally
        destroy!(diff)
    end
    isfinite(d) && return d <= max(atol, rtol * max(norm(v), norm(w)))
    return entrywise_isapprox(v, w; atol, rtol, nans)
end

# Base's default for `isapprox`: relative tolerance √eps, unless an absolute one is given
default_rtol(::Type{R}, atol) where {R} = atol > 0 ? zero(R) : sqrt(eps(R))

# Every owned entry of `v` ≈ the matching entry of `w`, reduced over the ranks
function entrywise_isapprox(v::AbstractPetscVec{PetscLib}, w::AbstractPetscVec{PetscLib}; kwargs...) where {PetscLib}
    ok = with_local_array!((v, w); write = false) do a, b
        all(isapprox.(a, b; kwargs...))
    end
    comm = LibPETSc.PetscObjectGetComm(getlib(PetscLib), v)
    return MPI.Allreduce(Int(ok), min, comm) == 1
end

# The communicator of `v`
comm_of(v::AbstractPetscVec{PetscLib}) where {PetscLib} =
    LibPETSc.PetscObjectGetComm(getlib(PetscLib), v)

# `f` applied to the local array of `v`, combined over the ranks with `op`
reduce_owned(f, op, v::AbstractPetscVec) =
    MPI.Allreduce(with_local_array!(f, v; write = false), op, comm_of(v))

# `cmp` over the owned entries of `v` and `w`, with the same answer on every rank.
# Vectors of one length but different layouts can't be compared entry by entry; the
# layout check goes through the reduction so that every rank throws, not only some.
function compare_owned(cmp, v::AbstractPetscVec, w::AbstractPetscVec)
    length(v) == length(w) || return false
    r = with_local_array!((v, w); write = false) do a, b
        length(a) == length(b) ? Int(cmp(a, b)) : -1
    end
    r = MPI.Allreduce(r, min, comm_of(v))
    r < 0 && throw(ArgumentError("the vectors have the same length but different parallel layouts"))
    return r == 1
end

"""
    ==(v::AbstractPetscVec, w::AbstractPetscVec)
    isequal(v::AbstractPetscVec, w::AbstractPetscVec)

Whether `v` and `w` have the same length and the same entries, compared as Julia
arrays compare them (so `NaN != NaN` and `0.0 == -0.0` under `==`). Collective: every
rank gets the same answer. Vectors of the same length with a different parallel layout
throw an `ArgumentError`.

`hash` of a PETSc vector throws, since its entries can change: key a `Dict` by the
handle with an `IdDict` instead.
"""
Base.:(==)(v::AbstractPetscVec{PetscLib}, w::AbstractPetscVec{PetscLib}) where {PetscLib} =
    compare_owned(==, v, w)
Base.isequal(v::AbstractPetscVec{PetscLib}, w::AbstractPetscVec{PetscLib}) where {PetscLib} =
    compare_owned(isequal, v, w)

Base.hash(::AbstractPetscVec, ::UInt) = throw(ArgumentError(
    "a PETSc vector has no hash: `==` compares its entries, which can change. " *
    "Key a Dict by the handle with an IdDict instead."
))

"""
    sum(v::AbstractPetscVec)
    maximum(v::AbstractPetscVec)
    minimum(v::AbstractPetscVec)

The sum, largest and smallest entry of `v`, from `VecSum`, `VecMax` and `VecMin`.
Collective. `maximum` and `minimum` need a real scalar type and a nonempty vector.

# External Links
$(doc_external("Vec/VecSum"))
$(doc_external("Vec/VecMax"))
$(doc_external("Vec/VecMin"))
"""
Base.sum(v::AbstractPetscVec{PetscLib}) where {PetscLib} = LibPETSc.VecSum(PetscLib, v)

function Base.maximum(v::AbstractPetscVec{PetscLib, <:Real}) where {PetscLib}
    check_nonempty(v)
    return last(LibPETSc.VecMax(PetscLib, v))
end

function Base.minimum(v::AbstractPetscVec{PetscLib, <:Real}) where {PetscLib}
    check_nonempty(v)
    return last(LibPETSc.VecMin(PetscLib, v))
end

check_nonempty(v) = isempty(v) && throw(ArgumentError("reducing over an empty collection is not allowed"))

"""
    mapreduce(f, op, v::AbstractPetscVec; [init])

`mapreduce` over all entries of `v`, collective: each rank reduces the entries it owns
and the partial results are combined in rank order, so every rank gets the same
answer, also for an `op` that is not commutative. This is what `sum(f, v)`,
`prod(v)`, `extrema(v)` and `reduce(op, v)` call. The results of `f` must be `isbits`.
`init` is applied once, to the combined result.
"""
function Base.mapreduce(f, op, v::AbstractPetscVec; init = Base._InitialValue())
    R = reduction_type(f, op, eltype(v))
    partial = with_local_array!(v; write = false) do a
        isempty(a) ? nothing : mapreduce(f, op, a)
    end
    parts = gather_partials(partial, R, comm_of(v))
    isempty(parts) && return reduce_empty_result(f, op, eltype(v), init)
    return with_init(op, init, foldl(op, parts))
end

# The type `op` gives when it combines results of `f` on entries of type `T`. The
# partial results are gathered over the ranks, so it must be a concrete isbits type.
function reduction_type(f, op, ::Type{T}) where {T}
    F = Base.promote_op(f, T)
    R = Base.promote_op(op, F, F)
    isconcretetype(R) && isbitstype(R) || throw(ArgumentError(
        "mapreduce over a PETSc object needs results of an isbits type, got $R"
    ))
    return R
end

# The partial results of the ranks that have one (`nothing` for none), in rank order,
# on every rank
function gather_partials(partial, ::Type{R}, comm) where {R}
    has = !isnothing(partial)
    flags = MPI.Allgather(Int[has], comm)
    vals = MPI.Allgather(R[has ? convert(R, partial) : only(Vector{R}(undef, 1))], comm)
    return vals[flags .== 1]
end

# What Base gives for a reduction over no entries: `init` if given, else the neutral
# element of `op`, or an error for an `op` without one, such as `max`
reduce_empty_result(f, op, ::Type{T}, init) where {T} =
    init isa Base._InitialValue ? Base.mapreduce_empty(f, op, T) : init

with_init(op, init, r) = init isa Base._InitialValue ? r : op(init, r)

"""
    any(f, v::AbstractPetscVec)
    all(f, v::AbstractPetscVec)
    count(f, v::AbstractPetscVec; init = 0)

Whether `f` holds for any or all entries of `v`, or for how many. Collective: `f` runs
on the entries each rank owns, and every rank gets the combined answer. On a device
vector `f` runs on the device array.
"""
Base.any(f::Function, v::AbstractPetscVec) = reduce_owned(a -> Int(any(f, a)), max, v) == 1
Base.all(f::Function, v::AbstractPetscVec) = reduce_owned(a -> Int(all(f, a)), min, v) == 1
Base.count(f, v::AbstractPetscVec; init = 0) = init + reduce_owned(a -> count(f, a), +, v)

"""
    collect(v::AbstractPetscVec)
    Vector(v::AbstractPetscVec)
    copyto!(dest::Vector, v::AbstractPetscVec)

Every entry of `v` as a Julia `Vector`. Collective on a communicator of more than one
rank: each rank receives the whole vector (`VecScatterCreateToAll`), so every rank must
make the call. On one rank this is `v[:]`. Needs memory for all `length(v)` entries on
each rank.

# External Links
$(doc_external("Vec/VecScatterCreateToAll"))
"""
function Base.collect(v::AbstractPetscVec{PetscLib}) where {PetscLib}
    MPI.Comm_size(comm_of(v)) == 1 && return v[:]
    petsclib = getlib(PetscLib)
    ctx, w = LibPETSc.VecScatterCreateToAll(petsclib, v)
    try
        LibPETSc.VecScatterBegin(petsclib, ctx, v, w, LibPETSc.INSERT_VALUES, LibPETSc.SCATTER_FORWARD)
        LibPETSc.VecScatterEnd(petsclib, ctx, v, w, LibPETSc.INSERT_VALUES, LibPETSc.SCATTER_FORWARD)
        return w[:]
    finally
        destroy!(w)
        LibPETSc.VecScatterDestroy(petsclib, ctx)
    end
end

function Base.copyto!(dest::Vector, src::AbstractPetscVec)
    length(dest) == length(src) || throw(DimensionMismatch(
        "destination has length $(length(dest)), source has length $(length(src))"))
    return copyto!(dest, collect(src))
end
