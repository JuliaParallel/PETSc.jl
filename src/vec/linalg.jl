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
    w = finalize_if_serial!(similar(x))
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
    w = finalize_if_serial!(similar(x))
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
