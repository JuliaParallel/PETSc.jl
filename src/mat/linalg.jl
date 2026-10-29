function Base.:(==)(
    A::PetscMat{PetscLib},
    B::PetscMat{PetscLib},
) where {PetscLib}
    return LibPETSc.MatEqual(PetscLib, A, B)
end

"""
    diagonal!(d::AbstractPetscVec, A::AbstractPetscMat)

Write the diagonal of `A` into `d` and return `d`. `d` needs the row layout of
`A`, as the left vector from `MatCreateVecs` has. `LinearAlgebra.diag` is not
extended, because it returns a new vector.

# External Links
$(doc_external("Mat/MatGetDiagonal"))
"""
function diagonal!(
    d::AbstractPetscVec{PetscLib},
    A::AbstractPetscMat{PetscLib},
) where {PetscLib}
    LibPETSc.MatGetDiagonal(PetscLib, A, d)
    return d
end


"""
    mul!(y::AbstractPetscVec, M::AbstractPetscMat, x::AbstractPetscVec)

Computes
    `y` = `M`*`x`

"""
function LinearAlgebra.mul!(y::AbstractPetscVec{PetscLib}, M::AbstractPetscMat{PetscLib}, x::AbstractPetscVec{PetscLib}) where {PetscLib}
    capture_callback_errors(() -> LibPETSc.MatMult(PetscLib, M, x, y))
    return y
end

function Base.:*(
    M::AbstractPetscMat{PetscLib},
    x::AbstractPetscVec{PetscLib},
) where {PetscLib}
    _, y = LibPETSc.MatCreateVecs(getlib(PetscLib), M)
    mul!(y, M, x)
    return y
end

# `A'` and `transpose(A)` are lazy wrappers, so `A' * x` and `mul!(y, A', x)` reach
# MatMultHermitianTranspose and MatMultTranspose without forming the transpose.
LinearAlgebra.adjoint(A::AbstractPetscMat) = Adjoint(A)
LinearAlgebra.transpose(A::AbstractPetscMat) = Transpose(A)

const TransposedPetscMat{PetscLib} = Union{
    Adjoint{<:Any, <:AbstractPetscMat{PetscLib}},
    Transpose{<:Any, <:AbstractPetscMat{PetscLib}},
}

"""
    mul!(y::AbstractPetscVec, A', x::AbstractPetscVec)
    mul!(y::AbstractPetscVec, transpose(A), x::AbstractPetscVec)

Compute `y = Aᴴx` or `y = Aᵀx` without forming the transpose, and return `y`.

# External Links
$(doc_external("Mat/MatMultHermitianTranspose"))
$(doc_external("Mat/MatMultTranspose"))
"""
function LinearAlgebra.mul!(
    y::AbstractPetscVec{PetscLib},
    M::TransposedPetscMat{PetscLib},
    x::AbstractPetscVec{PetscLib},
) where {PetscLib}
    multiply = M isa Adjoint ? LibPETSc.MatMultHermitianTranspose : LibPETSc.MatMultTranspose
    capture_callback_errors(() -> multiply(PetscLib, parent(M), x, y))
    return y
end

# The result has the layout of the columns of `A`, which MatCreateVecs returns first
function Base.:*(
    M::TransposedPetscMat{PetscLib},
    x::AbstractPetscVec{PetscLib},
) where {PetscLib}
    y, _ = LibPETSc.MatCreateVecs(getlib(PetscLib), parent(M))
    mul!(y, M, x)
    return y
end

function LinearAlgebra.issymmetric(A::PetscMat{PetscLib}; tol = 0.0) where {PetscLib}
    PetscReal = real(scalartype(PetscLib))
    return LibPETSc.MatIsSymmetric(PetscLib, A, PetscReal(tol))
end
function LinearAlgebra.ishermitian(A::PetscMat{PetscLib}; tol = 0.0) where {PetscLib}
    PetscReal = real(scalartype(PetscLib))
    return LibPETSc.MatIsHermitian(PetscLib, A, PetscReal(tol))
end

"""
    copy(A::PetscMat)

A new assembled `PetscMat` with the type, layout, nonzero pattern and entries of
`A`, which must be assembled. Collective. As from a constructor, the copy has a
finalizer only on a one-rank communicator; on more ranks call
[`destroy!`](@ref) on it on every rank.

# External Links
$(doc_external("Mat/MatDuplicate"))
"""
function Base.copy(A::PetscMat{PetscLib}) where {PetscLib}
    B = LibPETSc.MatDuplicate(getlib(PetscLib), A, LibPETSc.MAT_COPY_VALUES)
    return finalize_if_serial!(B)
end

"""
    A + B
    A - B
    -A
    α * A, A * α

Matrix arithmetic on assembled `PetscMat`s, each giving a new `PetscMat` of the
type and layout of `A`. Collective. The nonzero pattern of `A ± B` is the union
of the two patterns. The result has a finalizer only on a one-rank communicator,
as [`copy`](@ref) does.

# External Links
$(doc_external("Mat/MatAXPY"))
$(doc_external("Mat/MatScale"))
"""
function Base.:+(A::PetscMat{PetscLib}, B::PetscMat{PetscLib}) where {PetscLib}
    return add_scaled(A, one(scalartype(PetscLib)), B)
end

function Base.:-(A::PetscMat{PetscLib}, B::PetscMat{PetscLib}) where {PetscLib}
    return add_scaled(A, -one(scalartype(PetscLib)), B)
end

Base.:-(A::PetscMat) = scaled_copy(A, -1)
Base.:*(α::Number, A::PetscMat) = scaled_copy(A, α)
Base.:*(A::PetscMat, α::Number) = scaled_copy(A, α)

function scaled_copy(A::PetscMat{PetscLib}, α::Number) where {PetscLib}
    C = copy(A)
    LibPETSc.MatScale(PetscLib, C, scalartype(PetscLib)(α))
    return C
end

# A copy of A plus αB. UNKNOWN_NONZERO_PATTERN lets PETSc compare the two
# patterns and merge them when they differ.
function add_scaled(A::PetscMat{PetscLib}, α, B::PetscMat{PetscLib}) where {PetscLib}
    size(A) == size(B) ||
        throw(DimensionMismatch("matrices have sizes $(size(A)) and $(size(B))"))
    C = copy(A)
    LibPETSc.MatAXPY(PetscLib, C, α, B, LibPETSc.UNKNOWN_NONZERO_PATTERN)
    return C
end

"""
    norm(A::AbstractPetscMat, p::Real = 2)
    norm(A::AbstractPetscMat, normtype::NormType)

The norm of the entries of `A` taken as one vector, as `norm` means for a Julia matrix.
Only `p = 2`, the Frobenius norm, is available: PETSc computes no other entrywise norm.
For the norms induced by a vector norm use [`opnorm`](@ref). The second form passes a
PETSc `NormType` straight to `MatNorm`, where `NORM_1` and `NORM_INFINITY` are the
induced norms.

# External Links
$(doc_external("Mat/MatNorm"))
"""
function LinearAlgebra.norm(
    M::AbstractPetscMat{PetscLib},
    normtype::NormType = NORM_FROBENIUS,
) where {PetscLib}
    return LibPETSc.MatNorm(PetscLib, M, normtype)
end

function LinearAlgebra.norm(M::AbstractPetscMat, p::Real)
    p == 2 || throw(
        ArgumentError(
            "PETSc computes only the p = 2 (Frobenius) entrywise norm of a matrix, " *
            "got p = $p; use opnorm(A, 1) or opnorm(A, Inf) for the induced norms",
        ),
    )
    return norm(M, NORM_FROBENIUS)
end

"""
    opnorm(A::AbstractPetscMat, p::Real = 2)

The operator norm of `A` induced by the vector `p`-norm: the largest column sum of
`|A|` for `p = 1`, the largest row sum for `p = Inf`. PETSc does not compute the
induced 2-norm, so `p = 2` throws an `ArgumentError`, as does any other `p`.

# External Links
$(doc_external("Mat/MatNorm"))
"""
function LinearAlgebra.opnorm(M::AbstractPetscMat, p::Real = 2)
    p == 1 && return norm(M, NORM_1)
    p == Inf && return norm(M, NORM_INFINITY)
    throw(
        ArgumentError(
            "PETSc computes the induced matrix norm only for p = 1 and p = Inf, got p = $p",
        ),
    )
end
