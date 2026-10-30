"""
    ==(A::AbstractPetscMat, B::AbstractPetscMat)
    isequal(A::AbstractPetscMat, B::AbstractPetscMat)

Whether `A` and `B` have the same size and the same entries, compared as Julia matrices
compare them: an entry stored in only one of them is compared with zero, so an
explicitly stored zero makes no difference, and `NaN != NaN` under `==`. Collective:
every rank compares the rows it owns and gets the combined answer. Matrices of the same
size with a different row layout throw an `ArgumentError`. Needs `MatGetRow`, so a
`MatShell` can't be compared.

[`equal`](@ref) is PETSc's own test, `MatEqual`, which compares what is stored.
"""
Base.:(==)(A::AbstractPetscMat{PetscLib}, B::AbstractPetscMat{PetscLib}) where {PetscLib} =
    compare_rows(==, A, B)
Base.isequal(A::AbstractPetscMat{PetscLib}, B::AbstractPetscMat{PetscLib}) where {PetscLib} =
    compare_rows(isequal, A, B)

# `cmp` over the rows each rank owns, the same answer on every rank. The layout check
# goes through the reduction so that every rank throws, not only some.
function compare_rows(cmp, A::AbstractPetscMat, B::AbstractPetscMat)
    size(A) == size(B) || return false
    rows = ownership_range(A)
    r = rows == ownership_range(B) ? Int(all(i -> rows_match(cmp, A, B, i), rows)) : -1
    r = MPI.Allreduce(r, min, comm(A))
    r < 0 && throw(ArgumentError("the matrices have the same size but different row layouts"))
    return r == 1
end

# Row `i` (1-based) of `A` and `B` compared entry by entry. PETSc allows one row of a
# matrix out at a time, so a matrix compared with itself is read once.
function rows_match(cmp, A::AbstractPetscMat{PetscLib}, B::AbstractPetscMat{PetscLib}, i) where {PetscLib}
    row = inttype(PetscLib)(i - 1)
    with_row(A, row) do ca, va
        A.ptr == B.ptr && return sparse_rows_match(cmp, ca, va, ca, va)
        with_row(B, row) do cb, vb
            sparse_rows_match(cmp, ca, va, cb, vb)
        end
    end
end

# `f(cols, vals)` with row `row` (0-based) of `A`, handed back afterwards
function with_row(f, A::AbstractPetscMat{PetscLib}, row) where {PetscLib}
    n, cols, vals = LibPETSc.MatGetRow(PetscLib, A, row)
    try
        return f(cols, vals)
    finally
        LibPETSc.MatRestoreRow(PetscLib, A, row, n, cols, vals)
    end
end

# Two sparse rows, each given by its sorted column indices and values; a column present
# in only one row is compared with zero
function sparse_rows_match(cmp, ca, va, cb, vb)
    z = zero(eltype(va))
    i, j = firstindex(ca), firstindex(cb)
    while i <= lastindex(ca) || j <= lastindex(cb)
        if j > lastindex(cb) || (i <= lastindex(ca) && ca[i] < cb[j])
            cmp(va[i], z) || return false
            i += 1
        elseif i > lastindex(ca) || cb[j] < ca[i]
            cmp(z, vb[j]) || return false
            j += 1
        else
            cmp(va[i], vb[j]) || return false
            i += 1
            j += 1
        end
    end
    return true
end

"""
    equal(A::AbstractPetscMat, B::AbstractPetscMat)
    equal(v::AbstractPetscVec, w::AbstractPetscVec)

PETSc's own equality test, `MatEqual` or `VecEqual`: whether `A` and `B` have the same
size and the same stored entries. Collective and faster than `==`, but a different
test. `MatEqual` also compares which entries are stored, so a matrix with an explicitly
stored zero is not `equal` to the same matrix without it, and `VecEqual` compares the
bits of the entries, so `NaN` equals `NaN` and `-0.0` does not equal `0.0`. Use `==`
for the comparison Julia arrays make.

# External Links
$(doc_external("Mat/MatEqual"))
$(doc_external("Vec/VecEqual"))
"""
equal(A::AbstractPetscMat{PetscLib}, B::AbstractPetscMat{PetscLib}) where {PetscLib} =
    LibPETSc.MatEqual(PetscLib, A, B) == LibPETSc.PETSC_TRUE
equal(v::AbstractPetscVec{PetscLib}, w::AbstractPetscVec{PetscLib}) where {PetscLib} =
    LibPETSc.VecEqual(PetscLib, v, w) == LibPETSc.PETSC_TRUE

"""
    isapprox(A::PetscMat, B::PetscMat; atol = 0, rtol, nans = false)

Whether `A ≈ B` as for Julia matrices: ‖A − B‖ ≤ max(atol, rtol·max(‖A‖, ‖B‖)) in the
Frobenius norm, where `rtol` defaults to `√eps` of the real type when `atol` is zero.
Collective. If ‖A − B‖ is not finite, the entries are compared one by one, as Base
does, with an entry stored in only one matrix compared with zero.
"""
function Base.isapprox(
    A::PetscMat{PetscLib},
    B::PetscMat{PetscLib};
    atol::Real = 0,
    rtol::Real = default_rtol(real(scalartype(PetscLib)), atol),
    nans::Bool = false,
) where {PetscLib}
    D = A - B
    d = try
        norm(D)
    finally
        destroy!(D)
    end
    isfinite(d) && return d <= max(atol, rtol * max(norm(A), norm(B)))
    return compare_rows((a, b) -> isapprox(a, b; atol, rtol, nans), A, B)
end

"""
    copyto!(dst::AbstractPetscMat, src::AbstractPetscMat)

Copy the entries of `src` into `dst` and return `dst`. The two must have the same size,
which throws a `DimensionMismatch` otherwise; `dst` takes the nonzero pattern of `src`.
Collective.

# External Links
$(doc_external("Mat/MatCopy"))
"""
function Base.copyto!(dst::AbstractPetscMat{PetscLib}, src::AbstractPetscMat{PetscLib}) where {PetscLib}
    size(dst) == size(src) || throw(DimensionMismatch(
        "destination has size $(size(dst)), source has size $(size(src))"))
    LibPETSc.MatCopy(PetscLib, src, dst, LibPETSc.DIFFERENT_NONZERO_PATTERN)
    return dst
end

"""
    sum(A::AbstractPetscMat)

The sum of all entries of `A`, from `MatGetRowSum`. Collective.

# External Links
$(doc_external("Mat/MatGetRowSum"))
"""
function Base.sum(A::AbstractPetscMat{PetscLib}) where {PetscLib}
    right, left = LibPETSc.MatCreateVecs(getlib(PetscLib), A)
    try
        LibPETSc.MatGetRowSum(PetscLib, A, left)
        return sum(left)
    finally
        destroy!(right)
        destroy!(left)
    end
end

"""
    mapreduce(f, op, A::AbstractPetscMat; [init])

`mapreduce` over all entries of `A`, as for a sparse Julia matrix: `f` runs on the
stored entries of the rows each rank owns, and `f(0)` stands in for every entry that
is not stored. Collective, with the same answer on every rank. As for a
`SparseMatrixCSC`, the order of the entries is not defined, so `op` should be
associative and commutative. This is what `sum(f, A)`, `prod`, `maximum`, `minimum`,
`extrema`, `reduce`, `count`, `any` and `all` call. The results must be `isbits`.
Needs `MatGetRow`, so not for a `MatShell`.
"""
function Base.mapreduce(f, op, A::AbstractPetscMat{PetscLib}; dims::Colon = :, init = Base._InitialValue()) where {PetscLib}
    T = eltype(A)
    R = reduction_type(f, op, T)
    partial, nstored = nothing, 0
    for i in ownership_range(A)
        p, k = with_row(A, inttype(PetscLib)(i - 1)) do _, vals
            isempty(vals) ? (nothing, 0) : (convert(R, mapreduce(f, op, vals)), length(vals))
        end
        isnothing(p) && continue
        partial = isnothing(partial) ? p : convert(R, op(partial, p))
        nstored += k
    end
    c = comm(A)
    parts = gather_partials(partial, R, c)
    nzeros = prod(size(A)) - MPI.Allreduce(nstored, +, c)
    r = if !isempty(parts)
        reduce_zeros(f, op, T, nzeros, foldl(op, parts))
    elseif nzeros > 0
        reduce_zeros(f, op, T, nzeros - 1, f(zero(T)))
    else
        return reduce_empty_result(f, op, T, init)
    end
    return with_init(op, init, r)
end

# `v` combined with `f(0)` once for each of `k` entries that are not stored, as
# SparseArrays does it: a sum or product in one step, any other `op` until the value
# stops changing (`max`, `|` and the like stop after one step)
function reduce_zeros(f, op, ::Type{T}, k, v) where {T}
    z = f(zero(T))
    for _ in 1:k
        w = op(v, z)
        isequal(w, v) && return w
        v = w
    end
    return v
end
reduce_zeros(f, ::Union{typeof(+), typeof(Base.add_sum)}, ::Type{T}, k, v) where {T} =
    k == 0 ? v : v + f(zero(T)) * k
reduce_zeros(f, ::Union{typeof(*), typeof(Base.mul_prod)}, ::Type{T}, k, v) where {T} =
    k == 0 ? v : v * f(zero(T))^k

# Base's `any`, `all` and `count` on an array iterate it entry by entry
Base.any(f::Function, A::AbstractPetscMat) = mapreduce(x -> f(x)::Bool, |, A)
Base.all(f::Function, A::AbstractPetscMat) = mapreduce(x -> f(x)::Bool, &, A)
Base.count(f, A::AbstractPetscMat; init = 0) = mapreduce(x -> f(x)::Bool, Base.add_sum, A; init)

"""
    collect(A::AbstractPetscMat)
    Matrix(A::AbstractPetscMat)
    copyto!(dest::Matrix, A::AbstractPetscMat)

Every entry of `A` as a dense Julia `Matrix`. Collective on a communicator of more than
one rank: each rank reads the rows it owns and receives everyone else's, so every rank
must make the call and needs memory for all `size(A, 1) * size(A, 2)` entries. On one
rank this is `A[:, :]`.
"""
function Base.collect(A::AbstractPetscMat{PetscLib}) where {PetscLib}
    m, n = size(A)
    block = A[ownership_range(A), 1:n]
    c = comm(A)
    MPI.Comm_size(c) == 1 && return block
    # Each rank's rows are contiguous in the transposed (column-major) layout
    counts = MPI.Allgather(Int[length(block)], c)
    recv = Matrix{scalartype(PetscLib)}(undef, n, m)
    MPI.Allgatherv!(permutedims(block), MPI.VBuffer(recv, counts), c)
    return permutedims(recv)
end

function Base.copyto!(dest::Matrix, src::AbstractPetscMat)
    size(dest) == size(src) || throw(DimensionMismatch(
        "destination has size $(size(dest)), source has size $(size(src))"))
    return copyto!(dest, collect(src))
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
