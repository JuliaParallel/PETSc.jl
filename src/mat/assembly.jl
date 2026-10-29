"""
    assemble!(A::PetscMat) 

Assembles a PETSc matrix after setting values.
"""
function assemble!(A::AbstractPetscMat{PetscLib}) where {PetscLib}
    LibPETSc.MatAssemblyBegin(PetscLib, A, PETSc.MAT_FINAL_ASSEMBLY)
    LibPETSc.MatAssemblyEnd(PetscLib, A, PETSc.MAT_FINAL_ASSEMBLY)
    return A
end

"""
    isassembled(A::AbstractPetscMat)

Whether `A` has been assembled with [`assemble!`](@ref) since its last change.

# External Links
$(doc_external("Mat/MatAssembled"))
"""
isassembled(A::AbstractPetscMat{PetscLib}) where {PetscLib} =
    Bool(LibPETSc.MatAssembled(PetscLib, A))

"""
    fill!(A::AbstractPetscMat, 0)

Set every stored entry of `A` to zero and return `A`. The nonzero pattern
stays, so the matrix can be refilled without a new allocation. Only zero is
accepted: PETSc has no operation that sets every entry to another value, and
any other `x` throws an `ArgumentError`.

# External Links
$(doc_external("Mat/MatZeroEntries"))
"""
function Base.fill!(A::AbstractPetscMat{PetscLib}, x) where {PetscLib}
    iszero(x) || throw(ArgumentError("fill! on a PETSc matrix accepts only zero, got $x"))
    LibPETSc.MatZeroEntries(PetscLib, A)
    return A
end

"""
    zero_rows!(A::AbstractPetscMat, rows_0b; diag = 1, x = nothing, b = nothing)

Zero the rows `rows_0b` of the assembled matrix `A`, put `diag` on their
diagonal entries and return `A`. With `x` and `b` given, also set
`b[i] = diag * x[i]` for each zeroed row `i`, so that a solve keeps the values
of `x` there: the usual way to impose a Dirichlet condition.

`rows_0b` are 0-based global row numbers, as PETSc takes them (docs/src/man/naming.md
§12.1). Collective: every process calls it, each with its own rows, possibly none.
[`zero_rows_local!`](@ref) takes local numbers instead.

# External Links
$(doc_external("Mat/MatZeroRows"))
"""
function zero_rows!(
    A::AbstractPetscMat{PetscLib},
    rows_0b::AbstractVector{<:Integer};
    diag = 1,
    x = nothing,
    b = nothing,
) where {PetscLib}
    PetscInt = PetscLib.PetscInt
    xv, bv = _zero_rows_vecs(PetscLib, x, b)
    LibPETSc.MatZeroRows(PetscLib, A, PetscInt(length(rows_0b)), Vector{PetscInt}(rows_0b),
        PetscLib.PetscScalar(diag), xv, bv)
    return A
end

"""
    zero_rows_local!(A::AbstractPetscMat, rows_0b; diag = 1, x = nothing, b = nothing)

[`zero_rows!`](@ref) with 0-based local row numbers, translated through the
local-to-global mapping of `A`. A matrix from `DMCreateMatrix` has that mapping;
one built without it throws a `PetscError`.

# External Links
$(doc_external("Mat/MatZeroRowsLocal"))
"""
function zero_rows_local!(
    A::AbstractPetscMat{PetscLib},
    rows_0b::AbstractVector{<:Integer};
    diag = 1,
    x = nothing,
    b = nothing,
) where {PetscLib}
    PetscInt = PetscLib.PetscInt
    xv, bv = _zero_rows_vecs(PetscLib, x, b)
    LibPETSc.MatZeroRowsLocal(PetscLib, A, PetscInt(length(rows_0b)), Vector{PetscInt}(rows_0b),
        PetscLib.PetscScalar(diag), xv, bv)
    return A
end

# `x` and `b` of MatZeroRows go together; a missing one is passed as NULL
function _zero_rows_vecs(::Type{PetscLib}, x, b) where {PetscLib}
    isnothing(x) == isnothing(b) ||
        throw(ArgumentError("zero_rows! takes both x and b, or neither"))
    null_vec = LibPETSc.PetscVec{PetscLib}(C_NULL, 0; own = false)
    return something(x, null_vec), something(b, null_vec)
end

"""
    copyto!(M::AbstractPetscMat, S::SparseMatrixCSC)

Write the nonzeros of `S` into `M` and return `M`. Call [`assemble!`](@ref) afterwards.

`S` is indexed from the first row this rank owns: entry `S[i, j]` goes to global row
`r + i` and global column `r + j`, where `r` is the number of rows before this rank's.
In serial `r = 0`, so `S` is copied as it stands. In parallel `S` is the rank's diagonal
block, and couplings to other ranks' unknowns are not written. Throws a
`DimensionMismatch` when `S` reaches past this rank's rows or past the last column.
"""
function Base.copyto!(
    M::AbstractPetscMat{PetscLib},
    S::SparseMatrixCSC,
) where {PetscLib}
    row_rng = LibPETSc.MatGetOwnershipRange(PetscLib,M)
    PetscInt = PetscLib.PetscInt
    PetscScalar = PetscLib.PetscScalar
    row_start = row_rng[1]
    m, n = size(S)
    m <= row_rng[2] - row_start && row_start + n <= size(M, 2) || throw(
        DimensionMismatch(
            "a $(m)x$(n) block starting at row $(row_start + 1) does not fit the " *
            "$(row_rng[2] - row_start) rows this rank owns of a matrix with $(size(M, 2)) columns",
        ),
    )
    for j in 1:n
        for ii in S.colptr[j]:(S.colptr[j + 1] - 1)
            i = S.rowval[ii]
            M[PetscInt(i + row_start), PetscInt(j + row_start)] = PetscScalar(S.nzval[ii])
        end
    end
    return M
end

"""
    set_values!(
        M::AbstractPetscMat,
        rows_0b::AbstractVector{MatStencil},
        cols_0b::AbstractVector{MatStencil},
        rowvals::AbstractArray;
        mode::InsertMode = INSERT_VALUES,
        num_rows = length(rows_0b),
        num_cols = length(cols_0b),
    )

Set values of the matrix `M` with base-0 row and column indices `rows_0b` and
`cols_0b`, inserting the values `rowvals`, read row by row in linear order.
A `Vector{MatStencil}` and a `Vector` (or `Array`) of the library's scalar type 
go to PETSc as they are; anything else is copied first.

The `_0b` suffix says what the base is: a bulk index array handed to C keeps
PETSc's base rather than being rebuilt on a hot path (docs/src/man/naming.md
§12.1). `A[i, j] = v` is the 1-based route.

`mode` is `INSERT_VALUES` or `ADD_VALUES`. If the keyword arguments `num_rows` or
`num_cols` is specified then only the first `num_rows * num_cols` values of
`rowvals` will be used.

# External Links
$(doc_external("Mat/MatSetValuesStencil"))
"""
function set_values!(
    M::AbstractPetscMat{PetscLib},
    rows_0b::AbstractVector{MatStencil},
    cols_0b::AbstractVector{MatStencil},
    rowvals::AbstractArray;
    mode::InsertMode = INSERT_VALUES,
    num_rows = length(rows_0b),
    num_cols = length(cols_0b),
) where {PetscLib}
    PetscInt = PetscLib.PetscInt
    check_block(num_rows, num_cols, length(rows_0b), length(cols_0b), length(rowvals))
    LibPETSc.MatSetValuesStencil(
        PetscLib,
        M,
        PetscInt(num_rows),
        c_vector(MatStencil, rows_0b),
        PetscInt(num_cols),
        c_vector(MatStencil, cols_0b),
        c_vector(PetscLib.PetscScalar, rowvals),
        mode,
    )
    return M
end

# The array PETSc reads, of element type `T`: a `Vector{T}` as it is, any other
# `Array{T}` reshaped without a copy, anything else copied in linear order, so
# views, ranges and offset arrays all work.
c_vector(::Type{T}, v::Vector{T}) where {T} = v
c_vector(::Type{T}, v::Array{T}) where {T} = vec(v)
c_vector(::Type{T}, v::AbstractArray) where {T} = copyto!(Vector{T}(undef, length(v)), v)

# A num_rows x num_cols block needs that many indices and values
function check_block(num_rows, num_cols, nrows, ncols, nvals)
    num_rows <= nrows && num_cols <= ncols || throw(
        DimensionMismatch(
            "a $(num_rows)x$(num_cols) block needs $num_rows row and $num_cols column " *
            "indices, got $nrows and $ncols",
        ),
    )
    num_rows * num_cols <= nvals || throw(
        DimensionMismatch(
            "a $(num_rows)x$(num_cols) block needs $(num_rows * num_cols) values, " *
            "got $nvals",
        ),
    )
    return nothing
end

"""
    set_values!(
        M::AbstractPetscMat,
        rows_0b::AbstractVector{<:Integer},
        cols_0b::AbstractVector{<:Integer},
        rowvals::AbstractArray;
        mode::InsertMode = INSERT_VALUES,
        num_rows = length(rows_0b),
        num_cols = length(cols_0b),
    )

Set values of the matrix `M` with base-0 row and column indices `rows_0b` and
`cols_0b`, inserting the values `rowvals`, read row by row in linear order.
A `Vector` of the library's integer type for the indices, and a `Vector` (or `Array`)
of its scalar type for the values, go to PETSc as they are; anything else is
copied first.

The `_0b` suffix says what the base is: a bulk index array handed to C keeps
PETSc's base rather than being rebuilt on a hot path (docs/src/man/naming.md
§12.1). `A[i, j] = v` is the 1-based route.

`mode` is `INSERT_VALUES` or `ADD_VALUES`. If the keyword arguments `num_rows` or
`num_cols` is specified then only the first `num_rows * num_cols` values of
`rowvals` will be used.

# External Links
$(doc_external("Mat/MatSetValues"))
"""
function set_values!(
    M::AbstractPetscMat{PetscLib},
    rows_0b::AbstractVector{<:Integer},
    cols_0b::AbstractVector{<:Integer},
    rowvals::AbstractArray;
    mode::InsertMode = INSERT_VALUES,
    num_rows = length(rows_0b),
    num_cols = length(cols_0b),
) where {PetscLib}
    PetscInt = PetscLib.PetscInt
    check_block(num_rows, num_cols, length(rows_0b), length(cols_0b), length(rowvals))
    LibPETSc.MatSetValues(
        PetscLib,
        M,
        PetscInt(num_rows),
        c_vector(PetscInt, rows_0b),
        PetscInt(num_cols),
        c_vector(PetscInt, cols_0b),
        c_vector(PetscLib.PetscScalar, rowvals),
        mode,
    )
    return M
end
