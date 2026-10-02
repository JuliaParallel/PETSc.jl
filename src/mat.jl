
import .LibPETSc: AbstractPetscMat, PetscMat, CMat, MatStencil, InsertMode

# Custom display for REPL
function Base.show(io::IO, v::AbstractPetscMat{PetscLib}) where {PetscLib}
    if v.ptr == C_NULL
        print(io, "PETSc Mat (null pointer)")
        return
    end
    
    mat_type = type_name(v)
    if isnothing(mat_type)
        print(io, "PETSc Mat (type not set)")
    else
        print(io, "PETSc $(mat_type) Mat of size $(size(v))")
    end
    return nothing
end

"""
    MatPtr(petsclib, ptr::CMat, own::Bool)

Container type for a PETSc Mat that is just a raw pointer.

If `own` is `true` a finalizer is set on the matrix, but only on a serial
communicator, since `MatDestroy` is collective and a GC finalizer runs at an
arbitrary point. If `own` is `false` the handle belongs to PETSc and `destroy!`
is a no-op, leaving the wrapper usable.
"""
mutable struct MatPtr{PetscLib} <:
               AbstractPetscMat{PetscLib}
    ptr::CMat
    age::Int
    own::Bool
end
function MatPtr(
    petsclib::PetscLib,
    ptr::CMat,
    own,
) where {PetscLib <: PetscLibType}
    m = MatPtr{PetscLib}(ptr, petsclib.age, own)
    # Short-circuits on a borrowed handle, which is the hot path: callbacks wrap
    # PETSc-owned matrices on every invocation and never need the communicator.
    if own && MPI.Comm_size(LibPETSc.PetscObjectGetComm(getlib(PetscLib), m)) == 1
        finalizer(destroy!, m)
    end
    return m
end
MatPtr(::Type{PetscLib}, x...) where {PetscLib <: PetscLibType} =
    MatPtr(getlib(PetscLib), x...)

Base.size(m::AbstractPetscMat{PetscLib}) where {PetscLib} = LibPETSc.MatGetSize(PetscLib,m)
Base.length(m::AbstractPetscMat{PetscLib}) where {PetscLib} = prod(size(m))
Base.ndims(::Type{<:AbstractPetscMat}) = 2
Base.ndims(m::AbstractPetscMat) = ndims(typeof(m))

# As for any matrix: the global size along dimensions 1 and 2, and 1 beyond them
function Base.size(m::AbstractPetscMat, d::Integer)
    d >= 1 || throw(ArgumentError("dimension must be ≥ 1, got $d"))
    return d <= 2 ? Int(size(m)[d]) : 1
end
Base.axes(m::AbstractPetscMat) = map(Base.OneTo, size(m))
Base.axes(m::AbstractPetscMat, d::Integer) = Base.OneTo(size(m, d))
"""
    type_name(A::AbstractPetscMat)

The name PETSc knows this matrix's implementation by, as a `Symbol`
(`:seqaij`, `:mpiaij`, …), or `nothing` when the matrix has no type yet
(docs/src/man/naming.md §3.1). v0.4 answered with a `String`, using the
sentinel `"(not set)"`; that is a break with no shim (§16).

# External Links
$(doc_external("Mat/MatGetType"))
"""
type_name(m::AbstractPetscMat{PetscLib}) where {PetscLib} =
    type_name_symbol(LibPETSc.MatGetType(PetscLib, m))

"""
    set_type!(A::AbstractPetscMat, type::Symbol)

Set the matrix implementation, for example `:seqaij` or `:dense`.

# External Links
$(doc_external("Mat/MatSetType"))
"""
function set_type!(m::AbstractPetscMat{PetscLib}, type::Symbol) where {PetscLib}
    LibPETSc.MatSetType(getlib(PetscLib), m, String(type))
    return m
end

"""
    M::PetscMat = PetscMat(petsclib, S::SparseMatrixCSC; with_arrays = false)
    M::PetscMat = PetscMat(petsclib, comm, S::SparseMatrixCSC; with_arrays = false)

Creates a PetscMat object from a Julia SparseMatrixCSC `S` in sequential AIJ
format. `comm` defaults to `MPI.COMM_SELF`.

By default the entries are copied into storage PETSc allocates. With
`with_arrays = true` the CSR arrays converted from `S` are handed to PETSc and
borrowed rather than copied, which is what v0.4's `MatSeqAIJWithArrays` did; the
arrays are kept alive for as long as the matrix.

Replaces v0.4's `MatCreateSeqAIJ`: construction goes through the type
(docs/src/man/naming.md §5.1). The unrelated `MatSeqAIJ`, which allocated from
sizes, is now `PetscMat(petsclib, m, n, nnz)`.

# External Links
$(doc_external("Mat/MatCreateSeqAIJ"))
$(doc_external("Mat/MatCreateSeqAIJWithArrays"))
"""
LibPETSc.PetscMat(petsclib::PetscLibType, S::SparseMatrixCSC; kwargs...) =
    LibPETSc.PetscMat(petsclib, MPI.COMM_SELF, S; kwargs...)

function LibPETSc.PetscMat(
    petsclib::PetscLibType,
    comm::MPI.Comm,
    S::SparseMatrixCSC{PetscScalar};
    with_arrays::Bool = false,
) where {PetscScalar}
    check_initialized(petsclib)

    with_arrays && return mat_seqaij_with_arrays(petsclib, comm, S)

    PetscInt = petsclib.PetscInt

    # Set values from sparse matrix into PETSc Mat
    m, n = size(S)
    
    # Calculate non-zeros per row
    nnz = zeros(PetscInt, m)

    for r in S.rowval
        nnz[r] += 1
    end
    M = LibPETSc.MatCreateSeqAIJ(petsclib, comm, 
                PetscInt(m), 
                PetscInt(n), 
                PetscInt(0), nnz)

    for j in 1:n
        for ii in S.colptr[j]:(S.colptr[j + 1] - 1)
            i = S.rowval[ii]
            M[i, j] = S.nzval[ii]
        end
    end
    assemble!(M)
    # The garbage collector can only destroy it when no other rank takes part
    if MPI.Comm_size(comm) == 1
        finalizer(destroy!, M)
    end
    return M
end

"""
    mat = PetscMat(petsclib, num_rows, num_cols, nonzeros)

Create a PETSc serial sparse array using AIJ format (also known as a compressed
sparse row or CSR format) of size `num_rows X num_cols` with `nonzeros` per row

If `nonzeros` is an `Integer` the same number of non-zeros will be used for each
row, if `nonzeros` is a `Vector{PetscInt}` then one value must be specified for
each row.

Memory allocation is handled by PETSc and garbage collection can be used.

# External Links
$(doc_external("Mat/MatCreateSeqAIJ"))
"""
function LibPETSc.PetscMat(
    petsclib::PetscLib,
    num_rows::Integer,
    num_cols::Integer,
    nonzeros::Union{Integer, Vector},
) where {PetscLib <: PetscLibType}
    comm = MPI.COMM_SELF
    check_initialized(petsclib)
    PetscInt = petsclib.PetscInt
    if nonzeros isa Integer
        mat = LibPETSc.MatCreateSeqAIJ(petsclib, comm, 
                PetscInt(num_rows), 
                PetscInt(num_cols), 
                PetscInt(nonzeros), C_NULL)
    else
        eltype(nonzeros) === petsclib.PetscInt || throw(
            ArgumentError(
                "nonzeros has element type $(eltype(nonzeros)), " *
                "but the library uses $(petsclib.PetscInt)",
            ),
        )
        length(nonzeros) >= num_rows || throw(
            DimensionMismatch(
                "nonzeros has $(length(nonzeros)) entries, " *
                "but the matrix has $num_rows rows",
            ),
        )
        mat = LibPETSc.MatCreateSeqAIJ(petsclib, comm, 
                PetscInt(num_rows), 
                PetscInt(num_cols), 
                PetscInt(0), PetscInt.(nonzeros))    
    end

    finalizer(destroy!, mat)

    return mat
end

"""
    mat = PetscMat(petsclib, A::Matrix{PetscScalar})

PETSc dense array. This wraps a Julia `Matrix{PetscScalar}` object.

PETSc uses the memory of `A` as its storage, so `A` must be a `Matrix` of exactly the
library's scalar type; anything else throws an `ArgumentError`. To start from another
array, convert it first: `PetscMat(petsclib, Matrix{PetscScalar}(B))`.

Replaces v0.4's `MatSeqDense`: construction goes through the type
(docs/src/man/naming.md §5.1).

# External Links
$(doc_external("Mat/MatCreateSeqDense"))
"""
function LibPETSc.PetscMat(
    petsclib::PetscLib,
    A::Matrix{PetscScalar},
) where {PetscLib <: PetscLibType, PetscScalar}
    comm = MPI.COMM_SELF
    check_initialized(petsclib)
    PetscScalar === petsclib.PetscScalar || throw(
        ArgumentError(
            "matrix has element type $PetscScalar, " *
            "but the library uses $(petsclib.PetscScalar)",
        ),
    )
    
    PetscInt = petsclib.PetscInt
    # PETSc stores the data pointer directly without copying, so we must keep the
    # backing array alive for the entire lifetime of the PETSc Mat.  The finalizer
    # closure captures `data`, making it reachable (and therefore uncollectable) for
    # as long as `mat` is alive.
    data = vec(A)
    mat = LibPETSc.MatCreateSeqDense(petsclib, comm, PetscInt(size(A, 1)), PetscInt(size(A, 2)), data)

    finalizer(m -> (destroy!(m); data), mat)
    return mat
end

"""
    mat = PetscMat(petsclib, num_rows, num_cols; type = :seqaij)

An empty PETSc matrix of the given size.

`type` is the PETSc implementation name as a `Symbol`
(docs/src/man/naming.md §3.1): `:seqaij` (the default) preallocates nothing,
`:seqdense` (spelled `:dense` as well) allocates the dense storage. Flavour is a
keyword rather than a type parameter because it only affects construction, and
every variant returns a `PetscMat` (§9).

# External Links
$(doc_external("Mat/MatCreateSeqAIJ"))
$(doc_external("Mat/MatCreateSeqDense"))
"""
function LibPETSc.PetscMat(
    petsclib::PetscLib,
    num_rows::Integer,
    num_cols::Integer;
    type::Symbol = :seqaij,
) where {PetscLib <: PetscLibType}
    comm = MPI.COMM_SELF
    check_initialized(petsclib)
    PetscInt = petsclib.PetscInt
    PetscScalar = petsclib.PetscScalar
    if type === :dense || type === :seqdense
        data = zeros(PetscScalar, Int(num_rows) * Int(num_cols))
        mat = LibPETSc.MatCreateSeqDense(
            petsclib,
            comm,
            PetscInt(num_rows),
            PetscInt(num_cols),
            data,
        )
        finalizer(m -> (destroy!(m); data), mat)
        return mat
    elseif type === :seqaij || type === :aij
        return LibPETSc.PetscMat(petsclib, num_rows, num_cols, 0)
    else
        throw(
            ArgumentError(
                "unknown matrix type :$type, expected :seqaij or :dense",
            ),
        )
    end
end


# Matrix indexing - set single value
function Base.setindex!(m::AbstractPetscMat{PetscLib}, val, i::Integer, j::Integer) where {PetscLib}
    PetscInt = inttype(PetscLib)
    PetscScalar = scalartype(PetscLib)
    
    # Convert to 0-based indexing for PETSc (Julia uses 1-based)
    row = PetscInt(i - 1)
    col = PetscInt(j - 1)
    value = PetscScalar(val)
    
    # Use MatSetValues for single entry
    # MatSetValues(mat, m, idxm, n, idxn, y, INSERT_VALUES)
    LibPETSc.MatSetValues(PetscLib, m, PetscInt(1), [row], PetscInt(1), [col], [value], LibPETSc.INSERT_VALUES)
    
    return m
end

# Matrix indexing - set multiple values with vectors of indices  
function Base.setindex!(m::AbstractPetscMat{PetscLib}, vals::AbstractMatrix, rows::AbstractVector{<:Integer}, cols::AbstractVector{<:Integer}) where {PetscLib}
    PetscInt = inttype(PetscLib)
    PetscScalar = scalartype(PetscLib)
    
    # Convert to 0-based indexing for PETSc
    petsc_rows = PetscInt[r - 1 for r in rows]
    petsc_cols = PetscInt[c - 1 for c in cols]

    size(vals) == (length(rows), length(cols)) || throw(
        DimensionMismatch(
            "block is $(size(vals)) but the index ranges are " *
            "$(length(rows))x$(length(cols))",
        ),
    )

    # MatSetValues reads its value array row by row, so the block is transposed
    # before it is flattened: Julia lays a matrix out column by column.
    petsc_vals = vec(permutedims(PetscScalar.(vals)))

    nrows = PetscInt(length(petsc_rows))
    ncols = PetscInt(length(petsc_cols))

    LibPETSc.MatSetValues(PetscLib, m, nrows, petsc_rows, ncols, petsc_cols, petsc_vals, LibPETSc.INSERT_VALUES)

    return m
end

# Matrix indexing - set row with vector of values
function Base.setindex!(m::AbstractPetscMat{PetscLib}, vals::AbstractVector, i::Integer, cols::AbstractVector{<:Integer}) where {PetscLib}
    PetscInt = inttype(PetscLib)
    PetscScalar = scalartype(PetscLib)
    
    # Convert to 0-based indexing
    row = PetscInt(i - 1)
    petsc_cols = PetscInt[c - 1 for c in cols]
    petsc_vals = PetscScalar.(vals)
    
    ncols = PetscInt.(length(petsc_cols))
    LibPETSc.MatSetValues(PetscLib, m, PetscInt(1), [row], ncols, petsc_cols, petsc_vals, LibPETSc.INSERT_VALUES)
    
    return m
end

# Matrix indexing - set column with vector of values  
function Base.setindex!(m::AbstractPetscMat{PetscLib}, vals::AbstractVector, rows::AbstractVector{<:Integer}, j::Integer) where {PetscLib}
    PetscInt = inttype(PetscLib)
    PetscScalar = scalartype(PetscLib)
    
    # Convert to 0-based indexing
    petsc_rows = PetscInt[r - 1 for r in rows]
    col = PetscInt(j - 1)
    petsc_vals = PetscScalar.(vals)
    
    nrows = PetscInt.(length(petsc_rows))
    LibPETSc.MatSetValues(PetscLib, m, nrows, petsc_rows, PetscInt(1), [col], petsc_vals, LibPETSc.INSERT_VALUES)
    
    return m
end

# Matrix indexing - get single value
function Base.getindex(m::AbstractPetscMat{PetscLib}, i::Integer, j::Integer) where {PetscLib}
    PetscInt = inttype(PetscLib)
    PetscScalar = scalartype(PetscLib)
    
    # Convert to 0-based indexing for PETSc (Julia uses 1-based)
    row = PetscInt(i - 1)
    col = PetscInt(j - 1)
    
    # Use MatGetValues for single entry
    values = Vector{PetscScalar}(undef, 1)
    LibPETSc.MatGetValues(PetscLib, m, PetscInt(1), [row], PetscInt(1), [col], values)
    
    return values[1]
end

# Matrix indexing - get block of values
function Base.getindex(m::AbstractPetscMat{PetscLib}, rows::AbstractVector{<:Integer}, cols::AbstractVector{<:Integer}) where {PetscLib}
    PetscInt = inttype(PetscLib)
    PetscScalar = scalartype(PetscLib)
    
    # Convert to 0-based indexing for PETSc
    petsc_rows = PetscInt[r - 1 for r in rows]
    petsc_cols = PetscInt[c - 1 for c in cols]
    
    # Use MatGetValues for block of entries
    nrows = PetscInt.(length(petsc_rows))
    ncols = PetscInt.(length(petsc_cols))

    # PETSc returns values in row-major order
    values = Vector{PetscScalar}(undef, nrows * ncols)
    LibPETSc.MatGetValues(PetscLib, m, nrows, petsc_rows, ncols, petsc_cols, values)
    
    # Reshape to Julia matrix (column-major)
    result = Matrix{PetscScalar}(undef, nrows, ncols)
    for i in 1:nrows
        for j in 1:ncols
            # Convert from PETSc row-major to Julia column-major
            petsc_idx = (i-1) * ncols + (j-1) + 1
            result[i, j] = values[petsc_idx]
        end
    end
    
    return result
end

# Matrix indexing - get row values at specified columns
function Base.getindex(m::AbstractPetscMat{PetscLib}, i::Integer, cols::AbstractVector{<:Integer}) where {PetscLib}
    PetscInt = inttype(PetscLib)
    PetscScalar = scalartype(PetscLib)
    
    # Convert to 0-based indexing
    row = PetscInt(i - 1)
    petsc_cols = PetscInt[c - 1 for c in cols]

    ncols = PetscInt.(length(petsc_cols))
    values = Vector{PetscScalar}(undef, ncols)
    LibPETSc.MatGetValues(PetscLib, m, PetscInt(1), [row], ncols, petsc_cols, values)
    
    return values
end

# Matrix indexing - get column values at specified rows
function Base.getindex(m::AbstractPetscMat{PetscLib}, rows::AbstractVector{<:Integer}, j::Integer) where {PetscLib}
    PetscInt = inttype(PetscLib)
    PetscScalar = scalartype(PetscLib)
    
    # Convert to 0-based indexing
    petsc_rows = PetscInt[r - 1 for r in rows]
    col = PetscInt(j - 1)
    
    nrows = length(petsc_rows)
    values = Vector{PetscScalar}(undef, nrows)
    LibPETSc.MatGetValues(PetscLib, m, nrows, petsc_rows, PetscInt(1), [col], values)
    return values
end

# Matrix indexing - get entire row
function Base.getindex(m::AbstractPetscMat{PetscLib}, i::Integer, ::Colon) where {PetscLib}
    nrows, ncols = size(m)
    return getindex(m, i, 1:ncols)
end

# Matrix indexing - get entire column  
function Base.getindex(m::AbstractPetscMat{PetscLib}, ::Colon, j::Integer) where {PetscLib}
    nrows, ncols = size(m)
    return getindex(m, 1:nrows, j)
end

# Matrix indexing - get all values (use with caution for large matrices!)
function Base.getindex(m::AbstractPetscMat{PetscLib}, ::Colon, ::Colon) where {PetscLib}
    nrows, ncols = size(m)
    return getindex(m, 1:nrows, 1:ncols)
end

function Base.:(==)(
    A::PetscMat{PetscLib},
    B::PetscMat{PetscLib},
) where {PetscLib}
    return LibPETSc.MatEqual(PetscLib, A, B)
end

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
    zero_rows!(A::AbstractPetscMat, rows_0b, diag = 1; x = nothing, b = nothing)

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
    rows_0b::AbstractVector{<:Integer},
    diag = 1;
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
    zero_rows_local!(A::AbstractPetscMat, rows_0b, diag = 1; x = nothing, b = nothing)

[`zero_rows!`](@ref) with 0-based local row numbers, translated through the
local-to-global mapping of `A`. A matrix from `DMCreateMatrix` has that mapping;
one built without it throws a `PetscError`.

# External Links
$(doc_external("Mat/MatZeroRowsLocal"))
"""
function zero_rows_local!(
    A::AbstractPetscMat{PetscLib},
    rows_0b::AbstractVector{<:Integer},
    diag = 1;
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
    set_option!(A::AbstractPetscMat, option::LibPETSc.MatOption, flag::Bool)

Turn the matrix option `option` on or off and return `A`, e.g.
`set_option!(A, LibPETSc.MAT_NEW_NONZERO_ALLOCATION_ERR, false)` to allow
insertions outside the preallocated pattern.

# External Links
$(doc_external("Mat/MatSetOption"))
"""
function set_option!(
    A::AbstractPetscMat{PetscLib},
    option::LibPETSc.MatOption,
    flag::Bool,
) where {PetscLib}
    LibPETSc.MatSetOption(PetscLib, A, option, LibPETSc.PetscBool(flag))
    return A
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
    mul!(y::PetscVec{PetscLib}, M::AbstractPetscMat{PetscLib}, x::PetscVec{PetscLib})

Computes
    `y` = `M`*`x`

"""
function LinearAlgebra.mul!(y::PetscVec{PetscLib},M::AbstractPetscMat{PetscLib},x::PetscVec{PetscLib}) where {PetscLib} 
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
    mul!(y::PetscVec, A', x::PetscVec)
    mul!(y::PetscVec, transpose(A), x::PetscVec)

Compute `y = Aᴴx` or `y = Aᵀx` without forming the transpose, and return `y`.

# External Links
$(doc_external("Mat/MatMultHermitianTranspose"))
$(doc_external("Mat/MatMultTranspose"))
"""
function LinearAlgebra.mul!(
    y::PetscVec{PetscLib},
    M::TransposedPetscMat{PetscLib},
    x::PetscVec{PetscLib},
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
    setup!(mat::AbstractMat)

Set up the interal data for `mat`

# External Links
$(doc_external("Mat/MatSetUp"))
"""
function setup!(mat::PetscMat{PetscLib}) where {PetscLib}
    check_initialized(PetscLib)
    LibPETSc.MatSetUp(PetscLib, mat)
    return mat
end

"""
    rowptr, colval, nzval = csr_from_csc(petsclib, A::SparseMatrixCSC)

The CSR triple PETSc wants, built from Julia's CSC storage. `rowptr` and
`colval` are 0-based, as `MatCreateSeqAIJWithArrays` expects.
"""
function csr_from_csc(petsclib::PetscLibType, A::SparseMatrixCSC)
    PetscInt = PETSc.inttype(petsclib)
    PetscScalar = PETSc.scalartype(petsclib)
    
    m, n = size(A)
    
    # Convert Julia's CSC to CSR format manually
    # First, count non-zeros per row
    nnz_per_row = zeros(Int, m)
    for j in 1:n
        for k in A.colptr[j]:(A.colptr[j+1]-1)
            i = A.rowval[k]
            nnz_per_row[i] += 1
        end
    end
    
    # Build CSR arrays
    row_ptr = PetscInt[0]  # Start with 0 (PETSc uses 0-based indexing)
    for i in 1:m
        push!(row_ptr, row_ptr[end] + nnz_per_row[i])
    end
    
    # Pre-allocate arrays
    nnz_total = length(A.nzval)
    col_idx = Vector{PetscInt}(undef, nnz_total)
    values = Vector{PetscScalar}(undef, nnz_total)
    
    # Fill CSR arrays
    current_pos = copy(row_ptr[1:end-1]) .+ 1  # Track current position for each row
    
    for j in 1:n
        for k in A.colptr[j]:(A.colptr[j+1]-1)
            i = A.rowval[k]
            pos = current_pos[i]
            col_idx[pos] = j - 1  # Convert to 0-based indexing
            values[pos] = PetscScalar(A.nzval[k])
            current_pos[i] += 1
        end
    end
    
    return row_ptr, col_idx, values
end

"""
    B = PetscMat(petsclib, rowptr, colval, nzval; comm = MPI.COMM_SELF, ncols = …)

Create a PETSc SeqAIJ matrix on the CSR arrays `rowptr`, `colval` and `nzval`.
Arrays that are already `Vector`s of the library's integer and scalar types are
borrowed rather than copied; 
any other vector (another element type, a view, an offset array) is copied once.

`rowptr` and `colval` are 0-based, PETSc's own base for bulk index arrays
(docs/src/man/naming.md §12.1). The number of rows is `length(rowptr) - 1`;
`ncols` defaults to one past the largest column index.

The matrix keeps the arrays alive for as long as it exists, including when a
solver still holds it after `destroy!` on this handle.

Replaces v0.4's `MatSeqAIJWithArrays`, which took a `SparseMatrixCSC` and so
could not be told apart from `MatCreateSeqAIJ` (docs/src/man/naming.md §6).

# External Links
$(doc_external("Mat/MatCreateSeqAIJWithArrays"))
"""
function LibPETSc.PetscMat(
    petsclib::PetscLibType,
    rowptr::AbstractVector{<:Integer},
    colval::AbstractVector{<:Integer},
    nzval::AbstractVector;
    comm = MPI.COMM_SELF,
    ncols::Integer = isempty(colval) ? 0 : (maximum(colval) + 1),
)
    check_initialized(petsclib)
    PetscInt = inttype(petsclib)
    PetscScalar = scalartype(petsclib)

    row_ptr = c_vector(PetscInt, rowptr)
    col_idx = c_vector(PetscInt, colval)
    values = c_vector(PetscScalar, nzval)

    mat = LibPETSc.MatCreateSeqAIJWithArrays(
        petsclib,
        comm,
        PetscInt(length(row_ptr) - 1),
        PetscInt(ncols),
        row_ptr,
        col_idx,
        values,
    )

    keep_alive!(mat, (row_ptr, col_idx, values))
    # The garbage collector can only destroy it when no other rank takes part
    if MPI.Comm_size(comm) == 1
        finalizer(destroy!, mat)
    end
    return mat
end

"""
    mat_seqaij_with_arrays(petsclib, comm, A::SparseMatrixCSC)

The v0.4 `MatSeqAIJWithArrays` body, kept for its deprecation shim: converts
`A` to CSR and hands the arrays to the `PetscMat` constructor.
"""
function mat_seqaij_with_arrays(petsclib::PetscLibType, comm, A::SparseMatrixCSC)
    rowptr, colval, nzval = csr_from_csc(petsclib, A)
    return LibPETSc.PetscMat(
        petsclib,
        rowptr,
        colval,
        nzval;
        comm = comm,
        ncols = size(A, 2),
    )
end

"""
    destroy!(m::AbstractPetscMat)

Destroy a Mat (matrix) object and release associated resources.

This function is typically called automatically via finalizers when the object
is garbage collected, but can be called explicitly to free resources immediately.
Does nothing on a matrix that only borrows its handle: see [`owns`](@ref).

# External Links
$(doc_external("Mat/MatDestroy"))
"""
function destroy!(m::AbstractPetscMat{PetscLib}) where {PetscLib}
    owns(m) || return nothing
    if isdestroyable(m, PetscLib)
        LibPETSc.MatDestroy(PetscLib, m)
    end
    m.ptr = C_NULL
    return nothing
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
        rowvals::AbstractArray,
        insertmode::InsertMode = INSERT_VALUES;
        num_rows = length(rows_0b),
        num_cols = length(cols_0b)
    )

Set values of the matrix `M` with base-0 row and column indices `rows_0b` and
`cols_0b`, inserting the values `rowvals`, read row by row in linear order.
A `Vector{MatStencil}` and a `Vector` (or `Array`) of the library's scalar type 
go to PETSc as they are; anything else is copied first.

The `_0b` suffix says what the base is: a bulk index array handed to C keeps
PETSc's base rather than being rebuilt on a hot path (docs/src/man/naming.md
§12.1). `A[i, j] = v` is the 1-based route.

If the keyword arguments `num_rows` or `num_cols` is specified then only the
first `num_rows * num_cols` values of `rowvals` will be used.

# External Links
$(doc_external("Mat/MatSetValuesStencil"))
"""
function set_values!(
    M::AbstractPetscMat{PetscLib},
    rows_0b::AbstractVector{MatStencil},
    cols_0b::AbstractVector{MatStencil},
    rowvals::AbstractArray,
    insertmode::InsertMode = INSERT_VALUES;
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
        insertmode,
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


function Base.setindex!(
    M::AbstractPetscMat{PetscLib},
    val,
    i::CartesianIndex{N},
    j::CartesianIndex{N},
) where {PetscLib, N}
    PetscInt = PetscLib.PetscInt
    PetscScalar = PetscLib.PetscScalar
     ms_i = MatStencil(
        N < 3 ? PetscInt(0) : PetscInt(i[3] - 1),
        N < 2 ? PetscInt(0) : PetscInt(i[2] - 1),
        PetscInt(i[1] - 1),
        N < 4 ? PetscInt(0) : PetscInt(i[4] - 1),
    )
    ms_j = MatStencil(
        N < 3 ? PetscInt(0) : PetscInt(j[3] - 1),
        N < 2 ? PetscInt(0) : PetscInt(j[2] - 1),
        PetscInt(j[1] - 1),
        N < 4 ? PetscInt(0) : PetscInt(j[4] - 1),
    )
    set_values!(M, [ms_i], [ms_j], [PetscScalar(val)], INSERT_VALUES)
    return M
end

function add_index!(
    M::AbstractPetscMat{PetscLib},
    val,
    i::CartesianIndex{N},
    j::CartesianIndex{N},
) where {PetscLib, N}
    PetscInt = PetscLib.PetscInt
    PetscScalar = PetscLib.PetscScalar
    ms_i = MatStencil(
        N < 3 ? PetscInt(0) : PetscInt(i[3] - 1),
        N < 2 ? PetscInt(0) : PetscInt(i[2] - 1),
        PetscInt(i[1] - 1),
        N < 4 ? PetscInt(0) : PetscInt(i[4] - 1),
    )
    ms_j = MatStencil(
        N < 3 ? PetscInt(0) : PetscInt(j[3] - 1),
        N < 2 ? PetscInt(0) : PetscInt(j[2] - 1),
        PetscInt(j[1] - 1),
        N < 4 ? PetscInt(0) : PetscInt(j[4] - 1),
    )
    set_values!(M, [ms_i], [ms_j], [PetscScalar(val)], ADD_VALUES)
    return M
end


# ====

struct MatOp{PetscLib, Op} end

function (::MatOp{PetscLib, LibPETSc.MATOP_MULT})(
            M::CMat,
            cx::CVec,
            cy::CVec,
        ) where {PetscLib}
    return run_callback("MatShell multiply") do
        state = unsafe_pointer_to_objref(LibPETSc.MatShellGetContext(PetscLib, M))::MatShellState
        petsclib = getlib(PetscLib)
        _mul!(PetscVec(cy, petsclib; own = false), state.obj, PetscVec(cx, petsclib; own = false))
    end
end




# The macro is needed because of the @cfunction


"""
    MatShell(
        petsclib::PetscLib,
        obj::OType,
        comm::MPI.Comm,
        local_rows,
        local_cols,
        global_rows = LibPETSc.PETSC_DECIDE,
        global_cols = LibPETSc.PETSC_DECIDE,
    )

Create a `global_rows X global_cols` PETSc shell matrix object wrapping `obj`
with local size `local_rows X local_cols`.

The `obj` will be registered as an `MATOP_MULT` function and if if `obj` is a
`Function`, then the multiply action `obj(y,x)`; otherwise it calls `mul!(y,
obj, x)`.

if `comm == MPI.COMM_SELF` then the garbage connector can finalize the object,
otherwise the user is responsible for calling [`destroy!`](@ref).

$(doc_callback())

# External Links
$(doc_external("Mat/MatCreateShell"))
$(doc_external("Mat/MatShellSetOperation"))
$(doc_external("Mat/MATOP_MULT"))
"""
mutable struct MatShell{PetscLib, OType} <: AbstractPetscMat{PetscLib}
    ptr::CMat
    obj::OType
    age
    own::Bool
end

# The Julia side of a MatShell, kept with the PETSc object.
# It holds the operator, not the wrapper, so the wrapper can still be collected
# and its finalizer run.
mutable struct MatShellState <: ObjectState
    obj::Any
    alive::Bool
end
MatShellState() = MatShellState(nothing, true)
state_type(::Type{<:MatShell}) = MatShellState

LibPETSc.@for_petsc function MatShell(
    petsclib::$PetscLib,
    obj::OType,
    comm::MPI.Comm,
    local_rows,
    local_cols,
    global_rows = LibPETSc.PETSC_DECIDE,
    global_cols = LibPETSc.PETSC_DECIDE,
) where {OType}
    mat = MatShell{$PetscLib, OType}(C_NULL, obj, petsclib.age, true)

#=
    ccall(
        (:MatCreateShell, $petsc_library),
        LibPETSc.PetscErrorCode,
        (
            LibPETSc.MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Cvoid},
            Ptr{CMat},
        ),
        comm,
        local_rows,
        local_cols,
        global_rows,
        global_cols,
        pointer_from_objref(mat),
        A_,
    )
=#
    A_ = Ref{CMat}()
    ccall(
               (:MatCreateShell, $petsc_library),
               LibPETSc.PetscErrorCode,
               (LibPETSc.MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{Cvoid}, Ptr{CMat}),
               comm, local_rows, local_cols, global_rows, global_cols, C_NULL, A_,
              )

    mat.ptr = A_[]
    state = object_state!(mat)
    state.obj = obj
    LibPETSc.MatShellSetContext(petsclib, mat, pointer_from_objref(state))

    #=
     LibPETSc.MatCreateShell(
        petsclib,
        comm,
        local_rows,
        local_cols,
        global_rows,
        global_cols,
        pointer_from_objref(mat),
        mat,
    )


    mat = LibPETSc.MatCreateShell(
        petsclib,
        comm,
        local_rows,
        local_cols,
        global_rows,
        global_cols,
        pointer_from_objref(mat),
    )
    =#
  
    mulptr = @cfunction(
        MatOp{$PetscLib, LibPETSc.MATOP_MULT}(),
        LibPETSc.PetscErrorCode,
        (CMat, CVec, CVec)
    )

    LibPETSc.MatShellSetOperation(petsclib, mat, LibPETSc.MATOP_MULT, mulptr)

    # The garbage collector can only destroy it when no other rank takes part
    if MPI.Comm_size(comm) == 1
        finalizer(destroy!, mat)
    end

    return mat
end

# The operator of a MatShell applied to `x`: a function is called as `obj(y, x)`,
# anything else goes through `mul!(y, obj, x)`
_mul!(y, obj::Function, x) = obj(y, x)
_mul!(y, obj, x) = LinearAlgebra.mul!(y, obj, x)


# Matrix-vector multiplication for MatShell with Julia arrays
function Base.:*(M::MatShell{PetscLib}, x::AbstractVector) where {PetscLib}
    PetscScalar = scalartype(PetscLib)
    PetscInt = inttype(PetscLib)
    
    # Get matrix dimensions
    m, n = size(M)
    
    # Create PETSc vectors wrapping the Julia arrays
    petsclib = getlib(PetscLib)
    # PETSc works on this copy of `x` in place, so it must outlive the product
    x_copy = PetscScalar.(x)
    return GC.@preserve x_copy begin
        petsc_x = seq_vec_with_array(petsclib, MPI.COMM_SELF, PetscInt(1), PetscInt(n), x_copy)
        petsc_y = LibPETSc.VecCreateSeq(petsclib, MPI.COMM_SELF, PetscInt(m))
        mul!(petsc_y, M, petsc_x)
        result = petsc_y[:]
        destroy!(petsc_x)
        destroy!(petsc_y)
        result
    end
end

"""
    ownership_range(mat::AbstractPetscMat)

The range of row indices owned by this processor, assuming that the `mat` is
laid out with the first `n1` rows on the first processor, next `n2` rows on the
second, etc. For certain parallel layouts this range may not be well defined.

The range is **1-based**, always: an index into Julia data is 1-based
(docs/src/man/naming.md §12.1). v0.4 took `base_one::Bool` positionally and
made the convention a runtime choice; `ownership_range(A, false)` warns in
v0.5 and is a `MethodError` in v0.6.

!!! note

    unlike the C function, the range returned is inclusive (`idx_first:idx_last`)

# External Links
$(doc_external("Mat/MatGetOwnershipRange"))
"""
function ownership_range(mat::AbstractPetscMat{PetscLib}) where {PetscLib}
    PetscInt = PetscLib.PetscInt
    # The wrapper returns two plain integers, not `Ref`s.
    r_lo, r_hi = LibPETSc.MatGetOwnershipRange(PetscLib, mat)
    return (r_lo + PetscInt(1)):r_hi
end

"""
    set_values!(
        M::AbstractPetscMat,
        rows_0b::AbstractVector{<:Integer},
        cols_0b::AbstractVector{<:Integer},
        rowvals::AbstractArray,
        insertmode::InsertMode = INSERT_VALUES;
        num_rows = length(rows_0b),
        num_cols = length(cols_0b)
    )

Set values of the matrix `M` with base-0 row and column indices `rows_0b` and
`cols_0b`, inserting the values `rowvals`, read row by row in linear order.
A `Vector` of the library's integer type for the indices, and a `Vector` (or `Array`)
of its scalar type for the values, go to PETSc as they are; anything else is
copied first.

The `_0b` suffix says what the base is: a bulk index array handed to C keeps
PETSc's base rather than being rebuilt on a hot path (docs/src/man/naming.md
§12.1). `A[i, j] = v` is the 1-based route.

If the keyword arguments `num_rows` or `num_cols` is specified then only the
first `num_rows * num_cols` values of `rowvals` will be used.

# External Links
$(doc_external("Mat/MatSetValues"))
"""
function set_values!(
    M::AbstractPetscMat{PetscLib},
    rows_0b::AbstractVector{<:Integer},
    cols_0b::AbstractVector{<:Integer},
    rowvals::AbstractArray,
    insertmode::InsertMode = INSERT_VALUES;
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
        insertmode,
    )
    return M
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

# ====