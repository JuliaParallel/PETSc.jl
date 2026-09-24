
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
Base.ndims(m::AbstractPetscMat{PetscLib}) where {PetscLib} = length(LibPETSc.MatGetSize(PetscLib,m))
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
Base.axes(m::PetscMat{PetscLib}, i::Integer) where {PetscLib} = Base.OneTo(Base.size(m)[i])

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


LinearAlgebra.norm(M::PetscMat{PetscLib}, normtype::NormType = NORM_FROBENIUS) where {PetscLib} = LibPETSc.MatNorm(PetscLib, M, normtype)

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

function LinearAlgebra.mul!(
    y::PetscVec{PetscLib},
    M::Adjoint{AM},
    x::PetscVec{PetscLib},
) where {PetscLib, AM <: PetscMat{PetscLib}}
    LibPETSc.MatMultHermitianTranspose(PetscLib, parent(M), x, y)
    return y
end

function LinearAlgebra.mul!(
    y::PetscVec{PetscLib},
    M::Transpose{AM},
    x::PetscVec{PetscLib},
) where {PetscLib, AM <: PetscMat{PetscLib}}
    LibPETSc.MatMultTranspose(PetscLib, parent(M), x, y)
    return y
end

function LinearAlgebra.issymmetric(A::PetscMat{PetscLib}; tol = 0.0) where {PetscLib} 
    PetscReal = real(scalartype(PetscLib))        
    return LibPETSc.MatIsSymmetric(PetscLib, A, PetscReal(tol))
end
function LinearAlgebra.ishermitian(A::PetscMat{PetscLib}; tol = 0.0) where {PetscLib} 
    PetscReal = scalartype(PetscLib)
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

Create a PETSc SeqAIJ matrix directly on the CSR arrays `rowptr`, `colval` and
`nzval`, which PETSc borrows rather than copies.

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
    rowptr::Vector{<:Integer},
    colval::Vector{<:Integer},
    nzval::Vector;
    comm = MPI.COMM_SELF,
    ncols::Integer = isempty(colval) ? 0 : (maximum(colval) + 1),
)
    check_initialized(petsclib)
    PetscInt = inttype(petsclib)
    PetscScalar = scalartype(petsclib)

    row_ptr = convert(Vector{PetscInt}, rowptr)
    col_idx = convert(Vector{PetscInt}, colval)
    values = convert(Vector{PetscScalar}, nzval)

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

function LinearAlgebra.mul!(
    y::PetscVec{PetscLib},
    M::Transpose{PetscScalar, AM},
    x::PetscVec{PetscLib},
) where {PetscLib, PetscScalar, AM <: AbstractPetscMat{PetscLib}}
    LibPETSc.MatMultTranspose(PetscLib, parent(M), x, y)
    return y
end

function Base.copyto!(
    M::AbstractPetscMat{PetscLib},
    S::SparseMatrixCSC,
) where {PetscLib}
    row_rng = LibPETSc.MatGetOwnershipRange(PetscLib,M)
    PetscInt = PetscLib.PetscInt
    PetscScalar = PetscLib.PetscScalar
    row_start = row_rng[1]
    _, n = size(S)
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
        M::AbstractPetscMat{PetscLib},
        rows_0b::Vector{MatStencil},
        cols_0b::Vector{MatStencil},
        rowvals::Array{PetscScalar},
        insertmode::InsertMode = INSERT_VALUES;
        num_rows = length(rows_0b),
        num_cols = length(cols_0b)
    )

Set values of the matrix `M` with base-0 row and column indices `rows_0b` and
`cols_0b`, inserting the values `rowvals`.

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
    rows_0b::Vector{MatStencil},
    cols_0b::Vector{MatStencil},
    rowvals::Array{PetscScalar},
    insertmode::InsertMode = INSERT_VALUES;
    num_rows = length(rows_0b),
    num_cols = length(cols_0b),
) where {PetscLib, PetscScalar}
    PetscScalar === PetscLib.PetscScalar || throw(
        ArgumentError(
            "values have element type $PetscScalar, but the library uses " *
            "$(PetscLib.PetscScalar)",
        ),
    )
    num_rows * num_cols <= length(rowvals) || throw(
        DimensionMismatch(
            "a $(num_rows)x$(num_cols) block needs $(num_rows * num_cols) values, " *
            "got $(length(rowvals))",
        ),
    )
    LibPETSc.MatSetValuesStencil(
        PetscLib,
        M,
        num_rows,
        rows_0b,
        num_cols,
        cols_0b,
        rowvals,
        insertmode,
    )
    return M
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
        petsc_x = LibPETSc.VecCreateSeqWithArray(petsclib, MPI.COMM_SELF, PetscInt(1), PetscInt(n), x_copy)
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
        M::AbstractMat{PetscLib},
        rows_0b::Vector{PetscInt},
        cols_0b::Vector{PetscInt},
        rowvals::Array{PetscScalar},
        insertmode::InsertMode = INSERT_VALUES;
        num_rows = length(rows_0b),
        num_cols = length(cols_0b)
    )

Set values of the matrix `M` with base-0 row and column indices `rows_0b` and
`cols_0b`, inserting the values `rowvals`.

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
    rows_0b::Vector{PetscInt},
    cols_0b::Vector{PetscInt},
    rowvals::Array{PetscScalar},
    insertmode::InsertMode = INSERT_VALUES;
    num_rows = length(rows_0b),
    num_cols = length(cols_0b),
) where {PetscLib, PetscScalar, PetscInt}
    PetscScalar === PetscLib.PetscScalar || throw(
        ArgumentError(
            "values have element type $PetscScalar, " *
            "but the library uses $(PetscLib.PetscScalar)",
        ),
    )
    PetscInt === PetscLib.PetscInt || throw(
        ArgumentError(
            "indices have element type $PetscInt, " *
            "but the library uses $(PetscLib.PetscInt)",
        ),
    )
    num_rows * num_cols <= length(rowvals) || throw(
        DimensionMismatch(
            "a $(num_rows)x$(num_cols) block needs $(num_rows * num_cols) values, " *
            "got $(length(rowvals))",
        ),
    )
    LibPETSc.MatSetValues(
        PetscLib,
        M,
        PetscInt(num_rows),
        PetscInt.(rows_0b),
        PetscInt(num_cols),
        PetscInt.(cols_0b),
        rowvals,
        insertmode,
    )
    return M
end

function LinearAlgebra.norm(
    M::AbstractPetscMat{PetscLib},
    normtype::NormType = NORM_FROBENIUS,
) where {PetscLib}
    PetscReal = PetscLib.PetscReal
    #r_val = Ref{PetscReal}()
    r_val = LibPETSc.MatNorm(PetscLib, M, normtype)
    return r_val
end

# ====