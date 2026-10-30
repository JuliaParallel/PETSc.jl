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

# The REPL shows the handle, not the entries, which are one C call each
Base.show(io::IO, ::MIME"text/plain", A::AbstractPetscMat) = show(io, A)

# `==` on matrices compares entries, which can change
Base.hash(::AbstractPetscMat, ::UInt) = throw(ArgumentError(
    "a PETSc matrix has no hash: `==` compares its entries, which can change. " *
    "Key a Dict by the handle with an IdDict instead."
))

"""
    MatPtr(petsclib, ptr::CMat; own::Bool = false)

Container type for a PETSc Mat that is just a raw pointer.

If `own` is `true` a finalizer is set on the matrix, but only on a serial
communicator, since `MatDestroy` is collective and a GC finalizer runs at an
arbitrary point. If `own` is `false`, the default, the handle belongs to PETSc and `destroy!`
is a no-op, leaving the wrapper usable.

Calling `MatPtr` is deprecated and goes in v0.6: wrap a raw handle with
[`PetscMat(petsclib, ptr; own)`](@ref LibPETSc.PetscMat(::PetscLibType, ::CMat)).
"""
mutable struct MatPtr{PetscLib, T} <:
               AbstractPetscMat{PetscLib, T}
    ptr::CMat
    age::Int
    own::Bool

    MatPtr{PetscLib}(ptr::CMat, age::Int, own::Bool) where {PetscLib} =
        new{PetscLib, scalartype(PetscLib)}(ptr, age, own)
end
function MatPtr(petsclib::PetscLibType, ptr::CMat; own::Bool = false)
    @warn "MatPtr goes in v0.6: wrap a raw Mat with PetscMat(petsclib, ptr; own)" maxlog = 1
    return mat_from_ptr(petsclib, ptr, own)
end

function mat_from_ptr(petsclib::PetscLib, ptr::CMat, own::Bool) where {PetscLib <: PetscLibType}
    m = MatPtr{PetscLib}(ptr, petsclib.age, own)
    # Short-circuits on a borrowed handle, which is the hot path: callbacks wrap
    # PETSc-owned matrices on every invocation and never need the communicator.
    if own && MPI.Comm_size(LibPETSc.PetscObjectGetComm(getlib(PetscLib), m)) == 1
        finalizer(destroy!, m)
    end
    return m
end
MatPtr(::Type{PetscLib}, x...; kwargs...) where {PetscLib <: PetscLibType} =
    MatPtr(getlib(PetscLib), x...; kwargs...)

"""
    PetscMat(petsclib, ptr::CMat; own::Bool = false)

Wrap a raw PETSc `Mat` handle.

`own = false` is the borrowed case of docs/src/man/naming.md §3.3: no finalizer
is attached and `destroy!` is a no-op. The result is a [`MatPtr`](@ref), the
wrapper type that carries the ownership flag.
"""
LibPETSc.PetscMat(petsclib::PetscLibType, ptr::CMat; own::Bool = false) =
    mat_from_ptr(petsclib, ptr, own)

# Global size, as `Int`s, like any `AbstractMatrix`
Base.size(m::AbstractPetscMat{PetscLib}) where {PetscLib} = map(Int, LibPETSc.MatGetSize(PetscLib, m))
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
    # PETSc uses `A`'s memory without copying it, and a KSP can hold the Mat
    # after this handle is gone, so the array lives with the PETSc object
    data = vec(A)
    mat = LibPETSc.MatCreateSeqDense(petsclib, comm, PetscInt(size(A, 1)), PetscInt(size(A, 2)), data)
    keep_alive!(mat, data)

    finalizer(destroy!, mat)
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
    if type === :dense || type === :seqdense
        # Null data: PETSc allocates zeroed storage and frees it with the Mat
        mat = LibPETSc.MatCreateSeqDense(
            petsclib,
            comm,
            PetscInt(num_rows),
            PetscInt(num_cols),
            C_NULL,
        )
        finalizer(destroy!, mat)
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
    B = PetscMat(petsclib, rowptr, colval, nzval; ncols = …)
    B = PetscMat(petsclib, comm, rowptr, colval, nzval; ncols = …)

Create a PETSc SeqAIJ matrix on the CSR arrays `rowptr`, `colval` and `nzval`,
on `comm` (`MPI.COMM_SELF` when it is left out).
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
    comm = nothing,
    kwargs...,
)
    # TODO: remove the `comm` keyword in v0.6 with the other 0.5 deprecations.
    if comm === nothing
        comm = MPI.COMM_SELF
    else
        @warn "PetscMat(petsclib, rowptr, colval, nzval; comm) is deprecated, " *
              "use PetscMat(petsclib, comm, rowptr, colval, nzval)" maxlog = 1
    end
    return LibPETSc.PetscMat(petsclib, comm, rowptr, colval, nzval; kwargs...)
end

function LibPETSc.PetscMat(
    petsclib::PetscLibType,
    comm::MPI.Comm,
    rowptr::AbstractVector{<:Integer},
    colval::AbstractVector{<:Integer},
    nzval::AbstractVector;
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
    return LibPETSc.PetscMat(petsclib, comm, rowptr, colval, nzval; ncols = size(A, 2))
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
destroy!(m::AbstractPetscMat{PetscLib}) where {PetscLib} =
    destroy_handle!(LibPETSc.MatDestroy, m, PetscLib)

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
