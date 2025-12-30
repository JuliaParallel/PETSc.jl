
import .LibPETSc: AbstractPetscMat, PetscMat, CMat, MatStencil, InsertMode

# Custom display for REPL
function Base.show(io::IO, v::AbstractPetscMat{PetscLib}) where {PetscLib}
    if v.ptr == C_NULL
        print(io, "PETSc Mat (null pointer)")
        return
    else
        print(io, "PETSc $(type(v)) Mat of size $(size(v))")
    end
    return nothing
end

"""
    MatPtr(petsclib, mat::CMat)

Container type for a PETSc Mat that is just a raw pointer.
"""
mutable struct MatPtr{PetscLib} <:
               AbstractPetscMat{PetscLib}
    ptr::CMat
end

Base.size(m::AbstractPetscMat{PetscLib}) where {PetscLib} = LibPETSc.MatGetSize(PetscLib,m)
Base.length(m::AbstractPetscMat{PetscLib}) where {PetscLib} = prod(size(m))
Base.ndims(m::AbstractPetscMat{PetscLib}) where {PetscLib} = length(LibPETSc.MatGetSize(PetscLib,m))
type(m::AbstractPetscMat{PetscLib}) where {PetscLib} = LibPETSc.MatGetType(PetscLib,m)
Base.axes(m::PetscMat{PetscLib}, i::Integer) where {PetscLib} = Base.OneTo(Base.size(m)[i])


"""
    M::PetscMat = MatCreateSeqAIJ(petsclib, comm, S)

Creates a PetscMat object from a Julia SparseMatrixCSC `S` in sequential AIJ format.

"""
function MatCreateSeqAIJ(petsclib, comm, S::SparseMatrixCSC{PetscScalar})   where {PetscScalar}

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
    return M
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
    petsc_vals = PetscScalar.(vals)
    
    # Use MatSetValues for block of entries
    nrows = PetscInt.(length(petsc_rows))
    ncols = PetscInt.(length(petsc_cols))
    
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
end


LinearAlgebra.norm(M::PetscMat{PetscLib}, normtype::NormType = NORM_FROBENIUS) where {PetscLib} = LibPETSc.MatNorm(PetscLib, M, normtype)
function LinearAlgebra.mul!(y::PetscVec{PetscLib},M::AbstractPetscMat{PetscLib},x::PetscVec{PetscLib}) where {PetscLib} 
    LibPETSc.MatMult(PetscLib, M, x, y)
    return nothing
end

function Base.:*(
    M::AbstractPetscMat{PetscLib},
    x::AbstractPetscVec{PetscLib},
) where {PetscLib}
    y = LibPETSc.VecDuplicate(getlib(PetscLib), x)
    mul!(y, M, x)
    return y
end

function LinearAlgebra.mul!(
    y::PetscVec{PetscLib},
    M::Adjoint{AM},
    x::PetscVec{PetscLib},
) where {PetscLib, AM <: PetscMat{PetscLib}}
    LibPETSc.MatMultHermitianTranspose(PetscLib, parent(M), x, y)
    return nothing
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
$(_doc_external("Mat/MatSetUp"))
"""
function setup!(mat::PetscMat{PetscLib}) where {PetscLib}
    @assert initialized(PetscLib)
    LibPETSc.MatSetUp(PetscLib, mat)
    return mat
end

destroy(m::AbstractPetscMat{PetscLib}) where {PetscLib} = LibPETSc.MatDestroy(PetscLib,m)

"""
    B = MatSeqAIJWithArrays(petsclib, comm, A::SparseMatrixCSC)

Create a PETSc SeqAIJ matrix from a Julia SparseMatrixCSC.
Since Julia uses CSC and PETSc AIJ uses CSR, we convert the format properly.

# Arguments
- `petsclib`: PETSc library instance  
- `comm`: MPI communicator
- `A`: Julia sparse matrix in CSC format

# Returns
- `B`: PETSc matrix in AIJ (CSR) format
"""
function MatSeqAIJWithArrays(petsclib::PetscLibType, comm, A::SparseMatrixCSC{T}) where {T}
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
    
    # Create the PETSc matrix
    mat = LibPETSc.MatCreateSeqAIJWithArrays(
        petsclib,
        comm,
        PetscInt(m),
        PetscInt(n), 
        PetscInt.(row_ptr),
        PetscInt.(col_idx),
        PetscScalar.(values)
    )
    
    return mat
end

const MatAT{PetscLib, PetscScalar} = Union{
    PetscMat{PetscLib},
    Transpose{PetscScalar, <:PetscMat{PetscLib}},
    Adjoint{PetscScalar, <:PetscMat{PetscLib}},
}

function LinearAlgebra.mul!(
    y::PetscVec{PetscLib},
    M::Transpose{PetscScalar, AM},
    x::PetscVec{PetscLib},
) where {PetscLib, PetscScalar, AM <: AbstractPetscMat{PetscLib}}
    LibPETSc.MatMultTranspose(PetscLib, parent(M), x, y)
    return y
end

function Base.copyto!(
    M::PetscMat{PetscLib},
    S::SparseMatrixCSC,
) where {PetscLib}
    row_rng = LibPETSc.MatGetOwnershipRange(PetscLib,M)
    row_start = row_rng[1]
    _, n = size(S)
    for j in 1:n
        for ii in S.colptr[j]:(S.colptr[j + 1] - 1)
            i = S.rowval[ii]
            M[i + row_start, j + row_start] = S.nzval[ii]
        end
    end
end

"""
    setvalues!(
        M::AbstractPetscMat{PetscLib},
        row0idxs::Vector{MatStencil},
        col0idxs::Vector{MatStencil},
        rowvals::Array{PetscScalar},
        insertmode::InsertMode = INSERT_VALUES;
        num_rows = length(row0idxs),
        num_cols = length(col0idxs)
    )

Set values of the matrix `M` with base-0  row and column indices `row0idxs` and
`col0idxs` inserting the values `rowvals`.

If the keyword arguments `num_rows` or `num_cols` is specified then only the
first `num_rows * num_cols` values of `rowvals` will be used.

# External Links
$(_doc_external("Mat/MatSetValuesStencil"))
"""
function setvalues!(
    M::AbstractPetscMat{PetscLib},
    row0idxs::Vector{MatStencil},
    col0idxs::Vector{MatStencil},
    rowvals::Array{PetscScalar},
    insertmode::InsertMode = INSERT_VALUES;
    num_rows = length(row0idxs),
    num_cols = length(col0idxs),
) where {PetscLib, PetscScalar}
    @assert PetscScalar == PetscLib.PetscScalar
    @assert num_rows * num_cols <= length(rowvals)
    LibPETSc.MatSetValuesStencil(
        PetscLib,
        M,
        num_rows,
        row0idxs,
        num_cols,
        col0idxs,
        rowvals,
        insertmode,
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
    setvalues!(M, [ms_i], [ms_j], [PetscScalar(val)], INSERT_VALUES)
    return val
end

function addindex!(
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
    setvalues!(M, [ms_i], [ms_j], [PetscScalar(val)], ADD_VALUES)
    return val
end


# ====

struct MatOp{PetscLib, Op} end

function (::MatOp{PetscLib, LibPETSc.MATOP_MULT})(
            M::CMat,
            cx::CVec,
            cy::CVec,
        ) where {PetscLib}
    r_ctx = Ref{Ptr{Cvoid}}()
    LibPETSc.MatShellGetContext(PetscLib, M, r_ctx)
    ptr = r_ctx[]

    mat = unsafe_pointer_to_objref(ptr)

    PetscScalar = PetscLib.PetscScalar
    x = PetscVec(cx, getlib(PetscLib))
    y = PetscVec(cy, getlib(PetscLib))

    _mul!(y, mat, x)

    return 0
end




# NOTE: MatShell remains work in progress - doesn't function yet
# We have to use the macro here because of the @cfunction


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
otherwise the user is responsible for calling [`destroy`](@ref).

# External Links
$(_doc_external("Mat/MatCreateShell"))
$(_doc_external("Mat/MatShellSetOperation"))
$(_doc_external("Mat/MATOP_MULT"))
"""
mutable struct MatShell{PetscLib, OType} <: AbstractPetscMat{PetscLib}
    ptr::CMat
    obj::OType
    age
end

LibPETSc.@for_petsc function MatShell(
    petsclib::$PetscLib,
    obj::OType,
    comm::MPI.Comm,
    local_rows,
    local_cols,
    global_rows = LibPETSc.PETSC_DECIDE,
    global_cols = LibPETSc.PETSC_DECIDE,
) where {OType}
    mat = MatShell{$PetscLib, OType}(C_NULL, obj, 0)

    # we use the MatShell object itself
    ctx = pointer_from_objref(mat)

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
               comm, local_rows, local_cols, global_rows, global_cols, ctx, A_,
              )

    mat.ptr = A_[]  

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
        $PetscInt,
        (CMat, CVec, CVec)
    )

    LibPETSc.MatShellSetOperation(petsclib, mat, LibPETSc.MATOP_MULT, mulptr)

    #if MPI.Comm_size(comm) == 1
    #    finalizer(destroy, mat)
    #end
    
    return mat
end

function _mul!(
    y,
    mat::MatShell{PetscLib, F},
    x,
) where {PetscLib, F <: Function}
    mat.obj(y, x)
end

function _mul!(y, mat::MatShell, x)
    LinearAlgebra.mul!(y, mat.obj, x)
end


# Matrix-vector multiplication for MatShell with Julia arrays
function Base.:*(M::MatShell{PetscLib}, x::AbstractVector) where {PetscLib}
    PetscScalar = scalartype(PetscLib)
    PetscInt = inttype(PetscLib)
    
    # Get matrix dimensions
    m, n = size(M)
    
    # Create PETSc vectors wrapping the Julia arrays
    petsclib = getlib(PetscLib)
    petsc_x = LibPETSc.VecCreateSeqWithArray(petsclib, MPI.COMM_SELF, PetscInt(1), PetscInt(n), PetscScalar.(x))
    petsc_y = LibPETSc.VecCreateSeq(petsclib, MPI.COMM_SELF, PetscInt(m))
    
    # Perform matrix-vector multiplication
    LibPETSc.MatMult(petsclib, M, petsc_x, petsc_y)
    
    # Extract result to Julia array
    result = PetscScalar.(petsc_y[:])
    
    # Clean up
    PETSc.destroy(petsc_x)
    PETSc.destroy(petsc_y)
    
    return result
end

"""
    ownershiprange(mat::AbstractMat, [base_one = true])

The range of row indices owned by this processor, assuming that the `mat` is
laid out with the first `n1` rows on the first processor, next `n2` rows on the
second, etc. For certain parallel layouts this range may not be well defined.

If the optional argument `base_one == true` then base-1 indexing is used,
otherwise base-0 index is used.

!!! note

    unlike the C function, the range returned is inclusive (`idx_first:idx_last`)

# External Links
$(_doc_external("Mat/MatGetOwnershipRange"))
"""
function ownershiprange(
    mat::AbstractPetscMat{PetscLib},
    base_one::Bool = true,
) where {PetscLib}
    PetscInt = PetscLib.PetscInt
    r_lo, r_hi = LibPETSc.MatGetOwnershipRange(PetscLib, mat)
    return base_one ? ((r_lo[] + PetscInt(1)):(r_hi[])) :
           ((r_lo[]):(r_hi[] - PetscInt(1)))
end

"""
    setvalues!(
        M::AbstractMat{PetscLib},
        row0idxs::Vector{PetscInt},
        col0idxs::Vector{PetscInt},
        rowvals::Array{PetscScalar},
        insertmode::InsertMode = INSERT_VALUES;
        num_rows = length(row0idxs),
        num_cols = length(col0idxs)
    )

Set values of the matrix `M` with base-0  row and column indices `row0idxs` and
`col0idxs` inserting the values `rowvals`.

If the keyword arguments `num_rows` or `num_cols` is specified then only the
first `num_rows * num_cols` values of `rowvals` will be used.

# External Links
$(_doc_external("Mat/MatSetValues"))
"""
function setvalues!(
    M::AbstractPetscMat{PetscLib},
    row0idxs::Vector{PetscInt},
    col0idxs::Vector{PetscInt},
    rowvals::Array{PetscScalar},
    insertmode::InsertMode = INSERT_VALUES;
    num_rows = length(row0idxs),
    num_cols = length(col0idxs),
) where {PetscLib, PetscScalar, PetscInt}
    @assert PetscScalar == PetscLib.PetscScalar
    @assert PetscInt == PetscLib.PetscInt
    @assert num_rows * num_cols <= length(rowvals)
    LibPETSc.MatSetValues(
        PetscLib,
        M,
        PetscInt(num_rows),
        PetscInt.(row0idxs),
        PetscInt(num_cols),
        PetscInt.(col0idxs),
        rowvals,
        insertmode,
    )
    return nothing
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