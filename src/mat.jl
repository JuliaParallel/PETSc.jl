
import .LibPETSc: AbstractPetscMat, PetscMat, CMat

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
function assemble!(A::PetscMat{PetscLib}) where {PetscLib}
    LibPETSc.MatAssemblyBegin(PetscLib, A, PETSc.MAT_FINAL_ASSEMBLY)
    LibPETSc.MatAssemblyEnd(PetscLib, A, PETSc.MAT_FINAL_ASSEMBLY)
end


LinearAlgebra.norm(M::PetscMat{PetscLib}, normtype::NormType = NORM_FROBENIUS) where {PetscLib} = LibPETSc.MatNorm(PetscLib, M, normtype)
function LinearAlgebra.mul!(y::PetscVec{PetscLib},M::AbstractPetscMat{PetscLib},x::PetscVec{PetscLib}) where {PetscLib} 
    LibPETSc.MatMult(PetscLib, M, x, y)
    return nothing
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
        row_ptr,
        col_idx,
        values
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


# ====


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
mutable struct MatShell{PetscLib, OType} <:
               AbstractPetscMat{PetscLib}
    ptr::CMat
    obj::OType
end

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
    x = VecPtr(PetscLib, cx, false)
    y = VecPtr(PetscLib, cy, false)

    _mul!(y, mat, x)

    return PetscInt(0)
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


# NOTE: MatShell remains work in progress - doesn't function yet
# We have to use the macro here because of the @cfunction
LibPETSc.@for_petsc function MatShell(
    petsclib::$PetscLib,
    obj::OType,
    comm::MPI.Comm,
    local_rows,
    local_cols,
    global_rows = LibPETSc.PETSC_DECIDE,
    global_cols = LibPETSc.PETSC_DECIDE,
) where {OType}
    mat = MatShell{$PetscLib, OType}(C_NULL, obj)

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

    if MPI.Comm_size(comm) == 1
        finalizer(destroy, mat)
    end
    
    return mat
end


# ====