const CMat = Ptr{Cvoid}

abstract type AbstractMat{PetscLib, PetscScalar} <: AbstractMatrix{PetscScalar} end

Base.eltype(
    ::Type{V},
) where {
    V <: AbstractMat{PetscLib, PetscScalar},
} where {PetscLib, PetscScalar} = PetscScalar
Base.eltype(
    v::AbstractMat{PetscLib, PetscScalar},
) where {PetscLib, PetscScalar} = PetscScalar

function destroy(M::AbstractMat{PetscLib}) where {PetscLib}
    finalized(PetscLib) || LibPETSc.MatDestroy(PetscLib, M)
    return nothing
end

"""
    MatSeqAIJ(petsclib, num_rows, num_cols, nonzeros)

Create a PETSc serial sparse array using AIJ format (also known as a compressed
sparse row or CSR format) of size `num_rows X num_cols` with `nonzeros` per row

If `nonzeros` is an `Integer` the same number of non-zeros will be used for each
row, if `nonzeros` is a `Vector{PetscInt}` then one value must be specified for
each row.

Memory allocation is handled by PETSc and garbage collection can be used.

# External Links
$(_doc_external("Mat/MatCreateSeqAIJ"))
"""
mutable struct MatSeqAIJ{PetscLib, PetscScalar} <:
               AbstractMat{PetscLib, PetscScalar}
    ptr::CMat
end

function MatSeqAIJ(
    petsclib::PetscLib,
    num_rows::Integer,
    num_cols::Integer,
    nonzeros::Union{Integer, Vector},
) where {PetscLib <: PetscLibType}
    comm = MPI.COMM_SELF
    @assert initialized(petsclib)
    PetscScalar = petsclib.PetscScalar
    mat = MatSeqAIJ{PetscLib, PetscScalar}(C_NULL)
    if nonzeros isa Integer
        LibPETSc.MatCreateSeqAIJ(
            petsclib,
            comm,
            num_rows,
            num_cols,
            nonzeros,
            C_NULL,
            mat,
        )
    else
        @assert eltype(nonzeros) == petsclib.PetscInt
        @assert length(nonzeros) >= num_rows
        LibPETSc.MatCreateSeqAIJ(
            petsclib,
            comm,
            num_rows,
            num_cols,
            0,
            nonzeros,
            mat,
        )
    end

    finalizer(destroy, mat)

    return mat
end

"""
    MatSeqDense{PetscLib, PetscScalar}

PETSc dense array. This wraps a Julia `Matrix{PetscScalar}` object.

# External Links
$(_doc_external("Mat/MatCreateSeqDense"))
"""
mutable struct MatSeqDense{PetscLib, PetscScalar} <:
               AbstractMat{PetscLib, PetscScalar}
    ptr::CMat
    array::Matrix{PetscScalar}
end

function MatSeqDense(
    petsclib::PetscLib,
    A::Matrix{PetscScalar},
) where {PetscLib <: PetscLibType, PetscScalar}
    comm = MPI.COMM_SELF
    @assert initialized(petsclib)
    @assert PetscScalar == petsclib.PetscScalar
    mat = MatSeqDense{PetscLib, PetscScalar}(C_NULL, A)
    LibPETSc.MatCreateSeqDense(petsclib, comm, size(A, 1), size(A, 2), A, mat)
    finalizer(destroy, mat)
    return mat
end

"""
    MatAIJ(
        petsclib::PetscLib,
        comm::MPI.Comm,
        loc_num_rows::Integer,
        loc_num_cols::Integer,
        diag_nonzeros::Union{Integer, Vector},
        off_diag_nonzeros::Union{Integer, Vector};
        glo_num_rows = PETSC_DETERMINE,
        glo_num_cols = PETSC_DETERMINE,
        setup = true
    ) where {PetscLib <: PetscLibType}

Create an MPI PETSc sparse array on the `comm` using AIJ format (also known as a
compressed sparse row or CSR format) of size `glo_num_rows X glo_num_cols` with
local size `loc_num_rows X loc_num_cols`.

The diagonal block and off-diagonal block non-zeros are `diag_nonzeros` and
`off_diag_nonzeros` which can be either an integer (same for all rows) or a
Vector of `PetscInt`s with on entry per row.

Memory allocation is handled by PETSc and garbage collection can be used.

If `glo_num_rows isa Integer` or `glo_num_cols isa Integer` then the
corresponding local variable can be `PETSC_DECIDE`.

If `setup == true` then [`setup!`](@ref) is called

# External Links
$(_doc_external("Mat/MatCreateAIJ"))
$(_doc_external("Mat/MatSetUp"))

!!! note

    The user is responsible for calling `destroy(mat)` on the `MatAIJ` since
    this cannot be handled by the garbage collector do to the MPI nature of the
    object.
"""
mutable struct MatAIJ{PetscLib, PetscScalar} <:
               AbstractMat{PetscLib, PetscScalar}
    ptr::CMat
end

function MatAIJ(
    petsclib::PetscLib,
    comm::MPI.Comm,
    loc_num_rows::Integer,
    loc_num_cols::Integer,
    diag_nonzeros::Union{Integer, Vector},
    off_diag_nonzeros::Union{Integer, Vector};
    glo_num_rows = PETSC_DETERMINE,
    glo_num_cols = PETSC_DETERMINE,
    setup = true,
) where {PetscLib <: PetscLibType}
    @assert initialized(petsclib)
    PetscScalar = petsclib.PetscScalar
    mat = MatAIJ{PetscLib, PetscScalar}(C_NULL)
    if diag_nonzeros isa Integer
        diag_nonzero = diag_nonzeros
        diag_nonzeros = C_NULL
    else
        diag_nonzero = -1
    end
    if off_diag_nonzeros isa Integer
        off_diag_nonzero = off_diag_nonzeros
        off_diag_nonzeros = C_NULL
    else
        off_diag_nonzero = -1
    end

    LibPETSc.MatCreateAIJ(
        petsclib,
        comm,
        loc_num_rows,
        loc_num_cols,
        glo_num_rows,
        glo_num_cols,
        diag_nonzero,
        diag_nonzeros,
        off_diag_nonzero,
        off_diag_nonzeros,
        mat,
    )

    setup && setup!(mat)

    return mat
end

"""
    setup!(mat::AbstractMat)

Set up the interal data for `mat`

# External Links
$(_doc_external("Mat/MatSetUp"))
"""
function setup!(mat::AbstractMat{PetscLib}) where {PetscLib}
    @assert initialized(PetscLib)
    LibPETSc.MatSetUp(PetscLib, mat)
    return mat
end

"""
    assemble!(M::AbstractMat[, t::MatAssemblyType = MAT_FINAL_ASSEMBLY)

Assemble the matrix `M` with assembly type `t`.

For overlapping assembly see [`assemblybegin!`](@ref) and  [`assemblyend!`](@ref)

# External Links
$(_doc_external("Mat/MatAssemblyBegin"))
$(_doc_external("Mat/MatAssemblyEnd"))
"""
function assemble!(M::AbstractMat, t::MatAssemblyType = MAT_FINAL_ASSEMBLY)
    assemblybegin!(M, t)
    assemblyend!(M, t)
    return M
end

"""
    assemblybegin!(M::AbstractMat[, t::MatAssemblyType = MAT_FINAL_ASSEMBLY)

Begin assembly of the matrix `M` with assembly type `t`; finished with
[`assemblyend!`](@ref).

# External Links
$(_doc_external("Mat/MatAssemblyBegin"))
"""
function assemblybegin!(
    M::AbstractMat{PetscLib},
    t::MatAssemblyType = MAT_FINAL_ASSEMBLY,
) where {PetscLib}
    LibPETSc.MatAssemblyBegin(PetscLib, M, t)
    return M
end

"""
    assemblyend!(M::AbstractMat[, t::MatAssemblyType = MAT_FINAL_ASSEMBLY)

Finish assembly of the matrix `M` with assembly type `t`; start assembly with
[`assemblybegin!`](@ref).

# External Links
$(_doc_external("Mat/MatAssemblyEnd"))
"""
function assemblyend!(
    M::AbstractMat{PetscLib},
    t::MatAssemblyType = MAT_FINAL_ASSEMBLY,
) where {PetscLib}
    LibPETSc.MatAssemblyEnd(PetscLib, M, t)
    return M
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
    mat::AbstractMat{PetscLib},
    base_one::Bool = true,
) where {PetscLib}
    PetscInt = PetscLib.PetscInt
    r_lo = Ref{PetscInt}()
    r_hi = Ref{PetscInt}()
    LibPETSc.MatGetOwnershipRange(PetscLib, mat, r_lo, r_hi)
    return base_one ? ((r_lo[] + PetscInt(1)):(r_hi[])) :
           ((r_lo[]):(r_hi[] - PetscInt(1)))
end

function Base.size(A::AbstractMat{PetscLib}) where {PetscLib}
    m = Ref{PetscLib.PetscInt}()
    n = Ref{PetscLib.PetscInt}()
    LibPETSc.MatGetSize(PetscLib, A, m, n)
    return (m[], n[])
end

function Base.:(==)(
    A::AbstractMat{PetscLib},
    B::AbstractMat{PetscLib},
) where {PetscLib}
    fr = Ref{PetscBool}()
    LibPETSc.MatEqual(PetscLib, A, B, fr)
    return fr[] == PETSC_TRUE
end

function view(
    mat::AbstractMat{PetscLib},
    viewer = LibPETSc.PETSC_VIEWER_STDOUT_(PetscLib, getcomm(mat)),
) where {PetscLib}
    LibPETSc.MatView(PetscLib, mat, viewer)
    return nothing
end
Base.show(io::IO, mat::AbstractMat) = _show(io, mat)
Base.show(io::IO, ::MIME"text/plain", mat::AbstractMat) = _show(io, mat)

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
    M::AbstractMat{PetscLib},
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
        num_rows,
        row0idxs,
        num_cols,
        col0idxs,
        rowvals,
        insertmode,
    )
    return nothing
end

"""
    setvalues!(
        M::AbstractMat{PetscLib},
        row0idxs::Vector{MatStencil{PetscInt}},
        col0idxs::Vector{MatStencil{PetscInt}},
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
    M::AbstractMat{PetscLib},
    row0idxs::Vector{MatStencil{PetscInt}},
    col0idxs::Vector{MatStencil{PetscInt}},
    rowvals::Array{PetscScalar},
    insertmode::InsertMode = INSERT_VALUES;
    num_rows = length(row0idxs),
    num_cols = length(col0idxs),
) where {PetscLib, PetscScalar, PetscInt}
    @assert PetscScalar == PetscLib.PetscScalar
    @assert PetscInt == PetscLib.PetscInt
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
    M::AbstractMat{PetscLib},
    val,
    i::Integer,
    j::Integer,
) where {PetscLib}
    PetscInt = PetscLib.PetscInt
    PetscScalar = PetscLib.PetscScalar
    setvalues!(
        M,
        [PetscInt(i - 1)],
        [PetscInt(j - 1)],
        [PetscScalar(val)],
        INSERT_VALUES,
    )
    return val
end

function Base.setindex!(
    M::AbstractMat{PetscLib},
    val,
    i::CartesianIndex{N},
    j::CartesianIndex{N},
) where {PetscLib, N}
    PetscInt = PetscLib.PetscInt
    PetscScalar = PetscLib.PetscScalar
    ms_i = MatStencil{PetscInt}(
        N < 3 ? 0 : i[3] - 1,
        N < 2 ? 0 : i[2] - 1,
        i[1] - 1,
        N < 4 ? 0 : i[4] - 1,
    )
    ms_j = MatStencil{PetscInt}(
        N < 3 ? 0 : j[3] - 1,
        N < 2 ? 0 : j[2] - 1,
        j[1] - 1,
        N < 4 ? 0 : j[4] - 1,
    )
    setvalues!(M, [ms_i], [ms_j], [PetscScalar(val)], INSERT_VALUES)
    return val
end

function addindex!(
    M::AbstractMat{PetscLib},
    val,
    i::CartesianIndex{N},
    j::CartesianIndex{N},
) where {PetscLib, N}
    PetscInt = PetscLib.PetscInt
    PetscScalar = PetscLib.PetscScalar
    ms_i = MatStencil{PetscInt}(
        N < 3 ? 0 : i[3] - 1,
        N < 2 ? 0 : i[2] - 1,
        i[1] - 1,
        N < 4 ? 0 : i[4] - 1,
    )
    ms_j = MatStencil{PetscInt}(
        N < 3 ? 0 : j[3] - 1,
        N < 2 ? 0 : j[2] - 1,
        j[1] - 1,
        N < 4 ? 0 : j[4] - 1,
    )
    setvalues!(M, [ms_i], [ms_j], [PetscScalar(val)], ADD_VALUES)
    return val
end

"""
    getvalues!(
        rowvals::Array{PetscScalar};
        M::AbstractMat{PetscLib},
        row0idxs::Vector{PetscInt},
        col0idxs::Vector{PetscInt},
        num_rows = length(row0idxs),
        num_cols = length(col0idxs)
    )

get values of the matrix `M` with base-0  row and column indices `row0idxs` and
`col0idxs` inserting the values `rowvals`.

If the keyword arguments `num_rows` or `num_cols` is specified then only the
first `num_rows * num_cols` values of `rowvals` will be used.

# External Links
$(_doc_external("Mat/MatGetValues"))
"""
function getvalues!(
    rowvals::Array{PetscScalar},
    M::AbstractMat{PetscLib},
    row0idxs::Vector{PetscInt},
    col0idxs::Vector{PetscInt};
    num_rows = length(row0idxs),
    num_cols = length(col0idxs),
) where {PetscLib, PetscScalar, PetscInt}
    @assert PetscScalar == PetscLib.PetscScalar
    @assert PetscInt == PetscLib.PetscInt
    @assert num_rows * num_cols <= length(rowvals)
    LibPETSc.MatGetValues(
        PetscLib,
        M,
        num_rows,
        row0idxs,
        num_cols,
        col0idxs,
        rowvals,
    )
    return nothing
end

function Base.getindex(
    M::AbstractMat{PetscLib},
    i::Integer,
    j::Integer,
) where {PetscLib}
    PetscInt = PetscLib.PetscInt
    PetscScalar = PetscLib.PetscScalar
    v_array = [zero(PetscScalar)]
    getvalues!(v_array, M, [PetscInt(i - 1)], [PetscInt(j - 1)])
    return @inbounds v_array[1]
end

function LinearAlgebra.norm(
    M::AbstractMat{PetscLib},
    normtype::NormType = NORM_FROBENIUS,
) where {PetscLib}
    PetscReal = PetscLib.PetscReal
    r_val = Ref{PetscReal}()
    LibPETSc.MatNorm(PetscLib, M, normtype, r_val)
    return r_val[]
end

function LinearAlgebra.mul!(
    y::AbstractVec{PetscLib, PetscScalar},
    M::AbstractMat{PetscLib, PetscScalar},
    x::AbstractVec{PetscLib, PetscScalar},
) where {PetscLib, PetscScalar}
    LibPETSc.MatMult(PetscLib, M, x, y)
    return y
end

function LinearAlgebra.mul!(
    y::AbstractVec{PetscLib, PetscScalar},
    M::Adjoint{PetscScalar, AM},
    x::AbstractVec{PetscLib, PetscScalar},
) where {PetscLib, PetscScalar, AM <: AbstractMat{PetscLib, PetscScalar}}
    LibPETSc.MatMultHermitianTranspose(PetscLib, parent(M), x, y)
    return y
end

function LinearAlgebra.mul!(
    y::AbstractVec{PetscLib, PetscScalar},
    M::Transpose{PetscScalar, AM},
    x::AbstractVec{PetscLib, PetscScalar},
) where {PetscLib, PetscScalar, AM <: AbstractMat{PetscLib, PetscScalar}}
    LibPETSc.MatMultTranspose(PetscLib, parent(M), x, y)
    return y
end

const MatAT{PetscLib, PetscScalar} = Union{
    AbstractMat{PetscLib, PetscScalar},
    Transpose{PetscScalar, <:AbstractMat{PetscLib, PetscScalar}},
    Adjoint{PetscScalar, <:AbstractMat{PetscLib, PetscScalar}},
}

function LinearAlgebra.mul!(
    y::Vector{PetscScalar},
    M::MatAT{PetscLib, PetscScalar},
    x::Vector{PetscScalar},
) where {PetscScalar, PetscLib}
    LinearAlgebra.mul!(
        VecSeqWithArray(PetscLib, y),
        M,
        VecSeqWithArray(PetscLib, x),
    )
    return y
end

function LinearAlgebra.issymmetric(
    A::AbstractMat{PetscLib};
    tol = 0,
) where {PetscLib}
    fr = Ref{PetscBool}()
    LibPETSc.MatIsSymmetric(PetscLib, A, tol, fr)
    return fr[] == PETSC_TRUE
end

function LinearAlgebra.ishermitian(
    A::AbstractMat{PetscLib};
    tol = 0,
) where {PetscLib}
    fr = Ref{PetscBool}()
    LibPETSc.MatIsHermitian(PetscLib, A, tol, fr)
    return fr[] == PETSC_TRUE
end

#=
function MatSeqAIJ(S::SparseMatrixCSC{T}) where {T}
    PetscInt = inttype(T)
    m,n = size(S)
    nnz = zeros(PetscInt,m)
    for r in S.rowval
        nnz[r] += 1
    end
    M = MatSeqAIJ{T}(m, n, nnz)
    for j = 1:n
        for ii = S.colptr[j]:S.colptr[j+1]-1
            i = S.rowval[ii]
            M[i,j] = S.nzval[ii]
        end
    end
    assemble(M)
    return M
end

function Base.copyto!(M::PETSc.MatSeqAIJ{T}, S::SparseMatrixCSC{T}) where {T}
    for j = 1:size(S,2)
        for ii = S.colptr[j]:S.colptr[j+1]-1
            i = S.rowval[ii]
            M[i,j] = S.nzval[ii]
        end
    end
    assemble(M);  
end
=#

"""
    createvecs(
        M::AbstractMat{PetscLib},
    )

Returns vectors `V` which are compatible with `M`. A right compatible vectors is
`V.right` and a left compatible vector is `V.left`; positionally these are
returned as `(right, left)`

The created vectors are not garbage collected and should be destroyed with
[`destroy`](@ref).

# External Links
$(_doc_external("Mat/MatCreateVecs"))
"""
function createvecs(M::AbstractMat{PetscLib}) where {PetscLib}
    r_right = Ref{CVec}()
    r_left = Ref{CVec}()
    LibPETSc.MatCreateVecs(PetscLib, M, r_right, r_left)
    right = VecPtr(PetscLib, r_right[], true)
    left = VecPtr(PetscLib, r_left[], true)
    return (right = right, left = left)
end

"""
    getpetsctype(mat::AbstractMat)

return a string with the matrix type

# External Links
$(_doc_external("Mat/MatGetType"))
"""
function getpetsctype(mat::AbstractMat{PetscLib}) where {PetscLib}
    name_r = Ref{LibPETSc.MatType}()
    LibPETSc.MatGetType(PetscLib, mat, name_r)
    return unsafe_string(name_r[])
end
