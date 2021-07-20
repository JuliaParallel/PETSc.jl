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
    MatSeqAIJ{PetscLib, PetscScalar}

PETSc sparse array using AIJ format (also known as a compressed sparse row or
CSR format).

Memory allocation is handled by PETSc.

# External Links
$(_doc_external("Mat/MatCreateSeqAIJ"))
"""
mutable struct MatSeqAIJ{PetscLib, PetscScalar} <:
               AbstractMat{PetscLib, PetscScalar}
    ptr::CMat
end

"""
    MatSeqAIJ(petsclib, num_rows, num_cols, nonzeros)

Create a PETSc sparse array using AIJ format (also known as a compressed sparse
row or CSR format) of size `num_rows X num_cols` with `nonzeros` per row

If `nonzeros` is an `Integer` the same number of non-zeros will be used for each
row, if `nonzeros` is a `Vector{PetscInt}` then one value must be specified for
each row.

Memory allocation is handled by PETSc and garbage collection can be used.

# External Links
$(_doc_external("Mat/MatCreateSeqAIJ"))
"""
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
    assemble(M::AbstractMat[, t::MatAssemblyType = MAT_FINAL_ASSEMBLY)

Assemble the matrix `M` with assembly type `t`.

For overlapping assembly see [`assemblybegin`](@ref) and  [`assemblyend`](@ref)

# External Links
$(_doc_external("Mat/MatAssemblyBegin"))
$(_doc_external("Mat/MatAssemblyEnd"))
"""
function assemble(M::AbstractMat, t::MatAssemblyType = MAT_FINAL_ASSEMBLY)
    assemblybegin(M, t)
    assemblyend(M, t)
end

"""
    assemblybegin(M::AbstractMat[, t::MatAssemblyType = MAT_FINAL_ASSEMBLY)

Begin assembly of the matrix `M` with assembly type `t`; finished with
[`assemblyend`](@ref).

# External Links
$(_doc_external("Mat/MatAssemblyBegin"))
"""
function assemblybegin(
    M::AbstractMat{PetscLib},
    t::MatAssemblyType = MAT_FINAL_ASSEMBLY,
) where {PetscLib}
    LibPETSc.MatAssemblyBegin(PetscLib, M, t)
    return nothing
end

"""
    assemblyend(M::AbstractMat[, t::MatAssemblyType = MAT_FINAL_ASSEMBLY)

Finish assembly of the matrix `M` with assembly type `t`; start assembly with
[`assemblybegin`](@ref).

# External Links
$(_doc_external("Mat/MatAssemblyEnd"))
"""
function assemblyend(
    M::AbstractMat{PetscLib},
    t::MatAssemblyType = MAT_FINAL_ASSEMBLY,
) where {PetscLib}
    LibPETSc.MatAssemblyEnd(PetscLib, M, t)
    return nothing
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
        insertmode::InsertMode;
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
    insertmode::InsertMode;
    num_rows = length(row0idxs),
    num_cols = length(col0idxs),
) where {PetscLib, PetscScalar, PetscInt}
    @assert PetscScalar == PetscLib.PetscScalar
    @assert PetscInt == PetscLib.PetscInt
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
