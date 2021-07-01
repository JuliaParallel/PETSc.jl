const CMat = Ptr{Cvoid}

abstract type AbstractMat{T} <: AbstractMatrix{T} end
scalartype(::AbstractMat{T}) where {T} = T

# allows us to pass XXMat objects directly into CMat ccall signatures
Base.cconvert(::Type{CMat}, obj::AbstractMat) = obj.ptr
# allows us to pass XXMat objects directly into Ptr{CMat} ccall signatures
Base.unsafe_convert(::Type{Ptr{CMat}}, obj::AbstractMat) =
    convert(Ptr{CMat}, pointer_from_objref(obj))


Base.eltype(::Type{A}) where {A<:AbstractMat{T}} where {T} = T
Base.eltype(A::AbstractMat{T}) where {T} = T

"""
    MatSeqAIJ{T}

PETSc sparse array using AIJ format (also known as a compressed sparse row or CSR format).

Memory allocation is handled by PETSc.
"""
mutable struct MatSeqAIJ{T} <: AbstractMat{T}
    ptr::CMat
    comm::MPI.Comm
end

"""
    MatSeqDense{T}

PETSc dense array. This wraps a Julia `Matrix{T}` object.
"""
mutable struct MatSeqDense{T} <: AbstractMat{T}
    ptr::CMat
    comm::MPI.Comm
    array::Matrix{T}
end



@for_libpetsc begin
    function MatSeqAIJ{$PetscScalar}(m::Integer, n::Integer, nnz::Vector{$PetscInt})
        initialize($PetscScalar)
        comm = MPI.COMM_SELF
        mat = MatSeqAIJ{$PetscScalar}(C_NULL, comm)
        @chk ccall((:MatCreateSeqAIJ, $libpetsc), PetscErrorCode,
                (MPI.MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{CMat}),
                comm, m, n, 0, nnz, mat)
        finalizer(destroy, mat)
        return mat
    end
    function MatSeqDense(A::Matrix{$PetscScalar})
        initialize($PetscScalar)
        comm = MPI.COMM_SELF
        mat = MatSeqDense(C_NULL, comm, A)
        @chk ccall((:MatCreateSeqDense, $libpetsc), PetscErrorCode, 
            (MPI.MPI_Comm, $PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{CMat}),
            comm, size(A,1), size(A,2), A, mat)
        finalizer(destroy, mat)
        return mat
    end



    function destroy(M::AbstractMat{$PetscScalar})
        finalized($PetscScalar) ||
        @chk ccall((:MatDestroy, $libpetsc), PetscErrorCode, (Ptr{CMat},), M)
        return nothing
    end

    function setvalues!(M::AbstractMat{$PetscScalar}, row0idxs::Vector{$PetscInt}, col0idxs::Vector{$PetscInt}, rowvals::Array{$PetscScalar}, insertmode::InsertMode)
        @chk ccall((:MatSetValues, $libpetsc), PetscErrorCode, 
            (CMat, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar},InsertMode),
            M, length(row0idxs), row0idxs, length(col0idxs), col0idxs, rowvals, insertmode)
        return nothing
    end


    function Base.setindex!(M::AbstractMat{$PetscScalar}, val, i::Integer, j::Integer)    
        @chk ccall((:MatSetValues, $libpetsc), PetscErrorCode, 
            (CMat, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}, InsertMode),
            M, 1, Ref{$PetscInt}(i-1), 1, Ref{$PetscInt}(j-1), Ref{$PetscScalar}(val), INSERT_VALUES)
        return val
    end


    function assemblybegin(M::AbstractMat{$PetscScalar}, t::MatAssemblyType=MAT_FINAL_ASSEMBLY)
        @chk ccall((:MatAssemblyBegin, $libpetsc), PetscErrorCode, (CMat, MatAssemblyType), M, t)
        return nothing
    end
    function assemblyend(M::AbstractMat{$PetscScalar}, t::MatAssemblyType=MAT_FINAL_ASSEMBLY)
        @chk ccall((:MatAssemblyEnd, $libpetsc), PetscErrorCode, (CMat, MatAssemblyType), M, t)
        return nothing
    end
    function view(mat::AbstractMat{$PetscScalar}, viewer::Viewer{$PetscScalar}=ViewerStdout{$PetscScalar}(mat.comm))
        @chk ccall((:MatView, $libpetsc), PetscErrorCode, 
                    (CMat, CPetscViewer),
                mat, viewer);
        return nothing
    end


    function Base.size(A::AbstractMat{$PetscScalar})
        m = Ref{$PetscInt}()
        n = Ref{$PetscInt}()
        @chk ccall((:MatGetSize, $libpetsc), PetscErrorCode,
            (CMat, Ptr{$PetscInt}, Ptr{$PetscInt}), 
            A, m, n)
        return (m[], n[])
    end
    function Base.:(==)(A::AbstractMat{$PetscScalar}, B::AbstractMat{$PetscScalar})
        fr = Ref{PetscBool}()
        @chk ccall((:MatEqual, $libpetsc), PetscErrorCode, 
             (CMat, CMat, Ptr{PetscBool}), 
             A, B, fr)
        return fr[] == PETSC_TRUE
    end
    
    function LinearAlgebra.issymmetric(A::AbstractMat{$PetscScalar}; tol=zero($PetscReal))
        fr = Ref{PetscBool}()
        @chk ccall((:MatIsSymmetric, $libpetsc), PetscErrorCode,
            (CMat, $PetscReal, Ptr{PetscBool}),
            A, tol, fr)
        return fr[] == PETSC_TRUE
    end
    function LinearAlgebra.ishermitian(A::AbstractMat{$PetscScalar}; tol=zero($PetscReal))
        fr = Ref{PetscBool}()
        @chk ccall((:MatIsHermitian, $libpetsc), PetscErrorCode,
            (CMat, $PetscReal, Ptr{PetscBool}),
            A, tol, fr)
        return fr[] == PETSC_TRUE
    end
    function LinearAlgebra.norm(M::AbstractMat{$PetscScalar}, normtype::NormType=NORM_FROBENIUS)
        r_val = Ref{$PetscReal}()
        @chk ccall((:MatNorm, $libpetsc), PetscErrorCode, 
            (CMat, NormType, Ptr{$PetscReal}),
            M, normtype, r_val)
        return r_val[]
    end
    
    function LinearAlgebra.mul!(y::AbstractVec{$PetscScalar}, M::AbstractMat{$PetscScalar}, x::AbstractVec{$PetscScalar})
        @chk ccall((:MatMult, $libpetsc), PetscErrorCode, (CMat, CVec, CVec), M, x, y)
        return y
    end
    function LinearAlgebra.mul!(y::AbstractVec{$PetscScalar}, M::Adjoint{T,A}, x::AbstractVec{$PetscScalar}) where {T,A<:AbstractMat{$PetscScalar}}
        @chk ccall((:MatMultHermitianTranspose, $libpetsc), PetscErrorCode, (CMat, CVec, CVec), parent(M), x, y)
        return y
    end
    function LinearAlgebra.mul!(y::AbstractVec{$PetscScalar}, M::Transpose{T,A}, x::AbstractVec{$PetscScalar}) where {T,A<:AbstractMat{$PetscScalar}}
        @chk ccall((:MatMultTranspose, $libpetsc), PetscErrorCode, (CMat, CVec, CVec), parent(M), x, y)
        return y
    end

   

end    

function assemble(M::AbstractMat, t::MatAssemblyType=MAT_FINAL_ASSEMBLY)
    assemblybegin(M, t)
    assemblyend(M, t)
end

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

function Base.show(io::IO, ::MIME"text/plain", mat::AbstractMat)
    _show(io, mat)
end
AbstractMat(A::Matrix) = MatSeqDense(A)
AbstractMat(A::SparseMatrixCSC) = MatSeqAIJ(A)

const MatAT{T} = Union{AbstractMat{T}, Transpose{T, <:AbstractMat{T}}, Adjoint{T, <:AbstractMat{T}}}

LinearAlgebra.mul!(y::AbstractVector{T}, M::MatAT{T}, x::AbstractVector{T}) where {T} =
    parent(LinearAlgebra.mul!(AbstractVec(y), M, AbstractVec(x)))
