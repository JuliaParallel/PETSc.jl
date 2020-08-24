const CMat = Ptr{Cvoid}

abstract type AbstractMat <: AbstractMatrix{PetscScalar} end

# allows us to pass XXMat objects directly into CMat ccall signatures
function Base.cconvert(::Type{CMat}, obj::AbstractMat)
    obj.mat
end

# allows us to pass XXMat objects directly into Ptr{CMat} ccall signatures
function Base.unsafe_convert(::Type{Ptr{CMat}}, obj::AbstractMat)
    convert(Ptr{CMat}, pointer_from_objref(obj))
end

mutable struct SeqAIJMat <: AbstractMat
    mat::CMat
    comm::MPI.Comm
end

function SeqAIJMat(comm::MPI.Comm, m::Integer, n::Integer, nnz::Vector{PetscInt})
    mat = SeqAIJMat(C_NULL, comm)
    @chk ccall((:MatCreateSeqAIJ, libpetsc), PetscErrorCode,
               (MPI.MPI_Comm, PetscInt, PetscInt, PetscInt, Ptr{PetscInt}, Ptr{CMat}),
               comm, m, n, 0, nnz, mat)
    return mat
end

function destroy(M::AbstractMat)
    @chk ccall((:MatDestroy, libpetsc), PetscErrorCode, (Ptr{CMat},), M)
    return nothing
end


function Base.setindex!(M::AbstractMat, val, i::Integer, j::Integer)    
    @chk ccall((:MatSetValues, libpetsc), PetscErrorCode, 
        (CMat, PetscInt, Ptr{PetscInt}, PetscInt, Ptr{PetscInt}, Ptr{PetscScalar}, InsertMode),
        M, 1, Ref{PetscInt}(i-1), 1, Ref{PetscInt}(j-1), Ref{PetscScalar}(val), INSERT_VALUES)
    return val
end


function assemblybegin(M::AbstractMat, t::MatAssemblyType=MAT_FINAL_ASSEMBLY)
    @chk ccall((:MatAssemblyBegin, libpetsc), PetscErrorCode, (CMat, MatAssemblyType), M, t)
    return nothing
end

function assemblyend(M::AbstractMat, t::MatAssemblyType=MAT_FINAL_ASSEMBLY)
    @chk ccall((:MatAssemblyEnd, libpetsc), PetscErrorCode, (CMat, MatAssemblyType), M, t)
    return nothing
end

function SeqAIJMat(S::SparseMatrixCSC)
    m,n = size(S)
    nnz = zeros(PetscInt,m)
    for r in S.rowval
        nnz[r] += 1
    end
    M = SeqAIJMat(MPI.COMM_SELF, m, n, nnz)
    for j = 1:n
        for ii = S.colptr[j]:S.colptr[j+1]-1
            i = S.rowval[ii]
            M[i,j] = S.nzval[ii]
        end
    end
    assemblybegin(M)
    assemblyend(M)
    return M
end

function LinearAlgebra.norm(M::AbstractMat, normtype=NORM_FROBENIUS)
    r_val = Ref{PetscReal}()
    @chk ccall((:MatNorm, libpetsc), PetscErrorCode, 
        (CMat, NormType, Ptr{PetscReal}),
        M, normtype, r_val)
    return r_val[]
end

function LinearAlgebra.mul!(y::AbstractVec, M::AbstractMat, x::AbstractVec)
    @chk ccall((:MatMult, libpetsc), PetscErrorCode, (CMat, CVec, CVec), M, x, y)
    return y
end

function LinearAlgebra.mul!(y::AbstractVec, M::Adjoint{T,A}, x::AbstractVec) where {T,A<:AbstractMat}
    @chk ccall((:MatMultHermitianTranspose, libpetsc), PetscErrorCode, (CMat, CVec, CVec), parent(M), x, y)
    return y
end

function LinearAlgebra.mul!(y::AbstractVec, M::Transpose{T,A}, x::AbstractVec) where {T,A<:AbstractMat}
    @chk ccall((:MatMultTranspose, libpetsc), PetscErrorCode, (CMat, CVec, CVec), parent(M), x, y)
    return y
end
