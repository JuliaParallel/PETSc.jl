
const CVec = Ptr{Cvoid}

abstract type AbstractVec # <: AbstractVector{PetscScalar}
end

# allows us to pass XXVec objects directly into CVec ccall signatures
function Base.cconvert(::Type{CVec}, obj::AbstractVec)
    obj.vec
end
# allows us to pass XXVec objects directly into Ptr{CVec} ccall signatures
function Base.unsafe_convert(::Type{Ptr{CVec}}, obj::AbstractVec)
    convert(Ptr{CVec}, pointer_from_objref(obj))
end

mutable struct SeqVec <: AbstractVec
    vec::CVec
    comm::MPI.Comm
    array::Vector{PetscScalar}
end

function SeqVec(comm::MPI.Comm, X::Vector{PetscScalar}; blocksize=1)
    v = SeqVec(C_NULL, comm, X)
    @chk ccall((:VecCreateSeqWithArray, libpetsc), PetscErrorCode,
               (MPI.MPI_Comm, PetscInt, PetscInt, Ptr{PetscScalar}, Ptr{CVec}),
               comm, blocksize, length(X), X, v)
    return v
end

    
function destroy(v::AbstractVec)
    @chk ccall((:VecDestroy, libpetsc), PetscErrorCode, (Ptr{CVec},), v)
    return nothing
end



function LinearAlgebra.norm(v::AbstractVec, normtype::NormType=NORM_2)
    r_val = Ref{PetscReal}()
    @chk ccall((:VecNorm, libpetsc), PetscErrorCode,
               (CVec, NormType, Ptr{PetscReal}),
               v, normtype,r_val)
    return r_val[]
end
