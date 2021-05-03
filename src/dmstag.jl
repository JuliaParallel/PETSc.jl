# Attempt at include dmstag functions
const CDMStag = Ptr{Cvoid}
const CDMStagType = Cstring

mutable struct DMStag{T} <: Factorization{T}
    ptr::CDMStag
    comm::MPI.Comm
end

# allows us to pass XXMat objects directly into CMat ccall signatures
Base.cconvert(::Type{CDMStag}, obj::DMStag) = obj.ptr
# allows us to pass XXMat objects directly into Ptr{CMat} ccall signatures
Base.unsafe_convert(::Type{Ptr{CDMStag}}, obj::DMStag) =
    convert(Ptr{CDMStag}, pointer_from_objref(obj))

Base.eltype(::DMStag{T}) where {T} = T

@for_libpetsc begin

    function DMStagCreate1d{$PetscScalar}(comm::MPI.Comm, bndx::DMBoundaryType, M::Int32,dof0::Int32,dof1::Int32,stencilType::DMStagStencilType,stencilWidth::Int32,lx::Vector{Int32})
        dm = DMStag{$PetscScalar}(C_NULL, comm)
        #@chk ccall((:DMStagCreate1d, $libpetsc), PetscErrorCode,
        #        (MPI.MPI_Comm, DMBoundaryType, $PetscInt, $PetscInt, $PetscInt, DMStagStencilType, $PetscInt, $PetscInt, Ptr{CDMStag}),
        #        comm, bndx, M,dof0,dof1,stencilType,stencilWidth,lx,dm)
        return dm
    end



end

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------