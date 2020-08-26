
const CKSP = Ptr{Cvoid}


mutable struct KSP
    ptr::Ptr{Cvoid}
    comm::MPI.Comm
end

# allows us to pass XXMat objects directly into CMat ccall signatures
function Base.cconvert(::Type{CKSP}, obj::KSP)
  obj.ptr
end

# allows us to pass XXMat objects directly into Ptr{CMat} ccall signatures
function Base.unsafe_convert(::Type{Ptr{CKSP}}, obj::KSP)
    convert(Ptr{CKSP}, pointer_from_objref(obj))
end

function KSP(comm::MPI.Comm)
    ksp = KSP(C_NULL, comm)
    @chk ccall((:KSPCreate, libpetsc), PetscErrorCode, (MPI.MPI_Comm, Ptr{CKSP}), comm, ksp)
    return ksp
end

function destroy(ksp::KSP)
    @chk ccall((:KSPDestroy, libpetsc), PetscErrorCode, (Ptr{CKSP},), ksp)
    return nothing
end

function setoperators!(ksp::KSP, A::AbstractMat, P::AbstractMat)
    @chk ccall((:KSPSetOperators, libpetsc), PetscErrorCode, (CKSP, CMat, CMat), ksp, A, P)
    return nothing
end



const CPC = Ptr{Cvoid}
mutable struct PC
    ptr::CPC
end

Base.cconvert(::Type{CPC}, obj::PC) = obj.ptr
Base.unsafe_convert(::Type{Ptr{CPC}}, obj::PC) = 
    convert(Ptr{CPC}, pointer_from_objref(obj))


function PC(ksp::KSP)
    pc = PC(C_NULL)
    @chk ccall((:KSPGetPC, libpetsc), PetscErrorCode, (CKSP, Ptr{CPC}), ksp, pc)
    return pc
end
function settype!(pc::PC, pctype::String)
    @chk ccall((:PCSetType, libpetsc), PetscErrorCode, (CPC, Cstring), pc, pctype)
    return nothing
end
function settolerances!(ksp::KSP; rtol=PETSC_DEFAULT, atol=PETSC_DEFAULT, dtol=PETSC_DEFAULT, maxits=PETSC_DEFAULT)
    @chk ccall((:KSPSetTolerances, libpetsc), PetscErrorCode, 
                (CKSP, PetscReal, PetscReal, PetscReal, PetscInt),
                ksp, rtol, atol, dtol, maxits)
    return nothing
end

function solve!(x::AbstractVec, ksp::KSP, b::AbstractVec)
    @chk ccall((:KSPSolve, libpetsc), PetscErrorCode, 
      (CKSP, CVec, CVec), ksp, b, x)
end

function solve!(x::AbstractVec, ksp::Transpose{T,K}, b::AbstractVec) where {T,K <: KSP}
  @chk ccall((:KSPSolveTranspose, libpetsc), PetscErrorCode, 
    (CKSP, CVec, CVec), ksp, b, x)
end
