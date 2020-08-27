
const CKSP = Ptr{Cvoid}

mutable struct KSP{T}
    ptr::Ptr{Cvoid}
    comm::MPI.Comm
end

# allows us to pass XXMat objects directly into CMat ccall signatures
Base.cconvert(::Type{CKSP}, obj::KSP) = obj.ptr
# allows us to pass XXMat objects directly into Ptr{CMat} ccall signatures
Base.unsafe_convert(::Type{Ptr{CKSP}}, obj::KSP) =
    convert(Ptr{CKSP}, pointer_from_objref(obj))


const CPC = Ptr{Cvoid}
mutable struct PC{T}
    ptr::CPC
end

Base.cconvert(::Type{CPC}, obj::PC) = obj.ptr
Base.unsafe_convert(::Type{Ptr{CPC}}, obj::PC) = 
    convert(Ptr{CPC}, pointer_from_objref(obj))

@for_libpetsc begin
    function KSP{$PetscScalar}(comm::MPI.Comm)
        ksp = KSP{$PetscScalar}(C_NULL, comm)
        @chk ccall((:KSPCreate, $libpetsc), PetscErrorCode, (MPI.MPI_Comm, Ptr{CKSP}), comm, ksp)
        return ksp
    end
    function destroy(ksp::KSP{$PetscScalar})
        @chk ccall((:KSPDestroy, $libpetsc), PetscErrorCode, (Ptr{CKSP},), ksp)
        return nothing
    end
    function setoperators!(ksp::KSP{$PetscScalar}, A::AbstractMat{$PetscScalar}, P::AbstractMat{$PetscScalar})
        @chk ccall((:KSPSetOperators, $libpetsc), PetscErrorCode, (CKSP, CMat, CMat), ksp, A, P)
        return nothing
    end

    function PC(ksp::KSP{$PetscScalar})
        pc = PC{$PetscScalar}(C_NULL)
        @chk ccall((:KSPGetPC, $libpetsc), PetscErrorCode, (CKSP, Ptr{CPC}), ksp, pc)
        return pc
    end
    function settype!(pc::PC{$PetscScalar}, pctype::String)
        @chk ccall((:PCSetType, $libpetsc), PetscErrorCode, (CPC, Cstring), pc, pctype)
        return nothing
    end
    function settolerances!(ksp::KSP{$PetscScalar}; rtol=PETSC_DEFAULT, atol=PETSC_DEFAULT, dtol=PETSC_DEFAULT, maxits=PETSC_DEFAULT)
        @chk ccall((:KSPSetTolerances, $libpetsc), PetscErrorCode, 
                    (CKSP, $PetscReal, $PetscReal, $PetscReal, $PetscInt),
                    ksp, rtol, atol, dtol, maxits)
        return nothing
    end

    function iters(ksp::KSP{$PetscScalar})
        r_its = Ref{$PetscInt}()
        @chk ccall((:KSPGetIterationNumber, $libpetsc), PetscErrorCode, 
        (KSP, Ptr{$PetscInt}), ksp, r_its)
        return r_its[]
    end

    function resnorm(ksp::KSP{$PetscScalar})
        r_rnorm = Ref{$PetscReal}()
        @chk ccall((:KSPGetResidualNorm, $libpetsc), PetscErrorCode, 
        (KSP, Ptr{$PetscReal}), ksp, r_rnorm)
        return r_rnorm[]
    end
    function solve!(x::AbstractVec{$PetscScalar}, ksp::KSP{$PetscScalar}, b::AbstractVec{$PetscScalar})
        @chk ccall((:KSPSolve, $libpetsc), PetscErrorCode, 
        (CKSP, CVec, CVec), ksp, b, x)
    end
    function solve!(x::AbstractVec{$PetscScalar}, ksp::Transpose{T,K}, b::AbstractVec{$PetscScalar}) where {T,K <: KSP{$PetscScalar}}
        @chk ccall((:KSPSolveTranspose, $libpetsc), PetscErrorCode, 
        (CKSP, CVec, CVec), ksp, b, x)
    end
end

function KSP(A::AbstractMat{T}, P::AbstractMat{T}=A) where {T}
    ksp = KSP{T}(A.comm)
    setoperators!(ksp, A, P)
    return ksp
end

"""
    iters(ksp::KSP)

Gets the current iteration number; if the `solve!` is complete, returns the number of iterations used.

https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPGetIterationNumber.html
"""
iters


"""
    resnorm(ksp::KSP)

Gets the last (approximate preconditioned) residual norm that has been computed.

https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPGetResidualNorm.html
"""
resnorm

