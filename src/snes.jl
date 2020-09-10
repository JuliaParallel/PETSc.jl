
const CSNES = Ptr{Cvoid}
const CKSPType = Cstring


mutable struct SNES{T}
    ptr::CSNES
    comm::MPI.Comm
    fn
    fn_vec
    jac
    jac_A
    jac_P
end
scalartype(::SNES{T}) where {T} = T

Base.cconvert(::Type{CSNES}, obj::SNES) = obj.ptr
Base.unsafe_convert(::Type{Ptr{CSNES}}, obj::SNES) =
    convert(Ptr{CSNES}, pointer_from_objref(obj))



# How to handle Jacobians?
#  - https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESComputeJacobianDefault.html
#  - https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESComputeJacobianDefaultColor.html
#  - 

struct SNESFn{T}
end

struct SNESJac{T}
end

@for_libpetsc begin
    function SNES{$PetscScalar}(comm::MPI.Comm, vec, obj)
        snes = SNES{$PetscScalar}(C_NULL, comm, vec, obj)
        @chk ccall((:SNESCreate, $libpetsc), PetscErrorCode, (MPI.MPI_Comm, Ptr{CSNES}), comm, snes)
        
 
        return snes
    end


    function (::SNESFn{$PetscScalar})(csnes::CSNES, cx::CVec, cfx::CVec, ctx::Ptr{Cvoid})::$PetscInt
      snes = unsafe_pointer_to_objref(ctx)
      x = unsafe_localarray(::Type{$PetscScalar}, cx::CVec)
      fx = unsafe_localarray(::Type{$PetscScalar}, cfx::CVec)
      snes.fn(fx, x)
      Base.finalize(x)
      Base.finalize(fx)
      return $PetscInt(0)
    end

    function setfunction!(snes::SNES{$PetscScalar}, fn, vec::AbstractVec{$PetscScalar})
      ctx = pointer_from_objref(snes)
      fptr = @cfunction(SNESFn{$PetscScalar}(), $PetscInt, (CSNES, CVec, CVec, Ptr{Cvoid}))
      @chk ccall((:SNESSetFunction, $libpetsc), PetscErrorCode,
         (CSNES, CVec, Ptr{Cvoid}, Ptr{Cvoid}),
         snes, vec, fptr, ctx)
      snes.fn_vec = vec
      snes.fn = fn
      return nothing
    end



    function (::SNESJac{$PetscScalar})(csnes::CSNES, cx::CVec, cA::CMat, cP::CMat, ctx::Ptr{Cvoid})::$PetscInt
      snes = unsafe_pointer_to_objref(ctx)
      x = unsafe_localarray(::Type{$PetscScalar}, cx::CVec)      
      snes.jac(snes.A, snes.P, x)
      Base.finalize(x)
      return $PetscInt(0)
    end

    function setjacobian!(snes::SNES{$PetscScalar}, jacfn, A::AbstactMat{$PetscScalar}, P::AbstractMat{$PetscScalar})
      ctx = pointer_from_objref(snes)
      jacptr = @cfunction(SNESJac{$PetscScalar}(), $PetscInt, (CSNES, CVec, CMat, CMat, Ptr{Cvoid}))

      @chk ccall((:SNESSetJacobian, $libpetsc), PetscErrorCode,
        (CSNES, CMat, CMat, Ptr{Cvoid}, Ptr{Cvoid}),
        snes, A, P, jacptr, ctx)
    end
end