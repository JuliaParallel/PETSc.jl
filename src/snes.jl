
const CSNES = Ptr{Cvoid}
const CKSPType = Cstring


mutable struct SNES{T,V<:AbstractVec,A}
    ptr::CSNES
    comm::MPI.Comm
    vec::V
    obj::A
end
scalartype(::SNES{T}) where {T} = T

# allows us to pass XXMat objects directly into CMat ccall signatures
Base.cconvert(::Type{CSNES}, obj::SNES) = obj.ptr
# allows us to pass XXMat objects directly into Ptr{CMat} ccall signatures
Base.unsafe_convert(::Type{Ptr{CSNES}}, obj::SNES) =
    convert(Ptr{CSNES}, pointer_from_objref(obj))


function _snesfn(csnes::CSNES, cx::CVec, cfx::CVec, ctx::Ptr{Cvoid})
    snes = unsafe_pointer_to_objref(ctx)
    snes.obj(fx, x)
end


@for_libpetsc begin
    function SNES{$PetscScalar}(comm::MPI.Comm, vec, obj)
        snes = SNES{$PetscScalar}(C_NULL, comm, vec, obj)
        @chk ccall((:SNESCreate, $libpetsc), PetscErrorCode, (MPI.MPI_Comm, Ptr{CSNES}), comm, snes)

        ctx = pointer_from_objref(snes)

        @chk ccall((:SNESSetFunction, $libpetsc), PetscErrorCode,
           (CSNES, CVec, Ptr{Cvoid}, Ptr{Cvoid}),
           snes, vec, fptr, ctx)
        
        # jacptr
        # @chk ccall((:SNESSetJacobian, $libpetsc), PetscErrorCode,
        #    (CSNES, CMat, CMat, Ptr{Cvoid}, Ptr{Cvoid}),
        #    snes, Amat, Pmat, jacptr, ctx)

        return snes
    end

end