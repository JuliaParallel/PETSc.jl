
const CSNES = Ptr{Cvoid}
const CSNESType = Cstring


mutable struct SNES{T}
    ptr::CSNES
    comm::MPI.Comm
    Feval    # user-defined residual evaluation routine
    Jeval    # user-defined jacobian evaluation routine
    opts::Options{T}
end

scalartype(::SNES{T}) where {T} = T

# allows us to pass XXMat objects directly into CMat ccall signatures
Base.cconvert(::Type{CSNES}, obj::SNES) = obj.ptr
# allows us to pass XXMat objects directly into Ptr{CMat} ccall signatures
Base.unsafe_convert(::Type{Ptr{CSNES}}, obj::SNES) =
    convert(Ptr{CSNES}, pointer_from_objref(obj))

Base.eltype(::SNES{T}) where {T} = T

function _snesfn(csnes::CSNES, cx::CVec, cfx::CVec, ctx::Ptr{Cvoid})
    snes = unsafe_pointer_to_objref(ctx)
    snes.Feval(cfx, cx)
end

function _snesjac(csnes::CSNES, cx::CVec, cAmat::CMat, cPmat::CMat, ctx::Ptr{Cvoid})
    snes = unsafe_pointer_to_objref(ctx)
    snes.Jeval(cAmat, cPmat, cx)
end

@for_libpetsc begin

    function SNES{$PetscScalar}(comm::MPI.Comm; kwargs...)
        initialize($PetscScalar)
        opts = Options{$PetscScalar}(kwargs...)
        snes = SNES{$PetscScalar}(C_NULL, comm, nothing, nothing, opts)
        @chk ccall((:SNESCreate, $libpetsc), PetscErrorCode, (MPI.MPI_Comm, Ptr{CSNES}), comm, snes)
        if comm == MPI.COMM_SELF
            finalizer(destroy, snes)
        end
        return snes
    end

    function destroy(snes::SNES{$PetscScalar})
        finalized($PetscScalar) ||
        @chk ccall((:SNESDestroy, $libpetsc), PetscErrorCode, (Ptr{CSNES},), snes)
        return nothing
    end

    function setfromoptions!(snes::SNES{$PetscScalar})
        @chk ccall((:SNESSetFromOptions, $libpetsc), PetscErrorCode, (CSNES,), snes)
    end

    function gettype(snes::SNES{$PetscScalar})
        t_r = Ref{CSNESType}()
        @chk ccall((:SNESGetType, $libpetsc), PetscErrorCode, (CSNES, Ptr{CSNESType}), snes, t_r)
        return unsafe_string(t_r[])
    end

    function view(snes::SNES{$PetscScalar}, viewer::Viewer{$PetscScalar}=ViewerStdout{$PetscScalar}(snes.comm))
        @chk ccall((:SNESView, $libpetsc), PetscErrorCode,
                    (CSNES, CPetscViewer),
                snes, viewer);
        return nothing
    end

end