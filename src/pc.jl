
const CPC = Ptr{Cvoid}
const CPCType = Cstring


mutable struct PC{T}
    ptr::Ptr{Cvoid}
end

Base.cconvert(::Type{CPC}, obj::PC) = obj.ptr
Base.unsafe_convert(::Type{Ptr{CPC}}, obj::PC) =
    convert(Ptr{CPC}, pointer_from_objref(obj))

scalartype(::PC{T}) where {T} = T

@for_libpetsc begin

    function PC{$PetscScalar}(comm::MPI.Comm)
        pc = PC{$PetscScalar}(C_NULL)
        @chk ccall((:PCCreate, $libpetsc), PetscErrorCode, (MPI.MPI_Comm, Ptr{CPC}), comm, pc)
        finalizer(destroy, pc)
        return pc
    end

    function PC(ksp::KSP{$PetscScalar})
        pc = PC{$PetscScalar}(C_NULL)
        @chk ccall((:KSPGetPC, $libpetsc), PetscErrorCode, (CKSP, Ptr{CPC}), ksp, pc)
        incref(pc) # need to manually increment the reference counter
        finalizer(destroy, pc)
        return pc
    end

    function destroy(pc::PC{$PetscScalar})
        finalized($petsclib) ||
        @chk ccall((:PCDestroy, $libpetsc), PetscErrorCode, (Ptr{CPC},), pc)
        return nothing
    end

    function settype!(pc::PC{$PetscScalar}, pctype::String)
        @chk ccall((:PCSetType, $libpetsc), PetscErrorCode, (CPC, Cstring), pc, pctype)
        return nothing
    end

    function setpc!(ksp::KSP{$PetscScalar}, pc::PC{$PetscScalar})
        @chk ccall((:KSPSetPC, $libpetsc), PetscErrorCode, (CKSP, CPC), ksp, pc)
        return nothing
    end

    function gettype(pc::PC{$PetscScalar})
        t_r = Ref{CPCType}()
        @chk ccall((:PCGetType, $libpetsc), PetscErrorCode, (CPC, Ptr{CPCType}),  pc, t_r)
        return unsafe_string(t_r[])
    end

    function view(pc::PC{$PetscScalar}, viewer::AbstractViewer{$PetscLib}=ViewerStdout($petsclib, getcomm(pc)))
        @chk ccall((:PCView, $libpetsc), PetscErrorCode,
                    (CPC, CPetscViewer),
                pc, viewer);
        return nothing
    end

end


Base.show(io::IO, pc::PC) = _show(io, pc)
