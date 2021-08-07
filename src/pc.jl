
const CPC = Ptr{Cvoid}
const CPCType = Cstring

abstract type AbstractPC{T} end

mutable struct PC{T} <: AbstractPC{T}
    ptr::Ptr{Cvoid}
end

scalartype(::AbstractPC{T}) where {T} = T

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

    function destroy(pc::AbstractPC{$PetscScalar})
        if pc.age == getlib(PetscLib).age && !(finalized(PetscLib)) && pc.ptr != C_NULL
            @chk ccall((:PCDestroy, $libpetsc), PetscErrorCode, (Ptr{CPC},), pc)
        end
        pc.ptr = C_NULL
        return nothing
    end

    function settype!(pc::AbstractPC{$PetscScalar}, pctype::String)
        @chk ccall((:PCSetType, $libpetsc), PetscErrorCode, (CPC, Cstring), pc, pctype)
        return nothing
    end

    function setpc!(ksp::KSP{$PetscScalar}, pc::AbstractPC{$PetscScalar})
        @chk ccall((:KSPSetPC, $libpetsc), PetscErrorCode, (CKSP, CPC), ksp, pc)
        return nothing
    end

    function gettype(pc::AbstractPC{$PetscScalar})
        t_r = Ref{CPCType}()
        @chk ccall((:PCGetType, $libpetsc), PetscErrorCode, (CPC, Ptr{CPCType}),  pc, t_r)
        return unsafe_string(t_r[])
    end

    function view(pc::AbstractPC{$PetscScalar}, viewer::AbstractViewer{$PetscLib}=ViewerStdout($petsclib, getcomm(pc)))
        @chk ccall((:PCView, $libpetsc), PetscErrorCode,
                    (CPC, CPetscViewer),
                pc, viewer);
        return nothing
    end

end


Base.show(io::IO, pc::AbstractPC) = _show(io, pc)
