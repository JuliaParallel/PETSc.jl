
const CPC = Ptr{Cvoid}
const CPCType = Cstring


mutable struct PC{T}
    ptr::Ptr{Cvoid}
end

Base.cconvert(::Type{CPC}, obj::PC) = obj.ptr
Base.unsafe_convert(::Type{Ptr{CPC}}, obj::PC) =
    convert(Ptr{CPC}, pointer_from_objref(obj))

@for_libpetsc begin

    function PC(ksp::KSP{$PetscScalar})
        pc = PC{$PetscScalar}(C_NULL)
        @chk ccall((:KSPGetPC, $libpetsc), PetscErrorCode, (CKSP, Ptr{CPC}), ksp, pc)
        return pc
    end

    function destroy(pc::PC{$PetscScalar})
        finalized($PetscScalar) ||
        @chk ccall((:PCDestroy, $libpetsc), PetscErrorCode, (Ptr{CPC},), pc)
        return nothing
    end

    function settype!(pc::PC{$PetscScalar}, pctype::String)
        @chk ccall((:PCSetType, $libpetsc), PetscErrorCode, (CPC, Cstring), pc, pctype)
        return nothing
    end

    function gettype(pc::PC{$PetscScalar})
        t_r = Ref{CPCType}()
        @chk ccall((:PCGetType, $libpetsc), PetscErrorCode, (CPC, Ptr{CPCType}),  pc, t_r)
        return unsafe_string(t_r[])
    end

end