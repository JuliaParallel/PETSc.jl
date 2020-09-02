
@for_libpetsc begin

    function incref(::Type{$PetscScalar}, obj)
        @chk ccall((:PetscObjectReference, $libpetsc), PetscErrorCode, (Ptr{Cvoid},), obj)
    end

    function decref(::Type{$PetscScalar}, obj)
        @chk ccall((:PetscObjectDereference, $libpetsc), PetscErrorCode, (Ptr{Cvoid},), obj)
        if nrefs($PetscScalar, obj) == 0
            obj.ptr = C_NULL
        end
    end

    function nrefs(::Type{$PetscScalar}, obj)
        r_n = Ref{$PetscInt}()
        @chk ccall((:PetscObjectGetReference, $libpetsc), PetscErrorCode, (Ptr{Cvoid},Ptr{$PetscInt}), obj, r_n)
        return r_n[]
    end

end

incref(obj) = incref(scalartype(obj), obj)
