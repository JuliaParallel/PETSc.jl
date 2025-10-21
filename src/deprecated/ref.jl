
@for_libpetsc begin

    function incref(::Type{$PetscScalar}, obj)
        @chk ccall((:PetscObjectReference, $libpetsc), PetscErrorCode, (Ptr{Cvoid},), obj)
        return nothing
    end

    function decref(::Type{$PetscScalar}, obj)
        if !finalized($petsclib)
            @chk ccall((:PetscObjectDereference, $libpetsc), PetscErrorCode, (Ptr{Cvoid},), obj)
            obj.ptr = C_NULL
        end
        return nothing
    end

    function nrefs(::Type{$PetscScalar}, obj)
        r_n = Ref{$PetscInt}()
        @chk ccall((:PetscObjectGetReference, $libpetsc), PetscErrorCode, (Ptr{Cvoid},Ptr{$PetscInt}), obj, r_n)
        return r_n[]
    end

end

"""
    incref(obj)

Increment the reference counter fo `obj`.
This usually only needs to be called when accessing objects owned by other objects, e.g. via `KSPGetPC`.
"""
incref(obj) = incref(scalartype(obj), obj)

"""
    decref(obj)

Decrement the reference counter for `obj`.

In general we don't need to use this, as we can call `destroy` instead.
"""
decref(obj) = decref(scalartype(obj), obj)

"""
    nrefs(obj)

The current reference count for `obj`.
"""
nrefs(obj)  = nrefs(scalartype(obj), obj)
