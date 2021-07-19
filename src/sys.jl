const CPetscObject = Ptr{Cvoid}

const UnionPetscTypes = Union{Options}

# allows us to pass PETSc_XXX objects directly into CXXX ccall signatures
Base.cconvert(::Type{CPetscObject}, obj::UnionPetscTypes) = obj
Base.unsafe_convert(::Type{CPetscObject}, obj::UnionPetscTypes) = obj.ptr

# allows us to pass PETSc_XXX objects directly into Ptr{CXXX} ccall signatures
function Base.unsafe_convert(::Type{Ptr{CPetscObject}}, obj::UnionPetscTypes)
    convert(Ptr{CPetscObject}, pointer_from_objref(obj))
end
