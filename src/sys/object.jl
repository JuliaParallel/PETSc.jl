"""
    comm(obj)

The MPI communicator `obj` lives on: a vector, matrix, DM, `KSP`, `SNES`, `PC` or
index set (`TS` has a method of its own).

# External Links
$(doc_external("Sys/PetscObjectGetComm"))
"""
function comm(
    obj::Union{
        AbstractPetscVec{PetscLib},
        AbstractPetscMat{PetscLib},
        AbstractKSP{PetscLib},
        AbstractSNES{PetscLib},
        AbstractPetscDM{PetscLib},
        LibPETSc.AbstractPC{PetscLib},
        LibPETSc.AbstractIS{PetscLib},
    },
) where {PetscLib}
    comm = LibPETSc.PetscObjectGetComm(PetscLib, obj)
    return comm
end



"""
    set_name!(petsclib, obj, name)

Set the name of any PETSc object (DM, Vec, FE, …) to `name`.
`obj` can be any pointer type that is convertible to `Ptr{Cvoid}`.

Thin wrapper around `PetscObjectSetName` that avoids needing an explicit
`convert(Ptr{Cvoid}, obj)` at the call site.

# External Links
$(doc_external("Sys/PetscObjectSetName"))
"""
function set_name!(petsclib::LibPETSc.PetscLibType, obj, name::AbstractString)
    LibPETSc.PetscObjectSetName(petsclib, convert(Ptr{Cvoid}, obj), String(name))
    return obj
end
