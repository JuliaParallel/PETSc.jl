"""
    comm(obj)

The MPI communicator `obj` lives on: a vector, matrix, DM, `KSP`, `SNES`, `TS`, `PC`
or index set.

# External Links
$(doc_external("Sys/PetscObjectGetComm"))
"""
function comm(
    obj::Union{
        AbstractPetscVec{PetscLib},
        AbstractPetscMat{PetscLib},
        AbstractKSP{PetscLib},
        AbstractSNES{PetscLib},
        LibPETSc.AbstractTS{PetscLib},
        AbstractPetscDM{PetscLib},
        LibPETSc.AbstractPC{PetscLib},
        LibPETSc.AbstractIS{PetscLib},
    },
) where {PetscLib}
    comm = LibPETSc.PetscObjectGetComm(PetscLib, obj)
    return comm
end



"""
    set_name!(obj, name)

Set the name of the PETSc object `obj` (a Vec, Mat, DM, FE, solver, …) to `name`,
and return `obj`.

# External Links
$(doc_external("Sys/PetscObjectSetName"))
"""
function set_name!(
    obj::Union{
        AbstractPetscVec{PetscLib},
        AbstractPetscMat{PetscLib},
        AbstractKSP{PetscLib},
        AbstractSNES{PetscLib},
        LibPETSc.AbstractTS{PetscLib},
        AbstractPetscDM{PetscLib},
        LibPETSc.AbstractPC{PetscLib},
        LibPETSc.AbstractIS{PetscLib},
        LibPETSc.AbstractPetscFE{PetscLib},
        LibPETSc.AbstractDMLabel{PetscLib},
        LibPETSc.AbstractMatNullSpace{PetscLib},
    },
    name::AbstractString,
) where {PetscLib}
    LibPETSc.PetscObjectSetName(PetscLib, convert(Ptr{Cvoid}, obj), String(name))
    return obj
end
