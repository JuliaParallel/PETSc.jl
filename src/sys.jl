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
