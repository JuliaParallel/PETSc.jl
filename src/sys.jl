

"""
    comm = function comm(
                            obj::Union{
                                PetscVec{PetscLib},
                                PetscMat{PetscLib},
                                KSP{PetscLib},
                                #SNES{PetscLib},
                                #PetscDM{PetscLib},
                            },
                        ) where {PetscLib}

Gets the MPI communicator for any of the objects above                         

"""
function comm(
    obj::Union{
        AbstractPetscVec{PetscLib},
        AbstractPetscMat{PetscLib},
        AbstractKSP{PetscLib},
        AbstractSNES{PetscLib},
        AbstractPetscDM{PetscLib},
    },
) where {PetscLib}
    comm = LibPETSc.PetscObjectGetComm(PetscLib, obj)
    return comm
end
