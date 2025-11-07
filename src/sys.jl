

"""
    comm = function getcomm(
                            obj::Union{
                                PetscVec{PetscLib},
                                PetscMat{PetscLib},
                                PetscKSP{PetscLib},
                                #PetscSNES{PetscLib},
                                #PetscDM{PetscLib},
                            },
                        ) where {PetscLib}

Gets the MPI communicator for any of the objects above                         

"""
function getcomm(
    obj::Union{
        PetscVec{PetscLib},
        PetscMat{PetscLib},
        PetscKSP{PetscLib},
        PetscSNES{PetscLib},
        PetscDM{PetscLib},
    },
) where {PetscLib}
    comm = LibPETSc.PetscObjectGetComm(PetscLib, obj)
    return comm
end
