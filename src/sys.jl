const CPetscObject = Ptr{Cvoid}

@for_libpetsc begin
    function getcomm(
        obj::Union{
            AbstractKSP{$PetscScalar},
            AbstractMat{$PetscScalar},
            AbstractVec{$PetscScalar},
        },
    )
        comm = MPI.Comm()
        @chk ccall(
            (:PetscObjectGetComm, $libpetsc),
            PetscErrorCode,
            (CPetscObject, Ptr{MPI.MPI_Comm}),
            obj,
            comm,
        )
        return comm
    end
end
