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

        #XXX We should really increase the petsc reference counter.
        #    But for for some reason the PETSc says that this communicator is
        #    unknown
        #=
        # Call the PetscCommDuplicate to increase reference count
        @chk ccall(
            (:PetscCommDuplicate, $libpetsc),
            PetscErrorCode,
            (MPI.MPI_Comm, Ptr{MPI.MPI_Comm}, Ptr{Cvoid}),
            comm,
            comm,
            C_NULL,
        )

        # Register PetscCommDestroy to decriment the reference count
        finalizer(PetscCommDestroy, comm)
        =#

        return comm
    end
end

#=
#XXX Not sure why this doesn't work
@for_libpetsc begin
    function PetscCommDestroy(
        comm::MPI.Comm
    )
        @show comm.val
        @chk ccall(
            (:PetscCommDestroy, $libpetsc),
            PetscErrorCode,
            (Ptr{MPI.MPI_Comm},),
            comm,
        )
        return nothing
    end
end
=#
