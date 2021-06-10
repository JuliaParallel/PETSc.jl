const CPetscObject = Ptr{Cvoid}

const UnionPetscTypes = Union{
    AbstractVec,
    AbstractMat,
    AbstractKSP,
    AbstractViewer,
    AbstractSNES,
    AbstractPC,
    AbstractDM,
    Options,
}

# allows us to pass PETSc_XXX objects directly into CXXX ccall signatures
Base.cconvert(::Type{CPetscObject}, obj::UnionPetscTypes) = obj
Base.unsafe_convert(::Type{CPetscObject}, obj::UnionPetscTypes) = obj.ptr

# allows us to pass PETSc_XXX objects directly into Ptr{CXXX} ccall signatures
function Base.unsafe_convert(::Type{Ptr{CPetscObject}}, obj::UnionPetscTypes)
    convert(Ptr{CPetscObject}, pointer_from_objref(obj))
end

@for_petsc begin
    function getcomm(
        obj::Union{
            AbstractKSP{$PetscScalar},
            AbstractMat{$PetscScalar},
            AbstractVec{$PetscScalar},
            AbstractDM{$PetscScalar},
        },
    )
        comm = MPI.Comm()
        @chk ccall(
            (:PetscObjectGetComm, $petsc_library),
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
            (:PetscCommDuplicate, $petsc_library),
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
@for_petsc begin
    function PetscCommDestroy(
        comm::MPI.Comm
    )
        @show comm.val
        @chk ccall(
            (:PetscCommDestroy, $petsc_library),
            PetscErrorCode,
            (Ptr{MPI.MPI_Comm},),
            comm,
        )
        return nothing
    end
end
=#
