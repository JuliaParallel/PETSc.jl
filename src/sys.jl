const CPetscObject = Ptr{Cvoid}

const UnionPetscTypes = Union{Options, AbstractVec, AbstractMat, AbstractKSP, AbstractSNES, AbstractDM}

# allows us to pass PETSc_XXX objects directly into CXXX ccall signatures
Base.cconvert(::Type{CPetscObject}, obj::UnionPetscTypes) = obj
Base.unsafe_convert(::Type{CPetscObject}, obj::UnionPetscTypes) = obj.ptr

# allows us to pass PETSc_XXX objects directly into Ptr{CXXX} ccall signatures
function Base.unsafe_convert(::Type{Ptr{CPetscObject}}, obj::UnionPetscTypes)
    convert(Ptr{CPetscObject}, pointer_from_objref(obj))
end

function getcomm(
    obj::Union{
        AbstractVec{PetscLib},
        AbstractMat{PetscLib},
        AbstractKSP{PetscLib},
        AbstractSNES{PetscLib},
        AbstractDM{PetscLib},
    },
) where {PetscLib}
    comm = MPI.Comm()
    LibPETSc.PetscObjectGetComm(PetscLib, obj, comm)

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
