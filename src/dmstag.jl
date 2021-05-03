# Attempt at include dmstag functions
const CDMStag = Ptr{Cvoid}
const CDMStagType = Cstring

mutable struct DMStag{T} <: Factorization{T}
    ptr::CDMStag
    comm::MPI.Comm
end

# allows us to pass XXMat objects directly into CMat ccall signatures
Base.cconvert(::Type{CDMStag}, obj::DMStag) = obj.ptr
# allows us to pass XXMat objects directly into Ptr{CMat} ccall signatures
Base.unsafe_convert(::Type{Ptr{CDMStag}}, obj::DMStag) =
    convert(Ptr{CDMStag}, pointer_from_objref(obj))

#Base.eltype(::DMStag{T}) where {T} = T

@for_libpetsc begin

    function DMStagCreate1d(comm::MPI.Comm, bndx::DMBoundaryType, M, dof0 ,dof1 ,stencilType::DMStagStencilType,stencilWidth ,lx::Vector)
        
        dm = DMStag{$PetscScalar}(C_NULL, comm)

        @chk ccall((:DMStagCreate1d, $libpetsc), PetscErrorCode,
                (MPI.MPI_Comm, DMBoundaryType, $PetscInt, $PetscInt, $PetscInt, DMStagStencilType, $PetscInt,  Ptr{$PetscInt}, Ptr{CDMStag}),
                comm, bndx, M,dof0,dof1,stencilType,stencilWidth,lx, dm )
        
        @chk ccall((:DMSetUp, $libpetsc), PetscErrorCode, (Ptr{CDMStag}, ), dm )

       # if comm == MPI.COMM_SELF
       #     finalizer(destroy, dm)
       # end
        
        return dm
    end


    function DMStagGetGlobalSizes(dm::DMStag)

        M = Ref{$PetscInt}()
        N = Ref{$PetscInt}()
        P = Ref{$PetscInt}()

        @chk ccall((:DMStagGetGlobalSizes, $libpetsc), PetscErrorCode,
            (CDMStag, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}), 
            dm, M, N, P )

        return M[], N[], P[]    
    end








end

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------