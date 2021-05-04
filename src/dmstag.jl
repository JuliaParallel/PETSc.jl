# Attempt at include dmstag functions
const CDMStag = Ptr{Cvoid}
const CDMStagType = Cstring

mutable struct DMStag{T} <: Factorization{T}
    ptr::CDMStag
    comm::MPI.Comm
    dim::Int64
end

# allows us to pass XXMat objects directly into CMat ccall signatures
Base.cconvert(::Type{CDMStag}, obj::DMStag) = obj.ptr
# allows us to pass XXMat objects directly into Ptr{CMat} ccall signatures
Base.unsafe_convert(::Type{Ptr{CDMStag}}, obj::DMStag) =
    convert(Ptr{CDMStag}, pointer_from_objref(obj))

#Base.eltype(::DMStag{T}) where {T} = T

@for_libpetsc begin

    function DMStagCreate1d(comm::MPI.Comm, bndx::DMBoundaryType, M, dof0 ,dof1 ,stencilType::DMStagStencilType,stencilWidth ,lx::Vector)
        
        dm = DMStag{$PetscScalar}(C_NULL, comm, 1)

        @chk ccall((:DMStagCreate1d, $libpetsc), PetscErrorCode,
                (MPI.MPI_Comm, DMBoundaryType, $PetscInt, $PetscInt, $PetscInt, DMStagStencilType, $PetscInt,  Ptr{$PetscInt}, Ptr{CDMStag}),
                comm, bndx, M,dof0,dof1,stencilType,stencilWidth,lx, dm )
        
        @chk ccall((:DMSetUp, $libpetsc), PetscErrorCode, (Ptr{CDMStag}, ), dm )

        if comm == MPI.COMM_SELF
            finalizer(destroy, dm)
        end
        
        return dm
    end

    """
        Gets the global size of the DMStag object
            M,N,P = DMStagGetGlobalSizes(dm::DMStag)
    """
    function DMStagGetGlobalSizes(dm::DMStag)

        M = Ref{$PetscInt}()
        N = Ref{$PetscInt}()
        P = Ref{$PetscInt}()

        @chk ccall((:DMStagGetGlobalSizes, $libpetsc), PetscErrorCode,
            (CDMStag, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}), 
            dm, M, N, P )

        return M[], N[], P[]    
    end


    """
        Destroys the DMStag object
    """
    function destroy(dm::DMStag{$PetscScalar})
        finalized($PetscScalar) ||
            @chk ccall((:DMDestroy, $libpetsc), PetscErrorCode, (Ptr{CDMStag},), dm)
        return nothing
    end

    """
        Retrieves the Type of the DMStag object
    """
    function gettype(dm::DMStag{$PetscScalar})
        t_r = Ref{CDMStagType}()
        @chk ccall((:DMGetType, $libpetsc), PetscErrorCode, (CDMStag, Ptr{CDMStagType}), dm, t_r)
        return unsafe_string(t_r[])
    end

    function view(dm::DMStag{$PetscScalar}, viewer::Viewer{$PetscScalar}=ViewerStdout{$PetscScalar}(dm.comm))
        @chk ccall((:DMView, $libpetsc), PetscErrorCode, 
                    (CDMStag, CPetscViewer),
                dm, viewer);
        return nothing
    end

    """ 
        Gets the corners of the DMStag grid
            x,m,nExtrax = DMStagGetCorners(dm:DMStag)   in 1D

    """
    function  DMStagGetCorners(dm::DMStag)

        x = Ref{$PetscInt}()
        y = Ref{$PetscInt}()
        z = Ref{$PetscInt}()
        m = Ref{$PetscInt}()
        n = Ref{$PetscInt}()
        p = Ref{$PetscInt}()
        nExtrax = Ref{$PetscInt}()
        nExtray = Ref{$PetscInt}()
        nExtraz = Ref{$PetscInt}()
        
        @chk ccall((:DMStagGetCorners, $libpetsc), PetscErrorCode,
            (CDMStag,   Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt},
                        Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt},
                        Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}), 
            dm, x,y,z, m,n,p, nExtrax,nExtray,nExtraz )

            if dm.dim==1
                return x[], m[], nExtrax[]    
            elseif dm.dim==2
                return x[], y[], m[],n[], nExtrax[],nExtray[]    
            elseif dm.dim==3
                return x[], y[], z[], m[],n[],p[], nExtrax[],nExtray[],nExtraz[]    
            end
    end
   
    """
        returns the dimension of the DMStag object
    """
    function DMGetDimension(dm::DMStag)
        dim = Ref{$PetscInt}()

        @chk ccall((:DMStagGetCorners, $libpetsc), PetscErrorCode,
                    (CDMStag,Ptr{$PetscInt}), dm, dim )

        return dim[]
    end


   # DMStagGetCorners(DM dm,PetscInt *x,PetscInt *y,PetscInt *z,PetscInt *m,PetscInt *n,PetscInt *p,PetscInt *nExtrax,PetscInt *nExtray,PetscInt *nExtraz)


end

Base.show(io::IO, dm::DMStag) = _show(io, dm)
