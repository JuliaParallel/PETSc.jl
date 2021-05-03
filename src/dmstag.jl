# Attempt at include dmstag functions

mutable struct DM{T}
    ptr::CDM
    comm::MPI.Comm
end


@for_libpetsc begin

    function DMStagCreate1d(comm::MPI.Comm, X::Vector{$PetscScalar}; blocksize=1)
        v = VecSeq(C_NULL, comm, X)
        @chk ccall((:DMStagCreate1d, $libpetsc), PetscErrorCode,
                (MPI.MPI_Comm,DMBoundaryType,$PetscInt,$PetscInt,$PetscInt,DMStagStencilType,$PetscInt,const PetscInt,ptr{CDM}),
                comm, bndx, M,dof0,dof1,stencilType,stencilWidth,lx[],dm)
        return nothing
    end



end

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------





# AbstractVec
#   - VecSeq: wrap
#   - VecMPI (TODO)
#   - VecGhost (TODO)
# for the MPI variants we won't be able to attach finalizers, as destroy needs to be called collectively.

const CVec = Ptr{Cvoid}

abstract type AbstractVec{T} <: AbstractVector{T} end
scalartype(::AbstractVec{T}) where {T} = T

# allows us to pass XXVec objects directly into CVec ccall signatures
Base.cconvert(::Type{CVec}, obj::AbstractVec) = obj.ptr
# allows us to pass XXVec objects directly into Ptr{CVec} ccall signatures
Base.unsafe_convert(::Type{Ptr{CVec}}, obj::AbstractVec) =
    convert(Ptr{CVec}, pointer_from_objref(obj))

"""
    VecSeq(v::Vector)

A standard, sequentially-stored serial PETSc vector, wrapping the Julia vector `v`.

This reuses the array `v` as storage, and so `v` should not be `resize!`-ed or otherwise have its length modified while the PETSc object exists.

This should only be need to be called for more advanced uses, for most simple usecases, users should be able to pass `Vector`s directly and have the wrapping performed automatically

https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Vec/VecCreateSeqWithArray.html
"""
mutable struct VecSeq{T} <: AbstractVec{T}
    ptr::CVec
    comm::MPI.Comm
    array::Vector{T}
end

Base.eltype(::Type{V}) where {V<:AbstractVec{T}} where T = T
Base.eltype(v::AbstractVec{T}) where {T} = T
Base.size(v::AbstractVec) = (length(v),)
Base.parent(v::AbstractVec) = v.array

@for_libpetsc begin
    function VecSeq(comm::MPI.Comm, X::Vector{$PetscScalar}; blocksize=1)
        initialize($PetscScalar)
        v = VecSeq(C_NULL, comm, X)
        @chk ccall((:VecCreateSeqWithArray, $libpetsc), PetscErrorCode,
                (MPI.MPI_Comm, $PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{CVec}),
                comm, blocksize, length(X), X, v)
        finalizer(destroy, v)
        return v
    end
    function destroy(v::AbstractVec{$PetscScalar})
        finalized($PetscScalar) ||
        @chk ccall((:VecDestroy, $libpetsc), PetscErrorCode, (Ptr{CVec},), v)
        return nothing
    end
    function Base.length(v::AbstractVec{$PetscScalar})
        r_sz = Ref{$PetscInt}()
        @chk ccall((:VecGetSize, $libpetsc), PetscErrorCode,
          (CVec, Ptr{$PetscInt}), v, r_sz)
        return r_sz[]
    end
    function LinearAlgebra.norm(v::AbstractVec{$PetscScalar}, normtype::NormType=NORM_2)
        r_val = Ref{$PetscReal}()
        @chk ccall((:VecNorm, $libpetsc), PetscErrorCode,
                   (CVec, NormType, Ptr{$PetscReal}),
                   v, normtype,r_val)
        return r_val[]
    end

    function assemblybegin(V::AbstractVec{$PetscScalar})
        @chk ccall((:VecAssemblyBegin, $libpetsc), PetscErrorCode, (CVec,), V)
        return nothing
    end
    function assemblyend(V::AbstractVec{$PetscScalar})
        @chk ccall((:VecAssemblyEnd, $libpetsc), PetscErrorCode, (CVec,), V)
        return nothing
    end

    function ownershiprange(vec::AbstractVec{$PetscScalar})
        r_lo = Ref{$PetscInt}()
        r_hi = Ref{$PetscInt}()
        @chk ccall((:VecGetOwnershipRange, $libpetsc), PetscErrorCode,
          (CVec, Ptr{$PetscInt}, Ptr{$PetscInt}), vec, r_lo, r_hi)
        r_lo[]:(r_hi[]-$PetscInt(1))
    end

    function view(vec::AbstractVec{$PetscScalar}, viewer::Viewer{$PetscScalar}=ViewerStdout{$PetscScalar}(vec.comm))
        @chk ccall((:VecView, $libpetsc), PetscErrorCode,
                    (CVec, CPetscViewer),
                vec, viewer);
        return nothing
    end

    function localsize(v::AbstractVec{$PetscScalar})
        return r_sz[]
    end

    function unsafe_localarray(::Type{$PetscScalar}, cv::CVec; read::Bool=true, write::Bool=true)
        r_pv = Ref{Ptr{$PetscScalar}}()
        if write
            if read
                @chk ccall((:VecGetArray, $libpetsc), PetscErrorCode,
                    (CVec, Ptr{Ptr{$PetscScalar}}), cv, r_pv)
            else
                @chk ccall((:VecGetArrayWrite, $libpetsc), PetscErrorCode,
                    (CVec, Ptr{Ptr{$PetscScalar}}), cv, r_pv)
            end
        else
            @chk ccall((:VecGetArrayRead, $libpetsc), PetscErrorCode,
                (CVec, Ptr{Ptr{$PetscScalar}}), cv, r_pv)
        end
        r_sz = Ref{$PetscInt}()
        @chk ccall((:VecGetLocalSize, $libpetsc), PetscErrorCode,
            (CVec, Ptr{$PetscInt}), cv, r_sz)
        v = unsafe_wrap(Array, r_pv[], r_sz[]; own = false)

        if write
            if read
                finalizer(v) do v
                    @chk ccall((:VecRestoreArray, $libpetsc), PetscErrorCode, (CVec, Ptr{Ptr{$PetscScalar}}), cv, Ref(pointer(v)))
                    return nothing
                end
            else
                finalizer(v) do v
                    @chk ccall((:VecRestoreArrayWrite, $libpetsc), PetscErrorCode, (CVec, Ptr{Ptr{$PetscScalar}}), cv, Ref(pointer(v)))
                    return nothing
                end
            end
        else
            finalizer(v) do v
                @chk ccall((:VecRestoreArrayRead, $libpetsc), PetscErrorCode, (CVec, Ptr{Ptr{$PetscScalar}}), cv, Ref(pointer(v)))
                return nothing
            end
        end
        return v
    end
end

"""
    unsafe_localarray(PetscScalar, ptr:CVec; read=true, write=true)

Return an `Array{PetscScalar}` containing local portion of the PETSc data.

`finalize` should be called on the `Array` before the data can be used.

Use `read=false` if the array is write-only; `write=false` if read-only.
"""
unsafe_localarray

function Base.show(io::IO, ::MIME"text/plain", vec::AbstractVec)
    _show(io, vec)
end

VecSeq(X::Vector{T}; kwargs...) where {T} = VecSeq(MPI.COMM_SELF, X; kwargs...)
AbstractVec(X::AbstractVector) = VecSeq(X)


"""
    ownership_range(vec::AbstractVec)

The range of indices owned by this processor, assuming that the vectors are laid out with the first n1 elements on the first processor, next n2 elements on the second, etc. For certain parallel layouts this range may not be well defined.

Note: unlike the C function, the range returned is inclusive (`idx_first:idx_last`)

https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Vec/VecGetOwnershipRange.html
"""
ownershiprange

