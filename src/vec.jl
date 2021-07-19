# AbstractVec
#   - VecSeq: wrap
#   - VecMPI (TODO)
#   - VecGhost (TODO)
# for the MPI variants we won't be able to attach finalizers, as destroy needs to be called collectively.

const CVec = Ptr{Cvoid}

abstract type AbstractVec{PetscLib, PetscScalar} <: AbstractVector{PetscScalar} end

Base.eltype(
    ::Type{V},
) where {
    V <: AbstractVec{PetscLib, PetscScalar},
} where {PetscLib, PetscScalar} = PetscScalar
Base.eltype(
    v::AbstractVec{PetscLib, PetscScalar},
) where {PetscLib, PetscScalar} = PetscScalar
Base.size(v::AbstractVec) = (length(v),)

# allows us to pass XXVec objects directly into CVec ccall signatures
Base.cconvert(::Type{CVec}, obj::AbstractVec) = obj.ptr
# allows us to pass XXVec objects directly into Ptr{CVec} ccall signatures
Base.unsafe_convert(::Type{Ptr{CVec}}, obj::AbstractVec) =
    convert(Ptr{CVec}, pointer_from_objref(obj))


"""
    VecSeq(petsclib, v::Vector)

A standard, sequentially-stored serial PETSc vector, wrapping the Julia vector
`v`.

This reuses the array `v` as storage, and so `v` should not be `resize!`-ed or
otherwise have its length modified while the PETSc object exists.

This should only be need to be called for more advanced uses, for most simple
usecases, users should be able to pass `Vector`s directly and have the wrapping
performed automatically

# External Links
$(_doc_external("Vec/VecCreateSeqWithArray"))
"""
mutable struct VecSeq{PetscLib, PetscScalar} <:
               AbstractVec{PetscLib, PetscScalar}
    ptr::CVec
    array::Vector{PetscScalar}
end
Base.parent(v::VecSeq) = v.array

function VecSeq(
    petsclib::PetscLib,
    array::Vector{PetscScalar};
    blocksize = 1,
) where {PetscLib, PetscScalar}
    comm = MPI.COMM_SELF
    @assert initialized(petsclib)
    @assert PetscScalar == petsclib.PetscScalar
    v = VecSeq{PetscLib, PetscScalar}(C_NULL, array)
    LibPETSc.VecCreateSeqWithArray(
        petsclib,
        comm,
        blocksize,
        length(array),
        array,
        v,
    )
    finalizer(destroy, v)
    return v
end

function destroy(v::AbstractVec{PetscLib}) where {PetscLib}
    finalized(PetscLib) || LibPETSc.VecDestroy(PetscLib, v)
    return nothing
end

function Base.length(v::AbstractVec{PetscLib}) where {PetscLib}
    petsclib = getlib(PetscLib)
    PetscInt = petsclib.PetscInt
    r_sz = Ref{PetscInt}()
    LibPETSc.VecGetSize(PetscLib, v, r_sz)
    return r_sz[]
end

function locallength(v::AbstractVec{PetscLib}) where {PetscLib}
    petsclib = getlib(PetscLib)
    PetscInt = petsclib.PetscInt
    r_sz = Ref{PetscInt}()
    LibPETSc.VecGetLocalSize(PetscLib, v, r_sz)
    return r_sz[]
end

function LinearAlgebra.norm(
    v::AbstractVec{PetscLib},
    normtype::LibPETSc.NormType = LibPETSc.NORM_2,
) where {PetscLib}
    petsclib = getlib(PetscLib)
    PetscReal = petsclib.PetscReal
    r_val = Ref{PetscReal}()
    LibPETSc.VecNorm(PetscLib, v, normtype, r_val)
    return r_val[]
end

function view(
    vec::AbstractVec{PetscLib},
    viewer = LibPETSc.PETSC_VIEWER_STDOUT_(PetscLib, getcomm(vec)),
) where {PetscLib}
    LibPETSc.VecView(PetscLib, vec, viewer)
    return nothing
end
Base.show(io::IO, vec::AbstractVec) = _show(io, vec)
Base.show(io::IO, ::MIME"text/plain", vec::AbstractVec) = _show(io, vec)

"""
    ownership_range(vec::AbstractVec, [base_one = true])

The range of indices owned by this processor, assuming that the `vec` is laid
out with the first `n1` elements on the first processor, next `n2` elements on
the second, etc. For certain parallel layouts this range may not be well
defined.

If the optional argument `base_one == true` then base-1 indexing is used,
otherwise base-0 index is used.

!!! note

    unlike the C function, the range returned is inclusive (`idx_first:idx_last`)

# External Links
$(_doc_external("Vec/VecGetOwnershipRange"))
"""
function ownershiprange(
    vec::AbstractVec{PetscLib},
    base_one::Bool = true,
) where {PetscLib}
    petsclib = getlib(PetscLib)
    PetscInt = petsclib.PetscInt
    r_lo = Ref{PetscInt}()
    r_hi = Ref{PetscInt}()
    LibPETSc.VecGetOwnershipRange(PetscLib, vec, r_lo, r_hi)
    return base_one ? ((r_lo[] + PetscInt(1)):(r_hi[])) :
           ((r_lo[]):(r_hi[] - PetscInt(1)))
end

"""
    unsafe_localarray(vec::AbstractVec; read=true, write=true)

Return an `Array{PetscScalar}` containing local portion of the PETSc `vec`

Use `read=false` if the array is write-only; `write=false` if read-only.

!!! note
    `Base.finalize` should be called on the `Array` before the data can be used.

# External Links
$(_doc_external("Vec/VecGetArray"))
$(_doc_external("Vec/VecGetArrayWrite"))
$(_doc_external("Vec/VecGetArrayRead"))
$(_doc_external("Vec/VecRestoreArray"))
$(_doc_external("Vec/VecRestoreArrayWrite"))
$(_doc_external("Vec/VecRestoreArrayRead"))
"""
function unsafe_localarray(
    vec::AbstractVec{PetscLib};
    read::Bool = true,
    write::Bool = true,
) where {PetscLib}
    petsclib = getlib(PetscLib)
    PetscScalar = petsclib.PetscScalar
    r_pv = Ref{Ptr{PetscScalar}}()
    if write && read
        LibPETSc.VecGetArray(PetscLib, vec, r_pv)
    elseif write
        LibPETSc.VecGetArrayWrite(PetscLib, vec, r_pv)
    elseif read
        LibPETSc.VecGetArrayRead(PetscLib, vec, r_pv)
    else
        error("either read or write should be true")
    end
    sz = locallength(vec)
    v = unsafe_wrap(Array, r_pv[], sz; own = false)

    if write && read
        finalizer(v) do v
            LibPETSc.VecRestoreArray(PetscLib, vec, Ref(pointer(v)))
            return nothing
        end
    elseif write
        finalizer(v) do v
            LibPETSc.VecRestoreArrayWrite(PetscLib, vec, Ref(pointer(v)))
            return nothing
        end
    elseif read
        finalizer(v) do v
            LibPETSc.VecRestoreArrayRead(PetscLib, vec, Ref(pointer(v)))
            return nothing
        end
    end
    return v
end

"""
    with_unsafe_localarray!(
        f!,
        x::AbstractVec;
        read=true,
        write=true,
    )

Convert `x` to an `Array{PetscScalar}` using [`unsafe_localarray`](@ref) and
apply the function `f!`.

Use `read=false` if the array is write-only; `write=false` if read-only.

# Examples
```julia-repl
julia> map_unsafe_localarray(x; write=true) do x
   @. x .*= 2
end

!!! note
    `Base.finalize` is automatically called on the array.
"""
function with_unsafe_localarray!(f!, v::AbstractVec; kwargs...)
    array = unsafe_localarray(v; kwargs...)
    f!(array)
    Base.finalize(array)
end

#=
@for_libpetsc begin
    function assemblybegin(V::AbstractVec{$PetscScalar})
        @chk ccall((:VecAssemblyBegin, $libpetsc), PetscErrorCode, (CVec,), V)
        return nothing
    end
    function assemblyend(V::AbstractVec{$PetscScalar})
        @chk ccall((:VecAssemblyEnd, $libpetsc), PetscErrorCode, (CVec,), V)
        return nothing
    end
end
=#
