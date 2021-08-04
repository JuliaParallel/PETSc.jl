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

function destroy(v::AbstractVec{PetscLib}) where {PetscLib}
    finalized(PetscLib) || LibPETSc.VecDestroy(PetscLib, v)
    return nothing
end

"""
    VecPtr(petsclib, v::CVec, seq_finalize)

Container type for a PETSc Vec that is just a raw pointer.

If the `seq_finalize` and `v` points to a sequential vector, then finalizer will
be set, otherwise the user is responsible for calling `destroy` (e.g., MPI
vectors).
"""
mutable struct VecPtr{PetscLib, PetscScalar} <:
               AbstractVec{PetscLib, PetscScalar}
    ptr::CVec
end
function VecPtr(
    petsclib::PetscLib,
    ptr::CVec,
    seq_finalize = true,
) where {PetscLib <: PetscLibType}
    v = VecPtr{PetscLib, petsclib.PetscScalar}(ptr)
    if seq_finalize && occursin("seq", getpetsctype(v))
        finalizer(destroy, v)
    end
    return v
end
VecPtr(::Type{PetscLib}, x...) where {PetscLib <: PetscLibType} =
    VecPtr(getlib(PetscLib), x...)

"""
    VecSeqWithArray(petsclib, v::Vector)

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
mutable struct VecSeqWithArray{PetscLib, PetscScalar} <:
               AbstractVec{PetscLib, PetscScalar}
    ptr::CVec
    array::Vector{PetscScalar}
end
Base.parent(v::VecSeqWithArray) = v.array

function VecSeqWithArray(
    petsclib::PetscLib,
    array::Vector{PetscScalar};
    blocksize = 1,
) where {PetscLib <: PetscLibType, PetscScalar}
    comm = MPI.COMM_SELF
    @assert initialized(petsclib)
    @assert PetscScalar == petsclib.PetscScalar
    v = VecSeqWithArray{PetscLib, PetscScalar}(C_NULL, array)
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
function VecSeqWithArray(
    ::Type{PetscLib},
    x...;
    kw...,
) where {PetscLib <: PetscLibType}
    VecSeqWithArray(getlib(PetscLib), x...; kw...)
end

"""
    VecSeq(petsclib, n::Int)

A standard, sequentially-stored serial PETSc vector for `petsclib.PetscScalar`
of length `n`.

# External Links
$(_doc_external("Vec/VecCreateSeq"))
"""
mutable struct VecSeq{PetscLib, PetscScalar} <:
               AbstractVec{PetscLib, PetscScalar}
    ptr::CVec
end

function VecSeq(petsclib::PetscLib, n::Int) where {PetscLib <: PetscLibType}
    comm = MPI.COMM_SELF
    @assert initialized(petsclib)
    v = VecSeq{PetscLib, petsclib.PetscScalar}(C_NULL)
    LibPETSc.VecCreateSeq(petsclib, comm, n, v)
    finalizer(destroy, v)
    return v
end

"""
    VecMPI(
         petsclib,
         comm:MPI.Comm,
         local_length;
         global_length = PETSC_DETERMINE
    )

An sequentially-stored MPI PETSc vector for `petsclib.PetscScalar` of local
length `local_length` and global length `global_length` without ghost elements.

If `global_length isa Int` then `local_length` can be set to `PETSC_DECIDE` in
which case PETSc will decide the local_length.

# External Links
$(_doc_external("Vec/VecCreateMPI"))

!!! note

    The user is responsible for calling `destroy(vec)` on the `Vec` since
    this cannot be handled by the garbage collector do to the MPI nature of the
    object.
"""
mutable struct VecMPI{PetscLib, PetscScalar} <:
               AbstractVec{PetscLib, PetscScalar}
    ptr::CVec
end

function VecMPI(
    petsclib::PetscLib,
    comm::MPI.Comm,
    local_length;
    global_length = PETSC_DETERMINE,
) where {PetscLib <: PetscLibType}
    @assert initialized(petsclib)
    @assert local_length != PETSC_DECIDE || global_length != PETSC_DETERMINE
    v = VecMPI{PetscLib, petsclib.PetscScalar}(C_NULL)
    LibPETSc.VecCreateMPI(petsclib, comm, local_length, global_length, v)
    return v
end

"""
    VecGhost(
         petsclib,
         comm:MPI.Comm,
         local_length
         ghost::Vector{PetscInt};
         global_length = PETSC_DETERMINE,
         num_ghost = length(ghost),
    )

An sequentially-stored MPI PETSc vector for `petsclib.PetscScalar` of local
length `local_length` and global length `global_length` with ghost elements.

If `global_length isa Int` then `local_length` can be set to `PETSC_DECIDE`.

# External Links
$(_doc_external("Vec/VecCreateGhost"))

!!! note

    The user is responsible for calling `destroy(vec)` on the `Vec` since
    this cannot be handled by the garbage collector do to the MPI nature of the
    object.
"""
mutable struct VecGhost{PetscLib, PetscScalar} <:
               AbstractVec{PetscLib, PetscScalar}
    ptr::CVec
end

function VecGhost(
    petsclib::PetscLib,
    comm::MPI.Comm,
    local_length,
    ghost::Vector{PetscInt};
    global_length = PETSC_DETERMINE,
    num_ghost = length(ghost),
) where {PetscLib <: PetscLibType, PetscInt}
    @assert initialized(petsclib)
    @assert PetscInt == PetscLib.PetscInt
    @assert local_length != PETSC_DECIDE || global_length != PETSC_DETERMINE
    v = VecGhost{PetscLib, petsclib.PetscScalar}(C_NULL)
    LibPETSc.VecCreateGhost(
        petsclib,
        comm,
        local_length,
        global_length,
        num_ghost,
        ghost,
        v,
    )
    return v
end

function Base.length(v::AbstractVec{PetscLib}) where {PetscLib}
    PetscInt = PetscLib.PetscInt
    r_sz = Ref{PetscInt}()
    LibPETSc.VecGetSize(PetscLib, v, r_sz)
    return r_sz[]
end

function Base.getindex(v::AbstractVec{PetscLib}, i::Integer) where {PetscLib}
    vals = [PetscLib.PetscScalar(0)]
    LibPETSc.VecGetValues(PetscLib, v, 1, Ref{PetscLib.PetscInt}(i - 1), vals)
    return vals[1]
end

"""
    setvalues!(
        v::AbstractVec,
        indices::Vector{PetscInt},
        vals::Array{PetscScalar},
        insertmode::InsertMode,
    )

Assign the values `vals` in 0-based global `indices` of `vec`. The `insertmode`
can be `INSERT_VALUES` or `ADD_VALUES`.

!!! warning
    This function uses 0-based indexing!

# External Links
$(_doc_external("Vec/VecSetValues"))
"""
function setvalues!(
    v::AbstractVec{PetscLib},
    idxs0::Vector{PetscInt},
    vals::Array{PetscScalar},
    insertmode::InsertMode;
    num_idxs = length(idxs0),
) where {PetscLib, PetscInt, PetscScalar}
    @assert length(vals) >= num_idxs
    @assert PetscInt == PetscLib.PetscInt
    @assert PetscScalar == PetscLib.PetscScalar
    LibPETSc.VecSetValues(PetscLib, v, num_idxs, idxs0, vals, insertmode)
    return nothing
end

"""
    getvalues!(
        vals::Array{PetscScalar},
        v::AbstractVec,
        indices::Vector{PetscInt},
    )

Get the 0-based global `indices` of `vec` into the preallocated array `vals`.

!!! warning
    This function uses 0-based indexing!

# External Links
$(_doc_external("Vec/VecGetValues"))
"""
function getvalues!(
    vals::Array{PetscScalar},
    v::AbstractVec{PetscLib},
    idxs0::Vector{PetscInt};
    num_idxs = length(idxs0),
) where {PetscLib, PetscInt, PetscScalar}
    @assert length(vals) >= num_idxs
    @assert PetscInt == PetscLib.PetscInt
    @assert PetscScalar == PetscLib.PetscScalar
    LibPETSc.VecGetValues(PetscLib, v, num_idxs, idxs0, vals)
    return nothing
end

function Base.setindex!(
    v::AbstractVec{PetscLib},
    val,
    i::Integer,
) where {PetscLib}
    LibPETSc.VecSetValues(
        PetscLib,
        v,
        1,
        Ref{PetscLib.PetscInt}(i - 1),
        Ref{PetscLib.PetscScalar}(val),
        INSERT_VALUES,
    )

    return val
end

function locallength(v::AbstractVec{PetscLib}) where {PetscLib}
    PetscInt = PetscLib.PetscInt
    r_sz = Ref{PetscInt}()
    LibPETSc.VecGetLocalSize(PetscLib, v, r_sz)
    return r_sz[]
end

function LinearAlgebra.norm(
    v::AbstractVec{PetscLib},
    normtype::LibPETSc.NormType = LibPETSc.NORM_2,
) where {PetscLib}
    PetscReal = PetscLib.PetscReal
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
    ownershiprange(vec::AbstractVec, [base_one = true])

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
    PetscInt = PetscLib.PetscInt
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
    PetscScalar = PetscLib.PetscScalar
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
    withlocalarray!(f!, x::AbstractVec; read=true, write=true)
    withlocalarray!(f!, xs...; read=true, write=true)
    withlocalarray!(f!, xs::NTuple{N, AbstractVec}...; read=true, write=true)

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
function withlocalarray!(f!, vecs::NTuple{N, AbstractVec}; kwargs...) where {N}
    arrays = map(vecs) do v
        unsafe_localarray(v; kwargs...)
    end
    f!(arrays...)
    map(arrays) do array
        Base.finalize(array)
    end
end
withlocalarray!(f!, vecs...; kwargs...) = withlocalarray!(f!, vecs; kwargs...)

"""
    assemblybegin!(vec::AbstractVec)

Begin assembling `vec`

# External Links
$(_doc_external("Vec/VecAssemblyBegin"))
"""
function assemblybegin!(vec::AbstractVec{PetscLib}) where {PetscLib}
    LibPETSc.VecAssemblyBegin(PetscLib, vec)
    return nothing
end

"""
    assemblyend!(vec::AbstractVec)

Finish assembling `vec`

# External Links
$(_doc_external("Vec/VecAssemblyEnd"))
"""
function assemblyend!(vec::AbstractVec{PetscLib}) where {PetscLib}
    LibPETSc.VecAssemblyEnd(PetscLib, vec)
    return nothing
end

"""
    assemble!(v::AbstractVec)

Assemble the vector `v`.

For overlapping assembly see [`assemblybegin!`](@ref) and  [`assemblyend!`](@ref)

# External Links
$(_doc_external("Vec/VecAssemblyBegin"))
$(_doc_external("Vec/VecAssemblyEnd"))
"""
function assemble!(v::AbstractVec)
    assemblybegin!(v)
    assemblyend!(v)
    return v
end

mutable struct LocalVec{PetscLib, PetscScalar, GVec} <:
               AbstractVec{PetscLib, PetscScalar}
    ptr::CVec
    gvec::GVec
end
function LocalVec(gvec::AbstractVec{PetscLib}) where {PetscLib}
    GVec = typeof(gvec)
    PetscScalar = PetscLib.PetscScalar
    LocalVec{PetscLib, PetscScalar, GVec}(C_NULL, gvec)
end

"""
    getlocalform(vec::AbstractVec)

Obtains the local ghosted representation of a [`Vec`](@ref).

!!! note

    When done with the object the user should call [`restorelocalform!`](@ref)

# External Links
$(_doc_external("Vec/VecGhostGetLocalForm"))
"""
function getlocalform(gvec::AbstractVec{PetscLib}) where {PetscLib}
    lvec = LocalVec(gvec)
    LibPETSc.VecGhostGetLocalForm(PetscLib, gvec, lvec)
    if lvec.ptr == C_NULL
        restorelocalform!(lvec)
        throw(ArgumentError("no local form for vector"))
    end
    return lvec
end

"""
    restorelocalform!(local_vec::LocalVec)

Restore the `local_vec` to the associated global vector after a call to
[`getlocalform`](@ref).

# External Links
$(_doc_external("Vec/VecGhostRestoreLocalForm"))
"""
function restorelocalform!(lvec::LocalVec{PetscLib}) where {PetscLib}
    LibPETSc.VecGhostRestoreLocalForm(PetscLib, lvec.gvec, lvec)
    lvec.ptr = C_NULL
    return lvec.gvec
end

"""
    withlocalform(f::Function, vec::AbstractVec)

Convert `vec` to a `LocalVec` and apply the function `f!`.

```julia-repl
julia> withlocalform(vec) do l_vec
   # Do something with l_vec
end
```

!!! note

    This wrapper handles the calling of [`restorelocalform!`](@ref) before
    returning.
"""
function withlocalform(f!, vec::AbstractVec)
    lvec = getlocalform(vec)
    f!(lvec)
    restorelocalform!(lvec)
end

"""
    ghostupdatebegin!(
        vec::AbstractVec,
        insertmode = INSERT_VALUES,
        scattermode = SCATTER_FORWARD,
    )

Begins scattering `vec` to the local or global representations

# External Links
$(_doc_external("Vec/VecGhostUpdateBegin"))
"""
function ghostupdatebegin!(
    vec::AbstractVec{PetscLib},
    insertmode = INSERT_VALUES,
    scattermode = SCATTER_FORWARD,
) where {PetscLib}
    LibPETSc.VecGhostUpdateBegin(PetscLib, vec, insertmode, scattermode)
    return nothing
end

"""
    ghostupdateend!(
        vec::AbstractVec,
        insertmode = INSERT_VALUES,
        scattermode = SCATTER_FORWARD,
    )

Finishes scattering `vec` to the local or global representations

# External Links
$(_doc_external("Vec/VecGhostUpdateEnd"))
"""
function ghostupdateend!(
    vec::AbstractVec{PetscLib},
    insertmode = INSERT_VALUES,
    scattermode = SCATTER_FORWARD,
) where {PetscLib}
    LibPETSc.VecGhostUpdateEnd(PetscLib, vec, insertmode, scattermode)
    return nothing
end

"""
    getpetsctype(vec::AbstractVec)

return a string with the vector type

# External Links
$(_doc_external("Vec/VecGetType"))
"""
function getpetsctype(vec::AbstractVec{PetscLib}) where {PetscLib}
    name_r = Ref{LibPETSc.VecType}()
    LibPETSc.VecGetType(PetscLib, vec, name_r)
    return unsafe_string(name_r[])
end

function Base.similar(v::AbstractVec{PetscLib}) where {PetscLib}
    r_x = Ref{CVec}()
    LibPETSc.VecDuplicate(PetscLib, v, r_x)
    x = VecPtr(PetscLib, r_x[], true)
    return x
end
