# these are julia-specific PETSc functions, which makes PETSc a bit easier to use
# a (nearly) full wrapper is in wrapping/

import .LibPETSc: AbstractPetscVec, PetscVec, CVec

# Custom display for REPL
function Base.show(io::IO, v::AbstractPetscVec{PetscLib}) where {PetscLib}
    if v.ptr == C_NULL
        print(io, "PETSc Vec (null pointer)")
        return
    end
    # VecGetType internally calls VecInitializePackage which queries the PETSc
    # options database.  Calling it before PETSc is initialised causes a C-level
    # SIGSEGV that cannot be caught with try/catch.
    if !initialized(PetscLib)
        print(io, "PETSc Vec (PETSc not initialized)")
        return
    end
    try
        ty = LibPETSc.VecGetType(PetscLib, v)
        si = LibPETSc.VecGetSize(PetscLib, v)
        print(io, "PETSc $ty Vec; length=$si")
    catch
        print(io, "PETSc Vec (type not set)")
    end
    return nothing
end

"""
    VecPtr(petsclib, ptr::CVec, own::Bool)

Container type for a PETSc Vec that is just a raw pointer.

If `own` is `true`, the finalizer is set on the vector; calling `destroy` when
`own` is `false` is a no-op.
"""
mutable struct VecPtr{PetscLib} <:
               AbstractPetscVec{PetscLib}
    ptr::CVec
    own::Bool
end
function VecPtr(
    petsclib::PetscLib,
    ptr::CVec,
    own,
) where {PetscLib <: PetscLibType}
    v = VecPtr{PetscLib}(ptr, own)
    #comm = getcomm(v)

    #comm = MPI.Comm()
    comm = LibPETSc.PetscObjectGetComm(getlib(PetscLib), v)

    if own && MPI.Comm_size(comm) == 1
        finalizer(destroy, v)
    end
    return v
end
VecPtr(::Type{PetscLib}, x...) where {PetscLib <: PetscLibType} = VecPtr(getlib(PetscLib), x...)


"""
    VecSeq(petsclib, n::Int)

A standard, sequentially-stored serial PETSc vector for `petsclib.PetscScalar`
of length `n`.

# External Links
$(_doc_external("Vec/VecCreateSeq"))
"""
function VecSeq(petsclib::PetscLib, n::Int) where {PetscLib <: PetscLibType}
    comm = MPI.COMM_SELF
    @assert initialized(petsclib)
    v = LibPETSc.VecCreateSeq(petsclib, comm, n)
    finalizer(destroy, v)
    return v
end


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
function VecSeq(
    petsclib::PetscLib,
    array::Vector{PetscScalar};
    blocksize = 1,
) where {PetscLib <: PetscLibType, PetscScalar}
    comm = MPI.COMM_SELF
    @assert initialized(petsclib)
    @assert PetscScalar == petsclib.PetscScalar
    PetscInt = petsclib.PetscInt
    v = LibPETSc.VecCreateSeqWithArray(
        petsclib,
        comm,
        PetscInt(blocksize),
        PetscInt(length(array)),
        array,
    )
    finalizer(destroy, v)
    return v
end


# =============================================================================
# Multiple dispatch to make PetscVec behave like Julia Vector
# =============================================================================

# Treat PETSc vectors as 1-D array-like objects for broadcasting/indexing
Base.ndims(::Type{<:AbstractPetscVec}) = 1
Base.IndexStyle(::Type{<:AbstractPetscVec}) = IndexLinear()
Base.axes(v::AbstractPetscVec) = (Base.OneTo(length(v)),)
Base.BroadcastStyle(::Type{<:AbstractPetscVec}) = Broadcast.DefaultArrayStyle{1}()

# Array interface - size and length
Base.size(v::AbstractPetscVec{PetscLib}) where {PetscLib} = LibPETSc.VecGetSize(PetscLib,v)
Base.length(v::AbstractPetscVec{PetscLib}) where {PetscLib} = prod(size(v))
Base.lastindex(v::AbstractPetscVec{PetscLib}) where {PetscLib} = length(v)
Base.similar(v::AbstractPetscVec{PetscLib}) where {PetscLib} =  LibPETSc.VecDuplicate(getlib(PetscLib), v)
type(m::AbstractPetscVec{PetscLib}) where {PetscLib} = LibPETSc.VecGetType(PetscLib,m)

function Base.getindex(v::AbstractPetscVec{PetscLib}, i::Integer) where {PetscLib} 
    PetscInt = inttype(PetscLib)
    val = LibPETSc.VecGetValues(PetscLib,v, PetscInt(1), PetscInt.([i-1]))
    return val[1]
end

# Range indexing
function Base.getindex(v::AbstractPetscVec{PetscLib}, r::AbstractRange) where {PetscLib}
    PetscInt = inttype(PetscLib)
    val = LibPETSc.VecGetValues(PetscLib,v, PetscInt(length(r)), PetscInt.(Vector(r) .- 1))
    return val
end

# Get all values
function Base.getindex(v::AbstractPetscVec{PetscLib}, ::Colon) where {PetscLib}
    n  = length(v)
    return getindex(v, 1:n)
end

Base.isapprox(v::AbstractPetscVec{PetscLib}, w::AbstractPetscVec{PetscLib}; kwargs...) where {PetscLib} = all(v[:] .≈ w[:]; kwargs...)

function Base.setindex!(v::AbstractPetscVec{PetscLib}, val, i::Integer) where {PetscLib} 
     PetscInt = inttype(PetscLib)
     PetscScalar = PETSc.scalartype(PetscLib)
     LibPETSc.VecSetValues(PetscLib,v, PetscInt(1), PetscInt.([i-1]), [PetscScalar(val)], PETSc.INSERT_VALUES)
     return nothing
end

function Base.setindex!(v::AbstractPetscVec{PetscLib}, vals, r::AbstractRange) where {PetscLib} 
    PetscInt = inttype(PetscLib)
    LibPETSc.VecSetValues(PetscLib,v, PetscInt(length(r)), PetscInt.(Vector(r) .- 1), vals, PETSc.INSERT_VALUES)
    return nothing
end

Base.fill!(v::AbstractPetscVec{PetscLib}, val) where {PetscLib} = LibPETSc.VecSet(PetscLib,v, PetscLib.PetscScalar(val))

# Broadcasting assignment support (dest is not an AbstractArray)
function Base.copyto!(dest::AbstractPetscVec{PetscLib}, bc::Base.Broadcast.Broadcasted) where {PetscLib}
    # Evaluate the broadcasted RHS
    rhs = Base.materialize(bc)
    if rhs isa Number
        # Fast path for scalar RHS
        Base.fill!(dest, rhs)
        return dest
    end
    # Array-like RHS: shape must match
    axes(dest) == axes(rhs) || throw(DimensionMismatch("broadcast axes $(axes(rhs)) do not match destination axes $(axes(dest))"))
    @inbounds for i in eachindex(rhs)
        dest[i] = rhs[i]
    end
    return dest
end

Base.materialize!(dest::AbstractPetscVec, bc::Base.Broadcast.Broadcasted) = (Base.copyto!(dest, bc); dest)

# Iterator interface
Base.iterate(v::AbstractPetscVec{PetscLib})  where {PetscLib} = iterate(v, 1)
function Base.iterate(v::AbstractPetscVec{PetscLib}, state)  where {PetscLib}
    if state > length(v)
        return nothing
    end
    return (v[state], state + 1)
end

"""
    assemble!(A::PetscVec) 

Assembles a PETSc vector after setting values.
"""
function assemble!(A::AbstractPetscVec{PetscLib}) where {PetscLib}
    LibPETSc.VecAssemblyBegin(PetscLib, A)
    LibPETSc.VecAssemblyEnd(PetscLib, A)
end



destroy(m::AbstractPetscVec{PetscLib}) where {PetscLib} = LibPETSc.VecDestroy(PetscLib,m)


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
    vec::AbstractPetscVec{PetscLib};
    read::Bool = true,
    write::Bool = true,
) where {PetscLib}
    if write && read
        v = LibPETSc.VecGetArray(PetscLib, vec)
    elseif write
        v = LibPETSc.VecGetArrayWrite(PetscLib, vec)
    elseif read
        v = LibPETSc.VecGetArrayRead(PetscLib, vec)
    else
        error("either read or write should be true")
    end
    
    if write && read
        finalizer(v) do v
            LibPETSc.VecRestoreArray(PetscLib, vec, v)
            return nothing
        end
    elseif write
        finalizer(v) do v
            LibPETSc.VecRestoreArrayWrite(PetscLib, vec, v)         
            return nothing
        end
    elseif read
        finalizer(v) do v
            LibPETSc.VecRestoreArrayRead(PetscLib, vec,  v)
            return nothing
        end
    end
    return v
end


# ── Memory backend type hierarchy ─────────────────────────────────────────────
#
# Extensions add their own backend singletons (e.g. `CUDABackend`) and overload
# `_memtype_backend(::Val{PETSC_MEMTYPE_DEVICE})` to return them.  The base
# package handles only `PETSC_MEMTYPE_HOST` → `HostBackend`.

"""
    AbstractPETScMemBackend

Abstract supertype for PETSc memory backends.  The base package defines only
[`HostBackend`](@ref).  GPU extensions add their own (e.g. `CUDABackend`).
"""
abstract type AbstractPETScMemBackend end

"""
    HostBackend <: AbstractPETScMemBackend

Singleton dispatch type representing host (CPU) memory.
"""
struct HostBackend <: AbstractPETScMemBackend end

"""
    _memtype_backend(mtype::PetscMemType) → AbstractPETScMemBackend

Convert a `PetscMemType` runtime enum value to a singleton dispatch type.
GPU extensions overload `_memtype_backend(::Val{MT})` for their specific
`PetscMemType` values (e.g. `PETSC_MEMTYPE_DEVICE` for CUDA).
"""
_memtype_backend(::Val{LibPETSc.PETSC_MEMTYPE_HOST}) = HostBackend()
_memtype_backend(::Val{MT}) where {MT} =
    error("No GPU backend loaded for PetscMemType $MT — load CUDA.jl, AMDGPU.jl, …")
_memtype_backend(mt::LibPETSc.PetscMemType) = _memtype_backend(Val(mt))

# ── Device-aware local array access ───────────────────────────────────────────
#
# `_unsafe_localarray` is the unified entry point: it calls
# `VecGetArray*AndMemType`, converts the returned `PetscMemType` to a backend
# singleton via `_memtype_backend`, and dispatches to `_wrap_localarray`.
# GPU extensions add `_wrap_localarray` methods for their own backend types.
#
# The typed overload `_unsafe_localarray(::Type{A}, vec; ...)` additionally
# asserts that the returned array is of type `A`, giving a clear error when a
# Vec is on an unexpected device.

function _unsafe_localarray(
    vec::AbstractPetscVec{PetscLib};
    read::Bool = true,
    write::Bool = true,
) where {PetscLib}
    if write && read
        cpu_arr, mtype = LibPETSc.VecGetArrayAndMemType(PetscLib, vec)
    elseif write
        cpu_arr, mtype = LibPETSc.VecGetArrayWriteAndMemType(PetscLib, vec)
    else
        cpu_arr, mtype = LibPETSc.VecGetArrayReadAndMemType(PetscLib, vec)
    end
    return _wrap_localarray(cpu_arr, _memtype_backend(mtype), vec; read, write)
end

function _unsafe_localarray(
    ::Type{A},
    vec::AbstractPetscVec;
    read::Bool = true,
    write::Bool = true,
) where {A <: AbstractArray}
    arr = _unsafe_localarray(vec; read, write)
    arr isa A && return arr
    Base.finalize(arr)   # release the PETSc handle before throwing
    throw(ArgumentError(
        "expected array of type $A but Vec returned $(typeof(arr)). " *
        "Check that the Vec lives on the expected device."
    ))
end

function _wrap_localarray(
    cpu_arr, ::HostBackend, vec::AbstractPetscVec{PetscLib};
    read::Bool, write::Bool,
) where {PetscLib}
    finalizer(cpu_arr) do a
        if write && read
            LibPETSc.VecRestoreArrayAndMemType(PetscLib, vec, a)
        elseif write
            LibPETSc.VecRestoreArrayWriteAndMemType(PetscLib, vec, a)
        else
            LibPETSc.VecRestoreArrayReadAndMemType(PetscLib, vec, a)
        end
        return nothing
    end
    return cpu_arr
end

# Fallback: no backend loaded for this PetscMemType.
function _wrap_localarray(cpu_arr, b::AbstractPETScMemBackend, vec; kw...)
    error("_wrap_localarray not implemented for backend $(typeof(b)) — " *
          "load the corresponding GPU package (e.g. CUDA.jl)")
end

"""
    determine_memtype(vecs...) → Type{<:AbstractArray}

Query the `PetscMemType` of each Vec and return the corresponding array type.
Errors if the Vecs are on heterogeneous devices (different `PetscMemType`
values), since a single `withlocalarray!` call cannot handle mixed backends.
Returns `Vector` when all Vecs are host-resident.

Extensions overload `_array_type(::Val{MT})` for a `PetscMemType` enum value
`MT` to register the corresponding array type (e.g. `PETSC_MEMTYPE_DEVICE` →
`CuArray`).
"""
function determine_memtype(vecs::AbstractPetscVec...)
    mtypes = map(vecs) do v
        PetscLib = typeof(v).parameters[1]
        arr, mtype = LibPETSc.VecGetArrayReadAndMemType(PetscLib, v)
        LibPETSc.VecRestoreArrayReadAndMemType(PetscLib, v, arr)
        mtype
    end
    allequal(mtypes) || throw(ArgumentError(
        "Vecs are on heterogeneous devices: $(unique(mtypes)). " *
        "Use withlocalarray!(::Type{A}, ...) to handle each backend explicitly."
    ))
    return _array_type(Val(first(mtypes)))
end

_array_type(::Val{LibPETSc.PETSC_MEMTYPE_HOST}) = Vector
_array_type(::Val{MT}) where {MT} =
    error("No array type registered for PetscMemType $MT — load the corresponding GPU package (e.g. CUDA.jl)")
# GPU extensions add: _array_type(::Val{LibPETSc.PETSC_MEMTYPE_DEVICE}) = CuArray

"""
    withlocalarray!(
        f!,
        vecs::NTuple{N, AbstractVec};
        read::Union{Bool, NTuple{N, Bool}} = true,
        write::Union{Bool, NTuple{N, Bool}} = true,
    )
    withlocalarray!(::Type{A}, f!, vecs...; read, write) where {A <: AbstractArray}

Apply `f!` to local array views of `vecs`.

Uses `VecGetArray*AndMemType` internally.  When a GPU backend extension (e.g.
`PETScCUDAExt`) is loaded and a Vec lives on the device, `f!` receives a device
array (e.g. `CuArray`) — zero-copy, no host↔device transfer.  When all Vecs
are host-resident, `f!` receives plain `Array`s.

The optional `::Type{A}` first argument asserts that every array returned from
`VecGetArray*AndMemType` is of type `A`.  This is useful when Vecs are known to
be heterogeneous (e.g. some on host, some on device) and you need a type-stable
code path: passing `CuArray` will error immediately if any Vec is host-resident,
rather than silently returning a `Vector`.

Use `read=false` if the array is write-only; `write=false` if read-only.

!!! note
    Operations inside `f!` must be compatible with the actual array type.
    Scalar indexing is not supported on GPU arrays; use broadcasting or GPU
    kernels instead.

# Examples
```julia-repl
julia> withlocalarray!(x; write=true) do x
   @. x .*= 2
end

julia> withlocalarray!(
           x,
           y;
           read = (false, true),
           write = (true, false)
       ) do x, y
   @. x .= 2 .+ y
end
```

!!! note
    `Base.finalize` is automatically called on the array.
"""
function withlocalarray!(
    f!,
    vecs::NTuple{N, AbstractPetscVec};
    kwargs...,
) where {N}
    A = determine_memtype(vecs...)
    return withlocalarray!(A, f!, vecs; kwargs...)
end
withlocalarray!(f!, vecs...; kwargs...) = withlocalarray!(f!, vecs; kwargs...)

function withlocalarray!(
    ::Type{A},
    f!,
    vecs::NTuple{N, AbstractPetscVec};
    read::Union{Bool, NTuple{N, Bool}} = true,
    write::Union{Bool, NTuple{N, Bool}} = true,
) where {A <: AbstractArray, N}
    read isa NTuple{N, Bool} || (read = ntuple(_ -> read, N))
    write isa NTuple{N, Bool} || (write = ntuple(_ -> write, N))

    arrays = map(vecs, read, write) do v, r, w
        _unsafe_localarray(A, v; read = r, write = w)
    end
    val = f!(arrays...)
    map(Base.finalize, arrays)
    return val
end
withlocalarray!(::Type{A}, f!, vecs...; kwargs...) where {A <: AbstractArray} =
    withlocalarray!(A, f!, vecs; kwargs...)


"""
    ghostupdatebegin!(
        vec::AbstractPetscVec,
        insertmode = INSERT_VALUES,
        scattermode = SCATTER_FORWARD,
    )

Begins scattering `vec` to the local or global representations

# External Links
$(_doc_external("Vec/VecGhostUpdateBegin"))
"""
function ghostupdatebegin!(
    vec::AbstractPetscVec{PetscLib},
    insertmode = INSERT_VALUES,
    scattermode = SCATTER_FORWARD,
) where {PetscLib}
    LibPETSc.VecGhostUpdateBegin(PetscLib, vec, insertmode, scattermode)
    return nothing
end

"""
    ghostupdateend!(
        vec::AbstractPetscVec,
        insertmode = INSERT_VALUES,
        scattermode = SCATTER_FORWARD,
    )

Finishes scattering `vec` to the local or global representations

# External Links
$(_doc_external("Vec/VecGhostUpdateEnd"))
"""
function ghostupdateend!(
    vec::AbstractPetscVec{PetscLib},
    insertmode = INSERT_VALUES,
    scattermode = SCATTER_FORWARD,
) where {PetscLib}
    LibPETSc.VecGhostUpdateEnd(PetscLib, vec, insertmode, scattermode)
    return nothing
end

"""
    ghostupdate!(
        vec::AbstractPetscVec,
        insertmode = INSERT_VALUES,
        scattermode = SCATTER_FORWARD,
    )

Finishes scattering `vec` to the local or global representations

# External Links
$(_doc_external("Vec/VecGhostUpdateEnd"))
"""
function ghostupdate!(
    vec::AbstractPetscVec{PetscLib},
    insertmode = INSERT_VALUES,
    scattermode = SCATTER_FORWARD,
) where {PetscLib}
    ghostupdatebegin!(vec,insertmode,scattermode)
    ghostupdateend!(vec,insertmode,scattermode)
    return nothing
end

"""
    v = VecSeq(petsclib, comm, array)
Creates a sequential PETSc vector of length `n` given a julia array `array`` 
"""
function VecSeq(petsclib::PetscLib, comm, x::Vector) where {PetscLib <: PetscLibType}
    @assert initialized(petsclib)
    
    v = LibPETSc.VecCreateSeqWithArray(petsclib,comm, 1, length(x), x)    # solution vector
    finalizer(destroy, v)

    return v
end



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
    vec::AbstractPetscVec{PetscLib},
    base_one::Bool = true,
) where {PetscLib}
    PetscInt = PetscLib.PetscInt
    r_lo = Ref{PetscInt}()
    r_hi = Ref{PetscInt}()
    r_lo, r_hi = LibPETSc.VecGetOwnershipRange(PetscLib, vec)
    return base_one ? ((r_lo[] + PetscInt(1)):(r_hi[])) :
           ((r_lo[]):(r_hi[] - PetscInt(1)))
end

# Overload norm function
function LinearAlgebra.norm(
    v::AbstractPetscVec{PetscLib},
    normtype::LibPETSc.NormType = LibPETSc.NORM_2,
) where {PetscLib}
    PetscReal = PetscLib.PetscReal
    r_val = LibPETSc.VecNorm(PetscLib, v, normtype)
    return r_val
end

# ── GPU-aware array access helpers ────────────────────────────────────────────
#
# `get_petsc_arrays` calls `VecGetArrayAndMemType` on both Vecs, converts the
# returned `PetscMemType` values to backend singletons, and dispatches to
# `_get_petsc_arrays_impl`.  The base package handles the pure-CPU case
# (HostBackend × HostBackend).  GPU extensions add `_get_petsc_arrays_impl`
# methods for their backend combinations and a matching
# `_restore_petsc_arrays_impl` method dispatched by `restore_petsc_arrays`.
#
# Return tuple:  (fx, lx, fx_arr, lx_arr, fx_bounce)
#   CPU:  fx, lx are plain Arrays with VecRestore finalizers;
#         fx_arr = lx_arr = fx_bounce = nothing
#   GPU:  fx, lx are device arrays; fx_arr, lx_arr are raw PETSc arrays
#         (needed for restore); fx_bounce is a scratch device array or nothing.

"""
    get_petsc_arrays(petsclib, g_fx, l_x) -> (fx, lx, fx_arr, lx_arr, fx_bounce)

Return arrays for `g_fx` (read-write) and `l_x` (read-only) suitable for
passing to a compute kernel.  Dispatches on the memory location of each Vec
via `_memtype_backend`.

On the pure-CPU path (`HostBackend × HostBackend`) `fx`/`lx` are plain
`Array`s and `fx_arr = lx_arr = fx_bounce = nothing`.  When a GPU backend
extension is loaded and a Vec lives on the device the returned `fx`/`lx` are
device arrays.  An optional bounce buffer `fx_bounce` is allocated when `g_fx`
is host-resident; its contents must be written back by `restore_petsc_arrays`
after the kernel completes.

See also: [`restore_petsc_arrays`](@ref)
"""
function get_petsc_arrays(petsclib, g_fx, l_x)
    T = petsclib.PetscScalar
    fx_arr, fx_mtype = LibPETSc.VecGetArrayAndMemType(petsclib, g_fx)
    lx_arr, lx_mtype = LibPETSc.VecGetArrayReadAndMemType(petsclib, l_x)
    return _get_petsc_arrays_impl(
        petsclib, g_fx, l_x, T, fx_arr, lx_arr,
        _memtype_backend(fx_mtype), _memtype_backend(lx_mtype),
    )
end

# CPU base case: attach VecRestore finalizers and return the arrays directly.
function _get_petsc_arrays_impl(
    petsclib, g_fx, l_x, ::Type, fx_arr, lx_arr, ::HostBackend, ::HostBackend,
)
    finalizer(fx_arr) do a
        LibPETSc.VecRestoreArrayAndMemType(petsclib, g_fx, a)
    end
    finalizer(lx_arr) do a
        LibPETSc.VecRestoreArrayReadAndMemType(petsclib, l_x, a)
    end
    return fx_arr, lx_arr, nothing, nothing, nothing
end

"""
    restore_petsc_arrays(petsclib, g_fx, l_x, fx, lx, fx_arr, lx_arr, fx_bounce)

Restore PETSc Vecs after a kernel launched via [`get_petsc_arrays`](@ref).

Dispatches to `_restore_petsc_arrays_impl`.  On the CPU path (`fx_arr`,
`lx_arr`, `fx_bounce` all `nothing`) this simply finalizes `fx` and `lx`,
triggering the registered `VecRestoreArray*AndMemType` finalizers.  GPU backend
extensions add a `_restore_petsc_arrays_impl` method for their array types.
"""
function restore_petsc_arrays(petsclib, g_fx, l_x, fx, lx, fx_arr, lx_arr, fx_bounce)
    _restore_petsc_arrays_impl(petsclib, g_fx, l_x, fx, lx, fx_arr, lx_arr, fx_bounce)
end

# CPU base case: VecRestore finalizers on fx/lx do the work.
function _restore_petsc_arrays_impl(
    petsclib, g_fx, l_x, fx, lx, ::Nothing, ::Nothing, ::Nothing,
)
    Base.finalize(fx)
    Base.finalize(lx)
end
