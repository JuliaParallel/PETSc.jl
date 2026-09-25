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
    if !isinitialized(PetscLib)
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

If `own` is `true` a finalizer is set on the vector, but only on a serial
communicator, since `VecDestroy` is collective and a GC finalizer runs at an
arbitrary point. If `own` is `false` the handle belongs to PETSc and `destroy!`
is a no-op, leaving the wrapper usable.
"""
mutable struct VecPtr{PetscLib} <:
               AbstractPetscVec{PetscLib}
    ptr::CVec
    age::Int
    own::Bool
end
function VecPtr(
    petsclib::PetscLib,
    ptr::CVec,
    own,
) where {PetscLib <: PetscLibType}
    v = VecPtr{PetscLib}(ptr, petsclib.age, own)
    # Short-circuits on a borrowed handle, which is the hot path: callbacks wrap
    # PETSc-owned vectors on every invocation and never need the communicator.
    if own && MPI.Comm_size(LibPETSc.PetscObjectGetComm(getlib(PetscLib), v)) == 1
        finalizer(destroy!, v)
    end
    return v
end
VecPtr(::Type{PetscLib}, x...) where {PetscLib <: PetscLibType} = VecPtr(getlib(PetscLib), x...)


"""
    PetscVec(petsclib, ptr::CVec, own::Bool)

Wrap a raw PETSc `Vec` handle.

`own = false` is the borrowed case of docs/src/man/naming.md §3.3: no finalizer
is attached and `destroy!` is a no-op. The result is a [`VecPtr`](@ref), the
wrapper type that carries the ownership flag.
"""
LibPETSc.PetscVec(petsclib::PetscLibType, ptr::CVec, own::Bool) =
    VecPtr(petsclib, ptr, own)

"""
    PetscVec(v::AbstractPetscVec)

The plain `PetscVec` handle behind any high-level vector wrapper, borrowed from
`v`: `destroy!` on it does nothing (see [`owns`](@ref)).

The autowrapped `*AndMemType` routines are typed `x::PetscVec`, while
`AbstractPetscVec` also covers [`VecPtr`](@ref); this converts transparently.
"""
LibPETSc.PetscVec(v::AbstractPetscVec{PetscLib}) where {PetscLib} =
    LibPETSc.PetscVec{PetscLib}(v.ptr, v.age; own = false)

"""
    PetscVec(petsclib, n::Integer)

A standard, sequentially-stored serial PETSc vector for `petsclib.PetscScalar`
of length `n`.

Replaces v0.4's `VecSeq`: construction goes through the type
(docs/src/man/naming.md §5.1).

# External Links
$(doc_external("Vec/VecCreateSeq"))
"""
function LibPETSc.PetscVec(petsclib::PetscLib, n::Integer) where {PetscLib <: PetscLibType}
    comm = MPI.COMM_SELF
    check_initialized(petsclib)
    PetscInt = petsclib.PetscInt
    v = LibPETSc.VecCreateSeq(petsclib, comm, PetscInt(n))
    finalizer(destroy!, v)
    return v
end


"""
    PetscVec(petsclib, v::Vector)

A standard, sequentially-stored serial PETSc vector, wrapping the Julia vector
`v`.

This reuses the array `v` as storage, and so `v` should not be `resize!`-ed or
otherwise have its length modified while the PETSc object exists.
The vector keeps `v` alive, so `PetscVec(petsclib, [1.0, 2.0])` is safe.

This should only be need to be called for more advanced uses, for most simple
usecases, users should be able to pass `Vector`s directly and have the wrapping
performed automatically

# External Links
$(doc_external("Vec/VecCreateSeqWithArray"))
"""
function LibPETSc.PetscVec(
    petsclib::PetscLib,
    array::Vector{PetscScalar};
    blocksize = 1,
) where {PetscLib <: PetscLibType, PetscScalar}
    comm = MPI.COMM_SELF
    check_initialized(petsclib)
    PetscScalar === petsclib.PetscScalar || throw(
        ArgumentError(
            "array has element type $PetscScalar, " *
            "but the library uses $(petsclib.PetscScalar)",
        ),
    )
    PetscInt = petsclib.PetscInt
    v = LibPETSc.VecCreateSeqWithArray(
        petsclib,
        comm,
        PetscInt(blocksize),
        PetscInt(length(array)),
        array,
    )
    keep_alive!(v, array)
    finalizer(destroy!, v)
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

# The library a wrapper carries in its type parameter, as a value. Used where a
# call has to recover `petsclib` from an object rather than take it as an
# argument (§8), for instance `save_vtk(vecs, filename)`.
petsclib_of(::AbstractPetscVec{PetscLib}) where {PetscLib} = getlib(PetscLib)

# Array interface - size and length
Base.size(v::AbstractPetscVec{PetscLib}) where {PetscLib} = LibPETSc.VecGetSize(PetscLib,v)
Base.length(v::AbstractPetscVec{PetscLib}) where {PetscLib} = prod(size(v))
Base.lastindex(v::AbstractPetscVec{PetscLib}) where {PetscLib} = length(v)
Base.similar(v::AbstractPetscVec{PetscLib}) where {PetscLib} =  LibPETSc.VecDuplicate(getlib(PetscLib), v)
"""
    type_name(v::AbstractPetscVec)

The name PETSc knows this vector's implementation by, as a `Symbol` (`:seq`,
`:mpi`, …), or `nothing` when no type has been set yet (docs/src/man/naming.md
§3.1). v0.4 answered with a `String`; that is a break with no shim (§16).

# External Links
$(doc_external("Vec/VecGetType"))
"""
type_name(v::AbstractPetscVec{PetscLib}) where {PetscLib} =
    type_name_symbol(LibPETSc.VecGetType(PetscLib, v))

"""
    set_type!(v::AbstractPetscVec, type::Symbol)

Set the vector implementation, for example `:seq` or `:mpi`.

# External Links
$(doc_external("Vec/VecSetType"))
"""
function set_type!(v::AbstractPetscVec{PetscLib}, type::Symbol) where {PetscLib}
    LibPETSc.VecSetType(getlib(PetscLib), v, String(type))
    return v
end

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
     return v
end

function Base.setindex!(v::AbstractPetscVec{PetscLib}, vals, r::AbstractRange) where {PetscLib} 
    PetscInt = inttype(PetscLib)
    LibPETSc.VecSetValues(PetscLib,v, PetscInt(length(r)), PetscInt.(Vector(r) .- 1), vals, PETSc.INSERT_VALUES)
    return v
end

function Base.fill!(v::AbstractPetscVec{PetscLib}, val) where {PetscLib}
    LibPETSc.VecSet(PetscLib, v, PetscLib.PetscScalar(val))
    return v
end

"""
    copyto!(dst::AbstractPetscVec, src::AbstractPetscVec)

Copy the entries of `src` into `dst` and return `dst`. The two vectors must
have the same global length, which throws a `DimensionMismatch` otherwise, and
the same parallel layout, which PETSc checks.

# External Links
$(doc_external("Vec/VecCopy"))
"""
function Base.copyto!(
    dst::AbstractPetscVec{PetscLib},
    src::AbstractPetscVec{PetscLib},
) where {PetscLib}
    length(dst) == length(src) || throw(
        DimensionMismatch("destination has length $(length(dst)), source has length $(length(src))"),
    )
    LibPETSc.VecCopy(PetscLib, src, dst)
    return dst
end

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
    return A
end



"""
    destroy!(v::AbstractPetscVec)

Destroy a PETSc vector and release its resources.

Safe to call more than once, and safe to reach as a GC finalizer after the
library has been finalized or re-initialized: see [`isdestroyable`](@ref). Does
nothing on a vector that only borrows its handle: see [`owns`](@ref).

# External Links
$(doc_external("Vec/VecDestroy"))
"""
function destroy!(m::AbstractPetscVec{PetscLib}) where {PetscLib}
    owns(m) || return nothing
    if isdestroyable(m, PetscLib)
        LibPETSc.VecDestroy(PetscLib, m)
    end
    m.ptr = C_NULL
    return nothing
end


"""
    unsafe_local_array(vec::AbstractVec; read=true, write=true)

Return an `Array{PetscScalar}` containing local portion of the PETSc `vec`

Use `read=false` if the array is write-only; `write=false` if read-only.

!!! note
    `Base.finalize` should be called on the `Array` before the data can be used.

# External Links
$(doc_external("Vec/VecGetArray"))
$(doc_external("Vec/VecGetArrayWrite"))
$(doc_external("Vec/VecGetArrayRead"))
$(doc_external("Vec/VecRestoreArray"))
$(doc_external("Vec/VecRestoreArrayWrite"))
$(doc_external("Vec/VecRestoreArrayRead"))
"""
function unsafe_local_array(
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
    
    restore = write && read ? LibPETSc.VecRestoreArray :
              write ? LibPETSc.VecRestoreArrayWrite : LibPETSc.VecRestoreArrayRead
    # The finalizer may run after the library was finalized/re-initialized or after the
    # vector's owner destroyed it (borrowed handles): then there is nothing to restore.
    finalizer(v) do v
        (vec.ptr == C_NULL || !isdestroyable(vec, PetscLib)) && return nothing
        try
            restore(PetscLib, vec, v)
        catch err
            err isa LibPETSc.PetscError || rethrow()
        end
        return nothing
    end
    return v
end


# ── Memory backend type hierarchy ─────────────────────────────────────────────
#
# `memtype_backend` maps a `PetscMemType` to a dispatch tag.
# Host memory returns `nothing` (no KA backend — avoids confusion with
# KernelAbstractions.CPU()).  GPU extensions return their own singleton
# (e.g. `CUDAMemBackend`) by overloading `memtype_backend(::Val{MT})`.

"""
    AbstractPetscMemBackend

Abstract supertype for GPU memory backends used by PETSc extensions.
Host memory is represented by `nothing`, not a subtype of this.
GPU extensions define their own concrete subtype (e.g. `CUDAMemBackend`).
"""
abstract type AbstractPetscMemBackend end
const AbstractPETScMemBackend = AbstractPetscMemBackend   # deprecated spelling, remove in v0.6

"""
    memtype_backend(mtype::PetscMemType) → Nothing | AbstractPetscMemBackend

Convert a `PetscMemType` to a dispatch tag.  Returns `nothing` for host memory;
GPU extensions return their own singleton for device memory.
"""
memtype_backend(::Val{LibPETSc.PETSC_MEMTYPE_HOST}) = nothing
memtype_backend(::Val{MT}) where {MT} =
    error("No GPU backend loaded for PetscMemType $MT — load CUDA.jl, AMDGPU.jl, …")
memtype_backend(mt::LibPETSc.PetscMemType) = memtype_backend(Val(mt))

# ── Device-aware local array access ───────────────────────────────────────────
#
# `_unsafe_local_array` is the unified entry point: it calls
# `VecGetArray*AndMemType`, converts the returned `PetscMemType` to a backend
# singleton via `memtype_backend`, and dispatches to `wrap_local_array`.
# GPU extensions add `wrap_local_array` methods for their own backend types.
#
# The typed overload `_unsafe_local_array(::Type{A}, vec; ...)` additionally
# asserts that the returned array is of type `A`, giving a clear error when a
# Vec is on an unexpected device.

function _unsafe_local_array(
    vec::AbstractPetscVec{PetscLib};
    read::Bool = true,
    write::Bool = true,
) where {PetscLib}
    pv = as_petsc_vec(vec)
    if write && read
        cpu_arr, mtype = LibPETSc.VecGetArrayAndMemType(PetscLib, pv)
    elseif write
        cpu_arr, mtype = LibPETSc.VecGetArrayWriteAndMemType(PetscLib, pv)
    else
        cpu_arr, mtype = LibPETSc.VecGetArrayReadAndMemType(PetscLib, pv)
    end
    return wrap_local_array(cpu_arr, memtype_backend(mtype), vec; read, write)
end

function _unsafe_local_array(
    ::Type{A},
    vec::AbstractPetscVec;
    read::Bool = true,
    write::Bool = true,
) where {A <: AbstractArray}
    arr = _unsafe_local_array(vec; read, write)
    arr isa A && return arr
    Base.finalize(arr)   # release the PETSc handle before throwing
    throw(ArgumentError(
        "expected array of type $A but Vec returned $(typeof(arr)). " *
        "Check that the Vec lives on the expected device."
    ))
end

function wrap_local_array(
    cpu_arr, ::Nothing, vec::AbstractPetscVec{PetscLib};
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
function wrap_local_array(cpu_arr, b::AbstractPetscMemBackend, vec; kw...)
    error("wrap_local_array not implemented for backend $(typeof(b)) — " *
          "load the corresponding GPU package (e.g. CUDA.jl)")
end

# ── No-finalizer acquire/release ─────────────────────────────────────────────
#
# `with_local_array!` uses these instead of the finalizer-based `unsafe_local_array`
# to avoid a documented Julia pitfall: after `Base.finalize(x)` is called, if
# `x` later becomes unreachable GC may invoke the finalizer *again*, leading to
# a double VecRestore call on an already-freed Vec (→ SIGSEGV).
# `try/finally` provides deterministic, single-execution cleanup.

"""
    acquire_local_array(vec; read, write) -> (arr, cpu_arr, backend)

Get the local array from `vec` via `VecGetArray*AndMemType` without
registering a Julia finalizer.  Returns the user-visible array, the raw PETSc
cpu_arr needed for restore, and the backend singleton.
Extensions overload `make_local_array(cpu_arr, backend)` to wrap the raw
array for their device (e.g. `CUDAMemBackend` → `CuArray`).
"""
function acquire_local_array(
    vec::AbstractPetscVec{PLib}; read::Bool, write::Bool,
) where {PLib}
    pv = as_petsc_vec(vec)
    cpu_arr, mtype = if write && read
        LibPETSc.VecGetArrayAndMemType(PLib, pv)
    elseif write
        LibPETSc.VecGetArrayWriteAndMemType(PLib, pv)
    else
        LibPETSc.VecGetArrayReadAndMemType(PLib, pv)
    end
    backend = memtype_backend(mtype)
    arr = make_local_array(cpu_arr, backend)
    return arr, cpu_arr, backend
end

# CPU: the raw PETSc array is already a Vector — return it directly.
make_local_array(cpu_arr, ::Nothing) = cpu_arr
make_local_array(_, b::AbstractPetscMemBackend) =
    error("make_local_array not implemented for backend $(typeof(b)) — " *
          "load the corresponding GPU package (e.g. CUDA.jl)")

"""
    release_local_array(cpu_arr, backend, vec; read, write)

Restore a previously acquired local array.  Called in `finally` blocks by
`with_local_array!`.  Extensions overload this for GPU backends.
"""
function release_local_array(
    cpu_arr, ::Nothing, vec::AbstractPetscVec{PLib}; read::Bool, write::Bool,
) where {PLib}
    pv = as_petsc_vec(vec)
    if write && read
        LibPETSc.VecRestoreArrayAndMemType(PLib, pv, cpu_arr)
    elseif write
        LibPETSc.VecRestoreArrayWriteAndMemType(PLib, pv, cpu_arr)
    else
        LibPETSc.VecRestoreArrayReadAndMemType(PLib, pv, cpu_arr)
    end
    return nothing
end
release_local_array(cpu_arr, b::AbstractPetscMemBackend, vec; kw...) =
    error("release_local_array not implemented for backend $(typeof(b)) — " *
          "load the corresponding GPU package (e.g. CUDA.jl)")

# The auto-generated *AndMemType wrappers are typed `x::PetscVec`, but
# `AbstractPetscVec` also includes `VecPtr`.  Convert transparently.
as_petsc_vec(v::LibPETSc.PetscVec) = v
as_petsc_vec(v::AbstractPetscVec) = LibPETSc.PetscVec(v)

"""

Query the `PetscMemType` of each Vec and return the corresponding array type.
Errors if the Vecs are on heterogeneous devices (different `PetscMemType`
values), since a single `with_local_array!` call cannot handle mixed backends.
Returns `Vector` when all Vecs are host-resident.

Extensions overload `array_type(::Val{MT})` for a `PetscMemType` enum value
`MT` to register the corresponding array type (e.g. `PETSC_MEMTYPE_DEVICE` →
`CuArray`).
"""
function memtype(vecs::AbstractPetscVec...)
    mtypes = map(vecs) do v
        PetscLib = typeof(v).parameters[1]
        pv = as_petsc_vec(v)
        arr, mtype = LibPETSc.VecGetArrayReadAndMemType(PetscLib, pv)
        LibPETSc.VecRestoreArrayReadAndMemType(PetscLib, pv, arr)
        mtype
    end
    allequal(mtypes) || throw(ArgumentError(
        "Vecs are on heterogeneous devices: $(unique(mtypes)). " *
        "Use with_local_array!(f!, ::Type{A}, ...) to handle each backend explicitly."
    ))
    return array_type(Val(first(mtypes)))
end

array_type(::Val{LibPETSc.PETSC_MEMTYPE_HOST}) = Vector
array_type(::Val{MT}) where {MT} =
    error("No array type registered for PetscMemType $MT — load the corresponding GPU package (e.g. CUDA.jl)")
# GPU extensions add: array_type(::Val{LibPETSc.PETSC_MEMTYPE_DEVICE}) = CuArray

"""
    with_local_array!(
        f!,
        vecs::NTuple{N, AbstractVec};
        read::Union{Bool, NTuple{N, Bool}} = true,
        write::Union{Bool, NTuple{N, Bool}} = true,
    )
    with_local_array!(::Type{A}, f!, vecs...; read, write) where {A <: AbstractArray}

Apply `f!` to local array views of `vecs`.

The optional `::Type{A}` second argument (after the do-block function) asserts
that every array returned from `VecGetArray*AndMemType` is of type `A`.  Use
it with do-block syntax:

```julia
with_local_array!(Vector, petsc_x; write=true) do x
    x .= 1
end
```

"""
function with_local_array!(
    f!,
    vecs::NTuple{N, AbstractPetscVec};
    kwargs...,
) where {N}
    A = memtype(vecs...)
    return with_local_array!(f!, A, vecs; kwargs...)
end
with_local_array!(f!, vecs...; kwargs...) = with_local_array!(f!, vecs; kwargs...)

function with_local_array!(
    f!,
    ::Type{A},
    vecs::NTuple{N, AbstractPetscVec};
    read::Union{Bool, NTuple{N, Bool}} = true,
    write::Union{Bool, NTuple{N, Bool}} = true,
) where {A <: AbstractArray, N}
    read isa NTuple{N, Bool} || (read = ntuple(_ -> read, N))
    write isa NTuple{N, Bool} || (write = ntuple(_ -> write, N))
    # Acquire all arrays first (no finalizers), then use try/finally for release.
    # This avoids the Julia pitfall where Base.finalize + GC can both run the
    # finalizer if the object becomes unreachable again (double-restore → crash).
    acquired = map(vecs, read, write) do v, r, w
        acquire_local_array(v; read=r, write=w)
    end
    try
        # Type check inside try so finally still releases on mismatch.
        arrays = map(acquired) do (arr, cpu_arr, backend)
            arr isa A || throw(ArgumentError(
                "expected array of type $A but Vec returned $(typeof(arr)). " *
                "Check that the Vec lives on the expected device."
            ))
            arr
        end
        return f!(arrays...)
    finally
        foreach(vecs, acquired, read, write) do v, (_, cpu_arr, backend), r, w
            release_local_array(cpu_arr, backend, v; read=r, write=w)
        end
    end
end
with_local_array!(f!, ::Type{A}, vecs...; kwargs...) where {A <: AbstractArray} =
    with_local_array!(f!, A, vecs; kwargs...)


"""
    ghost_update_begin!(
        vec::AbstractPetscVec,
        insertmode = INSERT_VALUES,
        scattermode = SCATTER_FORWARD,
    )

Begins scattering `vec` to the local or global representations

# External Links
$(doc_external("Vec/VecGhostUpdateBegin"))
"""
function ghost_update_begin!(
    vec::AbstractPetscVec{PetscLib},
    insertmode = INSERT_VALUES,
    scattermode = SCATTER_FORWARD,
) where {PetscLib}
    LibPETSc.VecGhostUpdateBegin(PetscLib, vec, insertmode, scattermode)
    return vec
end

"""
    ghost_update_end!(
        vec::AbstractPetscVec,
        insertmode = INSERT_VALUES,
        scattermode = SCATTER_FORWARD,
    )

Finishes scattering `vec` to the local or global representations

# External Links
$(doc_external("Vec/VecGhostUpdateEnd"))
"""
function ghost_update_end!(
    vec::AbstractPetscVec{PetscLib},
    insertmode = INSERT_VALUES,
    scattermode = SCATTER_FORWARD,
) where {PetscLib}
    LibPETSc.VecGhostUpdateEnd(PetscLib, vec, insertmode, scattermode)
    return vec
end

"""
    ghost_update!(
        vec::AbstractPetscVec,
        insertmode = INSERT_VALUES,
        scattermode = SCATTER_FORWARD,
    )

Finishes scattering `vec` to the local or global representations

# External Links
$(doc_external("Vec/VecGhostUpdateEnd"))
"""
function ghost_update!(
    vec::AbstractPetscVec{PetscLib},
    insertmode = INSERT_VALUES,
    scattermode = SCATTER_FORWARD,
) where {PetscLib}
    ghost_update_begin!(vec,insertmode,scattermode)
    ghost_update_end!(vec,insertmode,scattermode)
    return vec
end

"""
    v = PetscVec(petsclib, comm, array)

Creates a sequential PETSc vector of length `n` given a julia array `array`,
on the communicator `comm`.
The vector uses `array` as its storage and keeps it alive.

# External Links
$(doc_external("Vec/VecCreateSeqWithArray"))
"""
function LibPETSc.PetscVec(petsclib::PetscLib, comm, x::Vector) where {PetscLib <: PetscLibType}
    check_initialized(petsclib)
    PetscInt = petsclib.PetscInt

    v = LibPETSc.VecCreateSeqWithArray(petsclib, comm, PetscInt(1), PetscInt(length(x)), x)
    keep_alive!(v, x)
    finalizer(destroy!, v)

    return v
end



"""
    ownership_range(vec::AbstractPetscVec)

The range of indices owned by this processor, assuming that the `vec` is laid
out with the first `n1` elements on the first processor, next `n2` elements on
the second, etc. For certain parallel layouts this range may not be well
defined.

The range is **1-based**, always: an index into Julia data is 1-based
(docs/src/man/naming.md §12.1). v0.4 took `base_one::Bool` positionally and
made the convention a runtime choice; `ownership_range(v, false)` warns in
v0.5 and is a `MethodError` in v0.6.

!!! note

    unlike the C function, the range returned is inclusive (`idx_first:idx_last`)

# External Links
$(doc_external("Vec/VecGetOwnershipRange"))
"""
function ownership_range(vec::AbstractPetscVec{PetscLib}) where {PetscLib}
    PetscInt = PetscLib.PetscInt
    # The wrapper returns two plain integers, not `Ref`s.
    r_lo, r_hi = LibPETSc.VecGetOwnershipRange(PetscLib, vec)
    return (r_lo + PetscInt(1)):r_hi
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
# `local_arrays` calls `VecGetArrayAndMemType` on both Vecs, converts the
# returned `PetscMemType` values to backend singletons, and dispatches to
# `_local_arrays`.  The base package handles the pure-CPU case
# (host × host (both backends nothing)).  GPU extensions add `_local_arrays`
# methods for their backend combinations and a matching
# `_restore_local_arrays!` method dispatched by `restore_local_arrays!`.
#
# Return tuple:  (fx, lx, fx_arr, lx_arr, fx_bounce)
#   CPU:  fx, lx are plain Arrays with VecRestore finalizers;
#         fx_arr = lx_arr = fx_bounce = nothing
#   GPU:  fx, lx are device arrays; fx_arr, lx_arr are raw PETSc arrays
#         (needed for restore); fx_bounce is a scratch device array or nothing.

"""
    local_arrays(petsclib, g_fx, l_x) -> (fx, lx, fx_arr, lx_arr, fx_bounce)

Return arrays for `g_fx` (read-write) and `l_x` (read-only) suitable for
passing to a compute kernel.  Dispatches on the memory location of each Vec
via `memtype_backend`.

On the pure-CPU path (`host × host (both backends nothing)`) `fx`/`lx` are plain
`Array`s and `fx_arr = lx_arr = fx_bounce = nothing`.  When a GPU backend
extension is loaded and a Vec lives on the device the returned `fx`/`lx` are
device arrays.  An optional bounce buffer `fx_bounce` is allocated when `g_fx`
is host-resident; its contents must be written back by `restore_local_arrays!`
after the kernel completes.

See also: [`restore_local_arrays!`](@ref)
"""
function local_arrays(petsclib, g_fx, l_x)
    T = petsclib.PetscScalar
    fx_arr, fx_mtype = LibPETSc.VecGetArrayAndMemType(petsclib, as_petsc_vec(g_fx))
    lx_arr, lx_mtype = LibPETSc.VecGetArrayReadAndMemType(petsclib, as_petsc_vec(l_x))
    return _local_arrays(
        petsclib, g_fx, l_x, T, fx_arr, lx_arr,
        memtype_backend(fx_mtype), memtype_backend(lx_mtype),
    )
end

# CPU base case: return arrays directly. restore_local_arrays! calls VecRestore
# explicitly — no finalizers to avoid the double-finalization crash.
function _local_arrays(
    petsclib, g_fx, l_x, ::Type, fx_arr, lx_arr, ::Nothing, ::Nothing,
)
    return fx_arr, lx_arr, nothing, nothing, nothing
end

"""
    restore_local_arrays!(petsclib, g_fx, l_x, fx, lx, fx_arr, lx_arr, fx_bounce)

Restore PETSc Vecs after a kernel launched via [`local_arrays`](@ref).

Dispatches to `_restore_local_arrays!`.  On the CPU path (`fx_arr`,
`lx_arr`, `fx_bounce` all `nothing`) this simply finalizes `fx` and `lx`,
triggering the registered `VecRestoreArray*AndMemType` finalizers.  GPU backend
extensions add a `_restore_local_arrays!` method for their array types.
"""
function restore_local_arrays!(petsclib, g_fx, l_x, fx, lx, fx_arr, lx_arr, fx_bounce)
    _restore_local_arrays!(petsclib, g_fx, l_x, fx, lx, fx_arr, lx_arr, fx_bounce)
end

# CPU base case: call VecRestore directly (no finalizers).
function _restore_local_arrays!(
    petsclib, g_fx, l_x, fx, lx, ::Nothing, ::Nothing, ::Nothing,
)
    LibPETSc.VecRestoreArrayAndMemType(petsclib, as_petsc_vec(g_fx), fx)
    LibPETSc.VecRestoreArrayReadAndMemType(petsclib, as_petsc_vec(l_x), lx)
end
