"""
    unsafe_local_array(vec::AbstractVec; read=true, write=true)

Return an `Array{PetscScalar}` containing local portion of the PETSc `vec`

Use `read=false` if the array is write-only; `write=false` if read-only.

The caller hands the array back: `vec` stays checked out, and the next
`VecGetArray` or assembly on it fails, until `Base.finalize` runs on the array
(or the garbage collector gets to it). Prefer [`with_local_array!`](@ref), which
hands it back when its block returns or throws; [`acquire_local_array`](@ref) and
[`release_local_array`](@ref) are the explicit pair with no finalizer.

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
    error("No GPU backend loaded for PetscMemType $MT: load CUDA.jl or AMDGPU.jl, " *
          "or pass memtype = :host for host arrays")
# `Val(mt)` on a value PETSc returned at run time is a dynamic dispatch, so host
# memory, the common case, answers before it: that keeps a checkout inferable and
# allocation free. A device goes on dispatching, so extensions are unaffected.
memtype_backend(mt::LibPETSc.PetscMemType) =
    mt === LibPETSc.PETSC_MEMTYPE_HOST ? nothing : memtype_backend(Val(mt))

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
    acquire_local_array(vec; read, write, memtype = :native) -> (arr, cpu_arr, backend)

Get the local array from `vec` without registering a Julia finalizer.  Returns
the user-visible array, the raw PETSc cpu_arr needed for restore, and the
backend tag that [`release_local_array`](@ref) takes.

With `memtype = :native` the array is where the Vec keeps its data, through
`VecGetArray*AndMemType`. Extensions overload `make_local_array(cpu_arr, backend)`
to wrap the raw array for their device (e.g. `CUDAMemBackend` → `CuArray`).
With `memtype = :host`, or while no GPU extension is loaded, it is the host array
of `VecGetArray*`, for any Vec type: for a device Vec PETSc copies the data to the
host, and back to the device when it next needs it there after a write.
"""
function acquire_local_array(
    vec::AbstractPetscVec; read::Bool, write::Bool, memtype::Symbol = :native,
)
    return array_getter(memtype)(vec; read, write)
end

function acquire_native_array(vec::AbstractPetscVec{PLib}; read::Bool, write::Bool) where {PLib}
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

Restore a local array taken with `acquire_local_array`, after waiting for the device to
finish with it.  Called in `finally` blocks by `with_local_array!`.
"""
function release_local_array(
    cpu_arr, backend::Union{Nothing, AbstractPetscMemBackend}, vec::AbstractPetscVec{PLib};
    read::Bool, write::Bool,
) where {PLib}
    synchronize_device(backend)
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

# ── Host access to any Vec ────────────────────────────────────────────────────
#
# PETSc decides the memory type of `VecGetArray*AndMemType` from the Vec type,
# and nothing tells it before a checkout. So the host path is chosen up front:
# plain `VecGetArray*` gives the host array of any Vec, and the matching
# `VecRestoreArray*` must hand it back, which the `HostCopy` tag records.

# The memory types a loaded GPU extension wraps; each extension adds its own
# when it loads.
const device_memtypes = Set{LibPETSc.PetscMemType}()

# Whether a checkout goes through plain `VecGetArray*`. `:native` does while no
# GPU extension is loaded, since nothing could wrap device memory then.
function host_access(memtype::Symbol)
    memtype === :host && return true
    memtype === :native ||
        throw(ArgumentError("memtype must be :native or :host, got :$memtype"))
    return isempty(device_memtypes)
end

# The checkout function for `memtype`. Callers choose it once and pass it on, so
# that Julia compiles the code after it for each getter, and a checkout stays
# inferable and allocation free.
array_getter(memtype::Symbol) = host_access(memtype) ? acquire_host_array : acquire_native_array

# Backend tag of a checkout through `VecGetArray*`
struct HostCopy end

function acquire_host_array(vec::AbstractPetscVec{PLib}; read::Bool, write::Bool) where {PLib}
    pv = as_petsc_vec(vec)
    arr = if write && read
        LibPETSc.VecGetArray(PLib, pv)
    elseif write
        LibPETSc.VecGetArrayWrite(PLib, pv)
    else
        LibPETSc.VecGetArrayRead(PLib, pv)
    end
    return arr, arr, HostCopy()
end

function release_local_array(
    cpu_arr, ::HostCopy, vec::AbstractPetscVec{PLib}; read::Bool, write::Bool,
) where {PLib}
    pv = as_petsc_vec(vec)
    if write && read
        LibPETSc.VecRestoreArray(PLib, pv, cpu_arr)
    elseif write
        LibPETSc.VecRestoreArrayWrite(PLib, pv, cpu_arr)
    else
        LibPETSc.VecRestoreArrayRead(PLib, pv, cpu_arr)
    end
    return nothing
end

# Wait until the kernels Julia launched on the arrays are done. The GPU packages
# launch on their own streams, not on PETSc's, so PETSc could otherwise read the
# memory before a kernel has written it. Extensions add a method per backend.
synchronize_device(::Union{Nothing, AbstractPetscMemBackend}) = nothing

# The auto-generated *AndMemType wrappers are typed `x::PetscVec`, but
# `AbstractPetscVec` covers every wrapper. Convert transparently.
as_petsc_vec(v::LibPETSc.PetscVec) = v
as_petsc_vec(v::AbstractPetscVec) = LibPETSc.PetscVec(v)

"""
    memtype(vecs::AbstractPetscVec...)

Query the `PetscMemType` of each Vec and return the corresponding array type.
Errors if the Vecs are on heterogeneous devices (different `PetscMemType`
values), since a single `with_local_array!` call cannot handle mixed backends.
Returns `Vector` when all Vecs are host-resident, and also while no GPU extension is
loaded, since [`with_local_array!`](@ref) then gives host arrays for every Vec.

Extensions overload `array_type(::Val{MT})` for a `PetscMemType` enum value
`MT` to register the corresponding array type (e.g. `PETSC_MEMTYPE_DEVICE` →
`CuArray`).
"""
function memtype(vecs::AbstractPetscVec...)
    host_access(:native) && return Vector
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
        memtype::Symbol = :native,
    )
    with_local_array!(::Type{A}, f!, vecs...; read, write, memtype) where {A <: AbstractArray}

Apply `f!` to local array views of `vecs`.

With `memtype = :native` each array is where its Vec keeps the data: a device
array for a Vec on a GPU, with the GPU package loaded. `memtype = :host` gives
host arrays for every Vec, so host code runs on device Vecs; each such checkout of
a device Vec copies its data to the host. Without a GPU package loaded, `:native`
gives host arrays too, so the type follows what is loaded in the session: code that
must stay on the host whatever else is loaded passes `memtype = :host`.

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
    read::Union{Bool, NTuple{N, Bool}} = true,
    write::Union{Bool, NTuple{N, Bool}} = true,
    memtype::Symbol = :native,
) where {N}
    return checkout_arrays(
        f!, nothing, vecs, flags(read, Val(N)), flags(write, Val(N)), array_getter(memtype),
    )
end
with_local_array!(f!, vecs...; kwargs...) = with_local_array!(f!, vecs; kwargs...)

function with_local_array!(
    f!,
    ::Type{A},
    vecs::NTuple{N, AbstractPetscVec};
    read::Union{Bool, NTuple{N, Bool}} = true,
    write::Union{Bool, NTuple{N, Bool}} = true,
    memtype::Symbol = :native,
) where {A <: AbstractArray, N}
    return checkout_arrays(
        f!, A, vecs, flags(read, Val(N)), flags(write, Val(N)), array_getter(memtype),
    )
end
with_local_array!(f!, ::Type{A}, vecs...; kwargs...) where {A <: AbstractArray} =
    with_local_array!(f!, A, vecs; kwargs...)

flags(b::Bool, ::Val{N}) where {N} = ntuple(_ -> b, Val(N))
flags(b::NTuple{N, Bool}, ::Val{N}) where {N} = b

# Check out every vector, call `f!` with the arrays, and hand them all back, 
# also when `f!` throws. `A` is the array type the caller asked for, or `nothing` to
# take whatever the vectors are on. The arrays are acquired without a finalizer,
# because after `Base.finalize` the GC may run the finalizer a second time once
# the array becomes unreachable, and restoring twice crashes PETSc.
function checkout_arrays(
    f!, A, vecs::NTuple{N, AbstractPetscVec}, read, write, acquire::G,
) where {N, G}
    acquired = map(vecs, read, write) do v, r, w
        acquire(v; read = r, write = w)
    end
    try
        # inside the try, so an array of the wrong type is still handed back
        arrays = map(a -> checked_array(a[1], A), acquired)
        return f!(arrays...)
    finally
        foreach(vecs, acquired, read, write) do v, (_, cpu_arr, backend), r, w
            release_local_array(cpu_arr, backend, v; read = r, write = w)
        end
    end
end

checked_array(arr, ::Nothing) = arr
checked_array(arr, ::Type{A}) where {A} =
    arr isa A ? arr : throw(ArgumentError(
        "expected array of type $A but Vec returned $(typeof(arr)). " *
        "Check that the Vec lives on the expected device."
    ))

# ── GPU-aware array access for a pair of Vecs ─────────────────────────────────

"""
    with_local_arrays!(f!, g_fx::AbstractPetscVec, l_x::AbstractPetscVec; memtype = :native)

Call `f!(fx, lx)` with the local array of `g_fx` (written) and of `l_x` (read),
as arrays a compute kernel can use, and hand both back to PETSc afterwards,
also when `f!` throws. Returns what `f!` returns.

Where each array lives is decided for the pair: with a GPU backend extension
loaded and both Vecs on the device, `fx` and `lx` wrap the device memory. When
`g_fx` is on the host, `fx` is a device scratch array whose contents are copied
back into `g_fx` on return, so the kernel can run on the device either way.
Without a GPU extension loaded, or with `memtype = :host`, both are host arrays,
also for device Vecs, whose data PETSc then copies to the host.

```julia
with_local_arrays!(g_fx, l_x) do fx, lx
    residual_kernel!(fx, lx)
end
```

# External Links
$(doc_external("Vec/VecGetArrayAndMemType"))
"""
function with_local_arrays!(
    f!,
    g_fx::AbstractPetscVec{PetscLib},
    l_x::AbstractPetscVec;
    memtype::Symbol = :native,
) where {PetscLib}
    host_access(memtype) &&
        return with_local_array!(f!, (g_fx, l_x); write = (true, false), memtype = :host)
    petsclib = getlib(PetscLib)
    fx, lx, fx_arr, lx_arr, fx_bounce = checkout_local_arrays(petsclib, g_fx, l_x)
    try
        return f!(fx, lx)
    finally
        _restore_local_arrays!(petsclib, g_fx, l_x, fx, lx, fx_arr, lx_arr, fx_bounce)
    end
end

"""
    local_arrays(petsclib, g_fx, l_x) -> (fx, lx, fx_arr, lx_arr, fx_bounce)

Check out the arrays of `g_fx` (read-write) and `l_x` (read-only) for a compute
kernel; [`restore_local_arrays!`](@ref) must hand them back.

Deprecated: use [`with_local_arrays!`](@ref), which restores the arrays also
when the kernel throws.
"""
function local_arrays(petsclib, g_fx, l_x)
    # TODO: remove in v0.6 with the other 0.5 deprecations.
    @warn "local_arrays is deprecated, use with_local_arrays!(f!, g_fx, l_x)" maxlog = 1
    return checkout_local_arrays(petsclib, g_fx, l_x)
end

# Checks out both arrays and picks where they live: `VecGetArrayAndMemType` gives
# each Vec's memory type, and `_local_arrays` dispatches on the pair of backends.
# The base package handles host × host; GPU extensions add the other cases and a
# matching `_restore_local_arrays!`.
#
# Returns (fx, lx, fx_arr, lx_arr, fx_bounce):
#   CPU:  fx, lx are the PETSc arrays; fx_arr = lx_arr = fx_bounce = nothing
#   GPU:  fx, lx are device arrays; fx_arr, lx_arr are the raw PETSc arrays
#         (needed for restore); fx_bounce is a scratch device array or nothing.
function checkout_local_arrays(petsclib, g_fx, l_x)
    T = petsclib.PetscScalar
    fx_arr, fx_mtype = LibPETSc.VecGetArrayAndMemType(petsclib, as_petsc_vec(g_fx))
    lx_arr, lx_mtype = LibPETSc.VecGetArrayReadAndMemType(petsclib, as_petsc_vec(l_x))
    return _local_arrays(
        petsclib, g_fx, l_x, T, fx_arr, lx_arr,
        memtype_backend(fx_mtype), memtype_backend(lx_mtype),
    )
end

# CPU base case: the PETSc arrays themselves, restored explicitly (no finalizers,
# which could restore a second time and crash PETSc).
function _local_arrays(
    petsclib, g_fx, l_x, ::Type, fx_arr, lx_arr, ::Nothing, ::Nothing,
)
    return fx_arr, lx_arr, nothing, nothing, nothing
end

"""
    restore_local_arrays!(petsclib, g_fx, l_x, fx, lx, fx_arr, lx_arr, fx_bounce)

Hand back the arrays [`local_arrays`](@ref) checked out.

Deprecated: use [`with_local_arrays!`](@ref).
"""
function restore_local_arrays!(petsclib, g_fx, l_x, fx, lx, fx_arr, lx_arr, fx_bounce)
    # TODO: remove in v0.6 with the other 0.5 deprecations.
    @warn "restore_local_arrays! is deprecated, use with_local_arrays!(f!, g_fx, l_x)" maxlog = 1
    _restore_local_arrays!(petsclib, g_fx, l_x, fx, lx, fx_arr, lx_arr, fx_bounce)
    return nothing
end

# CPU base case: call VecRestore directly (no finalizers).
function _restore_local_arrays!(
    petsclib, g_fx, l_x, fx, lx, ::Nothing, ::Nothing, ::Nothing,
)
    LibPETSc.VecRestoreArrayAndMemType(petsclib, as_petsc_vec(g_fx), fx)
    LibPETSc.VecRestoreArrayReadAndMemType(petsclib, as_petsc_vec(l_x), lx)
end
