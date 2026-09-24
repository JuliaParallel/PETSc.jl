# ============================================================================
#   Julia state of a PETSc object
# ============================================================================
#
# A wrapper can be rebuilt on every read (`pc(ksp)`, `snes(ts)`), so nothing
# but `ptr`, `age` and `own` can live on it. Everything else Julia keeps about a
# PETSc object (callbacks, the user context, the options applied at `solve!`)
# goes in one state object per PETSc object, attached to it through a
# PetscContainer: PETSc destroys the container with the object, and the
# container's destroy function marks the state dead. `object_states` roots
# every state for the GC. Dead states are dropped on the next attach, never in
# the destroy function, which can run from a GC finalizer; `finalize` drops a
# library's states once none of its objects can call back.

# Concrete subtypes are mutable, hold the Julia side of one kind of object, have
# a zero-argument constructor giving the defaults, and a field `alive::Bool`
abstract type ObjectState end

# The state type of a wrapper type, defined next to each wrapper
function state_type end

# PetscLib => (PETSc object address => its state)
const object_states = Dict{Any, Dict{Ptr{Cvoid}, ObjectState}}()
const object_states_lock = ReentrantLock()

# PetscCtxDestroyFn: PETSc passes the address of the container's pointer
function object_state_destroy(ctx::Ptr{Ptr{Cvoid}})
    state = unsafe_pointer_to_objref(unsafe_load(ctx))::ObjectState
    state.alive = false
    return LibPETSc.PetscErrorCode(0)
end

"""
    object_state(obj)

The live state of the PETSc object `obj` points at, or `nothing` if it has none.
"""
function object_state(obj::T) where {T}
    S = state_type(T)
    ptr = Ptr{Cvoid}(obj.ptr)
    ptr == C_NULL && return nothing
    state = lock(object_states_lock) do
        states = get(object_states, library_type(T), nothing)
        states === nothing ? nothing : get(states, ptr, nothing)
    end
    return (state isa S && state.alive) ? state : nothing
end

"""
    object_state!(obj)

The state of the PETSc object `obj` points at, created and attached to it on
first use.
"""
function object_state!(obj::T) where {T}
    state = object_state(obj)
    state === nothing || return state
    ptr = Ptr{Cvoid}(obj.ptr)
    ptr == C_NULL && throw(ArgumentError("$(nameof(T)) has no PETSc object (null pointer)"))
    PetscLib = library_type(T)
    petsclib = getlib(PetscLib)
    state = state_type(T)()
    lock(object_states_lock) do
        states = get!(Dict{Ptr{Cvoid}, ObjectState}, object_states, PetscLib)
        filter!(kv -> kv.second.alive, states)
        states[ptr] = state
    end
    container = LibPETSc.PetscContainerCreate(petsclib, LibPETSc.PetscObjectGetComm(petsclib, ptr))
    LibPETSc.PetscContainerSetPointer(petsclib, container, pointer_from_objref(state))
    destroy_fn = @cfunction(object_state_destroy, LibPETSc.PetscErrorCode, (Ptr{Ptr{Cvoid}},))
    LibPETSc.PetscContainerSetCtxDestroy(petsclib, container, destroy_fn)
    LibPETSc.PetscObjectCompose(petsclib, ptr, "PETSc.jl:state", container)
    LibPETSc.PetscContainerDestroy(petsclib, container)   # the object now holds the only reference
    return state
end

# The context pointer handed to PETSc with a callback: the state, never a wrapper
state_pointer(obj) = pointer_from_objref(object_state!(obj))

# Called by `finalize` after PetscFinalize, when no object of `PetscLib` can call back
function drop_object_states!(::Type{PetscLib}) where {PetscLib}
    lock(object_states_lock) do
        delete!(object_states, PetscLib)
    end
    return nothing
end

library_type(::Type{<:LibPETSc.AbstractKSP{PetscLib}})      where {PetscLib} = PetscLib
library_type(::Type{<:LibPETSc.AbstractSNES{PetscLib}})     where {PetscLib} = PetscLib
library_type(::Type{<:LibPETSc.AbstractTS{PetscLib}})       where {PetscLib} = PetscLib
library_type(::Type{<:LibPETSc.AbstractPC{PetscLib}})       where {PetscLib} = PetscLib
library_type(::Type{<:LibPETSc.AbstractPetscMat{PetscLib}}) where {PetscLib} = PetscLib
library_type(::Type{<:LibPETSc.AbstractPetscVec{PetscLib}}) where {PetscLib} = PetscLib

# Julia arrays a Vec or Mat uses as its storage without copying them
# (`VecCreateSeqWithArray`, `MatCreateSeqAIJWithArrays`). Kept with the PETSc
# object, so they live exactly as long as it does, whoever destroys it.
mutable struct WrappedArrays <: ObjectState
    arrays::Any
    alive::Bool
end
WrappedArrays() = WrappedArrays(nothing, true)
state_type(::Type{<:LibPETSc.AbstractPetscVec}) = WrappedArrays
state_type(::Type{<:LibPETSc.AbstractPetscMat}) = WrappedArrays

# Keep `arrays` alive for as long as the PETSc object `obj` holds exists
function keep_alive!(obj, arrays)
    object_state!(obj).arrays = arrays
    return obj
end

# The wrapper fields that moved into the state stay readable and
# writable as properties of the wrapper
# ----------------------------------------------------------------------------

const wrapper_fields = (:ptr, :age, :own)

function forward_getproperty(obj::T, name::Symbol) where {T}
    state = object_state(obj)
    return getfield(state === nothing ? state_type(T)() : state, name)
end

function forward_setproperty!(obj, name::Symbol, value)
    Base.setfield!(object_state!(obj), name, value)
    return value
end

forward_propertynames(::Type{T}) where {T} =
    (wrapper_fields..., filter(!=(:alive), fieldnames(state_type(T)))...)

# ============================================================================
#   Running Julia callbacks from PETSc
# ============================================================================

# "error in library called by PETSc", from `petscsystypes.h`.
const _PETSC_ERR_LIB = LibPETSc.PetscErrorCode(76)

# The task-local slot where a callback leaves its exception for the high-level
# call that started the work
const CALLBACK_ERROR = :PETSc_callback_error

# A callback's return value is ignored. A nonzero `Integer` is still the PETSc
# error code 0.5.0 read it as, and warns once per callback; v0.6 ignores it.
function errorcode(r, name)
    (r isa Integer && !iszero(r)) || return LibPETSc.PetscErrorCode(0)
    @warn "PETSc.jl: the $name callback returned $r, which fails the call as a PETSc error code. " *
          "From v0.6 a callback's return value is ignored: throw an exception to report a failure." maxlog = 1 _id = Symbol(name)
    return LibPETSc.PetscErrorCode(r)
end

# Call `f` for PETSc and turn the outcome into a PETSc error code.
#
# A Julia exception must not cross the `@cfunction` boundary: PETSc is C and
# cannot unwind a Julia frame. The exception is left for the high-level call
# that started the work (`capture_callback_errors`), which rethrows it once
# PETSc has unwound its own stack. Without one, it is logged, and the
# `LibPETSc` call raises `PetscError`.
function run_callback(f, name)
    try
        return errorcode(f(), name)
    catch e
        bt = catch_backtrace()
        # Nothing here may throw either
        try
            slot = get(task_local_storage(), CALLBACK_ERROR, nothing)
            if slot isa Base.RefValue{Any} && slot[] === nothing
                slot[] = (e, bt, name)
            else
                @error "PETSc.jl: the $name callback failed" exception = (e, bt)
            end
        catch
            Core.println("PETSc.jl: the ", name, " callback failed")
        end
        return _PETSC_ERR_LIB
    end
end

"""
    capture_callback_errors(f)

Run `f`, a high-level call that can make PETSc call back into Julia, and
rethrow the first exception a callback threw, as itself rather than as the
`PetscError` PETSc reports for it. The callback's own frames are logged at
debug level.
"""
function capture_callback_errors(f)
    tls = task_local_storage()
    outer = get(tls, CALLBACK_ERROR, nothing)
    slot = Ref{Any}(nothing)
    tls[CALLBACK_ERROR] = slot
    local result
    try
        result = f()
    catch
        slot[] === nothing && rethrow()
        rethrow_callback_error(slot[])
    finally
        outer === nothing ? delete!(tls, CALLBACK_ERROR) : (tls[CALLBACK_ERROR] = outer)
    end
    # PETSc absorbed the failure, but the callback still threw
    slot[] === nothing || rethrow_callback_error(slot[])
    return result
end

function rethrow_callback_error((e, bt, name))
    @debug "PETSc.jl: the $name callback threw; its frames:" exception = (e, bt)
    throw(e)
end
