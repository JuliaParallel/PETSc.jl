# these are julia-specific PETSc functions, which makes PETSc a bit easier to use
# a (nearly) full wrapper is in wrapping/

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
    VecPtr(petsclib, ptr::CVec; own::Bool = false)

Container type for a PETSc Vec that is just a raw pointer.

If `own` is `true` a finalizer is set on the vector, but only on a serial
communicator, since `VecDestroy` is collective and a GC finalizer runs at an
arbitrary point. If `own` is `false`, the default, the handle belongs to PETSc and `destroy!`
is a no-op, leaving the wrapper usable.

Calling `VecPtr` is deprecated and goes in v0.6: wrap a raw handle with
[`PetscVec(petsclib, ptr; own)`](@ref LibPETSc.PetscVec(::PetscLibType, ::CVec)).
"""
mutable struct VecPtr{PetscLib} <:
               AbstractPetscVec{PetscLib}
    ptr::CVec
    age::Int
    own::Bool
end
function VecPtr(petsclib::PetscLibType, ptr::CVec; own::Bool = false)
    @warn "VecPtr goes in v0.6: wrap a raw Vec with PetscVec(petsclib, ptr; own)" maxlog = 1
    return vec_from_ptr(petsclib, ptr, own)
end

function vec_from_ptr(petsclib::PetscLib, ptr::CVec, own::Bool) where {PetscLib <: PetscLibType}
    v = VecPtr{PetscLib}(ptr, petsclib.age, own)
    # Short-circuits on a borrowed handle, which is the hot path: callbacks wrap
    # PETSc-owned vectors on every invocation and never need the communicator.
    if own && MPI.Comm_size(LibPETSc.PetscObjectGetComm(getlib(PetscLib), v)) == 1
        finalizer(destroy!, v)
    end
    return v
end
VecPtr(::Type{PetscLib}, x...; kwargs...) where {PetscLib <: PetscLibType} =
    VecPtr(getlib(PetscLib), x...; kwargs...)


"""
    PetscVec(petsclib, ptr::CVec; own::Bool = false)

Wrap a raw PETSc `Vec` handle.

`own = false` is the borrowed case of docs/src/man/naming.md §3.3: no finalizer
is attached and `destroy!` is a no-op. The result is a [`VecPtr`](@ref), the
wrapper type that carries the ownership flag.
"""
LibPETSc.PetscVec(petsclib::PetscLibType, ptr::CVec; own::Bool = false) =
    vec_from_ptr(petsclib, ptr, own)

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
    PetscVec(petsclib, comm, n::Integer; local_size = PETSC_DECIDE)

A standard PETSc vector of global length `n` distributed over `comm`, with
`local_size` entries on this rank. Pass `n = PETSC_DETERMINE` to give only the
local sizes. On a communicator of more than one rank the vector has no
finalizer, since `VecDestroy` is collective: call [`destroy!`](@ref) on every
rank.

# External Links
$(doc_external("Vec/VecCreateMPI"))
"""
function LibPETSc.PetscVec(
    petsclib::PetscLib,
    comm::MPI.Comm,
    n::Integer;
    local_size::Integer = LibPETSc.PETSC_DECIDE,
) where {PetscLib <: PetscLibType}
    check_initialized(petsclib)
    PetscInt = petsclib.PetscInt
    v = LibPETSc.VecCreateMPI(petsclib, comm, PetscInt(local_size), PetscInt(n))
    if MPI.Comm_size(comm) == 1
        finalizer(destroy!, v)
    end
    return v
end


"""
    PetscVec(petsclib, v::Vector)

A standard, sequentially-stored serial PETSc vector, wrapping the Julia vector
`v`.

This reuses the array `v` as storage, and so `v` should not be `resize!`-ed or
otherwise have its length modified while the PETSc object exists. For the same
reason `v` must be a `Vector` of exactly the library's scalar type; anything else
throws an `ArgumentError`. Convert other arrays first: `PetscVec(petsclib, Vector{PetscScalar}(w))`.
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
    v = seq_vec_with_array(petsclib, comm, PetscInt(blocksize), PetscInt(length(array)), array)
    keep_alive!(v, array)
    finalizer(destroy!, v)
    return v
end


# A sequential Vec on the Julia array `array`, without a copy. PETSc starts it at
# object state 0, which is also the state the norm cache slots get when PETSc first
# allocates them, so after one norm every other norm type would read as cached and
# zero. Assembling once moves the Vec past that state; on one process it costs nothing.
function seq_vec_with_array(petsclib, comm, bs, n, array)
    v = LibPETSc.VecCreateSeqWithArray(petsclib, comm, bs, n, array)
    LibPETSc.VecAssemblyBegin(petsclib, v)
    LibPETSc.VecAssemblyEnd(petsclib, v)
    return v
end
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
destroy!(m::AbstractPetscVec{PetscLib}) where {PetscLib} =
    destroy_handle!(LibPETSc.VecDestroy, m, PetscLib)

"""
    v = PetscVec(petsclib, comm, array::Vector)

Creates a sequential PETSc vector on the communicator `comm` that uses `array` as its
storage and keeps it alive. The element type of `array` must be the library's scalar
type, which throws an `ArgumentError` otherwise: PETSc reads the memory as it is.

# External Links
$(doc_external("Vec/VecCreateSeqWithArray"))
"""
function LibPETSc.PetscVec(petsclib::PetscLib, comm, x::Vector) where {PetscLib <: PetscLibType}
    check_initialized(petsclib)
    eltype(x) === petsclib.PetscScalar || throw(
        ArgumentError(
            "array has element type $(eltype(x)), but the library uses $(petsclib.PetscScalar)",
        ),
    )
    PetscInt = petsclib.PetscInt

    v = seq_vec_with_array(petsclib, comm, PetscInt(1), PetscInt(length(x)), x)
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
