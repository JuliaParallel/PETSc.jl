# The generated wrappers are in src/autowrapped/.

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

# The REPL shows the handle, not the entries, which are one C call each
Base.show(io::IO, ::MIME"text/plain", v::AbstractPetscVec) = show(io, v)

function vec_from_ptr(petsclib::PetscLib, ptr::CVec, own::Bool) where {PetscLib <: PetscLibType}
    v = LibPETSc.PetscVec{PetscLib}(ptr, petsclib.age; own)
    # Short-circuits on a borrowed handle, which is the hot path: callbacks wrap
    # PETSc-owned vectors on every invocation and never need the communicator.
    own && finalize_if_serial!(v)
    return v
end

"""
    PetscVec(petsclib, ptr::CVec; own::Bool = false)

Wrap a raw PETSc `Vec` handle.

`own = false`, the default, borrows the handle: no
finalizer is attached and `destroy!` is a no-op. With `own = true` a finalizer
is attached, but only on a serial communicator, since `VecDestroy` is collective
and a GC finalizer runs at an arbitrary point.
"""
LibPETSc.PetscVec(petsclib::PetscLibType, ptr::CVec; own::Bool = false) =
    vec_from_ptr(petsclib, ptr, own)

"""
    PetscVec(v::AbstractPetscVec)

The plain `PetscVec` handle behind any high-level vector wrapper, borrowed from
`v`: `destroy!` on it does nothing (see [`owns`](@ref)).

The autowrapped `*AndMemType` routines are typed `x::PetscVec`, while
`AbstractPetscVec` covers every wrapper; this converts transparently.
"""
LibPETSc.PetscVec(v::AbstractPetscVec{PetscLib}) where {PetscLib} =
    LibPETSc.PetscVec{PetscLib}(v.ptr, v.age; own = false)

"""
    PetscVec(petsclib, n::Integer)

A standard, sequentially-stored serial PETSc vector for `petsclib.PetscScalar`
of length `n`.

# External Links
$(doc_external("Vec/VecCreateSeq"))
"""
function LibPETSc.PetscVec(petsclib::PetscLib, n::Integer) where {PetscLib <: PetscLibType}
    comm = MPI.COMM_SELF
    check_initialized(petsclib)
    PetscInt = petsclib.PetscInt
    v = LibPETSc.VecCreateSeq(petsclib, comm, PetscInt(n))
    return own_serial!(v, comm)
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
    return own_serial!(v, comm)
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
LibPETSc.PetscVec(petsclib::PetscLibType, array::Vector; blocksize = 1) =
    LibPETSc.PetscVec(petsclib, MPI.COMM_SELF, array; blocksize)


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
`:mpi`, …), or `nothing` when no type has been set yet.

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
    v = PetscVec(petsclib, comm, array::Vector; blocksize = 1)

Creates a sequential PETSc vector on the communicator `comm`, which must have one
rank, that uses `array` as its
storage and keeps it alive. The element type of `array` must be the library's scalar
type, which throws an `ArgumentError` otherwise: PETSc reads the memory as it is.

# External Links
$(doc_external("Vec/VecCreateSeqWithArray"))
"""
function LibPETSc.PetscVec(petsclib::PetscLibType, comm::MPI.Comm, x::Vector; blocksize = 1)
    check_initialized(petsclib)
    eltype(x) === petsclib.PetscScalar || throw(
        ArgumentError(
            "array has element type $(eltype(x)), but the library uses $(petsclib.PetscScalar)",
        ),
    )
    MPI.Comm_size(comm) == 1 || throw(ArgumentError(
        "a PetscVec on a Julia array is sequential: its communicator must have one rank, " *
        "this one has $(MPI.Comm_size(comm))"))
    PetscInt = petsclib.PetscInt
    v = seq_vec_with_array(petsclib, comm, PetscInt(blocksize), PetscInt(length(x)), x)
    keep_alive!(v, x)
    return own_serial!(v, comm)
end



"""
    ownership_range(vec::AbstractPetscVec)

The range of indices owned by this processor, assuming that the `vec` is laid
out with the first `n1` elements on the first processor, next `n2` elements on
the second, etc. For certain parallel layouts this range may not be well
defined.

The range is **1-based**, always: an index into Julia data is 1-based.

!!! note

    unlike the C function, the range returned is inclusive (`idx_first:idx_last`)

# External Links
$(doc_external("Vec/VecGetOwnershipRange"))
"""
function ownership_range(vec::AbstractPetscVec{PetscLib}) where {PetscLib}
    PetscInt = inttype(PetscLib)
    # The wrapper returns two plain integers, not `Ref`s.
    r_lo, r_hi = LibPETSc.VecGetOwnershipRange(PetscLib, vec)
    return (r_lo + PetscInt(1)):r_hi
end
