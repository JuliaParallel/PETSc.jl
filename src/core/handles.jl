# Ownership and release rules shared by every handle, and the docstrings of the
# abstract handle types, which src/autowrapped/ defines without one.

"""
    isdestroyable(obj, ::Type{PetscLib})

Whether `obj` still refers to a PETSc object this process is allowed to destroy.

It is not destroyable when the library is finalized, when the pointer is already
null, or when the object predates the current initialize/finalize cycle.
`initialize` and `finalize` both bump `petsclib.age`, and every object records
the age it was created under. `PetscFinalize` frees everything it owns,
including the inner communicator, so calling `xxxDestroy` on an object from an
earlier cycle reaches a communicator that no longer exists and aborts inside
MPI with "Invalid communicator". That happens from a GC finalizer, so it
surfaces at an arbitrary later point rather than where the object was dropped.
"""
function isdestroyable(obj, ::Type{PetscLib}) where {PetscLib}
    isfinalized(PetscLib) && return false
    obj.ptr == C_NULL && return false
    return obj.age == getlib(PetscLib).age
end

"""
    owns(obj)

Whether `destroy!` on `obj` destroys the PETSc object it holds, which is the
wrapper's `own` field.

A constructor returns an owning wrapper. A reader returns a borrowed one: the
object belongs to whatever it was read from, and `destroy!` on the borrowed
wrapper does nothing and leaves it usable. That covers the high-level readers
(`pc(ksp)`, `snes(ts)`, `dm(ksp)`, `solution(ksp)` and the rest), every
`LibPETSc` function with `Get` in its name except the few that hand out a new
reference (`MatGetFactor`, `DMLabelGetStratumIS`, ...: their manual pages say
the caller destroys the result), and the objects a callback receives.
"""
owns(obj) = obj.own

# What every `destroy!` does: release the PETSc object with `destroy_fn` (the
# class's `XDestroy`) if `obj` owns it and it belongs to the current
# initialize/finalize cycle, then null the handle. A borrowed handle is left as is.
function destroy_handle!(destroy_fn, obj, ::Type{PetscLib}) where {PetscLib}
    owns(obj) || return nothing
    isdestroyable(obj, PetscLib) && destroy_fn(PetscLib, obj)
    obj.ptr = C_NULL
    return nothing
end

# ----------------------------------------------------------------------------
# Type names (docs/src/man/naming.md §3.1)
#
# PETSc registers type names at runtime as strings, so they cannot become enums.
# They are `Symbol` at the Julia API and `String` at the C boundary, and the
# conversion happens here, once, where the call meets C.
#
# The generated `*GetType` wrappers answer `""` (and `MatGetType` answers
# `"(not set)"`) when PETSc has no type for the object yet. A reader that can
# decline to answer is documented rather than papered over with a default
# (§15.1), so both spellings come back as `nothing`.
# ----------------------------------------------------------------------------

type_name_symbol(s::AbstractString) =
    (isempty(s) || s == "(not set)") ? nothing : Symbol(s)

"""
    set_type!(obj, type::Symbol)

Set the PETSc implementation `obj` uses, for example `set_type!(ksp, :gmres)`.

Defined for `PetscVec`, `PetscMat`, `KSP`, `PC`, `SNES`, `TS` and `AbstractPetscDM`.
The `Symbol` is converted to a `String` at the C boundary (§3.1).

`set_type!(obj, "gmres")` still works in v0.5 and warns; it is a `MethodError`
in v0.6.
"""
function set_type! end


# Docstrings of the handle types
# ----------------------------------------------------------------------------

# Abstract handle types: the type to annotate with, so that every wrapper of
# that kind of object is accepted. The parameter is the library type.
for (T, object, wrappers) in (
    (:AbstractPetscVec, "vector (`Vec`)", "[`PetscVec`](@ref)"),
    (:AbstractPetscMat, "matrix (`Mat`)", "[`PetscMat`](@ref) and [`MatShell`](@ref)"),
    (:AbstractPetscDM, "DM", "[`DMDA`](@ref), [`DMStag`](@ref), [`DMPlex`](@ref) and [`PetscDM`](@ref)"),
    (:AbstractKSP, "linear solver (`KSP`)", "[`KSP`](@ref)"),
    (:AbstractPC, "preconditioner (`PC`)", "[`PC`](@ref)"),
    (:AbstractSNES, "nonlinear solver (`SNES`)", "[`SNES`](@ref)"),
    (:AbstractTS, "time stepper (`TS`)", "[`TS`](@ref)"),
    (:AbstractPetscOptions, "options database", "[`PetscOptions`](@ref)"),
    (:AbstractIS, "index set (`IS`)", "[`IS`](@ref)"),
    (:AbstractAO, "application ordering (`AO`)", "[`AO`](@ref)"),
    (:AbstractPF, "mathematical function (`PF`)", "[`PF`](@ref)"),
    (:AbstractTao, "optimization solver (`Tao`)", "[`Tao`](@ref)"),
)
    doc = """
        $T{PetscLib}

    Supertype of every handle to a PETSc $object: $wrappers. Annotate an
    argument with it to accept any of them. `PetscLib` is the library type.
    """
    @eval @doc $doc LibPETSc.$T
end
