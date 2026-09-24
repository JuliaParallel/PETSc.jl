import .LibPETSc: AbstractPC, CPC, PC, AbstractIS, CVec

# Custom display for REPL
function Base.show(io::IO, v::AbstractPC{PetscLib}) where {PetscLib}
    if v.ptr == C_NULL
        print(io, "PETSc PC (null pointer)")
    else
        print(io, "PETSc PC object")
    end
    return nothing
end

"""
    destroy!(p::AbstractPC)

Destroy the preconditioner `p` holds. Does nothing on the borrowed handle
[`pc`](@ref) hands back, which its `KSP` destroys: see [`owns`](@ref).

# External Links
$(doc_external("PC/PCDestroy"))
"""
function destroy!(p::AbstractPC{PetscLib}) where {PetscLib}
    owns(p) || return nothing
    if isdestroyable(p, PetscLib)
        LibPETSc.PCDestroy(PetscLib, p)
    end
    p.ptr = C_NULL
    return nothing
end

"""
    pc(ksp::AbstractKSP)

The preconditioner `ksp` applies. Inside the package, bind it as `p = pc(ksp)`,
not `pc = pc(ksp)` (docs/src/man/naming.md §3.2).

$(doc_borrowed())

# External Links
$(doc_external("KSP/KSPGetPC"))
"""
pc(ksp::AbstractKSP{PetscLib}) where {PetscLib} =
    LibPETSc.KSPGetPC(getlib(PetscLib), ksp)

"""
    type_name(p::AbstractPC)

The name PETSc knows this preconditioner by, as a `Symbol` (`:jacobi`, `:ilu`,
`:fieldsplit`, …), or `nothing` when no type has been set yet.

# External Links
$(doc_external("PC/PCGetType"))
"""
type_name(p::AbstractPC{PetscLib}) where {PetscLib} =
    type_name_symbol(LibPETSc.PCGetType(getlib(PetscLib), p))

"""
    set_type!(p::AbstractPC, type::Symbol)

Set the preconditioner, for example `:jacobi`, `:ilu`, `:gamg` or `:fieldsplit`.

# External Links
$(doc_external("PC/PCSetType"))
"""
function set_type!(p::AbstractPC{PetscLib}, type::Symbol) where {PetscLib}
    LibPETSc.PCSetType(getlib(PetscLib), p, String(type))
    return p
end

"""
    set_fieldsplit_is!(p::AbstractPC, name::AbstractString, is::AbstractIS)

Define the split called `name` of a `:fieldsplit` preconditioner as the rows
listed in `is`, which holds global, 0-based indices. Call it once per split, in
the order the splits should have; `name` is also the options prefix of the
split's solver (`-fieldsplit_<name>_ksp_type`).

Throws an `ArgumentError` unless `p` is already of type `:fieldsplit`: PETSc
would otherwise ignore the call without saying so.

# External Links
$(doc_external("PC/PCFieldSplitSetIS"))
"""
function set_fieldsplit_is!(
    p::AbstractPC{PetscLib},
    name::AbstractString,
    is::AbstractIS{PetscLib},
) where {PetscLib}
    t = type_name(p)
    t === :fieldsplit || throw(ArgumentError(
        "set_fieldsplit_is! needs a :fieldsplit preconditioner, got $(repr(t)); " *
        "call set_type!(p, :fieldsplit) first",
    ))
    LibPETSc.PCFieldSplitSetIS(getlib(PetscLib), p, String(name), is)
    return p
end

# ============================================================================
#   Shell preconditioner
# ============================================================================

# The Julia side of a PC, kept with the PETSc object (naming.md §18.3). PETSc
# also gets its address as the shell context, so the trampolines find it
# without a lookup.
mutable struct PCState <: ObjectState
    apply!::Any
    setup!::Any
    alive::Bool
end
PCState() = PCState(nothing, nothing, true)
state_type(::Type{<:PC}) = PCState

# The state of the `:shell` PC `p`, registered as its shell context
function shell_state(p::AbstractPC{PetscLib}, setter::Symbol) where {PetscLib}
    t = type_name(p)
    t === :shell || throw(ArgumentError(
        "$setter needs a :shell preconditioner, got $(repr(t)); " *
        "call set_type!(p, :shell) first",
    ))
    state = object_state!(p)
    LibPETSc.PCShellSetContext(getlib(PetscLib), p, pointer_from_objref(state))
    return state
end

# The PC PETSc is calling back from, borrowed, and its state
function shell_state(pc_ptr::CPC, ::Type{PetscLib}) where {PetscLib}
    petsclib = getlib(PetscLib)
    p = PC{PetscLib}(pc_ptr, petsclib.age; own = false)
    return p, unsafe_pointer_to_objref(LibPETSc.PCShellGetContext(petsclib, p))::PCState
end

"""
    set_shell_apply!(apply!, p::AbstractPC)

Make `apply!` the action of the `:shell` preconditioner `p`. 
It is called as `apply!(y, p, x)` and writes the preconditioned `x` into `y`.
Returns `p`.

Throws an `ArgumentError` unless `p` is already of type `:shell`: PETSc would
otherwise ignore the call without saying so.

The callback comes first (docs/src/man/naming.md §8.1), so `do` block syntax works.

$(doc_callback())

# External Links
$(doc_external("PC/PCShellSetApply"))
"""
function set_shell_apply! end

struct PCShellSetApplyFn{PetscLib} end
function (::PCShellSetApplyFn{PetscLib})(pc_ptr::CPC, x_ptr::CVec, y_ptr::CVec) where {PetscLib}
    return run_callback("shell apply!") do
        p, state = shell_state(pc_ptr, PetscLib)
        state.apply!(PetscVec{PetscLib}(y_ptr, p.age; own = false), p, PetscVec{PetscLib}(x_ptr, p.age; own = false))
    end
end

LibPETSc.@for_petsc function set_shell_apply!(apply!, p::AbstractPC{$PetscLib})
    state = shell_state(p, :set_shell_apply!)
    state.apply! = apply!
    fptr = @cfunction(
        PCShellSetApplyFn{$PetscLib}(),
        LibPETSc.PetscErrorCode,
        (CPC, CVec, CVec)
    )
    LibPETSc.PCShellSetApply($PetscLib, p, fptr)
    return p
end

"""
    set_shell_setup!(setup!, p::AbstractPC)

Make `setup!` the setup step of the `:shell` preconditioner `p`, called as
`setup!(p)` whenever PETSc sets the preconditioner up, for example after the
operator changes. The `:shell` requirement is as for [`set_shell_apply!`](@ref).
Returns `p`.

The callback comes first (docs/src/man/naming.md §8.1), so `do` block syntax works.

$(doc_callback())

# External Links
$(doc_external("PC/PCShellSetSetUp"))
"""
function set_shell_setup! end

struct PCShellSetSetUpFn{PetscLib} end
function (::PCShellSetSetUpFn{PetscLib})(pc_ptr::CPC) where {PetscLib}
    return run_callback("shell setup!") do
        p, state = shell_state(pc_ptr, PetscLib)
        state.setup!(p)
    end
end

LibPETSc.@for_petsc function set_shell_setup!(setup!, p::AbstractPC{$PetscLib})
    state = shell_state(p, :set_shell_setup!)
    state.setup! = setup!
    fptr = @cfunction(PCShellSetSetUpFn{$PetscLib}(), LibPETSc.PetscErrorCode, (CPC,))
    LibPETSc.PCShellSetSetUp($PetscLib, p, fptr)
    return p
end
