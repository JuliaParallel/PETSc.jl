import .LibPETSc: AbstractPC, CPC, PC, AbstractIS

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

Does nothing. The high-level layer only hands out a `PC` borrowed from its
`KSP` (see [`pc`](@ref)), and the `KSP` destroys it. A PC made with
`LibPETSc.PCCreate` is destroyed with `LibPETSc.PCDestroy`.
"""
destroy!(::AbstractPC) = nothing

# Borrowed from its KSP, see `destroy!` above
owns(::AbstractPC) = false

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
    return nothing
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
    return nothing
end
