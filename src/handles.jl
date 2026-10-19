# src/handles.jl
# Release for the handles that have no high-level layer of their own: IS, AO, PF
# and Tao. They carry `own` like every handle, so `destroy!` behaves as it does
# on a Vec or a KSP. An IS also has a `length`. The docstrings of the handle
# types, which src/autowrapped/ defines without one, are at the end.

import .LibPETSc: AbstractIS, AbstractAO, AbstractPF, AbstractTao, IS, AO, PF, Tao

"""
    destroy!(is::LibPETSc.AbstractIS)

Destroy the index set `is` holds. Safe to call more than once. Does nothing on
a borrowed index set, such as one `ISColoringGetIS` leaves with its coloring:
see [`owns`](@ref).

# External Links
$(doc_external("IS/ISDestroy"))
"""
function destroy!(is::AbstractIS{PetscLib}) where {PetscLib}
    owns(is) || return nothing
    if isdestroyable(is, PetscLib)
        LibPETSc.ISDestroy(PetscLib, is)
    end
    is.ptr = C_NULL
    return nothing
end

"""
    length(is::LibPETSc.AbstractIS)

The global number of indices in `is`, summed over its ranks, as `length` of a
`PetscVec` is its global length.

# External Links
$(doc_external("IS/ISGetSize"))
"""
Base.length(is::AbstractIS{PetscLib}) where {PetscLib} =
    Int(LibPETSc.ISGetSize(getlib(PetscLib), is))

"""
    destroy!(ao::LibPETSc.AbstractAO)

Destroy the application ordering `ao` holds. Safe to call more than once; does
nothing on a borrowed handle: see [`owns`](@ref).

# External Links
$(doc_external("AO/AODestroy"))
"""
function destroy!(ao::AbstractAO{PetscLib}) where {PetscLib}
    owns(ao) || return nothing
    if isdestroyable(ao, PetscLib)
        LibPETSc.AODestroy(PetscLib, ao)
    end
    ao.ptr = C_NULL
    return nothing
end

"""
    destroy!(pf::LibPETSc.AbstractPF)

Destroy the mathematical function `pf` holds. Safe to call more than once; does
nothing on a borrowed handle: see [`owns`](@ref).

# External Links
$(doc_external("PF/PFDestroy"))
"""
function destroy!(pf::AbstractPF{PetscLib}) where {PetscLib}
    owns(pf) || return nothing
    if isdestroyable(pf, PetscLib)
        LibPETSc.PFDestroy(PetscLib, pf)
    end
    pf.ptr = C_NULL
    return nothing
end

"""
    destroy!(tao::LibPETSc.AbstractTao)

Destroy the optimization solver `tao` holds. Safe to call more than once; does
nothing on a borrowed handle: see [`owns`](@ref).

# External Links
$(doc_external("Tao/TaoDestroy"))
"""
function destroy!(tao::AbstractTao{PetscLib}) where {PetscLib}
    owns(tao) || return nothing
    if isdestroyable(tao, PetscLib)
        LibPETSc.TaoDestroy(PetscLib, tao)
    end
    tao.ptr = C_NULL
    return nothing
end

# Docstrings of the handle types
# ----------------------------------------------------------------------------

# Abstract handle types: the type to annotate with, so that every wrapper of
# that kind of object is accepted. The parameter is the library type.
for (T, object, wrappers) in (
    (:AbstractPetscVec, "vector (`Vec`)", "[`PetscVec`](@ref) and [`VecPtr`](@ref)"),
    (:AbstractPetscMat, "matrix (`Mat`)", "[`PetscMat`](@ref), [`MatPtr`](@ref) and [`MatShell`](@ref)"),
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

"""
    PC{PetscLib}

Handle to a PETSc preconditioner. [`pc`](@ref) hands one back, borrowed from its
KSP; there is no high-level constructor, since a PC always belongs to a solver.
"""
LibPETSc.PC

"""
    IS{PetscLib}

Handle to a PETSc index set. Built through `LibPETSc` (`ISCreateGeneral`,
`ISCreateStride`, …) or from a `DMStag` layout with `IS(dm, …)`. It has a
`length` and a [`destroy!`](@ref); everything else is `LibPETSc`.
"""
LibPETSc.IS

"""
    AO{PetscLib}

Handle to a PETSc application ordering. Only `LibPETSc` works with it, apart
from [`destroy!`](@ref).
"""
LibPETSc.AO

"""
    PF{PetscLib}

Handle to a PETSc mathematical function. Only `LibPETSc` works with it, apart
from [`destroy!`](@ref).
"""
LibPETSc.PF

"""
    Tao{PetscLib}

Handle to a PETSc optimization solver. Only `LibPETSc` works with it, apart
from [`destroy!`](@ref): there is no high-level Tao layer.
"""
LibPETSc.Tao

"""
    PetscDM{PetscLib}

Low-level handle to a PETSc DM of any flavour, as `LibPETSc.DMCreate` and the
other generated functions return it. [`narrow`](@ref) gives the typed handle,
`DMDA`, `DMStag` or `DMPlex`, that the high-level functions dispatch on.
"""
LibPETSc.PetscDM
