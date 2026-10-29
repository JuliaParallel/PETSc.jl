"""
    destroy!(is::LibPETSc.AbstractIS)

Destroy the index set `is` holds. Safe to call more than once. Does nothing on
a borrowed index set, such as one `ISColoringGetIS` leaves with its coloring:
see [`owns`](@ref).

# External Links
$(doc_external("IS/ISDestroy"))
"""
destroy!(is::AbstractIS{PetscLib}) where {PetscLib} =
    destroy_handle!(LibPETSc.ISDestroy, is, PetscLib)

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
    IS{PetscLib}

Handle to a PETSc index set. Built through `LibPETSc` (`ISCreateGeneral`,
`ISCreateStride`, …) or from a `DMStag` layout with `IS(dm, …)`. It has a
`length` and a [`destroy!`](@ref); everything else is `LibPETSc`.
"""
LibPETSc.IS
