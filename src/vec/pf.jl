"""
    destroy!(pf::LibPETSc.AbstractPF)

Destroy the mathematical function `pf` holds. Safe to call more than once; does
nothing on a borrowed handle: see [`owns`](@ref).

# External Links
$(doc_external("PF/PFDestroy"))
"""
destroy!(pf::AbstractPF{PetscLib}) where {PetscLib} =
    destroy_handle!(LibPETSc.PFDestroy, pf, PetscLib)

"""
    PF{PetscLib}

Handle to a PETSc mathematical function. Only `LibPETSc` works with it, apart
from [`destroy!`](@ref).
"""
LibPETSc.PF
