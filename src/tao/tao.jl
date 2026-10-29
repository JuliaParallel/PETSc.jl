"""
    destroy!(tao::LibPETSc.AbstractTao)

Destroy the optimization solver `tao` holds. Safe to call more than once; does
nothing on a borrowed handle: see [`owns`](@ref).

# External Links
$(doc_external("Tao/TaoDestroy"))
"""
destroy!(tao::AbstractTao{PetscLib}) where {PetscLib} =
    destroy_handle!(LibPETSc.TaoDestroy, tao, PetscLib)

"""
    Tao{PetscLib}

Handle to a PETSc optimization solver. Only `LibPETSc` works with it, apart
from [`destroy!`](@ref): there is no high-level Tao layer.
"""
LibPETSc.Tao
