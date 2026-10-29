"""
    destroy!(ao::LibPETSc.AbstractAO)

Destroy the application ordering `ao` holds. Safe to call more than once; does
nothing on a borrowed handle: see [`owns`](@ref).

# External Links
$(doc_external("AO/AODestroy"))
"""
destroy!(ao::AbstractAO{PetscLib}) where {PetscLib} =
    destroy_handle!(LibPETSc.AODestroy, ao, PetscLib)

"""
    AO{PetscLib}

Handle to a PETSc application ordering. Only `LibPETSc` works with it, apart
from [`destroy!`](@ref).
"""
LibPETSc.AO
