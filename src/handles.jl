# src/handles.jl
# Release for the handles that have no high-level layer of their own: IS, AO, PF
# and Tao. They carry `own` like every handle, so `destroy!` behaves as it does
# on a Vec or a KSP.

import .LibPETSc: AbstractIS, AbstractAO, AbstractPF, AbstractTao

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
