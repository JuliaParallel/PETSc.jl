"""
    ghost_update_begin!(vec::AbstractPetscVec; mode = INSERT_VALUES, scatter = SCATTER_FORWARD)

Begin updating the ghost values of `vec` and return it. See [`ghost_update!`](@ref) for the keywords.

# External Links
$(doc_external("Vec/VecGhostUpdateBegin"))
"""
function ghost_update_begin!(
    vec::AbstractPetscVec{PetscLib};
    mode::InsertMode = INSERT_VALUES,
    scatter::LibPETSc.ScatterMode = SCATTER_FORWARD,
) where {PetscLib}
    LibPETSc.VecGhostUpdateBegin(PetscLib, vec, mode, scatter)
    return vec
end

"""
    ghost_update_end!(vec::AbstractPetscVec; mode = INSERT_VALUES, scatter = SCATTER_FORWARD)

Finish updating the ghost values of `vec` and return it. See [`ghost_update!`](@ref) for the keywords.

# External Links
$(doc_external("Vec/VecGhostUpdateEnd"))
"""
function ghost_update_end!(
    vec::AbstractPetscVec{PetscLib};
    mode::InsertMode = INSERT_VALUES,
    scatter::LibPETSc.ScatterMode = SCATTER_FORWARD,
) where {PetscLib}
    LibPETSc.VecGhostUpdateEnd(PetscLib, vec, mode, scatter)
    return vec
end

"""
    ghost_update!(vec::AbstractPetscVec; mode = INSERT_VALUES, scatter = SCATTER_FORWARD)

Update the ghost values of `vec` and return it: [`ghost_update_begin!`](@ref)
followed by [`ghost_update_end!`](@ref). Use the two halves separately to overlap
the communication with other work.

`scatter = SCATTER_FORWARD` copies the owned values into the ghost points;
`SCATTER_REVERSE` sends the ghost values back to their owners, combined by `mode`
(`INSERT_VALUES` or `ADD_VALUES`).

# External Links
$(doc_external("Vec/VecGhostUpdateBegin"))
$(doc_external("Vec/VecGhostUpdateEnd"))
"""
function ghost_update!(
    vec::AbstractPetscVec;
    mode::InsertMode = INSERT_VALUES,
    scatter::LibPETSc.ScatterMode = SCATTER_FORWARD,
)
    ghost_update_begin!(vec; mode, scatter)
    ghost_update_end!(vec; mode, scatter)
    return vec
end
