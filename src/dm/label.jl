"""
    label(dm::AbstractPetscDM, name::AbstractString) -> DMLabel or nothing

Look up the `DMLabel` named `name` on `dm`. Returns `nothing` if `dm` has no
label with that name.

$(doc_borrowed())

# External Links
$(doc_external("DM/DMGetLabel"))
"""
function label end

LibPETSc.@for_petsc function label(
    dm::AbstractPetscDM{$PetscLib},
    name::AbstractString,
)
    found = LibPETSc.DMGetLabel(getlib($PetscLib), dm, String(name))
    return found.ptr == C_NULL ? nothing : found
end

"""
    destroy!(label::DMLabel)

Destroy a `DMLabel` you created. Does nothing on a borrowed handle such as the
one [`label`](@ref) returns.

# External Links
$(doc_external("DM/DMLabelDestroy"))
"""
function destroy!(label::LibPETSc.AbstractDMLabel{PetscLib}) where {PetscLib}
    return destroy_handle!(LibPETSc.DMLabelDestroy, label, PetscLib)
end
