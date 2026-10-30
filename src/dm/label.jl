"""
    label(dm::AbstractPetscDM, name::AbstractString) -> Ptr{Cvoid}

Look up a `DMLabel` on `dm` by name and return its raw pointer (`Ptr{Cvoid}`).
Returns `C_NULL` if no label with that name exists.

# External Links
$(doc_external("DM/DMGetLabel"))
"""
function label end

LibPETSc.@for_petsc function label(
    dm::AbstractPetscDM{$PetscLib},
    name::AbstractString,
)
    return LibPETSc.DMGetLabel(getlib($PetscLib), dm, String(name))
end
