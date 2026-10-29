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

"""
    create_split_boundary_labels!(dm)

Split the single `"marker"` label on a 2D box mesh into four per-wall labels
(`"markerBottom"`, `"markerRight"`, `"markerTop"`, `"markerLeft"`, marker IDs 1–4)
so that corner DOFs can belong to multiple wall labels simultaneously.

Deprecated: it fits one example's mesh rather than DMPlex in general, so it
leaves PETSc.jl in v0.6. examples/dmplex/ex69.jl has a copy.
"""
function create_split_boundary_labels! end

LibPETSc.@for_petsc function create_split_boundary_labels!(
    dm::AbstractPetscDM{$PetscLib},
)
    # TODO: remove in v0.6 with the other 0.5 deprecations.
    @warn "create_split_boundary_labels! is deprecated and leaves PETSc.jl in v0.6; " *
          "copy it from examples/dmplex/ex69.jl" maxlog = 1
    petsclib = getlib($PetscLib)
    names = ("markerBottom", "markerRight", "markerTop", "markerLeft")
    for (name, id) in zip(names, $PetscInt[1, 2, 3, 4])
        LibPETSc.DMCreateLabel(petsclib, dm, name)
        is = LibPETSc.DMGetStratumIS(petsclib, dm, "marker", id)
        is.ptr == C_NULL && continue
        label = LibPETSc.DMGetLabel(petsclib, dm, name)
        LibPETSc.DMLabelInsertIS(petsclib, label, is, $PetscInt(1))
        LibPETSc.ISDestroy(petsclib, is)
    end
    return dm
end
