# override for DMStagVecGetValuesStencil; C signature: DMStagVecGetValuesStencil(DM dm, Vec vec, PetscInt n, DMStagStencil* pos, PetscScalar* val)
"""
	val::PetscScalar = DMStagVecGetValuesStencil(petsclib::PetscLibType,dm::AbstractPetscDM, vec::AbstractPetscVec, n::PetscInt, pos::DMStagStencil) 
get vector values using grid indexing

Not Collective

Input Parameters:
- `dm`  - the `DMSTAG` object
- `vec` - the vector object
- `n`   - the number of values to obtain
- `pos` - locations to obtain values from (as an array of `DMStagStencil` values)

Output Parameter:
- `val` - value at the point

Notes:
Accepts stencils which refer to global element numbers, but
only allows access to entries in the local representation (including ghosts).

This approach is not as efficient as getting values directly with `DMStagVecGetArray()`,
which is recommended for matrix-free operators.

Level: advanced

See also: 
=== 
`DMSTAG`, `DMStagStencil`, `DMStagStencilLocation`, `DMStagVecSetValuesStencil()`, `DMStagMatSetValuesStencil()`, `DMStagVecGetArray()`

# External Links
$(_doc_external("DMStag/DMStagVecGetValuesStencil"))
"""
function DMStagVecGetValuesStencil(petsclib::PetscLibType, dm::AbstractPetscDM, vec::AbstractPetscVec, n::PetscInt, pos::Vector{DMStagStencil}) end

@for_petsc function DMStagVecGetValuesStencil(petsclib::$UnionPetscLib, dm::AbstractPetscDM, vec::AbstractPetscVec, n::$PetscInt, pos::Vector{DMStagStencil})
    # Allocate the output array (PETSc writes into this buffer)
    @assert n>0

    vals = Vector{$PetscScalar}(undef, n)

    @chk ccall(
        (:DMStagVecGetValuesStencil, $petsc_library),
        PetscErrorCode,
        (CDM, CVec, $PetscInt, Ptr{DMStagStencil}, Ptr{$PetscScalar}),
        dm, vec, n, pos, vals,
    )

    return vals
end

