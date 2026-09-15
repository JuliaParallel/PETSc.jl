# override for DMStagVecRestoreArrayRead; C signature: DMStagVecRestoreArrayRead(DM dm, Vec vec, void* array)
"""
	DMStagVecRestoreArrayRead(petsclib::PetscLibType,dm::AbstractPetscDM, vec::AbstractPetscVec, array::PetscArray) 
restore read

Logically Collective

Input Parameters:
- `dm`  - the `DMSTAG` object
- `vec` - the Vec object

Output Parameter:
- `array` - the read-only array

Level: beginner

See also: 
=== 
`DMSTAG`, `DMStagVecGetArrayRead()`, `DMDAVecRestoreArrayRead()`, `DMDAVecRestoreArrayDOFRead()`

# External Links
$(_doc_external("DMStag/DMStagVecRestoreArrayRead"))
"""
function DMStagVecRestoreArrayRead(petsclib::PetscLibType, dm::AbstractPetscDM, vec::AbstractPetscVec, array::PetscArray) end

@for_petsc function DMStagVecRestoreArrayRead(petsclib::$UnionPetscLib, dm::AbstractPetscDM, vec::AbstractPetscVec, array::PetscArray{$PetscScalar, N} ) where N
    # Restore the read-only Vec array that was obtained from VecGetArrayRead
    if array.ptr[] !== nothing
        vec_array = array.ptr[]::Vector{$PetscScalar}
        VecRestoreArrayRead(petsclib, vec, vec_array)
        array.ptr[] = nothing
    end

    return nothing
end

