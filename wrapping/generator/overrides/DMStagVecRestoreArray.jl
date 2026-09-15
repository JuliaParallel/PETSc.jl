# override for DMStagVecRestoreArray; C signature: DMStagVecRestoreArray(DM dm, Vec vec, void* array)
"""
	DMStagVecRestoreArray(petsclib::PetscLibType,dm::AbstractPetscDM, vec::AbstractPetscVec, array::PetscArray) 
restore access to a raw array

Logically Collective

Input Parameters:
- `dm`  - the `DMSTAG` object
- `vec` - the `Vec` object

Output Parameter:
- `array` - the array

Level: beginner

See also: 
=== 
`DMSTAG`, `DMStagVecGetArray()`, `DMDAVecRestoreArray()`, `DMDAVecRestoreArrayDOF()`

# External Links
$(_doc_external("DMStag/DMStagVecRestoreArray"))
"""
function DMStagVecRestoreArray(petsclib::PetscLibType, dm::AbstractPetscDM, vec::AbstractPetscVec, array::PetscArray) end

@for_petsc function DMStagVecRestoreArray(petsclib::$UnionPetscLib, dm::AbstractPetscDM, vec::AbstractPetscVec, array::PetscArray{$PetscScalar, N} ) where N
    # Since we're using PermutedDimsArray (a view), modifications to the array
    # automatically affect the underlying vec_array. We just need to restore it.
    if !isnothing(array.ptr) && array.ptr[] !== nothing
        vec_array = array.ptr[]
        VecRestoreArray(petsclib, vec, vec_array)
        # Null out the reference to prevent double-restore
        array.ptr[] = nothing
    end
    
    return nothing
end

