# override for DMStagVecGetArrayRead; C signature: DMStagVecGetArrayRead(DM dm, Vec vec, void* array)
"""
	array::PetscArray = DMStagVecGetArrayRead(petsclib::PetscLibType,dm::AbstractPetscDM, vec::AbstractPetscVec) 
get read

Logically Collective

See the man page for `DMStagVecGetArray()` for more information.

Input Parameters:
- `dm`  - the `DMSTAG` object
- `vec` - the `Vec` object

Output Parameter:
- `array` - the read-only array

Level: beginner

Note:
`DMStagVecRestoreArrayRead()` must be called, once finished with the array

See also: 
=== 
`DMSTAG`, `DMStagGetLocationSlot()`, `DMGetLocalVector()`, `DMCreateLocalVector()`, `DMGetGlobalVector()`, `DMCreateGlobalVector()`, `DMDAVecGetArrayRead()`, `DMDAVecGetArrayDOFRead()`

# External Links
$(_doc_external("DMStag/DMStagVecGetArrayRead"))
"""
function DMStagVecGetArrayRead(petsclib::PetscLibType, dm::AbstractPetscDM, vec::AbstractPetscVec) end

@for_petsc function DMStagVecGetArrayRead(petsclib::$UnionPetscLib, dm::AbstractPetscDM, vec::AbstractPetscVec)
    # DMStagVecGetArrayRead in PETSc's C API returns a multi-dimensional pointer structure
    # that does not map directly to the Vec's contiguous storage. Instead, we use
    # the Vec's data directly (read-only) and reshape it according to DMStag's layout.
    
    xs,ys,zs,m,n,p = DMStagGetGhostCorners(petsclib, dm)
    dim         = DMGetDimension(petsclib, dm)
    q           = DMStagGetEntriesPerElement(petsclib, dm)

    # Get the underlying Vec data (read-only) - this needs to be restored later
    vec_array = VecGetArrayRead(petsclib, vec)
    
    # Reshape to (m, n[, p], q) layout where spatial dimensions come first,
    # then DOF index last for convenient slicing.
    # Use PermutedDimsArray to create a VIEW (not a copy) for read-only access.
    if dim==1
        mat = reshape(vec_array, (m, q))
        mat_oa = OffsetArray(mat, xs+1:xs+m, 1:q)
    elseif dim==2
        mat = PermutedDimsArray(reshape(vec_array, (q, m, n)), (2, 3, 1))
        mat_oa = OffsetArray(mat, xs+1:xs+m, ys+1:ys+n, 1:q)
    elseif dim==3
        mat = PermutedDimsArray(reshape(vec_array, (q, m, n, p)), (2, 3, 4, 1))
        mat_oa = OffsetArray(mat, xs+1:xs+m, ys+1:ys+n, zs+1:zs+p, 1:q)
    else
        error("Unsupported dimension: $dim")
    end

    # Create a PetscArray wrapper
    # Store the original vec_array in ptr field (as Ref{Any}) so we can restore it later
    arr = PetscArray(mat_oa, Ref{Any}(vec_array))

    return arr
end

