# override for DMStagVecGetArray; C signature: DMStagVecGetArray(DM dm, Vec vec, void* array)
"""
	array::PetscArray = DMStagVecGetArray(petsclib::PetscLibType,dm::AbstractPetscDM, vec::AbstractPetscVec) 
get access to local array

Logically Collective

Input Parameters:
- `dm`  - the `DMSTAG` object
- `vec` - the `Vec` object

Output Parameter:
- `array` - the array

Level: beginner

Note:
This function returns a (dim+1)-dimensional array for a dim-dimensional `DMSTAG`.

Also: in Julia we use OffsetArrays, such that index 1 (e.g. `array[1]` always refers to the true first point in the domain, 
and `array[0]` would be the ghost point (if available) . 

The first 1-3 dimensions indicate an element in the global
numbering, using the standard C ordering.

The final dimension in this array corresponds to a degree
of freedom with respect to this element, for example corresponding to
the element or one of its neighboring faces, edges, or vertices.

For example, for a 3D `DMSTAG`, indexing is `array[k][j][i][idx]`, where `k` is the
index in the z-direction, `j` is the index in the y-direction, and `i` is the
index in the x-direction.

`idx` is obtained with `DMStagGetLocationSlot()`, since the correct offset
into the (d+1)-dimensional C array for a d-dimensional `DMSTAG` depends on the grid size and the number
of DOF stored at each location.

`DMStagVecRestoreArray()` must be called, once finished with the array

See also: 
=== 
`DMSTAG`, `DMStagVecGetArrayRead()`, `DMStagGetLocationSlot()`, `DMGetLocalVector()`, `DMCreateLocalVector()`, `DMGetGlobalVector()`, `DMCreateGlobalVector()`, `DMDAVecGetArray()`, `DMDAVecGetArrayDOF()`

# External Links
$(_doc_external("DMStag/DMStagVecGetArray"))
"""
function DMStagVecGetArray(petsclib::PetscLibType, dm::AbstractPetscDM, vec::AbstractPetscVec) end

@for_petsc function DMStagVecGetArray(petsclib::$UnionPetscLib, dm::AbstractPetscDM, vec::AbstractPetscVec)
    # DMStagVecGetArray in PETSc's C API returns a multi-dimensional pointer structure
    # that does not map directly to the Vec's contiguous storage. Instead, we use
    # the Vec's data directly and reshape it according to DMStag's layout.
    
    xs,ys,zs,m,n,p = DMStagGetGhostCorners(petsclib, dm)
    dim         = DMGetDimension(petsclib, dm)
    q           = DMStagGetEntriesPerElement(petsclib, dm)

    # Get the underlying Vec data - this needs to be restored later
    vec_array = VecGetArray(petsclib, vec)
    
    # Reshape to (m, n[, p], q) layout where spatial dimensions come first,
    # then DOF index last for convenient slicing.
    # Use PermutedDimsArray to create a VIEW (not a copy) so modifications
    # to the array directly affect the underlying Vec data.
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

