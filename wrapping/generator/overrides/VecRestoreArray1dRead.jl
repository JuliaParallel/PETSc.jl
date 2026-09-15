# override for VecRestoreArray1dRead; C signature: VecRestoreArray1dRead(Vec x, PetscInt m, PetscInt mstart, PetscScalar* a[])
"""
	VecRestoreArray1dRead(petsclib::PetscLibType,x::AbstractPetscVec, m::PetscInt, mstart::PetscInt, a::PetscArray{PetscScalar, 1}) 
Restores a vector after `VecGetArray1dRead()` has been called.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of two dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `a`      - location of pointer to array obtained from `VecGetArray1dRead()`

Level: developer

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `VecRestoreArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray2d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray1dRead"))
"""
function VecRestoreArray1dRead(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, mstart::PetscInt, a::PetscArray{PetscScalar, 1}) end

@for_petsc function VecRestoreArray1dRead(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, mstart::$PetscInt, a::PetscArray{$PetscScalar, 1} )
	if a.ptr[]  != C_NULL 
        @chk ccall(
                (:VecRestoreArray1dRead, $petsc_library),
                PetscErrorCode,
                (CVec, $PetscInt, $PetscInt, Ptr{Ptr{$PetscScalar}}),
                x, m, mstart, a.ptr,
                )
	else
		error("The input array is already restored")
	end


	return nothing
end

