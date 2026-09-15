# override for VecRestoreArray1dWrite; C signature: VecRestoreArray1dWrite(Vec x, PetscInt m, PetscInt mstart, PetscScalar* a[])
"""
	VecRestoreArray1dWrite(petsclib::PetscLibType,x::AbstractPetscVec, m::PetscInt, mstart::PetscInt, a::PetscArray{PetscScalar, 1}) 
Restores a vector after `VecGetArray1dWrite()` has been called.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of two dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `a`      - location of pointer to array obtained from `VecGetArray1d()`

Level: developer

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `VecRestoreArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray2d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray1dWrite"))
"""
function VecRestoreArray1dWrite(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, mstart::PetscInt, a::PetscArray{PetscScalar, 1}) end

@for_petsc function VecRestoreArray1dWrite(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, mstart::$PetscInt, a::PetscArray{$PetscScalar, 1} )
	if a.ptr[]  != C_NULL 

    @chk ccall(
               (:VecRestoreArray1dWrite, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, Ref{Ptr{$PetscScalar}}),
               x, m, mstart, a.ptr,
              )

		
	else
		error("The input array is already restored")
	end


	return nothing
end

