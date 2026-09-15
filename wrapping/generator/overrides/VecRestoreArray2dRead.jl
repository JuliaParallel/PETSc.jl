# override for VecRestoreArray2dRead; C signature: VecRestoreArray2dRead(Vec x, PetscInt m, PetscInt n, PetscInt mstart, PetscInt nstart, PetscScalar** a[])
"""
	VecRestoreArray2dRead(petsclib::PetscLibType,x::AbstractPetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt, a::PetscArray{PetscScalar, 2}) 
Restores a vector after `VecGetArray2dRead()` has been called.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of two dimensional array
- `n`      - second dimension of the two dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `nstart` - first index in the second coordinate direction (often 0)
- `a`      - location of pointer to array obtained from VecGetArray2d()

Level: developer

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `VecRestoreArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray2dRead"))
"""
function VecRestoreArray2dRead(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt, a::PetscArray{PetscScalar, 2}) end

@for_petsc function VecRestoreArray2dRead(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, n::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, a::PetscArray{$PetscScalar, 2} )
	if a.ptr[]  != C_NULL  

        @chk ccall(
                (:VecRestoreArray2dRead, $petsc_library),
                PetscErrorCode,
                (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{$PetscScalar}}}),
                x, m, n, mstart, nstart, a.ptr,
                )

		
	else
		error("The input array is already restored")
	end


	return nothing
end

