# override for VecRestoreArray3d; C signature: VecRestoreArray3d(Vec x, PetscInt m, PetscInt n, PetscInt p, PetscInt mstart, PetscInt nstart, PetscInt pstart, PetscScalar** a[])
"""
	VecRestoreArray3d(petsclib::PetscLibType,x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, a::PetscArray{PetscScalar, 3}) 
Restores a vector after `VecGetArray3d()` has been called.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of three dimensional array
- `n`      - second dimension of the three dimensional array
- `p`      - third dimension of the three dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `nstart` - first index in the second coordinate direction (often 0)
- `pstart` - first index in the third coordinate direction (often 0)
- `a`      - location of pointer to array obtained from VecGetArray3d()

Level: developer

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray3d"))
"""
function VecRestoreArray3d(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, a::PetscArray{PetscScalar, 3}) end

@for_petsc function VecRestoreArray3d(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, n::$PetscInt, p::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, pstart::$PetscInt, a::PetscArray{$PetscScalar, 3} )
	if a.ptr[]  != C_NULL  

    @chk ccall(
               (:VecRestoreArray3d, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{$PetscScalar}}}),
               x, m, n, p, mstart, nstart, pstart, a.ptr,
              )

		
	else
		error("The input array is already restored")
	end


	return nothing
end

