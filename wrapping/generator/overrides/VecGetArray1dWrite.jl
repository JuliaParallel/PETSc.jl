# override for VecGetArray1dWrite; C signature: VecGetArray1dWrite(Vec x, PetscInt m, PetscInt mstart, PetscScalar* a[])
"""
	a::PetscArray = VecGetArray1dWrite(petsclib::PetscLibType,x::AbstractPetscVec, m::PetscInt, mstart::PetscInt) 
Returns a pointer to a 1d contiguous array that will contain this
processor's portion of the vector data.  You MUST call `VecRestoreArray1dWrite()`
when you no longer need access to the array.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of two dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)

Output Parameter:
- `a` - location to put pointer to the array

Level: developer

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrays()`, `VecPlaceArray()`,
`VecRestoreArray2d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`, `VecGetArray3d()`, `VecRestoreArray3d()`,
`VecGetArray2d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecGetArray1dWrite"))
"""
function VecGetArray1dWrite(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, mstart::PetscInt) end

@for_petsc function VecGetArray1dWrite(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, mstart::$PetscInt )
	a_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:VecGetArray1dWrite, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, Ref{Ptr{$PetscScalar}}),
               x, m, mstart, a_,
              )

	data_ptr = unsafe_load(a_[])
	mat = unsafe_wrap(Array, data_ptr, m) 
	a = PetscArray(mat,data_ptr) 

	return a
end

