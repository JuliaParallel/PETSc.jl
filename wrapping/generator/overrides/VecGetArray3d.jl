# override for VecGetArray3d; C signature: VecGetArray3d(Vec x, PetscInt m, PetscInt n, PetscInt p, PetscInt mstart, PetscInt nstart, PetscInt pstart, PetscScalar** a[])
"""
	a::PetscArray = VecGetArray3d(petsclib::PetscLibType,x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt) 
Returns a pointer to a 3d contiguous array that contains this
processor's portion of the vector data.  You MUST call `VecRestoreArray3d()`
when you no longer need access to the array.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of three dimensional array
- `n`      - second dimension of three dimensional array
- `p`      - third dimension of three dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `nstart` - first index in the second coordinate direction (often 0)
- `pstart` - first index in the third coordinate direction (often 0)

Output Parameter:
- `a` - location to put pointer to the array

Level: developer

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrays()`, `VecPlaceArray()`,
`VecRestoreArray2d()`, `DMDAVecGetarray()`, `DMDAVecRestoreArray()`, `VecRestoreArray3d()`,
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecGetArray3d"))
"""
function VecGetArray3d(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt) end

@for_petsc function VecGetArray3d(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, n::$PetscInt, p::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, pstart::$PetscInt )
	a_ = Ref{Ptr{Ptr{Ptr{$PetscScalar}}}}()

    @chk ccall(
               (:VecGetArray3d, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{Ptr{$PetscScalar}}}}),
               x, m, n, p, mstart, nstart, pstart, a_,
              )

	data_ptr = unsafe_load(a_[])
	mat = unsafe_wrap(Array, data_ptr, (m,n,p)) 
	a = PetscArray(mat,data_ptr) 

	return a
end

