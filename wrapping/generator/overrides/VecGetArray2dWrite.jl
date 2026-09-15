# override for VecGetArray2dWrite; C signature: VecGetArray2dWrite(Vec x, PetscInt m, PetscInt n, PetscInt mstart, PetscInt nstart, PetscScalar** a[])
"""
	a::Vector{PetscScalar} = VecGetArray2dWrite(petsclib::PetscLibType,x::AbstractPetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt) 
Returns a pointer to a 2d contiguous array that will contain this
processor's portion of the vector data.  You MUST call `VecRestoreArray2dWrite()`
when you no longer need access to the array.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of two dimensional array
- `n`      - second dimension of two dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `nstart` - first index in the second coordinate direction (often 0)

Output Parameter:
- `a` - location to put pointer to the array

Level: developer

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrays()`, `VecPlaceArray()`,
`VecRestoreArray2d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`, `VecGetArray3d()`, `VecRestoreArray3d()`,
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecGetArray2dWrite"))
"""
function VecGetArray2dWrite(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt) end

@for_petsc function VecGetArray2dWrite(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, n::$PetscInt, mstart::$PetscInt, nstart::$PetscInt )

    arr_ptr = Ref{Ptr{Ptr{$PetscScalar}}}()

    @chk ccall(
               (:VecGetArray2dWrite, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{$PetscScalar}}}),
               x, m, n, mstart, nstart, arr_ptr,
            )

    # Assume contiguous storage, use first row pointer
    data_ptr = unsafe_load(arr_ptr[])
    mat = unsafe_wrap(Array, data_ptr, (m, n))

	return PetscArray(mat,arr_ptr)            
end

