# override for VecGetArray2d; C signature: VecGetArray2d(Vec x, PetscInt m, PetscInt n, PetscInt mstart, PetscInt nstart, PetscScalar** a[])
"""
	a::Vector{PetscScalar} = VecGetArray2d(petsclib::PetscLibType,x::AbstractPetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt) 
Returns a pointer to a 2d contiguous array that contains this
processor's portion of the vector data.  You MUST call `VecRestoreArray2d()`
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
$(_doc_external("Vec/VecGetArray2d"))
"""
function VecGetArray2d(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt) end

@for_petsc function VecGetArray2d(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, n::$PetscInt, mstart::$PetscInt, nstart::$PetscInt )

    arr_ptr = Ref{Ptr{Ptr{$PetscScalar}}}()

    @chk ccall(
               (:VecGetArray2d, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{$PetscScalar}}}),
               x, m, n, mstart, nstart, arr_ptr,
              )

    # Assume contiguous storage, use first row pointer
    #mat = unsafe_wrap(Array, data_ptr, (m, n))
    
    # there is a difference in C vs julia storage of arrays
    data_ptr = unsafe_load(arr_ptr[])
    sz = (m,n)
    perm = (2,1)
    mat = unsafe_wrap(Array, data_ptr, sz)
    mat = PermutedDimsArray(mat, perm)
    #arr = PetscArray(mat, a_)



	return PetscArray(mat, arr_ptr)
end

