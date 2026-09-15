# override for VecGetArray4dWrite; C signature: VecGetArray4dWrite(Vec x, PetscInt m, PetscInt n, PetscInt p, PetscInt q, PetscInt mstart, PetscInt nstart, PetscInt pstart, PetscInt qstart, PetscScalar** a[])
"""
	a::PetscArray = VecGetArray4dWrite(petsclib::PetscLibType,x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, q::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, qstart::PetscInt) 
Returns a pointer to a 4d contiguous array that will contain this
processor's portion of the vector data.  You MUST call `VecRestoreArray4dWrite()`
when you no longer need access to the array.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of four dimensional array
- `n`      - second dimension of four dimensional array
- `p`      - third dimension of four dimensional array
- `q`      - fourth dimension of four dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `nstart` - first index in the second coordinate direction (often 0)
- `pstart` - first index in the third coordinate direction (often 0)
- `qstart` - first index in the fourth coordinate direction (often 0)

Output Parameter:
- `a` - location to put pointer to the array

Level: developer

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrays()`, `VecPlaceArray()`,
`VecRestoreArray2d()`, `DMDAVecGetarray()`, `DMDAVecRestoreArray()`, `VecGetArray3d()`, `VecRestoreArray3d()`,
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecGetArray4dWrite"))
"""
function VecGetArray4dWrite(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, q::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, qstart::PetscInt) end

@for_petsc function VecGetArray4dWrite(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, n::$PetscInt, p::$PetscInt, q::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, pstart::$PetscInt, qstart::$PetscInt )
	a_ = Ref{Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}}}()

    @chk ccall(
               (:VecGetArray4dWrite, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}}}),
               x, m, n, p, q, mstart, nstart, pstart, qstart, a_,
              )

	data_ptr = unsafe_load(a_[])
	mat = unsafe_wrap(Array, data_ptr, (m,n,p,q)) 
	a = PetscArray(mat,data_ptr) 

	return a
end

