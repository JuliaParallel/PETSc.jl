# override for VecRestoreArray4dWrite; C signature: VecRestoreArray4dWrite(Vec x, PetscInt m, PetscInt n, PetscInt p, PetscInt q, PetscInt mstart, PetscInt nstart, PetscInt pstart, PetscInt qstart, PetscScalar** a[])
"""
	VecRestoreArray4dWrite(petsclib::PetscLibType,x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, q::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, qstart::PetscInt, a::PetscArray{PetscScalar, 4}) 
Restores a vector after `VecGetArray4dWrite()` has been called.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of four dimensional array
- `n`      - second dimension of the four dimensional array
- `p`      - third dimension of the four dimensional array
- `q`      - fourth dimension of the four dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `nstart` - first index in the second coordinate direction (often 0)
- `pstart` - first index in the third coordinate direction (often 0)
- `qstart` - first index in the fourth coordinate direction (often 0)
- `a`      - location of pointer to array obtained from `VecGetArray4d()`

Level: developer

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `VecRestoreArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray4dWrite"))
"""
function VecRestoreArray4dWrite(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, q::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, qstart::PetscInt, a::PetscArray{PetscScalar, 4}) end

@for_petsc function VecRestoreArray4dWrite(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, n::$PetscInt, p::$PetscInt, q::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, pstart::$PetscInt, qstart::$PetscInt, a::PetscArray{$PetscScalar, 4} )
	if a.ptr[]  != C_NULL  

    @chk ccall(
               (:VecRestoreArray4dWrite, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}}}),
               x, m, n, p, q, mstart, nstart, pstart, qstart, a.ptr,
              )
	else
		error("The input array is already restored")
	end


	return nothing
end

