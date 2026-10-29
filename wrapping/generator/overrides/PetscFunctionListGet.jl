# override for PetscFunctionListGet; C signature: PetscFunctionListGet(PetscFunctionList list, const char*** array, int* n)
"""
	array::Vector{String} = PetscFunctionListGet(petsclib::PetscLibType, list::PetscFunctionList) 
Gets the names of the entries in a `PetscFunctionList`, this is used
by help etc.

Not Collective, No Fortran Support

Input Parameter:
- `list` - list of types

Output Parameter:
- `array` - the names

Level: developer

Note:
The names are copied into Julia `String`s. The C array that held them is freed here; the names themselves belong to the list.

See also: `PetscFunctionListAdd()`, `PetscFunctionList`

# External Links
$(_doc_external("Sys/PetscFunctionListGet"))
"""
function PetscFunctionListGet(petsclib::PetscLibType, list::PetscFunctionList) end

@for_petsc function PetscFunctionListGet(petsclib::$UnionPetscLib, list::PetscFunctionList )
	array_ = Ref{Ptr{Ptr{Cchar}}}(C_NULL)
	n_ = Ref{Cint}()

    @chk ccall(
               (:PetscFunctionListGet, $petsc_library),
               PetscErrorCode,
               (PetscFunctionList, Ptr{Ptr{Ptr{Cchar}}}, Ptr{Cint}),
               list, array_, n_,
              )

	# the array is the caller's to free, the names belong to the list
	array = String[unsafe_string(unsafe_load(array_[], i)) for i in 1:(array_[] == C_NULL ? 0 : n_[])]
	array_[] == C_NULL || PetscFree(petsclib, array_[])

	return array
end
