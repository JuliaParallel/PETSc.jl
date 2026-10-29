# override for PetscViewerBinaryReadStringArray; C signature: PetscViewerBinaryReadStringArray(PetscViewer viewer, char*** data)
"""
	data::Vector{String} = PetscViewerBinaryReadStringArray(petsclib::PetscLibType, viewer::PetscViewer) 
Reads a null terminated array of strings from a binary file.

Collective

Input Parameter:
- `viewer` - the binary viewer

Output Parameter:
- `data` - the strings

Level: intermediate

Note:
The strings are copied into Julia `String`s and the C memory is freed here. PETSc allocates the pointers and the
characters as one block, which is released with a single `PetscFree()`.

See also: `PetscViewerBinaryOpen()`, `PetscViewerBinaryRead()`, `PetscViewerBinaryWriteStringArray()`

# External Links
$(_doc_external("Viewer/PetscViewerBinaryReadStringArray"))
"""
function PetscViewerBinaryReadStringArray(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerBinaryReadStringArray(petsclib::$UnionPetscLib, viewer::PetscViewer )
	data_ = Ref{Ptr{Ptr{Cchar}}}(C_NULL)

    @chk ccall(
               (:PetscViewerBinaryReadStringArray, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Ptr{Ptr{Cchar}}}),
               viewer, data_,
              )

	# NULL-terminated array whose strings live in the same allocation
	data = String[]
	if data_[] != C_NULL
		i = 1
		while (str_ = unsafe_load(data_[], i)) != C_NULL
			push!(data, unsafe_string(str_))
			i += 1
		end
		PetscFree(petsclib, data_[])
	end

	return data
end
