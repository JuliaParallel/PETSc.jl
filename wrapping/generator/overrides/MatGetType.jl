# override for MatGetType; C signature: MatGetType(Mat mat, MatType* type)
"""
	type::MatType = MatGetType(petsclib::PetscLibType,mat::AbstractPetscMat) 
Gets the matrix type as a string from the matrix object.

Not Collective

Input Parameter:
- `mat` - the matrix

Output Parameter:
- `type` - name of matrix type

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatType`, `MatSetType()`

# External Links
$(_doc_external("Mat/MatGetType"))
"""
function MatGetType(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatGetType(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	type_ = Ref{MatType}()

    @chk ccall(
               (:MatGetType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatType}),
               mat, type_,
              )

	# Handle NULL type (matrix not yet set up)
	if type_[] == C_NULL
		return "(not set)"
	end
	type = unsafe_string(type_[])

	return type
end

