# override for MatDestroySubMatrices; C signature: MatDestroySubMatrices(PetscInt n, Mat* mat[])
"""
	MatDestroySubMatrices(petsclib::PetscLibType, n::PetscInt, mat::Union{Ptr, AbstractVector{<:AbstractPetscMat}})
Destroys a set of matrices obtained with `MatCreateSubMatrices()`.

Collective

Input Parameters:
- `n`   - the number of local matrices
- `mat` - the vector of matrices `MatCreateSubMatrices()` or `MatCreateSubMatricesMPI()` returned, or the raw pointer to
  PETSc's array of them

PETSc frees its own array of the matrices, so only a vector a PETSc function returned is accepted: any other vector
throws an `ArgumentError`. The matrices in the vector are left with a null pointer.

Level: advanced

See also: `Mat`, `MatCreateSubMatrices()`, `MatDestroyMatrices()`

# External Links
$(_doc_external("Mat/MatDestroySubMatrices"))
"""
function MatDestroySubMatrices(petsclib::PetscLibType, n::Integer, mat::Union{Ptr, AbstractVector{<:AbstractPetscMat}}) end

@for_petsc function MatDestroySubMatrices(petsclib::$UnionPetscLib, n::$PetscInt, mat::Union{Ptr, AbstractVector{<:AbstractPetscMat}})
	if mat isa AbstractVector
		n == length(mat) || throw(DimensionMismatch("n = $n, but the vector holds $(length(mat)) matrices"))
		mat_ = Ref{Ptr{CMat}}(Ptr{CMat}(handle_array(mat, "MatDestroySubMatrices"; take = true)))
	else
		mat_ = Ref{Ptr{CMat}}(mat)
	end

	@chk ccall(
		(:MatDestroySubMatrices, $petsc_library),
		PetscErrorCode,
		($PetscInt, Ptr{Ptr{CMat}}),
		n, mat_,
	)

	mat isa AbstractVector && foreach(m -> m.ptr = C_NULL, mat)
	return nothing
end
