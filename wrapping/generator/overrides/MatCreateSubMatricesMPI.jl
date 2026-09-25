# override for MatCreateSubMatricesMPI; C signature: MatCreateSubMatricesMPI(Mat mat, PetscInt n, const IS irow[], const IS icol[], MatReuse scall, Mat* submat[])
"""
	submat::Vector{PetscMat} = MatCreateSubMatricesMPI(petsclib::PetscLibType, mat::AbstractPetscMat, n::PetscInt, irow::Vector{<:AbstractIS}, icol::Vector{<:AbstractIS}, scall::MatReuse, submat = nothing)
Extracts MPI submatrices across a sub communicator of `mat` (by pairs of `IS` that may live on subcomms).

Collective

Input Parameters:
- `mat`    - the matrix
- `n`      - the number of submatrixes to be extracted
- `irow`   - index set of rows to extract
- `icol`   - index set of columns to extract
- `scall`  - either `MAT_INITIAL_MATRIX` or `MAT_REUSE_MATRIX`
- `submat` - with `MAT_REUSE_MATRIX`, the vector an earlier call returned, which is refilled and returned

Output Parameter:
- `submat` - the submatrices, owned by the caller

PETSc keeps the submatrices in an array it allocated, and `MatDestroySubMatrices()` needs that array back. Release the
submatrices by passing the returned vector to `MatDestroySubMatrices()`, not by destroying them one at a time.

Level: advanced

See also: `Mat`, `PCGASM`, `MatCreateSubMatrices()`, `MatDestroySubMatrices()`, `MatCreateSubMatrix()`, `MatReuse`

# External Links
$(_doc_external("Mat/MatCreateSubMatricesMPI"))
"""
function MatCreateSubMatricesMPI(petsclib::PetscLibType, mat::AbstractPetscMat, n::Integer, irow::Vector{<:AbstractIS}, icol::Vector{<:AbstractIS}, scall::MatReuse, submat = nothing) end

@for_petsc function MatCreateSubMatricesMPI(petsclib::$UnionPetscLib, mat::AbstractPetscMat, n::$PetscInt, irow::Vector{<:AbstractIS}, icol::Vector{<:AbstractIS}, scall::MatReuse, submat = nothing)
	reuse = scall == MAT_REUSE_MATRIX
	reuse && submat === nothing && throw(ArgumentError("MAT_REUSE_MATRIX needs the submatrices an earlier call returned"))
	submat_ = Ref{Ptr{CMat}}(reuse ? Ptr{CMat}(handle_array(submat, "MatCreateSubMatricesMPI")) : C_NULL)

	@chk ccall(
		(:MatCreateSubMatricesMPI, $petsc_library),
		PetscErrorCode,
		(CMat, $PetscInt, Ptr{CIS}, Ptr{CIS}, MatReuse, Ptr{Ptr{CMat}}),
		mat, n, irow, icol, scall, submat_,
	)

	reuse && return submat
	submat_[] == C_NULL && return PetscMat{$PetscLib}[]
	submat = [PetscMat(p, petsclib) for p in unsafe_wrap(Array, submat_[], n; own = false)]
	return record_handle_array!(submat, submat_[])
end
