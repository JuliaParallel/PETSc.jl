# override for MatCreateSubMatrices; C signature: MatCreateSubMatrices(Mat mat, PetscInt n, const IS irow[], const IS icol[], MatReuse scall, Mat* submat[])
"""
	submat::Vector{PetscMat} = MatCreateSubMatrices(petsclib::PetscLibType, mat::AbstractPetscMat, n::PetscInt, irow::Vector{<:AbstractIS}, icol::Vector{<:AbstractIS}, scall::MatReuse, submat = nothing)
Extracts several submatrices from a matrix.

Collective

Input Parameters:
- `mat`    - the matrix
- `n`      - the number of submatrixes to be extracted (on this processor, may be zero)
- `irow`   - index set of rows to extract
- `icol`   - index set of columns to extract
- `scall`  - either `MAT_INITIAL_MATRIX` or `MAT_REUSE_MATRIX`
- `submat` - with `MAT_REUSE_MATRIX`, the vector an earlier call returned, which is refilled and returned

Output Parameter:
- `submat` - the submatrices, owned by the caller

PETSc keeps the submatrices in an array it allocated, and `MatDestroySubMatrices()` needs that array back. Release the
submatrices by passing the returned vector to `MatDestroySubMatrices()`, not by destroying them one at a time.

Level: advanced

See also: `Mat`, `MatDestroySubMatrices()`, `MatCreateSubMatrix()`, `MatGetRow()`, `MatGetDiagonal()`, `MatReuse`

# External Links
$(_doc_external("Mat/MatCreateSubMatrices"))
"""
function MatCreateSubMatrices(petsclib::PetscLibType, mat::AbstractPetscMat, n::Integer, irow::Vector{<:AbstractIS}, icol::Vector{<:AbstractIS}, scall::MatReuse, submat = nothing) end

@for_petsc function MatCreateSubMatrices(petsclib::$UnionPetscLib, mat::AbstractPetscMat, n::$PetscInt, irow::Vector{<:AbstractIS}, icol::Vector{<:AbstractIS}, scall::MatReuse, submat = nothing)
	reuse = scall == MAT_REUSE_MATRIX
	reuse && submat === nothing && throw(ArgumentError("MAT_REUSE_MATRIX needs the submatrices an earlier call returned"))
	submat_ = Ref{Ptr{CMat}}(reuse ? Ptr{CMat}(handle_array(submat, "MatCreateSubMatrices")) : C_NULL)

	@chk ccall(
		(:MatCreateSubMatrices, $petsc_library),
		PetscErrorCode,
		(CMat, $PetscInt, Ptr{CIS}, Ptr{CIS}, MatReuse, Ptr{Ptr{CMat}}),
		mat, n, irow, icol, scall, submat_,
	)

	reuse && return submat
	submat_[] == C_NULL && return PetscMat{$PetscLib}[]
	submat = [PetscMat(p, petsclib) for p in unsafe_wrap(Array, submat_[], n; own = false)]
	return record_handle_array!(submat, submat_[])
end
