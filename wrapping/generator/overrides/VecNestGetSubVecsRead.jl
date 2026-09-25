# override for VecNestGetSubVecsRead; C signature: VecNestGetSubVecsRead(Vec X, PetscInt* N, const Vec* sx[])
"""
	N::PetscInt,sx::Vector{PetscVec} = VecNestGetSubVecsRead(petsclib::PetscLibType, X::AbstractPetscVec)
Access the subvecs of a `VECNEST` vector for read-only access

Logically collective

Input Parameter:
- `X` - nest vector

Output Parameters:
- `N`  - number of nested vecs
- `sx` - array of read-locked vectors, borrowed from `X`

`VecNestRestoreSubVecsRead()` checks that it gets back the array PETSc handed out, so pass it the returned vector.

Level: advanced

See also: `VECNEST`, `Vec`, `VecType`, `VecNestGetSize()`, `VecNestGetSubVec()`, `VecNestRestoreSubVecsRead()`

# External Links
$(_doc_external("Vec/VecNestGetSubVecsRead"))
"""
function VecNestGetSubVecsRead(petsclib::PetscLibType, X::AbstractPetscVec) end

@for_petsc function VecNestGetSubVecsRead(petsclib::$UnionPetscLib, X::AbstractPetscVec)
	N_ = Ref{$PetscInt}()
	sx_ = Ref{Ptr{CVec}}()

	@chk ccall(
		(:VecNestGetSubVecsRead, $petsc_library),
		PetscErrorCode,
		(CVec, Ptr{$PetscInt}, Ptr{Ptr{CVec}}),
		X, N_, sx_,
	)

	N = N_[]
	sx_[] == C_NULL && return N, PetscVec{$PetscLib}[]
	sx = [PetscVec(p, petsclib; own = false) for p in unsafe_wrap(Array, sx_[], N; own = false)]
	return N, record_handle_array!(sx, sx_[])
end
