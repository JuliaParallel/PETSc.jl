# override for VecNestRestoreSubVecsRead; C signature: VecNestRestoreSubVecsRead(Vec X, PetscInt* N, const Vec* sx[])
"""
	VecNestRestoreSubVecsRead(petsclib::PetscLibType, X::AbstractPetscVec, N::PetscInt, sx::Union{Ptr, AbstractVector{<:AbstractPetscVec}})
Restore access the subvecs of a `VECNEST` vector obtained with `VecNestGetSubVecsRead()`

Logically collective

Input Parameters:
- `X`  - nest vector
- `N`  - number of nested vecs
- `sx` - the vector `VecNestGetSubVecsRead()` returned, or the raw pointer to PETSc's array

PETSc checks that it gets back its own array, so only a vector `VecNestGetSubVecsRead()` returned is accepted: any
other vector throws an `ArgumentError`. The vectors in it are left with a null pointer.

Level: advanced

See also: `VECNEST`, `Vec`, `VecType`, `VecNestGetSize()`, `VecNestGetSubVec()`, `VecNestGetSubVecsRead()`

# External Links
$(_doc_external("Vec/VecNestRestoreSubVecsRead"))
"""
function VecNestRestoreSubVecsRead(petsclib::PetscLibType, X::AbstractPetscVec, N::Integer, sx::Union{Ptr, AbstractVector{<:AbstractPetscVec}}) end

@for_petsc function VecNestRestoreSubVecsRead(petsclib::$UnionPetscLib, X::AbstractPetscVec, N::$PetscInt, sx::Union{Ptr, AbstractVector{<:AbstractPetscVec}})
	if sx isa AbstractVector
		N == length(sx) || throw(DimensionMismatch("N = $N, but the vector holds $(length(sx)) vectors"))
		sx_ = Ref{Ptr{CVec}}(Ptr{CVec}(handle_array(sx, "VecNestRestoreSubVecsRead"; take = true)))
	else
		sx_ = Ref{Ptr{CVec}}(sx)
	end
	N_ = Ref{$PetscInt}(N)

	@chk ccall(
		(:VecNestRestoreSubVecsRead, $petsc_library),
		PetscErrorCode,
		(CVec, Ptr{$PetscInt}, Ptr{Ptr{CVec}}),
		X, N_, sx_,
	)

	sx isa AbstractVector && foreach(v -> v.ptr = C_NULL, sx)
	return nothing
end
