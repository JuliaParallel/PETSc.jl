# override for MatRestoreRowIJ; C signature: MatRestoreRowIJ(Mat mat, PetscInt shift, PetscBool symmetric, PetscBool inodecompressed, PetscInt* n, PetscInt* ia[], PetscInt* ja[], PetscBool* done)
"""
	done::PetscBool = MatRestoreRowIJ(petsclib::PetscLibType,mat::AbstractPetscMat, shift::PetscInt, symmetric::PetscBool, inodecompressed::PetscBool, ia::Vector{PetscInt}, ja::Vector{PetscInt}) 
Call after you are completed with the ia,ja indices obtained with `MatGetRowIJ()`.

Collective

Input Parameters:
- `mat`             - the matrix
- `shift`           - 1 or zero indicating we want the indices starting at 0 or 1
- `symmetric`       - `PETSC_TRUE` or `PETSC_FALSE` indicating the matrix data structure should be symmetrized
- `inodecompressed` - `PETSC_TRUE` or `PETSC_FALSE` indicating if the nonzero structure of the
inodes or the nonzero elements is wanted. For `MATBAIJ` matrices the compressed version is
always used.
- `n`               - size of (possibly compressed) matrix
- `ia`              - the row pointers
- `ja`              - the column indices

Output Parameter:
- `done` - `PETSC_TRUE` or `PETSC_FALSE` indicated that the values have been returned

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatGetRowIJ()`, `MatRestoreColumnIJ()`

# External Links
$(_doc_external("Mat/MatRestoreRowIJ"))
"""
function MatRestoreRowIJ(petsclib::PetscLibType, mat::AbstractPetscMat, shift::PetscInt, symmetric::PetscBool, inodecompressed::PetscBool, ia::Vector{PetscInt}, ja::Vector{PetscInt}) end

# `ia`/`ja` are the arrays obtained from MatGetRowIJ (their pointers are handed back to PETSc)
@for_petsc function MatRestoreRowIJ(petsclib::$UnionPetscLib, mat::AbstractPetscMat, shift::$PetscInt, symmetric::PetscBool, inodecompressed::PetscBool, ia::Vector{$PetscInt}, ja::Vector{$PetscInt} )
	n_ = Ref{$PetscInt}()
	ia_ = Ref{Ptr{$PetscInt}}(isempty(ia) ? C_NULL : pointer(ia))
	ja_ = Ref{Ptr{$PetscInt}}(isempty(ja) ? C_NULL : pointer(ja))
	done_ = Ref{PetscBool}()

    @chk ccall(
               (:MatRestoreRowIJ, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, PetscBool, PetscBool, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}, Ptr{PetscBool}),
               mat, shift, symmetric, inodecompressed, n_, ia_, ja_, done_,
              )

	return done_[]
end

