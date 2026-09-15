# override for MatGetColumnIJ; C signature: MatGetColumnIJ(Mat mat, PetscInt shift, PetscBool symmetric, PetscBool inodecompressed, PetscInt* n, PetscInt* ia[], PetscInt* ja[], PetscBool* done)
"""
	n::PetscInt,ia::Vector{PetscInt},ja::Vector{PetscInt},done::PetscBool = MatGetColumnIJ(petsclib::PetscLibType,mat::AbstractPetscMat, shift::PetscInt, symmetric::PetscBool, inodecompressed::PetscBool) 
Returns the compressed column storage i and j indices for sequential matrices.

Collective

Input Parameters:
- `mat`             - the matrix
- `shift`           - 1 or zero indicating we want the indices starting at 0 or 1
- `symmetric`       - `PETSC_TRUE` or `PETSC_FALSE` indicating the matrix data structure should be
symmetrized
- `inodecompressed` - `PETSC_TRUE` or `PETSC_FALSE` indicating if the nonzero structure of the
inodes or the nonzero elements is wanted. For `MATBAIJ` matrices the compressed version is
always used.

Output Parameters:
- `n`    - number of columns in the (possibly compressed) matrix
- `ia`   - the column pointers; that is ia[0] = 0, ia[col] = i[col-1] + number of elements in that col of the matrix
- `ja`   - the row indices
- `done` - `PETSC_TRUE` or `PETSC_FALSE`, indicating whether the values have been returned

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatGetRowIJ()`, `MatRestoreColumnIJ()`

# External Links
$(_doc_external("Mat/MatGetColumnIJ"))
"""
function MatGetColumnIJ(petsclib::PetscLibType, mat::AbstractPetscMat, shift::PetscInt, symmetric::PetscBool, inodecompressed::PetscBool) end

@for_petsc function MatGetColumnIJ(petsclib::$UnionPetscLib, mat::AbstractPetscMat, shift::$PetscInt, symmetric::PetscBool, inodecompressed::PetscBool )
	n_ = Ref{$PetscInt}()
	ia_ = Ref{Ptr{$PetscInt}}(C_NULL)
	ja_ = Ref{Ptr{$PetscInt}}(C_NULL)
	done_ = Ref{PetscBool}()

    @chk ccall(
               (:MatGetColumnIJ, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, PetscBool, PetscBool, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}, Ptr{PetscBool}),
               mat, shift, symmetric, inodecompressed, n_, ia_, ja_, done_,
              )

	n = n_[]
	done = done_[]
	# `ia` has n+1 entries (ia[1] == shift), `ja` has ia[n+1] - shift entries; both are
	# PETSc-owned and must be given back with MatRestoreColumnIJ. When PETSc cannot
	# provide them (`done == PETSC_FALSE`) empty arrays are returned.
	if done == PETSC_TRUE && ia_[] != C_NULL
		ia = unsafe_wrap(Array, ia_[], Int(n) + 1; own = false)
		nnz = Int(ia[end]) - Int(shift)
		ja = ja_[] == C_NULL ? $PetscInt[] : unsafe_wrap(Array, ja_[], nnz; own = false)
	else
		ia = $PetscInt[]
		ja = $PetscInt[]
	end

	return n,ia,ja,done
end

