"""
	A_loc::PetscMat = MatAIJGetLocalMat(petsclib::PetscLibType,A::AbstractPetscMat) 
Creates a `MATSEQAIJ` from a `MATAIJ` matrix.

Not Collective

Input Parameter:
- `A` - the matrix

Output Parameter:
- `A_loc` - the local sequential matrix generated

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatMPIAIJGetLocalMat()`

# External Links
$(_doc_external("Mat/MatAIJGetLocalMat"))
"""
function MatAIJGetLocalMat(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatAIJGetLocalMat(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	A_loc_ = Ref{CMat}()

    @chk ccall(
               (:MatAIJGetLocalMat, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, A_loc_,
              )

	A_loc = PetscMat(A_loc_[], petsclib)

	return A_loc
end 

"""
	MatAXPY(petsclib::PetscLibType,Y::AbstractPetscMat, a::PetscScalar, X::AbstractPetscMat, str::MatStructure) 
Computes Y = a*X + Y.

Logically Collective

Input Parameters:
- `a`   - the scalar multiplier
- `X`   - the first matrix
- `Y`   - the second matrix
- `str` - either `SAME_NONZERO_PATTERN`, `DIFFERENT_NONZERO_PATTERN`, `UNKNOWN_NONZERO_PATTERN`, or `SUBSET_NONZERO_PATTERN` (nonzeros of `X` is a subset of `Y`'s)

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatAYPX()`

# External Links
$(_doc_external("Mat/MatAXPY"))
"""
function MatAXPY(petsclib::PetscLibType, Y::AbstractPetscMat, a::PetscScalar, X::AbstractPetscMat, str::MatStructure) end

@for_petsc function MatAXPY(petsclib::$UnionPetscLib, Y::AbstractPetscMat, a::$PetscScalar, X::AbstractPetscMat, str::MatStructure )

    @chk ccall(
               (:MatAXPY, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscScalar, CMat, MatStructure),
               Y, a, X, str,
              )


	return nothing
end 

"""
	MatAYPX(petsclib::PetscLibType,Y::AbstractPetscMat, a::PetscScalar, X::AbstractPetscMat, str::MatStructure) 
Computes Y = a*Y + X.

Logically Collective

Input Parameters:
- `a`   - the `PetscScalar` multiplier
- `Y`   - the first matrix
- `X`   - the second matrix
- `str` - either `SAME_NONZERO_PATTERN`, `DIFFERENT_NONZERO_PATTERN`, `UNKNOWN_NONZERO_PATTERN`, or `SUBSET_NONZERO_PATTERN` (nonzeros of `X` is a subset of `Y`'s)

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatAXPY()`

# External Links
$(_doc_external("Mat/MatAYPX"))
"""
function MatAYPX(petsclib::PetscLibType, Y::AbstractPetscMat, a::PetscScalar, X::AbstractPetscMat, str::MatStructure) end

@for_petsc function MatAYPX(petsclib::$UnionPetscLib, Y::AbstractPetscMat, a::$PetscScalar, X::AbstractPetscMat, str::MatStructure )

    @chk ccall(
               (:MatAYPX, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscScalar, CMat, MatStructure),
               Y, a, X, str,
              )


	return nothing
end 

"""
	MatAppendOptionsPrefix(petsclib::PetscLibType,A::AbstractPetscMat, prefix::String) 
Appends to the prefix used for searching for all
matrix options in the database.

Logically Collective

Input Parameters:
- `A`      - the matrix
- `prefix` - the prefix to prepend to all option names

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatGetOptionsPrefix()`, `MatAppendOptionsPrefixFactor()`, `MatSetOptionsPrefix()`

# External Links
$(_doc_external("Mat/MatAppendOptionsPrefix"))
"""
function MatAppendOptionsPrefix(petsclib::PetscLibType, A::AbstractPetscMat, prefix::String) end

@for_petsc function MatAppendOptionsPrefix(petsclib::$UnionPetscLib, A::AbstractPetscMat, prefix::String )

    @chk ccall(
               (:MatAppendOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Cchar}),
               A, prefix,
              )


	return nothing
end 

"""
	MatAppendOptionsPrefixFactor(petsclib::PetscLibType,A::AbstractPetscMat, prefix::String) 
Appends to the prefix used for searching for all matrix factor options in the database for
for matrices created with `MatGetFactor()`

Logically Collective

Input Parameters:
- `A`      - the matrix
- `prefix` - the prefix to prepend to all option names for the factored matrix

Level: developer

-seealso: [](ch_matrices), `Mat`, [Matrix Factorization](sec_matfactor), `MatGetFactor()`, `PetscOptionsCreate()`, `PetscOptionsDestroy()`, `PetscObjectSetOptionsPrefix()`, `PetscObjectPrependOptionsPrefix()`,
`PetscObjectGetOptionsPrefix()`, `TSAppendOptionsPrefix()`, `SNESAppendOptionsPrefix()`, `KSPAppendOptionsPrefix()`, `MatSetOptionsPrefixFactor()`,
`MatSetOptionsPrefix()`

# External Links
$(_doc_external("Mat/MatAppendOptionsPrefixFactor"))
"""
function MatAppendOptionsPrefixFactor(petsclib::PetscLibType, A::AbstractPetscMat, prefix::String) end

@for_petsc function MatAppendOptionsPrefixFactor(petsclib::$UnionPetscLib, A::AbstractPetscMat, prefix::String )

    @chk ccall(
               (:MatAppendOptionsPrefixFactor, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Cchar}),
               A, prefix,
              )


	return nothing
end 

"""
	assembled::PetscBool = MatAssembled(petsclib::PetscLibType,mat::AbstractPetscMat) 
Indicates if a matrix has been assembled and is ready for
use; for example, in matrix-vector product.

Not Collective

Input Parameter:
- `mat` - the matrix

Output Parameter:
- `assembled` - `PETSC_TRUE` or `PETSC_FALSE`

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatAssemblyEnd()`, `MatSetValues()`, `MatAssemblyBegin()`

# External Links
$(_doc_external("Mat/MatAssembled"))
"""
function MatAssembled(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatAssembled(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	assembled_ = Ref{PetscBool}()

    @chk ccall(
               (:MatAssembled, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscBool}),
               mat, assembled_,
              )

	assembled = assembled_[]

	return assembled
end 

"""
	MatAssemblyBegin(petsclib::PetscLibType,mat::AbstractPetscMat, type::MatAssemblyType) 
Begins assembling the matrix.  This routine should
be called after completing all calls to `MatSetValues()`.

Collective

Input Parameters:
- `mat`  - the matrix
- `type` - type of assembly, either `MAT_FLUSH_ASSEMBLY` or `MAT_FINAL_ASSEMBLY`

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatAssemblyEnd()`, `MatSetValues()`, `MatAssembled()`

# External Links
$(_doc_external("Mat/MatAssemblyBegin"))
"""
function MatAssemblyBegin(petsclib::PetscLibType, mat::AbstractPetscMat, type::MatAssemblyType) end

@for_petsc function MatAssemblyBegin(petsclib::$UnionPetscLib, mat::AbstractPetscMat, type::MatAssemblyType )

    @chk ccall(
               (:MatAssemblyBegin, $petsc_library),
               PetscErrorCode,
               (CMat, MatAssemblyType),
               mat, type,
              )


	return nothing
end 

"""
	MatAssemblyEnd(petsclib::PetscLibType,mat::AbstractPetscMat, type::MatAssemblyType) 
Completes assembling the matrix.  This routine should
be called after `MatAssemblyBegin()`.

Collective

Input Parameters:
- `mat`  - the matrix
- `type` - type of assembly, either `MAT_FLUSH_ASSEMBLY` or `MAT_FINAL_ASSEMBLY`

Options Database Keys:
- `-mat_view ::ascii_info`             - Prints info on matrix at conclusion of `MatAssemblyEnd()`
- `-mat_view ::ascii_info_detail`      - Prints more detailed info
- `-mat_view`                          - Prints matrix in ASCII format
- `-mat_view ::ascii_matlab`           - Prints matrix in MATLAB format
- `-mat_view draw`                     - draws nonzero structure of matrix, using `MatView()` and `PetscDrawOpenX()`.
- `-display <name>`                    - Sets display name (default is host)
- `-draw_pause <sec>`                  - Sets number of seconds to pause after display
- `-mat_view socket`                   - Sends matrix to socket, can be accessed from MATLAB (See [Using MATLAB with PETSc](ch_matlab))
- `-viewer_socket_machine <machine>`   - Machine to use for socket
- `-viewer_socket_port <port>`         - Port number to use for socket
- `-mat_view binary:filename[:append]` - Save matrix to file in binary format

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatAssemblyBegin()`, `MatSetValues()`, `PetscDrawOpenX()`, `PetscDrawCreate()`, `MatView()`, `MatAssembled()`, `PetscViewerSocketOpen()`

# External Links
$(_doc_external("Mat/MatAssemblyEnd"))
"""
function MatAssemblyEnd(petsclib::PetscLibType, mat::AbstractPetscMat, type::MatAssemblyType) end

@for_petsc function MatAssemblyEnd(petsclib::$UnionPetscLib, mat::AbstractPetscMat, type::MatAssemblyType )

    @chk ccall(
               (:MatAssemblyEnd, $petsc_library),
               PetscErrorCode,
               (CMat, MatAssemblyType),
               mat, type,
              )


	return nothing
end 

"""
	MatBackwardSolve(petsclib::PetscLibType,mat::AbstractPetscMat, b::AbstractPetscVec, x::AbstractPetscVec) 
Solves U x = b, given a factored matrix, A = LU.
D^(1/2) U x = b, given a factored symmetric matrix, A = U^T*D*U,

Neighbor-wise Collective

Input Parameters:
- `mat` - the factored matrix
- `b`   - the right-hand-side vector

Output Parameter:
- `x` - the result vector

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatForwardSolve()`, `MatGetFactor()`, `MatSolve()`

# External Links
$(_doc_external("Mat/MatBackwardSolve"))
"""
function MatBackwardSolve(petsclib::PetscLibType, mat::AbstractPetscMat, b::AbstractPetscVec, x::AbstractPetscVec) end

@for_petsc function MatBackwardSolve(petsclib::$UnionPetscLib, mat::AbstractPetscMat, b::AbstractPetscVec, x::AbstractPetscVec )

    @chk ccall(
               (:MatBackwardSolve, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               mat, b, x,
              )


	return nothing
end 

"""
	MatBindToCPU(petsclib::PetscLibType,A::AbstractPetscMat, flg::PetscBool) 
marks a matrix to temporarily stay on the CPU and perform computations on the CPU

Logically Collective

Input Parameters:
- `A`   - the matrix
- `flg` - bind to the CPU if value of `PETSC_TRUE`

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatBoundToCPU()`

# External Links
$(_doc_external("Mat/MatBindToCPU"))
"""
function MatBindToCPU(petsclib::PetscLibType, A::AbstractPetscMat, flg::PetscBool) end

@for_petsc function MatBindToCPU(petsclib::$UnionPetscLib, A::AbstractPetscMat, flg::PetscBool )

    @chk ccall(
               (:MatBindToCPU, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               A, flg,
              )


	return nothing
end 

"""
	MatBlockMatSetPreallocation(petsclib::PetscLibType,B::AbstractPetscMat, bs::PetscInt, nz::PetscInt, nnz::Vector{PetscInt}) 
For good matrix assembly performance
the user should preallocate the matrix storage by setting the parameter nz
(or the array nnz).  By setting these parameters accurately, performance
during matrix assembly can be increased by more than a factor of 50.

Collective

Input Parameters:
- `B`   - The matrix
- `bs`  - size of each block in matrix
- `nz`  - number of nonzeros per block row (same for all rows)
- `nnz` - array containing the number of nonzeros in the various block rows
(possibly different for each row) or `NULL`

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `MatCreateBlockMat()`, `MatSetValues()`

# External Links
$(_doc_external("Mat/MatBlockMatSetPreallocation"))
"""
function MatBlockMatSetPreallocation(petsclib::PetscLibType, B::AbstractPetscMat, bs::PetscInt, nz::PetscInt, nnz::Vector{PetscInt}) end

@for_petsc function MatBlockMatSetPreallocation(petsclib::$UnionPetscLib, B::AbstractPetscMat, bs::$PetscInt, nz::$PetscInt, nnz::Vector{$PetscInt} )

    @chk ccall(
               (:MatBlockMatSetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, Ptr{$PetscInt}),
               B, bs, nz, nnz,
              )


	return nothing
end 

"""
	flg::PetscBool = MatBoundToCPU(petsclib::PetscLibType,A::AbstractPetscMat) 
query if a matrix is bound to the CPU

Input Parameter:
- `A` - the matrix

Output Parameter:
- `flg` - the logical flag

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatBindToCPU()`

# External Links
$(_doc_external("Mat/MatBoundToCPU"))
"""
function MatBoundToCPU(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatBoundToCPU(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatBoundToCPU, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscBool}),
               A, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	MatCholeskyFactor(petsclib::PetscLibType,mat::AbstractPetscMat, perm::AbstractIS, info::Vector{MatFactorInfo}) 
Performs in
symmetric matrix.

Collective

Input Parameters:
- `mat`  - the matrix
- `perm` - row and column permutations
- `info` - expected fill as ratio of original fill

Level: developer

-seealso: [](ch_matrices), `Mat`, [Matrix Factorization](sec_matfactor), `MatGetFactor()`, `MatFactorInfo`, `MatLUFactor()`, `MatCholeskyFactorSymbolic()`, `MatCholeskyFactorNumeric()`
`MatGetOrdering()`

# External Links
$(_doc_external("Mat/MatCholeskyFactor"))
"""
function MatCholeskyFactor(petsclib::PetscLibType, mat::AbstractPetscMat, perm::AbstractIS, info::Vector{MatFactorInfo}) end

@for_petsc function MatCholeskyFactor(petsclib::$UnionPetscLib, mat::AbstractPetscMat, perm::AbstractIS, info::Vector{MatFactorInfo} )

    @chk ccall(
               (:MatCholeskyFactor, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, Ptr{MatFactorInfo}),
               mat, perm, info,
              )


	return nothing
end 

"""
	MatCholeskyFactorNumeric(petsclib::PetscLibType,fact::AbstractPetscMat, mat::AbstractPetscMat, info::Vector{MatFactorInfo}) 
Performs numeric Cholesky factorization
of a symmetric matrix. Call this routine after first calling `MatGetFactor()` and
`MatCholeskyFactorSymbolic()`.

Collective

Input Parameters:
- `fact` - the factor matrix obtained with `MatGetFactor()`, where the factored values are stored
- `mat`  - the initial matrix that is to be factored
- `info` - options for factorization

Level: developer

-seealso: [](ch_matrices), `Mat`, [Matrix Factorization](sec_matfactor), `MatFactorInfo`, `MatGetFactor()`, `MatCholeskyFactorSymbolic()`, `MatCholeskyFactor()`, `MatLUFactorNumeric()`

# External Links
$(_doc_external("Mat/MatCholeskyFactorNumeric"))
"""
function MatCholeskyFactorNumeric(petsclib::PetscLibType, fact::AbstractPetscMat, mat::AbstractPetscMat, info::Vector{MatFactorInfo}) end

@for_petsc function MatCholeskyFactorNumeric(petsclib::$UnionPetscLib, fact::AbstractPetscMat, mat::AbstractPetscMat, info::Vector{MatFactorInfo} )

    @chk ccall(
               (:MatCholeskyFactorNumeric, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, Ptr{MatFactorInfo}),
               fact, mat, info,
              )


	return nothing
end 

"""
	MatCholeskyFactorSymbolic(petsclib::PetscLibType,fact::AbstractPetscMat, mat::AbstractPetscMat, perm::AbstractIS, info::Vector{MatFactorInfo}) 
Performs symbolic Cholesky factorization
of a symmetric matrix.

Collective

Input Parameters:
- `fact` - the factor matrix obtained with `MatGetFactor()`
- `mat`  - the matrix
- `perm` - row and column permutations
- `info` - options for factorization, includes
-seealso: [](ch_matrices), `Mat`, [Matrix Factorization](sec_matfactor), `MatFactorInfo`, `MatGetFactor()`, `MatLUFactorSymbolic()`, `MatCholeskyFactor()`, `MatCholeskyFactorNumeric()`
`MatGetOrdering()`

# External Links
$(_doc_external("Mat/MatCholeskyFactorSymbolic"))
"""
function MatCholeskyFactorSymbolic(petsclib::PetscLibType, fact::AbstractPetscMat, mat::AbstractPetscMat, perm::AbstractIS, info::Vector{MatFactorInfo}) end

@for_petsc function MatCholeskyFactorSymbolic(petsclib::$UnionPetscLib, fact::AbstractPetscMat, mat::AbstractPetscMat, perm::AbstractIS, info::Vector{MatFactorInfo} )

    @chk ccall(
               (:MatCholeskyFactorSymbolic, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CIS, Ptr{MatFactorInfo}),
               fact, mat, perm, info,
              )


	return nothing
end 

"""
	MatCompositeAddMat(petsclib::PetscLibType,mat::AbstractPetscMat, smat::AbstractPetscMat) 
Add another matrix to a composite matrix.

Collective

Input Parameters:
- `mat`  - the composite matrix
- `smat` - the partial matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatCreateComposite()`, `MatCompositeGetMat()`, `MATCOMPOSITE`

# External Links
$(_doc_external("Mat/MatCompositeAddMat"))
"""
function MatCompositeAddMat(petsclib::PetscLibType, mat::AbstractPetscMat, smat::AbstractPetscMat) end

@for_petsc function MatCompositeAddMat(petsclib::$UnionPetscLib, mat::AbstractPetscMat, smat::AbstractPetscMat )

    @chk ccall(
               (:MatCompositeAddMat, $petsc_library),
               PetscErrorCode,
               (CMat, CMat),
               mat, smat,
              )


	return nothing
end 

"""
	Ai::PetscMat = MatCompositeGetMat(petsclib::PetscLibType,mat::AbstractPetscMat, i::PetscInt) 
Returns the ith matrix from the composite matrix.

Logically Collective

Input Parameters:
- `mat` - the composite matrix
- `i`   - the number of requested matrix

Output Parameter:
- `Ai` - ith matrix in composite

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatCreateComposite()`, `MatCompositeGetNumberMat()`, `MatCompositeAddMat()`, `MATCOMPOSITE`

# External Links
$(_doc_external("Mat/MatCompositeGetMat"))
"""
function MatCompositeGetMat(petsclib::PetscLibType, mat::AbstractPetscMat, i::PetscInt) end

@for_petsc function MatCompositeGetMat(petsclib::$UnionPetscLib, mat::AbstractPetscMat, i::$PetscInt )
	Ai_ = Ref{CMat}()

    @chk ccall(
               (:MatCompositeGetMat, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CMat}),
               mat, i, Ai_,
              )

	Ai = PetscMat(Ai_[], petsclib)

	return Ai
end 

"""
	str::MatStructure = MatCompositeGetMatStructure(petsclib::PetscLibType,mat::AbstractPetscMat) 
Returns the structure of matrices in the composite matrix.

Not Collective

Input Parameter:
- `mat` - the composite matrix

Output Parameter:
- `str` - structure of the matrices

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatCreateComposite()`, `MatCompositeSetMatStructure()`, `MATCOMPOSITE`

# External Links
$(_doc_external("Mat/MatCompositeGetMatStructure"))
"""
function MatCompositeGetMatStructure(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatCompositeGetMatStructure(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	str_ = Ref{MatStructure}()

    @chk ccall(
               (:MatCompositeGetMatStructure, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatStructure}),
               mat, str_,
              )

	str = str_[]

	return str
end 

"""
	nmat::PetscInt = MatCompositeGetNumberMat(petsclib::PetscLibType,mat::AbstractPetscMat) 
Returns the number of matrices in the composite matrix.

Not Collective

Input Parameter:
- `mat` - the composite matrix

Output Parameter:
- `nmat` - number of matrices in the composite matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatCreateComposite()`, `MatCompositeGetMat()`, `MATCOMPOSITE`

# External Links
$(_doc_external("Mat/MatCompositeGetNumberMat"))
"""
function MatCompositeGetNumberMat(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatCompositeGetNumberMat(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	nmat_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatCompositeGetNumberMat, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}),
               mat, nmat_,
              )

	nmat = nmat_[]

	return nmat
end 

"""
	type::MatCompositeType = MatCompositeGetType(petsclib::PetscLibType,mat::AbstractPetscMat) 
Returns type of composite.

Not Collective

Input Parameter:
- `mat` - the composite matrix

Output Parameter:
- `type` - type of composite

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatCreateComposite()`, `MatCompositeSetType()`, `MATCOMPOSITE`, `MatCompositeType`

# External Links
$(_doc_external("Mat/MatCompositeGetType"))
"""
function MatCompositeGetType(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatCompositeGetType(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	type_ = Ref{MatCompositeType}()

    @chk ccall(
               (:MatCompositeGetType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatCompositeType}),
               mat, type_,
              )

	type = type_[]

	return type
end 

"""
	MatCompositeMerge(petsclib::PetscLibType,mat::AbstractPetscMat) 
Given a composite matrix, replaces it with a "regular" matrix
by summing or computing the product of all the matrices inside the composite matrix.

Collective

Input Parameter:
- `mat` - the composite matrix

Options Database Keys:
- `-mat_composite_merge`      - merge in `MatAssemblyEnd()`
- `-mat_composite_merge_type` - set merge direction

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatDestroy()`, `MatMult()`, `MatCompositeAddMat()`, `MatCreateComposite()`, `MatCompositeSetMatStructure()`, `MatCompositeSetMergeType()`, `MATCOMPOSITE`

# External Links
$(_doc_external("Mat/MatCompositeMerge"))
"""
function MatCompositeMerge(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatCompositeMerge(petsclib::$UnionPetscLib, mat::AbstractPetscMat )

    @chk ccall(
               (:MatCompositeMerge, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	MatCompositeSetMatStructure(petsclib::PetscLibType,mat::AbstractPetscMat, str::MatStructure) 
Indicates structure of matrices in the composite matrix.

Not Collective

Input Parameters:
- `mat` - the composite matrix
- `str` - either `SAME_NONZERO_PATTERN`, `DIFFERENT_NONZERO_PATTERN` (default) or `SUBSET_NONZERO_PATTERN`

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatAXPY()`, `MatCreateComposite()`, `MatCompositeMerge()` `MatCompositeGetMatStructure()`, `MATCOMPOSITE`

# External Links
$(_doc_external("Mat/MatCompositeSetMatStructure"))
"""
function MatCompositeSetMatStructure(petsclib::PetscLibType, mat::AbstractPetscMat, str::MatStructure) end

@for_petsc function MatCompositeSetMatStructure(petsclib::$UnionPetscLib, mat::AbstractPetscMat, str::MatStructure )

    @chk ccall(
               (:MatCompositeSetMatStructure, $petsc_library),
               PetscErrorCode,
               (CMat, MatStructure),
               mat, str,
              )


	return nothing
end 

"""
	MatCompositeSetMergeType(petsclib::PetscLibType,mat::AbstractPetscMat, type::MatCompositeMergeType) 
Sets order of `MatCompositeMerge()`.

Logically Collective

Input Parameters:
- `mat`  - the composite matrix
- `type` - `MAT_COMPOSITE_MERGE RIGHT` (default) to start merge from right with the first added matrix (mat[0]),
`MAT_COMPOSITE_MERGE_LEFT` to start merge from left with the last added matrix (mat[nmat-1])

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatCreateComposite()`, `MatCompositeMerge()`, `MATCOMPOSITE`

# External Links
$(_doc_external("Mat/MatCompositeSetMergeType"))
"""
function MatCompositeSetMergeType(petsclib::PetscLibType, mat::AbstractPetscMat, type::MatCompositeMergeType) end

@for_petsc function MatCompositeSetMergeType(petsclib::$UnionPetscLib, mat::AbstractPetscMat, type::MatCompositeMergeType )

    @chk ccall(
               (:MatCompositeSetMergeType, $petsc_library),
               PetscErrorCode,
               (CMat, MatCompositeMergeType),
               mat, type,
              )


	return nothing
end 

"""
	MatCompositeSetScalings(petsclib::PetscLibType,mat::AbstractPetscMat, scalings::Vector{PetscScalar}) 
Sets separate scaling factors for component matrices.

Logically Collective

Input Parameters:
- `mat`      - the composite matrix
- `scalings` - array of scaling factors with scalings[i] being factor of i-th matrix, for i in [0, nmat)

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatScale()`, `MatDiagonalScale()`, `MATCOMPOSITE`

# External Links
$(_doc_external("Mat/MatCompositeSetScalings"))
"""
function MatCompositeSetScalings(petsclib::PetscLibType, mat::AbstractPetscMat, scalings::Vector{PetscScalar}) end

@for_petsc function MatCompositeSetScalings(petsclib::$UnionPetscLib, mat::AbstractPetscMat, scalings::Vector{$PetscScalar} )

    @chk ccall(
               (:MatCompositeSetScalings, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscScalar}),
               mat, scalings,
              )


	return nothing
end 

"""
	MatCompositeSetType(petsclib::PetscLibType,mat::AbstractPetscMat, type::MatCompositeType) 
Indicates if the matrix is defined as the sum of a set of matrices or the product.

Logically Collective

Input Parameters:
- `mat`  - the composite matrix
- `type` - the `MatCompositeType` to use for the matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatDestroy()`, `MatMult()`, `MatCompositeAddMat()`, `MatCreateComposite()`, `MatCompositeGetType()`, `MATCOMPOSITE`,
`MatCompositeType`

# External Links
$(_doc_external("Mat/MatCompositeSetType"))
"""
function MatCompositeSetType(petsclib::PetscLibType, mat::AbstractPetscMat, type::MatCompositeType) end

@for_petsc function MatCompositeSetType(petsclib::$UnionPetscLib, mat::AbstractPetscMat, type::MatCompositeType )

    @chk ccall(
               (:MatCompositeSetType, $petsc_library),
               PetscErrorCode,
               (CMat, MatCompositeType),
               mat, type,
              )


	return nothing
end 

"""
	bw::PetscInt = MatComputeBandwidth(petsclib::PetscLibType,A::AbstractPetscMat, fraction::PetscReal) 
Calculate the full bandwidth of the matrix, meaning the width 2k+1 where k diagonals on either side are sufficient to contain all the matrix nonzeros.

Collective

Input Parameters:
- `A`        - The `Mat`
- `fraction` - An optional percentage of the Frobenius norm of the matrix that the bandwidth should enclose

Output Parameter:
- `bw` - The matrix bandwidth

Level: beginner

-seealso: `DMPlexCreate()`, `DMPlexSetConeSize()`, `DMPlexSetChart()`

# External Links
$(_doc_external("Mat/MatComputeBandwidth"))
"""
function MatComputeBandwidth(petsclib::PetscLibType, A::AbstractPetscMat, fraction::PetscReal) end

@for_petsc function MatComputeBandwidth(petsclib::$UnionPetscLib, A::AbstractPetscMat, fraction::$PetscReal )
	bw_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatComputeBandwidth, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscReal, Ptr{$PetscInt}),
               A, fraction, bw_,
              )

	bw = bw_[]

	return bw
end 

"""
	mat::PetscMat = MatComputeOperator(petsclib::PetscLibType,inmat::AbstractPetscMat, mattype::MatType) 
Computes the explicit matrix

Collective

Input Parameters:
- `inmat`   - the matrix
- `mattype` - the matrix type for the explicit operator

Output Parameter:
- `mat` - the explicit  operator

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatConvert()`, `MatMult()`, `MatComputeOperatorTranspose()`

# External Links
$(_doc_external("Mat/MatComputeOperator"))
"""
function MatComputeOperator(petsclib::PetscLibType, inmat::AbstractPetscMat, mattype::MatType) end

@for_petsc function MatComputeOperator(petsclib::$UnionPetscLib, inmat::AbstractPetscMat, mattype::MatType )
	mat_ = Ref{CMat}()

    @chk ccall(
               (:MatComputeOperator, $petsc_library),
               PetscErrorCode,
               (CMat, MatType, Ptr{CMat}),
               inmat, mattype, mat_,
              )

	mat = PetscMat(mat_[], petsclib)

	return mat
end 

"""
	mat::PetscMat = MatComputeOperatorTranspose(petsclib::PetscLibType,inmat::AbstractPetscMat, mattype::MatType) 
Computes the explicit matrix representation of
a give matrix that can apply `MatMultTranspose()`

Collective

Input Parameters:
- `inmat`   - the matrix
- `mattype` - the matrix type for the explicit operator

Output Parameter:
- `mat` - the explicit  operator transposed

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatConvert()`, `MatMult()`, `MatComputeOperator()`

# External Links
$(_doc_external("Mat/MatComputeOperatorTranspose"))
"""
function MatComputeOperatorTranspose(petsclib::PetscLibType, inmat::AbstractPetscMat, mattype::MatType) end

@for_petsc function MatComputeOperatorTranspose(petsclib::$UnionPetscLib, inmat::AbstractPetscMat, mattype::MatType )
	mat_ = Ref{CMat}()

    @chk ccall(
               (:MatComputeOperatorTranspose, $petsc_library),
               PetscErrorCode,
               (CMat, MatType, Ptr{CMat}),
               inmat, mattype, mat_,
              )

	mat = PetscMat(mat_[], petsclib)

	return mat
end 

"""
	MatComputeVariableBlockEnvelope(petsclib::PetscLibType,mat::AbstractPetscMat) 
Given a matrix whose nonzeros are in blocks along the diagonal this computes and stores
the sizes of these blocks in the matrix. An individual block may lie over several processes.

Collective

Input Parameter:
- `mat` - the matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatInvertVariableBlockEnvelope()`, `MatSetVariableBlockSizes()`

# External Links
$(_doc_external("Mat/MatComputeVariableBlockEnvelope"))
"""
function MatComputeVariableBlockEnvelope(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatComputeVariableBlockEnvelope(petsclib::$UnionPetscLib, mat::AbstractPetscMat )

    @chk ccall(
               (:MatComputeVariableBlockEnvelope, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	MatConjugate(petsclib::PetscLibType,mat::AbstractPetscMat) 
replaces the matrix values with their complex conjugates

Logically Collective

Input Parameter:
- `mat` - the matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatRealPart()`, `MatImaginaryPart()`, `VecConjugate()`, `MatTranspose()`

# External Links
$(_doc_external("Mat/MatConjugate"))
"""
function MatConjugate(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatConjugate(petsclib::$UnionPetscLib, mat::AbstractPetscMat )

    @chk ccall(
               (:MatConjugate, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	value::PetscScalar = MatConstantDiagonalGetConstant(petsclib::PetscLibType,mat::AbstractPetscMat) 
Get the scalar constant of a constant diagonal matrix

Not collective

Input Parameter:
- `mat` - a `MATCONSTANTDIAGONAL`

Output Parameter:
- `value` - the scalar value

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatDestroy()`, `MATCONSTANTDIAGONAL`

# External Links
$(_doc_external("Mat/MatConstantDiagonalGetConstant"))
"""
function MatConstantDiagonalGetConstant(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatConstantDiagonalGetConstant(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	value_ = Ref{$PetscScalar}()

    @chk ccall(
               (:MatConstantDiagonalGetConstant, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscScalar}),
               mat, value_,
              )

	value = value_[]

	return value
end 

"""
	M::PetscMat = MatConvert(petsclib::PetscLibType,mat::AbstractPetscMat, newtype::MatType, reuse::MatReuse) 
Converts a matrix to another matrix, either of the same
or different type.

Collective

Input Parameters:
- `mat`     - the matrix
- `newtype` - new matrix type.  Use `MATSAME` to create a new matrix of the
same type as the original matrix.
- `reuse`   - denotes if the destination matrix is to be created or reused.
Use `MAT_INPLACE_MATRIX` for inplace conversion (that is when you want the input `Mat` to be changed to contain the matrix in the new format), otherwise use
`MAT_INITIAL_MATRIX` or `MAT_REUSE_MATRIX` (can only be used after the first call was made with `MAT_INITIAL_MATRIX`, causes the matrix space in M to be reused).

Output Parameter:
- `M` - pointer to place new matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatCopy()`, `MatDuplicate()`, `MAT_INITIAL_MATRIX`, `MAT_REUSE_MATRIX`, `MAT_INPLACE_MATRIX`

# External Links
$(_doc_external("Mat/MatConvert"))
"""
function MatConvert(petsclib::PetscLibType, mat::AbstractPetscMat, newtype::MatType, reuse::MatReuse) end

@for_petsc function MatConvert(petsclib::$UnionPetscLib, mat::AbstractPetscMat, newtype::MatType, reuse::MatReuse )
	M_ = Ref{CMat}()

    @chk ccall(
               (:MatConvert, $petsc_library),
               PetscErrorCode,
               (CMat, MatType, MatReuse, Ptr{CMat}),
               mat, newtype, reuse, M_,
              )

	M = PetscMat(M_[], petsclib)

	return M
end 

"""
	MatCopy(petsclib::PetscLibType,A::AbstractPetscMat, B::AbstractPetscMat, str::MatStructure) 
Copies a matrix to another matrix.

Collective

Input Parameters:
- `A`   - the matrix
- `str` - `SAME_NONZERO_PATTERN` or `DIFFERENT_NONZERO_PATTERN`

Output Parameter:
- `B` - where the copy is put

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatConvert()`, `MatDuplicate()`

# External Links
$(_doc_external("Mat/MatCopy"))
"""
function MatCopy(petsclib::PetscLibType, A::AbstractPetscMat, B::AbstractPetscMat, str::MatStructure) end

@for_petsc function MatCopy(petsclib::$UnionPetscLib, A::AbstractPetscMat, B::AbstractPetscMat, str::MatStructure )

    @chk ccall(
               (:MatCopy, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, MatStructure),
               A, B, str,
              )


	return nothing
end 

"""
	MatCopyHashToXAIJ(petsclib::PetscLibType,A::AbstractPetscMat, B::AbstractPetscMat) 
copy hash table entries into an XAIJ matrix type

Logically Collective

Input Parameter:
- `A` - A matrix in unassembled, hash table form

Output Parameter:
- `B` - The XAIJ matrix. This can either be `A` or some matrix of equivalent size, e.g. obtained from `A` via `MatDuplicate()`

Example:
-seealso: [](ch_matrices), `Mat`, `MAT_USE_HASH_TABLE`

# External Links
$(_doc_external("Mat/MatCopyHashToXAIJ"))
"""
function MatCopyHashToXAIJ(petsclib::PetscLibType, A::AbstractPetscMat, B::AbstractPetscMat) end

@for_petsc function MatCopyHashToXAIJ(petsclib::$UnionPetscLib, A::AbstractPetscMat, B::AbstractPetscMat )

    @chk ccall(
               (:MatCopyHashToXAIJ, $petsc_library),
               PetscErrorCode,
               (CMat, CMat),
               A, B,
              )


	return nothing
end 

"""
	A::PetscMat = MatCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Creates a matrix where the type is determined
from either a call to `MatSetType()` or from the options database
with a call to `MatSetFromOptions()`.

Collective

Input Parameter:
- `comm` - MPI communicator

Output Parameter:
- `A` - the matrix

Options Database Keys:
- `-mat_type seqaij`   - `MATSEQAIJ` type, uses `MatCreateSeqAIJ()`
- `-mat_type mpiaij`   - `MATMPIAIJ` type, uses `MatCreateAIJ()`
- `-mat_type seqdense` - `MATSEQDENSE`, uses `MatCreateSeqDense()`
- `-mat_type mpidense` - `MATMPIDENSE` type, uses `MatCreateDense()`
- `-mat_type seqbaij`  - `MATSEQBAIJ` type, uses `MatCreateSeqBAIJ()`
- `-mat_type mpibaij`  - `MATMPIBAIJ` type, uses `MatCreateBAIJ()`

See the manpages for particular formats (e.g., `MATSEQAIJ`)
for additional format-specific options.

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatCreateSeqAIJ()`, `MatCreateAIJ()`,
`MatCreateSeqDense()`, `MatCreateDense()`,
`MatCreateSeqBAIJ()`, `MatCreateBAIJ()`,
`MatCreateSeqSBAIJ()`, `MatCreateSBAIJ()`,
`MatConvert()`

# External Links
$(_doc_external("Mat/MatCreate"))
"""
function MatCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function MatCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{CMat}),
               comm, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	A::PetscMat = MatCreateAIJ(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_nz::PetscInt, d_nnz::Union{Ptr, Vector{PetscInt}}, o_nz::PetscInt, o_nnz::Union{Ptr, Vector{PetscInt}}) 
Creates a sparse parallel matrix in `MATAIJ` format
(the default parallel PETSc format).  For good matrix assembly performance
the user should preallocate the matrix storage by setting the parameters
`d_nz` (or `d_nnz`) and `o_nz` (or `o_nnz`).

Collective

Input Parameters:
- `comm`  - MPI communicator
- `m`     - number of local rows (or `PETSC_DECIDE` to have calculated if M is given)
This value should be the same as the local size used in creating the
y vector for the matrix-vector product y = Ax.
- `n`     - This value should be the same as the local size used in creating the
x vector for the matrix-vector product y = Ax. (or `PETSC_DECIDE` to have
calculated if N is given) For square matrices n is almost always m.
- `M`     - number of global rows (or `PETSC_DETERMINE` to have calculated if m is given)
- `N`     - number of global columns (or `PETSC_DETERMINE` to have calculated if n is given)
- `d_nz`  - number of nonzeros per row in DIAGONAL portion of local submatrix
(same value is used for all local rows)
- `d_nnz` - array containing the number of nonzeros in the various rows of the
DIAGONAL portion of the local submatrix (possibly different for each row)
or `NULL`, if `d_nz` is used to specify the nonzero structure.
The size of this array is equal to the number of local rows, i.e 'm'.
- `o_nz`  - number of nonzeros per row in the OFF-DIAGONAL portion of local
submatrix (same value is used for all local rows).
- `o_nnz` - array containing the number of nonzeros in the various rows of the
OFF-DIAGONAL portion of the local submatrix (possibly different for
each row) or `NULL`, if `o_nz` is used to specify the nonzero
structure. The size of this array is equal to the number
of local rows, i.e 'm'.

Output Parameter:
- `A` - the matrix

Options Database Keys:
- `-mat_no_inode`                     - Do not use inodes
- `-mat_inode_limit <limit>`          - Sets inode limit (max limit=5)
- `-matmult_vecscatter_view <viewer>` - View the vecscatter (i.e., communication pattern) used in `MatMult()` of sparse parallel matrices.
See viewer types in manual of `MatView()`. Of them, ascii_matlab, draw or binary cause the `VecScatter`
to be viewed as a matrix. Entry (i,j) is the size of message (in bytes) rank i sends to rank j in one `MatMult()` call.

Level: intermediate

-seealso: [](ch_matrices), `Mat`, [Sparse Matrix Creation](sec_matsparse), `MatCreate()`, `MatCreateSeqAIJ()`, `MatSetValues()`, `MatMPIAIJSetPreallocation()`, `MatMPIAIJSetPreallocationCSR()`,
`MATMPIAIJ`, `MatCreateMPIAIJWithArrays()`, `MatGetOwnershipRange()`, `MatGetOwnershipRanges()`, `MatGetOwnershipRangeColumn()`,
`MatGetOwnershipRangesColumn()`, `PetscLayout`

# External Links
$(_doc_external("Mat/MatCreateAIJ"))
"""
function MatCreateAIJ(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_nz::PetscInt, d_nnz::Union{Ptr, Vector{PetscInt}}, o_nz::PetscInt, o_nnz::Union{Ptr, Vector{PetscInt}}) end

@for_petsc function MatCreateAIJ(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, d_nz::$PetscInt, d_nnz::Union{Ptr, Vector{$PetscInt}}, o_nz::$PetscInt, o_nnz::Union{Ptr, Vector{$PetscInt}} )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateAIJ, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, Ptr{CMat}),
               comm, m, n, M, N, d_nz, d_nnz, o_nz, o_nnz, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	A::PetscMat = MatCreateAIJKokkos(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_nz::PetscInt, d_nnz::Vector{PetscInt}, o_nz::PetscInt, o_nnz::Vector{PetscInt}) 

# External Links
$(_doc_external("Mat/MatCreateAIJKokkos"))
"""
function MatCreateAIJKokkos(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_nz::PetscInt, d_nnz::Vector{PetscInt}, o_nz::PetscInt, o_nnz::Vector{PetscInt}) end

@for_petsc function MatCreateAIJKokkos(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, d_nz::$PetscInt, d_nnz::Vector{$PetscInt}, o_nz::$PetscInt, o_nnz::Vector{$PetscInt} )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateAIJKokkos, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, Ptr{CMat}),
               comm, m, n, M, N, d_nz, d_nnz, o_nz, o_nnz, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	A::PetscMat = MatCreateAIJViennaCL(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_nz::PetscInt, d_nnz::Vector{PetscInt}, o_nz::PetscInt, o_nnz::Vector{PetscInt}) 

# External Links
$(_doc_external("Mat/MatCreateAIJViennaCL"))
"""
function MatCreateAIJViennaCL(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_nz::PetscInt, d_nnz::Vector{PetscInt}, o_nz::PetscInt, o_nnz::Vector{PetscInt}) end

@for_petsc function MatCreateAIJViennaCL(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, d_nz::$PetscInt, d_nnz::Vector{$PetscInt}, o_nz::$PetscInt, o_nnz::Vector{$PetscInt} )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateAIJViennaCL, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, Ptr{CMat}),
               comm, m, n, M, N, d_nz, d_nnz, o_nz, o_nnz, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	A::PetscMat = MatCreateBAIJ(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_nz::PetscInt, d_nnz::Union{Ptr, Vector{PetscInt}}, o_nz::PetscInt, o_nnz::Union{Ptr, Vector{PetscInt}}) 
Creates a sparse parallel matrix in `MATBAIJ` format
(block compressed row).

Collective

Input Parameters:
- `comm`  - MPI communicator
- `bs`    - size of block, the blocks are ALWAYS square. One can use `MatSetBlockSizes()` to set a different row and column blocksize but the row
blocksize always defines the size of the blocks. The column blocksize sets the blocksize of the vectors obtained with `MatCreateVecs()`
- `m`     - number of local rows (or `PETSC_DECIDE` to have calculated if M is given)
This value should be the same as the local size used in creating the
y vector for the matrix-vector product y = Ax.
- `n`     - number of local columns (or `PETSC_DECIDE` to have calculated if N is given)
This value should be the same as the local size used in creating the
x vector for the matrix-vector product y = Ax.
- `M`     - number of global rows (or `PETSC_DETERMINE` to have calculated if m is given)
- `N`     - number of global columns (or `PETSC_DETERMINE` to have calculated if n is given)
- `d_nz`  - number of nonzero blocks per block row in diagonal portion of local
submatrix  (same for all local rows)
- `d_nnz` - array containing the number of nonzero blocks in the various block rows
of the in diagonal portion of the local (possibly different for each block
row) or NULL.  If you plan to factor the matrix you must leave room for the diagonal entry
and set it even if it is zero.
- `o_nz`  - number of nonzero blocks per block row in the off-diagonal portion of local
submatrix (same for all local rows).
- `o_nnz` - array containing the number of nonzero blocks in the various block rows of the
off-diagonal portion of the local submatrix (possibly different for
each block row) or NULL.

Output Parameter:
- `A` - the matrix

Options Database Keys:
- `-mat_block_size`            - size of the blocks to use
- `-mat_use_hash_table <fact>` - set hash table factor

Level: intermediate

-seealso: `Mat`, `MatCreate()`, `MatCreateSeqBAIJ()`, `MatSetValues()`, `MatMPIBAIJSetPreallocation()`, `MatMPIBAIJSetPreallocationCSR()`,
`MatGetOwnershipRange()`,  `MatGetOwnershipRanges()`, `MatGetOwnershipRangeColumn()`, `MatGetOwnershipRangesColumn()`, `PetscLayout`

# External Links
$(_doc_external("Mat/MatCreateBAIJ"))
"""
function MatCreateBAIJ(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_nz::PetscInt, d_nnz::Union{Ptr, Vector{PetscInt}}, o_nz::PetscInt, o_nnz::Union{Ptr, Vector{PetscInt}}) end

@for_petsc function MatCreateBAIJ(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, d_nz::$PetscInt, d_nnz::Union{Ptr, Vector{$PetscInt}}, o_nz::$PetscInt, o_nnz::Union{Ptr, Vector{$PetscInt}} )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateBAIJ, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, Ptr{CMat}),
               comm, bs, m, n, M, N, d_nz, d_nnz, o_nz, o_nnz, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	A::PetscMat = MatCreateBAIJMKL(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_nz::PetscInt, d_nnz::Union{Ptr, Vector{PetscInt}}, o_nz::PetscInt, o_nnz::Union{Ptr, Vector{PetscInt}}) 
Creates a sparse parallel matrix in `MATBAIJMKL` format (block compressed row).

Collective

Input Parameters:
- `comm`  - MPI communicator
- `bs`    - size of block, the blocks are ALWAYS square. One can use `MatSetBlockSizes()` to set a different row and column blocksize but the row
blocksize always defines the size of the blocks. The column blocksize sets the blocksize of the vectors obtained with `MatCreateVecs()`
- `m`     - number of local rows (or `PETSC_DECIDE` to have calculated if `M` is given)
This value should be the same as the local size used in creating the
y vector for the matrix-vector product y = Ax.
- `n`     - number of local columns (or `PETSC_DECIDE` to have calculated if `N` is given)
This value should be the same as the local size used in creating the
x vector for the matrix-vector product y = Ax.
- `M`     - number of global rows (or `PETSC_DETERMINE` to have calculated if `m` is given)
- `N`     - number of global columns (or `PETSC_DETERMINE` to have calculated if `n` is given)
- `d_nz`  - number of nonzero blocks per block row in diagonal portion of local
submatrix  (same for all local rows)
- `d_nnz` - array containing the number of nonzero blocks in the various block rows
of the in diagonal portion of the local (possibly different for each block
row) or `NULL`.  If you plan to factor the matrix you must leave room for the diagonal entry
and set it even if it is zero.
- `o_nz`  - number of nonzero blocks per block row in the off-diagonal portion of local
submatrix (same for all local rows).
- `o_nnz` - array containing the number of nonzero blocks in the various block rows of the
off-diagonal portion of the local submatrix (possibly different for
each block row) or `NULL`.

Output Parameter:
- `A` - the matrix

Options Database Keys:
- `-mat_block_size`            - size of the blocks to use
- `-mat_use_hash_table <fact>` - set hash table factor

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATBAIJMKL`, `MATBAIJ`, `MatCreate()`, `MatCreateSeqBAIJMKL()`, `MatSetValues()`, `MatMPIBAIJSetPreallocation()`, `MatMPIBAIJSetPreallocationCSR()`

# External Links
$(_doc_external("Mat/MatCreateBAIJMKL"))
"""
function MatCreateBAIJMKL(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_nz::PetscInt, d_nnz::Union{Ptr, Vector{PetscInt}}, o_nz::PetscInt, o_nnz::Union{Ptr, Vector{PetscInt}}) end

@for_petsc function MatCreateBAIJMKL(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, d_nz::$PetscInt, d_nnz::Union{Ptr, Vector{$PetscInt}}, o_nz::$PetscInt, o_nnz::Union{Ptr, Vector{$PetscInt}} )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateBAIJMKL, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, Ptr{CMat}),
               comm, bs, m, n, M, N, d_nz, d_nnz, o_nz, o_nnz, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	nnz::PetscInt,A::PetscMat = MatCreateBlockMat(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, bs::PetscInt, nz::PetscInt) 
Creates a new matrix in which each block contains a uniform

Collective

Input Parameters:
- `comm` - MPI communicator
- `m`    - number of rows
- `n`    - number of columns
- `bs`   - size of each submatrix
- `nz`   - expected maximum number of nonzero blocks in row (use `PETSC_DEFAULT` if not known)
- `nnz`  - expected number of nonzers per block row if known (use `NULL` otherwise)

Output Parameter:
- `A` - the matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATBLOCKMAT`, `MatCreateNest()`

# External Links
$(_doc_external("Mat/MatCreateBlockMat"))
"""
function MatCreateBlockMat(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, bs::PetscInt, nz::PetscInt) end

@for_petsc function MatCreateBlockMat(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, bs::$PetscInt, nz::$PetscInt )
	nnz_ = Ref{$PetscInt}()
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateBlockMat, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{CMat}),
               comm, m, n, bs, nz, nnz_, A_,
              )

	nnz = nnz_[]
	A = PetscMat(A_[], petsclib)

	return nnz,A
end 

"""
	C::PetscMat = MatCreateCentering(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt, N::PetscInt) 
Creates a new matrix object that implements the (symmetric and idempotent) centering matrix,  I

Collective

Input Parameters:
- `comm` - MPI communicator
- `n`    - number of local rows (or `PETSC_DECIDE` to have calculated if `N` is given)
This value should be the same as the local size used in creating the
`y` vector for the matrix-vector product y = Ax.
- `N`    - number of global rows (or `PETSC_DETERMINE` to have calculated if `n` is given)

Output Parameter:
- `C` - the matrix

-seealso: [](ch_matrices), `Mat`, `MatCreateLRC()`, `MatCreateComposite()`

# External Links
$(_doc_external("Mat/MatCreateCentering"))
"""
function MatCreateCentering(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, N::PetscInt) end

@for_petsc function MatCreateCentering(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, N::$PetscInt )
	C_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateCentering, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{CMat}),
               comm, n, N, C_,
              )

	C = PetscMat(C_[], petsclib)

	return C
end 

"""
	mat::PetscMat = MatCreateComposite(petsclib::PetscLibType,comm::MPI_Comm, nmat::PetscInt, mats::AbstractPetscMat) 
Creates a matrix as the sum or product of one or more matrices

Collective

Input Parameters:
- `comm` - MPI communicator
- `nmat` - number of matrices to put in
- `mats` - the matrices

Output Parameter:
- `mat` - the matrix

Options Database Keys:
- `-mat_composite_merge`       - merge in `MatAssemblyEnd()`
- `-mat_composite_merge_mvctx` - merge Mvctx of component matrices to optimize communication in `MatMult()` for ADDITIVE matrices
- `-mat_composite_merge_type`  - set merge direction

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatDestroy()`, `MatMult()`, `MatCompositeAddMat()`, `MatCompositeGetMat()`, `MatCompositeMerge()`, `MatCompositeSetType()`,
`MATCOMPOSITE`, `MatCompositeType`

# External Links
$(_doc_external("Mat/MatCreateComposite"))
"""
function MatCreateComposite(petsclib::PetscLibType, comm::MPI_Comm, nmat::PetscInt, mats::AbstractPetscMat) end

@for_petsc function MatCreateComposite(petsclib::$UnionPetscLib, comm::MPI_Comm, nmat::$PetscInt, mats::AbstractPetscMat )
	mats_ = Ref(mats.ptr)
	mat_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateComposite, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, Ptr{CMat}, Ptr{CMat}),
               comm, nmat, mats_, mat_,
              )

	mats.ptr = mats_[]
	mat = PetscMat(mat_[], petsclib)

	return mat
end 

"""
	J::PetscMat = MatCreateConstantDiagonal(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, diag::PetscScalar) 
Creates a matrix with a uniform value along the diagonal

Collective

Input Parameters:
- `comm` - MPI communicator
- `m`    - number of local rows (or `PETSC_DECIDE` to have calculated if `M` is given)
This value should be the same as the local size used in creating the
y vector for the matrix-vector product y = Ax.
- `n`    - This value should be the same as the local size used in creating the
x vector for the matrix-vector product y = Ax. (or `PETSC_DECIDE` to have
calculated if `N` is given) For square matrices n is almost always `m`.
- `M`    - number of global rows (or `PETSC_DETERMINE` to have calculated if m is given)
- `N`    - number of global columns (or `PETSC_DETERMINE` to have calculated if n is given)
- `diag` - the diagonal value

Output Parameter:
- `J` - the diagonal matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatDestroy()`, `MATCONSTANTDIAGONAL`, `MatScale()`, `MatShift()`, `MatMult()`, `MatGetDiagonal()`, `MatGetFactor()`, `MatSolve()`

# External Links
$(_doc_external("Mat/MatCreateConstantDiagonal"))
"""
function MatCreateConstantDiagonal(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, diag::PetscScalar) end

@for_petsc function MatCreateConstantDiagonal(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, diag::$PetscScalar )
	J_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateConstantDiagonal, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscScalar, Ptr{CMat}),
               comm, m, n, M, N, diag, J_,
              )

	J = PetscMat(J_[], petsclib)

	return J
end 

"""
	A::PetscMat = MatCreateDense(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, data::Union{Ptr, Vector{PetscScalar}}) 
Creates a matrix in `MATDENSE` format.

Collective

Input Parameters:
- `comm` - MPI communicator
- `m`    - number of local rows (or `PETSC_DECIDE` to have calculated if `M` is given)
- `n`    - number of local columns (or `PETSC_DECIDE` to have calculated if `N` is given)
- `M`    - number of global rows (or `PETSC_DECIDE` to have calculated if `m` is given)
- `N`    - number of global columns (or `PETSC_DECIDE` to have calculated if `n` is given)
- `data` - optional location of matrix data.  Set data to `NULL` (`PETSC_NULL_SCALAR_ARRAY` for Fortran users) for PETSc
to control all matrix memory allocation.

Output Parameter:
- `A` - the matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MatCreate()`, `MatCreateSeqDense()`, `MatSetValues()`

# External Links
$(_doc_external("Mat/MatCreateDense"))
"""
function MatCreateDense(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, data::Union{Ptr, Vector{PetscScalar}}) end

@for_petsc function MatCreateDense(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, data::Union{Ptr, Vector{$PetscScalar}} )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateDense, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{CMat}),
               comm, m, n, M, N, data, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	data::PetscScalar,A::PetscMat = MatCreateDenseFromVecType(petsclib::PetscLibType,comm::MPI_Comm, vtype::VecType, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, lda::PetscInt) 
Create a matrix that matches the type of a Vec.

Collective

Input Parameters:
- `comm`  - the communicator
- `vtype` - the vector type
- `m`     - number of local rows (or `PETSC_DECIDE` to have calculated if `M` is given)
- `n`     - number of local columns (or `PETSC_DECIDE` to have calculated if `N` is given)
- `M`     - number of global rows (or `PETSC_DECIDE` to have calculated if `m` is given)
- `N`     - number of global columns (or `PETSC_DECIDE` to have calculated if `n` is given)
- `lda`   - optional leading dimension. Pass any non-positive number to use the default.
- `data`  - optional location of matrix data, which should have the same memory type as the vector. Pass `NULL` to have PETSc take care of matrix memory allocation.

Output Parameter:
- `A` - the dense matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatCreateDense()`, `MatCreateDenseCUDA()`, `MatCreateDenseHIP()`, `PetscMemType`

# External Links
$(_doc_external("Mat/MatCreateDenseFromVecType"))
"""
function MatCreateDenseFromVecType(petsclib::PetscLibType, comm::MPI_Comm, vtype::VecType, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, lda::PetscInt) end

@for_petsc function MatCreateDenseFromVecType(petsclib::$UnionPetscLib, comm::MPI_Comm, vtype::VecType, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, lda::$PetscInt )
	data_ = Ref{$PetscScalar}()
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateDenseFromVecType, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, VecType, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{CMat}),
               comm, vtype, m, n, M, N, lda, data_, A_,
              )

	data = data_[]
	A = PetscMat(A_[], petsclib)

	return data,A
end 

"""
	J::PetscMat = MatCreateDiagonal(petsclib::PetscLibType,diag::AbstractPetscVec) 
Creates a matrix defined by a given vector along its diagonal.

Collective

Input Parameter:
- `diag` - vector for the diagonal

Output Parameter:
- `J` - the diagonal matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatDestroy()`, `MATCONSTANTDIAGONAL`, `MatScale()`, `MatShift()`, `MatMult()`, `MatGetDiagonal()`, `MatSolve()`
`MatDiagonalRestoreInverseDiagonal()`, `MatDiagonalGetDiagonal()`, `MatDiagonalRestoreDiagonal()`, `MatDiagonalGetInverseDiagonal()`

# External Links
$(_doc_external("Mat/MatCreateDiagonal"))
"""
function MatCreateDiagonal(petsclib::PetscLibType, diag::AbstractPetscVec) end

@for_petsc function MatCreateDiagonal(petsclib::$UnionPetscLib, diag::AbstractPetscVec )
	J_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateDiagonal, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{CMat}),
               diag, J_,
              )

	J = PetscMat(J_[], petsclib)

	return J
end 

"""
	A::PetscMat = MatCreateFFT(petsclib::PetscLibType,comm::MPI_Comm, ndim::PetscInt, dim::Vector{PetscInt}, mattype::MatType) 
Creates a matrix object that provides FFT via an external package

Collective

Input Parameters:
- `comm`    - MPI communicator
- `ndim`    - the ndim-dimensional transform
- `dim`     - array of size ndim, dim[i] contains the vector length in the i-dimension
- `mattype` - package type, e.g., `MATFFTW` or `MATSEQCUFFT`

Output Parameter:
- `A` - the matrix

Options Database Key:
- `-mat_fft_type` - set FFT type fft or seqcufft

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATFFTW`, `MATSEQCUFFT`, `MatCreateVecsFFTW()`

# External Links
$(_doc_external("Mat/MatCreateFFT"))
"""
function MatCreateFFT(petsclib::PetscLibType, comm::MPI_Comm, ndim::PetscInt, dim::Vector{PetscInt}, mattype::MatType) end

@for_petsc function MatCreateFFT(petsclib::$UnionPetscLib, comm::MPI_Comm, ndim::$PetscInt, dim::Vector{$PetscInt}, mattype::MatType )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateFFT, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, Ptr{$PetscInt}, MatType, Ptr{CMat}),
               comm, ndim, dim, mattype, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	A::PetscMat = MatCreateFromOptions(petsclib::PetscLibType,comm::MPI_Comm, prefix::String, bs::PetscInt, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt) 
Creates a matrix whose type is set from the options database

Collective

Input Parameters:
- `comm`   - MPI communicator
- `prefix` - [optional] prefix for the options database
- `bs`     - the blocksize (commonly 1)
- `m`      - the local number of rows (or `PETSC_DECIDE`)
- `n`      - the local number of columns (or `PETSC_DECIDE` or `PETSC_DETERMINE`)
- `M`      - the global number of rows (or `PETSC_DETERMINE`)
- `N`      - the global number of columns (or `PETSC_DETERMINE`)

Output Parameter:
- `A` - the matrix

Options Database Key:
- `-mat_type` - see `MatType`, for example `aij`, `aijcusparse`, `baij`, `sbaij`, `dense`, defaults to `aij`

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatCreateSeqAIJ()`, `MatCreateAIJ()`,
`MatCreateSeqDense()`, `MatCreateDense()`,
`MatCreateSeqBAIJ()`, `MatCreateBAIJ()`,
`MatCreateSeqSBAIJ()`, `MatCreateSBAIJ()`,
`MatConvert()`, `MatCreate()`

# External Links
$(_doc_external("Mat/MatCreateFromOptions"))
"""
function MatCreateFromOptions(petsclib::PetscLibType, comm::MPI_Comm, prefix::String, bs::PetscInt, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt) end

@for_petsc function MatCreateFromOptions(petsclib::$UnionPetscLib, comm::MPI_Comm, prefix::String, bs::$PetscInt, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateFromOptions, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{CMat}),
               comm, prefix, bs, m, n, M, N, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	graph::PetscMat = MatCreateGraph(petsclib::PetscLibType,A::AbstractPetscMat, sym::PetscBool, scale::PetscBool, filter::PetscReal, num_idx::PetscInt, index::Vector{PetscInt}) 
create a scalar matrix (that is a matrix with one vertex for each block vertex in the original matrix), for use in graph algorithms
and possibly removes small values from the graph structure.

Collective

Input Parameters:
- `A`       - the matrix
- `sym`     - `PETSC_TRUE` indicates that the graph should be symmetrized
- `scale`   - `PETSC_TRUE` indicates that the graph edge weights should be symmetrically scaled with the diagonal entry
- `filter`  - filter value - < 0: does nothing; == 0: removes only 0.0 entries; otherwise: removes entries with abs(entries) <= value
- `num_idx` - size of 'index' array
- `index`   - array of block indices to use for graph strength of connection weight

Output Parameter:
- `graph` - the resulting graph

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `PCGAMG`

# External Links
$(_doc_external("Mat/MatCreateGraph"))
"""
function MatCreateGraph(petsclib::PetscLibType, A::AbstractPetscMat, sym::PetscBool, scale::PetscBool, filter::PetscReal, num_idx::PetscInt, index::Vector{PetscInt}) end

@for_petsc function MatCreateGraph(petsclib::$UnionPetscLib, A::AbstractPetscMat, sym::PetscBool, scale::PetscBool, filter::$PetscReal, num_idx::$PetscInt, index::Vector{$PetscInt} )
	graph_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateGraph, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool, PetscBool, $PetscReal, $PetscInt, Ptr{$PetscInt}, Ptr{CMat}),
               A, sym, scale, filter, num_idx, index, graph_,
              )

	graph = PetscMat(graph_[], petsclib)

	return graph
end 

"""
	N::PetscMat = MatCreateHermitianTranspose(petsclib::PetscLibType,A::AbstractPetscMat) 
Creates a new matrix object of `MatType` `MATHERMITIANTRANSPOSEVIRTUAL` that behaves like A'*

Collective

Input Parameter:
- `A` - the (possibly rectangular) matrix

Output Parameter:
- `N` - the matrix that represents A'*

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatCreateNormal()`, `MatMult()`, `MatMultHermitianTranspose()`, `MatCreate()`,
`MATTRANSPOSEVIRTUAL`, `MatCreateTranspose()`, `MatHermitianTransposeGetMat()`, `MATNORMAL`, `MATNORMALHERMITIAN`

# External Links
$(_doc_external("Mat/MatCreateHermitianTranspose"))
"""
function MatCreateHermitianTranspose(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatCreateHermitianTranspose(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	N_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateHermitianTranspose, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, N_,
              )

	N = PetscMat(N_[], petsclib)

	return N
end 

"""
	B::PetscMat = MatCreateHtoolFromKernel(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, spacedim::PetscInt, coords_target::Vector{PetscReal}, coords_source::Vector{PetscReal}, kernel::Ptr{Cvoid}, kernelctx::Ptr{Cvoid}) 

# External Links
$(_doc_external("Mat/MatCreateHtoolFromKernel"))
"""
function MatCreateHtoolFromKernel(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, spacedim::PetscInt, coords_target::Vector{PetscReal}, coords_source::Vector{PetscReal}, kernel::Ptr{Cvoid}, kernelctx::Ptr{Cvoid}) end

@for_petsc function MatCreateHtoolFromKernel(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, spacedim::$PetscInt, coords_target::Vector{$PetscReal}, coords_source::Vector{$PetscReal}, kernel::Ptr{Cvoid}, kernelctx::Ptr{Cvoid} )
	B_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateHtoolFromKernel, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{CMat}),
               comm, m, n, M, N, spacedim, coords_target, coords_source, kernel, kernelctx, B_,
              )

	B = PetscMat(B_[], petsclib)

	return B
end 

"""
	A::PetscMat = MatCreateIS(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, rmap::ISLocalToGlobalMapping, cmap::ISLocalToGlobalMapping) 
Creates a "process" unassembled matrix.

Collective.

Input Parameters:
- `comm` - MPI communicator that will share the matrix
- `bs`   - block size of the matrix
- `m`    - local size of left vector used in matrix vector products
- `n`    - local size of right vector used in matrix vector products
- `M`    - global size of left vector used in matrix vector products
- `N`    - global size of right vector used in matrix vector products
- `rmap` - local to global map for rows
- `cmap` - local to global map for cols

Output Parameter:
- `A` - the resulting matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATIS`, `MatSetLocalToGlobalMapping()`

# External Links
$(_doc_external("Mat/MatCreateIS"))
"""
function MatCreateIS(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, rmap::ISLocalToGlobalMapping, cmap::ISLocalToGlobalMapping) end

@for_petsc function MatCreateIS(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, rmap::ISLocalToGlobalMapping, cmap::ISLocalToGlobalMapping )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateIS, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, ISLocalToGlobalMapping, ISLocalToGlobalMapping, Ptr{CMat}),
               comm, bs, m, n, M, N, rmap, cmap, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	kaij::PetscMat = MatCreateKAIJ(petsclib::PetscLibType,A::AbstractPetscMat, p::PetscInt, q::PetscInt, S::Union{Ptr, Vector{PetscScalar}}, T::Union{Ptr, Vector{PetscScalar}}) 
Creates a matrix of type `MATKAIJ`.

Collective

Input Parameters:
- `A` - the `MATAIJ` matrix
- `p` - number of rows in `S` and `T`
- `q` - number of columns in `S` and `T`
- `S` - the `S` matrix (can be `NULL`), stored as a `PetscScalar` array (column-major)
- `T` - the `T` matrix (can be `NULL`), stored as a `PetscScalar` array (column-major)

Output Parameter:
- `kaij` - the new `MATKAIJ` matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatKAIJSetAIJ()`, `MatKAIJSetS()`, `MatKAIJSetT()`, `MatKAIJGetAIJ()`, `MatKAIJGetS()`, `MatKAIJGetT()`, `MATKAIJ`

# External Links
$(_doc_external("Mat/MatCreateKAIJ"))
"""
function MatCreateKAIJ(petsclib::PetscLibType, A::AbstractPetscMat, p::PetscInt, q::PetscInt, S::Union{Ptr, Vector{PetscScalar}}, T::Union{Ptr, Vector{PetscScalar}}) end

@for_petsc function MatCreateKAIJ(petsclib::$UnionPetscLib, A::AbstractPetscMat, p::$PetscInt, q::$PetscInt, S::Union{Ptr, Vector{$PetscScalar}}, T::Union{Ptr, Vector{$PetscScalar}} )
	kaij_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateKAIJ, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{$PetscScalar}, Ptr{CMat}),
               A, p, q, S, T, kaij_,
              )

	kaij = PetscMat(kaij_[], petsclib)

	return kaij
end 

"""
	B::PetscMat = MatCreateLMVMBFGS(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt, N::PetscInt) 
Creates a limited
matrix used for approximating Jacobians. L-BFGS is symmetric positive-definite by
construction, and is commonly used to approximate Hessians in optimization
problems.

To use the L-BFGS matrix with other vector types, the matrix must be
created using `MatCreate()` and `MatSetType()`, followed by `MatLMVMAllocate()`.
This ensures that the internal storage and work vectors are duplicated from the
correct type of vector.

Collective

Input Parameters:
- `comm` - MPI communicator
- `n`    - number of local rows for storage vectors
- `N`    - global size of the storage vectors

Output Parameter:
- `B` - the matrix

Options Database Keys:
- `-mat_lmvm_scale_type` - (developer) type of scaling applied to J0 (none, scalar, diagonal)
- `-mat_lmvm_theta`      - (developer) convex ratio between BFGS and DFP components of the diagonal J0 scaling
- `-mat_lmvm_rho`        - (developer) update limiter for the J0 scaling
- `-mat_lmvm_alpha`      - (developer) coefficient factor for the quadratic subproblem in J0 scaling
- `-mat_lmvm_beta`       - (developer) exponential factor for the diagonal J0 scaling
- `-mat_lmvm_sigma_hist` - (developer) number of past updates to use in J0 scaling

Level: intermediate

-seealso: [](ch_ksp), `MatCreate()`, `MATLMVM`, `MATLMVMBFGS`, `MatCreateLMVMDFP()`, `MatCreateLMVMSR1()`,
`MatCreateLMVMBroyden()`, `MatCreateLMVMBadBroyden()`, `MatCreateLMVMSymBroyden()`

# External Links
$(_doc_external("KSP/MatCreateLMVMBFGS"))
"""
function MatCreateLMVMBFGS(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, N::PetscInt) end

@for_petsc function MatCreateLMVMBFGS(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, N::$PetscInt )
	B_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateLMVMBFGS, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{CMat}),
               comm, n, N, B_,
              )

	B = PetscMat(B_[], petsclib)

	return B
end 

"""
	B::PetscMat = MatCreateLMVMBadBroyden(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt, N::PetscInt) 
Creates a limited
approximation matrix used for a Jacobian. L-BadBrdn is not guaranteed to be
symmetric or positive-definite.

To use the L-BadBrdn matrix with other vector types, the matrix must be
created using `MatCreate()` and `MatSetType()`, followed by `MatLMVMAllocate()`.
This ensures that the internal storage and work vectors are duplicated from the
correct type of vector.

Collective

Input Parameters:
- `comm` - MPI communicator
- `n`    - number of local rows for storage vectors
- `N`    - global size of the storage vectors

Output Parameter:
- `B` - the matrix

Options Database Keys:
- `-mat_lmvm_hist_size`         - the number of history vectors to keep
- `-mat_lmvm_mult_algorithm`    - the algorithm to use for multiplication (recursive, dense, compact_dense)
- `-mat_lmvm_cache_J0_products` - whether products between the base Jacobian J0 and history vectors should be cached or recomputed
- `-mat_lmvm_debug`             - (developer) perform internal debugging checks

Level: intermediate

-seealso: [](ch_ksp), `MatCreate()`, `MATLMVM`, `MATLMVMBADBRDN`, `MatCreateLMVMDFP()`, `MatCreateLMVMSR1()`,
`MatCreateLMVMBFGS()`, `MatCreateLMVMBroyden()`, `MatCreateLMVMSymBroyden()`

# External Links
$(_doc_external("KSP/MatCreateLMVMBadBroyden"))
"""
function MatCreateLMVMBadBroyden(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, N::PetscInt) end

@for_petsc function MatCreateLMVMBadBroyden(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, N::$PetscInt )
	B_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateLMVMBadBroyden, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{CMat}),
               comm, n, N, B_,
              )

	B = PetscMat(B_[], petsclib)

	return B
end 

"""
	B::PetscMat = MatCreateLMVMBroyden(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt, N::PetscInt) 
Creates a limited
matrix used for a Jacobian. L-Brdn is not guaranteed to be symmetric or
positive-definite.

To use the L-Brdn matrix with other vector types, the matrix must be
created using `MatCreate()` and `MatSetType()`, followed by `MatLMVMAllocate()`.
This ensures that the internal storage and work vectors are duplicated from the
correct type of vector.

Collective

Input Parameters:
- `comm` - MPI communicator
- `n`    - number of local rows for storage vectors
- `N`    - global size of the storage vectors

Output Parameter:
- `B` - the matrix

Options Database Keys:
- `-mat_lmvm_hist_size`         - the number of history vectors to keep
- `-mat_lmvm_mult_algorithm`    - the algorithm to use for multiplication (recursive, dense, compact_dense)
- `-mat_lmvm_cache_J0_products` - whether products between the base Jacobian J0 and history vectors should be cached or recomputed
- `-mat_lmvm_debug`             - (developer) perform internal debugging checks

Level: intermediate

-seealso: [](ch_ksp), `MatCreate()`, `MATLMVM`, `MATLMVMBRDN`, `MatCreateLMVMDFP()`, `MatCreateLMVMSR1()`,
`MatCreateLMVMBFGS()`, `MatCreateLMVMBadBroyden()`, `MatCreateLMVMSymBroyden()`

# External Links
$(_doc_external("KSP/MatCreateLMVMBroyden"))
"""
function MatCreateLMVMBroyden(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, N::PetscInt) end

@for_petsc function MatCreateLMVMBroyden(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, N::$PetscInt )
	B_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateLMVMBroyden, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{CMat}),
               comm, n, N, B_,
              )

	B = PetscMat(B_[], petsclib)

	return B
end 

"""
	B::PetscMat = MatCreateLMVMDBFGS(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt, N::PetscInt) 
Creates a dense representation of the limited
Broyden-Fletcher-Goldfarb-Shanno (BFGS) approximation to a Hessian.

Collective

Input Parameters:
- `comm` - MPI communicator
- `n`    - number of local rows for storage vectors
- `N`    - global size of the storage vectors

Output Parameter:
- `B` - the matrix

Level: advanced

-seealso: `MatCreate()`, `MATLMVM`, `MATLMVMDBFGS`, `MatCreateLMVMBFGS()`

# External Links
$(_doc_external("KSP/MatCreateLMVMDBFGS"))
"""
function MatCreateLMVMDBFGS(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, N::PetscInt) end

@for_petsc function MatCreateLMVMDBFGS(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, N::$PetscInt )
	B_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateLMVMDBFGS, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{CMat}),
               comm, n, N, B_,
              )

	B = PetscMat(B_[], petsclib)

	return B
end 

"""
	B::PetscMat = MatCreateLMVMDDFP(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt, N::PetscInt) 
Creates a dense representation of the limited
Davidon-Fletcher-Powell (DFP) approximation to a Hessian.

Collective

Input Parameters:
- `comm` - MPI communicator
- `n`    - number of local rows for storage vectors
- `N`    - global size of the storage vectors

Output Parameter:
- `B` - the matrix

Level: advanced

-seealso: `MatCreate()`, `MATLMVM`, `MATLMVMDDFP`, `MatCreateLMVMDFP()`

# External Links
$(_doc_external("KSP/MatCreateLMVMDDFP"))
"""
function MatCreateLMVMDDFP(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, N::PetscInt) end

@for_petsc function MatCreateLMVMDDFP(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, N::$PetscInt )
	B_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateLMVMDDFP, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{CMat}),
               comm, n, N, B_,
              )

	B = PetscMat(B_[], petsclib)

	return B
end 

"""
	B::PetscMat = MatCreateLMVMDFP(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt, N::PetscInt) 
Creates a limited
used for approximating Jacobians. L-DFP is symmetric positive-definite by
construction, and is the dual of L-BFGS where Y and S vectors swap roles.

To use the L-DFP matrix with other vector types, the matrix must be
created using `MatCreate()` and `MatSetType()`, followed by `MatLMVMAllocate()`.
This ensures that the internal storage and work vectors are duplicated from the
correct type of vector.

Collective

Input Parameters:
- `comm` - MPI communicator
- `n`    - number of local rows for storage vectors
- `N`    - global size of the storage vectors

Output Parameter:
- `B` - the matrix

Options Database Keys:
- `-mat_lmvm_scale_type` - (developer) type of scaling applied to J0 (none, scalar, diagonal)
- `-mat_lmvm_theta`      - (developer) convex ratio between BFGS and DFP components of the diagonal J0 scaling
- `-mat_lmvm_rho`        - (developer) update limiter for the J0 scaling
- `-mat_lmvm_alpha`      - (developer) coefficient factor for the quadratic subproblem in J0 scaling
- `-mat_lmvm_beta`       - (developer) exponential factor for the diagonal J0 scaling
- `-mat_lmvm_sigma_hist` - (developer) number of past updates to use in J0 scaling

Level: intermediate

-seealso: [](ch_ksp), `MatCreate()`, `MATLMVM`, `MATLMVMDFP`, `MatCreateLMVMBFGS()`, `MatCreateLMVMSR1()`,
`MatCreateLMVMBroyden()`, `MatCreateLMVMBadBroyden()`, `MatCreateLMVMSymBroyden()`

# External Links
$(_doc_external("KSP/MatCreateLMVMDFP"))
"""
function MatCreateLMVMDFP(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, N::PetscInt) end

@for_petsc function MatCreateLMVMDFP(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, N::$PetscInt )
	B_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateLMVMDFP, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{CMat}),
               comm, n, N, B_,
              )

	B = PetscMat(B_[], petsclib)

	return B
end 

"""
	B::PetscMat = MatCreateLMVMDQN(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt, N::PetscInt) 
Creates a dense representation of the limited
Quasi-Newton approximation to a Hessian.

Collective

Input Parameters:
- `comm` - MPI communicator
- `n`    - number of local rows for storage vectors
- `N`    - global size of the storage vectors

Output Parameter:
- `B` - the matrix

Level: advanced

-seealso: `MatCreate()`, `MATLMVM`, `MATLMVMDBFGS`, `MATLMVMDDFP`, `MatCreateLMVMDDFP()`, `MatCreateLMVMDBFGS()`

# External Links
$(_doc_external("KSP/MatCreateLMVMDQN"))
"""
function MatCreateLMVMDQN(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, N::PetscInt) end

@for_petsc function MatCreateLMVMDQN(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, N::$PetscInt )
	B_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateLMVMDQN, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{CMat}),
               comm, n, N, B_,
              )

	B = PetscMat(B_[], petsclib)

	return B
end 

"""
	B::PetscMat = MatCreateLMVMDiagBroyden(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt, N::PetscInt) 
DiagBrdn creates a symmetric Broyden
for approximating Hessians.

Collective

Input Parameters:
- `comm` - MPI communicator
- `n`    - number of local rows for storage vectors
- `N`    - global size of the storage vectors

Output Parameter:
- `B` - the matrix

Options Database Keys:
- `-mat_lmvm_theta`      - (developer) convex ratio between BFGS and DFP components of the diagonal J0 scaling
- `-mat_lmvm_rho`        - (developer) update limiter for the J0 scaling
- `-mat_lmvm_alpha`      - (developer) coefficient factor for the quadratic subproblem in J0 scaling
- `-mat_lmvm_beta`       - (developer) exponential factor for the diagonal J0 scaling
- `-mat_lmvm_sigma_hist` - (developer) number of past updates to use in J0 scaling.
- `-mat_lmvm_tol`        - (developer) tolerance for bounding the denominator of the rescaling away from 0.
- `-mat_lmvm_forward`    - (developer) whether or not to use the forward or backward Broyden update to the diagonal

Level: intermediate

-seealso: [](ch_ksp), `MatCreate()`, `MATLMVM`, `MATLMVMDIAGBRDN`, `MatCreateLMVMDFP()`, `MatCreateLMVMSR1()`,
`MatCreateLMVMBFGS()`, `MatCreateLMVMBroyden()`, `MatCreateLMVMSymBroyden()`

# External Links
$(_doc_external("KSP/MatCreateLMVMDiagBroyden"))
"""
function MatCreateLMVMDiagBroyden(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, N::PetscInt) end

@for_petsc function MatCreateLMVMDiagBroyden(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, N::$PetscInt )
	B_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateLMVMDiagBroyden, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{CMat}),
               comm, n, N, B_,
              )

	B = PetscMat(B_[], petsclib)

	return B
end 

"""
	B::PetscMat = MatCreateLMVMSR1(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt, N::PetscInt) 
Creates a limited
matrix used for a Jacobian. L-SR1 is symmetric by construction, but is not
guaranteed to be positive-definite.

To use the L-SR1 matrix with other vector types, the matrix must be
created using `MatCreate()` and `MatSetType()`, followed by `MatLMVMAllocate()`.
This ensures that the internal storage and work vectors are duplicated from the
correct type of vector.

Collective

Input Parameters:
- `comm` - MPI communicator
- `n`    - number of local rows for storage vectors
- `N`    - global size of the storage vectors

Output Parameter:
- `B` - the matrix

Options Database Keys:
- `-mat_lmvm_hist_size`         - the number of history vectors to keep
- `-mat_lmvm_mult_algorithm`    - the algorithm to use for multiplication (recursive, dense, compact_dense)
- `-mat_lmvm_cache_J0_products` - whether products between the base Jacobian J0 and history vectors should be cached or recomputed
- `-mat_lmvm_eps`               - (developer) numerical zero tolerance for testing when an update should be skipped
- `-mat_lmvm_debug`             - (developer) perform internal debugging checks

Level: intermediate

-seealso: [](ch_ksp), `MatCreate()`, `MATLMVM`, `MATLMVMSR1`, `MatCreateLMVMBFGS()`, `MatCreateLMVMDFP()`,
`MatCreateLMVMBroyden()`, `MatCreateLMVMBadBroyden()`, `MatCreateLMVMSymBroyden()`

# External Links
$(_doc_external("KSP/MatCreateLMVMSR1"))
"""
function MatCreateLMVMSR1(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, N::PetscInt) end

@for_petsc function MatCreateLMVMSR1(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, N::$PetscInt )
	B_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateLMVMSR1, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{CMat}),
               comm, n, N, B_,
              )

	B = PetscMat(B_[], petsclib)

	return B
end 

"""
	B::PetscMat = MatCreateLMVMSymBadBroyden(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt, N::PetscInt) 
Creates a limited
for approximating Jacobians.

Collective

Input Parameters:
- `comm` - MPI communicator
- `n`    - number of local rows for storage vectors
- `N`    - global size of the storage vectors

Output Parameter:
- `B` - the matrix

Options Database Keys:
- `-mat_lmvm_hist_size`         - the number of history vectors to keep
- `-mat_lmvm_psi`               - convex ratio between BFGS and DFP components of the update
- `-mat_lmvm_scale_type`        - type of scaling applied to J0 (none, scalar, diagonal)
- `-mat_lmvm_mult_algorithm`    - the algorithm to use for multiplication (recursive, dense, compact_dense)
- `-mat_lmvm_cache_J0_products` - whether products between the base Jacobian J0 and history vectors should be cached or recomputed
- `-mat_lmvm_eps`               - (developer) numerical zero tolerance for testing when an update should be skipped
- `-mat_lmvm_debug`             - (developer) perform internal debugging checks
- `-mat_lmvm_theta`             - (developer) convex ratio between BFGS and DFP components of the diagonal J0 scaling
- `-mat_lmvm_rho`               - (developer) update limiter for the J0 scaling
- `-mat_lmvm_alpha`             - (developer) coefficient factor for the quadratic subproblem in J0 scaling
- `-mat_lmvm_beta`              - (developer) exponential factor for the diagonal J0 scaling
- `-mat_lmvm_sigma_hist`        - (developer) number of past updates to use in J0 scaling

Level: intermediate

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MatCreate()`, `MATLMVM`, `MATLMVMSYMBROYDEN`, `MatCreateLMVMDFP()`, `MatCreateLMVMSR1()`,
`MatCreateLMVMBFGS()`, `MatCreateLMVMBroyden()`, `MatCreateLMVMBadBroyden()`

# External Links
$(_doc_external("KSP/MatCreateLMVMSymBadBroyden"))
"""
function MatCreateLMVMSymBadBroyden(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, N::PetscInt) end

@for_petsc function MatCreateLMVMSymBadBroyden(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, N::$PetscInt )
	B_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateLMVMSymBadBroyden, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{CMat}),
               comm, n, N, B_,
              )

	B = PetscMat(B_[], petsclib)

	return B
end 

"""
	B::PetscMat = MatCreateLMVMSymBroyden(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt, N::PetscInt) 
Creates a limited
for approximating Jacobians.

Collective

Input Parameters:
- `comm` - MPI communicator, set to `PETSC_COMM_SELF`
- `n`    - number of local rows for storage vectors
- `N`    - global size of the storage vectors

Output Parameter:
- `B` - the matrix

Options Database Keys:
- `-mat_lmvm_hist_size`         - the number of history vectors to keep
- `-mat_lmvm_phi`               - convex ratio between BFGS and DFP components of the update
- `-mat_lmvm_scale_type`        - type of scaling applied to J0 (none, scalar, diagonal)
- `-mat_lmvm_mult_algorithm`    - the algorithm to use for multiplication (recursive, dense, compact_dense)
- `-mat_lmvm_cache_J0_products` - whether products between the base Jacobian J0 and history vectors should be cached or recomputed
- `-mat_lmvm_eps`               - (developer) numerical zero tolerance for testing when an update should be skipped
- `-mat_lmvm_debug`             - (developer) perform internal debugging checks
- `-mat_lmvm_theta`             - (developer) convex ratio between BFGS and DFP components of the diagonal J0 scaling
- `-mat_lmvm_rho`               - (developer) update limiter for the J0 scaling
- `-mat_lmvm_alpha`             - (developer) coefficient factor for the quadratic subproblem in J0 scaling
- `-mat_lmvm_beta`              - (developer) exponential factor for the diagonal J0 scaling
- `-mat_lmvm_sigma_hist`        - (developer) number of past updates to use in J0 scaling

Level: intermediate

-seealso: [](ch_ksp), `MatCreate()`, `MATLMVM`, `MATLMVMSYMBROYDEN`, `MatCreateLMVMDFP()`, `MatCreateLMVMSR1()`,
`MatCreateLMVMBFGS()`, `MatCreateLMVMBroyden()`, `MatCreateLMVMBadBroyden()`

# External Links
$(_doc_external("KSP/MatCreateLMVMSymBroyden"))
"""
function MatCreateLMVMSymBroyden(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, N::PetscInt) end

@for_petsc function MatCreateLMVMSymBroyden(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, N::$PetscInt )
	B_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateLMVMSymBroyden, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{CMat}),
               comm, n, N, B_,
              )

	B = PetscMat(B_[], petsclib)

	return B
end 

"""
	N::PetscMat = MatCreateLRC(petsclib::PetscLibType,A::AbstractPetscMat, U::AbstractPetscMat, c::Union{Ptr, AbstractPetscVec}, V::AbstractPetscMat) 
Creates a new matrix object that behaves like A + U*C*V' of type `MATLRC`

Collective

Input Parameters:
- `A` - the (sparse) matrix (can be `NULL`)
- `U` - dense rectangular (tall and skinny) matrix
- `V` - dense rectangular (tall and skinny) matrix
- `c` - a vector containing the diagonal of C (can be `NULL`)

Output Parameter:
- `N` - the matrix that represents A + U*C*V'

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATLRC`, `MatLRCGetMats()`

# External Links
$(_doc_external("Mat/MatCreateLRC"))
"""
function MatCreateLRC(petsclib::PetscLibType, A::AbstractPetscMat, U::AbstractPetscMat, c::Union{Ptr, AbstractPetscVec}, V::AbstractPetscMat) end

@for_petsc function MatCreateLRC(petsclib::$UnionPetscLib, A::AbstractPetscMat, U::AbstractPetscMat, c::Union{Ptr, AbstractPetscVec}, V::AbstractPetscMat )
	N_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateLRC, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CVec, CMat, Ptr{CMat}),
               A, U, c, V, N_,
              )

	N = PetscMat(N_[], petsclib)

	return N
end 

"""
	L::PetscMat = MatCreateLaplacian(petsclib::PetscLibType,A::AbstractPetscMat, tol::PetscReal, weighted::PetscBool) 
Create the matrix Laplacian, with all values in the matrix less than the tolerance set to zero

Input Parameters:
- `A`        - The matrix
- `tol`      - The zero tolerance
- `weighted` - Flag for using edge weights

Output Parameter:
- `L` - The graph Laplacian matrix

Level: intermediate

-seealso: `MatFilter()`, `MatGetGraph()`

# External Links
$(_doc_external("MatGraphOperations/MatCreateLaplacian"))
"""
function MatCreateLaplacian(petsclib::PetscLibType, A::AbstractPetscMat, tol::PetscReal, weighted::PetscBool) end

@for_petsc function MatCreateLaplacian(petsclib::$UnionPetscLib, A::AbstractPetscMat, tol::$PetscReal, weighted::PetscBool )
	L_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateLaplacian, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscReal, PetscBool, Ptr{CMat}),
               A, tol, weighted, L_,
              )

	L = PetscMat(L_[], petsclib)

	return L
end 

"""
	newmat::PetscMat = MatCreateLocalRef(petsclib::PetscLibType,A::AbstractPetscMat, isrow::AbstractIS, iscol::AbstractIS) 
Gets a logical reference to a local submatrix, for use in assembly, that is to set values into the matrix

Not Collective

Input Parameters:
- `A`     - full matrix, generally parallel
- `isrow` - Local index set for the rows
- `iscol` - Local index set for the columns

Output Parameter:
- `newmat` - new serial `Mat`

Level: developer

-seealso: [](ch_matrices), `Mat`, `MATSUBMATRIX`, `MatCreateSubMatrixVirtual()`, `MatSetValuesLocal()`, `MatSetValuesBlockedLocal()`, `MatGetLocalSubMatrix()`, `MatCreateSubMatrix()`

# External Links
$(_doc_external("Mat/MatCreateLocalRef"))
"""
function MatCreateLocalRef(petsclib::PetscLibType, A::AbstractPetscMat, isrow::AbstractIS, iscol::AbstractIS) end

@for_petsc function MatCreateLocalRef(petsclib::$UnionPetscLib, A::AbstractPetscMat, isrow::AbstractIS, iscol::AbstractIS )
	newmat_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateLocalRef, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, CIS, Ptr{CMat}),
               A, isrow, iscol, newmat_,
              )

	newmat = PetscMat(newmat_[], petsclib)

	return newmat
end 

"""
	maij::PetscMat = MatCreateMAIJ(petsclib::PetscLibType,A::AbstractPetscMat, dof::PetscInt) 
Creates a matrix type providing restriction and interpolation
operations for multicomponent problems.  It interpolates each component the same
way independently.  The matrix type is based on `MATSEQAIJ` for sequential matrices,
and `MATMPIAIJ` for distributed matrices.

Collective

Input Parameters:
- `A`   - the `MATAIJ` matrix describing the action on blocks
- `dof` - the block size (number of components per node)

Output Parameter:
- `maij` - the new `MATMAIJ` matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATAIJ`, `MATMAIJ`, `MatMAIJGetAIJ()`, `MatMAIJRedimension()`

# External Links
$(_doc_external("Mat/MatCreateMAIJ"))
"""
function MatCreateMAIJ(petsclib::PetscLibType, A::AbstractPetscMat, dof::PetscInt) end

@for_petsc function MatCreateMAIJ(petsclib::$UnionPetscLib, A::AbstractPetscMat, dof::$PetscInt )
	maij_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateMAIJ, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CMat}),
               A, dof, maij_,
              )

	maij = PetscMat(maij_[], petsclib)

	return maij
end 

"""
	J::PetscMat = MatCreateMFFD(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt) 
Creates a matrix
approximately multiply a vector by the matrix (Jacobian) . See also `MatCreateSNESMF()`

Collective

Input Parameters:
- `comm` - MPI communicator
- `m`    - number of local rows (or `PETSC_DECIDE` to have calculated if `M` is given)
This value should be the same as the local size used in creating the
y vector for the matrix-vector product y = Ax.
- `n`    - This value should be the same as the local size used in creating the
x vector for the matrix-vector product y = Ax. (or `PETSC_DECIDE` to have
calculated if `N` is given) For square matrices `n` is almost always `m`.
- `M`    - number of global rows (or `PETSC_DETERMINE` to have calculated if `m` is given)
- `N`    - number of global columns (or `PETSC_DETERMINE` to have calculated if `n` is given)

Output Parameter:
- `J` - the matrix-free matrix

Options Database Keys:
- `-mat_mffd_type`             - wp or ds (see `MATMFFD_WP` or `MATMFFD_DS`)
- `-mat_mffd_err`              - square root of estimated relative error in function evaluation
- `-mat_mffd_period`           - how often h is recomputed, defaults to 1, every time
- `-mat_mffd_check_positivity` - possibly decrease `h` until U + h*a has only positive values
- `-mat_mffd_umin <umin>`      - Sets umin (for default PETSc routine that computes h only)
- `-mat_mffd_complex`          - use the Lyness trick with complex numbers to compute the matrix-vector product instead of differencing
(requires real valued functions but that PETSc be configured for complex numbers)
- `-snes_mf`                   - use the finite difference based matrix-free matrix with `SNESSolve()` and no preconditioner
- `-snes_mf_operator`          - use the finite difference based matrix-free matrix with `SNESSolve()` but construct a preconditioner
using the matrix passed as `pmat` to `SNESSetJacobian()`.

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATMFFD`, `MatDestroy()`, `MatMFFDSetFunctionError()`, `MatMFFDDSSetUmin()`, `MatMFFDSetFunction()`
`MatMFFDSetHHistory()`, `MatMFFDResetHHistory()`, `MatCreateSNESMF()`, `MatCreateShell()`, `MATSHELL`,
`MatMFFDGetH()`, `MatMFFDRegister()`, `MatMFFDComputeJacobian()`

# External Links
$(_doc_external("Mat/MatCreateMFFD"))
"""
function MatCreateMFFD(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt) end

@for_petsc function MatCreateMFFD(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt )
	J_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateMFFD, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{CMat}),
               comm, m, n, M, N, J_,
              )

	J = PetscMat(J_[], petsclib)

	return J
end 

"""
	A::PetscMat = MatCreateMPIAIJCRL(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Vector{PetscInt}, onz::PetscInt, onnz::Vector{PetscInt}) 
Creates a sparse matrix of type `MATMPIAIJCRL`.

Collective

Input Parameters:
- `comm` - MPI communicator, set to `PETSC_COMM_SELF`
- `m`    - number of rows
- `n`    - number of columns
- `nz`   - number of nonzeros per row (same for all rows), for the "diagonal" submatrix
- `nnz`  - array containing the number of nonzeros in the various rows (possibly different for each row) or `NULL`, for the "diagonal" submatrix
- `onz`  - number of nonzeros per row (same for all rows), for the "off-diagonal" submatrix
- `onnz` - array containing the number of nonzeros in the various rows (possibly different for each row) or `NULL`, for the "off-diagonal" submatrix

Output Parameter:
- `A` - the matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, [Sparse Matrix Creation](sec_matsparse), `MATAIJ`, `MATAIJSELL`, `MATAIJPERM`, `MATAIJMKL`, `MatCreate()`, `MatCreateMPIAIJPERM()`, `MatSetValues()`

# External Links
$(_doc_external("Mat/MatCreateMPIAIJCRL"))
"""
function MatCreateMPIAIJCRL(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Vector{PetscInt}, onz::PetscInt, onnz::Vector{PetscInt}) end

@for_petsc function MatCreateMPIAIJCRL(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, nz::$PetscInt, nnz::Vector{$PetscInt}, onz::$PetscInt, onnz::Vector{$PetscInt} )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateMPIAIJCRL, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, Ptr{CMat}),
               comm, m, n, nz, nnz, onz, onnz, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	A::PetscMat = MatCreateMPIAIJMKL(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_nz::PetscInt, d_nnz::Vector{PetscInt}, o_nz::PetscInt, o_nnz::Vector{PetscInt}) 
Creates a sparse parallel matrix whose local
portions are stored as `MATSEQAIJMKL` matrices (a matrix class that inherits
from `MATSEQAIJ` but uses some operations provided by Intel MKL).

Collective

Input Parameters:
- `comm`  - MPI communicator
- `m`     - number of local rows (or `PETSC_DECIDE` to have calculated if `M` is given)
This value should be the same as the local size used in creating the
y vector for the matrix-vector product y = Ax.
- `n`     - This value should be the same as the local size used in creating the
x vector for the matrix-vector product y = Ax. (or `PETSC_DECIDE` to have
calculated if N is given) For square matrices n is almost always `m`.
- `M`     - number of global rows (or `PETSC_DETERMINE` to have calculated if `m` is given)
- `N`     - number of global columns (or `PETSC_DETERMINE` to have calculated if `n` is given)
- `d_nz`  - number of nonzeros per row in DIAGONAL portion of local submatrix
(same value is used for all local rows)
- `d_nnz` - array containing the number of nonzeros in the various rows of the
DIAGONAL portion of the local submatrix (possibly different for each row)
or `NULL`, if `d_nz` is used to specify the nonzero structure.
The size of this array is equal to the number of local rows, i.e `m`.
For matrices you plan to factor you must leave room for the diagonal entry and
put in the entry even if it is zero.
- `o_nz`  - number of nonzeros per row in the OFF-DIAGONAL portion of local
submatrix (same value is used for all local rows).
- `o_nnz` - array containing the number of nonzeros in the various rows of the
OFF-DIAGONAL portion of the local submatrix (possibly different for
each row) or `NULL`, if `o_nz` is used to specify the nonzero
structure. The size of this array is equal to the number
of local rows, i.e `m`.

Output Parameter:
- `A` - the matrix

Options Database Key:
- `-mat_aijmkl_no_spmv2` - disables use of the SpMV2 inspector-executor routines

Level: intermediate

-seealso: [](ch_matrices), `Mat`, [Sparse Matrix Creation](sec_matsparse), `MATMPIAIJMKL`, `MatCreate()`, `MatCreateSeqAIJMKL()`,
`MatSetValues()`, `MatGetOwnershipRange()`, `MatGetOwnershipRanges()`, `MatGetOwnershipRangeColumn()`,
`MatGetOwnershipRangesColumn()`, `PetscLayout`

# External Links
$(_doc_external("Mat/MatCreateMPIAIJMKL"))
"""
function MatCreateMPIAIJMKL(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_nz::PetscInt, d_nnz::Vector{PetscInt}, o_nz::PetscInt, o_nnz::Vector{PetscInt}) end

@for_petsc function MatCreateMPIAIJMKL(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, d_nz::$PetscInt, d_nnz::Vector{$PetscInt}, o_nz::$PetscInt, o_nnz::Vector{$PetscInt} )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateMPIAIJMKL, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, Ptr{CMat}),
               comm, m, n, M, N, d_nz, d_nnz, o_nz, o_nnz, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	A::PetscMat = MatCreateMPIAIJPERM(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_nz::PetscInt, d_nnz::Union{Ptr, Vector{PetscInt}}, o_nz::PetscInt, o_nnz::Union{Ptr, Vector{PetscInt}}) 
Creates a sparse parallel matrix whose local
portions are stored as `MATSEQAIJPERM` matrices (a matrix class that inherits
from SEQAIJ but includes some optimizations to allow more effective
vectorization).

Collective

Input Parameters:
- `comm`  - MPI communicator
- `m`     - number of local rows (or `PETSC_DECIDE` to have calculated if `M` is given)
This value should be the same as the local size used in creating the
y vector for the matrix-vector product y = Ax.
- `n`     - This value should be the same as the local size used in creating the
x vector for the matrix-vector product y = Ax. (or PETSC_DECIDE to have
calculated if `N` is given) For square matrices `n` is almost always `m`.
- `M`     - number of global rows (or `PETSC_DETERMINE` to have calculated if `m` is given)
- `N`     - number of global columns (or `PETSC_DETERMINE` to have calculated if `n` is given)
- `d_nz`  - number of nonzeros per row in DIAGONAL portion of local submatrix
(same value is used for all local rows)
- `d_nnz` - array containing the number of nonzeros in the various rows of the
DIAGONAL portion of the local submatrix (possibly different for each row)
or `NULL`, if `d_nz` is used to specify the nonzero structure.
The size of this array is equal to the number of local rows, i.e `m`.
For matrices you plan to factor you must leave room for the diagonal entry and
put in the entry even if it is zero.
- `o_nz`  - number of nonzeros per row in the OFF-DIAGONAL portion of local
submatrix (same value is used for all local rows).
- `o_nnz` - array containing the number of nonzeros in the various rows of the
OFF-DIAGONAL portion of the local submatrix (possibly different for
each row) or `NULL`, if `o_nz` is used to specify the nonzero
structure. The size of this array is equal to the number
of local rows, i.e `m`.

Output Parameter:
- `A` - the matrix

Options Database Keys:
- `-mat_no_inode`            - Do not use inodes
- `-mat_inode_limit <limit>` - Sets inode limit (max limit=5)

Level: intermediate

-seealso: [](ch_matrices), `Mat`, [Sparse Matrix Creation](sec_matsparse), `MATMPIAIJPERM`, `MatCreate()`, `MatCreateSeqAIJPERM()`, `MatSetValues()`

# External Links
$(_doc_external("Mat/MatCreateMPIAIJPERM"))
"""
function MatCreateMPIAIJPERM(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_nz::PetscInt, d_nnz::Union{Ptr, Vector{PetscInt}}, o_nz::PetscInt, o_nnz::Union{Ptr, Vector{PetscInt}}) end

@for_petsc function MatCreateMPIAIJPERM(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, d_nz::$PetscInt, d_nnz::Union{Ptr, Vector{$PetscInt}}, o_nz::$PetscInt, o_nnz::Union{Ptr, Vector{$PetscInt}} )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateMPIAIJPERM, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, Ptr{CMat}),
               comm, m, n, M, N, d_nz, d_nnz, o_nz, o_nnz, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	A::PetscMat = MatCreateMPIAIJSELL(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_nz::PetscInt, d_nnz::Vector{PetscInt}, o_nz::PetscInt, o_nnz::Vector{PetscInt}) 
Creates a sparse parallel matrix whose local
portions are stored as `MATSEQAIJSELL` matrices (a matrix class that inherits
from SEQAIJ but performs some operations in SELL format).

Collective

Input Parameters:
- `comm`  - MPI communicator
- `m`     - number of local rows (or `PETSC_DECIDE` to have calculated if `M` is given)
This value should be the same as the local size used in creating the
y vector for the matrix-vector product y = Ax.
- `n`     - This value should be the same as the local size used in creating the
x vector for the matrix-vector product y = Ax. (or `PETSC_DECIDE` to have
calculated if `N` is given) For square matrices `n` is almost always `m`.
- `M`     - number of global rows (or `PETSC_DETERMINE` to have calculated if `m` is given)
- `N`     - number of global columns (or `PETSC_DETERMINE` to have calculated if `n` is given)
- `d_nz`  - number of nonzeros per row in DIAGONAL portion of local submatrix
(same value is used for all local rows)
- `d_nnz` - array containing the number of nonzeros in the various rows of the
DIAGONAL portion of the local submatrix (possibly different for each row)
or `NULL`, if `d_nz` is used to specify the nonzero structure.
The size of this array is equal to the number of local rows, i.e `m`.
For matrices you plan to factor you must leave room for the diagonal entry and
put in the entry even if it is zero.
- `o_nz`  - number of nonzeros per row in the OFF-DIAGONAL portion of local
submatrix (same value is used for all local rows).
- `o_nnz` - array containing the number of nonzeros in the various rows of the
OFF-DIAGONAL portion of the local submatrix (possibly different for
each row) or `NULL`, if `o_nz` is used to specify the nonzero
structure. The size of this array is equal to the number
of local rows, i.e `m`.

Output Parameter:
- `A` - the matrix

Options Database Key:
- `-mat_aijsell_eager_shadow` - Construct shadow matrix upon matrix assembly; default is to take a "lazy" approach, performing this step the first
time the matrix is applied

Level: intermediate

-seealso: [](ch_matrices), `Mat`, [Sparse Matrix Creation](sec_matsparse), `MATSEQAIJSELL`, `MATMPIAIJSELL`, `MATAIJSELL`, `MatCreate()`, `MatCreateSeqAIJSELL()`, `MatSetValues()`

# External Links
$(_doc_external("Mat/MatCreateMPIAIJSELL"))
"""
function MatCreateMPIAIJSELL(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_nz::PetscInt, d_nnz::Vector{PetscInt}, o_nz::PetscInt, o_nnz::Vector{PetscInt}) end

@for_petsc function MatCreateMPIAIJSELL(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, d_nz::$PetscInt, d_nnz::Vector{$PetscInt}, o_nz::$PetscInt, o_nnz::Vector{$PetscInt} )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateMPIAIJSELL, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, Ptr{CMat}),
               comm, m, n, M, N, d_nz, d_nnz, o_nz, o_nnz, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	mpimat::PetscMat = MatCreateMPIAIJSumSeqAIJ(petsclib::PetscLibType,comm::MPI_Comm, seqmat::AbstractPetscMat, m::PetscInt, n::PetscInt, scall::MatReuse) 
Creates a `MATMPIAIJ` matrix by adding sequential
matrices from each processor

Collective

Input Parameters:
- `comm`   - the communicators the parallel matrix will live on
- `seqmat` - the input sequential matrices
- `m`      - number of local rows (or `PETSC_DECIDE`)
- `n`      - number of local columns (or `PETSC_DECIDE`)
- `scall`  - either `MAT_INITIAL_MATRIX` or `MAT_REUSE_MATRIX`

Output Parameter:
- `mpimat` - the parallel matrix generated

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatCreateAIJ()`

# External Links
$(_doc_external("Mat/MatCreateMPIAIJSumSeqAIJ"))
"""
function MatCreateMPIAIJSumSeqAIJ(petsclib::PetscLibType, comm::MPI_Comm, seqmat::AbstractPetscMat, m::PetscInt, n::PetscInt, scall::MatReuse) end

@for_petsc function MatCreateMPIAIJSumSeqAIJ(petsclib::$UnionPetscLib, comm::MPI_Comm, seqmat::AbstractPetscMat, m::$PetscInt, n::$PetscInt, scall::MatReuse )
	mpimat_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateMPIAIJSumSeqAIJ, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, CMat, $PetscInt, $PetscInt, MatReuse, Ptr{CMat}),
               comm, seqmat, m, n, scall, mpimat_,
              )

	mpimat = PetscMat(mpimat_[], petsclib)

	return mpimat
end 

"""
	MatCreateMPIAIJSumSeqAIJNumeric(petsclib::PetscLibType,seqmat::AbstractPetscMat, mpimat::AbstractPetscMat) 

# External Links
$(_doc_external("Mat/MatCreateMPIAIJSumSeqAIJNumeric"))
"""
function MatCreateMPIAIJSumSeqAIJNumeric(petsclib::PetscLibType, seqmat::AbstractPetscMat, mpimat::AbstractPetscMat) end

@for_petsc function MatCreateMPIAIJSumSeqAIJNumeric(petsclib::$UnionPetscLib, seqmat::AbstractPetscMat, mpimat::AbstractPetscMat )

    @chk ccall(
               (:MatCreateMPIAIJSumSeqAIJNumeric, $petsc_library),
               PetscErrorCode,
               (CMat, CMat),
               seqmat, mpimat,
              )


	return nothing
end 

"""
	mpimat::PetscMat = MatCreateMPIAIJSumSeqAIJSymbolic(petsclib::PetscLibType,comm::MPI_Comm, seqmat::AbstractPetscMat, m::PetscInt, n::PetscInt) 

# External Links
$(_doc_external("Mat/MatCreateMPIAIJSumSeqAIJSymbolic"))
"""
function MatCreateMPIAIJSumSeqAIJSymbolic(petsclib::PetscLibType, comm::MPI_Comm, seqmat::AbstractPetscMat, m::PetscInt, n::PetscInt) end

@for_petsc function MatCreateMPIAIJSumSeqAIJSymbolic(petsclib::$UnionPetscLib, comm::MPI_Comm, seqmat::AbstractPetscMat, m::$PetscInt, n::$PetscInt )
	mpimat_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateMPIAIJSumSeqAIJSymbolic, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, CMat, $PetscInt, $PetscInt, Ptr{CMat}),
               comm, seqmat, m, n, mpimat_,
              )

	mpimat = PetscMat(mpimat_[], petsclib)

	return mpimat
end 

"""
	mat::PetscMat = MatCreateMPIAIJWithArrays(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, a::Vector{PetscScalar}) 
creates a `MATMPIAIJ` matrix using arrays that contain in standard
CSR format for the local rows.

Collective

Input Parameters:
- `comm` - MPI communicator
- `m`    - number of local rows (Cannot be `PETSC_DECIDE`)
- `n`    - This value should be the same as the local size used in creating the
x vector for the matrix-vector product  y = Ax. (or `PETSC_DECIDE` to have
calculated if `N` is given) For square matrices n is almost always `m`.
- `M`    - number of global rows (or `PETSC_DETERMINE` to have calculated if `m` is given)
- `N`    - number of global columns (or `PETSC_DETERMINE` to have calculated if `n` is given)
- `i`    - row indices (of length m+1); that is i[0] = 0, i[row] = i[row-1] + number of elements in that row of the matrix
- `j`    - global column indices
- `a`    - optional matrix values

Output Parameter:
- `mat` - the matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `MatCreateSeqAIJ()`, `MatSetValues()`, `MatMPIAIJSetPreallocation()`, `MatMPIAIJSetPreallocationCSR()`,
`MATMPIAIJ`, `MatCreateAIJ()`, `MatCreateMPIAIJWithSplitArrays()`, `MatUpdateMPIAIJWithArray()`, `MatSetPreallocationCOO()`, `MatSetValuesCOO()`

# External Links
$(_doc_external("Mat/MatCreateMPIAIJWithArrays"))
"""
function MatCreateMPIAIJWithArrays(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, a::Vector{PetscScalar}) end

@for_petsc function MatCreateMPIAIJWithArrays(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, i::Vector{$PetscInt}, j::Vector{$PetscInt}, a::Vector{$PetscScalar} )
	mat_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateMPIAIJWithArrays, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}, Ptr{CMat}),
               comm, m, n, M, N, i, j, a, mat_,
              )

	mat = PetscMat(mat_[], petsclib)

	return mat
end 

"""
	garray::PetscInt,mat::PetscMat = MatCreateMPIAIJWithSeqAIJ(petsclib::PetscLibType,comm::MPI_Comm, M::PetscInt, N::PetscInt, A::AbstractPetscMat, B::AbstractPetscMat) 
creates a `MATMPIAIJ` matrix using `MATSEQAIJ` matrices that contain the "diagonal"
and "off-diagonal" part of the matrix in CSR format.

Collective

Input Parameters:
- `comm`   - MPI communicator
- `M`      - the global row size
- `N`      - the global column size
- `A`      - "diagonal" portion of matrix
- `B`      - if garray is `NULL`, B should be the offdiag matrix using global col ids and of size N - if garray is not `NULL`, B should be the offdiag matrix using local col ids and of size garray
- `garray` - either `NULL` or the global index of `B` columns. If not `NULL`, it should be allocated by `PetscMalloc1()` and will be owned by `mat` thereafter.

Output Parameter:
- `mat` - the matrix, with input `A` as its local diagonal matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATMPIAIJ`, `MATSEQAIJ`, `MatCreateMPIAIJWithSplitArrays()`

# External Links
$(_doc_external("Mat/MatCreateMPIAIJWithSeqAIJ"))
"""
function MatCreateMPIAIJWithSeqAIJ(petsclib::PetscLibType, comm::MPI_Comm, M::PetscInt, N::PetscInt, A::AbstractPetscMat, B::AbstractPetscMat) end

@for_petsc function MatCreateMPIAIJWithSeqAIJ(petsclib::$UnionPetscLib, comm::MPI_Comm, M::$PetscInt, N::$PetscInt, A::AbstractPetscMat, B::AbstractPetscMat )
	garray_ = Ref{$PetscInt}()
	mat_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateMPIAIJWithSeqAIJ, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, CMat, CMat, Ptr{$PetscInt}, Ptr{CMat}),
               comm, M, N, A, B, garray_, mat_,
              )

	garray = garray_[]
	mat = PetscMat(mat_[], petsclib)

	return garray,mat
end 

"""
	mat::PetscMat = MatCreateMPIAIJWithSplitArrays(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, a::Vector{PetscScalar}, oi::Vector{PetscInt}, oj::Vector{PetscInt}, oa::Vector{PetscScalar}) 
creates a `MATMPIAIJ` matrix using arrays that contain the "diagonal"
and "off-diagonal" part of the matrix in CSR format.

Collective

Input Parameters:
- `comm` - MPI communicator
- `m`    - number of local rows (Cannot be `PETSC_DECIDE`)
- `n`    - This value should be the same as the local size used in creating the
x vector for the matrix-vector product y = Ax. (or `PETSC_DECIDE` to have
calculated if `N` is given) For square matrices `n` is almost always `m`.
- `M`    - number of global rows (or `PETSC_DETERMINE` to have calculated if `m` is given)
- `N`    - number of global columns (or `PETSC_DETERMINE` to have calculated if `n` is given)
- `i`    - row indices for "diagonal" portion of matrix; that is i[0] = 0, i[row] = i[row-1] + number of elements in that row of the matrix
- `j`    - column indices, which must be local, i.e., based off the start column of the diagonal portion
- `a`    - matrix values
- `oi`   - row indices for "off-diagonal" portion of matrix; that is oi[0] = 0, oi[row] = oi[row-1] + number of elements in that row of the matrix
- `oj`   - column indices, which must be global, representing global columns in the `MATMPIAIJ` matrix
- `oa`   - matrix values

Output Parameter:
- `mat` - the matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `MatCreateSeqAIJ()`, `MatSetValues()`, `MatMPIAIJSetPreallocation()`, `MatMPIAIJSetPreallocationCSR()`,
`MATMPIAIJ`, `MatCreateAIJ()`, `MatCreateMPIAIJWithArrays()`

# External Links
$(_doc_external("Mat/MatCreateMPIAIJWithSplitArrays"))
"""
function MatCreateMPIAIJWithSplitArrays(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, a::Vector{PetscScalar}, oi::Vector{PetscInt}, oj::Vector{PetscInt}, oa::Vector{PetscScalar}) end

@for_petsc function MatCreateMPIAIJWithSplitArrays(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, i::Vector{$PetscInt}, j::Vector{$PetscInt}, a::Vector{$PetscScalar}, oi::Vector{$PetscInt}, oj::Vector{$PetscInt}, oa::Vector{$PetscScalar} )
	mat_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateMPIAIJWithSplitArrays, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}, Ptr{CMat}),
               comm, m, n, M, N, i, j, a, oi, oj, oa, mat_,
              )

	mat = PetscMat(mat_[], petsclib)

	return mat
end 

"""
	A::PetscMat = MatCreateMPIAdj(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, N::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, values::Union{Ptr, Vector{PetscInt}}) 
Creates a sparse matrix representing an adjacency list.
The matrix need not have numerical values associated with it, it is
intended for ordering (to reduce bandwidth etc) and partitioning.

Collective

Input Parameters:
- `comm`   - MPI communicator
- `m`      - number of local rows
- `N`      - number of global columns
- `i`      - the indices into `j` for the start of each row
- `j`      - the column indices for each row (sorted for each row).
- `values` - the values, optional, use `NULL` if not provided

Output Parameter:
- `A` - the matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `MatConvert()`, `MatGetOrdering()`, `MATMPIADJ`, `MatMPIAdjSetPreallocation()`

# External Links
$(_doc_external("Mat/MatCreateMPIAdj"))
"""
function MatCreateMPIAdj(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, N::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, values::Union{Ptr, Vector{PetscInt}}) end

@for_petsc function MatCreateMPIAdj(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, N::$PetscInt, i::Vector{$PetscInt}, j::Vector{$PetscInt}, values::Union{Ptr, Vector{$PetscInt}} )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateMPIAdj, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{CMat}),
               comm, m, N, i, j, values, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	mat::PetscMat = MatCreateMPIBAIJWithArrays(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, a::Vector{PetscScalar}) 
creates a `MATMPIBAIJ` matrix using arrays that contain in standard block CSR format for the local rows.

Collective

Input Parameters:
- `comm` - MPI communicator
- `bs`   - the block size, only a block size of 1 is supported
- `m`    - number of local rows (Cannot be `PETSC_DECIDE`)
- `n`    - This value should be the same as the local size used in creating the
x vector for the matrix-vector product  y = Ax . (or `PETSC_DECIDE` to have
calculated if `N` is given) For square matrices `n` is almost always `m`.
- `M`    - number of global rows (or `PETSC_DETERMINE` to have calculated if `m` is given)
- `N`    - number of global columns (or `PETSC_DETERMINE` to have calculated if `n` is given)
- `i`    - row indices; that is i[0] = 0, i[row] = i[row-1] + number of block elements in that rowth block row of the matrix
- `j`    - column indices
- `a`    - matrix values

Output Parameter:
- `mat` - the matrix

Level: intermediate

-seealso: `Mat`, `MatCreate()`, `MatCreateSeqAIJ()`, `MatSetValues()`, `MatMPIAIJSetPreallocation()`, `MatMPIAIJSetPreallocationCSR()`,
`MATMPIAIJ`, `MatCreateAIJ()`, `MatCreateMPIAIJWithSplitArrays()`

# External Links
$(_doc_external("Mat/MatCreateMPIBAIJWithArrays"))
"""
function MatCreateMPIBAIJWithArrays(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, a::Vector{PetscScalar}) end

@for_petsc function MatCreateMPIBAIJWithArrays(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, i::Vector{$PetscInt}, j::Vector{$PetscInt}, a::Vector{$PetscScalar} )
	mat_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateMPIBAIJWithArrays, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}, Ptr{CMat}),
               comm, bs, m, n, M, N, i, j, a, mat_,
              )

	mat = PetscMat(mat_[], petsclib)

	return mat
end 

"""
	mpimat::PetscMat = MatCreateMPIMatConcatenateSeqMat(petsclib::PetscLibType,comm::MPI_Comm, seqmat::AbstractPetscMat, n::PetscInt, reuse::MatReuse) 
Creates a single large PETSc matrix by concatenating sequential
matrices from each processor

Collective

Input Parameters:
- `comm`   - the communicators the parallel matrix will live on
- `seqmat` - the input sequential matrices
- `n`      - number of local columns (or `PETSC_DECIDE`)
- `reuse`  - either `MAT_INITIAL_MATRIX` or `MAT_REUSE_MATRIX`

Output Parameter:
- `mpimat` - the parallel matrix generated

Level: developer

-seealso: [](ch_matrices), `Mat`

# External Links
$(_doc_external("Mat/MatCreateMPIMatConcatenateSeqMat"))
"""
function MatCreateMPIMatConcatenateSeqMat(petsclib::PetscLibType, comm::MPI_Comm, seqmat::AbstractPetscMat, n::PetscInt, reuse::MatReuse) end

@for_petsc function MatCreateMPIMatConcatenateSeqMat(petsclib::$UnionPetscLib, comm::MPI_Comm, seqmat::AbstractPetscMat, n::$PetscInt, reuse::MatReuse )
	mpimat_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateMPIMatConcatenateSeqMat, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, CMat, $PetscInt, MatReuse, Ptr{CMat}),
               comm, seqmat, n, reuse, mpimat_,
              )

	mpimat = PetscMat(mpimat_[], petsclib)

	return mpimat
end 

"""
	mat::PetscMat = MatCreateMPISBAIJWithArrays(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, a::Vector{PetscScalar}) 
creates a `MATMPISBAIJ` matrix using arrays that contain in standard CSR format for the local rows.

Collective

Input Parameters:
- `comm` - MPI communicator
- `bs`   - the block size, only a block size of 1 is supported
- `m`    - number of local rows (Cannot be `PETSC_DECIDE`)
- `n`    - This value should be the same as the local size used in creating the
x vector for the matrix-vector product  y = Ax . (or `PETSC_DECIDE` to have
calculated if `N` is given) For square matrices `n` is almost always `m`.
- `M`    - number of global rows (or `PETSC_DETERMINE` to have calculated if `m` is given)
- `N`    - number of global columns (or `PETSC_DETERMINE` to have calculated if `n` is given)
- `i`    - row indices; that is i[0] = 0, i[row] = i[row-1] + number of block elements in that row block row of the matrix
- `j`    - column indices
- `a`    - matrix values

Output Parameter:
- `mat` - the matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATMPISBAIJ`, `MatCreate()`, `MatCreateSeqAIJ()`, `MatSetValues()`, `MatMPIAIJSetPreallocation()`, `MatMPIAIJSetPreallocationCSR()`,
`MATMPIAIJ`, `MatCreateAIJ()`, `MatCreateMPIAIJWithSplitArrays()`, `MatMPISBAIJSetPreallocationCSR()`

# External Links
$(_doc_external("Mat/MatCreateMPISBAIJWithArrays"))
"""
function MatCreateMPISBAIJWithArrays(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, a::Vector{PetscScalar}) end

@for_petsc function MatCreateMPISBAIJWithArrays(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, i::Vector{$PetscInt}, j::Vector{$PetscInt}, a::Vector{$PetscScalar} )
	mat_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateMPISBAIJWithArrays, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}, Ptr{CMat}),
               comm, bs, m, n, M, N, i, j, a, mat_,
              )

	mat = PetscMat(mat_[], petsclib)

	return mat
end 

"""
	B::PetscMat = MatCreateNest(petsclib::PetscLibType,comm::MPI_Comm, nr::PetscInt, is_row::Vector{<:AbstractIS}, nc::PetscInt, is_col::Vector{<:AbstractIS}, a::Vector{<:AbstractPetscMat}) 
Creates a new `MATNEST` matrix containing several nested submatrices, each stored separately

Collective

Input Parameters:
- `comm`   - Communicator for the new `MATNEST`
- `nr`     - number of nested row blocks
- `is_row` - index sets for each nested row block, or `NULL` to make contiguous
- `nc`     - number of nested column blocks
- `is_col` - index sets for each nested column block, or `NULL` to make contiguous
- `a`      - array of nr \\times nc submatrices, empty submatrices can be passed using `NULL`

Output Parameter:
- `B` - new matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATNEST`, `MatCreate()`, `VecCreateNest()`, `DMCreateMatrix()`, `MatNestSetSubMat()`,
`MatNestGetSubMat()`, `MatNestGetLocalISs()`, `MatNestGetSize()`,
`MatNestGetISs()`, `MatNestSetSubMats()`, `MatNestGetSubMats()`

# External Links
$(_doc_external("Mat/MatCreateNest"))
"""
function MatCreateNest(petsclib::PetscLibType, comm::MPI_Comm, nr::PetscInt, is_row::Vector{<:AbstractIS}, nc::PetscInt, is_col::Vector{<:AbstractIS}, a::Vector{<:AbstractPetscMat}) end

@for_petsc function MatCreateNest(petsclib::$UnionPetscLib, comm::MPI_Comm, nr::$PetscInt, is_row::Vector{<:AbstractIS}, nc::$PetscInt, is_col::Vector{<:AbstractIS}, a::Vector{<:AbstractPetscMat} )
	B_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateNest, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, Ptr{CIS}, $PetscInt, Ptr{CIS}, Ptr{CMat}, Ptr{CMat}),
               comm, nr, is_row, nc, is_col, a, B_,
              )

	B = PetscMat(B_[], petsclib)

	return B
end 

"""
	N::PetscMat = MatCreateNormal(petsclib::PetscLibType,A::AbstractPetscMat) 
Creates a new `MATNORMAL` matrix object that behaves like A^T A.

Collective

Input Parameter:
- `A` - the (possibly rectangular) matrix

Output Parameter:
- `N` - the matrix that represents A^T A 

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATNORMAL`, `MatMult()`, `MatNormalGetMat()`, `MATNORMALHERMITIAN`, `MatCreateNormalHermitian()`

# External Links
$(_doc_external("Mat/MatCreateNormal"))
"""
function MatCreateNormal(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatCreateNormal(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	N_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateNormal, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, N_,
              )

	N = PetscMat(N_[], petsclib)

	return N
end 

"""
	N::PetscMat = MatCreateNormalHermitian(petsclib::PetscLibType,A::AbstractPetscMat) 
Creates a new matrix object `MATNORMALHERMITIAN` that behaves like A^* A.

Collective

Input Parameter:
- `A` - the (possibly rectangular complex) matrix

Output Parameter:
- `N` - the matrix that represents  A^* A

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATNORMAL`, `MATNORMALHERMITIAN`, `MatNormalHermitianGetMat()`

# External Links
$(_doc_external("Mat/MatCreateNormalHermitian"))
"""
function MatCreateNormalHermitian(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatCreateNormalHermitian(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	N_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateNormalHermitian, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, N_,
              )

	N = PetscMat(N_[], petsclib)

	return N
end 

"""
	matredundant::PetscMat = MatCreateRedundantMatrix(petsclib::PetscLibType,mat::AbstractPetscMat, nsubcomm::PetscInt, subcomm::MPI_Comm, reuse::MatReuse) 
Create redundant matrices and put them into processors of subcommunicators.

Collective

Input Parameters:
- `mat`      - the matrix
- `nsubcomm` - the number of subcommunicators (= number of redundant parallel or sequential matrices)
- `subcomm`  - MPI communicator split from the communicator where mat resides in (or `MPI_COMM_NULL` if nsubcomm is used)
- `reuse`    - either `MAT_INITIAL_MATRIX` or `MAT_REUSE_MATRIX`

Output Parameter:
- `matredundant` - redundant matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatDestroy()`, `PetscSubcommCreate()`, `PetscSubcomm`

# External Links
$(_doc_external("Mat/MatCreateRedundantMatrix"))
"""
function MatCreateRedundantMatrix(petsclib::PetscLibType, mat::AbstractPetscMat, nsubcomm::PetscInt, subcomm::MPI_Comm, reuse::MatReuse) end

@for_petsc function MatCreateRedundantMatrix(petsclib::$UnionPetscLib, mat::AbstractPetscMat, nsubcomm::$PetscInt, subcomm::MPI_Comm, reuse::MatReuse )
	matredundant_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateRedundantMatrix, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, MPI_Comm, MatReuse, Ptr{CMat}),
               mat, nsubcomm, subcomm, reuse, matredundant_,
              )

	matredundant = PetscMat(matredundant_[], petsclib)

	return matredundant
end 

"""
	A::PetscMat = MatCreateSBAIJ(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_nz::PetscInt, d_nnz::Union{Ptr, Vector{PetscInt}}, o_nz::PetscInt, o_nnz::Union{Ptr, Vector{PetscInt}}) 
Creates a sparse parallel matrix in symmetric block AIJ format, `MATSBAIJ`,
(block compressed row).  For good matrix assembly performance
the user should preallocate the matrix storage by setting the parameters
`d_nz` (or `d_nnz`) and `o_nz` (or `o_nnz`).

Collective

Input Parameters:
- `comm`  - MPI communicator
- `bs`    - size of block, the blocks are ALWAYS square. One can use `MatSetBlockSizes()` to set a different row and column blocksize but the row
blocksize always defines the size of the blocks. The column blocksize sets the blocksize of the vectors obtained with `MatCreateVecs()`
- `m`     - number of local rows (or `PETSC_DECIDE` to have calculated if `M` is given)
This value should be the same as the local size used in creating the
y vector for the matrix-vector product y = Ax.
- `n`     - number of local columns (or `PETSC_DECIDE` to have calculated if `N` is given)
This value should be the same as the local size used in creating the
x vector for the matrix-vector product y = Ax.
- `M`     - number of global rows (or `PETSC_DETERMINE` to have calculated if `m` is given)
- `N`     - number of global columns (or `PETSC_DETERMINE` to have calculated if `n` is given)
- `d_nz`  - number of block nonzeros per block row in diagonal portion of local
submatrix (same for all local rows)
- `d_nnz` - array containing the number of block nonzeros in the various block rows
in the upper triangular portion of the in diagonal portion of the local
(possibly different for each block block row) or `NULL`.
If you plan to factor the matrix you must leave room for the diagonal entry and
set its value even if it is zero.
- `o_nz`  - number of block nonzeros per block row in the off-diagonal portion of local
submatrix (same for all local rows).
- `o_nnz` - array containing the number of nonzeros in the various block rows of the
off-diagonal portion of the local submatrix (possibly different for
each block row) or `NULL`.

Output Parameter:
- `A` - the matrix

Options Database Keys:
- `-mat_no_unroll`  - uses code that does not unroll the loops in the
block calculations (much slower)
- `-mat_block_size` - size of the blocks to use
- `-mat_mpi`        - use the parallel matrix data structures even on one processor
(defaults to using SeqBAIJ format on one processor)

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATSBAIJ`, `MatCreate()`, `MatCreateSeqSBAIJ()`, `MatSetValues()`, `MatCreateBAIJ()`,
`MatGetOwnershipRange()`,  `MatGetOwnershipRanges()`, `MatGetOwnershipRangeColumn()`, `MatGetOwnershipRangesColumn()`, `PetscLayout`

# External Links
$(_doc_external("Mat/MatCreateSBAIJ"))
"""
function MatCreateSBAIJ(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_nz::PetscInt, d_nnz::Union{Ptr, Vector{PetscInt}}, o_nz::PetscInt, o_nnz::Union{Ptr, Vector{PetscInt}}) end

@for_petsc function MatCreateSBAIJ(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, d_nz::$PetscInt, d_nnz::Union{Ptr, Vector{$PetscInt}}, o_nz::$PetscInt, o_nnz::Union{Ptr, Vector{$PetscInt}} )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateSBAIJ, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, Ptr{CMat}),
               comm, bs, m, n, M, N, d_nz, d_nnz, o_nz, o_nnz, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	A::PetscMat = MatCreateSELL(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_rlenmax::PetscInt, d_rlen::Union{Ptr, Vector{PetscInt}}, o_rlenmax::PetscInt, o_rlen::Union{Ptr, Vector{PetscInt}}) 
Creates a sparse parallel matrix in `MATSELL` format.

Collective

Input Parameters:
- `comm`      - MPI communicator
- `m`         - number of local rows (or `PETSC_DECIDE` to have calculated if M is given)
This value should be the same as the local size used in creating the
y vector for the matrix-vector product y = Ax.
- `n`         - This value should be the same as the local size used in creating the
x vector for the matrix-vector product y = Ax. (or `PETSC_DECIDE` to have
calculated if `N` is given) For square matrices n is almost always `m`.
- `M`         - number of global rows (or `PETSC_DETERMINE` to have calculated if `m` is given)
- `N`         - number of global columns (or `PETSC_DETERMINE` to have calculated if `n` is given)
- `d_rlenmax` - max number of nonzeros per row in DIAGONAL portion of local submatrix
(same value is used for all local rows)
- `d_rlen`    - array containing the number of nonzeros in the various rows of the
DIAGONAL portion of the local submatrix (possibly different for each row)
or `NULL`, if d_rlenmax is used to specify the nonzero structure.
The size of this array is equal to the number of local rows, i.e `m`.
- `o_rlenmax` - max number of nonzeros per row in the OFF-DIAGONAL portion of local
submatrix (same value is used for all local rows).
- `o_rlen`    - array containing the number of nonzeros in the various rows of the
OFF-DIAGONAL portion of the local submatrix (possibly different for
each row) or `NULL`, if `o_rlenmax` is used to specify the nonzero
structure. The size of this array is equal to the number
of local rows, i.e `m`.

Output Parameter:
- `A` - the matrix

Options Database Key:
- `-mat_sell_oneindex` - Internally use indexing starting at 1
rather than 0.  When calling `MatSetValues()`,
the user still MUST index entries starting at 0!

Example:
Consider the following 8x8 matrix with 34 non-zero values, that is
assembled across 3 processors. Lets assume that proc0 owns 3 rows,
proc1 owns 3 rows, proc2 owns 2 rows. This division can be shown
as follows

-seealso: `Mat`, `MATSELL`, `MatCreate()`, `MatCreateSeqSELL()`, `MatSetValues()`, `MatMPISELLSetPreallocation()`, `MATMPISELL`

# External Links
$(_doc_external("Mat/MatCreateSELL"))
"""
function MatCreateSELL(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_rlenmax::PetscInt, d_rlen::Union{Ptr, Vector{PetscInt}}, o_rlenmax::PetscInt, o_rlen::Union{Ptr, Vector{PetscInt}}) end

@for_petsc function MatCreateSELL(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, d_rlenmax::$PetscInt, d_rlen::Union{Ptr, Vector{$PetscInt}}, o_rlenmax::$PetscInt, o_rlen::Union{Ptr, Vector{$PetscInt}} )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateSELL, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, Ptr{CMat}),
               comm, m, n, M, N, d_rlenmax, d_rlen, o_rlenmax, o_rlen, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	J::PetscMat = MatCreateSNESMF(petsclib::PetscLibType,snes::AbstractPetscSNES) 
Creates a finite differencing based matrix
a `SNES` solver.  This matrix can be used as the Jacobian argument for
the routine `SNESSetJacobian()`. See `MatCreateMFFD()` for details on how
the finite difference computation is done.

Collective

Input Parameters:
- `snes` - the `SNES` context

Output Parameter:
- `J` - the matrix-free matrix which is of type `MATMFFD`

Level: advanced

-seealso: [](ch_snes), `SNES`, `MATMFFD`, `MatDestroy()`, `MatMFFDSetFunction()`, `MatMFFDSetFunctionError()`, `MatMFFDDSSetUmin()`
`MatMFFDSetHHistory()`, `MatMFFDResetHHistory()`, `MatCreateMFFD()`, `MatCreateShell()`,
`MatMFFDGetH()`, `MatMFFDRegister()`, `MatMFFDComputeJacobian()`, `MatSNESMFSetReuseBase()`, `MatSNESMFGetReuseBase()`

# External Links
$(_doc_external("SNES/MatCreateSNESMF"))
"""
function MatCreateSNESMF(petsclib::PetscLibType, snes::AbstractPetscSNES) end

@for_petsc function MatCreateSNESMF(petsclib::$UnionPetscLib, snes::AbstractPetscSNES )
	J_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateSNESMF, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CMat}),
               snes, J_,
              )

	J = PetscMat(J_[], petsclib)

	return J
end 

"""
	J::PetscMat = MatCreateSNESMFMore(petsclib::PetscLibType,snes::AbstractPetscSNES, x::AbstractPetscVec) 
Creates a matrix
context for use with a `SNES` solver that uses the More method to compute an optimal h based on the noise of the function.  This matrix can be used as
the Jacobian argument for the routine `SNESSetJacobian()`.

Input Parameters:
- `snes` - the `SNES` context
- `x`    - vector where `SNES` solution is to be stored.

Output Parameter:
- `J` - the matrix-free matrix

Options Database Keys:
- `-snes_mf_err <error_rel>` - see `MatCreateSNESMF()`
- `-snes_mf_umin <umin>`     - see `MatCreateSNESMF()`
- `-snes_mf_compute_err`     - compute the square root or relative error in function
- `-snes_mf_freq_err <freq>` - set the frequency to recompute the parameters
- `-snes_mf_jorge`           - use the method of Jorge More

Level: advanced

-seealso: [](ch_snes), `SNESCreateMF()`, `MatCreateMFFD()`, `MatDestroy()`, `MatMFFDSetFunctionError()`

# External Links
$(_doc_external("SNES/MatCreateSNESMFMore"))
"""
function MatCreateSNESMFMore(petsclib::PetscLibType, snes::AbstractPetscSNES, x::AbstractPetscVec) end

@for_petsc function MatCreateSNESMFMore(petsclib::$UnionPetscLib, snes::AbstractPetscSNES, x::AbstractPetscVec )
	J_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateSNESMFMore, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, Ptr{CMat}),
               snes, x, J_,
              )

	J = PetscMat(J_[], petsclib)

	return J
end 

"""
	A::PetscMat = MatCreateScaLAPACK(petsclib::PetscLibType,comm::MPI_Comm, mb::PetscInt, nb::PetscInt, M::PetscInt, N::PetscInt, rsrc::PetscInt, csrc::PetscInt) 
Creates a dense parallel matrix in ScaLAPACK format
(2D block cyclic distribution) for a `MATSCALAPACK` matrix

Collective

Input Parameters:
- `comm` - MPI communicator
- `mb`   - row block size (or `PETSC_DECIDE` to have it set)
- `nb`   - column block size (or `PETSC_DECIDE` to have it set)
- `M`    - number of global rows
- `N`    - number of global columns
- `rsrc` - coordinate of process that owns the first row of the distributed matrix
- `csrc` - coordinate of process that owns the first column of the distributed matrix

Output Parameter:
- `A` - the matrix

Options Database Key:
- `-mat_scalapack_block_sizes` - size of the blocks to use (one or two integers separated by comma)

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATSCALAPACK`, `MATDENSE`, `MATELEMENTAL`, `MatCreate()`, `MatCreateDense()`, `MatSetValues()`

# External Links
$(_doc_external("Mat/MatCreateScaLAPACK"))
"""
function MatCreateScaLAPACK(petsclib::PetscLibType, comm::MPI_Comm, mb::PetscInt, nb::PetscInt, M::PetscInt, N::PetscInt, rsrc::PetscInt, csrc::PetscInt) end

@for_petsc function MatCreateScaLAPACK(petsclib::$UnionPetscLib, comm::MPI_Comm, mb::$PetscInt, nb::$PetscInt, M::$PetscInt, N::$PetscInt, rsrc::$PetscInt, csrc::$PetscInt )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateScaLAPACK, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{CMat}),
               comm, mb, nb, M, N, rsrc, csrc, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	A::PetscMat = MatCreateScatter(petsclib::PetscLibType,comm::MPI_Comm, scatter::VecScatter) 
Creates a new matrix of `MatType` `MATSCATTER`, based on a VecScatter

Collective

Input Parameters:
- `comm`    - MPI communicator
- `scatter` - a `VecScatter`

Output Parameter:
- `A` - the matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatScatterSetVecScatter()`, `MatScatterGetVecScatter()`, `MATSCATTER`

# External Links
$(_doc_external("Mat/MatCreateScatter"))
"""
function MatCreateScatter(petsclib::PetscLibType, comm::MPI_Comm, scatter::VecScatter) end

@for_petsc function MatCreateScatter(petsclib::$UnionPetscLib, comm::MPI_Comm, scatter::VecScatter )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateScatter, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, VecScatter, Ptr{CMat}),
               comm, scatter, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	S::PetscMat = MatCreateSchurComplement(petsclib::PetscLibType,A00::AbstractPetscMat, Ap00::AbstractPetscMat, A01::AbstractPetscMat, A10::AbstractPetscMat, A11::AbstractPetscMat) 
Creates a new `Mat` that behaves like the Schur complement of a matrix

Collective

Input Parameters:
- `A00`  - the upper-left block of the original matrix A = [A00 A01; A10 A11]
- `Ap00` - matrix from which the preconditioner is constructed for use in ksp(A00,Ap00) to approximate the action of A00^{-1}
- `A01`  - the upper-right block of the original matrix A = [A00 A01; A10 A11]
- `A10`  - the lower-left block of the original matrix A = [A00 A01; A10 A11]
- `A11`  - (optional) the lower-right block of the original matrix A = [A00 A01; A10 A11]

Output Parameter:
- `S` - the matrix that behaves as the Schur complement S = A11 - A10 ksp(A00,Ap00) A01

Level: intermediate

-seealso: [](ch_ksp), `MatCreateNormal()`, `MatMult()`, `MatCreate()`, `MatSchurComplementGetKSP()`, `MatSchurComplementUpdateSubMatrices()`, `MatCreateTranspose()`, `MatGetSchurComplement()`,
`MatSchurComplementGetPmat()`, `MatSchurComplementSetSubMatrices()`

# External Links
$(_doc_external("KSP/MatCreateSchurComplement"))
"""
function MatCreateSchurComplement(petsclib::PetscLibType, A00::AbstractPetscMat, Ap00::AbstractPetscMat, A01::AbstractPetscMat, A10::AbstractPetscMat, A11::AbstractPetscMat) end

@for_petsc function MatCreateSchurComplement(petsclib::$UnionPetscLib, A00::AbstractPetscMat, Ap00::AbstractPetscMat, A01::AbstractPetscMat, A10::AbstractPetscMat, A11::AbstractPetscMat )
	S_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateSchurComplement, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat, CMat, CMat, Ptr{CMat}),
               A00, Ap00, A01, A10, A11, S_,
              )

	S = PetscMat(S_[], petsclib)

	return S
end 

"""
	Sp::PetscMat = MatCreateSchurComplementPmat(petsclib::PetscLibType,A00::AbstractPetscMat, A01::AbstractPetscMat, A10::AbstractPetscMat, A11::AbstractPetscMat, ainvtype::MatSchurComplementAinvType, preuse::MatReuse) 
create a matrix for preconditioning the Schur complement by explicitly assembling the sparse matrix
Sp = A11 - A10 inv(DIAGFORM(A00)) A01

Collective

Input Parameters:
- `A00`      - the upper-left part of the original matrix A = [A00 A01; A10 A11]
- `A01`      - (optional) the upper-right part of the original matrix A = [A00 A01; A10 A11]
- `A10`      - (optional) the lower-left part of the original matrix A = [A00 A01; A10 A11]
- `A11`      - (optional) the lower-right part of the original matrix A = [A00 A01; A10 A11]
- `ainvtype` - type of approximation for DIAGFORM(A00) used when forming Sp = A11 - A10 inv(DIAGFORM(A00)) A01. See `MatSchurComplementAinvType`.
- `preuse`   - `MAT_INITIAL_MATRIX` for a new `Sp`, or `MAT_REUSE_MATRIX` to reuse an existing `Sp`, or `MAT_IGNORE_MATRIX` to put nothing in `Sp`

Output Parameter:
- `Sp` - approximate Schur complement suitable for constructing a preconditioner for the true Schur complement S = A11 - A10 inv(A00) A01

Level: advanced

-seealso: [](ch_ksp), `MatCreateSchurComplement()`, `MatGetSchurComplement()`, `MatSchurComplementGetPmat()`, `MatSchurComplementAinvType`

# External Links
$(_doc_external("KSP/MatCreateSchurComplementPmat"))
"""
function MatCreateSchurComplementPmat(petsclib::PetscLibType, A00::AbstractPetscMat, A01::AbstractPetscMat, A10::AbstractPetscMat, A11::AbstractPetscMat, ainvtype::MatSchurComplementAinvType, preuse::MatReuse) end

@for_petsc function MatCreateSchurComplementPmat(petsclib::$UnionPetscLib, A00::AbstractPetscMat, A01::AbstractPetscMat, A10::AbstractPetscMat, A11::AbstractPetscMat, ainvtype::MatSchurComplementAinvType, preuse::MatReuse )
	Sp_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateSchurComplementPmat, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat, CMat, MatSchurComplementAinvType, MatReuse, Ptr{CMat}),
               A00, A01, A10, A11, ainvtype, preuse, Sp_,
              )

	Sp = PetscMat(Sp_[], petsclib)

	return Sp
end 

"""
	A::PetscMat = MatCreateSeqAIJ(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Union{Ptr, Vector{PetscInt}}) 
Creates a sparse matrix in `MATSEQAIJ` (compressed row) format
(the default parallel PETSc format).  For good matrix assembly performance
the user should preallocate the matrix storage by setting the parameter `nz`
(or the array `nnz`).

Collective

Input Parameters:
- `comm` - MPI communicator, set to `PETSC_COMM_SELF`
- `m`    - number of rows
- `n`    - number of columns
- `nz`   - number of nonzeros per row (same for all rows)
- `nnz`  - array containing the number of nonzeros in the various rows
(possibly different for each row) or NULL

Output Parameter:
- `A` - the matrix

Options Database Keys:
- `-mat_no_inode`            - Do not use inodes
- `-mat_inode_limit <limit>` - Sets inode limit (max limit=5)

Level: intermediate

-seealso: [](ch_matrices), `Mat`, [Sparse Matrix Creation](sec_matsparse), `MatCreate()`, `MatCreateAIJ()`, `MatSetValues()`, `MatSeqAIJSetColumnIndices()`, `MatCreateSeqAIJWithArrays()`

# External Links
$(_doc_external("Mat/MatCreateSeqAIJ"))
"""
function MatCreateSeqAIJ(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Union{Ptr, Vector{PetscInt}}) end

@for_petsc function MatCreateSeqAIJ(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, nz::$PetscInt, nnz::Union{Ptr, Vector{$PetscInt}} )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateSeqAIJ, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{CMat}),
               comm, m, n, nz, nnz, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	A::PetscMat = MatCreateSeqAIJCRL(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Union{Ptr, Vector{PetscInt}}) 
Creates a sparse matrix of type `MATSEQAIJCRL`.

Collective

Input Parameters:
- `comm` - MPI communicator, set to `PETSC_COMM_SELF`
- `m`    - number of rows
- `n`    - number of columns
- `nz`   - number of nonzeros per row (same for all rows), ignored if `nnz` is given
- `nnz`  - array containing the number of nonzeros in the various rows
(possibly different for each row) or `NULL`

Output Parameter:
- `A` - the matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `MatCreateMPIAIJPERM()`, `MatSetValues()`

# External Links
$(_doc_external("Mat/MatCreateSeqAIJCRL"))
"""
function MatCreateSeqAIJCRL(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Union{Ptr, Vector{PetscInt}}) end

@for_petsc function MatCreateSeqAIJCRL(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, nz::$PetscInt, nnz::Union{Ptr, Vector{$PetscInt}} )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateSeqAIJCRL, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{CMat}),
               comm, m, n, nz, nnz, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	mat::PetscMat = MatCreateSeqAIJFromTriple(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, a::Vector{PetscScalar}, nz::PetscCount, idx::PetscBool) 
Creates an sequential `MATSEQAIJ` matrix using matrix elements (in COO format)
provided by the user.

Collective

Input Parameters:
- `comm` - must be an MPI communicator of size 1
- `m`    - number of rows
- `n`    - number of columns
- `i`    - row indices
- `j`    - column indices
- `a`    - matrix values
- `nz`   - number of nonzeros
- `idx`  - if the `i` and `j` indices start with 1 use `PETSC_TRUE` otherwise use `PETSC_FALSE`

Output Parameter:
- `mat` - the matrix

Level: intermediate

Example:
For the following matrix, the input data expected is as shown (using 0 based indexing)
-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `MatCreateAIJ()`, `MatCreateSeqAIJ()`, `MatCreateSeqAIJWithArrays()`, `MatMPIAIJSetPreallocationCSR()`, `MatSetValuesCOO()`, `MatSetPreallocationCOO()`

# External Links
$(_doc_external("Mat/MatCreateSeqAIJFromTriple"))
"""
function MatCreateSeqAIJFromTriple(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, a::Vector{PetscScalar}, nz::PetscCount, idx::PetscBool) end

@for_petsc function MatCreateSeqAIJFromTriple(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, i::Vector{$PetscInt}, j::Vector{$PetscInt}, a::Vector{$PetscScalar}, nz::PetscCount, idx::PetscBool )
	mat_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateSeqAIJFromTriple, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}, Ptr{CMat}, PetscCount, PetscBool),
               comm, m, n, i, j, a, mat_, nz, idx,
              )

	mat = PetscMat(mat_[], petsclib)

	return mat
end 

"""
	A::PetscMat = MatCreateSeqAIJKokkos(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Union{Ptr, Vector{PetscInt}}) 

# External Links
$(_doc_external("Mat/MatCreateSeqAIJKokkos"))
"""
function MatCreateSeqAIJKokkos(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Union{Ptr, Vector{PetscInt}}) end

@for_petsc function MatCreateSeqAIJKokkos(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, nz::$PetscInt, nnz::Union{Ptr, Vector{$PetscInt}} )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateSeqAIJKokkos, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{CMat}),
               comm, m, n, nz, nnz, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	A::PetscMat = MatCreateSeqAIJMKL(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Vector{PetscInt}) 
Creates a sparse matrix of type `MATSEQAIJMKL`.

Collective

Input Parameters:
- `comm` - MPI communicator, set to `PETSC_COMM_SELF`
- `m`    - number of rows
- `n`    - number of columns
- `nz`   - number of nonzeros per row (same for all rows)
- `nnz`  - array containing the number of nonzeros in the various rows
(possibly different for each row) or `NULL`

Output Parameter:
- `A` - the matrix

Options Database Keys:
- `-mat_aijmkl_no_spmv2`         - disable use of the SpMV2 inspector-executor routines
- `-mat_aijmkl_eager_inspection` - perform MKL "inspection" phase upon matrix assembly; default is to do "lazy" inspection,
performing this step the first time the matrix is applied

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `MatCreateMPIAIJMKL()`, `MatSetValues()`

# External Links
$(_doc_external("Mat/MatCreateSeqAIJMKL"))
"""
function MatCreateSeqAIJMKL(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Vector{PetscInt}) end

@for_petsc function MatCreateSeqAIJMKL(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, nz::$PetscInt, nnz::Vector{$PetscInt} )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateSeqAIJMKL, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{CMat}),
               comm, m, n, nz, nnz, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	A::PetscMat = MatCreateSeqAIJPERM(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Vector{PetscInt}) 
Creates a sparse matrix of type `MATSEQAIJPERM`.

Collective

Input Parameters:
- `comm` - MPI communicator, set to `PETSC_COMM_SELF`
- `m`    - number of rows
- `n`    - number of columns
- `nz`   - number of nonzeros per row (same for all rows), ignored if `nnz` is given
- `nnz`  - array containing the number of nonzeros in the various rows (possibly different for each row) or `NULL`

Output Parameter:
- `A` - the matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `MatCreateMPIAIJPERM()`, `MatSetValues()`

# External Links
$(_doc_external("Mat/MatCreateSeqAIJPERM"))
"""
function MatCreateSeqAIJPERM(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Vector{PetscInt}) end

@for_petsc function MatCreateSeqAIJPERM(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, nz::$PetscInt, nnz::Vector{$PetscInt} )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateSeqAIJPERM, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{CMat}),
               comm, m, n, nz, nnz, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	A::PetscMat = MatCreateSeqAIJSELL(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Union{Ptr, Vector{PetscInt}}) 
Creates a sparse matrix of type `MATSEQAIJSELL`.

Collective

Input Parameters:
- `comm` - MPI communicator, set to `PETSC_COMM_SELF`
- `m`    - number of rows
- `n`    - number of columns
- `nz`   - number of nonzeros per row (same for all rows)
- `nnz`  - array containing the number of nonzeros in the various rows
(possibly different for each row) or `NULL`

Output Parameter:
- `A` - the matrix

Options Database Keys:
- `-mat_aijsell_eager_shadow` - Construct shadow matrix upon matrix assembly; default is to take a "lazy" approach,
performing this step the first time the matrix is applied

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `MatCreateMPIAIJSELL()`, `MatSetValues()`

# External Links
$(_doc_external("Mat/MatCreateSeqAIJSELL"))
"""
function MatCreateSeqAIJSELL(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Union{Ptr, Vector{PetscInt}}) end

@for_petsc function MatCreateSeqAIJSELL(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, nz::$PetscInt, nnz::Union{Ptr, Vector{$PetscInt}} )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateSeqAIJSELL, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{CMat}),
               comm, m, n, nz, nnz, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	A::PetscMat = MatCreateSeqAIJViennaCL(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Union{Ptr, Vector{PetscInt}}) 

# External Links
$(_doc_external("Mat/MatCreateSeqAIJViennaCL"))
"""
function MatCreateSeqAIJViennaCL(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Union{Ptr, Vector{PetscInt}}) end

@for_petsc function MatCreateSeqAIJViennaCL(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, nz::$PetscInt, nnz::Union{Ptr, Vector{$PetscInt}} )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateSeqAIJViennaCL, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{CMat}),
               comm, m, n, nz, nnz, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	mat::PetscMat = MatCreateSeqAIJWithArrays(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, a::Vector{PetscScalar}) 
Creates an sequential `MATSEQAIJ` matrix using matrix elements (in CSR format)
provided by the user.

Collective

Input Parameters:
- `comm` - must be an MPI communicator of size 1
- `m`    - number of rows
- `n`    - number of columns
- `i`    - row indices; that is i[0] = 0, i[row] = i[row-1] + number of elements in that row of the matrix
- `j`    - column indices
- `a`    - matrix values

Output Parameter:
- `mat` - the matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `MatCreateAIJ()`, `MatCreateSeqAIJ()`, `MatCreateMPIAIJWithArrays()`, `MatMPIAIJSetPreallocationCSR()`

# External Links
$(_doc_external("Mat/MatCreateSeqAIJWithArrays"))
"""
function MatCreateSeqAIJWithArrays(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, a::Vector{PetscScalar}) end

@for_petsc function MatCreateSeqAIJWithArrays(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, i::Vector{$PetscInt}, j::Vector{$PetscInt}, a::Vector{$PetscScalar} )
	mat_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateSeqAIJWithArrays, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}, Ptr{CMat}),
               comm, m, n, i, j, a, mat_,
              )

	mat = PetscMat(mat_[], petsclib)

	return mat
end 

"""
	A::PetscMat = MatCreateSeqBAIJ(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Union{Ptr, Vector{PetscInt}}) 
Creates a sparse matrix in `MATSEQAIJ` (block
compressed row) format.  For good matrix assembly performance the
user should preallocate the matrix storage by setting the parameter `nz`
(or the array `nnz`).

Collective

Input Parameters:
- `comm` - MPI communicator, set to `PETSC_COMM_SELF`
- `bs`   - size of block, the blocks are ALWAYS square. One can use `MatSetBlockSizes()` to set a different row and column blocksize but the row
blocksize always defines the size of the blocks. The column blocksize sets the blocksize of the vectors obtained with `MatCreateVecs()`
- `m`    - number of rows
- `n`    - number of columns
- `nz`   - number of nonzero blocks  per block row (same for all rows)
- `nnz`  - array containing the number of nonzero blocks in the various block rows
(possibly different for each block row) or `NULL`

Output Parameter:
- `A` - the matrix

Options Database Keys:
- `-mat_no_unroll`  - uses code that does not unroll the loops in the block calculations (much slower)
- `-mat_block_size` - size of the blocks to use

Level: intermediate

-seealso: [](ch_matrices), `Mat`, [Sparse Matrices](sec_matsparse), `MatCreate()`, `MatCreateSeqAIJ()`, `MatSetValues()`, `MatCreateBAIJ()`

# External Links
$(_doc_external("Mat/MatCreateSeqBAIJ"))
"""
function MatCreateSeqBAIJ(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Union{Ptr, Vector{PetscInt}}) end

@for_petsc function MatCreateSeqBAIJ(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, m::$PetscInt, n::$PetscInt, nz::$PetscInt, nnz::Union{Ptr, Vector{$PetscInt}} )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateSeqBAIJ, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{CMat}),
               comm, bs, m, n, nz, nnz, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	A::PetscMat = MatCreateSeqBAIJMKL(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Union{Ptr, Vector{PetscInt}}) 
Creates a sparse matrix of type `MATSEQBAIJMKL`.
This type inherits from `MATSEQBAIJ` and is largely identical, but uses sparse BLAS
routines from Intel MKL whenever possible.

Input Parameters:
- `comm` - MPI communicator, set to `PETSC_COMM_SELF`
- `bs`   - size of block, the blocks are ALWAYS square. One can use `MatSetBlockSizes()` to set a different row and column blocksize but the row
blocksize always defines the size of the blocks. The column blocksize sets the blocksize of the vectors obtained with `MatCreateVecs()`
- `m`    - number of rows
- `n`    - number of columns
- `nz`   - number of nonzero blocks  per block row (same for all rows)
- `nnz`  - array containing the number of nonzero blocks in the various block rows
(possibly different for each block row) or `NULL`

Output Parameter:
- `A` - the matrix

It is recommended that one use the `MatCreate()`, `MatSetType()` and/or `MatSetFromOptions()`,
MatXXXXSetPreallocation() paradigm instead of this routine directly.
[MatXXXXSetPreallocation() is, for example, `MatSeqBAIJSetPreallocation()`]

Options Database Keys:
- `-mat_no_unroll`  - uses code that does not unroll the loops in the block calculations (much slower)
- `-mat_block_size` - size of the blocks to use

Level: intermediate

-seealso: [Sparse Matrices](sec_matsparse), `MatCreate()`, `MatCreateSeqAIJ()`, `MatSetValues()`, `MatCreateBAIJ()`

# External Links
$(_doc_external("Mat/MatCreateSeqBAIJMKL"))
"""
function MatCreateSeqBAIJMKL(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Union{Ptr, Vector{PetscInt}}) end

@for_petsc function MatCreateSeqBAIJMKL(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, m::$PetscInt, n::$PetscInt, nz::$PetscInt, nnz::Union{Ptr, Vector{$PetscInt}} )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateSeqBAIJMKL, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{CMat}),
               comm, bs, m, n, nz, nnz, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	mat::PetscMat = MatCreateSeqBAIJWithArrays(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, a::Vector{PetscScalar}) 
Creates a `MATSEQBAIJ` matrix using matrix elements provided by the user.

Collective

Input Parameters:
- `comm` - must be an MPI communicator of size 1
- `bs`   - size of block
- `m`    - number of rows
- `n`    - number of columns
- `i`    - row indices; that is i[0] = 0, i[row] = i[row-1] + number of elements in that row block row of the matrix
- `j`    - column indices
- `a`    - matrix values

Output Parameter:
- `mat` - the matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `MatCreateBAIJ()`, `MatCreateSeqBAIJ()`

# External Links
$(_doc_external("Mat/MatCreateSeqBAIJWithArrays"))
"""
function MatCreateSeqBAIJWithArrays(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, a::Vector{PetscScalar}) end

@for_petsc function MatCreateSeqBAIJWithArrays(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, m::$PetscInt, n::$PetscInt, i::Vector{$PetscInt}, j::Vector{$PetscInt}, a::Vector{$PetscScalar} )
	mat_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateSeqBAIJWithArrays, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}, Ptr{CMat}),
               comm, bs, m, n, i, j, a, mat_,
              )

	mat = PetscMat(mat_[], petsclib)

	return mat
end 

"""
	A::PetscMat = MatCreateSeqDense(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, data::Union{Ptr, Vector{PetscScalar}}) 
Creates a `MATSEQDENSE` that
is stored in column major order (the usual Fortran format).

Collective

Input Parameters:
- `comm` - MPI communicator, set to `PETSC_COMM_SELF`
- `m`    - number of rows
- `n`    - number of columns
- `data` - optional location of matrix data in column major order.  Use `NULL` for PETSc
to control all matrix memory allocation.

Output Parameter:
- `A` - the matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATSEQDENSE`, `MatCreate()`, `MatCreateDense()`, `MatSetValues()`

# External Links
$(_doc_external("Mat/MatCreateSeqDense"))
"""
function MatCreateSeqDense(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, data::Union{Ptr, Vector{PetscScalar}}) end

@for_petsc function MatCreateSeqDense(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, data::Union{Ptr, Vector{$PetscScalar}} )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateSeqDense, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{CMat}),
               comm, m, n, data, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	A::PetscMat = MatCreateSeqSBAIJ(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Union{Ptr, Vector{PetscInt}}) 
Creates a sparse symmetric matrix in (block
compressed row) `MATSEQSBAIJ` format.  For good matrix assembly performance the
user should preallocate the matrix storage by setting the parameter `nz`
(or the array `nnz`).

Collective

Input Parameters:
- `comm` - MPI communicator, set to `PETSC_COMM_SELF`
- `bs`   - size of block, the blocks are ALWAYS square. One can use `MatSetBlockSizes()` to set a different row and column blocksize but the row
blocksize always defines the size of the blocks. The column blocksize sets the blocksize of the vectors obtained with MatCreateVecs()
- `m`    - number of rows
- `n`    - number of columns
- `nz`   - number of block nonzeros per block row (same for all rows)
- `nnz`  - array containing the number of block nonzeros in the upper triangular plus
diagonal portion of each block (possibly different for each block row) or `NULL`

Output Parameter:
- `A` - the symmetric matrix

Options Database Keys:
- `-mat_no_unroll`  - uses code that does not unroll the loops in the block calculations (much slower)
- `-mat_block_size` - size of the blocks to use

Level: intermediate

-seealso: [](ch_matrices), `Mat`, [Sparse Matrices](sec_matsparse), `MATSEQSBAIJ`, `MatCreate()`, `MatCreateSeqAIJ()`, `MatSetValues()`, `MatCreateSBAIJ()`

# External Links
$(_doc_external("Mat/MatCreateSeqSBAIJ"))
"""
function MatCreateSeqSBAIJ(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Union{Ptr, Vector{PetscInt}}) end

@for_petsc function MatCreateSeqSBAIJ(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, m::$PetscInt, n::$PetscInt, nz::$PetscInt, nnz::Union{Ptr, Vector{$PetscInt}} )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateSeqSBAIJ, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{CMat}),
               comm, bs, m, n, nz, nnz, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	mat::PetscMat = MatCreateSeqSBAIJWithArrays(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, a::Vector{PetscScalar}) 
Creates an sequential `MATSEQSBAIJ` matrix using matrix elements
(upper triangular entries in CSR format) provided by the user.

Collective

Input Parameters:
- `comm` - must be an MPI communicator of size 1
- `bs`   - size of block
- `m`    - number of rows
- `n`    - number of columns
- `i`    - row indices; that is i[0] = 0, i[row] = i[row-1] + number of block elements in that row block row of the matrix
- `j`    - column indices
- `a`    - matrix values

Output Parameter:
- `mat` - the matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATSEQSBAIJ`, `MatCreate()`, `MatCreateSBAIJ()`, `MatCreateSeqSBAIJ()`

# External Links
$(_doc_external("Mat/MatCreateSeqSBAIJWithArrays"))
"""
function MatCreateSeqSBAIJWithArrays(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, a::Vector{PetscScalar}) end

@for_petsc function MatCreateSeqSBAIJWithArrays(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, m::$PetscInt, n::$PetscInt, i::Vector{$PetscInt}, j::Vector{$PetscInt}, a::Vector{$PetscScalar} )
	mat_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateSeqSBAIJWithArrays, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}, Ptr{CMat}),
               comm, bs, m, n, i, j, a, mat_,
              )

	mat = PetscMat(mat_[], petsclib)

	return mat
end 

"""
	A::PetscMat = MatCreateSeqSELL(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, rlenmax::PetscInt, rlen::Union{Ptr, Vector{PetscInt}}) 
Creates a sparse matrix in `MATSEQSELL` format.

Collective

Input Parameters:
- `comm`    - MPI communicator, set to `PETSC_COMM_SELF`
- `m`       - number of rows
- `n`       - number of columns
- `rlenmax` - maximum number of nonzeros in a row, ignored if `rlen` is provided
- `rlen`    - array containing the number of nonzeros in the various rows (possibly different for each row) or NULL

Output Parameter:
- `A` - the matrix

Level: intermediate

-seealso: `Mat`, `MATSEQSELL`, `MatCreate()`, `MatCreateSELL()`, `MatSetValues()`, `MatSeqSELLSetPreallocation()`, `MATSELL`, `MATMPISELL`

# External Links
$(_doc_external("Mat/MatCreateSeqSELL"))
"""
function MatCreateSeqSELL(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, rlenmax::PetscInt, rlen::Union{Ptr, Vector{PetscInt}}) end

@for_petsc function MatCreateSeqSELL(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, rlenmax::$PetscInt, rlen::Union{Ptr, Vector{$PetscInt}} )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateSeqSELL, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{CMat}),
               comm, m, n, rlenmax, rlen, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	A::PetscMat = MatCreateShell(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, ctx::Ptr{Cvoid}) 
Creates a new matrix of `MatType` `MATSHELL` for use with a user
private matrix data storage format.

Collective

Input Parameters:
- `comm` - MPI communicator
- `m`    - number of local rows (or `PETSC_DECIDE` to have calculated if `M` is given)
- `n`    - number of local columns (or `PETSC_DECIDE` to have calculated if `N` is given)
- `M`    - number of global rows (may be `PETSC_DETERMINE` to have calculated if `m` is given)
- `N`    - number of global columns (may be `PETSC_DETERMINE` to have calculated if `n` is given)
- `ctx`  - pointer to data needed by the shell matrix routines

Output Parameter:
- `A` - the matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATSHELL`, `MatShellSetOperation()`, `MatHasOperation()`, `MatShellGetContext()`, `MatShellSetContext()`, `MatShellSetManageScalingShifts()`, `MatShellSetMatProductOperation()`

# External Links
$(_doc_external("Mat/MatCreateShell"))
"""
function MatCreateShell(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, ctx::Ptr{Cvoid}) end

@for_petsc function MatCreateShell(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, ctx::Ptr{Cvoid} )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateShell, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{Cvoid}, Ptr{CMat}),
               comm, m, n, M, N, ctx, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	submat::Ptr{PetscMat} = MatCreateSubMatrices(petsclib::PetscLibType,mat::AbstractPetscMat, n::PetscInt, irow::Vector{<:AbstractIS}, icol::Vector{<:AbstractIS}, scall::MatReuse) 
Extracts several submatrices from a matrix. If submat
points to an array of valid matrices, they may be reused to store the new
submatrices.

Collective

Input Parameters:
- `mat`   - the matrix
- `n`     - the number of submatrixes to be extracted (on this processor, may be zero)
- `irow`  - index set of rows to extract
- `icol`  - index set of columns to extract
- `scall` - either `MAT_INITIAL_MATRIX` or `MAT_REUSE_MATRIX`

Output Parameter:
- `submat` - the array of submatrices

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatDestroySubMatrices()`, `MatCreateSubMatrix()`, `MatGetRow()`, `MatGetDiagonal()`, `MatReuse`

# External Links
$(_doc_external("Mat/MatCreateSubMatrices"))
"""
function MatCreateSubMatrices(petsclib::PetscLibType, mat::AbstractPetscMat, n::PetscInt, irow::Vector{<:AbstractIS}, icol::Vector{<:AbstractIS}, scall::MatReuse) end

@for_petsc function MatCreateSubMatrices(petsclib::$UnionPetscLib, mat::AbstractPetscMat, n::$PetscInt, irow::Vector{<:AbstractIS}, icol::Vector{<:AbstractIS}, scall::MatReuse )
	submat_ = Ref{Ptr{PetscMat}}()

    @chk ccall(
               (:MatCreateSubMatrices, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CIS}, Ptr{CIS}, MatReuse, Ptr{Ptr{CMat}}),
               mat, n, irow, icol, scall, submat_,
              )

	submat = submat_[]

	return submat
end 

"""
	submat::Ptr{PetscMat} = MatCreateSubMatricesMPI(petsclib::PetscLibType,mat::AbstractPetscMat, n::PetscInt, irow::Vector{<:AbstractIS}, icol::Vector{<:AbstractIS}, scall::MatReuse) 
Extracts MPI submatrices across a sub communicator of `mat` (by pairs of `IS` that may live on subcomms).

Collective

Input Parameters:
- `mat`   - the matrix
- `n`     - the number of submatrixes to be extracted
- `irow`  - index set of rows to extract
- `icol`  - index set of columns to extract
- `scall` - either `MAT_INITIAL_MATRIX` or `MAT_REUSE_MATRIX`

Output Parameter:
- `submat` - the array of submatrices

Level: advanced

-seealso: [](ch_matrices), `Mat`, `PCGASM`, `MatCreateSubMatrices()`, `MatCreateSubMatrix()`, `MatGetRow()`, `MatGetDiagonal()`, `MatReuse`

# External Links
$(_doc_external("Mat/MatCreateSubMatricesMPI"))
"""
function MatCreateSubMatricesMPI(petsclib::PetscLibType, mat::AbstractPetscMat, n::PetscInt, irow::Vector{<:AbstractIS}, icol::Vector{<:AbstractIS}, scall::MatReuse) end

@for_petsc function MatCreateSubMatricesMPI(petsclib::$UnionPetscLib, mat::AbstractPetscMat, n::$PetscInt, irow::Vector{<:AbstractIS}, icol::Vector{<:AbstractIS}, scall::MatReuse )
	submat_ = Ref{Ptr{PetscMat}}()

    @chk ccall(
               (:MatCreateSubMatricesMPI, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CIS}, Ptr{CIS}, MatReuse, Ptr{Ptr{CMat}}),
               mat, n, irow, icol, scall, submat_,
              )

	submat = submat_[]

	return submat
end 

"""
	newmat::PetscMat = MatCreateSubMatrix(petsclib::PetscLibType,mat::AbstractPetscMat, isrow::AbstractIS, iscol::AbstractIS, cll::MatReuse) 
Gets a single submatrix on the same number of processors
as the original matrix.

Collective

Input Parameters:
- `mat`   - the original matrix
- `isrow` - parallel `IS` containing the rows this processor should obtain
- `iscol` - parallel `IS` containing all columns you wish to keep. Each process should list the columns that will be in IT's "diagonal part" in the new matrix.
- `cll`   - either `MAT_INITIAL_MATRIX` or `MAT_REUSE_MATRIX`

Output Parameter:
- `newmat` - the new submatrix, of the same type as the original matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatCreateSubMatrices()`, `MatCreateSubMatricesMPI()`, `MatCreateSubMatrixVirtual()`, `MatSubMatrixVirtualUpdate()`

# External Links
$(_doc_external("Mat/MatCreateSubMatrix"))
"""
function MatCreateSubMatrix(petsclib::PetscLibType, mat::AbstractPetscMat, isrow::AbstractIS, iscol::AbstractIS, cll::MatReuse) end

@for_petsc function MatCreateSubMatrix(petsclib::$UnionPetscLib, mat::AbstractPetscMat, isrow::AbstractIS, iscol::AbstractIS, cll::MatReuse )
	newmat_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateSubMatrix, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, CIS, MatReuse, Ptr{CMat}),
               mat, isrow, iscol, cll, newmat_,
              )

	newmat = PetscMat(newmat_[], petsclib)

	return newmat
end 

"""
	J::PetscMat = MatCreateSubMatrixFree(petsclib::PetscLibType,mat::AbstractPetscMat, Rows::AbstractIS, Cols::AbstractIS) 
Creates a reduced matrix by masking a
full matrix.

Collective

Input Parameters:
- `mat`  - matrix of arbitrary type
- `Rows` - the rows that will be in the submatrix
- `Cols` - the columns that will be in the submatrix

Output Parameter:
- `J` - New matrix

Level: developer

-seealso: `MatCreate()`

# External Links
$(_doc_external("Tao/MatCreateSubMatrixFree"))
"""
function MatCreateSubMatrixFree(petsclib::PetscLibType, mat::AbstractPetscMat, Rows::AbstractIS, Cols::AbstractIS) end

@for_petsc function MatCreateSubMatrixFree(petsclib::$UnionPetscLib, mat::AbstractPetscMat, Rows::AbstractIS, Cols::AbstractIS )
	J_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateSubMatrixFree, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, CIS, Ptr{CMat}),
               mat, Rows, Cols, J_,
              )

	J = PetscMat(J_[], petsclib)

	return J
end 

"""
	newmat::PetscMat = MatCreateSubMatrixVirtual(petsclib::PetscLibType,A::AbstractPetscMat, isrow::AbstractIS, iscol::AbstractIS) 
Creates a virtual matrix `MATSUBMATRIX` that acts as a submatrix

Collective

Input Parameters:
- `A`     - matrix that we will extract a submatrix of
- `isrow` - rows to be present in the submatrix
- `iscol` - columns to be present in the submatrix

Output Parameter:
- `newmat` - new matrix

Level: developer

-seealso: [](ch_matrices), `Mat`, `MATSUBMATRIX`, `MATLOCALREF`, `MatCreateLocalRef()`, `MatCreateSubMatrix()`, `MatSubMatrixVirtualUpdate()`

# External Links
$(_doc_external("Mat/MatCreateSubMatrixVirtual"))
"""
function MatCreateSubMatrixVirtual(petsclib::PetscLibType, A::AbstractPetscMat, isrow::AbstractIS, iscol::AbstractIS) end

@for_petsc function MatCreateSubMatrixVirtual(petsclib::$UnionPetscLib, A::AbstractPetscMat, isrow::AbstractIS, iscol::AbstractIS )
	newmat_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateSubMatrixVirtual, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, CIS, Ptr{CMat}),
               A, isrow, iscol, newmat_,
              )

	newmat = PetscMat(newmat_[], petsclib)

	return newmat
end 

"""
	N::PetscMat = MatCreateTranspose(petsclib::PetscLibType,A::AbstractPetscMat) 
Creates a new matrix `MATTRANSPOSEVIRTUAL` object that behaves like A'

Collective

Input Parameter:
- `A` - the (possibly rectangular) matrix

Output Parameter:
- `N` - the matrix that represents A'

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATTRANSPOSEVIRTUAL`, `MatCreateNormal()`, `MatMult()`, `MatMultTranspose()`, `MatCreate()`,
`MATNORMALHERMITIAN`

# External Links
$(_doc_external("Mat/MatCreateTranspose"))
"""
function MatCreateTranspose(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatCreateTranspose(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	N_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateTranspose, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, N_,
              )

	N = PetscMat(N_[], petsclib)

	return N
end 

"""
	right::PetscVec,left::PetscVec = MatCreateVecs(petsclib::PetscLibType,mat::AbstractPetscMat) 
Get vector(s) compatible with the matrix, i.e. with the same
parallel layout, `PetscLayout` for rows and columns

Collective

Input Parameter:
- `mat` - the matrix

Output Parameters:
- `right` - (optional) vector that the matrix can be multiplied against
- `left`  - (optional) vector that the matrix vector product can be stored in

Level: advanced

-seealso: [](ch_matrices), `Mat`, `Vec`, `VecCreate()`, `VecDestroy()`, `DMCreateGlobalVector()`

# External Links
$(_doc_external("Mat/MatCreateVecs"))
"""
function MatCreateVecs(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatCreateVecs(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	right_ = Ref{CVec}()
	left_ = Ref{CVec}()

    @chk ccall(
               (:MatCreateVecs, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CVec}, Ptr{CVec}),
               mat, right_, left_,
              )

	right = PetscVec(right_[], petsclib)
	left = PetscVec(left_[], petsclib)

	return right,left
end 

"""
	x::PetscVec,y::PetscVec,z::PetscVec = MatCreateVecsFFTW(petsclib::PetscLibType,A::AbstractPetscMat) 
Get vector(s) compatible with the matrix, i.e. with the
parallel layout determined by `MATFFTW`

Collective

Input Parameter:
- `A` - the matrix

Output Parameters:
- `x` - (optional) input vector of forward FFTW
- `y` - (optional) output vector of forward FFTW
- `z` - (optional) output vector of backward FFTW

Options Database Key:
- `-mat_fftw_plannerflags` - set FFTW planner flags

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATFFTW`, `MatCreateFFT()`, `MatCreateVecs()`

# External Links
$(_doc_external("Mat/MatCreateVecsFFTW"))
"""
function MatCreateVecsFFTW(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatCreateVecsFFTW(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	x_ = Ref{CVec}()
	y_ = Ref{CVec}()
	z_ = Ref{CVec}()

    @chk ccall(
               (:MatCreateVecsFFTW, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CVec}, Ptr{CVec}, Ptr{CVec}),
               A, x_, y_, z_,
              )

	x = PetscVec(x_[], petsclib)
	y = PetscVec(y_[], petsclib)
	z = PetscVec(z_[], petsclib)

	return x,y,z
end 

"""
	MatDFischer(petsclib::PetscLibType,jac::AbstractPetscMat, X::AbstractPetscVec, Con::AbstractPetscVec, XL::AbstractPetscVec, XU::AbstractPetscVec, T1::AbstractPetscVec, T2::AbstractPetscVec, Da::AbstractPetscVec, Db::AbstractPetscVec) 
Calculates an element of the B
Fischer-Burmeister function for complementarity problems.

Collective

Input Parameters:
- `jac` - the jacobian of `f` at `X`
- `X`   - current point
- `Con` - constraints function evaluated at `X`
- `XL`  - lower bounds
- `XU`  - upper bounds
- `T1`  - work vector
- `T2`  - work vector

Output Parameters:
- `Da` - diagonal perturbation component of the result
- `Db` - row scaling component of the result

Level: developer

-seealso: `Mat`, `VecFischer()`, `VecSFischer()`, `MatDSFischer()`

# External Links
$(_doc_external("Tao/MatDFischer"))
"""
function MatDFischer(petsclib::PetscLibType, jac::AbstractPetscMat, X::AbstractPetscVec, Con::AbstractPetscVec, XL::AbstractPetscVec, XU::AbstractPetscVec, T1::AbstractPetscVec, T2::AbstractPetscVec, Da::AbstractPetscVec, Db::AbstractPetscVec) end

@for_petsc function MatDFischer(petsclib::$UnionPetscLib, jac::AbstractPetscMat, X::AbstractPetscVec, Con::AbstractPetscVec, XL::AbstractPetscVec, XU::AbstractPetscVec, T1::AbstractPetscVec, T2::AbstractPetscVec, Da::AbstractPetscVec, Db::AbstractPetscVec )

    @chk ccall(
               (:MatDFischer, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec, CVec, CVec, CVec, CVec, CVec, CVec),
               jac, X, Con, XL, XU, T1, T2, Da, Db,
              )


	return nothing
end 

"""
	MatDSFischer(petsclib::PetscLibType,jac::AbstractPetscMat, X::AbstractPetscVec, Con::AbstractPetscVec, XL::AbstractPetscVec, XU::AbstractPetscVec, mu::PetscReal, T1::AbstractPetscVec, T2::AbstractPetscVec, Da::AbstractPetscVec, Db::AbstractPetscVec, Dm::AbstractPetscVec) 
Calculates an element of the B
smoothed Fischer-Burmeister function for complementarity problems.

Collective

Input Parameters:
- `jac` - the jacobian of f at X
- `X`   - current point
- `Con` - constraint function evaluated at X
- `XL`  - lower bounds
- `XU`  - upper bounds
- `mu`  - smoothing parameter
- `T1`  - work vector
- `T2`  - work vector

Output Parameters:
- `Da` - diagonal perturbation component of the result
- `Db` - row scaling component of the result
- `Dm` - derivative with respect to scaling parameter

Level: developer

-seealso: `Mat`, `VecFischer()`, `VecDFischer()`, `MatDFischer()`

# External Links
$(_doc_external("Tao/MatDSFischer"))
"""
function MatDSFischer(petsclib::PetscLibType, jac::AbstractPetscMat, X::AbstractPetscVec, Con::AbstractPetscVec, XL::AbstractPetscVec, XU::AbstractPetscVec, mu::PetscReal, T1::AbstractPetscVec, T2::AbstractPetscVec, Da::AbstractPetscVec, Db::AbstractPetscVec, Dm::AbstractPetscVec) end

@for_petsc function MatDSFischer(petsclib::$UnionPetscLib, jac::AbstractPetscMat, X::AbstractPetscVec, Con::AbstractPetscVec, XL::AbstractPetscVec, XU::AbstractPetscVec, mu::$PetscReal, T1::AbstractPetscVec, T2::AbstractPetscVec, Da::AbstractPetscVec, Db::AbstractPetscVec, Dm::AbstractPetscVec )

    @chk ccall(
               (:MatDSFischer, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec, CVec, CVec, $PetscReal, CVec, CVec, CVec, CVec, CVec),
               jac, X, Con, XL, XU, mu, T1, T2, Da, Db, Dm,
              )


	return nothing
end 

"""
	array::Vector{PetscScalar} = MatDenseGetArray(petsclib::PetscLibType,A::AbstractPetscMat) 
gives read

Logically Collective

Input Parameter:
- `A` - a dense matrix

Output Parameter:
- `array` - pointer to the data

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MatDenseRestoreArray()`, `MatDenseGetArrayRead()`, `MatDenseRestoreArrayRead()`, `MatDenseGetArrayWrite()`, `MatDenseRestoreArrayWrite()`

# External Links
$(_doc_external("Mat/MatDenseGetArray"))
"""
function MatDenseGetArray(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatDenseGetArray(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatDenseGetArray, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

	m, n = MatGetLocalSize(petsclib, A)
	array = unsafe_wrap(Array, array_[], (Int(m), Int(n)); own = false)

	return array
end 

"""
	array::Vector{PetscScalar},mtype::PetscMemType = MatDenseGetArrayAndMemType(petsclib::PetscLibType,A::AbstractPetscMat) 
gives read

Logically Collective

Input Parameter:
- `A` - a dense matrix

Output Parameters:
- `array` - pointer to the data
- `mtype` - memory type of the returned pointer

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MatDenseRestoreArrayAndMemType()`, `MatDenseGetArrayReadAndMemType()`, `MatDenseGetArrayWriteAndMemType()`, `MatDenseGetArrayRead()`,
`MatDenseRestoreArrayRead()`, `MatDenseGetArrayWrite()`, `MatDenseRestoreArrayWrite()`, `MatSeqAIJGetCSRAndMemType()`

# External Links
$(_doc_external("Mat/MatDenseGetArrayAndMemType"))
"""
function MatDenseGetArrayAndMemType(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatDenseGetArrayAndMemType(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()
	mtype_ = Ref{PetscMemType}()

    @chk ccall(
               (:MatDenseGetArrayAndMemType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}, Ptr{PetscMemType}),
               A, array_, mtype_,
              )

	m, n = MatGetLocalSize(petsclib, A)
	array = unsafe_wrap(Array, array_[], (Int(m), Int(n)); own = false)
	mtype = mtype_[]

	return array,mtype
end 

"""
	array::Vector{PetscScalar} = MatDenseGetArrayRead(petsclib::PetscLibType,A::AbstractPetscMat) 
gives read

Not Collective

Input Parameter:
- `A` - a dense matrix

Output Parameter:
- `array` - pointer to the data

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MatDenseRestoreArrayRead()`, `MatDenseGetArray()`, `MatDenseRestoreArray()`, `MatDenseGetArrayWrite()`, `MatDenseRestoreArrayWrite()`

# External Links
$(_doc_external("Mat/MatDenseGetArrayRead"))
"""
function MatDenseGetArrayRead(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatDenseGetArrayRead(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatDenseGetArrayRead, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

	m, n = MatGetLocalSize(petsclib, A)
	array = unsafe_wrap(Array, array_[], (Int(m), Int(n)); own = false)

	return array
end 

"""
	array::Vector{PetscScalar},mtype::PetscMemType = MatDenseGetArrayReadAndMemType(petsclib::PetscLibType,A::AbstractPetscMat) 
gives read

Logically Collective

Input Parameter:
- `A` - a dense matrix

Output Parameters:
- `array` - pointer to the data
- `mtype` - memory type of the returned pointer

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MatDenseRestoreArrayReadAndMemType()`, `MatDenseGetArrayWriteAndMemType()`,
`MatDenseGetArrayRead()`, `MatDenseRestoreArrayRead()`, `MatDenseGetArrayWrite()`, `MatDenseRestoreArrayWrite()`, `MatSeqAIJGetCSRAndMemType()`

# External Links
$(_doc_external("Mat/MatDenseGetArrayReadAndMemType"))
"""
function MatDenseGetArrayReadAndMemType(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatDenseGetArrayReadAndMemType(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()
	mtype_ = Ref{PetscMemType}()

    @chk ccall(
               (:MatDenseGetArrayReadAndMemType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}, Ptr{PetscMemType}),
               A, array_, mtype_,
              )

	m, n = MatGetLocalSize(petsclib, A)
	array = unsafe_wrap(Array, array_[], (Int(m), Int(n)); own = false)
	mtype = mtype_[]

	return array,mtype
end 

"""
	array::Vector{PetscScalar} = MatDenseGetArrayWrite(petsclib::PetscLibType,A::AbstractPetscMat) 
gives write

Not Collective

Input Parameter:
- `A` - a dense matrix

Output Parameter:
- `array` - pointer to the data

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MatDenseRestoreArrayWrite()`, `MatDenseGetArray()`, `MatDenseRestoreArray()`, `MatDenseGetArrayRead()`, `MatDenseRestoreArrayRead()`

# External Links
$(_doc_external("Mat/MatDenseGetArrayWrite"))
"""
function MatDenseGetArrayWrite(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatDenseGetArrayWrite(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatDenseGetArrayWrite, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

	m, n = MatGetLocalSize(petsclib, A)
	array = unsafe_wrap(Array, array_[], (Int(m), Int(n)); own = false)

	return array
end 

"""
	array::Vector{PetscScalar},mtype::PetscMemType = MatDenseGetArrayWriteAndMemType(petsclib::PetscLibType,A::AbstractPetscMat) 
gives write

Logically Collective

Input Parameter:
- `A` - a dense matrix

Output Parameters:
- `array` - pointer to the data
- `mtype` - memory type of the returned pointer

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MatDenseRestoreArrayWriteAndMemType()`, `MatDenseGetArrayReadAndMemType()`, `MatDenseGetArrayRead()`,
`MatDenseRestoreArrayRead()`, `MatDenseGetArrayWrite()`, `MatDenseRestoreArrayWrite()`, `MatSeqAIJGetCSRAndMemType()`

# External Links
$(_doc_external("Mat/MatDenseGetArrayWriteAndMemType"))
"""
function MatDenseGetArrayWriteAndMemType(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatDenseGetArrayWriteAndMemType(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()
	mtype_ = Ref{PetscMemType}()

    @chk ccall(
               (:MatDenseGetArrayWriteAndMemType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}, Ptr{PetscMemType}),
               A, array_, mtype_,
              )

	m, n = MatGetLocalSize(petsclib, A)
	array = unsafe_wrap(Array, array_[], (Int(m), Int(n)); own = false)
	mtype = mtype_[]

	return array,mtype
end 

"""
	vals::Ptr{PetscScalar} = MatDenseGetColumn(petsclib::PetscLibType,A::AbstractPetscMat, col::PetscInt) 
gives access to a column of a dense matrix. This is only the local part of the column. You MUST call `MatDenseRestoreColumn()` to avoid memory bleeding.

Not Collective

Input Parameters:
- `A`   - a `MATSEQDENSE` or `MATMPIDENSE` matrix
- `col` - column index

Output Parameter:
- `vals` - pointer to the data

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MatDenseRestoreColumn()`, `MatDenseGetColumnVec()`

# External Links
$(_doc_external("Mat/MatDenseGetColumn"))
"""
function MatDenseGetColumn(petsclib::PetscLibType, A::AbstractPetscMat, col::PetscInt) end

@for_petsc function MatDenseGetColumn(petsclib::$UnionPetscLib, A::AbstractPetscMat, col::$PetscInt )
	vals_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatDenseGetColumn, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{Ptr{$PetscScalar}}),
               A, col, vals_,
              )

	vals = vals_[]

	return vals
end 

"""
	v::PetscVec = MatDenseGetColumnVec(petsclib::PetscLibType,A::AbstractPetscMat, col::PetscInt) 
Gives read

Collective

Input Parameters:
- `A`   - the `Mat` object
- `col` - the column index

Output Parameter:
- `v` - the vector

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MATDENSECUDA`, `MATDENSEHIP`, `MatDenseGetColumnVecRead()`, `MatDenseGetColumnVecWrite()`, `MatDenseRestoreColumnVec()`, `MatDenseRestoreColumnVecRead()`, `MatDenseRestoreColumnVecWrite()`, `MatDenseGetColumn()`

# External Links
$(_doc_external("Mat/MatDenseGetColumnVec"))
"""
function MatDenseGetColumnVec(petsclib::PetscLibType, A::AbstractPetscMat, col::PetscInt) end

@for_petsc function MatDenseGetColumnVec(petsclib::$UnionPetscLib, A::AbstractPetscMat, col::$PetscInt )
	v_ = Ref{CVec}()

    @chk ccall(
               (:MatDenseGetColumnVec, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CVec}),
               A, col, v_,
              )

	v = PetscVec(v_[], petsclib)

	return v
end 

"""
	v::PetscVec = MatDenseGetColumnVecRead(petsclib::PetscLibType,A::AbstractPetscMat, col::PetscInt) 
Gives read

Collective

Input Parameters:
- `A`   - the `Mat` object
- `col` - the column index

Output Parameter:
- `v` - the vector

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MATDENSECUDA`, `MATDENSEHIP`, `MatDenseGetColumnVec()`, `MatDenseGetColumnVecWrite()`, `MatDenseRestoreColumnVec()`, `MatDenseRestoreColumnVecRead()`, `MatDenseRestoreColumnVecWrite()`

# External Links
$(_doc_external("Mat/MatDenseGetColumnVecRead"))
"""
function MatDenseGetColumnVecRead(petsclib::PetscLibType, A::AbstractPetscMat, col::PetscInt) end

@for_petsc function MatDenseGetColumnVecRead(petsclib::$UnionPetscLib, A::AbstractPetscMat, col::$PetscInt )
	v_ = Ref{CVec}()

    @chk ccall(
               (:MatDenseGetColumnVecRead, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CVec}),
               A, col, v_,
              )

	v = PetscVec(v_[], petsclib)

	return v
end 

"""
	v::PetscVec = MatDenseGetColumnVecWrite(petsclib::PetscLibType,A::AbstractPetscMat, col::PetscInt) 
Gives write

Collective

Input Parameters:
- `A`   - the `Mat` object
- `col` - the column index

Output Parameter:
- `v` - the vector

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MATDENSECUDA`, `MATDENSEHIP`, `MatDenseGetColumnVec()`, `MatDenseGetColumnVecRead()`, `MatDenseRestoreColumnVec()`, `MatDenseRestoreColumnVecRead()`, `MatDenseRestoreColumnVecWrite()`

# External Links
$(_doc_external("Mat/MatDenseGetColumnVecWrite"))
"""
function MatDenseGetColumnVecWrite(petsclib::PetscLibType, A::AbstractPetscMat, col::PetscInt) end

@for_petsc function MatDenseGetColumnVecWrite(petsclib::$UnionPetscLib, A::AbstractPetscMat, col::$PetscInt )
	v_ = Ref{CVec}()

    @chk ccall(
               (:MatDenseGetColumnVecWrite, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CVec}),
               A, col, v_,
              )

	v = PetscVec(v_[], petsclib)

	return v
end 

"""
	lda::PetscInt = MatDenseGetLDA(petsclib::PetscLibType,A::AbstractPetscMat) 
gets the leading dimension of the array returned from `MatDenseGetArray()`

Not Collective

Input Parameter:
- `A` - a `MATDENSE` or `MATDENSECUDA` matrix

Output Parameter:
- `lda` - the leading dimension

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MATDENSECUDA`, `MatDenseGetArray()`, `MatDenseRestoreArray()`, `MatDenseGetArrayRead()`, `MatDenseRestoreArrayRead()`, `MatDenseSetLDA()`

# External Links
$(_doc_external("Mat/MatDenseGetLDA"))
"""
function MatDenseGetLDA(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatDenseGetLDA(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	lda_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatDenseGetLDA, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}),
               A, lda_,
              )

	lda = lda_[]

	return lda
end 

"""
	B::PetscMat = MatDenseGetLocalMatrix(petsclib::PetscLibType,A::AbstractPetscMat) 
For a `MATMPIDENSE` or `MATSEQDENSE` matrix returns the sequential
matrix that represents the operator. For sequential matrices it returns itself.

Input Parameter:
- `A` - the sequential or MPI `MATDENSE` matrix

Output Parameter:
- `B` - the inner matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MATMPIDENSE`, `MATSEQDENSE`

# External Links
$(_doc_external("Mat/MatDenseGetLocalMatrix"))
"""
function MatDenseGetLocalMatrix(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatDenseGetLocalMatrix(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	B_ = Ref{CMat}()

    @chk ccall(
               (:MatDenseGetLocalMatrix, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, B_,
              )

	B = PetscMat(B_[], petsclib)

	return B
end 

"""
	v::PetscMat = MatDenseGetSubMatrix(petsclib::PetscLibType,A::AbstractPetscMat, rbegin::PetscInt, rend::PetscInt, cbegin::PetscInt, cend::PetscInt) 
Gives access to a block of rows and columns of a dense matrix, represented as a `Mat`.

Collective

Input Parameters:
- `A`      - the `Mat` object
- `rbegin` - the first global row index in the block (if `PETSC_DECIDE`, is 0)
- `rend`   - the global row index past the last one in the block (if `PETSC_DECIDE`, is `M`)
- `cbegin` - the first global column index in the block (if `PETSC_DECIDE`, is 0)
- `cend`   - the global column index past the last one in the block (if `PETSC_DECIDE`, is `N`)

Output Parameter:
- `v` - the matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MATDENSECUDA`, `MATDENSEHIP`, `MatDenseGetColumnVec()`, `MatDenseRestoreColumnVec()`, `MatDenseRestoreSubMatrix()`

# External Links
$(_doc_external("Mat/MatDenseGetSubMatrix"))
"""
function MatDenseGetSubMatrix(petsclib::PetscLibType, A::AbstractPetscMat, rbegin::PetscInt, rend::PetscInt, cbegin::PetscInt, cend::PetscInt) end

@for_petsc function MatDenseGetSubMatrix(petsclib::$UnionPetscLib, A::AbstractPetscMat, rbegin::$PetscInt, rend::$PetscInt, cbegin::$PetscInt, cend::$PetscInt )
	v_ = Ref{CMat}()

    @chk ccall(
               (:MatDenseGetSubMatrix, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{CMat}),
               A, rbegin, rend, cbegin, cend, v_,
              )

	v = PetscMat(v_[], petsclib)

	return v
end 

"""
	MatDensePlaceArray(petsclib::PetscLibType,mat::AbstractPetscMat, array::Vector{PetscScalar}) 
Allows one to replace the array in a `MATDENSE` matrix with an
array provided by the user. This is useful to avoid copying an array
into a matrix

Not Collective

Input Parameters:
- `mat`   - the matrix
- `array` - the array in column major order

Level: developer

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MatDenseGetArray()`, `MatDenseResetArray()`, `VecPlaceArray()`, `VecGetArray()`, `VecRestoreArray()`, `VecReplaceArray()`, `VecResetArray()`,
`MatDenseReplaceArray()`

# External Links
$(_doc_external("Mat/MatDensePlaceArray"))
"""
function MatDensePlaceArray(petsclib::PetscLibType, mat::AbstractPetscMat, array::Vector{PetscScalar}) end

@for_petsc function MatDensePlaceArray(petsclib::$UnionPetscLib, mat::AbstractPetscMat, array::Vector{$PetscScalar} )

    @chk ccall(
               (:MatDensePlaceArray, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscScalar}),
               mat, array,
              )


	return nothing
end 

"""
	MatDenseReplaceArray(petsclib::PetscLibType,mat::AbstractPetscMat, array::Vector{PetscScalar}) 
Allows one to replace the array in a dense matrix with an
array provided by the user. This is useful to avoid copying an array
into a matrix

Not Collective

Input Parameters:
- `mat`   - the matrix
- `array` - the array in column major order

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatDensePlaceArray()`, `MatDenseGetArray()`, `VecReplaceArray()`

# External Links
$(_doc_external("Mat/MatDenseReplaceArray"))
"""
function MatDenseReplaceArray(petsclib::PetscLibType, mat::AbstractPetscMat, array::Vector{PetscScalar}) end

@for_petsc function MatDenseReplaceArray(petsclib::$UnionPetscLib, mat::AbstractPetscMat, array::Vector{$PetscScalar} )

    @chk ccall(
               (:MatDenseReplaceArray, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscScalar}),
               mat, array,
              )


	return nothing
end 

"""
	MatDenseResetArray(petsclib::PetscLibType,mat::AbstractPetscMat) 
Resets the matrix array to that it previously had before the call to `MatDensePlaceArray()`

Not Collective

Input Parameter:
- `mat` - the matrix

Level: developer

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MatDenseGetArray()`, `MatDensePlaceArray()`, `VecPlaceArray()`, `VecGetArray()`, `VecRestoreArray()`, `VecReplaceArray()`, `VecResetArray()`

# External Links
$(_doc_external("Mat/MatDenseResetArray"))
"""
function MatDenseResetArray(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatDenseResetArray(petsclib::$UnionPetscLib, mat::AbstractPetscMat )

    @chk ccall(
               (:MatDenseResetArray, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	MatDenseRestoreArray(petsclib::PetscLibType,A::AbstractPetscMat, array::AbstractArray{PetscScalar}) 
returns access to the array where the data for a `MATDENSE` matrix is stored obtained by `MatDenseGetArray()`

Logically Collective

Input Parameters:
- `A`     - a dense matrix
- `array` - pointer to the data (may be `NULL`)

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MatDenseGetArray()`, `MatDenseGetArrayRead()`, `MatDenseRestoreArrayRead()`, `MatDenseGetArrayWrite()`, `MatDenseRestoreArrayWrite()`

# External Links
$(_doc_external("Mat/MatDenseRestoreArray"))
"""
function MatDenseRestoreArray(petsclib::PetscLibType, A::AbstractPetscMat, array::AbstractArray{PetscScalar}) end

@for_petsc function MatDenseRestoreArray(petsclib::$UnionPetscLib, A::AbstractPetscMat, array::AbstractArray{$PetscScalar} )
	array_ = Ref(pointer(array))

    @chk ccall(
               (:MatDenseRestoreArray, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )


	return nothing
end 

"""
	MatDenseRestoreArrayAndMemType(petsclib::PetscLibType,A::AbstractPetscMat, array::AbstractArray{PetscScalar}) 
returns access to the array that is obtained by `MatDenseGetArrayAndMemType()`

Logically Collective

Input Parameters:
- `A`     - a dense matrix
- `array` - pointer to the data

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MatDenseGetArrayAndMemType()`, `MatDenseGetArray()`, `MatDenseGetArrayRead()`, `MatDenseRestoreArrayRead()`, `MatDenseGetArrayWrite()`, `MatDenseRestoreArrayWrite()`

# External Links
$(_doc_external("Mat/MatDenseRestoreArrayAndMemType"))
"""
function MatDenseRestoreArrayAndMemType(petsclib::PetscLibType, A::AbstractPetscMat, array::AbstractArray{PetscScalar}) end

@for_petsc function MatDenseRestoreArrayAndMemType(petsclib::$UnionPetscLib, A::AbstractPetscMat, array::AbstractArray{$PetscScalar} )
	array_ = Ref(pointer(array))

    @chk ccall(
               (:MatDenseRestoreArrayAndMemType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )


	return nothing
end 

"""
	MatDenseRestoreArrayRead(petsclib::PetscLibType,A::AbstractPetscMat, array::AbstractArray{PetscScalar}) 
returns access to the array where the data for a `MATDENSE` matrix is stored obtained by `MatDenseGetArrayRead()`

Not Collective

Input Parameters:
- `A`     - a dense matrix
- `array` - pointer to the data (may be `NULL`)

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MatDenseGetArrayRead()`, `MatDenseGetArray()`, `MatDenseRestoreArray()`, `MatDenseGetArrayWrite()`, `MatDenseRestoreArrayWrite()`

# External Links
$(_doc_external("Mat/MatDenseRestoreArrayRead"))
"""
function MatDenseRestoreArrayRead(petsclib::PetscLibType, A::AbstractPetscMat, array::AbstractArray{PetscScalar}) end

@for_petsc function MatDenseRestoreArrayRead(petsclib::$UnionPetscLib, A::AbstractPetscMat, array::AbstractArray{$PetscScalar} )
	array_ = Ref(pointer(array))

    @chk ccall(
               (:MatDenseRestoreArrayRead, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )


	return nothing
end 

"""
	MatDenseRestoreArrayReadAndMemType(petsclib::PetscLibType,A::AbstractPetscMat, array::AbstractArray{PetscScalar}) 
returns access to the array that is obtained by `MatDenseGetArrayReadAndMemType()`

Logically Collective

Input Parameters:
- `A`     - a dense matrix
- `array` - pointer to the data

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MatDenseGetArrayReadAndMemType()`, `MatDenseGetArray()`, `MatDenseGetArrayRead()`, `MatDenseRestoreArrayRead()`, `MatDenseGetArrayWrite()`, `MatDenseRestoreArrayWrite()`

# External Links
$(_doc_external("Mat/MatDenseRestoreArrayReadAndMemType"))
"""
function MatDenseRestoreArrayReadAndMemType(petsclib::PetscLibType, A::AbstractPetscMat, array::AbstractArray{PetscScalar}) end

@for_petsc function MatDenseRestoreArrayReadAndMemType(petsclib::$UnionPetscLib, A::AbstractPetscMat, array::AbstractArray{$PetscScalar} )
	array_ = Ref(pointer(array))

    @chk ccall(
               (:MatDenseRestoreArrayReadAndMemType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )


	return nothing
end 

"""
	MatDenseRestoreArrayWrite(petsclib::PetscLibType,A::AbstractPetscMat, array::AbstractArray{PetscScalar}) 
returns access to the array where the data for a `MATDENSE` matrix is stored obtained by `MatDenseGetArrayWrite()`

Not Collective

Input Parameters:
- `A`     - a dense matrix
- `array` - pointer to the data (may be `NULL`)

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MatDenseGetArrayWrite()`, `MatDenseGetArray()`, `MatDenseRestoreArray()`, `MatDenseGetArrayRead()`, `MatDenseRestoreArrayRead()`

# External Links
$(_doc_external("Mat/MatDenseRestoreArrayWrite"))
"""
function MatDenseRestoreArrayWrite(petsclib::PetscLibType, A::AbstractPetscMat, array::AbstractArray{PetscScalar}) end

@for_petsc function MatDenseRestoreArrayWrite(petsclib::$UnionPetscLib, A::AbstractPetscMat, array::AbstractArray{$PetscScalar} )
	array_ = Ref(pointer(array))

    @chk ccall(
               (:MatDenseRestoreArrayWrite, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )


	return nothing
end 

"""
	MatDenseRestoreArrayWriteAndMemType(petsclib::PetscLibType,A::AbstractPetscMat, array::AbstractArray{PetscScalar}) 
returns access to the array that is obtained by `MatDenseGetArrayReadAndMemType()`

Logically Collective

Input Parameters:
- `A`     - a dense matrix
- `array` - pointer to the data

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MatDenseGetArrayWriteAndMemType()`, `MatDenseGetArray()`, `MatDenseGetArrayRead()`, `MatDenseRestoreArrayRead()`, `MatDenseGetArrayWrite()`, `MatDenseRestoreArrayWrite()`

# External Links
$(_doc_external("Mat/MatDenseRestoreArrayWriteAndMemType"))
"""
function MatDenseRestoreArrayWriteAndMemType(petsclib::PetscLibType, A::AbstractPetscMat, array::AbstractArray{PetscScalar}) end

@for_petsc function MatDenseRestoreArrayWriteAndMemType(petsclib::$UnionPetscLib, A::AbstractPetscMat, array::AbstractArray{$PetscScalar} )
	array_ = Ref(pointer(array))

    @chk ccall(
               (:MatDenseRestoreArrayWriteAndMemType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )


	return nothing
end 

"""
	MatDenseRestoreColumn(petsclib::PetscLibType,A::AbstractPetscMat, vals::AbstractArray{PetscScalar}) 
returns access to a column of a `MATDENSE` matrix which is returned by `MatDenseGetColumn()`.

Not Collective

Input Parameters:
- `A`    - a `MATSEQDENSE` or `MATMPIDENSE` matrix
- `vals` - pointer to the data (may be `NULL`)

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MatDenseGetColumn()`

# External Links
$(_doc_external("Mat/MatDenseRestoreColumn"))
"""
function MatDenseRestoreColumn(petsclib::PetscLibType, A::AbstractPetscMat, vals::AbstractArray{PetscScalar}) end

@for_petsc function MatDenseRestoreColumn(petsclib::$UnionPetscLib, A::AbstractPetscMat, vals::AbstractArray{$PetscScalar} )
	vals_ = Ref(pointer(vals))

    @chk ccall(
               (:MatDenseRestoreColumn, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, vals_,
              )


	return nothing
end 

"""
	MatDenseRestoreColumnVec(petsclib::PetscLibType,A::AbstractPetscMat, col::PetscInt, v::AbstractPetscVec) 
Returns access to a column of a dense matrix obtained from `MatDenseGetColumnVec()`.

Collective

Input Parameters:
- `A`   - the `Mat` object
- `col` - the column index
- `v`   - the `Vec` object (may be `NULL`)

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MATDENSECUDA`, `MATDENSEHIP`, `MatDenseGetColumnVec()`, `MatDenseGetColumnVecRead()`, `MatDenseGetColumnVecWrite()`, `MatDenseRestoreColumnVecRead()`, `MatDenseRestoreColumnVecWrite()`

# External Links
$(_doc_external("Mat/MatDenseRestoreColumnVec"))
"""
function MatDenseRestoreColumnVec(petsclib::PetscLibType, A::AbstractPetscMat, col::PetscInt, v::AbstractPetscVec) end

@for_petsc function MatDenseRestoreColumnVec(petsclib::$UnionPetscLib, A::AbstractPetscMat, col::$PetscInt, v::AbstractPetscVec )
	v_ = Ref(v.ptr)

    @chk ccall(
               (:MatDenseRestoreColumnVec, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CVec}),
               A, col, v_,
              )

	v.ptr = v_[]

	return nothing
end 

"""
	MatDenseRestoreColumnVecRead(petsclib::PetscLibType,A::AbstractPetscMat, col::PetscInt, v::AbstractPetscVec) 
Returns access to a column of a dense matrix obtained from `MatDenseGetColumnVecRead()`.

Collective

Input Parameters:
- `A`   - the `Mat` object
- `col` - the column index
- `v`   - the `Vec` object (may be `NULL`)

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MATDENSECUDA`, `MATDENSEHIP`, `MatDenseGetColumnVec()`, `MatDenseGetColumnVecRead()`, `MatDenseGetColumnVecWrite()`, `MatDenseRestoreColumnVec()`, `MatDenseRestoreColumnVecWrite()`

# External Links
$(_doc_external("Mat/MatDenseRestoreColumnVecRead"))
"""
function MatDenseRestoreColumnVecRead(petsclib::PetscLibType, A::AbstractPetscMat, col::PetscInt, v::AbstractPetscVec) end

@for_petsc function MatDenseRestoreColumnVecRead(petsclib::$UnionPetscLib, A::AbstractPetscMat, col::$PetscInt, v::AbstractPetscVec )
	v_ = Ref(v.ptr)

    @chk ccall(
               (:MatDenseRestoreColumnVecRead, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CVec}),
               A, col, v_,
              )

	v.ptr = v_[]

	return nothing
end 

"""
	MatDenseRestoreColumnVecWrite(petsclib::PetscLibType,A::AbstractPetscMat, col::PetscInt, v::AbstractPetscVec) 
Returns access to a column of a dense matrix obtained from `MatDenseGetColumnVecWrite()`.

Collective

Input Parameters:
- `A`   - the `Mat` object
- `col` - the column index
- `v`   - the `Vec` object (may be `NULL`)

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MATDENSECUDA`, `MATDENSEHIP`, `MatDenseGetColumnVec()`, `MatDenseGetColumnVecRead()`, `MatDenseGetColumnVecWrite()`, `MatDenseRestoreColumnVec()`, `MatDenseRestoreColumnVecRead()`

# External Links
$(_doc_external("Mat/MatDenseRestoreColumnVecWrite"))
"""
function MatDenseRestoreColumnVecWrite(petsclib::PetscLibType, A::AbstractPetscMat, col::PetscInt, v::AbstractPetscVec) end

@for_petsc function MatDenseRestoreColumnVecWrite(petsclib::$UnionPetscLib, A::AbstractPetscMat, col::$PetscInt, v::AbstractPetscVec )
	v_ = Ref(v.ptr)

    @chk ccall(
               (:MatDenseRestoreColumnVecWrite, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CVec}),
               A, col, v_,
              )

	v.ptr = v_[]

	return nothing
end 

"""
	MatDenseRestoreSubMatrix(petsclib::PetscLibType,A::AbstractPetscMat, v::AbstractPetscMat) 
Returns access to a block of columns of a dense matrix obtained from `MatDenseGetSubMatrix()`.

Collective

Input Parameters:
- `A` - the `Mat` object
- `v` - the `Mat` object (may be `NULL`)

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MATDENSECUDA`, `MATDENSEHIP`, `MatDenseGetColumnVec()`, `MatDenseRestoreColumnVec()`, `MatDenseGetSubMatrix()`

# External Links
$(_doc_external("Mat/MatDenseRestoreSubMatrix"))
"""
function MatDenseRestoreSubMatrix(petsclib::PetscLibType, A::AbstractPetscMat, v::AbstractPetscMat) end

@for_petsc function MatDenseRestoreSubMatrix(petsclib::$UnionPetscLib, A::AbstractPetscMat, v::AbstractPetscMat )
	v_ = Ref(v.ptr)

    @chk ccall(
               (:MatDenseRestoreSubMatrix, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, v_,
              )

	v.ptr = v_[]

	return nothing
end 

"""
	MatDenseSetLDA(petsclib::PetscLibType,A::AbstractPetscMat, lda::PetscInt) 
Sets the leading dimension of the array used by the `MATDENSE` matrix

Collective if the matrix layouts have not yet been setup

Input Parameters:
- `A`   - a `MATDENSE` or `MATDENSECUDA` matrix
- `lda` - the leading dimension

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MATDENSECUDA`, `MatDenseGetArray()`, `MatDenseRestoreArray()`, `MatDenseGetArrayRead()`, `MatDenseRestoreArrayRead()`, `MatDenseGetLDA()`

# External Links
$(_doc_external("Mat/MatDenseSetLDA"))
"""
function MatDenseSetLDA(petsclib::PetscLibType, A::AbstractPetscMat, lda::PetscInt) end

@for_petsc function MatDenseSetLDA(petsclib::$UnionPetscLib, A::AbstractPetscMat, lda::$PetscInt )

    @chk ccall(
               (:MatDenseSetLDA, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt),
               A, lda,
              )


	return nothing
end 

"""
	MatDestroy(petsclib::PetscLibType,A::AbstractPetscMat) 
Frees space taken by a matrix.

Collective

Input Parameter:
- `A` - the matrix

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatCreate()`

# External Links
$(_doc_external("Mat/MatDestroy"))
"""
function MatDestroy(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatDestroy(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	A_ = Ref(A.ptr)

    @chk ccall(
               (:MatDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{CMat},),
               A_,
              )

	A.ptr = C_NULL

	return nothing
end 

"""
	MatDestroyMatrices(petsclib::PetscLibType,n::PetscInt, mat::AbstractArray{PetscMat}) 
Destroys an array of matrices

Collective

Input Parameters:
- `n`   - the number of local matrices
- `mat` - the matrices (this is a pointer to the array of matrices)

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatCreateSubMatrices()`, `MatDestroySubMatrices()`

# External Links
$(_doc_external("Mat/MatDestroyMatrices"))
"""
function MatDestroyMatrices(petsclib::PetscLibType, n::PetscInt, mat::AbstractArray{PetscMat}) end

@for_petsc function MatDestroyMatrices(petsclib::$UnionPetscLib, n::$PetscInt, mat::AbstractArray{PetscMat} )
	mat_ = Ref(pointer(mat))

    @chk ccall(
               (:MatDestroyMatrices, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{Ptr{CMat}}),
               n, mat_,
              )


	return nothing
end 

"""
	MatDestroySeqNonzeroStructure(petsclib::PetscLibType,mat::AbstractPetscMat) 
Destroys matrix obtained with `MatGetSeqNonzeroStructure()`.

Collective

Input Parameter:
- `mat` - the matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatGetSeqNonzeroStructure()`

# External Links
$(_doc_external("Mat/MatDestroySeqNonzeroStructure"))
"""
function MatDestroySeqNonzeroStructure(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatDestroySeqNonzeroStructure(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	mat_ = Ref(mat.ptr)

    @chk ccall(
               (:MatDestroySeqNonzeroStructure, $petsc_library),
               PetscErrorCode,
               (Ptr{CMat},),
               mat_,
              )

	mat.ptr = mat_[]

	return nothing
end 

"""
	MatDestroySubMatrices(petsclib::PetscLibType,n::PetscInt, mat::AbstractArray{PetscMat}) 
Destroys a set of matrices obtained with `MatCreateSubMatrices()`.

Collective

Input Parameters:
- `n`   - the number of local matrices
- `mat` - the matrices (this is a pointer to the array of matrices, to match the calling sequence of `MatCreateSubMatrices()`)

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatCreateSubMatrices()`, `MatDestroyMatrices()`

# External Links
$(_doc_external("Mat/MatDestroySubMatrices"))
"""
function MatDestroySubMatrices(petsclib::PetscLibType, n::PetscInt, mat::AbstractArray{PetscMat}) end

@for_petsc function MatDestroySubMatrices(petsclib::$UnionPetscLib, n::$PetscInt, mat::AbstractArray{PetscMat} )
	mat_ = Ref(pointer(mat))

    @chk ccall(
               (:MatDestroySubMatrices, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{Ptr{CMat}}),
               n, mat_,
              )


	return nothing
end 

"""
	diag::PetscVec = MatDiagonalGetDiagonal(petsclib::PetscLibType,A::AbstractPetscMat) 
Get the diagonal of a `MATDIAGONAL`

Input Parameter:
- `A` - the `MATDIAGONAL`

Output Parameter:
- `diag` - the `Vec` that defines the diagonal

Level: developer

-seealso: [](ch_matrices), `MATDIAGONAL`, `MatCreateDiagonal()`, `MatDiagonalRestoreDiagonal()`, `MatDiagonalGetInverseDiagonal()`, `MatGetDiagonal()`

# External Links
$(_doc_external("Mat/MatDiagonalGetDiagonal"))
"""
function MatDiagonalGetDiagonal(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatDiagonalGetDiagonal(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	diag_ = Ref{CVec}()

    @chk ccall(
               (:MatDiagonalGetDiagonal, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CVec}),
               A, diag_,
              )

	diag = PetscVec(diag_[], petsclib)

	return diag
end 

"""
	inv_diag::PetscVec = MatDiagonalGetInverseDiagonal(petsclib::PetscLibType,A::AbstractPetscMat) 
Get the inverse diagonal of a `MATDIAGONAL`

Input Parameter:
- `A` - the `MATDIAGONAL`

Output Parameter:
- `inv_diag` - the `Vec` that defines the inverse diagonal

Level: developer

-seealso: [](ch_matrices), `MATDIAGONAL`, `MatCreateDiagonal()`, `MatDiagonalRestoreInverseDiagonal()`, `MatDiagonalGetDiagonal()`, `MATLMVMBROYDEN`, `MatSolve()`

# External Links
$(_doc_external("Mat/MatDiagonalGetInverseDiagonal"))
"""
function MatDiagonalGetInverseDiagonal(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatDiagonalGetInverseDiagonal(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	inv_diag_ = Ref{CVec}()

    @chk ccall(
               (:MatDiagonalGetInverseDiagonal, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CVec}),
               A, inv_diag_,
              )

	inv_diag = PetscVec(inv_diag_[], petsclib)

	return inv_diag
end 

"""
	MatDiagonalRestoreDiagonal(petsclib::PetscLibType,A::AbstractPetscMat, diag::AbstractPetscVec) 
Restore the diagonal of a `MATDIAGONAL`

Input Parameters:
- `A`    - the `MATDIAGONAL`
- `diag` - the `Vec` obtained from `MatDiagonalGetDiagonal()`

Level: developer

-seealso: [](ch_matrices), `MATDIAGONAL`, `MatCreateDiagonal()`, `MatDiagonalGetDiagonal()`

# External Links
$(_doc_external("Mat/MatDiagonalRestoreDiagonal"))
"""
function MatDiagonalRestoreDiagonal(petsclib::PetscLibType, A::AbstractPetscMat, diag::AbstractPetscVec) end

@for_petsc function MatDiagonalRestoreDiagonal(petsclib::$UnionPetscLib, A::AbstractPetscMat, diag::AbstractPetscVec )
	diag_ = Ref(diag.ptr)

    @chk ccall(
               (:MatDiagonalRestoreDiagonal, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CVec}),
               A, diag_,
              )

	diag.ptr = diag_[]

	return nothing
end 

"""
	MatDiagonalRestoreInverseDiagonal(petsclib::PetscLibType,A::AbstractPetscMat, inv_diag::AbstractPetscVec) 
Restore the inverse diagonal of a `MATDIAGONAL`

Input Parameters:
- `A`        - the `MATDIAGONAL`
- `inv_diag` - the `Vec` obtained from `MatDiagonalGetInverseDiagonal()`

Level: developer

-seealso: [](ch_matrices), `MATDIAGONAL`, `MatCreateDiagonal()`, `MatDiagonalGetInverseDiagonal()`

# External Links
$(_doc_external("Mat/MatDiagonalRestoreInverseDiagonal"))
"""
function MatDiagonalRestoreInverseDiagonal(petsclib::PetscLibType, A::AbstractPetscMat, inv_diag::AbstractPetscVec) end

@for_petsc function MatDiagonalRestoreInverseDiagonal(petsclib::$UnionPetscLib, A::AbstractPetscMat, inv_diag::AbstractPetscVec )
	inv_diag_ = Ref(inv_diag.ptr)

    @chk ccall(
               (:MatDiagonalRestoreInverseDiagonal, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CVec}),
               A, inv_diag_,
              )

	inv_diag.ptr = inv_diag_[]

	return nothing
end 

"""
	MatDiagonalScale(petsclib::PetscLibType,mat::AbstractPetscMat, l::Union{Ptr, AbstractPetscVec}, r::Union{Ptr, AbstractPetscVec}) 
Scales a matrix on the left and right by diagonal
matrices that are stored as vectors.  Either of the two scaling
matrices can be `NULL`.

Collective

Input Parameters:
- `mat` - the matrix to be scaled
- `l`   - the left scaling vector (or `NULL`)
- `r`   - the right scaling vector (or `NULL`)

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatScale()`, `MatShift()`, `MatDiagonalSet()`

# External Links
$(_doc_external("Mat/MatDiagonalScale"))
"""
function MatDiagonalScale(petsclib::PetscLibType, mat::AbstractPetscMat, l::Union{Ptr, AbstractPetscVec}, r::Union{Ptr, AbstractPetscVec}) end

@for_petsc function MatDiagonalScale(petsclib::$UnionPetscLib, mat::AbstractPetscMat, l::Union{Ptr, AbstractPetscVec}, r::Union{Ptr, AbstractPetscVec} )

    @chk ccall(
               (:MatDiagonalScale, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               mat, l, r,
              )


	return nothing
end 

"""
	MatDiagonalScaleLocal(petsclib::PetscLibType,mat::AbstractPetscMat, diag::AbstractPetscVec) 
Scales columns of a matrix given the scaling values including the
ghosted ones.

Not Collective

Input Parameters:
- `mat`  - the matrix
- `diag` - the diagonal values, including ghost ones

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatDiagonalScale()`

# External Links
$(_doc_external("Mat/MatDiagonalScaleLocal"))
"""
function MatDiagonalScaleLocal(petsclib::PetscLibType, mat::AbstractPetscMat, diag::AbstractPetscVec) end

@for_petsc function MatDiagonalScaleLocal(petsclib::$UnionPetscLib, mat::AbstractPetscMat, diag::AbstractPetscVec )

    @chk ccall(
               (:MatDiagonalScaleLocal, $petsc_library),
               PetscErrorCode,
               (CMat, CVec),
               mat, diag,
              )


	return nothing
end 

"""
	MatDiagonalSet(petsclib::PetscLibType,Y::AbstractPetscMat, D::AbstractPetscVec, is::InsertMode) 
Computes `Y` = `Y` + `D`, where `D` is a diagonal matrix
that is represented as a vector. Or Y[i,i] = D[i] if `InsertMode` is
`INSERT_VALUES`.

Neighbor-wise Collective

Input Parameters:
- `Y`  - the input matrix
- `D`  - the diagonal matrix, represented as a vector
- `is` - `INSERT_VALUES` or `ADD_VALUES`

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatShift()`, `MatScale()`, `MatDiagonalScale()`

# External Links
$(_doc_external("Mat/MatDiagonalSet"))
"""
function MatDiagonalSet(petsclib::PetscLibType, Y::AbstractPetscMat, D::AbstractPetscVec, is::InsertMode) end

@for_petsc function MatDiagonalSet(petsclib::$UnionPetscLib, Y::AbstractPetscMat, D::AbstractPetscVec, is::InsertMode )

    @chk ccall(
               (:MatDiagonalSet, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, InsertMode),
               Y, D, is,
              )


	return nothing
end 

"""
	M::PetscMat = MatDuplicate(petsclib::PetscLibType,mat::AbstractPetscMat, op::MatDuplicateOption) 
Duplicates a matrix including the non

Collective

Input Parameters:
- `mat` - the matrix
- `op`  - One of `MAT_DO_NOT_COPY_VALUES`, `MAT_COPY_VALUES`, or `MAT_SHARE_NONZERO_PATTERN`.
See the manual page for `MatDuplicateOption()` for an explanation of these options.

Output Parameter:
- `M` - pointer to place new matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatCopy()`, `MatConvert()`, `MatDuplicateOption`

# External Links
$(_doc_external("Mat/MatDuplicate"))
"""
function MatDuplicate(petsclib::PetscLibType, mat::AbstractPetscMat, op::MatDuplicateOption) end

@for_petsc function MatDuplicate(petsclib::$UnionPetscLib, mat::AbstractPetscMat, op::MatDuplicateOption )
	M_ = Ref{CMat}()

    @chk ccall(
               (:MatDuplicate, $petsc_library),
               PetscErrorCode,
               (CMat, MatDuplicateOption, Ptr{CMat}),
               mat, op, M_,
              )

	M = PetscMat(M_[], petsclib)

	return M
end 

"""
	MatEliminateZeros(petsclib::PetscLibType,A::AbstractPetscMat, keep::PetscBool) 
eliminate the nondiagonal zero entries in place from the nonzero structure of a sparse `Mat` in place,
meaning the same memory is used for the matrix, and no new memory is allocated.

Collective

Input Parameters:
- `A`    - the matrix
- `keep` - if for a given row of `A`, the diagonal coefficient is zero, indicates whether it should be left in the structure or eliminated as well

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `MatCreateGraph()`, `MatFilter()`

# External Links
$(_doc_external("Mat/MatEliminateZeros"))
"""
function MatEliminateZeros(petsclib::PetscLibType, A::AbstractPetscMat, keep::PetscBool) end

@for_petsc function MatEliminateZeros(petsclib::$UnionPetscLib, A::AbstractPetscMat, keep::PetscBool )

    @chk ccall(
               (:MatEliminateZeros, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               A, keep,
              )


	return nothing
end 

"""
	flg::PetscBool = MatEqual(petsclib::PetscLibType,A::AbstractPetscMat, B::AbstractPetscMat) 
Compares two matrices.

Collective

Input Parameters:
- `A` - the first matrix
- `B` - the second matrix

Output Parameter:
- `flg` - `PETSC_TRUE` if the matrices are equal; `PETSC_FALSE` otherwise.

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatMultEqual()`

# External Links
$(_doc_external("Mat/MatEqual"))
"""
function MatEqual(petsclib::PetscLibType, A::AbstractPetscMat, B::AbstractPetscMat) end

@for_petsc function MatEqual(petsclib::$UnionPetscLib, A::AbstractPetscMat, B::AbstractPetscMat )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatEqual, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, Ptr{PetscBool}),
               A, B, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	MatFactorClearError(petsclib::PetscLibType,mat::AbstractPetscMat) 
clears the error code in a factorization

Logically Collective

Input Parameter:
- `mat` - the factored matrix

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatZeroEntries()`, `MatFactor()`, `MatGetFactor()`, `MatLUFactorSymbolic()`, `MatCholeskyFactorSymbolic()`, `MatFactorGetError()`, `MatFactorGetErrorZeroPivot()`,
`MatGetErrorCode()`, `MatFactorError`

# External Links
$(_doc_external("Mat/MatFactorClearError"))
"""
function MatFactorClearError(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatFactorClearError(petsclib::$UnionPetscLib, mat::AbstractPetscMat )

    @chk ccall(
               (:MatFactorClearError, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	status::MatFactorSchurStatus = MatFactorCreateSchurComplement(petsclib::PetscLibType,F::AbstractPetscMat, S::AbstractPetscMat) 
Create a Schur complement matrix object using Schur data computed during the factorization step

Logically Collective

Input Parameters:
- `F`      - the factored matrix obtained by calling `MatGetFactor()`
- `S`      - location where to return the Schur complement, can be `NULL`
- `status` - the status of the Schur complement matrix, can be `NULL`

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatGetFactor()`, `MatFactorSetSchurIS()`, `MatFactorGetSchurComplement()`, `MatFactorSchurStatus`, `MATSOLVERMUMPS`, `MATSOLVERMKL_PARDISO`

# External Links
$(_doc_external("Mat/MatFactorCreateSchurComplement"))
"""
function MatFactorCreateSchurComplement(petsclib::PetscLibType, F::AbstractPetscMat, S::AbstractPetscMat) end

@for_petsc function MatFactorCreateSchurComplement(petsclib::$UnionPetscLib, F::AbstractPetscMat, S::AbstractPetscMat )
	S_ = Ref(S.ptr)
	status_ = Ref{MatFactorSchurStatus}()

    @chk ccall(
               (:MatFactorCreateSchurComplement, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}, Ptr{MatFactorSchurStatus}),
               F, S_, status_,
              )

	S.ptr = S_[]
	status = status_[]

	return status
end 

"""
	MatFactorFactorizeSchurComplement(petsclib::PetscLibType,F::AbstractPetscMat) 
Factorize the Schur complement matrix computed during the factorization step

Logically Collective

Input Parameter:
- `F` - the factored matrix obtained by calling `MatGetFactor()`

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatGetFactor()`, `MatFactorSetSchurIS()`, `MatFactorInvertSchurComplement()`

# External Links
$(_doc_external("Mat/MatFactorFactorizeSchurComplement"))
"""
function MatFactorFactorizeSchurComplement(petsclib::PetscLibType, F::AbstractPetscMat) end

@for_petsc function MatFactorFactorizeSchurComplement(petsclib::$UnionPetscLib, F::AbstractPetscMat )

    @chk ccall(
               (:MatFactorFactorizeSchurComplement, $petsc_library),
               PetscErrorCode,
               (CMat,),
               F,
              )


	return nothing
end 

"""
	flg::PetscBool = MatFactorGetCanUseOrdering(petsclib::PetscLibType,mat::AbstractPetscMat) 
Indicates if the factorization can use the ordering provided in `MatLUFactorSymbolic()`, `MatCholeskyFactorSymbolic()`

Logically Collective

Input Parameter:
- `mat` - the matrix

Output Parameter:
- `flg` - `PETSC_TRUE` if uses the ordering

Level: developer

-seealso: [](ch_matrices), `Mat`, [Matrix Factorization](sec_matfactor), `MatCopy()`, `MatDuplicate()`, `MatGetFactorAvailable()`, `MatGetFactor()`, `MatLUFactorSymbolic()`, `MatCholeskyFactorSymbolic()`

# External Links
$(_doc_external("Mat/MatFactorGetCanUseOrdering"))
"""
function MatFactorGetCanUseOrdering(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatFactorGetCanUseOrdering(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatFactorGetCanUseOrdering, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscBool}),
               mat, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	err::MatFactorError = MatFactorGetError(petsclib::PetscLibType,mat::AbstractPetscMat) 
gets the error code from a factorization

Logically Collective

Input Parameter:
- `mat` - the factored matrix

Output Parameter:
- `err` - the error code

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatZeroEntries()`, `MatFactor()`, `MatGetFactor()`, `MatLUFactorSymbolic()`, `MatCholeskyFactorSymbolic()`,
`MatFactorClearError()`, `MatFactorGetErrorZeroPivot()`, `MatFactorError`

# External Links
$(_doc_external("Mat/MatFactorGetError"))
"""
function MatFactorGetError(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatFactorGetError(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	err_ = Ref{MatFactorError}()

    @chk ccall(
               (:MatFactorGetError, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatFactorError}),
               mat, err_,
              )

	err = err_[]

	return err
end 

"""
	pivot::PetscReal,row::PetscInt = MatFactorGetErrorZeroPivot(petsclib::PetscLibType,mat::AbstractPetscMat) 
returns the pivot value that was determined to be zero and the row it occurred in

Logically Collective

Input Parameter:
- `mat` - the factored matrix

Output Parameters:
- `pivot` - the pivot value computed
- `row`   - the row that the zero pivot occurred. This row value must be interpreted carefully due to row reorderings and which processes
the share the matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatZeroEntries()`, `MatFactor()`, `MatGetFactor()`,
`MatLUFactorSymbolic()`, `MatCholeskyFactorSymbolic()`, `MatFactorClearError()`,
`MAT_FACTOR_NUMERIC_ZEROPIVOT`

# External Links
$(_doc_external("Mat/MatFactorGetErrorZeroPivot"))
"""
function MatFactorGetErrorZeroPivot(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatFactorGetErrorZeroPivot(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	pivot_ = Ref{$PetscReal}()
	row_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatFactorGetErrorZeroPivot, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscReal}, Ptr{$PetscInt}),
               mat, pivot_, row_,
              )

	pivot = pivot_[]
	row = row_[]

	return pivot,row
end 

"""
	otype::MatOrderingType = MatFactorGetPreferredOrdering(petsclib::PetscLibType,mat::AbstractPetscMat, ftype::MatFactorType) 
The preferred ordering for a particular matrix factor object

Logically Collective

Input Parameters:
- `mat`   - the matrix obtained with `MatGetFactor()`
- `ftype` - the factorization type to be used

Output Parameter:
- `otype` - the preferred ordering type

Level: developer

-seealso: [](ch_matrices), `Mat`, [Matrix Factorization](sec_matfactor), `MatFactorType`, `MatOrderingType`, `MatCopy()`, `MatDuplicate()`, `MatGetFactorAvailable()`, `MatGetFactor()`, `MatLUFactorSymbolic()`, `MatCholeskyFactorSymbolic()`

# External Links
$(_doc_external("Mat/MatFactorGetPreferredOrdering"))
"""
function MatFactorGetPreferredOrdering(petsclib::PetscLibType, mat::AbstractPetscMat, ftype::MatFactorType) end

@for_petsc function MatFactorGetPreferredOrdering(petsclib::$UnionPetscLib, mat::AbstractPetscMat, ftype::MatFactorType )
	otype_ = Ref{MatOrderingType}()

    @chk ccall(
               (:MatFactorGetPreferredOrdering, $petsc_library),
               PetscErrorCode,
               (CMat, MatFactorType, Ptr{MatOrderingType}),
               mat, ftype, otype_,
              )

	otype = otype_[] == C_NULL ? "" : unsafe_string(otype_[])

	return otype
end 

"""
	status::MatFactorSchurStatus = MatFactorGetSchurComplement(petsclib::PetscLibType,F::AbstractPetscMat, S::Union{Ptr, AbstractPetscMat}) 
Gets access to a Schur complement matrix using the current Schur data within a factored matrix

Logically Collective

Input Parameters:
- `F`      - the factored matrix obtained by calling `MatGetFactor()`
- `S`      - location where to return the Schur complement, can be `NULL`
- `status` - the status of the Schur complement matrix, can be `NULL`

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatGetFactor()`, `MatFactorSetSchurIS()`, `MatFactorRestoreSchurComplement()`, `MatFactorCreateSchurComplement()`, `MatFactorSchurStatus`

# External Links
$(_doc_external("Mat/MatFactorGetSchurComplement"))
"""
function MatFactorGetSchurComplement(petsclib::PetscLibType, F::AbstractPetscMat, S::Union{Ptr, AbstractPetscMat}) end

@for_petsc function MatFactorGetSchurComplement(petsclib::$UnionPetscLib, F::AbstractPetscMat, S::Union{Ptr, AbstractPetscMat} )
	S_ = Ref(S.ptr)
	status_ = Ref{MatFactorSchurStatus}()

    @chk ccall(
               (:MatFactorGetSchurComplement, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}, Ptr{MatFactorSchurStatus}),
               F, S_, status_,
              )

	S.ptr = S_[]
	status = status_[]

	return status
end 

"""
	type::MatSolverType = MatFactorGetSolverType(petsclib::PetscLibType,mat::AbstractPetscMat) 
Returns name of the package providing the factorization routines

Not Collective

Input Parameter:
- `mat` - the matrix, must be a factored matrix

Output Parameter:
- `type` - the string name of the package (do not free this string)

Level: intermediate

-seealso: [](ch_matrices), `Mat`, [Matrix Factorization](sec_matfactor), `MatGetFactor()`, `MatSolverType`, `MatCopy()`, `MatDuplicate()`, `MatGetFactorAvailable()`

# External Links
$(_doc_external("Mat/MatFactorGetSolverType"))
"""
function MatFactorGetSolverType(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatFactorGetSolverType(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	type_ = Ref{MatSolverType}()

    @chk ccall(
               (:MatFactorGetSolverType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatSolverType}),
               mat, type_,
              )

	type = type_[] == C_NULL ? "" : unsafe_string(type_[])

	return type
end 

"""
	MatFactorInfoInitialize(petsclib::PetscLibType,info::Vector{MatFactorInfo}) 
Initializes a `MatFactorInfo` data structure
with default values.

Not Collective

Input Parameter:
- `info` - the `MatFactorInfo` data structure

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatGetFactor()`, `MatFactorInfo`

# External Links
$(_doc_external("Mat/MatFactorInfoInitialize"))
"""
function MatFactorInfoInitialize(petsclib::PetscLibType, info::Vector{MatFactorInfo}) end

@for_petsc function MatFactorInfoInitialize(petsclib::$UnionPetscLib, info::Vector{MatFactorInfo} )

    @chk ccall(
               (:MatFactorInfoInitialize, $petsc_library),
               PetscErrorCode,
               (Ptr{MatFactorInfo},),
               info,
              )


	return nothing
end 

"""
	MatFactorInvertSchurComplement(petsclib::PetscLibType,F::AbstractPetscMat) 
Invert the Schur complement matrix computed during the factorization step

Logically Collective

Input Parameter:
- `F` - the factored matrix obtained by calling `MatGetFactor()`

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatGetFactor()`, `MatFactorSetSchurIS()`, `MatFactorGetSchurComplement()`, `MatFactorCreateSchurComplement()`

# External Links
$(_doc_external("Mat/MatFactorInvertSchurComplement"))
"""
function MatFactorInvertSchurComplement(petsclib::PetscLibType, F::AbstractPetscMat) end

@for_petsc function MatFactorInvertSchurComplement(petsclib::$UnionPetscLib, F::AbstractPetscMat )

    @chk ccall(
               (:MatFactorInvertSchurComplement, $petsc_library),
               PetscErrorCode,
               (CMat,),
               F,
              )


	return nothing
end 

"""
	MatFactorRestoreSchurComplement(petsclib::PetscLibType,F::AbstractPetscMat, S::AbstractPetscMat, status::MatFactorSchurStatus) 
Restore the Schur complement matrix object obtained from a call to `MatFactorGetSchurComplement()`

Logically Collective

Input Parameters:
- `F`      - the factored matrix obtained by calling `MatGetFactor()`
- `S`      - location where the Schur complement is stored
- `status` - the status of the Schur complement matrix (see `MatFactorSchurStatus`)

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatGetFactor()`, `MatFactorSetSchurIS()`, `MatFactorCreateSchurComplement()`, `MatFactorSchurStatus`

# External Links
$(_doc_external("Mat/MatFactorRestoreSchurComplement"))
"""
function MatFactorRestoreSchurComplement(petsclib::PetscLibType, F::AbstractPetscMat, S::AbstractPetscMat, status::MatFactorSchurStatus) end

@for_petsc function MatFactorRestoreSchurComplement(petsclib::$UnionPetscLib, F::AbstractPetscMat, S::AbstractPetscMat, status::MatFactorSchurStatus )
	S_ = Ref(S.ptr)

    @chk ccall(
               (:MatFactorRestoreSchurComplement, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}, MatFactorSchurStatus),
               F, S_, status,
              )

	S.ptr = S_[]

	return nothing
end 

"""
	MatFactorSetSchurIS(petsclib::PetscLibType,mat::AbstractPetscMat, is::AbstractIS) 
Set indices corresponding to the Schur complement you wish to have computed

Collective

Input Parameters:
- `mat` - the factored matrix
- `is`  - the index set defining the Schur indices (0-based)

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatGetFactor()`, `MatFactorGetSchurComplement()`, `MatFactorRestoreSchurComplement()`, `MatFactorCreateSchurComplement()`, `MatFactorSolveSchurComplement()`,
`MatFactorSolveSchurComplementTranspose()`, `MATSOLVERMUMPS`, `MATSOLVERMKL_PARDISO`

# External Links
$(_doc_external("Mat/MatFactorSetSchurIS"))
"""
function MatFactorSetSchurIS(petsclib::PetscLibType, mat::AbstractPetscMat, is::AbstractIS) end

@for_petsc function MatFactorSetSchurIS(petsclib::$UnionPetscLib, mat::AbstractPetscMat, is::AbstractIS )

    @chk ccall(
               (:MatFactorSetSchurIS, $petsc_library),
               PetscErrorCode,
               (CMat, CIS),
               mat, is,
              )


	return nothing
end 

"""
	MatFactorSolveSchurComplement(petsclib::PetscLibType,F::AbstractPetscMat, rhs::AbstractPetscVec, sol::AbstractPetscVec) 
Solve the Schur complement system computed during the factorization step

Logically Collective

Input Parameters:
- `F`   - the factored matrix obtained by calling `MatGetFactor()`
- `rhs` - location where the right-hand side of the Schur complement system is stored
- `sol` - location where the solution of the Schur complement system has to be returned

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatGetFactor()`, `MatFactorSetSchurIS()`, `MatFactorSolveSchurComplementTranspose()`

# External Links
$(_doc_external("Mat/MatFactorSolveSchurComplement"))
"""
function MatFactorSolveSchurComplement(petsclib::PetscLibType, F::AbstractPetscMat, rhs::AbstractPetscVec, sol::AbstractPetscVec) end

@for_petsc function MatFactorSolveSchurComplement(petsclib::$UnionPetscLib, F::AbstractPetscMat, rhs::AbstractPetscVec, sol::AbstractPetscVec )

    @chk ccall(
               (:MatFactorSolveSchurComplement, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               F, rhs, sol,
              )


	return nothing
end 

"""
	MatFactorSolveSchurComplementTranspose(petsclib::PetscLibType,F::AbstractPetscMat, rhs::AbstractPetscVec, sol::AbstractPetscVec) 
Solve the transpose of the Schur complement system computed during the factorization step

Logically Collective

Input Parameters:
- `F`   - the factored matrix obtained by calling `MatGetFactor()`
- `rhs` - location where the right-hand side of the Schur complement system is stored
- `sol` - location where the solution of the Schur complement system has to be returned

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatGetFactor()`, `MatFactorSetSchurIS()`, `MatFactorSolveSchurComplement()`

# External Links
$(_doc_external("Mat/MatFactorSolveSchurComplementTranspose"))
"""
function MatFactorSolveSchurComplementTranspose(petsclib::PetscLibType, F::AbstractPetscMat, rhs::AbstractPetscVec, sol::AbstractPetscVec) end

@for_petsc function MatFactorSolveSchurComplementTranspose(petsclib::$UnionPetscLib, F::AbstractPetscMat, rhs::AbstractPetscVec, sol::AbstractPetscVec )

    @chk ccall(
               (:MatFactorSolveSchurComplementTranspose, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               F, rhs, sol,
              )


	return nothing
end 

"""
	MatFilter(petsclib::PetscLibType,A::AbstractPetscMat, tol::PetscReal, compress::PetscBool, keep::PetscBool) 
Set all values in the matrix with an absolute value less than or equal to the tolerance to zero, and optionally compress the underlying storage

Input Parameters:
- `A`        - The matrix
- `tol`      - The zero tolerance
- `compress` - Whether the storage from the input matrix `A` should be compressed once values less than or equal to `tol` are set to zero
- `keep`     - If `compress` is true and for a given row of `A`, the diagonal coefficient is less than or equal to `tol`, indicates whether it should be left in the structure or eliminated as well

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `MatZeroEntries()`, `MatEliminateZeros()`, `VecFilter()`

# External Links
$(_doc_external("Mat/MatFilter"))
"""
function MatFilter(petsclib::PetscLibType, A::AbstractPetscMat, tol::PetscReal, compress::PetscBool, keep::PetscBool) end

@for_petsc function MatFilter(petsclib::$UnionPetscLib, A::AbstractPetscMat, tol::$PetscReal, compress::PetscBool, keep::PetscBool )

    @chk ccall(
               (:MatFilter, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscReal, PetscBool, PetscBool),
               A, tol, compress, keep,
              )


	return nothing
end 

"""
	MatFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the PETSc interface to the `Mat`
package. It is called from `PetscFinalize()`.

Level: developer

-seealso: `Mat`, `PetscFinalize()`, `MatInitializePackage()`

# External Links
$(_doc_external("Mat/MatFinalizePackage"))
"""
function MatFinalizePackage(petsclib::PetscLibType) end

@for_petsc function MatFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:MatFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	keptrows::IS = MatFindNonzeroRows(petsclib::PetscLibType,mat::AbstractPetscMat) 
Locate all rows that are not completely zero in the matrix

Input Parameter:
- `mat` - the matrix

Output Parameter:
- `keptrows` - the rows that are not completely zero

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatFindZeroRows()`

# External Links
$(_doc_external("Mat/MatFindNonzeroRows"))
"""
function MatFindNonzeroRows(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatFindNonzeroRows(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	keptrows_ = Ref{CIS}()

    @chk ccall(
               (:MatFindNonzeroRows, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CIS}),
               mat, keptrows_,
              )

	keptrows = IS(keptrows_[], petsclib)

	return keptrows
end 

"""
	is::IS = MatFindOffBlockDiagonalEntries(petsclib::PetscLibType,mat::AbstractPetscMat) 
Finds all the rows of a matrix that have entries outside of the main diagonal block (defined by the matrix block size)

Collective

Input Parameter:
- `mat` - the matrix

Output Parameter:
- `is` - contains the list of rows with off block diagonal entries

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatMultTranspose()`, `MatMultAdd()`, `MatMultTransposeAdd()`

# External Links
$(_doc_external("Mat/MatFindOffBlockDiagonalEntries"))
"""
function MatFindOffBlockDiagonalEntries(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatFindOffBlockDiagonalEntries(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	is_ = Ref{CIS}()

    @chk ccall(
               (:MatFindOffBlockDiagonalEntries, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CIS}),
               mat, is_,
              )

	is = IS(is_[], petsclib)

	return is
end 

"""
	is::IS = MatFindZeroDiagonals(petsclib::PetscLibType,mat::AbstractPetscMat) 
Finds all the rows of a matrix that have zero or no diagonal entry in the matrix

Collective

Input Parameter:
- `mat` - the matrix

Output Parameter:
- `is` - if any rows have zero diagonals this contains the list of them

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatMultTranspose()`, `MatMultAdd()`, `MatMultTransposeAdd()`

# External Links
$(_doc_external("Mat/MatFindZeroDiagonals"))
"""
function MatFindZeroDiagonals(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatFindZeroDiagonals(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	is_ = Ref{CIS}()

    @chk ccall(
               (:MatFindZeroDiagonals, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CIS}),
               mat, is_,
              )

	is = IS(is_[], petsclib)

	return is
end 

"""
	zerorows::IS = MatFindZeroRows(petsclib::PetscLibType,mat::AbstractPetscMat) 
Locate all rows that are completely zero in the matrix

Input Parameter:
- `mat` - the matrix

Output Parameter:
- `zerorows` - the rows that are completely zero

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatFindNonzeroRows()`

# External Links
$(_doc_external("Mat/MatFindZeroRows"))
"""
function MatFindZeroRows(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatFindZeroRows(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	zerorows_ = Ref{CIS}()

    @chk ccall(
               (:MatFindZeroRows, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CIS}),
               mat, zerorows_,
              )

	zerorows = IS(zerorows_[], petsclib)

	return zerorows
end 

"""
	MatForwardSolve(petsclib::PetscLibType,mat::AbstractPetscMat, b::AbstractPetscVec, x::AbstractPetscVec) 
Solves  L x = b , given a factored matrix, A = LU , or
U^T*D^(1/2) x = b, given a factored symmetric matrix, A = U^T*D*U,

Neighbor-wise Collective

Input Parameters:
- `mat` - the factored matrix
- `b`   - the right-hand-side vector

Output Parameter:
- `x` - the result vector

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatBackwardSolve()`, `MatGetFactor()`, `MatSolve()`

# External Links
$(_doc_external("Mat/MatForwardSolve"))
"""
function MatForwardSolve(petsclib::PetscLibType, mat::AbstractPetscMat, b::AbstractPetscVec, x::AbstractPetscVec) end

@for_petsc function MatForwardSolve(petsclib::$UnionPetscLib, mat::AbstractPetscMat, b::AbstractPetscVec, x::AbstractPetscVec )

    @chk ccall(
               (:MatForwardSolve, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               mat, b, x,
              )


	return nothing
end 

"""
	A::PetscMat = MatGalerkin(petsclib::PetscLibType,restrct::AbstractPetscMat, dA::AbstractPetscMat, interpolate::AbstractPetscMat, reuse::MatReuse, fill::PetscReal) 
Constructs the coarse grid problem matrix via Galerkin projection.

If the interpolation and restriction operators are the same, uses `MatPtAP()`.
If they are not the same, uses `MatMatMatMult()`.

Once the coarse grid problem is constructed, correct for interpolation operators
that are not of full rank, which can legitimately happen in the case of non-nested
geometric multigrid.

Input Parameters:
- `restrct`     - restriction operator
- `dA`          - fine grid matrix
- `interpolate` - interpolation operator
- `reuse`       - either `MAT_INITIAL_MATRIX` or `MAT_REUSE_MATRIX`
- `fill`        - expected fill, use `PETSC_DETERMINE` or `PETSC_DETERMINE` if you do not have a good estimate

Output Parameter:
- `A` - the Galerkin coarse matrix

Options Database Key:
- `-pc_mg_galerkin <both,pmat,mat,none>` - for what matrices the Galerkin process should be used

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatPtAP()`, `MatMatMatMult()`

# External Links
$(_doc_external("Mat/MatGalerkin"))
"""
function MatGalerkin(petsclib::PetscLibType, restrct::AbstractPetscMat, dA::AbstractPetscMat, interpolate::AbstractPetscMat, reuse::MatReuse, fill::PetscReal) end

@for_petsc function MatGalerkin(petsclib::$UnionPetscLib, restrct::AbstractPetscMat, dA::AbstractPetscMat, interpolate::AbstractPetscMat, reuse::MatReuse, fill::$PetscReal )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatGalerkin, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat, MatReuse, $PetscReal, Ptr{CMat}),
               restrct, dA, interpolate, reuse, fill, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	flg::PetscBool = MatGetBindingPropagates(petsclib::PetscLibType,A::AbstractPetscMat) 
Gets whether the state of being bound to the CPU for a GPU matrix type propagates to child and some other associated objects

Input Parameter:
- `A` - the matrix

Output Parameter:
- `flg` - flag indicating whether the boundtocpu flag will be propagated

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatSetBindingPropagates()`

# External Links
$(_doc_external("Mat/MatGetBindingPropagates"))
"""
function MatGetBindingPropagates(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatGetBindingPropagates(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatGetBindingPropagates, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscBool}),
               A, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	bs::PetscInt = MatGetBlockSize(petsclib::PetscLibType,mat::AbstractPetscMat) 
Returns the matrix block size.

Not Collective

Input Parameter:
- `mat` - the matrix

Output Parameter:
- `bs` - block size

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATBAIJ`, `MATSBAIJ`, `MatCreateSeqBAIJ()`, `MatCreateBAIJ()`, `MatGetBlockSizes()`

# External Links
$(_doc_external("Mat/MatGetBlockSize"))
"""
function MatGetBlockSize(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatGetBlockSize(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	bs_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatGetBlockSize, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}),
               mat, bs_,
              )

	bs = bs_[]

	return bs
end 

"""
	rbs::PetscInt,cbs::PetscInt = MatGetBlockSizes(petsclib::PetscLibType,mat::AbstractPetscMat) 
Returns the matrix block row and column sizes.

Not Collective

Input Parameter:
- `mat` - the matrix

Output Parameters:
- `rbs` - row block size
- `cbs` - column block size

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATBAIJ`, `MATSBAIJ`, `MatCreateSeqBAIJ()`, `MatCreateBAIJ()`, `MatGetBlockSize()`, `MatSetBlockSize()`, `MatSetBlockSizes()`

# External Links
$(_doc_external("Mat/MatGetBlockSizes"))
"""
function MatGetBlockSizes(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatGetBlockSizes(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	rbs_ = Ref{$PetscInt}()
	cbs_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatGetBlockSizes, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}, Ptr{$PetscInt}),
               mat, rbs_, cbs_,
              )

	rbs = rbs_[]
	cbs = cbs_[]

	return rbs,cbs
end 

"""
	rowb::IS,colb::IS,B_seq::PetscMat = MatGetBrowsOfAcols(petsclib::PetscLibType,A::AbstractPetscMat, B::AbstractPetscMat, scall::MatReuse) 
Returns `IS` that contain rows of `B` that equal to nonzero columns of local `A`

Collective

Input Parameters:
- `A`     - the first matrix in `MATMPIAIJ` format
- `B`     - the second matrix in `MATMPIAIJ` format
- `scall` - either `MAT_INITIAL_MATRIX` or `MAT_REUSE_MATRIX`

Output Parameters:
- `rowb`  - On input index sets of rows of B to extract (or `NULL`), modified on output
- `colb`  - On input index sets of columns of B to extract (or `NULL`), modified on output
- `B_seq` - the sequential matrix generated

Level: developer

-seealso: `Mat`, `MATMPIAIJ`, `IS`, `MatReuse`

# External Links
$(_doc_external("Mat/MatGetBrowsOfAcols"))
"""
function MatGetBrowsOfAcols(petsclib::PetscLibType, A::AbstractPetscMat, B::AbstractPetscMat, scall::MatReuse) end

@for_petsc function MatGetBrowsOfAcols(petsclib::$UnionPetscLib, A::AbstractPetscMat, B::AbstractPetscMat, scall::MatReuse )
	rowb_ = Ref{CIS}()
	colb_ = Ref{CIS}()
	B_seq_ = Ref{CMat}()

    @chk ccall(
               (:MatGetBrowsOfAcols, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, MatReuse, Ptr{CIS}, Ptr{CIS}, Ptr{CMat}),
               A, B, scall, rowb_, colb_, B_seq_,
              )

	rowb = IS(rowb_[], petsclib)
	colb = IS(colb_[], petsclib)
	B_seq = PetscMat(B_seq_[], petsclib)

	return rowb,colb,B_seq
end 

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

"""
	MatGetColumnMeans(petsclib::PetscLibType,A::AbstractPetscMat, means::Vector{PetscScalar}) 
Gets the arithmetic means of each column of a sparse or dense matrix.

Input Parameter:
- `A` - the matrix

Output Parameter:
- `means` - an array as large as the TOTAL number of columns in the matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `VecSum()`, `MatGetColumnSums()`, `MatGetColumnNorms()`, `MatGetColumnReductions()`

# External Links
$(_doc_external("Mat/MatGetColumnMeans"))
"""
function MatGetColumnMeans(petsclib::PetscLibType, A::AbstractPetscMat, means::Vector{PetscScalar}) end

@for_petsc function MatGetColumnMeans(petsclib::$UnionPetscLib, A::AbstractPetscMat, means::Vector{$PetscScalar} )

    @chk ccall(
               (:MatGetColumnMeans, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscScalar}),
               A, means,
              )


	return nothing
end 

"""
	MatGetColumnMeansImaginaryPart(petsclib::PetscLibType,A::AbstractPetscMat, means::Vector{PetscReal}) 
Gets the arithmetic means of the imaginary part of each column of a sparse or dense matrix.

Input Parameter:
- `A` - the matrix

Output Parameter:
- `means` - an array as large as the TOTAL number of columns in the matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatGetColumnMeansRealPart()`, `VecSum()`, `MatGetColumnSums()`, `MatGetColumnNorms()`, `MatGetColumnReductions()`

# External Links
$(_doc_external("Mat/MatGetColumnMeansImaginaryPart"))
"""
function MatGetColumnMeansImaginaryPart(petsclib::PetscLibType, A::AbstractPetscMat, means::Vector{PetscReal}) end

@for_petsc function MatGetColumnMeansImaginaryPart(petsclib::$UnionPetscLib, A::AbstractPetscMat, means::Vector{$PetscReal} )

    @chk ccall(
               (:MatGetColumnMeansImaginaryPart, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscReal}),
               A, means,
              )


	return nothing
end 

"""
	MatGetColumnMeansRealPart(petsclib::PetscLibType,A::AbstractPetscMat, means::Vector{PetscReal}) 
Gets the arithmetic means of the real part of each column of a sparse or dense matrix.

Input Parameter:
- `A` - the matrix

Output Parameter:
- `means` - an array as large as the TOTAL number of columns in the matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatGetColumnMeansImaginaryPart()`, `VecSum()`, `MatGetColumnSums()`, `MatGetColumnNorms()`, `MatGetColumnReductions()`

# External Links
$(_doc_external("Mat/MatGetColumnMeansRealPart"))
"""
function MatGetColumnMeansRealPart(petsclib::PetscLibType, A::AbstractPetscMat, means::Vector{PetscReal}) end

@for_petsc function MatGetColumnMeansRealPart(petsclib::$UnionPetscLib, A::AbstractPetscMat, means::Vector{$PetscReal} )

    @chk ccall(
               (:MatGetColumnMeansRealPart, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscReal}),
               A, means,
              )


	return nothing
end 

"""
	MatGetColumnNorms(petsclib::PetscLibType,A::AbstractPetscMat, type::NormType, norms::Vector{PetscReal}) 
Gets the norms of each column of a sparse or dense matrix.

Input Parameters:
- `A`    - the matrix
- `type` - `NORM_2`, `NORM_1` or `NORM_INFINITY`

Output Parameter:
- `norms` - an array as large as the TOTAL number of columns in the matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `NormType`, `MatNorm()`

# External Links
$(_doc_external("Mat/MatGetColumnNorms"))
"""
function MatGetColumnNorms(petsclib::PetscLibType, A::AbstractPetscMat, type::NormType, norms::Vector{PetscReal}) end

@for_petsc function MatGetColumnNorms(petsclib::$UnionPetscLib, A::AbstractPetscMat, type::NormType, norms::Vector{$PetscReal} )

    @chk ccall(
               (:MatGetColumnNorms, $petsc_library),
               PetscErrorCode,
               (CMat, NormType, Ptr{$PetscReal}),
               A, type, norms,
              )


	return nothing
end 

"""
	MatGetColumnReductions(petsclib::PetscLibType,A::AbstractPetscMat, type::PetscInt, reductions::Vector{PetscReal}) 
Gets the reductions of each column of a sparse or dense matrix.

Input Parameters:
- `A`    - the matrix
- `type` - A constant defined in `NormType` or `ReductionType`: `NORM_2`, `NORM_1`, `NORM_INFINITY`, `REDUCTION_SUM_REALPART`,
`REDUCTION_SUM_IMAGINARYPART`, `REDUCTION_MEAN_REALPART`, `REDUCTION_MEAN_IMAGINARYPART`

Output Parameter:
- `reductions` - an array as large as the TOTAL number of columns in the matrix

Level: developer

-seealso: [](ch_matrices), `Mat`, `ReductionType`, `NormType`, `MatGetColumnNorms()`, `MatGetColumnSums()`, `MatGetColumnMeans()`

# External Links
$(_doc_external("Mat/MatGetColumnReductions"))
"""
function MatGetColumnReductions(petsclib::PetscLibType, A::AbstractPetscMat, type::PetscInt, reductions::Vector{PetscReal}) end

@for_petsc function MatGetColumnReductions(petsclib::$UnionPetscLib, A::AbstractPetscMat, type::$PetscInt, reductions::Vector{$PetscReal} )

    @chk ccall(
               (:MatGetColumnReductions, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscReal}),
               A, type, reductions,
              )


	return nothing
end 

"""
	MatGetColumnSums(petsclib::PetscLibType,A::AbstractPetscMat, sums::Vector{PetscScalar}) 
Gets the sums of each column of a sparse or dense matrix.

Input Parameter:
- `A` - the matrix

Output Parameter:
- `sums` - an array as large as the TOTAL number of columns in the matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `VecSum()`, `MatGetColumnMeans()`, `MatGetColumnNorms()`, `MatGetColumnReductions()`

# External Links
$(_doc_external("Mat/MatGetColumnSums"))
"""
function MatGetColumnSums(petsclib::PetscLibType, A::AbstractPetscMat, sums::Vector{PetscScalar}) end

@for_petsc function MatGetColumnSums(petsclib::$UnionPetscLib, A::AbstractPetscMat, sums::Vector{$PetscScalar} )

    @chk ccall(
               (:MatGetColumnSums, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscScalar}),
               A, sums,
              )


	return nothing
end 

"""
	MatGetColumnSumsImaginaryPart(petsclib::PetscLibType,A::AbstractPetscMat, sums::Vector{PetscReal}) 
Gets the sums of the imaginary part of each column of a sparse or dense matrix.

Input Parameter:
- `A` - the matrix

Output Parameter:
- `sums` - an array as large as the TOTAL number of columns in the matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatGetColumnSumsRealPart()`, `VecSum()`, `MatGetColumnMeans()`, `MatGetColumnNorms()`, `MatGetColumnReductions()`

# External Links
$(_doc_external("Mat/MatGetColumnSumsImaginaryPart"))
"""
function MatGetColumnSumsImaginaryPart(petsclib::PetscLibType, A::AbstractPetscMat, sums::Vector{PetscReal}) end

@for_petsc function MatGetColumnSumsImaginaryPart(petsclib::$UnionPetscLib, A::AbstractPetscMat, sums::Vector{$PetscReal} )

    @chk ccall(
               (:MatGetColumnSumsImaginaryPart, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscReal}),
               A, sums,
              )


	return nothing
end 

"""
	MatGetColumnSumsRealPart(petsclib::PetscLibType,A::AbstractPetscMat, sums::Vector{PetscReal}) 
Gets the sums of the real part of each column of a sparse or dense matrix.

Input Parameter:
- `A` - the matrix

Output Parameter:
- `sums` - an array as large as the TOTAL number of columns in the matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatGetColumnSumsImaginaryPart()`, `VecSum()`, `MatGetColumnMeans()`, `MatGetColumnNorms()`, `MatGetColumnReductions()`

# External Links
$(_doc_external("Mat/MatGetColumnSumsRealPart"))
"""
function MatGetColumnSumsRealPart(petsclib::PetscLibType, A::AbstractPetscMat, sums::Vector{PetscReal}) end

@for_petsc function MatGetColumnSumsRealPart(petsclib::$UnionPetscLib, A::AbstractPetscMat, sums::Vector{$PetscReal} )

    @chk ccall(
               (:MatGetColumnSumsRealPart, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscReal}),
               A, sums,
              )


	return nothing
end 

"""
	MatGetColumnVector(petsclib::PetscLibType,A::AbstractPetscMat, yy::AbstractPetscVec, col::PetscInt) 
Gets the values from a given column of a matrix.

Not Collective

Input Parameters:
- `A`   - the matrix
- `yy`  - the vector
- `col` - the column requested (in global numbering)

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatGetRow()`, `MatGetDiagonal()`, `MatMult()`

# External Links
$(_doc_external("Mat/MatGetColumnVector"))
"""
function MatGetColumnVector(petsclib::PetscLibType, A::AbstractPetscMat, yy::AbstractPetscVec, col::PetscInt) end

@for_petsc function MatGetColumnVector(petsclib::$UnionPetscLib, A::AbstractPetscMat, yy::AbstractPetscVec, col::$PetscInt )

    @chk ccall(
               (:MatGetColumnVector, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, $PetscInt),
               A, yy, col,
              )


	return nothing
end 

"""
	m::PetscMemType = MatGetCurrentMemType(petsclib::PetscLibType,A::AbstractPetscMat) 
Get the memory location of the matrix

Not Collective, but the result will be the same on all MPI processes

Input Parameter:
- `A` - the matrix whose memory type we are checking

Output Parameter:
- `m` - the memory type

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatBoundToCPU()`, `PetscMemType`

# External Links
$(_doc_external("Mat/MatGetCurrentMemType"))
"""
function MatGetCurrentMemType(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatGetCurrentMemType(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	m_ = Ref{PetscMemType}()

    @chk ccall(
               (:MatGetCurrentMemType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscMemType}),
               A, m_,
              )

	m = m_[]

	return m
end 

"""
	dm::PetscDM = MatGetDM(petsclib::PetscLibType,A::AbstractPetscMat) 
Gets the `DM` defining the data layout of the matrix

Not Collective

Input Parameter:
- `A` - The `Mat`

Output Parameter:
- `dm` - The `DM`

Level: intermediate

Note:
A matrix may not have a `DM` associated with it

Developer Note:
Since the `Mat` class doesn't know about the `DM` class the `DM` object is associated with the `Mat` through a `PetscObjectCompose()` operation

See also: 
=== 
`DM`, `MatSetDM()`, `DMCreateMatrix()`, `DMSetMatType()`

# External Links
$(_doc_external("DM/MatGetDM"))
"""
function MatGetDM(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatGetDM(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	dm_ = Ref{CDM}()

    @chk ccall(
               (:MatGetDM, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CDM}),
               A, dm_,
              )

	dm = PetscDM(dm_[], petsclib)

	return dm
end 

"""
	MatGetDiagonal(petsclib::PetscLibType,mat::AbstractPetscMat, v::AbstractPetscVec) 
Gets the diagonal of a matrix as a `Vec`

Logically Collective

Input Parameter:
- `mat` - the matrix

Output Parameter:
- `v` - the diagonal of the matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `Vec`, `MatGetRow()`, `MatCreateSubMatrices()`, `MatCreateSubMatrix()`, `MatGetRowMaxAbs()`

# External Links
$(_doc_external("Mat/MatGetDiagonal"))
"""
function MatGetDiagonal(petsclib::PetscLibType, mat::AbstractPetscMat, v::AbstractPetscVec) end

@for_petsc function MatGetDiagonal(petsclib::$UnionPetscLib, mat::AbstractPetscMat, v::AbstractPetscVec )

    @chk ccall(
               (:MatGetDiagonal, $petsc_library),
               PetscErrorCode,
               (CMat, CVec),
               mat, v,
              )


	return nothing
end 

"""
	a::PetscMat = MatGetDiagonalBlock(petsclib::PetscLibType,A::AbstractPetscMat) 
Returns the part of the matrix associated with the on

Not Collective

Input Parameter:
- `A` - the matrix

Output Parameter:
- `a` - the diagonal part (which is a SEQUENTIAL matrix)

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatCreateAIJ()`, `MATAIJ`, `MATBAIJ`, `MATSBAIJ`

# External Links
$(_doc_external("Mat/MatGetDiagonalBlock"))
"""
function MatGetDiagonalBlock(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatGetDiagonalBlock(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	a_ = Ref{CMat}()

    @chk ccall(
               (:MatGetDiagonalBlock, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, a_,
              )

	a = PetscMat(a_[], petsclib)

	return a
end 

"""
	f::PetscMat = MatGetFactor(petsclib::PetscLibType,mat::AbstractPetscMat, type::MatSolverType, ftype::MatFactorType) 
Returns a matrix suitable to calls to MatXXFactorSymbolic,Numeric()

Collective

Input Parameters:
- `mat`   - the matrix
- `type`  - name of solver type, for example, `superlu`, `petsc` (to use PETSc's solver if it is available), if this is 'NULL', then the first result that satisfies
the other criteria is returned
- `ftype` - factor type, `MAT_FACTOR_LU`, `MAT_FACTOR_CHOLESKY`, `MAT_FACTOR_ICC`, `MAT_FACTOR_ILU`, `MAT_FACTOR_QR`

Output Parameter:
- `f` - the factor matrix used with MatXXFactorSymbolic,Numeric() calls. Can be `NULL` in some cases, see notes below.

Options Database Keys:
- `-pc_factor_mat_solver_type <type>`    - choose the type at run time. When using `KSP` solvers
- `-pc_factor_mat_factor_on_host <bool>` - do mat factorization on host (with device matrices). Default is doing it on device
- `-pc_factor_mat_solve_on_host <bool>`  - do mat solve on host (with device matrices). Default is doing it on device

Level: intermediate

-seealso: [](ch_matrices), `Mat`, [Matrix Factorization](sec_matfactor), `KSP`, `MatSolverType`, `MatFactorType`, `MatCopy()`, `MatDuplicate()`,
`MatGetFactorAvailable()`, `MatFactorGetCanUseOrdering()`, `MatSolverTypeRegister()`, `MatSolverTypeGet()`
`MAT_FACTOR_LU`, `MAT_FACTOR_CHOLESKY`, `MAT_FACTOR_ICC`, `MAT_FACTOR_ILU`, `MAT_FACTOR_QR`, `MatInitializePackage()`

# External Links
$(_doc_external("Mat/MatGetFactor"))
"""
function MatGetFactor(petsclib::PetscLibType, mat::AbstractPetscMat, type::MatSolverType, ftype::MatFactorType) end

@for_petsc function MatGetFactor(petsclib::$UnionPetscLib, mat::AbstractPetscMat, type::MatSolverType, ftype::MatFactorType )
	f_ = Ref{CMat}()

    @chk ccall(
               (:MatGetFactor, $petsc_library),
               PetscErrorCode,
               (CMat, MatSolverType, MatFactorType, Ptr{CMat}),
               mat, type, ftype, f_,
              )

	f = PetscMat(f_[], petsclib)

	return f
end 

"""
	flg::PetscBool = MatGetFactorAvailable(petsclib::PetscLibType,mat::AbstractPetscMat, type::MatSolverType, ftype::MatFactorType) 
Returns a flag if matrix supports particular type and factor type

Not Collective

Input Parameters:
- `mat`   - the matrix
- `type`  - name of solver type, for example, `superlu`, `petsc` (to use PETSc's default)
- `ftype` - factor type, `MAT_FACTOR_LU`, `MAT_FACTOR_CHOLESKY`, `MAT_FACTOR_ICC`, `MAT_FACTOR_ILU`, `MAT_FACTOR_QR`

Output Parameter:
- `flg` - PETSC_TRUE if the factorization is available

Level: intermediate

-seealso: [](ch_matrices), `Mat`, [Matrix Factorization](sec_matfactor), `MatSolverType`, `MatFactorType`, `MatGetFactor()`, `MatCopy()`, `MatDuplicate()`, `MatSolverTypeRegister()`,
`MAT_FACTOR_LU`, `MAT_FACTOR_CHOLESKY`, `MAT_FACTOR_ICC`, `MAT_FACTOR_ILU`, `MAT_FACTOR_QR`, `MatSolverTypeGet()`

# External Links
$(_doc_external("Mat/MatGetFactorAvailable"))
"""
function MatGetFactorAvailable(petsclib::PetscLibType, mat::AbstractPetscMat, type::MatSolverType, ftype::MatFactorType) end

@for_petsc function MatGetFactorAvailable(petsclib::$UnionPetscLib, mat::AbstractPetscMat, type::MatSolverType, ftype::MatFactorType )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatGetFactorAvailable, $petsc_library),
               PetscErrorCode,
               (CMat, MatSolverType, MatFactorType, Ptr{PetscBool}),
               mat, type, ftype, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	t::MatFactorType = MatGetFactorType(petsclib::PetscLibType,mat::AbstractPetscMat) 
gets the type of factorization a matrix is

Not Collective

Input Parameter:
- `mat` - the matrix

Output Parameter:
- `t` - the type, one of `MAT_FACTOR_NONE`, `MAT_FACTOR_LU`, `MAT_FACTOR_CHOLESKY`, `MAT_FACTOR_ILU`, `MAT_FACTOR_ICC,MAT_FACTOR_ILUDT`, `MAT_FACTOR_QR`

Level: intermediate

-seealso: [](ch_matrices), `Mat`, [Matrix Factorization](sec_matfactor), `MatFactorType`, `MatGetFactor()`, `MatSetFactorType()`, `MAT_FACTOR_NONE`, `MAT_FACTOR_LU`, `MAT_FACTOR_CHOLESKY`, `MAT_FACTOR_ILU`,
`MAT_FACTOR_ICC`,`MAT_FACTOR_ILUDT`, `MAT_FACTOR_QR`

# External Links
$(_doc_external("Mat/MatGetFactorType"))
"""
function MatGetFactorType(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatGetFactorType(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	t_ = Ref{MatFactorType}()

    @chk ccall(
               (:MatGetFactorType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatFactorType}),
               mat, t_,
              )

	t = t_[]

	return t
end 

"""
	nghosts::PetscInt,ghosts::Vector{PetscInt} = MatGetGhosts(petsclib::PetscLibType,mat::AbstractPetscMat) 
Get the global indices of all ghost nodes defined by the sparse matrix

Collective

Input Parameter:
- `mat` - the matrix

Output Parameters:
- `nghosts` - number of ghosts (for `MATBAIJ` and `MATSBAIJ` matrices there is one ghost for each matrix block)
- `ghosts`  - the global indices of the ghost points

Level: advanced

-seealso: [](ch_matrices), `Mat`, `VecCreateGhost()`, `VecCreateGhostBlock()`

# External Links
$(_doc_external("Mat/MatGetGhosts"))
"""
function MatGetGhosts(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatGetGhosts(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	nghosts_ = Ref{$PetscInt}()
	ghosts_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:MatGetGhosts, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
               mat, nghosts_, ghosts_,
              )

	nghosts = nghosts_[]
	nghosts = nghosts_[]
	ghosts = unsafe_wrap(Array, ghosts_[], Int(nghosts); own = false)

	return nghosts,ghosts
end 

"""
	nneg::PetscInt,nzero::PetscInt,npos::PetscInt = MatGetInertia(petsclib::PetscLibType,mat::AbstractPetscMat) 
Gets the inertia from a factored matrix

Collective

Input Parameter:
- `mat` - the matrix

Output Parameters:
- `nneg`  - number of negative eigenvalues
- `nzero` - number of zero eigenvalues
- `npos`  - number of positive eigenvalues

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatGetFactor()`, `MatCholeskyFactor()`

# External Links
$(_doc_external("Mat/MatGetInertia"))
"""
function MatGetInertia(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatGetInertia(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	nneg_ = Ref{$PetscInt}()
	nzero_ = Ref{$PetscInt}()
	npos_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatGetInertia, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
               mat, nneg_, nzero_, npos_,
              )

	nneg = nneg_[]
	nzero = nzero_[]
	npos = npos_[]

	return nneg,nzero,npos
end 

"""
	info::MatInfo = MatGetInfo(petsclib::PetscLibType,mat::AbstractPetscMat, flag::MatInfoType) 
Returns information about matrix storage (number of
nonzeros, memory, etc.).

Collective if `MAT_GLOBAL_MAX` or `MAT_GLOBAL_SUM` is used as the flag

Input Parameters:
- `mat`  - the matrix
- `flag` - flag indicating the type of parameters to be returned (`MAT_LOCAL` - local matrix, `MAT_GLOBAL_MAX` - maximum over all processors, `MAT_GLOBAL_SUM` - sum over all processors)

Output Parameter:
- `info` - matrix information context

Options Database Key:
- `-mat_view ::ascii_info` - print matrix info to `PETSC_STDOUT`

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatInfo`, `MatStashGetInfo()`

# External Links
$(_doc_external("Mat/MatGetInfo"))
"""
function MatGetInfo(petsclib::PetscLibType, mat::AbstractPetscMat, flag::MatInfoType) end

@for_petsc function MatGetInfo(petsclib::$UnionPetscLib, mat::AbstractPetscMat, flag::MatInfoType )
	info_ = Ref{MatInfo}()

    @chk ccall(
               (:MatGetInfo, $petsc_library),
               PetscErrorCode,
               (CMat, MatInfoType, Ptr{MatInfo}),
               mat, flag, info_,
              )

	info = info_[]

	return info
end 

"""
	rmap::PetscLayout,cmap::PetscLayout = MatGetLayouts(petsclib::PetscLibType,A::AbstractPetscMat) 
Gets the `PetscLayout` objects for rows and columns

Not Collective

Input Parameter:
- `A` - the matrix

Output Parameters:
- `rmap` - row layout
- `cmap` - column layout

Level: advanced

-seealso: [](ch_matrices), `Mat`, [Matrix Layouts](sec_matlayout), `PetscLayout`, `MatCreateVecs()`, `MatGetLocalToGlobalMapping()`, `MatSetLayouts()`

# External Links
$(_doc_external("Mat/MatGetLayouts"))
"""
function MatGetLayouts(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatGetLayouts(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	rmap_ = Ref{PetscLayout}()
	cmap_ = Ref{PetscLayout}()

    @chk ccall(
               (:MatGetLayouts, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscLayout}, Ptr{PetscLayout}),
               A, rmap_, cmap_,
              )

	rmap = rmap_[]
	cmap = cmap_[]

	return rmap,cmap
end 

"""
	m::PetscInt,n::PetscInt = MatGetLocalSize(petsclib::PetscLibType,mat::AbstractPetscMat) 
For most matrix formats, excluding `MATELEMENTAL` and `MATSCALAPACK`, Returns the number of local rows and local columns
of a matrix. For all matrices this is the local size of the left and right vectors as returned by `MatCreateVecs()`.

Not Collective

Input Parameter:
- `mat` - the matrix

Output Parameters:
- `m` - the number of local rows, use `NULL` to not obtain this value
- `n` - the number of local columns, use `NULL` to not obtain this value

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatSetSizes()`, `MatGetSize()`

# External Links
$(_doc_external("Mat/MatGetLocalSize"))
"""
function MatGetLocalSize(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatGetLocalSize(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	m_ = Ref{$PetscInt}()
	n_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatGetLocalSize, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}, Ptr{$PetscInt}),
               mat, m_, n_,
              )

	m = m_[]
	n = n_[]

	return m,n
end 

"""
	submat::PetscMat = MatGetLocalSubMatrix(petsclib::PetscLibType,mat::AbstractPetscMat, isrow::AbstractIS, iscol::AbstractIS) 
Gets a reference to a submatrix specified in local numbering

Not Collective

Input Parameters:
- `mat`   - matrix to extract local submatrix from
- `isrow` - local row indices for submatrix
- `iscol` - local column indices for submatrix

Output Parameter:
- `submat` - the submatrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatRestoreLocalSubMatrix()`, `MatCreateLocalRef()`, `MatSetLocalToGlobalMapping()`

# External Links
$(_doc_external("Mat/MatGetLocalSubMatrix"))
"""
function MatGetLocalSubMatrix(petsclib::PetscLibType, mat::AbstractPetscMat, isrow::AbstractIS, iscol::AbstractIS) end

@for_petsc function MatGetLocalSubMatrix(petsclib::$UnionPetscLib, mat::AbstractPetscMat, isrow::AbstractIS, iscol::AbstractIS )
	submat_ = Ref{CMat}()

    @chk ccall(
               (:MatGetLocalSubMatrix, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, CIS, Ptr{CMat}),
               mat, isrow, iscol, submat_,
              )

	submat = PetscMat(submat_[], petsclib)

	return submat
end 

"""
	rmapping::ISLocalToGlobalMapping,cmapping::ISLocalToGlobalMapping = MatGetLocalToGlobalMapping(petsclib::PetscLibType,A::AbstractPetscMat) 
Gets the local

Not Collective

Input Parameter:
- `A` - the matrix

Output Parameters:
- `rmapping` - row mapping
- `cmapping` - column mapping

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatSetLocalToGlobalMapping()`, `MatSetValuesLocal()`

# External Links
$(_doc_external("Mat/MatGetLocalToGlobalMapping"))
"""
function MatGetLocalToGlobalMapping(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatGetLocalToGlobalMapping(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	rmapping_ = Ref{ISLocalToGlobalMapping}()
	cmapping_ = Ref{ISLocalToGlobalMapping}()

    @chk ccall(
               (:MatGetLocalToGlobalMapping, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{ISLocalToGlobalMapping}, Ptr{ISLocalToGlobalMapping}),
               A, rmapping_, cmapping_,
              )

	rmapping = rmapping_[]
	cmapping = cmapping_[]

	return rmapping,cmapping
end 

"""
	subMat::PetscMat = MatGetMultiProcBlock(petsclib::PetscLibType,mat::AbstractPetscMat, subComm::MPI_Comm, scall::MatReuse) 
Create multiple 'parallel submatrices' from
a given `Mat`. Each submatrix can span multiple procs.

Collective

Input Parameters:
- `mat`     - the matrix
- `subComm` - the sub communicator obtained as if by `MPI_Comm_split(PetscObjectComm((PetscObject)mat))`
- `scall`   - either `MAT_INITIAL_MATRIX` or `MAT_REUSE_MATRIX`

Output Parameter:
- `subMat` - parallel sub-matrices each spanning a given `subcomm`

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatCreateRedundantMatrix()`, `MatCreateSubMatrices()`, `PCBJACOBI`

# External Links
$(_doc_external("Mat/MatGetMultiProcBlock"))
"""
function MatGetMultiProcBlock(petsclib::PetscLibType, mat::AbstractPetscMat, subComm::MPI_Comm, scall::MatReuse) end

@for_petsc function MatGetMultiProcBlock(petsclib::$UnionPetscLib, mat::AbstractPetscMat, subComm::MPI_Comm, scall::MatReuse )
	subMat_ = Ref{CMat}()

    @chk ccall(
               (:MatGetMultiProcBlock, $petsc_library),
               PetscErrorCode,
               (CMat, MPI_Comm, MatReuse, Ptr{CMat}),
               mat, subComm, scall, subMat_,
              )

	subMat = PetscMat(subMat_[], petsclib)

	return subMat
end 

"""
	nullsp::MatNullSpace = MatGetNearNullSpace(petsclib::PetscLibType,mat::AbstractPetscMat) 
Get null space attached with `MatSetNearNullSpace()`

Not Collective

Input Parameter:
- `mat` - the matrix

Output Parameter:
- `nullsp` - the null space object, `NULL` if not set

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatNullSpace`, `MatSetNearNullSpace()`, `MatGetNullSpace()`, `MatNullSpaceCreate()`

# External Links
$(_doc_external("Mat/MatGetNearNullSpace"))
"""
function MatGetNearNullSpace(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatGetNearNullSpace(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	nullsp_ = Ref{MatNullSpace}()

    @chk ccall(
               (:MatGetNearNullSpace, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatNullSpace}),
               mat, nullsp_,
              )

	nullsp = nullsp_[]

	return nullsp
end 

"""
	state::PetscObjectState = MatGetNonzeroState(petsclib::PetscLibType,mat::AbstractPetscMat) 
Returns a 64
matrix has had new nonzero locations added to (or removed from) the matrix since the previous call, the value will be larger.

Not Collective

Input Parameter:
- `mat` - the matrix

Output Parameter:
- `state` - the current state

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `PetscObjectStateGet()`, `PetscObjectGetId()`

# External Links
$(_doc_external("Mat/MatGetNonzeroState"))
"""
function MatGetNonzeroState(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatGetNonzeroState(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	state_ = Ref{PetscObjectState}()

    @chk ccall(
               (:MatGetNonzeroState, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscObjectState}),
               mat, state_,
              )

	state = state_[]

	return state
end 

"""
	MatGetNullSpace(petsclib::PetscLibType,mat::AbstractPetscMat, nullsp::MatNullSpace) 
retrieves the null space of a matrix.

Logically Collective

Input Parameters:
- `mat`    - the matrix
- `nullsp` - the null space object

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `MatNullSpaceCreate()`, `MatSetNearNullSpace()`, `MatSetNullSpace()`, `MatNullSpace`

# External Links
$(_doc_external("Mat/MatGetNullSpace"))
"""
function MatGetNullSpace(petsclib::PetscLibType, mat::AbstractPetscMat, nullsp::MatNullSpace) end

@for_petsc function MatGetNullSpace(petsclib::$UnionPetscLib, mat::AbstractPetscMat, nullsp::MatNullSpace )

    @chk ccall(
               (:MatGetNullSpace, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatNullSpace}),
               mat, nullsp,
              )


	return nothing
end 

"""
	nullsp::Ptr{MatNullSpace} = MatGetNullSpaces(petsclib::PetscLibType,n::PetscInt, mat::Vector{<:AbstractPetscMat}) 
gets the null spaces, transpose null spaces, and near null spaces from an array of matrices

Logically Collective

Input Parameters:
- `n`   - the number of matrices
- `mat` - the array of matrices

Output Parameters:
- `nullsp` - an array of null spaces, `NULL` for each matrix that does not have a null space, length 3 * `n`

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `MatNullSpaceCreate()`, `MatSetNearNullSpace()`, `MatGetNullSpace()`, `MatSetTransposeNullSpace()`, `MatGetTransposeNullSpace()`,
`MatNullSpaceRemove()`, `MatRestoreNullSpaces()`

# External Links
$(_doc_external("Mat/MatGetNullSpaces"))
"""
function MatGetNullSpaces(petsclib::PetscLibType, n::PetscInt, mat::Vector{<:AbstractPetscMat}) end

@for_petsc function MatGetNullSpaces(petsclib::$UnionPetscLib, n::$PetscInt, mat::Vector{<:AbstractPetscMat} )
	nullsp_ = Ref{Ptr{MatNullSpace}}()

    @chk ccall(
               (:MatGetNullSpaces, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{CMat}, Ptr{Ptr{MatNullSpace}}),
               n, mat, nullsp_,
              )

	nullsp = nullsp_[]

	return nullsp
end 

"""
	flg::PetscBool = MatGetOption(petsclib::PetscLibType,mat::AbstractPetscMat, op::MatOption) 
Gets a parameter option that has been set for a matrix.

Logically Collective

Input Parameters:
- `mat` - the matrix
- `op`  - the option, this only responds to certain options, check the code for which ones

Output Parameter:
- `flg` - turn the option on (`PETSC_TRUE`) or off (`PETSC_FALSE`)

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatOption`, `MatSetOption()`, `MatIsSymmetric()`, `MatIsHermitian()`, `MatIsStructurallySymmetric()`,
`MatIsSymmetricKnown()`, `MatIsHermitianKnown()`, `MatIsStructurallySymmetricKnown()`

# External Links
$(_doc_external("Mat/MatGetOption"))
"""
function MatGetOption(petsclib::PetscLibType, mat::AbstractPetscMat, op::MatOption) end

@for_petsc function MatGetOption(petsclib::$UnionPetscLib, mat::AbstractPetscMat, op::MatOption )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatGetOption, $petsc_library),
               PetscErrorCode,
               (CMat, MatOption, Ptr{PetscBool}),
               mat, op, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	prefix::Ptr{Cchar} = MatGetOptionsPrefix(petsclib::PetscLibType,A::AbstractPetscMat) 
Gets the prefix used for searching for all
matrix options in the database.

Not Collective

Input Parameter:
- `A` - the matrix

Output Parameter:
- `prefix` - pointer to the prefix string used

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatAppendOptionsPrefix()`, `MatSetOptionsPrefix()`, `MatAppendOptionsPrefixFactor()`, `MatSetOptionsPrefixFactor()`

# External Links
$(_doc_external("Mat/MatGetOptionsPrefix"))
"""
function MatGetOptionsPrefix(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatGetOptionsPrefix(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	prefix_ = Ref{Ptr{Cchar}}()

    @chk ccall(
               (:MatGetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{Cchar}}),
               A, prefix_,
              )

	prefix = prefix_[]

	return prefix
end 

"""
	rperm::IS,cperm::IS = MatGetOrdering(petsclib::PetscLibType,mat::AbstractPetscMat, type::MatOrderingType) 
Gets a reordering for a matrix to reduce fill or to
improve numerical stability of LU factorization.

Collective

Input Parameters:
- `mat`  - the matrix
- `type` - type of reordering, one of the following
-seealso: `MatOrderingRegister()`, `PCFactorSetMatOrderingType()`, `MatColoring`, `MatColoringCreate()`, `MatOrderingType`, `Mat`

# External Links
$(_doc_external("MatGraphOperations/MatGetOrdering"))
"""
function MatGetOrdering(petsclib::PetscLibType, mat::AbstractPetscMat, type::MatOrderingType) end

@for_petsc function MatGetOrdering(petsclib::$UnionPetscLib, mat::AbstractPetscMat, type::MatOrderingType )
	rperm_ = Ref{CIS}()
	cperm_ = Ref{CIS}()

    @chk ccall(
               (:MatGetOrdering, $petsc_library),
               PetscErrorCode,
               (CMat, MatOrderingType, Ptr{CIS}, Ptr{CIS}),
               mat, type, rperm_, cperm_,
              )

	rperm = IS(rperm_[], petsclib)
	cperm = IS(cperm_[], petsclib)

	return rperm,cperm
end 

"""
	MatGetOrderingList(petsclib::PetscLibType,list::PetscFunctionList) 

# External Links
$(_doc_external("MatGraphOperations/MatGetOrderingList"))
"""
function MatGetOrderingList(petsclib::PetscLibType, list::PetscFunctionList) end

@for_petsc function MatGetOrderingList(petsclib::$UnionPetscLib, list::PetscFunctionList )

    @chk ccall(
               (:MatGetOrderingList, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscFunctionList},),
               list,
              )


	return nothing
end 

"""
	rows::IS,cols::IS = MatGetOwnershipIS(petsclib::PetscLibType,A::AbstractPetscMat) 
Get row and column ownership of a matrices' values as index sets.

Not Collective

Input Parameter:
- `A` - matrix

Output Parameters:
- `rows` - rows in which this process owns elements, , use `NULL` to not obtain this value
- `cols` - columns in which this process owns elements, use `NULL` to not obtain this value

Level: intermediate

-seealso: [](ch_matrices), `IS`, `Mat`, `MatGetOwnershipRanges()`, `MatSetValues()`, `MATELEMENTAL`, `MATSCALAPACK`

# External Links
$(_doc_external("Mat/MatGetOwnershipIS"))
"""
function MatGetOwnershipIS(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatGetOwnershipIS(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	rows_ = Ref{CIS}()
	cols_ = Ref{CIS}()

    @chk ccall(
               (:MatGetOwnershipIS, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CIS}, Ptr{CIS}),
               A, rows_, cols_,
              )

	rows = IS(rows_[], petsclib)
	cols = IS(cols_[], petsclib)

	return rows,cols
end 

"""
	m::PetscInt,n::PetscInt = MatGetOwnershipRange(petsclib::PetscLibType,mat::AbstractPetscMat) 
For matrices that own values by row, excludes `MATELEMENTAL` and `MATSCALAPACK`, returns the range of matrix rows owned by
this MPI process.

Not Collective

Input Parameter:
- `mat` - the matrix

Output Parameters:
- `m` - the global index of the first local row, use `NULL` to not obtain this value
- `n` - one more than the global index of the last local row, use `NULL` to not obtain this value

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatGetOwnershipRanges()`, `MatGetOwnershipRangeColumn()`, `MatGetOwnershipRangesColumn()`, `PetscSplitOwnership()`,
`PetscSplitOwnershipBlock()`, `PetscLayout`, `MatSetSizes()`, `MatCreateAIJ()`, `DMDAGetGhostCorners()`, `DM`

# External Links
$(_doc_external("Mat/MatGetOwnershipRange"))
"""
function MatGetOwnershipRange(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatGetOwnershipRange(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	m_ = Ref{$PetscInt}()
	n_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatGetOwnershipRange, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}, Ptr{$PetscInt}),
               mat, m_, n_,
              )

	m = m_[]
	n = n_[]

	return m,n
end 

"""
	m::PetscInt,n::PetscInt = MatGetOwnershipRangeColumn(petsclib::PetscLibType,mat::AbstractPetscMat) 
Returns the range of matrix columns associated with rows of a
vector one multiplies this matrix by that are owned by this processor.

Not Collective, unless matrix has not been allocated, then collective

Input Parameter:
- `mat` - the matrix

Output Parameters:
- `m` - the global index of the first local column, use `NULL` to not obtain this value
- `n` - one more than the global index of the last local column, use `NULL` to not obtain this value

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatGetOwnershipRange()`, `MatGetOwnershipRanges()`, `MatGetOwnershipRangesColumn()`, `PetscLayout`,
`MatSetSizes()`, `MatCreateAIJ()`, `DMDAGetGhostCorners()`, `DM`

# External Links
$(_doc_external("Mat/MatGetOwnershipRangeColumn"))
"""
function MatGetOwnershipRangeColumn(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatGetOwnershipRangeColumn(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	m_ = Ref{$PetscInt}()
	n_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatGetOwnershipRangeColumn, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}, Ptr{$PetscInt}),
               mat, m_, n_,
              )

	m = m_[]
	n = n_[]

	return m,n
end 

"""
	ranges::Vector{PetscInt} = MatGetOwnershipRanges(petsclib::PetscLibType,mat::AbstractPetscMat) 
For matrices that own values by row, excludes `MATELEMENTAL` and
`MATSCALAPACK`, returns the range of matrix rows owned by each process.

Not Collective, unless matrix has not been allocated

Input Parameter:
- `mat` - the matrix

Output Parameter:
- `ranges` - start of each processors portion plus one more than the total length at the end, of length `size` + 1
where `size` is the number of MPI processes used by `mat`

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatGetOwnershipRange()`, `MatGetOwnershipRangeColumn()`, `MatGetOwnershipRangesColumn()`, `PetscLayout`,
`PetscSplitOwnership()`, `PetscSplitOwnershipBlock()`, `MatSetSizes()`, `MatCreateAIJ()`,
`DMDAGetGhostCorners()`, `DM`

# External Links
$(_doc_external("Mat/MatGetOwnershipRanges"))
"""
function MatGetOwnershipRanges(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatGetOwnershipRanges(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	ranges_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:MatGetOwnershipRanges, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscInt}}),
               mat, ranges_,
              )

	comm_ref = Ref{MPI.MPI_Comm}()
	comm = MPI.Comm(comm_ref[])
	nproc = MPI.Comm_size(comm)
	ranges = unsafe_wrap(Array, ranges_[], nproc + 1; own = false)

	return ranges
end 

"""
	ranges::Vector{PetscInt} = MatGetOwnershipRangesColumn(petsclib::PetscLibType,mat::AbstractPetscMat) 
Returns the ranges of matrix columns associated with rows of a
vector one multiplies this vector by that are owned by each processor.

Not Collective, unless matrix has not been allocated

Input Parameter:
- `mat` - the matrix

Output Parameter:
- `ranges` - start of each processors portion plus one more than the total length at the end

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatGetOwnershipRange()`, `MatGetOwnershipRangeColumn()`, `MatGetOwnershipRanges()`,
`PetscSplitOwnership()`, `PetscSplitOwnershipBlock()`, `PetscLayout`, `MatSetSizes()`, `MatCreateAIJ()`,
`DMDAGetGhostCorners()`, `DM`

# External Links
$(_doc_external("Mat/MatGetOwnershipRangesColumn"))
"""
function MatGetOwnershipRangesColumn(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatGetOwnershipRangesColumn(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	ranges_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:MatGetOwnershipRangesColumn, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscInt}}),
               mat, ranges_,
              )

	comm_ref = Ref{MPI.MPI_Comm}()
	comm = MPI.Comm(comm_ref[])
	nproc = MPI.Comm_size(comm)
	ranges = unsafe_wrap(Array, ranges_[], nproc + 1; own = false)

	return ranges
end 

"""
	ncols::PetscInt,cols::Vector{PetscInt},vals::Vector{PetscScalar} = MatGetRow(petsclib::PetscLibType,mat::AbstractPetscMat, row::PetscInt) 
Gets a row of a matrix.  You MUST call `MatRestoreRow()`
for each row that you get to ensure that your application does
not bleed memory.

Not Collective

Input Parameters:
- `mat` - the matrix
- `row` - the row to get

Output Parameters:
- `ncols` - if not `NULL`, the number of nonzeros in `row`
- `cols`  - if not `NULL`, the column numbers
- `vals`  - if not `NULL`, the numerical values

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatRestoreRow()`, `MatSetValues()`, `MatGetValues()`, `MatCreateSubMatrices()`, `MatGetDiagonal()`, `MatGetRowIJ()`, `MatRestoreRowIJ()`

# External Links
$(_doc_external("Mat/MatGetRow"))
"""
function MatGetRow(petsclib::PetscLibType, mat::AbstractPetscMat, row::PetscInt) end

@for_petsc function MatGetRow(petsclib::$UnionPetscLib, mat::AbstractPetscMat, row::$PetscInt )
	ncols_ = Ref{$PetscInt}()
	cols_ = Ref{Ptr{$PetscInt}}()
	vals_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatGetRow, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscScalar}}),
               mat, row, ncols_, cols_, vals_,
              )

	ncols = ncols_[]
	ncols = ncols_[]
	cols = unsafe_wrap(Array, cols_[], Int(ncols); own = false)
	ncols = ncols_[]
	vals = unsafe_wrap(Array, vals_[], Int(ncols); own = false)

	return ncols,cols,vals
end 

# override for MatGetRowIJ; C signature: MatGetRowIJ(Mat mat, PetscInt shift, PetscBool symmetric, PetscBool inodecompressed, PetscInt* n, PetscInt* ia[], PetscInt* ja[], PetscBool* done)
"""
	n::PetscInt,ia::Vector{PetscInt},ja::Vector{PetscInt},done::PetscBool = MatGetRowIJ(petsclib::PetscLibType,mat::AbstractPetscMat, shift::PetscInt, symmetric::PetscBool, inodecompressed::PetscBool) 
Returns the compressed row storage i and j indices for the local rows of a sparse matrix

Collective

Input Parameters:
- `mat`             - the matrix
- `shift`           - 0 or 1 indicating we want the indices starting at 0 or 1
- `symmetric`       - `PETSC_TRUE` or `PETSC_FALSE` indicating the matrix data structure should be symmetrized
- `inodecompressed` - `PETSC_TRUE` or `PETSC_FALSE`  indicating if the nonzero structure of the
inodes or the nonzero elements is wanted. For `MATBAIJ` matrices the compressed version is
always used.

Output Parameters:
- `n`    - number of local rows in the (possibly compressed) matrix, use `NULL` if not needed
- `ia`   - the row pointers; that is ia[0] = 0, ia[row] = ia[row-1] + number of elements in that row of the matrix, use `NULL` if not needed
- `ja`   - the column indices, use `NULL` if not needed
- `done` - indicates if the routine actually worked and returned appropriate ia[] and ja[] arrays; callers
are responsible for handling the case when done == `PETSC_FALSE` and ia and ja are not set

Level: developer

-seealso: [](ch_matrices), `Mat`, `MATAIJ`, `MatGetColumnIJ()`, `MatRestoreRowIJ()`, `MatSeqAIJGetArray()`

# External Links
$(_doc_external("Mat/MatGetRowIJ"))
"""
function MatGetRowIJ(petsclib::PetscLibType, mat::AbstractPetscMat, shift::PetscInt, symmetric::PetscBool, inodecompressed::PetscBool) end

@for_petsc function MatGetRowIJ(petsclib::$UnionPetscLib, mat::AbstractPetscMat, shift::$PetscInt, symmetric::PetscBool, inodecompressed::PetscBool )
	n_ = Ref{$PetscInt}()
	ia_ = Ref{Ptr{$PetscInt}}(C_NULL)
	ja_ = Ref{Ptr{$PetscInt}}(C_NULL)
	done_ = Ref{PetscBool}()

    @chk ccall(
               (:MatGetRowIJ, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, PetscBool, PetscBool, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}, Ptr{PetscBool}),
               mat, shift, symmetric, inodecompressed, n_, ia_, ja_, done_,
              )

	n = n_[]
	done = done_[]
	# `ia` has n+1 entries (ia[1] == shift), `ja` has ia[n+1] - shift entries; both are
	# PETSc-owned and must be given back with MatRestoreRowIJ. When PETSc cannot
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

"""
	MatGetRowMax(petsclib::PetscLibType,mat::AbstractPetscMat, v::AbstractPetscVec, idx::Vector{PetscInt}) 
Gets the maximum value (of the real part) of each
row of the matrix

Logically Collective

Input Parameter:
- `mat` - the matrix

Output Parameters:
- `v`   - the vector for storing the maximums
- `idx` - the indices of the column found for each row (optional, otherwise pass `NULL`)

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatGetDiagonal()`, `MatCreateSubMatrices()`, `MatCreateSubMatrix()`, `MatGetRowMaxAbs()`, `MatGetRowMin()`, `MatGetRowMinAbs()`

# External Links
$(_doc_external("Mat/MatGetRowMax"))
"""
function MatGetRowMax(petsclib::PetscLibType, mat::AbstractPetscMat, v::AbstractPetscVec, idx::Vector{PetscInt}) end

@for_petsc function MatGetRowMax(petsclib::$UnionPetscLib, mat::AbstractPetscMat, v::AbstractPetscVec, idx::Vector{$PetscInt} )

    @chk ccall(
               (:MatGetRowMax, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, Ptr{$PetscInt}),
               mat, v, idx,
              )


	return nothing
end 

"""
	MatGetRowMaxAbs(petsclib::PetscLibType,mat::AbstractPetscMat, v::AbstractPetscVec, idx::Vector{PetscInt}) 
Gets the maximum value (in absolute value) of each
row of the matrix

Logically Collective

Input Parameter:
- `mat` - the matrix

Output Parameters:
- `v`   - the vector for storing the maximums
- `idx` - the indices of the column found for each row (or `NULL` if not needed)

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatGetDiagonal()`, `MatCreateSubMatrices()`, `MatCreateSubMatrix()`, `MatGetRowSum()`, `MatGetRowMin()`, `MatGetRowMinAbs()`

# External Links
$(_doc_external("Mat/MatGetRowMaxAbs"))
"""
function MatGetRowMaxAbs(petsclib::PetscLibType, mat::AbstractPetscMat, v::AbstractPetscVec, idx::Vector{PetscInt}) end

@for_petsc function MatGetRowMaxAbs(petsclib::$UnionPetscLib, mat::AbstractPetscMat, v::AbstractPetscVec, idx::Vector{$PetscInt} )

    @chk ccall(
               (:MatGetRowMaxAbs, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, Ptr{$PetscInt}),
               mat, v, idx,
              )


	return nothing
end 

"""
	MatGetRowMin(petsclib::PetscLibType,mat::AbstractPetscMat, v::AbstractPetscVec, idx::Vector{PetscInt}) 
Gets the minimum value (of the real part) of each
row of the matrix

Logically Collective

Input Parameter:
- `mat` - the matrix

Output Parameters:
- `v`   - the vector for storing the maximums
- `idx` - the indices of the column found for each row (optional, pass `NULL` if not needed)

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatGetDiagonal()`, `MatCreateSubMatrices()`, `MatCreateSubMatrix()`, `MatGetRowMaxAbs()`, `MatGetRowMinAbs()`,
`MatGetRowMax()`

# External Links
$(_doc_external("Mat/MatGetRowMin"))
"""
function MatGetRowMin(petsclib::PetscLibType, mat::AbstractPetscMat, v::AbstractPetscVec, idx::Vector{PetscInt}) end

@for_petsc function MatGetRowMin(petsclib::$UnionPetscLib, mat::AbstractPetscMat, v::AbstractPetscVec, idx::Vector{$PetscInt} )

    @chk ccall(
               (:MatGetRowMin, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, Ptr{$PetscInt}),
               mat, v, idx,
              )


	return nothing
end 

"""
	MatGetRowMinAbs(petsclib::PetscLibType,mat::AbstractPetscMat, v::AbstractPetscVec, idx::Vector{PetscInt}) 
Gets the minimum value (in absolute value) of each
row of the matrix

Logically Collective

Input Parameter:
- `mat` - the matrix

Output Parameters:
- `v`   - the vector for storing the minimums
- `idx` - the indices of the column found for each row (or `NULL` if not needed)

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatGetDiagonal()`, `MatCreateSubMatrices()`, `MatCreateSubMatrix()`, `MatGetRowMax()`, `MatGetRowMaxAbs()`, `MatGetRowMin()`

# External Links
$(_doc_external("Mat/MatGetRowMinAbs"))
"""
function MatGetRowMinAbs(petsclib::PetscLibType, mat::AbstractPetscMat, v::AbstractPetscVec, idx::Vector{PetscInt}) end

@for_petsc function MatGetRowMinAbs(petsclib::$UnionPetscLib, mat::AbstractPetscMat, v::AbstractPetscVec, idx::Vector{$PetscInt} )

    @chk ccall(
               (:MatGetRowMinAbs, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, Ptr{$PetscInt}),
               mat, v, idx,
              )


	return nothing
end 

"""
	MatGetRowSum(petsclib::PetscLibType,mat::AbstractPetscMat, v::AbstractPetscVec) 
Gets the sum of each row of the matrix

Logically or Neighborhood Collective

Input Parameter:
- `mat` - the matrix

Output Parameter:
- `v` - the vector for storing the sum of rows

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatGetDiagonal()`, `MatCreateSubMatrices()`, `MatCreateSubMatrix()`, `MatGetRowMax()`, `MatGetRowMin()`, `MatGetRowMaxAbs()`, `MatGetRowMinAbs()`, `MatGetRowSumAbs()`

# External Links
$(_doc_external("Mat/MatGetRowSum"))
"""
function MatGetRowSum(petsclib::PetscLibType, mat::AbstractPetscMat, v::AbstractPetscVec) end

@for_petsc function MatGetRowSum(petsclib::$UnionPetscLib, mat::AbstractPetscMat, v::AbstractPetscVec )

    @chk ccall(
               (:MatGetRowSum, $petsc_library),
               PetscErrorCode,
               (CMat, CVec),
               mat, v,
              )


	return nothing
end 

"""
	MatGetRowSumAbs(petsclib::PetscLibType,mat::AbstractPetscMat, v::AbstractPetscVec) 
Gets the sum value (in absolute value) of each row of the matrix

Logically Collective

Input Parameter:
- `mat` - the matrix

Output Parameter:
- `v` - the vector for storing the sum

Level: intermediate

This code is only implemented for a couple of matrix formats.

-seealso: [](ch_matrices), `Mat`, `MatGetDiagonal()`, `MatCreateSubMatrices()`, `MatCreateSubMatrix()`, `MatGetRowMax()`, `MatGetRowMin()`, `MatGetRowMinAbs()`

# External Links
$(_doc_external("Mat/MatGetRowSumAbs"))
"""
function MatGetRowSumAbs(petsclib::PetscLibType, mat::AbstractPetscMat, v::AbstractPetscVec) end

@for_petsc function MatGetRowSumAbs(petsclib::$UnionPetscLib, mat::AbstractPetscMat, v::AbstractPetscVec )

    @chk ccall(
               (:MatGetRowSumAbs, $petsc_library),
               PetscErrorCode,
               (CMat, CVec),
               mat, v,
              )


	return nothing
end 

"""
	MatGetRowUpperTriangular(petsclib::PetscLibType,mat::AbstractPetscMat) 
Sets a flag to enable calls to `MatGetRow()` for matrix in `MATSBAIJ` format.
You should call `MatRestoreRowUpperTriangular()` after calling` MatGetRow()` and `MatRestoreRow()` to disable the flag.

Not Collective

Input Parameter:
- `mat` - the matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATSBAIJ`, `MatRestoreRowUpperTriangular()`

# External Links
$(_doc_external("Mat/MatGetRowUpperTriangular"))
"""
function MatGetRowUpperTriangular(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatGetRowUpperTriangular(petsclib::$UnionPetscLib, mat::AbstractPetscMat )

    @chk ccall(
               (:MatGetRowUpperTriangular, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	S::PetscMat,Sp::PetscMat = MatGetSchurComplement(petsclib::PetscLibType,A::AbstractPetscMat, isrow0::AbstractIS, iscol0::AbstractIS, isrow1::AbstractIS, iscol1::AbstractIS, mreuse::MatReuse, ainvtype::MatSchurComplementAinvType, preuse::MatReuse) 
Obtain the Schur complement from eliminating part of the matrix in another part.

Collective

Input Parameters:
- `A`        - matrix in which the complement is to be taken
- `isrow0`   - rows to eliminate
- `iscol0`   - columns to eliminate, (isrow0,iscol0) should be square and nonsingular
- `isrow1`   - rows in which the Schur complement is formed
- `iscol1`   - columns in which the Schur complement is formed
- `mreuse`   - `MAT_INITIAL_MATRIX` or `MAT_REUSE_MATRIX`, use `MAT_IGNORE_MATRIX` to put nothing in `S`
- `ainvtype` - the type of approximation used for the inverse of the (0,0) block used in forming `Sp`:
`MAT_SCHUR_COMPLEMENT_AINV_DIAG`, `MAT_SCHUR_COMPLEMENT_AINV_LUMP`, `MAT_SCHUR_COMPLEMENT_AINV_BLOCK_DIAG`, or `MAT_SCHUR_COMPLEMENT_AINV_FULL`
- `preuse`   - `MAT_INITIAL_MATRIX` or `MAT_REUSE_MATRIX`, use `MAT_IGNORE_MATRIX` to put nothing in `Sp`

Output Parameters:
- `S`  - exact Schur complement, often of type `MATSCHURCOMPLEMENT` which is difficult to use for preconditioning
- `Sp` - approximate Schur complement from which a preconditioner can be built A11 - A10 inv(DIAGFORM(A00)) A01

Level: advanced

-seealso: [](ch_ksp), `MatCreateSubMatrix()`, `PCFIELDSPLIT`, `MatCreateSchurComplement()`, `MatSchurComplementAinvType`

# External Links
$(_doc_external("KSP/MatGetSchurComplement"))
"""
function MatGetSchurComplement(petsclib::PetscLibType, A::AbstractPetscMat, isrow0::AbstractIS, iscol0::AbstractIS, isrow1::AbstractIS, iscol1::AbstractIS, mreuse::MatReuse, ainvtype::MatSchurComplementAinvType, preuse::MatReuse) end

@for_petsc function MatGetSchurComplement(petsclib::$UnionPetscLib, A::AbstractPetscMat, isrow0::AbstractIS, iscol0::AbstractIS, isrow1::AbstractIS, iscol1::AbstractIS, mreuse::MatReuse, ainvtype::MatSchurComplementAinvType, preuse::MatReuse )
	S_ = Ref{CMat}()
	Sp_ = Ref{CMat}()

    @chk ccall(
               (:MatGetSchurComplement, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, CIS, CIS, CIS, MatReuse, Ptr{CMat}, MatSchurComplementAinvType, MatReuse, Ptr{CMat}),
               A, isrow0, iscol0, isrow1, iscol1, mreuse, S_, ainvtype, preuse, Sp_,
              )

	S = PetscMat(S_[], petsclib)
	Sp = PetscMat(Sp_[], petsclib)

	return S,Sp
end 

"""
	matstruct::PetscMat = MatGetSeqNonzeroStructure(petsclib::PetscLibType,mat::AbstractPetscMat) 
Extracts the nonzero structure from a matrix and stores it, in its entirety, on each process

Collective

Input Parameter:
- `mat` - the matrix

Output Parameter:
- `matstruct` - the sequential matrix with the nonzero structure of `mat`

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatDestroySeqNonzeroStructure()`, `MatCreateSubMatrices()`, `MatDestroyMatrices()`

# External Links
$(_doc_external("Mat/MatGetSeqNonzeroStructure"))
"""
function MatGetSeqNonzeroStructure(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatGetSeqNonzeroStructure(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	matstruct_ = Ref{CMat}()

    @chk ccall(
               (:MatGetSeqNonzeroStructure, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               mat, matstruct_,
              )

	matstruct = PetscMat(matstruct_[], petsclib)

	return matstruct
end 

"""
	m::PetscInt,n::PetscInt = MatGetSize(petsclib::PetscLibType,mat::AbstractPetscMat) 
Returns the numbers of rows and columns in a matrix.

Not Collective

Input Parameter:
- `mat` - the matrix

Output Parameters:
- `m` - the number of global rows
- `n` - the number of global columns

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatSetSizes()`, `MatGetLocalSize()`

# External Links
$(_doc_external("Mat/MatGetSize"))
"""
function MatGetSize(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatGetSize(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	m_ = Ref{$PetscInt}()
	n_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatGetSize, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}, Ptr{$PetscInt}),
               mat, m_, n_,
              )

	m = m_[]
	n = n_[]

	return m,n
end 

"""
	state::PetscObjectState = MatGetState(petsclib::PetscLibType,A::AbstractPetscMat) 
Gets the state of a `Mat`. Same value as returned by `PetscObjectStateGet()`

Not Collective

Input Parameter:
- `A` - the matrix

Output Parameter:
- `state` - the object state

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `PetscObjectStateGet()`, `MatGetNonzeroState()`

# External Links
$(_doc_external("Mat/MatGetState"))
"""
function MatGetState(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatGetState(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	state_ = Ref{PetscObjectState}()

    @chk ccall(
               (:MatGetState, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscObjectState}),
               A, state_,
              )

	state = state_[]

	return state
end 

"""
	trace::PetscScalar = MatGetTrace(petsclib::PetscLibType,mat::AbstractPetscMat) 
Gets the trace of a matrix. The sum of the diagonal entries.

Collective

Input Parameter:
- `mat` - the matrix

Output Parameter:
- `trace` - the sum of the diagonal entries

Level: advanced

-seealso: [](ch_matrices), `Mat`

# External Links
$(_doc_external("Mat/MatGetTrace"))
"""
function MatGetTrace(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatGetTrace(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	trace_ = Ref{$PetscScalar}()

    @chk ccall(
               (:MatGetTrace, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscScalar}),
               mat, trace_,
              )

	trace = trace_[]

	return trace
end 

"""
	MatGetTransposeNullSpace(petsclib::PetscLibType,mat::AbstractPetscMat, nullsp::MatNullSpace) 
retrieves the null space of the transpose of a matrix.

Logically Collective

Input Parameters:
- `mat`    - the matrix
- `nullsp` - the null space object

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatNullSpace`, `MatCreate()`, `MatNullSpaceCreate()`, `MatSetNearNullSpace()`, `MatSetTransposeNullSpace()`, `MatSetNullSpace()`, `MatGetNullSpace()`

# External Links
$(_doc_external("Mat/MatGetTransposeNullSpace"))
"""
function MatGetTransposeNullSpace(petsclib::PetscLibType, mat::AbstractPetscMat, nullsp::MatNullSpace) end

@for_petsc function MatGetTransposeNullSpace(petsclib::$UnionPetscLib, mat::AbstractPetscMat, nullsp::MatNullSpace )

    @chk ccall(
               (:MatGetTransposeNullSpace, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatNullSpace}),
               mat, nullsp,
              )


	return nothing
end 

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

"""
	va::PetscScalar = MatGetValue(petsclib::PetscLibType,mat::AbstractPetscMat, row::PetscInt, col::PetscInt) 

# External Links
$(_doc_external("Mat/MatGetValue"))
"""
function MatGetValue(petsclib::PetscLibType, mat::AbstractPetscMat, row::PetscInt, col::PetscInt) end

@for_petsc function MatGetValue(petsclib::$UnionPetscLib, mat::AbstractPetscMat, row::$PetscInt, col::$PetscInt )
	va_ = Ref{$PetscScalar}()

    @chk ccall(
               (:MatGetValue, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, Ptr{$PetscScalar}),
               mat, row, col, va_,
              )

	va = va_[]

	return va
end 

"""
	MatGetValues(petsclib::PetscLibType,mat::AbstractPetscMat, m::PetscInt, idxm::Vector{PetscInt}, n::PetscInt, idxn::Vector{PetscInt}, v::Vector{PetscScalar}) 
Gets a block of local values from a matrix.

Not Collective; can only return values that are owned by the give process

Input Parameters:
- `mat`  - the matrix
- `v`    - a logically two-dimensional array for storing the values
- `m`    - the number of rows
- `idxm` - the  global indices of the rows
- `n`    - the number of columns
- `idxn` - the global indices of the columns

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatGetRow()`, `MatCreateSubMatrices()`, `MatSetValues()`, `MatGetOwnershipRange()`, `MatGetValuesLocal()`, `MatGetValue()`

# External Links
$(_doc_external("Mat/MatGetValues"))
"""
function MatGetValues(petsclib::PetscLibType, mat::AbstractPetscMat, m::PetscInt, idxm::Vector{PetscInt}, n::PetscInt, idxn::Vector{PetscInt}, v::Vector{PetscScalar}) end

@for_petsc function MatGetValues(petsclib::$UnionPetscLib, mat::AbstractPetscMat, m::$PetscInt, idxm::Vector{$PetscInt}, n::$PetscInt, idxn::Vector{$PetscInt}, v::Vector{$PetscScalar} )

    @chk ccall(
               (:MatGetValues, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}),
               mat, m, idxm, n, idxn, v,
              )


	return nothing
end 

"""
	MatGetValuesLocal(petsclib::PetscLibType,mat::AbstractPetscMat, nrow::PetscInt, irow::Vector{PetscInt}, ncol::PetscInt, icol::Vector{PetscInt}, y::Vector{PetscScalar}) 
retrieves values from certain locations in a matrix using the local numbering of the indices
defined previously by `MatSetLocalToGlobalMapping()`

Not Collective

Input Parameters:
- `mat`  - the matrix
- `nrow` - number of rows
- `irow` - the row local indices
- `ncol` - number of columns
- `icol` - the column local indices

Output Parameter:
- `y` - a one-dimensional array that contains the values implicitly stored as a two-dimensional array, by default in row-major order.
See `MAT_ROW_ORIENTED` in `MatSetOption()` for how to use column-major order.

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatAssemblyBegin()`, `MatAssemblyEnd()`, `MatSetValues()`, `MatSetLocalToGlobalMapping()`,
`MatSetValuesLocal()`, `MatGetValues()`

# External Links
$(_doc_external("Mat/MatGetValuesLocal"))
"""
function MatGetValuesLocal(petsclib::PetscLibType, mat::AbstractPetscMat, nrow::PetscInt, irow::Vector{PetscInt}, ncol::PetscInt, icol::Vector{PetscInt}, y::Vector{PetscScalar}) end

@for_petsc function MatGetValuesLocal(petsclib::$UnionPetscLib, mat::AbstractPetscMat, nrow::$PetscInt, irow::Vector{$PetscInt}, ncol::$PetscInt, icol::Vector{$PetscInt}, y::Vector{$PetscScalar} )

    @chk ccall(
               (:MatGetValuesLocal, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}),
               mat, nrow, irow, ncol, icol, y,
              )


	return nothing
end 

"""
	nblocks::PetscInt,bsizes::Ptr{PetscInt} = MatGetVariableBlockSizes(petsclib::PetscLibType,mat::AbstractPetscMat) 
Gets a diagonal blocks of the matrix that need not be of the same size

Not Collective; No Fortran Support

Input Parameter:
- `mat` - the matrix

Output Parameters:
- `nblocks` - the number of blocks on this process
- `bsizes`  - the block sizes

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatCreateSeqBAIJ()`, `MatCreateBAIJ()`, `MatGetBlockSize()`, `MatSetBlockSizes()`, `MatGetBlockSizes()`, `MatSetVariableBlockSizes()`, `MatComputeVariableBlockEnvelope()`

# External Links
$(_doc_external("Mat/MatGetVariableBlockSizes"))
"""
function MatGetVariableBlockSizes(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatGetVariableBlockSizes(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	nblocks_ = Ref{$PetscInt}()
	bsizes_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:MatGetVariableBlockSizes, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
               mat, nblocks_, bsizes_,
              )

	nblocks = nblocks_[]
	bsizes = bsizes_[]

	return nblocks,bsizes
end 

"""
	vtype::VecType = MatGetVecType(petsclib::PetscLibType,mat::AbstractPetscMat) 
Gets the vector type the matrix will return with `MatCreateVecs()`

Not Collective

Input Parameter:
- `mat` - the matrix

Output Parameter:
- `vtype` - name of vector type

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatType`, `MatSetVecType()`, `VecType`

# External Links
$(_doc_external("Mat/MatGetVecType"))
"""
function MatGetVecType(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatGetVecType(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	vtype_ = Ref{VecType}()

    @chk ccall(
               (:MatGetVecType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{VecType}),
               mat, vtype_,
              )

	vtype = vtype_[] == C_NULL ? "" : unsafe_string(vtype_[])

	return vtype
end 

"""
	parcsr::Ptr{hypre_ParCSRMatrix} = MatHYPREGetParCSR(petsclib::PetscLibType,A::AbstractPetscMat) 
Gets the pointer to the ParCSR matrix

Not Collective, No Fortran Support

Input Parameter:
- `A` - the `MATHYPRE` object

Output Parameter:
- `parcsr` - the pointer to the `hypre_ParCSRMatrix`

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATHYPRE`, `PetscCopyMode`

# External Links
$(_doc_external("Mat/MatHYPREGetParCSR"))
"""
function MatHYPREGetParCSR(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatHYPREGetParCSR(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	parcsr_ = Ref{Ptr{hypre_ParCSRMatrix}}()

    @chk ccall(
               (:MatHYPREGetParCSR, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{hypre_ParCSRMatrix}}),
               A, parcsr_,
              )

	parcsr = parcsr_[]

	return parcsr
end 

"""
	MatHYPRESetPreallocation(petsclib::PetscLibType,A::AbstractPetscMat, dnz::PetscInt, dnnz::Vector{PetscInt}, onz::PetscInt, onnz::Vector{PetscInt}) 
Preallocates memory for a sparse parallel matrix in HYPRE IJ format

Collective

Input Parameters:
- `A`    - the matrix
- `dnz`  - number of nonzeros per row in DIAGONAL portion of local submatrix
(same value is used for all local rows)
- `dnnz` - array containing the number of nonzeros in the various rows of the
DIAGONAL portion of the local submatrix (possibly different for each row)
or `NULL` (`PETSC_NULL_INTEGER` in Fortran), if `d_nz` is used to specify the nonzero structure.
The size of this array is equal to the number of local rows, i.e `m`.
For matrices that will be factored, you must leave room for (and set)
the diagonal entry even if it is zero.
- `onz`  - number of nonzeros per row in the OFF-DIAGONAL portion of local
submatrix (same value is used for all local rows).
- `onnz` - array containing the number of nonzeros in the various rows of the
OFF-DIAGONAL portion of the local submatrix (possibly different for
each row) or `NULL` (`PETSC_NULL_INTEGER` in Fortran), if `o_nz` is used to specify the nonzero
structure. The size of this array is equal to the number
of local rows, i.e `m`.

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `MatMPIAIJSetPreallocation()`, `MATHYPRE`, `MATAIJ`

# External Links
$(_doc_external("Mat/MatHYPRESetPreallocation"))
"""
function MatHYPRESetPreallocation(petsclib::PetscLibType, A::AbstractPetscMat, dnz::PetscInt, dnnz::Vector{PetscInt}, onz::PetscInt, onnz::Vector{PetscInt}) end

@for_petsc function MatHYPRESetPreallocation(petsclib::$UnionPetscLib, A::AbstractPetscMat, dnz::$PetscInt, dnnz::Vector{$PetscInt}, onz::$PetscInt, onnz::Vector{$PetscInt} )

    @chk ccall(
               (:MatHYPRESetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}),
               A, dnz, dnnz, onz, onnz,
              )


	return nothing
end 

"""
	cong::PetscBool = MatHasCongruentLayouts(petsclib::PetscLibType,mat::AbstractPetscMat) 
Determines whether the rows and columns layouts of the matrix are congruent

Collective

Input Parameter:
- `mat` - the matrix

Output Parameter:
- `cong` - either `PETSC_TRUE` or `PETSC_FALSE`

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `MatSetSizes()`, `PetscLayout`

# External Links
$(_doc_external("Mat/MatHasCongruentLayouts"))
"""
function MatHasCongruentLayouts(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatHasCongruentLayouts(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	cong_ = Ref{PetscBool}()

    @chk ccall(
               (:MatHasCongruentLayouts, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscBool}),
               mat, cong_,
              )

	cong = cong_[]

	return cong
end 

"""
	has::PetscBool = MatHasOperation(petsclib::PetscLibType,mat::AbstractPetscMat, op::MatOperation) 
Determines whether the given matrix supports the particular operation.

Not Collective

Input Parameters:
- `mat` - the matrix
- `op`  - the operation, for example, `MATOP_GET_DIAGONAL`

Output Parameter:
- `has` - either `PETSC_TRUE` or `PETSC_FALSE`

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatCreateShell()`, `MatGetOperation()`, `MatSetOperation()`

# External Links
$(_doc_external("Mat/MatHasOperation"))
"""
function MatHasOperation(petsclib::PetscLibType, mat::AbstractPetscMat, op::MatOperation) end

@for_petsc function MatHasOperation(petsclib::$UnionPetscLib, mat::AbstractPetscMat, op::MatOperation )
	has_ = Ref{PetscBool}()

    @chk ccall(
               (:MatHasOperation, $petsc_library),
               PetscErrorCode,
               (CMat, MatOperation, Ptr{PetscBool}),
               mat, op, has_,
              )

	has = has_[]

	return has
end 

"""
	MatHeaderMerge(petsclib::PetscLibType,A::AbstractPetscMat, C::AbstractPetscMat) 
Merges some information from the header of `C` to `A`; the `C` object is then destroyed

Collective, No Fortran Support

Input Parameters:
- `A` - a `Mat` being merged into
- `C` - the `Mat` providing the merge information

Level: developer

-seealso: `Mat`, `MatHeaderReplace()`

# External Links
$(_doc_external("Mat/MatHeaderMerge"))
"""
function MatHeaderMerge(petsclib::PetscLibType, A::AbstractPetscMat, C::AbstractPetscMat) end

@for_petsc function MatHeaderMerge(petsclib::$UnionPetscLib, A::AbstractPetscMat, C::AbstractPetscMat )
	C_ = Ref(C.ptr)

    @chk ccall(
               (:MatHeaderMerge, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, C_,
              )

	C.ptr = C_[]

	return nothing
end 

"""
	MatHeaderReplace(petsclib::PetscLibType,A::AbstractPetscMat, C::AbstractPetscMat) 
Replaces the internal data of matrix `A` by the internal data of matrix `C` while deleting the outer wrapper of `C`

Input Parameters:
- `A` - a `Mat` whose internal data is to be replaced
- `C` - the `Mat` providing new internal data for `A`

Level: advanced

-seealso: `Mat`, `MatHeaderMerge()`

# External Links
$(_doc_external("Mat/MatHeaderReplace"))
"""
function MatHeaderReplace(petsclib::PetscLibType, A::AbstractPetscMat, C::AbstractPetscMat) end

@for_petsc function MatHeaderReplace(petsclib::$UnionPetscLib, A::AbstractPetscMat, C::AbstractPetscMat )
	C_ = Ref(C.ptr)

    @chk ccall(
               (:MatHeaderReplace, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, C_,
              )

	C.ptr = C_[]

	return nothing
end 

"""
	B::PetscMat = MatHermitianTranspose(petsclib::PetscLibType,mat::AbstractPetscMat, reuse::MatReuse) 
Computes an in

Collective

Input Parameters:
- `mat`   - the matrix to transpose and complex conjugate
- `reuse` - either `MAT_INITIAL_MATRIX`, `MAT_REUSE_MATRIX`, or `MAT_INPLACE_MATRIX`

Output Parameter:
- `B` - the Hermitian transpose

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatTranspose()`, `MatMultTranspose()`, `MatMultTransposeAdd()`, `MatIsTranspose()`, `MatReuse`

# External Links
$(_doc_external("Mat/MatHermitianTranspose"))
"""
function MatHermitianTranspose(petsclib::PetscLibType, mat::AbstractPetscMat, reuse::MatReuse) end

@for_petsc function MatHermitianTranspose(petsclib::$UnionPetscLib, mat::AbstractPetscMat, reuse::MatReuse )
	B_ = Ref{CMat}()

    @chk ccall(
               (:MatHermitianTranspose, $petsc_library),
               PetscErrorCode,
               (CMat, MatReuse, Ptr{CMat}),
               mat, reuse, B_,
              )

	B = PetscMat(B_[], petsclib)

	return B
end 

"""
	M::PetscMat = MatHermitianTransposeGetMat(petsclib::PetscLibType,A::AbstractPetscMat) 
Gets the `Mat` object stored inside a `MATHERMITIANTRANSPOSEVIRTUAL`

Logically Collective

Input Parameter:
- `A` - the `MATHERMITIANTRANSPOSEVIRTUAL` matrix

Output Parameter:
- `M` - the matrix object stored inside A

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATHERMITIANTRANSPOSEVIRTUAL`, `MatCreateHermitianTranspose()`

# External Links
$(_doc_external("Mat/MatHermitianTransposeGetMat"))
"""
function MatHermitianTransposeGetMat(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatHermitianTransposeGetMat(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	M_ = Ref{CMat}()

    @chk ccall(
               (:MatHermitianTransposeGetMat, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, M_,
              )

	M = PetscMat(M_[], petsclib)

	return M
end 

"""
	MatHtoolGetPermutationSource(petsclib::PetscLibType,A::AbstractPetscMat, is::AbstractIS) 

# External Links
$(_doc_external("Mat/MatHtoolGetPermutationSource"))
"""
function MatHtoolGetPermutationSource(petsclib::PetscLibType, A::AbstractPetscMat, is::AbstractIS) end

@for_petsc function MatHtoolGetPermutationSource(petsclib::$UnionPetscLib, A::AbstractPetscMat, is::AbstractIS )
	is_ = Ref(is.ptr)

    @chk ccall(
               (:MatHtoolGetPermutationSource, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CIS}),
               A, is_,
              )

	is.ptr = is_[]

	return nothing
end 

"""
	MatHtoolGetPermutationTarget(petsclib::PetscLibType,A::AbstractPetscMat, is::AbstractIS) 

# External Links
$(_doc_external("Mat/MatHtoolGetPermutationTarget"))
"""
function MatHtoolGetPermutationTarget(petsclib::PetscLibType, A::AbstractPetscMat, is::AbstractIS) end

@for_petsc function MatHtoolGetPermutationTarget(petsclib::$UnionPetscLib, A::AbstractPetscMat, is::AbstractIS )
	is_ = Ref(is.ptr)

    @chk ccall(
               (:MatHtoolGetPermutationTarget, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CIS}),
               A, is_,
              )

	is.ptr = is_[]

	return nothing
end 

"""
	MatHtoolSetKernel(petsclib::PetscLibType,A::AbstractPetscMat, kernel::Ptr{Cvoid}, kernelctx::Ptr{Cvoid}) 

# External Links
$(_doc_external("Mat/MatHtoolSetKernel"))
"""
function MatHtoolSetKernel(petsclib::PetscLibType, A::AbstractPetscMat, kernel::Ptr{Cvoid}, kernelctx::Ptr{Cvoid}) end

@for_petsc function MatHtoolSetKernel(petsclib::$UnionPetscLib, A::AbstractPetscMat, kernel::Ptr{Cvoid}, kernelctx::Ptr{Cvoid} )

    @chk ccall(
               (:MatHtoolSetKernel, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Cvoid}, Ptr{Cvoid}),
               A, kernel, kernelctx,
              )


	return nothing
end 

"""
	MatHtoolUsePermutation(petsclib::PetscLibType,A::AbstractPetscMat, use::PetscBool) 

# External Links
$(_doc_external("Mat/MatHtoolUsePermutation"))
"""
function MatHtoolUsePermutation(petsclib::PetscLibType, A::AbstractPetscMat, use::PetscBool) end

@for_petsc function MatHtoolUsePermutation(petsclib::$UnionPetscLib, A::AbstractPetscMat, use::PetscBool )

    @chk ccall(
               (:MatHtoolUsePermutation, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               A, use,
              )


	return nothing
end 

"""
	MatHtoolUseRecompression(petsclib::PetscLibType,A::AbstractPetscMat, use::PetscBool) 

# External Links
$(_doc_external("Mat/MatHtoolUseRecompression"))
"""
function MatHtoolUseRecompression(petsclib::PetscLibType, A::AbstractPetscMat, use::PetscBool) end

@for_petsc function MatHtoolUseRecompression(petsclib::$UnionPetscLib, A::AbstractPetscMat, use::PetscBool )

    @chk ccall(
               (:MatHtoolUseRecompression, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               A, use,
              )


	return nothing
end 

"""
	MatICCFactor(petsclib::PetscLibType,mat::AbstractPetscMat, row::AbstractIS, info::Vector{MatFactorInfo}) 
Performs in

Collective

Input Parameters:
- `mat`  - the matrix
- `row`  - row/column permutation
- `info` - information on desired factorization process

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatFactorInfo`, `MatGetFactor()`, `MatICCFactorSymbolic()`, `MatLUFactorNumeric()`, `MatCholeskyFactor()`

# External Links
$(_doc_external("Mat/MatICCFactor"))
"""
function MatICCFactor(petsclib::PetscLibType, mat::AbstractPetscMat, row::AbstractIS, info::Vector{MatFactorInfo}) end

@for_petsc function MatICCFactor(petsclib::$UnionPetscLib, mat::AbstractPetscMat, row::AbstractIS, info::Vector{MatFactorInfo} )

    @chk ccall(
               (:MatICCFactor, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, Ptr{MatFactorInfo}),
               mat, row, info,
              )


	return nothing
end 

"""
	MatICCFactorSymbolic(petsclib::PetscLibType,fact::AbstractPetscMat, mat::AbstractPetscMat, perm::AbstractIS, info::Vector{MatFactorInfo}) 
Performs symbolic incomplete
Cholesky factorization for a symmetric matrix.  Use
`MatCholeskyFactorNumeric()` to complete the factorization.

Collective

Input Parameters:
- `fact` - the factorized matrix obtained with `MatGetFactor()`
- `mat`  - the matrix to be factored
- `perm` - row and column permutation
- `info` - structure containing
-seealso: [](ch_matrices), `Mat`, `MatGetFactor()`, `MatCholeskyFactorNumeric()`, `MatCholeskyFactor()`, `MatFactorInfo`

# External Links
$(_doc_external("Mat/MatICCFactorSymbolic"))
"""
function MatICCFactorSymbolic(petsclib::PetscLibType, fact::AbstractPetscMat, mat::AbstractPetscMat, perm::AbstractIS, info::Vector{MatFactorInfo}) end

@for_petsc function MatICCFactorSymbolic(petsclib::$UnionPetscLib, fact::AbstractPetscMat, mat::AbstractPetscMat, perm::AbstractIS, info::Vector{MatFactorInfo} )

    @chk ccall(
               (:MatICCFactorSymbolic, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CIS, Ptr{MatFactorInfo}),
               fact, mat, perm, info,
              )


	return nothing
end 

"""
	MatILUFactor(petsclib::PetscLibType,mat::AbstractPetscMat, row::AbstractIS, col::AbstractIS, info::Vector{MatFactorInfo}) 
Performs in

Collective

Input Parameters:
- `mat`  - the matrix
- `row`  - row permutation
- `col`  - column permutation
- `info` - structure containing
-seealso: [](ch_matrices), `Mat`, [Matrix Factorization](sec_matfactor), `MatILUFactorSymbolic()`, `MatLUFactorNumeric()`, `MatCholeskyFactor()`, `MatFactorInfo`

# External Links
$(_doc_external("Mat/MatILUFactor"))
"""
function MatILUFactor(petsclib::PetscLibType, mat::AbstractPetscMat, row::AbstractIS, col::AbstractIS, info::Vector{MatFactorInfo}) end

@for_petsc function MatILUFactor(petsclib::$UnionPetscLib, mat::AbstractPetscMat, row::AbstractIS, col::AbstractIS, info::Vector{MatFactorInfo} )

    @chk ccall(
               (:MatILUFactor, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, CIS, Ptr{MatFactorInfo}),
               mat, row, col, info,
              )


	return nothing
end 

"""
	MatILUFactorSymbolic(petsclib::PetscLibType,fact::AbstractPetscMat, mat::AbstractPetscMat, row::AbstractIS, col::AbstractIS, info::Vector{MatFactorInfo}) 
Performs symbolic ILU factorization of a matrix obtained with `MatGetFactor()`
Uses levels of fill only, not drop tolerance. Use `MatLUFactorNumeric()`
to complete the factorization.

Collective

Input Parameters:
- `fact` - the factorized matrix obtained with `MatGetFactor()`
- `mat`  - the matrix
- `row`  - row permutation
- `col`  - column permutation
- `info` - structure containing
-seealso: [](ch_matrices), `Mat`, [Matrix Factorization](sec_matfactor), `MatGetFactor()`, `MatLUFactorSymbolic()`, `MatLUFactorNumeric()`, `MatCholeskyFactor()`
`MatGetOrdering()`, `MatFactorInfo`

# External Links
$(_doc_external("Mat/MatILUFactorSymbolic"))
"""
function MatILUFactorSymbolic(petsclib::PetscLibType, fact::AbstractPetscMat, mat::AbstractPetscMat, row::AbstractIS, col::AbstractIS, info::Vector{MatFactorInfo}) end

@for_petsc function MatILUFactorSymbolic(petsclib::$UnionPetscLib, fact::AbstractPetscMat, mat::AbstractPetscMat, row::AbstractIS, col::AbstractIS, info::Vector{MatFactorInfo} )

    @chk ccall(
               (:MatILUFactorSymbolic, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CIS, CIS, Ptr{MatFactorInfo}),
               fact, mat, row, col, info,
              )


	return nothing
end 

"""
	MatISFixLocalEmpty(petsclib::PetscLibType,A::AbstractPetscMat, fix::PetscBool) 
Compress out zero local rows from the local matrices

Logically Collective

Input Parameters:
- `A`   - the matrix
- `fix` - the boolean flag

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATIS`, `MatCreate()`, `MatCreateIS()`, `MatISSetPreallocation()`, `MatAssemblyEnd()`, `MAT_FINAL_ASSEMBLY`

# External Links
$(_doc_external("Mat/MatISFixLocalEmpty"))
"""
function MatISFixLocalEmpty(petsclib::PetscLibType, A::AbstractPetscMat, fix::PetscBool) end

@for_petsc function MatISFixLocalEmpty(petsclib::$UnionPetscLib, A::AbstractPetscMat, fix::PetscBool )

    @chk ccall(
               (:MatISFixLocalEmpty, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               A, fix,
              )


	return nothing
end 

"""
	flg::PetscBool = MatISGetAllowRepeated(petsclib::PetscLibType,A::AbstractPetscMat) 
Get the flag to allow repeated entries in the local to global map

Not Collective

Input Parameter:
- `A` - the matrix

Output Parameter:
- `flg` - the boolean flag

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `MatCreateIS()`, `MatSetLocalToGlobalMapping()`, `MatISSetAllowRepeated()`

# External Links
$(_doc_external("Mat/MatISGetAllowRepeated"))
"""
function MatISGetAllowRepeated(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatISGetAllowRepeated(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatISGetAllowRepeated, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscBool}),
               A, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	MatISGetLocalMat(petsclib::PetscLibType,mat::AbstractPetscMat, loc::AbstractPetscMat) 
Gets the local matrix stored inside a `MATIS` matrix.

Not Collective.

Input Parameter:
- `mat` - the matrix

Output Parameter:
- `local` - the local matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATIS`, `MatISRestoreLocalMat()`

# External Links
$(_doc_external("Mat/MatISGetLocalMat"))
"""
function MatISGetLocalMat(petsclib::PetscLibType, mat::AbstractPetscMat, loc::AbstractPetscMat) end

@for_petsc function MatISGetLocalMat(petsclib::$UnionPetscLib, mat::AbstractPetscMat, loc::AbstractPetscMat )
	loc_ = Ref(loc.ptr)

    @chk ccall(
               (:MatISGetLocalMat, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               mat, loc_,
              )

	loc.ptr = loc_[]

	return nothing
end 

"""
	rmapping::ISLocalToGlobalMapping,cmapping::ISLocalToGlobalMapping = MatISGetLocalToGlobalMapping(petsclib::PetscLibType,A::AbstractPetscMat) 
Gets the local

Not Collective

Input Parameter:
- `A` - the matrix

Output Parameters:
- `rmapping` - row mapping
- `cmapping` - column mapping

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATIS`, `MatSetLocalToGlobalMapping()`

# External Links
$(_doc_external("Mat/MatISGetLocalToGlobalMapping"))
"""
function MatISGetLocalToGlobalMapping(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatISGetLocalToGlobalMapping(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	rmapping_ = Ref{ISLocalToGlobalMapping}()
	cmapping_ = Ref{ISLocalToGlobalMapping}()

    @chk ccall(
               (:MatISGetLocalToGlobalMapping, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{ISLocalToGlobalMapping}, Ptr{ISLocalToGlobalMapping}),
               A, rmapping_, cmapping_,
              )

	rmapping = rmapping_[]
	cmapping = cmapping_[]

	return rmapping,cmapping
end 

"""
	MatISRestoreLocalMat(petsclib::PetscLibType,mat::AbstractPetscMat, loc::AbstractPetscMat) 
Restores the local matrix obtained with `MatISGetLocalMat()`

Not Collective.

Input Parameters:
- `mat`   - the matrix
- `local` - the local matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATIS`, `MatISGetLocalMat()`

# External Links
$(_doc_external("Mat/MatISRestoreLocalMat"))
"""
function MatISRestoreLocalMat(petsclib::PetscLibType, mat::AbstractPetscMat, loc::AbstractPetscMat) end

@for_petsc function MatISRestoreLocalMat(petsclib::$UnionPetscLib, mat::AbstractPetscMat, loc::AbstractPetscMat )
	loc_ = Ref(loc.ptr)

    @chk ccall(
               (:MatISRestoreLocalMat, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               mat, loc_,
              )

	loc.ptr = loc_[]

	return nothing
end 

"""
	MatISSetAllowRepeated(petsclib::PetscLibType,A::AbstractPetscMat, flg::PetscBool) 
Set the flag to allow repeated entries in the local to global map

Logically Collective

Input Parameters:
- `A`   - the matrix
- `flg` - the boolean flag

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `MatCreateIS()`, `MatSetLocalToGlobalMapping()`, `MatISGetAllowRepeated()`

# External Links
$(_doc_external("Mat/MatISSetAllowRepeated"))
"""
function MatISSetAllowRepeated(petsclib::PetscLibType, A::AbstractPetscMat, flg::PetscBool) end

@for_petsc function MatISSetAllowRepeated(petsclib::$UnionPetscLib, A::AbstractPetscMat, flg::PetscBool )

    @chk ccall(
               (:MatISSetAllowRepeated, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               A, flg,
              )


	return nothing
end 

"""
	MatISSetLocalMat(petsclib::PetscLibType,mat::AbstractPetscMat, loc::AbstractPetscMat) 
Replace the local matrix stored inside a `MATIS` object.

Not Collective

Input Parameters:
- `mat`   - the matrix
- `local` - the local matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATIS`, `MatISSetLocalMatType`, `MatISGetLocalMat()`

# External Links
$(_doc_external("Mat/MatISSetLocalMat"))
"""
function MatISSetLocalMat(petsclib::PetscLibType, mat::AbstractPetscMat, loc::AbstractPetscMat) end

@for_petsc function MatISSetLocalMat(petsclib::$UnionPetscLib, mat::AbstractPetscMat, loc::AbstractPetscMat )

    @chk ccall(
               (:MatISSetLocalMat, $petsc_library),
               PetscErrorCode,
               (CMat, CMat),
               mat, loc,
              )


	return nothing
end 

"""
	MatISSetLocalMatType(petsclib::PetscLibType,mat::AbstractPetscMat, mtype::MatType) 
Specifies the type of local matrix inside the `MATIS`

Logically Collective.

Input Parameters:
- `mat`   - the matrix
- `mtype` - the local matrix type

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATIS`, `MatSetType()`, `MatType`

# External Links
$(_doc_external("Mat/MatISSetLocalMatType"))
"""
function MatISSetLocalMatType(petsclib::PetscLibType, mat::AbstractPetscMat, mtype::MatType) end

@for_petsc function MatISSetLocalMatType(petsclib::$UnionPetscLib, mat::AbstractPetscMat, mtype::MatType )

    @chk ccall(
               (:MatISSetLocalMatType, $petsc_library),
               PetscErrorCode,
               (CMat, MatType),
               mat, mtype,
              )


	return nothing
end 

"""
	MatISSetPreallocation(petsclib::PetscLibType,B::AbstractPetscMat, d_nz::PetscInt, d_nnz::Union{Ptr, Vector{PetscInt}}, o_nz::PetscInt, o_nnz::Union{Ptr, Vector{PetscInt}}) 
Preallocates memory for a `MATIS` parallel matrix.

Collective

Input Parameters:
- `B`     - the matrix
- `d_nz`  - number of nonzeros per row in DIAGONAL portion of local submatrix
(same value is used for all local rows)
- `d_nnz` - array containing the number of nonzeros in the various rows of the
DIAGONAL portion of the local submatrix (possibly different for each row)
or `NULL`, if `d_nz` is used to specify the nonzero structure.
The size of this array is equal to the number of local rows, i.e `m`.
For matrices that will be factored, you must leave room for (and set)
the diagonal entry even if it is zero.
- `o_nz`  - number of nonzeros per row in the OFF-DIAGONAL portion of local
submatrix (same value is used for all local rows).
- `o_nnz` - array containing the number of nonzeros in the various rows of the
OFF-DIAGONAL portion of the local submatrix (possibly different for
each row) or `NULL`, if `o_nz` is used to specify the nonzero
structure. The size of this array is equal to the number
of local rows, i.e `m`.

If the *_nnz parameter is given then the *_nz parameter is ignored

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `MatCreateIS()`, `MatMPIAIJSetPreallocation()`, `MatISGetLocalMat()`, `MATIS`

# External Links
$(_doc_external("Mat/MatISSetPreallocation"))
"""
function MatISSetPreallocation(petsclib::PetscLibType, B::AbstractPetscMat, d_nz::PetscInt, d_nnz::Union{Ptr, Vector{PetscInt}}, o_nz::PetscInt, o_nnz::Union{Ptr, Vector{PetscInt}}) end

@for_petsc function MatISSetPreallocation(petsclib::$UnionPetscLib, B::AbstractPetscMat, d_nz::$PetscInt, d_nnz::Union{Ptr, Vector{$PetscInt}}, o_nz::$PetscInt, o_nnz::Union{Ptr, Vector{$PetscInt}} )

    @chk ccall(
               (:MatISSetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}),
               B, d_nz, d_nnz, o_nz, o_nnz,
              )


	return nothing
end 

"""
	MatISStoreL2L(petsclib::PetscLibType,A::AbstractPetscMat, store::PetscBool) 
Store local

Logically Collective

Input Parameters:
- `A`     - the matrix
- `store` - the boolean flag

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `MatCreateIS()`, `MatISSetPreallocation()`, `MatPtAP()`

# External Links
$(_doc_external("Mat/MatISStoreL2L"))
"""
function MatISStoreL2L(petsclib::PetscLibType, A::AbstractPetscMat, store::PetscBool) end

@for_petsc function MatISStoreL2L(petsclib::$UnionPetscLib, A::AbstractPetscMat, store::PetscBool )

    @chk ccall(
               (:MatISStoreL2L, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               A, store,
              )


	return nothing
end 

"""
	MatImaginaryPart(petsclib::PetscLibType,mat::AbstractPetscMat) 
Moves the imaginary part of the matrix to the real part and zeros the imaginary part

Logically Collective

Input Parameter:
- `mat` - the matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatRealPart()`

# External Links
$(_doc_external("Mat/MatImaginaryPart"))
"""
function MatImaginaryPart(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatImaginaryPart(petsclib::$UnionPetscLib, mat::AbstractPetscMat )

    @chk ccall(
               (:MatImaginaryPart, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	MatIncreaseOverlap(petsclib::PetscLibType,mat::AbstractPetscMat, n::PetscInt, is::Vector{<:AbstractIS}, ov::PetscInt) 
Given a set of submatrices indicated by index sets,
replaces the index sets by larger ones that represent submatrices with
additional overlap.

Collective

Input Parameters:
- `mat` - the matrix
- `n`   - the number of index sets
- `is`  - the array of index sets (these index sets will changed during the call)
- `ov`  - the additional overlap requested

Options Database Key:
- `-mat_increase_overlap_scalable` - use a scalable algorithm to compute the overlap (supported by MPIAIJ matrix)

Level: developer

-seealso: [](ch_matrices), `Mat`, `PCASM`, `MatSetBlockSize()`, `MatIncreaseOverlapSplit()`, `MatCreateSubMatrices()`

# External Links
$(_doc_external("Mat/MatIncreaseOverlap"))
"""
function MatIncreaseOverlap(petsclib::PetscLibType, mat::AbstractPetscMat, n::PetscInt, is::Vector{<:AbstractIS}, ov::PetscInt) end

@for_petsc function MatIncreaseOverlap(petsclib::$UnionPetscLib, mat::AbstractPetscMat, n::$PetscInt, is::Vector{<:AbstractIS}, ov::$PetscInt )

    @chk ccall(
               (:MatIncreaseOverlap, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CIS}, $PetscInt),
               mat, n, is, ov,
              )


	return nothing
end 

"""
	MatIncreaseOverlapSplit(petsclib::PetscLibType,mat::AbstractPetscMat, n::PetscInt, is::Vector{<:AbstractIS}, ov::PetscInt) 
Given a set of submatrices indicated by index sets across
a sub communicator, replaces the index sets by larger ones that represent submatrices with
additional overlap.

Collective

Input Parameters:
- `mat` - the matrix
- `n`   - the number of index sets
- `is`  - the array of index sets (these index sets will changed during the call)
- `ov`  - the additional overlap requested

`   Options Database Key:
- `-mat_increase_overlap_scalable` - use a scalable algorithm to compute the overlap (supported by MPIAIJ matrix)

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatCreateSubMatrices()`, `MatIncreaseOverlap()`

# External Links
$(_doc_external("Mat/MatIncreaseOverlapSplit"))
"""
function MatIncreaseOverlapSplit(petsclib::PetscLibType, mat::AbstractPetscMat, n::PetscInt, is::Vector{<:AbstractIS}, ov::PetscInt) end

@for_petsc function MatIncreaseOverlapSplit(petsclib::$UnionPetscLib, mat::AbstractPetscMat, n::$PetscInt, is::Vector{<:AbstractIS}, ov::$PetscInt )

    @chk ccall(
               (:MatIncreaseOverlapSplit, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CIS}, $PetscInt),
               mat, n, is, ov,
              )


	return nothing
end 

"""
	MatInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `Mat` package. It is called
from `PetscDLLibraryRegister_petscmat()` when using dynamic libraries, and on the first call to `MatCreate()`
when using shared or static libraries.

Level: developer

-seealso: [](ch_matrices), `Mat`, `PetscInitialize()`, `MatFinalizePackage()`

# External Links
$(_doc_external("Mat/MatInitializePackage"))
"""
function MatInitializePackage(petsclib::PetscLibType) end

@for_petsc function MatInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:MatInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	MatInodeAdjustForInodes(petsclib::PetscLibType,A::AbstractPetscMat, rperm::AbstractIS, cperm::AbstractIS) 

# External Links
$(_doc_external("Mat/MatInodeAdjustForInodes"))
"""
function MatInodeAdjustForInodes(petsclib::PetscLibType, A::AbstractPetscMat, rperm::AbstractIS, cperm::AbstractIS) end

@for_petsc function MatInodeAdjustForInodes(petsclib::$UnionPetscLib, A::AbstractPetscMat, rperm::AbstractIS, cperm::AbstractIS )
	rperm_ = Ref(rperm.ptr)
	cperm_ = Ref(cperm.ptr)

    @chk ccall(
               (:MatInodeAdjustForInodes, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CIS}, Ptr{CIS}),
               A, rperm_, cperm_,
              )

	rperm.ptr = rperm_[]
	cperm.ptr = cperm_[]

	return nothing
end 

"""
	node_count::PetscInt,sizes::Ptr{PetscInt},limit::PetscInt = MatInodeGetInodeSizes(petsclib::PetscLibType,A::AbstractPetscMat) 
Returns the inode information of a matrix with inodes

Not Collective

Input Parameter:
- `A` - the Inode matrix or matrix derived from the Inode class -- e.g., `MATSEQAIJ`

Output Parameters:
- `node_count` - no of inodes present in the matrix.
- `sizes`      - an array of size `node_count`, with the sizes of each inode.
- `limit`      - the max size used to generate the inodes.

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatGetInfo()`

# External Links
$(_doc_external("Mat/MatInodeGetInodeSizes"))
"""
function MatInodeGetInodeSizes(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatInodeGetInodeSizes(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	node_count_ = Ref{$PetscInt}()
	sizes_ = Ref{Ptr{$PetscInt}}()
	limit_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatInodeGetInodeSizes, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{$PetscInt}),
               A, node_count_, sizes_, limit_,
              )

	node_count = node_count_[]
	sizes = sizes_[]
	limit = limit_[]

	return node_count,sizes,limit
end 

"""
	MatInterpolate(petsclib::PetscLibType,A::AbstractPetscMat, x::AbstractPetscVec, y::AbstractPetscVec) 
y = A*x or A^T*x depending on the shape of
the matrix

Neighbor-wise Collective

Input Parameters:
- `A` - the matrix
- `x` - the vector to be interpolated

Output Parameter:
- `y` - the resulting vector

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatMultAdd()`, `MatMultTransposeAdd()`, `MatRestrict()`, `PCMG`

# External Links
$(_doc_external("Mat/MatInterpolate"))
"""
function MatInterpolate(petsclib::PetscLibType, A::AbstractPetscMat, x::AbstractPetscVec, y::AbstractPetscVec) end

@for_petsc function MatInterpolate(petsclib::$UnionPetscLib, A::AbstractPetscMat, x::AbstractPetscVec, y::AbstractPetscVec )

    @chk ccall(
               (:MatInterpolate, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               A, x, y,
              )


	return nothing
end 

"""
	MatInterpolateAdd(petsclib::PetscLibType,A::AbstractPetscMat, x::AbstractPetscVec, y::AbstractPetscVec, w::AbstractPetscVec) 
w = y + A*x or A^T*x depending on the shape of
the matrix

Neighbor-wise Collective

Input Parameters:
- `A` - the matrix
- `x` - the vector to be multiplied by the interpolation operator
- `y` - the vector to be added to the result

Output Parameter:
- `w` - the resulting vector

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatMultAdd()`, `MatMultTransposeAdd()`, `MatRestrict()`, `PCMG`

# External Links
$(_doc_external("Mat/MatInterpolateAdd"))
"""
function MatInterpolateAdd(petsclib::PetscLibType, A::AbstractPetscMat, x::AbstractPetscVec, y::AbstractPetscVec, w::AbstractPetscVec) end

@for_petsc function MatInterpolateAdd(petsclib::$UnionPetscLib, A::AbstractPetscMat, x::AbstractPetscVec, y::AbstractPetscVec, w::AbstractPetscVec )

    @chk ccall(
               (:MatInterpolateAdd, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec, CVec),
               A, x, y, w,
              )


	return nothing
end 

"""
	values::Ptr{PetscScalar} = MatInvertBlockDiagonal(petsclib::PetscLibType,mat::AbstractPetscMat) 
Inverts the block diagonal entries.

Collective; No Fortran Support

Input Parameter:
- `mat` - the matrix

Output Parameter:
- `values` - the block inverses in column major order (FORTRAN-like)

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatInvertVariableBlockEnvelope()`, `MatInvertBlockDiagonalMat()`

# External Links
$(_doc_external("Mat/MatInvertBlockDiagonal"))
"""
function MatInvertBlockDiagonal(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatInvertBlockDiagonal(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	values_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatInvertBlockDiagonal, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               mat, values_,
              )

	values = values_[]

	return values
end 

"""
	MatInvertBlockDiagonalMat(petsclib::PetscLibType,A::AbstractPetscMat, C::AbstractPetscMat) 
set the values of matrix C to be the inverted block diagonal of matrix A

Collective

Input Parameters:
- `A` - the matrix
- `C` - matrix with inverted block diagonal of `A`.  This matrix should be created and may have its type set.

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatInvertBlockDiagonal()`

# External Links
$(_doc_external("Mat/MatInvertBlockDiagonalMat"))
"""
function MatInvertBlockDiagonalMat(petsclib::PetscLibType, A::AbstractPetscMat, C::AbstractPetscMat) end

@for_petsc function MatInvertBlockDiagonalMat(petsclib::$UnionPetscLib, A::AbstractPetscMat, C::AbstractPetscMat )

    @chk ccall(
               (:MatInvertBlockDiagonalMat, $petsc_library),
               PetscErrorCode,
               (CMat, CMat),
               A, C,
              )


	return nothing
end 

"""
	MatInvertVariableBlockDiagonal(petsclib::PetscLibType,mat::AbstractPetscMat, nblocks::PetscInt, bsizes::Vector{PetscInt}, values::Vector{PetscScalar}) 
Inverts the point block diagonal entries.

Collective; No Fortran Support

Input Parameters:
- `mat`     - the matrix
- `nblocks` - the number of blocks on the process, set with `MatSetVariableBlockSizes()`
- `bsizes`  - the size of each block on the process, set with `MatSetVariableBlockSizes()`

Output Parameter:
- `values` - the block inverses in column major order (FORTRAN-like)

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatInvertBlockDiagonal()`, `MatSetVariableBlockSizes()`, `MatInvertVariableBlockEnvelope()`

# External Links
$(_doc_external("Mat/MatInvertVariableBlockDiagonal"))
"""
function MatInvertVariableBlockDiagonal(petsclib::PetscLibType, mat::AbstractPetscMat, nblocks::PetscInt, bsizes::Vector{PetscInt}, values::Vector{PetscScalar}) end

@for_petsc function MatInvertVariableBlockDiagonal(petsclib::$UnionPetscLib, mat::AbstractPetscMat, nblocks::$PetscInt, bsizes::Vector{$PetscInt}, values::Vector{$PetscScalar} )

    @chk ccall(
               (:MatInvertVariableBlockDiagonal, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}),
               mat, nblocks, bsizes, values,
              )


	return nothing
end 

"""
	C::PetscMat = MatInvertVariableBlockEnvelope(petsclib::PetscLibType,A::AbstractPetscMat, reuse::MatReuse) 
set matrix C to be the inverted block diagonal of matrix A

Collective

Input Parameters:
- `A`     - the matrix
- `reuse` - indicates if the `C` matrix was obtained from a previous call to this routine

Output Parameter:
- `C` - matrix with inverted block diagonal of `A`

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatInvertBlockDiagonal()`, `MatComputeBlockDiagonal()`

# External Links
$(_doc_external("Mat/MatInvertVariableBlockEnvelope"))
"""
function MatInvertVariableBlockEnvelope(petsclib::PetscLibType, A::AbstractPetscMat, reuse::MatReuse) end

@for_petsc function MatInvertVariableBlockEnvelope(petsclib::$UnionPetscLib, A::AbstractPetscMat, reuse::MatReuse )
	C_ = Ref{CMat}()

    @chk ccall(
               (:MatInvertVariableBlockEnvelope, $petsc_library),
               PetscErrorCode,
               (CMat, MatReuse, Ptr{CMat}),
               A, reuse, C_,
              )

	C = PetscMat(C_[], petsclib)

	return C
end 

"""
	flg::PetscBool = MatIsHermitian(petsclib::PetscLibType,A::AbstractPetscMat, tol::PetscReal) 
Test whether a matrix is Hermitian

Collective

Input Parameters:
- `A`   - the matrix to test
- `tol` - difference between value and its transpose less than this amount counts as equal (use 0.0 for exact Hermitian)

Output Parameter:
- `flg` - the result

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatTranspose()`, `MatIsTranspose()`, `MatIsHermitianKnown()`, `MatIsStructurallySymmetric()`, `MatSetOption()`,
`MatIsSymmetricKnown()`, `MatIsSymmetric()`, `MAT_HERMITIAN`, `MAT_SYMMETRY_ETERNAL`

# External Links
$(_doc_external("Mat/MatIsHermitian"))
"""
function MatIsHermitian(petsclib::PetscLibType, A::AbstractPetscMat, tol::PetscReal) end

@for_petsc function MatIsHermitian(petsclib::$UnionPetscLib, A::AbstractPetscMat, tol::$PetscReal )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatIsHermitian, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscReal, Ptr{PetscBool}),
               A, tol, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	set::PetscBool,flg::PetscBool = MatIsHermitianKnown(petsclib::PetscLibType,A::AbstractPetscMat) 
Checks if a matrix knows if it is Hermitian or not and its Hermitian state

Not Collective

Input Parameter:
- `A` - the matrix to check

Output Parameters:
- `set` - `PETSC_TRUE` if the matrix knows its Hermitian state (this tells you if the next flag is valid)
- `flg` - the result (only valid if set is `PETSC_TRUE`)

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MAT_SYMMETRY_ETERNAL`, `MAT_HERMITIAN`, `MatTranspose()`, `MatIsTranspose()`, `MatIsHermitian()`, `MatIsStructurallySymmetric()`, `MatSetOption()`, `MatIsSymmetric()`

# External Links
$(_doc_external("Mat/MatIsHermitianKnown"))
"""
function MatIsHermitianKnown(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatIsHermitianKnown(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	set_ = Ref{PetscBool}()
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatIsHermitianKnown, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscBool}, Ptr{PetscBool}),
               A, set_, flg_,
              )

	set = set_[]
	flg = flg_[]

	return set,flg
end 

"""
	flg::PetscBool = MatIsHermitianTranspose(petsclib::PetscLibType,A::AbstractPetscMat, B::AbstractPetscMat, tol::PetscReal) 
Test whether a matrix is another one's Hermitian transpose,

Collective

Input Parameters:
- `A`   - the matrix to test
- `B`   - the matrix to test against, this can equal the first parameter
- `tol` - tolerance, differences between entries smaller than this are counted as zero

Output Parameter:
- `flg` - the result

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatTranspose()`, `MatIsSymmetric()`, `MatIsHermitian()`, `MatIsTranspose()`

# External Links
$(_doc_external("Mat/MatIsHermitianTranspose"))
"""
function MatIsHermitianTranspose(petsclib::PetscLibType, A::AbstractPetscMat, B::AbstractPetscMat, tol::PetscReal) end

@for_petsc function MatIsHermitianTranspose(petsclib::$UnionPetscLib, A::AbstractPetscMat, B::AbstractPetscMat, tol::$PetscReal )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatIsHermitianTranspose, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, $PetscReal, Ptr{PetscBool}),
               A, B, tol, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	flg::PetscBool = MatIsLinear(petsclib::PetscLibType,A::AbstractPetscMat, n::PetscInt) 
Check if a shell matrix `A` is a linear operator.

Collective

Input Parameters:
- `A` - the shell matrix
- `n` - number of random vectors to be tested

Output Parameter:
- `flg` - `PETSC_TRUE` if the shell matrix is linear; `PETSC_FALSE` otherwise.

Level: intermediate

-seealso: `Mat`, `MatMatMultEqual()`, `MatMultEqual()`, `MatMultAddEqual()`, `MatMultTransposeEqual()`

# External Links
$(_doc_external("Mat/MatIsLinear"))
"""
function MatIsLinear(petsclib::PetscLibType, A::AbstractPetscMat, n::PetscInt) end

@for_petsc function MatIsLinear(petsclib::$UnionPetscLib, A::AbstractPetscMat, n::$PetscInt )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatIsLinear, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{PetscBool}),
               A, n, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	set::PetscBool,flg::PetscBool = MatIsSPDKnown(petsclib::PetscLibType,A::AbstractPetscMat) 
Checks if a matrix knows if it is symmetric positive definite or not and its symmetric positive definite state

Not Collective

Input Parameter:
- `A` - the matrix to check

Output Parameters:
- `set` - `PETSC_TRUE` if the matrix knows its symmetric positive definite state (this tells you if the next flag is valid)
- `flg` - the result (only valid if set is `PETSC_TRUE`)

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MAT_SPD_ETERNAL`, `MAT_SPD`, `MatTranspose()`, `MatIsTranspose()`, `MatIsHermitian()`, `MatIsStructurallySymmetric()`, `MatSetOption()`, `MatIsSymmetric()`, `MatIsHermitianKnown()`

# External Links
$(_doc_external("Mat/MatIsSPDKnown"))
"""
function MatIsSPDKnown(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatIsSPDKnown(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	set_ = Ref{PetscBool}()
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatIsSPDKnown, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscBool}, Ptr{PetscBool}),
               A, set_, flg_,
              )

	set = set_[]
	flg = flg_[]

	return set,flg
end 

"""
	flg::PetscBool = MatIsShell(petsclib::PetscLibType,mat::AbstractPetscMat) 
Inquires if a matrix is derived from `MATSHELL`

Input Parameter:
- `mat` - the matrix

Output Parameter:
- `flg` - the Boolean value

Level: developer

-seealso: [](ch_matrices), `Mat`, `MATSHELL`, `MATMFFD`, `MatCreateShell()`, `MATTRANSPOSEVIRTUAL`, `MATSCHURCOMPLEMENT`

# External Links
$(_doc_external("Mat/MatIsShell"))
"""
function MatIsShell(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatIsShell(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatIsShell, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscBool}),
               mat, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	flg::PetscBool = MatIsStructurallySymmetric(petsclib::PetscLibType,A::AbstractPetscMat) 
Test whether a matrix is structurally symmetric

Collective

Input Parameter:
- `A` - the matrix to test

Output Parameter:
- `flg` - the result

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MAT_STRUCTURALLY_SYMMETRIC`, `MAT_STRUCTURAL_SYMMETRY_ETERNAL`, `MatTranspose()`, `MatIsTranspose()`, `MatIsHermitian()`, `MatIsSymmetric()`, `MatSetOption()`, `MatIsStructurallySymmetricKnown()`

# External Links
$(_doc_external("Mat/MatIsStructurallySymmetric"))
"""
function MatIsStructurallySymmetric(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatIsStructurallySymmetric(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatIsStructurallySymmetric, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscBool}),
               A, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	set::PetscBool,flg::PetscBool = MatIsStructurallySymmetricKnown(petsclib::PetscLibType,A::AbstractPetscMat) 
Checks if a matrix knows if it is structurally symmetric or not and its structurally symmetric state

Not Collective

Input Parameter:
- `A` - the matrix to check

Output Parameters:
- `set` - PETSC_TRUE if the matrix knows its structurally symmetric state (this tells you if the next flag is valid)
- `flg` - the result (only valid if set is PETSC_TRUE)

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MAT_STRUCTURALLY_SYMMETRIC`, `MatTranspose()`, `MatIsTranspose()`, `MatIsHermitian()`, `MatIsStructurallySymmetric()`, `MatSetOption()`, `MatIsSymmetric()`, `MatIsHermitianKnown()`

# External Links
$(_doc_external("Mat/MatIsStructurallySymmetricKnown"))
"""
function MatIsStructurallySymmetricKnown(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatIsStructurallySymmetricKnown(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	set_ = Ref{PetscBool}()
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatIsStructurallySymmetricKnown, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscBool}, Ptr{PetscBool}),
               A, set_, flg_,
              )

	set = set_[]
	flg = flg_[]

	return set,flg
end 

"""
	flg::PetscBool = MatIsSymmetric(petsclib::PetscLibType,A::AbstractPetscMat, tol::PetscReal) 
Test whether a matrix is symmetric

Collective

Input Parameters:
- `A`   - the matrix to test
- `tol` - difference between value and its transpose less than this amount counts as equal (use 0.0 for exact transpose)

Output Parameter:
- `flg` - the result

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatTranspose()`, `MatIsTranspose()`, `MatIsHermitian()`, `MatIsStructurallySymmetric()`, `MatSetOption()`, `MatIsSymmetricKnown()`,
`MAT_SYMMETRIC`, `MAT_SYMMETRY_ETERNAL`

# External Links
$(_doc_external("Mat/MatIsSymmetric"))
"""
function MatIsSymmetric(petsclib::PetscLibType, A::AbstractPetscMat, tol::PetscReal) end

@for_petsc function MatIsSymmetric(petsclib::$UnionPetscLib, A::AbstractPetscMat, tol::$PetscReal )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatIsSymmetric, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscReal, Ptr{PetscBool}),
               A, tol, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	set::PetscBool,flg::PetscBool = MatIsSymmetricKnown(petsclib::PetscLibType,A::AbstractPetscMat) 
Checks if a matrix knows if it is symmetric or not and its symmetric state

Not Collective

Input Parameter:
- `A` - the matrix to check

Output Parameters:
- `set` - `PETSC_TRUE` if the matrix knows its symmetry state (this tells you if the next flag is valid)
- `flg` - the result (only valid if set is `PETSC_TRUE`)

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MAT_SYMMETRY_ETERNAL`, `MatTranspose()`, `MatIsTranspose()`, `MatIsHermitian()`, `MatIsStructurallySymmetric()`, `MatSetOption()`, `MatIsSymmetric()`, `MatIsHermitianKnown()`

# External Links
$(_doc_external("Mat/MatIsSymmetricKnown"))
"""
function MatIsSymmetricKnown(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatIsSymmetricKnown(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	set_ = Ref{PetscBool}()
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatIsSymmetricKnown, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscBool}, Ptr{PetscBool}),
               A, set_, flg_,
              )

	set = set_[]
	flg = flg_[]

	return set,flg
end 

"""
	flg::PetscBool = MatIsTranspose(petsclib::PetscLibType,A::AbstractPetscMat, B::AbstractPetscMat, tol::PetscReal) 
Test whether a matrix is another one's transpose,
or its own, in which case it tests symmetry.

Collective

Input Parameters:
- `A`   - the matrix to test
- `B`   - the matrix to test against, this can equal the first parameter
- `tol` - tolerance, differences between entries smaller than this are counted as zero

Output Parameter:
- `flg` - the result

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatTranspose()`, `MatIsSymmetric()`, `MatIsHermitian()`

# External Links
$(_doc_external("Mat/MatIsTranspose"))
"""
function MatIsTranspose(petsclib::PetscLibType, A::AbstractPetscMat, B::AbstractPetscMat, tol::PetscReal) end

@for_petsc function MatIsTranspose(petsclib::$UnionPetscLib, A::AbstractPetscMat, B::AbstractPetscMat, tol::$PetscReal )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatIsTranspose, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, $PetscReal, Ptr{PetscBool}),
               A, B, tol, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	B::PetscMat = MatKAIJGetAIJ(petsclib::PetscLibType,A::AbstractPetscMat) 
Get the `MATAIJ` matrix describing the blockwise action of the `MATKAIJ` matrix

Not Collective, but if the `MATKAIJ` matrix is parallel, the `MATAIJ` matrix is also parallel

Input Parameter:
- `A` - the `MATKAIJ` matrix

Output Parameter:
- `B` - the `MATAIJ` matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatCreateKAIJ()`, `MATKAIJ`, `MATAIJ`

# External Links
$(_doc_external("Mat/MatKAIJGetAIJ"))
"""
function MatKAIJGetAIJ(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatKAIJGetAIJ(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	B_ = Ref{CMat}()

    @chk ccall(
               (:MatKAIJGetAIJ, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, B_,
              )

	B = PetscMat(B_[], petsclib)

	return B
end 

"""
	m::PetscInt,n::PetscInt,S::Ptr{PetscScalar} = MatKAIJGetS(petsclib::PetscLibType,A::AbstractPetscMat) 
Get the `S` matrix describing the shift action of the `MATKAIJ` matrix

Not Collective; the entire `S` is stored and returned independently on all processes.

Input Parameter:
- `A` - the `MATKAIJ` matrix

Output Parameters:
- `m` - the number of rows in `S`
- `n` - the number of columns in `S`
- `S` - the S matrix, in form of a scalar array in column-major format

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATKAIJ`, `MatCreateKAIJ()`, `MatGetBlockSizes()`

# External Links
$(_doc_external("Mat/MatKAIJGetS"))
"""
function MatKAIJGetS(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatKAIJGetS(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	m_ = Ref{$PetscInt}()
	n_ = Ref{$PetscInt}()
	S_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatKAIJGetS, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscScalar}}),
               A, m_, n_, S_,
              )

	m = m_[]
	n = n_[]
	S = S_[]

	return m,n,S
end 

"""
	m::PetscInt,n::PetscInt,S::Ptr{PetscScalar} = MatKAIJGetSRead(petsclib::PetscLibType,A::AbstractPetscMat) 
Get a read

Not Collective; the entire `S` is stored and returned independently on all processes.

Input Parameter:
- `A` - the `MATKAIJ` matrix

Output Parameters:
- `m` - the number of rows in `S`
- `n` - the number of columns in `S`
- `S` - the S matrix, in form of a scalar array in column-major format

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATKAIJ`, `MatCreateKAIJ()`, `MatGetBlockSizes()`

# External Links
$(_doc_external("Mat/MatKAIJGetSRead"))
"""
function MatKAIJGetSRead(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatKAIJGetSRead(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	m_ = Ref{$PetscInt}()
	n_ = Ref{$PetscInt}()
	S_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatKAIJGetSRead, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscScalar}}),
               A, m_, n_, S_,
              )

	m = m_[]
	n = n_[]
	S = S_[]

	return m,n,S
end 

"""
	identity::PetscBool = MatKAIJGetScaledIdentity(petsclib::PetscLibType,A::AbstractPetscMat) 
Check if both `S` and `T` are scaled identities.

Logically Collective.

Input Parameter:
- `A` - the `MATKAIJ` matrix

Output Parameter:
- `identity` - the Boolean value

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATKAIJ`, `MatKAIJGetS()`, `MatKAIJGetT()`

# External Links
$(_doc_external("Mat/MatKAIJGetScaledIdentity"))
"""
function MatKAIJGetScaledIdentity(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatKAIJGetScaledIdentity(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	identity_ = Ref{PetscBool}()

    @chk ccall(
               (:MatKAIJGetScaledIdentity, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscBool}),
               A, identity_,
              )

	identity = identity_[]

	return identity
end 

"""
	m::PetscInt,n::PetscInt,T::Ptr{PetscScalar} = MatKAIJGetT(petsclib::PetscLibType,A::AbstractPetscMat) 
Get the transformation matrix `T` associated with the `MATKAIJ` matrix

Not Collective; the entire `T` is stored and returned independently on all processes

Input Parameter:
- `A` - the `MATKAIJ` matrix

Output Parameters:
- `m` - the number of rows in `T`
- `n` - the number of columns in `T`
- `T` - the T matrix, in form of a scalar array in column-major format

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATKAIJ`, `MatCreateKAIJ()`, `MatGetBlockSizes()`

# External Links
$(_doc_external("Mat/MatKAIJGetT"))
"""
function MatKAIJGetT(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatKAIJGetT(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	m_ = Ref{$PetscInt}()
	n_ = Ref{$PetscInt}()
	T_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatKAIJGetT, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscScalar}}),
               A, m_, n_, T_,
              )

	m = m_[]
	n = n_[]
	T = T_[]

	return m,n,T
end 

"""
	m::PetscInt,n::PetscInt,T::Ptr{PetscScalar} = MatKAIJGetTRead(petsclib::PetscLibType,A::AbstractPetscMat) 
Get a read

Not Collective; the entire `T` is stored and returned independently on all processes

Input Parameter:
- `A` - the `MATKAIJ` matrix

Output Parameters:
- `m` - the number of rows in `T`
- `n` - the number of columns in `T`
- `T` - the T matrix, in form of a scalar array in column-major format

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATKAIJ`, `MatCreateKAIJ()`, `MatGetBlockSizes()`

# External Links
$(_doc_external("Mat/MatKAIJGetTRead"))
"""
function MatKAIJGetTRead(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatKAIJGetTRead(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	m_ = Ref{$PetscInt}()
	n_ = Ref{$PetscInt}()
	T_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatKAIJGetTRead, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscScalar}}),
               A, m_, n_, T_,
              )

	m = m_[]
	n = n_[]
	T = T_[]

	return m,n,T
end 

"""
	MatKAIJRestoreS(petsclib::PetscLibType,A::AbstractPetscMat, S::AbstractArray{PetscScalar}) 
Restore array obtained with `MatKAIJGetS()`

Not Collective

Input Parameters:
- `A` - the `MATKAIJ` matrix
- `S` - location of pointer to array obtained with `MatKAIJGetS()`

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATKAIJ`, `MatKAIJGetS()`, `MatKAIJGetSRead()`, `MatKAIJRestoreSRead()`

# External Links
$(_doc_external("Mat/MatKAIJRestoreS"))
"""
function MatKAIJRestoreS(petsclib::PetscLibType, A::AbstractPetscMat, S::AbstractArray{PetscScalar}) end

@for_petsc function MatKAIJRestoreS(petsclib::$UnionPetscLib, A::AbstractPetscMat, S::AbstractArray{$PetscScalar} )
	S_ = Ref(pointer(S))

    @chk ccall(
               (:MatKAIJRestoreS, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, S_,
              )


	return nothing
end 

"""
	MatKAIJRestoreSRead(petsclib::PetscLibType,A::AbstractPetscMat, S::AbstractArray{PetscScalar}) 
Restore array obtained with `MatKAIJGetSRead()`

Not Collective

Input Parameters:
- `A` - the `MATKAIJ` matrix
- `S` - location of pointer to array obtained with `MatKAIJGetS()`

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATKAIJ`, `MatKAIJGetS()`, `MatKAIJGetSRead()`

# External Links
$(_doc_external("Mat/MatKAIJRestoreSRead"))
"""
function MatKAIJRestoreSRead(petsclib::PetscLibType, A::AbstractPetscMat, S::AbstractArray{PetscScalar}) end

@for_petsc function MatKAIJRestoreSRead(petsclib::$UnionPetscLib, A::AbstractPetscMat, S::AbstractArray{$PetscScalar} )
	S_ = Ref(pointer(S))

    @chk ccall(
               (:MatKAIJRestoreSRead, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, S_,
              )


	return nothing
end 

"""
	MatKAIJRestoreT(petsclib::PetscLibType,A::AbstractPetscMat, T::AbstractArray{PetscScalar}) 
Restore array obtained with `MatKAIJGetT()`

Not Collective

Input Parameters:
- `A` - the `MATKAIJ` matrix
- `T` - location of pointer to array obtained with `MatKAIJGetS()`

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATKAIJ`, `MatKAIJGetT()`, `MatKAIJGetTRead()`, `MatKAIJRestoreTRead()`

# External Links
$(_doc_external("Mat/MatKAIJRestoreT"))
"""
function MatKAIJRestoreT(petsclib::PetscLibType, A::AbstractPetscMat, T::AbstractArray{PetscScalar}) end

@for_petsc function MatKAIJRestoreT(petsclib::$UnionPetscLib, A::AbstractPetscMat, T::AbstractArray{$PetscScalar} )
	T_ = Ref(pointer(T))

    @chk ccall(
               (:MatKAIJRestoreT, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, T_,
              )


	return nothing
end 

"""
	MatKAIJRestoreTRead(petsclib::PetscLibType,A::AbstractPetscMat, T::AbstractArray{PetscScalar}) 
Restore array obtained with `MatKAIJGetTRead()`

Not Collective

Input Parameters:
- `A` - the `MATKAIJ` matrix
- `T` - location of pointer to array obtained with `MatKAIJGetS()`

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATKAIJ`, `MatKAIJGetT()`, `MatKAIJGetTRead()`

# External Links
$(_doc_external("Mat/MatKAIJRestoreTRead"))
"""
function MatKAIJRestoreTRead(petsclib::PetscLibType, A::AbstractPetscMat, T::AbstractArray{PetscScalar}) end

@for_petsc function MatKAIJRestoreTRead(petsclib::$UnionPetscLib, A::AbstractPetscMat, T::AbstractArray{$PetscScalar} )
	T_ = Ref(pointer(T))

    @chk ccall(
               (:MatKAIJRestoreTRead, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, T_,
              )


	return nothing
end 

"""
	MatKAIJSetAIJ(petsclib::PetscLibType,A::AbstractPetscMat, B::AbstractPetscMat) 
Set the `MATAIJ` matrix describing the blockwise action of the `MATKAIJ` matrix

Logically Collective; if the `MATAIJ` matrix is parallel, the `MATKAIJ` matrix is also parallel

Input Parameters:
- `A` - the `MATKAIJ` matrix
- `B` - the `MATAIJ` matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATKAIJ`, `MatKAIJGetAIJ()`, `MatKAIJSetS()`, `MatKAIJSetT()`

# External Links
$(_doc_external("Mat/MatKAIJSetAIJ"))
"""
function MatKAIJSetAIJ(petsclib::PetscLibType, A::AbstractPetscMat, B::AbstractPetscMat) end

@for_petsc function MatKAIJSetAIJ(petsclib::$UnionPetscLib, A::AbstractPetscMat, B::AbstractPetscMat )

    @chk ccall(
               (:MatKAIJSetAIJ, $petsc_library),
               PetscErrorCode,
               (CMat, CMat),
               A, B,
              )


	return nothing
end 

"""
	MatKAIJSetS(petsclib::PetscLibType,A::AbstractPetscMat, p::PetscInt, q::PetscInt, S::Vector{PetscScalar}) 
Set the `S` matrix describing the shift action of the `MATKAIJ` matrix

Logically Collective; the entire `S` is stored independently on all processes.

Input Parameters:
- `A` - the `MATKAIJ` matrix
- `p` - the number of rows in `S`
- `q` - the number of columns in `S`
- `S` - the S matrix, in form of a scalar array in column-major format

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATKAIJ`, `MatKAIJGetS()`, `MatKAIJSetT()`, `MatKAIJSetAIJ()`

# External Links
$(_doc_external("Mat/MatKAIJSetS"))
"""
function MatKAIJSetS(petsclib::PetscLibType, A::AbstractPetscMat, p::PetscInt, q::PetscInt, S::Vector{PetscScalar}) end

@for_petsc function MatKAIJSetS(petsclib::$UnionPetscLib, A::AbstractPetscMat, p::$PetscInt, q::$PetscInt, S::Vector{$PetscScalar} )

    @chk ccall(
               (:MatKAIJSetS, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, Ptr{$PetscScalar}),
               A, p, q, S,
              )


	return nothing
end 

"""
	MatKAIJSetT(petsclib::PetscLibType,A::AbstractPetscMat, p::PetscInt, q::PetscInt, T::Vector{PetscScalar}) 
Set the transformation matrix `T` associated with the `MATKAIJ` matrix

Logically Collective; the entire `T` is stored independently on all processes.

Input Parameters:
- `A` - the `MATKAIJ` matrix
- `p` - the number of rows in `S`
- `q` - the number of columns in `S`
- `T` - the `T` matrix, in form of a scalar array in column-major format

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATKAIJ`, `MatKAIJGetT()`, `MatKAIJSetS()`, `MatKAIJSetAIJ()`

# External Links
$(_doc_external("Mat/MatKAIJSetT"))
"""
function MatKAIJSetT(petsclib::PetscLibType, A::AbstractPetscMat, p::PetscInt, q::PetscInt, T::Vector{PetscScalar}) end

@for_petsc function MatKAIJSetT(petsclib::$UnionPetscLib, A::AbstractPetscMat, p::$PetscInt, q::$PetscInt, T::Vector{$PetscScalar} )

    @chk ccall(
               (:MatKAIJSetT, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, Ptr{$PetscScalar}),
               A, p, q, T,
              )


	return nothing
end 

"""
	MatLMVMAllocate(petsclib::PetscLibType,B::AbstractPetscMat, X::AbstractPetscVec, F::AbstractPetscVec) 
Produces all necessary common memory for
LMVM approximations based on the solution and function vectors
provided.

Input Parameters:
- `B` - A `MATLMVM` matrix
- `X` - Solution vector
- `F` - Function vector

Level: intermediate

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMReset()`, `MatLMVMUpdate()`

# External Links
$(_doc_external("KSP/MatLMVMAllocate"))
"""
function MatLMVMAllocate(petsclib::PetscLibType, B::AbstractPetscMat, X::AbstractPetscVec, F::AbstractPetscVec) end

@for_petsc function MatLMVMAllocate(petsclib::$UnionPetscLib, B::AbstractPetscMat, X::AbstractPetscVec, F::AbstractPetscVec )

    @chk ccall(
               (:MatLMVMAllocate, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               B, X, F,
              )


	return nothing
end 

"""
	MatLMVMApplyJ0Fwd(petsclib::PetscLibType,B::AbstractPetscMat, X::AbstractPetscVec, Y::AbstractPetscVec) 
Applies an approximation of the forward
matrix-vector product with the initial Jacobian.

Input Parameters:
- `B` - A `MATLMVM` matrix
- `X` - vector to multiply with J0

Output Parameter:
- `Y` - resulting vector for the operation

Level: advanced

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMSetJ0()`, `MatLMVMSetJ0Scale()`, `MatLMVMSetJ0ScaleDiag()`,
`MatLMVMSetJ0PC()`, `MatLMVMSetJ0KSP()`, `MatLMVMApplyJ0Inv()`

# External Links
$(_doc_external("KSP/MatLMVMApplyJ0Fwd"))
"""
function MatLMVMApplyJ0Fwd(petsclib::PetscLibType, B::AbstractPetscMat, X::AbstractPetscVec, Y::AbstractPetscVec) end

@for_petsc function MatLMVMApplyJ0Fwd(petsclib::$UnionPetscLib, B::AbstractPetscMat, X::AbstractPetscVec, Y::AbstractPetscVec )

    @chk ccall(
               (:MatLMVMApplyJ0Fwd, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               B, X, Y,
              )


	return nothing
end 

"""
	MatLMVMApplyJ0Inv(petsclib::PetscLibType,B::AbstractPetscMat, X::AbstractPetscVec, Y::AbstractPetscVec) 
Applies some estimation of the initial Jacobian
inverse to the given vector.

Input Parameters:
- `B` - A `MATLMVM` matrix
- `X` - vector to "multiply" with J0^{-1}

Output Parameter:
- `Y` - resulting vector for the operation

Level: advanced

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMSetJ0()`, `MatLMVMSetJ0Scale()`, `MatLMVMSetJ0ScaleDiag()`,
`MatLMVMSetJ0PC()`, `MatLMVMSetJ0KSP()`, `MatLMVMApplyJ0Fwd()`

# External Links
$(_doc_external("KSP/MatLMVMApplyJ0Inv"))
"""
function MatLMVMApplyJ0Inv(petsclib::PetscLibType, B::AbstractPetscMat, X::AbstractPetscVec, Y::AbstractPetscVec) end

@for_petsc function MatLMVMApplyJ0Inv(petsclib::$UnionPetscLib, B::AbstractPetscMat, X::AbstractPetscVec, Y::AbstractPetscVec )

    @chk ccall(
               (:MatLMVMApplyJ0Inv, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               B, X, Y,
              )


	return nothing
end 

"""
	MatLMVMClearJ0(petsclib::PetscLibType,B::AbstractPetscMat) 
Removes all definitions of J0 and reverts to
an identity matrix (scale = 1.0).

Input Parameter:
- `B` - A `MATLMVM` matrix

Level: advanced

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMSetJ0()`

# External Links
$(_doc_external("KSP/MatLMVMClearJ0"))
"""
function MatLMVMClearJ0(petsclib::PetscLibType, B::AbstractPetscMat) end

@for_petsc function MatLMVMClearJ0(petsclib::$UnionPetscLib, B::AbstractPetscMat )

    @chk ccall(
               (:MatLMVMClearJ0, $petsc_library),
               PetscErrorCode,
               (CMat,),
               B,
              )


	return nothing
end 

"""
	MatLMVMDenseSetType(petsclib::PetscLibType,B::AbstractPetscMat, type::MatLMVMDenseType) 
Sets the memory storage type for dense `MATLMVM`

Input Parameters:
- `B`    - the `MATLMVM` matrix
- `type` - scale type, see `MatLMVMDenseSetType`

Options Database Keys:
- `-mat_lqn_type   <reorder,inplace>` - set the strategy
- `-mat_lbfgs_type <reorder,inplace>` - set the strategy
- `-mat_ldfp_type  <reorder,inplace>` - set the strategy

Level: intermediate

MatLMVMDenseTypes:
- `MAT_LMVM_DENSE_REORDER`   - reorders memory to minimize kernel launch
- `MAT_LMVM_DENSE_INPLACE`   - launches kernel inplace to minimize memory movement

-seealso: [](ch_ksp), `MATLMVMDQN`, `MATLMVMDBFGS`, `MATLMVMDDFP`, `MatLMVMDenseType`

# External Links
$(_doc_external("KSP/MatLMVMDenseSetType"))
"""
function MatLMVMDenseSetType(petsclib::PetscLibType, B::AbstractPetscMat, type::MatLMVMDenseType) end

@for_petsc function MatLMVMDenseSetType(petsclib::$UnionPetscLib, B::AbstractPetscMat, type::MatLMVMDenseType )

    @chk ccall(
               (:MatLMVMDenseSetType, $petsc_library),
               PetscErrorCode,
               (CMat, MatLMVMDenseType),
               B, type,
              )


	return nothing
end 

"""
	hist_size::PetscInt = MatLMVMGetHistorySize(petsclib::PetscLibType,B::AbstractPetscMat) 

# External Links
$(_doc_external("KSP/MatLMVMGetHistorySize"))
"""
function MatLMVMGetHistorySize(petsclib::PetscLibType, B::AbstractPetscMat) end

@for_petsc function MatLMVMGetHistorySize(petsclib::$UnionPetscLib, B::AbstractPetscMat )
	hist_size_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatLMVMGetHistorySize, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}),
               B, hist_size_,
              )

	hist_size = hist_size_[]

	return hist_size
end 

"""
	J0::PetscMat = MatLMVMGetJ0(petsclib::PetscLibType,B::AbstractPetscMat) 
Returns a pointer to the internal `J0` matrix.

Input Parameter:
- `B` - A `MATLMVM` matrix

Output Parameter:
- `J0` - `Mat` object for defining the initial Jacobian

Level: advanced

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMSetJ0()`

# External Links
$(_doc_external("KSP/MatLMVMGetJ0"))
"""
function MatLMVMGetJ0(petsclib::PetscLibType, B::AbstractPetscMat) end

@for_petsc function MatLMVMGetJ0(petsclib::$UnionPetscLib, B::AbstractPetscMat )
	J0_ = Ref{CMat}()

    @chk ccall(
               (:MatLMVMGetJ0, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               B, J0_,
              )

	J0 = PetscMat(J0_[], petsclib)

	return J0
end 

"""
	J0ksp::PetscKSP = MatLMVMGetJ0KSP(petsclib::PetscLibType,B::AbstractPetscMat) 
Returns a pointer to the internal `KSP` solver
associated with the initial Jacobian.

Input Parameter:
- `B` - A `MATLMVM` matrix

Output Parameter:
- `J0ksp` - `KSP` solver for defining the initial inverse-Jacobian

Level: advanced

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMSetJ0KSP()`

# External Links
$(_doc_external("KSP/MatLMVMGetJ0KSP"))
"""
function MatLMVMGetJ0KSP(petsclib::PetscLibType, B::AbstractPetscMat) end

@for_petsc function MatLMVMGetJ0KSP(petsclib::$UnionPetscLib, B::AbstractPetscMat )
	J0ksp_ = Ref{CKSP}()

    @chk ccall(
               (:MatLMVMGetJ0KSP, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CKSP}),
               B, J0ksp_,
              )

	J0ksp = PetscKSP(J0ksp_[], petsclib)

	return J0ksp
end 

"""
	J0pc::PC = MatLMVMGetJ0PC(petsclib::PetscLibType,B::AbstractPetscMat) 
Returns a pointer to the internal `PC` object
associated with the initial Jacobian.

Input Parameter:
- `B` - A `MATLMVM` matrix

Output Parameter:
- `J0pc` - `PC` object for defining the initial inverse-Jacobian

Level: advanced

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMSetJ0PC()`

# External Links
$(_doc_external("KSP/MatLMVMGetJ0PC"))
"""
function MatLMVMGetJ0PC(petsclib::PetscLibType, B::AbstractPetscMat) end

@for_petsc function MatLMVMGetJ0PC(petsclib::$UnionPetscLib, B::AbstractPetscMat )
	J0pc_ = Ref{PC}()

    @chk ccall(
               (:MatLMVMGetJ0PC, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PC}),
               B, J0pc_,
              )

	J0pc = J0pc_[]

	return J0pc
end 

"""
	x_prev::PetscVec,f_prev::PetscVec = MatLMVMGetLastUpdate(petsclib::PetscLibType,B::AbstractPetscMat) 
Get the last vectors passed to `MatLMVMUpdate()`

Not collective

Input Parameter:
- `B` - a `MatLMVM` matrix

Output Parameters:
- `x_prev` - the last solution vector
- `f_prev` - the last function vector

Level: intermediate

-seealso: [](ch_matrices), `MatLMVM`, `MatLMVMUpdate()`

# External Links
$(_doc_external("KSP/MatLMVMGetLastUpdate"))
"""
function MatLMVMGetLastUpdate(petsclib::PetscLibType, B::AbstractPetscMat) end

@for_petsc function MatLMVMGetLastUpdate(petsclib::$UnionPetscLib, B::AbstractPetscMat )
	x_prev_ = Ref{CVec}()
	f_prev_ = Ref{CVec}()

    @chk ccall(
               (:MatLMVMGetLastUpdate, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CVec}, Ptr{CVec}),
               B, x_prev_, f_prev_,
              )

	x_prev = PetscVec(x_prev_[], petsclib)
	f_prev = PetscVec(f_prev_[], petsclib)

	return x_prev,f_prev
end 

"""
	alg::MatLMVMMultAlgorithm = MatLMVMGetMultAlgorithm(petsclib::PetscLibType,B::AbstractPetscMat) 
Get the algorithm used by a `MatLMVM` for products

Not collective

Input Parameter:
- `B` - a `MatLMVM` matrix

Output Parameter:
- `alg` - one of the algorithm classes (`MAT_LMVM_MULT_RECURSIVE`, `MAT_LMVM_MULT_DENSE`, `MAT_LMVM_MULT_COMPACT_DENSE`)

Level: advanced

-seealso: [](ch_matrices), `MatLMVM`, `MatLMVMMultAlgorithm`, `MatLMVMSetMultAlgorithm()`

# External Links
$(_doc_external("KSP/MatLMVMGetMultAlgorithm"))
"""
function MatLMVMGetMultAlgorithm(petsclib::PetscLibType, B::AbstractPetscMat) end

@for_petsc function MatLMVMGetMultAlgorithm(petsclib::$UnionPetscLib, B::AbstractPetscMat )
	alg_ = Ref{MatLMVMMultAlgorithm}()

    @chk ccall(
               (:MatLMVMGetMultAlgorithm, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatLMVMMultAlgorithm}),
               B, alg_,
              )

	alg = alg_[]

	return alg
end 

"""
	nrejects::PetscInt = MatLMVMGetRejectCount(petsclib::PetscLibType,B::AbstractPetscMat) 
Returns the number of rejected updates.
The counters are reset when `MatLMVMReset()` is called.

Input Parameter:
- `B` - A `MATLMVM` matrix

Output Parameter:
- `nrejects` - number of rejected updates

Level: intermediate

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMReset()`

# External Links
$(_doc_external("KSP/MatLMVMGetRejectCount"))
"""
function MatLMVMGetRejectCount(petsclib::PetscLibType, B::AbstractPetscMat) end

@for_petsc function MatLMVMGetRejectCount(petsclib::$UnionPetscLib, B::AbstractPetscMat )
	nrejects_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatLMVMGetRejectCount, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}),
               B, nrejects_,
              )

	nrejects = nrejects_[]

	return nrejects
end 

"""
	nupdates::PetscInt = MatLMVMGetUpdateCount(petsclib::PetscLibType,B::AbstractPetscMat) 
Returns the number of accepted updates.

Input Parameter:
- `B` - A `MATLMVM` matrix

Output Parameter:
- `nupdates` - number of accepted updates

Level: intermediate

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMGetRejectCount()`, `MatLMVMReset()`

# External Links
$(_doc_external("KSP/MatLMVMGetUpdateCount"))
"""
function MatLMVMGetUpdateCount(petsclib::PetscLibType, B::AbstractPetscMat) end

@for_petsc function MatLMVMGetUpdateCount(petsclib::$UnionPetscLib, B::AbstractPetscMat )
	nupdates_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatLMVMGetUpdateCount, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}),
               B, nupdates_,
              )

	nupdates = nupdates_[]

	return nupdates
end 

"""
	flg::PetscBool = MatLMVMIsAllocated(petsclib::PetscLibType,B::AbstractPetscMat) 
Returns a boolean flag that shows whether
the necessary data structures for the underlying matrix is allocated.

Input Parameter:
- `B` - A `MATLMVM` matrix

Output Parameter:
- `flg` - `PETSC_TRUE` if allocated, `PETSC_FALSE` otherwise

Level: intermediate

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMAllocate()`, `MatLMVMReset()`

# External Links
$(_doc_external("KSP/MatLMVMIsAllocated"))
"""
function MatLMVMIsAllocated(petsclib::PetscLibType, B::AbstractPetscMat) end

@for_petsc function MatLMVMIsAllocated(petsclib::$UnionPetscLib, B::AbstractPetscMat )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatLMVMIsAllocated, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscBool}),
               B, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	MatLMVMReset(petsclib::PetscLibType,B::AbstractPetscMat, destructive::PetscBool) 
Flushes all of the accumulated updates out of
the `MATLMVM` approximation.

Input Parameters:
- `B`           - A `MATLMVM` matrix
- `destructive` - flag for enabling destruction of data structures

Level: intermediate

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMAllocate()`, `MatLMVMUpdate()`

# External Links
$(_doc_external("KSP/MatLMVMReset"))
"""
function MatLMVMReset(petsclib::PetscLibType, B::AbstractPetscMat, destructive::PetscBool) end

@for_petsc function MatLMVMReset(petsclib::$UnionPetscLib, B::AbstractPetscMat, destructive::PetscBool )

    @chk ccall(
               (:MatLMVMReset, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               B, destructive,
              )


	return nothing
end 

"""
	MatLMVMResetShift(petsclib::PetscLibType,B::AbstractPetscMat) 
Zero the shift factor for a `MATLMVM`.

Input Parameter:
- `B` - A `MATLMVM` matrix

Level: intermediate

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMAllocate()`, `MatLMVMUpdate()`

# External Links
$(_doc_external("KSP/MatLMVMResetShift"))
"""
function MatLMVMResetShift(petsclib::PetscLibType, B::AbstractPetscMat) end

@for_petsc function MatLMVMResetShift(petsclib::$UnionPetscLib, B::AbstractPetscMat )

    @chk ccall(
               (:MatLMVMResetShift, $petsc_library),
               PetscErrorCode,
               (CMat,),
               B,
              )


	return nothing
end 

"""
	MatLMVMSetHistorySize(petsclib::PetscLibType,B::AbstractPetscMat, hist_size::PetscInt) 
Set the number of past iterates to be
stored for the construction of the limited-memory quasi-Newton update.

Input Parameters:
- `B`         - A `MATLMVM` matrix
- `hist_size` - number of past iterates (default 5)

Options Database Key:
- `-mat_lmvm_hist_size <m>` - set number of past iterates

Level: beginner

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMGetUpdateCount()`

# External Links
$(_doc_external("KSP/MatLMVMSetHistorySize"))
"""
function MatLMVMSetHistorySize(petsclib::PetscLibType, B::AbstractPetscMat, hist_size::PetscInt) end

@for_petsc function MatLMVMSetHistorySize(petsclib::$UnionPetscLib, B::AbstractPetscMat, hist_size::$PetscInt )

    @chk ccall(
               (:MatLMVMSetHistorySize, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt),
               B, hist_size,
              )


	return nothing
end 

"""
	MatLMVMSetJ0(petsclib::PetscLibType,B::AbstractPetscMat, J0::AbstractPetscMat) 
Allows the user to define the initial Jacobian matrix from which the LMVM
up.

Input Parameters:
- `B`  - An LMVM-type matrix
- `J0` - The initial Jacobian matrix, will be referenced by B.

Level: advanced

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMSetJ0PC()`, `MatLMVMSetJ0KSP()`

# External Links
$(_doc_external("KSP/MatLMVMSetJ0"))
"""
function MatLMVMSetJ0(petsclib::PetscLibType, B::AbstractPetscMat, J0::AbstractPetscMat) end

@for_petsc function MatLMVMSetJ0(petsclib::$UnionPetscLib, B::AbstractPetscMat, J0::AbstractPetscMat )

    @chk ccall(
               (:MatLMVMSetJ0, $petsc_library),
               PetscErrorCode,
               (CMat, CMat),
               B, J0,
              )


	return nothing
end 

"""
	MatLMVMSetJ0Diag(petsclib::PetscLibType,B::AbstractPetscMat, V::AbstractPetscVec) 
Allows the user to define a vector
V such that J0 = diag(V).

Input Parameters:
- `B` - An LMVM-type matrix
- `V` - Vector that defines the diagonal of the initial Jacobian: values are copied, V is not referenced

Level: advanced

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMSetScale()`, `MatLMVMSetJ0()`

# External Links
$(_doc_external("KSP/MatLMVMSetJ0Diag"))
"""
function MatLMVMSetJ0Diag(petsclib::PetscLibType, B::AbstractPetscMat, V::AbstractPetscVec) end

@for_petsc function MatLMVMSetJ0Diag(petsclib::$UnionPetscLib, B::AbstractPetscMat, V::AbstractPetscVec )

    @chk ccall(
               (:MatLMVMSetJ0Diag, $petsc_library),
               PetscErrorCode,
               (CMat, CVec),
               B, V,
              )


	return nothing
end 

"""
	MatLMVMSetJ0KSP(petsclib::PetscLibType,B::AbstractPetscMat, J0ksp::AbstractPetscKSP) 
Allows the user to provide a pre
approximation.

Input Parameters:
- `B`     - A `MATLMVM` matrix
- `J0ksp` - `KSP` solver for the initial inverse-Jacobian application

Level: advanced

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMGetJ0KSP()`

# External Links
$(_doc_external("KSP/MatLMVMSetJ0KSP"))
"""
function MatLMVMSetJ0KSP(petsclib::PetscLibType, B::AbstractPetscMat, J0ksp::AbstractPetscKSP) end

@for_petsc function MatLMVMSetJ0KSP(petsclib::$UnionPetscLib, B::AbstractPetscMat, J0ksp::AbstractPetscKSP )

    @chk ccall(
               (:MatLMVMSetJ0KSP, $petsc_library),
               PetscErrorCode,
               (CMat, CKSP),
               B, J0ksp,
              )


	return nothing
end 

"""
	MatLMVMSetJ0PC(petsclib::PetscLibType,B::AbstractPetscMat, J0pc::PC) 
Allows the user to define a `PC` object that acts as the initial inverse

Input Parameters:
- `B`    - A `MATLMVM` matrix
- `J0pc` - `PC` object where `PCApply()` defines an inverse application for J0

Level: advanced

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMGetJ0PC()`

# External Links
$(_doc_external("KSP/MatLMVMSetJ0PC"))
"""
function MatLMVMSetJ0PC(petsclib::PetscLibType, B::AbstractPetscMat, J0pc::PC) end

@for_petsc function MatLMVMSetJ0PC(petsclib::$UnionPetscLib, B::AbstractPetscMat, J0pc::PC )

    @chk ccall(
               (:MatLMVMSetJ0PC, $petsc_library),
               PetscErrorCode,
               (CMat, PC),
               B, J0pc,
              )


	return nothing
end 

"""
	MatLMVMSetJ0Scale(petsclib::PetscLibType,B::AbstractPetscMat, scale::PetscReal) 
Allows the user to define a scalar value
mu such that J0 = mu*I.

Input Parameters:
- `B`     - A `MATLMVM` matrix
- `scale` - Scalar value mu that defines the initial Jacobian

Level: advanced

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMSetDiagScale()`, `MatLMVMSetJ0()`

# External Links
$(_doc_external("KSP/MatLMVMSetJ0Scale"))
"""
function MatLMVMSetJ0Scale(petsclib::PetscLibType, B::AbstractPetscMat, scale::PetscReal) end

@for_petsc function MatLMVMSetJ0Scale(petsclib::$UnionPetscLib, B::AbstractPetscMat, scale::$PetscReal )

    @chk ccall(
               (:MatLMVMSetJ0Scale, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscReal),
               B, scale,
              )


	return nothing
end 

"""
	MatLMVMSetMultAlgorithm(petsclib::PetscLibType,B::AbstractPetscMat, alg::MatLMVMMultAlgorithm) 
Set the algorithm used by a `MatLMVM` for products

Logically collective

Input Parameters:
- `B`   - a `MatLMVM` matrix
- `alg` - one of the algorithm classes (`MAT_LMVM_MULT_RECURSIVE`, `MAT_LMVM_MULT_DENSE`, `MAT_LMVM_MULT_COMPACT_DENSE`)

Level: advanced

-seealso: [](ch_matrices), `MatLMVM`, `MatLMVMMultAlgorithm`, `MatLMVMGetMultAlgorithm()`

# External Links
$(_doc_external("KSP/MatLMVMSetMultAlgorithm"))
"""
function MatLMVMSetMultAlgorithm(petsclib::PetscLibType, B::AbstractPetscMat, alg::MatLMVMMultAlgorithm) end

@for_petsc function MatLMVMSetMultAlgorithm(petsclib::$UnionPetscLib, B::AbstractPetscMat, alg::MatLMVMMultAlgorithm )

    @chk ccall(
               (:MatLMVMSetMultAlgorithm, $petsc_library),
               PetscErrorCode,
               (CMat, MatLMVMMultAlgorithm),
               B, alg,
              )


	return nothing
end 

"""
	psi::PetscReal = MatLMVMSymBadBroydenGetPsi(petsclib::PetscLibType,B::AbstractPetscMat) 
Get the psi parameter for a Broyden class quasi

Input Parameter:
- `B` - The matrix

Output Parameter:
- `psi` - a number defining an update that is an affine combination of the BFGS update (psi = 1) and DFP update (psi = 0)

Level: advanced

-seealso: [](ch_ksp),
`MATLMVMSYMBROYDEN`, `MATLMVMSYMBADBROYDEN`,
`MATLMVMDFP`, `MATLMVMBFGS`,
`MatLMVMSymBadBroydenSetPsi()`,
`MatLMVMSymBroydenGetPhi()`, `MatLMVMSymBroydenSetPhi()`

# External Links
$(_doc_external("KSP/MatLMVMSymBadBroydenGetPsi"))
"""
function MatLMVMSymBadBroydenGetPsi(petsclib::PetscLibType, B::AbstractPetscMat) end

@for_petsc function MatLMVMSymBadBroydenGetPsi(petsclib::$UnionPetscLib, B::AbstractPetscMat )
	psi_ = Ref{$PetscReal}()

    @chk ccall(
               (:MatLMVMSymBadBroydenGetPsi, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscReal}),
               B, psi_,
              )

	psi = psi_[]

	return psi
end 

"""
	MatLMVMSymBadBroydenSetPsi(petsclib::PetscLibType,B::AbstractPetscMat, psi::PetscReal) 
Get the psi parameter for a Broyden class quasi

Input Parameters:
- `B`   - The matrix
- `psi` - a number defining an update that is a convex combination of the BFGS update (psi = 1) and DFP update (psi = 0)

Level: developer

-seealso: [](ch_ksp),
`MATLMVMSYMBROYDEN`, `MATLMVMSYMBADBROYDEN`,
`MATLMVMDFP`, `MATLMVMBFGS`,
`MatLMVMSymBadBroydenGetPsi()`,
`MatLMVMSymBroydenGetPhi()`, `MatLMVMSymBroydenSetPhi()`

# External Links
$(_doc_external("KSP/MatLMVMSymBadBroydenSetPsi"))
"""
function MatLMVMSymBadBroydenSetPsi(petsclib::PetscLibType, B::AbstractPetscMat, psi::PetscReal) end

@for_petsc function MatLMVMSymBadBroydenSetPsi(petsclib::$UnionPetscLib, B::AbstractPetscMat, psi::$PetscReal )

    @chk ccall(
               (:MatLMVMSymBadBroydenSetPsi, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscReal),
               B, psi,
              )


	return nothing
end 

"""
	phi::PetscReal = MatLMVMSymBroydenGetPhi(petsclib::PetscLibType,B::AbstractPetscMat) 
Get the phi parameter for a Broyden class quasi

Input Parameter:
- `B` - The matrix

Output Parameter:
- `phi` - a number defining an update that is an affine combination of the BFGS update (phi = 0) and DFP update (phi = 1)

Level: advanced

-seealso: [](ch_ksp),
`MATLMVMSYMBROYDEN`, `MATLMVMSYMBADBROYDEN`,
`MATLMVMDFP`, `MATLMVMBFGS`,
`MatLMVMSymBroydenSetPhi()`,
`MatLMVMSymBadBroydenGetPsi()`, `MatLMVMSymBadBroydenSetPsi()`

# External Links
$(_doc_external("KSP/MatLMVMSymBroydenGetPhi"))
"""
function MatLMVMSymBroydenGetPhi(petsclib::PetscLibType, B::AbstractPetscMat) end

@for_petsc function MatLMVMSymBroydenGetPhi(petsclib::$UnionPetscLib, B::AbstractPetscMat )
	phi_ = Ref{$PetscReal}()

    @chk ccall(
               (:MatLMVMSymBroydenGetPhi, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscReal}),
               B, phi_,
              )

	phi = phi_[]

	return phi
end 

"""
	MatLMVMSymBroydenSetDelta(petsclib::PetscLibType,B::AbstractPetscMat, delta::PetscScalar) 
Sets the starting value for the diagonal scaling vector computed
in the SymBrdn approximations (also works for BFGS and DFP).

Input Parameters:
- `B`     - `MATLMVM` matrix
- `delta` - initial value for diagonal scaling

Level: intermediate

-seealso: [](ch_ksp), `MATLMVMSYMBROYDEN`

# External Links
$(_doc_external("KSP/MatLMVMSymBroydenSetDelta"))
"""
function MatLMVMSymBroydenSetDelta(petsclib::PetscLibType, B::AbstractPetscMat, delta::PetscScalar) end

@for_petsc function MatLMVMSymBroydenSetDelta(petsclib::$UnionPetscLib, B::AbstractPetscMat, delta::$PetscScalar )

    @chk ccall(
               (:MatLMVMSymBroydenSetDelta, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscScalar),
               B, delta,
              )


	return nothing
end 

"""
	MatLMVMSymBroydenSetPhi(petsclib::PetscLibType,B::AbstractPetscMat, phi::PetscReal) 
Get the phi parameter for a Broyden class quasi

Input Parameters:
- `B`   - The matrix
- `phi` - a number defining an update that is a convex combination of the BFGS update (phi = 0) and DFP update (phi = 1)

Level: advanced

-seealso: [](ch_ksp),
`MATLMVMSYMBROYDEN`, `MATLMVMSYMBADBROYDEN`,
`MATLMVMDFP`, `MATLMVMBFGS`,
`MatLMVMSymBroydenGetPhi()`,
`MatLMVMSymBadBroydenGetPsi()`, `MatLMVMSymBadBroydenSetPsi()`

# External Links
$(_doc_external("KSP/MatLMVMSymBroydenSetPhi"))
"""
function MatLMVMSymBroydenSetPhi(petsclib::PetscLibType, B::AbstractPetscMat, phi::PetscReal) end

@for_petsc function MatLMVMSymBroydenSetPhi(petsclib::$UnionPetscLib, B::AbstractPetscMat, phi::$PetscReal )

    @chk ccall(
               (:MatLMVMSymBroydenSetPhi, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscReal),
               B, phi,
              )


	return nothing
end 

"""
	MatLMVMSymBroydenSetScaleType(petsclib::PetscLibType,B::AbstractPetscMat, stype::MatLMVMSymBroydenScaleType) 
Sets the scale type for symmetric Broyden

Input Parameters:
- `B`     - the `MATLMVM` matrix
- `stype` - scale type, see `MatLMVMSymBroydenScaleType`

Options Database Key:
- `-mat_lmvm_scale_type <none,scalar,diagonal>` - set the scaling type

Level: intermediate

MatLMVMSymBrdnScaleTypes:
- `MAT_LMVM_SYMBROYDEN_SCALE_NONE`       - use whatever initial Hessian is already there (will be the identity if the user does nothing)
- `MAT_LMVM_SYMBROYDEN_SCALE_SCALAR`     - use the Shanno scalar as the initial Hessian
- `MAT_LMVM_SYMBROYDEN_SCALE_DIAGONAL`   - use a diagonalized BFGS update as the initial Hessian
- `MAT_LMVM_SYMBROYDEN_SCALE_USER`       - same as `MAT_LMVM_SYMBROYDEN_NONE`
- `MAT_LMVM_SYMBROYDEN_SCALE_DECIDE`     - let PETSc decide

-seealso: [](ch_ksp), `MATLMVMSYMBROYDEN`, `MatCreateLMVMSymBroyden()`, `MatLMVMSymBroydenScaleType`

# External Links
$(_doc_external("KSP/MatLMVMSymBroydenSetScaleType"))
"""
function MatLMVMSymBroydenSetScaleType(petsclib::PetscLibType, B::AbstractPetscMat, stype::MatLMVMSymBroydenScaleType) end

@for_petsc function MatLMVMSymBroydenSetScaleType(petsclib::$UnionPetscLib, B::AbstractPetscMat, stype::MatLMVMSymBroydenScaleType )

    @chk ccall(
               (:MatLMVMSymBroydenSetScaleType, $petsc_library),
               PetscErrorCode,
               (CMat, MatLMVMSymBroydenScaleType),
               B, stype,
              )


	return nothing
end 

"""
	MatLMVMUpdate(petsclib::PetscLibType,B::AbstractPetscMat, X::AbstractPetscVec, F::AbstractPetscVec) 
Adds (X

Input Parameters:
- `B` - A `MATLMVM` matrix
- `X` - Solution vector
- `F` - Function vector

Level: intermediate

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMReset()`, `MatLMVMAllocate()`

# External Links
$(_doc_external("KSP/MatLMVMUpdate"))
"""
function MatLMVMUpdate(petsclib::PetscLibType, B::AbstractPetscMat, X::AbstractPetscVec, F::AbstractPetscVec) end

@for_petsc function MatLMVMUpdate(petsclib::$UnionPetscLib, B::AbstractPetscMat, X::AbstractPetscVec, F::AbstractPetscVec )

    @chk ccall(
               (:MatLMVMUpdate, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               B, X, F,
              )


	return nothing
end 

"""
	A::PetscMat,U::PetscMat,c::PetscVec,V::PetscMat = MatLRCGetMats(petsclib::PetscLibType,N::AbstractPetscMat) 
Returns the constituents of an LRC matrix

Not collective

Input Parameter:
- `N` - matrix of type `MATLRC`

Output Parameters:
- `A` - the (sparse) matrix
- `U` - first dense rectangular (tall and skinny) matrix
- `c` - a sequential vector containing the diagonal of C
- `V` - second dense rectangular (tall and skinny) matrix

Level: intermediate

-seealso: [](ch_matrices), `MatLRCSetMats()`, `Mat`, `MATLRC`, `MatCreateLRC()`

# External Links
$(_doc_external("Mat/MatLRCGetMats"))
"""
function MatLRCGetMats(petsclib::PetscLibType, N::AbstractPetscMat) end

@for_petsc function MatLRCGetMats(petsclib::$UnionPetscLib, N::AbstractPetscMat )
	A_ = Ref{CMat}()
	U_ = Ref{CMat}()
	c_ = Ref{CVec}()
	V_ = Ref{CMat}()

    @chk ccall(
               (:MatLRCGetMats, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}, Ptr{CMat}, Ptr{CVec}, Ptr{CMat}),
               N, A_, U_, c_, V_,
              )

	A = PetscMat(A_[], petsclib)
	U = PetscMat(U_[], petsclib)
	c = PetscVec(c_[], petsclib)
	V = PetscMat(V_[], petsclib)

	return A,U,c,V
end 

"""
	MatLRCSetMats(petsclib::PetscLibType,N::AbstractPetscMat, A::AbstractPetscMat, U::AbstractPetscMat, c::AbstractPetscVec, V::AbstractPetscMat) 
Sets the constituents of an LRC matrix

Logically collective

Input Parameters:
- `N` - matrix of type `MATLRC`
- `A` - the (sparse) matrix
- `U` - first dense rectangular (tall and skinny) matrix
- `c` - a sequential vector containing the diagonal of C
- `V` - second dense rectangular (tall and skinny) matrix

Level: intermediate

-seealso: [](ch_matrices), `MatLRCGetMats()`, `Mat`, `MATLRC`, `MatCreateLRC()`

# External Links
$(_doc_external("Mat/MatLRCSetMats"))
"""
function MatLRCSetMats(petsclib::PetscLibType, N::AbstractPetscMat, A::AbstractPetscMat, U::AbstractPetscMat, c::AbstractPetscVec, V::AbstractPetscMat) end

@for_petsc function MatLRCSetMats(petsclib::$UnionPetscLib, N::AbstractPetscMat, A::AbstractPetscMat, U::AbstractPetscMat, c::AbstractPetscVec, V::AbstractPetscMat )

    @chk ccall(
               (:MatLRCSetMats, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat, CVec, CMat),
               N, A, U, c, V,
              )


	return nothing
end 

"""
	MatLUFactor(petsclib::PetscLibType,mat::AbstractPetscMat, row::AbstractIS, col::AbstractIS, info::Vector{MatFactorInfo}) 
Performs in

Collective

Input Parameters:
- `mat`  - the matrix
- `row`  - row permutation
- `col`  - column permutation
- `info` - options for factorization, includes
-seealso: [](ch_matrices), [Matrix Factorization](sec_matfactor), `Mat`, `MatFactorType`, `MatLUFactorSymbolic()`, `MatLUFactorNumeric()`, `MatCholeskyFactor()`,
`MatGetOrdering()`, `MatSetUnfactored()`, `MatFactorInfo`, `MatGetFactor()`

# External Links
$(_doc_external("Mat/MatLUFactor"))
"""
function MatLUFactor(petsclib::PetscLibType, mat::AbstractPetscMat, row::AbstractIS, col::AbstractIS, info::Vector{MatFactorInfo}) end

@for_petsc function MatLUFactor(petsclib::$UnionPetscLib, mat::AbstractPetscMat, row::AbstractIS, col::AbstractIS, info::Vector{MatFactorInfo} )

    @chk ccall(
               (:MatLUFactor, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, CIS, Ptr{MatFactorInfo}),
               mat, row, col, info,
              )


	return nothing
end 

"""
	MatLUFactorNumeric(petsclib::PetscLibType,fact::AbstractPetscMat, mat::AbstractPetscMat, info::Vector{MatFactorInfo}) 
Performs numeric LU factorization of a matrix.
Call this routine after first calling `MatLUFactorSymbolic()` and `MatGetFactor()`.

Collective

Input Parameters:
- `fact` - the factor matrix obtained with `MatGetFactor()`
- `mat`  - the matrix
- `info` - options for factorization

Level: developer

-seealso: [](ch_matrices), `Mat`, [Matrix Factorization](sec_matfactor), `MatGetFactor()`, `MatFactorInfo`, `MatLUFactorSymbolic()`, `MatLUFactor()`, `MatCholeskyFactor()`

# External Links
$(_doc_external("Mat/MatLUFactorNumeric"))
"""
function MatLUFactorNumeric(petsclib::PetscLibType, fact::AbstractPetscMat, mat::AbstractPetscMat, info::Vector{MatFactorInfo}) end

@for_petsc function MatLUFactorNumeric(petsclib::$UnionPetscLib, fact::AbstractPetscMat, mat::AbstractPetscMat, info::Vector{MatFactorInfo} )

    @chk ccall(
               (:MatLUFactorNumeric, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, Ptr{MatFactorInfo}),
               fact, mat, info,
              )


	return nothing
end 

"""
	MatLUFactorSymbolic(petsclib::PetscLibType,fact::AbstractPetscMat, mat::AbstractPetscMat, row::AbstractIS, col::AbstractIS, info::Vector{MatFactorInfo}) 
Performs symbolic LU factorization of matrix.
Call this routine before calling `MatLUFactorNumeric()` and after `MatGetFactor()`.

Collective

Input Parameters:
- `fact` - the factor matrix obtained with `MatGetFactor()`
- `mat`  - the matrix
- `row`  - the row permutation
- `col`  - the column permutation
- `info` - options for factorization, includes
-seealso: [](ch_matrices), `Mat`, [Matrix Factorization](sec_matfactor), `MatGetFactor()`, `MatLUFactor()`, `MatLUFactorNumeric()`, `MatCholeskyFactor()`, `MatFactorInfo`, `MatFactorInfoInitialize()`

# External Links
$(_doc_external("Mat/MatLUFactorSymbolic"))
"""
function MatLUFactorSymbolic(petsclib::PetscLibType, fact::AbstractPetscMat, mat::AbstractPetscMat, row::AbstractIS, col::AbstractIS, info::Vector{MatFactorInfo}) end

@for_petsc function MatLUFactorSymbolic(petsclib::$UnionPetscLib, fact::AbstractPetscMat, mat::AbstractPetscMat, row::AbstractIS, col::AbstractIS, info::Vector{MatFactorInfo} )

    @chk ccall(
               (:MatLUFactorSymbolic, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CIS, CIS, Ptr{MatFactorInfo}),
               fact, mat, row, col, info,
              )


	return nothing
end 

"""
	MatLoad(petsclib::PetscLibType,mat::AbstractPetscMat, viewer::PetscViewer) 
Loads a matrix that has been stored in binary/HDF5 format
with `MatView()`.  The matrix format is determined from the options database.
Generates a parallel MPI matrix if the communicator has more than one
processor.  The default matrix type is `MATAIJ`.

Collective

Input Parameters:
- `mat`    - the newly loaded matrix, this needs to have been created with `MatCreate()`
or some related function before a call to `MatLoad()`
- `viewer` - `PETSCVIEWERBINARY`/`PETSCVIEWERHDF5` file viewer

Options Database Key:
- `-matload_block_size <bs>` - set block size

Level: beginner

-seealso: [](ch_matrices), `Mat`, `PetscViewerBinaryOpen()`, `PetscViewerSetType()`, `MatView()`, `VecLoad()`

# External Links
$(_doc_external("Mat/MatLoad"))
"""
function MatLoad(petsclib::PetscLibType, mat::AbstractPetscMat, viewer::PetscViewer) end

@for_petsc function MatLoad(petsclib::$UnionPetscLib, mat::AbstractPetscMat, viewer::PetscViewer )

    @chk ccall(
               (:MatLoad, $petsc_library),
               PetscErrorCode,
               (CMat, PetscViewer),
               mat, viewer,
              )


	return nothing
end 

"""
	B::PetscMat = MatMAIJGetAIJ(petsclib::PetscLibType,A::AbstractPetscMat) 
Get the `MATAIJ` matrix describing the blockwise action of the `MATMAIJ` matrix

Not Collective, but if the `MATMAIJ` matrix is parallel, the `MATAIJ` matrix is also parallel

Input Parameter:
- `A` - the `MATMAIJ` matrix

Output Parameter:
- `B` - the `MATAIJ` matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATMAIJ`, `MATAIJ`, `MatCreateMAIJ()`

# External Links
$(_doc_external("Mat/MatMAIJGetAIJ"))
"""
function MatMAIJGetAIJ(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatMAIJGetAIJ(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	B_ = Ref{CMat}()

    @chk ccall(
               (:MatMAIJGetAIJ, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, B_,
              )

	B = PetscMat(B_[], petsclib)

	return B
end 

"""
	B::PetscMat = MatMAIJRedimension(petsclib::PetscLibType,A::AbstractPetscMat, dof::PetscInt) 
Get a new `MATMAIJ` matrix with the same action, but for a different block size

Logically Collective

Input Parameters:
- `A`   - the `MATMAIJ` matrix
- `dof` - the block size for the new matrix

Output Parameter:
- `B` - the new `MATMAIJ` matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATMAIJ`, `MatCreateMAIJ()`

# External Links
$(_doc_external("Mat/MatMAIJRedimension"))
"""
function MatMAIJRedimension(petsclib::PetscLibType, A::AbstractPetscMat, dof::PetscInt) end

@for_petsc function MatMAIJRedimension(petsclib::$UnionPetscLib, A::AbstractPetscMat, dof::$PetscInt )
	B_ = Ref{CMat}()

    @chk ccall(
               (:MatMAIJRedimension, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CMat}),
               A, dof, B_,
              )

	B = PetscMat(B_[], petsclib)

	return B
end 

"""
	A_loc::PetscMat = MatMPIAIJGetLocalMat(petsclib::PetscLibType,A::AbstractPetscMat, scall::MatReuse) 
Creates a `MATSEQAIJ` from a `MATMPIAIJ` matrix.

Not Collective

Input Parameters:
- `A`     - the matrix
- `scall` - either `MAT_INITIAL_MATRIX` or `MAT_REUSE_MATRIX`

Output Parameter:
- `A_loc` - the local sequential matrix generated

Level: developer

-seealso: [](ch_matrices), `Mat`, `MATMPIAIJ`, `MatGetOwnershipRange()`, `MatMPIAIJGetLocalMatCondensed()`, `MatMPIAIJGetLocalMatMerge()`

# External Links
$(_doc_external("Mat/MatMPIAIJGetLocalMat"))
"""
function MatMPIAIJGetLocalMat(petsclib::PetscLibType, A::AbstractPetscMat, scall::MatReuse) end

@for_petsc function MatMPIAIJGetLocalMat(petsclib::$UnionPetscLib, A::AbstractPetscMat, scall::MatReuse )
	A_loc_ = Ref{CMat}()

    @chk ccall(
               (:MatMPIAIJGetLocalMat, $petsc_library),
               PetscErrorCode,
               (CMat, MatReuse, Ptr{CMat}),
               A, scall, A_loc_,
              )

	A_loc = PetscMat(A_loc_[], petsclib)

	return A_loc
end 

"""
	A_loc::PetscMat = MatMPIAIJGetLocalMatCondensed(petsclib::PetscLibType,A::AbstractPetscMat, scall::MatReuse, row::AbstractIS, col::AbstractIS) 
Creates a `MATSEQAIJ` matrix from an `MATMPIAIJ` matrix by taking all its local rows and NON

Not Collective

Input Parameters:
- `A`     - the matrix
- `scall` - either `MAT_INITIAL_MATRIX` or `MAT_REUSE_MATRIX`
- `row`   - index set of rows to extract (or `NULL`)
- `col`   - index set of columns to extract (or `NULL`)

Output Parameter:
- `A_loc` - the local sequential matrix generated

Level: developer

-seealso: [](ch_matrices), `Mat`, `MATMPIAIJ`, `MatGetOwnershipRange()`, `MatMPIAIJGetLocalMat()`

# External Links
$(_doc_external("Mat/MatMPIAIJGetLocalMatCondensed"))
"""
function MatMPIAIJGetLocalMatCondensed(petsclib::PetscLibType, A::AbstractPetscMat, scall::MatReuse, row::AbstractIS, col::AbstractIS) end

@for_petsc function MatMPIAIJGetLocalMatCondensed(petsclib::$UnionPetscLib, A::AbstractPetscMat, scall::MatReuse, row::AbstractIS, col::AbstractIS )
	row_ = Ref(row.ptr)
	col_ = Ref(col.ptr)
	A_loc_ = Ref{CMat}()

    @chk ccall(
               (:MatMPIAIJGetLocalMatCondensed, $petsc_library),
               PetscErrorCode,
               (CMat, MatReuse, Ptr{CIS}, Ptr{CIS}, Ptr{CMat}),
               A, scall, row_, col_, A_loc_,
              )

	row.ptr = row_[]
	col.ptr = col_[]
	A_loc = PetscMat(A_loc_[], petsclib)

	return A_loc
end 

"""
	glob::IS,A_loc::PetscMat = MatMPIAIJGetLocalMatMerge(petsclib::PetscLibType,A::AbstractPetscMat, scall::MatReuse) 
Creates a `MATSEQAIJ` from a `MATMPIAIJ` matrix by taking all its local rows and putting them into a sequential matrix with
mlocal rows and n columns. Where n is the sum of the number of columns of the diagonal and off-diagonal part

Not Collective

Input Parameters:
- `A`     - the matrix
- `scall` - either `MAT_INITIAL_MATRIX` or `MAT_REUSE_MATRIX`

Output Parameters:
- `glob`  - sequential `IS` with global indices associated with the columns of the local sequential matrix generated (can be `NULL`)
- `A_loc` - the local sequential matrix generated

Level: developer

-seealso: [](ch_matrices), `Mat`, `MATMPIAIJ`, `MatGetOwnershipRange()`, `MatMPIAIJGetLocalMat()`, `MatMPIAIJGetLocalMatCondensed()`

# External Links
$(_doc_external("Mat/MatMPIAIJGetLocalMatMerge"))
"""
function MatMPIAIJGetLocalMatMerge(petsclib::PetscLibType, A::AbstractPetscMat, scall::MatReuse) end

@for_petsc function MatMPIAIJGetLocalMatMerge(petsclib::$UnionPetscLib, A::AbstractPetscMat, scall::MatReuse )
	glob_ = Ref{CIS}()
	A_loc_ = Ref{CMat}()

    @chk ccall(
               (:MatMPIAIJGetLocalMatMerge, $petsc_library),
               PetscErrorCode,
               (CMat, MatReuse, Ptr{CIS}, Ptr{CMat}),
               A, scall, glob_, A_loc_,
              )

	glob = IS(glob_[], petsclib)
	A_loc = PetscMat(A_loc_[], petsclib)

	return glob,A_loc
end 

"""
	nz::PetscCount = MatMPIAIJGetNumberNonzeros(petsclib::PetscLibType,A::AbstractPetscMat) 
gets the number of nonzeros in the matrix on this MPI rank

Not Collective

Input Parameter:
- `A` - the matrix

Output Parameter:
- `nz` - the number of nonzeros

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATMPIAIJ`

# External Links
$(_doc_external("Mat/MatMPIAIJGetNumberNonzeros"))
"""
function MatMPIAIJGetNumberNonzeros(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatMPIAIJGetNumberNonzeros(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	nz_ = Ref{PetscCount}()

    @chk ccall(
               (:MatMPIAIJGetNumberNonzeros, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscCount}),
               A, nz_,
              )

	nz = nz_[]

	return nz
end 

# override for MatMPIAIJGetSeqAIJ; C signature: MatMPIAIJGetSeqAIJ(Mat A, Mat* Ad, Mat* Ao, PetscInt* colmap[])
"""
	Ad::PetscMat,Ao::PetscMat,colmap::Vector{PetscInt} = MatMPIAIJGetSeqAIJ(petsclib::PetscLibType,A::AbstractPetscMat) 
Returns the local pieces of this distributed matrix

Not Collective

Input Parameter:
- `A` - The `MATMPIAIJ` matrix

Output Parameters:
- `Ad`     - The local diagonal block as a `MATSEQAIJ` matrix
- `Ao`     - The local off-diagonal block as a `MATSEQAIJ` matrix
- `colmap` - An array mapping local column numbers of `Ao` to global column numbers of the parallel matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATMPIAIJ`, `MatMPIAIJGetLocalMat()`, `MatMPIAIJGetLocalMatCondensed()`, `MatCreateAIJ()`, `MATSEQAIJ`

# External Links
$(_doc_external("Mat/MatMPIAIJGetSeqAIJ"))
"""
function MatMPIAIJGetSeqAIJ(petsclib::PetscLibType, A::AbstractPetscMat) end

# `Ad`, `Ao` and `colmap` are owned by `A` (do not destroy them); `colmap` has one entry per
# column of `Ao`
@for_petsc function MatMPIAIJGetSeqAIJ(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	Ad_ = Ref{CMat}()
	Ao_ = Ref{CMat}()
	colmap_ = Ref{Ptr{$PetscInt}}(C_NULL)

    @chk ccall(
               (:MatMPIAIJGetSeqAIJ, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}, Ptr{CMat}, Ptr{Ptr{$PetscInt}}),
               A, Ad_, Ao_, colmap_,
              )

	Ad = PetscMat(Ad_[], petsclib)
	Ao = PetscMat(Ao_[], petsclib)
	_, ncols_o = MatGetLocalSize(petsclib, Ao)
	colmap = colmap_[] == C_NULL ? $PetscInt[] : unsafe_wrap(Array, colmap_[], Int(ncols_o); own = false)

	return Ad,Ao,colmap
end

"""
	MatMPIAIJSetPreallocation(petsclib::PetscLibType,B::AbstractPetscMat, d_nz::PetscInt, d_nnz::Union{Ptr, Vector{PetscInt}}, o_nz::PetscInt, o_nnz::Union{Ptr, Vector{PetscInt}}) 
Preallocates memory for a sparse parallel matrix in `MATMPIAIJ` format
(the default parallel PETSc format).  For good matrix assembly performance
the user should preallocate the matrix storage by setting the parameters
`d_nz` (or `d_nnz`) and `o_nz` (or `o_nnz`).

Collective

Input Parameters:
- `B`     - the matrix
- `d_nz`  - number of nonzeros per row in DIAGONAL portion of local submatrix
(same value is used for all local rows)
- `d_nnz` - array containing the number of nonzeros in the various rows of the
DIAGONAL portion of the local submatrix (possibly different for each row)
or `NULL` (`PETSC_NULL_INTEGER` in Fortran), if `d_nz` is used to specify the nonzero structure.
The size of this array is equal to the number of local rows, i.e 'm'.
For matrices that will be factored, you must leave room for (and set)
the diagonal entry even if it is zero.
- `o_nz`  - number of nonzeros per row in the OFF-DIAGONAL portion of local
submatrix (same value is used for all local rows).
- `o_nnz` - array containing the number of nonzeros in the various rows of the
OFF-DIAGONAL portion of the local submatrix (possibly different for
each row) or `NULL` (`PETSC_NULL_INTEGER` in Fortran), if `o_nz` is used to specify the nonzero
structure. The size of this array is equal to the number
of local rows, i.e 'm'.

-seealso: [](ch_matrices), `Mat`, [Sparse Matrices](sec_matsparse), `MATMPIAIJ`, `MATAIJ`, `MatCreate()`, `MatCreateSeqAIJ()`, `MatSetValues()`, `MatCreateAIJ()`, `MatMPIAIJSetPreallocationCSR()`,
`MatGetInfo()`, `PetscSplitOwnership()`, `MatSetPreallocationCOO()`, `MatSetValuesCOO()`

# External Links
$(_doc_external("Mat/MatMPIAIJSetPreallocation"))
"""
function MatMPIAIJSetPreallocation(petsclib::PetscLibType, B::AbstractPetscMat, d_nz::PetscInt, d_nnz::Union{Ptr, Vector{PetscInt}}, o_nz::PetscInt, o_nnz::Union{Ptr, Vector{PetscInt}}) end

@for_petsc function MatMPIAIJSetPreallocation(petsclib::$UnionPetscLib, B::AbstractPetscMat, d_nz::$PetscInt, d_nnz::Union{Ptr, Vector{$PetscInt}}, o_nz::$PetscInt, o_nnz::Union{Ptr, Vector{$PetscInt}} )

    @chk ccall(
               (:MatMPIAIJSetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}),
               B, d_nz, d_nnz, o_nz, o_nnz,
              )


	return nothing
end 

"""
	MatMPIAIJSetPreallocationCSR(petsclib::PetscLibType,B::AbstractPetscMat, i::Vector{PetscInt}, j::Vector{PetscInt}, v::Vector{PetscScalar}) 
Allocates memory for a sparse parallel matrix in `MATAIJ` format
(the default parallel PETSc format).

Collective

Input Parameters:
- `B` - the matrix
- `i` - the indices into `j` for the start of each local row (indices start with zero)
- `j` - the column indices for each local row (indices start with zero)
- `v` - optional values in the matrix

Level: developer

-seealso: [](ch_matrices), `Mat`, `MATMPIAIJ`, `MatCreate()`, `MatCreateSeqAIJ()`, `MatSetValues()`, `MatMPIAIJSetPreallocation()`, `MatCreateAIJ()`,
`MatCreateSeqAIJWithArrays()`, `MatCreateMPIAIJWithSplitArrays()`, `MatCreateMPIAIJWithArrays()`, `MatSetPreallocationCOO()`, `MatSetValuesCOO()`

# External Links
$(_doc_external("Mat/MatMPIAIJSetPreallocationCSR"))
"""
function MatMPIAIJSetPreallocationCSR(petsclib::PetscLibType, B::AbstractPetscMat, i::Vector{PetscInt}, j::Vector{PetscInt}, v::Vector{PetscScalar}) end

@for_petsc function MatMPIAIJSetPreallocationCSR(petsclib::$UnionPetscLib, B::AbstractPetscMat, i::Vector{$PetscInt}, j::Vector{$PetscInt}, v::Vector{$PetscScalar} )

    @chk ccall(
               (:MatMPIAIJSetPreallocationCSR, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}),
               B, i, j, v,
              )


	return nothing
end 

"""
	MatMPIAIJSetUseScalableIncreaseOverlap(petsclib::PetscLibType,A::AbstractPetscMat, sc::PetscBool) 
Determine if the matrix uses a scalable algorithm to compute the overlap

Collective

Input Parameters:
- `A`  - the matrix
- `sc` - `PETSC_TRUE` indicates use the scalable algorithm (default is not to use the scalable algorithm)

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATMPIAIJ`

# External Links
$(_doc_external("Mat/MatMPIAIJSetUseScalableIncreaseOverlap"))
"""
function MatMPIAIJSetUseScalableIncreaseOverlap(petsclib::PetscLibType, A::AbstractPetscMat, sc::PetscBool) end

@for_petsc function MatMPIAIJSetUseScalableIncreaseOverlap(petsclib::$UnionPetscLib, A::AbstractPetscMat, sc::PetscBool )

    @chk ccall(
               (:MatMPIAIJSetUseScalableIncreaseOverlap, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               A, sc,
              )


	return nothing
end 

"""
	B::PetscMat = MatMPIAdjCreateNonemptySubcommMat(petsclib::PetscLibType,A::AbstractPetscMat) 
create the same `MATMPIADJ` matrix on a subcommunicator containing only processes owning a positive number of rows

Collective

Input Parameter:
- `A` - original `MATMPIADJ` matrix

Output Parameter:
- `B` - matrix on subcommunicator, `NULL` on MPI processes that own zero rows of `A`

Level: developer

-seealso: [](ch_matrices), `Mat`, `MATMPIADJ`, `MatCreateMPIAdj()`

# External Links
$(_doc_external("Mat/MatMPIAdjCreateNonemptySubcommMat"))
"""
function MatMPIAdjCreateNonemptySubcommMat(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatMPIAdjCreateNonemptySubcommMat(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	B_ = Ref{CMat}()

    @chk ccall(
               (:MatMPIAdjCreateNonemptySubcommMat, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, B_,
              )

	B = PetscMat(B_[], petsclib)

	return B
end 

"""
	i::PetscInt,j::PetscInt,values::PetscInt = MatMPIAdjSetPreallocation(petsclib::PetscLibType,B::AbstractPetscMat) 
Sets the array used for storing the matrix elements

Logically Collective

Input Parameters:
- `B`      - the matrix
- `i`      - the indices into `j` for the start of each row
- `j`      - the column indices for each row (sorted for each row).
The indices in `i` and `j` start with zero (NOT with one).
- `values` - [use `NULL` if not provided] edge weights

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `MatCreateMPIAdj()`, `MatSetValues()`, `MATMPIADJ`

# External Links
$(_doc_external("Mat/MatMPIAdjSetPreallocation"))
"""
function MatMPIAdjSetPreallocation(petsclib::PetscLibType, B::AbstractPetscMat) end

@for_petsc function MatMPIAdjSetPreallocation(petsclib::$UnionPetscLib, B::AbstractPetscMat )
	i_ = Ref{$PetscInt}()
	j_ = Ref{$PetscInt}()
	values_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatMPIAdjSetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
               B, i_, j_, values_,
              )

	i = i_[]
	j = j_[]
	values = values_[]

	return i,j,values
end 

"""
	B::PetscMat = MatMPIAdjToSeq(petsclib::PetscLibType,A::AbstractPetscMat) 
Converts an parallel `MATMPIADJ` matrix to complete `MATMPIADJ` on each process (needed by sequential partitioners)

Logically Collective

Input Parameter:
- `A` - the matrix

Output Parameter:
- `B` - the same matrix on all processes

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATMPIADJ`, `MatCreate()`, `MatCreateMPIAdj()`, `MatSetValues()`, `MatMPIAdjToSeqRankZero()`

# External Links
$(_doc_external("Mat/MatMPIAdjToSeq"))
"""
function MatMPIAdjToSeq(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatMPIAdjToSeq(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	B_ = Ref{CMat}()

    @chk ccall(
               (:MatMPIAdjToSeq, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, B_,
              )

	B = PetscMat(B_[], petsclib)

	return B
end 

"""
	B::PetscMat = MatMPIAdjToSeqRankZero(petsclib::PetscLibType,A::AbstractPetscMat) 
Converts an parallel `MATMPIADJ` matrix to complete `MATMPIADJ` on rank zero (needed by sequential partitioners)

Logically Collective

Input Parameter:
- `A` - the matrix

Output Parameter:
- `B` - the same matrix on rank zero, not set on other ranks

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATMPIADJ`, `MatCreate()`, `MatCreateMPIAdj()`, `MatSetValues()`, `MatMPIAdjToSeq()`

# External Links
$(_doc_external("Mat/MatMPIAdjToSeqRankZero"))
"""
function MatMPIAdjToSeqRankZero(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatMPIAdjToSeqRankZero(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	B_ = Ref{CMat}()

    @chk ccall(
               (:MatMPIAdjToSeqRankZero, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, B_,
              )

	B = PetscMat(B_[], petsclib)

	return B
end 

"""
	colmap::Ptr{PetscInt} = MatMPIBAIJGetSeqBAIJ(petsclib::PetscLibType,A::AbstractPetscMat, Ad::AbstractPetscMat, Ao::AbstractPetscMat) 

# External Links
$(_doc_external("Mat/MatMPIBAIJGetSeqBAIJ"))
"""
function MatMPIBAIJGetSeqBAIJ(petsclib::PetscLibType, A::AbstractPetscMat, Ad::AbstractPetscMat, Ao::AbstractPetscMat) end

@for_petsc function MatMPIBAIJGetSeqBAIJ(petsclib::$UnionPetscLib, A::AbstractPetscMat, Ad::AbstractPetscMat, Ao::AbstractPetscMat )
	Ad_ = Ref(Ad.ptr)
	Ao_ = Ref(Ao.ptr)
	colmap_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:MatMPIBAIJGetSeqBAIJ, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}, Ptr{CMat}, Ptr{Ptr{$PetscInt}}),
               A, Ad_, Ao_, colmap_,
              )

	Ad.ptr = Ad_[]
	Ao.ptr = Ao_[]
	colmap = colmap_[]

	return colmap
end 

"""
	MatMPIBAIJSetHashTableFactor(petsclib::PetscLibType,mat::AbstractPetscMat, fact::PetscReal) 
Sets the factor required to compute the size of the matrices hash table

Input Parameters:
- `mat`  - the matrix
- `fact` - factor

Options Database Key:
- `-mat_use_hash_table <fact>` - provide the factor

Level: advanced

-seealso: `Mat`, `MATMPIBAIJ`, `MatSetOption()`

# External Links
$(_doc_external("Mat/MatMPIBAIJSetHashTableFactor"))
"""
function MatMPIBAIJSetHashTableFactor(petsclib::PetscLibType, mat::AbstractPetscMat, fact::PetscReal) end

@for_petsc function MatMPIBAIJSetHashTableFactor(petsclib::$UnionPetscLib, mat::AbstractPetscMat, fact::$PetscReal )

    @chk ccall(
               (:MatMPIBAIJSetHashTableFactor, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscReal),
               mat, fact,
              )


	return nothing
end 

"""
	MatMPIBAIJSetPreallocation(petsclib::PetscLibType,B::AbstractPetscMat, bs::PetscInt, d_nz::PetscInt, d_nnz::Union{Ptr, Vector{PetscInt}}, o_nz::PetscInt, o_nnz::Union{Ptr, Vector{PetscInt}}) 
Allocates memory for a sparse parallel matrix in `MATMPIBAIJ` format
(block compressed row).

Collective

Input Parameters:
- `B`     - the matrix
- `bs`    - size of block, the blocks are ALWAYS square. One can use `MatSetBlockSizes()` to set a different row and column blocksize but the row
blocksize always defines the size of the blocks. The column blocksize sets the blocksize of the vectors obtained with `MatCreateVecs()`
- `d_nz`  - number of block nonzeros per block row in diagonal portion of local
submatrix  (same for all local rows)
- `d_nnz` - array containing the number of block nonzeros in the various block rows
of the in diagonal portion of the local (possibly different for each block
row) or `NULL`.  If you plan to factor the matrix you must leave room for the diagonal entry and
set it even if it is zero.
- `o_nz`  - number of block nonzeros per block row in the off-diagonal portion of local
submatrix (same for all local rows).
- `o_nnz` - array containing the number of nonzeros in the various block rows of the
off-diagonal portion of the local submatrix (possibly different for
each block row) or `NULL`.

If the *_nnz parameter is given then the *_nz parameter is ignored

Options Database Keys:
- `-mat_block_size`            - size of the blocks to use
- `-mat_use_hash_table <fact>` - set hash table factor

Level: intermediate

-seealso: `Mat`, `MATMPIBAIJ`, `MatCreate()`, `MatCreateSeqBAIJ()`, `MatSetValues()`, `MatCreateBAIJ()`, `MatMPIBAIJSetPreallocationCSR()`, `PetscSplitOwnership()`

# External Links
$(_doc_external("Mat/MatMPIBAIJSetPreallocation"))
"""
function MatMPIBAIJSetPreallocation(petsclib::PetscLibType, B::AbstractPetscMat, bs::PetscInt, d_nz::PetscInt, d_nnz::Union{Ptr, Vector{PetscInt}}, o_nz::PetscInt, o_nnz::Union{Ptr, Vector{PetscInt}}) end

@for_petsc function MatMPIBAIJSetPreallocation(petsclib::$UnionPetscLib, B::AbstractPetscMat, bs::$PetscInt, d_nz::$PetscInt, d_nnz::Union{Ptr, Vector{$PetscInt}}, o_nz::$PetscInt, o_nnz::Union{Ptr, Vector{$PetscInt}} )

    @chk ccall(
               (:MatMPIBAIJSetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}),
               B, bs, d_nz, d_nnz, o_nz, o_nnz,
              )


	return nothing
end 

"""
	MatMPIBAIJSetPreallocationCSR(petsclib::PetscLibType,B::AbstractPetscMat, bs::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, v::Union{Ptr, Vector{PetscScalar}}) 
Creates a sparse parallel matrix in `MATBAIJ` format using the given nonzero structure and (optional) numerical values

Collective

Input Parameters:
- `B`  - the matrix
- `bs` - the block size
- `i`  - the indices into `j` for the start of each local row (starts with zero)
- `j`  - the column indices for each local row (starts with zero) these must be sorted for each row
- `v`  - optional values in the matrix, use `NULL` if not provided

Level: advanced

-seealso: `Mat`, `MatCreate()`, `MatCreateSeqAIJ()`, `MatSetValues()`, `MatMPIBAIJSetPreallocation()`, `MatCreateAIJ()`, `MATMPIAIJ`, `MatCreateMPIBAIJWithArrays()`, `MATMPIBAIJ`

# External Links
$(_doc_external("Mat/MatMPIBAIJSetPreallocationCSR"))
"""
function MatMPIBAIJSetPreallocationCSR(petsclib::PetscLibType, B::AbstractPetscMat, bs::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, v::Union{Ptr, Vector{PetscScalar}}) end

@for_petsc function MatMPIBAIJSetPreallocationCSR(petsclib::$UnionPetscLib, B::AbstractPetscMat, bs::$PetscInt, i::Vector{$PetscInt}, j::Vector{$PetscInt}, v::Union{Ptr, Vector{$PetscScalar}} )

    @chk ccall(
               (:MatMPIBAIJSetPreallocationCSR, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}),
               B, bs, i, j, v,
              )


	return nothing
end 

"""
	data::PetscScalar = MatMPIDenseSetPreallocation(petsclib::PetscLibType,B::AbstractPetscMat) 
Sets the array used to store the matrix entries

Collective

Input Parameters:
- `B`    - the matrix
- `data` - optional location of matrix data.  Set to `NULL` for PETSc
to control all matrix memory allocation.

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATMPIDENSE`, `MatCreate()`, `MatCreateSeqDense()`, `MatSetValues()`

# External Links
$(_doc_external("Mat/MatMPIDenseSetPreallocation"))
"""
function MatMPIDenseSetPreallocation(petsclib::PetscLibType, B::AbstractPetscMat) end

@for_petsc function MatMPIDenseSetPreallocation(petsclib::$UnionPetscLib, B::AbstractPetscMat )
	data_ = Ref{$PetscScalar}()

    @chk ccall(
               (:MatMPIDenseSetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscScalar}),
               B, data_,
              )

	data = data_[]

	return data
end 

"""
	MatMPISBAIJSetPreallocation(petsclib::PetscLibType,B::AbstractPetscMat, bs::PetscInt, d_nz::PetscInt, d_nnz::Vector{PetscInt}, o_nz::PetscInt, o_nnz::Vector{PetscInt}) 
For good matrix assembly performance
the user should preallocate the matrix storage by setting the parameters
d_nz (or d_nnz) and o_nz (or o_nnz).  By setting these parameters accurately,
performance can be increased by more than a factor of 50.

Collective

Input Parameters:
- `B`     - the matrix
- `bs`    - size of block, the blocks are ALWAYS square. One can use MatSetBlockSizes() to set a different row and column blocksize but the row
blocksize always defines the size of the blocks. The column blocksize sets the blocksize of the vectors obtained with MatCreateVecs()
- `d_nz`  - number of block nonzeros per block row in diagonal portion of local
submatrix  (same for all local rows)
- `d_nnz` - array containing the number of block nonzeros in the various block rows
in the upper triangular and diagonal part of the in diagonal portion of the local
(possibly different for each block row) or `NULL`.  If you plan to factor the matrix you must leave room
for the diagonal entry and set a value even if it is zero.
- `o_nz`  - number of block nonzeros per block row in the off-diagonal portion of local
submatrix (same for all local rows).
- `o_nnz` - array containing the number of nonzeros in the various block rows of the
off-diagonal portion of the local submatrix that is right of the diagonal
(possibly different for each block row) or `NULL`.

Options Database Keys:
- `-mat_no_unroll`  - uses code that does not unroll the loops in the
block calculations (much slower)
- `-mat_block_size` - size of the blocks to use

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATMPISBAIJ`, `MATSBAIJ`, `MatCreate()`, `MatCreateSeqSBAIJ()`, `MatSetValues()`, `MatCreateBAIJ()`, `PetscSplitOwnership()`

# External Links
$(_doc_external("Mat/MatMPISBAIJSetPreallocation"))
"""
function MatMPISBAIJSetPreallocation(petsclib::PetscLibType, B::AbstractPetscMat, bs::PetscInt, d_nz::PetscInt, d_nnz::Vector{PetscInt}, o_nz::PetscInt, o_nnz::Vector{PetscInt}) end

@for_petsc function MatMPISBAIJSetPreallocation(petsclib::$UnionPetscLib, B::AbstractPetscMat, bs::$PetscInt, d_nz::$PetscInt, d_nnz::Vector{$PetscInt}, o_nz::$PetscInt, o_nnz::Vector{$PetscInt} )

    @chk ccall(
               (:MatMPISBAIJSetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}),
               B, bs, d_nz, d_nnz, o_nz, o_nnz,
              )


	return nothing
end 

"""
	MatMPISBAIJSetPreallocationCSR(petsclib::PetscLibType,B::AbstractPetscMat, bs::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, v::Vector{PetscScalar}) 
Creates a sparse parallel matrix in `MATMPISBAIJ` format using the given nonzero structure and (optional) numerical values

Collective

Input Parameters:
- `B`  - the matrix
- `bs` - the block size
- `i`  - the indices into `j` for the start of each local row (indices start with zero)
- `j`  - the column indices for each local row (indices start with zero) these must be sorted for each row
- `v`  - optional values in the matrix, pass `NULL` if not provided

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATMPISBAIJ`, `MatCreate()`, `MatCreateSeqAIJ()`, `MatSetValues()`, `MatMPIBAIJSetPreallocation()`, `MatCreateAIJ()`, `MATMPIAIJ`,
`MatCreateMPISBAIJWithArrays()`

# External Links
$(_doc_external("Mat/MatMPISBAIJSetPreallocationCSR"))
"""
function MatMPISBAIJSetPreallocationCSR(petsclib::PetscLibType, B::AbstractPetscMat, bs::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, v::Vector{PetscScalar}) end

@for_petsc function MatMPISBAIJSetPreallocationCSR(petsclib::$UnionPetscLib, B::AbstractPetscMat, bs::$PetscInt, i::Vector{$PetscInt}, j::Vector{$PetscInt}, v::Vector{$PetscScalar} )

    @chk ccall(
               (:MatMPISBAIJSetPreallocationCSR, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}),
               B, bs, i, j, v,
              )


	return nothing
end 

"""
	A_loc::PetscMat = MatMPISELLGetLocalMatCondensed(petsclib::PetscLibType,A::AbstractPetscMat, scall::MatReuse, row::Union{Ptr, AbstractIS}, col::Union{Ptr, AbstractIS}) 
Creates a `MATSEQSELL` matrix from an `MATMPISELL` matrix by
taking all its local rows and NON-ZERO columns

Not Collective

Input Parameters:
- `A`     - the matrix
- `scall` - either `MAT_INITIAL_MATRIX` or `MAT_REUSE_MATRIX`
- `row`   - index sets of rows to extract (or `NULL`)
- `col`   - index sets of columns to extract (or `NULL`)

Output Parameter:
- `A_loc` - the local sequential matrix generated

Level: advanced

-seealso: `Mat`, `MATSEQSELL`, `MATMPISELL`, `MatGetOwnershipRange()`, `MatMPISELLGetLocalMat()`

# External Links
$(_doc_external("Mat/MatMPISELLGetLocalMatCondensed"))
"""
function MatMPISELLGetLocalMatCondensed(petsclib::PetscLibType, A::AbstractPetscMat, scall::MatReuse, row::Union{Ptr, AbstractIS}, col::Union{Ptr, AbstractIS}) end

@for_petsc function MatMPISELLGetLocalMatCondensed(petsclib::$UnionPetscLib, A::AbstractPetscMat, scall::MatReuse, row::Union{Ptr, AbstractIS}, col::Union{Ptr, AbstractIS} )
	row_ = Ref(row.ptr)
	col_ = Ref(col.ptr)
	A_loc_ = Ref{CMat}()

    @chk ccall(
               (:MatMPISELLGetLocalMatCondensed, $petsc_library),
               PetscErrorCode,
               (CMat, MatReuse, Ptr{CIS}, Ptr{CIS}, Ptr{CMat}),
               A, scall, row_, col_, A_loc_,
              )

	row.ptr = row_[]
	col.ptr = col_[]
	A_loc = PetscMat(A_loc_[], petsclib)

	return A_loc
end 

"""
	Ad::PetscMat,Ao::PetscMat,colmap::Ptr{PetscInt} = MatMPISELLGetSeqSELL(petsclib::PetscLibType,A::AbstractPetscMat) 
Returns the local pieces of this distributed matrix

Not Collective

Input Parameter:
- `A` - the `MATMPISELL` matrix

Output Parameters:
- `Ad`     - The diagonal portion of `A`
- `Ao`     - The off-diagonal portion of `A`
- `colmap` - An array mapping local column numbers of `Ao` to global column numbers of the parallel matrix

Level: advanced

-seealso: `Mat`, `MATSEQSELL`, `MATMPISELL`

# External Links
$(_doc_external("Mat/MatMPISELLGetSeqSELL"))
"""
function MatMPISELLGetSeqSELL(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatMPISELLGetSeqSELL(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	Ad_ = Ref{CMat}()
	Ao_ = Ref{CMat}()
	colmap_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:MatMPISELLGetSeqSELL, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}, Ptr{CMat}, Ptr{Ptr{$PetscInt}}),
               A, Ad_, Ao_, colmap_,
              )

	Ad = PetscMat(Ad_[], petsclib)
	Ao = PetscMat(Ao_[], petsclib)
	colmap = colmap_[]

	return Ad,Ao,colmap
end 

"""
	MatMPISELLSetPreallocation(petsclib::PetscLibType,B::AbstractPetscMat, d_nz::PetscInt, d_nnz::Vector{PetscInt}, o_nz::PetscInt, o_nnz::Vector{PetscInt}) 
Preallocates memory for a `MATMPISELL` sparse parallel matrix in sell format.
For good matrix assembly performance the user should preallocate the matrix storage by
setting the parameters `d_nz` (or `d_nnz`) and `o_nz` (or `o_nnz`).

Collective

Input Parameters:
- `B`     - the matrix
- `d_nz`  - number of nonzeros per row in DIAGONAL portion of local submatrix
(same value is used for all local rows)
- `d_nnz` - array containing the number of nonzeros in the various rows of the
DIAGONAL portion of the local submatrix (possibly different for each row)
or NULL (`PETSC_NULL_INTEGER` in Fortran), if `d_nz` is used to specify the nonzero structure.
The size of this array is equal to the number of local rows, i.e 'm'.
For matrices that will be factored, you must leave room for (and set)
the diagonal entry even if it is zero.
- `o_nz`  - number of nonzeros per row in the OFF-DIAGONAL portion of local
submatrix (same value is used for all local rows).
- `o_nnz` - array containing the number of nonzeros in the various rows of the
OFF-DIAGONAL portion of the local submatrix (possibly different for
each row) or NULL (`PETSC_NULL_INTEGER` in Fortran), if `o_nz` is used to specify the nonzero
structure. The size of this array is equal to the number
of local rows, i.e 'm'.

Example usage:
Consider the following 8x8 matrix with 34 non-zero values, that is
assembled across 3 processors. Lets assume that proc0 owns 3 rows,
proc1 owns 3 rows, proc2 owns 2 rows. This division can be shown
as follows

-seealso: `Mat`, `MatCreate()`, `MatCreateSeqSELL()`, `MatSetValues()`, `MatCreateSELL()`,
`MATMPISELL`, `MatGetInfo()`, `PetscSplitOwnership()`, `MATSELL`

# External Links
$(_doc_external("Mat/MatMPISELLSetPreallocation"))
"""
function MatMPISELLSetPreallocation(petsclib::PetscLibType, B::AbstractPetscMat, d_nz::PetscInt, d_nnz::Vector{PetscInt}, o_nz::PetscInt, o_nnz::Vector{PetscInt}) end

@for_petsc function MatMPISELLSetPreallocation(petsclib::$UnionPetscLib, B::AbstractPetscMat, d_nz::$PetscInt, d_nnz::Vector{$PetscInt}, o_nz::$PetscInt, o_nnz::Vector{$PetscInt} )

    @chk ccall(
               (:MatMPISELLSetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}),
               B, d_nz, d_nnz, o_nz, o_nnz,
              )


	return nothing
end 

"""
	y::PetscMat = MatMatInterpolate(petsclib::PetscLibType,A::AbstractPetscMat, x::AbstractPetscMat) 
Y = A*X or A^T*X depending on the shape of `A`

Neighbor-wise Collective

Input Parameters:
- `A` - the matrix
- `x` - the input dense matrix

Output Parameter:
- `y` - the output dense matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatInterpolate()`, `MatRestrict()`, `MatMatRestrict()`, `PCMG`

# External Links
$(_doc_external("Mat/MatMatInterpolate"))
"""
function MatMatInterpolate(petsclib::PetscLibType, A::AbstractPetscMat, x::AbstractPetscMat) end

@for_petsc function MatMatInterpolate(petsclib::$UnionPetscLib, A::AbstractPetscMat, x::AbstractPetscMat )
	y_ = Ref{CMat}()

    @chk ccall(
               (:MatMatInterpolate, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, Ptr{CMat}),
               A, x, y_,
              )

	y = PetscMat(y_[], petsclib)

	return y
end 

"""
	y::PetscMat = MatMatInterpolateAdd(petsclib::PetscLibType,A::AbstractPetscMat, x::AbstractPetscMat, w::AbstractPetscMat) 
Y = W + A*X or W + A^T*X depending on the shape of `A`

Neighbor-wise Collective

Input Parameters:
- `A` - the matrix
- `x` - the input dense matrix to be multiplied
- `w` - the input dense matrix to be added to the result

Output Parameter:
- `y` - the output dense matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatInterpolateAdd()`, `MatMatInterpolate()`, `MatMatRestrict()`, `PCMG`

# External Links
$(_doc_external("Mat/MatMatInterpolateAdd"))
"""
function MatMatInterpolateAdd(petsclib::PetscLibType, A::AbstractPetscMat, x::AbstractPetscMat, w::AbstractPetscMat) end

@for_petsc function MatMatInterpolateAdd(petsclib::$UnionPetscLib, A::AbstractPetscMat, x::AbstractPetscMat, w::AbstractPetscMat )
	y_ = Ref{CMat}()

    @chk ccall(
               (:MatMatInterpolateAdd, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat, Ptr{CMat}),
               A, x, w, y_,
              )

	y = PetscMat(y_[], petsclib)

	return y
end 

"""
	D::PetscMat = MatMatMatMult(petsclib::PetscLibType,A::AbstractPetscMat, B::AbstractPetscMat, C::AbstractPetscMat, scall::MatReuse, fill::PetscReal) 
Performs matrix

Neighbor-wise Collective

Input Parameters:
- `A`     - the left matrix
- `B`     - the middle matrix
- `C`     - the right matrix
- `scall` - either `MAT_INITIAL_MATRIX` or `MAT_REUSE_MATRIX`
- `fill`  - expected fill as ratio of nnz(D)/(nnz(A) + nnz(B)+nnz(C)), use `PETSC_DETERMINE` or `PETSC_CURRENT` if you do not have a good estimate
if the result is a dense matrix this is irrelevant

Output Parameter:
- `D` - the product matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatProductCreate()`, `MATPRODUCT_ABC`, `MatMatMult`, `MatPtAP()`, `MatMatTransposeMult()`, `MatTransposeMatMult()`

# External Links
$(_doc_external("Mat/MatMatMatMult"))
"""
function MatMatMatMult(petsclib::PetscLibType, A::AbstractPetscMat, B::AbstractPetscMat, C::AbstractPetscMat, scall::MatReuse, fill::PetscReal) end

@for_petsc function MatMatMatMult(petsclib::$UnionPetscLib, A::AbstractPetscMat, B::AbstractPetscMat, C::AbstractPetscMat, scall::MatReuse, fill::$PetscReal )
	D_ = Ref{CMat}()

    @chk ccall(
               (:MatMatMatMult, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat, MatReuse, $PetscReal, Ptr{CMat}),
               A, B, C, scall, fill, D_,
              )

	D = PetscMat(D_[], petsclib)

	return D
end 

"""
	C::PetscMat = MatMatMult(petsclib::PetscLibType,A::AbstractPetscMat, B::AbstractPetscMat, scall::MatReuse, fill::PetscReal) 
Performs matrix

Neighbor-wise Collective

Input Parameters:
- `A`     - the left matrix
- `B`     - the right matrix
- `scall` - either `MAT_INITIAL_MATRIX` or `MAT_REUSE_MATRIX`
- `fill`  - expected fill as ratio of nnz(C)/(nnz(A) + nnz(B)), use `PETSC_DETERMINE` or `PETSC_CURRENT` if you do not have a good estimate
if the result is a dense matrix this is irrelevant

Output Parameter:
- `C` - the product matrix

-seealso: [](ch_matrices), `Mat`, `MatProductType`, `MATPRODUCT_AB`, `MatTransposeMatMult()`, `MatMatTransposeMult()`, `MatPtAP()`, `MatProductCreate()`, `MatProductSymbolic()`, `MatProductReplaceMats()`, `MatProductNumeric()`

# External Links
$(_doc_external("Mat/MatMatMult"))
"""
function MatMatMult(petsclib::PetscLibType, A::AbstractPetscMat, B::AbstractPetscMat, scall::MatReuse, fill::PetscReal) end

@for_petsc function MatMatMult(petsclib::$UnionPetscLib, A::AbstractPetscMat, B::AbstractPetscMat, scall::MatReuse, fill::$PetscReal )
	C_ = Ref{CMat}()

    @chk ccall(
               (:MatMatMult, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, MatReuse, $PetscReal, Ptr{CMat}),
               A, B, scall, fill, C_,
              )

	C = PetscMat(C_[], petsclib)

	return C
end 

"""
	flg::PetscBool = MatMatMultEqual(petsclib::PetscLibType,A::AbstractPetscMat, B::AbstractPetscMat, C::AbstractPetscMat, n::PetscInt) 
Test A*B*x = C*x for n random vector x

Collective

Input Parameters:
- `A` - the first matrix
- `B` - the second matrix
- `C` - the third matrix
- `n` - number of random vectors to be tested

Output Parameter:
- `flg` - `PETSC_TRUE` if the products are equal; `PETSC_FALSE` otherwise.

Level: intermediate

-seealso: `Mat`, `MatMultEqual()`, `MatMultAddEqual()`, `MatMultTransposeEqual()`

# External Links
$(_doc_external("Mat/MatMatMultEqual"))
"""
function MatMatMultEqual(petsclib::PetscLibType, A::AbstractPetscMat, B::AbstractPetscMat, C::AbstractPetscMat, n::PetscInt) end

@for_petsc function MatMatMultEqual(petsclib::$UnionPetscLib, A::AbstractPetscMat, B::AbstractPetscMat, C::AbstractPetscMat, n::$PetscInt )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatMatMultEqual, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat, $PetscInt, Ptr{PetscBool}),
               A, B, C, n, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	y::PetscMat = MatMatRestrict(petsclib::PetscLibType,A::AbstractPetscMat, x::AbstractPetscMat) 
Y = A*X or A^T*X depending on the shape of `A`

Neighbor-wise Collective

Input Parameters:
- `A` - the matrix
- `x` - the input dense matrix

Output Parameter:
- `y` - the output dense matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatRestrict()`, `MatInterpolate()`, `MatMatInterpolate()`, `PCMG`

# External Links
$(_doc_external("Mat/MatMatRestrict"))
"""
function MatMatRestrict(petsclib::PetscLibType, A::AbstractPetscMat, x::AbstractPetscMat) end

@for_petsc function MatMatRestrict(petsclib::$UnionPetscLib, A::AbstractPetscMat, x::AbstractPetscMat )
	y_ = Ref{CMat}()

    @chk ccall(
               (:MatMatRestrict, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, Ptr{CMat}),
               A, x, y_,
              )

	y = PetscMat(y_[], petsclib)

	return y
end 

"""
	MatMatSolve(petsclib::PetscLibType,A::AbstractPetscMat, B::AbstractPetscMat, X::AbstractPetscMat) 
Solves A X = B, given a factored matrix.

Neighbor-wise Collective

Input Parameters:
- `A` - the factored matrix
- `B` - the right-hand-side matrix `MATDENSE` (or sparse `MATAIJ`-- when using MUMPS)

Output Parameter:
- `X` - the result matrix (dense matrix)

Level: developer

-seealso: [](ch_matrices), `Mat`, [Matrix Factorization](sec_matfactor), `MatGetFactor()`, `MatSolve()`, `MatMatSolveTranspose()`, `MatLUFactor()`, `MatCholeskyFactor()`

# External Links
$(_doc_external("Mat/MatMatSolve"))
"""
function MatMatSolve(petsclib::PetscLibType, A::AbstractPetscMat, B::AbstractPetscMat, X::AbstractPetscMat) end

@for_petsc function MatMatSolve(petsclib::$UnionPetscLib, A::AbstractPetscMat, B::AbstractPetscMat, X::AbstractPetscMat )

    @chk ccall(
               (:MatMatSolve, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat),
               A, B, X,
              )


	return nothing
end 

"""
	MatMatSolveTranspose(petsclib::PetscLibType,A::AbstractPetscMat, B::AbstractPetscMat, X::AbstractPetscMat) 
Solves A^T X = B , given a factored matrix.

Neighbor-wise Collective

Input Parameters:
- `A` - the factored matrix
- `B` - the right-hand-side matrix  (`MATDENSE` matrix)

Output Parameter:
- `X` - the result matrix (dense matrix)

Level: developer

-seealso: [](ch_matrices), `Mat`, [Matrix Factorization](sec_matfactor), `MatGetFactor()`, `MatSolveTranspose()`, `MatMatSolve()`, `MatLUFactor()`, `MatCholeskyFactor()`

# External Links
$(_doc_external("Mat/MatMatSolveTranspose"))
"""
function MatMatSolveTranspose(petsclib::PetscLibType, A::AbstractPetscMat, B::AbstractPetscMat, X::AbstractPetscMat) end

@for_petsc function MatMatSolveTranspose(petsclib::$UnionPetscLib, A::AbstractPetscMat, B::AbstractPetscMat, X::AbstractPetscMat )

    @chk ccall(
               (:MatMatSolveTranspose, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat),
               A, B, X,
              )


	return nothing
end 

"""
	C::PetscMat = MatMatTransposeMult(petsclib::PetscLibType,A::AbstractPetscMat, B::AbstractPetscMat, scall::MatReuse, fill::PetscReal) 
Performs matrix

Neighbor-wise Collective

Input Parameters:
- `A`     - the left matrix
- `B`     - the right matrix
- `scall` - either `MAT_INITIAL_MATRIX` or `MAT_REUSE_MATRIX`
- `fill`  - expected fill as ratio of nnz(C)/(nnz(A) + nnz(B)), use `PETSC_DETERMINE` or `PETSC_CURRENT` if not known

Output Parameter:
- `C` - the product matrix

Options Database Key:
- `-matmattransmult_mpidense_mpidense_via {allgatherv,cyclic}` - Choose between algorithms for `MATMPIDENSE` matrices: the
first redundantly copies the transposed `B` matrix on each process and requires O(log P) communication complexity;
the second never stores more than one portion of the `B` matrix at a time but requires O(P) communication complexity.

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatProductCreate()`, `MATPRODUCT_ABt`, `MatMatMult()`, `MatTransposeMatMult()` `MatPtAP()`, `MatProductAlgorithm`, `MatProductType`

# External Links
$(_doc_external("Mat/MatMatTransposeMult"))
"""
function MatMatTransposeMult(petsclib::PetscLibType, A::AbstractPetscMat, B::AbstractPetscMat, scall::MatReuse, fill::PetscReal) end

@for_petsc function MatMatTransposeMult(petsclib::$UnionPetscLib, A::AbstractPetscMat, B::AbstractPetscMat, scall::MatReuse, fill::$PetscReal )
	C_ = Ref{CMat}()

    @chk ccall(
               (:MatMatTransposeMult, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, MatReuse, $PetscReal, Ptr{CMat}),
               A, B, scall, fill, C_,
              )

	C = PetscMat(C_[], petsclib)

	return C
end 

"""
	flg::PetscBool = MatMatTransposeMultEqual(petsclib::PetscLibType,A::AbstractPetscMat, B::AbstractPetscMat, C::AbstractPetscMat, n::PetscInt) 
Test A*B^T*x = C*x for n random vector x

Collective

Input Parameters:
- `A` - the first matrix
- `B` - the second matrix
- `C` - the third matrix
- `n` - number of random vectors to be tested

Output Parameter:
- `flg` - `PETSC_TRUE` if the products are equal; `PETSC_FALSE` otherwise.

Level: intermediate

-seealso: `Mat`, `MatMatMultEqual()`, `MatMultEqual()`, `MatMultAddEqual()`, `MatMultTransposeEqual()`

# External Links
$(_doc_external("Mat/MatMatTransposeMultEqual"))
"""
function MatMatTransposeMultEqual(petsclib::PetscLibType, A::AbstractPetscMat, B::AbstractPetscMat, C::AbstractPetscMat, n::PetscInt) end

@for_petsc function MatMatTransposeMultEqual(petsclib::$UnionPetscLib, A::AbstractPetscMat, B::AbstractPetscMat, C::AbstractPetscMat, n::$PetscInt )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatMatTransposeMultEqual, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat, $PetscInt, Ptr{PetscBool}),
               A, B, C, n, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	MatMatTransposeSolve(petsclib::PetscLibType,A::AbstractPetscMat, Bt::AbstractPetscMat, X::AbstractPetscMat) 
Solves A X = B^T, given a factored matrix.

Neighbor-wise Collective

Input Parameters:
- `A`  - the factored matrix
- `Bt` - the transpose of right-hand-side matrix as a `MATDENSE`

Output Parameter:
- `X` - the result matrix (dense matrix)

Level: developer

-seealso: [](ch_matrices), `Mat`, [Matrix Factorization](sec_matfactor), `MatMatSolve()`, `MatMatSolveTranspose()`, `MatLUFactor()`, `MatCholeskyFactor()`

# External Links
$(_doc_external("Mat/MatMatTransposeSolve"))
"""
function MatMatTransposeSolve(petsclib::PetscLibType, A::AbstractPetscMat, Bt::AbstractPetscMat, X::AbstractPetscMat) end

@for_petsc function MatMatTransposeSolve(petsclib::$UnionPetscLib, A::AbstractPetscMat, Bt::AbstractPetscMat, X::AbstractPetscMat )

    @chk ccall(
               (:MatMatTransposeSolve, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat),
               A, Bt, X,
              )


	return nothing
end 

"""
	dual::PetscMat = MatMeshToCellGraph(petsclib::PetscLibType,mesh::AbstractPetscMat, ncommonnodes::PetscInt) 
Convert a mesh to a cell graph.

Collective

Input Parameters:
- `mesh`         - the graph that represents the coupling of the vertices of the mesh
- `ncommonnodes` - mesh elements that share this number of common nodes are considered neighbors, use 2 for triangles and
quadrilaterials, 3 for tetrahedrals and 4 for hexahedrals

Output Parameter:
- `dual` - the dual graph

Level: advanced

-seealso: `MatCreateMPIAdj()`, `MatPartitioningCreate()`

# External Links
$(_doc_external("MatGraphOperations/MatMeshToCellGraph"))
"""
function MatMeshToCellGraph(petsclib::PetscLibType, mesh::AbstractPetscMat, ncommonnodes::PetscInt) end

@for_petsc function MatMeshToCellGraph(petsclib::$UnionPetscLib, mesh::AbstractPetscMat, ncommonnodes::$PetscInt )
	dual_ = Ref{CMat}()

    @chk ccall(
               (:MatMeshToCellGraph, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CMat}),
               mesh, ncommonnodes, dual_,
              )

	dual = PetscMat(dual_[], petsclib)

	return dual
end 

"""
	missing::PetscBool,dd::PetscInt = MatMissingDiagonal(petsclib::PetscLibType,mat::AbstractPetscMat) 
Determine if sparse matrix is missing a diagonal entry (or block entry for `MATBAIJ` and `MATSBAIJ` matrices) in the nonzero structure

Not Collective

Input Parameter:
- `mat` - the matrix

Output Parameters:
- `missing` - is any diagonal entry missing
- `dd`      - first diagonal entry that is missing (optional) on this process

Level: advanced

-seealso: [](ch_matrices), `Mat`

# External Links
$(_doc_external("Mat/MatMissingDiagonal"))
"""
function MatMissingDiagonal(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatMissingDiagonal(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	missing_ = Ref{PetscBool}()
	dd_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatMissingDiagonal, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscBool}, Ptr{$PetscInt}),
               mat, missing_, dd_,
              )

	missing = missing_[]
	dd = dd_[]

	return missing,dd
end 

"""
	MatMult(petsclib::PetscLibType,mat::AbstractPetscMat, x::AbstractPetscVec, y::AbstractPetscVec) 
Computes the matrix

Neighbor-wise Collective

Input Parameters:
- `mat` - the matrix
- `x`   - the vector to be multiplied

Output Parameter:
- `y` - the result

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatMultTranspose()`, `MatMultAdd()`, `MatMultTransposeAdd()`

# External Links
$(_doc_external("Mat/MatMult"))
"""
function MatMult(petsclib::PetscLibType, mat::AbstractPetscMat, x::AbstractPetscVec, y::AbstractPetscVec) end

@for_petsc function MatMult(petsclib::$UnionPetscLib, mat::AbstractPetscMat, x::AbstractPetscVec, y::AbstractPetscVec )

    @chk ccall(
               (:MatMult, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               mat, x, y,
              )


	return nothing
end 

"""
	MatMultAdd(petsclib::PetscLibType,mat::AbstractPetscMat, v1::AbstractPetscVec, v2::AbstractPetscVec, v3::AbstractPetscVec) 
Computes v3 = v2 + A * v1.

Neighbor-wise Collective

Input Parameters:
- `mat` - the matrix
- `v1`  - the vector to be multiplied by `mat`
- `v2`  - the vector to be added to the result

Output Parameter:
- `v3` - the result

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatMultTranspose()`, `MatMult()`, `MatMultTransposeAdd()`

# External Links
$(_doc_external("Mat/MatMultAdd"))
"""
function MatMultAdd(petsclib::PetscLibType, mat::AbstractPetscMat, v1::AbstractPetscVec, v2::AbstractPetscVec, v3::AbstractPetscVec) end

@for_petsc function MatMultAdd(petsclib::$UnionPetscLib, mat::AbstractPetscMat, v1::AbstractPetscVec, v2::AbstractPetscVec, v3::AbstractPetscVec )

    @chk ccall(
               (:MatMultAdd, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec, CVec),
               mat, v1, v2, v3,
              )


	return nothing
end 

"""
	flg::PetscBool = MatMultAddEqual(petsclib::PetscLibType,A::AbstractPetscMat, B::AbstractPetscMat, n::PetscInt) 
Compares matrix

Collective

Input Parameters:
- `A` - the first matrix
- `B` - the second matrix
- `n` - number of random vectors to be tested

Output Parameter:
- `flg` - `PETSC_TRUE` if the products are equal; `PETSC_FALSE` otherwise.

Level: intermediate

-seealso: `Mat`, `MatMultEqual()`, `MatMultTransposeEqual()`, `MatMultTransposeAddEqual()`

# External Links
$(_doc_external("Mat/MatMultAddEqual"))
"""
function MatMultAddEqual(petsclib::PetscLibType, A::AbstractPetscMat, B::AbstractPetscMat, n::PetscInt) end

@for_petsc function MatMultAddEqual(petsclib::$UnionPetscLib, A::AbstractPetscMat, B::AbstractPetscMat, n::$PetscInt )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatMultAddEqual, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, $PetscInt, Ptr{PetscBool}),
               A, B, n, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	MatMultDiagonalBlock(petsclib::PetscLibType,mat::AbstractPetscMat, x::AbstractPetscVec, y::AbstractPetscVec) 
Computes the matrix

Collective

Input Parameters:
- `mat` - the matrix
- `x`   - the vector to be multiplied

Output Parameter:
- `y` - the result

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatMult()`, `MatMultTranspose()`, `MatMultAdd()`, `MatMultTransposeAdd()`

# External Links
$(_doc_external("Mat/MatMultDiagonalBlock"))
"""
function MatMultDiagonalBlock(petsclib::PetscLibType, mat::AbstractPetscMat, x::AbstractPetscVec, y::AbstractPetscVec) end

@for_petsc function MatMultDiagonalBlock(petsclib::$UnionPetscLib, mat::AbstractPetscMat, x::AbstractPetscVec, y::AbstractPetscVec )

    @chk ccall(
               (:MatMultDiagonalBlock, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               mat, x, y,
              )


	return nothing
end 

"""
	flg::PetscBool = MatMultEqual(petsclib::PetscLibType,A::AbstractPetscMat, B::AbstractPetscMat, n::PetscInt) 
Compares matrix

Collective

Input Parameters:
- `A` - the first matrix
- `B` - the second matrix
- `n` - number of random vectors to be tested

Output Parameter:
- `flg` - `PETSC_TRUE` if the products are equal; `PETSC_FALSE` otherwise.

Level: intermediate

-seealso: `Mat`, `MatMultAddEqual()`, `MatMultTransposeEqual()`, `MatMultTransposeAddEqual()`, `MatIsLinear()`, `MatEqual()`

# External Links
$(_doc_external("Mat/MatMultEqual"))
"""
function MatMultEqual(petsclib::PetscLibType, A::AbstractPetscMat, B::AbstractPetscMat, n::PetscInt) end

@for_petsc function MatMultEqual(petsclib::$UnionPetscLib, A::AbstractPetscMat, B::AbstractPetscMat, n::$PetscInt )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatMultEqual, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, $PetscInt, Ptr{PetscBool}),
               A, B, n, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	MatMultHermitianTranspose(petsclib::PetscLibType,mat::AbstractPetscMat, x::AbstractPetscVec, y::AbstractPetscVec) 
Computes matrix Hermitian

Neighbor-wise Collective

Input Parameters:
- `mat` - the matrix
- `x`   - the vector to be multiplied

Output Parameter:
- `y` - the result

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatMult()`, `MatMultAdd()`, `MatMultHermitianTransposeAdd()`, `MatMultTranspose()`

# External Links
$(_doc_external("Mat/MatMultHermitianTranspose"))
"""
function MatMultHermitianTranspose(petsclib::PetscLibType, mat::AbstractPetscMat, x::AbstractPetscVec, y::AbstractPetscVec) end

@for_petsc function MatMultHermitianTranspose(petsclib::$UnionPetscLib, mat::AbstractPetscMat, x::AbstractPetscVec, y::AbstractPetscVec )

    @chk ccall(
               (:MatMultHermitianTranspose, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               mat, x, y,
              )


	return nothing
end 

"""
	MatMultHermitianTransposeAdd(petsclib::PetscLibType,mat::AbstractPetscMat, v1::AbstractPetscVec, v2::AbstractPetscVec, v3::AbstractPetscVec) 
Computes v3 = v2 + A^H * v1.

Neighbor-wise Collective

Input Parameters:
- `mat` - the matrix
- `v1`  - the vector to be multiplied by the Hermitian transpose
- `v2`  - the vector to be added to the result

Output Parameter:
- `v3` - the result

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatMultHermitianTranspose()`, `MatMultTranspose()`, `MatMultAdd()`, `MatMult()`

# External Links
$(_doc_external("Mat/MatMultHermitianTransposeAdd"))
"""
function MatMultHermitianTransposeAdd(petsclib::PetscLibType, mat::AbstractPetscMat, v1::AbstractPetscVec, v2::AbstractPetscVec, v3::AbstractPetscVec) end

@for_petsc function MatMultHermitianTransposeAdd(petsclib::$UnionPetscLib, mat::AbstractPetscMat, v1::AbstractPetscVec, v2::AbstractPetscVec, v3::AbstractPetscVec )

    @chk ccall(
               (:MatMultHermitianTransposeAdd, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec, CVec),
               mat, v1, v2, v3,
              )


	return nothing
end 

"""
	flg::PetscBool = MatMultHermitianTransposeAddEqual(petsclib::PetscLibType,A::AbstractPetscMat, B::AbstractPetscMat, n::PetscInt) 
Compares matrix

Collective

Input Parameters:
- `A` - the first matrix
- `B` - the second matrix
- `n` - number of random vectors to be tested

Output Parameter:
- `flg` - `PETSC_TRUE` if the products are equal; `PETSC_FALSE` otherwise.

Level: intermediate

-seealso: `Mat`, `MatMatMultEqual()`, `MatMultEqual()`, `MatMultAddEqual()`, `MatMultTransposeEqual()`

# External Links
$(_doc_external("Mat/MatMultHermitianTransposeAddEqual"))
"""
function MatMultHermitianTransposeAddEqual(petsclib::PetscLibType, A::AbstractPetscMat, B::AbstractPetscMat, n::PetscInt) end

@for_petsc function MatMultHermitianTransposeAddEqual(petsclib::$UnionPetscLib, A::AbstractPetscMat, B::AbstractPetscMat, n::$PetscInt )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatMultHermitianTransposeAddEqual, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, $PetscInt, Ptr{PetscBool}),
               A, B, n, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	flg::PetscBool = MatMultHermitianTransposeEqual(petsclib::PetscLibType,A::AbstractPetscMat, B::AbstractPetscMat, n::PetscInt) 
Compares matrix

Collective

Input Parameters:
- `A` - the first matrix
- `B` - the second matrix
- `n` - number of random vectors to be tested

Output Parameter:
- `flg` - `PETSC_TRUE` if the products are equal; `PETSC_FALSE` otherwise.

Level: intermediate

-seealso: `Mat`, `MatMatMultEqual()`, `MatMultEqual()`, `MatMultAddEqual()`, `MatMultTransposeEqual()`

# External Links
$(_doc_external("Mat/MatMultHermitianTransposeEqual"))
"""
function MatMultHermitianTransposeEqual(petsclib::PetscLibType, A::AbstractPetscMat, B::AbstractPetscMat, n::PetscInt) end

@for_petsc function MatMultHermitianTransposeEqual(petsclib::$UnionPetscLib, A::AbstractPetscMat, B::AbstractPetscMat, n::$PetscInt )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatMultHermitianTransposeEqual, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, $PetscInt, Ptr{PetscBool}),
               A, B, n, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	MatMultTranspose(petsclib::PetscLibType,mat::AbstractPetscMat, x::AbstractPetscVec, y::AbstractPetscVec) 
Computes matrix transpose times a vector y = A^T * x.

Neighbor-wise Collective

Input Parameters:
- `mat` - the matrix
- `x`   - the vector to be multiplied

Output Parameter:
- `y` - the result

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatMult()`, `MatMultAdd()`, `MatMultTransposeAdd()`, `MatMultHermitianTranspose()`, `MatTranspose()`

# External Links
$(_doc_external("Mat/MatMultTranspose"))
"""
function MatMultTranspose(petsclib::PetscLibType, mat::AbstractPetscMat, x::AbstractPetscVec, y::AbstractPetscVec) end

@for_petsc function MatMultTranspose(petsclib::$UnionPetscLib, mat::AbstractPetscMat, x::AbstractPetscVec, y::AbstractPetscVec )

    @chk ccall(
               (:MatMultTranspose, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               mat, x, y,
              )


	return nothing
end 

"""
	MatMultTransposeAdd(petsclib::PetscLibType,mat::AbstractPetscMat, v1::AbstractPetscVec, v2::AbstractPetscVec, v3::AbstractPetscVec) 
Computes v3 = v2 + A^T * v1.

Neighbor-wise Collective

Input Parameters:
- `mat` - the matrix
- `v1`  - the vector to be multiplied by the transpose of the matrix
- `v2`  - the vector to be added to the result

Output Parameter:
- `v3` - the result

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatMultTranspose()`, `MatMultAdd()`, `MatMult()`

# External Links
$(_doc_external("Mat/MatMultTransposeAdd"))
"""
function MatMultTransposeAdd(petsclib::PetscLibType, mat::AbstractPetscMat, v1::AbstractPetscVec, v2::AbstractPetscVec, v3::AbstractPetscVec) end

@for_petsc function MatMultTransposeAdd(petsclib::$UnionPetscLib, mat::AbstractPetscMat, v1::AbstractPetscVec, v2::AbstractPetscVec, v3::AbstractPetscVec )

    @chk ccall(
               (:MatMultTransposeAdd, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec, CVec),
               mat, v1, v2, v3,
              )


	return nothing
end 

"""
	flg::PetscBool = MatMultTransposeAddEqual(petsclib::PetscLibType,A::AbstractPetscMat, B::AbstractPetscMat, n::PetscInt) 
Compares matrix

Collective

Input Parameters:
- `A` - the first matrix
- `B` - the second matrix
- `n` - number of random vectors to be tested

Output Parameter:
- `flg` - `PETSC_TRUE` if the products are equal; `PETSC_FALSE` otherwise.

Level: intermediate

-seealso: `Mat`, `MatMultEqual()`, `MatMultAddEqual()`, `MatMultTransposeEqual()`

# External Links
$(_doc_external("Mat/MatMultTransposeAddEqual"))
"""
function MatMultTransposeAddEqual(petsclib::PetscLibType, A::AbstractPetscMat, B::AbstractPetscMat, n::PetscInt) end

@for_petsc function MatMultTransposeAddEqual(petsclib::$UnionPetscLib, A::AbstractPetscMat, B::AbstractPetscMat, n::$PetscInt )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatMultTransposeAddEqual, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, $PetscInt, Ptr{PetscBool}),
               A, B, n, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	flg::PetscBool = MatMultTransposeEqual(petsclib::PetscLibType,A::AbstractPetscMat, B::AbstractPetscMat, n::PetscInt) 
Compares matrix

Collective

Input Parameters:
- `A` - the first matrix
- `B` - the second matrix
- `n` - number of random vectors to be tested

Output Parameter:
- `flg` - `PETSC_TRUE` if the products are equal; `PETSC_FALSE` otherwise.

Level: intermediate

-seealso: `Mat`, `MatMultEqual()`, `MatMultAddEqual()`, `MatMultTransposeAddEqual()`

# External Links
$(_doc_external("Mat/MatMultTransposeEqual"))
"""
function MatMultTransposeEqual(petsclib::PetscLibType, A::AbstractPetscMat, B::AbstractPetscMat, n::PetscInt) end

@for_petsc function MatMultTransposeEqual(petsclib::$UnionPetscLib, A::AbstractPetscMat, B::AbstractPetscMat, n::$PetscInt )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatMultTransposeEqual, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, $PetscInt, Ptr{PetscBool}),
               A, B, n, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	val::PetscReal = MatMumpsGetCntl(petsclib::PetscLibType,F::AbstractPetscMat, icntl::PetscInt) 
Get MUMPS parameter CNTL() <https://mumps

Logically Collective

Input Parameters:
- `F`     - the factored matrix obtained by calling `MatGetFactor()` with a `MatSolverType` of `MATSOLVERMUMPS` and a `MatFactorType` of `MAT_FACTOR_LU` or `MAT_FACTOR_CHOLESKY`
- `icntl` - index of MUMPS parameter array CNTL()

Output Parameter:
- `val` - value of MUMPS CNTL(icntl)

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatGetFactor()`, `MatMumpsSetIcntl()`, `MatMumpsGetIcntl()`, `MatMumpsSetCntl()`, `MatMumpsGetInfo()`, `MatMumpsGetInfog()`, `MatMumpsGetRinfo()`, `MatMumpsGetRinfog()`

# External Links
$(_doc_external("Mat/MatMumpsGetCntl"))
"""
function MatMumpsGetCntl(petsclib::PetscLibType, F::AbstractPetscMat, icntl::PetscInt) end

@for_petsc function MatMumpsGetCntl(petsclib::$UnionPetscLib, F::AbstractPetscMat, icntl::$PetscInt )
	val_ = Ref{$PetscReal}()

    @chk ccall(
               (:MatMumpsGetCntl, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscReal}),
               F, icntl, val_,
              )

	val = val_[]

	return val
end 

"""
	ival::PetscInt = MatMumpsGetIcntl(petsclib::PetscLibType,F::AbstractPetscMat, icntl::PetscInt) 
Get MUMPS parameter ICNTL() <https://mumps

Logically Collective

Input Parameters:
- `F`     - the factored matrix obtained by calling `MatGetFactor()` with a `MatSolverType` of `MATSOLVERMUMPS` and a `MatFactorType` of `MAT_FACTOR_LU` or `MAT_FACTOR_CHOLESKY`
- `icntl` - index of MUMPS parameter array ICNTL()

Output Parameter:
- `ival` - value of MUMPS ICNTL(icntl)

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatGetFactor()`, `MatMumpsSetIcntl()`, `MatMumpsSetCntl()`, `MatMumpsGetCntl()`, `MatMumpsGetInfo()`, `MatMumpsGetInfog()`, `MatMumpsGetRinfo()`, `MatMumpsGetRinfog()`

# External Links
$(_doc_external("Mat/MatMumpsGetIcntl"))
"""
function MatMumpsGetIcntl(petsclib::PetscLibType, F::AbstractPetscMat, icntl::PetscInt) end

@for_petsc function MatMumpsGetIcntl(petsclib::$UnionPetscLib, F::AbstractPetscMat, icntl::$PetscInt )
	ival_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatMumpsGetIcntl, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}),
               F, icntl, ival_,
              )

	ival = ival_[]

	return ival
end 

"""
	ival::PetscInt = MatMumpsGetInfo(petsclib::PetscLibType,F::AbstractPetscMat, icntl::PetscInt) 
Get MUMPS parameter INFO() <https://mumps

Logically Collective

Input Parameters:
- `F`     - the factored matrix obtained by calling `MatGetFactor()` with a `MatSolverType` of `MATSOLVERMUMPS` and a `MatFactorType` of `MAT_FACTOR_LU` or `MAT_FACTOR_CHOLESKY`
- `icntl` - index of MUMPS parameter array INFO()

Output Parameter:
- `ival` - value of MUMPS INFO(icntl)

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatGetFactor()`, `MatMumpsSetIcntl()`, `MatMumpsGetIcntl()`, `MatMumpsSetCntl()`, `MatMumpsGetCntl()`, `MatMumpsGetInfog()`, `MatMumpsGetRinfo()`, `MatMumpsGetRinfog()`

# External Links
$(_doc_external("Mat/MatMumpsGetInfo"))
"""
function MatMumpsGetInfo(petsclib::PetscLibType, F::AbstractPetscMat, icntl::PetscInt) end

@for_petsc function MatMumpsGetInfo(petsclib::$UnionPetscLib, F::AbstractPetscMat, icntl::$PetscInt )
	ival_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatMumpsGetInfo, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}),
               F, icntl, ival_,
              )

	ival = ival_[]

	return ival
end 

"""
	ival::PetscInt = MatMumpsGetInfog(petsclib::PetscLibType,F::AbstractPetscMat, icntl::PetscInt) 
Get MUMPS parameter INFOG() <https://mumps

Logically Collective

Input Parameters:
- `F`     - the factored matrix obtained by calling `MatGetFactor()` with a `MatSolverType` of `MATSOLVERMUMPS` and a `MatFactorType` of `MAT_FACTOR_LU` or `MAT_FACTOR_CHOLESKY`
- `icntl` - index of MUMPS parameter array INFOG()

Output Parameter:
- `ival` - value of MUMPS INFOG(icntl)

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatGetFactor()`, `MatMumpsSetIcntl()`, `MatMumpsGetIcntl()`, `MatMumpsSetCntl()`, `MatMumpsGetCntl()`, `MatMumpsGetInfo()`, `MatMumpsGetRinfo()`, `MatMumpsGetRinfog()`

# External Links
$(_doc_external("Mat/MatMumpsGetInfog"))
"""
function MatMumpsGetInfog(petsclib::PetscLibType, F::AbstractPetscMat, icntl::PetscInt) end

@for_petsc function MatMumpsGetInfog(petsclib::$UnionPetscLib, F::AbstractPetscMat, icntl::$PetscInt )
	ival_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatMumpsGetInfog, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}),
               F, icntl, ival_,
              )

	ival = ival_[]

	return ival
end 

"""
	MatMumpsGetInverse(petsclib::PetscLibType,F::AbstractPetscMat, spRHS::AbstractPetscMat) 
Get user

Logically Collective

Input Parameter:
- `F` - the factored matrix obtained by calling `MatGetFactor()` with a `MatSolverType` of `MATSOLVERMUMPS` and a `MatFactorType` of `MAT_FACTOR_LU` or `MAT_FACTOR_CHOLESKY`

Output Parameter:
- `spRHS` - sequential sparse matrix in `MATTRANSPOSEVIRTUAL` format with requested entries of inverse of `A`

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatGetFactor()`, `MatCreateTranspose()`

# External Links
$(_doc_external("Mat/MatMumpsGetInverse"))
"""
function MatMumpsGetInverse(petsclib::PetscLibType, F::AbstractPetscMat, spRHS::AbstractPetscMat) end

@for_petsc function MatMumpsGetInverse(petsclib::$UnionPetscLib, F::AbstractPetscMat, spRHS::AbstractPetscMat )

    @chk ccall(
               (:MatMumpsGetInverse, $petsc_library),
               PetscErrorCode,
               (CMat, CMat),
               F, spRHS,
              )


	return nothing
end 

"""
	MatMumpsGetInverseTranspose(petsclib::PetscLibType,F::AbstractPetscMat, spRHST::AbstractPetscMat) 
Get user

Logically Collective

Input Parameter:
- `F` - the factored matrix of A obtained by calling `MatGetFactor()` with a `MatSolverType` of `MATSOLVERMUMPS` and a `MatFactorType` of `MAT_FACTOR_LU` or `MAT_FACTOR_CHOLESKY`

Output Parameter:
- `spRHST` - sequential sparse matrix in `MATAIJ` format containing the requested entries of inverse of `A`^T

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatGetFactor()`, `MatCreateTranspose()`, `MatMumpsGetInverse()`

# External Links
$(_doc_external("Mat/MatMumpsGetInverseTranspose"))
"""
function MatMumpsGetInverseTranspose(petsclib::PetscLibType, F::AbstractPetscMat, spRHST::AbstractPetscMat) end

@for_petsc function MatMumpsGetInverseTranspose(petsclib::$UnionPetscLib, F::AbstractPetscMat, spRHST::AbstractPetscMat )

    @chk ccall(
               (:MatMumpsGetInverseTranspose, $petsc_library),
               PetscErrorCode,
               (CMat, CMat),
               F, spRHST,
              )


	return nothing
end 

"""
	size::PetscInt,array::Ptr{PetscInt} = MatMumpsGetNullPivots(petsclib::PetscLibType,F::AbstractPetscMat) 
Get MUMPS parameter PIVNUL_LIST() <https://mumps

Logically Collective

Input Parameter:
- `F` - the factored matrix obtained by calling `MatGetFactor()` with a `MatSolverType` of `MATSOLVERMUMPS` and a `MatFactorType` of `MAT_FACTOR_LU` or `MAT_FACTOR_CHOLESKY`

Output Parameters:
- `size`  - local size of the array. The size of the array is non-zero only on MPI rank 0
- `array` - array of rows with null pivot, these rows follow 0-based indexing. The array gets allocated within the function and the user is responsible
for freeing this array.

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatGetFactor()`, `MatMumpsSetIcntl()`, `MatMumpsGetIcntl()`, `MatMumpsSetCntl()`, `MatMumpsGetCntl()`, `MatMumpsGetInfo()`, `MatMumpsGetInfog()`, `MatMumpsGetRinfo()`

# External Links
$(_doc_external("Mat/MatMumpsGetNullPivots"))
"""
function MatMumpsGetNullPivots(petsclib::PetscLibType, F::AbstractPetscMat) end

@for_petsc function MatMumpsGetNullPivots(petsclib::$UnionPetscLib, F::AbstractPetscMat )
	size_ = Ref{$PetscInt}()
	array_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:MatMumpsGetNullPivots, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
               F, size_, array_,
              )

	size = size_[]
	array = array_[]

	return size,array
end 

"""
	val::PetscReal = MatMumpsGetRinfo(petsclib::PetscLibType,F::AbstractPetscMat, icntl::PetscInt) 
Get MUMPS parameter RINFO() <https://mumps

Logically Collective

Input Parameters:
- `F`     - the factored matrix obtained by calling `MatGetFactor()` with a `MatSolverType` of `MATSOLVERMUMPS` and a `MatFactorType` of `MAT_FACTOR_LU` or `MAT_FACTOR_CHOLESKY`
- `icntl` - index of MUMPS parameter array RINFO()

Output Parameter:
- `val` - value of MUMPS RINFO(icntl)

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatGetFactor()`, `MatMumpsSetIcntl()`, `MatMumpsGetIcntl()`, `MatMumpsSetCntl()`, `MatMumpsGetCntl()`, `MatMumpsGetInfo()`, `MatMumpsGetInfog()`, `MatMumpsGetRinfog()`

# External Links
$(_doc_external("Mat/MatMumpsGetRinfo"))
"""
function MatMumpsGetRinfo(petsclib::PetscLibType, F::AbstractPetscMat, icntl::PetscInt) end

@for_petsc function MatMumpsGetRinfo(petsclib::$UnionPetscLib, F::AbstractPetscMat, icntl::$PetscInt )
	val_ = Ref{$PetscReal}()

    @chk ccall(
               (:MatMumpsGetRinfo, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscReal}),
               F, icntl, val_,
              )

	val = val_[]

	return val
end 

"""
	val::PetscReal = MatMumpsGetRinfog(petsclib::PetscLibType,F::AbstractPetscMat, icntl::PetscInt) 
Get MUMPS parameter RINFOG() <https://mumps

Logically Collective

Input Parameters:
- `F`     - the factored matrix obtained by calling `MatGetFactor()` with a `MatSolverType` of `MATSOLVERMUMPS` and a `MatFactorType` of `MAT_FACTOR_LU` or `MAT_FACTOR_CHOLESKY`
- `icntl` - index of MUMPS parameter array RINFOG()

Output Parameter:
- `val` - value of MUMPS RINFOG(icntl)

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatGetFactor()`, `MatMumpsSetIcntl()`, `MatMumpsGetIcntl()`, `MatMumpsSetCntl()`, `MatMumpsGetCntl()`, `MatMumpsGetInfo()`, `MatMumpsGetInfog()`, `MatMumpsGetRinfo()`

# External Links
$(_doc_external("Mat/MatMumpsGetRinfog"))
"""
function MatMumpsGetRinfog(petsclib::PetscLibType, F::AbstractPetscMat, icntl::PetscInt) end

@for_petsc function MatMumpsGetRinfog(petsclib::$UnionPetscLib, F::AbstractPetscMat, icntl::$PetscInt )
	val_ = Ref{$PetscReal}()

    @chk ccall(
               (:MatMumpsGetRinfog, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscReal}),
               F, icntl, val_,
              )

	val = val_[]

	return val
end 

"""
	MatMumpsSetBlk(petsclib::PetscLibType,F::AbstractPetscMat, nblk::PetscInt, blkvar::Vector{PetscInt}, blkptr::Vector{PetscInt}) 
Set user

Not collective, only relevant on the first process of the MPI communicator

Input Parameters:
- `F`      - the factored matrix of A obtained by calling `MatGetFactor()` with a `MatSolverType` of `MATSOLVERMUMPS` and a `MatFactorType` of `MAT_FACTOR_LU` or `MAT_FACTOR_CHOLESKY`
- `nblk`   - the number of blocks
- `blkvar` - see MUMPS documentation, `blkvar(blkptr(iblk):blkptr(iblk+1)-1)`, (`iblk=1, nblk`) holds the variables associated to block `iblk`
- `blkptr` - array starting at 1 and of size `nblk + 1` storing the prefix sum of all blocks

Level: advanced

-seealso: [](ch_matrices), `MATSOLVERMUMPS`, `Mat`, `MatGetFactor()`, `MatMumpsSetIcntl()`, `MatSetVariableBlockSizes()`

# External Links
$(_doc_external("Mat/MatMumpsSetBlk"))
"""
function MatMumpsSetBlk(petsclib::PetscLibType, F::AbstractPetscMat, nblk::PetscInt, blkvar::Vector{PetscInt}, blkptr::Vector{PetscInt}) end

@for_petsc function MatMumpsSetBlk(petsclib::$UnionPetscLib, F::AbstractPetscMat, nblk::$PetscInt, blkvar::Vector{$PetscInt}, blkptr::Vector{$PetscInt} )

    @chk ccall(
               (:MatMumpsSetBlk, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
               F, nblk, blkvar, blkptr,
              )


	return nothing
end 

"""
	MatMumpsSetCntl(petsclib::PetscLibType,F::AbstractPetscMat, icntl::PetscInt, val::PetscReal) 
Set MUMPS parameter CNTL() <https://mumps

Logically Collective

Input Parameters:
- `F`     - the factored matrix obtained by calling `MatGetFactor()` with a `MatSolverType` of `MATSOLVERMUMPS` and a `MatFactorType` of `MAT_FACTOR_LU` or `MAT_FACTOR_CHOLESKY`
- `icntl` - index of MUMPS parameter array `CNTL()`
- `val`   - value of MUMPS `CNTL(icntl)`

Options Database Key:
- `-mat_mumps_cntl_<icntl> <val>` - change the option numbered icntl to ival

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatGetFactor()`, `MatMumpsSetIcntl()`, `MatMumpsGetIcntl()`, `MatMumpsGetCntl()`, `MatMumpsGetInfo()`, `MatMumpsGetInfog()`, `MatMumpsGetRinfo()`, `MatMumpsGetRinfog()`

# External Links
$(_doc_external("Mat/MatMumpsSetCntl"))
"""
function MatMumpsSetCntl(petsclib::PetscLibType, F::AbstractPetscMat, icntl::PetscInt, val::PetscReal) end

@for_petsc function MatMumpsSetCntl(petsclib::$UnionPetscLib, F::AbstractPetscMat, icntl::$PetscInt, val::$PetscReal )

    @chk ccall(
               (:MatMumpsSetCntl, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscReal),
               F, icntl, val,
              )


	return nothing
end 

"""
	MatMumpsSetIcntl(petsclib::PetscLibType,F::AbstractPetscMat, icntl::PetscInt, ival::PetscInt) 
Set MUMPS parameter ICNTL() <https://mumps

Logically Collective

Input Parameters:
- `F`     - the factored matrix obtained by calling `MatGetFactor()` with a `MatSolverType` of `MATSOLVERMUMPS` and a `MatFactorType` of `MAT_FACTOR_LU` or `MAT_FACTOR_CHOLESKY`
- `icntl` - index of MUMPS parameter array `ICNTL()`
- `ival`  - value of MUMPS `ICNTL(icntl)`

Options Database Key:
- `-mat_mumps_icntl_<icntl> <ival>` - change the option numbered `icntl` to `ival`

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatGetFactor()`, `MatMumpsGetIcntl()`, `MatMumpsSetCntl()`, `MatMumpsGetCntl()`, `MatMumpsGetInfo()`, `MatMumpsGetInfog()`, `MatMumpsGetRinfo()`, `MatMumpsGetRinfog()`

# External Links
$(_doc_external("Mat/MatMumpsSetIcntl"))
"""
function MatMumpsSetIcntl(petsclib::PetscLibType, F::AbstractPetscMat, icntl::PetscInt, ival::PetscInt) end

@for_petsc function MatMumpsSetIcntl(petsclib::$UnionPetscLib, F::AbstractPetscMat, icntl::$PetscInt, ival::$PetscInt )

    @chk ccall(
               (:MatMumpsSetIcntl, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt),
               F, icntl, ival,
              )


	return nothing
end 

"""
	MatNestGetISs(petsclib::PetscLibType,A::AbstractPetscMat, rows::Vector{<:AbstractIS}, cols::Vector{<:AbstractIS}) 
Returns the index sets partitioning the row and column spaces of a `MATNEST`

Not Collective

Input Parameter:
- `A` - `MATNEST` matrix

Output Parameters:
- `rows` - array of row index sets (pass `NULL` to ignore)
- `cols` - array of column index sets (pass `NULL` to ignore)

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATNEST`, `MatNestGetSubMat()`, `MatNestGetSubMats()`, `MatNestGetSize()`, `MatNestGetLocalISs()`,
`MatCreateNest()`, `MatNestSetSubMats()`

# External Links
$(_doc_external("Mat/MatNestGetISs"))
"""
function MatNestGetISs(petsclib::PetscLibType, A::AbstractPetscMat, rows::Vector{<:AbstractIS}, cols::Vector{<:AbstractIS}) end

@for_petsc function MatNestGetISs(petsclib::$UnionPetscLib, A::AbstractPetscMat, rows::Vector{<:AbstractIS}, cols::Vector{<:AbstractIS} )

    @chk ccall(
               (:MatNestGetISs, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CIS}, Ptr{CIS}),
               A, rows, cols,
              )


	return nothing
end 

"""
	MatNestGetLocalISs(petsclib::PetscLibType,A::AbstractPetscMat, rows::Vector{<:AbstractIS}, cols::Vector{<:AbstractIS}) 
Returns the index sets partitioning the row and column spaces of a `MATNEST`

Not Collective

Input Parameter:
- `A` - `MATNEST` matrix

Output Parameters:
- `rows` - array of row index sets (pass `NULL` to ignore)
- `cols` - array of column index sets (pass `NULL` to ignore)

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATNEST`, `MatNestGetSubMat()`, `MatNestGetSubMats()`, `MatNestGetSize()`, `MatNestGetISs()`, `MatCreateNest()`,
`MatNestSetSubMats()`, `MatNestSetSubMat()`

# External Links
$(_doc_external("Mat/MatNestGetLocalISs"))
"""
function MatNestGetLocalISs(petsclib::PetscLibType, A::AbstractPetscMat, rows::Vector{<:AbstractIS}, cols::Vector{<:AbstractIS}) end

@for_petsc function MatNestGetLocalISs(petsclib::$UnionPetscLib, A::AbstractPetscMat, rows::Vector{<:AbstractIS}, cols::Vector{<:AbstractIS} )

    @chk ccall(
               (:MatNestGetLocalISs, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CIS}, Ptr{CIS}),
               A, rows, cols,
              )


	return nothing
end 

"""
	M::PetscInt,N::PetscInt = MatNestGetSize(petsclib::PetscLibType,A::AbstractPetscMat) 
Returns the size of the `MATNEST` matrix.

Not Collective

Input Parameter:
- `A` - `MATNEST` matrix

Output Parameters:
- `M` - number of rows in the nested mat
- `N` - number of cols in the nested mat

Level: developer

-seealso: [](ch_matrices), `Mat`, `MATNEST`, `MatNestGetSubMat()`, `MatNestGetSubMats()`, `MatCreateNest()`, `MatNestGetLocalISs()`,
`MatNestGetISs()`

# External Links
$(_doc_external("Mat/MatNestGetSize"))
"""
function MatNestGetSize(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatNestGetSize(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	M_ = Ref{$PetscInt}()
	N_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatNestGetSize, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}, Ptr{$PetscInt}),
               A, M_, N_,
              )

	M = M_[]
	N = N_[]

	return M,N
end 

"""
	sub::PetscMat = MatNestGetSubMat(petsclib::PetscLibType,A::AbstractPetscMat, idxm::PetscInt, jdxm::PetscInt) 
Returns a single, sub

Not Collective

Input Parameters:
- `A`    - `MATNEST` matrix
- `idxm` - index of the matrix within the nest matrix
- `jdxm` - index of the matrix within the nest matrix

Output Parameter:
- `sub` - matrix at index `idxm`, `jdxm` within the nest matrix

Level: developer

-seealso: [](ch_matrices), `Mat`, `MATNEST`, `MatNestGetSize()`, `MatNestGetSubMats()`, `MatCreateNest()`, `MatNestSetSubMat()`,
`MatNestGetLocalISs()`, `MatNestGetISs()`

# External Links
$(_doc_external("Mat/MatNestGetSubMat"))
"""
function MatNestGetSubMat(petsclib::PetscLibType, A::AbstractPetscMat, idxm::PetscInt, jdxm::PetscInt) end

@for_petsc function MatNestGetSubMat(petsclib::$UnionPetscLib, A::AbstractPetscMat, idxm::$PetscInt, jdxm::$PetscInt )
	sub_ = Ref{CMat}()

    @chk ccall(
               (:MatNestGetSubMat, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, Ptr{CMat}),
               A, idxm, jdxm, sub_,
              )

	sub = PetscMat(sub_[], petsclib)

	return sub
end 

"""
	M::PetscInt,N::PetscInt,mat::PetscMat = MatNestGetSubMats(petsclib::PetscLibType,A::AbstractPetscMat) 
Returns the entire two dimensional array of matrices defining a `MATNEST` matrix.

Not Collective

Input Parameter:
- `A` - nest matrix

Output Parameters:
- `M`   - number of submatrix rows in the nest matrix
- `N`   - number of submatrix columns in the nest matrix
- `mat` - array of matrices

Level: developer

-seealso: [](ch_matrices), `Mat`, `MATNEST`, `MatNestGetSize()`, `MatNestGetSubMat()`, `MatNestGetLocalISs()`, `MatCreateNest()`,
`MatNestSetSubMats()`, `MatNestGetISs()`, `MatNestSetSubMat()`

# External Links
$(_doc_external("Mat/MatNestGetSubMats"))
"""
function MatNestGetSubMats(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatNestGetSubMats(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	M_ = Ref{$PetscInt}()
	N_ = Ref{$PetscInt}()
	mat_ = Ref{CMat}()

    @chk ccall(
               (:MatNestGetSubMats, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{CMat}}),
               A, M_, N_, mat_,
              )

	M = M_[]
	N = N_[]
	mat = PetscMat(mat_[], petsclib)

	return M,N,mat
end 

"""
	MatNestSetSubMat(petsclib::PetscLibType,A::AbstractPetscMat, idxm::PetscInt, jdxm::PetscInt, sub::AbstractPetscMat) 
Set a single submatrix in the `MATNEST`

Logically Collective

Input Parameters:
- `A`    - `MATNEST` matrix
- `idxm` - index of the matrix within the nest matrix
- `jdxm` - index of the matrix within the nest matrix
- `sub`  - matrix at index `idxm`, `jdxm` within the nest matrix

Level: developer

-seealso: [](ch_matrices), `Mat`, `MATNEST`, `MatNestSetSubMats()`, `MatNestGetSubMats()`, `MatNestGetLocalISs()`, `MatCreateNest()`,
`MatNestGetSubMat()`, `MatNestGetISs()`, `MatNestGetSize()`

# External Links
$(_doc_external("Mat/MatNestSetSubMat"))
"""
function MatNestSetSubMat(petsclib::PetscLibType, A::AbstractPetscMat, idxm::PetscInt, jdxm::PetscInt, sub::AbstractPetscMat) end

@for_petsc function MatNestSetSubMat(petsclib::$UnionPetscLib, A::AbstractPetscMat, idxm::$PetscInt, jdxm::$PetscInt, sub::AbstractPetscMat )

    @chk ccall(
               (:MatNestSetSubMat, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, CMat),
               A, idxm, jdxm, sub,
              )


	return nothing
end 

"""
	MatNestSetSubMats(petsclib::PetscLibType,A::AbstractPetscMat, nr::PetscInt, is_row::Vector{<:AbstractIS}, nc::PetscInt, is_col::Vector{<:AbstractIS}, a::Vector{<:AbstractPetscMat}) 
Sets the nested submatrices in a `MATNEST`

Collective

Input Parameters:
- `A`      - `MATNEST` matrix
- `nr`     - number of nested row blocks
- `is_row` - index sets for each nested row block, or `NULL` to make contiguous
- `nc`     - number of nested column blocks
- `is_col` - index sets for each nested column block, or `NULL` to make contiguous
- `a`      - array of  nr \\times nc submatrices, or `NULL`

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATNEST`, `MatCreateNest()`, `MatNestSetSubMat()`, `MatNestGetSubMat()`, `MatNestGetSubMats()`

# External Links
$(_doc_external("Mat/MatNestSetSubMats"))
"""
function MatNestSetSubMats(petsclib::PetscLibType, A::AbstractPetscMat, nr::PetscInt, is_row::Vector{<:AbstractIS}, nc::PetscInt, is_col::Vector{<:AbstractIS}, a::Vector{<:AbstractPetscMat}) end

@for_petsc function MatNestSetSubMats(petsclib::$UnionPetscLib, A::AbstractPetscMat, nr::$PetscInt, is_row::Vector{<:AbstractIS}, nc::$PetscInt, is_col::Vector{<:AbstractIS}, a::Vector{<:AbstractPetscMat} )

    @chk ccall(
               (:MatNestSetSubMats, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CIS}, $PetscInt, Ptr{CIS}, Ptr{CMat}),
               A, nr, is_row, nc, is_col, a,
              )


	return nothing
end 

"""
	MatNestSetVecType(petsclib::PetscLibType,A::AbstractPetscMat, vtype::VecType) 
Sets the type of `Vec` returned by `MatCreateVecs()`

Not Collective

Input Parameters:
- `A`     - `MATNEST` matrix
- `vtype` - `VecType` to use for creating vectors

Level: developer

-seealso: [](ch_matrices), `Mat`, `MATNEST`, `MatCreateVecs()`, `MatCreateNest()`, `VecType`

# External Links
$(_doc_external("Mat/MatNestSetVecType"))
"""
function MatNestSetVecType(petsclib::PetscLibType, A::AbstractPetscMat, vtype::VecType) end

@for_petsc function MatNestSetVecType(petsclib::$UnionPetscLib, A::AbstractPetscMat, vtype::VecType )

    @chk ccall(
               (:MatNestSetVecType, $petsc_library),
               PetscErrorCode,
               (CMat, VecType),
               A, vtype,
              )


	return nothing
end 

"""
	nrm::PetscReal = MatNorm(petsclib::PetscLibType,mat::AbstractPetscMat, type::NormType) 
Calculates various norms of a matrix.

Collective

Input Parameters:
- `mat`  - the matrix
- `type` - the type of norm, `NORM_1`, `NORM_FROBENIUS`, `NORM_INFINITY`

Output Parameter:
- `nrm` - the resulting norm

Level: intermediate

-seealso: [](ch_matrices), `Mat`

# External Links
$(_doc_external("Mat/MatNorm"))
"""
function MatNorm(petsclib::PetscLibType, mat::AbstractPetscMat, type::NormType) end

@for_petsc function MatNorm(petsclib::$UnionPetscLib, mat::AbstractPetscMat, type::NormType )
	nrm_ = Ref{$PetscReal}()

    @chk ccall(
               (:MatNorm, $petsc_library),
               PetscErrorCode,
               (CMat, NormType, Ptr{$PetscReal}),
               mat, type, nrm_,
              )

	nrm = nrm_[]

	return nrm
end 

"""
	M::PetscMat = MatNormalGetMat(petsclib::PetscLibType,A::AbstractPetscMat) 
Gets the `Mat` object stored inside a `MATNORMAL`

Logically Collective

Input Parameter:
- `A` - the `MATNORMAL` matrix

Output Parameter:
- `M` - the matrix object stored inside `A`

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATNORMAL`, `MATNORMALHERMITIAN`, `MatCreateNormal()`

# External Links
$(_doc_external("Mat/MatNormalGetMat"))
"""
function MatNormalGetMat(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatNormalGetMat(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	M_ = Ref{CMat}()

    @chk ccall(
               (:MatNormalGetMat, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, M_,
              )

	M = PetscMat(M_[], petsclib)

	return M
end 

"""
	M::PetscMat = MatNormalHermitianGetMat(petsclib::PetscLibType,A::AbstractPetscMat) 
Gets the `Mat` object stored inside a `MATNORMALHERMITIAN`

Logically Collective

Input Parameter:
- `A` - the `MATNORMALHERMITIAN` matrix

Output Parameter:
- `M` - the matrix object stored inside `A`

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATNORMALHERMITIAN`, `MatCreateNormalHermitian()`

# External Links
$(_doc_external("Mat/MatNormalHermitianGetMat"))
"""
function MatNormalHermitianGetMat(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatNormalHermitianGetMat(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	M_ = Ref{CMat}()

    @chk ccall(
               (:MatNormalHermitianGetMat, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, M_,
              )

	M = PetscMat(M_[], petsclib)

	return M
end 

"""
	MatOrderingRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Adds a new sparse matrix ordering to the matrix package.

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - name of ordering (for example `MATORDERINGND`)
- `function` - function pointer that creates the ordering

Level: developer

-seealso: `Mat`, `MatOrderingType`, `MatOrderingRegisterAll()`, `MatGetOrdering()`

# External Links
$(_doc_external("MatGraphOperations/MatOrderingRegister"))
"""
function MatOrderingRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function MatOrderingRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:MatOrderingRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	B::PetscMat = MatPermute(petsclib::PetscLibType,mat::AbstractPetscMat, row::AbstractIS, col::AbstractIS) 
Creates a new matrix with rows and columns permuted from the
original.

Collective

Input Parameters:
- `mat` - the matrix to permute
- `row` - row permutation, each processor supplies only the permutation for its rows
- `col` - column permutation, each processor supplies only the permutation for its columns

Output Parameter:
- `B` - the permuted matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatGetOrdering()`, `ISAllGather()`, `MatCreateSubMatrix()`

# External Links
$(_doc_external("Mat/MatPermute"))
"""
function MatPermute(petsclib::PetscLibType, mat::AbstractPetscMat, row::AbstractIS, col::AbstractIS) end

@for_petsc function MatPermute(petsclib::$UnionPetscLib, mat::AbstractPetscMat, row::AbstractIS, col::AbstractIS )
	B_ = Ref{CMat}()

    @chk ccall(
               (:MatPermute, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, CIS, Ptr{CMat}),
               mat, row, col, B_,
              )

	B = PetscMat(B_[], petsclib)

	return B
end 

"""
	MatPreallocatorPreallocate(petsclib::PetscLibType,mat::AbstractPetscMat, fill::PetscBool, A::AbstractPetscMat) 
Preallocates the A matrix, using information from a `MATPREALLOCATOR` mat, optionally filling A with zeros

Input Parameters:
- `mat`  - the `MATPREALLOCATOR` preallocator matrix
- `fill` - fill the matrix with zeros
- `A`    - the matrix to be preallocated

-seealso: `MATPREALLOCATOR`, `MatXAIJSetPreallocation()`

# External Links
$(_doc_external("Mat/MatPreallocatorPreallocate"))
"""
function MatPreallocatorPreallocate(petsclib::PetscLibType, mat::AbstractPetscMat, fill::PetscBool, A::AbstractPetscMat) end

@for_petsc function MatPreallocatorPreallocate(petsclib::$UnionPetscLib, mat::AbstractPetscMat, fill::PetscBool, A::AbstractPetscMat )

    @chk ccall(
               (:MatPreallocatorPreallocate, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool, CMat),
               mat, fill, A,
              )


	return nothing
end 

"""
	MatProductClear(petsclib::PetscLibType,mat::AbstractPetscMat) 
Clears from the matrix any internal data structures related to the computation of the values of the matrix from matrix

Collective

Input Parameter:
- `mat` - the matrix whose values are to be computed via a matrix-matrix product operation

Options Database Key:
- `-mat_product_clear` - Clear intermediate data structures after `MatProductNumeric()` has been called

Level: intermediate

-seealso: [](ch_matrices), `MatProduct`, `Mat`, `MatProductCreate()`

# External Links
$(_doc_external("Mat/MatProductClear"))
"""
function MatProductClear(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatProductClear(petsclib::$UnionPetscLib, mat::AbstractPetscMat )

    @chk ccall(
               (:MatProductClear, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	D::PetscMat = MatProductCreate(petsclib::PetscLibType,A::AbstractPetscMat, B::AbstractPetscMat, C::Union{Ptr, AbstractPetscMat}) 
create a matrix to hold the result of a matrix

Collective

Input Parameters:
- `A` - the first matrix
- `B` - the second matrix
- `C` - the third matrix (or `NULL`)

Output Parameter:
- `D` - the matrix whose values are to be computed via a matrix-matrix product operation

Level: intermediate

Example:
-seealso: [](ch_matrices), `MatProduct`, `Mat`, `MatProductCreateWithMat()`, `MatProductSetType()`, `MatProductSetAlgorithm()`, `MatProductClear()`

# External Links
$(_doc_external("Mat/MatProductCreate"))
"""
function MatProductCreate(petsclib::PetscLibType, A::AbstractPetscMat, B::AbstractPetscMat, C::Union{Ptr, AbstractPetscMat}) end

@for_petsc function MatProductCreate(petsclib::$UnionPetscLib, A::AbstractPetscMat, B::AbstractPetscMat, C::Union{Ptr, AbstractPetscMat} )
	D_ = Ref{CMat}()

    @chk ccall(
               (:MatProductCreate, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat, Ptr{CMat}),
               A, B, C, D_,
              )

	D = PetscMat(D_[], petsclib)

	return D
end 

"""
	MatProductCreateWithMat(petsclib::PetscLibType,A::AbstractPetscMat, B::AbstractPetscMat, C::Union{Ptr, AbstractPetscMat}, D::AbstractPetscMat) 
Set a given matrix to have its values computed via matrix

Collective

Input Parameters:
- `A` - the first matrix
- `B` - the second matrix
- `C` - the third matrix (optional, use `NULL` if not needed)
- `D` - the matrix whose values are to be computed via a matrix-matrix product operation

Level: intermediate

-seealso: [](ch_matrices), `MatProduct`, `Mat`, `MatProductType`, `MatProductSetType()`, `MatProductAlgorithm`,
`MatProductSetAlgorithm`, `MatProductCreate()`, `MatProductClear()`

# External Links
$(_doc_external("Mat/MatProductCreateWithMat"))
"""
function MatProductCreateWithMat(petsclib::PetscLibType, A::AbstractPetscMat, B::AbstractPetscMat, C::Union{Ptr, AbstractPetscMat}, D::AbstractPetscMat) end

@for_petsc function MatProductCreateWithMat(petsclib::$UnionPetscLib, A::AbstractPetscMat, B::AbstractPetscMat, C::Union{Ptr, AbstractPetscMat}, D::AbstractPetscMat )

    @chk ccall(
               (:MatProductCreateWithMat, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat, CMat),
               A, B, C, D,
              )


	return nothing
end 

"""
	alg::MatProductAlgorithm = MatProductGetAlgorithm(petsclib::PetscLibType,mat::AbstractPetscMat) 
Returns the selected algorithm for a matrix

Not Collective

Input Parameter:
- `mat` - the matrix whose values are computed via a matrix-matrix product operation

Output Parameter:
- `alg` - the selected algorithm of the matrix product, e.g., `MATPRODUCTALGORITHMDEFAULT`.

Level: intermediate

-seealso: [](ch_matrices), `MatProduct`, `Mat`, `MatProductSetAlgorithm()`

# External Links
$(_doc_external("Mat/MatProductGetAlgorithm"))
"""
function MatProductGetAlgorithm(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatProductGetAlgorithm(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	alg_ = Ref{MatProductAlgorithm}()

    @chk ccall(
               (:MatProductGetAlgorithm, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatProductAlgorithm}),
               mat, alg_,
              )

	alg = alg_[] == C_NULL ? "" : unsafe_string(alg_[])

	return alg
end 

"""
	A::PetscMat,B::PetscMat,C::PetscMat = MatProductGetMats(petsclib::PetscLibType,mat::AbstractPetscMat) 
Returns the matrices associated with the matrix

Not Collective

Input Parameter:
- `mat` - the matrix whose values are to be computed via a matrix-matrix product operation

Output Parameters:
- `A` - the first matrix
- `B` - the second matrix
- `C` - the third matrix (may be `NULL` for some `MatProductType`)

Level: intermediate

-seealso: [](ch_matrices), `MatProduct`, `Mat`, `MatProductCreateWithMat()`, `MatProductSetType()`, `MatProductSetAlgorithm()`, `MatProductCreate()`

# External Links
$(_doc_external("Mat/MatProductGetMats"))
"""
function MatProductGetMats(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatProductGetMats(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	A_ = Ref{CMat}()
	B_ = Ref{CMat}()
	C_ = Ref{CMat}()

    @chk ccall(
               (:MatProductGetMats, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}, Ptr{CMat}, Ptr{CMat}),
               mat, A_, B_, C_,
              )

	A = PetscMat(A_[], petsclib)
	B = PetscMat(B_[], petsclib)
	C = PetscMat(C_[], petsclib)

	return A,B,C
end 

"""
	mtype::MatProductType = MatProductGetType(petsclib::PetscLibType,mat::AbstractPetscMat) 
Returns the type of matrix

Not Collective

Input Parameter:
- `mat` - the matrix whose values are to be computed via a matrix-matrix product operation

Output Parameter:
- `mtype` - the `MatProductType`

Level: intermediate

-seealso: [](ch_matrices), `MatProduct`, `Mat`, `MatProductCreateWithMat()`, `MatProductSetType()`, `MatProductCreate()`, `MatProductType`, `MatProductAlgorithm`

# External Links
$(_doc_external("Mat/MatProductGetType"))
"""
function MatProductGetType(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatProductGetType(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	mtype_ = Ref{MatProductType}()

    @chk ccall(
               (:MatProductGetType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatProductType}),
               mat, mtype_,
              )

	mtype = mtype_[]

	return mtype
end 

"""
	MatProductNumeric(petsclib::PetscLibType,mat::AbstractPetscMat) 
Compute a matrix

Collective

Input/Output Parameter:
- `mat` - the matrix whose values are computed via a matrix-matrix product operation

Level: intermediate

-seealso: [](ch_matrices), `MatProduct`, `Mat`, `MatProductSetAlgorithm()`, `MatProductSetType()`, `MatProductCreate()`, `MatSetType()`, `MatProductSymbolic()`

# External Links
$(_doc_external("Mat/MatProductNumeric"))
"""
function MatProductNumeric(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatProductNumeric(petsclib::$UnionPetscLib, mat::AbstractPetscMat )

    @chk ccall(
               (:MatProductNumeric, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	MatProductReplaceMats(petsclib::PetscLibType,A::Union{Ptr, AbstractPetscMat}, B::Union{Ptr, AbstractPetscMat}, C::Union{Ptr, AbstractPetscMat}, D::Union{Ptr, AbstractPetscMat}) 
Replace the input matrices for the matrix

Collective

Input Parameters:
- `A` - the matrix or `NULL` if not being replaced
- `B` - the matrix or `NULL` if not being replaced
- `C` - the matrix or `NULL` if not being replaced
- `D` - the matrix whose values are computed via a matrix-matrix product operation

Level: intermediate

-seealso: [](ch_matrices), `MatProduct`, `Mat`, `MatProductCreate()`, `MatProductSetFromOptions()`, `MatProductSymbolic()`, `MatProductClear()`

# External Links
$(_doc_external("Mat/MatProductReplaceMats"))
"""
function MatProductReplaceMats(petsclib::PetscLibType, A::Union{Ptr, AbstractPetscMat}, B::Union{Ptr, AbstractPetscMat}, C::Union{Ptr, AbstractPetscMat}, D::Union{Ptr, AbstractPetscMat}) end

@for_petsc function MatProductReplaceMats(petsclib::$UnionPetscLib, A::Union{Ptr, AbstractPetscMat}, B::Union{Ptr, AbstractPetscMat}, C::Union{Ptr, AbstractPetscMat}, D::Union{Ptr, AbstractPetscMat} )

    @chk ccall(
               (:MatProductReplaceMats, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat, CMat),
               A, B, C, D,
              )


	return nothing
end 

"""
	MatProductSetAlgorithm(petsclib::PetscLibType,mat::AbstractPetscMat, alg::MatProductAlgorithm) 
Requests a particular algorithm for a matrix

Collective

Input Parameters:
- `mat` - the matrix whose values are computed via a matrix-matrix product operation
- `alg` - particular implementation algorithm of the matrix product, e.g., `MATPRODUCTALGORITHMDEFAULT`.

Options Database Key:
- `-mat_product_algorithm <algorithm>` - Sets the algorithm, see `MatProductAlgorithm`

Level: intermediate

-seealso: [](ch_matrices), `MatProduct`, `Mat`, `MatProductClear()`, `MatProductSetType()`, `MatProductSetFill()`, `MatProductCreate()`, `MatProductAlgorithm`, `MatProductType`, `MatProductGetAlgorithm()`

# External Links
$(_doc_external("Mat/MatProductSetAlgorithm"))
"""
function MatProductSetAlgorithm(petsclib::PetscLibType, mat::AbstractPetscMat, alg::MatProductAlgorithm) end

@for_petsc function MatProductSetAlgorithm(petsclib::$UnionPetscLib, mat::AbstractPetscMat, alg::MatProductAlgorithm )

    @chk ccall(
               (:MatProductSetAlgorithm, $petsc_library),
               PetscErrorCode,
               (CMat, MatProductAlgorithm),
               mat, alg,
              )


	return nothing
end 

"""
	MatProductSetFill(petsclib::PetscLibType,mat::AbstractPetscMat, fill::PetscReal) 
Set an expected fill of the matrix whose values are computed via a matrix

Collective

Input Parameters:
- `mat`  - the matrix whose values are to be computed via a matrix-matrix product operation
- `fill` - expected fill as ratio of nnz(mat)/(nnz(A) + nnz(B) + nnz(C)); use `PETSC_DETERMINE` or `PETSC_CURRENT` if you do not have a good estimate.
If the product is a dense matrix, this value is not used.

Level: intermediate

-seealso: [](ch_matrices), `MatProduct`, `PETSC_DETERMINE`, `Mat`, `MatProductSetFromOptions()`, `MatProductSetType()`, `MatProductSetAlgorithm()`, `MatProductCreate()`

# External Links
$(_doc_external("Mat/MatProductSetFill"))
"""
function MatProductSetFill(petsclib::PetscLibType, mat::AbstractPetscMat, fill::PetscReal) end

@for_petsc function MatProductSetFill(petsclib::$UnionPetscLib, mat::AbstractPetscMat, fill::$PetscReal )

    @chk ccall(
               (:MatProductSetFill, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscReal),
               mat, fill,
              )


	return nothing
end 

"""
	MatProductSetFromOptions(petsclib::PetscLibType,mat::AbstractPetscMat) 
Sets the options for the computation of a matrix
the algorithm etc are determined from the options database.

Logically Collective

Input Parameter:
- `mat` - the matrix whose values are computed via a matrix-matrix product operation

Options Database Keys:
- `-mat_product_clear`                 - Clear intermediate data structures after `MatProductNumeric()` has been called
- `-mat_product_algorithm <algorithm>` - Sets the algorithm, see `MatProductAlgorithm` for possible values
- `-mat_product_algorithm_backend_cpu` - Use the CPU to perform the computation even if the matrix is a GPU matrix

Level: intermediate

-seealso: [](ch_matrices), `MatProduct`, `Mat`, `MatSetFromOptions()`, `MatProductCreate()`, `MatProductCreateWithMat()`, `MatProductNumeric()`,
`MatProductSetType()`, `MatProductSetAlgorithm()`, `MatProductAlgorithm`

# External Links
$(_doc_external("Mat/MatProductSetFromOptions"))
"""
function MatProductSetFromOptions(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatProductSetFromOptions(petsclib::$UnionPetscLib, mat::AbstractPetscMat )

    @chk ccall(
               (:MatProductSetFromOptions, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	MatProductSetType(petsclib::PetscLibType,mat::AbstractPetscMat, productype::MatProductType) 
Sets a particular matrix

Collective

Input Parameters:
- `mat`        - the matrix whose values are computed via a matrix-matrix product operation
- `productype` - matrix product type, e.g., `MATPRODUCT_AB`,`MATPRODUCT_AtB`,`MATPRODUCT_ABt`,`MATPRODUCT_PtAP`,`MATPRODUCT_RARt`,`MATPRODUCT_ABC`,
see `MatProductType`

Level: intermediate

-seealso: [](ch_matrices), `MatProduct`, `Mat`, `MatProductCreate()`, `MatProductType`,
`MATPRODUCT_AB`, `MATPRODUCT_AtB`, `MATPRODUCT_ABt`, `MATPRODUCT_PtAP`, `MATPRODUCT_RARt`, `MATPRODUCT_ABC`

# External Links
$(_doc_external("Mat/MatProductSetType"))
"""
function MatProductSetType(petsclib::PetscLibType, mat::AbstractPetscMat, productype::MatProductType) end

@for_petsc function MatProductSetType(petsclib::$UnionPetscLib, mat::AbstractPetscMat, productype::MatProductType )

    @chk ccall(
               (:MatProductSetType, $petsc_library),
               PetscErrorCode,
               (CMat, MatProductType),
               mat, productype,
              )


	return nothing
end 

"""
	MatProductSymbolic(petsclib::PetscLibType,mat::AbstractPetscMat) 
Perform the symbolic portion of a matrix
product to be done with `MatProductNumeric()`

Collective

Input/Output Parameter:
- `mat` - the matrix whose values are to be computed via a matrix-matrix product operation

Level: intermediate

-seealso: [](ch_matrices), `MatProduct`, `Mat`, `MatProductCreate()`, `MatProductCreateWithMat()`, `MatProductSetFromOptions()`, `MatProductNumeric()`, `MatProductSetType()`, `MatProductSetAlgorithm()`

# External Links
$(_doc_external("Mat/MatProductSymbolic"))
"""
function MatProductSymbolic(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatProductSymbolic(petsclib::$UnionPetscLib, mat::AbstractPetscMat )

    @chk ccall(
               (:MatProductSymbolic, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	MatProductView(petsclib::PetscLibType,mat::AbstractPetscMat, viewer::PetscViewer) 
View the private matrix

Logically Collective

Input Parameters:
- `mat`    - the matrix obtained with `MatProductCreate()` or `MatProductCreateWithMat()`
- `viewer` - where the information on the matrix-matrix algorithm of `mat` should be reviewed

Level: intermediate

-seealso: [](ch_matrices), `MatProductType`, `Mat`, `MatProductSetFromOptions()`, `MatView()`, `MatProductCreate()`, `MatProductCreateWithMat()`

# External Links
$(_doc_external("Mat/MatProductView"))
"""
function MatProductView(petsclib::PetscLibType, mat::AbstractPetscMat, viewer::PetscViewer) end

@for_petsc function MatProductView(petsclib::$UnionPetscLib, mat::AbstractPetscMat, viewer::PetscViewer )

    @chk ccall(
               (:MatProductView, $petsc_library),
               PetscErrorCode,
               (CMat, PetscViewer),
               mat, viewer,
              )


	return nothing
end 

"""
	MatPropagateSymmetryOptions(petsclib::PetscLibType,A::AbstractPetscMat, B::AbstractPetscMat) 
Propagates symmetry options set on a matrix to another matrix

Not Collective

Input Parameters:
- `A` - the matrix we wish to propagate options from
- `B` - the matrix we wish to propagate options to

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatSetOption()`, `MatIsSymmetricKnown()`, `MatIsSPDKnown()`, `MatIsHermitianKnown()`, `MatIsStructurallySymmetricKnown()`

# External Links
$(_doc_external("Mat/MatPropagateSymmetryOptions"))
"""
function MatPropagateSymmetryOptions(petsclib::PetscLibType, A::AbstractPetscMat, B::AbstractPetscMat) end

@for_petsc function MatPropagateSymmetryOptions(petsclib::$UnionPetscLib, A::AbstractPetscMat, B::AbstractPetscMat )

    @chk ccall(
               (:MatPropagateSymmetryOptions, $petsc_library),
               PetscErrorCode,
               (CMat, CMat),
               A, B,
              )


	return nothing
end 

"""
	C::PetscMat = MatPtAP(petsclib::PetscLibType,A::AbstractPetscMat, P::AbstractPetscMat, scall::MatReuse, fill::PetscReal) 
Creates the matrix product C = P^T * A * P

Neighbor-wise Collective

Input Parameters:
- `A`     - the matrix
- `P`     - the projection matrix
- `scall` - either `MAT_INITIAL_MATRIX` or `MAT_REUSE_MATRIX`
- `fill`  - expected fill as ratio of nnz(C)/(nnz(A) + nnz(P)), use `PETSC_DETERMINE` or `PETSC_CURRENT` if you do not have a good estimate
if the result is a dense matrix this is irrelevant

Output Parameter:
- `C` - the product matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatProductCreate()`, `MatMatMult()`, `MatRARt()`

# External Links
$(_doc_external("Mat/MatPtAP"))
"""
function MatPtAP(petsclib::PetscLibType, A::AbstractPetscMat, P::AbstractPetscMat, scall::MatReuse, fill::PetscReal) end

@for_petsc function MatPtAP(petsclib::$UnionPetscLib, A::AbstractPetscMat, P::AbstractPetscMat, scall::MatReuse, fill::$PetscReal )
	C_ = Ref{CMat}()

    @chk ccall(
               (:MatPtAP, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, MatReuse, $PetscReal, Ptr{CMat}),
               A, P, scall, fill, C_,
              )

	C = PetscMat(C_[], petsclib)

	return C
end 

"""
	flg::PetscBool = MatPtAPMultEqual(petsclib::PetscLibType,A::AbstractPetscMat, B::AbstractPetscMat, C::AbstractPetscMat, n::PetscInt) 
Compares matrix

Collective

Input Parameters:
- `A` - the first matrix
- `B` - the second matrix
- `C` - the third matrix
- `n` - number of random vectors to be tested

Output Parameter:
- `flg` - `PETSC_TRUE` if the products are equal; `PETSC_FALSE` otherwise.

Level: intermediate

-seealso: `Mat`, `MatMatMultEqual()`, `MatMultEqual()`, `MatMultAddEqual()`, `MatMultTransposeEqual()`

# External Links
$(_doc_external("Mat/MatPtAPMultEqual"))
"""
function MatPtAPMultEqual(petsclib::PetscLibType, A::AbstractPetscMat, B::AbstractPetscMat, C::AbstractPetscMat, n::PetscInt) end

@for_petsc function MatPtAPMultEqual(petsclib::$UnionPetscLib, A::AbstractPetscMat, B::AbstractPetscMat, C::AbstractPetscMat, n::$PetscInt )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatPtAPMultEqual, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat, $PetscInt, Ptr{PetscBool}),
               A, B, C, n, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	A::PetscMat = MatPythonCreate(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, pyname::String) 
Create a `Mat` object implemented in Python.

Collective

Input Parameters:
- `comm`  - MPI communicator
- `m`  - number of local rows (or `PETSC_DECIDE` to have calculated if `M` is given)
- `n`  - number of local columns (or `PETSC_DECIDE` to have calculated if `N` is given)
- `M`  - number of global rows (or `PETSC_DECIDE` to have calculated if `m` is given)
- `N`  - number of global columns (or `PETSC_DECIDE` to have calculated if `n` is given)
- `pyname`  - full dotted Python name [package].module[.{class|function}]

Output Parameter:
- `A`  - the matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatType`, `MATPYTHON`, `MatPythonSetType()`, `PetscPythonInitialize()`

# External Links
$(_doc_external("Mat/MatPythonCreate"))
"""
function MatPythonCreate(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, pyname::String) end

@for_petsc function MatPythonCreate(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, pyname::String )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatPythonCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{Cchar}, Ptr{CMat}),
               comm, m, n, M, N, pyname, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	pyname::Ptr{Cchar} = MatPythonGetType(petsclib::PetscLibType,mat::AbstractPetscMat) 
Get the Python name of a `Mat` object implemented in Python.

Not Collective

Input Parameter:
- `mat`  - the matrix

Output Parameter:
- `pyname`  - full dotted Python name [package].module[.{class|function}]

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatType`, `MatCreate()`, `MatSetType()`, `MATPYTHON`, `PetscPythonInitialize()`, `MatPythonSetType()`

# External Links
$(_doc_external("Mat/MatPythonGetType"))
"""
function MatPythonGetType(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatPythonGetType(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	pyname_ = Ref{Ptr{Cchar}}()

    @chk ccall(
               (:MatPythonGetType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{Cchar}}),
               mat, pyname_,
              )

	pyname = pyname_[]

	return pyname
end 

"""
	MatPythonSetType(petsclib::PetscLibType,mat::AbstractPetscMat, pyname::String) 
Initialize a `Mat` object implemented in Python.

Collective

Input Parameters:
- `mat`  - the matrix object.
- `pyname`  - full dotted Python name [package].module[.{class|function}]

Options Database Key:
- `-mat_python_type <pyname>`  - python class

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatType`, `MatCreate()`, `MatSetType()`, `MATPYTHON`, `PetscPythonInitialize()`

# External Links
$(_doc_external("Mat/MatPythonSetType"))
"""
function MatPythonSetType(petsclib::PetscLibType, mat::AbstractPetscMat, pyname::String) end

@for_petsc function MatPythonSetType(petsclib::$UnionPetscLib, mat::AbstractPetscMat, pyname::String )

    @chk ccall(
               (:MatPythonSetType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Cchar}),
               mat, pyname,
              )


	return nothing
end 

"""
	MatQRFactor(petsclib::PetscLibType,mat::AbstractPetscMat, col::AbstractIS, info::Vector{MatFactorInfo}) 
Performs in

Collective

Input Parameters:
- `mat`  - the matrix
- `col`  - column permutation
- `info` - options for factorization, includes
-seealso: [](ch_matrices), `Mat`, [Matrix Factorization](sec_matfactor), `MatFactorInfo`, `MatGetFactor()`, `MatQRFactorSymbolic()`, `MatQRFactorNumeric()`, `MatLUFactor()`,
`MatSetUnfactored()`

# External Links
$(_doc_external("Mat/MatQRFactor"))
"""
function MatQRFactor(petsclib::PetscLibType, mat::AbstractPetscMat, col::AbstractIS, info::Vector{MatFactorInfo}) end

@for_petsc function MatQRFactor(petsclib::$UnionPetscLib, mat::AbstractPetscMat, col::AbstractIS, info::Vector{MatFactorInfo} )

    @chk ccall(
               (:MatQRFactor, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, Ptr{MatFactorInfo}),
               mat, col, info,
              )


	return nothing
end 

"""
	MatQRFactorNumeric(petsclib::PetscLibType,fact::AbstractPetscMat, mat::AbstractPetscMat, info::Vector{MatFactorInfo}) 
Performs numeric QR factorization of a matrix.
Call this routine after first calling `MatGetFactor()`, and `MatQRFactorSymbolic()`.

Collective

Input Parameters:
- `fact` - the factor matrix obtained with `MatGetFactor()`
- `mat`  - the matrix
- `info` - options for factorization

Level: developer

-seealso: [](ch_matrices), `Mat`, [Matrix Factorization](sec_matfactor), `MatFactorInfo`, `MatGetFactor()`, `MatQRFactor()`, `MatQRFactorSymbolic()`, `MatLUFactor()`

# External Links
$(_doc_external("Mat/MatQRFactorNumeric"))
"""
function MatQRFactorNumeric(petsclib::PetscLibType, fact::AbstractPetscMat, mat::AbstractPetscMat, info::Vector{MatFactorInfo}) end

@for_petsc function MatQRFactorNumeric(petsclib::$UnionPetscLib, fact::AbstractPetscMat, mat::AbstractPetscMat, info::Vector{MatFactorInfo} )

    @chk ccall(
               (:MatQRFactorNumeric, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, Ptr{MatFactorInfo}),
               fact, mat, info,
              )


	return nothing
end 

"""
	MatQRFactorSymbolic(petsclib::PetscLibType,fact::AbstractPetscMat, mat::AbstractPetscMat, col::AbstractIS, info::Vector{MatFactorInfo}) 
Performs symbolic QR factorization of matrix.
Call this routine after `MatGetFactor()` but before calling `MatQRFactorNumeric()`.

Collective

Input Parameters:
- `fact` - the factor matrix obtained with `MatGetFactor()`
- `mat`  - the matrix
- `col`  - column permutation
- `info` - options for factorization, includes
-seealso: [](ch_matrices), `Mat`, [Matrix Factorization](sec_matfactor), `MatGetFactor()`, `MatFactorInfo`, `MatQRFactor()`, `MatQRFactorNumeric()`, `MatLUFactor()`, `MatFactorInfoInitialize()`

# External Links
$(_doc_external("Mat/MatQRFactorSymbolic"))
"""
function MatQRFactorSymbolic(petsclib::PetscLibType, fact::AbstractPetscMat, mat::AbstractPetscMat, col::AbstractIS, info::Vector{MatFactorInfo}) end

@for_petsc function MatQRFactorSymbolic(petsclib::$UnionPetscLib, fact::AbstractPetscMat, mat::AbstractPetscMat, col::AbstractIS, info::Vector{MatFactorInfo} )

    @chk ccall(
               (:MatQRFactorSymbolic, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CIS, Ptr{MatFactorInfo}),
               fact, mat, col, info,
              )


	return nothing
end 

"""
	C::PetscMat = MatRARt(petsclib::PetscLibType,A::AbstractPetscMat, R::AbstractPetscMat, scall::MatReuse, fill::PetscReal) 
Creates the matrix product C = R * A * R^T

Neighbor-wise Collective

Input Parameters:
- `A`     - the matrix
- `R`     - the projection matrix
- `scall` - either `MAT_INITIAL_MATRIX` or `MAT_REUSE_MATRIX`
- `fill`  - expected fill as ratio of nnz(C)/nnz(A), use `PETSC_DETERMINE` or `PETSC_CURRENT` if you do not have a good estimate
if the result is a dense matrix this is irrelevant

Output Parameter:
- `C` - the product matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatProductCreate()`, `MatMatMult()`, `MatPtAP()`

# External Links
$(_doc_external("Mat/MatRARt"))
"""
function MatRARt(petsclib::PetscLibType, A::AbstractPetscMat, R::AbstractPetscMat, scall::MatReuse, fill::PetscReal) end

@for_petsc function MatRARt(petsclib::$UnionPetscLib, A::AbstractPetscMat, R::AbstractPetscMat, scall::MatReuse, fill::$PetscReal )
	C_ = Ref{CMat}()

    @chk ccall(
               (:MatRARt, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, MatReuse, $PetscReal, Ptr{CMat}),
               A, R, scall, fill, C_,
              )

	C = PetscMat(C_[], petsclib)

	return C
end 

"""
	flg::PetscBool = MatRARtMultEqual(petsclib::PetscLibType,A::AbstractPetscMat, B::AbstractPetscMat, C::AbstractPetscMat, n::PetscInt) 
Compares matrix

Collective

Input Parameters:
- `A` - the first matrix
- `B` - the second matrix
- `C` - the third matrix
- `n` - number of random vectors to be tested

Output Parameter:
- `flg` - `PETSC_TRUE` if the products are equal; `PETSC_FALSE` otherwise.

Level: intermediate

-seealso: `Mat`, `MatMatMultEqual()`, `MatMultEqual()`, `MatMultAddEqual()`, `MatMultTransposeEqual()`

# External Links
$(_doc_external("Mat/MatRARtMultEqual"))
"""
function MatRARtMultEqual(petsclib::PetscLibType, A::AbstractPetscMat, B::AbstractPetscMat, C::AbstractPetscMat, n::PetscInt) end

@for_petsc function MatRARtMultEqual(petsclib::$UnionPetscLib, A::AbstractPetscMat, B::AbstractPetscMat, C::AbstractPetscMat, n::$PetscInt )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatRARtMultEqual, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat, $PetscInt, Ptr{PetscBool}),
               A, B, C, n, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	MatRealPart(petsclib::PetscLibType,mat::AbstractPetscMat) 
Zeros out the imaginary part of the matrix

Logically Collective

Input Parameter:
- `mat` - the matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatImaginaryPart()`

# External Links
$(_doc_external("Mat/MatRealPart"))
"""
function MatRealPart(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatRealPart(petsclib::$UnionPetscLib, mat::AbstractPetscMat )

    @chk ccall(
               (:MatRealPart, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	MatRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Adds a new matrix type implementation that is usable as a `Mat` in PETSc

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - name of a new user-defined matrix type
- `function` - routine to create method context

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatType`, `MatSetType()`, `MatRegisterAll()`

# External Links
$(_doc_external("Mat/MatRegister"))
"""
function MatRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function MatRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:MatRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	MatRegisterRootName(petsclib::PetscLibType,rname::String, sname::String, mname::String) 
Registers a name that can be used for either a sequential or its corresponding parallel matrix type.

Input Parameters:
- `rname` - the rootname, for example, `MATAIJ`
- `sname` - the name of the sequential matrix type, for example, `MATSEQAIJ`
- `mname` - the name of the parallel matrix type, for example, `MATMPIAIJ`

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatType`, `PetscObjectBaseTypeCompare()`

# External Links
$(_doc_external("Mat/MatRegisterRootName"))
"""
function MatRegisterRootName(petsclib::PetscLibType, rname::String, sname::String, mname::String) end

@for_petsc function MatRegisterRootName(petsclib::$UnionPetscLib, rname::String, sname::String, mname::String )

    @chk ccall(
               (:MatRegisterRootName, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}),
               rname, sname, mname,
              )


	return nothing
end 

"""
	MatReorderForNonzeroDiagonal(petsclib::PetscLibType,mat::AbstractPetscMat, abstol::PetscReal, ris::AbstractIS, cis::AbstractIS) 
Changes matrix ordering to remove
zeros from diagonal. This may help in the `PCLU` factorization to
prevent a zero pivot.

Collective

Input Parameters:
- `mat`    - matrix to reorder
- `abstol` - absolute tolerance, it attempts to move all values smaller off the diagonal
- `ris`    - the row reordering
- `cis`    - the column reordering; this may be changed

Level: intermediate

Options Database Key:
- `-pc_factor_nonzeros_along_diagonal` - Reorder to remove zeros from diagonal

-seealso: `Mat`, `MatGetFactor()`, `MatGetOrdering()`

# External Links
$(_doc_external("Mat/MatReorderForNonzeroDiagonal"))
"""
function MatReorderForNonzeroDiagonal(petsclib::PetscLibType, mat::AbstractPetscMat, abstol::PetscReal, ris::AbstractIS, cis::AbstractIS) end

@for_petsc function MatReorderForNonzeroDiagonal(petsclib::$UnionPetscLib, mat::AbstractPetscMat, abstol::$PetscReal, ris::AbstractIS, cis::AbstractIS )

    @chk ccall(
               (:MatReorderForNonzeroDiagonal, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscReal, CIS, CIS),
               mat, abstol, ris, cis,
              )


	return nothing
end 

"""
	MatReorderingSeqSBAIJ(petsclib::PetscLibType,A::AbstractPetscMat, perm::AbstractIS) 

# External Links
$(_doc_external("Mat/MatReorderingSeqSBAIJ"))
"""
function MatReorderingSeqSBAIJ(petsclib::PetscLibType, A::AbstractPetscMat, perm::AbstractIS) end

@for_petsc function MatReorderingSeqSBAIJ(petsclib::$UnionPetscLib, A::AbstractPetscMat, perm::AbstractIS )

    @chk ccall(
               (:MatReorderingSeqSBAIJ, $petsc_library),
               PetscErrorCode,
               (CMat, CIS),
               A, perm,
              )


	return nothing
end 

"""
	MatResetHash(petsclib::PetscLibType,A::AbstractPetscMat) 
Reset the matrix so that it will use a hash table for the next round of `MatSetValues()` and `MatAssemblyBegin()`/`MatAssemblyEnd()`.

Collective

Input Parameter:
- `A` - the matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatResetPreallocation()`

# External Links
$(_doc_external("Mat/MatResetHash"))
"""
function MatResetHash(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatResetHash(petsclib::$UnionPetscLib, A::AbstractPetscMat )

    @chk ccall(
               (:MatResetHash, $petsc_library),
               PetscErrorCode,
               (CMat,),
               A,
              )


	return nothing
end 

"""
	MatResetPreallocation(petsclib::PetscLibType,A::AbstractPetscMat) 
Reset matrix to use the original preallocation values provided by the user, for example with `MatXAIJSetPreallocation()`

Collective

Input Parameter:
- `A` - the matrix

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatSeqAIJSetPreallocation()`, `MatMPIAIJSetPreallocation()`, `MatXAIJSetPreallocation()`

# External Links
$(_doc_external("Mat/MatResetPreallocation"))
"""
function MatResetPreallocation(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatResetPreallocation(petsclib::$UnionPetscLib, A::AbstractPetscMat )

    @chk ccall(
               (:MatResetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat,),
               A,
              )


	return nothing
end 

"""
	MatResidual(petsclib::PetscLibType,mat::AbstractPetscMat, b::AbstractPetscVec, x::AbstractPetscVec, r::AbstractPetscVec) 
Default routine to calculate the residual r = b

Collective

Input Parameters:
- `mat` - the matrix
- `b`   - the right-hand-side
- `x`   - the approximate solution

Output Parameter:
- `r` - location to store the residual

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatMult()`, `MatMultAdd()`, `PCMGSetResidual()`

# External Links
$(_doc_external("Mat/MatResidual"))
"""
function MatResidual(petsclib::PetscLibType, mat::AbstractPetscMat, b::AbstractPetscVec, x::AbstractPetscVec, r::AbstractPetscVec) end

@for_petsc function MatResidual(petsclib::$UnionPetscLib, mat::AbstractPetscMat, b::AbstractPetscVec, x::AbstractPetscVec, r::AbstractPetscVec )

    @chk ccall(
               (:MatResidual, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec, CVec),
               mat, b, x, r,
              )


	return nothing
end 

# override for MatRestoreColumnIJ; C signature: MatRestoreColumnIJ(Mat mat, PetscInt shift, PetscBool symmetric, PetscBool inodecompressed, PetscInt* n, PetscInt* ia[], PetscInt* ja[], PetscBool* done)
"""
	done::PetscBool = MatRestoreColumnIJ(petsclib::PetscLibType,mat::AbstractPetscMat, shift::PetscInt, symmetric::PetscBool, inodecompressed::PetscBool, ia::Vector{PetscInt}, ja::Vector{PetscInt}) 
Call after you are completed with the ia,ja indices obtained with `MatGetColumnIJ()`.

Collective

Input Parameters:
- `mat`             - the matrix
- `shift`           - 1 or zero indicating we want the indices starting at 0 or 1
- `symmetric`       - `PETSC_TRUE` or `PETSC_FALSE` indicating the matrix data structure should be symmetrized
- `inodecompressed` - `PETSC_TRUE` or `PETSC_FALSE` indicating if the nonzero structure of the
inodes or the nonzero elements is wanted. For `MATBAIJ` matrices the compressed version is
always used.

Output Parameters:
- `n`    - size of (possibly compressed) matrix
- `ia`   - the column pointers
- `ja`   - the row indices
- `done` - `PETSC_TRUE` or `PETSC_FALSE` indicated that the values have been returned

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatGetColumnIJ()`, `MatRestoreRowIJ()`

# External Links
$(_doc_external("Mat/MatRestoreColumnIJ"))
"""
function MatRestoreColumnIJ(petsclib::PetscLibType, mat::AbstractPetscMat, shift::PetscInt, symmetric::PetscBool, inodecompressed::PetscBool, ia::Vector{PetscInt}, ja::Vector{PetscInt}) end

# `ia`/`ja` are the arrays obtained from MatGetColumnIJ (their pointers are handed back to PETSc)
@for_petsc function MatRestoreColumnIJ(petsclib::$UnionPetscLib, mat::AbstractPetscMat, shift::$PetscInt, symmetric::PetscBool, inodecompressed::PetscBool, ia::Vector{$PetscInt}, ja::Vector{$PetscInt} )
	n_ = Ref{$PetscInt}()
	ia_ = Ref{Ptr{$PetscInt}}(isempty(ia) ? C_NULL : pointer(ia))
	ja_ = Ref{Ptr{$PetscInt}}(isempty(ja) ? C_NULL : pointer(ja))
	done_ = Ref{PetscBool}()

    @chk ccall(
               (:MatRestoreColumnIJ, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, PetscBool, PetscBool, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}, Ptr{PetscBool}),
               mat, shift, symmetric, inodecompressed, n_, ia_, ja_, done_,
              )

	return done_[]
end

"""
	MatRestoreLocalSubMatrix(petsclib::PetscLibType,mat::AbstractPetscMat, isrow::AbstractIS, iscol::AbstractIS, submat::AbstractPetscMat) 
Restores a reference to a submatrix specified in local numbering obtained with `MatGetLocalSubMatrix()`

Not Collective

Input Parameters:
- `mat`    - matrix to extract local submatrix from
- `isrow`  - local row indices for submatrix
- `iscol`  - local column indices for submatrix
- `submat` - the submatrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatGetLocalSubMatrix()`

# External Links
$(_doc_external("Mat/MatRestoreLocalSubMatrix"))
"""
function MatRestoreLocalSubMatrix(petsclib::PetscLibType, mat::AbstractPetscMat, isrow::AbstractIS, iscol::AbstractIS, submat::AbstractPetscMat) end

@for_petsc function MatRestoreLocalSubMatrix(petsclib::$UnionPetscLib, mat::AbstractPetscMat, isrow::AbstractIS, iscol::AbstractIS, submat::AbstractPetscMat )
	submat_ = Ref(submat.ptr)

    @chk ccall(
               (:MatRestoreLocalSubMatrix, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, CIS, Ptr{CMat}),
               mat, isrow, iscol, submat_,
              )

	submat.ptr = submat_[]

	return nothing
end 

"""
	MatRestoreNullSpaces(petsclib::PetscLibType,n::PetscInt, mat::Vector{<:AbstractPetscMat}, nullsp::AbstractArray{MatNullSpace}) 
sets the null spaces, transpose null spaces, and near null spaces obtained with `MatGetNullSpaces()` for an array of matrices

Logically Collective

Input Parameters:
- `n`      - the number of matrices
- `mat`    - the array of matrices
- `nullsp` - an array of null spaces

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `MatNullSpaceCreate()`, `MatSetNearNullSpace()`, `MatGetNullSpace()`, `MatSetTransposeNullSpace()`, `MatGetTransposeNullSpace()`,
`MatNullSpaceRemove()`, `MatGetNullSpaces()`

# External Links
$(_doc_external("Mat/MatRestoreNullSpaces"))
"""
function MatRestoreNullSpaces(petsclib::PetscLibType, n::PetscInt, mat::Vector{<:AbstractPetscMat}, nullsp::AbstractArray{MatNullSpace}) end

@for_petsc function MatRestoreNullSpaces(petsclib::$UnionPetscLib, n::$PetscInt, mat::Vector{<:AbstractPetscMat}, nullsp::AbstractArray{MatNullSpace} )
	nullsp_ = Ref(pointer(nullsp))

    @chk ccall(
               (:MatRestoreNullSpaces, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{CMat}, Ptr{Ptr{MatNullSpace}}),
               n, mat, nullsp_,
              )


	return nothing
end 

"""
	MatRestoreRow(petsclib::PetscLibType,mat::AbstractPetscMat, row::PetscInt, ncols::PetscInt, cols::AbstractArray{PetscInt}, vals::AbstractArray{PetscScalar}) 
Frees any temporary space allocated by `MatGetRow()`.

Not Collective

Input Parameters:
- `mat`   - the matrix
- `row`   - the row to get
- `ncols` - the number of nonzeros
- `cols`  - the columns of the nonzeros
- `vals`  - if nonzero the column values

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatGetRow()`

# External Links
$(_doc_external("Mat/MatRestoreRow"))
"""
function MatRestoreRow(petsclib::PetscLibType, mat::AbstractPetscMat, row::PetscInt, ncols::PetscInt, cols::AbstractArray{PetscInt}, vals::AbstractArray{PetscScalar}) end

@for_petsc function MatRestoreRow(petsclib::$UnionPetscLib, mat::AbstractPetscMat, row::$PetscInt, ncols::$PetscInt, cols::AbstractArray{$PetscInt}, vals::AbstractArray{$PetscScalar} )
	ncols_ = Ref{$PetscInt}(ncols)
	cols_ = Ref(pointer(cols))
	vals_ = Ref(pointer(vals))

    @chk ccall(
               (:MatRestoreRow, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscScalar}}),
               mat, row, ncols_, cols_, vals_,
              )


	return nothing
end 

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

"""
	MatRestoreRowUpperTriangular(petsclib::PetscLibType,mat::AbstractPetscMat) 
Disable calls to `MatGetRow()` for matrix in `MATSBAIJ` format.

Not Collective

Input Parameter:
- `mat` - the matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATSBAIJ`, `MatGetRowUpperTriangular()`

# External Links
$(_doc_external("Mat/MatRestoreRowUpperTriangular"))
"""
function MatRestoreRowUpperTriangular(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatRestoreRowUpperTriangular(petsclib::$UnionPetscLib, mat::AbstractPetscMat )

    @chk ccall(
               (:MatRestoreRowUpperTriangular, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	MatRestrict(petsclib::PetscLibType,A::AbstractPetscMat, x::AbstractPetscVec, y::AbstractPetscVec) 
y = A*x or A^T*x

Neighbor-wise Collective

Input Parameters:
- `A` - the matrix
- `x` - the vector to be restricted

Output Parameter:
- `y` - the resulting vector

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatMultAdd()`, `MatMultTransposeAdd()`, `MatInterpolate()`, `PCMG`

# External Links
$(_doc_external("Mat/MatRestrict"))
"""
function MatRestrict(petsclib::PetscLibType, A::AbstractPetscMat, x::AbstractPetscVec, y::AbstractPetscVec) end

@for_petsc function MatRestrict(petsclib::$UnionPetscLib, A::AbstractPetscMat, x::AbstractPetscVec, y::AbstractPetscVec )

    @chk ccall(
               (:MatRestrict, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               A, x, y,
              )


	return nothing
end 

"""
	MatRetrieveValues(petsclib::PetscLibType,mat::AbstractPetscMat) 
Retrieves the copy of the matrix values that was stored with `MatStoreValues()`

Logically Collect

Input Parameter:
- `mat` - the matrix (currently only `MATAIJ` matrices support this option)

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatStoreValues()`

# External Links
$(_doc_external("Mat/MatRetrieveValues"))
"""
function MatRetrieveValues(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatRetrieveValues(petsclib::$UnionPetscLib, mat::AbstractPetscMat )

    @chk ccall(
               (:MatRetrieveValues, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	use::PetscBool = MatSNESMFGetReuseBase(petsclib::PetscLibType,J::AbstractPetscMat) 
Determines if the base vector is to be used for differencing even if the function provided to `SNESSetFunction()` is not the
same as that provided to `MatMFFDSetFunction()`.

Logically Collective

Input Parameter:
- `J` - the `MATMFFD` matrix

Output Parameter:
- `use` - if true always reuse the base vector instead of recomputing f(u) even if the function in the `MATMFFD` is
not `SNESComputeFunction()`

Level: advanced

-seealso: [](ch_snes), `Mat`, `SNES`, `MatSNESMFSetReuseBase()`, `MatCreateSNESMF()`

# External Links
$(_doc_external("SNES/MatSNESMFGetReuseBase"))
"""
function MatSNESMFGetReuseBase(petsclib::PetscLibType, J::AbstractPetscMat) end

@for_petsc function MatSNESMFGetReuseBase(petsclib::$UnionPetscLib, J::AbstractPetscMat )
	use_ = Ref{PetscBool}()

    @chk ccall(
               (:MatSNESMFGetReuseBase, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscBool}),
               J, use_,
              )

	use = use_[]

	return use
end 

"""
	snes::PetscSNES = MatSNESMFGetSNES(petsclib::PetscLibType,J::AbstractPetscMat) 
returns the `SNES` associated with a matrix created with `MatCreateSNESMF()`

Not Collective

Input Parameter:
- `J` - the matrix

Output Parameter:
- `snes` - the `SNES` object

Level: advanced

-seealso: [](ch_snes), `Mat`, `SNES`, `MatCreateSNESMF()`

# External Links
$(_doc_external("SNES/MatSNESMFGetSNES"))
"""
function MatSNESMFGetSNES(petsclib::PetscLibType, J::AbstractPetscMat) end

@for_petsc function MatSNESMFGetSNES(petsclib::$UnionPetscLib, J::AbstractPetscMat )
	snes_ = Ref{CSNES}()

    @chk ccall(
               (:MatSNESMFGetSNES, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CSNES}),
               J, snes_,
              )

	snes = PetscSNES(snes_[], petsclib)

	return snes
end 

"""
	MatSNESMFMoreSetParameters(petsclib::PetscLibType,mat::AbstractPetscMat, error::PetscReal, umin::PetscReal, h::PetscReal) 
Sets the parameters for the approximation of
matrix-vector products using finite differences, see  `MatCreateSNESMFMore()`

Input Parameters:
- `mat`   - the matrix
- `error` - relative error (should be set to the square root of the relative error in the function evaluations)
- `umin`  - minimum allowable u-value
- `h`     - differencing parameter

Options Database Keys:
- `-snes_mf_err <error_rel>` - see `MatCreateSNESMF()`
- `-snes_mf_umin <umin>`     - see `MatCreateSNESMF()`
- `-snes_mf_compute_err`     - compute the square root or relative error in function
- `-snes_mf_freq_err <freq>` - set the frequency to recompute the parameters
- `-snes_mf_jorge`           - use the method of Jorge More

Level: advanced

-seealso: [](ch_snes), `SNES`, `MatCreateSNESMF()`, `MatCreateSNESMFMore()`

# External Links
$(_doc_external("SNES/MatSNESMFMoreSetParameters"))
"""
function MatSNESMFMoreSetParameters(petsclib::PetscLibType, mat::AbstractPetscMat, error::PetscReal, umin::PetscReal, h::PetscReal) end

@for_petsc function MatSNESMFMoreSetParameters(petsclib::$UnionPetscLib, mat::AbstractPetscMat, error::$PetscReal, umin::$PetscReal, h::$PetscReal )

    @chk ccall(
               (:MatSNESMFMoreSetParameters, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscReal, $PetscReal, $PetscReal),
               mat, error, umin, h,
              )


	return nothing
end 

"""
	MatSNESMFSetReuseBase(petsclib::PetscLibType,J::AbstractPetscMat, use::PetscBool) 
Causes the base vector to be used for differencing even if the function provided to `SNESSetFunction()` is not the
same as that provided to `MatMFFDSetFunction()`.

Logically Collective

Input Parameters:
- `J`   - the `MATMFFD` matrix
- `use` - if true always reuse the base vector instead of recomputing f(u) even if the function in the `MATMFFD` is
not `SNESComputeFunction()`

Level: advanced

-seealso: [](ch_snes), `SNES`, `MATMFFD`, `MatMFFDSetFunction()`, `SNESSetFunction()`, `MatCreateSNESMF()`, `MatSNESMFGetReuseBase()`

# External Links
$(_doc_external("SNES/MatSNESMFSetReuseBase"))
"""
function MatSNESMFSetReuseBase(petsclib::PetscLibType, J::AbstractPetscMat, use::PetscBool) end

@for_petsc function MatSNESMFSetReuseBase(petsclib::$UnionPetscLib, J::AbstractPetscMat, use::PetscBool )

    @chk ccall(
               (:MatSNESMFSetReuseBase, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               J, use,
              )


	return nothing
end 

"""
	MatSOR(petsclib::PetscLibType,mat::AbstractPetscMat, b::AbstractPetscVec, omega::PetscReal, flag::MatSORType, shift::PetscReal, its::PetscInt, lits::PetscInt, x::AbstractPetscVec) 
Computes relaxation (SOR, Gauss

Neighbor-wise Collective

Input Parameters:
- `mat`   - the matrix
- `b`     - the right-hand side
- `omega` - the relaxation factor
- `flag`  - flag indicating the type of SOR (see below)
- `shift` - diagonal shift
- `its`   - the number of iterations
- `lits`  - the number of local iterations

Output Parameter:
- `x` - the solution (can contain an initial guess, use option `SOR_ZERO_INITIAL_GUESS` to indicate no guess)

SOR Flags:
- `SOR_FORWARD_SWEEP`     - forward SOR
- `SOR_BACKWARD_SWEEP`     - backward SOR
- `SOR_SYMMETRIC_SWEEP`     - SSOR (symmetric SOR)
- `SOR_LOCAL_FORWARD_SWEEP`     - local forward SOR
- `SOR_LOCAL_BACKWARD_SWEEP`     - local forward SOR
- `SOR_LOCAL_SYMMETRIC_SWEEP`     - local SSOR
- `SOR_EISENSTAT`     - SOR with Eisenstat trick
- `SOR_APPLY_UPPER`, `SOR_APPLY_LOWER`     - applies
upper/lower triangular part of matrix to
vector (with omega)
- `SOR_ZERO_INITIAL_GUESS`     - zero initial guess

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatMult()`, `KSP`, `PC`, `MatGetFactor()`

# External Links
$(_doc_external("Mat/MatSOR"))
"""
function MatSOR(petsclib::PetscLibType, mat::AbstractPetscMat, b::AbstractPetscVec, omega::PetscReal, flag::MatSORType, shift::PetscReal, its::PetscInt, lits::PetscInt, x::AbstractPetscVec) end

@for_petsc function MatSOR(petsclib::$UnionPetscLib, mat::AbstractPetscMat, b::AbstractPetscVec, omega::$PetscReal, flag::MatSORType, shift::$PetscReal, its::$PetscInt, lits::$PetscInt, x::AbstractPetscVec )

    @chk ccall(
               (:MatSOR, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, $PetscReal, MatSORType, $PetscReal, $PetscInt, $PetscInt, CVec),
               mat, b, omega, flag, shift, its, lits, x,
              )


	return nothing
end 

"""
	cperm::PetscBool = MatSTRUMPACKGetColPerm(petsclib::PetscLibType,F::AbstractPetscMat) 
Get whether STRUMPACK <https://portal.nersc.gov/project/sparse/strumpack/master>
will try to permute the columns of the matrix in order to get a nonzero diagonal

Logically Collective

Input Parameters:
- `F` - the factored matrix obtained by calling `MatGetFactor()`

Output Parameter:
- `cperm` - Indicates whether STRUMPACK will permute columns

Level: intermediate

-seealso: `MATSOLVERSTRUMPACK`, [](ch_matrices), `MatSTRUMPACKSetReordering()`, `Mat`, `MatGetFactor()`, `MatSTRUMPACKSetColPerm()`

# External Links
$(_doc_external("Mat/MatSTRUMPACKGetColPerm"))
"""
function MatSTRUMPACKGetColPerm(petsclib::PetscLibType, F::AbstractPetscMat) end

@for_petsc function MatSTRUMPACKGetColPerm(petsclib::$UnionPetscLib, F::AbstractPetscMat )
	cperm_ = Ref{PetscBool}()

    @chk ccall(
               (:MatSTRUMPACKGetColPerm, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscBool}),
               F, cperm_,
              )

	cperm = cperm_[]

	return cperm
end 

"""
	atol::PetscReal = MatSTRUMPACKGetCompAbsTol(petsclib::PetscLibType,F::AbstractPetscMat) 
Get STRUMPACK <https://portal.nersc.gov/project/sparse/strumpack/master> absolute tolerance for compression

Logically Collective

Input Parameters:
- `F` - the factored matrix obtained by calling `MatGetFactor()`

Output Parameter:
- `atol` - absolute compression tolerance

Level: intermediate

-seealso: `MATSOLVERSTRUMPACK`, [](ch_matrices), `Mat`, `MatGetFactor()`, `MatSTRUMPACKSetCompAbsTol()`, `MatSTRUMPACKSetReordering()`, `MatSTRUMPACKSetColPerm()`

# External Links
$(_doc_external("Mat/MatSTRUMPACKGetCompAbsTol"))
"""
function MatSTRUMPACKGetCompAbsTol(petsclib::PetscLibType, F::AbstractPetscMat) end

@for_petsc function MatSTRUMPACKGetCompAbsTol(petsclib::$UnionPetscLib, F::AbstractPetscMat )
	atol_ = Ref{$PetscReal}()

    @chk ccall(
               (:MatSTRUMPACKGetCompAbsTol, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscReal}),
               F, atol_,
              )

	atol = atol_[]

	return atol
end 

"""
	bfly_lvls::PetscInt = MatSTRUMPACKGetCompButterflyLevels(petsclib::PetscLibType,F::AbstractPetscMat) 
Get STRUMPACK <https://portal.nersc.gov/project/sparse/strumpack/master>
number of butterfly levels in HODLR compression (requires ButterflyPACK support)

Logically Collective

Input Parameters:
- `F` - the factored matrix obtained by calling `MatGetFactor()` from PETSc-STRUMPACK interface

Output Parameter:
- `bfly_lvls` - Number of levels of butterfly compression in HODLR compression

Level: intermediate

-seealso: `MATSOLVERSTRUMPACK`, `MatGetFactor()`, `MatSTRUMPACKSetCompButterflyLevels()`

# External Links
$(_doc_external("Mat/MatSTRUMPACKGetCompButterflyLevels"))
"""
function MatSTRUMPACKGetCompButterflyLevels(petsclib::PetscLibType, F::AbstractPetscMat) end

@for_petsc function MatSTRUMPACKGetCompButterflyLevels(petsclib::$UnionPetscLib, F::AbstractPetscMat )
	bfly_lvls_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatSTRUMPACKGetCompButterflyLevels, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}),
               F, bfly_lvls_,
              )

	bfly_lvls = bfly_lvls_[]

	return bfly_lvls
end 

"""
	leaf_size::PetscInt = MatSTRUMPACKGetCompLeafSize(petsclib::PetscLibType,F::AbstractPetscMat) 
Get STRUMPACK <https://portal.nersc.gov/project/sparse/strumpack/master> leaf size for HSS, BLR, HODLR...

Logically Collective

Input Parameters:
- `F` - the factored matrix obtained by calling `MatGetFactor()` from PETSc-STRUMPACK interface

Output Parameter:
- `leaf_size` - Size of diagonal blocks in rank-structured approximation

Level: intermediate

-seealso: `MATSOLVERSTRUMPACK`, [](ch_matrices), `Mat`, `MatGetFactor()`, `MatSTRUMPACKSetCompLeafSize()`, `MatSTRUMPACKSetReordering()`, `MatSTRUMPACKSetColPerm()`

# External Links
$(_doc_external("Mat/MatSTRUMPACKGetCompLeafSize"))
"""
function MatSTRUMPACKGetCompLeafSize(petsclib::PetscLibType, F::AbstractPetscMat) end

@for_petsc function MatSTRUMPACKGetCompLeafSize(petsclib::$UnionPetscLib, F::AbstractPetscMat )
	leaf_size_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatSTRUMPACKGetCompLeafSize, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}),
               F, leaf_size_,
              )

	leaf_size = leaf_size_[]

	return leaf_size
end 

"""
	lossy_prec::PetscInt = MatSTRUMPACKGetCompLossyPrecision(petsclib::PetscLibType,F::AbstractPetscMat) 
Get STRUMPACK <https://portal.nersc.gov/project/sparse/strumpack/master> precision for lossy compression (requires ZFP support)

Logically Collective

Input Parameters:
- `F` - the factored matrix obtained by calling `MatGetFactor()` from PETSc-STRUMPACK interface

Output Parameter:
- `lossy_prec` - Number of bitplanes to use in lossy compression

Level: intermediate

-seealso: `MATSOLVERSTRUMPACK`, `MatGetFactor()`, `MatSTRUMPACKSetCompLossyPrecision()`

# External Links
$(_doc_external("Mat/MatSTRUMPACKGetCompLossyPrecision"))
"""
function MatSTRUMPACKGetCompLossyPrecision(petsclib::PetscLibType, F::AbstractPetscMat) end

@for_petsc function MatSTRUMPACKGetCompLossyPrecision(petsclib::$UnionPetscLib, F::AbstractPetscMat )
	lossy_prec_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatSTRUMPACKGetCompLossyPrecision, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}),
               F, lossy_prec_,
              )

	lossy_prec = lossy_prec_[]

	return lossy_prec
end 

"""
	min_sep_size::PetscInt = MatSTRUMPACKGetCompMinSepSize(petsclib::PetscLibType,F::AbstractPetscMat) 
Get STRUMPACK <https://portal.nersc.gov/project/sparse/strumpack/master> minimum separator size for low

Logically Collective

Input Parameters:
- `F` - the factored matrix obtained by calling `MatGetFactor()` from PETSc-STRUMPACK interface

Output Parameter:
- `min_sep_size` - minimum dense matrix size for low-rank approximation

Level: intermediate

-seealso: `MATSOLVERSTRUMPACK`, `MatGetFactor()`, `MatSTRUMPACKSetCompMinSepSize()`

# External Links
$(_doc_external("Mat/MatSTRUMPACKGetCompMinSepSize"))
"""
function MatSTRUMPACKGetCompMinSepSize(petsclib::PetscLibType, F::AbstractPetscMat) end

@for_petsc function MatSTRUMPACKGetCompMinSepSize(petsclib::$UnionPetscLib, F::AbstractPetscMat )
	min_sep_size_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatSTRUMPACKGetCompMinSepSize, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}),
               F, min_sep_size_,
              )

	min_sep_size = min_sep_size_[]

	return min_sep_size
end 

"""
	rtol::PetscReal = MatSTRUMPACKGetCompRelTol(petsclib::PetscLibType,F::AbstractPetscMat) 
Get STRUMPACK <https://portal.nersc.gov/project/sparse/strumpack/master> relative tolerance for compression

Logically Collective

Input Parameters:
- `F` - the factored matrix obtained by calling `MatGetFactor()`

Output Parameter:
- `rtol` - relative compression tolerance

Level: intermediate

-seealso: `MATSOLVERSTRUMPACK`, [](ch_matrices), `Mat`, `MatGetFactor()`, `MatSTRUMPACKSetCompRelTol()`, `MatSTRUMPACKSetReordering()`, `MatSTRUMPACKSetColPerm()`

# External Links
$(_doc_external("Mat/MatSTRUMPACKGetCompRelTol"))
"""
function MatSTRUMPACKGetCompRelTol(petsclib::PetscLibType, F::AbstractPetscMat) end

@for_petsc function MatSTRUMPACKGetCompRelTol(petsclib::$UnionPetscLib, F::AbstractPetscMat )
	rtol_ = Ref{$PetscReal}()

    @chk ccall(
               (:MatSTRUMPACKGetCompRelTol, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscReal}),
               F, rtol_,
              )

	rtol = rtol_[]

	return rtol
end 

"""
	comp::MatSTRUMPACKCompressionType = MatSTRUMPACKGetCompression(petsclib::PetscLibType,F::AbstractPetscMat) 
Get STRUMPACK <https://portal.nersc.gov/project/sparse/strumpack/master> compression type

Input Parameters:
- `F` - the factored matrix obtained by calling `MatGetFactor()` from PETSc-STRUMPACK interface

Output Parameter:
- `comp` - Type of compression to be used in the approximate sparse factorization

Level: intermediate

-seealso: `MATSOLVERSTRUMPACK`, `MatGetFactor()`, `MatSTRUMPACKCompressionType`, `MatSTRUMPACKSetCompression()`

# External Links
$(_doc_external("Mat/MatSTRUMPACKGetCompression"))
"""
function MatSTRUMPACKGetCompression(petsclib::PetscLibType, F::AbstractPetscMat) end

@for_petsc function MatSTRUMPACKGetCompression(petsclib::$UnionPetscLib, F::AbstractPetscMat )
	comp_ = Ref{MatSTRUMPACKCompressionType}()

    @chk ccall(
               (:MatSTRUMPACKGetCompression, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatSTRUMPACKCompressionType}),
               F, comp_,
              )

	comp = comp_[]

	return comp
end 

"""
	gpu::PetscBool = MatSTRUMPACKGetGPU(petsclib::PetscLibType,F::AbstractPetscMat) 
Get whether STRUMPACK <https://portal.nersc.gov/project/sparse/strumpack/master>
will try to use GPU acceleration (not supported for all compression types)

Logically Collective

Input Parameters:
- `F` - the factored matrix obtained by calling `MatGetFactor()` from PETSc-STRUMPACK interface

Output Parameter:
- `gpu` - whether or not STRUMPACK will try to use GPU acceleration

Level: intermediate

-seealso: `MATSOLVERSTRUMPACK`, `MatGetFactor()`, `MatSTRUMPACKSetGPU()`

# External Links
$(_doc_external("Mat/MatSTRUMPACKGetGPU"))
"""
function MatSTRUMPACKGetGPU(petsclib::PetscLibType, F::AbstractPetscMat) end

@for_petsc function MatSTRUMPACKGetGPU(petsclib::$UnionPetscLib, F::AbstractPetscMat )
	gpu_ = Ref{PetscBool}()

    @chk ccall(
               (:MatSTRUMPACKGetGPU, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscBool}),
               F, gpu_,
              )

	gpu = gpu_[]

	return gpu
end 

"""
	reordering::MatSTRUMPACKReordering = MatSTRUMPACKGetReordering(petsclib::PetscLibType,F::AbstractPetscMat) 
Get STRUMPACK <https://portal.nersc.gov/project/sparse/strumpack/master> fill

Logically Collective

Input Parameters:
- `F` - the factored matrix obtained by calling `MatGetFactor()` from PETSc-STRUMPACK interface

Output Parameter:
- `reordering` - the code to be used to find the fill-reducing reordering

Level: intermediate

-seealso: `MATSOLVERSTRUMPACK`, [](ch_matrices), `Mat`, `MatSTRUMPACKReordering`, `MatGetFactor()`, `MatSTRUMPACKSetReordering()`, `MatSTRUMPACKSetColPerm()`

# External Links
$(_doc_external("Mat/MatSTRUMPACKGetReordering"))
"""
function MatSTRUMPACKGetReordering(petsclib::PetscLibType, F::AbstractPetscMat) end

@for_petsc function MatSTRUMPACKGetReordering(petsclib::$UnionPetscLib, F::AbstractPetscMat )
	reordering_ = Ref{MatSTRUMPACKReordering}()

    @chk ccall(
               (:MatSTRUMPACKGetReordering, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatSTRUMPACKReordering}),
               F, reordering_,
              )

	reordering = reordering_[]

	return reordering
end 

"""
	MatSTRUMPACKSetColPerm(petsclib::PetscLibType,F::AbstractPetscMat, cperm::PetscBool) 
Set whether STRUMPACK <https://portal.nersc.gov/project/sparse/strumpack/master>
should try to permute the columns of the matrix in order to get a nonzero diagonal

Logically Collective

Input Parameters:
- `F`     - the factored matrix obtained by calling `MatGetFactor()`
- `cperm` - `PETSC_TRUE` to permute (internally) the columns of the matrix

Options Database Key:
- `-mat_strumpack_colperm <cperm>` - true to use the permutation

Level: intermediate

-seealso: `MATSOLVERSTRUMPACK`, [](ch_matrices), `MatSTRUMPACKSetReordering()`, `Mat`, `MatGetFactor()`, `MatSTRUMPACKGetColPerm()`

# External Links
$(_doc_external("Mat/MatSTRUMPACKSetColPerm"))
"""
function MatSTRUMPACKSetColPerm(petsclib::PetscLibType, F::AbstractPetscMat, cperm::PetscBool) end

@for_petsc function MatSTRUMPACKSetColPerm(petsclib::$UnionPetscLib, F::AbstractPetscMat, cperm::PetscBool )

    @chk ccall(
               (:MatSTRUMPACKSetColPerm, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               F, cperm,
              )


	return nothing
end 

"""
	MatSTRUMPACKSetCompAbsTol(petsclib::PetscLibType,F::AbstractPetscMat, atol::PetscReal) 
Set STRUMPACK <https://portal.nersc.gov/project/sparse/strumpack/master> absolute tolerance for compression

Logically Collective

Input Parameters:
- `F`    - the factored matrix obtained by calling `MatGetFactor()`
- `atol` - absolute compression tolerance

Options Database Key:
- `-mat_strumpack_compression_abs_tol <1e-10>` - Absolute compression tolerance, when using `-pctype ilu`

Level: intermediate

-seealso: `MATSOLVERSTRUMPACK`, [](ch_matrices), `Mat`, `MatGetFactor()`, `MatSTRUMPACKGetCompAbsTol()`, `MatSTRUMPACKSetReordering()`, `MatSTRUMPACKSetColPerm()`

# External Links
$(_doc_external("Mat/MatSTRUMPACKSetCompAbsTol"))
"""
function MatSTRUMPACKSetCompAbsTol(petsclib::PetscLibType, F::AbstractPetscMat, atol::PetscReal) end

@for_petsc function MatSTRUMPACKSetCompAbsTol(petsclib::$UnionPetscLib, F::AbstractPetscMat, atol::$PetscReal )

    @chk ccall(
               (:MatSTRUMPACKSetCompAbsTol, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscReal),
               F, atol,
              )


	return nothing
end 

"""
	MatSTRUMPACKSetCompButterflyLevels(petsclib::PetscLibType,F::AbstractPetscMat, bfly_lvls::PetscInt) 
Set STRUMPACK <https://portal.nersc.gov/project/sparse/strumpack/master>
number of butterfly levels in HODLR compression (requires ButterflyPACK support)

Logically Collective

Input Parameters:
- `F`         - the factored matrix obtained by calling `MatGetFactor()` from PETSc-STRUMPACK interface
- `bfly_lvls` - Number of levels of butterfly compression in HODLR compression

Options Database Key:
- `-mat_strumpack_compression_butterfly_levels <bfly_lvls>` - Number of levels in the hierarchically off-diagonal matrix for which to use butterfly,
when using `-pctype ilu`, (BLR_)HODLR compression

Level: intermediate

-seealso: `MATSOLVERSTRUMPACK`, `MatGetFactor()`, `MatSTRUMPACKGetCompButterflyLevels()`

# External Links
$(_doc_external("Mat/MatSTRUMPACKSetCompButterflyLevels"))
"""
function MatSTRUMPACKSetCompButterflyLevels(petsclib::PetscLibType, F::AbstractPetscMat, bfly_lvls::PetscInt) end

@for_petsc function MatSTRUMPACKSetCompButterflyLevels(petsclib::$UnionPetscLib, F::AbstractPetscMat, bfly_lvls::$PetscInt )

    @chk ccall(
               (:MatSTRUMPACKSetCompButterflyLevels, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt),
               F, bfly_lvls,
              )


	return nothing
end 

"""
	MatSTRUMPACKSetCompLeafSize(petsclib::PetscLibType,F::AbstractPetscMat, leaf_size::PetscInt) 
Set STRUMPACK <https://portal.nersc.gov/project/sparse/strumpack/master> leaf size for HSS, BLR, HODLR...

Logically Collective

Input Parameters:
- `F`         - the factored matrix obtained by calling `MatGetFactor()` from PETSc-STRUMPACK interface
- `leaf_size` - Size of diagonal blocks in rank-structured approximation

Options Database Key:
- `-mat_strumpack_compression_leaf_size` - Size of diagonal blocks in rank-structured approximation, when using `-pctype ilu`

Level: intermediate

-seealso: `MATSOLVERSTRUMPACK`, [](ch_matrices), `Mat`, `MatGetFactor()`, `MatSTRUMPACKGetCompLeafSize()`, `MatSTRUMPACKSetReordering()`, `MatSTRUMPACKSetColPerm()`

# External Links
$(_doc_external("Mat/MatSTRUMPACKSetCompLeafSize"))
"""
function MatSTRUMPACKSetCompLeafSize(petsclib::PetscLibType, F::AbstractPetscMat, leaf_size::PetscInt) end

@for_petsc function MatSTRUMPACKSetCompLeafSize(petsclib::$UnionPetscLib, F::AbstractPetscMat, leaf_size::$PetscInt )

    @chk ccall(
               (:MatSTRUMPACKSetCompLeafSize, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt),
               F, leaf_size,
              )


	return nothing
end 

"""
	MatSTRUMPACKSetCompLossyPrecision(petsclib::PetscLibType,F::AbstractPetscMat, lossy_prec::PetscInt) 
Set STRUMPACK <https://portal.nersc.gov/project/sparse/strumpack/master> precision for lossy compression (requires ZFP support)

Logically Collective

Input Parameters:
- `F`          - the factored matrix obtained by calling `MatGetFactor()` from PETSc-STRUMPACK interface
- `lossy_prec` - Number of bitplanes to use in lossy compression

Options Database Key:
- `-mat_strumpack_compression_lossy_precision <lossy_prec>` - Precision when using lossy compression [1-64], when using `-pctype ilu -mat_strumpack_compression MAT_STRUMPACK_COMPRESSION_TYPE_LOSSY`

Level: intermediate

-seealso: `MATSOLVERSTRUMPACK`, `MatGetFactor()`, `MatSTRUMPACKGetCompLossyPrecision()`

# External Links
$(_doc_external("Mat/MatSTRUMPACKSetCompLossyPrecision"))
"""
function MatSTRUMPACKSetCompLossyPrecision(petsclib::PetscLibType, F::AbstractPetscMat, lossy_prec::PetscInt) end

@for_petsc function MatSTRUMPACKSetCompLossyPrecision(petsclib::$UnionPetscLib, F::AbstractPetscMat, lossy_prec::$PetscInt )

    @chk ccall(
               (:MatSTRUMPACKSetCompLossyPrecision, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt),
               F, lossy_prec,
              )


	return nothing
end 

"""
	MatSTRUMPACKSetCompMinSepSize(petsclib::PetscLibType,F::AbstractPetscMat, min_sep_size::PetscInt) 
Set STRUMPACK <https://portal.nersc.gov/project/sparse/strumpack/master> minimum separator size for low

Logically Collective

Input Parameters:
- `F`            - the factored matrix obtained by calling `MatGetFactor()` from PETSc-STRUMPACK interface
- `min_sep_size` - minimum dense matrix size for low-rank approximation

Options Database Key:
- `-mat_strumpack_compression_min_sep_size <min_sep_size>` - Minimum size of dense sub-block for low-rank compression

Level: intermediate

-seealso: `MATSOLVERSTRUMPACK`, `MatGetFactor()`, `MatSTRUMPACKGetCompMinSepSize()`

# External Links
$(_doc_external("Mat/MatSTRUMPACKSetCompMinSepSize"))
"""
function MatSTRUMPACKSetCompMinSepSize(petsclib::PetscLibType, F::AbstractPetscMat, min_sep_size::PetscInt) end

@for_petsc function MatSTRUMPACKSetCompMinSepSize(petsclib::$UnionPetscLib, F::AbstractPetscMat, min_sep_size::$PetscInt )

    @chk ccall(
               (:MatSTRUMPACKSetCompMinSepSize, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt),
               F, min_sep_size,
              )


	return nothing
end 

"""
	MatSTRUMPACKSetCompRelTol(petsclib::PetscLibType,F::AbstractPetscMat, rtol::PetscReal) 
Set STRUMPACK <https://portal.nersc.gov/project/sparse/strumpack/master> relative tolerance for compression

Logically Collective

Input Parameters:
- `F`    - the factored matrix obtained by calling `MatGetFactor()`
- `rtol` - relative compression tolerance

Options Database Key:
- `-mat_strumpack_compression_rel_tol <1e-4>` - Relative compression tolerance, when using `-pctype ilu`

Level: intermediate

-seealso: `MATSOLVERSTRUMPACK`, [](ch_matrices), `Mat`, `MatGetFactor()`, `MatSTRUMPACKGetCompRelTol()`, `MatSTRUMPACKSetReordering()`, `MatSTRUMPACKSetColPerm()`

# External Links
$(_doc_external("Mat/MatSTRUMPACKSetCompRelTol"))
"""
function MatSTRUMPACKSetCompRelTol(petsclib::PetscLibType, F::AbstractPetscMat, rtol::PetscReal) end

@for_petsc function MatSTRUMPACKSetCompRelTol(petsclib::$UnionPetscLib, F::AbstractPetscMat, rtol::$PetscReal )

    @chk ccall(
               (:MatSTRUMPACKSetCompRelTol, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscReal),
               F, rtol,
              )


	return nothing
end 

"""
	MatSTRUMPACKSetCompression(petsclib::PetscLibType,F::AbstractPetscMat, comp::MatSTRUMPACKCompressionType) 
Set STRUMPACK <https://portal.nersc.gov/project/sparse/strumpack/master> compression type

Input Parameters:
- `F`    - the factored matrix obtained by calling `MatGetFactor()` from PETSc-STRUMPACK interface
- `comp` - Type of compression to be used in the approximate sparse factorization

Options Database Key:
- `-mat_strumpack_compression <NONE>` - Type of rank-structured compression in sparse LU factors (choose one of) NONE HSS BLR HODLR BLR_HODLR ZFP_BLR_HODLR LOSSLESS LOSSY

Level: intermediate

-seealso: `MATSOLVERSTRUMPACK`, `MatGetFactor()`, `MatSTRUMPACKCompressionType`, `MatSTRUMPACKGetCompression()`

# External Links
$(_doc_external("Mat/MatSTRUMPACKSetCompression"))
"""
function MatSTRUMPACKSetCompression(petsclib::PetscLibType, F::AbstractPetscMat, comp::MatSTRUMPACKCompressionType) end

@for_petsc function MatSTRUMPACKSetCompression(petsclib::$UnionPetscLib, F::AbstractPetscMat, comp::MatSTRUMPACKCompressionType )

    @chk ccall(
               (:MatSTRUMPACKSetCompression, $petsc_library),
               PetscErrorCode,
               (CMat, MatSTRUMPACKCompressionType),
               F, comp,
              )


	return nothing
end 

"""
	MatSTRUMPACKSetGPU(petsclib::PetscLibType,F::AbstractPetscMat, gpu::PetscBool) 
Set whether STRUMPACK <https://portal.nersc.gov/project/sparse/strumpack/master>
should enable GPU acceleration (not supported for all compression types)

Logically Collective

Input Parameters:
- `F`   - the factored matrix obtained by calling `MatGetFactor()` from PETSc-STRUMPACK interface
- `gpu` - whether or not to use GPU acceleration

Options Database Key:
- `-mat_strumpack_gpu <gpu>` - true to use gpu offload

Level: intermediate

-seealso: `MATSOLVERSTRUMPACK`, `MatGetFactor()`, `MatSTRUMPACKGetGPU()`

# External Links
$(_doc_external("Mat/MatSTRUMPACKSetGPU"))
"""
function MatSTRUMPACKSetGPU(petsclib::PetscLibType, F::AbstractPetscMat, gpu::PetscBool) end

@for_petsc function MatSTRUMPACKSetGPU(petsclib::$UnionPetscLib, F::AbstractPetscMat, gpu::PetscBool )

    @chk ccall(
               (:MatSTRUMPACKSetGPU, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               F, gpu,
              )


	return nothing
end 

"""
	MatSTRUMPACKSetGeometricComponents(petsclib::PetscLibType,F::AbstractPetscMat, nc::PetscInt) 
Set STRUMPACK <https://portal.nersc.gov/project/sparse/strumpack/master>
number of degrees of freedom per mesh point, for use with GEOMETRIC ordering.

Logically Collective

Input Parameters:
- `F`  - the factored matrix obtained by calling `MatGetFactor()` from PETSc-STRUMPACK interface
- `nc` - Number of components/dof's per grid point

Options Database Key:
- `-mat_strumpack_geometric_components <1>` - Number of components per mesh point, for geometric nested dissection ordering

Level: intermediate

-seealso: `MATSOLVERSTRUMPACK`, `MatGetFactor()`

# External Links
$(_doc_external("Mat/MatSTRUMPACKSetGeometricComponents"))
"""
function MatSTRUMPACKSetGeometricComponents(petsclib::PetscLibType, F::AbstractPetscMat, nc::PetscInt) end

@for_petsc function MatSTRUMPACKSetGeometricComponents(petsclib::$UnionPetscLib, F::AbstractPetscMat, nc::$PetscInt )

    @chk ccall(
               (:MatSTRUMPACKSetGeometricComponents, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt),
               F, nc,
              )


	return nothing
end 

"""
	MatSTRUMPACKSetGeometricNxyz(petsclib::PetscLibType,F::AbstractPetscMat, nx::PetscInt, ny::PetscInt, nz::PetscInt) 
Set STRUMPACK <https://portal.nersc.gov/project/sparse/strumpack/master> mesh x, y and z dimensions, for use with GEOMETRIC ordering.

Logically Collective

Input Parameters:
- `F`  - the factored matrix obtained by calling `MatGetFactor()` from PETSc-STRUMPACK interface
- `nx` - x dimension of the mesh
- `ny` - y dimension of the mesh
- `nz` - z dimension of the mesh

Level: intermediate

-seealso: `MATSOLVERSTRUMPACK`, `MatGetFactor()`

# External Links
$(_doc_external("Mat/MatSTRUMPACKSetGeometricNxyz"))
"""
function MatSTRUMPACKSetGeometricNxyz(petsclib::PetscLibType, F::AbstractPetscMat, nx::PetscInt, ny::PetscInt, nz::PetscInt) end

@for_petsc function MatSTRUMPACKSetGeometricNxyz(petsclib::$UnionPetscLib, F::AbstractPetscMat, nx::$PetscInt, ny::$PetscInt, nz::$PetscInt )

    @chk ccall(
               (:MatSTRUMPACKSetGeometricNxyz, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, $PetscInt),
               F, nx, ny, nz,
              )


	return nothing
end 

"""
	MatSTRUMPACKSetGeometricWidth(petsclib::PetscLibType,F::AbstractPetscMat, w::PetscInt) 
Set STRUMPACK <https://portal.nersc.gov/project/sparse/strumpack/master> width of the separator, for use with GEOMETRIC ordering.

Logically Collective

Input Parameters:
- `F` - the factored matrix obtained by calling `MatGetFactor()` from PETSc-STRUMPACK interface
- `w` - width of the separator

Options Database Key:
- `-mat_strumpack_geometric_width <1>` - Width of the separator of the mesh, for geometric nested dissection ordering

Level: intermediate

-seealso: `MATSOLVERSTRUMPACK`, `MatGetFactor()`

# External Links
$(_doc_external("Mat/MatSTRUMPACKSetGeometricWidth"))
"""
function MatSTRUMPACKSetGeometricWidth(petsclib::PetscLibType, F::AbstractPetscMat, w::PetscInt) end

@for_petsc function MatSTRUMPACKSetGeometricWidth(petsclib::$UnionPetscLib, F::AbstractPetscMat, w::$PetscInt )

    @chk ccall(
               (:MatSTRUMPACKSetGeometricWidth, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt),
               F, w,
              )


	return nothing
end 

"""
	MatSTRUMPACKSetReordering(petsclib::PetscLibType,F::AbstractPetscMat, reordering::MatSTRUMPACKReordering) 
Set STRUMPACK <https://portal.nersc.gov/project/sparse/strumpack/master> fill

Logically Collective

Input Parameters:
- `F`          - the factored matrix obtained by calling `MatGetFactor()` from PETSc-STRUMPACK interface
- `reordering` - the code to be used to find the fill-reducing reordering

Options Database Key:
- `-mat_strumpack_reordering <METIS>` - Sparsity reducing matrix reordering, see `MatSTRUMPACKReordering`

Level: intermediate

-seealso: `MATSOLVERSTRUMPACK`, [](ch_matrices), `Mat`, `MatSTRUMPACKReordering`, `MatGetFactor()`, `MatSTRUMPACKSetColPerm()`, `MatSTRUMPACKGetReordering()`

# External Links
$(_doc_external("Mat/MatSTRUMPACKSetReordering"))
"""
function MatSTRUMPACKSetReordering(petsclib::PetscLibType, F::AbstractPetscMat, reordering::MatSTRUMPACKReordering) end

@for_petsc function MatSTRUMPACKSetReordering(petsclib::$UnionPetscLib, F::AbstractPetscMat, reordering::MatSTRUMPACKReordering )

    @chk ccall(
               (:MatSTRUMPACKSetReordering, $petsc_library),
               PetscErrorCode,
               (CMat, MatSTRUMPACKReordering),
               F, reordering,
              )


	return nothing
end 

"""
	mb::PetscInt,nb::PetscInt = MatScaLAPACKGetBlockSizes(petsclib::PetscLibType,A::AbstractPetscMat) 
Gets the block sizes used in the distribution of
the `MATSCALAPACK` matrix

Not Collective

Input Parameter:
- `A` - a `MATSCALAPACK` matrix

Output Parameters:
- `mb` - the row block size
- `nb` - the column block size

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATSCALAPACK`, `MatCreateScaLAPACK()`, `MatScaLAPACKSetBlockSizes()`

# External Links
$(_doc_external("Mat/MatScaLAPACKGetBlockSizes"))
"""
function MatScaLAPACKGetBlockSizes(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatScaLAPACKGetBlockSizes(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	mb_ = Ref{$PetscInt}()
	nb_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatScaLAPACKGetBlockSizes, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}, Ptr{$PetscInt}),
               A, mb_, nb_,
              )

	mb = mb_[]
	nb = nb_[]

	return mb,nb
end 

"""
	MatScaLAPACKSetBlockSizes(petsclib::PetscLibType,A::AbstractPetscMat, mb::PetscInt, nb::PetscInt) 
Sets the block sizes to be used for the distribution of
the `MATSCALAPACK` matrix

Logically Collective

Input Parameters:
- `A`  - a `MATSCALAPACK` matrix
- `mb` - the row block size
- `nb` - the column block size

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATSCALAPACK`, `MatCreateScaLAPACK()`, `MatScaLAPACKGetBlockSizes()`

# External Links
$(_doc_external("Mat/MatScaLAPACKSetBlockSizes"))
"""
function MatScaLAPACKSetBlockSizes(petsclib::PetscLibType, A::AbstractPetscMat, mb::PetscInt, nb::PetscInt) end

@for_petsc function MatScaLAPACKSetBlockSizes(petsclib::$UnionPetscLib, A::AbstractPetscMat, mb::$PetscInt, nb::$PetscInt )

    @chk ccall(
               (:MatScaLAPACKSetBlockSizes, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt),
               A, mb, nb,
              )


	return nothing
end 

"""
	MatScale(petsclib::PetscLibType,mat::AbstractPetscMat, a::PetscScalar) 
Scales all elements of a matrix by a given number.

Logically Collective

Input Parameters:
- `mat` - the matrix to be scaled
- `a`   - the scaling value

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatDiagonalScale()`

# External Links
$(_doc_external("Mat/MatScale"))
"""
function MatScale(petsclib::PetscLibType, mat::AbstractPetscMat, a::PetscScalar) end

@for_petsc function MatScale(petsclib::$UnionPetscLib, mat::AbstractPetscMat, a::$PetscScalar )

    @chk ccall(
               (:MatScale, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscScalar),
               mat, a,
              )


	return nothing
end 

"""
	scatter::VecScatter = MatScatterGetVecScatter(petsclib::PetscLibType,mat::AbstractPetscMat) 
Returns the user

Logically Collective

Input Parameter:
- `mat` - the matrix, should have been created with MatCreateScatter() or have type `MATSCATTER`

Output Parameter:
- `scatter` - the scatter context

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATSCATTER`, `MatCreateScatter()`, `MatScatterSetVecScatter()`

# External Links
$(_doc_external("Mat/MatScatterGetVecScatter"))
"""
function MatScatterGetVecScatter(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatScatterGetVecScatter(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	scatter_ = Ref{VecScatter}()

    @chk ccall(
               (:MatScatterGetVecScatter, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{VecScatter}),
               mat, scatter_,
              )

	scatter = scatter_[]

	return scatter
end 

"""
	MatScatterSetVecScatter(petsclib::PetscLibType,mat::AbstractPetscMat, scatter::VecScatter) 
sets the scatter that the matrix is to apply as its linear operator in a `MATSCATTER`

Logically Collective

Input Parameters:
- `mat`     - the `MATSCATTER` matrix
- `scatter` - the scatter context create with `VecScatterCreate()`

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATSCATTER`, `MatCreateScatter()`

# External Links
$(_doc_external("Mat/MatScatterSetVecScatter"))
"""
function MatScatterSetVecScatter(petsclib::PetscLibType, mat::AbstractPetscMat, scatter::VecScatter) end

@for_petsc function MatScatterSetVecScatter(petsclib::$UnionPetscLib, mat::AbstractPetscMat, scatter::VecScatter )

    @chk ccall(
               (:MatScatterSetVecScatter, $petsc_library),
               PetscErrorCode,
               (CMat, VecScatter),
               mat, scatter,
              )


	return nothing
end 

"""
	S::PetscMat = MatSchurComplementComputeExplicitOperator(petsclib::PetscLibType,A::AbstractPetscMat) 
Compute the Schur complement matrix explicitly

Collective

Input Parameter:
- `A` - the matrix obtained with `MatCreateSchurComplement()`

Output Parameter:
- `S` - the Schur complement matrix

Level: advanced

-seealso: [](ch_ksp), `MatCreateSchurComplement()`, `MatSchurComplementUpdateSubMatrices()`, `MatSchurComplementGetPmat()`

# External Links
$(_doc_external("KSP/MatSchurComplementComputeExplicitOperator"))
"""
function MatSchurComplementComputeExplicitOperator(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatSchurComplementComputeExplicitOperator(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	S_ = Ref{CMat}()

    @chk ccall(
               (:MatSchurComplementComputeExplicitOperator, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, S_,
              )

	S = PetscMat(S_[], petsclib)

	return S
end 

"""
	ainvtype::MatSchurComplementAinvType = MatSchurComplementGetAinvType(petsclib::PetscLibType,S::AbstractPetscMat) 
get the type of approximation for the inverse of the (0,0) block used in forming `Sp` in `MatSchurComplementGetPmat()`

Not Collective

Input Parameter:
- `S` - matrix obtained with `MatCreateSchurComplement()` (or equivalent) and implementing the action of A11 - A10 ksp(A00,Ap00) A01

Output Parameter:
- `ainvtype` - type of approximation used to form approximate Schur complement Sp = A11 - A10 inv(DIAGFORM(A00)) A01:
`MAT_SCHUR_COMPLEMENT_AINV_DIAG`, `MAT_SCHUR_COMPLEMENT_AINV_LUMP`, `MAT_SCHUR_COMPLEMENT_AINV_BLOCK_DIAG`, or `MAT_SCHUR_COMPLEMENT_AINV_FULL`

Level: advanced

-seealso: [](ch_ksp), `MatSchurComplementAinvType`, `MatCreateSchurComplement()`, `MatGetSchurComplement()`, `MatSchurComplementGetPmat()`, `MatSchurComplementSetAinvType()`

# External Links
$(_doc_external("KSP/MatSchurComplementGetAinvType"))
"""
function MatSchurComplementGetAinvType(petsclib::PetscLibType, S::AbstractPetscMat) end

@for_petsc function MatSchurComplementGetAinvType(petsclib::$UnionPetscLib, S::AbstractPetscMat )
	ainvtype_ = Ref{MatSchurComplementAinvType}()

    @chk ccall(
               (:MatSchurComplementGetAinvType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatSchurComplementAinvType}),
               S, ainvtype_,
              )

	ainvtype = ainvtype_[]

	return ainvtype
end 

"""
	ksp::PetscKSP = MatSchurComplementGetKSP(petsclib::PetscLibType,S::AbstractPetscMat) 
Gets the `KSP` object that is used to solve with `A00` in the Schur complement matrix S = A11

Not Collective

Input Parameter:
- `S` - matrix obtained with `MatCreateSchurComplement()` (or equivalent) and implementing the action of  A11 - A10 ksp(A00,Ap00) A01 

Output Parameter:
- `ksp` - the linear solver object

Options Database Key:
- `-fieldsplit_<splitname_0>_XXX` - sets `KSP` and `PC` options for the 0-split solver inside the Schur complement used in `PCFIELDSPLIT`; default <splitname_0> is 0.

Level: intermediate

-seealso: [](ch_ksp), `Mat`, `MatSchurComplementSetKSP()`, `MatCreateSchurComplement()`, `MatCreateNormal()`, `MatMult()`, `MatCreate()`

# External Links
$(_doc_external("KSP/MatSchurComplementGetKSP"))
"""
function MatSchurComplementGetKSP(petsclib::PetscLibType, S::AbstractPetscMat) end

@for_petsc function MatSchurComplementGetKSP(petsclib::$UnionPetscLib, S::AbstractPetscMat )
	ksp_ = Ref{CKSP}()

    @chk ccall(
               (:MatSchurComplementGetKSP, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CKSP}),
               S, ksp_,
              )

	ksp = PetscKSP(ksp_[], petsclib)

	return ksp
end 

"""
	Sp::PetscMat = MatSchurComplementGetPmat(petsclib::PetscLibType,S::AbstractPetscMat, preuse::MatReuse) 
Obtain a matrix for preconditioning the Schur complement by assembling Sp = A11

Collective

Input Parameters:
- `S`      - matrix obtained with MatCreateSchurComplement() (or equivalent) that implements the action of A11 - A10 ksp(A00,Ap00) A01
- `preuse` - `MAT_INITIAL_MATRIX` for a new `Sp`, or `MAT_REUSE_MATRIX` to reuse an existing `Sp`, or `MAT_IGNORE_MATRIX` to put nothing in `Sp`

Output Parameter:
- `Sp` - approximate Schur complement suitable for preconditioning the exact Schur complement S = A11 - A10 inv(A00) A01

Level: advanced

-seealso: [](ch_ksp), `MatCreateSubMatrix()`, `PCFIELDSPLIT`, `MatGetSchurComplement()`, `MatCreateSchurComplement()`, `MatSchurComplementSetAinvType()`

# External Links
$(_doc_external("KSP/MatSchurComplementGetPmat"))
"""
function MatSchurComplementGetPmat(petsclib::PetscLibType, S::AbstractPetscMat, preuse::MatReuse) end

@for_petsc function MatSchurComplementGetPmat(petsclib::$UnionPetscLib, S::AbstractPetscMat, preuse::MatReuse )
	Sp_ = Ref{CMat}()

    @chk ccall(
               (:MatSchurComplementGetPmat, $petsc_library),
               PetscErrorCode,
               (CMat, MatReuse, Ptr{CMat}),
               S, preuse, Sp_,
              )

	Sp = PetscMat(Sp_[], petsclib)

	return Sp
end 

"""
	A00::PetscMat,Ap00::PetscMat,A01::PetscMat,A10::PetscMat,A11::PetscMat = MatSchurComplementGetSubMatrices(petsclib::PetscLibType,S::AbstractPetscMat) 
Get the individual submatrices in the Schur complement

Collective

Input Parameter:
- `S` - matrix obtained with `MatCreateSchurComplement()` (or equivalent) and implementing the action of A11 - A10 ksp(A00,Ap00) A01

Output Parameters:
- `A00`  - the upper-left block of the original matrix A = [A00 A01; A10 A11]
- `Ap00` - matrix from which the preconditioner is constructed for use in ksp(A00,Ap00) to approximate the action of A^{-1}
- `A01`  - the upper-right block of the original matrix A = [A00 A01; A10 A11]
- `A10`  - the lower-left block of the original matrix A = [A00 A01; A10 A11]
- `A11`  - (optional) the lower-right block of the original matrix A = [A00 A01; A10 A11]

Level: intermediate

-seealso: [](ch_ksp), `MatCreateNormal()`, `MatMult()`, `MatCreate()`, `MatSchurComplementGetKSP()`, `MatCreateSchurComplement()`, `MatSchurComplementUpdateSubMatrices()`

# External Links
$(_doc_external("KSP/MatSchurComplementGetSubMatrices"))
"""
function MatSchurComplementGetSubMatrices(petsclib::PetscLibType, S::AbstractPetscMat) end

@for_petsc function MatSchurComplementGetSubMatrices(petsclib::$UnionPetscLib, S::AbstractPetscMat )
	A00_ = Ref{CMat}()
	Ap00_ = Ref{CMat}()
	A01_ = Ref{CMat}()
	A10_ = Ref{CMat}()
	A11_ = Ref{CMat}()

    @chk ccall(
               (:MatSchurComplementGetSubMatrices, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}, Ptr{CMat}, Ptr{CMat}, Ptr{CMat}, Ptr{CMat}),
               S, A00_, Ap00_, A01_, A10_, A11_,
              )

	A00 = PetscMat(A00_[], petsclib)
	Ap00 = PetscMat(Ap00_[], petsclib)
	A01 = PetscMat(A01_[], petsclib)
	A10 = PetscMat(A10_[], petsclib)
	A11 = PetscMat(A11_[], petsclib)

	return A00,Ap00,A01,A10,A11
end 

"""
	MatSchurComplementSetAinvType(petsclib::PetscLibType,S::AbstractPetscMat, ainvtype::MatSchurComplementAinvType) 
set the type of approximation used for the inverse of the (0,0) block used in forming `Sp` in `MatSchurComplementGetPmat()`

Not Collective

Input Parameters:
- `S`        - matrix obtained with `MatCreateSchurComplement()` (or equivalent) and implementing the action of A11 - A10 ksp(A00,Ap00) A01
- `ainvtype` - type of approximation to be used to form approximate Schur complement Sp = A11 - A10 inv(DIAGFORM(A00)) A01:
`MAT_SCHUR_COMPLEMENT_AINV_DIAG`, `MAT_SCHUR_COMPLEMENT_AINV_LUMP`, `MAT_SCHUR_COMPLEMENT_AINV_BLOCK_DIAG`, or `MAT_SCHUR_COMPLEMENT_AINV_FULL`

Options Database Key:
- `-mat_schur_complement_ainv_type diag | lump | blockdiag | full` - set schur complement type

Level: advanced

-seealso: [](ch_ksp), `MatSchurComplementAinvType`, `MatCreateSchurComplement()`, `MatGetSchurComplement()`, `MatSchurComplementGetPmat()`, `MatSchurComplementGetAinvType()`

# External Links
$(_doc_external("KSP/MatSchurComplementSetAinvType"))
"""
function MatSchurComplementSetAinvType(petsclib::PetscLibType, S::AbstractPetscMat, ainvtype::MatSchurComplementAinvType) end

@for_petsc function MatSchurComplementSetAinvType(petsclib::$UnionPetscLib, S::AbstractPetscMat, ainvtype::MatSchurComplementAinvType )

    @chk ccall(
               (:MatSchurComplementSetAinvType, $petsc_library),
               PetscErrorCode,
               (CMat, MatSchurComplementAinvType),
               S, ainvtype,
              )


	return nothing
end 

"""
	MatSchurComplementSetKSP(petsclib::PetscLibType,S::AbstractPetscMat, ksp::AbstractPetscKSP) 
Sets the `KSP` object that is used to solve with `A00` in the Schur complement matrix  S = A11

Not Collective

Input Parameters:
- `S`   - matrix created with `MatCreateSchurComplement()`
- `ksp` - the linear solver object

Level: developer

-seealso: [](ch_ksp), `Mat`, `MatSchurComplementGetKSP()`, `MatCreateSchurComplement()`, `MatCreateNormal()`, `MatMult()`, `MatCreate()`, `MATSCHURCOMPLEMENT`

# External Links
$(_doc_external("KSP/MatSchurComplementSetKSP"))
"""
function MatSchurComplementSetKSP(petsclib::PetscLibType, S::AbstractPetscMat, ksp::AbstractPetscKSP) end

@for_petsc function MatSchurComplementSetKSP(petsclib::$UnionPetscLib, S::AbstractPetscMat, ksp::AbstractPetscKSP )

    @chk ccall(
               (:MatSchurComplementSetKSP, $petsc_library),
               PetscErrorCode,
               (CMat, CKSP),
               S, ksp,
              )


	return nothing
end 

"""
	MatSchurComplementSetSubMatrices(petsclib::PetscLibType,S::AbstractPetscMat, A00::AbstractPetscMat, Ap00::AbstractPetscMat, A01::AbstractPetscMat, A10::AbstractPetscMat, A11::AbstractPetscMat) 
Sets the matrices that define the Schur complement

Collective

Input Parameters:
- `S`    - matrix obtained with `MatSetType`(S,`MATSCHURCOMPLEMENT`)
- `A00`  - the upper-left block of the original matrix A = [A00 A01; A10 A11]
- `Ap00` - matrix from which the preconditioner is constructed for use in ksp(A00,Ap00) to approximate the action of A00^{-1}
- `A01`  - the upper-right block of the original matrix A = [A00 A01; A10 A11]
- `A10`  - the lower-left block of the original matrix A = [A00 A01; A10 A11]
- `A11`  - (optional) the lower-right block of the original matrix A = [A00 A01; A10 A11]

Level: intermediate

-seealso: [](ch_ksp), `Mat`, `MatCreateNormal()`, `MatMult()`, `MatCreate()`, `MatSchurComplementGetKSP()`, `MatSchurComplementUpdateSubMatrices()`, `MatCreateTranspose()`, `MatCreateSchurComplement()`, `MatGetSchurComplement()`

# External Links
$(_doc_external("KSP/MatSchurComplementSetSubMatrices"))
"""
function MatSchurComplementSetSubMatrices(petsclib::PetscLibType, S::AbstractPetscMat, A00::AbstractPetscMat, Ap00::AbstractPetscMat, A01::AbstractPetscMat, A10::AbstractPetscMat, A11::AbstractPetscMat) end

@for_petsc function MatSchurComplementSetSubMatrices(petsclib::$UnionPetscLib, S::AbstractPetscMat, A00::AbstractPetscMat, Ap00::AbstractPetscMat, A01::AbstractPetscMat, A10::AbstractPetscMat, A11::AbstractPetscMat )

    @chk ccall(
               (:MatSchurComplementSetSubMatrices, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat, CMat, CMat, CMat),
               S, A00, Ap00, A01, A10, A11,
              )


	return nothing
end 

"""
	MatSchurComplementUpdateSubMatrices(petsclib::PetscLibType,S::AbstractPetscMat, A00::AbstractPetscMat, Ap00::AbstractPetscMat, A01::AbstractPetscMat, A10::AbstractPetscMat, A11::AbstractPetscMat) 
Updates the Schur complement matrix object with new submatrices

Collective

Input Parameters:
- `S`    - matrix obtained with `MatCreateSchurComplement()` (or `MatSchurSetSubMatrices()`) and implementing the action of A11 - A10 ksp(A00,Ap00) A01
- `A00`  - the upper-left block of the original matrix A = [A00 A01; A10 A11]
- `Ap00` - matrix from which the preconditioner is constructed for use in ksp(A00,Ap00) to approximate the action of A00^{-1}
- `A01`  - the upper-right block of the original matrix A = [A00 A01; A10 A11]
- `A10`  - the lower-left block of the original matrix A = [A00 A01; A10 A11]
- `A11`  - (optional) the lower-right block of the original matrix A = [A00 A01; A10 A11]

Level: intermediate

-seealso: [](ch_ksp), `Mat`, `MatCreateNormal()`, `MatMult()`, `MatCreate()`, `MatSchurComplementGetKSP()`, `MatCreateSchurComplement()`

# External Links
$(_doc_external("KSP/MatSchurComplementUpdateSubMatrices"))
"""
function MatSchurComplementUpdateSubMatrices(petsclib::PetscLibType, S::AbstractPetscMat, A00::AbstractPetscMat, Ap00::AbstractPetscMat, A01::AbstractPetscMat, A10::AbstractPetscMat, A11::AbstractPetscMat) end

@for_petsc function MatSchurComplementUpdateSubMatrices(petsclib::$UnionPetscLib, S::AbstractPetscMat, A00::AbstractPetscMat, Ap00::AbstractPetscMat, A01::AbstractPetscMat, A10::AbstractPetscMat, A11::AbstractPetscMat )

    @chk ccall(
               (:MatSchurComplementUpdateSubMatrices, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat, CMat, CMat, CMat),
               S, A00, Ap00, A01, A10, A11,
              )


	return nothing
end 

"""
	MatSelectVariableBlockSizes(petsclib::PetscLibType,subA::AbstractPetscMat, A::AbstractPetscMat, isrow::AbstractIS) 
When creating a submatrix, pass on the variable block sizes

Not Collective

Input Parameter:
- `subA`  - the submatrix
- `A`     - the original matrix
- `isrow` - The `IS` of selected rows for the submatrix, must be sorted

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatSetVariableBlockSizes()`, `MatComputeVariableBlockEnvelope()`

# External Links
$(_doc_external("Mat/MatSelectVariableBlockSizes"))
"""
function MatSelectVariableBlockSizes(petsclib::PetscLibType, subA::AbstractPetscMat, A::AbstractPetscMat, isrow::AbstractIS) end

@for_petsc function MatSelectVariableBlockSizes(petsclib::$UnionPetscLib, subA::AbstractPetscMat, A::AbstractPetscMat, isrow::AbstractIS )

    @chk ccall(
               (:MatSelectVariableBlockSizes, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CIS),
               subA, A, isrow,
              )


	return nothing
end 

# override for MatSeqAIJGetArray; C signature: MatSeqAIJGetArray(Mat A, PetscScalar* array[])
"""
    array::Vector{PetscScalar} = MatSeqAIJGetArray(petsclib::PetscLibType, A::AbstractPetscMat)

Returns a 1D Julia Array (Vector) providing direct read/write access to the internal
numerical values of a `MATSEQAIJ` (Sparse) matrix.

Note: This only provides access to the values array. It does not allow changing
the sparsity pattern (row pointers or column indices).

Not Collective

Input Parameter:
- `A` - a `MATSEQAIJ` matrix

Output Parameter:
- `array` - A `Vector` view of the non-zero values.

Level: intermediate

See also: `MatSeqAIJRestoreArray()`, `MatDenseGetArray()`

# External Links
$(_doc_external("Mat/MatSeqAIJGetArray"))
"""
function MatSeqAIJGetArray(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatSeqAIJGetArray(petsclib::$UnionPetscLib, A::AbstractPetscMat)
    array_ = Ref{Ptr{$PetscScalar}}()

    # 1. Get the pointer to the numerical values
    @chk ccall(
               (:MatSeqAIJGetArray, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

    # 2. Get the matrix info to determine the number of nonzeros
    info_ref = Ref{LibPETSc.MatInfo}()
    LibPETSc.MatGetInfo(petsclib, A, LibPETSc.MAT_LOCAL, info_ref)

    # 3. Extract the nnz count using the [] syntax to unwrap the Ref
    nnz = Int(info_ref[].nz_used)

    return unsafe_wrap(Array, array_[], nnz; own = false)
end

"""
	array::Ptr{PetscScalar} = MatSeqAIJGetArrayRead(petsclib::PetscLibType,A::AbstractPetscMat) 
gives read

Not Collective; No Fortran Support

Input Parameter:
- `A` - a `MATSEQAIJ` matrix

Output Parameter:
- `array` - pointer to the data

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatSeqAIJGetArray()`, `MatSeqAIJRestoreArrayRead()`

# External Links
$(_doc_external("Mat/MatSeqAIJGetArrayRead"))
"""
function MatSeqAIJGetArrayRead(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatSeqAIJGetArrayRead(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatSeqAIJGetArrayRead, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

	array = array_[]

	return array
end 

"""
	array::Ptr{PetscScalar} = MatSeqAIJGetArrayWrite(petsclib::PetscLibType,A::AbstractPetscMat) 
gives write

Not Collective; No Fortran Support

Input Parameter:
- `A` - a `MATSEQAIJ` matrix

Output Parameter:
- `array` - pointer to the data

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatSeqAIJGetArray()`, `MatSeqAIJRestoreArrayRead()`

# External Links
$(_doc_external("Mat/MatSeqAIJGetArrayWrite"))
"""
function MatSeqAIJGetArrayWrite(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatSeqAIJGetArrayWrite(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatSeqAIJGetArrayWrite, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

	array = array_[]

	return array
end 

"""
	i::Vector{PetscInt},j::Vector{PetscInt},a::Vector{PetscScalar},mtype::PetscMemType = MatSeqAIJGetCSRAndMemType(petsclib::PetscLibType,mat::AbstractPetscMat) 
Get the CSR arrays and the memory type of the `MATSEQAIJ` matrix

Not Collective; No Fortran Support

Input Parameter:
- `mat` - a matrix of type `MATSEQAIJ` or its subclasses

Output Parameters:
- `i`     - row map array of the matrix
- `j`     - column index array of the matrix
- `a`     - data array of the matrix
- `mtype` - memory type of the arrays

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatSeqAIJGetArray()`, `MatSeqAIJGetArrayRead()`

# External Links
$(_doc_external("Mat/MatSeqAIJGetCSRAndMemType"))
"""
function MatSeqAIJGetCSRAndMemType(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatSeqAIJGetCSRAndMemType(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	i_ = Ref{Ptr{$PetscInt}}()
	j_ = Ref{Ptr{$PetscInt}}()
	a_ = Ref{Ptr{$PetscScalar}}()
	mtype_ = Ref{PetscMemType}()

    @chk ccall(
               (:MatSeqAIJGetCSRAndMemType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscScalar}}, Ptr{PetscMemType}),
               mat, i_, j_, a_, mtype_,
              )

	m, _ = MatGetLocalSize(petsclib, mat)
	i = unsafe_wrap(Array, i_[], Int(m) + 1; own = false)
	m, _ = MatGetLocalSize(petsclib, mat)
	nnz = Int(i[end])
	j = unsafe_wrap(Array, j_[], nnz; own = false)
	m, _ = MatGetLocalSize(petsclib, mat)
	nnz = Int(i[end])
	a = unsafe_wrap(Array, a_[], nnz; own = false)
	mtype = mtype_[]

	return i,j,a,mtype
end 

"""
	nz::PetscInt = MatSeqAIJGetMaxRowNonzeros(petsclib::PetscLibType,A::AbstractPetscMat) 
returns the maximum number of nonzeros in any row

Not Collective

Input Parameter:
- `A` - a `MATSEQAIJ` matrix

Output Parameter:
- `nz` - the maximum number of nonzeros in any row

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatSeqAIJRestoreArray()`

# External Links
$(_doc_external("Mat/MatSeqAIJGetMaxRowNonzeros"))
"""
function MatSeqAIJGetMaxRowNonzeros(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatSeqAIJGetMaxRowNonzeros(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	nz_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatSeqAIJGetMaxRowNonzeros, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}),
               A, nz_,
              )

	nz = nz_[]

	return nz
end 

"""
	C::PetscMat = MatSeqAIJKron(petsclib::PetscLibType,A::AbstractPetscMat, B::AbstractPetscMat, reuse::MatReuse) 
Computes `C`, the Kronecker product of `A` and `B`.

Input Parameters:
- `A`     - left-hand side matrix
- `B`     - right-hand side matrix
- `reuse` - either `MAT_INITIAL_MATRIX` or `MAT_REUSE_MATRIX`

Output Parameter:
- `C` - Kronecker product of `A` and `B`

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatCreateSeqAIJ()`, `MATSEQAIJ`, `MATKAIJ`, `MatReuse`

# External Links
$(_doc_external("Mat/MatSeqAIJKron"))
"""
function MatSeqAIJKron(petsclib::PetscLibType, A::AbstractPetscMat, B::AbstractPetscMat, reuse::MatReuse) end

@for_petsc function MatSeqAIJKron(petsclib::$UnionPetscLib, A::AbstractPetscMat, B::AbstractPetscMat, reuse::MatReuse )
	C_ = Ref{CMat}()

    @chk ccall(
               (:MatSeqAIJKron, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, MatReuse, Ptr{CMat}),
               A, B, reuse, C_,
              )

	C = PetscMat(C_[], petsclib)

	return C
end 

"""
	MatSeqAIJRegister(petsclib::PetscLibType,sname::String, fnc::external) 


Not Collective, No Fortran Support

Input Parameters:
- `sname`    - name of a new user-defined matrix type, for example `MATSEQAIJCRL`
- `function` - routine to convert to subtype

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatSeqAIJRegisterAll()`

# External Links
$(_doc_external("Mat/MatSeqAIJRegister"))
"""
function MatSeqAIJRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function MatSeqAIJRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:MatSeqAIJRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	MatSeqAIJRestoreArray(petsclib::PetscLibType,A::AbstractPetscMat, array::AbstractArray{PetscScalar}) 
returns access to the array where the data for a `MATSEQAIJ` matrix is stored obtained by `MatSeqAIJGetArray()`

Not Collective

Input Parameters:
- `A`     - a `MATSEQAIJ` matrix
- `array` - pointer to the data

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatSeqAIJGetArray()`

# External Links
$(_doc_external("Mat/MatSeqAIJRestoreArray"))
"""
function MatSeqAIJRestoreArray(petsclib::PetscLibType, A::AbstractPetscMat, array::AbstractArray{PetscScalar}) end

@for_petsc function MatSeqAIJRestoreArray(petsclib::$UnionPetscLib, A::AbstractPetscMat, array::AbstractArray{$PetscScalar} )
	array_ = Ref(pointer(array))

    @chk ccall(
               (:MatSeqAIJRestoreArray, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )


	return nothing
end 

"""
	array::Ptr{PetscScalar} = MatSeqAIJRestoreArrayRead(petsclib::PetscLibType,A::AbstractPetscMat) 
restore the read

Not Collective; No Fortran Support

Input Parameter:
- `A` - a `MATSEQAIJ` matrix

Output Parameter:
- `array` - pointer to the data

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatSeqAIJGetArray()`, `MatSeqAIJGetArrayRead()`

# External Links
$(_doc_external("Mat/MatSeqAIJRestoreArrayRead"))
"""
function MatSeqAIJRestoreArrayRead(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatSeqAIJRestoreArrayRead(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatSeqAIJRestoreArrayRead, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

	array = array_[]

	return array
end 

"""
	array::Ptr{PetscScalar} = MatSeqAIJRestoreArrayWrite(petsclib::PetscLibType,A::AbstractPetscMat) 
restore the read

Not Collective; No Fortran Support

Input Parameter:
- `A` - a MATSEQAIJ matrix

Output Parameter:
- `array` - pointer to the data

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatSeqAIJGetArray()`, `MatSeqAIJGetArrayRead()`

# External Links
$(_doc_external("Mat/MatSeqAIJRestoreArrayWrite"))
"""
function MatSeqAIJRestoreArrayWrite(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatSeqAIJRestoreArrayWrite(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatSeqAIJRestoreArrayWrite, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

	array = array_[]

	return array
end 

"""
	indices::PetscInt = MatSeqAIJSetColumnIndices(petsclib::PetscLibType,mat::AbstractPetscMat) 
Set the column indices for all the rows
in the matrix.

Input Parameters:
- `mat`     - the `MATSEQAIJ` matrix
- `indices` - the column indices

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATSEQAIJ`

# External Links
$(_doc_external("Mat/MatSeqAIJSetColumnIndices"))
"""
function MatSeqAIJSetColumnIndices(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatSeqAIJSetColumnIndices(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	indices_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatSeqAIJSetColumnIndices, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}),
               mat, indices_,
              )

	indices = indices_[]

	return indices
end 

"""
	MatSeqAIJSetPreallocation(petsclib::PetscLibType,B::AbstractPetscMat, nz::PetscInt, nnz::Union{Ptr, Vector{PetscInt}}) 
For good matrix assembly performance
the user should preallocate the matrix storage by setting the parameter nz
(or the array nnz).  By setting these parameters accurately, performance
during matrix assembly can be increased by more than a factor of 50.

Collective

Input Parameters:
- `B`   - The matrix
- `nz`  - number of nonzeros per row (same for all rows)
- `nnz` - array containing the number of nonzeros in the various rows
(possibly different for each row) or NULL

Options Database Keys:
- `-mat_no_inode`            - Do not use inodes
- `-mat_inode_limit <limit>` - Sets inode limit (max limit=5)

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `MatCreateAIJ()`, `MatSetValues()`, `MatSeqAIJSetColumnIndices()`, `MatCreateSeqAIJWithArrays()`, `MatGetInfo()`,
`MatSeqAIJSetTotalPreallocation()`

# External Links
$(_doc_external("Mat/MatSeqAIJSetPreallocation"))
"""
function MatSeqAIJSetPreallocation(petsclib::PetscLibType, B::AbstractPetscMat, nz::PetscInt, nnz::Union{Ptr, Vector{PetscInt}}) end

@for_petsc function MatSeqAIJSetPreallocation(petsclib::$UnionPetscLib, B::AbstractPetscMat, nz::$PetscInt, nnz::Union{Ptr, Vector{$PetscInt}} )

    @chk ccall(
               (:MatSeqAIJSetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}),
               B, nz, nnz,
              )


	return nothing
end 

"""
	MatSeqAIJSetPreallocationCSR(petsclib::PetscLibType,B::AbstractPetscMat, i::Vector{PetscInt}, j::Vector{PetscInt}, v::Vector{PetscScalar}) 
Allocates memory for a sparse sequential matrix in `MATSEQAIJ` format.

Input Parameters:
- `B` - the matrix
- `i` - the indices into `j` for the start of each row (indices start with zero)
- `j` - the column indices for each row (indices start with zero) these must be sorted for each row
- `v` - optional values in the matrix, use `NULL` if not provided

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `MatCreateSeqAIJ()`, `MatSetValues()`, `MatSeqAIJSetPreallocation()`, `MATSEQAIJ`, `MatResetPreallocation()`

# External Links
$(_doc_external("Mat/MatSeqAIJSetPreallocationCSR"))
"""
function MatSeqAIJSetPreallocationCSR(petsclib::PetscLibType, B::AbstractPetscMat, i::Vector{PetscInt}, j::Vector{PetscInt}, v::Vector{PetscScalar}) end

@for_petsc function MatSeqAIJSetPreallocationCSR(petsclib::$UnionPetscLib, B::AbstractPetscMat, i::Vector{$PetscInt}, j::Vector{$PetscInt}, v::Vector{$PetscScalar} )

    @chk ccall(
               (:MatSeqAIJSetPreallocationCSR, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}),
               B, i, j, v,
              )


	return nothing
end 

"""
	MatSeqAIJSetTotalPreallocation(petsclib::PetscLibType,A::AbstractPetscMat, nztotal::PetscInt) 
Sets an upper bound on the total number of expected nonzeros in the matrix.

Input Parameters:
- `A`       - the `MATSEQAIJ` matrix
- `nztotal` - bound on the number of nonzeros

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatSetOption()`, `MAT_SORTED_FULL`, `MatSetValues()`, `MatSeqAIJSetPreallocation()`

# External Links
$(_doc_external("Mat/MatSeqAIJSetTotalPreallocation"))
"""
function MatSeqAIJSetTotalPreallocation(petsclib::PetscLibType, A::AbstractPetscMat, nztotal::PetscInt) end

@for_petsc function MatSeqAIJSetTotalPreallocation(petsclib::$UnionPetscLib, A::AbstractPetscMat, nztotal::$PetscInt )

    @chk ccall(
               (:MatSeqAIJSetTotalPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt),
               A, nztotal,
              )


	return nothing
end 

"""
	MatSeqAIJSetType(petsclib::PetscLibType,mat::AbstractPetscMat, matype::MatType) 
Converts a `MATSEQAIJ` matrix to a subtype

Collective

Input Parameters:
- `mat`    - the matrix object
- `matype` - matrix type

Options Database Key:
- `-mat_seqaij_type  <method>` - for example seqaijcrl

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `PCSetType()`, `VecSetType()`, `MatCreate()`, `MatType`

# External Links
$(_doc_external("Mat/MatSeqAIJSetType"))
"""
function MatSeqAIJSetType(petsclib::PetscLibType, mat::AbstractPetscMat, matype::MatType) end

@for_petsc function MatSeqAIJSetType(petsclib::$UnionPetscLib, mat::AbstractPetscMat, matype::MatType )

    @chk ccall(
               (:MatSeqAIJSetType, $petsc_library),
               PetscErrorCode,
               (CMat, MatType),
               mat, matype,
              )


	return nothing
end 

"""
	MatSeqAIJSetValuesLocalFast(petsclib::PetscLibType,A::AbstractPetscMat, m::PetscInt, im::Vector{PetscInt}, n::PetscInt, in::Vector{PetscInt}, v::Vector{PetscScalar}, is::InsertMode) 

# External Links
$(_doc_external("Mat/MatSeqAIJSetValuesLocalFast"))
"""
function MatSeqAIJSetValuesLocalFast(petsclib::PetscLibType, A::AbstractPetscMat, m::PetscInt, im::Vector{PetscInt}, n::PetscInt, in::Vector{PetscInt}, v::Vector{PetscScalar}, is::InsertMode) end

@for_petsc function MatSeqAIJSetValuesLocalFast(petsclib::$UnionPetscLib, A::AbstractPetscMat, m::$PetscInt, im::Vector{$PetscInt}, n::$PetscInt, in::Vector{$PetscInt}, v::Vector{$PetscScalar}, is::InsertMode )

    @chk ccall(
               (:MatSeqAIJSetValuesLocalFast, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}, InsertMode),
               A, m, im, n, in, v, is,
              )


	return nothing
end 

"""
	array::Ptr{PetscScalar} = MatSeqBAIJGetArray(petsclib::PetscLibType,A::AbstractPetscMat) 
gives read/write access to the array where the data for a `MATSEQBAIJ` matrix is stored

Not Collective

Input Parameter:
- `A` - a `MATSEQBAIJ` matrix

Output Parameter:
- `array` - pointer to the data

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATSEQBAIJ`, `MatSeqBAIJRestoreArray()`, `MatSeqAIJGetArray()`, `MatSeqAIJRestoreArray()`

# External Links
$(_doc_external("Mat/MatSeqBAIJGetArray"))
"""
function MatSeqBAIJGetArray(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatSeqBAIJGetArray(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatSeqBAIJGetArray, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

	array = array_[]

	return array
end 

"""
	MatSeqBAIJRestoreArray(petsclib::PetscLibType,A::AbstractPetscMat, array::AbstractArray{PetscScalar}) 
returns access to the array where the data for a `MATSEQBAIJ` matrix is stored obtained by `MatSeqBAIJGetArray()`

Not Collective

Input Parameters:
- `A`     - a `MATSEQBAIJ` matrix
- `array` - pointer to the data

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatSeqBAIJGetArray()`, `MatSeqAIJGetArray()`, `MatSeqAIJRestoreArray()`

# External Links
$(_doc_external("Mat/MatSeqBAIJRestoreArray"))
"""
function MatSeqBAIJRestoreArray(petsclib::PetscLibType, A::AbstractPetscMat, array::AbstractArray{PetscScalar}) end

@for_petsc function MatSeqBAIJRestoreArray(petsclib::$UnionPetscLib, A::AbstractPetscMat, array::AbstractArray{$PetscScalar} )
	array_ = Ref(pointer(array))

    @chk ccall(
               (:MatSeqBAIJRestoreArray, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )


	return nothing
end 

"""
	indices::PetscInt = MatSeqBAIJSetColumnIndices(petsclib::PetscLibType,mat::AbstractPetscMat) 
Set the column indices for all the block rows in the matrix.

Input Parameters:
- `mat`     - the `MATSEQBAIJ` matrix
- `indices` - the block column indices

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATSEQBAIJ`, `MatSetValues()`

# External Links
$(_doc_external("Mat/MatSeqBAIJSetColumnIndices"))
"""
function MatSeqBAIJSetColumnIndices(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatSeqBAIJSetColumnIndices(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	indices_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatSeqBAIJSetColumnIndices, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}),
               mat, indices_,
              )

	indices = indices_[]

	return indices
end 

"""
	MatSeqBAIJSetPreallocation(petsclib::PetscLibType,B::AbstractPetscMat, bs::PetscInt, nz::PetscInt, nnz::Union{Ptr, Vector{PetscInt}}) 
Sets the block size and expected nonzeros
per row in the matrix. For good matrix assembly performance the
user should preallocate the matrix storage by setting the parameter `nz`
(or the array `nnz`).

Collective

Input Parameters:
- `B`   - the matrix
- `bs`  - size of block, the blocks are ALWAYS square. One can use `MatSetBlockSizes()` to set a different row and column blocksize but the row
blocksize always defines the size of the blocks. The column blocksize sets the blocksize of the vectors obtained with `MatCreateVecs()`
- `nz`  - number of block nonzeros per block row (same for all rows)
- `nnz` - array containing the number of block nonzeros in the various block rows
(possibly different for each block row) or `NULL`

Options Database Keys:
- `-mat_no_unroll`  - uses code that does not unroll the loops in the block calculations (much slower)
- `-mat_block_size` - size of the blocks to use

Level: intermediate

-seealso: [](ch_matrices), `Mat`, [Sparse Matrices](sec_matsparse), `MatCreate()`, `MatCreateSeqAIJ()`, `MatSetValues()`, `MatCreateBAIJ()`, `MatGetInfo()`

# External Links
$(_doc_external("Mat/MatSeqBAIJSetPreallocation"))
"""
function MatSeqBAIJSetPreallocation(petsclib::PetscLibType, B::AbstractPetscMat, bs::PetscInt, nz::PetscInt, nnz::Union{Ptr, Vector{PetscInt}}) end

@for_petsc function MatSeqBAIJSetPreallocation(petsclib::$UnionPetscLib, B::AbstractPetscMat, bs::$PetscInt, nz::$PetscInt, nnz::Union{Ptr, Vector{$PetscInt}} )

    @chk ccall(
               (:MatSeqBAIJSetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, Ptr{$PetscInt}),
               B, bs, nz, nnz,
              )


	return nothing
end 

"""
	MatSeqBAIJSetPreallocationCSR(petsclib::PetscLibType,B::AbstractPetscMat, bs::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, v::Union{Ptr, Vector{PetscScalar}}) 
Creates a sparse sequential matrix in `MATSEQBAIJ` format using the given nonzero structure and (optional) numerical values

Collective

Input Parameters:
- `B`  - the matrix
- `bs` - the blocksize
- `i`  - the indices into `j` for the start of each local row (indices start with zero)
- `j`  - the column indices for each local row (indices start with zero) these must be sorted for each row
- `v`  - optional values in the matrix, use `NULL` if not provided

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `MatCreateSeqBAIJ()`, `MatSetValues()`, `MatSeqBAIJSetPreallocation()`, `MATSEQBAIJ`

# External Links
$(_doc_external("Mat/MatSeqBAIJSetPreallocationCSR"))
"""
function MatSeqBAIJSetPreallocationCSR(petsclib::PetscLibType, B::AbstractPetscMat, bs::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, v::Union{Ptr, Vector{PetscScalar}}) end

@for_petsc function MatSeqBAIJSetPreallocationCSR(petsclib::$UnionPetscLib, B::AbstractPetscMat, bs::$PetscInt, i::Vector{$PetscInt}, j::Vector{$PetscInt}, v::Union{Ptr, Vector{$PetscScalar}} )

    @chk ccall(
               (:MatSeqBAIJSetPreallocationCSR, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}),
               B, bs, i, j, v,
              )


	return nothing
end 

"""
	MatSeqDenseInvert(petsclib::PetscLibType,A::AbstractPetscMat) 

# External Links
$(_doc_external("Mat/MatSeqDenseInvert"))
"""
function MatSeqDenseInvert(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatSeqDenseInvert(petsclib::$UnionPetscLib, A::AbstractPetscMat )

    @chk ccall(
               (:MatSeqDenseInvert, $petsc_library),
               PetscErrorCode,
               (CMat,),
               A,
              )


	return nothing
end 

"""
	MatSeqDenseSetPreallocation(petsclib::PetscLibType,B::AbstractPetscMat, data::Vector{PetscScalar}) 
Sets the array used for storing the matrix elements of a `MATSEQDENSE` matrix

Collective

Input Parameters:
- `B`    - the matrix
- `data` - the array (or `NULL`)

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATSEQDENSE`, `MatCreate()`, `MatCreateDense()`, `MatSetValues()`, `MatDenseSetLDA()`

# External Links
$(_doc_external("Mat/MatSeqDenseSetPreallocation"))
"""
function MatSeqDenseSetPreallocation(petsclib::PetscLibType, B::AbstractPetscMat, data::Vector{PetscScalar}) end

@for_petsc function MatSeqDenseSetPreallocation(petsclib::$UnionPetscLib, B::AbstractPetscMat, data::Vector{$PetscScalar} )

    @chk ccall(
               (:MatSeqDenseSetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscScalar}),
               B, data,
              )


	return nothing
end 

"""
	array::Ptr{PetscScalar} = MatSeqSBAIJGetArray(petsclib::PetscLibType,A::AbstractPetscMat) 
gives access to the array where the numerical data for a `MATSEQSBAIJ` matrix is stored

Not Collective

Input Parameter:
- `A` - a `MATSEQSBAIJ` matrix

Output Parameter:
- `array` - pointer to the data

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATSEQSBAIJ`, `MatSeqSBAIJRestoreArray()`, `MatSeqAIJGetArray()`, `MatSeqAIJRestoreArray()`

# External Links
$(_doc_external("Mat/MatSeqSBAIJGetArray"))
"""
function MatSeqSBAIJGetArray(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatSeqSBAIJGetArray(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatSeqSBAIJGetArray, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

	array = array_[]

	return array
end 

"""
	MatSeqSBAIJRestoreArray(petsclib::PetscLibType,A::AbstractPetscMat, array::AbstractArray{PetscScalar}) 
returns access to the array where the numerical data for a `MATSEQSBAIJ` matrix is stored obtained by `MatSeqSBAIJGetArray()`

Not Collective

Input Parameters:
- `A`     - a `MATSEQSBAIJ` matrix
- `array` - pointer to the data

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATSEQSBAIJ`, `MatSeqSBAIJGetArray()`, `MatSeqAIJGetArray()`, `MatSeqAIJRestoreArray()`

# External Links
$(_doc_external("Mat/MatSeqSBAIJRestoreArray"))
"""
function MatSeqSBAIJRestoreArray(petsclib::PetscLibType, A::AbstractPetscMat, array::AbstractArray{PetscScalar}) end

@for_petsc function MatSeqSBAIJRestoreArray(petsclib::$UnionPetscLib, A::AbstractPetscMat, array::AbstractArray{$PetscScalar} )
	array_ = Ref(pointer(array))

    @chk ccall(
               (:MatSeqSBAIJRestoreArray, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )


	return nothing
end 

"""
	indices::PetscInt = MatSeqSBAIJSetColumnIndices(petsclib::PetscLibType,mat::AbstractPetscMat) 
Set the column indices for all the rows
in a `MATSEQSBAIJ` matrix.

Input Parameters:
- `mat`     - the `MATSEQSBAIJ` matrix
- `indices` - the column indices

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATSEQSBAIJ`, `MatCreateSeqSBAIJ`

# External Links
$(_doc_external("Mat/MatSeqSBAIJSetColumnIndices"))
"""
function MatSeqSBAIJSetColumnIndices(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatSeqSBAIJSetColumnIndices(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	indices_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatSeqSBAIJSetColumnIndices, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}),
               mat, indices_,
              )

	indices = indices_[]

	return indices
end 

"""
	MatSeqSBAIJSetPreallocation(petsclib::PetscLibType,B::AbstractPetscMat, bs::PetscInt, nz::PetscInt, nnz::Vector{PetscInt}) 
Creates a sparse symmetric matrix in block AIJ (block
compressed row) `MATSEQSBAIJ` format.  For good matrix assembly performance the
user should preallocate the matrix storage by setting the parameter `nz`
(or the array `nnz`).

Collective

Input Parameters:
- `B`   - the symmetric matrix
- `bs`  - size of block, the blocks are ALWAYS square. One can use `MatSetBlockSizes()` to set a different row and column blocksize but the row
blocksize always defines the size of the blocks. The column blocksize sets the blocksize of the vectors obtained with `MatCreateVecs()`
- `nz`  - number of block nonzeros per block row (same for all rows)
- `nnz` - array containing the number of block nonzeros in the upper triangular plus
diagonal portion of each block (possibly different for each block row) or `NULL`

Options Database Keys:
- `-mat_no_unroll`  - uses code that does not unroll the loops in the block calculations (much slower)
- `-mat_block_size` - size of the blocks to use (only works if a negative bs is passed in

Level: intermediate

-seealso: [](ch_matrices), `Mat`, [Sparse Matrices](sec_matsparse), `MATSEQSBAIJ`, `MatCreate()`, `MatCreateSeqAIJ()`, `MatSetValues()`, `MatCreateSBAIJ()`

# External Links
$(_doc_external("Mat/MatSeqSBAIJSetPreallocation"))
"""
function MatSeqSBAIJSetPreallocation(petsclib::PetscLibType, B::AbstractPetscMat, bs::PetscInt, nz::PetscInt, nnz::Vector{PetscInt}) end

@for_petsc function MatSeqSBAIJSetPreallocation(petsclib::$UnionPetscLib, B::AbstractPetscMat, bs::$PetscInt, nz::$PetscInt, nnz::Vector{$PetscInt} )

    @chk ccall(
               (:MatSeqSBAIJSetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, Ptr{$PetscInt}),
               B, bs, nz, nnz,
              )


	return nothing
end 

"""
	MatSeqSBAIJSetPreallocationCSR(petsclib::PetscLibType,B::AbstractPetscMat, bs::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, v::Vector{PetscScalar}) 
Creates a sparse parallel matrix in `MATSEQSBAIJ` format using the given nonzero structure and (optional) numerical values

Input Parameters:
- `B`  - the matrix
- `bs` - size of block, the blocks are ALWAYS square.
- `i`  - the indices into `j` for the start of each local row (indices start with zero)
- `j`  - the column indices for each local row (indices start with zero) these must be sorted for each row
- `v`  - optional values in the matrix, use `NULL` if not provided

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATSEQSBAIJ`, `MatCreate()`, `MatCreateSeqSBAIJ()`, `MatSetValuesBlocked()`, `MatSeqSBAIJSetPreallocation()`

# External Links
$(_doc_external("Mat/MatSeqSBAIJSetPreallocationCSR"))
"""
function MatSeqSBAIJSetPreallocationCSR(petsclib::PetscLibType, B::AbstractPetscMat, bs::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, v::Vector{PetscScalar}) end

@for_petsc function MatSeqSBAIJSetPreallocationCSR(petsclib::$UnionPetscLib, B::AbstractPetscMat, bs::$PetscInt, i::Vector{$PetscInt}, j::Vector{$PetscInt}, v::Vector{$PetscScalar} )

    @chk ccall(
               (:MatSeqSBAIJSetPreallocationCSR, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}),
               B, bs, i, j, v,
              )


	return nothing
end 

"""
	slicewidth::PetscReal = MatSeqSELLGetAvgSliceWidth(petsclib::PetscLibType,A::AbstractPetscMat) 
returns the average slice width.

Not Collective

Input Parameter:
- `A` - a MATSEQSELL matrix

Output Parameter:
- `slicewidth` - average slice width

Level: intermediate

-seealso: `MATSEQSELL`, `MatSeqSELLGetMaxSliceWidth()`

# External Links
$(_doc_external("Mat/MatSeqSELLGetAvgSliceWidth"))
"""
function MatSeqSELLGetAvgSliceWidth(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatSeqSELLGetAvgSliceWidth(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	slicewidth_ = Ref{$PetscReal}()

    @chk ccall(
               (:MatSeqSELLGetAvgSliceWidth, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscReal}),
               A, slicewidth_,
              )

	slicewidth = slicewidth_[]

	return slicewidth
end 

"""
	ratio::PetscReal = MatSeqSELLGetFillRatio(petsclib::PetscLibType,A::AbstractPetscMat) 
returns a ratio that indicates the irregularity of the matrix.

Not Collective

Input Parameter:
- `A` - a MATSEQSELL matrix

Output Parameter:
- `ratio` - ratio of number of padded zeros to number of allocated elements

Level: intermediate

-seealso: `MATSEQSELL`, `MatSeqSELLGetAvgSliceWidth()`

# External Links
$(_doc_external("Mat/MatSeqSELLGetFillRatio"))
"""
function MatSeqSELLGetFillRatio(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatSeqSELLGetFillRatio(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	ratio_ = Ref{$PetscReal}()

    @chk ccall(
               (:MatSeqSELLGetFillRatio, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscReal}),
               A, ratio_,
              )

	ratio = ratio_[]

	return ratio
end 

"""
	slicewidth::PetscInt = MatSeqSELLGetMaxSliceWidth(petsclib::PetscLibType,A::AbstractPetscMat) 
returns the maximum slice width.

Not Collective

Input Parameter:
- `A` - a MATSEQSELL matrix

Output Parameter:
- `slicewidth` - maximum slice width

Level: intermediate

-seealso: `MATSEQSELL`, `MatSeqSELLGetAvgSliceWidth()`

# External Links
$(_doc_external("Mat/MatSeqSELLGetMaxSliceWidth"))
"""
function MatSeqSELLGetMaxSliceWidth(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatSeqSELLGetMaxSliceWidth(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	slicewidth_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatSeqSELLGetMaxSliceWidth, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}),
               A, slicewidth_,
              )

	slicewidth = slicewidth_[]

	return slicewidth
end 

"""
	variance::PetscReal = MatSeqSELLGetVarSliceSize(petsclib::PetscLibType,A::AbstractPetscMat) 
returns the variance of the slice size.

Not Collective

Input Parameter:
- `A` - a MATSEQSELL matrix

Output Parameter:
- `variance` - variance of the slice size

Level: intermediate

-seealso: `MATSEQSELL`, `MatSeqSELLSetSliceHeight()`

# External Links
$(_doc_external("Mat/MatSeqSELLGetVarSliceSize"))
"""
function MatSeqSELLGetVarSliceSize(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatSeqSELLGetVarSliceSize(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	variance_ = Ref{$PetscReal}()

    @chk ccall(
               (:MatSeqSELLGetVarSliceSize, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscReal}),
               A, variance_,
              )

	variance = variance_[]

	return variance
end 

"""
	MatSeqSELLSetPreallocation(petsclib::PetscLibType,B::AbstractPetscMat, rlenmax::PetscInt, rlen::Vector{PetscInt}) 
For good matrix assembly performance
the user should preallocate the matrix storage by setting the parameter `nz`
(or the array `nnz`).

Collective

Input Parameters:
- `B`       - The `MATSEQSELL` matrix
- `rlenmax` - number of nonzeros per row (same for all rows), ignored if `rlen` is provided
- `rlen`    - array containing the number of nonzeros in the various rows (possibly different for each row) or `NULL`

Level: intermediate

-seealso: `Mat`, `MATSEQSELL`, `MATSELL`, `MatCreate()`, `MatCreateSELL()`, `MatSetValues()`, `MatGetInfo()`

# External Links
$(_doc_external("Mat/MatSeqSELLSetPreallocation"))
"""
function MatSeqSELLSetPreallocation(petsclib::PetscLibType, B::AbstractPetscMat, rlenmax::PetscInt, rlen::Vector{PetscInt}) end

@for_petsc function MatSeqSELLSetPreallocation(petsclib::$UnionPetscLib, B::AbstractPetscMat, rlenmax::$PetscInt, rlen::Vector{$PetscInt} )

    @chk ccall(
               (:MatSeqSELLSetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}),
               B, rlenmax, rlen,
              )


	return nothing
end 

"""
	MatSeqSELLSetSliceHeight(petsclib::PetscLibType,A::AbstractPetscMat, sliceheight::PetscInt) 
sets the slice height.

Not Collective

Input Parameters:
- `A`           - a MATSEQSELL matrix
- `sliceheight` - slice height

-seealso: `MATSEQSELL`, `MatSeqSELLGetVarSliceSize()`

# External Links
$(_doc_external("Mat/MatSeqSELLSetSliceHeight"))
"""
function MatSeqSELLSetSliceHeight(petsclib::PetscLibType, A::AbstractPetscMat, sliceheight::PetscInt) end

@for_petsc function MatSeqSELLSetSliceHeight(petsclib::$UnionPetscLib, A::AbstractPetscMat, sliceheight::$PetscInt )

    @chk ccall(
               (:MatSeqSELLSetSliceHeight, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt),
               A, sliceheight,
              )


	return nothing
end 

"""
	MatSetBindingPropagates(petsclib::PetscLibType,A::AbstractPetscMat, flg::PetscBool) 
Sets whether the state of being bound to the CPU for a GPU matrix type propagates to child and some other associated objects

Input Parameters:
- `A`   - the matrix
- `flg` - flag indicating whether the boundtocpu flag should be propagated

Level: developer

-seealso: [](ch_matrices), `Mat`, `VecSetBindingPropagates()`, `MatGetBindingPropagates()`

# External Links
$(_doc_external("Mat/MatSetBindingPropagates"))
"""
function MatSetBindingPropagates(petsclib::PetscLibType, A::AbstractPetscMat, flg::PetscBool) end

@for_petsc function MatSetBindingPropagates(petsclib::$UnionPetscLib, A::AbstractPetscMat, flg::PetscBool )

    @chk ccall(
               (:MatSetBindingPropagates, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               A, flg,
              )


	return nothing
end 

"""
	MatSetBlockSize(petsclib::PetscLibType,mat::AbstractPetscMat, bs::PetscInt) 
Sets the matrix block size.

Logically Collective

Input Parameters:
- `mat` - the matrix
- `bs`  - block size

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATBAIJ`, `MATSBAIJ`, `MATAIJ`, `MatCreateSeqBAIJ()`, `MatCreateBAIJ()`, `MatGetBlockSize()`, `MatSetBlockSizes()`, `MatGetBlockSizes()`

# External Links
$(_doc_external("Mat/MatSetBlockSize"))
"""
function MatSetBlockSize(petsclib::PetscLibType, mat::AbstractPetscMat, bs::PetscInt) end

@for_petsc function MatSetBlockSize(petsclib::$UnionPetscLib, mat::AbstractPetscMat, bs::$PetscInt )

    @chk ccall(
               (:MatSetBlockSize, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt),
               mat, bs,
              )


	return nothing
end 

"""
	MatSetBlockSizes(petsclib::PetscLibType,mat::AbstractPetscMat, rbs::PetscInt, cbs::PetscInt) 
Sets the matrix block row and column sizes.

Logically Collective

Input Parameters:
- `mat` - the matrix
- `rbs` - row block size
- `cbs` - column block size

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatCreateSeqBAIJ()`, `MatCreateBAIJ()`, `MatGetBlockSize()`, `MatSetBlockSize()`, `MatGetBlockSizes()`

# External Links
$(_doc_external("Mat/MatSetBlockSizes"))
"""
function MatSetBlockSizes(petsclib::PetscLibType, mat::AbstractPetscMat, rbs::PetscInt, cbs::PetscInt) end

@for_petsc function MatSetBlockSizes(petsclib::$UnionPetscLib, mat::AbstractPetscMat, rbs::$PetscInt, cbs::$PetscInt )

    @chk ccall(
               (:MatSetBlockSizes, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt),
               mat, rbs, cbs,
              )


	return nothing
end 

"""
	MatSetBlockSizesFromMats(petsclib::PetscLibType,mat::AbstractPetscMat, fromRow::AbstractPetscMat, fromCol::AbstractPetscMat) 
Sets the matrix block row and column sizes to match a pair of matrices

Logically Collective

Input Parameters:
- `mat`     - the matrix
- `fromRow` - matrix from which to copy row block size
- `fromCol` - matrix from which to copy column block size (can be same as fromRow)

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatCreateSeqBAIJ()`, `MatCreateBAIJ()`, `MatGetBlockSize()`, `MatSetBlockSizes()`

# External Links
$(_doc_external("Mat/MatSetBlockSizesFromMats"))
"""
function MatSetBlockSizesFromMats(petsclib::PetscLibType, mat::AbstractPetscMat, fromRow::AbstractPetscMat, fromCol::AbstractPetscMat) end

@for_petsc function MatSetBlockSizesFromMats(petsclib::$UnionPetscLib, mat::AbstractPetscMat, fromRow::AbstractPetscMat, fromCol::AbstractPetscMat )

    @chk ccall(
               (:MatSetBlockSizesFromMats, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat),
               mat, fromRow, fromCol,
              )


	return nothing
end 

"""
	MatSetDM(petsclib::PetscLibType,A::AbstractPetscMat, dm::AbstractPetscDM) 
Sets the `DM` defining the data layout of the matrix

Not Collective

Input Parameters:
- `A`  - The `Mat`
- `dm` - The `DM`

Level: developer

Note:
This is rarely used in practice, rather `DMCreateMatrix()` is used to create a matrix associated with a particular `DM`

Developer Note:
Since the `Mat` class doesn't know about the `DM` class the `DM` object is associated with
the `Mat` through a `PetscObjectCompose()` operation

See also: 
=== 
`DM`, `MatGetDM()`, `DMCreateMatrix()`, `DMSetMatType()`

# External Links
$(_doc_external("DM/MatSetDM"))
"""
function MatSetDM(petsclib::PetscLibType, A::AbstractPetscMat, dm::AbstractPetscDM) end

@for_petsc function MatSetDM(petsclib::$UnionPetscLib, A::AbstractPetscMat, dm::AbstractPetscDM )

    @chk ccall(
               (:MatSetDM, $petsc_library),
               PetscErrorCode,
               (CMat, CDM),
               A, dm,
              )


	return nothing
end 

"""
	MatSetErrorIfFailure(petsclib::PetscLibType,mat::AbstractPetscMat, flg::PetscBool) 
Causes `Mat` to generate an immediate error, for example a zero pivot, is detected.

Logically Collective

Input Parameters:
- `mat` - matrix obtained from `MatCreate()`
- `flg` - `PETSC_TRUE` indicates you want the error generated

Level: advanced

-seealso: [](ch_matrices), `Mat`, `PCSetErrorIfFailure()`, `KSPConvergedReason`, `SNESConvergedReason`

# External Links
$(_doc_external("Mat/MatSetErrorIfFailure"))
"""
function MatSetErrorIfFailure(petsclib::PetscLibType, mat::AbstractPetscMat, flg::PetscBool) end

@for_petsc function MatSetErrorIfFailure(petsclib::$UnionPetscLib, mat::AbstractPetscMat, flg::PetscBool )

    @chk ccall(
               (:MatSetErrorIfFailure, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               mat, flg,
              )


	return nothing
end 

"""
	MatSetFactorType(petsclib::PetscLibType,mat::AbstractPetscMat, t::MatFactorType) 
sets the type of factorization a matrix is

Logically Collective

Input Parameters:
- `mat` - the matrix
- `t`   - the type, one of `MAT_FACTOR_NONE`, `MAT_FACTOR_LU`, `MAT_FACTOR_CHOLESKY`, `MAT_FACTOR_ILU`, `MAT_FACTOR_ICC,MAT_FACTOR_ILUDT`, `MAT_FACTOR_QR`

Level: intermediate

-seealso: [](ch_matrices), `Mat`, [Matrix Factorization](sec_matfactor), `MatFactorType`, `MatGetFactor()`, `MatGetFactorType()`, `MAT_FACTOR_NONE`, `MAT_FACTOR_LU`, `MAT_FACTOR_CHOLESKY`, `MAT_FACTOR_ILU`,
`MAT_FACTOR_ICC`,`MAT_FACTOR_ILUDT`, `MAT_FACTOR_QR`

# External Links
$(_doc_external("Mat/MatSetFactorType"))
"""
function MatSetFactorType(petsclib::PetscLibType, mat::AbstractPetscMat, t::MatFactorType) end

@for_petsc function MatSetFactorType(petsclib::$UnionPetscLib, mat::AbstractPetscMat, t::MatFactorType )

    @chk ccall(
               (:MatSetFactorType, $petsc_library),
               PetscErrorCode,
               (CMat, MatFactorType),
               mat, t,
              )


	return nothing
end 

"""
	MatSetFromOptions(petsclib::PetscLibType,B::AbstractPetscMat) 
Creates a matrix where the type is determined
from the options database.

Collective

Input Parameter:
- `B` - the matrix

Options Database Keys:
- `-mat_type seqaij`   - `MATSEQAIJ` type, uses `MatCreateSeqAIJ()`
- `-mat_type mpiaij`   - `MATMPIAIJ` type, uses `MatCreateAIJ()`
- `-mat_type seqdense` - `MATSEQDENSE` type, uses `MatCreateSeqDense()`
- `-mat_type mpidense` - `MATMPIDENSE`, uses `MatCreateDense()`
- `-mat_type seqbaij`  - `MATSEQBAIJ`, uses `MatCreateSeqBAIJ()`
- `-mat_type mpibaij`  - `MATMPIBAIJ`, uses `MatCreateBAIJ()`

See the manpages for particular formats (e.g., `MATSEQAIJ`)
for additional format-specific options.

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatCreateSeqAIJ()`, `MatCreateAIJ()`,
`MatCreateSeqDense()`, `MatCreateDense()`,
`MatCreateSeqBAIJ()`, `MatCreateBAIJ()`,
`MatCreateSeqSBAIJ()`, `MatCreateSBAIJ()`,
`MatConvert()`

# External Links
$(_doc_external("Mat/MatSetFromOptions"))
"""
function MatSetFromOptions(petsclib::PetscLibType, B::AbstractPetscMat) end

@for_petsc function MatSetFromOptions(petsclib::$UnionPetscLib, B::AbstractPetscMat )

    @chk ccall(
               (:MatSetFromOptions, $petsc_library),
               PetscErrorCode,
               (CMat,),
               B,
              )


	return nothing
end 

"""
	MatSetHPL(petsclib::PetscLibType,A::AbstractPetscMat, iseed::Cint) 
fills a `MATSEQDENSE` matrix using the HPL 2.3 random matrix generation routine

Collective

Input Parameters:
- `A`     - the matrix
- `iseed` - the random number seed

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatCreate()`

# External Links
$(_doc_external("Mat/MatSetHPL"))
"""
function MatSetHPL(petsclib::PetscLibType, A::AbstractPetscMat, iseed::Cint) end

@for_petsc function MatSetHPL(petsclib::$UnionPetscLib, A::AbstractPetscMat, iseed::Cint )

    @chk ccall(
               (:MatSetHPL, $petsc_library),
               PetscErrorCode,
               (CMat, Cint),
               A, iseed,
              )


	return nothing
end 

"""
	MatSetInf(petsclib::PetscLibType,A::AbstractPetscMat) 

# External Links
$(_doc_external("Mat/MatSetInf"))
"""
function MatSetInf(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatSetInf(petsclib::$UnionPetscLib, A::AbstractPetscMat )

    @chk ccall(
               (:MatSetInf, $petsc_library),
               PetscErrorCode,
               (CMat,),
               A,
              )


	return nothing
end 

"""
	MatSetLayouts(petsclib::PetscLibType,A::AbstractPetscMat, rmap::PetscLayout, cmap::PetscLayout) 
Sets the `PetscLayout` objects for rows and columns of a matrix

Logically Collective

Input Parameters:
- `A`    - the matrix
- `rmap` - row layout
- `cmap` - column layout

Level: advanced

-seealso: [](ch_matrices), `Mat`, `PetscLayout`, `MatCreateVecs()`, `MatGetLocalToGlobalMapping()`, `MatGetLayouts()`

# External Links
$(_doc_external("Mat/MatSetLayouts"))
"""
function MatSetLayouts(petsclib::PetscLibType, A::AbstractPetscMat, rmap::PetscLayout, cmap::PetscLayout) end

@for_petsc function MatSetLayouts(petsclib::$UnionPetscLib, A::AbstractPetscMat, rmap::PetscLayout, cmap::PetscLayout )

    @chk ccall(
               (:MatSetLayouts, $petsc_library),
               PetscErrorCode,
               (CMat, PetscLayout, PetscLayout),
               A, rmap, cmap,
              )


	return nothing
end 

"""
	MatSetLocalToGlobalMapping(petsclib::PetscLibType,x::AbstractPetscMat, rmapping::ISLocalToGlobalMapping, cmapping::ISLocalToGlobalMapping) 
Sets a local
the routine `MatSetValuesLocal()` to allow users to insert matrix entries
using a local (per-processor) numbering.

Not Collective

Input Parameters:
- `x`        - the matrix
- `rmapping` - row mapping created with `ISLocalToGlobalMappingCreate()` or `ISLocalToGlobalMappingCreateIS()`
- `cmapping` - column mapping

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `DM`, `DMCreateMatrix()`, `MatGetLocalToGlobalMapping()`, `MatAssemblyBegin()`, `MatAssemblyEnd()`, `MatSetValues()`, `MatSetValuesLocal()`, `MatGetValuesLocal()`

# External Links
$(_doc_external("Mat/MatSetLocalToGlobalMapping"))
"""
function MatSetLocalToGlobalMapping(petsclib::PetscLibType, x::AbstractPetscMat, rmapping::ISLocalToGlobalMapping, cmapping::ISLocalToGlobalMapping) end

@for_petsc function MatSetLocalToGlobalMapping(petsclib::$UnionPetscLib, x::AbstractPetscMat, rmapping::ISLocalToGlobalMapping, cmapping::ISLocalToGlobalMapping )

    @chk ccall(
               (:MatSetLocalToGlobalMapping, $petsc_library),
               PetscErrorCode,
               (CMat, ISLocalToGlobalMapping, ISLocalToGlobalMapping),
               x, rmapping, cmapping,
              )


	return nothing
end 

"""
	MatSetNearNullSpace(petsclib::PetscLibType,mat::AbstractPetscMat, nullsp::MatNullSpace) 
attaches a null space to a matrix, which is often the null space (rigid body modes) of the operator without boundary conditions
This null space will be used to provide near null space vectors to a multigrid preconditioner built from this matrix.

Logically Collective

Input Parameters:
- `mat`    - the matrix
- `nullsp` - the null space object

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatNullSpace`, `MatCreate()`, `MatNullSpaceCreate()`, `MatSetNullSpace()`, `MatNullSpaceCreateRigidBody()`, `MatGetNearNullSpace()`

# External Links
$(_doc_external("Mat/MatSetNearNullSpace"))
"""
function MatSetNearNullSpace(petsclib::PetscLibType, mat::AbstractPetscMat, nullsp::MatNullSpace) end

@for_petsc function MatSetNearNullSpace(petsclib::$UnionPetscLib, mat::AbstractPetscMat, nullsp::MatNullSpace )

    @chk ccall(
               (:MatSetNearNullSpace, $petsc_library),
               PetscErrorCode,
               (CMat, MatNullSpace),
               mat, nullsp,
              )


	return nothing
end 

"""
	MatSetNullSpace(petsclib::PetscLibType,mat::AbstractPetscMat, nullsp::MatNullSpace) 
attaches a null space to a matrix.

Logically Collective

Input Parameters:
- `mat`    - the matrix
- `nullsp` - the null space object

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `MatNullSpaceCreate()`, `MatSetNearNullSpace()`, `MatGetNullSpace()`, `MatSetTransposeNullSpace()`, `MatGetTransposeNullSpace()`, `MatNullSpaceRemove()`,
`KSPSetPCSide()`

# External Links
$(_doc_external("Mat/MatSetNullSpace"))
"""
function MatSetNullSpace(petsclib::PetscLibType, mat::AbstractPetscMat, nullsp::MatNullSpace) end

@for_petsc function MatSetNullSpace(petsclib::$UnionPetscLib, mat::AbstractPetscMat, nullsp::MatNullSpace )

    @chk ccall(
               (:MatSetNullSpace, $petsc_library),
               PetscErrorCode,
               (CMat, MatNullSpace),
               mat, nullsp,
              )


	return nothing
end 

"""
	MatSetOption(petsclib::PetscLibType,mat::AbstractPetscMat, op::MatOption, flg::PetscBool) 
Sets a parameter option for a matrix. Some options
may be specific to certain storage formats.  Some options
determine how values will be inserted (or added). Sorted,
row-oriented input will generally assemble the fastest. The default
is row-oriented.

Logically Collective for certain operations, such as `MAT_SPD`, not collective for `MAT_ROW_ORIENTED`, see `MatOption`

Input Parameters:
- `mat` - the matrix
- `op`  - the option, one of those listed below (and possibly others),
- `flg` - turn the option on (`PETSC_TRUE`) or off (`PETSC_FALSE`)

Options Describing Matrix Structure:
- `MAT_SPD`                         - symmetric positive definite
- `MAT_SYMMETRIC`                   - symmetric in terms of both structure and value
- `MAT_HERMITIAN`                   - transpose is the complex conjugation
- `MAT_STRUCTURALLY_SYMMETRIC`      - symmetric nonzero structure
- `MAT_SYMMETRY_ETERNAL`            - indicates the symmetry (or Hermitian structure) or its absence will persist through any changes to the matrix
- `MAT_STRUCTURAL_SYMMETRY_ETERNAL` - indicates the structural symmetry or its absence will persist through any changes to the matrix
- `MAT_SPD_ETERNAL`                 - indicates the value of `MAT_SPD` (true or false) will persist through any changes to the matrix

These are not really options of the matrix, they are knowledge about the structure of the matrix that users may provide so that they
do not need to be computed (usually at a high cost)

Options For Use with `MatSetValues()`:
Insert a logically dense subblock, which can be
- `MAT_ROW_ORIENTED`                - row-oriented (default)

These options reflect the data you pass in with `MatSetValues()`; it has
nothing to do with how the data is stored internally in the matrix
data structure.

When (re)assembling a matrix, we can restrict the input for
efficiency/debugging purposes.  These options include
- `MAT_NEW_NONZERO_LOCATIONS`       - additional insertions will be allowed if they generate a new nonzero (slow)
- `MAT_FORCE_DIAGONAL_ENTRIES`      - forces diagonal entries to be allocated
- `MAT_IGNORE_OFF_PROC_ENTRIES`     - drops off-processor entries
- `MAT_NEW_NONZERO_LOCATION_ERR`    - generates an error for new matrix entry
- `MAT_USE_HASH_TABLE`              - uses a hash table to speed up matrix assembly
- `MAT_NO_OFF_PROC_ENTRIES`         - you know each process will only set values for its own rows, will generate an error if
any process sets values for another process. This avoids all reductions in the MatAssembly routines and thus improves
performance for very large process counts.
- `MAT_SUBSET_OFF_PROC_ENTRIES`     - you know that the first assembly after setting this flag will set a superset
of the off-process entries required for all subsequent assemblies. This avoids a rendezvous step in the MatAssembly
functions, instead sending only neighbor messages.

Level: intermediate

-seealso: [](ch_matrices), `MatOption`, `Mat`, `MatGetOption()`

# External Links
$(_doc_external("Mat/MatSetOption"))
"""
function MatSetOption(petsclib::PetscLibType, mat::AbstractPetscMat, op::MatOption, flg::PetscBool) end

@for_petsc function MatSetOption(petsclib::$UnionPetscLib, mat::AbstractPetscMat, op::MatOption, flg::PetscBool )

    @chk ccall(
               (:MatSetOption, $petsc_library),
               PetscErrorCode,
               (CMat, MatOption, PetscBool),
               mat, op, flg,
              )


	return nothing
end 

"""
	MatSetOptionsPrefix(petsclib::PetscLibType,A::AbstractPetscMat, prefix::String) 
Sets the prefix used for searching for all
`Mat` options in the database.

Logically Collective

Input Parameters:
- `A`      - the matrix
- `prefix` - the prefix to prepend to all option names

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatSetFromOptions()`, `MatSetOptionsPrefixFactor()`

# External Links
$(_doc_external("Mat/MatSetOptionsPrefix"))
"""
function MatSetOptionsPrefix(petsclib::PetscLibType, A::AbstractPetscMat, prefix::String) end

@for_petsc function MatSetOptionsPrefix(petsclib::$UnionPetscLib, A::AbstractPetscMat, prefix::String )

    @chk ccall(
               (:MatSetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Cchar}),
               A, prefix,
              )


	return nothing
end 

"""
	MatSetOptionsPrefixFactor(petsclib::PetscLibType,A::AbstractPetscMat, prefix::String) 
Sets the prefix used for searching for all matrix factor options in the database for
for matrices created with `MatGetFactor()`

Logically Collective

Input Parameters:
- `A`      - the matrix
- `prefix` - the prefix to prepend to all option names for the factored matrix

Level: developer

-seealso: [](ch_matrices), `Mat`,   [Matrix Factorization](sec_matfactor), `MatGetFactor()`, `MatSetFromOptions()`, `MatSetOptionsPrefix()`, `MatAppendOptionsPrefixFactor()`

# External Links
$(_doc_external("Mat/MatSetOptionsPrefixFactor"))
"""
function MatSetOptionsPrefixFactor(petsclib::PetscLibType, A::AbstractPetscMat, prefix::String) end

@for_petsc function MatSetOptionsPrefixFactor(petsclib::$UnionPetscLib, A::AbstractPetscMat, prefix::String )

    @chk ccall(
               (:MatSetOptionsPrefixFactor, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Cchar}),
               A, prefix,
              )


	return nothing
end 

"""
	MatSetPreallocationCOO(petsclib::PetscLibType,A::AbstractPetscMat, ncoo::PetscCount, coo_i::Vector{PetscInt}, coo_j::Vector{PetscInt}) 
set preallocation for matrices using a coordinate format of the entries with global indices

Collective

Input Parameters:
- `A`     - matrix being preallocated
- `ncoo`  - number of entries
- `coo_i` - row indices
- `coo_j` - column indices

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatSetValuesCOO()`, `MatSeqAIJSetPreallocation()`, `MatMPIAIJSetPreallocation()`, `MatSeqBAIJSetPreallocation()`,
`MatMPIBAIJSetPreallocation()`, `MatSeqSBAIJSetPreallocation()`, `MatMPISBAIJSetPreallocation()`, `MatSetPreallocationCOOLocal()`,
`DMSetMatrixPreallocateSkip()`, `MatCreateSeqAIJFromTriple()`

# External Links
$(_doc_external("Mat/MatSetPreallocationCOO"))
"""
function MatSetPreallocationCOO(petsclib::PetscLibType, A::AbstractPetscMat, ncoo::PetscCount, coo_i::Vector{PetscInt}, coo_j::Vector{PetscInt}) end

@for_petsc function MatSetPreallocationCOO(petsclib::$UnionPetscLib, A::AbstractPetscMat, ncoo::PetscCount, coo_i::Vector{$PetscInt}, coo_j::Vector{$PetscInt} )

    @chk ccall(
               (:MatSetPreallocationCOO, $petsc_library),
               PetscErrorCode,
               (CMat, PetscCount, Ptr{$PetscInt}, Ptr{$PetscInt}),
               A, ncoo, coo_i, coo_j,
              )


	return nothing
end 

"""
	MatSetPreallocationCOOLocal(petsclib::PetscLibType,A::AbstractPetscMat, ncoo::PetscCount, coo_i::Vector{PetscInt}, coo_j::Vector{PetscInt}) 
set preallocation for matrices using a coordinate format of the entries with local indices

Collective

Input Parameters:
- `A`     - matrix being preallocated
- `ncoo`  - number of entries
- `coo_i` - row indices (local numbering; may be modified)
- `coo_j` - column indices (local numbering; may be modified)

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatSetValuesCOO()`, `MatSeqAIJSetPreallocation()`, `MatMPIAIJSetPreallocation()`, `MatSeqBAIJSetPreallocation()`,
`MatMPIBAIJSetPreallocation()`, `MatSeqSBAIJSetPreallocation()`, `MatMPISBAIJSetPreallocation()`, `MatSetPreallocationCOO()`,
`DMSetMatrixPreallocateSkip()`

# External Links
$(_doc_external("Mat/MatSetPreallocationCOOLocal"))
"""
function MatSetPreallocationCOOLocal(petsclib::PetscLibType, A::AbstractPetscMat, ncoo::PetscCount, coo_i::Vector{PetscInt}, coo_j::Vector{PetscInt}) end

@for_petsc function MatSetPreallocationCOOLocal(petsclib::$UnionPetscLib, A::AbstractPetscMat, ncoo::PetscCount, coo_i::Vector{$PetscInt}, coo_j::Vector{$PetscInt} )

    @chk ccall(
               (:MatSetPreallocationCOOLocal, $petsc_library),
               PetscErrorCode,
               (CMat, PetscCount, Ptr{$PetscInt}, Ptr{$PetscInt}),
               A, ncoo, coo_i, coo_j,
              )


	return nothing
end 

"""
	MatSetRandom(petsclib::PetscLibType,x::AbstractPetscMat, rctx::Union{Ptr, PetscRandom}) 
Sets all components of a matrix to random numbers.

Logically Collective

Input Parameters:
- `x`    - the matrix
- `rctx` - the `PetscRandom` object, formed by `PetscRandomCreate()`, or `NULL` and
it will create one internally.

Example:
-seealso: [](ch_matrices), `Mat`, `PetscRandom`, `PetscRandomCreate()`, `MatZeroEntries()`, `MatSetValues()`, `PetscRandomDestroy()`

# External Links
$(_doc_external("Mat/MatSetRandom"))
"""
function MatSetRandom(petsclib::PetscLibType, x::AbstractPetscMat, rctx::Union{Ptr, PetscRandom}) end

@for_petsc function MatSetRandom(petsclib::$UnionPetscLib, x::AbstractPetscMat, rctx::Union{Ptr, PetscRandom} )

    @chk ccall(
               (:MatSetRandom, $petsc_library),
               PetscErrorCode,
               (CMat, PetscRandom),
               x, rctx,
              )


	return nothing
end 

"""
	MatSetSizes(petsclib::PetscLibType,A::AbstractPetscMat, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt) 
Sets the local and global sizes, and checks to determine compatibility

Collective

Input Parameters:
- `A` - the matrix
- `m` - number of local rows (or `PETSC_DECIDE`)
- `n` - number of local columns (or `PETSC_DECIDE`)
- `M` - number of global rows (or `PETSC_DETERMINE`)
- `N` - number of global columns (or `PETSC_DETERMINE`)

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatGetSize()`, `PetscSplitOwnership()`, `MatGetOwnershipRange()`, `MatGetOwnershipRanges()`,
`MatGetOwnershipRangeColumn()`, `MatGetOwnershipRangesColumn()`, `PetscLayout`, `VecSetSizes()`

# External Links
$(_doc_external("Mat/MatSetSizes"))
"""
function MatSetSizes(petsclib::PetscLibType, A::AbstractPetscMat, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt) end

@for_petsc function MatSetSizes(petsclib::$UnionPetscLib, A::AbstractPetscMat, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt )

    @chk ccall(
               (:MatSetSizes, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, $PetscInt, $PetscInt),
               A, m, n, M, N,
              )


	return nothing
end 

"""
	MatSetStencil(petsclib::PetscLibType,mat::AbstractPetscMat, dim::PetscInt, dims::Vector{PetscInt}, starts::Vector{PetscInt}, dof::PetscInt) 
Sets the grid information for setting values into a matrix via
`MatSetValuesStencil()`

Not Collective

Input Parameters:
- `mat`    - the matrix
- `dim`    - dimension of the grid 1, 2, or 3
- `dims`   - number of grid points in x, y, and z direction, including ghost points on your processor
- `starts` - starting point of ghost nodes on your processor in x, y, and z direction
- `dof`    - number of degrees of freedom per node

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatStencil`, `MatSetOption()`, `MatAssemblyBegin()`, `MatAssemblyEnd()`, `MatSetValuesBlocked()`, `MatSetValuesLocal()`
`MatSetValues()`, `MatSetValuesBlockedStencil()`, `MatSetValuesStencil()`

# External Links
$(_doc_external("Mat/MatSetStencil"))
"""
function MatSetStencil(petsclib::PetscLibType, mat::AbstractPetscMat, dim::PetscInt, dims::Vector{PetscInt}, starts::Vector{PetscInt}, dof::PetscInt) end

@for_petsc function MatSetStencil(petsclib::$UnionPetscLib, mat::AbstractPetscMat, dim::$PetscInt, dims::Vector{$PetscInt}, starts::Vector{$PetscInt}, dof::$PetscInt )

    @chk ccall(
               (:MatSetStencil, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, $PetscInt),
               mat, dim, dims, starts, dof,
              )


	return nothing
end 

"""
	MatSetTransposeNullSpace(petsclib::PetscLibType,mat::AbstractPetscMat, nullsp::MatNullSpace) 
attaches the null space of a transpose of a matrix to the matrix

Logically Collective

Input Parameters:
- `mat`    - the matrix
- `nullsp` - the null space object

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatNullSpace`, `MatCreate()`, `MatNullSpaceCreate()`, `MatSetNearNullSpace()`, `MatGetNullSpace()`, `MatSetNullSpace()`, `MatGetTransposeNullSpace()`, `MatNullSpaceRemove()`, `KSPSetPCSide()`

# External Links
$(_doc_external("Mat/MatSetTransposeNullSpace"))
"""
function MatSetTransposeNullSpace(petsclib::PetscLibType, mat::AbstractPetscMat, nullsp::MatNullSpace) end

@for_petsc function MatSetTransposeNullSpace(petsclib::$UnionPetscLib, mat::AbstractPetscMat, nullsp::MatNullSpace )

    @chk ccall(
               (:MatSetTransposeNullSpace, $petsc_library),
               PetscErrorCode,
               (CMat, MatNullSpace),
               mat, nullsp,
              )


	return nothing
end 

"""
	MatSetType(petsclib::PetscLibType,mat::AbstractPetscMat, matype::MatType) 
Builds matrix object for a particular matrix type

Collective

Input Parameters:
- `mat`    - the matrix object
- `matype` - matrix type

Options Database Key:
- `-mat_type  <method>` - Sets the type; see `MatType`

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `PCSetType()`, `VecSetType()`, `MatCreate()`, `MatType`

# External Links
$(_doc_external("Mat/MatSetType"))
"""
function MatSetType(petsclib::PetscLibType, mat::AbstractPetscMat, matype::MatType) end

@for_petsc function MatSetType(petsclib::$UnionPetscLib, mat::AbstractPetscMat, matype::MatType )

    @chk ccall(
               (:MatSetType, $petsc_library),
               PetscErrorCode,
               (CMat, MatType),
               mat, matype,
              )


	return nothing
end 

"""
	MatSetUnfactored(petsclib::PetscLibType,mat::AbstractPetscMat) 
Resets a factored matrix to be treated as unfactored.

Logically Collective

Input Parameter:
- `mat` - the factored matrix to be reset

Level: developer

-seealso: [](ch_matrices), `Mat`, `PCFactorSetUseInPlace()`, `PCFactorGetUseInPlace()`

# External Links
$(_doc_external("Mat/MatSetUnfactored"))
"""
function MatSetUnfactored(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatSetUnfactored(petsclib::$UnionPetscLib, mat::AbstractPetscMat )

    @chk ccall(
               (:MatSetUnfactored, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	MatSetUp(petsclib::PetscLibType,A::AbstractPetscMat) 
Sets up the internal matrix data structures for later use by the matrix

Collective

Input Parameter:
- `A` - the matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatMult()`, `MatCreate()`, `MatDestroy()`, `MatXAIJSetPreallocation()`

# External Links
$(_doc_external("Mat/MatSetUp"))
"""
function MatSetUp(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatSetUp(petsclib::$UnionPetscLib, A::AbstractPetscMat )

    @chk ccall(
               (:MatSetUp, $petsc_library),
               PetscErrorCode,
               (CMat,),
               A,
              )


	return nothing
end 

"""
	MatSetValue(petsclib::PetscLibType,mat::AbstractPetscMat, i::PetscInt, j::PetscInt, va::PetscScalar, mode::InsertMode) 

# External Links
$(_doc_external("Mat/MatSetValue"))
"""
function MatSetValue(petsclib::PetscLibType, mat::AbstractPetscMat, i::PetscInt, j::PetscInt, va::PetscScalar, mode::InsertMode) end

@for_petsc function MatSetValue(petsclib::$UnionPetscLib, mat::AbstractPetscMat, i::$PetscInt, j::$PetscInt, va::$PetscScalar, mode::InsertMode )

    @chk ccall(
               (:MatSetValue, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, $PetscScalar, InsertMode),
               mat, i, j, va, mode,
              )


	return nothing
end 

"""
	MatSetValueLocal(petsclib::PetscLibType,mat::AbstractPetscMat, i::PetscInt, j::PetscInt, va::PetscScalar, mode::InsertMode) 

# External Links
$(_doc_external("Mat/MatSetValueLocal"))
"""
function MatSetValueLocal(petsclib::PetscLibType, mat::AbstractPetscMat, i::PetscInt, j::PetscInt, va::PetscScalar, mode::InsertMode) end

@for_petsc function MatSetValueLocal(petsclib::$UnionPetscLib, mat::AbstractPetscMat, i::$PetscInt, j::$PetscInt, va::$PetscScalar, mode::InsertMode )

    @chk ccall(
               (:MatSetValueLocal, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, $PetscScalar, InsertMode),
               mat, i, j, va, mode,
              )


	return nothing
end 

"""
	MatSetValues(petsclib::PetscLibType,mat::AbstractPetscMat, m::PetscInt, idxm::Vector{PetscInt}, n::PetscInt, idxn::Vector{PetscInt}, v::Vector{PetscScalar}, addv::InsertMode) 
Inserts or adds a block of values into a matrix.
These values may be cached, so `MatAssemblyBegin()` and `MatAssemblyEnd()`
MUST be called after all calls to `MatSetValues()` have been completed.

Not Collective

Input Parameters:
- `mat`  - the matrix
- `m`    - the number of rows
- `idxm` - the global indices of the rows
- `n`    - the number of columns
- `idxn` - the global indices of the columns
- `v`    - a one-dimensional array that contains the values implicitly stored as a two-dimensional array, by default in row-major order.
See `MAT_ROW_ORIENTED` in `MatSetOption()` for how to use column-major order.
- `addv` - either `ADD_VALUES` to add values to any existing entries, or `INSERT_VALUES` to replace existing entries with new values

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatSetOption()`, `MatAssemblyBegin()`, `MatAssemblyEnd()`, `MatSetValuesBlocked()`, `MatSetValuesLocal()`,
`InsertMode`, `INSERT_VALUES`, `ADD_VALUES`

# External Links
$(_doc_external("Mat/MatSetValues"))
"""
function MatSetValues(petsclib::PetscLibType, mat::AbstractPetscMat, m::PetscInt, idxm::Vector{PetscInt}, n::PetscInt, idxn::Vector{PetscInt}, v::Vector{PetscScalar}, addv::InsertMode) end

@for_petsc function MatSetValues(petsclib::$UnionPetscLib, mat::AbstractPetscMat, m::$PetscInt, idxm::Vector{$PetscInt}, n::$PetscInt, idxn::Vector{$PetscInt}, v::Vector{$PetscScalar}, addv::InsertMode )

    @chk ccall(
               (:MatSetValues, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}, InsertMode),
               mat, m, idxm, n, idxn, v, addv,
              )


	return nothing
end 

"""
	MatSetValuesBatch(petsclib::PetscLibType,mat::AbstractPetscMat, nb::PetscInt, bs::PetscInt, rows::Vector{PetscInt}, v::Vector{PetscScalar}) 
Adds (`ADD_VALUES`) many blocks of values into a matrix at once. The blocks must all be square and
the same size. Currently, this can only be called once and creates the given matrix.

Not Collective

Input Parameters:
- `mat`  - the matrix
- `nb`   - the number of blocks
- `bs`   - the number of rows (and columns) in each block
- `rows` - a concatenation of the rows for each block
- `v`    - a concatenation of logically two-dimensional arrays of values

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatSetOption()`, `MatAssemblyBegin()`, `MatAssemblyEnd()`, `MatSetValuesBlocked()`, `MatSetValuesLocal()`,
`InsertMode`, `INSERT_VALUES`, `ADD_VALUES`, `MatSetValues()`, `MatSetPreallocationCOO()`, `MatSetValuesCOO()`

# External Links
$(_doc_external("Mat/MatSetValuesBatch"))
"""
function MatSetValuesBatch(petsclib::PetscLibType, mat::AbstractPetscMat, nb::PetscInt, bs::PetscInt, rows::Vector{PetscInt}, v::Vector{PetscScalar}) end

@for_petsc function MatSetValuesBatch(petsclib::$UnionPetscLib, mat::AbstractPetscMat, nb::$PetscInt, bs::$PetscInt, rows::Vector{$PetscInt}, v::Vector{$PetscScalar} )

    @chk ccall(
               (:MatSetValuesBatch, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}),
               mat, nb, bs, rows, v,
              )


	return nothing
end 

"""
	MatSetValuesBlocked(petsclib::PetscLibType,mat::AbstractPetscMat, m::PetscInt, idxm::Vector{PetscInt}, n::PetscInt, idxn::Vector{PetscInt}, v::Vector{PetscScalar}, addv::InsertMode) 
Inserts or adds a block of values into a matrix.

Not Collective

Input Parameters:
- `mat`  - the matrix
- `m`    - the number of block rows
- `idxm` - the global block indices
- `n`    - the number of block columns
- `idxn` - the global block indices
- `v`    - a one-dimensional array that contains the values implicitly stored as a two-dimensional array, by default in row-major order.
See `MAT_ROW_ORIENTED` in `MatSetOption()` for how to use column-major order.
- `addv` - either `ADD_VALUES` to add values to any existing entries, or `INSERT_VALUES` replaces existing entries with new values

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatSetBlockSize()`, `MatSetOption()`, `MatAssemblyBegin()`, `MatAssemblyEnd()`, `MatSetValues()`, `MatSetValuesBlockedLocal()`

# External Links
$(_doc_external("Mat/MatSetValuesBlocked"))
"""
function MatSetValuesBlocked(petsclib::PetscLibType, mat::AbstractPetscMat, m::PetscInt, idxm::Vector{PetscInt}, n::PetscInt, idxn::Vector{PetscInt}, v::Vector{PetscScalar}, addv::InsertMode) end

@for_petsc function MatSetValuesBlocked(petsclib::$UnionPetscLib, mat::AbstractPetscMat, m::$PetscInt, idxm::Vector{$PetscInt}, n::$PetscInt, idxn::Vector{$PetscInt}, v::Vector{$PetscScalar}, addv::InsertMode )

    @chk ccall(
               (:MatSetValuesBlocked, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}, InsertMode),
               mat, m, idxm, n, idxn, v, addv,
              )


	return nothing
end 

"""
	MatSetValuesBlockedLocal(petsclib::PetscLibType,mat::AbstractPetscMat, nrow::PetscInt, irow::Vector{PetscInt}, ncol::PetscInt, icol::Vector{PetscInt}, y::Vector{PetscScalar}, addv::InsertMode) 
Inserts or adds values into certain locations of a matrix,
using a local ordering of the nodes a block at a time.

Not Collective

Input Parameters:
- `mat`  - the matrix
- `nrow` - number of rows
- `irow` - the row local indices
- `ncol` - number of columns
- `icol` - the column local indices
- `y`    - a one-dimensional array that contains the values implicitly stored as a two-dimensional array, by default in row-major order.
See `MAT_ROW_ORIENTED` in `MatSetOption()` for how to use column-major order.
- `addv` - either `ADD_VALUES` to add values to any existing entries, or `INSERT_VALUES` to replace existing entries with new values

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatSetBlockSize()`, `MatSetLocalToGlobalMapping()`, `MatAssemblyBegin()`, `MatAssemblyEnd()`,
`MatSetValuesLocal()`, `MatSetValuesBlocked()`

# External Links
$(_doc_external("Mat/MatSetValuesBlockedLocal"))
"""
function MatSetValuesBlockedLocal(petsclib::PetscLibType, mat::AbstractPetscMat, nrow::PetscInt, irow::Vector{PetscInt}, ncol::PetscInt, icol::Vector{PetscInt}, y::Vector{PetscScalar}, addv::InsertMode) end

@for_petsc function MatSetValuesBlockedLocal(petsclib::$UnionPetscLib, mat::AbstractPetscMat, nrow::$PetscInt, irow::Vector{$PetscInt}, ncol::$PetscInt, icol::Vector{$PetscInt}, y::Vector{$PetscScalar}, addv::InsertMode )

    @chk ccall(
               (:MatSetValuesBlockedLocal, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}, InsertMode),
               mat, nrow, irow, ncol, icol, y, addv,
              )


	return nothing
end 

"""
	MatSetValuesBlockedStencil(petsclib::PetscLibType,mat::AbstractPetscMat, m::PetscInt, idxm::Vector{MatStencil}, n::PetscInt, idxn::Vector{MatStencil}, v::Vector{PetscScalar}, addv::InsertMode) 
Inserts or adds a block of values into a matrix.
Using structured grid indexing

Not Collective

Input Parameters:
- `mat`  - the matrix
- `m`    - number of rows being entered
- `idxm` - grid coordinates for matrix rows being entered
- `n`    - number of columns being entered
- `idxn` - grid coordinates for matrix columns being entered
- `v`    - a one-dimensional array that contains the values implicitly stored as a two-dimensional array, by default in row-major order.
See `MAT_ROW_ORIENTED` in `MatSetOption()` for how to use column-major order.
- `addv` - either `ADD_VALUES` to add to existing entries or `INSERT_VALUES` to replace existing entries with new values

Level: beginner

-seealso: [](ch_matrices), `Mat`, `DMDA`, `MatSetOption()`, `MatAssemblyBegin()`, `MatAssemblyEnd()`, `MatSetValuesBlocked()`, `MatSetValuesLocal()`
`MatSetValues()`, `MatSetValuesStencil()`, `MatSetStencil()`, `DMCreateMatrix()`, `DMDAVecGetArray()`, `MatStencil`,
`MatSetBlockSize()`, `MatSetLocalToGlobalMapping()`

# External Links
$(_doc_external("Mat/MatSetValuesBlockedStencil"))
"""
function MatSetValuesBlockedStencil(petsclib::PetscLibType, mat::AbstractPetscMat, m::PetscInt, idxm::Vector{MatStencil}, n::PetscInt, idxn::Vector{MatStencil}, v::Vector{PetscScalar}, addv::InsertMode) end

@for_petsc function MatSetValuesBlockedStencil(petsclib::$UnionPetscLib, mat::AbstractPetscMat, m::$PetscInt, idxm::Vector{MatStencil}, n::$PetscInt, idxn::Vector{MatStencil}, v::Vector{$PetscScalar}, addv::InsertMode )

    @chk ccall(
               (:MatSetValuesBlockedStencil, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{MatStencil}, $PetscInt, Ptr{MatStencil}, Ptr{$PetscScalar}, InsertMode),
               mat, m, idxm, n, idxn, v, addv,
              )


	return nothing
end 

"""
	MatSetValuesCOO(petsclib::PetscLibType,A::AbstractPetscMat, coo_v::Vector{PetscScalar}, imode::InsertMode) 
set values at once in a matrix preallocated using `MatSetPreallocationCOO()`

Collective

Input Parameters:
- `A`     - matrix being preallocated
- `coo_v` - the matrix values (can be `NULL`)
- `imode` - the insert mode

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatSetPreallocationCOO()`, `MatSetPreallocationCOOLocal()`, `InsertMode`, `INSERT_VALUES`, `ADD_VALUES`

# External Links
$(_doc_external("Mat/MatSetValuesCOO"))
"""
function MatSetValuesCOO(petsclib::PetscLibType, A::AbstractPetscMat, coo_v::Vector{PetscScalar}, imode::InsertMode) end

@for_petsc function MatSetValuesCOO(petsclib::$UnionPetscLib, A::AbstractPetscMat, coo_v::Vector{$PetscScalar}, imode::InsertMode )

    @chk ccall(
               (:MatSetValuesCOO, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscScalar}, InsertMode),
               A, coo_v, imode,
              )


	return nothing
end 

"""
	MatSetValuesIS(petsclib::PetscLibType,mat::AbstractPetscMat, ism::AbstractIS, isn::AbstractIS, v::Vector{PetscScalar}, addv::InsertMode) 
Inserts or adds a block of values into a matrix using an `IS` to indicate the rows and columns
These values may be cached, so `MatAssemblyBegin()` and `MatAssemblyEnd()`
MUST be called after all calls to `MatSetValues()` have been completed.

Not Collective

Input Parameters:
- `mat`  - the matrix
- `ism`  - the rows to provide
- `isn`  - the columns to provide
- `v`    - a one-dimensional array that contains the values implicitly stored as a two-dimensional array, by default in row-major order.
See `MAT_ROW_ORIENTED` in `MatSetOption()` for how to use column-major order.
- `addv` - either `ADD_VALUES` to add values to any existing entries, or `INSERT_VALUES` to replace existing entries with new values

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatSetOption()`, `MatSetValues()`, `MatAssemblyBegin()`, `MatAssemblyEnd()`, `MatSetValuesBlocked()`, `MatSetValuesLocal()`,
`InsertMode`, `INSERT_VALUES`, `ADD_VALUES`

# External Links
$(_doc_external("Mat/MatSetValuesIS"))
"""
function MatSetValuesIS(petsclib::PetscLibType, mat::AbstractPetscMat, ism::AbstractIS, isn::AbstractIS, v::Vector{PetscScalar}, addv::InsertMode) end

@for_petsc function MatSetValuesIS(petsclib::$UnionPetscLib, mat::AbstractPetscMat, ism::AbstractIS, isn::AbstractIS, v::Vector{$PetscScalar}, addv::InsertMode )

    @chk ccall(
               (:MatSetValuesIS, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, CIS, Ptr{$PetscScalar}, InsertMode),
               mat, ism, isn, v, addv,
              )


	return nothing
end 

"""
	MatSetValuesLocal(petsclib::PetscLibType,mat::AbstractPetscMat, nrow::PetscInt, irow::Vector{PetscInt}, ncol::PetscInt, icol::Vector{PetscInt}, y::Vector{PetscScalar}, addv::InsertMode) 
Inserts or adds values into certain locations of a matrix,
using a local numbering of the rows and columns.

Not Collective

Input Parameters:
- `mat`  - the matrix
- `nrow` - number of rows
- `irow` - the row local indices
- `ncol` - number of columns
- `icol` - the column local indices
- `y`    - a one-dimensional array that contains the values implicitly stored as a two-dimensional array, by default in row-major order.
See `MAT_ROW_ORIENTED` in `MatSetOption()` for how to use column-major order.
- `addv` - either `ADD_VALUES` to add values to any existing entries, or `INSERT_VALUES` to replace existing entries with new values

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatAssemblyBegin()`, `MatAssemblyEnd()`, `MatSetValues()`, `MatSetLocalToGlobalMapping()`,
`MatGetValuesLocal()`

# External Links
$(_doc_external("Mat/MatSetValuesLocal"))
"""
function MatSetValuesLocal(petsclib::PetscLibType, mat::AbstractPetscMat, nrow::PetscInt, irow::Vector{PetscInt}, ncol::PetscInt, icol::Vector{PetscInt}, y::Vector{PetscScalar}, addv::InsertMode) end

@for_petsc function MatSetValuesLocal(petsclib::$UnionPetscLib, mat::AbstractPetscMat, nrow::$PetscInt, irow::Vector{$PetscInt}, ncol::$PetscInt, icol::Vector{$PetscInt}, y::Vector{$PetscScalar}, addv::InsertMode )

    @chk ccall(
               (:MatSetValuesLocal, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}, InsertMode),
               mat, nrow, irow, ncol, icol, y, addv,
              )


	return nothing
end 

"""
	MatSetValuesRow(petsclib::PetscLibType,mat::AbstractPetscMat, row::PetscInt, v::Vector{PetscScalar}) 
Inserts a row (block row for `MATBAIJ` matrices) of nonzero
values into a matrix

Not Collective

Input Parameters:
- `mat` - the matrix
- `row` - the (block) row to set
- `v`   - a logically two-dimensional (column major) array of values for  block matrices with blocksize larger than one, otherwise a one dimensional array of values

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatSetValues()`, `MatSetOption()`, `MatAssemblyBegin()`, `MatAssemblyEnd()`, `MatSetValuesBlocked()`, `MatSetValuesLocal()`,
`InsertMode`, `INSERT_VALUES`, `ADD_VALUES`

# External Links
$(_doc_external("Mat/MatSetValuesRow"))
"""
function MatSetValuesRow(petsclib::PetscLibType, mat::AbstractPetscMat, row::PetscInt, v::Vector{PetscScalar}) end

@for_petsc function MatSetValuesRow(petsclib::$UnionPetscLib, mat::AbstractPetscMat, row::$PetscInt, v::Vector{$PetscScalar} )

    @chk ccall(
               (:MatSetValuesRow, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscScalar}),
               mat, row, v,
              )


	return nothing
end 

"""
	MatSetValuesRowLocal(petsclib::PetscLibType,mat::AbstractPetscMat, row::PetscInt, v::Vector{PetscScalar}) 
Inserts a row (block row for `MATBAIJ` matrices) of nonzero
values into a matrix

Not Collective

Input Parameters:
- `mat` - the matrix
- `row` - the (block) row to set
- `v`   - a one-dimensional array that contains the values. For `MATBAIJ` they are implicitly stored as a two-dimensional array, by default in row-major order.
See `MAT_ROW_ORIENTED` in `MatSetOption()` for how to use column-major order.

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatSetOption()`, `MatAssemblyBegin()`, `MatAssemblyEnd()`, `MatSetValuesBlocked()`, `MatSetValuesLocal()`,
`InsertMode`, `INSERT_VALUES`, `ADD_VALUES`, `MatSetValues()`, `MatSetValuesRow()`, `MatSetLocalToGlobalMapping()`

# External Links
$(_doc_external("Mat/MatSetValuesRowLocal"))
"""
function MatSetValuesRowLocal(petsclib::PetscLibType, mat::AbstractPetscMat, row::PetscInt, v::Vector{PetscScalar}) end

@for_petsc function MatSetValuesRowLocal(petsclib::$UnionPetscLib, mat::AbstractPetscMat, row::$PetscInt, v::Vector{$PetscScalar} )

    @chk ccall(
               (:MatSetValuesRowLocal, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscScalar}),
               mat, row, v,
              )


	return nothing
end 

"""
	MatSetValuesStencil(petsclib::PetscLibType,mat::AbstractPetscMat, m::PetscInt, idxm::Vector{MatStencil}, n::PetscInt, idxn::Vector{MatStencil}, v::Vector{PetscScalar}, addv::InsertMode) 
Inserts or adds a block of values into a matrix.
Using structured grid indexing

Not Collective

Input Parameters:
- `mat`  - the matrix
- `m`    - number of rows being entered
- `idxm` - grid coordinates (and component number when dof > 1) for matrix rows being entered
- `n`    - number of columns being entered
- `idxn` - grid coordinates (and component number when dof > 1) for matrix columns being entered
- `v`    - a one-dimensional array that contains the values implicitly stored as a two-dimensional array, by default in row-major order.
See `MAT_ROW_ORIENTED` in `MatSetOption()` for how to use column-major order.
- `addv` - either `ADD_VALUES` to add to existing entries at that location or `INSERT_VALUES` to replace existing entries with new values

Level: beginner

-seealso: [](ch_matrices), `Mat`, `DMDA`, `MatSetOption()`, `MatAssemblyBegin()`, `MatAssemblyEnd()`, `MatSetValuesBlocked()`, `MatSetValuesLocal()`
`MatSetValues()`, `MatSetValuesBlockedStencil()`, `MatSetStencil()`, `DMCreateMatrix()`, `DMDAVecGetArray()`, `MatStencil`

# External Links
$(_doc_external("Mat/MatSetValuesStencil"))
"""
function MatSetValuesStencil(petsclib::PetscLibType, mat::AbstractPetscMat, m::PetscInt, idxm::Vector{MatStencil}, n::PetscInt, idxn::Vector{MatStencil}, v::Vector{PetscScalar}, addv::InsertMode) end

@for_petsc function MatSetValuesStencil(petsclib::$UnionPetscLib, mat::AbstractPetscMat, m::$PetscInt, idxm::Vector{MatStencil}, n::$PetscInt, idxn::Vector{MatStencil}, v::Vector{$PetscScalar}, addv::InsertMode )

    @chk ccall(
               (:MatSetValuesStencil, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{MatStencil}, $PetscInt, Ptr{MatStencil}, Ptr{$PetscScalar}, InsertMode),
               mat, m, idxm, n, idxn, v, addv,
              )


	return nothing
end 

"""
	MatSetVariableBlockSizes(petsclib::PetscLibType,mat::AbstractPetscMat, nblocks::PetscInt, bsizes::Vector{PetscInt}) 
Sets diagonal point

Not Collective

Input Parameters:
- `mat`     - the matrix
- `nblocks` - the number of blocks on this process, each block can only exist on a single process
- `bsizes`  - the block sizes

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatCreateSeqBAIJ()`, `MatCreateBAIJ()`, `MatGetBlockSize()`, `MatSetBlockSizes()`, `MatGetBlockSizes()`, `MatGetVariableBlockSizes()`,
`MatComputeVariableBlockEnvelope()`, `PCVPBJACOBI`

# External Links
$(_doc_external("Mat/MatSetVariableBlockSizes"))
"""
function MatSetVariableBlockSizes(petsclib::PetscLibType, mat::AbstractPetscMat, nblocks::PetscInt, bsizes::Vector{PetscInt}) end

@for_petsc function MatSetVariableBlockSizes(petsclib::$UnionPetscLib, mat::AbstractPetscMat, nblocks::$PetscInt, bsizes::Vector{$PetscInt} )

    @chk ccall(
               (:MatSetVariableBlockSizes, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}),
               mat, nblocks, bsizes,
              )


	return nothing
end 

"""
	MatSetVecType(petsclib::PetscLibType,mat::AbstractPetscMat, vtype::VecType) 
Set the vector type the matrix will return with `MatCreateVecs()`

Collective

Input Parameters:
- `mat`   - the matrix object
- `vtype` - vector type

Level: advanced

-seealso: [](ch_matrices), `Mat`, `VecType`, `VecSetType()`, `MatGetVecType()`

# External Links
$(_doc_external("Mat/MatSetVecType"))
"""
function MatSetVecType(petsclib::PetscLibType, mat::AbstractPetscMat, vtype::VecType) end

@for_petsc function MatSetVecType(petsclib::$UnionPetscLib, mat::AbstractPetscMat, vtype::VecType )

    @chk ccall(
               (:MatSetVecType, $petsc_library),
               PetscErrorCode,
               (CMat, VecType),
               mat, vtype,
              )


	return nothing
end 

"""
	MatShellGetContext(petsclib::PetscLibType,mat::Union{Ptr, AbstractPetscMat}, ctx::Ptr{Cvoid}) 
Returns the user

Not Collective

Input Parameter:
- `mat` - the matrix, should have been created with `MatCreateShell()`

Output Parameter:
- `ctx` - the user provided context

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATSHELL`, `MatCreateShell()`, `MatShellSetOperation()`, `MatShellSetContext()`

# External Links
$(_doc_external("Mat/MatShellGetContext"))
"""
function MatShellGetContext(petsclib::PetscLibType, mat::Union{Ptr, AbstractPetscMat}, ctx::Ptr{Cvoid}) end

@for_petsc function MatShellGetContext(petsclib::$UnionPetscLib, mat::Union{Ptr, AbstractPetscMat}, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:MatShellGetContext, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Cvoid}),
               mat, ctx,
              )


	return nothing
end 

"""
	MatShellSetContext(petsclib::PetscLibType,mat::AbstractPetscMat, ctx::Ptr{Cvoid}) 
sets the context for a `MATSHELL` shell matrix

Logically Collective

Input Parameters:
- `mat` - the `MATSHELL` shell matrix
- `ctx` - the context

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATSHELL`, `MatCreateShell()`, `MatShellGetContext()`, `MatShellGetOperation()`

# External Links
$(_doc_external("Mat/MatShellSetContext"))
"""
function MatShellSetContext(petsclib::PetscLibType, mat::AbstractPetscMat, ctx::Ptr{Cvoid}) end

@for_petsc function MatShellSetContext(petsclib::$UnionPetscLib, mat::AbstractPetscMat, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:MatShellSetContext, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Cvoid}),
               mat, ctx,
              )


	return nothing
end 

"""
	MatShellSetContextDestroy(petsclib::PetscLibType,mat::AbstractPetscMat, f::Ptr{Cvoid}) 
sets the destroy function for a `MATSHELL` shell matrix context

Logically Collective

Input Parameters:
- `mat` - the shell matrix
- `f`   - the context destroy function, see `PetscCtxDestroyFn` for calling sequence

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATSHELL`, `MatCreateShell()`, `MatShellSetContext()`,
`PetscCtxDestroyFn`

# External Links
$(_doc_external("Mat/MatShellSetContextDestroy"))
"""
function MatShellSetContextDestroy(petsclib::PetscLibType, mat::AbstractPetscMat, f::Ptr{Cvoid}) end

@for_petsc function MatShellSetContextDestroy(petsclib::$UnionPetscLib, mat::AbstractPetscMat, f::Ptr{Cvoid} )

    @chk ccall(
               (:MatShellSetContextDestroy, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Cvoid}),
               mat, f,
              )


	return nothing
end 

"""
	MatShellSetManageScalingShifts(petsclib::PetscLibType,A::AbstractPetscMat) 
Allows the user to control the scaling and shift operations of the `MATSHELL`. Must be called immediately
after `MatCreateShell()`

Logically Collective

Input Parameter:
- `A` - the `MATSHELL` shell matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATSHELL`, `MatCreateShell()`, `MatShellGetContext()`, `MatShellGetOperation()`, `MatShellSetContext()`, `MatShellSetOperation()`

# External Links
$(_doc_external("Mat/MatShellSetManageScalingShifts"))
"""
function MatShellSetManageScalingShifts(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatShellSetManageScalingShifts(petsclib::$UnionPetscLib, A::AbstractPetscMat )

    @chk ccall(
               (:MatShellSetManageScalingShifts, $petsc_library),
               PetscErrorCode,
               (CMat,),
               A,
              )


	return nothing
end 

"""
	MatShellSetMatProductOperation(petsclib::PetscLibType,A::AbstractPetscMat, ptype::MatProductType, symbolic::Union{Ptr, external}, numeric::Union{Ptr, external}, destroy::Union{Ptr, external}, Btype::MatType, Ctype::Union{Ptr, MatType}) 
Allows user to set a matrix matrix operation for a `MATSHELL` shell matrix.

Logically Collective; No Fortran Support

Input Parameters:
- `A`        - the `MATSHELL` shell matrix
- `ptype`    - the product type
- `symbolic` - the function for the symbolic phase (can be `NULL`)
- `numeric`  - the function for the numerical phase
- `destroy`  - the function for the destruction of the needed data generated during the symbolic phase (can be `NULL`)
- `Btype`    - the matrix type for the matrix to be multiplied against
- `Ctype`    - the matrix type for the result (can be `NULL`)

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATSHELL`, `MatCreateShell()`, `MatShellGetContext()`, `MatShellGetOperation()`, `MatShellSetContext()`, `MatSetOperation()`, `MatProductType`, `MatType`, `MatSetUp()`

# External Links
$(_doc_external("Mat/MatShellSetMatProductOperation"))
"""
function MatShellSetMatProductOperation(petsclib::PetscLibType, A::AbstractPetscMat, ptype::MatProductType, symbolic::Union{Ptr, external}, numeric::Union{Ptr, external}, destroy::Union{Ptr, external}, Btype::MatType, Ctype::Union{Ptr, MatType}) end

@for_petsc function MatShellSetMatProductOperation(petsclib::$UnionPetscLib, A::AbstractPetscMat, ptype::MatProductType, symbolic::Union{Ptr, external}, numeric::Union{Ptr, external}, destroy::Union{Ptr, external}, Btype::MatType, Ctype::Union{Ptr, MatType} )

    @chk ccall(
               (:MatShellSetMatProductOperation, $petsc_library),
               PetscErrorCode,
               (CMat, MatProductType, external, external, external, MatType, MatType),
               A, ptype, symbolic, numeric, destroy, Btype, Ctype,
              )


	return nothing
end 

"""
	MatShellSetVecType(petsclib::PetscLibType,mat::AbstractPetscMat, vtype::VecType) 
Sets the `VecType` of `Vec` returned by `MatCreateVecs()`

Logically Collective

Input Parameters:
- `mat`   - the `MATSHELL` shell matrix
- `vtype` - type to use for creating vectors

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATSHELL`, `MatCreateVecs()`

# External Links
$(_doc_external("Mat/MatShellSetVecType"))
"""
function MatShellSetVecType(petsclib::PetscLibType, mat::AbstractPetscMat, vtype::VecType) end

@for_petsc function MatShellSetVecType(petsclib::$UnionPetscLib, mat::AbstractPetscMat, vtype::VecType )

    @chk ccall(
               (:MatShellSetVecType, $petsc_library),
               PetscErrorCode,
               (CMat, VecType),
               mat, vtype,
              )


	return nothing
end 

"""
	flg::PetscBool = MatShellTestMult(petsclib::PetscLibType,mat::AbstractPetscMat, f::external, base::AbstractPetscVec, ctx::Ptr{Cvoid}) 
Compares the multiply routine provided to the `MATSHELL` with differencing on a given function.

Logically Collective; No Fortran Support

Input Parameters:
- `mat`  - the `MATSHELL` shell matrix
- `f`    - the function
- `base` - differences are computed around this vector, see `MatMFFDSetBase()`, for Jacobians this is the point at which the Jacobian is being evaluated
- `ctx`  - an optional context for the function

Output Parameter:
- `flg` - `PETSC_TRUE` if the multiply is likely correct

Options Database Key:
- `-mat_shell_test_mult_view` - print if any differences are detected between the products and print the difference

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATSHELL`, `MatCreateShell()`, `MatShellGetContext()`, `MatShellGetOperation()`, `MatShellTestMultTranspose()`

# External Links
$(_doc_external("Mat/MatShellTestMult"))
"""
function MatShellTestMult(petsclib::PetscLibType, mat::AbstractPetscMat, f::external, base::AbstractPetscVec, ctx::Ptr{Cvoid}) end

@for_petsc function MatShellTestMult(petsclib::$UnionPetscLib, mat::AbstractPetscMat, f::external, base::AbstractPetscVec, ctx::Ptr{Cvoid} )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatShellTestMult, $petsc_library),
               PetscErrorCode,
               (CMat, external, CVec, Ptr{Cvoid}, Ptr{PetscBool}),
               mat, f, base, ctx, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	flg::PetscBool = MatShellTestMultTranspose(petsclib::PetscLibType,mat::AbstractPetscMat, f::external, base::AbstractPetscVec, ctx::Ptr{Cvoid}) 
Compares the multiply transpose routine provided to the `MATSHELL` with differencing on a given function.

Logically Collective; No Fortran Support

Input Parameters:
- `mat`  - the `MATSHELL` shell matrix
- `f`    - the function
- `base` - differences are computed around this vector, see `MatMFFDSetBase()`, for Jacobians this is the point at which the Jacobian is being evaluated
- `ctx`  - an optional context for the function

Output Parameter:
- `flg` - `PETSC_TRUE` if the multiply is likely correct

Options Database Key:
- `-mat_shell_test_mult_view` - print if any differences are detected between the products and print the difference

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATSHELL`, `MatCreateShell()`, `MatShellGetContext()`, `MatShellGetOperation()`, `MatShellTestMult()`

# External Links
$(_doc_external("Mat/MatShellTestMultTranspose"))
"""
function MatShellTestMultTranspose(petsclib::PetscLibType, mat::AbstractPetscMat, f::external, base::AbstractPetscVec, ctx::Ptr{Cvoid}) end

@for_petsc function MatShellTestMultTranspose(petsclib::$UnionPetscLib, mat::AbstractPetscMat, f::external, base::AbstractPetscVec, ctx::Ptr{Cvoid} )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatShellTestMultTranspose, $petsc_library),
               PetscErrorCode,
               (CMat, external, CVec, Ptr{Cvoid}, Ptr{PetscBool}),
               mat, f, base, ctx, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	MatShift(petsclib::PetscLibType,Y::AbstractPetscMat, a::PetscScalar) 
Computes `Y =  Y + a I`, where `a` is a `PetscScalar`

Neighbor-wise Collective

Input Parameters:
- `Y` - the matrix
- `a` - the `PetscScalar`

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatDiagonalSet()`, `MatScale()`, `MatDiagonalScale()`

# External Links
$(_doc_external("Mat/MatShift"))
"""
function MatShift(petsclib::PetscLibType, Y::AbstractPetscMat, a::PetscScalar) end

@for_petsc function MatShift(petsclib::$UnionPetscLib, Y::AbstractPetscMat, a::$PetscScalar )

    @chk ccall(
               (:MatShift, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscScalar),
               Y, a,
              )


	return nothing
end 

"""
	MatSolve(petsclib::PetscLibType,mat::AbstractPetscMat, b::AbstractPetscVec, x::AbstractPetscVec) 
Solves A x = b, given a factored matrix.

Neighbor-wise Collective

Input Parameters:
- `mat` - the factored matrix
- `b`   - the right-hand-side vector

Output Parameter:
- `x` - the result vector

Level: developer

-seealso: [](ch_matrices), `Mat`, [Matrix Factorization](sec_matfactor), `MatGetFactor()`, `MatLUFactor()`, `MatSolveAdd()`, `MatSolveTranspose()`, `MatSolveTransposeAdd()`

# External Links
$(_doc_external("Mat/MatSolve"))
"""
function MatSolve(petsclib::PetscLibType, mat::AbstractPetscMat, b::AbstractPetscVec, x::AbstractPetscVec) end

@for_petsc function MatSolve(petsclib::$UnionPetscLib, mat::AbstractPetscMat, b::AbstractPetscVec, x::AbstractPetscVec )

    @chk ccall(
               (:MatSolve, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               mat, b, x,
              )


	return nothing
end 

"""
	MatSolveAdd(petsclib::PetscLibType,mat::AbstractPetscMat, b::AbstractPetscVec, y::AbstractPetscVec, x::AbstractPetscVec) 
Computes x = y + A^{

Neighbor-wise Collective

Input Parameters:
- `mat` - the factored matrix
- `b`   - the right-hand-side vector
- `y`   - the vector to be added to

Output Parameter:
- `x` - the result vector

Level: developer

-seealso: [](ch_matrices), `Mat`, [Matrix Factorization](sec_matfactor), `MatSolve()`, `MatGetFactor()`, `MatSolveTranspose()`, `MatSolveTransposeAdd()`

# External Links
$(_doc_external("Mat/MatSolveAdd"))
"""
function MatSolveAdd(petsclib::PetscLibType, mat::AbstractPetscMat, b::AbstractPetscVec, y::AbstractPetscVec, x::AbstractPetscVec) end

@for_petsc function MatSolveAdd(petsclib::$UnionPetscLib, mat::AbstractPetscMat, b::AbstractPetscVec, y::AbstractPetscVec, x::AbstractPetscVec )

    @chk ccall(
               (:MatSolveAdd, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec, CVec),
               mat, b, y, x,
              )


	return nothing
end 

"""
	MatSolveTranspose(petsclib::PetscLibType,mat::AbstractPetscMat, b::AbstractPetscVec, x::AbstractPetscVec) 
Solves A^T x = b, given a factored matrix.

Neighbor-wise Collective

Input Parameters:
- `mat` - the factored matrix
- `b`   - the right-hand-side vector

Output Parameter:
- `x` - the result vector

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatGetFactor()`, `KSP`, `MatSolve()`, `MatSolveAdd()`, `MatSolveTransposeAdd()`

# External Links
$(_doc_external("Mat/MatSolveTranspose"))
"""
function MatSolveTranspose(petsclib::PetscLibType, mat::AbstractPetscMat, b::AbstractPetscVec, x::AbstractPetscVec) end

@for_petsc function MatSolveTranspose(petsclib::$UnionPetscLib, mat::AbstractPetscMat, b::AbstractPetscVec, x::AbstractPetscVec )

    @chk ccall(
               (:MatSolveTranspose, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               mat, b, x,
              )


	return nothing
end 

"""
	MatSolveTransposeAdd(petsclib::PetscLibType,mat::AbstractPetscMat, b::AbstractPetscVec, y::AbstractPetscVec, x::AbstractPetscVec) 
Computes x = y + A^{
factored matrix.

Neighbor-wise Collective

Input Parameters:
- `mat` - the factored matrix
- `b`   - the right-hand-side vector
- `y`   - the vector to be added to

Output Parameter:
- `x` - the result vector

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatGetFactor()`, `MatSolve()`, `MatSolveAdd()`, `MatSolveTranspose()`

# External Links
$(_doc_external("Mat/MatSolveTransposeAdd"))
"""
function MatSolveTransposeAdd(petsclib::PetscLibType, mat::AbstractPetscMat, b::AbstractPetscVec, y::AbstractPetscVec, x::AbstractPetscVec) end

@for_petsc function MatSolveTransposeAdd(petsclib::$UnionPetscLib, mat::AbstractPetscMat, b::AbstractPetscVec, y::AbstractPetscVec, x::AbstractPetscVec )

    @chk ccall(
               (:MatSolveTransposeAdd, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec, CVec),
               mat, b, y, x,
              )


	return nothing
end 

"""
	MatSolverTypeRegister(petsclib::PetscLibType,package::MatSolverType, mtype::MatType, ftype::MatFactorType, createfactor::external) 
Registers a `MatSolverType` that works for a particular matrix type

Logically Collective, No Fortran Support

Input Parameters:
- `package`      - name of the package, for example `petsc` or `superlu`
- `mtype`        - the matrix type that works with this package
- `ftype`        - the type of factorization supported by the package
- `createfactor` - routine that will create the factored matrix ready to be used

Level: developer

-seealso: [](ch_matrices), `Mat`, [Matrix Factorization](sec_matfactor), `MatFactorGetSolverType()`, `MatCopy()`, `MatDuplicate()`, `MatGetFactorAvailable()`,
`MatGetFactor()`

# External Links
$(_doc_external("Mat/MatSolverTypeRegister"))
"""
function MatSolverTypeRegister(petsclib::PetscLibType, package::MatSolverType, mtype::MatType, ftype::MatFactorType, createfactor::external) end

@for_petsc function MatSolverTypeRegister(petsclib::$UnionPetscLib, package::MatSolverType, mtype::MatType, ftype::MatFactorType, createfactor::external )

    @chk ccall(
               (:MatSolverTypeRegister, $petsc_library),
               PetscErrorCode,
               (MatSolverType, MatType, MatFactorType, external),
               package, mtype, ftype, createfactor,
              )


	return nothing
end 

"""
	MatSolves(petsclib::PetscLibType,mat::AbstractPetscMat, b::Vecs, x::Vecs) 
Solves A x = b, given a factored matrix, for a collection of vectors

Neighbor-wise Collective

Input Parameters:
- `mat` - the factored matrix obtained with `MatGetFactor()`
- `b`   - the right-hand-side vectors

Output Parameter:
- `x` - the result vectors

Level: developer

-seealso: [](ch_matrices), `Mat`, `Vecs`, `MatSolveAdd()`, `MatSolveTranspose()`, `MatSolveTransposeAdd()`, `MatSolve()`

# External Links
$(_doc_external("Mat/MatSolves"))
"""
function MatSolves(petsclib::PetscLibType, mat::AbstractPetscMat, b::Vecs, x::Vecs) end

@for_petsc function MatSolves(petsclib::$UnionPetscLib, mat::AbstractPetscMat, b::Vecs, x::Vecs )

    @chk ccall(
               (:MatSolves, $petsc_library),
               PetscErrorCode,
               (CMat, Vecs, Vecs),
               mat, b, x,
              )


	return nothing
end 

"""
	nstash::PetscInt,reallocs::PetscInt,bnstash::PetscInt,breallocs::PetscInt = MatStashGetInfo(petsclib::PetscLibType,mat::AbstractPetscMat) 
Gets how many values are currently in the matrix stash, i.e. need
to be communicated to other processors during the `MatAssemblyBegin()`/`MatAssemblyEnd()` process

Not Collective

Input Parameter:
- `mat` - the matrix

Output Parameters:
- `nstash`    - the size of the stash
- `reallocs`  - the number of additional mallocs incurred.
- `bnstash`   - the size of the block stash
- `breallocs` - the number of additional mallocs incurred.in the block stash

Level: advanced

-seealso: [](ch_matrices), `MatAssemblyBegin()`, `MatAssemblyEnd()`, `Mat`, `MatStashSetInitialSize()`

# External Links
$(_doc_external("Mat/MatStashGetInfo"))
"""
function MatStashGetInfo(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatStashGetInfo(petsclib::$UnionPetscLib, mat::AbstractPetscMat )
	nstash_ = Ref{$PetscInt}()
	reallocs_ = Ref{$PetscInt}()
	bnstash_ = Ref{$PetscInt}()
	breallocs_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatStashGetInfo, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
               mat, nstash_, reallocs_, bnstash_, breallocs_,
              )

	nstash = nstash_[]
	reallocs = reallocs_[]
	bnstash = bnstash_[]
	breallocs = breallocs_[]

	return nstash,reallocs,bnstash,breallocs
end 

"""
	MatStashSetInitialSize(petsclib::PetscLibType,mat::AbstractPetscMat, size::PetscInt, bsize::PetscInt) 
sets the sizes of the matrix stash, that is
used during the assembly process to store values that belong to
other processors.

Not Collective

Input Parameters:
- `mat`   - the matrix
- `size`  - the initial size of the stash.
- `bsize` - the initial size of the block-stash(if used).

Options Database Keys:
- `-matstash_initial_size <size> or <size0,size1,...sizep-1>`            - set initial size
- `-matstash_block_initial_size <bsize>  or <bsize0,bsize1,...bsizep-1>` - set initial block size

Level: intermediate

-seealso: [](ch_matrices), `MatAssemblyBegin()`, `MatAssemblyEnd()`, `Mat`, `MatStashGetInfo()`

# External Links
$(_doc_external("Mat/MatStashSetInitialSize"))
"""
function MatStashSetInitialSize(petsclib::PetscLibType, mat::AbstractPetscMat, size::PetscInt, bsize::PetscInt) end

@for_petsc function MatStashSetInitialSize(petsclib::$UnionPetscLib, mat::AbstractPetscMat, size::$PetscInt, bsize::$PetscInt )

    @chk ccall(
               (:MatStashSetInitialSize, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt),
               mat, size, bsize,
              )


	return nothing
end 

"""
	MatStoreValues(petsclib::PetscLibType,mat::AbstractPetscMat) 
Stashes a copy of the matrix values; this allows reusing of the linear part of a Jacobian, while recomputing only the
nonlinear portion.

Logically Collect

Input Parameter:
- `mat` - the matrix (currently only `MATAIJ` matrices support this option)

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatRetrieveValues()`

# External Links
$(_doc_external("Mat/MatStoreValues"))
"""
function MatStoreValues(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatStoreValues(petsclib::$UnionPetscLib, mat::AbstractPetscMat )

    @chk ccall(
               (:MatStoreValues, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	MatSubMatrixVirtualUpdate(petsclib::PetscLibType,N::AbstractPetscMat, A::AbstractPetscMat, isrow::AbstractIS, iscol::AbstractIS) 
Updates a `MATSUBMATRIX` virtual submatrix

Collective

Input Parameters:
- `N`     - submatrix to update
- `A`     - full matrix in the submatrix
- `isrow` - rows in the update (same as the first time the submatrix was created)
- `iscol` - columns in the update (same as the first time the submatrix was created)

Level: developer

-seealso: [](ch_matrices), `Mat`, `MATSUBMATRIX`, `MatCreateSubMatrixVirtual()`

# External Links
$(_doc_external("Mat/MatSubMatrixVirtualUpdate"))
"""
function MatSubMatrixVirtualUpdate(petsclib::PetscLibType, N::AbstractPetscMat, A::AbstractPetscMat, isrow::AbstractIS, iscol::AbstractIS) end

@for_petsc function MatSubMatrixVirtualUpdate(petsclib::$UnionPetscLib, N::AbstractPetscMat, A::AbstractPetscMat, isrow::AbstractIS, iscol::AbstractIS )

    @chk ccall(
               (:MatSubMatrixVirtualUpdate, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CIS, CIS),
               N, A, isrow, iscol,
              )


	return nothing
end 

"""
	n::PetscInt,iss::Ptr{IS} = MatSubdomainsCreateCoalesce(petsclib::PetscLibType,A::AbstractPetscMat, N::PetscInt) 
Creates index subdomains by coalescing adjacent MPI processes' ownership ranges.

Collective

Input Parameters:
- `A` - the matrix to create subdomains from
- `N` - requested number of subdomains

Output Parameters:
- `n`   - number of subdomains resulting on this MPI process
- `iss` - `IS` list with indices of subdomains on this MPI process

Level: advanced

-seealso: [](ch_matrices), `Mat`, `IS`

# External Links
$(_doc_external("Mat/MatSubdomainsCreateCoalesce"))
"""
function MatSubdomainsCreateCoalesce(petsclib::PetscLibType, A::AbstractPetscMat, N::PetscInt) end

@for_petsc function MatSubdomainsCreateCoalesce(petsclib::$UnionPetscLib, A::AbstractPetscMat, N::$PetscInt )
	n_ = Ref{$PetscInt}()
	iss_ = Ref{Ptr{IS}}()

    @chk ccall(
               (:MatSubdomainsCreateCoalesce, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, Ptr{Ptr{CIS}}),
               A, N, n_, iss_,
              )

	n = n_[]
	iss = iss_[]

	return n,iss
end 

"""
	diagU::PetscScalar = MatSuperluDistGetDiagU(petsclib::PetscLibType,F::AbstractPetscMat) 

# External Links
$(_doc_external("Mat/MatSuperluDistGetDiagU"))
"""
function MatSuperluDistGetDiagU(petsclib::PetscLibType, F::AbstractPetscMat) end

@for_petsc function MatSuperluDistGetDiagU(petsclib::$UnionPetscLib, F::AbstractPetscMat )
	diagU_ = Ref{$PetscScalar}()

    @chk ccall(
               (:MatSuperluDistGetDiagU, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscScalar}),
               F, diagU_,
              )

	diagU = diagU_[]

	return diagU
end 

"""
	MatSuperluSetILUDropTol(petsclib::PetscLibType,F::AbstractPetscMat, dtol::PetscReal) 
Set SuperLU <https://portal.nersc.gov/project/sparse/superlu/superlu_ug.pdf> ILU drop tolerance

Logically Collective

Input Parameters:
- `F`    - the factored matrix obtained by calling `MatGetFactor()`
- `dtol` - drop tolerance

Options Database Key:
- `-mat_superlu_ilu_droptol <dtol>` - the drop tolerance

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatGetFactor()`, `MATSOLVERSUPERLU`

# External Links
$(_doc_external("Mat/MatSuperluSetILUDropTol"))
"""
function MatSuperluSetILUDropTol(petsclib::PetscLibType, F::AbstractPetscMat, dtol::PetscReal) end

@for_petsc function MatSuperluSetILUDropTol(petsclib::$UnionPetscLib, F::AbstractPetscMat, dtol::$PetscReal )

    @chk ccall(
               (:MatSuperluSetILUDropTol, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscReal),
               F, dtol,
              )


	return nothing
end 

"""
	MatTransColoringApplyDenToSp(petsclib::PetscLibType,matcoloring::MatTransposeColoring, Cden::AbstractPetscMat, Csp::AbstractPetscMat) 
Given a symbolic matrix product C_{sp} = A*B^T for which
a `MatTransposeColoring` context has been created and a dense matrix C_{den} = A*B^T_{dense}
in which `B^T_{dens}` is obtained from `MatTransColoringApplySpToDen()`, recover sparse matrix
C_{sp} from C_{den}.

Collective

Input Parameters:
- `matcoloring` - coloring context created with `MatTransposeColoringCreate()`
- `Cden`        - matrix product of a sparse matrix and a dense matrix Btdense

Output Parameter:
- `Csp` - sparse matrix

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatTransposeColoringCreate()`, `MatTransposeColoringDestroy()`, `MatTransColoringApplySpToDen()`

# External Links
$(_doc_external("Mat/MatTransColoringApplyDenToSp"))
"""
function MatTransColoringApplyDenToSp(petsclib::PetscLibType, matcoloring::MatTransposeColoring, Cden::AbstractPetscMat, Csp::AbstractPetscMat) end

@for_petsc function MatTransColoringApplyDenToSp(petsclib::$UnionPetscLib, matcoloring::MatTransposeColoring, Cden::AbstractPetscMat, Csp::AbstractPetscMat )

    @chk ccall(
               (:MatTransColoringApplyDenToSp, $petsc_library),
               PetscErrorCode,
               (MatTransposeColoring, CMat, CMat),
               matcoloring, Cden, Csp,
              )


	return nothing
end 

"""
	MatTransColoringApplySpToDen(petsclib::PetscLibType,coloring::MatTransposeColoring, B::AbstractPetscMat, Btdense::AbstractPetscMat) 
Given a symbolic matrix product C = A*B^T for which
a `MatTransposeColoring` context has been created, computes a dense B^T by applying
`MatTransposeColoring` to sparse `B`.

Collective

Input Parameters:
- `coloring` - coloring context created with `MatTransposeColoringCreate()`
- `B`        - sparse matrix

Output Parameter:
- `Btdense` - dense matrix B^T

Level: developer

-seealso: [](ch_matrices), `Mat`, `MatTransposeColoringCreate()`, `MatTransposeColoringDestroy()`, `MatTransColoringApplyDenToSp()`

# External Links
$(_doc_external("Mat/MatTransColoringApplySpToDen"))
"""
function MatTransColoringApplySpToDen(petsclib::PetscLibType, coloring::MatTransposeColoring, B::AbstractPetscMat, Btdense::AbstractPetscMat) end

@for_petsc function MatTransColoringApplySpToDen(petsclib::$UnionPetscLib, coloring::MatTransposeColoring, B::AbstractPetscMat, Btdense::AbstractPetscMat )

    @chk ccall(
               (:MatTransColoringApplySpToDen, $petsc_library),
               PetscErrorCode,
               (MatTransposeColoring, CMat, CMat),
               coloring, B, Btdense,
              )


	return nothing
end 

"""
	B::PetscMat = MatTranspose(petsclib::PetscLibType,mat::AbstractPetscMat, reuse::MatReuse) 
Computes the transpose of a matrix, either in

Collective

Input Parameters:
- `mat`   - the matrix to transpose
- `reuse` - either `MAT_INITIAL_MATRIX`, `MAT_REUSE_MATRIX`, or `MAT_INPLACE_MATRIX`

Output Parameter:
- `B` - the transpose of the matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatTransposeSetPrecursor()`, `MatMultTranspose()`, `MatMultTransposeAdd()`, `MatIsTranspose()`, `MatReuse`, `MAT_INITIAL_MATRIX`, `MAT_REUSE_MATRIX`, `MAT_INPLACE_MATRIX`,
`MatTransposeSymbolic()`, `MatCreateTranspose()`

# External Links
$(_doc_external("Mat/MatTranspose"))
"""
function MatTranspose(petsclib::PetscLibType, mat::AbstractPetscMat, reuse::MatReuse) end

@for_petsc function MatTranspose(petsclib::$UnionPetscLib, mat::AbstractPetscMat, reuse::MatReuse )
	B_ = Ref{CMat}()

    @chk ccall(
               (:MatTranspose, $petsc_library),
               PetscErrorCode,
               (CMat, MatReuse, Ptr{CMat}),
               mat, reuse, B_,
              )

	B = PetscMat(B_[], petsclib)

	return B
end 

"""
	M::PetscMat = MatTransposeGetMat(petsclib::PetscLibType,A::AbstractPetscMat) 
Gets the `Mat` object stored inside a `MATTRANSPOSEVIRTUAL`

Logically Collective

Input Parameter:
- `A` - the `MATTRANSPOSEVIRTUAL` matrix

Output Parameter:
- `M` - the matrix object stored inside `A`

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATTRANSPOSEVIRTUAL`, `MatCreateTranspose()`

# External Links
$(_doc_external("Mat/MatTransposeGetMat"))
"""
function MatTransposeGetMat(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatTransposeGetMat(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	M_ = Ref{CMat}()

    @chk ccall(
               (:MatTransposeGetMat, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, M_,
              )

	M = PetscMat(M_[], petsclib)

	return M
end 

"""
	C::PetscMat = MatTransposeMatMult(petsclib::PetscLibType,A::AbstractPetscMat, B::AbstractPetscMat, scall::MatReuse, fill::PetscReal) 
Performs matrix

Neighbor-wise Collective

Input Parameters:
- `A`     - the left matrix
- `B`     - the right matrix
- `scall` - either `MAT_INITIAL_MATRIX` or `MAT_REUSE_MATRIX`
- `fill`  - expected fill as ratio of nnz(C)/(nnz(A) + nnz(B)), use `PETSC_DETERMINE` or `PETSC_CURRENT` if not known

Output Parameter:
- `C` - the product matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatProductCreate()`, `MATPRODUCT_AtB`, `MatMatMult()`, `MatMatTransposeMult()`, `MatPtAP()`

# External Links
$(_doc_external("Mat/MatTransposeMatMult"))
"""
function MatTransposeMatMult(petsclib::PetscLibType, A::AbstractPetscMat, B::AbstractPetscMat, scall::MatReuse, fill::PetscReal) end

@for_petsc function MatTransposeMatMult(petsclib::$UnionPetscLib, A::AbstractPetscMat, B::AbstractPetscMat, scall::MatReuse, fill::$PetscReal )
	C_ = Ref{CMat}()

    @chk ccall(
               (:MatTransposeMatMult, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, MatReuse, $PetscReal, Ptr{CMat}),
               A, B, scall, fill, C_,
              )

	C = PetscMat(C_[], petsclib)

	return C
end 

"""
	flg::PetscBool = MatTransposeMatMultEqual(petsclib::PetscLibType,A::AbstractPetscMat, B::AbstractPetscMat, C::AbstractPetscMat, n::PetscInt) 
Test A^T*B*x = C*x for n random vector x

Collective

Input Parameters:
- `A` - the first matrix
- `B` - the second matrix
- `C` - the third matrix
- `n` - number of random vectors to be tested

Output Parameter:
- `flg` - `PETSC_TRUE` if the products are equal; `PETSC_FALSE` otherwise.

Level: intermediate

-seealso: `Mat`, `MatMatMultEqual()`, `MatMultEqual()`, `MatMultAddEqual()`, `MatMultTransposeEqual()`

# External Links
$(_doc_external("Mat/MatTransposeMatMultEqual"))
"""
function MatTransposeMatMultEqual(petsclib::PetscLibType, A::AbstractPetscMat, B::AbstractPetscMat, C::AbstractPetscMat, n::PetscInt) end

@for_petsc function MatTransposeMatMultEqual(petsclib::$UnionPetscLib, A::AbstractPetscMat, B::AbstractPetscMat, C::AbstractPetscMat, n::$PetscInt )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatTransposeMatMultEqual, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat, $PetscInt, Ptr{PetscBool}),
               A, B, C, n, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	MatTransposeSetPrecursor(petsclib::PetscLibType,mat::AbstractPetscMat, B::AbstractPetscMat) 
Set the matrix from which the second matrix will receive numerical transpose data with a call to `MatTranspose`(A,`MAT_REUSE_MATRIX`,&B)
when B was not obtained with `MatTranspose`(A,`MAT_INITIAL_MATRIX`,&B)

Collective

Input Parameter:
- `mat` - the matrix to provide the transpose

Output Parameter:
- `B` - the matrix to contain the transpose; it MUST have the nonzero structure of the transpose of A or the code will crash or generate incorrect results

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatTransposeSymbolic()`, `MatTranspose()`, `MatMultTranspose()`, `MatMultTransposeAdd()`, `MatIsTranspose()`, `MatReuse`, `MAT_INITIAL_MATRIX`, `MAT_REUSE_MATRIX`, `MAT_INPLACE_MATRIX`

# External Links
$(_doc_external("Mat/MatTransposeSetPrecursor"))
"""
function MatTransposeSetPrecursor(petsclib::PetscLibType, mat::AbstractPetscMat, B::AbstractPetscMat) end

@for_petsc function MatTransposeSetPrecursor(petsclib::$UnionPetscLib, mat::AbstractPetscMat, B::AbstractPetscMat )

    @chk ccall(
               (:MatTransposeSetPrecursor, $petsc_library),
               PetscErrorCode,
               (CMat, CMat),
               mat, B,
              )


	return nothing
end 

"""
	B::PetscMat = MatTransposeSymbolic(petsclib::PetscLibType,A::AbstractPetscMat) 
Computes the symbolic part of the transpose of a matrix.

Collective

Input Parameter:
- `A` - the matrix to transpose

Output Parameter:
- `B` - the transpose. This is a complete matrix but the numerical portion is invalid. One can call `MatTranspose`(A,`MAT_REUSE_MATRIX`,&B) to compute the
numerical portion.

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatTransposeSetPrecursor()`, `MatTranspose()`, `MatMultTranspose()`, `MatMultTransposeAdd()`, `MatIsTranspose()`, `MatReuse`, `MAT_INITIAL_MATRIX`, `MAT_REUSE_MATRIX`, `MAT_INPLACE_MATRIX`

# External Links
$(_doc_external("Mat/MatTransposeSymbolic"))
"""
function MatTransposeSymbolic(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatTransposeSymbolic(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	B_ = Ref{CMat}()

    @chk ccall(
               (:MatTransposeSymbolic, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, B_,
              )

	B = PetscMat(B_[], petsclib)

	return B
end 

"""
	MatUpdateMPIAIJWithArray(petsclib::PetscLibType,mat::AbstractPetscMat, v::Vector{PetscScalar}) 
updates an `MATMPIAIJ` matrix using an array that contains the nonzero values

Collective

Input Parameters:
- `mat` - the matrix
- `v`   - matrix values, stored by row

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatCreate()`, `MatCreateSeqAIJ()`, `MatSetValues()`, `MatMPIAIJSetPreallocation()`, `MatMPIAIJSetPreallocationCSR()`,
`MATMPIAIJ`, `MatCreateAIJ()`, `MatCreateMPIAIJWithSplitArrays()`, `MatUpdateMPIAIJWithArrays()`, `MatSetPreallocationCOO()`, `MatSetValuesCOO()`

# External Links
$(_doc_external("Mat/MatUpdateMPIAIJWithArray"))
"""
function MatUpdateMPIAIJWithArray(petsclib::PetscLibType, mat::AbstractPetscMat, v::Vector{PetscScalar}) end

@for_petsc function MatUpdateMPIAIJWithArray(petsclib::$UnionPetscLib, mat::AbstractPetscMat, v::Vector{$PetscScalar} )

    @chk ccall(
               (:MatUpdateMPIAIJWithArray, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscScalar}),
               mat, v,
              )


	return nothing
end 

"""
	MatUpdateMPIAIJWithArrays(petsclib::PetscLibType,mat::AbstractPetscMat, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, Ii::Vector{PetscInt}, J::Vector{PetscInt}, v::Vector{PetscScalar}) 
updates a `MATMPIAIJ` matrix using arrays that contain in standard
CSR format for the local rows. Only the numerical values are updated the other arrays must be identical to what was passed
from `MatCreateMPIAIJWithArrays()`

Deprecated: Use `MatUpdateMPIAIJWithArray()`

Collective

Input Parameters:
- `mat` - the matrix
- `m`   - number of local rows (Cannot be `PETSC_DECIDE`)
- `n`   - This value should be the same as the local size used in creating the
x vector for the matrix-vector product y = Ax. (or `PETSC_DECIDE` to have
calculated if N is given) For square matrices n is almost always m.
- `M`   - number of global rows (or `PETSC_DETERMINE` to have calculated if m is given)
- `N`   - number of global columns (or `PETSC_DETERMINE` to have calculated if n is given)
- `Ii`  - row indices; that is Ii[0] = 0, Ii[row] = Ii[row-1] + number of elements in that row of the matrix
- `J`   - column indices
- `v`   - matrix values

Level: deprecated

-seealso: [](ch_matrices), `Mat`, `MATMPIAIJ`, `MatCreate()`, `MatCreateSeqAIJ()`, `MatSetValues()`, `MatMPIAIJSetPreallocation()`, `MatMPIAIJSetPreallocationCSR()`,
`MatCreateAIJ()`, `MatCreateMPIAIJWithSplitArrays()`, `MatUpdateMPIAIJWithArray()`, `MatSetPreallocationCOO()`, `MatSetValuesCOO()`

# External Links
$(_doc_external("Mat/MatUpdateMPIAIJWithArrays"))
"""
function MatUpdateMPIAIJWithArrays(petsclib::PetscLibType, mat::AbstractPetscMat, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, Ii::Vector{PetscInt}, J::Vector{PetscInt}, v::Vector{PetscScalar}) end

@for_petsc function MatUpdateMPIAIJWithArrays(petsclib::$UnionPetscLib, mat::AbstractPetscMat, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, Ii::Vector{$PetscInt}, J::Vector{$PetscInt}, v::Vector{$PetscScalar} )

    @chk ccall(
               (:MatUpdateMPIAIJWithArrays, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}),
               mat, m, n, M, N, Ii, J, v,
              )


	return nothing
end 

"""
	MatView(petsclib::PetscLibType,mat::AbstractPetscMat, viewer::PetscViewer) 
display information about a matrix in a variety ways

Collective on viewer

Input Parameters:
- `mat`    - the matrix
- `viewer` - visualization context

Options Database Keys:
- `-mat_view ::ascii_info`           - Prints info on matrix at conclusion of `MatAssemblyEnd()`
- `-mat_view ::ascii_info_detail`    - Prints more detailed info
- `-mat_view`                        - Prints matrix in ASCII format
- `-mat_view ::ascii_matlab`         - Prints matrix in MATLAB format
- `-mat_view draw`                   - PetscDraws nonzero structure of matrix, using `MatView()` and `PetscDrawOpenX()`.
- `-display <name>`                  - Sets display name (default is host)
- `-draw_pause <sec>`                - Sets number of seconds to pause after display
- `-mat_view socket`                 - Sends matrix to socket, can be accessed from MATLAB (see Users-Manual: ch_matlab for details)
- `-viewer_socket_machine <machine>` - -
- `-viewer_socket_port <port>`       - -
- `-mat_view binary`                 - save matrix to file in binary format
- `-viewer_binary_filename <name>`   - -

Level: beginner

-seealso: [](ch_matrices), `Mat`, `PetscViewerPushFormat()`, `PetscViewerASCIIOpen()`, `PetscViewerDrawOpen()`, `PetscViewer`,
`PetscViewerSocketOpen()`, `PetscViewerBinaryOpen()`, `MatLoad()`, `MatViewFromOptions()`

# External Links
$(_doc_external("Mat/MatView"))
"""
function MatView(petsclib::PetscLibType, mat::AbstractPetscMat, viewer::PetscViewer) end

@for_petsc function MatView(petsclib::$UnionPetscLib, mat::AbstractPetscMat, viewer::PetscViewer )

    @chk ccall(
               (:MatView, $petsc_library),
               PetscErrorCode,
               (CMat, PetscViewer),
               mat, viewer,
              )


	return nothing
end 

"""
	MatViewFromOptions(petsclib::PetscLibType,A::AbstractPetscMat, obj, name::String) 
View properties of the matrix based on options set in the options database

Collective

Input Parameters:
- `A`    - the matrix
- `obj`  - optional additional object that provides the options prefix to use
- `name` - command line option

Options Database Key:
- `-mat_view [viewertype]:...` - the viewer and its options

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatView()`, `PetscObjectViewFromOptions()`, `MatCreate()`

# External Links
$(_doc_external("Mat/MatViewFromOptions"))
"""
function MatViewFromOptions(petsclib::PetscLibType, A::AbstractPetscMat, obj, name::String) end

@for_petsc function MatViewFromOptions(petsclib::$UnionPetscLib, A::AbstractPetscMat, obj, name::String )

    @chk ccall(
               (:MatViewFromOptions, $petsc_library),
               PetscErrorCode,
               (CMat, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

"""
	MatXAIJSetPreallocation(petsclib::PetscLibType,A::AbstractPetscMat, bs::PetscInt, dnnz::Union{Ptr, Vector{PetscInt}}, onnz::Union{Ptr, Vector{PetscInt}}, dnnzu::Union{Ptr, Vector{PetscInt}}, onnzu::Union{Ptr, Vector{PetscInt}}) 
set preallocation for serial and parallel `MATAIJ`, `MATBAIJ`, and `MATSBAIJ` matrices and their unassembled versions.

Collective

Input Parameters:
- `A`     - matrix being preallocated
- `bs`    - block size
- `dnnz`  - number of nonzero column blocks per block row of diagonal part of parallel matrix
- `onnz`  - number of nonzero column blocks per block row of off-diagonal part of parallel matrix
- `dnnzu` - number of nonzero column blocks per block row of upper-triangular part of diagonal part of parallel matrix
- `onnzu` - number of nonzero column blocks per block row of upper-triangular part of off-diagonal part of parallel matrix

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatSeqAIJSetPreallocation()`, `MatMPIAIJSetPreallocation()`, `MatSeqBAIJSetPreallocation()`, `MatMPIBAIJSetPreallocation()`,
`MatSeqSBAIJSetPreallocation()`, `MatMPISBAIJSetPreallocation()`,
`PetscSplitOwnership()`

# External Links
$(_doc_external("Mat/MatXAIJSetPreallocation"))
"""
function MatXAIJSetPreallocation(petsclib::PetscLibType, A::AbstractPetscMat, bs::PetscInt, dnnz::Union{Ptr, Vector{PetscInt}}, onnz::Union{Ptr, Vector{PetscInt}}, dnnzu::Union{Ptr, Vector{PetscInt}}, onnzu::Union{Ptr, Vector{PetscInt}}) end

@for_petsc function MatXAIJSetPreallocation(petsclib::$UnionPetscLib, A::AbstractPetscMat, bs::$PetscInt, dnnz::Union{Ptr, Vector{$PetscInt}}, onnz::Union{Ptr, Vector{$PetscInt}}, dnnzu::Union{Ptr, Vector{$PetscInt}}, onnzu::Union{Ptr, Vector{$PetscInt}} )

    @chk ccall(
               (:MatXAIJSetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
               A, bs, dnnz, onnz, dnnzu, onnzu,
              )


	return nothing
end 

"""
	MatZeroEntries(petsclib::PetscLibType,mat::AbstractPetscMat) 
Zeros all entries of a matrix.  For sparse matrices
this routine retains the old nonzero structure.

Logically Collective

Input Parameter:
- `mat` - the matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatZeroRows()`, `MatZeroRowsColumns()`

# External Links
$(_doc_external("Mat/MatZeroEntries"))
"""
function MatZeroEntries(petsclib::PetscLibType, mat::AbstractPetscMat) end

@for_petsc function MatZeroEntries(petsclib::$UnionPetscLib, mat::AbstractPetscMat )

    @chk ccall(
               (:MatZeroEntries, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	MatZeroRows(petsclib::PetscLibType,mat::AbstractPetscMat, numRows::PetscInt, rows::Vector{PetscInt}, diag::PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec) 
Zeros all entries (except possibly the main diagonal)
of a set of rows of a matrix.

Collective

Input Parameters:
- `mat`     - the matrix
- `numRows` - the number of rows to zero
- `rows`    - the global row indices
- `diag`    - value put in the diagonal of the zeroed rows
- `x`       - optional vector of solutions for zeroed rows (other entries in vector are not used), these must be set before this call
- `b`       - optional vector of right-hand side, that will be adjusted by provided solution entries

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatZeroRowsIS()`, `MatZeroRowsColumns()`, `MatZeroRowsLocalIS()`, `MatZeroRowsStencil()`, `MatZeroEntries()`, `MatZeroRowsLocal()`, `MatSetOption()`,
`MatZeroRowsColumnsLocal()`, `MatZeroRowsColumnsLocalIS()`, `MatZeroRowsColumnsIS()`, `MatZeroRowsColumnsStencil()`, `PCREDISTRIBUTE`, `MAT_KEEP_NONZERO_PATTERN`

# External Links
$(_doc_external("Mat/MatZeroRows"))
"""
function MatZeroRows(petsclib::PetscLibType, mat::AbstractPetscMat, numRows::PetscInt, rows::Vector{PetscInt}, diag::PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec) end

@for_petsc function MatZeroRows(petsclib::$UnionPetscLib, mat::AbstractPetscMat, numRows::$PetscInt, rows::Vector{$PetscInt}, diag::$PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec )

    @chk ccall(
               (:MatZeroRows, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscScalar, CVec, CVec),
               mat, numRows, rows, diag, x, b,
              )


	return nothing
end 

"""
	MatZeroRowsColumns(petsclib::PetscLibType,mat::AbstractPetscMat, numRows::PetscInt, rows::Vector{PetscInt}, diag::PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec) 
Zeros all entries (except possibly the main diagonal)
of a set of rows and columns of a matrix.

Collective

Input Parameters:
- `mat`     - the matrix
- `numRows` - the number of rows/columns to zero
- `rows`    - the global row indices
- `diag`    - value put in the diagonal of the eliminated rows
- `x`       - optional vector of the solution for zeroed rows (other entries in vector are not used), these must be set before this call
- `b`       - optional vector of the right-hand side, that will be adjusted by provided solution entries

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatZeroRowsIS()`, `MatZeroRows()`, `MatZeroRowsLocalIS()`, `MatZeroRowsStencil()`, `MatZeroEntries()`, `MatZeroRowsLocal()`, `MatSetOption()`,
`MatZeroRowsColumnsLocal()`, `MatZeroRowsColumnsLocalIS()`, `MatZeroRowsColumnsIS()`, `MatZeroRowsColumnsStencil()`

# External Links
$(_doc_external("Mat/MatZeroRowsColumns"))
"""
function MatZeroRowsColumns(petsclib::PetscLibType, mat::AbstractPetscMat, numRows::PetscInt, rows::Vector{PetscInt}, diag::PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec) end

@for_petsc function MatZeroRowsColumns(petsclib::$UnionPetscLib, mat::AbstractPetscMat, numRows::$PetscInt, rows::Vector{$PetscInt}, diag::$PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec )

    @chk ccall(
               (:MatZeroRowsColumns, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscScalar, CVec, CVec),
               mat, numRows, rows, diag, x, b,
              )


	return nothing
end 

"""
	MatZeroRowsColumnsIS(petsclib::PetscLibType,mat::AbstractPetscMat, is::AbstractIS, diag::PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec) 
Zeros all entries (except possibly the main diagonal)
of a set of rows and columns of a matrix.

Collective

Input Parameters:
- `mat`  - the matrix
- `is`   - the rows to zero
- `diag` - value put in all diagonals of eliminated rows (0.0 will even eliminate diagonal entry)
- `x`    - optional vector of solutions for zeroed rows (other entries in vector are not used)
- `b`    - optional vector of right-hand side, that will be adjusted by provided solution

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatZeroRowsIS()`, `MatZeroRowsColumns()`, `MatZeroRowsLocalIS()`, `MatZeroRowsStencil()`, `MatZeroEntries()`, `MatZeroRowsLocal()`, `MatSetOption()`,
`MatZeroRowsColumnsLocal()`, `MatZeroRowsColumnsLocalIS()`, `MatZeroRows()`, `MatZeroRowsColumnsStencil()`

# External Links
$(_doc_external("Mat/MatZeroRowsColumnsIS"))
"""
function MatZeroRowsColumnsIS(petsclib::PetscLibType, mat::AbstractPetscMat, is::AbstractIS, diag::PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec) end

@for_petsc function MatZeroRowsColumnsIS(petsclib::$UnionPetscLib, mat::AbstractPetscMat, is::AbstractIS, diag::$PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec )

    @chk ccall(
               (:MatZeroRowsColumnsIS, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, $PetscScalar, CVec, CVec),
               mat, is, diag, x, b,
              )


	return nothing
end 

"""
	MatZeroRowsColumnsLocal(petsclib::PetscLibType,mat::AbstractPetscMat, numRows::PetscInt, rows::Vector{PetscInt}, diag::PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec) 
Zeros all entries (except possibly the main diagonal)
of a set of rows and columns of a matrix; using local numbering of rows.

Collective

Input Parameters:
- `mat`     - the matrix
- `numRows` - the number of rows to remove
- `rows`    - the global row indices
- `diag`    - value put in all diagonals of eliminated rows
- `x`       - optional vector of solutions for zeroed rows (other entries in vector are not used)
- `b`       - optional vector of right-hand side, that will be adjusted by provided solution

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatZeroRowsIS()`, `MatZeroRowsColumns()`, `MatZeroRowsLocalIS()`, `MatZeroRowsStencil()`, `MatZeroEntries()`, `MatZeroRowsLocal()`, `MatSetOption()`,
`MatZeroRows()`, `MatZeroRowsColumnsLocalIS()`, `MatZeroRowsColumnsIS()`, `MatZeroRowsColumnsStencil()`

# External Links
$(_doc_external("Mat/MatZeroRowsColumnsLocal"))
"""
function MatZeroRowsColumnsLocal(petsclib::PetscLibType, mat::AbstractPetscMat, numRows::PetscInt, rows::Vector{PetscInt}, diag::PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec) end

@for_petsc function MatZeroRowsColumnsLocal(petsclib::$UnionPetscLib, mat::AbstractPetscMat, numRows::$PetscInt, rows::Vector{$PetscInt}, diag::$PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec )

    @chk ccall(
               (:MatZeroRowsColumnsLocal, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscScalar, CVec, CVec),
               mat, numRows, rows, diag, x, b,
              )


	return nothing
end 

"""
	MatZeroRowsColumnsLocalIS(petsclib::PetscLibType,mat::AbstractPetscMat, is::AbstractIS, diag::PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec) 
Zeros all entries (except possibly the main diagonal)
of a set of rows and columns of a matrix; using local numbering of rows.

Collective

Input Parameters:
- `mat`  - the matrix
- `is`   - index set of rows to remove
- `diag` - value put in all diagonals of eliminated rows
- `x`    - optional vector of solutions for zeroed rows (other entries in vector are not used)
- `b`    - optional vector of right-hand side, that will be adjusted by provided solution

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatZeroRowsIS()`, `MatZeroRowsColumns()`, `MatZeroRowsLocalIS()`, `MatZeroRowsStencil()`, `MatZeroEntries()`, `MatZeroRowsLocal()`, `MatSetOption()`,
`MatZeroRowsColumnsLocal()`, `MatZeroRows()`, `MatZeroRowsColumnsIS()`, `MatZeroRowsColumnsStencil()`

# External Links
$(_doc_external("Mat/MatZeroRowsColumnsLocalIS"))
"""
function MatZeroRowsColumnsLocalIS(petsclib::PetscLibType, mat::AbstractPetscMat, is::AbstractIS, diag::PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec) end

@for_petsc function MatZeroRowsColumnsLocalIS(petsclib::$UnionPetscLib, mat::AbstractPetscMat, is::AbstractIS, diag::$PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec )

    @chk ccall(
               (:MatZeroRowsColumnsLocalIS, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, $PetscScalar, CVec, CVec),
               mat, is, diag, x, b,
              )


	return nothing
end 

"""
	MatZeroRowsColumnsStencil(petsclib::PetscLibType,mat::AbstractPetscMat, numRows::PetscInt, rows::Vector{MatStencil}, diag::PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec) 
Zeros all row and column entries (except possibly the main diagonal)
of a set of rows and columns of a matrix.

Collective

Input Parameters:
- `mat`     - the matrix
- `numRows` - the number of rows/columns to remove
- `rows`    - the grid coordinates (and component number when dof > 1) for matrix rows
- `diag`    - value put in all diagonals of eliminated rows (0.0 will even eliminate diagonal entry)
- `x`       - optional vector of solutions for zeroed rows (other entries in vector are not used)
- `b`       - optional vector of right-hand side, that will be adjusted by provided solution

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatZeroRowsIS()`, `MatZeroRowsColumns()`, `MatZeroRowsLocalIS()`, `MatZeroRowsStencil()`, `MatZeroEntries()`, `MatZeroRowsLocal()`, `MatSetOption()`,
`MatZeroRowsColumnsLocal()`, `MatZeroRowsColumnsLocalIS()`, `MatZeroRowsColumnsIS()`, `MatZeroRows()`

# External Links
$(_doc_external("Mat/MatZeroRowsColumnsStencil"))
"""
function MatZeroRowsColumnsStencil(petsclib::PetscLibType, mat::AbstractPetscMat, numRows::PetscInt, rows::Vector{MatStencil}, diag::PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec) end

@for_petsc function MatZeroRowsColumnsStencil(petsclib::$UnionPetscLib, mat::AbstractPetscMat, numRows::$PetscInt, rows::Vector{MatStencil}, diag::$PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec )

    @chk ccall(
               (:MatZeroRowsColumnsStencil, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{MatStencil}, $PetscScalar, CVec, CVec),
               mat, numRows, rows, diag, x, b,
              )


	return nothing
end 

"""
	MatZeroRowsIS(petsclib::PetscLibType,mat::AbstractPetscMat, is::Union{Ptr, AbstractIS}, diag::PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec) 
Zeros all entries (except possibly the main diagonal)
of a set of rows of a matrix indicated by an `IS`

Collective

Input Parameters:
- `mat`  - the matrix
- `is`   - index set, `IS`, of rows to remove (if `NULL` then no row is removed)
- `diag` - value put in all diagonals of eliminated rows
- `x`    - optional vector of solutions for zeroed rows (other entries in vector are not used)
- `b`    - optional vector of right-hand side, that will be adjusted by provided solution

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatZeroRows()`, `MatZeroRowsColumns()`, `MatZeroRowsLocalIS()`, `MatZeroRowsStencil()`, `MatZeroEntries()`, `MatZeroRowsLocal()`, `MatSetOption()`,
`MatZeroRowsColumnsLocal()`, `MatZeroRowsColumnsLocalIS()`, `MatZeroRowsColumnsIS()`, `MatZeroRowsColumnsStencil()`, `IS`

# External Links
$(_doc_external("Mat/MatZeroRowsIS"))
"""
function MatZeroRowsIS(petsclib::PetscLibType, mat::AbstractPetscMat, is::Union{Ptr, AbstractIS}, diag::PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec) end

@for_petsc function MatZeroRowsIS(petsclib::$UnionPetscLib, mat::AbstractPetscMat, is::Union{Ptr, AbstractIS}, diag::$PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec )

    @chk ccall(
               (:MatZeroRowsIS, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, $PetscScalar, CVec, CVec),
               mat, is, diag, x, b,
              )


	return nothing
end 

"""
	MatZeroRowsLocal(petsclib::PetscLibType,mat::AbstractPetscMat, numRows::PetscInt, rows::Vector{PetscInt}, diag::PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec) 
Zeros all entries (except possibly the main diagonal)
of a set of rows of a matrix; using local numbering of rows.

Collective

Input Parameters:
- `mat`     - the matrix
- `numRows` - the number of rows to remove
- `rows`    - the local row indices
- `diag`    - value put in all diagonals of eliminated rows
- `x`       - optional vector of solutions for zeroed rows (other entries in vector are not used)
- `b`       - optional vector of right-hand side, that will be adjusted by provided solution

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatZeroRowsIS()`, `MatZeroRowsColumns()`, `MatZeroRowsLocalIS()`, `MatZeroRowsStencil()`, `MatZeroEntries()`, `MatZeroRows()`, `MatSetOption()`,
`MatZeroRowsColumnsLocal()`, `MatZeroRowsColumnsLocalIS()`, `MatZeroRowsColumnsIS()`, `MatZeroRowsColumnsStencil()`

# External Links
$(_doc_external("Mat/MatZeroRowsLocal"))
"""
function MatZeroRowsLocal(petsclib::PetscLibType, mat::AbstractPetscMat, numRows::PetscInt, rows::Vector{PetscInt}, diag::PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec) end

@for_petsc function MatZeroRowsLocal(petsclib::$UnionPetscLib, mat::AbstractPetscMat, numRows::$PetscInt, rows::Vector{$PetscInt}, diag::$PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec )

    @chk ccall(
               (:MatZeroRowsLocal, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscScalar, CVec, CVec),
               mat, numRows, rows, diag, x, b,
              )


	return nothing
end 

"""
	MatZeroRowsLocalIS(petsclib::PetscLibType,mat::AbstractPetscMat, is::AbstractIS, diag::PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec) 
Zeros all entries (except possibly the main diagonal)
of a set of rows of a matrix; using local numbering of rows.

Collective

Input Parameters:
- `mat`  - the matrix
- `is`   - index set of rows to remove
- `diag` - value put in all diagonals of eliminated rows
- `x`    - optional vector of solutions for zeroed rows (other entries in vector are not used)
- `b`    - optional vector of right-hand side, that will be adjusted by provided solution

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatZeroRowsIS()`, `MatZeroRowsColumns()`, `MatZeroRows()`, `MatZeroRowsStencil()`, `MatZeroEntries()`, `MatZeroRowsLocal()`, `MatSetOption()`,
`MatZeroRowsColumnsLocal()`, `MatZeroRowsColumnsLocalIS()`, `MatZeroRowsColumnsIS()`, `MatZeroRowsColumnsStencil()`

# External Links
$(_doc_external("Mat/MatZeroRowsLocalIS"))
"""
function MatZeroRowsLocalIS(petsclib::PetscLibType, mat::AbstractPetscMat, is::AbstractIS, diag::PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec) end

@for_petsc function MatZeroRowsLocalIS(petsclib::$UnionPetscLib, mat::AbstractPetscMat, is::AbstractIS, diag::$PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec )

    @chk ccall(
               (:MatZeroRowsLocalIS, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, $PetscScalar, CVec, CVec),
               mat, is, diag, x, b,
              )


	return nothing
end 

"""
	MatZeroRowsStencil(petsclib::PetscLibType,mat::AbstractPetscMat, numRows::PetscInt, rows::Vector{MatStencil}, diag::PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec) 
Zeros all entries (except possibly the main diagonal)
of a set of rows of a matrix indicated by a `MatStencil`. These rows must be local to the process.

Collective

Input Parameters:
- `mat`     - the matrix
- `numRows` - the number of rows to remove
- `rows`    - the grid coordinates (and component number when dof > 1) for matrix rows indicated by an array of `MatStencil`
- `diag`    - value put in all diagonals of eliminated rows (0.0 will even eliminate diagonal entry)
- `x`       - optional vector of solutions for zeroed rows (other entries in vector are not used)
- `b`       - optional vector of right-hand side, that will be adjusted by provided solution

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatStencil`, `MatZeroRowsIS()`, `MatZeroRowsColumns()`, `MatZeroRowsLocalIS()`, `MatZeroRows()`, `MatZeroEntries()`, `MatZeroRowsLocal()`, `MatSetOption()`,
`MatZeroRowsColumnsLocal()`, `MatZeroRowsColumnsLocalIS()`, `MatZeroRowsColumnsIS()`, `MatZeroRowsColumnsStencil()`

# External Links
$(_doc_external("Mat/MatZeroRowsStencil"))
"""
function MatZeroRowsStencil(petsclib::PetscLibType, mat::AbstractPetscMat, numRows::PetscInt, rows::Vector{MatStencil}, diag::PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec) end

@for_petsc function MatZeroRowsStencil(petsclib::$UnionPetscLib, mat::AbstractPetscMat, numRows::$PetscInt, rows::Vector{MatStencil}, diag::$PetscScalar, x::AbstractPetscVec, b::AbstractPetscVec )

    @chk ccall(
               (:MatZeroRowsStencil, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{MatStencil}, $PetscScalar, CVec, CVec),
               mat, numRows, rows, diag, x, b,
              )


	return nothing
end 

