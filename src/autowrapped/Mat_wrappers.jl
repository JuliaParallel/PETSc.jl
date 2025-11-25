# autodefined type arguments for class Mat ------
mutable struct _n_MatNullSpace end
const MatNullSpace = Ptr{_n_MatNullSpace}

mutable struct _n_MatTransposeColoring end
const MatTransposeColoring = Ptr{_n_MatTransposeColoring}

mutable struct _n_hypre_ParCSRMatrix end
const hypre_ParCSRMatrix = Ptr{_n_hypre_ParCSRMatrix}

mutable struct _n_PetscFunctionList end
const PetscFunctionList = Ptr{_n_PetscFunctionList}

mutable struct MatHtoolKernelFn end

# -------------------------------------------------------
"""
	MatSetType(petsclib::PetscLibType,mat::PetscMat, matype::MatType) 
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
function MatSetType(petsclib::PetscLibType, mat::PetscMat, matype::MatType) end

@for_petsc function MatSetType(petsclib::$UnionPetscLib, mat::PetscMat, matype::MatType )

    @chk ccall(
               (:MatSetType, $petsc_library),
               PetscErrorCode,
               (CMat, MatType),
               mat, matype,
              )


	return nothing
end 

"""
	type::MatType = MatGetType(petsclib::PetscLibType,mat::PetscMat) 
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

	type = unsafe_string(type_[])

	return type
end 

"""
	vtype::VecType = MatGetVecType(petsclib::PetscLibType,mat::PetscMat) 
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
function MatGetVecType(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatGetVecType(petsclib::$UnionPetscLib, mat::PetscMat )
	vtype_ = Ref{VecType}()

    @chk ccall(
               (:MatGetVecType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{VecType}),
               mat, vtype_,
              )

	vtype = unsafe_string(vtype_[])

	return vtype
end 

"""
	MatSetVecType(petsclib::PetscLibType,mat::PetscMat, vtype::VecType) 
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
function MatSetVecType(petsclib::PetscLibType, mat::PetscMat, vtype::VecType) end

@for_petsc function MatSetVecType(petsclib::$UnionPetscLib, mat::PetscMat, vtype::VecType )

    @chk ccall(
               (:MatSetVecType, $petsc_library),
               PetscErrorCode,
               (CMat, VecType),
               mat, vtype,
              )


	return nothing
end 

"""
	MatRegister(petsclib::PetscLibType,sname::Vector{Cchar}, fnc::external) 
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
function MatRegister(petsclib::PetscLibType, sname::Vector{Cchar}, fnc::external) end

@for_petsc function MatRegister(petsclib::$UnionPetscLib, sname::Vector{Cchar}, fnc::external )

    @chk ccall(
               (:MatRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	MatRegisterRootName(petsclib::PetscLibType,rname::Vector{Cchar}, sname::Vector{Cchar}, mname::Vector{Cchar}) 
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
function MatRegisterRootName(petsclib::PetscLibType, rname::Vector{Cchar}, sname::Vector{Cchar}, mname::Vector{Cchar}) end

@for_petsc function MatRegisterRootName(petsclib::$UnionPetscLib, rname::Vector{Cchar}, sname::Vector{Cchar}, mname::Vector{Cchar} )

    @chk ccall(
               (:MatRegisterRootName, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}),
               rname, sname, mname,
              )


	return nothing
end 

"""
	MatProductReplaceMats(petsclib::PetscLibType,A::PetscMat, B::PetscMat, C::PetscMat, D::PetscMat) 
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
function MatProductReplaceMats(petsclib::PetscLibType, A::Union{Ptr,PetscMat}, B::Union{Ptr,PetscMat}, C::Union{Ptr,PetscMat}, D::Union{Ptr,PetscMat}) end

@for_petsc function MatProductReplaceMats(petsclib::$UnionPetscLib, A::Union{Ptr,PetscMat}, B::Union{Ptr,PetscMat}, C::Union{Ptr,PetscMat}, D::Union{Ptr,PetscMat})

    @chk ccall(
               (:MatProductReplaceMats, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat, CMat),
               A, B, C, D,
              )


	return nothing
end 

"""
	MatProductSetFromOptions(petsclib::PetscLibType,mat::PetscMat) 
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
function MatProductSetFromOptions(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatProductSetFromOptions(petsclib::$UnionPetscLib, mat::PetscMat )

    @chk ccall(
               (:MatProductSetFromOptions, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	MatProductView(petsclib::PetscLibType,mat::PetscMat, viewer::PetscViewer) 
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
function MatProductView(petsclib::PetscLibType, mat::PetscMat, viewer::PetscViewer) end

@for_petsc function MatProductView(petsclib::$UnionPetscLib, mat::PetscMat, viewer::PetscViewer )

    @chk ccall(
               (:MatProductView, $petsc_library),
               PetscErrorCode,
               (CMat, PetscViewer),
               mat, viewer,
              )


	return nothing
end 

"""
	MatProductNumeric(petsclib::PetscLibType,mat::PetscMat) 
Compute a matrix

Collective

Input/Output Parameter:
- `mat` - the matrix whose values are computed via a matrix-matrix product operation

Level: intermediate

-seealso: [](ch_matrices), `MatProduct`, `Mat`, `MatProductSetAlgorithm()`, `MatProductSetType()`, `MatProductCreate()`, `MatSetType()`, `MatProductSymbolic()`

# External Links
$(_doc_external("Mat/MatProductNumeric"))
"""
function MatProductNumeric(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatProductNumeric(petsclib::$UnionPetscLib, mat::PetscMat )

    @chk ccall(
               (:MatProductNumeric, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	MatProductSymbolic(petsclib::PetscLibType,mat::PetscMat) 
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
function MatProductSymbolic(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatProductSymbolic(petsclib::$UnionPetscLib, mat::PetscMat )

    @chk ccall(
               (:MatProductSymbolic, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	MatProductSetFill(petsclib::PetscLibType,mat::PetscMat, fill::PetscReal) 
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
function MatProductSetFill(petsclib::PetscLibType, mat::PetscMat, fill::PetscReal) end

@for_petsc function MatProductSetFill(petsclib::$UnionPetscLib, mat::PetscMat, fill::$PetscReal )

    @chk ccall(
               (:MatProductSetFill, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscReal),
               mat, fill,
              )


	return nothing
end 

"""
	MatProductSetAlgorithm(petsclib::PetscLibType,mat::PetscMat, alg::MatProductAlgorithm) 
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
function MatProductSetAlgorithm(petsclib::PetscLibType, mat::PetscMat, alg::MatProductAlgorithm) end

@for_petsc function MatProductSetAlgorithm(petsclib::$UnionPetscLib, mat::PetscMat, alg::MatProductAlgorithm )

    @chk ccall(
               (:MatProductSetAlgorithm, $petsc_library),
               PetscErrorCode,
               (CMat, MatProductAlgorithm),
               mat, alg,
              )


	return nothing
end 

"""
	MatProductGetAlgorithm(petsclib::PetscLibType,mat::PetscMat, alg::MatProductAlgorithm) 
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
function MatProductGetAlgorithm(petsclib::PetscLibType, mat::PetscMat, alg::MatProductAlgorithm) end

@for_petsc function MatProductGetAlgorithm(petsclib::$UnionPetscLib, mat::PetscMat, alg::MatProductAlgorithm )

    @chk ccall(
               (:MatProductGetAlgorithm, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatProductAlgorithm}),
               mat, alg,
              )


	return nothing
end 

"""
	MatProductSetType(petsclib::PetscLibType,mat::PetscMat, productype::MatProductType) 
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
function MatProductSetType(petsclib::PetscLibType, mat::PetscMat, productype::MatProductType) end

@for_petsc function MatProductSetType(petsclib::$UnionPetscLib, mat::PetscMat, productype::MatProductType )

    @chk ccall(
               (:MatProductSetType, $petsc_library),
               PetscErrorCode,
               (CMat, MatProductType),
               mat, productype,
              )


	return nothing
end 

"""
	MatProductClear(petsclib::PetscLibType,mat::PetscMat) 
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
function MatProductClear(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatProductClear(petsclib::$UnionPetscLib, mat::PetscMat )

    @chk ccall(
               (:MatProductClear, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	MatProductCreateWithMat(petsclib::PetscLibType,A::PetscMat, B::PetscMat, C::PetscMat, D::PetscMat) 
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
function MatProductCreateWithMat(petsclib::PetscLibType, A::PetscMat, B::PetscMat, C::Union{Ptr,PetscMat}, D::PetscMat) end

@for_petsc function MatProductCreateWithMat(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat, C::Union{Ptr,PetscMat}, D::PetscMat )

    @chk ccall(
               (:MatProductCreateWithMat, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat, CMat),
               A, B, C, D,
              )


	return nothing
end 

"""
	D::PetscMat = MatProductCreate(petsclib::PetscLibType,A::PetscMat, B::PetscMat, C::PetscMat) 
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
function MatProductCreate(petsclib::PetscLibType, A::PetscMat, B::PetscMat, C::Union{Ptr,PetscMat}) end

@for_petsc function MatProductCreate(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat, C::Union{Ptr,PetscMat})
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
	mtype::MatProductType = MatProductGetType(petsclib::PetscLibType,mat::PetscMat) 
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
function MatProductGetType(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatProductGetType(petsclib::$UnionPetscLib, mat::PetscMat )
	mtype_ = Ref{MatProductType}()

    @chk ccall(
               (:MatProductGetType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatProductType}),
               mat, mtype_,
              )

	mtype = unsafe_string(mtype_[])

	return mtype
end 

"""
	MatProductGetMats(petsclib::PetscLibType,mat::PetscMat, A::PetscMat, B::PetscMat, C::PetscMat) 
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
function MatProductGetMats(petsclib::PetscLibType, mat::PetscMat, A::PetscMat, B::PetscMat, C::Union{Ptr,PetscMat}) end

@for_petsc function MatProductGetMats(petsclib::$UnionPetscLib, mat::PetscMat, A::PetscMat, B::PetscMat, C::Union{Ptr,PetscMat})
	A_ = Ref(A.ptr)
	B_ = Ref(B.ptr)
	C_ = Ref(C.ptr)

    @chk ccall(
               (:MatProductGetMats, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}, Ptr{CMat}, Ptr{CMat}),
               mat, A_, B_, C_,
              )

	A.ptr = C_NULL
	B.ptr = C_NULL
	C.ptr = C_NULL

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
	MatSetRandom(petsclib::PetscLibType,x::PetscMat, rctx::PetscRandom) 
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
function MatSetRandom(petsclib::PetscLibType, x::PetscMat, rctx::Union{Ptr,PetscRandom}) end

@for_petsc function MatSetRandom(petsclib::$UnionPetscLib, x::PetscMat, rctx::Union{Ptr,PetscRandom})

    @chk ccall(
               (:MatSetRandom, $petsc_library),
               PetscErrorCode,
               (CMat, PetscRandom),
               x, rctx,
              )


	return nothing
end 

"""
	MatCopyHashToXAIJ(petsclib::PetscLibType,A::PetscMat, B::PetscMat) 
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
function MatCopyHashToXAIJ(petsclib::PetscLibType, A::PetscMat, B::PetscMat) end

@for_petsc function MatCopyHashToXAIJ(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat )

    @chk ccall(
               (:MatCopyHashToXAIJ, $petsc_library),
               PetscErrorCode,
               (CMat, CMat),
               A, B,
              )


	return nothing
end 

"""
	pivot::PetscReal,row::PetscInt = MatFactorGetErrorZeroPivot(petsclib::PetscLibType,mat::PetscMat) 
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
function MatFactorGetErrorZeroPivot(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatFactorGetErrorZeroPivot(petsclib::$UnionPetscLib, mat::PetscMat )
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
	MatFactorGetError(petsclib::PetscLibType,mat::PetscMat, err::MatFactorError) 
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
function MatFactorGetError(petsclib::PetscLibType, mat::PetscMat, err::MatFactorError) end

@for_petsc function MatFactorGetError(petsclib::$UnionPetscLib, mat::PetscMat, err::MatFactorError )

    @chk ccall(
               (:MatFactorGetError, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatFactorError}),
               mat, err,
              )


	return nothing
end 

"""
	MatFactorClearError(petsclib::PetscLibType,mat::PetscMat) 
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
function MatFactorClearError(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatFactorClearError(petsclib::$UnionPetscLib, mat::PetscMat )

    @chk ccall(
               (:MatFactorClearError, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	MatFindNonzeroRows(petsclib::PetscLibType,mat::PetscMat, keptrows::IS) 
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
function MatFindNonzeroRows(petsclib::PetscLibType, mat::PetscMat, keptrows::IS) end

@for_petsc function MatFindNonzeroRows(petsclib::$UnionPetscLib, mat::PetscMat, keptrows::IS )

    @chk ccall(
               (:MatFindNonzeroRows, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CIS}),
               mat, keptrows,
              )


	return nothing
end 

"""
	MatFindZeroRows(petsclib::PetscLibType,mat::PetscMat, zerorows::IS) 
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
function MatFindZeroRows(petsclib::PetscLibType, mat::PetscMat, zerorows::IS) end

@for_petsc function MatFindZeroRows(petsclib::$UnionPetscLib, mat::PetscMat, zerorows::IS )

    @chk ccall(
               (:MatFindZeroRows, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CIS}),
               mat, zerorows,
              )


	return nothing
end 

"""
	MatGetDiagonalBlock(petsclib::PetscLibType,A::PetscMat, a::PetscMat) 
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
function MatGetDiagonalBlock(petsclib::PetscLibType, A::PetscMat, a::PetscMat) end

@for_petsc function MatGetDiagonalBlock(petsclib::$UnionPetscLib, A::PetscMat, a::PetscMat )
	a_ = Ref(a.ptr)

    @chk ccall(
               (:MatGetDiagonalBlock, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, a_,
              )

	a.ptr = C_NULL

	return nothing
end 

"""
	trace::PetscScalar = MatGetTrace(petsclib::PetscLibType,mat::PetscMat) 
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
function MatGetTrace(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatGetTrace(petsclib::$UnionPetscLib, mat::PetscMat )
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
	MatRealPart(petsclib::PetscLibType,mat::PetscMat) 
Zeros out the imaginary part of the matrix

Logically Collective

Input Parameter:
- `mat` - the matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatImaginaryPart()`

# External Links
$(_doc_external("Mat/MatRealPart"))
"""
function MatRealPart(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatRealPart(petsclib::$UnionPetscLib, mat::PetscMat )

    @chk ccall(
               (:MatRealPart, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	nghosts::PetscInt,ghosts::Vector{PetscInt} = MatGetGhosts(petsclib::PetscLibType,mat::PetscMat) 
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
function MatGetGhosts(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatGetGhosts(petsclib::$UnionPetscLib, mat::PetscMat )
	nghosts_ = Ref{$PetscInt}()
	ghosts_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:MatGetGhosts, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
               mat, nghosts_, ghosts_,
              )

	nghosts = nghosts_[]
	ghosts = unsafe_wrap(Array, ghosts_[], VecGetLocalSize(petsclib, x); own = false)

	return nghosts,ghosts
end 

"""
	MatImaginaryPart(petsclib::PetscLibType,mat::PetscMat) 
Moves the imaginary part of the matrix to the real part and zeros the imaginary part

Logically Collective

Input Parameter:
- `mat` - the matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatRealPart()`

# External Links
$(_doc_external("Mat/MatImaginaryPart"))
"""
function MatImaginaryPart(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatImaginaryPart(petsclib::$UnionPetscLib, mat::PetscMat )

    @chk ccall(
               (:MatImaginaryPart, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	missing::PetscBool,dd::PetscInt = MatMissingDiagonal(petsclib::PetscLibType,mat::PetscMat) 
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
function MatMissingDiagonal(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatMissingDiagonal(petsclib::$UnionPetscLib, mat::PetscMat )
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
	ncols::PetscInt,cols::Vector{PetscInt},vals::Vector{PetscScalar} = MatGetRow(petsclib::PetscLibType,mat::PetscMat, row::PetscInt) 
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
function MatGetRow(petsclib::PetscLibType, mat::PetscMat, row::PetscInt) end

@for_petsc function MatGetRow(petsclib::$UnionPetscLib, mat::PetscMat, row::$PetscInt )
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
	cols = unsafe_wrap(Array, cols_[], VecGetLocalSize(petsclib, x); own = false)
	vals = unsafe_wrap(Array, vals_[], VecGetLocalSize(petsclib, x); own = false)

	return ncols,cols,vals
end 

"""
	MatConjugate(petsclib::PetscLibType,mat::PetscMat) 
replaces the matrix values with their complex conjugates

Logically Collective

Input Parameter:
- `mat` - the matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatRealPart()`, `MatImaginaryPart()`, `VecConjugate()`, `MatTranspose()`

# External Links
$(_doc_external("Mat/MatConjugate"))
"""
function MatConjugate(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatConjugate(petsclib::$UnionPetscLib, mat::PetscMat )

    @chk ccall(
               (:MatConjugate, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	ncols::PetscInt,cols::Vector{PetscInt},vals::Vector{PetscScalar} = MatRestoreRow(petsclib::PetscLibType,mat::PetscMat, row::PetscInt) 
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
function MatRestoreRow(petsclib::PetscLibType, mat::PetscMat, row::PetscInt) end

@for_petsc function MatRestoreRow(petsclib::$UnionPetscLib, mat::PetscMat, row::$PetscInt )
	ncols_ = Ref{$PetscInt}()
	cols_ = Ref{Ptr{$PetscInt}}()
	vals_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatRestoreRow, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscScalar}}),
               mat, row, ncols_, cols_, vals_,
              )

	ncols = ncols_[]
	cols = unsafe_wrap(Array, cols_[], VecGetLocalSize(petsclib, x); own = false)
	vals = unsafe_wrap(Array, vals_[], VecGetLocalSize(petsclib, x); own = false)

	return ncols,cols,vals
end 

"""
	MatGetRowUpperTriangular(petsclib::PetscLibType,mat::PetscMat) 
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
function MatGetRowUpperTriangular(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatGetRowUpperTriangular(petsclib::$UnionPetscLib, mat::PetscMat )

    @chk ccall(
               (:MatGetRowUpperTriangular, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	MatRestoreRowUpperTriangular(petsclib::PetscLibType,mat::PetscMat) 
Disable calls to `MatGetRow()` for matrix in `MATSBAIJ` format.

Not Collective

Input Parameter:
- `mat` - the matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATSBAIJ`, `MatGetRowUpperTriangular()`

# External Links
$(_doc_external("Mat/MatRestoreRowUpperTriangular"))
"""
function MatRestoreRowUpperTriangular(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatRestoreRowUpperTriangular(petsclib::$UnionPetscLib, mat::PetscMat )

    @chk ccall(
               (:MatRestoreRowUpperTriangular, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	MatSetOptionsPrefix(petsclib::PetscLibType,A::PetscMat, prefix::Vector{Cchar}) 
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
function MatSetOptionsPrefix(petsclib::PetscLibType, A::PetscMat, prefix::Vector{Cchar}) end

@for_petsc function MatSetOptionsPrefix(petsclib::$UnionPetscLib, A::PetscMat, prefix::Vector{Cchar} )

    @chk ccall(
               (:MatSetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Cchar}),
               A, prefix,
              )


	return nothing
end 

"""
	MatSetOptionsPrefixFactor(petsclib::PetscLibType,A::PetscMat, prefix::Vector{Cchar}) 
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
function MatSetOptionsPrefixFactor(petsclib::PetscLibType, A::PetscMat, prefix::Vector{Cchar}) end

@for_petsc function MatSetOptionsPrefixFactor(petsclib::$UnionPetscLib, A::PetscMat, prefix::Vector{Cchar} )

    @chk ccall(
               (:MatSetOptionsPrefixFactor, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Cchar}),
               A, prefix,
              )


	return nothing
end 

"""
	MatAppendOptionsPrefixFactor(petsclib::PetscLibType,A::PetscMat, prefix::Vector{Cchar}) 
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
function MatAppendOptionsPrefixFactor(petsclib::PetscLibType, A::PetscMat, prefix::Vector{Cchar}) end

@for_petsc function MatAppendOptionsPrefixFactor(petsclib::$UnionPetscLib, A::PetscMat, prefix::Vector{Cchar} )

    @chk ccall(
               (:MatAppendOptionsPrefixFactor, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Cchar}),
               A, prefix,
              )


	return nothing
end 

"""
	MatAppendOptionsPrefix(petsclib::PetscLibType,A::PetscMat, prefix::Vector{Cchar}) 
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
function MatAppendOptionsPrefix(petsclib::PetscLibType, A::PetscMat, prefix::Vector{Cchar}) end

@for_petsc function MatAppendOptionsPrefix(petsclib::$UnionPetscLib, A::PetscMat, prefix::Vector{Cchar} )

    @chk ccall(
               (:MatAppendOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Cchar}),
               A, prefix,
              )


	return nothing
end 

"""
	MatGetOptionsPrefix(petsclib::PetscLibType,A::PetscMat, prefix::Vector{Cchar}) 
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
function MatGetOptionsPrefix(petsclib::PetscLibType, A::PetscMat, prefix::Vector{Cchar}) end

@for_petsc function MatGetOptionsPrefix(petsclib::$UnionPetscLib, A::PetscMat, prefix::Vector{Cchar} )
	prefix_ = Ref(pointer(prefix))

    @chk ccall(
               (:MatGetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{Cchar}}),
               A, prefix_,
              )


	return nothing
end 

"""
	MatGetState(petsclib::PetscLibType,A::PetscMat, state::PetscObjectState) 
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
function MatGetState(petsclib::PetscLibType, A::PetscMat, state::PetscObjectState) end

@for_petsc function MatGetState(petsclib::$UnionPetscLib, A::PetscMat, state::PetscObjectState )

    @chk ccall(
               (:MatGetState, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscObjectState}),
               A, state,
              )


	return nothing
end 

"""
	MatResetPreallocation(petsclib::PetscLibType,A::PetscMat) 
Reset matrix to use the original preallocation values provided by the user, for example with `MatXAIJSetPreallocation()`

Collective

Input Parameter:
- `A` - the matrix

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatSeqAIJSetPreallocation()`, `MatMPIAIJSetPreallocation()`, `MatXAIJSetPreallocation()`

# External Links
$(_doc_external("Mat/MatResetPreallocation"))
"""
function MatResetPreallocation(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatResetPreallocation(petsclib::$UnionPetscLib, A::PetscMat )

    @chk ccall(
               (:MatResetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat,),
               A,
              )


	return nothing
end 

"""
	MatResetHash(petsclib::PetscLibType,A::PetscMat) 
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
	MatSetUp(petsclib::PetscLibType,A::PetscMat) 
Sets up the internal matrix data structures for later use by the matrix

Collective

Input Parameter:
- `A` - the matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatMult()`, `MatCreate()`, `MatDestroy()`, `MatXAIJSetPreallocation()`

# External Links
$(_doc_external("Mat/MatSetUp"))
"""
function MatSetUp(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatSetUp(petsclib::$UnionPetscLib, A::PetscMat )

    @chk ccall(
               (:MatSetUp, $petsc_library),
               PetscErrorCode,
               (CMat,),
               A,
              )


	return nothing
end 

"""
	MatViewFromOptions(petsclib::PetscLibType,A::PetscMat, obj::PetscObject, name::Vector{Cchar}) 
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
function MatViewFromOptions(petsclib::PetscLibType, A::PetscMat, obj::PetscObject, name::Vector{Cchar}) end

@for_petsc function MatViewFromOptions(petsclib::$UnionPetscLib, A::PetscMat, obj::PetscObject, name::Vector{Cchar} )

    @chk ccall(
               (:MatViewFromOptions, $petsc_library),
               PetscErrorCode,
               (CMat, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

"""
	MatView(petsclib::PetscLibType,mat::PetscMat, viewer::PetscViewer) 
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
function MatView(petsclib::PetscLibType, mat::PetscMat, viewer::PetscViewer) end

@for_petsc function MatView(petsclib::$UnionPetscLib, mat::PetscMat, viewer::PetscViewer )

    @chk ccall(
               (:MatView, $petsc_library),
               PetscErrorCode,
               (CMat, PetscViewer),
               mat, viewer,
              )


	return nothing
end 

"""
	MatLoad(petsclib::PetscLibType,mat::PetscMat, viewer::PetscViewer) 
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
function MatLoad(petsclib::PetscLibType, mat::PetscMat, viewer::PetscViewer) end

@for_petsc function MatLoad(petsclib::$UnionPetscLib, mat::PetscMat, viewer::PetscViewer )

    @chk ccall(
               (:MatLoad, $petsc_library),
               PetscErrorCode,
               (CMat, PetscViewer),
               mat, viewer,
              )


	return nothing
end 

"""
	MatDestroy(petsclib::PetscLibType,A::PetscMat) 
Frees space taken by a matrix.

Collective

Input Parameter:
- `A` - the matrix

Level: beginner

-seealso: [](ch_matrices), `Mat`, `MatCreate()`

# External Links
$(_doc_external("Mat/MatDestroy"))
"""
function MatDestroy(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatDestroy(petsclib::$UnionPetscLib, A::PetscMat )
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
	MatSetValuesIS(petsclib::PetscLibType,mat::PetscMat, ism::IS, isn::IS, v::Vector{PetscScalar}, addv::InsertMode) 
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
function MatSetValuesIS(petsclib::PetscLibType, mat::PetscMat, ism::IS, isn::IS, v::Vector{PetscScalar}, addv::InsertMode) end

@for_petsc function MatSetValuesIS(petsclib::$UnionPetscLib, mat::PetscMat, ism::IS, isn::IS, v::Vector{$PetscScalar}, addv::InsertMode )

    @chk ccall(
               (:MatSetValuesIS, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, CIS, Ptr{$PetscScalar}, InsertMode),
               mat, ism, isn, v, addv,
              )


	return nothing
end 

"""
	MatSetValuesRowLocal(petsclib::PetscLibType,mat::PetscMat, row::PetscInt, v::Vector{PetscScalar}) 
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
function MatSetValuesRowLocal(petsclib::PetscLibType, mat::PetscMat, row::PetscInt, v::Vector{PetscScalar}) end

@for_petsc function MatSetValuesRowLocal(petsclib::$UnionPetscLib, mat::PetscMat, row::$PetscInt, v::Vector{$PetscScalar} )

    @chk ccall(
               (:MatSetValuesRowLocal, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscScalar}),
               mat, row, v,
              )


	return nothing
end 

"""
	MatSetValuesRow(petsclib::PetscLibType,mat::PetscMat, row::PetscInt, v::Vector{PetscScalar}) 
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
function MatSetValuesRow(petsclib::PetscLibType, mat::PetscMat, row::PetscInt, v::Vector{PetscScalar}) end

@for_petsc function MatSetValuesRow(petsclib::$UnionPetscLib, mat::PetscMat, row::$PetscInt, v::Vector{$PetscScalar} )

    @chk ccall(
               (:MatSetValuesRow, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscScalar}),
               mat, row, v,
              )


	return nothing
end 

"""
	MatSetValuesStencil(petsclib::PetscLibType,mat::PetscMat, m::PetscInt, idxm::Vector{MatStencil}, n::PetscInt, idxn::Vector{MatStencil}, v::Vector{PetscScalar}, addv::InsertMode) 
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
function MatSetValuesStencil(petsclib::PetscLibType, mat::PetscMat, m::PetscInt, idxm::Vector{MatStencil}, n::PetscInt, idxn::Vector{MatStencil}, v::Vector{PetscScalar}, addv::InsertMode) end

@for_petsc function MatSetValuesStencil(petsclib::$UnionPetscLib, mat::PetscMat, m::$PetscInt, idxm::Vector{MatStencil}, n::$PetscInt, idxn::Vector{MatStencil}, v::Vector{$PetscScalar}, addv::InsertMode )

    @chk ccall(
               (:MatSetValuesStencil, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{MatStencil}, $PetscInt, Ptr{MatStencil}, Ptr{$PetscScalar}, InsertMode),
               mat, m, idxm, n, idxn, v, addv,
              )


	return nothing
end 

"""
	MatSetValuesBlockedStencil(petsclib::PetscLibType,mat::PetscMat, m::PetscInt, idxm::Vector{MatStencil}, n::PetscInt, idxn::Vector{MatStencil}, v::Vector{PetscScalar}, addv::InsertMode) 
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
function MatSetValuesBlockedStencil(petsclib::PetscLibType, mat::PetscMat, m::PetscInt, idxm::Vector{MatStencil}, n::PetscInt, idxn::Vector{MatStencil}, v::Vector{PetscScalar}, addv::InsertMode) end

@for_petsc function MatSetValuesBlockedStencil(petsclib::$UnionPetscLib, mat::PetscMat, m::$PetscInt, idxm::Vector{MatStencil}, n::$PetscInt, idxn::Vector{MatStencil}, v::Vector{$PetscScalar}, addv::InsertMode )

    @chk ccall(
               (:MatSetValuesBlockedStencil, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{MatStencil}, $PetscInt, Ptr{MatStencil}, Ptr{$PetscScalar}, InsertMode),
               mat, m, idxm, n, idxn, v, addv,
              )


	return nothing
end 

"""
	MatSetStencil(petsclib::PetscLibType,mat::PetscMat, dim::PetscInt, dims::Vector{PetscInt}, starts::Vector{PetscInt}, dof::PetscInt) 
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
function MatSetStencil(petsclib::PetscLibType, mat::PetscMat, dim::PetscInt, dims::Vector{PetscInt}, starts::Vector{PetscInt}, dof::PetscInt) end

@for_petsc function MatSetStencil(petsclib::$UnionPetscLib, mat::PetscMat, dim::$PetscInt, dims::Vector{$PetscInt}, starts::Vector{$PetscInt}, dof::$PetscInt )

    @chk ccall(
               (:MatSetStencil, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, $PetscInt),
               mat, dim, dims, starts, dof,
              )


	return nothing
end 

"""
	MatSetValuesBlocked(petsclib::PetscLibType,mat::PetscMat, m::PetscInt, idxm::Vector{PetscInt}, n::PetscInt, idxn::Vector{PetscInt}, v::Vector{PetscScalar}, addv::InsertMode) 
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
function MatSetValuesBlocked(petsclib::PetscLibType, mat::PetscMat, m::PetscInt, idxm::Vector{PetscInt}, n::PetscInt, idxn::Vector{PetscInt}, v::Vector{PetscScalar}, addv::InsertMode) end

@for_petsc function MatSetValuesBlocked(petsclib::$UnionPetscLib, mat::PetscMat, m::$PetscInt, idxm::Vector{$PetscInt}, n::$PetscInt, idxn::Vector{$PetscInt}, v::Vector{$PetscScalar}, addv::InsertMode )

    @chk ccall(
               (:MatSetValuesBlocked, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}, InsertMode),
               mat, m, idxm, n, idxn, v, addv,
              )


	return nothing
end 

"""
	MatGetValues(petsclib::PetscLibType,mat::PetscMat, m::PetscInt, idxm::Vector{PetscInt}, n::PetscInt, idxn::Vector{PetscInt}, v::Vector{PetscScalar}) 
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
function MatGetValues(petsclib::PetscLibType, mat::PetscMat, m::PetscInt, idxm::Vector{PetscInt}, n::PetscInt, idxn::Vector{PetscInt}, v::Vector{PetscScalar}) end

@for_petsc function MatGetValues(petsclib::$UnionPetscLib, mat::PetscMat, m::$PetscInt, idxm::Vector{$PetscInt}, n::$PetscInt, idxn::Vector{$PetscInt}, v::Vector{$PetscScalar} )

    @chk ccall(
               (:MatGetValues, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}),
               mat, m, idxm, n, idxn, v,
              )


	return nothing
end 

"""
	MatGetValuesLocal(petsclib::PetscLibType,mat::PetscMat, nrow::PetscInt, irow::Vector{PetscInt}, ncol::PetscInt, icol::Vector{PetscInt}, y::Vector{PetscScalar}) 
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
function MatGetValuesLocal(petsclib::PetscLibType, mat::PetscMat, nrow::PetscInt, irow::Vector{PetscInt}, ncol::PetscInt, icol::Vector{PetscInt}, y::Vector{PetscScalar}) end

@for_petsc function MatGetValuesLocal(petsclib::$UnionPetscLib, mat::PetscMat, nrow::$PetscInt, irow::Vector{$PetscInt}, ncol::$PetscInt, icol::Vector{$PetscInt}, y::Vector{$PetscScalar} )

    @chk ccall(
               (:MatGetValuesLocal, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}),
               mat, nrow, irow, ncol, icol, y,
              )


	return nothing
end 

"""
	MatSetValuesBatch(petsclib::PetscLibType,mat::PetscMat, nb::PetscInt, bs::PetscInt, rows::Vector{PetscInt}, v::Vector{PetscScalar}) 
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
function MatSetValuesBatch(petsclib::PetscLibType, mat::PetscMat, nb::PetscInt, bs::PetscInt, rows::Vector{PetscInt}, v::Vector{PetscScalar}) end

@for_petsc function MatSetValuesBatch(petsclib::$UnionPetscLib, mat::PetscMat, nb::$PetscInt, bs::$PetscInt, rows::Vector{$PetscInt}, v::Vector{$PetscScalar} )

    @chk ccall(
               (:MatSetValuesBatch, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}),
               mat, nb, bs, rows, v,
              )


	return nothing
end 

"""
	MatSetLocalToGlobalMapping(petsclib::PetscLibType,x::PetscMat, rmapping::ISLocalToGlobalMapping, cmapping::ISLocalToGlobalMapping) 
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
function MatSetLocalToGlobalMapping(petsclib::PetscLibType, x::PetscMat, rmapping::ISLocalToGlobalMapping, cmapping::ISLocalToGlobalMapping) end

@for_petsc function MatSetLocalToGlobalMapping(petsclib::$UnionPetscLib, x::PetscMat, rmapping::ISLocalToGlobalMapping, cmapping::ISLocalToGlobalMapping )

    @chk ccall(
               (:MatSetLocalToGlobalMapping, $petsc_library),
               PetscErrorCode,
               (CMat, ISLocalToGlobalMapping, ISLocalToGlobalMapping),
               x, rmapping, cmapping,
              )


	return nothing
end 

"""
	MatGetLocalToGlobalMapping(petsclib::PetscLibType,A::PetscMat, rmapping::ISLocalToGlobalMapping, cmapping::ISLocalToGlobalMapping) 
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
function MatGetLocalToGlobalMapping(petsclib::PetscLibType, A::PetscMat, rmapping::ISLocalToGlobalMapping, cmapping::ISLocalToGlobalMapping) end

@for_petsc function MatGetLocalToGlobalMapping(petsclib::$UnionPetscLib, A::PetscMat, rmapping::ISLocalToGlobalMapping, cmapping::ISLocalToGlobalMapping )

    @chk ccall(
               (:MatGetLocalToGlobalMapping, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{ISLocalToGlobalMapping}, Ptr{ISLocalToGlobalMapping}),
               A, rmapping, cmapping,
              )


	return nothing
end 

"""
	MatSetLayouts(petsclib::PetscLibType,A::PetscMat, rmap::PetscLayout, cmap::PetscLayout) 
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
function MatSetLayouts(petsclib::PetscLibType, A::PetscMat, rmap::PetscLayout, cmap::PetscLayout) end

@for_petsc function MatSetLayouts(petsclib::$UnionPetscLib, A::PetscMat, rmap::PetscLayout, cmap::PetscLayout )

    @chk ccall(
               (:MatSetLayouts, $petsc_library),
               PetscErrorCode,
               (CMat, PetscLayout, PetscLayout),
               A, rmap, cmap,
              )


	return nothing
end 

"""
	MatGetLayouts(petsclib::PetscLibType,A::PetscMat, rmap::PetscLayout, cmap::PetscLayout) 
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
function MatGetLayouts(petsclib::PetscLibType, A::PetscMat, rmap::PetscLayout, cmap::PetscLayout) end

@for_petsc function MatGetLayouts(petsclib::$UnionPetscLib, A::PetscMat, rmap::PetscLayout, cmap::PetscLayout )

    @chk ccall(
               (:MatGetLayouts, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscLayout}, Ptr{PetscLayout}),
               A, rmap, cmap,
              )


	return nothing
end 

"""
	MatSetValuesLocal(petsclib::PetscLibType,mat::PetscMat, nrow::PetscInt, irow::Vector{PetscInt}, ncol::PetscInt, icol::Vector{PetscInt}, y::Vector{PetscScalar}, addv::InsertMode) 
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
function MatSetValuesLocal(petsclib::PetscLibType, mat::PetscMat, nrow::PetscInt, irow::Vector{PetscInt}, ncol::PetscInt, icol::Vector{PetscInt}, y::Vector{PetscScalar}, addv::InsertMode) end

@for_petsc function MatSetValuesLocal(petsclib::$UnionPetscLib, mat::PetscMat, nrow::$PetscInt, irow::Vector{$PetscInt}, ncol::$PetscInt, icol::Vector{$PetscInt}, y::Vector{$PetscScalar}, addv::InsertMode )

    @chk ccall(
               (:MatSetValuesLocal, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}, InsertMode),
               mat, nrow, irow, ncol, icol, y, addv,
              )


	return nothing
end 

"""
	MatSetValuesBlockedLocal(petsclib::PetscLibType,mat::PetscMat, nrow::PetscInt, irow::Vector{PetscInt}, ncol::PetscInt, icol::Vector{PetscInt}, y::Vector{PetscScalar}, addv::InsertMode) 
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
function MatSetValuesBlockedLocal(petsclib::PetscLibType, mat::PetscMat, nrow::PetscInt, irow::Vector{PetscInt}, ncol::PetscInt, icol::Vector{PetscInt}, y::Vector{PetscScalar}, addv::InsertMode) end

@for_petsc function MatSetValuesBlockedLocal(petsclib::$UnionPetscLib, mat::PetscMat, nrow::$PetscInt, irow::Vector{$PetscInt}, ncol::$PetscInt, icol::Vector{$PetscInt}, y::Vector{$PetscScalar}, addv::InsertMode )

    @chk ccall(
               (:MatSetValuesBlockedLocal, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}, InsertMode),
               mat, nrow, irow, ncol, icol, y, addv,
              )


	return nothing
end 

"""
	MatMultDiagonalBlock(petsclib::PetscLibType,mat::PetscMat, x::PetscVec, y::PetscVec) 
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
function MatMultDiagonalBlock(petsclib::PetscLibType, mat::PetscMat, x::PetscVec, y::PetscVec) end

@for_petsc function MatMultDiagonalBlock(petsclib::$UnionPetscLib, mat::PetscMat, x::PetscVec, y::PetscVec )

    @chk ccall(
               (:MatMultDiagonalBlock, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               mat, x, y,
              )


	return nothing
end 

"""
	MatMult(petsclib::PetscLibType,mat::PetscMat, x::PetscVec, y::PetscVec) 
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
function MatMult(petsclib::PetscLibType, mat::AbstractPetscMat, x::PetscVec, y::PetscVec) end

@for_petsc function MatMult(petsclib::$UnionPetscLib, mat::AbstractPetscMat, x::PetscVec, y::PetscVec )

    @chk ccall(
               (:MatMult, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               mat, x, y,
              )


	return nothing
end 

"""
	MatMultTranspose(petsclib::PetscLibType,mat::PetscMat, x::PetscVec, y::PetscVec) 
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
function MatMultTranspose(petsclib::PetscLibType, mat::PetscMat, x::PetscVec, y::PetscVec) end

@for_petsc function MatMultTranspose(petsclib::$UnionPetscLib, mat::PetscMat, x::PetscVec, y::PetscVec )

    @chk ccall(
               (:MatMultTranspose, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               mat, x, y,
              )


	return nothing
end 

"""
	MatMultHermitianTranspose(petsclib::PetscLibType,mat::PetscMat, x::PetscVec, y::PetscVec) 
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
function MatMultHermitianTranspose(petsclib::PetscLibType, mat::PetscMat, x::PetscVec, y::PetscVec) end

@for_petsc function MatMultHermitianTranspose(petsclib::$UnionPetscLib, mat::PetscMat, x::PetscVec, y::PetscVec )

    @chk ccall(
               (:MatMultHermitianTranspose, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               mat, x, y,
              )


	return nothing
end 

"""
	MatMultAdd(petsclib::PetscLibType,mat::PetscMat, v1::PetscVec, v2::PetscVec, v3::PetscVec) 
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
function MatMultAdd(petsclib::PetscLibType, mat::PetscMat, v1::PetscVec, v2::PetscVec, v3::PetscVec) end

@for_petsc function MatMultAdd(petsclib::$UnionPetscLib, mat::PetscMat, v1::PetscVec, v2::PetscVec, v3::PetscVec )

    @chk ccall(
               (:MatMultAdd, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec, CVec),
               mat, v1, v2, v3,
              )


	return nothing
end 

"""
	MatMultTransposeAdd(petsclib::PetscLibType,mat::PetscMat, v1::PetscVec, v2::PetscVec, v3::PetscVec) 
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
function MatMultTransposeAdd(petsclib::PetscLibType, mat::PetscMat, v1::PetscVec, v2::PetscVec, v3::PetscVec) end

@for_petsc function MatMultTransposeAdd(petsclib::$UnionPetscLib, mat::PetscMat, v1::PetscVec, v2::PetscVec, v3::PetscVec )

    @chk ccall(
               (:MatMultTransposeAdd, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec, CVec),
               mat, v1, v2, v3,
              )


	return nothing
end 

"""
	MatMultHermitianTransposeAdd(petsclib::PetscLibType,mat::PetscMat, v1::PetscVec, v2::PetscVec, v3::PetscVec) 
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
function MatMultHermitianTransposeAdd(petsclib::PetscLibType, mat::PetscMat, v1::PetscVec, v2::PetscVec, v3::PetscVec) end

@for_petsc function MatMultHermitianTransposeAdd(petsclib::$UnionPetscLib, mat::PetscMat, v1::PetscVec, v2::PetscVec, v3::PetscVec )

    @chk ccall(
               (:MatMultHermitianTransposeAdd, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec, CVec),
               mat, v1, v2, v3,
              )


	return nothing
end 

"""
	t::MatFactorType = MatGetFactorType(petsclib::PetscLibType,mat::PetscMat) 
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
function MatGetFactorType(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatGetFactorType(petsclib::$UnionPetscLib, mat::PetscMat )
	t_ = Ref{MatFactorType}()

    @chk ccall(
               (:MatGetFactorType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatFactorType}),
               mat, t_,
              )

	t = unsafe_string(t_[])

	return t
end 

"""
	MatSetFactorType(petsclib::PetscLibType,mat::PetscMat, t::MatFactorType) 
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
function MatSetFactorType(petsclib::PetscLibType, mat::PetscMat, t::MatFactorType) end

@for_petsc function MatSetFactorType(petsclib::$UnionPetscLib, mat::PetscMat, t::MatFactorType )

    @chk ccall(
               (:MatSetFactorType, $petsc_library),
               PetscErrorCode,
               (CMat, MatFactorType),
               mat, t,
              )


	return nothing
end 

"""
	MatGetInfo(petsclib::PetscLibType,mat::PetscMat, flag::MatInfoType, info::MatInfo) 
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
function MatGetInfo(petsclib::PetscLibType, mat::PetscMat, flag::MatInfoType, info::MatInfo) end

@for_petsc function MatGetInfo(petsclib::$UnionPetscLib, mat::PetscMat, flag::MatInfoType, info::MatInfo )

    @chk ccall(
               (:MatGetInfo, $petsc_library),
               PetscErrorCode,
               (CMat, MatInfoType, Ptr{MatInfo}),
               mat, flag, info,
              )


	return nothing
end 

"""
	MatLUFactor(petsclib::PetscLibType,mat::PetscMat, row::IS, col::IS, info::MatFactorInfo) 
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
function MatLUFactor(petsclib::PetscLibType, mat::PetscMat, row::IS, col::IS, info::MatFactorInfo) end

@for_petsc function MatLUFactor(petsclib::$UnionPetscLib, mat::PetscMat, row::IS, col::IS, info::MatFactorInfo )

    @chk ccall(
               (:MatLUFactor, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, CIS, Ptr{MatFactorInfo}),
               mat, row, col, info,
              )


	return nothing
end 

"""
	MatILUFactor(petsclib::PetscLibType,mat::PetscMat, row::IS, col::IS, info::MatFactorInfo) 
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
function MatILUFactor(petsclib::PetscLibType, mat::PetscMat, row::IS, col::IS, info::MatFactorInfo) end

@for_petsc function MatILUFactor(petsclib::$UnionPetscLib, mat::PetscMat, row::IS, col::IS, info::MatFactorInfo )

    @chk ccall(
               (:MatILUFactor, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, CIS, Ptr{MatFactorInfo}),
               mat, row, col, info,
              )


	return nothing
end 

"""
	MatLUFactorSymbolic(petsclib::PetscLibType,fact::PetscMat, mat::PetscMat, row::IS, col::IS, info::MatFactorInfo) 
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
function MatLUFactorSymbolic(petsclib::PetscLibType, fact::PetscMat, mat::PetscMat, row::IS, col::IS, info::MatFactorInfo) end

@for_petsc function MatLUFactorSymbolic(petsclib::$UnionPetscLib, fact::PetscMat, mat::PetscMat, row::IS, col::IS, info::MatFactorInfo )

    @chk ccall(
               (:MatLUFactorSymbolic, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CIS, CIS, Ptr{MatFactorInfo}),
               fact, mat, row, col, info,
              )


	return nothing
end 

"""
	MatLUFactorNumeric(petsclib::PetscLibType,fact::PetscMat, mat::PetscMat, info::MatFactorInfo) 
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
function MatLUFactorNumeric(petsclib::PetscLibType, fact::PetscMat, mat::PetscMat, info::MatFactorInfo) end

@for_petsc function MatLUFactorNumeric(petsclib::$UnionPetscLib, fact::PetscMat, mat::PetscMat, info::MatFactorInfo )

    @chk ccall(
               (:MatLUFactorNumeric, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, Ptr{MatFactorInfo}),
               fact, mat, info,
              )


	return nothing
end 

"""
	MatCholeskyFactor(petsclib::PetscLibType,mat::PetscMat, perm::IS, info::MatFactorInfo) 
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
function MatCholeskyFactor(petsclib::PetscLibType, mat::PetscMat, perm::IS, info::MatFactorInfo) end

@for_petsc function MatCholeskyFactor(petsclib::$UnionPetscLib, mat::PetscMat, perm::IS, info::MatFactorInfo )

    @chk ccall(
               (:MatCholeskyFactor, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, Ptr{MatFactorInfo}),
               mat, perm, info,
              )


	return nothing
end 

"""
	MatCholeskyFactorSymbolic(petsclib::PetscLibType,fact::PetscMat, mat::PetscMat, perm::IS, info::MatFactorInfo) 
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
function MatCholeskyFactorSymbolic(petsclib::PetscLibType, fact::PetscMat, mat::PetscMat, perm::IS, info::MatFactorInfo) end

@for_petsc function MatCholeskyFactorSymbolic(petsclib::$UnionPetscLib, fact::PetscMat, mat::PetscMat, perm::IS, info::MatFactorInfo )

    @chk ccall(
               (:MatCholeskyFactorSymbolic, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CIS, Ptr{MatFactorInfo}),
               fact, mat, perm, info,
              )


	return nothing
end 

"""
	MatCholeskyFactorNumeric(petsclib::PetscLibType,fact::PetscMat, mat::PetscMat, info::MatFactorInfo) 
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
function MatCholeskyFactorNumeric(petsclib::PetscLibType, fact::PetscMat, mat::PetscMat, info::MatFactorInfo) end

@for_petsc function MatCholeskyFactorNumeric(petsclib::$UnionPetscLib, fact::PetscMat, mat::PetscMat, info::MatFactorInfo )

    @chk ccall(
               (:MatCholeskyFactorNumeric, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, Ptr{MatFactorInfo}),
               fact, mat, info,
              )


	return nothing
end 

"""
	MatQRFactor(petsclib::PetscLibType,mat::PetscMat, col::IS, info::MatFactorInfo) 
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
function MatQRFactor(petsclib::PetscLibType, mat::PetscMat, col::IS, info::MatFactorInfo) end

@for_petsc function MatQRFactor(petsclib::$UnionPetscLib, mat::PetscMat, col::IS, info::MatFactorInfo )

    @chk ccall(
               (:MatQRFactor, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, Ptr{MatFactorInfo}),
               mat, col, info,
              )


	return nothing
end 

"""
	MatQRFactorSymbolic(petsclib::PetscLibType,fact::PetscMat, mat::PetscMat, col::IS, info::MatFactorInfo) 
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
function MatQRFactorSymbolic(petsclib::PetscLibType, fact::PetscMat, mat::PetscMat, col::IS, info::MatFactorInfo) end

@for_petsc function MatQRFactorSymbolic(petsclib::$UnionPetscLib, fact::PetscMat, mat::PetscMat, col::IS, info::MatFactorInfo )

    @chk ccall(
               (:MatQRFactorSymbolic, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CIS, Ptr{MatFactorInfo}),
               fact, mat, col, info,
              )


	return nothing
end 

"""
	MatQRFactorNumeric(petsclib::PetscLibType,fact::PetscMat, mat::PetscMat, info::MatFactorInfo) 
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
function MatQRFactorNumeric(petsclib::PetscLibType, fact::PetscMat, mat::PetscMat, info::MatFactorInfo) end

@for_petsc function MatQRFactorNumeric(petsclib::$UnionPetscLib, fact::PetscMat, mat::PetscMat, info::MatFactorInfo )

    @chk ccall(
               (:MatQRFactorNumeric, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, Ptr{MatFactorInfo}),
               fact, mat, info,
              )


	return nothing
end 

"""
	MatSolve(petsclib::PetscLibType,mat::PetscMat, b::PetscVec, x::PetscVec) 
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
function MatSolve(petsclib::PetscLibType, mat::PetscMat, b::PetscVec, x::PetscVec) end

@for_petsc function MatSolve(petsclib::$UnionPetscLib, mat::PetscMat, b::PetscVec, x::PetscVec )

    @chk ccall(
               (:MatSolve, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               mat, b, x,
              )


	return nothing
end 

"""
	MatMatSolve(petsclib::PetscLibType,A::PetscMat, B::PetscMat, X::PetscMat) 
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
function MatMatSolve(petsclib::PetscLibType, A::PetscMat, B::PetscMat, X::PetscMat) end

@for_petsc function MatMatSolve(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat, X::PetscMat )

    @chk ccall(
               (:MatMatSolve, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat),
               A, B, X,
              )


	return nothing
end 

"""
	MatMatSolveTranspose(petsclib::PetscLibType,A::PetscMat, B::PetscMat, X::PetscMat) 
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
function MatMatSolveTranspose(petsclib::PetscLibType, A::PetscMat, B::PetscMat, X::PetscMat) end

@for_petsc function MatMatSolveTranspose(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat, X::PetscMat )

    @chk ccall(
               (:MatMatSolveTranspose, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat),
               A, B, X,
              )


	return nothing
end 

"""
	MatMatTransposeSolve(petsclib::PetscLibType,A::PetscMat, Bt::PetscMat, X::PetscMat) 
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
function MatMatTransposeSolve(petsclib::PetscLibType, A::PetscMat, Bt::PetscMat, X::PetscMat) end

@for_petsc function MatMatTransposeSolve(petsclib::$UnionPetscLib, A::PetscMat, Bt::PetscMat, X::PetscMat )

    @chk ccall(
               (:MatMatTransposeSolve, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat),
               A, Bt, X,
              )


	return nothing
end 

"""
	MatForwardSolve(petsclib::PetscLibType,mat::PetscMat, b::PetscVec, x::PetscVec) 
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
function MatForwardSolve(petsclib::PetscLibType, mat::PetscMat, b::PetscVec, x::PetscVec) end

@for_petsc function MatForwardSolve(petsclib::$UnionPetscLib, mat::PetscMat, b::PetscVec, x::PetscVec )

    @chk ccall(
               (:MatForwardSolve, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               mat, b, x,
              )


	return nothing
end 

"""
	MatBackwardSolve(petsclib::PetscLibType,mat::PetscMat, b::PetscVec, x::PetscVec) 
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
function MatBackwardSolve(petsclib::PetscLibType, mat::PetscMat, b::PetscVec, x::PetscVec) end

@for_petsc function MatBackwardSolve(petsclib::$UnionPetscLib, mat::PetscMat, b::PetscVec, x::PetscVec )

    @chk ccall(
               (:MatBackwardSolve, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               mat, b, x,
              )


	return nothing
end 

"""
	MatSolveAdd(petsclib::PetscLibType,mat::PetscMat, b::PetscVec, y::PetscVec, x::PetscVec) 
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
function MatSolveAdd(petsclib::PetscLibType, mat::PetscMat, b::PetscVec, y::PetscVec, x::PetscVec) end

@for_petsc function MatSolveAdd(petsclib::$UnionPetscLib, mat::PetscMat, b::PetscVec, y::PetscVec, x::PetscVec )

    @chk ccall(
               (:MatSolveAdd, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec, CVec),
               mat, b, y, x,
              )


	return nothing
end 

"""
	MatSolveTranspose(petsclib::PetscLibType,mat::PetscMat, b::PetscVec, x::PetscVec) 
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
function MatSolveTranspose(petsclib::PetscLibType, mat::PetscMat, b::PetscVec, x::PetscVec) end

@for_petsc function MatSolveTranspose(petsclib::$UnionPetscLib, mat::PetscMat, b::PetscVec, x::PetscVec )

    @chk ccall(
               (:MatSolveTranspose, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               mat, b, x,
              )


	return nothing
end 

"""
	MatSolveTransposeAdd(petsclib::PetscLibType,mat::PetscMat, b::PetscVec, y::PetscVec, x::PetscVec) 
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
function MatSolveTransposeAdd(petsclib::PetscLibType, mat::PetscMat, b::PetscVec, y::PetscVec, x::PetscVec) end

@for_petsc function MatSolveTransposeAdd(petsclib::$UnionPetscLib, mat::PetscMat, b::PetscVec, y::PetscVec, x::PetscVec )

    @chk ccall(
               (:MatSolveTransposeAdd, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec, CVec),
               mat, b, y, x,
              )


	return nothing
end 

"""
	MatSOR(petsclib::PetscLibType,mat::PetscMat, b::PetscVec, omega::PetscReal, flag::MatSORType, shift::PetscReal, its::PetscInt, lits::PetscInt, x::PetscVec) 
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
function MatSOR(petsclib::PetscLibType, mat::PetscMat, b::PetscVec, omega::PetscReal, flag::MatSORType, shift::PetscReal, its::PetscInt, lits::PetscInt, x::PetscVec) end

@for_petsc function MatSOR(petsclib::$UnionPetscLib, mat::PetscMat, b::PetscVec, omega::$PetscReal, flag::MatSORType, shift::$PetscReal, its::$PetscInt, lits::$PetscInt, x::PetscVec )

    @chk ccall(
               (:MatSOR, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, $PetscReal, MatSORType, $PetscReal, $PetscInt, $PetscInt, CVec),
               mat, b, omega, flag, shift, its, lits, x,
              )


	return nothing
end 

"""
	MatCopy(petsclib::PetscLibType,A::PetscMat, B::PetscMat, str::MatStructure) 
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
function MatCopy(petsclib::PetscLibType, A::PetscMat, B::PetscMat, str::MatStructure) end

@for_petsc function MatCopy(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat, str::MatStructure )

    @chk ccall(
               (:MatCopy, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, MatStructure),
               A, B, str,
              )


	return nothing
end 

"""
	MatConvert(petsclib::PetscLibType,mat::PetscMat, newtype::MatType, reuse::MatReuse, M::PetscMat) 
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
function MatConvert(petsclib::PetscLibType, mat::PetscMat, newtype::MatType, reuse::MatReuse, M::PetscMat) end

@for_petsc function MatConvert(petsclib::$UnionPetscLib, mat::PetscMat, newtype::MatType, reuse::MatReuse, M::PetscMat )
	M_ = Ref(M.ptr)

    @chk ccall(
               (:MatConvert, $petsc_library),
               PetscErrorCode,
               (CMat, MatType, MatReuse, Ptr{CMat}),
               mat, newtype, reuse, M_,
              )

	M.ptr = C_NULL

	return nothing
end 

"""
	type::MatSolverType = MatFactorGetSolverType(petsclib::PetscLibType,mat::PetscMat) 
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
function MatFactorGetSolverType(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatFactorGetSolverType(petsclib::$UnionPetscLib, mat::PetscMat )
	type_ = Ref{MatSolverType}()

    @chk ccall(
               (:MatFactorGetSolverType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatSolverType}),
               mat, type_,
              )

	type = unsafe_string(type_[])

	return type
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
	flg::PetscBool = MatFactorGetCanUseOrdering(petsclib::PetscLibType,mat::PetscMat) 
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
function MatFactorGetCanUseOrdering(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatFactorGetCanUseOrdering(petsclib::$UnionPetscLib, mat::PetscMat )
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
	MatFactorGetPreferredOrdering(petsclib::PetscLibType,mat::PetscMat, ftype::MatFactorType, otype::MatOrderingType) 
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
function MatFactorGetPreferredOrdering(petsclib::PetscLibType, mat::PetscMat, ftype::MatFactorType, otype::MatOrderingType) end

@for_petsc function MatFactorGetPreferredOrdering(petsclib::$UnionPetscLib, mat::PetscMat, ftype::MatFactorType, otype::MatOrderingType )

    @chk ccall(
               (:MatFactorGetPreferredOrdering, $petsc_library),
               PetscErrorCode,
               (CMat, MatFactorType, Ptr{MatOrderingType}),
               mat, ftype, otype,
              )


	return nothing
end 

"""
	MatGetFactor(petsclib::PetscLibType,mat::PetscMat, type::MatSolverType, ftype::MatFactorType, f::PetscMat) 
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
function MatGetFactor(petsclib::PetscLibType, mat::PetscMat, type::MatSolverType, ftype::MatFactorType, f::PetscMat) end

@for_petsc function MatGetFactor(petsclib::$UnionPetscLib, mat::PetscMat, type::MatSolverType, ftype::MatFactorType, f::PetscMat )
	f_ = Ref(f.ptr)

    @chk ccall(
               (:MatGetFactor, $petsc_library),
               PetscErrorCode,
               (CMat, MatSolverType, MatFactorType, Ptr{CMat}),
               mat, type, ftype, f_,
              )

	f.ptr = C_NULL

	return nothing
end 

"""
	flg::PetscBool = MatGetFactorAvailable(petsclib::PetscLibType,mat::PetscMat, type::MatSolverType, ftype::MatFactorType) 
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
function MatGetFactorAvailable(petsclib::PetscLibType, mat::PetscMat, type::MatSolverType, ftype::MatFactorType) end

@for_petsc function MatGetFactorAvailable(petsclib::$UnionPetscLib, mat::PetscMat, type::MatSolverType, ftype::MatFactorType )
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
	MatDuplicate(petsclib::PetscLibType,mat::PetscMat, op::MatDuplicateOption, M::PetscMat) 
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
function MatDuplicate(petsclib::PetscLibType, mat::PetscMat, op::MatDuplicateOption, M::PetscMat) end

@for_petsc function MatDuplicate(petsclib::$UnionPetscLib, mat::PetscMat, op::MatDuplicateOption, M::PetscMat )
	M_ = Ref(M.ptr)

    @chk ccall(
               (:MatDuplicate, $petsc_library),
               PetscErrorCode,
               (CMat, MatDuplicateOption, Ptr{CMat}),
               mat, op, M_,
              )

	M.ptr = C_NULL

	return nothing
end 

"""
	MatGetDiagonal(petsclib::PetscLibType,mat::PetscMat, v::PetscVec) 
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
function MatGetDiagonal(petsclib::PetscLibType, mat::PetscMat, v::PetscVec) end

@for_petsc function MatGetDiagonal(petsclib::$UnionPetscLib, mat::PetscMat, v::PetscVec )

    @chk ccall(
               (:MatGetDiagonal, $petsc_library),
               PetscErrorCode,
               (CMat, CVec),
               mat, v,
              )


	return nothing
end 

"""
	MatGetRowMin(petsclib::PetscLibType,mat::PetscMat, v::PetscVec, idx::Vector{PetscInt}) 
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
function MatGetRowMin(petsclib::PetscLibType, mat::PetscMat, v::PetscVec, idx::Vector{PetscInt}) end

@for_petsc function MatGetRowMin(petsclib::$UnionPetscLib, mat::PetscMat, v::PetscVec, idx::Vector{$PetscInt} )

    @chk ccall(
               (:MatGetRowMin, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, Ptr{$PetscInt}),
               mat, v, idx,
              )


	return nothing
end 

"""
	MatGetRowMinAbs(petsclib::PetscLibType,mat::PetscMat, v::PetscVec, idx::Vector{PetscInt}) 
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
function MatGetRowMinAbs(petsclib::PetscLibType, mat::PetscMat, v::PetscVec, idx::Vector{PetscInt}) end

@for_petsc function MatGetRowMinAbs(petsclib::$UnionPetscLib, mat::PetscMat, v::PetscVec, idx::Vector{$PetscInt} )

    @chk ccall(
               (:MatGetRowMinAbs, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, Ptr{$PetscInt}),
               mat, v, idx,
              )


	return nothing
end 

"""
	MatGetRowMax(petsclib::PetscLibType,mat::PetscMat, v::PetscVec, idx::Vector{PetscInt}) 
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
function MatGetRowMax(petsclib::PetscLibType, mat::PetscMat, v::PetscVec, idx::Vector{PetscInt}) end

@for_petsc function MatGetRowMax(petsclib::$UnionPetscLib, mat::PetscMat, v::PetscVec, idx::Vector{$PetscInt} )

    @chk ccall(
               (:MatGetRowMax, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, Ptr{$PetscInt}),
               mat, v, idx,
              )


	return nothing
end 

"""
	MatGetRowMaxAbs(petsclib::PetscLibType,mat::PetscMat, v::PetscVec, idx::Vector{PetscInt}) 
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
function MatGetRowMaxAbs(petsclib::PetscLibType, mat::PetscMat, v::PetscVec, idx::Vector{PetscInt}) end

@for_petsc function MatGetRowMaxAbs(petsclib::$UnionPetscLib, mat::PetscMat, v::PetscVec, idx::Vector{$PetscInt} )

    @chk ccall(
               (:MatGetRowMaxAbs, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, Ptr{$PetscInt}),
               mat, v, idx,
              )


	return nothing
end 

"""
	MatGetRowSumAbs(petsclib::PetscLibType,mat::PetscMat, v::PetscVec) 
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
function MatGetRowSumAbs(petsclib::PetscLibType, mat::PetscMat, v::PetscVec) end

@for_petsc function MatGetRowSumAbs(petsclib::$UnionPetscLib, mat::PetscMat, v::PetscVec )

    @chk ccall(
               (:MatGetRowSumAbs, $petsc_library),
               PetscErrorCode,
               (CMat, CVec),
               mat, v,
              )


	return nothing
end 

"""
	MatGetRowSum(petsclib::PetscLibType,mat::PetscMat, v::PetscVec) 
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
function MatGetRowSum(petsclib::PetscLibType, mat::PetscMat, v::PetscVec) end

@for_petsc function MatGetRowSum(petsclib::$UnionPetscLib, mat::PetscMat, v::PetscVec )

    @chk ccall(
               (:MatGetRowSum, $petsc_library),
               PetscErrorCode,
               (CMat, CVec),
               mat, v,
              )


	return nothing
end 

"""
	MatTransposeSetPrecursor(petsclib::PetscLibType,mat::PetscMat, B::PetscMat) 
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
function MatTransposeSetPrecursor(petsclib::PetscLibType, mat::PetscMat, B::PetscMat) end

@for_petsc function MatTransposeSetPrecursor(petsclib::$UnionPetscLib, mat::PetscMat, B::PetscMat )

    @chk ccall(
               (:MatTransposeSetPrecursor, $petsc_library),
               PetscErrorCode,
               (CMat, CMat),
               mat, B,
              )


	return nothing
end 

"""
	MatTranspose(petsclib::PetscLibType,mat::PetscMat, reuse::MatReuse, B::PetscMat) 
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
function MatTranspose(petsclib::PetscLibType, mat::PetscMat, reuse::MatReuse, B::PetscMat) end

@for_petsc function MatTranspose(petsclib::$UnionPetscLib, mat::PetscMat, reuse::MatReuse, B::PetscMat )
	B_ = Ref(B.ptr)

    @chk ccall(
               (:MatTranspose, $petsc_library),
               PetscErrorCode,
               (CMat, MatReuse, Ptr{CMat}),
               mat, reuse, B_,
              )

	B.ptr = C_NULL

	return nothing
end 

"""
	MatTransposeSymbolic(petsclib::PetscLibType,A::PetscMat, B::PetscMat) 
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
function MatTransposeSymbolic(petsclib::PetscLibType, A::PetscMat, B::PetscMat) end

@for_petsc function MatTransposeSymbolic(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat )
	B_ = Ref(B.ptr)

    @chk ccall(
               (:MatTransposeSymbolic, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, B_,
              )

	B.ptr = C_NULL

	return nothing
end 

"""
	flg::PetscBool = MatIsTranspose(petsclib::PetscLibType,A::PetscMat, B::PetscMat, tol::PetscReal) 
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
function MatIsTranspose(petsclib::PetscLibType, A::PetscMat, B::PetscMat, tol::PetscReal) end

@for_petsc function MatIsTranspose(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat, tol::$PetscReal )
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
	MatHermitianTranspose(petsclib::PetscLibType,mat::PetscMat, reuse::MatReuse, B::PetscMat) 
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
function MatHermitianTranspose(petsclib::PetscLibType, mat::PetscMat, reuse::MatReuse, B::PetscMat) end

@for_petsc function MatHermitianTranspose(petsclib::$UnionPetscLib, mat::PetscMat, reuse::MatReuse, B::PetscMat )
	B_ = Ref(B.ptr)

    @chk ccall(
               (:MatHermitianTranspose, $petsc_library),
               PetscErrorCode,
               (CMat, MatReuse, Ptr{CMat}),
               mat, reuse, B_,
              )

	B.ptr = C_NULL

	return nothing
end 

"""
	flg::PetscBool = MatIsHermitianTranspose(petsclib::PetscLibType,A::PetscMat, B::PetscMat, tol::PetscReal) 
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
function MatIsHermitianTranspose(petsclib::PetscLibType, A::PetscMat, B::PetscMat, tol::PetscReal) end

@for_petsc function MatIsHermitianTranspose(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat, tol::$PetscReal )
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
	MatPermute(petsclib::PetscLibType,mat::PetscMat, row::IS, col::IS, B::PetscMat) 
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
function MatPermute(petsclib::PetscLibType, mat::PetscMat, row::IS, col::IS, B::PetscMat) end

@for_petsc function MatPermute(petsclib::$UnionPetscLib, mat::PetscMat, row::IS, col::IS, B::PetscMat )
	B_ = Ref(B.ptr)

    @chk ccall(
               (:MatPermute, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, CIS, Ptr{CMat}),
               mat, row, col, B_,
              )

	B.ptr = C_NULL

	return nothing
end 

"""
	flg::PetscBool = MatEqual(petsclib::PetscLibType,A::PetscMat, B::PetscMat) 
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
function MatEqual(petsclib::PetscLibType, A::PetscMat, B::PetscMat) end

@for_petsc function MatEqual(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat )
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
	MatDiagonalScale(petsclib::PetscLibType,mat::PetscMat, l::PetscVec, r::PetscVec) 
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
function MatDiagonalScale(petsclib::PetscLibType, mat::PetscMat, l::Union{Ptr,PetscVec}, r::Union{Ptr,PetscVec}) end

@for_petsc function MatDiagonalScale(petsclib::$UnionPetscLib, mat::PetscMat, l::Union{Ptr,PetscVec}, r::Union{Ptr,PetscVec})

    @chk ccall(
               (:MatDiagonalScale, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               mat, l, r,
              )


	return nothing
end 

"""
	MatScale(petsclib::PetscLibType,mat::PetscMat, a::PetscScalar) 
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
function MatScale(petsclib::PetscLibType, mat::PetscMat, a::PetscScalar) end

@for_petsc function MatScale(petsclib::$UnionPetscLib, mat::PetscMat, a::$PetscScalar )

    @chk ccall(
               (:MatScale, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscScalar),
               mat, a,
              )


	return nothing
end 

"""
	nrm::PetscReal = MatNorm(petsclib::PetscLibType,mat::PetscMat, type::NormType) 
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
function MatNorm(petsclib::PetscLibType, mat::PetscMat, type::NormType) end

@for_petsc function MatNorm(petsclib::$UnionPetscLib, mat::PetscMat, type::NormType )
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
	MatAssemblyBegin(petsclib::PetscLibType,mat::PetscMat, type::MatAssemblyType) 
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
	assembled::PetscBool = MatAssembled(petsclib::PetscLibType,mat::PetscMat) 
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
function MatAssembled(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatAssembled(petsclib::$UnionPetscLib, mat::PetscMat )
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
	MatAssemblyEnd(petsclib::PetscLibType,mat::PetscMat, type::MatAssemblyType) 
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
	MatSetOption(petsclib::PetscLibType,mat::PetscMat, op::MatOption, flg::PetscBool) 
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
function MatSetOption(petsclib::PetscLibType, mat::PetscMat, op::MatOption, flg::PetscBool) end

@for_petsc function MatSetOption(petsclib::$UnionPetscLib, mat::PetscMat, op::MatOption, flg::PetscBool )

    @chk ccall(
               (:MatSetOption, $petsc_library),
               PetscErrorCode,
               (CMat, MatOption, PetscBool),
               mat, op, flg,
              )


	return nothing
end 

"""
	flg::PetscBool = MatGetOption(petsclib::PetscLibType,mat::PetscMat, op::MatOption) 
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
function MatGetOption(petsclib::PetscLibType, mat::PetscMat, op::MatOption) end

@for_petsc function MatGetOption(petsclib::$UnionPetscLib, mat::PetscMat, op::MatOption )
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
	MatZeroEntries(petsclib::PetscLibType,mat::PetscMat) 
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
function MatZeroEntries(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatZeroEntries(petsclib::$UnionPetscLib, mat::PetscMat )

    @chk ccall(
               (:MatZeroEntries, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	MatZeroRowsColumns(petsclib::PetscLibType,mat::PetscMat, numRows::PetscInt, rows::Vector{PetscInt}, diag::PetscScalar, x::PetscVec, b::PetscVec) 
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
function MatZeroRowsColumns(petsclib::PetscLibType, mat::PetscMat, numRows::PetscInt, rows::Vector{PetscInt}, diag::PetscScalar, x::PetscVec, b::PetscVec) end

@for_petsc function MatZeroRowsColumns(petsclib::$UnionPetscLib, mat::PetscMat, numRows::$PetscInt, rows::Vector{$PetscInt}, diag::$PetscScalar, x::PetscVec, b::PetscVec )

    @chk ccall(
               (:MatZeroRowsColumns, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscScalar, CVec, CVec),
               mat, numRows, rows, diag, x, b,
              )


	return nothing
end 

"""
	MatZeroRowsColumnsIS(petsclib::PetscLibType,mat::PetscMat, is::IS, diag::PetscScalar, x::PetscVec, b::PetscVec) 
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
function MatZeroRowsColumnsIS(petsclib::PetscLibType, mat::PetscMat, is::IS, diag::PetscScalar, x::PetscVec, b::PetscVec) end

@for_petsc function MatZeroRowsColumnsIS(petsclib::$UnionPetscLib, mat::PetscMat, is::IS, diag::$PetscScalar, x::PetscVec, b::PetscVec )

    @chk ccall(
               (:MatZeroRowsColumnsIS, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, $PetscScalar, CVec, CVec),
               mat, is, diag, x, b,
              )


	return nothing
end 

"""
	MatZeroRows(petsclib::PetscLibType,mat::PetscMat, numRows::PetscInt, rows::Vector{PetscInt}, diag::PetscScalar, x::PetscVec, b::PetscVec) 
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
function MatZeroRows(petsclib::PetscLibType, mat::PetscMat, numRows::PetscInt, rows::Vector{PetscInt}, diag::PetscScalar, x::PetscVec, b::PetscVec) end

@for_petsc function MatZeroRows(petsclib::$UnionPetscLib, mat::PetscMat, numRows::$PetscInt, rows::Vector{$PetscInt}, diag::$PetscScalar, x::PetscVec, b::PetscVec )

    @chk ccall(
               (:MatZeroRows, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscScalar, CVec, CVec),
               mat, numRows, rows, diag, x, b,
              )


	return nothing
end 

"""
	MatZeroRowsIS(petsclib::PetscLibType,mat::PetscMat, is::IS, diag::PetscScalar, x::PetscVec, b::PetscVec) 
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
function MatZeroRowsIS(petsclib::PetscLibType, mat::PetscMat, is::Union{Ptr,IS}, diag::PetscScalar, x::PetscVec, b::PetscVec) end

@for_petsc function MatZeroRowsIS(petsclib::$UnionPetscLib, mat::PetscMat, is::Union{Ptr,IS}, diag::$PetscScalar, x::PetscVec, b::PetscVec )

    @chk ccall(
               (:MatZeroRowsIS, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, $PetscScalar, CVec, CVec),
               mat, is, diag, x, b,
              )


	return nothing
end 

"""
	MatZeroRowsStencil(petsclib::PetscLibType,mat::PetscMat, numRows::PetscInt, rows::Vector{MatStencil}, diag::PetscScalar, x::PetscVec, b::PetscVec) 
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
function MatZeroRowsStencil(petsclib::PetscLibType, mat::PetscMat, numRows::PetscInt, rows::Vector{MatStencil}, diag::PetscScalar, x::PetscVec, b::PetscVec) end

@for_petsc function MatZeroRowsStencil(petsclib::$UnionPetscLib, mat::PetscMat, numRows::$PetscInt, rows::Vector{MatStencil}, diag::$PetscScalar, x::PetscVec, b::PetscVec )

    @chk ccall(
               (:MatZeroRowsStencil, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{MatStencil}, $PetscScalar, CVec, CVec),
               mat, numRows, rows, diag, x, b,
              )


	return nothing
end 

"""
	MatZeroRowsColumnsStencil(petsclib::PetscLibType,mat::PetscMat, numRows::PetscInt, rows::Vector{MatStencil}, diag::PetscScalar, x::PetscVec, b::PetscVec) 
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
function MatZeroRowsColumnsStencil(petsclib::PetscLibType, mat::PetscMat, numRows::PetscInt, rows::Vector{MatStencil}, diag::PetscScalar, x::PetscVec, b::PetscVec) end

@for_petsc function MatZeroRowsColumnsStencil(petsclib::$UnionPetscLib, mat::PetscMat, numRows::$PetscInt, rows::Vector{MatStencil}, diag::$PetscScalar, x::PetscVec, b::PetscVec )

    @chk ccall(
               (:MatZeroRowsColumnsStencil, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{MatStencil}, $PetscScalar, CVec, CVec),
               mat, numRows, rows, diag, x, b,
              )


	return nothing
end 

"""
	MatZeroRowsLocal(petsclib::PetscLibType,mat::PetscMat, numRows::PetscInt, rows::Vector{PetscInt}, diag::PetscScalar, x::PetscVec, b::PetscVec) 
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
function MatZeroRowsLocal(petsclib::PetscLibType, mat::PetscMat, numRows::PetscInt, rows::Vector{PetscInt}, diag::PetscScalar, x::PetscVec, b::PetscVec) end

@for_petsc function MatZeroRowsLocal(petsclib::$UnionPetscLib, mat::PetscMat, numRows::$PetscInt, rows::Vector{$PetscInt}, diag::$PetscScalar, x::PetscVec, b::PetscVec )

    @chk ccall(
               (:MatZeroRowsLocal, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscScalar, CVec, CVec),
               mat, numRows, rows, diag, x, b,
              )


	return nothing
end 

"""
	MatZeroRowsLocalIS(petsclib::PetscLibType,mat::PetscMat, is::IS, diag::PetscScalar, x::PetscVec, b::PetscVec) 
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
function MatZeroRowsLocalIS(petsclib::PetscLibType, mat::PetscMat, is::IS, diag::PetscScalar, x::PetscVec, b::PetscVec) end

@for_petsc function MatZeroRowsLocalIS(petsclib::$UnionPetscLib, mat::PetscMat, is::IS, diag::$PetscScalar, x::PetscVec, b::PetscVec )

    @chk ccall(
               (:MatZeroRowsLocalIS, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, $PetscScalar, CVec, CVec),
               mat, is, diag, x, b,
              )


	return nothing
end 

"""
	MatZeroRowsColumnsLocal(petsclib::PetscLibType,mat::PetscMat, numRows::PetscInt, rows::Vector{PetscInt}, diag::PetscScalar, x::PetscVec, b::PetscVec) 
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
function MatZeroRowsColumnsLocal(petsclib::PetscLibType, mat::PetscMat, numRows::PetscInt, rows::Vector{PetscInt}, diag::PetscScalar, x::PetscVec, b::PetscVec) end

@for_petsc function MatZeroRowsColumnsLocal(petsclib::$UnionPetscLib, mat::PetscMat, numRows::$PetscInt, rows::Vector{$PetscInt}, diag::$PetscScalar, x::PetscVec, b::PetscVec )

    @chk ccall(
               (:MatZeroRowsColumnsLocal, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscScalar, CVec, CVec),
               mat, numRows, rows, diag, x, b,
              )


	return nothing
end 

"""
	MatZeroRowsColumnsLocalIS(petsclib::PetscLibType,mat::PetscMat, is::IS, diag::PetscScalar, x::PetscVec, b::PetscVec) 
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
function MatZeroRowsColumnsLocalIS(petsclib::PetscLibType, mat::PetscMat, is::IS, diag::PetscScalar, x::PetscVec, b::PetscVec) end

@for_petsc function MatZeroRowsColumnsLocalIS(petsclib::$UnionPetscLib, mat::PetscMat, is::IS, diag::$PetscScalar, x::PetscVec, b::PetscVec )

    @chk ccall(
               (:MatZeroRowsColumnsLocalIS, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, $PetscScalar, CVec, CVec),
               mat, is, diag, x, b,
              )


	return nothing
end 

"""
	m::PetscInt,n::PetscInt = MatGetSize(petsclib::PetscLibType,mat::PetscMat) 
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
	m::PetscInt,n::PetscInt = MatGetLocalSize(petsclib::PetscLibType,mat::PetscMat) 
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
function MatGetLocalSize(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatGetLocalSize(petsclib::$UnionPetscLib, mat::PetscMat )
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
	m::PetscInt,n::PetscInt = MatGetOwnershipRangeColumn(petsclib::PetscLibType,mat::PetscMat) 
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
function MatGetOwnershipRangeColumn(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatGetOwnershipRangeColumn(petsclib::$UnionPetscLib, mat::PetscMat )
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
	m::PetscInt,n::PetscInt = MatGetOwnershipRange(petsclib::PetscLibType,mat::PetscMat) 
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
function MatGetOwnershipRange(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatGetOwnershipRange(petsclib::$UnionPetscLib, mat::PetscMat )
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
	ranges::Vector{PetscInt} = MatGetOwnershipRanges(petsclib::PetscLibType,mat::PetscMat) 
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
function MatGetOwnershipRanges(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatGetOwnershipRanges(petsclib::$UnionPetscLib, mat::PetscMat )
	ranges_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:MatGetOwnershipRanges, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscInt}}),
               mat, ranges_,
              )

	ranges = unsafe_wrap(Array, ranges_[], VecGetLocalSize(petsclib, x); own = false)

	return ranges
end 

"""
	ranges::Vector{PetscInt} = MatGetOwnershipRangesColumn(petsclib::PetscLibType,mat::PetscMat) 
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
function MatGetOwnershipRangesColumn(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatGetOwnershipRangesColumn(petsclib::$UnionPetscLib, mat::PetscMat )
	ranges_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:MatGetOwnershipRangesColumn, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscInt}}),
               mat, ranges_,
              )

	ranges = unsafe_wrap(Array, ranges_[], VecGetLocalSize(petsclib, x); own = false)

	return ranges
end 

"""
	MatGetOwnershipIS(petsclib::PetscLibType,A::PetscMat, rows::IS, cols::IS) 
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
function MatGetOwnershipIS(petsclib::PetscLibType, A::PetscMat, rows::IS, cols::IS) end

@for_petsc function MatGetOwnershipIS(petsclib::$UnionPetscLib, A::PetscMat, rows::IS, cols::IS )

    @chk ccall(
               (:MatGetOwnershipIS, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CIS}, Ptr{CIS}),
               A, rows, cols,
              )


	return nothing
end 

"""
	MatILUFactorSymbolic(petsclib::PetscLibType,fact::PetscMat, mat::PetscMat, row::IS, col::IS, info::MatFactorInfo) 
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
function MatILUFactorSymbolic(petsclib::PetscLibType, fact::PetscMat, mat::PetscMat, row::IS, col::IS, info::MatFactorInfo) end

@for_petsc function MatILUFactorSymbolic(petsclib::$UnionPetscLib, fact::PetscMat, mat::PetscMat, row::IS, col::IS, info::MatFactorInfo )

    @chk ccall(
               (:MatILUFactorSymbolic, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CIS, CIS, Ptr{MatFactorInfo}),
               fact, mat, row, col, info,
              )


	return nothing
end 

"""
	MatICCFactorSymbolic(petsclib::PetscLibType,fact::PetscMat, mat::PetscMat, perm::IS, info::MatFactorInfo) 
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
function MatICCFactorSymbolic(petsclib::PetscLibType, fact::PetscMat, mat::PetscMat, perm::IS, info::MatFactorInfo) end

@for_petsc function MatICCFactorSymbolic(petsclib::$UnionPetscLib, fact::PetscMat, mat::PetscMat, perm::IS, info::MatFactorInfo )

    @chk ccall(
               (:MatICCFactorSymbolic, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CIS, Ptr{MatFactorInfo}),
               fact, mat, perm, info,
              )


	return nothing
end 

"""
	submat::Vector{PetscMat} = MatCreateSubMatrices(petsclib::PetscLibType,mat::PetscMat, n::PetscInt, irow::Vector{IS}, icol::Vector{IS}, scall::MatReuse) 
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
function MatCreateSubMatrices(petsclib::PetscLibType, mat::PetscMat, n::PetscInt, irow::Vector{IS}, icol::Vector{IS}, scall::MatReuse) end

@for_petsc function MatCreateSubMatrices(petsclib::$UnionPetscLib, mat::PetscMat, n::$PetscInt, irow::Vector{IS}, icol::Vector{IS}, scall::MatReuse )
	submat_ = Ref{Ptr{PetscMat}}()

    @chk ccall(
               (:MatCreateSubMatrices, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CIS}, Ptr{CIS}, MatReuse, Ptr{Ptr{CMat}}),
               mat, n, irow, icol, scall, submat_,
              )

	submat = unsafe_wrap(Array, submat_[], VecGetLocalSize(petsclib, x); own = false)

	return submat
end 

"""
	submat::Vector{PetscMat} = MatCreateSubMatricesMPI(petsclib::PetscLibType,mat::PetscMat, n::PetscInt, irow::Vector{IS}, icol::Vector{IS}, scall::MatReuse) 
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
function MatCreateSubMatricesMPI(petsclib::PetscLibType, mat::PetscMat, n::PetscInt, irow::Vector{IS}, icol::Vector{IS}, scall::MatReuse) end

@for_petsc function MatCreateSubMatricesMPI(petsclib::$UnionPetscLib, mat::PetscMat, n::$PetscInt, irow::Vector{IS}, icol::Vector{IS}, scall::MatReuse )
	submat_ = Ref{Ptr{PetscMat}}()

    @chk ccall(
               (:MatCreateSubMatricesMPI, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CIS}, Ptr{CIS}, MatReuse, Ptr{Ptr{CMat}}),
               mat, n, irow, icol, scall, submat_,
              )

	submat = unsafe_wrap(Array, submat_[], VecGetLocalSize(petsclib, x); own = false)

	return submat
end 

"""
	MatDestroyMatrices(petsclib::PetscLibType,n::PetscInt, mat::Vector{PetscMat}) 
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
function MatDestroyMatrices(petsclib::PetscLibType, n::PetscInt, mat::Vector{PetscMat}) end

@for_petsc function MatDestroyMatrices(petsclib::$UnionPetscLib, n::$PetscInt, mat::Vector{PetscMat} )
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
	MatDestroySubMatrices(petsclib::PetscLibType,n::PetscInt, mat::Vector{PetscMat}) 
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
function MatDestroySubMatrices(petsclib::PetscLibType, n::PetscInt, mat::Vector{PetscMat}) end

@for_petsc function MatDestroySubMatrices(petsclib::$UnionPetscLib, n::$PetscInt, mat::Vector{PetscMat} )
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
	MatGetSeqNonzeroStructure(petsclib::PetscLibType,mat::PetscMat, matstruct::PetscMat) 
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
function MatGetSeqNonzeroStructure(petsclib::PetscLibType, mat::PetscMat, matstruct::PetscMat) end

@for_petsc function MatGetSeqNonzeroStructure(petsclib::$UnionPetscLib, mat::PetscMat, matstruct::PetscMat )
	matstruct_ = Ref(matstruct.ptr)

    @chk ccall(
               (:MatGetSeqNonzeroStructure, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               mat, matstruct_,
              )

	matstruct.ptr = C_NULL

	return nothing
end 

"""
	MatDestroySeqNonzeroStructure(petsclib::PetscLibType,mat::PetscMat) 
Destroys matrix obtained with `MatGetSeqNonzeroStructure()`.

Collective

Input Parameter:
- `mat` - the matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatGetSeqNonzeroStructure()`

# External Links
$(_doc_external("Mat/MatDestroySeqNonzeroStructure"))
"""
function MatDestroySeqNonzeroStructure(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatDestroySeqNonzeroStructure(petsclib::$UnionPetscLib, mat::PetscMat )
	mat_ = Ref(mat.ptr)

    @chk ccall(
               (:MatDestroySeqNonzeroStructure, $petsc_library),
               PetscErrorCode,
               (Ptr{CMat},),
               mat_,
              )

	mat.ptr = C_NULL

	return nothing
end 

"""
	MatIncreaseOverlap(petsclib::PetscLibType,mat::PetscMat, n::PetscInt, is::Vector{IS}, ov::PetscInt) 
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
function MatIncreaseOverlap(petsclib::PetscLibType, mat::PetscMat, n::PetscInt, is::Vector{IS}, ov::PetscInt) end

@for_petsc function MatIncreaseOverlap(petsclib::$UnionPetscLib, mat::PetscMat, n::$PetscInt, is::Vector{IS}, ov::$PetscInt )

    @chk ccall(
               (:MatIncreaseOverlap, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CIS}, $PetscInt),
               mat, n, is, ov,
              )


	return nothing
end 

"""
	MatIncreaseOverlapSplit(petsclib::PetscLibType,mat::PetscMat, n::PetscInt, is::Vector{IS}, ov::PetscInt) 
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
function MatIncreaseOverlapSplit(petsclib::PetscLibType, mat::PetscMat, n::PetscInt, is::Vector{IS}, ov::PetscInt) end

@for_petsc function MatIncreaseOverlapSplit(petsclib::$UnionPetscLib, mat::PetscMat, n::$PetscInt, is::Vector{IS}, ov::$PetscInt )

    @chk ccall(
               (:MatIncreaseOverlapSplit, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CIS}, $PetscInt),
               mat, n, is, ov,
              )


	return nothing
end 

"""
	bs::PetscInt = MatGetBlockSize(petsclib::PetscLibType,mat::PetscMat) 
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
function MatGetBlockSize(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatGetBlockSize(petsclib::$UnionPetscLib, mat::PetscMat )
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
	rbs::PetscInt,cbs::PetscInt = MatGetBlockSizes(petsclib::PetscLibType,mat::PetscMat) 
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
function MatGetBlockSizes(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatGetBlockSizes(petsclib::$UnionPetscLib, mat::PetscMat )
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
	MatSetBlockSize(petsclib::PetscLibType,mat::PetscMat, bs::PetscInt) 
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
function MatSetBlockSize(petsclib::PetscLibType, mat::PetscMat, bs::PetscInt) end

@for_petsc function MatSetBlockSize(petsclib::$UnionPetscLib, mat::PetscMat, bs::$PetscInt )

    @chk ccall(
               (:MatSetBlockSize, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt),
               mat, bs,
              )


	return nothing
end 

"""
	MatComputeVariableBlockEnvelope(petsclib::PetscLibType,mat::PetscMat) 
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
function MatComputeVariableBlockEnvelope(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatComputeVariableBlockEnvelope(petsclib::$UnionPetscLib, mat::PetscMat )

    @chk ccall(
               (:MatComputeVariableBlockEnvelope, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	MatInvertVariableBlockEnvelope(petsclib::PetscLibType,A::PetscMat, reuse::MatReuse, C::PetscMat) 
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
function MatInvertVariableBlockEnvelope(petsclib::PetscLibType, A::PetscMat, reuse::MatReuse, C::PetscMat) end

@for_petsc function MatInvertVariableBlockEnvelope(petsclib::$UnionPetscLib, A::PetscMat, reuse::MatReuse, C::PetscMat )
	C_ = Ref(C.ptr)

    @chk ccall(
               (:MatInvertVariableBlockEnvelope, $petsc_library),
               PetscErrorCode,
               (CMat, MatReuse, Ptr{CMat}),
               A, reuse, C_,
              )

	C.ptr = C_NULL

	return nothing
end 

"""
	MatSetVariableBlockSizes(petsclib::PetscLibType,mat::PetscMat, nblocks::PetscInt, bsizes::Vector{PetscInt}) 
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
function MatSetVariableBlockSizes(petsclib::PetscLibType, mat::PetscMat, nblocks::PetscInt, bsizes::Vector{PetscInt}) end

@for_petsc function MatSetVariableBlockSizes(petsclib::$UnionPetscLib, mat::PetscMat, nblocks::$PetscInt, bsizes::Vector{$PetscInt} )

    @chk ccall(
               (:MatSetVariableBlockSizes, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}),
               mat, nblocks, bsizes,
              )


	return nothing
end 

"""
	nblocks::PetscInt,bsizes::Vector{PetscInt} = MatGetVariableBlockSizes(petsclib::PetscLibType,mat::PetscMat) 
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
function MatGetVariableBlockSizes(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatGetVariableBlockSizes(petsclib::$UnionPetscLib, mat::PetscMat )
	nblocks_ = Ref{$PetscInt}()
	bsizes_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:MatGetVariableBlockSizes, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
               mat, nblocks_, bsizes_,
              )

	nblocks = nblocks_[]
	bsizes = unsafe_wrap(Array, bsizes_[], VecGetLocalSize(petsclib, x); own = false)

	return nblocks,bsizes
end 

"""
	MatSelectVariableBlockSizes(petsclib::PetscLibType,subA::PetscMat, A::PetscMat, isrow::IS) 
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
function MatSelectVariableBlockSizes(petsclib::PetscLibType, subA::PetscMat, A::PetscMat, isrow::IS) end

@for_petsc function MatSelectVariableBlockSizes(petsclib::$UnionPetscLib, subA::PetscMat, A::PetscMat, isrow::IS )

    @chk ccall(
               (:MatSelectVariableBlockSizes, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CIS),
               subA, A, isrow,
              )


	return nothing
end 

"""
	MatSetBlockSizes(petsclib::PetscLibType,mat::PetscMat, rbs::PetscInt, cbs::PetscInt) 
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
function MatSetBlockSizes(petsclib::PetscLibType, mat::PetscMat, rbs::PetscInt, cbs::PetscInt) end

@for_petsc function MatSetBlockSizes(petsclib::$UnionPetscLib, mat::PetscMat, rbs::$PetscInt, cbs::$PetscInt )

    @chk ccall(
               (:MatSetBlockSizes, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt),
               mat, rbs, cbs,
              )


	return nothing
end 

"""
	MatSetBlockSizesFromMats(petsclib::PetscLibType,mat::PetscMat, fromRow::PetscMat, fromCol::PetscMat) 
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
function MatSetBlockSizesFromMats(petsclib::PetscLibType, mat::PetscMat, fromRow::PetscMat, fromCol::PetscMat) end

@for_petsc function MatSetBlockSizesFromMats(petsclib::$UnionPetscLib, mat::PetscMat, fromRow::PetscMat, fromCol::PetscMat )

    @chk ccall(
               (:MatSetBlockSizesFromMats, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat),
               mat, fromRow, fromCol,
              )


	return nothing
end 

"""
	MatResidual(petsclib::PetscLibType,mat::PetscMat, b::PetscVec, x::PetscVec, r::PetscVec) 
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
function MatResidual(petsclib::PetscLibType, mat::PetscMat, b::PetscVec, x::PetscVec, r::PetscVec) end

@for_petsc function MatResidual(petsclib::$UnionPetscLib, mat::PetscMat, b::PetscVec, x::PetscVec, r::PetscVec )

    @chk ccall(
               (:MatResidual, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec, CVec),
               mat, b, x, r,
              )


	return nothing
end 

"""
	n::PetscInt,ia::Vector{PetscInt},ja::Vector{PetscInt},done::PetscBool = MatGetRowIJ(petsclib::PetscLibType,mat::PetscMat, shift::PetscInt, symmetric::PetscBool, inodecompressed::PetscBool) 
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
function MatGetRowIJ(petsclib::PetscLibType, mat::PetscMat, shift::PetscInt, symmetric::PetscBool, inodecompressed::PetscBool) end

@for_petsc function MatGetRowIJ(petsclib::$UnionPetscLib, mat::PetscMat, shift::$PetscInt, symmetric::PetscBool, inodecompressed::PetscBool )
	n_ = Ref{$PetscInt}()
	ia_ = Ref{Ptr{$PetscInt}}()
	ja_ = Ref{Ptr{$PetscInt}}()
	done_ = Ref{PetscBool}()

    @chk ccall(
               (:MatGetRowIJ, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, PetscBool, PetscBool, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}, Ptr{PetscBool}),
               mat, shift, symmetric, inodecompressed, n_, ia_, ja_, done_,
              )

	n = n_[]
	ia = unsafe_wrap(Array, ia_[], VecGetLocalSize(petsclib, x); own = false)
	ja = unsafe_wrap(Array, ja_[], VecGetLocalSize(petsclib, x); own = false)
	done = done_[]

	return n,ia,ja,done
end 

"""
	n::PetscInt,ia::Vector{PetscInt},ja::Vector{PetscInt},done::PetscBool = MatGetColumnIJ(petsclib::PetscLibType,mat::PetscMat, shift::PetscInt, symmetric::PetscBool, inodecompressed::PetscBool) 
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
function MatGetColumnIJ(petsclib::PetscLibType, mat::PetscMat, shift::PetscInt, symmetric::PetscBool, inodecompressed::PetscBool) end

@for_petsc function MatGetColumnIJ(petsclib::$UnionPetscLib, mat::PetscMat, shift::$PetscInt, symmetric::PetscBool, inodecompressed::PetscBool )
	n_ = Ref{$PetscInt}()
	ia_ = Ref{Ptr{$PetscInt}}()
	ja_ = Ref{Ptr{$PetscInt}}()
	done_ = Ref{PetscBool}()

    @chk ccall(
               (:MatGetColumnIJ, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, PetscBool, PetscBool, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}, Ptr{PetscBool}),
               mat, shift, symmetric, inodecompressed, n_, ia_, ja_, done_,
              )

	n = n_[]
	ia = unsafe_wrap(Array, ia_[], VecGetLocalSize(petsclib, x); own = false)
	ja = unsafe_wrap(Array, ja_[], VecGetLocalSize(petsclib, x); own = false)
	done = done_[]

	return n,ia,ja,done
end 

"""
	n::PetscInt,ia::Vector{PetscInt},ja::Vector{PetscInt},done::PetscBool = MatRestoreRowIJ(petsclib::PetscLibType,mat::PetscMat, shift::PetscInt, symmetric::PetscBool, inodecompressed::PetscBool) 
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
function MatRestoreRowIJ(petsclib::PetscLibType, mat::PetscMat, shift::PetscInt, symmetric::PetscBool, inodecompressed::PetscBool) end

@for_petsc function MatRestoreRowIJ(petsclib::$UnionPetscLib, mat::PetscMat, shift::$PetscInt, symmetric::PetscBool, inodecompressed::PetscBool )
	n_ = Ref{$PetscInt}()
	ia_ = Ref{Ptr{$PetscInt}}()
	ja_ = Ref{Ptr{$PetscInt}}()
	done_ = Ref{PetscBool}()

    @chk ccall(
               (:MatRestoreRowIJ, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, PetscBool, PetscBool, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}, Ptr{PetscBool}),
               mat, shift, symmetric, inodecompressed, n_, ia_, ja_, done_,
              )

	n = n_[]
	ia = unsafe_wrap(Array, ia_[], VecGetLocalSize(petsclib, x); own = false)
	ja = unsafe_wrap(Array, ja_[], VecGetLocalSize(petsclib, x); own = false)
	done = done_[]

	return n,ia,ja,done
end 

"""
	n::PetscInt,ia::Vector{PetscInt},ja::Vector{PetscInt},done::PetscBool = MatRestoreColumnIJ(petsclib::PetscLibType,mat::PetscMat, shift::PetscInt, symmetric::PetscBool, inodecompressed::PetscBool) 
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
function MatRestoreColumnIJ(petsclib::PetscLibType, mat::PetscMat, shift::PetscInt, symmetric::PetscBool, inodecompressed::PetscBool) end

@for_petsc function MatRestoreColumnIJ(petsclib::$UnionPetscLib, mat::PetscMat, shift::$PetscInt, symmetric::PetscBool, inodecompressed::PetscBool )
	n_ = Ref{$PetscInt}()
	ia_ = Ref{Ptr{$PetscInt}}()
	ja_ = Ref{Ptr{$PetscInt}}()
	done_ = Ref{PetscBool}()

    @chk ccall(
               (:MatRestoreColumnIJ, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, PetscBool, PetscBool, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}, Ptr{PetscBool}),
               mat, shift, symmetric, inodecompressed, n_, ia_, ja_, done_,
              )

	n = n_[]
	ia = unsafe_wrap(Array, ia_[], VecGetLocalSize(petsclib, x); own = false)
	ja = unsafe_wrap(Array, ja_[], VecGetLocalSize(petsclib, x); own = false)
	done = done_[]

	return n,ia,ja,done
end 

"""
	MatSetUnfactored(petsclib::PetscLibType,mat::PetscMat) 
Resets a factored matrix to be treated as unfactored.

Logically Collective

Input Parameter:
- `mat` - the factored matrix to be reset

Level: developer

-seealso: [](ch_matrices), `Mat`, `PCFactorSetUseInPlace()`, `PCFactorGetUseInPlace()`

# External Links
$(_doc_external("Mat/MatSetUnfactored"))
"""
function MatSetUnfactored(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatSetUnfactored(petsclib::$UnionPetscLib, mat::PetscMat )

    @chk ccall(
               (:MatSetUnfactored, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	newmat::PetscMat = MatCreateSubMatrix(petsclib::PetscLibType,mat::PetscMat, isrow::IS, iscol::IS, cll::MatReuse) 
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
function MatCreateSubMatrix(petsclib::PetscLibType, mat::PetscMat, isrow::IS, iscol::IS, cll::MatReuse) end

@for_petsc function MatCreateSubMatrix(petsclib::$UnionPetscLib, mat::PetscMat, isrow::IS, iscol::IS, cll::MatReuse )
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
	MatPropagateSymmetryOptions(petsclib::PetscLibType,A::PetscMat, B::PetscMat) 
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
function MatPropagateSymmetryOptions(petsclib::PetscLibType, A::PetscMat, B::PetscMat) end

@for_petsc function MatPropagateSymmetryOptions(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat )

    @chk ccall(
               (:MatPropagateSymmetryOptions, $petsc_library),
               PetscErrorCode,
               (CMat, CMat),
               A, B,
              )


	return nothing
end 

"""
	MatStashSetInitialSize(petsclib::PetscLibType,mat::PetscMat, size::PetscInt, bsize::PetscInt) 
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
function MatStashSetInitialSize(petsclib::PetscLibType, mat::PetscMat, size::PetscInt, bsize::PetscInt) end

@for_petsc function MatStashSetInitialSize(petsclib::$UnionPetscLib, mat::PetscMat, size::$PetscInt, bsize::$PetscInt )

    @chk ccall(
               (:MatStashSetInitialSize, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt),
               mat, size, bsize,
              )


	return nothing
end 

"""
	MatInterpolateAdd(petsclib::PetscLibType,A::PetscMat, x::PetscVec, y::PetscVec, w::PetscVec) 
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
function MatInterpolateAdd(petsclib::PetscLibType, A::PetscMat, x::PetscVec, y::PetscVec, w::PetscVec) end

@for_petsc function MatInterpolateAdd(petsclib::$UnionPetscLib, A::PetscMat, x::PetscVec, y::PetscVec, w::PetscVec )

    @chk ccall(
               (:MatInterpolateAdd, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec, CVec),
               A, x, y, w,
              )


	return nothing
end 

"""
	MatInterpolate(petsclib::PetscLibType,A::PetscMat, x::PetscVec, y::PetscVec) 
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
function MatInterpolate(petsclib::PetscLibType, A::PetscMat, x::PetscVec, y::PetscVec) end

@for_petsc function MatInterpolate(petsclib::$UnionPetscLib, A::PetscMat, x::PetscVec, y::PetscVec )

    @chk ccall(
               (:MatInterpolate, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               A, x, y,
              )


	return nothing
end 

"""
	MatRestrict(petsclib::PetscLibType,A::PetscMat, x::PetscVec, y::PetscVec) 
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
function MatRestrict(petsclib::PetscLibType, A::PetscMat, x::PetscVec, y::PetscVec) end

@for_petsc function MatRestrict(petsclib::$UnionPetscLib, A::PetscMat, x::PetscVec, y::PetscVec )

    @chk ccall(
               (:MatRestrict, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               A, x, y,
              )


	return nothing
end 

"""
	MatMatInterpolateAdd(petsclib::PetscLibType,A::PetscMat, x::PetscMat, w::PetscMat, y::PetscMat) 
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
function MatMatInterpolateAdd(petsclib::PetscLibType, A::PetscMat, x::PetscMat, w::PetscMat, y::PetscMat) end

@for_petsc function MatMatInterpolateAdd(petsclib::$UnionPetscLib, A::PetscMat, x::PetscMat, w::PetscMat, y::PetscMat )
	y_ = Ref(y.ptr)

    @chk ccall(
               (:MatMatInterpolateAdd, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat, Ptr{CMat}),
               A, x, w, y_,
              )

	y.ptr = C_NULL

	return nothing
end 

"""
	MatMatInterpolate(petsclib::PetscLibType,A::PetscMat, x::PetscMat, y::PetscMat) 
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
function MatMatInterpolate(petsclib::PetscLibType, A::PetscMat, x::PetscMat, y::PetscMat) end

@for_petsc function MatMatInterpolate(petsclib::$UnionPetscLib, A::PetscMat, x::PetscMat, y::PetscMat )
	y_ = Ref(y.ptr)

    @chk ccall(
               (:MatMatInterpolate, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, Ptr{CMat}),
               A, x, y_,
              )

	y.ptr = C_NULL

	return nothing
end 

"""
	MatMatRestrict(petsclib::PetscLibType,A::PetscMat, x::PetscMat, y::PetscMat) 
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
function MatMatRestrict(petsclib::PetscLibType, A::PetscMat, x::PetscMat, y::PetscMat) end

@for_petsc function MatMatRestrict(petsclib::$UnionPetscLib, A::PetscMat, x::PetscMat, y::PetscMat )
	y_ = Ref(y.ptr)

    @chk ccall(
               (:MatMatRestrict, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, Ptr{CMat}),
               A, x, y_,
              )

	y.ptr = C_NULL

	return nothing
end 

"""
	MatGetNullSpace(petsclib::PetscLibType,mat::PetscMat, nullsp::MatNullSpace) 
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
function MatGetNullSpace(petsclib::PetscLibType, mat::PetscMat, nullsp::MatNullSpace) end

@for_petsc function MatGetNullSpace(petsclib::$UnionPetscLib, mat::PetscMat, nullsp::MatNullSpace )

    @chk ccall(
               (:MatGetNullSpace, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatNullSpace}),
               mat, nullsp,
              )


	return nothing
end 

"""
	MatGetNullSpaces(petsclib::PetscLibType,n::PetscInt, mat::Vector{PetscMat}, nullsp::Vector{MatNullSpace}) 
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
function MatGetNullSpaces(petsclib::PetscLibType, n::PetscInt, mat::Vector{PetscMat}, nullsp::Vector{MatNullSpace}) end

@for_petsc function MatGetNullSpaces(petsclib::$UnionPetscLib, n::$PetscInt, mat::Vector{PetscMat}, nullsp::Vector{MatNullSpace} )
	nullsp_ = Ref(pointer(nullsp))

    @chk ccall(
               (:MatGetNullSpaces, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{CMat}, Ptr{Ptr{MatNullSpace}}),
               n, mat, nullsp_,
              )


	return nothing
end 

"""
	MatRestoreNullSpaces(petsclib::PetscLibType,n::PetscInt, mat::Vector{PetscMat}, nullsp::Vector{MatNullSpace}) 
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
function MatRestoreNullSpaces(petsclib::PetscLibType, n::PetscInt, mat::Vector{PetscMat}, nullsp::Vector{MatNullSpace}) end

@for_petsc function MatRestoreNullSpaces(petsclib::$UnionPetscLib, n::$PetscInt, mat::Vector{PetscMat}, nullsp::Vector{MatNullSpace} )
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
	MatSetNullSpace(petsclib::PetscLibType,mat::PetscMat, nullsp::MatNullSpace) 
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
function MatSetNullSpace(petsclib::PetscLibType, mat::PetscMat, nullsp::MatNullSpace) end

@for_petsc function MatSetNullSpace(petsclib::$UnionPetscLib, mat::PetscMat, nullsp::MatNullSpace )

    @chk ccall(
               (:MatSetNullSpace, $petsc_library),
               PetscErrorCode,
               (CMat, MatNullSpace),
               mat, nullsp,
              )


	return nothing
end 

"""
	MatGetTransposeNullSpace(petsclib::PetscLibType,mat::PetscMat, nullsp::MatNullSpace) 
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
function MatGetTransposeNullSpace(petsclib::PetscLibType, mat::PetscMat, nullsp::MatNullSpace) end

@for_petsc function MatGetTransposeNullSpace(petsclib::$UnionPetscLib, mat::PetscMat, nullsp::MatNullSpace )

    @chk ccall(
               (:MatGetTransposeNullSpace, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatNullSpace}),
               mat, nullsp,
              )


	return nothing
end 

"""
	MatSetTransposeNullSpace(petsclib::PetscLibType,mat::PetscMat, nullsp::MatNullSpace) 
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
function MatSetTransposeNullSpace(petsclib::PetscLibType, mat::PetscMat, nullsp::MatNullSpace) end

@for_petsc function MatSetTransposeNullSpace(petsclib::$UnionPetscLib, mat::PetscMat, nullsp::MatNullSpace )

    @chk ccall(
               (:MatSetTransposeNullSpace, $petsc_library),
               PetscErrorCode,
               (CMat, MatNullSpace),
               mat, nullsp,
              )


	return nothing
end 

"""
	MatSetNearNullSpace(petsclib::PetscLibType,mat::PetscMat, nullsp::MatNullSpace) 
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
function MatSetNearNullSpace(petsclib::PetscLibType, mat::PetscMat, nullsp::MatNullSpace) end

@for_petsc function MatSetNearNullSpace(petsclib::$UnionPetscLib, mat::PetscMat, nullsp::MatNullSpace )

    @chk ccall(
               (:MatSetNearNullSpace, $petsc_library),
               PetscErrorCode,
               (CMat, MatNullSpace),
               mat, nullsp,
              )


	return nothing
end 

"""
	MatGetNearNullSpace(petsclib::PetscLibType,mat::PetscMat, nullsp::MatNullSpace) 
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
function MatGetNearNullSpace(petsclib::PetscLibType, mat::PetscMat, nullsp::MatNullSpace) end

@for_petsc function MatGetNearNullSpace(petsclib::$UnionPetscLib, mat::PetscMat, nullsp::MatNullSpace )

    @chk ccall(
               (:MatGetNearNullSpace, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatNullSpace}),
               mat, nullsp,
              )


	return nothing
end 

"""
	MatICCFactor(petsclib::PetscLibType,mat::PetscMat, row::IS, info::MatFactorInfo) 
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
function MatICCFactor(petsclib::PetscLibType, mat::PetscMat, row::IS, info::MatFactorInfo) end

@for_petsc function MatICCFactor(petsclib::$UnionPetscLib, mat::PetscMat, row::IS, info::MatFactorInfo )

    @chk ccall(
               (:MatICCFactor, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, Ptr{MatFactorInfo}),
               mat, row, info,
              )


	return nothing
end 

"""
	MatDiagonalScaleLocal(petsclib::PetscLibType,mat::PetscMat, diag::PetscVec) 
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
function MatDiagonalScaleLocal(petsclib::PetscLibType, mat::PetscMat, diag::PetscVec) end

@for_petsc function MatDiagonalScaleLocal(petsclib::$UnionPetscLib, mat::PetscMat, diag::PetscVec )

    @chk ccall(
               (:MatDiagonalScaleLocal, $petsc_library),
               PetscErrorCode,
               (CMat, CVec),
               mat, diag,
              )


	return nothing
end 

"""
	nneg::PetscInt,nzero::PetscInt,npos::PetscInt = MatGetInertia(petsclib::PetscLibType,mat::PetscMat) 
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
function MatGetInertia(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatGetInertia(petsclib::$UnionPetscLib, mat::PetscMat )
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
	MatSolves(petsclib::PetscLibType,mat::PetscMat, b::Vecs, x::Vecs) 
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
function MatSolves(petsclib::PetscLibType, mat::PetscMat, b::Vecs, x::Vecs) end

@for_petsc function MatSolves(petsclib::$UnionPetscLib, mat::PetscMat, b::Vecs, x::Vecs )

    @chk ccall(
               (:MatSolves, $petsc_library),
               PetscErrorCode,
               (CMat, Vecs, Vecs),
               mat, b, x,
              )


	return nothing
end 

"""
	flg::PetscBool = MatIsSymmetric(petsclib::PetscLibType,A::PetscMat, tol::PetscReal) 
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
function MatIsSymmetric(petsclib::PetscLibType, A::PetscMat, tol::PetscReal) end

@for_petsc function MatIsSymmetric(petsclib::$UnionPetscLib, A::PetscMat, tol::$PetscReal )
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
	flg::PetscBool = MatIsHermitian(petsclib::PetscLibType,A::PetscMat, tol::PetscReal) 
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
function MatIsHermitian(petsclib::PetscLibType, A::PetscMat, tol::PetscReal) end

@for_petsc function MatIsHermitian(petsclib::$UnionPetscLib, A::PetscMat, tol::$PetscScalar )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:MatIsHermitian, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscScalar, Ptr{PetscBool}),
               A, tol, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	set::PetscBool,flg::PetscBool = MatIsSymmetricKnown(petsclib::PetscLibType,A::PetscMat) 
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
function MatIsSymmetricKnown(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatIsSymmetricKnown(petsclib::$UnionPetscLib, A::PetscMat )
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
	set::PetscBool,flg::PetscBool = MatIsSPDKnown(petsclib::PetscLibType,A::PetscMat) 
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
function MatIsSPDKnown(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatIsSPDKnown(petsclib::$UnionPetscLib, A::PetscMat )
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
	set::PetscBool,flg::PetscBool = MatIsHermitianKnown(petsclib::PetscLibType,A::PetscMat) 
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
function MatIsHermitianKnown(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatIsHermitianKnown(petsclib::$UnionPetscLib, A::PetscMat )
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
	flg::PetscBool = MatIsStructurallySymmetric(petsclib::PetscLibType,A::PetscMat) 
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
function MatIsStructurallySymmetric(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatIsStructurallySymmetric(petsclib::$UnionPetscLib, A::PetscMat )
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
	set::PetscBool,flg::PetscBool = MatIsStructurallySymmetricKnown(petsclib::PetscLibType,A::PetscMat) 
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
function MatIsStructurallySymmetricKnown(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatIsStructurallySymmetricKnown(petsclib::$UnionPetscLib, A::PetscMat )
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
	nstash::PetscInt,reallocs::PetscInt,bnstash::PetscInt,breallocs::PetscInt = MatStashGetInfo(petsclib::PetscLibType,mat::PetscMat) 
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
function MatStashGetInfo(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatStashGetInfo(petsclib::$UnionPetscLib, mat::PetscMat )
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
	right::PetscVec,left::PetscVec = MatCreateVecs(petsclib::PetscLibType,mat::PetscMat) 
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
function MatCreateVecs(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatCreateVecs(petsclib::$UnionPetscLib, mat::PetscMat )
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
	MatFactorInfoInitialize(petsclib::PetscLibType,info::MatFactorInfo) 
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
function MatFactorInfoInitialize(petsclib::PetscLibType, info::MatFactorInfo) end

@for_petsc function MatFactorInfoInitialize(petsclib::$UnionPetscLib, info::MatFactorInfo )

    @chk ccall(
               (:MatFactorInfoInitialize, $petsc_library),
               PetscErrorCode,
               (Ptr{MatFactorInfo},),
               info,
              )


	return nothing
end 

"""
	MatFactorSetSchurIS(petsclib::PetscLibType,mat::PetscMat, is::IS) 
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
function MatFactorSetSchurIS(petsclib::PetscLibType, mat::PetscMat, is::IS) end

@for_petsc function MatFactorSetSchurIS(petsclib::$UnionPetscLib, mat::PetscMat, is::IS )

    @chk ccall(
               (:MatFactorSetSchurIS, $petsc_library),
               PetscErrorCode,
               (CMat, CIS),
               mat, is,
              )


	return nothing
end 

"""
	S::PetscMat,status::MatFactorSchurStatus = MatFactorCreateSchurComplement(petsclib::PetscLibType,F::PetscMat) 
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
function MatFactorCreateSchurComplement(petsclib::PetscLibType, F::PetscMat) end

@for_petsc function MatFactorCreateSchurComplement(petsclib::$UnionPetscLib, F::PetscMat )
	S_ = Ref{CMat}()
	status_ = Ref{MatFactorSchurStatus}()

    @chk ccall(
               (:MatFactorCreateSchurComplement, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}, Ptr{MatFactorSchurStatus}),
               F, S_, status_,
              )

	S = PetscMat(S_[], petsclib)
	status = status_[]

	return S,status
end 

"""
	MatFactorGetSchurComplement(petsclib::PetscLibType,F::PetscMat, S::PetscMat, status::MatFactorSchurStatus) 
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
function MatFactorGetSchurComplement(petsclib::PetscLibType, F::PetscMat, S::Union{Ptr,PetscMat}, status::Union{Ptr,MatFactorSchurStatus}) end

@for_petsc function MatFactorGetSchurComplement(petsclib::$UnionPetscLib, F::PetscMat, S::Union{Ptr,PetscMat}, status::Union{Ptr,MatFactorSchurStatus})
	S_ = Ref(S.ptr)

    @chk ccall(
               (:MatFactorGetSchurComplement, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}, Ptr{MatFactorSchurStatus}),
               F, S_, status,
              )

	S.ptr = C_NULL

	return nothing
end 

"""
	MatFactorRestoreSchurComplement(petsclib::PetscLibType,F::PetscMat, S::PetscMat, status::MatFactorSchurStatus) 
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
function MatFactorRestoreSchurComplement(petsclib::PetscLibType, F::PetscMat, S::PetscMat, status::MatFactorSchurStatus) end

@for_petsc function MatFactorRestoreSchurComplement(petsclib::$UnionPetscLib, F::PetscMat, S::PetscMat, status::MatFactorSchurStatus )
	S_ = Ref(S.ptr)

    @chk ccall(
               (:MatFactorRestoreSchurComplement, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}, MatFactorSchurStatus),
               F, S_, status,
              )

	S.ptr = C_NULL

	return nothing
end 

"""
	MatFactorSolveSchurComplementTranspose(petsclib::PetscLibType,F::PetscMat, rhs::PetscVec, sol::PetscVec) 
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
function MatFactorSolveSchurComplementTranspose(petsclib::PetscLibType, F::PetscMat, rhs::PetscVec, sol::PetscVec) end

@for_petsc function MatFactorSolveSchurComplementTranspose(petsclib::$UnionPetscLib, F::PetscMat, rhs::PetscVec, sol::PetscVec )

    @chk ccall(
               (:MatFactorSolveSchurComplementTranspose, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               F, rhs, sol,
              )


	return nothing
end 

"""
	MatFactorSolveSchurComplement(petsclib::PetscLibType,F::PetscMat, rhs::PetscVec, sol::PetscVec) 
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
function MatFactorSolveSchurComplement(petsclib::PetscLibType, F::PetscMat, rhs::PetscVec, sol::PetscVec) end

@for_petsc function MatFactorSolveSchurComplement(petsclib::$UnionPetscLib, F::PetscMat, rhs::PetscVec, sol::PetscVec )

    @chk ccall(
               (:MatFactorSolveSchurComplement, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               F, rhs, sol,
              )


	return nothing
end 

"""
	MatFactorInvertSchurComplement(petsclib::PetscLibType,F::PetscMat) 
Invert the Schur complement matrix computed during the factorization step

Logically Collective

Input Parameter:
- `F` - the factored matrix obtained by calling `MatGetFactor()`

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatGetFactor()`, `MatFactorSetSchurIS()`, `MatFactorGetSchurComplement()`, `MatFactorCreateSchurComplement()`

# External Links
$(_doc_external("Mat/MatFactorInvertSchurComplement"))
"""
function MatFactorInvertSchurComplement(petsclib::PetscLibType, F::PetscMat) end

@for_petsc function MatFactorInvertSchurComplement(petsclib::$UnionPetscLib, F::PetscMat )

    @chk ccall(
               (:MatFactorInvertSchurComplement, $petsc_library),
               PetscErrorCode,
               (CMat,),
               F,
              )


	return nothing
end 

"""
	MatFactorFactorizeSchurComplement(petsclib::PetscLibType,F::PetscMat) 
Factorize the Schur complement matrix computed during the factorization step

Logically Collective

Input Parameter:
- `F` - the factored matrix obtained by calling `MatGetFactor()`

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatGetFactor()`, `MatFactorSetSchurIS()`, `MatFactorInvertSchurComplement()`

# External Links
$(_doc_external("Mat/MatFactorFactorizeSchurComplement"))
"""
function MatFactorFactorizeSchurComplement(petsclib::PetscLibType, F::PetscMat) end

@for_petsc function MatFactorFactorizeSchurComplement(petsclib::$UnionPetscLib, F::PetscMat )

    @chk ccall(
               (:MatFactorFactorizeSchurComplement, $petsc_library),
               PetscErrorCode,
               (CMat,),
               F,
              )


	return nothing
end 

"""
	MatPtAP(petsclib::PetscLibType,A::PetscMat, P::PetscMat, scall::MatReuse, fill::PetscReal, C::PetscMat) 
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
function MatPtAP(petsclib::PetscLibType, A::PetscMat, P::PetscMat, scall::MatReuse, fill::PetscReal, C::PetscMat) end

@for_petsc function MatPtAP(petsclib::$UnionPetscLib, A::PetscMat, P::PetscMat, scall::MatReuse, fill::$PetscReal, C::PetscMat )
	C_ = Ref(C.ptr)

    @chk ccall(
               (:MatPtAP, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, MatReuse, $PetscReal, Ptr{CMat}),
               A, P, scall, fill, C_,
              )

	C.ptr = C_NULL

	return nothing
end 

"""
	MatRARt(petsclib::PetscLibType,A::PetscMat, R::PetscMat, scall::MatReuse, fill::PetscReal, C::PetscMat) 
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
function MatRARt(petsclib::PetscLibType, A::PetscMat, R::PetscMat, scall::MatReuse, fill::PetscReal, C::PetscMat) end

@for_petsc function MatRARt(petsclib::$UnionPetscLib, A::PetscMat, R::PetscMat, scall::MatReuse, fill::$PetscReal, C::PetscMat )
	C_ = Ref(C.ptr)

    @chk ccall(
               (:MatRARt, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, MatReuse, $PetscReal, Ptr{CMat}),
               A, R, scall, fill, C_,
              )

	C.ptr = C_NULL

	return nothing
end 

"""
	MatMatMult(petsclib::PetscLibType,A::PetscMat, B::PetscMat, scall::MatReuse, fill::PetscReal, C::PetscMat) 
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
function MatMatMult(petsclib::PetscLibType, A::PetscMat, B::PetscMat, scall::MatReuse, fill::PetscReal, C::PetscMat) end

@for_petsc function MatMatMult(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat, scall::MatReuse, fill::$PetscReal, C::PetscMat )
	C_ = Ref(C.ptr)

    @chk ccall(
               (:MatMatMult, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, MatReuse, $PetscReal, Ptr{CMat}),
               A, B, scall, fill, C_,
              )

	C.ptr = C_NULL

	return nothing
end 

"""
	MatMatTransposeMult(petsclib::PetscLibType,A::PetscMat, B::PetscMat, scall::MatReuse, fill::PetscReal, C::PetscMat) 
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
function MatMatTransposeMult(petsclib::PetscLibType, A::PetscMat, B::PetscMat, scall::MatReuse, fill::PetscReal, C::PetscMat) end

@for_petsc function MatMatTransposeMult(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat, scall::MatReuse, fill::$PetscReal, C::PetscMat )
	C_ = Ref(C.ptr)

    @chk ccall(
               (:MatMatTransposeMult, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, MatReuse, $PetscReal, Ptr{CMat}),
               A, B, scall, fill, C_,
              )

	C.ptr = C_NULL

	return nothing
end 

"""
	MatTransposeMatMult(petsclib::PetscLibType,A::PetscMat, B::PetscMat, scall::MatReuse, fill::PetscReal, C::PetscMat) 
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
function MatTransposeMatMult(petsclib::PetscLibType, A::PetscMat, B::PetscMat, scall::MatReuse, fill::PetscReal, C::PetscMat) end

@for_petsc function MatTransposeMatMult(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat, scall::MatReuse, fill::$PetscReal, C::PetscMat )
	C_ = Ref(C.ptr)

    @chk ccall(
               (:MatTransposeMatMult, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, MatReuse, $PetscReal, Ptr{CMat}),
               A, B, scall, fill, C_,
              )

	C.ptr = C_NULL

	return nothing
end 

"""
	MatMatMatMult(petsclib::PetscLibType,A::PetscMat, B::PetscMat, C::PetscMat, scall::MatReuse, fill::PetscReal, D::PetscMat) 
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
function MatMatMatMult(petsclib::PetscLibType, A::PetscMat, B::PetscMat, C::PetscMat, scall::MatReuse, fill::PetscReal, D::PetscMat) end

@for_petsc function MatMatMatMult(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat, C::PetscMat, scall::MatReuse, fill::$PetscReal, D::PetscMat )
	D_ = Ref(D.ptr)

    @chk ccall(
               (:MatMatMatMult, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat, MatReuse, $PetscReal, Ptr{CMat}),
               A, B, C, scall, fill, D_,
              )

	D.ptr = C_NULL

	return nothing
end 

"""
	matredundant::PetscMat = MatCreateRedundantMatrix(petsclib::PetscLibType,mat::PetscMat, nsubcomm::PetscInt, subcomm::MPI_Comm, reuse::MatReuse) 
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
function MatCreateRedundantMatrix(petsclib::PetscLibType, mat::PetscMat, nsubcomm::PetscInt, subcomm::MPI_Comm, reuse::MatReuse) end

@for_petsc function MatCreateRedundantMatrix(petsclib::$UnionPetscLib, mat::PetscMat, nsubcomm::$PetscInt, subcomm::MPI_Comm, reuse::MatReuse )
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
	MatGetMultiProcBlock(petsclib::PetscLibType,mat::PetscMat, subComm::MPI_Comm, scall::MatReuse, subMat::PetscMat) 
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
function MatGetMultiProcBlock(petsclib::PetscLibType, mat::PetscMat, subComm::MPI_Comm, scall::MatReuse, subMat::PetscMat) end

@for_petsc function MatGetMultiProcBlock(petsclib::$UnionPetscLib, mat::PetscMat, subComm::MPI_Comm, scall::MatReuse, subMat::PetscMat )
	subMat_ = Ref(subMat.ptr)

    @chk ccall(
               (:MatGetMultiProcBlock, $petsc_library),
               PetscErrorCode,
               (CMat, MPI_Comm, MatReuse, Ptr{CMat}),
               mat, subComm, scall, subMat_,
              )

	subMat.ptr = C_NULL

	return nothing
end 

"""
	MatGetLocalSubMatrix(petsclib::PetscLibType,mat::PetscMat, isrow::IS, iscol::IS, submat::PetscMat) 
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
function MatGetLocalSubMatrix(petsclib::PetscLibType, mat::PetscMat, isrow::IS, iscol::IS, submat::PetscMat) end

@for_petsc function MatGetLocalSubMatrix(petsclib::$UnionPetscLib, mat::PetscMat, isrow::IS, iscol::IS, submat::PetscMat )
	submat_ = Ref(submat.ptr)

    @chk ccall(
               (:MatGetLocalSubMatrix, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, CIS, Ptr{CMat}),
               mat, isrow, iscol, submat_,
              )

	submat.ptr = C_NULL

	return nothing
end 

"""
	MatRestoreLocalSubMatrix(petsclib::PetscLibType,mat::PetscMat, isrow::IS, iscol::IS, submat::PetscMat) 
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
function MatRestoreLocalSubMatrix(petsclib::PetscLibType, mat::PetscMat, isrow::IS, iscol::IS, submat::PetscMat) end

@for_petsc function MatRestoreLocalSubMatrix(petsclib::$UnionPetscLib, mat::PetscMat, isrow::IS, iscol::IS, submat::PetscMat )
	submat_ = Ref(submat.ptr)

    @chk ccall(
               (:MatRestoreLocalSubMatrix, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, CIS, Ptr{CMat}),
               mat, isrow, iscol, submat_,
              )

	submat.ptr = C_NULL

	return nothing
end 

"""
	MatFindZeroDiagonals(petsclib::PetscLibType,mat::PetscMat, is::IS) 
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
function MatFindZeroDiagonals(petsclib::PetscLibType, mat::PetscMat, is::IS) end

@for_petsc function MatFindZeroDiagonals(petsclib::$UnionPetscLib, mat::PetscMat, is::IS )

    @chk ccall(
               (:MatFindZeroDiagonals, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CIS}),
               mat, is,
              )


	return nothing
end 

"""
	MatFindOffBlockDiagonalEntries(petsclib::PetscLibType,mat::PetscMat, is::IS) 
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
function MatFindOffBlockDiagonalEntries(petsclib::PetscLibType, mat::PetscMat, is::IS) end

@for_petsc function MatFindOffBlockDiagonalEntries(petsclib::$UnionPetscLib, mat::PetscMat, is::IS )

    @chk ccall(
               (:MatFindOffBlockDiagonalEntries, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CIS}),
               mat, is,
              )


	return nothing
end 

"""
	values::Vector{PetscScalar} = MatInvertBlockDiagonal(petsclib::PetscLibType,mat::PetscMat) 
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
function MatInvertBlockDiagonal(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatInvertBlockDiagonal(petsclib::$UnionPetscLib, mat::PetscMat )
	values_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatInvertBlockDiagonal, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               mat, values_,
              )

	values = unsafe_wrap(Array, values_[], VecGetLocalSize(petsclib, x); own = false)

	return values
end 

"""
	MatInvertVariableBlockDiagonal(petsclib::PetscLibType,mat::PetscMat, nblocks::PetscInt, bsizes::Vector{PetscInt}, values::Vector{PetscScalar}) 
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
function MatInvertVariableBlockDiagonal(petsclib::PetscLibType, mat::PetscMat, nblocks::PetscInt, bsizes::Vector{PetscInt}, values::Vector{PetscScalar}) end

@for_petsc function MatInvertVariableBlockDiagonal(petsclib::$UnionPetscLib, mat::PetscMat, nblocks::$PetscInt, bsizes::Vector{$PetscInt}, values::Vector{$PetscScalar} )

    @chk ccall(
               (:MatInvertVariableBlockDiagonal, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}),
               mat, nblocks, bsizes, values,
              )


	return nothing
end 

"""
	MatInvertBlockDiagonalMat(petsclib::PetscLibType,A::PetscMat, C::PetscMat) 
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
function MatInvertBlockDiagonalMat(petsclib::PetscLibType, A::PetscMat, C::PetscMat) end

@for_petsc function MatInvertBlockDiagonalMat(petsclib::$UnionPetscLib, A::PetscMat, C::PetscMat )

    @chk ccall(
               (:MatInvertBlockDiagonalMat, $petsc_library),
               PetscErrorCode,
               (CMat, CMat),
               A, C,
              )


	return nothing
end 

"""
	MatTransColoringApplySpToDen(petsclib::PetscLibType,coloring::MatTransposeColoring, B::PetscMat, Btdense::PetscMat) 
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
function MatTransColoringApplySpToDen(petsclib::PetscLibType, coloring::MatTransposeColoring, B::PetscMat, Btdense::PetscMat) end

@for_petsc function MatTransColoringApplySpToDen(petsclib::$UnionPetscLib, coloring::MatTransposeColoring, B::PetscMat, Btdense::PetscMat )

    @chk ccall(
               (:MatTransColoringApplySpToDen, $petsc_library),
               PetscErrorCode,
               (MatTransposeColoring, CMat, CMat),
               coloring, B, Btdense,
              )


	return nothing
end 

"""
	MatTransColoringApplyDenToSp(petsclib::PetscLibType,matcoloring::MatTransposeColoring, Cden::PetscMat, Csp::PetscMat) 
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
function MatTransColoringApplyDenToSp(petsclib::PetscLibType, matcoloring::MatTransposeColoring, Cden::PetscMat, Csp::PetscMat) end

@for_petsc function MatTransColoringApplyDenToSp(petsclib::$UnionPetscLib, matcoloring::MatTransposeColoring, Cden::PetscMat, Csp::PetscMat )

    @chk ccall(
               (:MatTransColoringApplyDenToSp, $petsc_library),
               PetscErrorCode,
               (MatTransposeColoring, CMat, CMat),
               matcoloring, Cden, Csp,
              )


	return nothing
end 

"""
	MatGetNonzeroState(petsclib::PetscLibType,mat::PetscMat, state::PetscObjectState) 
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
function MatGetNonzeroState(petsclib::PetscLibType, mat::PetscMat, state::PetscObjectState) end

@for_petsc function MatGetNonzeroState(petsclib::$UnionPetscLib, mat::PetscMat, state::PetscObjectState )

    @chk ccall(
               (:MatGetNonzeroState, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscObjectState}),
               mat, state,
              )


	return nothing
end 

"""
	mpimat::PetscMat = MatCreateMPIMatConcatenateSeqMat(petsclib::PetscLibType,comm::MPI_Comm, seqmat::PetscMat, n::PetscInt, reuse::MatReuse) 
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
function MatCreateMPIMatConcatenateSeqMat(petsclib::PetscLibType, comm::MPI_Comm, seqmat::PetscMat, n::PetscInt, reuse::MatReuse) end

@for_petsc function MatCreateMPIMatConcatenateSeqMat(petsclib::$UnionPetscLib, comm::MPI_Comm, seqmat::PetscMat, n::$PetscInt, reuse::MatReuse )
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
	n::PetscInt,iss::Vector{IS} = MatSubdomainsCreateCoalesce(petsclib::PetscLibType,A::PetscMat, N::PetscInt) 
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
function MatSubdomainsCreateCoalesce(petsclib::PetscLibType, A::PetscMat, N::PetscInt) end

@for_petsc function MatSubdomainsCreateCoalesce(petsclib::$UnionPetscLib, A::PetscMat, N::$PetscInt )
	n_ = Ref{$PetscInt}()
	iss_ = Ref{Ptr{CIS}}()

    @chk ccall(
               (:MatSubdomainsCreateCoalesce, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, Ptr{Ptr{CIS}}),
               A, N, n_, iss_,
              )

	n = n_[]
	iss = unsafe_wrap(Array, iss_[], VecGetLocalSize(petsclib, x); own = false)

	return n,iss
end 

"""
	MatGalerkin(petsclib::PetscLibType,restrct::PetscMat, dA::PetscMat, interpolate::PetscMat, reuse::MatReuse, fill::PetscReal, A::PetscMat) 
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
function MatGalerkin(petsclib::PetscLibType, restrct::PetscMat, dA::PetscMat, interpolate::PetscMat, reuse::MatReuse, fill::PetscReal, A::PetscMat) end

@for_petsc function MatGalerkin(petsclib::$UnionPetscLib, restrct::PetscMat, dA::PetscMat, interpolate::PetscMat, reuse::MatReuse, fill::$PetscReal, A::PetscMat )
	A_ = Ref(A.ptr)

    @chk ccall(
               (:MatGalerkin, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat, MatReuse, $PetscReal, Ptr{CMat}),
               restrct, dA, interpolate, reuse, fill, A_,
              )

	A.ptr = C_NULL

	return nothing
end 

"""
	has::PetscBool = MatHasOperation(petsclib::PetscLibType,mat::PetscMat, op::MatOperation) 
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
function MatHasOperation(petsclib::PetscLibType, mat::PetscMat, op::MatOperation) end

@for_petsc function MatHasOperation(petsclib::$UnionPetscLib, mat::PetscMat, op::MatOperation )
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
	cong::PetscBool = MatHasCongruentLayouts(petsclib::PetscLibType,mat::PetscMat) 
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
function MatHasCongruentLayouts(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatHasCongruentLayouts(petsclib::$UnionPetscLib, mat::PetscMat )
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
	MatSetInf(petsclib::PetscLibType,A::PetscMat) 

# External Links
$(_doc_external("Mat/MatSetInf"))
"""
function MatSetInf(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatSetInf(petsclib::$UnionPetscLib, A::PetscMat )

    @chk ccall(
               (:MatSetInf, $petsc_library),
               PetscErrorCode,
               (CMat,),
               A,
              )


	return nothing
end 

"""
	graph::PetscMat = MatCreateGraph(petsclib::PetscLibType,A::PetscMat, sym::PetscBool, scale::PetscBool, filter::PetscReal, num_idx::PetscInt, index::Vector{PetscInt}) 
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
function MatCreateGraph(petsclib::PetscLibType, A::PetscMat, sym::PetscBool, scale::PetscBool, filter::PetscReal, num_idx::PetscInt, index::Vector{PetscInt}) end

@for_petsc function MatCreateGraph(petsclib::$UnionPetscLib, A::PetscMat, sym::PetscBool, scale::PetscBool, filter::$PetscReal, num_idx::$PetscInt, index::Vector{$PetscInt} )
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
	MatEliminateZeros(petsclib::PetscLibType,A::PetscMat, keep::PetscBool) 
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
function MatEliminateZeros(petsclib::PetscLibType, A::PetscMat, keep::PetscBool) end

@for_petsc function MatEliminateZeros(petsclib::$UnionPetscLib, A::PetscMat, keep::PetscBool )

    @chk ccall(
               (:MatEliminateZeros, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               A, keep,
              )


	return nothing
end 

"""
	m::PetscMemType = MatGetCurrentMemType(petsclib::PetscLibType,A::PetscMat) 
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
function MatGetCurrentMemType(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatGetCurrentMemType(petsclib::$UnionPetscLib, A::PetscMat )
	m_ = Ref{PetscMemType}()

    @chk ccall(
               (:MatGetCurrentMemType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscMemType}),
               A, m_,
              )

	m = unsafe_string(m_[])

	return m
end 

"""
	MatScaLAPACKSetBlockSizes(petsclib::PetscLibType,A::PetscMat, mb::PetscInt, nb::PetscInt) 
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
function MatScaLAPACKSetBlockSizes(petsclib::PetscLibType, A::PetscMat, mb::PetscInt, nb::PetscInt) end

@for_petsc function MatScaLAPACKSetBlockSizes(petsclib::$UnionPetscLib, A::PetscMat, mb::$PetscInt, nb::$PetscInt )

    @chk ccall(
               (:MatScaLAPACKSetBlockSizes, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt),
               A, mb, nb,
              )


	return nothing
end 

"""
	mb::PetscInt,nb::PetscInt = MatScaLAPACKGetBlockSizes(petsclib::PetscLibType,A::PetscMat) 
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
function MatScaLAPACKGetBlockSizes(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatScaLAPACKGetBlockSizes(petsclib::$UnionPetscLib, A::PetscMat )
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
	MatMPISELLSetPreallocation(petsclib::PetscLibType,B::PetscMat, d_nz::PetscInt, d_nnz::Vector{PetscInt}, o_nz::PetscInt, o_nnz::Vector{PetscInt}) 
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
function MatMPISELLSetPreallocation(petsclib::PetscLibType, B::PetscMat, d_nz::PetscInt, d_nnz::Vector{PetscInt}, o_nz::PetscInt, o_nnz::Vector{PetscInt}) end

@for_petsc function MatMPISELLSetPreallocation(petsclib::$UnionPetscLib, B::PetscMat, d_nz::$PetscInt, d_nnz::Vector{$PetscInt}, o_nz::$PetscInt, o_nnz::Vector{$PetscInt} )

    @chk ccall(
               (:MatMPISELLSetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}),
               B, d_nz, d_nnz, o_nz, o_nnz,
              )


	return nothing
end 

"""
	A::PetscMat = MatCreateSELL(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_rlenmax::PetscInt, d_rlen::Vector{PetscInt}, o_rlenmax::PetscInt, o_rlen::Vector{PetscInt}) 
Creates a sparse parallel matrix in `MATSELL` format.

Collective

Input Parameters:
- `comm`      - MPI communicator
- `m`         - number of local rows (or `PETSC_DECIDE` to have calculated if M is given). This value should be the same as the local size used in creating the y vector for the matrix-vector product y = Ax.
- `n`         - This value should be the same as the local size used in creating the x vector for the matrix-vector product y = Ax. (or `PETSC_DECIDE` to have calculated if `N` is given) For square matrices n is almost always `m`.
- `M`         - number of global rows (or `PETSC_DETERMINE` to have calculated if `m` is given)
- `N`         - number of global columns (or `PETSC_DETERMINE` to have calculated if `n` is given)
- `d_rlenmax` - max number of nonzeros per row in DIAGONAL portion of local submatrix (same value is used for all local rows)
- `d_rlen`    - array containing the number of nonzeros in the various rows of the DIAGONAL portion of the local submatrix (possibly different for each row) or `NULL`, if d_rlenmax is used to specify the nonzero structure. The size of this array is equal to the number of local rows, i.e `m`.
- `o_rlenmax` - max number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix (same value is used for all local rows).
- `o_rlen`    - array containing the number of nonzeros in the various rows of the OFF-DIAGONAL portion of the local submatrix (possibly different for each row) or `NULL`, if `o_rlenmax` is used to specify the nonzero structure. The size of this array is equal to the number of local rows, i.e `m`.

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
function MatCreateSELL(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_rlenmax::PetscInt, d_rlen::Union{Ptr,Vector{PetscInt}}, o_rlenmax::PetscInt, o_rlen::Union{Ptr,Vector{PetscInt}}) end

@for_petsc function MatCreateSELL(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, d_rlenmax::$PetscInt, d_rlen::Union{Ptr,Vector{$PetscInt}}, o_rlenmax::$PetscInt, o_rlen::Union{Ptr,Vector{$PetscInt}} )
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
	colmap::Vector{PetscInt} = MatMPISELLGetSeqSELL(petsclib::PetscLibType,A::PetscMat, Ad::PetscMat, Ao::PetscMat) 
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
function MatMPISELLGetSeqSELL(petsclib::PetscLibType, A::PetscMat, Ad::PetscMat, Ao::PetscMat) end

@for_petsc function MatMPISELLGetSeqSELL(petsclib::$UnionPetscLib, A::PetscMat, Ad::PetscMat, Ao::PetscMat )
	Ad_ = Ref(Ad.ptr)
	Ao_ = Ref(Ao.ptr)
	colmap_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:MatMPISELLGetSeqSELL, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}, Ptr{CMat}, Ptr{Ptr{$PetscInt}}),
               A, Ad_, Ao_, colmap_,
              )

	Ad.ptr = C_NULL
	Ao.ptr = C_NULL
	colmap = unsafe_wrap(Array, colmap_[], VecGetLocalSize(petsclib, x); own = false)

	return colmap
end 

"""
	MatMPISELLGetLocalMatCondensed(petsclib::PetscLibType,A::PetscMat, scall::MatReuse, row::IS, col::IS, A_loc::PetscMat) 
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
function MatMPISELLGetLocalMatCondensed(petsclib::PetscLibType, A::PetscMat, scall::MatReuse, row::Union{Ptr,IS}, col::Union{Ptr,IS}, A_loc::PetscMat) end

@for_petsc function MatMPISELLGetLocalMatCondensed(petsclib::$UnionPetscLib, A::PetscMat, scall::MatReuse, row::Union{Ptr,IS}, col::Union{Ptr,IS}, A_loc::PetscMat )
	A_loc_ = Ref(A_loc.ptr)

    @chk ccall(
               (:MatMPISELLGetLocalMatCondensed, $petsc_library),
               PetscErrorCode,
               (CMat, MatReuse, Ptr{CIS}, Ptr{CIS}, Ptr{CMat}),
               A, scall, row, col, A_loc_,
              )

	A_loc.ptr = C_NULL

	return nothing
end 

"""
	MatSeqSELLSetPreallocation(petsclib::PetscLibType,B::PetscMat, rlenmax::PetscInt, rlen::Vector{PetscInt}) 
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
function MatSeqSELLSetPreallocation(petsclib::PetscLibType, B::PetscMat, rlenmax::PetscInt, rlen::Vector{PetscInt}) end

@for_petsc function MatSeqSELLSetPreallocation(petsclib::$UnionPetscLib, B::PetscMat, rlenmax::$PetscInt, rlen::Vector{$PetscInt} )

    @chk ccall(
               (:MatSeqSELLSetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}),
               B, rlenmax, rlen,
              )


	return nothing
end 

"""
	ratio::PetscReal = MatSeqSELLGetFillRatio(petsclib::PetscLibType,A::PetscMat) 
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
function MatSeqSELLGetFillRatio(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatSeqSELLGetFillRatio(petsclib::$UnionPetscLib, A::PetscMat )
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
	slicewidth::PetscInt = MatSeqSELLGetMaxSliceWidth(petsclib::PetscLibType,A::PetscMat) 
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
function MatSeqSELLGetMaxSliceWidth(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatSeqSELLGetMaxSliceWidth(petsclib::$UnionPetscLib, A::PetscMat )
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
	slicewidth::PetscReal = MatSeqSELLGetAvgSliceWidth(petsclib::PetscLibType,A::PetscMat) 
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
function MatSeqSELLGetAvgSliceWidth(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatSeqSELLGetAvgSliceWidth(petsclib::$UnionPetscLib, A::PetscMat )
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
	MatSeqSELLSetSliceHeight(petsclib::PetscLibType,A::PetscMat, sliceheight::PetscInt) 
sets the slice height.

Not Collective

Input Parameters:
- `A`           - a MATSEQSELL matrix
- `sliceheight` - slice height

-seealso: `MATSEQSELL`, `MatSeqSELLGetVarSliceSize()`

# External Links
$(_doc_external("Mat/MatSeqSELLSetSliceHeight"))
"""
function MatSeqSELLSetSliceHeight(petsclib::PetscLibType, A::PetscMat, sliceheight::PetscInt) end

@for_petsc function MatSeqSELLSetSliceHeight(petsclib::$UnionPetscLib, A::PetscMat, sliceheight::$PetscInt )

    @chk ccall(
               (:MatSeqSELLSetSliceHeight, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt),
               A, sliceheight,
              )


	return nothing
end 

"""
	variance::PetscReal = MatSeqSELLGetVarSliceSize(petsclib::PetscLibType,A::PetscMat) 
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
function MatSeqSELLGetVarSliceSize(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatSeqSELLGetVarSliceSize(petsclib::$UnionPetscLib, A::PetscMat )
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
	A::PetscMat = MatCreateSeqSELL(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, rlenmax::PetscInt, rlen::Union{Ptr,Vector{PetscInt}}) 
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
	MatPreallocatorPreallocate(petsclib::PetscLibType,mat::PetscMat, fill::PetscBool, A::PetscMat) 
Preallocates the A matrix, using information from a `MATPREALLOCATOR` mat, optionally filling A with zeros

Input Parameters:
- `mat`  - the `MATPREALLOCATOR` preallocator matrix
- `fill` - fill the matrix with zeros
- `A`    - the matrix to be preallocated

-seealso: `MATPREALLOCATOR`, `MatXAIJSetPreallocation()`

# External Links
$(_doc_external("Mat/MatPreallocatorPreallocate"))
"""
function MatPreallocatorPreallocate(petsclib::PetscLibType, mat::PetscMat, fill::PetscBool, A::PetscMat) end

@for_petsc function MatPreallocatorPreallocate(petsclib::$UnionPetscLib, mat::PetscMat, fill::PetscBool, A::PetscMat )

    @chk ccall(
               (:MatPreallocatorPreallocate, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool, CMat),
               mat, fill, A,
              )


	return nothing
end 

"""
	J::PetscMat = MatCreateConstantDiagonal(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, diag::PetscScalar) 
Creates a matrix with a uniform value along the diagonal

Collective

Input Parameters:
- `comm` - MPI communicator
- `m`    - number of local rows (or `PETSC_DECIDE` to have calculated if `M` is given). This value should be the same as the local size used in creating the y vector for the matrix-vector product y = Ax.
- `n`    - This value should be the same as the local size used in creating the x vector for the matrix-vector product y = Ax. (or `PETSC_DECIDE` to have calculated if `N` is given) For square matrices n is almost always `m`.
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
	value::PetscScalar = MatConstantDiagonalGetConstant(petsclib::PetscLibType,mat::PetscMat) 
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
function MatConstantDiagonalGetConstant(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatConstantDiagonalGetConstant(petsclib::$UnionPetscLib, mat::PetscMat )
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
	MatMPIBAIJSetPreallocationCSR(petsclib::PetscLibType,B::PetscMat, bs::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, v::Vector{PetscScalar}) 
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
function MatMPIBAIJSetPreallocationCSR(petsclib::PetscLibType, B::PetscMat, bs::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, v::Union{Ptr,Vector{PetscScalar}}) end

@for_petsc function MatMPIBAIJSetPreallocationCSR(petsclib::$UnionPetscLib, B::PetscMat, bs::$PetscInt, i::Vector{$PetscInt}, j::Vector{$PetscInt}, v::Union{Ptr,Vector{$PetscScalar}})

    @chk ccall(
               (:MatMPIBAIJSetPreallocationCSR, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}),
               B, bs, i, j, v,
              )


	return nothing
end 

"""
	MatMPIBAIJSetPreallocation(petsclib::PetscLibType,B::PetscMat, bs::PetscInt, d_nz::PetscInt, d_nnz::Vector{PetscInt}, o_nz::PetscInt, o_nnz::Vector{PetscInt}) 
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
function MatMPIBAIJSetPreallocation(petsclib::PetscLibType, B::PetscMat, bs::PetscInt, d_nz::PetscInt, d_nnz::Union{Ptr,Vector{PetscInt}}, o_nz::PetscInt, o_nnz::Union{Ptr,Vector{PetscInt}}) end

@for_petsc function MatMPIBAIJSetPreallocation(petsclib::$UnionPetscLib, B::PetscMat, bs::$PetscInt, d_nz::$PetscInt, d_nnz::Union{Ptr,Vector{$PetscInt}}, o_nz::$PetscInt, o_nnz::Union{Ptr,Vector{$PetscInt}})

    @chk ccall(
               (:MatMPIBAIJSetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}),
               B, bs, d_nz, d_nnz, o_nz, o_nnz,
              )


	return nothing
end 

"""
	A::PetscMat = MatCreateBAIJ(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_nz::PetscInt, d_nnz::Union{Ptr,Vector{PetscInt}}, o_nz::PetscInt, o_nnz::Union{Ptr,Vector{PetscInt}}) 
Creates a sparse parallel matrix in `MATBAIJ` format
(block compressed row).

Collective

Input Parameters:
- `comm`  - MPI communicator
- `bs`    - size of block, the blocks are ALWAYS square. One can use `MatSetBlockSizes()` to set a different row and column blocksize but the row blocksize always defines the size of the blocks. The column blocksize sets the blocksize of the vectors obtained with `MatCreateVecs()`
- `m`     - number of local rows (or `PETSC_DECIDE` to have calculated if M is given). This value should be the same as the local size used in creating the y vector for the matrix-vector product y = Ax.
- `n`     - number of local columns (or `PETSC_DECIDE` to have calculated if N is given). This value should be the same as the local size used in creating the x vector for the matrix-vector product y = Ax.
- `M`     - number of global rows (or `PETSC_DETERMINE` to have calculated if m is given)
- `N`     - number of global columns (or `PETSC_DETERMINE` to have calculated if n is given)
- `d_nz`  - number of nonzero blocks per block row in diagonal portion of local submatrix  (same for all local rows)
- `d_nnz` - array containing the number of nonzero blocks in the various block rows of the in diagonal portion of the local (possibly different for each block row) or NULL.  If you plan to factor the matrix you must leave room for the diagonal entry and set it even if it is zero.
- `o_nz`  - number of nonzero blocks per block row in the off-diagonal portion of local submatrix (same for all local rows).
- `o_nnz` - array containing the number of nonzero blocks in the various block rows of the off-diagonal portion of the local submatrix (possibly different for each block row) or NULL.

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
function MatCreateBAIJ(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_nz::PetscInt, d_nnz::Union{Ptr,Vector{PetscInt}}, o_nz::PetscInt, o_nnz::Union{Ptr,Vector{PetscInt}}) end

@for_petsc function MatCreateBAIJ(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, d_nz::$PetscInt, d_nnz::Union{Ptr,Vector{$PetscInt}}, o_nz::$PetscInt, o_nnz::Union{Ptr,Vector{$PetscInt}} )
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
	MatMPIBAIJSetHashTableFactor(petsclib::PetscLibType,mat::PetscMat, fact::PetscReal) 
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
function MatMPIBAIJSetHashTableFactor(petsclib::PetscLibType, mat::PetscMat, fact::PetscReal) end

@for_petsc function MatMPIBAIJSetHashTableFactor(petsclib::$UnionPetscLib, mat::PetscMat, fact::$PetscReal )

    @chk ccall(
               (:MatMPIBAIJSetHashTableFactor, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscReal),
               mat, fact,
              )


	return nothing
end 

"""
	colmap::Vector{PetscInt} = MatMPIBAIJGetSeqBAIJ(petsclib::PetscLibType,A::PetscMat, Ad::PetscMat, Ao::PetscMat) 

# External Links
$(_doc_external("Mat/MatMPIBAIJGetSeqBAIJ"))
"""
function MatMPIBAIJGetSeqBAIJ(petsclib::PetscLibType, A::PetscMat, Ad::PetscMat, Ao::PetscMat) end

@for_petsc function MatMPIBAIJGetSeqBAIJ(petsclib::$UnionPetscLib, A::PetscMat, Ad::PetscMat, Ao::PetscMat )
	Ad_ = Ref(Ad.ptr)
	Ao_ = Ref(Ao.ptr)
	colmap_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:MatMPIBAIJGetSeqBAIJ, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}, Ptr{CMat}, Ptr{Ptr{$PetscInt}}),
               A, Ad_, Ao_, colmap_,
              )

	Ad.ptr = C_NULL
	Ao.ptr = C_NULL
	colmap = unsafe_wrap(Array, colmap_[], VecGetLocalSize(petsclib, x); own = false)

	return colmap
end 

"""
	mat::PetscMat = MatCreateMPIBAIJWithArrays(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, a::Vector{PetscScalar}) 
creates a `MATMPIBAIJ` matrix using arrays that contain in standard block CSR format for the local rows.

Collective

Input Parameters:
- `comm` - MPI communicator
- `bs`   - the block size, only a block size of 1 is supported
- `m`    - number of local rows (Cannot be `PETSC_DECIDE`)
- `n`    - This value should be the same as the local size used in creating the x vector for the matrix-vector product  y = Ax . (or `PETSC_DECIDE` to have calculated if `N` is given) For square matrices `n` is almost always `m`.
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
	A::PetscMat = MatCreateBAIJMKL(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_nz::PetscInt, d_nnz::Vector{PetscInt}, o_nz::PetscInt, o_nnz::Vector{PetscInt}) 
Creates a sparse parallel matrix in `MATBAIJMKL` format (block compressed row).

Collective

Input Parameters:
- `comm`  - MPI communicator
- `bs`    - size of block, the blocks are ALWAYS square. One can use `MatSetBlockSizes()` to set a different row and column blocksize but the row blocksize always defines the size of the blocks. The column blocksize sets the blocksize of the vectors obtained with `MatCreateVecs()`
- `m`     - number of local rows (or `PETSC_DECIDE` to have calculated if `M` is given) This value should be the same as the local size used in creating the y vector for the matrix-vector product y = Ax.
- `n`     - number of local columns (or `PETSC_DECIDE` to have calculated if `N` is given) This value should be the same as the local size used in creating the x vector for the matrix-vector product y = Ax.
- `M`     - number of global rows (or `PETSC_DETERMINE` to have calculated if `m` is given)
- `N`     - number of global columns (or `PETSC_DETERMINE` to have calculated if `n` is given)
- `d_nz`  - number of nonzero blocks per block row in diagonal portion of local submatrix  (same for all local rows)
- `d_nnz` - array containing the number of nonzero blocks in the various block rows of the in diagonal portion of the local (possibly different for each block row) or `NULL`.  If you plan to factor the matrix you must leave room for the diagonal entry and set it even if it is zero.
- `o_nz`  - number of nonzero blocks per block row in the off-diagonal portion of local submatrix (same for all local rows).
- `o_nnz` - array containing the number of nonzero blocks in the various block rows of the off-diagonal portion of the local submatrix (possibly different for each block row) or `NULL`.

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
function MatCreateBAIJMKL(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_nz::PetscInt, d_nnz::Union{Ptr,Vector{PetscInt}}, o_nz::PetscInt, o_nnz::Union{Ptr,Vector{PetscInt}}) end

@for_petsc function MatCreateBAIJMKL(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, d_nz::$PetscInt, d_nnz::Union{Ptr,Vector{$PetscInt}}, o_nz::$PetscInt, o_nnz::Union{Ptr,Vector{$PetscInt}} )
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
	indices::PetscInt = MatSeqBAIJSetColumnIndices(petsclib::PetscLibType,mat::PetscMat) 
Set the column indices for all the block rows in the matrix.

Input Parameters:
- `mat`     - the `MATSEQBAIJ` matrix
- `indices` - the block column indices

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATSEQBAIJ`, `MatSetValues()`

# External Links
$(_doc_external("Mat/MatSeqBAIJSetColumnIndices"))
"""
function MatSeqBAIJSetColumnIndices(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatSeqBAIJSetColumnIndices(petsclib::$UnionPetscLib, mat::PetscMat )
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
	array::Vector{PetscScalar} = MatSeqBAIJGetArray(petsclib::PetscLibType,A::PetscMat) 
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
function MatSeqBAIJGetArray(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatSeqBAIJGetArray(petsclib::$UnionPetscLib, A::PetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatSeqBAIJGetArray, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

	array = unsafe_wrap(Array, array_[], VecGetLocalSize(petsclib, x); own = false)

	return array
end 

"""
	array::Vector{PetscScalar} = MatSeqBAIJRestoreArray(petsclib::PetscLibType,A::PetscMat) 
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
function MatSeqBAIJRestoreArray(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatSeqBAIJRestoreArray(petsclib::$UnionPetscLib, A::PetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatSeqBAIJRestoreArray, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

	array = unsafe_wrap(Array, array_[], VecGetLocalSize(petsclib, x); own = false)

	return array
end 

"""
	A::PetscMat = MatCreateSeqBAIJ(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Vector{PetscInt}) 
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
function MatCreateSeqBAIJ(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Union{Ptr,Vector{PetscInt}}) end

@for_petsc function MatCreateSeqBAIJ(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, m::$PetscInt, n::$PetscInt, nz::$PetscInt, nnz::Union{Ptr,Vector{$PetscInt}} )
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
	MatSeqBAIJSetPreallocation(petsclib::PetscLibType,B::PetscMat, bs::PetscInt, nz::PetscInt, nnz::Vector{PetscInt}) 
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
function MatSeqBAIJSetPreallocation(petsclib::PetscLibType, B::PetscMat, bs::PetscInt, nz::PetscInt, nnz::Union{Ptr,Vector{PetscInt}}) end

@for_petsc function MatSeqBAIJSetPreallocation(petsclib::$UnionPetscLib, B::PetscMat, bs::$PetscInt, nz::$PetscInt, nnz::Union{Ptr,Vector{$PetscInt}})

    @chk ccall(
               (:MatSeqBAIJSetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, Ptr{$PetscInt}),
               B, bs, nz, nnz,
              )


	return nothing
end 

"""
	MatSeqBAIJSetPreallocationCSR(petsclib::PetscLibType,B::PetscMat, bs::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, v::Vector{PetscScalar}) 
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
function MatSeqBAIJSetPreallocationCSR(petsclib::PetscLibType, B::PetscMat, bs::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, v::Union{Ptr,Vector{PetscScalar}}) end

@for_petsc function MatSeqBAIJSetPreallocationCSR(petsclib::$UnionPetscLib, B::PetscMat, bs::$PetscInt, i::Vector{$PetscInt}, j::Vector{$PetscInt}, v::Union{Ptr,Vector{$PetscScalar}})

    @chk ccall(
               (:MatSeqBAIJSetPreallocationCSR, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}),
               B, bs, i, j, v,
              )


	return nothing
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
	A::PetscMat = MatCreateSeqBAIJMKL(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Vector{PetscInt}) 
Creates a sparse matrix of type `MATSEQBAIJMKL`.
This type inherits from `MATSEQBAIJ` and is largely identical, but uses sparse BLAS
routines from Intel MKL whenever possible.

Input Parameters:
- `comm` - MPI communicator, set to `PETSC_COMM_SELF`
- `bs`   - size of block, the blocks are ALWAYS square. One can use `MatSetBlockSizes()` to set a different row and column blocksize but the row blocksize always defines the size of the blocks. The column blocksize sets the blocksize of the vectors obtained with `MatCreateVecs()`
- `m`    - number of rows
- `n`    - number of columns
- `nz`   - number of nonzero blocks  per block row (same for all rows)
- `nnz`  - array containing the number of nonzero blocks in the various block rows (possibly different for each block row) or `NULL`

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
function MatCreateSeqBAIJMKL(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Union{Ptr,Vector{PetscInt}}) end

@for_petsc function MatCreateSeqBAIJMKL(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, m::$PetscInt, n::$PetscInt, nz::$PetscInt, nnz::Union{Ptr,Vector{$PetscInt}} )
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
	MatPythonSetType(petsclib::PetscLibType,mat::PetscMat, pyname::Vector{Cchar}) 
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
function MatPythonSetType(petsclib::PetscLibType, mat::PetscMat, pyname::Vector{Cchar}) end

@for_petsc function MatPythonSetType(petsclib::$UnionPetscLib, mat::PetscMat, pyname::Vector{Cchar} )

    @chk ccall(
               (:MatPythonSetType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Cchar}),
               mat, pyname,
              )


	return nothing
end 

"""
	pyname::Vector{Cchar} = MatPythonGetType(petsclib::PetscLibType,mat::PetscMat) 
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
function MatPythonGetType(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatPythonGetType(petsclib::$UnionPetscLib, mat::PetscMat )
	pyname_ = Ref{Ptr{Cchar}}()

    @chk ccall(
               (:MatPythonGetType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{Cchar}}),
               mat, pyname_,
              )

	pyname = unsafe_wrap(Array, pyname_[], VecGetLocalSize(petsclib, x); own = false)

	return pyname
end 

"""
	A::PetscMat = MatPythonCreate(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, pyname::Vector{Cchar}) 
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
function MatPythonCreate(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, pyname::Vector{Cchar}) end

@for_petsc function MatPythonCreate(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, pyname::Vector{Cchar} )
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
	newmat::PetscMat = MatCreateLocalRef(petsclib::PetscLibType,A::PetscMat, isrow::IS, iscol::IS) 
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
function MatCreateLocalRef(petsclib::PetscLibType, A::PetscMat, isrow::IS, iscol::IS) end

@for_petsc function MatCreateLocalRef(petsclib::$UnionPetscLib, A::PetscMat, isrow::IS, iscol::IS )
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
	MatShellGetContext(petsclib::PetscLibType,mat::PetscMat, ctx::Cvoid) 
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
function MatShellGetContext(petsclib::PetscLibType, mat::AbstractPetscMat, ctx::Union{Cvoid,Ptr}) end

@for_petsc function MatShellGetContext(petsclib::$UnionPetscLib, mat::AbstractPetscMat, ctx::Union{Cvoid,Ptr})

    @chk ccall(
               (:MatShellGetContext, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Cvoid}),
               mat, ctx,
              )

	return nothing
end 

"""
	MatShellSetMatProductOperation(petsclib::PetscLibType,A::PetscMat, ptype::MatProductType, symbolic::external, numeric::external, destroy::external, Btype::MatType, Ctype::MatType) 
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
function MatShellSetMatProductOperation(petsclib::PetscLibType, A::PetscMat, ptype::MatProductType, symbolic::Union{Ptr,external}, numeric::Union{Ptr,external}, destroy::Union{Ptr,external}, Btype::MatType, Ctype::Union{Ptr,MatType}) end

@for_petsc function MatShellSetMatProductOperation(petsclib::$UnionPetscLib, A::PetscMat, ptype::MatProductType, symbolic::Union{Ptr,external}, numeric::external, destroy::Union{Ptr,external}, Btype::MatType, Ctype::Union{Ptr,MatType})

    @chk ccall(
               (:MatShellSetMatProductOperation, $petsc_library),
               PetscErrorCode,
               (CMat, MatProductType, external, external, external, MatType, MatType),
               A, ptype, symbolic, numeric, destroy, Btype, Ctype,
              )


	return nothing
end 

"""
	A::PetscMat = MatCreateShell(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, ctx::Ptr) 
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
function MatCreateShell(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt) end

@for_petsc function MatCreateShell(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, ctx::Ptr )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateShell, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{Cvoid}, Ptr{CMat}),
               comm, m, n, M, N, ctx, A_,
              )

	#A = PetscMat(A_[], petsclib)
    #A = MatShell{$PetscLib, OType}(A_, ctx)
    

	return A
end 

"""
	MatShellSetContext(petsclib::PetscLibType,mat::PetscMat, ctx::Cvoid) 
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
function MatShellSetContext(petsclib::PetscLibType, mat::PetscMat, ctx::Cvoid) end

@for_petsc function MatShellSetContext(petsclib::$UnionPetscLib, mat::PetscMat, ctx::Cvoid )

    @chk ccall(
               (:MatShellSetContext, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Cvoid}),
               mat, ctx,
              )


	return nothing
end 

"""
	MatShellSetContextDestroy(petsclib::PetscLibType,mat::PetscMat, f::PetscCtxDestroyFn) 
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
function MatShellSetContextDestroy(petsclib::PetscLibType, mat::PetscMat, f::PetscCtxDestroyFn) end

@for_petsc function MatShellSetContextDestroy(petsclib::$UnionPetscLib, mat::PetscMat, f::PetscCtxDestroyFn )

    @chk ccall(
               (:MatShellSetContextDestroy, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscCtxDestroyFn}),
               mat, f,
              )


	return nothing
end 

"""
	MatShellSetVecType(petsclib::PetscLibType,mat::PetscMat, vtype::VecType) 
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
function MatShellSetVecType(petsclib::PetscLibType, mat::PetscMat, vtype::VecType) end

@for_petsc function MatShellSetVecType(petsclib::$UnionPetscLib, mat::PetscMat, vtype::VecType )

    @chk ccall(
               (:MatShellSetVecType, $petsc_library),
               PetscErrorCode,
               (CMat, VecType),
               mat, vtype,
              )


	return nothing
end 

"""
	MatShellSetManageScalingShifts(petsclib::PetscLibType,A::PetscMat) 
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
function MatShellSetManageScalingShifts(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatShellSetManageScalingShifts(petsclib::$UnionPetscLib, A::PetscMat )

    @chk ccall(
               (:MatShellSetManageScalingShifts, $petsc_library),
               PetscErrorCode,
               (CMat,),
               A,
              )


	return nothing
end 

"""
	flg::PetscBool = MatShellTestMult(petsclib::PetscLibType,mat::PetscMat, f::external, base::PetscVec, ctx::Cvoid) 
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
function MatShellTestMult(petsclib::PetscLibType, mat::PetscMat, f::external, base::PetscVec, ctx::Cvoid) end

@for_petsc function MatShellTestMult(petsclib::$UnionPetscLib, mat::PetscMat, f::external, base::PetscVec, ctx::Cvoid )
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
	flg::PetscBool = MatShellTestMultTranspose(petsclib::PetscLibType,mat::PetscMat, f::external, base::PetscVec, ctx::Cvoid) 
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
function MatShellTestMultTranspose(petsclib::PetscLibType, mat::PetscMat, f::external, base::PetscVec, ctx::Cvoid) end

@for_petsc function MatShellTestMultTranspose(petsclib::$UnionPetscLib, mat::PetscMat, f::external, base::PetscVec, ctx::Cvoid )
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
	flg::PetscBool = MatIsShell(petsclib::PetscLibType,mat::PetscMat) 
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
function MatIsShell(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatIsShell(petsclib::$UnionPetscLib, mat::PetscMat )
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
	MatMAIJGetAIJ(petsclib::PetscLibType,A::PetscMat, B::PetscMat) 
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
function MatMAIJGetAIJ(petsclib::PetscLibType, A::PetscMat, B::PetscMat) end

@for_petsc function MatMAIJGetAIJ(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat )
	B_ = Ref(B.ptr)

    @chk ccall(
               (:MatMAIJGetAIJ, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, B_,
              )

	B.ptr = C_NULL

	return nothing
end 

"""
	MatMAIJRedimension(petsclib::PetscLibType,A::PetscMat, dof::PetscInt, B::PetscMat) 
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
function MatMAIJRedimension(petsclib::PetscLibType, A::PetscMat, dof::PetscInt, B::PetscMat) end

@for_petsc function MatMAIJRedimension(petsclib::$UnionPetscLib, A::PetscMat, dof::$PetscInt, B::PetscMat )
	B_ = Ref(B.ptr)

    @chk ccall(
               (:MatMAIJRedimension, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CMat}),
               A, dof, B_,
              )

	B.ptr = C_NULL

	return nothing
end 

"""
	maij::PetscMat = MatCreateMAIJ(petsclib::PetscLibType,A::PetscMat, dof::PetscInt) 
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
function MatCreateMAIJ(petsclib::PetscLibType, A::PetscMat, dof::PetscInt) end

@for_petsc function MatCreateMAIJ(petsclib::$UnionPetscLib, A::PetscMat, dof::$PetscInt )
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
	flg::PetscBool = MatISGetAllowRepeated(petsclib::PetscLibType,A::PetscMat) 
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
function MatISGetAllowRepeated(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatISGetAllowRepeated(petsclib::$UnionPetscLib, A::PetscMat )
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
	MatISSetAllowRepeated(petsclib::PetscLibType,A::PetscMat, flg::PetscBool) 
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
function MatISSetAllowRepeated(petsclib::PetscLibType, A::PetscMat, flg::PetscBool) end

@for_petsc function MatISSetAllowRepeated(petsclib::$UnionPetscLib, A::PetscMat, flg::PetscBool )

    @chk ccall(
               (:MatISSetAllowRepeated, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               A, flg,
              )


	return nothing
end 

"""
	MatISStoreL2L(petsclib::PetscLibType,A::PetscMat, store::PetscBool) 
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
function MatISStoreL2L(petsclib::PetscLibType, A::PetscMat, store::PetscBool) end

@for_petsc function MatISStoreL2L(petsclib::$UnionPetscLib, A::PetscMat, store::PetscBool )

    @chk ccall(
               (:MatISStoreL2L, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               A, store,
              )


	return nothing
end 

"""
	MatISFixLocalEmpty(petsclib::PetscLibType,A::PetscMat, fix::PetscBool) 
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
function MatISFixLocalEmpty(petsclib::PetscLibType, A::PetscMat, fix::PetscBool) end

@for_petsc function MatISFixLocalEmpty(petsclib::$UnionPetscLib, A::PetscMat, fix::PetscBool )

    @chk ccall(
               (:MatISFixLocalEmpty, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               A, fix,
              )


	return nothing
end 

"""
	MatISSetPreallocation(petsclib::PetscLibType,B::PetscMat, d_nz::PetscInt, d_nnz::Vector{PetscInt}, o_nz::PetscInt, o_nnz::Vector{PetscInt}) 
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
function MatISSetPreallocation(petsclib::PetscLibType, B::PetscMat, d_nz::PetscInt, d_nnz::Union{Ptr,Vector{PetscInt}}, o_nz::PetscInt, o_nnz::Union{Ptr,Vector{PetscInt}}) end

@for_petsc function MatISSetPreallocation(petsclib::$UnionPetscLib, B::PetscMat, d_nz::$PetscInt, d_nnz::Union{Ptr,Vector{$PetscInt}}, o_nz::$PetscInt, o_nnz::Union{Ptr,Vector{$PetscInt}})

    @chk ccall(
               (:MatISSetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}),
               B, d_nz, d_nnz, o_nz, o_nnz,
              )


	return nothing
end 

"""
	MatISGetLocalMat(petsclib::PetscLibType,mat::PetscMat, loc::PetscMat) 
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
function MatISGetLocalMat(petsclib::PetscLibType, mat::PetscMat, loc::PetscMat) end

@for_petsc function MatISGetLocalMat(petsclib::$UnionPetscLib, mat::PetscMat, loc::PetscMat )
	loc_ = Ref(loc.ptr)

    @chk ccall(
               (:MatISGetLocalMat, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               mat, loc_,
              )

	loc.ptr = C_NULL

	return nothing
end 

"""
	MatISRestoreLocalMat(petsclib::PetscLibType,mat::PetscMat, loc::PetscMat) 
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
function MatISRestoreLocalMat(petsclib::PetscLibType, mat::PetscMat, loc::PetscMat) end

@for_petsc function MatISRestoreLocalMat(petsclib::$UnionPetscLib, mat::PetscMat, loc::PetscMat )
	loc_ = Ref(loc.ptr)

    @chk ccall(
               (:MatISRestoreLocalMat, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               mat, loc_,
              )

	loc.ptr = C_NULL

	return nothing
end 

"""
	MatISSetLocalMatType(petsclib::PetscLibType,mat::PetscMat, mtype::MatType) 
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
function MatISSetLocalMatType(petsclib::PetscLibType, mat::PetscMat, mtype::MatType) end

@for_petsc function MatISSetLocalMatType(petsclib::$UnionPetscLib, mat::PetscMat, mtype::MatType )

    @chk ccall(
               (:MatISSetLocalMatType, $petsc_library),
               PetscErrorCode,
               (CMat, MatType),
               mat, mtype,
              )


	return nothing
end 

"""
	MatISSetLocalMat(petsclib::PetscLibType,mat::PetscMat, loc::PetscMat) 
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
function MatISSetLocalMat(petsclib::PetscLibType, mat::PetscMat, loc::PetscMat) end

@for_petsc function MatISSetLocalMat(petsclib::$UnionPetscLib, mat::PetscMat, loc::PetscMat )

    @chk ccall(
               (:MatISSetLocalMat, $petsc_library),
               PetscErrorCode,
               (CMat, CMat),
               mat, loc,
              )


	return nothing
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
	MatISGetLocalToGlobalMapping(petsclib::PetscLibType,A::PetscMat, rmapping::ISLocalToGlobalMapping, cmapping::ISLocalToGlobalMapping) 
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
function MatISGetLocalToGlobalMapping(petsclib::PetscLibType, A::PetscMat, rmapping::ISLocalToGlobalMapping, cmapping::ISLocalToGlobalMapping) end

@for_petsc function MatISGetLocalToGlobalMapping(petsclib::$UnionPetscLib, A::PetscMat, rmapping::ISLocalToGlobalMapping, cmapping::ISLocalToGlobalMapping )

    @chk ccall(
               (:MatISGetLocalToGlobalMapping, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{ISLocalToGlobalMapping}, Ptr{ISLocalToGlobalMapping}),
               A, rmapping, cmapping,
              )


	return nothing
end 

"""
	MatHermitianTransposeGetMat(petsclib::PetscLibType,A::PetscMat, M::PetscMat) 
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
function MatHermitianTransposeGetMat(petsclib::PetscLibType, A::PetscMat, M::PetscMat) end

@for_petsc function MatHermitianTransposeGetMat(petsclib::$UnionPetscLib, A::PetscMat, M::PetscMat )
	M_ = Ref(M.ptr)

    @chk ccall(
               (:MatHermitianTransposeGetMat, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, M_,
              )

	M.ptr = C_NULL

	return nothing
end 

"""
	N::PetscMat = MatCreateHermitianTranspose(petsclib::PetscLibType,A::PetscMat) 
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
function MatCreateHermitianTranspose(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatCreateHermitianTranspose(petsclib::$UnionPetscLib, A::PetscMat )
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
	MatTransposeGetMat(petsclib::PetscLibType,A::PetscMat, M::PetscMat) 
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
function MatTransposeGetMat(petsclib::PetscLibType, A::PetscMat, M::PetscMat) end

@for_petsc function MatTransposeGetMat(petsclib::$UnionPetscLib, A::PetscMat, M::PetscMat )
	M_ = Ref(M.ptr)

    @chk ccall(
               (:MatTransposeGetMat, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, M_,
              )

	M.ptr = C_NULL

	return nothing
end 

"""
	N::PetscMat = MatCreateTranspose(petsclib::PetscLibType,A::PetscMat) 
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
function MatCreateTranspose(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatCreateTranspose(petsclib::$UnionPetscLib, A::PetscMat )
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
	MatDenseGetLocalMatrix(petsclib::PetscLibType,A::PetscMat, B::PetscMat) 
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
function MatDenseGetLocalMatrix(petsclib::PetscLibType, A::PetscMat, B::PetscMat) end

@for_petsc function MatDenseGetLocalMatrix(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat )
	B_ = Ref(B.ptr)

    @chk ccall(
               (:MatDenseGetLocalMatrix, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, B_,
              )

	B.ptr = C_NULL

	return nothing
end 

"""
	MatMPIDenseSetPreallocation(petsclib::PetscLibType,B::PetscMat, data::PetscScalar) 
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
function MatMPIDenseSetPreallocation(petsclib::PetscLibType, B::PetscMat) end

@for_petsc function MatMPIDenseSetPreallocation(petsclib::$UnionPetscLib, B::PetscMat, data::Union{Ptr,$PetscScalar} )

    @chk ccall(
               (:MatMPIDenseSetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscScalar),
               B, data,
              )

	return nothing
end 

"""
	array::PetscScalar = MatDensePlaceArray(petsclib::PetscLibType,mat::PetscMat) 
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
function MatDensePlaceArray(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatDensePlaceArray(petsclib::$UnionPetscLib, mat::PetscMat )
	array_ = Ref{$PetscScalar}()

    @chk ccall(
               (:MatDensePlaceArray, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscScalar}),
               mat, array_,
              )

	array = array_[]

	return array
end 

"""
	MatDenseResetArray(petsclib::PetscLibType,mat::PetscMat) 
Resets the matrix array to that it previously had before the call to `MatDensePlaceArray()`

Not Collective

Input Parameter:
- `mat` - the matrix

Level: developer

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MatDenseGetArray()`, `MatDensePlaceArray()`, `VecPlaceArray()`, `VecGetArray()`, `VecRestoreArray()`, `VecReplaceArray()`, `VecResetArray()`

# External Links
$(_doc_external("Mat/MatDenseResetArray"))
"""
function MatDenseResetArray(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatDenseResetArray(petsclib::$UnionPetscLib, mat::PetscMat )

    @chk ccall(
               (:MatDenseResetArray, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	array::PetscScalar = MatDenseReplaceArray(petsclib::PetscLibType,mat::PetscMat) 
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
function MatDenseReplaceArray(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatDenseReplaceArray(petsclib::$UnionPetscLib, mat::PetscMat )
	array_ = Ref{$PetscScalar}()

    @chk ccall(
               (:MatDenseReplaceArray, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscScalar}),
               mat, array_,
              )

	array = array_[]

	return array
end 

"""
	A::PetscMat = MatCreateDense(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, data::Union{Ptr,Vector{PetscScalar}}) 
Creates a matrix in `MATDENSE` format.

Collective

Input Parameters:
- `comm` - MPI communicator
- `m`    - number of local rows (or `PETSC_DECIDE` to have calculated if `M` is given)
- `n`    - number of local columns (or `PETSC_DECIDE` to have calculated if `N` is given)
- `M`    - number of global rows (or `PETSC_DECIDE` to have calculated if `m` is given)
- `N`    - number of global columns (or `PETSC_DECIDE` to have calculated if `n` is given)
- `data` - optional location of matrix data.  Set data to `NULL` (`PETSC_NULL_SCALAR_ARRAY` for Fortran users) for PETSc to control all matrix memory allocation.

Output Parameter:
- `A` - the matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATDENSE`, `MatCreate()`, `MatCreateSeqDense()`, `MatSetValues()`

# External Links
$(_doc_external("Mat/MatCreateDense"))
"""
function MatCreateDense(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, data::Union{Ptr,Vector{PetscScalar}}) end

@for_petsc function MatCreateDense(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, data::Union{Ptr,Vector{$PetscScalar}} )
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
	lda::PetscInt = MatDenseGetLDA(petsclib::PetscLibType,A::PetscMat) 
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
function MatDenseGetLDA(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatDenseGetLDA(petsclib::$UnionPetscLib, A::PetscMat )
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
	MatDenseSetLDA(petsclib::PetscLibType,A::PetscMat, lda::PetscInt) 
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
function MatDenseSetLDA(petsclib::PetscLibType, A::PetscMat, lda::PetscInt) end

@for_petsc function MatDenseSetLDA(petsclib::$UnionPetscLib, A::PetscMat, lda::$PetscInt )

    @chk ccall(
               (:MatDenseSetLDA, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt),
               A, lda,
              )


	return nothing
end 

"""
	array::Vector{PetscScalar} = MatDenseGetArray(petsclib::PetscLibType,A::PetscMat) 
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
function MatDenseGetArray(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatDenseGetArray(petsclib::$UnionPetscLib, A::PetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatDenseGetArray, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

	array = unsafe_wrap(Array, array_[], VecGetLocalSize(petsclib, x); own = false)

	return array
end 

"""
	array::Vector{PetscScalar} = MatDenseRestoreArray(petsclib::PetscLibType,A::PetscMat) 
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
function MatDenseRestoreArray(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatDenseRestoreArray(petsclib::$UnionPetscLib, A::PetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatDenseRestoreArray, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

	array = unsafe_wrap(Array, array_[], VecGetLocalSize(petsclib, x); own = false)

	return array
end 

"""
	array::Vector{PetscScalar} = MatDenseGetArrayRead(petsclib::PetscLibType,A::PetscMat) 
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
function MatDenseGetArrayRead(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatDenseGetArrayRead(petsclib::$UnionPetscLib, A::PetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatDenseGetArrayRead, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

	array = unsafe_wrap(Array, array_[], VecGetLocalSize(petsclib, x); own = false)

	return array
end 

"""
	array::Vector{PetscScalar} = MatDenseRestoreArrayRead(petsclib::PetscLibType,A::PetscMat) 
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
function MatDenseRestoreArrayRead(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatDenseRestoreArrayRead(petsclib::$UnionPetscLib, A::PetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatDenseRestoreArrayRead, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

	array = unsafe_wrap(Array, array_[], VecGetLocalSize(petsclib, x); own = false)

	return array
end 

"""
	array::Vector{PetscScalar} = MatDenseGetArrayWrite(petsclib::PetscLibType,A::PetscMat) 
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
function MatDenseGetArrayWrite(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatDenseGetArrayWrite(petsclib::$UnionPetscLib, A::PetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatDenseGetArrayWrite, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

	array = unsafe_wrap(Array, array_[], VecGetLocalSize(petsclib, x); own = false)

	return array
end 

"""
	array::Vector{PetscScalar} = MatDenseRestoreArrayWrite(petsclib::PetscLibType,A::PetscMat) 
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
function MatDenseRestoreArrayWrite(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatDenseRestoreArrayWrite(petsclib::$UnionPetscLib, A::PetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatDenseRestoreArrayWrite, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

	array = unsafe_wrap(Array, array_[], VecGetLocalSize(petsclib, x); own = false)

	return array
end 

"""
	array::Vector{PetscScalar},mtype::PetscMemType = MatDenseGetArrayAndMemType(petsclib::PetscLibType,A::PetscMat) 
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
function MatDenseGetArrayAndMemType(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatDenseGetArrayAndMemType(petsclib::$UnionPetscLib, A::PetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()
	mtype_ = Ref{PetscMemType}()

    @chk ccall(
               (:MatDenseGetArrayAndMemType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}, Ptr{PetscMemType}),
               A, array_, mtype_,
              )

	array = unsafe_wrap(Array, array_[], VecGetLocalSize(petsclib, x); own = false)
	mtype = unsafe_string(mtype_[])

	return array,mtype
end 

"""
	array::Vector{PetscScalar} = MatDenseRestoreArrayAndMemType(petsclib::PetscLibType,A::PetscMat) 
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
function MatDenseRestoreArrayAndMemType(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatDenseRestoreArrayAndMemType(petsclib::$UnionPetscLib, A::PetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatDenseRestoreArrayAndMemType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

	array = unsafe_wrap(Array, array_[], VecGetLocalSize(petsclib, x); own = false)

	return array
end 

"""
	array::Vector{PetscScalar},mtype::PetscMemType = MatDenseGetArrayReadAndMemType(petsclib::PetscLibType,A::PetscMat) 
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
function MatDenseGetArrayReadAndMemType(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatDenseGetArrayReadAndMemType(petsclib::$UnionPetscLib, A::PetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()
	mtype_ = Ref{PetscMemType}()

    @chk ccall(
               (:MatDenseGetArrayReadAndMemType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}, Ptr{PetscMemType}),
               A, array_, mtype_,
              )

	array = unsafe_wrap(Array, array_[], VecGetLocalSize(petsclib, x); own = false)
	mtype = unsafe_string(mtype_[])

	return array,mtype
end 

"""
	array::Vector{PetscScalar} = MatDenseRestoreArrayReadAndMemType(petsclib::PetscLibType,A::PetscMat) 
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
function MatDenseRestoreArrayReadAndMemType(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatDenseRestoreArrayReadAndMemType(petsclib::$UnionPetscLib, A::PetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatDenseRestoreArrayReadAndMemType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

	array = unsafe_wrap(Array, array_[], VecGetLocalSize(petsclib, x); own = false)

	return array
end 

"""
	array::Vector{PetscScalar},mtype::PetscMemType = MatDenseGetArrayWriteAndMemType(petsclib::PetscLibType,A::PetscMat) 
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
function MatDenseGetArrayWriteAndMemType(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatDenseGetArrayWriteAndMemType(petsclib::$UnionPetscLib, A::PetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()
	mtype_ = Ref{PetscMemType}()

    @chk ccall(
               (:MatDenseGetArrayWriteAndMemType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}, Ptr{PetscMemType}),
               A, array_, mtype_,
              )

	array = unsafe_wrap(Array, array_[], VecGetLocalSize(petsclib, x); own = false)
	mtype = unsafe_string(mtype_[])

	return array,mtype
end 

"""
	array::Vector{PetscScalar} = MatDenseRestoreArrayWriteAndMemType(petsclib::PetscLibType,A::PetscMat) 
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
function MatDenseRestoreArrayWriteAndMemType(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatDenseRestoreArrayWriteAndMemType(petsclib::$UnionPetscLib, A::PetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatDenseRestoreArrayWriteAndMemType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

	array = unsafe_wrap(Array, array_[], VecGetLocalSize(petsclib, x); own = false)

	return array
end 

"""
	A::PetscMat = MatCreateSeqDense(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, data::Vector{PetscScalar}) 
Creates a `MATSEQDENSE` that
is stored in column major order (the usual Fortran format).

Collective

Input Parameters:
- `comm` - MPI communicator, set to `PETSC_COMM_SELF`
- `m`    - number of rows
- `n`    - number of columns
- `data` - optional location of matrix data in column major order.  Use `NULL` for PETSc to control all matrix memory allocation.

Output Parameter:
- `A` - the matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATSEQDENSE`, `MatCreate()`, `MatCreateDense()`, `MatSetValues()`

# External Links
$(_doc_external("Mat/MatCreateSeqDense"))
"""
function MatCreateSeqDense(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, data::Union{Ptr,Vector{PetscScalar}}) end

@for_petsc function MatCreateSeqDense(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, data::Union{Ptr,Vector{$PetscScalar}})
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
	MatSeqDenseSetPreallocation(petsclib::PetscLibType,B::PetscMat, data::Vector{PetscScalar}) 
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
function MatSeqDenseSetPreallocation(petsclib::PetscLibType, B::PetscMat, data::Vector{PetscScalar}) end

@for_petsc function MatSeqDenseSetPreallocation(petsclib::$UnionPetscLib, B::PetscMat, data::Vector{$PetscScalar} )

    @chk ccall(
               (:MatSeqDenseSetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscScalar}),
               B, data,
              )


	return nothing
end 

"""
	vals::Vector{PetscScalar} = MatDenseGetColumn(petsclib::PetscLibType,A::PetscMat, col::PetscInt) 
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
function MatDenseGetColumn(petsclib::PetscLibType, A::PetscMat, col::PetscInt) end

@for_petsc function MatDenseGetColumn(petsclib::$UnionPetscLib, A::PetscMat, col::$PetscInt )
	vals_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatDenseGetColumn, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{Ptr{$PetscScalar}}),
               A, col, vals_,
              )

	vals = unsafe_wrap(Array, vals_[], VecGetLocalSize(petsclib, x); own = false)

	return vals
end 

"""
	vals::Vector{PetscScalar} = MatDenseRestoreColumn(petsclib::PetscLibType,A::PetscMat) 
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
function MatDenseRestoreColumn(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatDenseRestoreColumn(petsclib::$UnionPetscLib, A::PetscMat )
	vals_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatDenseRestoreColumn, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, vals_,
              )

	vals = unsafe_wrap(Array, vals_[], VecGetLocalSize(petsclib, x); own = false)

	return vals
end 

"""
	MatDenseGetColumnVec(petsclib::PetscLibType,A::PetscMat, col::PetscInt, v::PetscVec) 
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
function MatDenseGetColumnVec(petsclib::PetscLibType, A::PetscMat, col::PetscInt, v::PetscVec) end

@for_petsc function MatDenseGetColumnVec(petsclib::$UnionPetscLib, A::PetscMat, col::$PetscInt, v::PetscVec )
	v_ = Ref(v.ptr)

    @chk ccall(
               (:MatDenseGetColumnVec, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CVec}),
               A, col, v_,
              )

	v.ptr = C_NULL

	return nothing
end 

"""
	MatDenseRestoreColumnVec(petsclib::PetscLibType,A::PetscMat, col::PetscInt, v::PetscVec) 
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
function MatDenseRestoreColumnVec(petsclib::PetscLibType, A::PetscMat, col::PetscInt, v::PetscVec) end

@for_petsc function MatDenseRestoreColumnVec(petsclib::$UnionPetscLib, A::PetscMat, col::$PetscInt, v::PetscVec )
	v_ = Ref(v.ptr)

    @chk ccall(
               (:MatDenseRestoreColumnVec, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CVec}),
               A, col, v_,
              )

	v.ptr = C_NULL

	return nothing
end 

"""
	MatDenseGetColumnVecRead(petsclib::PetscLibType,A::PetscMat, col::PetscInt, v::PetscVec) 
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
function MatDenseGetColumnVecRead(petsclib::PetscLibType, A::PetscMat, col::PetscInt, v::PetscVec) end

@for_petsc function MatDenseGetColumnVecRead(petsclib::$UnionPetscLib, A::PetscMat, col::$PetscInt, v::PetscVec )
	v_ = Ref(v.ptr)

    @chk ccall(
               (:MatDenseGetColumnVecRead, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CVec}),
               A, col, v_,
              )

	v.ptr = C_NULL

	return nothing
end 

"""
	MatDenseRestoreColumnVecRead(petsclib::PetscLibType,A::PetscMat, col::PetscInt, v::PetscVec) 
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
function MatDenseRestoreColumnVecRead(petsclib::PetscLibType, A::PetscMat, col::PetscInt, v::PetscVec) end

@for_petsc function MatDenseRestoreColumnVecRead(petsclib::$UnionPetscLib, A::PetscMat, col::$PetscInt, v::PetscVec )
	v_ = Ref(v.ptr)

    @chk ccall(
               (:MatDenseRestoreColumnVecRead, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CVec}),
               A, col, v_,
              )

	v.ptr = C_NULL

	return nothing
end 

"""
	MatDenseGetColumnVecWrite(petsclib::PetscLibType,A::PetscMat, col::PetscInt, v::PetscVec) 
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
function MatDenseGetColumnVecWrite(petsclib::PetscLibType, A::PetscMat, col::PetscInt, v::PetscVec) end

@for_petsc function MatDenseGetColumnVecWrite(petsclib::$UnionPetscLib, A::PetscMat, col::$PetscInt, v::PetscVec )
	v_ = Ref(v.ptr)

    @chk ccall(
               (:MatDenseGetColumnVecWrite, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CVec}),
               A, col, v_,
              )

	v.ptr = C_NULL

	return nothing
end 

"""
	MatDenseRestoreColumnVecWrite(petsclib::PetscLibType,A::PetscMat, col::PetscInt, v::PetscVec) 
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
function MatDenseRestoreColumnVecWrite(petsclib::PetscLibType, A::PetscMat, col::PetscInt, v::PetscVec) end

@for_petsc function MatDenseRestoreColumnVecWrite(petsclib::$UnionPetscLib, A::PetscMat, col::$PetscInt, v::PetscVec )
	v_ = Ref(v.ptr)

    @chk ccall(
               (:MatDenseRestoreColumnVecWrite, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CVec}),
               A, col, v_,
              )

	v.ptr = C_NULL

	return nothing
end 

"""
	MatDenseGetSubMatrix(petsclib::PetscLibType,A::PetscMat, rbegin::PetscInt, rend::PetscInt, cbegin::PetscInt, cend::PetscInt, v::PetscMat) 
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
function MatDenseGetSubMatrix(petsclib::PetscLibType, A::PetscMat, rbegin::PetscInt, rend::PetscInt, cbegin::PetscInt, cend::PetscInt, v::PetscMat) end

@for_petsc function MatDenseGetSubMatrix(petsclib::$UnionPetscLib, A::PetscMat, rbegin::$PetscInt, rend::$PetscInt, cbegin::$PetscInt, cend::$PetscInt, v::PetscMat )
	v_ = Ref(v.ptr)

    @chk ccall(
               (:MatDenseGetSubMatrix, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{CMat}),
               A, rbegin, rend, cbegin, cend, v_,
              )

	v.ptr = C_NULL

	return nothing
end 

"""
	MatDenseRestoreSubMatrix(petsclib::PetscLibType,A::PetscMat, v::PetscMat) 
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
function MatDenseRestoreSubMatrix(petsclib::PetscLibType, A::PetscMat, v::PetscMat) end

@for_petsc function MatDenseRestoreSubMatrix(petsclib::$UnionPetscLib, A::PetscMat, v::PetscMat )
	v_ = Ref(v.ptr)

    @chk ccall(
               (:MatDenseRestoreSubMatrix, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, v_,
              )

	v.ptr = C_NULL

	return nothing
end 

"""
	MatSeqDenseInvert(petsclib::PetscLibType,A::PetscMat) 

# External Links
$(_doc_external("Mat/MatSeqDenseInvert"))
"""
function MatSeqDenseInvert(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatSeqDenseInvert(petsclib::$UnionPetscLib, A::PetscMat )

    @chk ccall(
               (:MatSeqDenseInvert, $petsc_library),
               PetscErrorCode,
               (CMat,),
               A,
              )


	return nothing
end 

"""
	mats::PetscMat,mat::PetscMat = MatCreateComposite(petsclib::PetscLibType,comm::MPI_Comm, nmat::PetscInt) 
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
function MatCreateComposite(petsclib::PetscLibType, comm::MPI_Comm, nmat::PetscInt) end

@for_petsc function MatCreateComposite(petsclib::$UnionPetscLib, comm::MPI_Comm, nmat::$PetscInt )
	mats_ = Ref{CMat}()
	mat_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateComposite, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, Ptr{CMat}, Ptr{CMat}),
               comm, nmat, mats_, mat_,
              )

	mats = PetscMat(mats_[], petsclib)
	mat = PetscMat(mat_[], petsclib)

	return mats,mat
end 

"""
	MatCompositeAddMat(petsclib::PetscLibType,mat::PetscMat, smat::PetscMat) 
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
function MatCompositeAddMat(petsclib::PetscLibType, mat::PetscMat, smat::PetscMat) end

@for_petsc function MatCompositeAddMat(petsclib::$UnionPetscLib, mat::PetscMat, smat::PetscMat )

    @chk ccall(
               (:MatCompositeAddMat, $petsc_library),
               PetscErrorCode,
               (CMat, CMat),
               mat, smat,
              )


	return nothing
end 

"""
	MatCompositeSetType(petsclib::PetscLibType,mat::PetscMat, type::MatCompositeType) 
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
function MatCompositeSetType(petsclib::PetscLibType, mat::PetscMat, type::MatCompositeType) end

@for_petsc function MatCompositeSetType(petsclib::$UnionPetscLib, mat::PetscMat, type::MatCompositeType )

    @chk ccall(
               (:MatCompositeSetType, $petsc_library),
               PetscErrorCode,
               (CMat, MatCompositeType),
               mat, type,
              )


	return nothing
end 

"""
	type::MatCompositeType = MatCompositeGetType(petsclib::PetscLibType,mat::PetscMat) 
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
function MatCompositeGetType(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatCompositeGetType(petsclib::$UnionPetscLib, mat::PetscMat )
	type_ = Ref{MatCompositeType}()

    @chk ccall(
               (:MatCompositeGetType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatCompositeType}),
               mat, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	MatCompositeSetMatStructure(petsclib::PetscLibType,mat::PetscMat, str::MatStructure) 
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
function MatCompositeSetMatStructure(petsclib::PetscLibType, mat::PetscMat, str::MatStructure) end

@for_petsc function MatCompositeSetMatStructure(petsclib::$UnionPetscLib, mat::PetscMat, str::MatStructure )

    @chk ccall(
               (:MatCompositeSetMatStructure, $petsc_library),
               PetscErrorCode,
               (CMat, MatStructure),
               mat, str,
              )


	return nothing
end 

"""
	MatCompositeGetMatStructure(petsclib::PetscLibType,mat::PetscMat, str::MatStructure) 
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
function MatCompositeGetMatStructure(petsclib::PetscLibType, mat::PetscMat, str::MatStructure) end

@for_petsc function MatCompositeGetMatStructure(petsclib::$UnionPetscLib, mat::PetscMat, str::MatStructure )

    @chk ccall(
               (:MatCompositeGetMatStructure, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatStructure}),
               mat, str,
              )


	return nothing
end 

"""
	MatCompositeSetMergeType(petsclib::PetscLibType,mat::PetscMat, type::MatCompositeMergeType) 
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
function MatCompositeSetMergeType(petsclib::PetscLibType, mat::PetscMat, type::MatCompositeMergeType) end

@for_petsc function MatCompositeSetMergeType(petsclib::$UnionPetscLib, mat::PetscMat, type::MatCompositeMergeType )

    @chk ccall(
               (:MatCompositeSetMergeType, $petsc_library),
               PetscErrorCode,
               (CMat, MatCompositeMergeType),
               mat, type,
              )


	return nothing
end 

"""
	MatCompositeMerge(petsclib::PetscLibType,mat::PetscMat) 
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
function MatCompositeMerge(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatCompositeMerge(petsclib::$UnionPetscLib, mat::PetscMat )

    @chk ccall(
               (:MatCompositeMerge, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	nmat::PetscInt = MatCompositeGetNumberMat(petsclib::PetscLibType,mat::PetscMat) 
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
function MatCompositeGetNumberMat(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatCompositeGetNumberMat(petsclib::$UnionPetscLib, mat::PetscMat )
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
	MatCompositeGetMat(petsclib::PetscLibType,mat::PetscMat, i::PetscInt, Ai::PetscMat) 
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
function MatCompositeGetMat(petsclib::PetscLibType, mat::PetscMat, i::PetscInt, Ai::PetscMat) end

@for_petsc function MatCompositeGetMat(petsclib::$UnionPetscLib, mat::PetscMat, i::$PetscInt, Ai::PetscMat )
	Ai_ = Ref(Ai.ptr)

    @chk ccall(
               (:MatCompositeGetMat, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CMat}),
               mat, i, Ai_,
              )

	Ai.ptr = C_NULL

	return nothing
end 

"""
	scalings::PetscScalar = MatCompositeSetScalings(petsclib::PetscLibType,mat::PetscMat) 
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
function MatCompositeSetScalings(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatCompositeSetScalings(petsclib::$UnionPetscLib, mat::PetscMat )
	scalings_ = Ref{$PetscScalar}()

    @chk ccall(
               (:MatCompositeSetScalings, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscScalar}),
               mat, scalings_,
              )

	scalings = scalings_[]

	return scalings
end 

"""
	MatLRCGetMats(petsclib::PetscLibType,N::PetscMat, A::PetscMat, U::PetscMat, c::PetscVec, V::PetscMat) 
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
function MatLRCGetMats(petsclib::PetscLibType, N::PetscMat, A::PetscMat, U::PetscMat, c::PetscVec, V::PetscMat) end

@for_petsc function MatLRCGetMats(petsclib::$UnionPetscLib, N::PetscMat, A::PetscMat, U::PetscMat, c::PetscVec, V::PetscMat )
	A_ = Ref(A.ptr)
	U_ = Ref(U.ptr)
	c_ = Ref(c.ptr)
	V_ = Ref(V.ptr)

    @chk ccall(
               (:MatLRCGetMats, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}, Ptr{CMat}, Ptr{CVec}, Ptr{CMat}),
               N, A_, U_, c_, V_,
              )

	A.ptr = C_NULL
	U.ptr = C_NULL
	c.ptr = C_NULL
	V.ptr = C_NULL

	return nothing
end 

"""
	MatLRCSetMats(petsclib::PetscLibType,N::PetscMat, A::PetscMat, U::PetscMat, c::PetscVec, V::PetscMat) 
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
function MatLRCSetMats(petsclib::PetscLibType, N::PetscMat, A::PetscMat, U::PetscMat, c::PetscVec, V::PetscMat) end

@for_petsc function MatLRCSetMats(petsclib::$UnionPetscLib, N::PetscMat, A::PetscMat, U::PetscMat, c::PetscVec, V::PetscMat )

    @chk ccall(
               (:MatLRCSetMats, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat, CVec, CMat),
               N, A, U, c, V,
              )


	return nothing
end 

"""
	N::PetscMat = MatCreateLRC(petsclib::PetscLibType,A::PetscMat, U::PetscMat, c::PetscVec, V::PetscMat) 
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
function MatCreateLRC(petsclib::PetscLibType, A::PetscMat, U::PetscMat, c::Union{Ptr,PetscVec}, V::PetscMat) end

@for_petsc function MatCreateLRC(petsclib::$UnionPetscLib, A::PetscMat, U::PetscMat, c::Union{Ptr,PetscVec}, V::PetscMat )
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
	J::PetscMat = MatCreateMFFD(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt) 
Creates a matrix
approximately multiply a vector by the matrix (Jacobian) . See also `MatCreateSNESMF()`

Collective

Input Parameters:
- `comm` - MPI communicator
- `m`    - number of local rows (or `PETSC_DECIDE` to have calculated if `M` is given). This value should be the same as the local size used in creating the y vector for the matrix-vector product y = Ax.
- `n`    - This value should be the same as the local size used in creating the x vector for the matrix-vector product y = Ax. (or `PETSC_DECIDE` to have calculated if `N` is given) For square matrices `n` is almost always `m`.
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
- `-mat_mffd_complex`          - use the Lyness trick with complex numbers to compute the matrix-vector product instead of differencing (requires real valued functions but that PETSc be configured for complex numbers)
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
	MatMPISBAIJSetPreallocation(petsclib::PetscLibType,B::PetscMat, bs::PetscInt, d_nz::PetscInt, d_nnz::Vector{PetscInt}, o_nz::PetscInt, o_nnz::Vector{PetscInt}) 
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
function MatMPISBAIJSetPreallocation(petsclib::PetscLibType, B::PetscMat, bs::PetscInt, d_nz::PetscInt, d_nnz::Vector{PetscInt}, o_nz::PetscInt, o_nnz::Vector{PetscInt}) end

@for_petsc function MatMPISBAIJSetPreallocation(petsclib::$UnionPetscLib, B::PetscMat, bs::$PetscInt, d_nz::$PetscInt, d_nnz::Vector{$PetscInt}, o_nz::$PetscInt, o_nnz::Vector{$PetscInt} )

    @chk ccall(
               (:MatMPISBAIJSetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}),
               B, bs, d_nz, d_nnz, o_nz, o_nnz,
              )


	return nothing
end 

"""
	A::PetscMat = MatCreateSBAIJ(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_nz::PetscInt, d_nnz::Vector{PetscInt}, o_nz::PetscInt, o_nnz::Vector{PetscInt}) 
Creates a sparse parallel matrix in symmetric block AIJ format, `MATSBAIJ`,
(block compressed row).  For good matrix assembly performance
the user should preallocate the matrix storage by setting the parameters
`d_nz` (or `d_nnz`) and `o_nz` (or `o_nnz`).

Collective

Input Parameters:
- `comm`  - MPI communicator
- `bs`    - size of block, the blocks are ALWAYS square. One can use `MatSetBlockSizes()` to set a different row and column blocksize but the row blocksize always defines the size of the blocks. The column blocksize sets the blocksize of the vectors obtained with `MatCreateVecs()`
- `m`     - number of local rows (or `PETSC_DECIDE` to have calculated if `M` is given) This value should be the same as the local size used in creating the y vector for the matrix-vector product y = Ax.
- `n`     - number of local columns (or `PETSC_DECIDE` to have calculated if `N` is given) This value should be the same as the local size used in creating the x vector for the matrix-vector product y = Ax.
- `M`     - number of global rows (or `PETSC_DETERMINE` to have calculated if `m` is given)
- `N`     - number of global columns (or `PETSC_DETERMINE` to have calculated if `n` is given)
- `d_nz`  - number of block nonzeros per block row in diagonal portion of local submatrix (same for all local rows)
- `d_nnz` - array containing the number of block nonzeros in the various block rows in the upper triangular portion of the in diagonal portion of the local (possibly different for each block block row) or `NULL`. If you plan to factor the matrix you must leave room for the diagonal entry and set its value even if it is zero.
- `o_nz`  - number of block nonzeros per block row in the off-diagonal portion of local submatrix (same for all local rows).
- `o_nnz` - array containing the number of nonzeros in the various block rows of the off-diagonal portion of the local submatrix (possibly different for each block row) or `NULL`.

Output Parameter:
- `A` - the matrix

Options Database Keys:
- `-mat_no_unroll`  - uses code that does not unroll the loops in the block calculations (much slower)
- `-mat_block_size` - size of the blocks to use
- `-mat_mpi`        - use the parallel matrix data structures even on one processor (defaults to using SeqBAIJ format on one processor)

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATSBAIJ`, `MatCreate()`, `MatCreateSeqSBAIJ()`, `MatSetValues()`, `MatCreateBAIJ()`,
`MatGetOwnershipRange()`,  `MatGetOwnershipRanges()`, `MatGetOwnershipRangeColumn()`, `MatGetOwnershipRangesColumn()`, `PetscLayout`

# External Links
$(_doc_external("Mat/MatCreateSBAIJ"))
"""
function MatCreateSBAIJ(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_nz::PetscInt, d_nnz::Union{Ptr,Vector{PetscInt}}, o_nz::PetscInt, o_nnz::Union{Ptr,Vector{PetscInt}}) end

@for_petsc function MatCreateSBAIJ(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, d_nz::$PetscInt, d_nnz::Union{Ptr,Vector{$PetscInt}}, o_nz::$PetscInt, o_nnz::Union{Ptr,Vector{$PetscInt}} )
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
	mat::PetscMat = MatCreateMPISBAIJWithArrays(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, a::Vector{PetscScalar}) 
creates a `MATMPISBAIJ` matrix using arrays that contain in standard CSR format for the local rows.

Collective

Input Parameters:
- `comm` - MPI communicator
- `bs`   - the block size, only a block size of 1 is supported
- `m`    - number of local rows (Cannot be `PETSC_DECIDE`)
- `n`    - This value should be the same as the local size used in creating the x vector for the matrix-vector product  y = Ax . (or `PETSC_DECIDE` to have calculated if `N` is given) For square matrices `n` is almost always `m`.
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
	MatMPISBAIJSetPreallocationCSR(petsclib::PetscLibType,B::PetscMat, bs::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, v::Vector{PetscScalar}) 
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
function MatMPISBAIJSetPreallocationCSR(petsclib::PetscLibType, B::PetscMat, bs::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, v::Vector{PetscScalar}) end

@for_petsc function MatMPISBAIJSetPreallocationCSR(petsclib::$UnionPetscLib, B::PetscMat, bs::$PetscInt, i::Vector{$PetscInt}, j::Vector{$PetscInt}, v::Vector{$PetscScalar} )

    @chk ccall(
               (:MatMPISBAIJSetPreallocationCSR, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}),
               B, bs, i, j, v,
              )


	return nothing
end 

"""
	indices::PetscInt = MatSeqSBAIJSetColumnIndices(petsclib::PetscLibType,mat::PetscMat) 
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
function MatSeqSBAIJSetColumnIndices(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatSeqSBAIJSetColumnIndices(petsclib::$UnionPetscLib, mat::PetscMat )
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
	array::Vector{PetscScalar} = MatSeqSBAIJGetArray(petsclib::PetscLibType,A::PetscMat) 
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
function MatSeqSBAIJGetArray(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatSeqSBAIJGetArray(petsclib::$UnionPetscLib, A::PetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatSeqSBAIJGetArray, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

	array = unsafe_wrap(Array, array_[], VecGetLocalSize(petsclib, x); own = false)

	return array
end 

"""
	array::Vector{PetscScalar} = MatSeqSBAIJRestoreArray(petsclib::PetscLibType,A::PetscMat) 
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
function MatSeqSBAIJRestoreArray(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatSeqSBAIJRestoreArray(petsclib::$UnionPetscLib, A::PetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatSeqSBAIJRestoreArray, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

	array = unsafe_wrap(Array, array_[], VecGetLocalSize(petsclib, x); own = false)

	return array
end 

"""
	MatSeqSBAIJSetPreallocation(petsclib::PetscLibType,B::PetscMat, bs::PetscInt, nz::PetscInt, nnz::Vector{PetscInt}) 
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
function MatSeqSBAIJSetPreallocation(petsclib::PetscLibType, B::PetscMat, bs::PetscInt, nz::PetscInt, nnz::Vector{PetscInt}) end

@for_petsc function MatSeqSBAIJSetPreallocation(petsclib::$UnionPetscLib, B::PetscMat, bs::$PetscInt, nz::$PetscInt, nnz::Vector{$PetscInt} )

    @chk ccall(
               (:MatSeqSBAIJSetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, Ptr{$PetscInt}),
               B, bs, nz, nnz,
              )


	return nothing
end 

"""
	MatSeqSBAIJSetPreallocationCSR(petsclib::PetscLibType,B::PetscMat, bs::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, v::Vector{PetscScalar}) 
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
function MatSeqSBAIJSetPreallocationCSR(petsclib::PetscLibType, B::PetscMat, bs::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, v::Vector{PetscScalar}) end

@for_petsc function MatSeqSBAIJSetPreallocationCSR(petsclib::$UnionPetscLib, B::PetscMat, bs::$PetscInt, i::Vector{$PetscInt}, j::Vector{$PetscInt}, v::Vector{$PetscScalar} )

    @chk ccall(
               (:MatSeqSBAIJSetPreallocationCSR, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}),
               B, bs, i, j, v,
              )


	return nothing
end 

"""
	A::PetscMat = MatCreateSeqSBAIJ(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Vector{PetscInt}) 
Creates a sparse symmetric matrix in (block
compressed row) `MATSEQSBAIJ` format.  For good matrix assembly performance the
user should preallocate the matrix storage by setting the parameter `nz`
(or the array `nnz`).

Collective

Input Parameters:
- `comm` - MPI communicator, set to `PETSC_COMM_SELF`
- `bs`   - size of block, the blocks are ALWAYS square. One can use `MatSetBlockSizes()` to set a different row and column blocksize but the row blocksize always defines the size of the blocks. The column blocksize sets the blocksize of the vectors obtained with MatCreateVecs()
- `m`    - number of rows
- `n`    - number of columns
- `nz`   - number of block nonzeros per block row (same for all rows)
- `nnz`  - array containing the number of block nonzeros in the upper triangular plus diagonal portion of each block (possibly different for each block row) or `NULL`

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
function MatCreateSeqSBAIJ(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Union{Ptr,Vector{PetscInt}}) end

@for_petsc function MatCreateSeqSBAIJ(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, m::$PetscInt, n::$PetscInt, nz::$PetscInt, nnz::Union{Ptr,Vector{$PetscInt}} )
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
	MatReorderingSeqSBAIJ(petsclib::PetscLibType,A::PetscMat, perm::IS) 

# External Links
$(_doc_external("Mat/MatReorderingSeqSBAIJ"))
"""
function MatReorderingSeqSBAIJ(petsclib::PetscLibType, A::PetscMat, perm::IS) end

@for_petsc function MatReorderingSeqSBAIJ(petsclib::$UnionPetscLib, A::PetscMat, perm::IS )

    @chk ccall(
               (:MatReorderingSeqSBAIJ, $petsc_library),
               PetscErrorCode,
               (CMat, CIS),
               A, perm,
              )


	return nothing
end 

"""
	MatKAIJGetAIJ(petsclib::PetscLibType,A::PetscMat, B::PetscMat) 
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
function MatKAIJGetAIJ(petsclib::PetscLibType, A::PetscMat, B::PetscMat) end

@for_petsc function MatKAIJGetAIJ(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat )
	B_ = Ref(B.ptr)

    @chk ccall(
               (:MatKAIJGetAIJ, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, B_,
              )

	B.ptr = C_NULL

	return nothing
end 

"""
	m::PetscInt,n::PetscInt,S::Vector{PetscScalar} = MatKAIJGetS(petsclib::PetscLibType,A::PetscMat) 
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
function MatKAIJGetS(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatKAIJGetS(petsclib::$UnionPetscLib, A::PetscMat )
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
	S = unsafe_wrap(Array, S_[], VecGetLocalSize(petsclib, x); own = false)

	return m,n,S
end 

"""
	m::PetscInt,n::PetscInt,S::Vector{PetscScalar} = MatKAIJGetSRead(petsclib::PetscLibType,A::PetscMat) 
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
function MatKAIJGetSRead(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatKAIJGetSRead(petsclib::$UnionPetscLib, A::PetscMat )
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
	S = unsafe_wrap(Array, S_[], VecGetLocalSize(petsclib, x); own = false)

	return m,n,S
end 

"""
	S::Vector{PetscScalar} = MatKAIJRestoreS(petsclib::PetscLibType,A::PetscMat) 
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
function MatKAIJRestoreS(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatKAIJRestoreS(petsclib::$UnionPetscLib, A::PetscMat )
	S_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatKAIJRestoreS, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, S_,
              )

	S = unsafe_wrap(Array, S_[], VecGetLocalSize(petsclib, x); own = false)

	return S
end 

"""
	S::Vector{PetscScalar} = MatKAIJRestoreSRead(petsclib::PetscLibType,A::PetscMat) 
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
function MatKAIJRestoreSRead(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatKAIJRestoreSRead(petsclib::$UnionPetscLib, A::PetscMat )
	S_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatKAIJRestoreSRead, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, S_,
              )

	S = unsafe_wrap(Array, S_[], VecGetLocalSize(petsclib, x); own = false)

	return S
end 

"""
	m::PetscInt,n::PetscInt,T::Vector{PetscScalar} = MatKAIJGetT(petsclib::PetscLibType,A::PetscMat) 
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
function MatKAIJGetT(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatKAIJGetT(petsclib::$UnionPetscLib, A::PetscMat )
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
	T = unsafe_wrap(Array, T_[], VecGetLocalSize(petsclib, x); own = false)

	return m,n,T
end 

"""
	m::PetscInt,n::PetscInt,T::Vector{PetscScalar} = MatKAIJGetTRead(petsclib::PetscLibType,A::PetscMat) 
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
function MatKAIJGetTRead(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatKAIJGetTRead(petsclib::$UnionPetscLib, A::PetscMat )
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
	T = unsafe_wrap(Array, T_[], VecGetLocalSize(petsclib, x); own = false)

	return m,n,T
end 

"""
	T::Vector{PetscScalar} = MatKAIJRestoreT(petsclib::PetscLibType,A::PetscMat) 
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
function MatKAIJRestoreT(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatKAIJRestoreT(petsclib::$UnionPetscLib, A::PetscMat )
	T_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatKAIJRestoreT, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, T_,
              )

	T = unsafe_wrap(Array, T_[], VecGetLocalSize(petsclib, x); own = false)

	return T
end 

"""
	T::Vector{PetscScalar} = MatKAIJRestoreTRead(petsclib::PetscLibType,A::PetscMat) 
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
function MatKAIJRestoreTRead(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatKAIJRestoreTRead(petsclib::$UnionPetscLib, A::PetscMat )
	T_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatKAIJRestoreTRead, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, T_,
              )

	T = unsafe_wrap(Array, T_[], VecGetLocalSize(petsclib, x); own = false)

	return T
end 

"""
	MatKAIJSetAIJ(petsclib::PetscLibType,A::PetscMat, B::PetscMat) 
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
function MatKAIJSetAIJ(petsclib::PetscLibType, A::PetscMat, B::PetscMat) end

@for_petsc function MatKAIJSetAIJ(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat )

    @chk ccall(
               (:MatKAIJSetAIJ, $petsc_library),
               PetscErrorCode,
               (CMat, CMat),
               A, B,
              )


	return nothing
end 

"""
	MatKAIJSetS(petsclib::PetscLibType,A::PetscMat, p::PetscInt, q::PetscInt, S::Vector{PetscScalar}) 
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
function MatKAIJSetS(petsclib::PetscLibType, A::PetscMat, p::PetscInt, q::PetscInt, S::Vector{PetscScalar}) end

@for_petsc function MatKAIJSetS(petsclib::$UnionPetscLib, A::PetscMat, p::$PetscInt, q::$PetscInt, S::Vector{$PetscScalar} )

    @chk ccall(
               (:MatKAIJSetS, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, Ptr{$PetscScalar}),
               A, p, q, S,
              )


	return nothing
end 

"""
	identity::PetscBool = MatKAIJGetScaledIdentity(petsclib::PetscLibType,A::PetscMat) 
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
function MatKAIJGetScaledIdentity(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatKAIJGetScaledIdentity(petsclib::$UnionPetscLib, A::PetscMat )
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
	MatKAIJSetT(petsclib::PetscLibType,A::PetscMat, p::PetscInt, q::PetscInt, T::Vector{PetscScalar}) 
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
function MatKAIJSetT(petsclib::PetscLibType, A::PetscMat, p::PetscInt, q::PetscInt, T::Vector{PetscScalar}) end

@for_petsc function MatKAIJSetT(petsclib::$UnionPetscLib, A::PetscMat, p::$PetscInt, q::$PetscInt, T::Vector{$PetscScalar} )

    @chk ccall(
               (:MatKAIJSetT, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, Ptr{$PetscScalar}),
               A, p, q, T,
              )


	return nothing
end 

"""
	kaij::PetscMat = MatCreateKAIJ(petsclib::PetscLibType,A::PetscMat, p::PetscInt, q::PetscInt, S::Vector{PetscScalar}, T::Vector{PetscScalar}) 
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
function MatCreateKAIJ(petsclib::PetscLibType, A::PetscMat, p::PetscInt, q::PetscInt, S::Union{Ptr,Vector{PetscScalar}}, T::Union{Ptr,Vector{PetscScalar}}) end

@for_petsc function MatCreateKAIJ(petsclib::$UnionPetscLib, A::PetscMat, p::$PetscInt, q::$PetscInt, S::Union{Ptr,Vector{$PetscScalar}}, T::Union{Ptr,Vector{$PetscScalar}} )
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
	B::PetscMat = MatMPIAdjCreateNonemptySubcommMat(petsclib::PetscLibType,A::PetscMat) 
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
function MatMPIAdjCreateNonemptySubcommMat(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatMPIAdjCreateNonemptySubcommMat(petsclib::$UnionPetscLib, A::PetscMat )
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
	MatMPIAdjToSeq(petsclib::PetscLibType,A::PetscMat, B::PetscMat) 
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
function MatMPIAdjToSeq(petsclib::PetscLibType, A::PetscMat, B::PetscMat) end

@for_petsc function MatMPIAdjToSeq(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat )
	B_ = Ref(B.ptr)

    @chk ccall(
               (:MatMPIAdjToSeq, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, B_,
              )

	B.ptr = C_NULL

	return nothing
end 

"""
	MatMPIAdjToSeqRankZero(petsclib::PetscLibType,A::PetscMat, B::PetscMat) 
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
function MatMPIAdjToSeqRankZero(petsclib::PetscLibType, A::PetscMat, B::PetscMat) end

@for_petsc function MatMPIAdjToSeqRankZero(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat )
	B_ = Ref(B.ptr)

    @chk ccall(
               (:MatMPIAdjToSeqRankZero, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, B_,
              )

	B.ptr = C_NULL

	return nothing
end 

"""
	i::PetscInt,j::PetscInt,values::PetscInt = MatMPIAdjSetPreallocation(petsclib::PetscLibType,B::PetscMat) 
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
function MatMPIAdjSetPreallocation(petsclib::PetscLibType, B::PetscMat) end

@for_petsc function MatMPIAdjSetPreallocation(petsclib::$UnionPetscLib, B::PetscMat )
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
	A::PetscMat = MatCreateMPIAdj(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, N::PetscInt, i::Vector{PetscInt}, j::Vector{PetscInt}, values::Vector{PetscInt}) 
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
	MatScatterGetVecScatter(petsclib::PetscLibType,mat::PetscMat, scatter::VecScatter) 
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
function MatScatterGetVecScatter(petsclib::PetscLibType, mat::PetscMat, scatter::VecScatter) end

@for_petsc function MatScatterGetVecScatter(petsclib::$UnionPetscLib, mat::PetscMat, scatter::VecScatter )

    @chk ccall(
               (:MatScatterGetVecScatter, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{VecScatter}),
               mat, scatter,
              )


	return nothing
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
	MatScatterSetVecScatter(petsclib::PetscLibType,mat::PetscMat, scatter::VecScatter) 
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
function MatScatterSetVecScatter(petsclib::PetscLibType, mat::PetscMat, scatter::VecScatter) end

@for_petsc function MatScatterSetVecScatter(petsclib::$UnionPetscLib, mat::PetscMat, scatter::VecScatter )

    @chk ccall(
               (:MatScatterSetVecScatter, $petsc_library),
               PetscErrorCode,
               (CMat, VecScatter),
               mat, scatter,
              )


	return nothing
end 

"""
	C::PetscMat = MatCreateCentering(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt, N::PetscInt) 
Creates a new matrix object that implements the (symmetric and idempotent) centering matrix,  I

Collective

Input Parameters:
- `comm` - MPI communicator
- `n`    - number of local rows (or `PETSC_DECIDE` to have calculated if `N` is given) This value should be the same as the local size used in creating the `y` vector for the matrix-vector product y = Ax.
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
	MatDiagonalGetDiagonal(petsclib::PetscLibType,A::PetscMat, diag::PetscVec) 
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
function MatDiagonalGetDiagonal(petsclib::PetscLibType, A::PetscMat, diag::PetscVec) end

@for_petsc function MatDiagonalGetDiagonal(petsclib::$UnionPetscLib, A::PetscMat, diag::PetscVec )
	diag_ = Ref(diag.ptr)

    @chk ccall(
               (:MatDiagonalGetDiagonal, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CVec}),
               A, diag_,
              )

	diag.ptr = C_NULL

	return nothing
end 

"""
	MatDiagonalRestoreDiagonal(petsclib::PetscLibType,A::PetscMat, diag::PetscVec) 
Restore the diagonal of a `MATDIAGONAL`

Input Parameters:
- `A`    - the `MATDIAGONAL`
- `diag` - the `Vec` obtained from `MatDiagonalGetDiagonal()`

Level: developer

-seealso: [](ch_matrices), `MATDIAGONAL`, `MatCreateDiagonal()`, `MatDiagonalGetDiagonal()`

# External Links
$(_doc_external("Mat/MatDiagonalRestoreDiagonal"))
"""
function MatDiagonalRestoreDiagonal(petsclib::PetscLibType, A::PetscMat, diag::PetscVec) end

@for_petsc function MatDiagonalRestoreDiagonal(petsclib::$UnionPetscLib, A::PetscMat, diag::PetscVec )
	diag_ = Ref(diag.ptr)

    @chk ccall(
               (:MatDiagonalRestoreDiagonal, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CVec}),
               A, diag_,
              )

	diag.ptr = C_NULL

	return nothing
end 

"""
	MatDiagonalGetInverseDiagonal(petsclib::PetscLibType,A::PetscMat, inv_diag::PetscVec) 
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
function MatDiagonalGetInverseDiagonal(petsclib::PetscLibType, A::PetscMat, inv_diag::PetscVec) end

@for_petsc function MatDiagonalGetInverseDiagonal(petsclib::$UnionPetscLib, A::PetscMat, inv_diag::PetscVec )
	inv_diag_ = Ref(inv_diag.ptr)

    @chk ccall(
               (:MatDiagonalGetInverseDiagonal, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CVec}),
               A, inv_diag_,
              )

	inv_diag.ptr = C_NULL

	return nothing
end 

"""
	MatDiagonalRestoreInverseDiagonal(petsclib::PetscLibType,A::PetscMat, inv_diag::PetscVec) 
Restore the inverse diagonal of a `MATDIAGONAL`

Input Parameters:
- `A`        - the `MATDIAGONAL`
- `inv_diag` - the `Vec` obtained from `MatDiagonalGetInverseDiagonal()`

Level: developer

-seealso: [](ch_matrices), `MATDIAGONAL`, `MatCreateDiagonal()`, `MatDiagonalGetInverseDiagonal()`

# External Links
$(_doc_external("Mat/MatDiagonalRestoreInverseDiagonal"))
"""
function MatDiagonalRestoreInverseDiagonal(petsclib::PetscLibType, A::PetscMat, inv_diag::PetscVec) end

@for_petsc function MatDiagonalRestoreInverseDiagonal(petsclib::$UnionPetscLib, A::PetscMat, inv_diag::PetscVec )
	inv_diag_ = Ref(inv_diag.ptr)

    @chk ccall(
               (:MatDiagonalRestoreInverseDiagonal, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CVec}),
               A, inv_diag_,
              )

	inv_diag.ptr = C_NULL

	return nothing
end 

"""
	J::PetscMat = MatCreateDiagonal(petsclib::PetscLibType,diag::PetscVec) 
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
function MatCreateDiagonal(petsclib::PetscLibType, diag::PetscVec) end

@for_petsc function MatCreateDiagonal(petsclib::$UnionPetscLib, diag::PetscVec )
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
	MatNormalGetMat(petsclib::PetscLibType,A::PetscMat, M::PetscMat) 
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
function MatNormalGetMat(petsclib::PetscLibType, A::PetscMat, M::PetscMat) end

@for_petsc function MatNormalGetMat(petsclib::$UnionPetscLib, A::PetscMat, M::PetscMat )
	M_ = Ref(M.ptr)

    @chk ccall(
               (:MatNormalGetMat, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, M_,
              )

	M.ptr = C_NULL

	return nothing
end 

"""
	N::PetscMat = MatCreateNormal(petsclib::PetscLibType,A::PetscMat) 
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
function MatCreateNormal(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatCreateNormal(petsclib::$UnionPetscLib, A::PetscMat )
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
	MatNormalHermitianGetMat(petsclib::PetscLibType,A::PetscMat, M::PetscMat) 
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
function MatNormalHermitianGetMat(petsclib::PetscLibType, A::PetscMat, M::PetscMat) end

@for_petsc function MatNormalHermitianGetMat(petsclib::$UnionPetscLib, A::PetscMat, M::PetscMat )
	M_ = Ref(M.ptr)

    @chk ccall(
               (:MatNormalHermitianGetMat, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, M_,
              )

	M.ptr = C_NULL

	return nothing
end 

"""
	N::PetscMat = MatCreateNormalHermitian(petsclib::PetscLibType,A::PetscMat) 
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
function MatCreateNormalHermitian(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatCreateNormalHermitian(petsclib::$UnionPetscLib, A::PetscMat )
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
	MatBlockMatSetPreallocation(petsclib::PetscLibType,B::PetscMat, bs::PetscInt, nz::PetscInt, nnz::Vector{PetscInt}) 
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
function MatBlockMatSetPreallocation(petsclib::PetscLibType, B::PetscMat, bs::PetscInt, nz::PetscInt, nnz::Vector{PetscInt}) end

@for_petsc function MatBlockMatSetPreallocation(petsclib::$UnionPetscLib, B::PetscMat, bs::$PetscInt, nz::$PetscInt, nnz::Vector{$PetscInt} )

    @chk ccall(
               (:MatBlockMatSetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, Ptr{$PetscInt}),
               B, bs, nz, nnz,
              )


	return nothing
end 

"""
	A::PetscMat = MatCreateBlockMat(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, bs::PetscInt, nz::PetscInt, nnz::Union{Ptr,Vector{PetscInt}}) 
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
function MatCreateBlockMat(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, bs::PetscInt, nz::PetscInt, nnz::Union{Ptr,Vector{PetscInt}}) end

@for_petsc function MatCreateBlockMat(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, bs::$PetscInt, nz::$PetscInt, nnz::Union{Ptr,Vector{$PetscInt}} )
	A_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateBlockMat, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{CMat}),
               comm, m, n, bs, nz, nnz, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
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
	x::PetscVec,y::PetscVec,z::PetscVec = MatCreateVecsFFTW(petsclib::PetscLibType,A::PetscMat) 
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
function MatCreateVecsFFTW(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatCreateVecsFFTW(petsclib::$UnionPetscLib, A::PetscMat )
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
	MatHtoolSetKernel(petsclib::PetscLibType,A::PetscMat, kernel::MatHtoolKernelFn, kernelctx::Cvoid) 

# External Links
$(_doc_external("Mat/MatHtoolSetKernel"))
"""
function MatHtoolSetKernel(petsclib::PetscLibType, A::PetscMat, kernel::MatHtoolKernelFn, kernelctx::Cvoid) end

@for_petsc function MatHtoolSetKernel(petsclib::$UnionPetscLib, A::PetscMat, kernel::MatHtoolKernelFn, kernelctx::Cvoid )

    @chk ccall(
               (:MatHtoolSetKernel, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatHtoolKernelFn}, Ptr{Cvoid}),
               A, kernel, kernelctx,
              )


	return nothing
end 

"""
	MatHtoolGetPermutationSource(petsclib::PetscLibType,A::PetscMat, is::IS) 

# External Links
$(_doc_external("Mat/MatHtoolGetPermutationSource"))
"""
function MatHtoolGetPermutationSource(petsclib::PetscLibType, A::PetscMat, is::IS) end

@for_petsc function MatHtoolGetPermutationSource(petsclib::$UnionPetscLib, A::PetscMat, is::IS )

    @chk ccall(
               (:MatHtoolGetPermutationSource, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CIS}),
               A, is,
              )


	return nothing
end 

"""
	MatHtoolGetPermutationTarget(petsclib::PetscLibType,A::PetscMat, is::IS) 

# External Links
$(_doc_external("Mat/MatHtoolGetPermutationTarget"))
"""
function MatHtoolGetPermutationTarget(petsclib::PetscLibType, A::PetscMat, is::IS) end

@for_petsc function MatHtoolGetPermutationTarget(petsclib::$UnionPetscLib, A::PetscMat, is::IS )

    @chk ccall(
               (:MatHtoolGetPermutationTarget, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CIS}),
               A, is,
              )


	return nothing
end 

"""
	MatHtoolUsePermutation(petsclib::PetscLibType,A::PetscMat, use::PetscBool) 

# External Links
$(_doc_external("Mat/MatHtoolUsePermutation"))
"""
function MatHtoolUsePermutation(petsclib::PetscLibType, A::PetscMat, use::PetscBool) end

@for_petsc function MatHtoolUsePermutation(petsclib::$UnionPetscLib, A::PetscMat, use::PetscBool )

    @chk ccall(
               (:MatHtoolUsePermutation, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               A, use,
              )


	return nothing
end 

"""
	MatHtoolUseRecompression(petsclib::PetscLibType,A::PetscMat, use::PetscBool) 

# External Links
$(_doc_external("Mat/MatHtoolUseRecompression"))
"""
function MatHtoolUseRecompression(petsclib::PetscLibType, A::PetscMat, use::PetscBool) end

@for_petsc function MatHtoolUseRecompression(petsclib::$UnionPetscLib, A::PetscMat, use::PetscBool )

    @chk ccall(
               (:MatHtoolUseRecompression, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               A, use,
              )


	return nothing
end 

"""
	kernel::MatHtoolKernelFn,kernelctx::Cvoid,B::PetscMat = MatCreateHtoolFromKernel(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, spacedim::PetscInt, coords_target::Vector{PetscReal}, coords_source::Vector{PetscReal}) 

# External Links
$(_doc_external("Mat/MatCreateHtoolFromKernel"))
"""
function MatCreateHtoolFromKernel(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, spacedim::PetscInt, coords_target::Vector{PetscReal}, coords_source::Vector{PetscReal}) end

@for_petsc function MatCreateHtoolFromKernel(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, spacedim::$PetscInt, coords_target::Vector{$PetscReal}, coords_source::Vector{$PetscReal} )
	kernel_ = Ref{MatHtoolKernelFn}()
	kernelctx_ = Ref{Cvoid}()
	B_ = Ref{CMat}()

    @chk ccall(
               (:MatCreateHtoolFromKernel, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{MatHtoolKernelFn}, Ptr{Cvoid}, Ptr{CMat}),
               comm, m, n, M, N, spacedim, coords_target, coords_source, kernel_, kernelctx_, B_,
              )

	kernel = kernel_[]
	kernelctx = kernelctx_[]
	B = PetscMat(B_[], petsclib)

	return kernel,kernelctx,B
end 

"""
	MatHYPRESetPreallocation(petsclib::PetscLibType,A::PetscMat, dnz::PetscInt, dnnz::Vector{PetscInt}, onz::PetscInt, onnz::Vector{PetscInt}) 
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
function MatHYPRESetPreallocation(petsclib::PetscLibType, A::PetscMat, dnz::PetscInt, dnnz::Vector{PetscInt}, onz::PetscInt, onnz::Vector{PetscInt}) end

@for_petsc function MatHYPRESetPreallocation(petsclib::$UnionPetscLib, A::PetscMat, dnz::$PetscInt, dnnz::Vector{$PetscInt}, onz::$PetscInt, onnz::Vector{$PetscInt} )

    @chk ccall(
               (:MatHYPRESetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}),
               A, dnz, dnnz, onz, onnz,
              )


	return nothing
end 

"""
	MatHYPREGetParCSR(petsclib::PetscLibType,A::PetscMat, parcsr::hypre_ParCSRMatrix) 
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
function MatHYPREGetParCSR(petsclib::PetscLibType, A::PetscMat, parcsr::hypre_ParCSRMatrix) end

@for_petsc function MatHYPREGetParCSR(petsclib::$UnionPetscLib, A::PetscMat, parcsr::hypre_ParCSRMatrix )

    @chk ccall(
               (:MatHYPREGetParCSR, $petsc_library),
               PetscErrorCode,
               (CMat, hypre_ParCSRMatrix),
               A, parcsr,
              )


	return nothing
end 

"""
	newmat::PetscMat = MatCreateSubMatrixVirtual(petsclib::PetscLibType,A::PetscMat, isrow::IS, iscol::IS) 
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
function MatCreateSubMatrixVirtual(petsclib::PetscLibType, A::PetscMat, isrow::IS, iscol::IS) end

@for_petsc function MatCreateSubMatrixVirtual(petsclib::$UnionPetscLib, A::PetscMat, isrow::IS, iscol::IS )
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
	MatSubMatrixVirtualUpdate(petsclib::PetscLibType,N::PetscMat, A::PetscMat, isrow::IS, iscol::IS) 
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
function MatSubMatrixVirtualUpdate(petsclib::PetscLibType, N::PetscMat, A::PetscMat, isrow::IS, iscol::IS) end

@for_petsc function MatSubMatrixVirtualUpdate(petsclib::$UnionPetscLib, N::PetscMat, A::PetscMat, isrow::IS, iscol::IS )

    @chk ccall(
               (:MatSubMatrixVirtualUpdate, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CIS, CIS),
               N, A, isrow, iscol,
              )


	return nothing
end 

"""
	MatNestGetSubMat(petsclib::PetscLibType,A::PetscMat, idxm::PetscInt, jdxm::PetscInt, sub::PetscMat) 
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
function MatNestGetSubMat(petsclib::PetscLibType, A::PetscMat, idxm::PetscInt, jdxm::PetscInt, sub::PetscMat) end

@for_petsc function MatNestGetSubMat(petsclib::$UnionPetscLib, A::PetscMat, idxm::$PetscInt, jdxm::$PetscInt, sub::PetscMat )
	sub_ = Ref(sub.ptr)

    @chk ccall(
               (:MatNestGetSubMat, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, Ptr{CMat}),
               A, idxm, jdxm, sub_,
              )

	sub.ptr = C_NULL

	return nothing
end 

"""
	MatNestSetSubMat(petsclib::PetscLibType,A::PetscMat, idxm::PetscInt, jdxm::PetscInt, sub::PetscMat) 
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
function MatNestSetSubMat(petsclib::PetscLibType, A::PetscMat, idxm::PetscInt, jdxm::PetscInt, sub::PetscMat) end

@for_petsc function MatNestSetSubMat(petsclib::$UnionPetscLib, A::PetscMat, idxm::$PetscInt, jdxm::$PetscInt, sub::PetscMat )

    @chk ccall(
               (:MatNestSetSubMat, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, CMat),
               A, idxm, jdxm, sub,
              )


	return nothing
end 

"""
	M::PetscInt,N::PetscInt = MatNestGetSubMats(petsclib::PetscLibType,A::PetscMat, mat::PetscMat) 
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
function MatNestGetSubMats(petsclib::PetscLibType, A::PetscMat, mat::PetscMat) end

@for_petsc function MatNestGetSubMats(petsclib::$UnionPetscLib, A::PetscMat, mat::PetscMat )
	M_ = Ref{$PetscInt}()
	N_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatNestGetSubMats, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}, Ptr{$PetscInt}, CMat),
               A, M_, N_, mat,
              )

	M = M_[]
	N = N_[]

	return M,N
end 

"""
	M::PetscInt,N::PetscInt = MatNestGetSize(petsclib::PetscLibType,A::PetscMat) 
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
function MatNestGetSize(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatNestGetSize(petsclib::$UnionPetscLib, A::PetscMat )
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
	MatNestGetISs(petsclib::PetscLibType,A::PetscMat, rows::Vector{IS}, cols::Vector{IS}) 
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
function MatNestGetISs(petsclib::PetscLibType, A::PetscMat, rows::Vector{IS}, cols::Vector{IS}) end

@for_petsc function MatNestGetISs(petsclib::$UnionPetscLib, A::PetscMat, rows::Vector{IS}, cols::Vector{IS} )

    @chk ccall(
               (:MatNestGetISs, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CIS}, Ptr{CIS}),
               A, rows, cols,
              )


	return nothing
end 

"""
	MatNestGetLocalISs(petsclib::PetscLibType,A::PetscMat, rows::Vector{IS}, cols::Vector{IS}) 
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
function MatNestGetLocalISs(petsclib::PetscLibType, A::PetscMat, rows::Vector{IS}, cols::Vector{IS}) end

@for_petsc function MatNestGetLocalISs(petsclib::$UnionPetscLib, A::PetscMat, rows::Vector{IS}, cols::Vector{IS} )

    @chk ccall(
               (:MatNestGetLocalISs, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CIS}, Ptr{CIS}),
               A, rows, cols,
              )


	return nothing
end 

"""
	MatNestSetVecType(petsclib::PetscLibType,A::PetscMat, vtype::VecType) 
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
function MatNestSetVecType(petsclib::PetscLibType, A::PetscMat, vtype::VecType) end

@for_petsc function MatNestSetVecType(petsclib::$UnionPetscLib, A::PetscMat, vtype::VecType )

    @chk ccall(
               (:MatNestSetVecType, $petsc_library),
               PetscErrorCode,
               (CMat, VecType),
               A, vtype,
              )


	return nothing
end 

"""
	MatNestSetSubMats(petsclib::PetscLibType,A::PetscMat, nr::PetscInt, is_row::Vector{IS}, nc::PetscInt, is_col::Vector{IS}, a::Vector{PetscMat}) 
Sets the nested submatrices in a `MATNEST`

Collective

Input Parameters:
- `A`      - `MATNEST` matrix
- `nr`     - number of nested row blocks
- `is_row` - index sets for each nested row block, or `NULL` to make contiguous
- `nc`     - number of nested column blocks
- `is_col` - index sets for each nested column block, or `NULL` to make contiguous
- `a`      - array of  nr \times nc submatrices, or `NULL`

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATNEST`, `MatCreateNest()`, `MatNestSetSubMat()`, `MatNestGetSubMat()`, `MatNestGetSubMats()`

# External Links
$(_doc_external("Mat/MatNestSetSubMats"))
"""
function MatNestSetSubMats(petsclib::PetscLibType, A::PetscMat, nr::PetscInt, is_row::Vector{IS}, nc::PetscInt, is_col::Vector{IS}, a::Vector{PetscMat}) end

@for_petsc function MatNestSetSubMats(petsclib::$UnionPetscLib, A::PetscMat, nr::$PetscInt, is_row::Vector{IS}, nc::$PetscInt, is_col::Vector{IS}, a::Vector{PetscMat} )

    @chk ccall(
               (:MatNestSetSubMats, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CIS}, $PetscInt, Ptr{CIS}, Ptr{CMat}),
               A, nr, is_row, nc, is_col, a,
              )


	return nothing
end 

"""
	B::PetscMat = MatCreateNest(petsclib::PetscLibType,comm::MPI_Comm, nr::PetscInt, is_row::Vector{IS}, nc::PetscInt, is_col::Vector{IS}, a::Vector{PetscMat}) 
Creates a new `MATNEST` matrix containing several nested submatrices, each stored separately

Collective

Input Parameters:
- `comm`   - Communicator for the new `MATNEST`
- `nr`     - number of nested row blocks
- `is_row` - index sets for each nested row block, or `NULL` to make contiguous
- `nc`     - number of nested column blocks
- `is_col` - index sets for each nested column block, or `NULL` to make contiguous
- `a`      - array of nr \times nc submatrices, empty submatrices can be passed using `NULL`

Output Parameter:
- `B` - new matrix

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MATNEST`, `MatCreate()`, `VecCreateNest()`, `DMCreateMatrix()`, `MatNestSetSubMat()`,
`MatNestGetSubMat()`, `MatNestGetLocalISs()`, `MatNestGetSize()`,
`MatNestGetISs()`, `MatNestSetSubMats()`, `MatNestGetSubMats()`

# External Links
$(_doc_external("Mat/MatCreateNest"))
"""
function MatCreateNest(petsclib::PetscLibType, comm::MPI_Comm, nr::PetscInt, is_row::Vector{IS}, nc::PetscInt, is_col::Vector{IS}, a::Vector{PetscMat}) end

@for_petsc function MatCreateNest(petsclib::$UnionPetscLib, comm::MPI_Comm, nr::$PetscInt, is_row::Vector{IS}, nc::$PetscInt, is_col::Vector{IS}, a::Vector{PetscMat} )
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
	MatMPIAIJGetNumberNonzeros(petsclib::PetscLibType,A::PetscMat, nz::PetscCount) 
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
function MatMPIAIJGetNumberNonzeros(petsclib::PetscLibType, A::PetscMat, nz::PetscCount) end

@for_petsc function MatMPIAIJGetNumberNonzeros(petsclib::$UnionPetscLib, A::PetscMat, nz::PetscCount )

    @chk ccall(
               (:MatMPIAIJGetNumberNonzeros, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PetscCount}),
               A, nz,
              )


	return nothing
end 

"""
	MatMPIAIJSetUseScalableIncreaseOverlap(petsclib::PetscLibType,A::PetscMat, sc::PetscBool) 
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
function MatMPIAIJSetUseScalableIncreaseOverlap(petsclib::PetscLibType, A::PetscMat, sc::PetscBool) end

@for_petsc function MatMPIAIJSetUseScalableIncreaseOverlap(petsclib::$UnionPetscLib, A::PetscMat, sc::PetscBool )

    @chk ccall(
               (:MatMPIAIJSetUseScalableIncreaseOverlap, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               A, sc,
              )


	return nothing
end 

"""
	garray::PetscInt,mat::PetscMat = MatCreateMPIAIJWithSeqAIJ(petsclib::PetscLibType,comm::MPI_Comm, M::PetscInt, N::PetscInt, A::PetscMat, B::PetscMat) 
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
function MatCreateMPIAIJWithSeqAIJ(petsclib::PetscLibType, comm::MPI_Comm, M::PetscInt, N::PetscInt, A::PetscMat, B::PetscMat) end

@for_petsc function MatCreateMPIAIJWithSeqAIJ(petsclib::$UnionPetscLib, comm::MPI_Comm, M::$PetscInt, N::$PetscInt, A::PetscMat, B::PetscMat )
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
	MatMPIAIJSetPreallocationCSR(petsclib::PetscLibType,B::PetscMat, i::Vector{PetscInt}, j::Vector{PetscInt}, v::Vector{PetscScalar}) 
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
function MatMPIAIJSetPreallocationCSR(petsclib::PetscLibType, B::PetscMat, i::Vector{PetscInt}, j::Vector{PetscInt}, v::Vector{PetscScalar}) end

@for_petsc function MatMPIAIJSetPreallocationCSR(petsclib::$UnionPetscLib, B::PetscMat, i::Vector{$PetscInt}, j::Vector{$PetscInt}, v::Vector{$PetscScalar} )

    @chk ccall(
               (:MatMPIAIJSetPreallocationCSR, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}),
               B, i, j, v,
              )


	return nothing
end 

"""
	MatMPIAIJSetPreallocation(petsclib::PetscLibType,B::PetscMat, d_nz::PetscInt, d_nnz::Vector{PetscInt}, o_nz::PetscInt, o_nnz::Vector{PetscInt}) 
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
function MatMPIAIJSetPreallocation(petsclib::PetscLibType, B::PetscMat, d_nz::PetscInt, d_nnz::Vector{PetscInt}, o_nz::PetscInt, o_nnz::Vector{PetscInt}) end

@for_petsc function MatMPIAIJSetPreallocation(petsclib::$UnionPetscLib, B::PetscMat, d_nz::$PetscInt, d_nnz::Vector{$PetscInt}, o_nz::$PetscInt, o_nnz::Vector{$PetscInt} )

    @chk ccall(
               (:MatMPIAIJSetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}),
               B, d_nz, d_nnz, o_nz, o_nnz,
              )


	return nothing
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
	MatUpdateMPIAIJWithArrays(petsclib::PetscLibType,mat::PetscMat, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, Ii::Vector{PetscInt}, J::Vector{PetscInt}, v::Vector{PetscScalar}) 
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
function MatUpdateMPIAIJWithArrays(petsclib::PetscLibType, mat::PetscMat, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, Ii::Vector{PetscInt}, J::Vector{PetscInt}, v::Vector{PetscScalar}) end

@for_petsc function MatUpdateMPIAIJWithArrays(petsclib::$UnionPetscLib, mat::PetscMat, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, Ii::Vector{$PetscInt}, J::Vector{$PetscInt}, v::Vector{$PetscScalar} )

    @chk ccall(
               (:MatUpdateMPIAIJWithArrays, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}),
               mat, m, n, M, N, Ii, J, v,
              )


	return nothing
end 

"""
	MatUpdateMPIAIJWithArray(petsclib::PetscLibType,mat::PetscMat, v::Vector{PetscScalar}) 
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
function MatUpdateMPIAIJWithArray(petsclib::PetscLibType, mat::PetscMat, v::Vector{PetscScalar}) end

@for_petsc function MatUpdateMPIAIJWithArray(petsclib::$UnionPetscLib, mat::PetscMat, v::Vector{$PetscScalar} )

    @chk ccall(
               (:MatUpdateMPIAIJWithArray, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscScalar}),
               mat, v,
              )


	return nothing
end 

"""
	A::PetscMat = MatCreateAIJ(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_nz::PetscInt, d_nnz::Vector{PetscInt}, o_nz::PetscInt, o_nnz::Vector{PetscInt}) 
Creates a sparse parallel matrix in `MATAIJ` format
(the default parallel PETSc format).  For good matrix assembly performance
the user should preallocate the matrix storage by setting the parameters
`d_nz` (or `d_nnz`) and `o_nz` (or `o_nnz`).

Collective

Input Parameters:
- `comm`  - MPI communicator
- `m`     - number of local rows (or `PETSC_DECIDE` to have calculated if M is given). This value should be the same as the local size used in creating the y vector for the matrix-vector product y = Ax.
- `n`     - This value should be the same as the local size used in creating the x vector for the matrix-vector product y = Ax. (or `PETSC_DECIDE` to have calculated if N is given) For square matrices n is almost always m.
- `M`     - number of global rows (or `PETSC_DETERMINE` to have calculated if m is given)
- `N`     - number of global columns (or `PETSC_DETERMINE` to have calculated if n is given)
- `d_nz`  - number of nonzeros per row in DIAGONAL portion of local submatrix (same value is used for all local rows)
- `d_nnz` - array containing the number of nonzeros in the various rows of the DIAGONAL portion of the local submatrix (possibly different for each row) or `NULL`, if `d_nz` is used to specify the nonzero structure. The size of this array is equal to the number of local rows, i.e 'm'.
- `o_nz`  - number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix (same value is used for all local rows).
- `o_nnz` - array containing the number of nonzeros in the various rows of the OFF-DIAGONAL portion of the local submatrix (possibly different for each row) or `NULL`, if `o_nz` is used to specify the nonzero structure. The size of this array is equal to the number of local rows, i.e 'm'.

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
function MatCreateAIJ(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_nz::PetscInt, d_nnz::Union{Ptr,Vector{PetscInt}}, o_nz::PetscInt, o_nnz::Union{Ptr,Vector{PetscInt}}) end

@for_petsc function MatCreateAIJ(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, d_nz::$PetscInt, d_nnz::Union{Ptr,Vector{$PetscInt}}, o_nz::$PetscInt, o_nnz::Union{Ptr,Vector{$PetscInt}} )
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
	colmap::Vector{PetscInt} = MatMPIAIJGetSeqAIJ(petsclib::PetscLibType,A::PetscMat, Ad::PetscMat, Ao::PetscMat) 
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
function MatMPIAIJGetSeqAIJ(petsclib::PetscLibType, A::PetscMat, Ad::PetscMat, Ao::PetscMat) end

@for_petsc function MatMPIAIJGetSeqAIJ(petsclib::$UnionPetscLib, A::PetscMat, Ad::PetscMat, Ao::PetscMat )
	Ad_ = Ref(Ad.ptr)
	Ao_ = Ref(Ao.ptr)
	colmap_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:MatMPIAIJGetSeqAIJ, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}, Ptr{CMat}, Ptr{Ptr{$PetscInt}}),
               A, Ad_, Ao_, colmap_,
              )

	Ad.ptr = C_NULL
	Ao.ptr = C_NULL
	colmap = unsafe_wrap(Array, colmap_[], VecGetLocalSize(petsclib, x); own = false)

	return colmap
end 

"""
	MatCreateMPIAIJSumSeqAIJNumeric(petsclib::PetscLibType,seqmat::PetscMat, mpimat::PetscMat) 

# External Links
$(_doc_external("Mat/MatCreateMPIAIJSumSeqAIJNumeric"))
"""
function MatCreateMPIAIJSumSeqAIJNumeric(petsclib::PetscLibType, seqmat::PetscMat, mpimat::PetscMat) end

@for_petsc function MatCreateMPIAIJSumSeqAIJNumeric(petsclib::$UnionPetscLib, seqmat::PetscMat, mpimat::PetscMat )

    @chk ccall(
               (:MatCreateMPIAIJSumSeqAIJNumeric, $petsc_library),
               PetscErrorCode,
               (CMat, CMat),
               seqmat, mpimat,
              )


	return nothing
end 

"""
	mpimat::PetscMat = MatCreateMPIAIJSumSeqAIJSymbolic(petsclib::PetscLibType,comm::MPI_Comm, seqmat::PetscMat, m::PetscInt, n::PetscInt) 

# External Links
$(_doc_external("Mat/MatCreateMPIAIJSumSeqAIJSymbolic"))
"""
function MatCreateMPIAIJSumSeqAIJSymbolic(petsclib::PetscLibType, comm::MPI_Comm, seqmat::PetscMat, m::PetscInt, n::PetscInt) end

@for_petsc function MatCreateMPIAIJSumSeqAIJSymbolic(petsclib::$UnionPetscLib, comm::MPI_Comm, seqmat::PetscMat, m::$PetscInt, n::$PetscInt )
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
	mpimat::PetscMat = MatCreateMPIAIJSumSeqAIJ(petsclib::PetscLibType,comm::MPI_Comm, seqmat::PetscMat, m::PetscInt, n::PetscInt, scall::MatReuse) 
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
function MatCreateMPIAIJSumSeqAIJ(petsclib::PetscLibType, comm::MPI_Comm, seqmat::PetscMat, m::PetscInt, n::PetscInt, scall::MatReuse) end

@for_petsc function MatCreateMPIAIJSumSeqAIJ(petsclib::$UnionPetscLib, comm::MPI_Comm, seqmat::PetscMat, m::$PetscInt, n::$PetscInt, scall::MatReuse )
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
	MatAIJGetLocalMat(petsclib::PetscLibType,A::PetscMat, A_loc::PetscMat) 
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
function MatAIJGetLocalMat(petsclib::PetscLibType, A::PetscMat, A_loc::PetscMat) end

@for_petsc function MatAIJGetLocalMat(petsclib::$UnionPetscLib, A::PetscMat, A_loc::PetscMat )
	A_loc_ = Ref(A_loc.ptr)

    @chk ccall(
               (:MatAIJGetLocalMat, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, A_loc_,
              )

	A_loc.ptr = C_NULL

	return nothing
end 

"""
	MatMPIAIJGetLocalMat(petsclib::PetscLibType,A::PetscMat, scall::MatReuse, A_loc::PetscMat) 
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
function MatMPIAIJGetLocalMat(petsclib::PetscLibType, A::PetscMat, scall::MatReuse, A_loc::PetscMat) end

@for_petsc function MatMPIAIJGetLocalMat(petsclib::$UnionPetscLib, A::PetscMat, scall::MatReuse, A_loc::PetscMat )
	A_loc_ = Ref(A_loc.ptr)

    @chk ccall(
               (:MatMPIAIJGetLocalMat, $petsc_library),
               PetscErrorCode,
               (CMat, MatReuse, Ptr{CMat}),
               A, scall, A_loc_,
              )

	A_loc.ptr = C_NULL

	return nothing
end 

"""
	MatMPIAIJGetLocalMatMerge(petsclib::PetscLibType,A::PetscMat, scall::MatReuse, glob::IS, A_loc::PetscMat) 
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
function MatMPIAIJGetLocalMatMerge(petsclib::PetscLibType, A::PetscMat, scall::MatReuse, glob::IS, A_loc::PetscMat) end

@for_petsc function MatMPIAIJGetLocalMatMerge(petsclib::$UnionPetscLib, A::PetscMat, scall::MatReuse, glob::IS, A_loc::PetscMat )
	A_loc_ = Ref(A_loc.ptr)

    @chk ccall(
               (:MatMPIAIJGetLocalMatMerge, $petsc_library),
               PetscErrorCode,
               (CMat, MatReuse, Ptr{CIS}, Ptr{CMat}),
               A, scall, glob, A_loc_,
              )

	A_loc.ptr = C_NULL

	return nothing
end 

"""
	MatMPIAIJGetLocalMatCondensed(petsclib::PetscLibType,A::PetscMat, scall::MatReuse, row::IS, col::IS, A_loc::PetscMat) 
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
function MatMPIAIJGetLocalMatCondensed(petsclib::PetscLibType, A::PetscMat, scall::MatReuse, row::IS, col::IS, A_loc::PetscMat) end

@for_petsc function MatMPIAIJGetLocalMatCondensed(petsclib::$UnionPetscLib, A::PetscMat, scall::MatReuse, row::IS, col::IS, A_loc::PetscMat )
	A_loc_ = Ref(A_loc.ptr)

    @chk ccall(
               (:MatMPIAIJGetLocalMatCondensed, $petsc_library),
               PetscErrorCode,
               (CMat, MatReuse, Ptr{CIS}, Ptr{CIS}, Ptr{CMat}),
               A, scall, row, col, A_loc_,
              )

	A_loc.ptr = C_NULL

	return nothing
end 

"""
	MatGetBrowsOfAcols(petsclib::PetscLibType,A::PetscMat, B::PetscMat, scall::MatReuse, rowb::IS, colb::IS, B_seq::PetscMat) 
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
function MatGetBrowsOfAcols(petsclib::PetscLibType, A::PetscMat, B::PetscMat, scall::MatReuse, rowb::IS, colb::IS, B_seq::PetscMat) end

@for_petsc function MatGetBrowsOfAcols(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat, scall::MatReuse, rowb::IS, colb::IS, B_seq::PetscMat )
	B_seq_ = Ref(B_seq.ptr)

    @chk ccall(
               (:MatGetBrowsOfAcols, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, MatReuse, Ptr{CIS}, Ptr{CIS}, Ptr{CMat}),
               A, B, scall, rowb, colb, B_seq_,
              )

	B_seq.ptr = C_NULL

	return nothing
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
	MatSTRUMPACKSetReordering(petsclib::PetscLibType,F::PetscMat, reordering::MatSTRUMPACKReordering) 
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
function MatSTRUMPACKSetReordering(petsclib::PetscLibType, F::PetscMat, reordering::MatSTRUMPACKReordering) end

@for_petsc function MatSTRUMPACKSetReordering(petsclib::$UnionPetscLib, F::PetscMat, reordering::MatSTRUMPACKReordering )

    @chk ccall(
               (:MatSTRUMPACKSetReordering, $petsc_library),
               PetscErrorCode,
               (CMat, MatSTRUMPACKReordering),
               F, reordering,
              )


	return nothing
end 

"""
	MatSTRUMPACKGetReordering(petsclib::PetscLibType,F::PetscMat, reordering::MatSTRUMPACKReordering) 
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
function MatSTRUMPACKGetReordering(petsclib::PetscLibType, F::PetscMat, reordering::MatSTRUMPACKReordering) end

@for_petsc function MatSTRUMPACKGetReordering(petsclib::$UnionPetscLib, F::PetscMat, reordering::MatSTRUMPACKReordering )

    @chk ccall(
               (:MatSTRUMPACKGetReordering, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatSTRUMPACKReordering}),
               F, reordering,
              )


	return nothing
end 

"""
	MatSTRUMPACKSetColPerm(petsclib::PetscLibType,F::PetscMat, cperm::PetscBool) 
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
function MatSTRUMPACKSetColPerm(petsclib::PetscLibType, F::PetscMat, cperm::PetscBool) end

@for_petsc function MatSTRUMPACKSetColPerm(petsclib::$UnionPetscLib, F::PetscMat, cperm::PetscBool )

    @chk ccall(
               (:MatSTRUMPACKSetColPerm, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               F, cperm,
              )


	return nothing
end 

"""
	cperm::PetscBool = MatSTRUMPACKGetColPerm(petsclib::PetscLibType,F::PetscMat) 
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
function MatSTRUMPACKGetColPerm(petsclib::PetscLibType, F::PetscMat) end

@for_petsc function MatSTRUMPACKGetColPerm(petsclib::$UnionPetscLib, F::PetscMat )
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
	MatSTRUMPACKSetGPU(petsclib::PetscLibType,F::PetscMat, gpu::PetscBool) 
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
function MatSTRUMPACKSetGPU(petsclib::PetscLibType, F::PetscMat, gpu::PetscBool) end

@for_petsc function MatSTRUMPACKSetGPU(petsclib::$UnionPetscLib, F::PetscMat, gpu::PetscBool )

    @chk ccall(
               (:MatSTRUMPACKSetGPU, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               F, gpu,
              )


	return nothing
end 

"""
	gpu::PetscBool = MatSTRUMPACKGetGPU(petsclib::PetscLibType,F::PetscMat) 
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
function MatSTRUMPACKGetGPU(petsclib::PetscLibType, F::PetscMat) end

@for_petsc function MatSTRUMPACKGetGPU(petsclib::$UnionPetscLib, F::PetscMat )
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
	MatSTRUMPACKSetCompression(petsclib::PetscLibType,F::PetscMat, comp::MatSTRUMPACKCompressionType) 
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
function MatSTRUMPACKSetCompression(petsclib::PetscLibType, F::PetscMat, comp::MatSTRUMPACKCompressionType) end

@for_petsc function MatSTRUMPACKSetCompression(petsclib::$UnionPetscLib, F::PetscMat, comp::MatSTRUMPACKCompressionType )

    @chk ccall(
               (:MatSTRUMPACKSetCompression, $petsc_library),
               PetscErrorCode,
               (CMat, MatSTRUMPACKCompressionType),
               F, comp,
              )


	return nothing
end 

"""
	MatSTRUMPACKGetCompression(petsclib::PetscLibType,F::PetscMat, comp::MatSTRUMPACKCompressionType) 
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
function MatSTRUMPACKGetCompression(petsclib::PetscLibType, F::PetscMat, comp::MatSTRUMPACKCompressionType) end

@for_petsc function MatSTRUMPACKGetCompression(petsclib::$UnionPetscLib, F::PetscMat, comp::MatSTRUMPACKCompressionType )

    @chk ccall(
               (:MatSTRUMPACKGetCompression, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatSTRUMPACKCompressionType}),
               F, comp,
              )


	return nothing
end 

"""
	MatSTRUMPACKSetCompRelTol(petsclib::PetscLibType,F::PetscMat, rtol::PetscReal) 
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
function MatSTRUMPACKSetCompRelTol(petsclib::PetscLibType, F::PetscMat, rtol::PetscReal) end

@for_petsc function MatSTRUMPACKSetCompRelTol(petsclib::$UnionPetscLib, F::PetscMat, rtol::$PetscReal )

    @chk ccall(
               (:MatSTRUMPACKSetCompRelTol, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscReal),
               F, rtol,
              )


	return nothing
end 

"""
	rtol::PetscReal = MatSTRUMPACKGetCompRelTol(petsclib::PetscLibType,F::PetscMat) 
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
function MatSTRUMPACKGetCompRelTol(petsclib::PetscLibType, F::PetscMat) end

@for_petsc function MatSTRUMPACKGetCompRelTol(petsclib::$UnionPetscLib, F::PetscMat )
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
	MatSTRUMPACKSetCompAbsTol(petsclib::PetscLibType,F::PetscMat, atol::PetscReal) 
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
function MatSTRUMPACKSetCompAbsTol(petsclib::PetscLibType, F::PetscMat, atol::PetscReal) end

@for_petsc function MatSTRUMPACKSetCompAbsTol(petsclib::$UnionPetscLib, F::PetscMat, atol::$PetscReal )

    @chk ccall(
               (:MatSTRUMPACKSetCompAbsTol, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscReal),
               F, atol,
              )


	return nothing
end 

"""
	atol::PetscReal = MatSTRUMPACKGetCompAbsTol(petsclib::PetscLibType,F::PetscMat) 
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
function MatSTRUMPACKGetCompAbsTol(petsclib::PetscLibType, F::PetscMat) end

@for_petsc function MatSTRUMPACKGetCompAbsTol(petsclib::$UnionPetscLib, F::PetscMat )
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
	MatSTRUMPACKSetCompLeafSize(petsclib::PetscLibType,F::PetscMat, leaf_size::PetscInt) 
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
function MatSTRUMPACKSetCompLeafSize(petsclib::PetscLibType, F::PetscMat, leaf_size::PetscInt) end

@for_petsc function MatSTRUMPACKSetCompLeafSize(petsclib::$UnionPetscLib, F::PetscMat, leaf_size::$PetscInt )

    @chk ccall(
               (:MatSTRUMPACKSetCompLeafSize, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt),
               F, leaf_size,
              )


	return nothing
end 

"""
	leaf_size::PetscInt = MatSTRUMPACKGetCompLeafSize(petsclib::PetscLibType,F::PetscMat) 
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
function MatSTRUMPACKGetCompLeafSize(petsclib::PetscLibType, F::PetscMat) end

@for_petsc function MatSTRUMPACKGetCompLeafSize(petsclib::$UnionPetscLib, F::PetscMat )
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
	MatSTRUMPACKSetGeometricNxyz(petsclib::PetscLibType,F::PetscMat, nx::PetscInt, ny::PetscInt, nz::PetscInt) 
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
function MatSTRUMPACKSetGeometricNxyz(petsclib::PetscLibType, F::PetscMat, nx::PetscInt, ny::PetscInt, nz::PetscInt) end

@for_petsc function MatSTRUMPACKSetGeometricNxyz(petsclib::$UnionPetscLib, F::PetscMat, nx::$PetscInt, ny::$PetscInt, nz::$PetscInt )

    @chk ccall(
               (:MatSTRUMPACKSetGeometricNxyz, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, $PetscInt),
               F, nx, ny, nz,
              )


	return nothing
end 

"""
	MatSTRUMPACKSetGeometricComponents(petsclib::PetscLibType,F::PetscMat, nc::PetscInt) 
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
function MatSTRUMPACKSetGeometricComponents(petsclib::PetscLibType, F::PetscMat, nc::PetscInt) end

@for_petsc function MatSTRUMPACKSetGeometricComponents(petsclib::$UnionPetscLib, F::PetscMat, nc::$PetscInt )

    @chk ccall(
               (:MatSTRUMPACKSetGeometricComponents, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt),
               F, nc,
              )


	return nothing
end 

"""
	MatSTRUMPACKSetGeometricWidth(petsclib::PetscLibType,F::PetscMat, w::PetscInt) 
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
function MatSTRUMPACKSetGeometricWidth(petsclib::PetscLibType, F::PetscMat, w::PetscInt) end

@for_petsc function MatSTRUMPACKSetGeometricWidth(petsclib::$UnionPetscLib, F::PetscMat, w::$PetscInt )

    @chk ccall(
               (:MatSTRUMPACKSetGeometricWidth, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt),
               F, w,
              )


	return nothing
end 

"""
	MatSTRUMPACKSetCompMinSepSize(petsclib::PetscLibType,F::PetscMat, min_sep_size::PetscInt) 
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
function MatSTRUMPACKSetCompMinSepSize(petsclib::PetscLibType, F::PetscMat, min_sep_size::PetscInt) end

@for_petsc function MatSTRUMPACKSetCompMinSepSize(petsclib::$UnionPetscLib, F::PetscMat, min_sep_size::$PetscInt )

    @chk ccall(
               (:MatSTRUMPACKSetCompMinSepSize, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt),
               F, min_sep_size,
              )


	return nothing
end 

"""
	min_sep_size::PetscInt = MatSTRUMPACKGetCompMinSepSize(petsclib::PetscLibType,F::PetscMat) 
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
function MatSTRUMPACKGetCompMinSepSize(petsclib::PetscLibType, F::PetscMat) end

@for_petsc function MatSTRUMPACKGetCompMinSepSize(petsclib::$UnionPetscLib, F::PetscMat )
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
	MatSTRUMPACKSetCompLossyPrecision(petsclib::PetscLibType,F::PetscMat, lossy_prec::PetscInt) 
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
function MatSTRUMPACKSetCompLossyPrecision(petsclib::PetscLibType, F::PetscMat, lossy_prec::PetscInt) end

@for_petsc function MatSTRUMPACKSetCompLossyPrecision(petsclib::$UnionPetscLib, F::PetscMat, lossy_prec::$PetscInt )

    @chk ccall(
               (:MatSTRUMPACKSetCompLossyPrecision, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt),
               F, lossy_prec,
              )


	return nothing
end 

"""
	lossy_prec::PetscInt = MatSTRUMPACKGetCompLossyPrecision(petsclib::PetscLibType,F::PetscMat) 
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
function MatSTRUMPACKGetCompLossyPrecision(petsclib::PetscLibType, F::PetscMat) end

@for_petsc function MatSTRUMPACKGetCompLossyPrecision(petsclib::$UnionPetscLib, F::PetscMat )
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
	MatSTRUMPACKSetCompButterflyLevels(petsclib::PetscLibType,F::PetscMat, bfly_lvls::PetscInt) 
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
function MatSTRUMPACKSetCompButterflyLevels(petsclib::PetscLibType, F::PetscMat, bfly_lvls::PetscInt) end

@for_petsc function MatSTRUMPACKSetCompButterflyLevels(petsclib::$UnionPetscLib, F::PetscMat, bfly_lvls::$PetscInt )

    @chk ccall(
               (:MatSTRUMPACKSetCompButterflyLevels, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt),
               F, bfly_lvls,
              )


	return nothing
end 

"""
	bfly_lvls::PetscInt = MatSTRUMPACKGetCompButterflyLevels(petsclib::PetscLibType,F::PetscMat) 
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
function MatSTRUMPACKGetCompButterflyLevels(petsclib::PetscLibType, F::PetscMat) end

@for_petsc function MatSTRUMPACKGetCompButterflyLevels(petsclib::$UnionPetscLib, F::PetscMat )
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
	diagU::PetscScalar = MatSuperluDistGetDiagU(petsclib::PetscLibType,F::PetscMat) 

# External Links
$(_doc_external("Mat/MatSuperluDistGetDiagU"))
"""
function MatSuperluDistGetDiagU(petsclib::PetscLibType, F::PetscMat) end

@for_petsc function MatSuperluDistGetDiagU(petsclib::$UnionPetscLib, F::PetscMat )
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
	A::PetscMat = MatCreateMPIAIJPERM(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_nz::PetscInt, d_nnz::Vector{PetscInt}, o_nz::PetscInt, o_nnz::Vector{PetscInt}) 
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
function MatCreateMPIAIJPERM(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, d_nz::PetscInt, d_nnz::Union{Ptr,Vector{PetscInt}}, o_nz::PetscInt, o_nnz::Union{Ptr,Vector{PetscInt}}) end

@for_petsc function MatCreateMPIAIJPERM(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, d_nz::$PetscInt, d_nnz::Union{Ptr,Vector{$PetscInt}}, o_nz::$PetscInt, o_nnz::Union{Ptr,Vector{$PetscInt}} )
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
	MatMumpsSetIcntl(petsclib::PetscLibType,F::PetscMat, icntl::PetscInt, ival::PetscInt) 
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
function MatMumpsSetIcntl(petsclib::PetscLibType, F::PetscMat, icntl::PetscInt, ival::PetscInt) end

@for_petsc function MatMumpsSetIcntl(petsclib::$UnionPetscLib, F::PetscMat, icntl::$PetscInt, ival::$PetscInt )

    @chk ccall(
               (:MatMumpsSetIcntl, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt),
               F, icntl, ival,
              )


	return nothing
end 

"""
	ival::PetscInt = MatMumpsGetIcntl(petsclib::PetscLibType,F::PetscMat, icntl::PetscInt) 
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
function MatMumpsGetIcntl(petsclib::PetscLibType, F::PetscMat, icntl::PetscInt) end

@for_petsc function MatMumpsGetIcntl(petsclib::$UnionPetscLib, F::PetscMat, icntl::$PetscInt )
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
	MatMumpsSetCntl(petsclib::PetscLibType,F::PetscMat, icntl::PetscInt, val::PetscReal) 
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
function MatMumpsSetCntl(petsclib::PetscLibType, F::PetscMat, icntl::PetscInt, val::PetscReal) end

@for_petsc function MatMumpsSetCntl(petsclib::$UnionPetscLib, F::PetscMat, icntl::$PetscInt, val::$PetscReal )

    @chk ccall(
               (:MatMumpsSetCntl, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscReal),
               F, icntl, val,
              )


	return nothing
end 

"""
	val::PetscReal = MatMumpsGetCntl(petsclib::PetscLibType,F::PetscMat, icntl::PetscInt) 
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
function MatMumpsGetCntl(petsclib::PetscLibType, F::PetscMat, icntl::PetscInt) end

@for_petsc function MatMumpsGetCntl(petsclib::$UnionPetscLib, F::PetscMat, icntl::$PetscInt )
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
	MatMumpsGetInverse(petsclib::PetscLibType,F::PetscMat, spRHS::PetscMat) 
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
function MatMumpsGetInverse(petsclib::PetscLibType, F::PetscMat, spRHS::PetscMat) end

@for_petsc function MatMumpsGetInverse(petsclib::$UnionPetscLib, F::PetscMat, spRHS::PetscMat )

    @chk ccall(
               (:MatMumpsGetInverse, $petsc_library),
               PetscErrorCode,
               (CMat, CMat),
               F, spRHS,
              )


	return nothing
end 

"""
	MatMumpsGetInverseTranspose(petsclib::PetscLibType,F::PetscMat, spRHST::PetscMat) 
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
function MatMumpsGetInverseTranspose(petsclib::PetscLibType, F::PetscMat, spRHST::PetscMat) end

@for_petsc function MatMumpsGetInverseTranspose(petsclib::$UnionPetscLib, F::PetscMat, spRHST::PetscMat )

    @chk ccall(
               (:MatMumpsGetInverseTranspose, $petsc_library),
               PetscErrorCode,
               (CMat, CMat),
               F, spRHST,
              )


	return nothing
end 

"""
	MatMumpsSetBlk(petsclib::PetscLibType,F::PetscMat, nblk::PetscInt, blkvar::Vector{PetscInt}, blkptr::Vector{PetscInt}) 
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
function MatMumpsSetBlk(petsclib::PetscLibType, F::PetscMat, nblk::PetscInt, blkvar::Vector{PetscInt}, blkptr::Vector{PetscInt}) end

@for_petsc function MatMumpsSetBlk(petsclib::$UnionPetscLib, F::PetscMat, nblk::$PetscInt, blkvar::Vector{$PetscInt}, blkptr::Vector{$PetscInt} )

    @chk ccall(
               (:MatMumpsSetBlk, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
               F, nblk, blkvar, blkptr,
              )


	return nothing
end 

"""
	ival::PetscInt = MatMumpsGetInfo(petsclib::PetscLibType,F::PetscMat, icntl::PetscInt) 
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
function MatMumpsGetInfo(petsclib::PetscLibType, F::PetscMat, icntl::PetscInt) end

@for_petsc function MatMumpsGetInfo(petsclib::$UnionPetscLib, F::PetscMat, icntl::$PetscInt )
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
	ival::PetscInt = MatMumpsGetInfog(petsclib::PetscLibType,F::PetscMat, icntl::PetscInt) 
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
function MatMumpsGetInfog(petsclib::PetscLibType, F::PetscMat, icntl::PetscInt) end

@for_petsc function MatMumpsGetInfog(petsclib::$UnionPetscLib, F::PetscMat, icntl::$PetscInt )
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
	val::PetscReal = MatMumpsGetRinfo(petsclib::PetscLibType,F::PetscMat, icntl::PetscInt) 
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
function MatMumpsGetRinfo(petsclib::PetscLibType, F::PetscMat, icntl::PetscInt) end

@for_petsc function MatMumpsGetRinfo(petsclib::$UnionPetscLib, F::PetscMat, icntl::$PetscInt )
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
	val::PetscReal = MatMumpsGetRinfog(petsclib::PetscLibType,F::PetscMat, icntl::PetscInt) 
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
function MatMumpsGetRinfog(petsclib::PetscLibType, F::PetscMat, icntl::PetscInt) end

@for_petsc function MatMumpsGetRinfog(petsclib::$UnionPetscLib, F::PetscMat, icntl::$PetscInt )
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
	size::PetscInt = MatMumpsGetNullPivots(petsclib::PetscLibType,F::PetscMat, array::PetscInt) 
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
function MatMumpsGetNullPivots(petsclib::PetscLibType, F::PetscMat, array::PetscInt) end

@for_petsc function MatMumpsGetNullPivots(petsclib::$UnionPetscLib, F::PetscMat, array::$PetscInt )
	size_ = Ref{$PetscInt}()

    @chk ccall(
               (:MatMumpsGetNullPivots, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}, $PetscInt),
               F, size_, array,
              )

	size = size_[]

	return size
end 

"""
	MatSeqAIJSetValuesLocalFast(petsclib::PetscLibType,A::PetscMat, m::PetscInt, im::Vector{PetscInt}, n::PetscInt, in::Vector{PetscInt}, v::Vector{PetscScalar}, is::InsertMode) 

# External Links
$(_doc_external("Mat/MatSeqAIJSetValuesLocalFast"))
"""
function MatSeqAIJSetValuesLocalFast(petsclib::PetscLibType, A::PetscMat, m::PetscInt, im::Vector{PetscInt}, n::PetscInt, in::Vector{PetscInt}, v::Vector{PetscScalar}, is::InsertMode) end

@for_petsc function MatSeqAIJSetValuesLocalFast(petsclib::$UnionPetscLib, A::PetscMat, m::$PetscInt, im::Vector{$PetscInt}, n::$PetscInt, in::Vector{$PetscInt}, v::Vector{$PetscScalar}, is::InsertMode )

    @chk ccall(
               (:MatSeqAIJSetValuesLocalFast, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}, InsertMode),
               A, m, im, n, in, v, is,
              )


	return nothing
end 

"""
	MatSeqAIJSetTotalPreallocation(petsclib::PetscLibType,A::PetscMat, nztotal::PetscInt) 
Sets an upper bound on the total number of expected nonzeros in the matrix.

Input Parameters:
- `A`       - the `MATSEQAIJ` matrix
- `nztotal` - bound on the number of nonzeros

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatSetOption()`, `MAT_SORTED_FULL`, `MatSetValues()`, `MatSeqAIJSetPreallocation()`

# External Links
$(_doc_external("Mat/MatSeqAIJSetTotalPreallocation"))
"""
function MatSeqAIJSetTotalPreallocation(petsclib::PetscLibType, A::PetscMat, nztotal::PetscInt) end

@for_petsc function MatSeqAIJSetTotalPreallocation(petsclib::$UnionPetscLib, A::PetscMat, nztotal::$PetscInt )

    @chk ccall(
               (:MatSeqAIJSetTotalPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt),
               A, nztotal,
              )


	return nothing
end 

"""
	indices::PetscInt = MatSeqAIJSetColumnIndices(petsclib::PetscLibType,mat::PetscMat) 
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
function MatSeqAIJSetColumnIndices(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatSeqAIJSetColumnIndices(petsclib::$UnionPetscLib, mat::PetscMat )
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
	MatStoreValues(petsclib::PetscLibType,mat::PetscMat) 
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
function MatStoreValues(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatStoreValues(petsclib::$UnionPetscLib, mat::PetscMat )

    @chk ccall(
               (:MatStoreValues, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	MatRetrieveValues(petsclib::PetscLibType,mat::PetscMat) 
Retrieves the copy of the matrix values that was stored with `MatStoreValues()`

Logically Collect

Input Parameter:
- `mat` - the matrix (currently only `MATAIJ` matrices support this option)

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatStoreValues()`

# External Links
$(_doc_external("Mat/MatRetrieveValues"))
"""
function MatRetrieveValues(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatRetrieveValues(petsclib::$UnionPetscLib, mat::PetscMat )

    @chk ccall(
               (:MatRetrieveValues, $petsc_library),
               PetscErrorCode,
               (CMat,),
               mat,
              )


	return nothing
end 

"""
	A::PetscMat = MatCreateSeqAIJ(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Vector{PetscInt}) 
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
function MatCreateSeqAIJ(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Union{Ptr,Vector{PetscInt}}) end

@for_petsc function MatCreateSeqAIJ(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, nz::$PetscInt, nnz::Union{Ptr,Vector{$PetscInt}} )
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
	MatSeqAIJSetPreallocation(petsclib::PetscLibType,B::PetscMat, nz::PetscInt, nnz::Vector{PetscInt}) 
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
function MatSeqAIJSetPreallocation(petsclib::PetscLibType, B::PetscMat, nz::PetscInt, nnz::Vector{PetscInt}) end

@for_petsc function MatSeqAIJSetPreallocation(petsclib::$UnionPetscLib, B::PetscMat, nz::$PetscInt, nnz::Vector{$PetscInt} )

    @chk ccall(
               (:MatSeqAIJSetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}),
               B, nz, nnz,
              )


	return nothing
end 

"""
	MatSeqAIJSetPreallocationCSR(petsclib::PetscLibType,B::PetscMat, i::Vector{PetscInt}, j::Vector{PetscInt}, v::Vector{PetscScalar}) 
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
function MatSeqAIJSetPreallocationCSR(petsclib::PetscLibType, B::PetscMat, i::Vector{PetscInt}, j::Vector{PetscInt}, v::Vector{PetscScalar}) end

@for_petsc function MatSeqAIJSetPreallocationCSR(petsclib::$UnionPetscLib, B::PetscMat, i::Vector{$PetscInt}, j::Vector{$PetscInt}, v::Vector{$PetscScalar} )

    @chk ccall(
               (:MatSeqAIJSetPreallocationCSR, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscScalar}),
               B, i, j, v,
              )


	return nothing
end 

"""
	MatSeqAIJKron(petsclib::PetscLibType,A::PetscMat, B::PetscMat, reuse::MatReuse, C::PetscMat) 
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
function MatSeqAIJKron(petsclib::PetscLibType, A::PetscMat, B::PetscMat, reuse::MatReuse, C::PetscMat) end

@for_petsc function MatSeqAIJKron(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat, reuse::MatReuse, C::PetscMat )
	C_ = Ref(C.ptr)

    @chk ccall(
               (:MatSeqAIJKron, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, MatReuse, Ptr{CMat}),
               A, B, reuse, C_,
              )

	C.ptr = C_NULL

	return nothing
end 

"""
	array::Vector{PetscScalar} = MatSeqAIJGetArray(petsclib::PetscLibType,A::PetscMat) 
gives read/write access to the array where the data for a `MATSEQAIJ` matrix is stored

Not Collective

Input Parameter:
- `A` - a `MATSEQAIJ` matrix

Output Parameter:
- `array` - pointer to the data

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MatSeqAIJRestoreArray()`

# External Links
$(_doc_external("Mat/MatSeqAIJGetArray"))
"""
function MatSeqAIJGetArray(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatSeqAIJGetArray(petsclib::$UnionPetscLib, A::PetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatSeqAIJGetArray, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

	array = unsafe_wrap(Array, array_[], VecGetLocalSize(petsclib, x); own = false)

	return array
end 

"""
	array::Vector{PetscScalar} = MatSeqAIJRestoreArray(petsclib::PetscLibType,A::PetscMat) 
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
function MatSeqAIJRestoreArray(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatSeqAIJRestoreArray(petsclib::$UnionPetscLib, A::PetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatSeqAIJRestoreArray, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

	array = unsafe_wrap(Array, array_[], VecGetLocalSize(petsclib, x); own = false)

	return array
end 

"""
	array::Vector{PetscScalar} = MatSeqAIJGetArrayRead(petsclib::PetscLibType,A::PetscMat) 
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
function MatSeqAIJGetArrayRead(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatSeqAIJGetArrayRead(petsclib::$UnionPetscLib, A::PetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatSeqAIJGetArrayRead, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

	array = unsafe_wrap(Array, array_[], VecGetLocalSize(petsclib, x); own = false)

	return array
end 

"""
	array::Vector{PetscScalar} = MatSeqAIJRestoreArrayRead(petsclib::PetscLibType,A::PetscMat) 
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
function MatSeqAIJRestoreArrayRead(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatSeqAIJRestoreArrayRead(petsclib::$UnionPetscLib, A::PetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatSeqAIJRestoreArrayRead, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

	array = unsafe_wrap(Array, array_[], VecGetLocalSize(petsclib, x); own = false)

	return array
end 

"""
	array::Vector{PetscScalar} = MatSeqAIJGetArrayWrite(petsclib::PetscLibType,A::PetscMat) 
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
function MatSeqAIJGetArrayWrite(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatSeqAIJGetArrayWrite(petsclib::$UnionPetscLib, A::PetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatSeqAIJGetArrayWrite, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

	array = unsafe_wrap(Array, array_[], VecGetLocalSize(petsclib, x); own = false)

	return array
end 

"""
	array::Vector{PetscScalar} = MatSeqAIJRestoreArrayWrite(petsclib::PetscLibType,A::PetscMat) 
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
function MatSeqAIJRestoreArrayWrite(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatSeqAIJRestoreArrayWrite(petsclib::$UnionPetscLib, A::PetscMat )
	array_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:MatSeqAIJRestoreArrayWrite, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

	array = unsafe_wrap(Array, array_[], VecGetLocalSize(petsclib, x); own = false)

	return array
end 

"""
	i::Vector{PetscInt},j::Vector{PetscInt},a::Vector{PetscScalar},mtype::PetscMemType = MatSeqAIJGetCSRAndMemType(petsclib::PetscLibType,mat::PetscMat) 
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
function MatSeqAIJGetCSRAndMemType(petsclib::PetscLibType, mat::PetscMat) end

@for_petsc function MatSeqAIJGetCSRAndMemType(petsclib::$UnionPetscLib, mat::PetscMat )
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

	i = unsafe_wrap(Array, i_[], VecGetLocalSize(petsclib, x); own = false)
	j = unsafe_wrap(Array, j_[], VecGetLocalSize(petsclib, x); own = false)
	a = unsafe_wrap(Array, a_[], VecGetLocalSize(petsclib, x); own = false)
	mtype = unsafe_string(mtype_[])

	return i,j,a,mtype
end 

"""
	nz::PetscInt = MatSeqAIJGetMaxRowNonzeros(petsclib::PetscLibType,A::PetscMat) 
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
function MatSeqAIJGetMaxRowNonzeros(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatSeqAIJGetMaxRowNonzeros(petsclib::$UnionPetscLib, A::PetscMat )
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
	MatSeqAIJSetType(petsclib::PetscLibType,mat::PetscMat, matype::MatType) 
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
function MatSeqAIJSetType(petsclib::PetscLibType, mat::PetscMat, matype::MatType) end

@for_petsc function MatSeqAIJSetType(petsclib::$UnionPetscLib, mat::PetscMat, matype::MatType )

    @chk ccall(
               (:MatSeqAIJSetType, $petsc_library),
               PetscErrorCode,
               (CMat, MatType),
               mat, matype,
              )


	return nothing
end 

"""
	MatSeqAIJRegister(petsclib::PetscLibType,sname::Vector{Cchar}, fnc::external) 


Not Collective, No Fortran Support

Input Parameters:
- `sname`    - name of a new user-defined matrix type, for example `MATSEQAIJCRL`
- `function` - routine to convert to subtype

Level: advanced

-seealso: [](ch_matrices), `Mat`, `MatSeqAIJRegisterAll()`

# External Links
$(_doc_external("Mat/MatSeqAIJRegister"))
"""
function MatSeqAIJRegister(petsclib::PetscLibType, sname::Vector{Cchar}, fnc::external) end

@for_petsc function MatSeqAIJRegister(petsclib::$UnionPetscLib, sname::Vector{Cchar}, fnc::external )

    @chk ccall(
               (:MatSeqAIJRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	MatInodeAdjustForInodes(petsclib::PetscLibType,A::PetscMat, rperm::IS, cperm::IS) 

# External Links
$(_doc_external("Mat/MatInodeAdjustForInodes"))
"""
function MatInodeAdjustForInodes(petsclib::PetscLibType, A::PetscMat, rperm::IS, cperm::IS) end

@for_petsc function MatInodeAdjustForInodes(petsclib::$UnionPetscLib, A::PetscMat, rperm::IS, cperm::IS )

    @chk ccall(
               (:MatInodeAdjustForInodes, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CIS}, Ptr{CIS}),
               A, rperm, cperm,
              )


	return nothing
end 

"""
	node_count::PetscInt,sizes::Vector{PetscInt},limit::PetscInt = MatInodeGetInodeSizes(petsclib::PetscLibType,A::PetscMat) 
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
function MatInodeGetInodeSizes(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatInodeGetInodeSizes(petsclib::$UnionPetscLib, A::PetscMat )
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
	sizes = unsafe_wrap(Array, sizes_[], VecGetLocalSize(petsclib, x); own = false)
	limit = limit_[]

	return node_count,sizes,limit
end 


"""
	A::PetscMat = MatCreateSeqAIJCRL(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Vector{PetscInt}) 
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
function MatCreateSeqAIJCRL(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Union{Ptr,Vector{PetscInt}}) end

@for_petsc function MatCreateSeqAIJCRL(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, nz::$PetscInt, nnz::Vector{$PetscInt} )
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
	A::PetscMat = MatCreateSeqAIJKokkos(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Vector{PetscInt}) 

# External Links
$(_doc_external("Mat/MatCreateSeqAIJKokkos"))
"""
function MatCreateSeqAIJKokkos(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Union{Ptr,Vector{PetscInt}}) end

@for_petsc function MatCreateSeqAIJKokkos(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, nz::$PetscInt, nnz::Union{Ptr,Vector{$PetscInt}} )
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
	MatSuperluSetILUDropTol(petsclib::PetscLibType,F::PetscMat, dtol::PetscReal) 
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
function MatSuperluSetILUDropTol(petsclib::PetscLibType, F::PetscMat, dtol::PetscReal) end

@for_petsc function MatSuperluSetILUDropTol(petsclib::$UnionPetscLib, F::PetscMat, dtol::$PetscReal )

    @chk ccall(
               (:MatSuperluSetILUDropTol, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscReal),
               F, dtol,
              )


	return nothing
end 

"""
	A::PetscMat = MatCreateSeqAIJViennaCL(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Vector{PetscInt}) 

# External Links
$(_doc_external("Mat/MatCreateSeqAIJViennaCL"))
"""
function MatCreateSeqAIJViennaCL(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Union{Ptr,Vector{PetscInt}}) end

@for_petsc function MatCreateSeqAIJViennaCL(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, nz::$PetscInt, nnz::Union{Ptr,Vector{$PetscInt}} )
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
	A::PetscMat = MatCreateSeqAIJSELL(petsclib::PetscLibType,comm::MPI_Comm, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Vector{PetscInt}) 
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
function MatCreateSeqAIJSELL(petsclib::PetscLibType, comm::MPI_Comm, m::PetscInt, n::PetscInt, nz::PetscInt, nnz::Union{Ptr,Vector{PetscInt}}) end

@for_petsc function MatCreateSeqAIJSELL(petsclib::$UnionPetscLib, comm::MPI_Comm, m::$PetscInt, n::$PetscInt, nz::$PetscInt, nnz::Union{Ptr,Vector{$PetscInt}} )
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
	MatGetColumnVector(petsclib::PetscLibType,A::PetscMat, yy::PetscVec, col::PetscInt) 
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
function MatGetColumnVector(petsclib::PetscLibType, A::PetscMat, yy::PetscVec, col::PetscInt) end

@for_petsc function MatGetColumnVector(petsclib::$UnionPetscLib, A::PetscMat, yy::PetscVec, col::$PetscInt )

    @chk ccall(
               (:MatGetColumnVector, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, $PetscInt),
               A, yy, col,
              )


	return nothing
end 

"""
	MatGetColumnNorms(petsclib::PetscLibType,A::PetscMat, type::NormType, norms::Vector{PetscReal}) 
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
function MatGetColumnNorms(petsclib::PetscLibType, A::PetscMat, type::NormType, norms::Vector{PetscReal}) end

@for_petsc function MatGetColumnNorms(petsclib::$UnionPetscLib, A::PetscMat, type::NormType, norms::Vector{$PetscReal} )

    @chk ccall(
               (:MatGetColumnNorms, $petsc_library),
               PetscErrorCode,
               (CMat, NormType, Ptr{$PetscReal}),
               A, type, norms,
              )


	return nothing
end 

"""
	MatGetColumnSumsRealPart(petsclib::PetscLibType,A::PetscMat, sums::Vector{PetscReal}) 
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
function MatGetColumnSumsRealPart(petsclib::PetscLibType, A::PetscMat, sums::Vector{PetscReal}) end

@for_petsc function MatGetColumnSumsRealPart(petsclib::$UnionPetscLib, A::PetscMat, sums::Vector{$PetscReal} )

    @chk ccall(
               (:MatGetColumnSumsRealPart, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscReal}),
               A, sums,
              )


	return nothing
end 

"""
	MatGetColumnSumsImaginaryPart(petsclib::PetscLibType,A::PetscMat, sums::Vector{PetscReal}) 
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
function MatGetColumnSumsImaginaryPart(petsclib::PetscLibType, A::PetscMat, sums::Vector{PetscReal}) end

@for_petsc function MatGetColumnSumsImaginaryPart(petsclib::$UnionPetscLib, A::PetscMat, sums::Vector{$PetscReal} )

    @chk ccall(
               (:MatGetColumnSumsImaginaryPart, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscReal}),
               A, sums,
              )


	return nothing
end 

"""
	MatGetColumnSums(petsclib::PetscLibType,A::PetscMat, sums::Vector{PetscScalar}) 
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
function MatGetColumnSums(petsclib::PetscLibType, A::PetscMat, sums::Vector{PetscScalar}) end

@for_petsc function MatGetColumnSums(petsclib::$UnionPetscLib, A::PetscMat, sums::Vector{$PetscScalar} )

    @chk ccall(
               (:MatGetColumnSums, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscScalar}),
               A, sums,
              )


	return nothing
end 

"""
	MatGetColumnMeansRealPart(petsclib::PetscLibType,A::PetscMat, means::Vector{PetscReal}) 
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
function MatGetColumnMeansRealPart(petsclib::PetscLibType, A::PetscMat, means::Vector{PetscReal}) end

@for_petsc function MatGetColumnMeansRealPart(petsclib::$UnionPetscLib, A::PetscMat, means::Vector{$PetscReal} )

    @chk ccall(
               (:MatGetColumnMeansRealPart, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscReal}),
               A, means,
              )


	return nothing
end 

"""
	MatGetColumnMeansImaginaryPart(petsclib::PetscLibType,A::PetscMat, means::Vector{PetscReal}) 
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
function MatGetColumnMeansImaginaryPart(petsclib::PetscLibType, A::PetscMat, means::Vector{PetscReal}) end

@for_petsc function MatGetColumnMeansImaginaryPart(petsclib::$UnionPetscLib, A::PetscMat, means::Vector{$PetscReal} )

    @chk ccall(
               (:MatGetColumnMeansImaginaryPart, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscReal}),
               A, means,
              )


	return nothing
end 

"""
	MatGetColumnMeans(petsclib::PetscLibType,A::PetscMat, means::Vector{PetscScalar}) 
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
function MatGetColumnMeans(petsclib::PetscLibType, A::PetscMat, means::Vector{PetscScalar}) end

@for_petsc function MatGetColumnMeans(petsclib::$UnionPetscLib, A::PetscMat, means::Vector{$PetscScalar} )

    @chk ccall(
               (:MatGetColumnMeans, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscScalar}),
               A, means,
              )


	return nothing
end 

"""
	MatGetColumnReductions(petsclib::PetscLibType,A::PetscMat, type::PetscInt, reductions::Vector{PetscReal}) 
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
function MatGetColumnReductions(petsclib::PetscLibType, A::PetscMat, type::PetscInt, reductions::Vector{PetscReal}) end

@for_petsc function MatGetColumnReductions(petsclib::$UnionPetscLib, A::PetscMat, type::$PetscInt, reductions::Vector{$PetscReal} )

    @chk ccall(
               (:MatGetColumnReductions, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscReal}),
               A, type, reductions,
              )


	return nothing
end 

"""
	MatReorderForNonzeroDiagonal(petsclib::PetscLibType,mat::PetscMat, abstol::PetscReal, ris::IS, cis::IS) 
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
function MatReorderForNonzeroDiagonal(petsclib::PetscLibType, mat::PetscMat, abstol::PetscReal, ris::IS, cis::IS) end

@for_petsc function MatReorderForNonzeroDiagonal(petsclib::$UnionPetscLib, mat::PetscMat, abstol::$PetscReal, ris::IS, cis::IS )

    @chk ccall(
               (:MatReorderForNonzeroDiagonal, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscReal, CIS, CIS),
               mat, abstol, ris, cis,
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
	A::PetscMat = MatCreateFromOptions(petsclib::PetscLibType,comm::MPI_Comm, prefix::Vector{Cchar}, bs::PetscInt, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt) 
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
function MatCreateFromOptions(petsclib::PetscLibType, comm::MPI_Comm, prefix::Vector{Cchar}, bs::PetscInt, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt) end

@for_petsc function MatCreateFromOptions(petsclib::$UnionPetscLib, comm::MPI_Comm, prefix::Vector{Cchar}, bs::$PetscInt, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt )
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
	MatSetErrorIfFailure(petsclib::PetscLibType,mat::PetscMat, flg::PetscBool) 
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
function MatSetErrorIfFailure(petsclib::PetscLibType, mat::PetscMat, flg::PetscBool) end

@for_petsc function MatSetErrorIfFailure(petsclib::$UnionPetscLib, mat::PetscMat, flg::PetscBool )

    @chk ccall(
               (:MatSetErrorIfFailure, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               mat, flg,
              )


	return nothing
end 

"""
	MatSetSizes(petsclib::PetscLibType,A::PetscMat, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt) 
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
function MatSetSizes(petsclib::PetscLibType, A::PetscMat, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt) end

@for_petsc function MatSetSizes(petsclib::$UnionPetscLib, A::PetscMat, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt )

    @chk ccall(
               (:MatSetSizes, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, $PetscInt, $PetscInt),
               A, m, n, M, N,
              )


	return nothing
end 

"""
	MatSetFromOptions(petsclib::PetscLibType,B::PetscMat) 
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
function MatSetFromOptions(petsclib::PetscLibType, B::PetscMat) end

@for_petsc function MatSetFromOptions(petsclib::$UnionPetscLib, B::PetscMat )

    @chk ccall(
               (:MatSetFromOptions, $petsc_library),
               PetscErrorCode,
               (CMat,),
               B,
              )


	return nothing
end 

"""
	MatXAIJSetPreallocation(petsclib::PetscLibType,A::PetscMat, bs::PetscInt, dnnz::Vector{PetscInt}, onnz::Vector{PetscInt}, dnnzu::Vector{PetscInt}, onnzu::Vector{PetscInt}) 
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
function MatXAIJSetPreallocation(petsclib::PetscLibType, A::PetscMat, bs::PetscInt, dnnz::Vector{PetscInt}, onnz::Vector{PetscInt}, dnnzu::Vector{PetscInt}, onnzu::Vector{PetscInt}) end

@for_petsc function MatXAIJSetPreallocation(petsclib::$UnionPetscLib, A::PetscMat, bs::$PetscInt, dnnz::Vector{$PetscInt}, onnz::Vector{$PetscInt}, dnnzu::Vector{$PetscInt}, onnzu::Vector{$PetscInt} )

    @chk ccall(
               (:MatXAIJSetPreallocation, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
               A, bs, dnnz, onnz, dnnzu, onnzu,
              )


	return nothing
end 

"""
	MatHeaderMerge(petsclib::PetscLibType,A::PetscMat, C::PetscMat) 
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
function MatHeaderMerge(petsclib::PetscLibType, A::PetscMat, C::PetscMat) end

@for_petsc function MatHeaderMerge(petsclib::$UnionPetscLib, A::PetscMat, C::PetscMat )
	C_ = Ref(C.ptr)

    @chk ccall(
               (:MatHeaderMerge, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, C_,
              )

	C.ptr = C_NULL

	return nothing
end 

"""
	MatHeaderReplace(petsclib::PetscLibType,A::PetscMat, C::PetscMat) 
Replaces the internal data of matrix `A` by the internal data of matrix `C` while deleting the outer wrapper of `C`

Input Parameters:
- `A` - a `Mat` whose internal data is to be replaced
- `C` - the `Mat` providing new internal data for `A`

Level: advanced

-seealso: `Mat`, `MatHeaderMerge()`

# External Links
$(_doc_external("Mat/MatHeaderReplace"))
"""
function MatHeaderReplace(petsclib::PetscLibType, A::PetscMat, C::PetscMat) end

@for_petsc function MatHeaderReplace(petsclib::$UnionPetscLib, A::PetscMat, C::PetscMat )
	C_ = Ref(C.ptr)

    @chk ccall(
               (:MatHeaderReplace, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, C_,
              )

	C.ptr = C_NULL

	return nothing
end 

"""
	MatBindToCPU(petsclib::PetscLibType,A::PetscMat, flg::PetscBool) 
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
function MatBindToCPU(petsclib::PetscLibType, A::PetscMat, flg::PetscBool) end

@for_petsc function MatBindToCPU(petsclib::$UnionPetscLib, A::PetscMat, flg::PetscBool )

    @chk ccall(
               (:MatBindToCPU, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               A, flg,
              )


	return nothing
end 

"""
	flg::PetscBool = MatBoundToCPU(petsclib::PetscLibType,A::PetscMat) 
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
function MatBoundToCPU(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatBoundToCPU(petsclib::$UnionPetscLib, A::PetscMat )
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
	MatSetPreallocationCOO(petsclib::PetscLibType,A::PetscMat, ncoo::PetscCount, coo_i::Vector{PetscInt}, coo_j::Vector{PetscInt}) 
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
function MatSetPreallocationCOO(petsclib::PetscLibType, A::PetscMat, ncoo::PetscCount, coo_i::Vector{PetscInt}, coo_j::Vector{PetscInt}) end

@for_petsc function MatSetPreallocationCOO(petsclib::$UnionPetscLib, A::PetscMat, ncoo::PetscCount, coo_i::Vector{$PetscInt}, coo_j::Vector{$PetscInt} )

    @chk ccall(
               (:MatSetPreallocationCOO, $petsc_library),
               PetscErrorCode,
               (CMat, PetscCount, Ptr{$PetscInt}, Ptr{$PetscInt}),
               A, ncoo, coo_i, coo_j,
              )


	return nothing
end 

"""
	MatSetPreallocationCOOLocal(petsclib::PetscLibType,A::PetscMat, ncoo::PetscCount, coo_i::Vector{PetscInt}, coo_j::Vector{PetscInt}) 
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
function MatSetPreallocationCOOLocal(petsclib::PetscLibType, A::PetscMat, ncoo::PetscCount, coo_i::Vector{PetscInt}, coo_j::Vector{PetscInt}) end

@for_petsc function MatSetPreallocationCOOLocal(petsclib::$UnionPetscLib, A::PetscMat, ncoo::PetscCount, coo_i::Vector{$PetscInt}, coo_j::Vector{$PetscInt} )

    @chk ccall(
               (:MatSetPreallocationCOOLocal, $petsc_library),
               PetscErrorCode,
               (CMat, PetscCount, Ptr{$PetscInt}, Ptr{$PetscInt}),
               A, ncoo, coo_i, coo_j,
              )


	return nothing
end 

"""
	MatSetValuesCOO(petsclib::PetscLibType,A::PetscMat, coo_v::Vector{PetscScalar}, imode::InsertMode) 
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
function MatSetValuesCOO(petsclib::PetscLibType, A::PetscMat, coo_v::Vector{PetscScalar}, imode::InsertMode) end

@for_petsc function MatSetValuesCOO(petsclib::$UnionPetscLib, A::PetscMat, coo_v::Vector{$PetscScalar}, imode::InsertMode )

    @chk ccall(
               (:MatSetValuesCOO, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{$PetscScalar}, InsertMode),
               A, coo_v, imode,
              )


	return nothing
end 

"""
	MatSetBindingPropagates(petsclib::PetscLibType,A::PetscMat, flg::PetscBool) 
Sets whether the state of being bound to the CPU for a GPU matrix type propagates to child and some other associated objects

Input Parameters:
- `A`   - the matrix
- `flg` - flag indicating whether the boundtocpu flag should be propagated

Level: developer

-seealso: [](ch_matrices), `Mat`, `VecSetBindingPropagates()`, `MatGetBindingPropagates()`

# External Links
$(_doc_external("Mat/MatSetBindingPropagates"))
"""
function MatSetBindingPropagates(petsclib::PetscLibType, A::PetscMat, flg::PetscBool) end

@for_petsc function MatSetBindingPropagates(petsclib::$UnionPetscLib, A::PetscMat, flg::PetscBool )

    @chk ccall(
               (:MatSetBindingPropagates, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               A, flg,
              )


	return nothing
end 

"""
	flg::PetscBool = MatGetBindingPropagates(petsclib::PetscLibType,A::PetscMat) 
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
function MatGetBindingPropagates(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatGetBindingPropagates(petsclib::$UnionPetscLib, A::PetscMat )
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
	bw::PetscInt = MatComputeBandwidth(petsclib::PetscLibType,A::PetscMat, fraction::PetscReal) 
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
function MatComputeBandwidth(petsclib::PetscLibType, A::PetscMat, fraction::PetscReal) end

@for_petsc function MatComputeBandwidth(petsclib::$UnionPetscLib, A::PetscMat, fraction::$PetscReal )
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
	MatAXPY(petsclib::PetscLibType,Y::PetscMat, a::PetscScalar, X::PetscMat, str::MatStructure) 
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
function MatAXPY(petsclib::PetscLibType, Y::PetscMat, a::PetscScalar, X::PetscMat, str::MatStructure) end

@for_petsc function MatAXPY(petsclib::$UnionPetscLib, Y::PetscMat, a::$PetscScalar, X::PetscMat, str::MatStructure )

    @chk ccall(
               (:MatAXPY, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscScalar, CMat, MatStructure),
               Y, a, X, str,
              )


	return nothing
end 

"""
	MatShift(petsclib::PetscLibType,Y::PetscMat, a::PetscScalar) 
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
function MatShift(petsclib::PetscLibType, Y::PetscMat, a::PetscScalar) end

@for_petsc function MatShift(petsclib::$UnionPetscLib, Y::PetscMat, a::$PetscScalar )

    @chk ccall(
               (:MatShift, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscScalar),
               Y, a,
              )


	return nothing
end 

"""
	MatDiagonalSet(petsclib::PetscLibType,Y::PetscMat, D::PetscVec, is::InsertMode) 
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
function MatDiagonalSet(petsclib::PetscLibType, Y::PetscMat, D::PetscVec, is::InsertMode) end

@for_petsc function MatDiagonalSet(petsclib::$UnionPetscLib, Y::PetscMat, D::PetscVec, is::InsertMode )

    @chk ccall(
               (:MatDiagonalSet, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, InsertMode),
               Y, D, is,
              )


	return nothing
end 

"""
	MatAYPX(petsclib::PetscLibType,Y::PetscMat, a::PetscScalar, X::PetscMat, str::MatStructure) 
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
function MatAYPX(petsclib::PetscLibType, Y::PetscMat, a::PetscScalar, X::PetscMat, str::MatStructure) end

@for_petsc function MatAYPX(petsclib::$UnionPetscLib, Y::PetscMat, a::$PetscScalar, X::PetscMat, str::MatStructure )

    @chk ccall(
               (:MatAYPX, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscScalar, CMat, MatStructure),
               Y, a, X, str,
              )


	return nothing
end 

"""
	MatComputeOperator(petsclib::PetscLibType,inmat::PetscMat, mattype::MatType, mat::PetscMat) 
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
function MatComputeOperator(petsclib::PetscLibType, inmat::PetscMat, mattype::MatType, mat::PetscMat) end

@for_petsc function MatComputeOperator(petsclib::$UnionPetscLib, inmat::PetscMat, mattype::MatType, mat::PetscMat )
	mat_ = Ref(mat.ptr)

    @chk ccall(
               (:MatComputeOperator, $petsc_library),
               PetscErrorCode,
               (CMat, MatType, Ptr{CMat}),
               inmat, mattype, mat_,
              )

	mat.ptr = C_NULL

	return nothing
end 

"""
	MatComputeOperatorTranspose(petsclib::PetscLibType,inmat::PetscMat, mattype::MatType, mat::PetscMat) 
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
function MatComputeOperatorTranspose(petsclib::PetscLibType, inmat::PetscMat, mattype::MatType, mat::PetscMat) end

@for_petsc function MatComputeOperatorTranspose(petsclib::$UnionPetscLib, inmat::PetscMat, mattype::MatType, mat::PetscMat )
	mat_ = Ref(mat.ptr)

    @chk ccall(
               (:MatComputeOperatorTranspose, $petsc_library),
               PetscErrorCode,
               (CMat, MatType, Ptr{CMat}),
               inmat, mattype, mat_,
              )

	mat.ptr = C_NULL

	return nothing
end 

"""
	MatFilter(petsclib::PetscLibType,A::PetscMat, tol::PetscReal, compress::PetscBool, keep::PetscBool) 
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
function MatFilter(petsclib::PetscLibType, A::PetscMat, tol::PetscReal, compress::PetscBool, keep::PetscBool) end

@for_petsc function MatFilter(petsclib::$UnionPetscLib, A::PetscMat, tol::$PetscReal, compress::PetscBool, keep::PetscBool )

    @chk ccall(
               (:MatFilter, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscReal, PetscBool, PetscBool),
               A, tol, compress, keep,
              )


	return nothing
end 

"""
	flg::PetscBool = MatMultEqual(petsclib::PetscLibType,A::PetscMat, B::PetscMat, n::PetscInt) 
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
function MatMultEqual(petsclib::PetscLibType, A::PetscMat, B::PetscMat, n::PetscInt) end

@for_petsc function MatMultEqual(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat, n::$PetscInt )
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
	flg::PetscBool = MatMultAddEqual(petsclib::PetscLibType,A::PetscMat, B::PetscMat, n::PetscInt) 
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
function MatMultAddEqual(petsclib::PetscLibType, A::PetscMat, B::PetscMat, n::PetscInt) end

@for_petsc function MatMultAddEqual(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat, n::$PetscInt )
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
	flg::PetscBool = MatMultTransposeEqual(petsclib::PetscLibType,A::PetscMat, B::PetscMat, n::PetscInt) 
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
function MatMultTransposeEqual(petsclib::PetscLibType, A::PetscMat, B::PetscMat, n::PetscInt) end

@for_petsc function MatMultTransposeEqual(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat, n::$PetscInt )
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
	flg::PetscBool = MatMultTransposeAddEqual(petsclib::PetscLibType,A::PetscMat, B::PetscMat, n::PetscInt) 
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
function MatMultTransposeAddEqual(petsclib::PetscLibType, A::PetscMat, B::PetscMat, n::PetscInt) end

@for_petsc function MatMultTransposeAddEqual(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat, n::$PetscInt )
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
	flg::PetscBool = MatMultHermitianTransposeEqual(petsclib::PetscLibType,A::PetscMat, B::PetscMat, n::PetscInt) 
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
function MatMultHermitianTransposeEqual(petsclib::PetscLibType, A::PetscMat, B::PetscMat, n::PetscInt) end

@for_petsc function MatMultHermitianTransposeEqual(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat, n::$PetscInt )
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
	flg::PetscBool = MatMultHermitianTransposeAddEqual(petsclib::PetscLibType,A::PetscMat, B::PetscMat, n::PetscInt) 
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
function MatMultHermitianTransposeAddEqual(petsclib::PetscLibType, A::PetscMat, B::PetscMat, n::PetscInt) end

@for_petsc function MatMultHermitianTransposeAddEqual(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat, n::$PetscInt )
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
	flg::PetscBool = MatMatMultEqual(petsclib::PetscLibType,A::PetscMat, B::PetscMat, C::PetscMat, n::PetscInt) 
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
function MatMatMultEqual(petsclib::PetscLibType, A::PetscMat, B::PetscMat, C::PetscMat, n::PetscInt) end

@for_petsc function MatMatMultEqual(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat, C::PetscMat, n::$PetscInt )
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
	flg::PetscBool = MatTransposeMatMultEqual(petsclib::PetscLibType,A::PetscMat, B::PetscMat, C::PetscMat, n::PetscInt) 
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
function MatTransposeMatMultEqual(petsclib::PetscLibType, A::PetscMat, B::PetscMat, C::PetscMat, n::PetscInt) end

@for_petsc function MatTransposeMatMultEqual(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat, C::PetscMat, n::$PetscInt )
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
	flg::PetscBool = MatMatTransposeMultEqual(petsclib::PetscLibType,A::PetscMat, B::PetscMat, C::PetscMat, n::PetscInt) 
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
function MatMatTransposeMultEqual(petsclib::PetscLibType, A::PetscMat, B::PetscMat, C::PetscMat, n::PetscInt) end

@for_petsc function MatMatTransposeMultEqual(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat, C::PetscMat, n::$PetscInt )
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
	flg::PetscBool = MatPtAPMultEqual(petsclib::PetscLibType,A::PetscMat, B::PetscMat, C::PetscMat, n::PetscInt) 
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
function MatPtAPMultEqual(petsclib::PetscLibType, A::PetscMat, B::PetscMat, C::PetscMat, n::PetscInt) end

@for_petsc function MatPtAPMultEqual(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat, C::PetscMat, n::$PetscInt )
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
	flg::PetscBool = MatRARtMultEqual(petsclib::PetscLibType,A::PetscMat, B::PetscMat, C::PetscMat, n::PetscInt) 
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
function MatRARtMultEqual(petsclib::PetscLibType, A::PetscMat, B::PetscMat, C::PetscMat, n::PetscInt) end

@for_petsc function MatRARtMultEqual(petsclib::$UnionPetscLib, A::PetscMat, B::PetscMat, C::PetscMat, n::$PetscInt )
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
	flg::PetscBool = MatIsLinear(petsclib::PetscLibType,A::PetscMat, n::PetscInt) 
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
function MatIsLinear(petsclib::PetscLibType, A::PetscMat, n::PetscInt) end

@for_petsc function MatIsLinear(petsclib::$UnionPetscLib, A::PetscMat, n::$PetscInt )
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
	MatSetHPL(petsclib::PetscLibType,A::PetscMat, iseed::Cint) 
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
function MatSetHPL(petsclib::PetscLibType, A::PetscMat, iseed::Cint) end

@for_petsc function MatSetHPL(petsclib::$UnionPetscLib, A::PetscMat, iseed::Cint )

    @chk ccall(
               (:MatSetHPL, $petsc_library),
               PetscErrorCode,
               (CMat, Cint),
               A, iseed,
              )


	return nothing
end 

"""
	L::PetscMat = MatCreateLaplacian(petsclib::PetscLibType,A::PetscMat, tol::PetscReal, weighted::PetscBool) 
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
$(_doc_external("Mat/MatCreateLaplacian"))
"""
function MatCreateLaplacian(petsclib::PetscLibType, A::PetscMat, tol::PetscReal, weighted::PetscBool) end

@for_petsc function MatCreateLaplacian(petsclib::$UnionPetscLib, A::PetscMat, tol::$PetscReal, weighted::PetscBool )
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
	MatOrderingRegister(petsclib::PetscLibType,sname::Vector{Cchar}, fnc::external) 
Adds a new sparse matrix ordering to the matrix package.

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - name of ordering (for example `MATORDERINGND`)
- `function` - function pointer that creates the ordering

Level: developer

-seealso: `Mat`, `MatOrderingType`, `MatOrderingRegisterAll()`, `MatGetOrdering()`

# External Links
$(_doc_external("Mat/MatOrderingRegister"))
"""
function MatOrderingRegister(petsclib::PetscLibType, sname::Vector{Cchar}, fnc::external) end

@for_petsc function MatOrderingRegister(petsclib::$UnionPetscLib, sname::Vector{Cchar}, fnc::external )

    @chk ccall(
               (:MatOrderingRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	MatGetOrdering(petsclib::PetscLibType,mat::PetscMat, type::MatOrderingType, rperm::IS, cperm::IS) 
Gets a reordering for a matrix to reduce fill or to
improve numerical stability of LU factorization.

Collective

Input Parameters:
- `mat`  - the matrix
- `type` - type of reordering, one of the following
-seealso: `MatOrderingRegister()`, `PCFactorSetMatOrderingType()`, `MatColoring`, `MatColoringCreate()`, `MatOrderingType`, `Mat`

# External Links
$(_doc_external("Mat/MatGetOrdering"))
"""
function MatGetOrdering(petsclib::PetscLibType, mat::PetscMat, type::MatOrderingType, rperm::IS, cperm::IS) end

@for_petsc function MatGetOrdering(petsclib::$UnionPetscLib, mat::PetscMat, type::MatOrderingType, rperm::IS, cperm::IS )

    @chk ccall(
               (:MatGetOrdering, $petsc_library),
               PetscErrorCode,
               (CMat, MatOrderingType, Ptr{CIS}, Ptr{CIS}),
               mat, type, rperm, cperm,
              )


	return nothing
end 

"""
	MatGetOrderingList(petsclib::PetscLibType,list::PetscFunctionList) 

# External Links
$(_doc_external("Mat/MatGetOrderingList"))
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
	MatMeshToCellGraph(petsclib::PetscLibType,mesh::PetscMat, ncommonnodes::PetscInt, dual::PetscMat) 
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
$(_doc_external("Mat/MatMeshToCellGraph"))
"""
function MatMeshToCellGraph(petsclib::PetscLibType, mesh::PetscMat, ncommonnodes::PetscInt, dual::PetscMat) end

@for_petsc function MatMeshToCellGraph(petsclib::$UnionPetscLib, mesh::PetscMat, ncommonnodes::$PetscInt, dual::PetscMat )
	dual_ = Ref(dual.ptr)

    @chk ccall(
               (:MatMeshToCellGraph, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{CMat}),
               mesh, ncommonnodes, dual_,
              )

	dual.ptr = C_NULL

	return nothing
end 

"""
	MatLMVMSetMultAlgorithm(petsclib::PetscLibType,B::PetscMat, alg::MatLMVMMultAlgorithm) 
Set the algorithm used by a `MatLMVM` for products

Logically collective

Input Parameters:
- `B`   - a `MatLMVM` matrix
- `alg` - one of the algorithm classes (`MAT_LMVM_MULT_RECURSIVE`, `MAT_LMVM_MULT_DENSE`, `MAT_LMVM_MULT_COMPACT_DENSE`)

Level: advanced

-seealso: [](ch_matrices), `MatLMVM`, `MatLMVMMultAlgorithm`, `MatLMVMGetMultAlgorithm()`

# External Links
$(_doc_external("Ksp/MatLMVMSetMultAlgorithm"))
"""
function MatLMVMSetMultAlgorithm(petsclib::PetscLibType, B::PetscMat, alg::MatLMVMMultAlgorithm) end

@for_petsc function MatLMVMSetMultAlgorithm(petsclib::$UnionPetscLib, B::PetscMat, alg::MatLMVMMultAlgorithm )

    @chk ccall(
               (:MatLMVMSetMultAlgorithm, $petsc_library),
               PetscErrorCode,
               (CMat, MatLMVMMultAlgorithm),
               B, alg,
              )


	return nothing
end 

"""
	MatLMVMGetMultAlgorithm(petsclib::PetscLibType,B::PetscMat, alg::MatLMVMMultAlgorithm) 
Get the algorithm used by a `MatLMVM` for products

Not collective

Input Parameter:
- `B` - a `MatLMVM` matrix

Output Parameter:
- `alg` - one of the algorithm classes (`MAT_LMVM_MULT_RECURSIVE`, `MAT_LMVM_MULT_DENSE`, `MAT_LMVM_MULT_COMPACT_DENSE`)

Level: advanced

-seealso: [](ch_matrices), `MatLMVM`, `MatLMVMMultAlgorithm`, `MatLMVMSetMultAlgorithm()`

# External Links
$(_doc_external("Ksp/MatLMVMGetMultAlgorithm"))
"""
function MatLMVMGetMultAlgorithm(petsclib::PetscLibType, B::PetscMat, alg::MatLMVMMultAlgorithm) end

@for_petsc function MatLMVMGetMultAlgorithm(petsclib::$UnionPetscLib, B::PetscMat, alg::MatLMVMMultAlgorithm )

    @chk ccall(
               (:MatLMVMGetMultAlgorithm, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatLMVMMultAlgorithm}),
               B, alg,
              )


	return nothing
end 

"""
	MatLMVMGetLastUpdate(petsclib::PetscLibType,B::PetscMat, x_prev::PetscVec, f_prev::PetscVec) 
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
$(_doc_external("Ksp/MatLMVMGetLastUpdate"))
"""
function MatLMVMGetLastUpdate(petsclib::PetscLibType, B::PetscMat, x_prev::PetscVec, f_prev::PetscVec) end

@for_petsc function MatLMVMGetLastUpdate(petsclib::$UnionPetscLib, B::PetscMat, x_prev::PetscVec, f_prev::PetscVec )
	x_prev_ = Ref(x_prev.ptr)
	f_prev_ = Ref(f_prev.ptr)

    @chk ccall(
               (:MatLMVMGetLastUpdate, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CVec}, Ptr{CVec}),
               B, x_prev_, f_prev_,
              )

	x_prev.ptr = C_NULL
	f_prev.ptr = C_NULL

	return nothing
end 

"""
	MatLMVMUpdate(petsclib::PetscLibType,B::PetscMat, X::PetscVec, F::PetscVec) 
Adds (X

Input Parameters:
- `B` - A `MATLMVM` matrix
- `X` - Solution vector
- `F` - Function vector

Level: intermediate

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMReset()`, `MatLMVMAllocate()`

# External Links
$(_doc_external("Ksp/MatLMVMUpdate"))
"""
function MatLMVMUpdate(petsclib::PetscLibType, B::PetscMat, X::PetscVec, F::PetscVec) end

@for_petsc function MatLMVMUpdate(petsclib::$UnionPetscLib, B::PetscMat, X::PetscVec, F::PetscVec )

    @chk ccall(
               (:MatLMVMUpdate, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               B, X, F,
              )


	return nothing
end 

"""
	MatLMVMClearJ0(petsclib::PetscLibType,B::PetscMat) 
Removes all definitions of J0 and reverts to
an identity matrix (scale = 1.0).

Input Parameter:
- `B` - A `MATLMVM` matrix

Level: advanced

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMSetJ0()`

# External Links
$(_doc_external("Ksp/MatLMVMClearJ0"))
"""
function MatLMVMClearJ0(petsclib::PetscLibType, B::PetscMat) end

@for_petsc function MatLMVMClearJ0(petsclib::$UnionPetscLib, B::PetscMat )

    @chk ccall(
               (:MatLMVMClearJ0, $petsc_library),
               PetscErrorCode,
               (CMat,),
               B,
              )


	return nothing
end 

"""
	MatLMVMSetJ0Scale(petsclib::PetscLibType,B::PetscMat, scale::PetscReal) 
Allows the user to define a scalar value
mu such that J0 = mu*I.

Input Parameters:
- `B`     - A `MATLMVM` matrix
- `scale` - Scalar value mu that defines the initial Jacobian

Level: advanced

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMSetDiagScale()`, `MatLMVMSetJ0()`

# External Links
$(_doc_external("Ksp/MatLMVMSetJ0Scale"))
"""
function MatLMVMSetJ0Scale(petsclib::PetscLibType, B::PetscMat, scale::PetscReal) end

@for_petsc function MatLMVMSetJ0Scale(petsclib::$UnionPetscLib, B::PetscMat, scale::$PetscReal )

    @chk ccall(
               (:MatLMVMSetJ0Scale, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscReal),
               B, scale,
              )


	return nothing
end 

"""
	MatLMVMSetJ0Diag(petsclib::PetscLibType,B::PetscMat, V::PetscVec) 
Allows the user to define a vector
V such that J0 = diag(V).

Input Parameters:
- `B` - An LMVM-type matrix
- `V` - Vector that defines the diagonal of the initial Jacobian: values are copied, V is not referenced

Level: advanced

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMSetScale()`, `MatLMVMSetJ0()`

# External Links
$(_doc_external("Ksp/MatLMVMSetJ0Diag"))
"""
function MatLMVMSetJ0Diag(petsclib::PetscLibType, B::PetscMat, V::PetscVec) end

@for_petsc function MatLMVMSetJ0Diag(petsclib::$UnionPetscLib, B::PetscMat, V::PetscVec )

    @chk ccall(
               (:MatLMVMSetJ0Diag, $petsc_library),
               PetscErrorCode,
               (CMat, CVec),
               B, V,
              )


	return nothing
end 

"""
	MatLMVMSetJ0(petsclib::PetscLibType,B::PetscMat, J0::PetscMat) 
Allows the user to define the initial Jacobian matrix from which the LMVM
up.

Input Parameters:
- `B`  - An LMVM-type matrix
- `J0` - The initial Jacobian matrix, will be referenced by B.

Level: advanced

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMSetJ0PC()`, `MatLMVMSetJ0KSP()`

# External Links
$(_doc_external("Ksp/MatLMVMSetJ0"))
"""
function MatLMVMSetJ0(petsclib::PetscLibType, B::PetscMat, J0::PetscMat) end

@for_petsc function MatLMVMSetJ0(petsclib::$UnionPetscLib, B::PetscMat, J0::PetscMat )

    @chk ccall(
               (:MatLMVMSetJ0, $petsc_library),
               PetscErrorCode,
               (CMat, CMat),
               B, J0,
              )


	return nothing
end 

"""
	MatLMVMSetJ0PC(petsclib::PetscLibType,B::PetscMat, J0pc::PC) 
Allows the user to define a `PC` object that acts as the initial inverse

Input Parameters:
- `B`    - A `MATLMVM` matrix
- `J0pc` - `PC` object where `PCApply()` defines an inverse application for J0

Level: advanced

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMGetJ0PC()`

# External Links
$(_doc_external("Ksp/MatLMVMSetJ0PC"))
"""
function MatLMVMSetJ0PC(petsclib::PetscLibType, B::PetscMat, J0pc::PC) end

@for_petsc function MatLMVMSetJ0PC(petsclib::$UnionPetscLib, B::PetscMat, J0pc::PC )

    @chk ccall(
               (:MatLMVMSetJ0PC, $petsc_library),
               PetscErrorCode,
               (CMat, PC),
               B, J0pc,
              )


	return nothing
end 

"""
	MatLMVMSetJ0KSP(petsclib::PetscLibType,B::PetscMat, J0ksp::PetscKSP) 
Allows the user to provide a pre
approximation.

Input Parameters:
- `B`     - A `MATLMVM` matrix
- `J0ksp` - `KSP` solver for the initial inverse-Jacobian application

Level: advanced

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMGetJ0KSP()`

# External Links
$(_doc_external("Ksp/MatLMVMSetJ0KSP"))
"""
function MatLMVMSetJ0KSP(petsclib::PetscLibType, B::PetscMat, J0ksp::PetscKSP) end

@for_petsc function MatLMVMSetJ0KSP(petsclib::$UnionPetscLib, B::PetscMat, J0ksp::PetscKSP )

    @chk ccall(
               (:MatLMVMSetJ0KSP, $petsc_library),
               PetscErrorCode,
               (CMat, CKSP),
               B, J0ksp,
              )


	return nothing
end 

"""
	MatLMVMGetJ0(petsclib::PetscLibType,B::PetscMat, J0::PetscMat) 
Returns a pointer to the internal `J0` matrix.

Input Parameter:
- `B` - A `MATLMVM` matrix

Output Parameter:
- `J0` - `Mat` object for defining the initial Jacobian

Level: advanced

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMSetJ0()`

# External Links
$(_doc_external("Ksp/MatLMVMGetJ0"))
"""
function MatLMVMGetJ0(petsclib::PetscLibType, B::PetscMat, J0::PetscMat) end

@for_petsc function MatLMVMGetJ0(petsclib::$UnionPetscLib, B::PetscMat, J0::PetscMat )
	J0_ = Ref(J0.ptr)

    @chk ccall(
               (:MatLMVMGetJ0, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               B, J0_,
              )

	J0.ptr = C_NULL

	return nothing
end 

"""
	MatLMVMGetJ0PC(petsclib::PetscLibType,B::PetscMat, J0pc::PC) 
Returns a pointer to the internal `PC` object
associated with the initial Jacobian.

Input Parameter:
- `B` - A `MATLMVM` matrix

Output Parameter:
- `J0pc` - `PC` object for defining the initial inverse-Jacobian

Level: advanced

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMSetJ0PC()`

# External Links
$(_doc_external("Ksp/MatLMVMGetJ0PC"))
"""
function MatLMVMGetJ0PC(petsclib::PetscLibType, B::PetscMat, J0pc::PC) end

@for_petsc function MatLMVMGetJ0PC(petsclib::$UnionPetscLib, B::PetscMat, J0pc::PC )

    @chk ccall(
               (:MatLMVMGetJ0PC, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{PC}),
               B, J0pc,
              )


	return nothing
end 

"""
	MatLMVMGetJ0KSP(petsclib::PetscLibType,B::PetscMat, J0ksp::PetscKSP) 
Returns a pointer to the internal `KSP` solver
associated with the initial Jacobian.

Input Parameter:
- `B` - A `MATLMVM` matrix

Output Parameter:
- `J0ksp` - `KSP` solver for defining the initial inverse-Jacobian

Level: advanced

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMSetJ0KSP()`

# External Links
$(_doc_external("Ksp/MatLMVMGetJ0KSP"))
"""
function MatLMVMGetJ0KSP(petsclib::PetscLibType, B::PetscMat, J0ksp::PetscKSP) end

@for_petsc function MatLMVMGetJ0KSP(petsclib::$UnionPetscLib, B::PetscMat, J0ksp::PetscKSP )
	J0ksp_ = Ref(J0ksp.ptr)

    @chk ccall(
               (:MatLMVMGetJ0KSP, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CKSP}),
               B, J0ksp_,
              )

	J0ksp.ptr = C_NULL

	return nothing
end 

"""
	MatLMVMApplyJ0Fwd(petsclib::PetscLibType,B::PetscMat, X::PetscVec, Y::PetscVec) 
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
$(_doc_external("Ksp/MatLMVMApplyJ0Fwd"))
"""
function MatLMVMApplyJ0Fwd(petsclib::PetscLibType, B::PetscMat, X::PetscVec, Y::PetscVec) end

@for_petsc function MatLMVMApplyJ0Fwd(petsclib::$UnionPetscLib, B::PetscMat, X::PetscVec, Y::PetscVec )

    @chk ccall(
               (:MatLMVMApplyJ0Fwd, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               B, X, Y,
              )


	return nothing
end 

"""
	MatLMVMApplyJ0Inv(petsclib::PetscLibType,B::PetscMat, X::PetscVec, Y::PetscVec) 
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
$(_doc_external("Ksp/MatLMVMApplyJ0Inv"))
"""
function MatLMVMApplyJ0Inv(petsclib::PetscLibType, B::PetscMat, X::PetscVec, Y::PetscVec) end

@for_petsc function MatLMVMApplyJ0Inv(petsclib::$UnionPetscLib, B::PetscMat, X::PetscVec, Y::PetscVec )

    @chk ccall(
               (:MatLMVMApplyJ0Inv, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               B, X, Y,
              )


	return nothing
end 

"""
	flg::PetscBool = MatLMVMIsAllocated(petsclib::PetscLibType,B::PetscMat) 
Returns a boolean flag that shows whether
the necessary data structures for the underlying matrix is allocated.

Input Parameter:
- `B` - A `MATLMVM` matrix

Output Parameter:
- `flg` - `PETSC_TRUE` if allocated, `PETSC_FALSE` otherwise

Level: intermediate

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMAllocate()`, `MatLMVMReset()`

# External Links
$(_doc_external("Ksp/MatLMVMIsAllocated"))
"""
function MatLMVMIsAllocated(petsclib::PetscLibType, B::PetscMat) end

@for_petsc function MatLMVMIsAllocated(petsclib::$UnionPetscLib, B::PetscMat )
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
	MatLMVMAllocate(petsclib::PetscLibType,B::PetscMat, X::PetscVec, F::PetscVec) 
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
$(_doc_external("Ksp/MatLMVMAllocate"))
"""
function MatLMVMAllocate(petsclib::PetscLibType, B::PetscMat, X::PetscVec, F::PetscVec) end

@for_petsc function MatLMVMAllocate(petsclib::$UnionPetscLib, B::PetscMat, X::PetscVec, F::PetscVec )

    @chk ccall(
               (:MatLMVMAllocate, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               B, X, F,
              )


	return nothing
end 

"""
	MatLMVMResetShift(petsclib::PetscLibType,B::PetscMat) 
Zero the shift factor for a `MATLMVM`.

Input Parameter:
- `B` - A `MATLMVM` matrix

Level: intermediate

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMAllocate()`, `MatLMVMUpdate()`

# External Links
$(_doc_external("Ksp/MatLMVMResetShift"))
"""
function MatLMVMResetShift(petsclib::PetscLibType, B::PetscMat) end

@for_petsc function MatLMVMResetShift(petsclib::$UnionPetscLib, B::PetscMat )

    @chk ccall(
               (:MatLMVMResetShift, $petsc_library),
               PetscErrorCode,
               (CMat,),
               B,
              )


	return nothing
end 

"""
	MatLMVMReset(petsclib::PetscLibType,B::PetscMat, destructive::PetscBool) 
Flushes all of the accumulated updates out of
the `MATLMVM` approximation.

Input Parameters:
- `B`           - A `MATLMVM` matrix
- `destructive` - flag for enabling destruction of data structures

Level: intermediate

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMAllocate()`, `MatLMVMUpdate()`

# External Links
$(_doc_external("Ksp/MatLMVMReset"))
"""
function MatLMVMReset(petsclib::PetscLibType, B::PetscMat, destructive::PetscBool) end

@for_petsc function MatLMVMReset(petsclib::$UnionPetscLib, B::PetscMat, destructive::PetscBool )

    @chk ccall(
               (:MatLMVMReset, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               B, destructive,
              )


	return nothing
end 

"""
	MatLMVMSetHistorySize(petsclib::PetscLibType,B::PetscMat, hist_size::PetscInt) 
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
$(_doc_external("Ksp/MatLMVMSetHistorySize"))
"""
function MatLMVMSetHistorySize(petsclib::PetscLibType, B::PetscMat, hist_size::PetscInt) end

@for_petsc function MatLMVMSetHistorySize(petsclib::$UnionPetscLib, B::PetscMat, hist_size::$PetscInt )

    @chk ccall(
               (:MatLMVMSetHistorySize, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt),
               B, hist_size,
              )


	return nothing
end 

"""
	hist_size::PetscInt = MatLMVMGetHistorySize(petsclib::PetscLibType,B::PetscMat) 

# External Links
$(_doc_external("Ksp/MatLMVMGetHistorySize"))
"""
function MatLMVMGetHistorySize(petsclib::PetscLibType, B::PetscMat) end

@for_petsc function MatLMVMGetHistorySize(petsclib::$UnionPetscLib, B::PetscMat )
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
	nupdates::PetscInt = MatLMVMGetUpdateCount(petsclib::PetscLibType,B::PetscMat) 
Returns the number of accepted updates.

Input Parameter:
- `B` - A `MATLMVM` matrix

Output Parameter:
- `nupdates` - number of accepted updates

Level: intermediate

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMGetRejectCount()`, `MatLMVMReset()`

# External Links
$(_doc_external("Ksp/MatLMVMGetUpdateCount"))
"""
function MatLMVMGetUpdateCount(petsclib::PetscLibType, B::PetscMat) end

@for_petsc function MatLMVMGetUpdateCount(petsclib::$UnionPetscLib, B::PetscMat )
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
	nrejects::PetscInt = MatLMVMGetRejectCount(petsclib::PetscLibType,B::PetscMat) 
Returns the number of rejected updates.
The counters are reset when `MatLMVMReset()` is called.

Input Parameter:
- `B` - A `MATLMVM` matrix

Output Parameter:
- `nrejects` - number of rejected updates

Level: intermediate

-seealso: [](ch_ksp), [LMVM Matrices](sec_matlmvm), `MATLMVM`, `MatLMVMReset()`

# External Links
$(_doc_external("Ksp/MatLMVMGetRejectCount"))
"""
function MatLMVMGetRejectCount(petsclib::PetscLibType, B::PetscMat) end

@for_petsc function MatLMVMGetRejectCount(petsclib::$UnionPetscLib, B::PetscMat )
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
$(_doc_external("Ksp/MatCreateLMVMDiagBroyden"))
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
$(_doc_external("Ksp/MatCreateLMVMBadBroyden"))
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
$(_doc_external("Ksp/MatCreateLMVMBroyden"))
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
$(_doc_external("Ksp/MatCreateLMVMDQN"))
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
$(_doc_external("Ksp/MatCreateLMVMDBFGS"))
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
$(_doc_external("Ksp/MatCreateLMVMDDFP"))
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
	MatLMVMDenseSetType(petsclib::PetscLibType,B::PetscMat, type::MatLMVMDenseType) 
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
$(_doc_external("Ksp/MatLMVMDenseSetType"))
"""
function MatLMVMDenseSetType(petsclib::PetscLibType, B::PetscMat, type::MatLMVMDenseType) end

@for_petsc function MatLMVMDenseSetType(petsclib::$UnionPetscLib, B::PetscMat, type::MatLMVMDenseType )

    @chk ccall(
               (:MatLMVMDenseSetType, $petsc_library),
               PetscErrorCode,
               (CMat, MatLMVMDenseType),
               B, type,
              )


	return nothing
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
$(_doc_external("Ksp/MatCreateLMVMSR1"))
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
$(_doc_external("Ksp/MatCreateLMVMBFGS"))
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
	phi::PetscReal = MatLMVMSymBroydenGetPhi(petsclib::PetscLibType,B::PetscMat) 
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
$(_doc_external("Ksp/MatLMVMSymBroydenGetPhi"))
"""
function MatLMVMSymBroydenGetPhi(petsclib::PetscLibType, B::PetscMat) end

@for_petsc function MatLMVMSymBroydenGetPhi(petsclib::$UnionPetscLib, B::PetscMat )
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
	MatLMVMSymBroydenSetPhi(petsclib::PetscLibType,B::PetscMat, phi::PetscReal) 
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
$(_doc_external("Ksp/MatLMVMSymBroydenSetPhi"))
"""
function MatLMVMSymBroydenSetPhi(petsclib::PetscLibType, B::PetscMat, phi::PetscReal) end

@for_petsc function MatLMVMSymBroydenSetPhi(petsclib::$UnionPetscLib, B::PetscMat, phi::$PetscReal )

    @chk ccall(
               (:MatLMVMSymBroydenSetPhi, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscReal),
               B, phi,
              )


	return nothing
end 

"""
	psi::PetscReal = MatLMVMSymBadBroydenGetPsi(petsclib::PetscLibType,B::PetscMat) 
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
$(_doc_external("Ksp/MatLMVMSymBadBroydenGetPsi"))
"""
function MatLMVMSymBadBroydenGetPsi(petsclib::PetscLibType, B::PetscMat) end

@for_petsc function MatLMVMSymBadBroydenGetPsi(petsclib::$UnionPetscLib, B::PetscMat )
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
	MatLMVMSymBadBroydenSetPsi(petsclib::PetscLibType,B::PetscMat, psi::PetscReal) 
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
$(_doc_external("Ksp/MatLMVMSymBadBroydenSetPsi"))
"""
function MatLMVMSymBadBroydenSetPsi(petsclib::PetscLibType, B::PetscMat, psi::PetscReal) end

@for_petsc function MatLMVMSymBadBroydenSetPsi(petsclib::$UnionPetscLib, B::PetscMat, psi::$PetscReal )

    @chk ccall(
               (:MatLMVMSymBadBroydenSetPsi, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscReal),
               B, psi,
              )


	return nothing
end 

"""
	MatLMVMSymBroydenSetDelta(petsclib::PetscLibType,B::PetscMat, delta::PetscScalar) 
Sets the starting value for the diagonal scaling vector computed
in the SymBrdn approximations (also works for BFGS and DFP).

Input Parameters:
- `B`     - `MATLMVM` matrix
- `delta` - initial value for diagonal scaling

Level: intermediate

-seealso: [](ch_ksp), `MATLMVMSYMBROYDEN`

# External Links
$(_doc_external("Ksp/MatLMVMSymBroydenSetDelta"))
"""
function MatLMVMSymBroydenSetDelta(petsclib::PetscLibType, B::PetscMat, delta::PetscScalar) end

@for_petsc function MatLMVMSymBroydenSetDelta(petsclib::$UnionPetscLib, B::PetscMat, delta::$PetscScalar )

    @chk ccall(
               (:MatLMVMSymBroydenSetDelta, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscScalar),
               B, delta,
              )


	return nothing
end 

"""
	MatLMVMSymBroydenSetScaleType(petsclib::PetscLibType,B::PetscMat, stype::MatLMVMSymBroydenScaleType) 
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
$(_doc_external("Ksp/MatLMVMSymBroydenSetScaleType"))
"""
function MatLMVMSymBroydenSetScaleType(petsclib::PetscLibType, B::PetscMat, stype::MatLMVMSymBroydenScaleType) end

@for_petsc function MatLMVMSymBroydenSetScaleType(petsclib::$UnionPetscLib, B::PetscMat, stype::MatLMVMSymBroydenScaleType )

    @chk ccall(
               (:MatLMVMSymBroydenSetScaleType, $petsc_library),
               PetscErrorCode,
               (CMat, MatLMVMSymBroydenScaleType),
               B, stype,
              )


	return nothing
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
$(_doc_external("Ksp/MatCreateLMVMSymBroyden"))
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
$(_doc_external("Ksp/MatCreateLMVMSymBadBroyden"))
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
$(_doc_external("Ksp/MatCreateLMVMDFP"))
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
	S::PetscMat = MatCreateSchurComplement(petsclib::PetscLibType,A00::PetscMat, Ap00::PetscMat, A01::PetscMat, A10::PetscMat, A11::PetscMat) 
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
$(_doc_external("Ksp/MatCreateSchurComplement"))
"""
function MatCreateSchurComplement(petsclib::PetscLibType, A00::PetscMat, Ap00::PetscMat, A01::PetscMat, A10::PetscMat, A11::PetscMat) end

@for_petsc function MatCreateSchurComplement(petsclib::$UnionPetscLib, A00::PetscMat, Ap00::PetscMat, A01::PetscMat, A10::PetscMat, A11::PetscMat )
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
	MatSchurComplementSetSubMatrices(petsclib::PetscLibType,S::PetscMat, A00::PetscMat, Ap00::PetscMat, A01::PetscMat, A10::PetscMat, A11::PetscMat) 
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
$(_doc_external("Ksp/MatSchurComplementSetSubMatrices"))
"""
function MatSchurComplementSetSubMatrices(petsclib::PetscLibType, S::PetscMat, A00::PetscMat, Ap00::PetscMat, A01::PetscMat, A10::PetscMat, A11::PetscMat) end

@for_petsc function MatSchurComplementSetSubMatrices(petsclib::$UnionPetscLib, S::PetscMat, A00::PetscMat, Ap00::PetscMat, A01::PetscMat, A10::PetscMat, A11::PetscMat )

    @chk ccall(
               (:MatSchurComplementSetSubMatrices, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat, CMat, CMat, CMat),
               S, A00, Ap00, A01, A10, A11,
              )


	return nothing
end 

"""
	MatSchurComplementGetKSP(petsclib::PetscLibType,S::PetscMat, ksp::PetscKSP) 
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
$(_doc_external("Ksp/MatSchurComplementGetKSP"))
"""
function MatSchurComplementGetKSP(petsclib::PetscLibType, S::PetscMat, ksp::PetscKSP) end

@for_petsc function MatSchurComplementGetKSP(petsclib::$UnionPetscLib, S::PetscMat, ksp::PetscKSP )
	ksp_ = Ref(ksp.ptr)

    @chk ccall(
               (:MatSchurComplementGetKSP, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CKSP}),
               S, ksp_,
              )

	ksp.ptr = C_NULL

	return nothing
end 

"""
	MatSchurComplementSetKSP(petsclib::PetscLibType,S::PetscMat, ksp::PetscKSP) 
Sets the `KSP` object that is used to solve with `A00` in the Schur complement matrix  S = A11

Not Collective

Input Parameters:
- `S`   - matrix created with `MatCreateSchurComplement()`
- `ksp` - the linear solver object

Level: developer

-seealso: [](ch_ksp), `Mat`, `MatSchurComplementGetKSP()`, `MatCreateSchurComplement()`, `MatCreateNormal()`, `MatMult()`, `MatCreate()`, `MATSCHURCOMPLEMENT`

# External Links
$(_doc_external("Ksp/MatSchurComplementSetKSP"))
"""
function MatSchurComplementSetKSP(petsclib::PetscLibType, S::PetscMat, ksp::PetscKSP) end

@for_petsc function MatSchurComplementSetKSP(petsclib::$UnionPetscLib, S::PetscMat, ksp::PetscKSP )

    @chk ccall(
               (:MatSchurComplementSetKSP, $petsc_library),
               PetscErrorCode,
               (CMat, CKSP),
               S, ksp,
              )


	return nothing
end 

"""
	MatSchurComplementUpdateSubMatrices(petsclib::PetscLibType,S::PetscMat, A00::PetscMat, Ap00::PetscMat, A01::PetscMat, A10::PetscMat, A11::PetscMat) 
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
$(_doc_external("Ksp/MatSchurComplementUpdateSubMatrices"))
"""
function MatSchurComplementUpdateSubMatrices(petsclib::PetscLibType, S::PetscMat, A00::PetscMat, Ap00::PetscMat, A01::PetscMat, A10::PetscMat, A11::PetscMat) end

@for_petsc function MatSchurComplementUpdateSubMatrices(petsclib::$UnionPetscLib, S::PetscMat, A00::PetscMat, Ap00::PetscMat, A01::PetscMat, A10::PetscMat, A11::PetscMat )

    @chk ccall(
               (:MatSchurComplementUpdateSubMatrices, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat, CMat, CMat, CMat),
               S, A00, Ap00, A01, A10, A11,
              )


	return nothing
end 

"""
	MatSchurComplementGetSubMatrices(petsclib::PetscLibType,S::PetscMat, A00::PetscMat, Ap00::PetscMat, A01::PetscMat, A10::PetscMat, A11::PetscMat) 
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
$(_doc_external("Ksp/MatSchurComplementGetSubMatrices"))
"""
function MatSchurComplementGetSubMatrices(petsclib::PetscLibType, S::PetscMat, A00::PetscMat, Ap00::PetscMat, A01::PetscMat, A10::PetscMat, A11::PetscMat) end

@for_petsc function MatSchurComplementGetSubMatrices(petsclib::$UnionPetscLib, S::PetscMat, A00::PetscMat, Ap00::PetscMat, A01::PetscMat, A10::PetscMat, A11::PetscMat )
	A00_ = Ref(A00.ptr)
	Ap00_ = Ref(Ap00.ptr)
	A01_ = Ref(A01.ptr)
	A10_ = Ref(A10.ptr)
	A11_ = Ref(A11.ptr)

    @chk ccall(
               (:MatSchurComplementGetSubMatrices, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}, Ptr{CMat}, Ptr{CMat}, Ptr{CMat}, Ptr{CMat}),
               S, A00_, Ap00_, A01_, A10_, A11_,
              )

	A00.ptr = C_NULL
	Ap00.ptr = C_NULL
	A01.ptr = C_NULL
	A10.ptr = C_NULL
	A11.ptr = C_NULL

	return nothing
end 

"""
	MatSchurComplementComputeExplicitOperator(petsclib::PetscLibType,A::PetscMat, S::PetscMat) 
Compute the Schur complement matrix explicitly

Collective

Input Parameter:
- `A` - the matrix obtained with `MatCreateSchurComplement()`

Output Parameter:
- `S` - the Schur complement matrix

Level: advanced

-seealso: [](ch_ksp), `MatCreateSchurComplement()`, `MatSchurComplementUpdateSubMatrices()`, `MatSchurComplementGetPmat()`

# External Links
$(_doc_external("Ksp/MatSchurComplementComputeExplicitOperator"))
"""
function MatSchurComplementComputeExplicitOperator(petsclib::PetscLibType, A::PetscMat, S::PetscMat) end

@for_petsc function MatSchurComplementComputeExplicitOperator(petsclib::$UnionPetscLib, A::PetscMat, S::PetscMat )
	S_ = Ref(S.ptr)

    @chk ccall(
               (:MatSchurComplementComputeExplicitOperator, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}),
               A, S_,
              )

	S.ptr = C_NULL

	return nothing
end 

"""
	MatGetSchurComplement(petsclib::PetscLibType,A::PetscMat, isrow0::IS, iscol0::IS, isrow1::IS, iscol1::IS, mreuse::MatReuse, S::PetscMat, ainvtype::MatSchurComplementAinvType, preuse::MatReuse, Sp::PetscMat) 
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
$(_doc_external("Ksp/MatGetSchurComplement"))
"""
function MatGetSchurComplement(petsclib::PetscLibType, A::PetscMat, isrow0::IS, iscol0::IS, isrow1::IS, iscol1::IS, mreuse::MatReuse, S::PetscMat, ainvtype::MatSchurComplementAinvType, preuse::MatReuse, Sp::PetscMat) end

@for_petsc function MatGetSchurComplement(petsclib::$UnionPetscLib, A::PetscMat, isrow0::IS, iscol0::IS, isrow1::IS, iscol1::IS, mreuse::MatReuse, S::PetscMat, ainvtype::MatSchurComplementAinvType, preuse::MatReuse, Sp::PetscMat )
	S_ = Ref(S.ptr)
	Sp_ = Ref(Sp.ptr)

    @chk ccall(
               (:MatGetSchurComplement, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, CIS, CIS, CIS, MatReuse, Ptr{CMat}, MatSchurComplementAinvType, MatReuse, Ptr{CMat}),
               A, isrow0, iscol0, isrow1, iscol1, mreuse, S_, ainvtype, preuse, Sp_,
              )

	S.ptr = C_NULL
	Sp.ptr = C_NULL

	return nothing
end 

"""
	MatSchurComplementSetAinvType(petsclib::PetscLibType,S::PetscMat, ainvtype::MatSchurComplementAinvType) 
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
$(_doc_external("Ksp/MatSchurComplementSetAinvType"))
"""
function MatSchurComplementSetAinvType(petsclib::PetscLibType, S::PetscMat, ainvtype::MatSchurComplementAinvType) end

@for_petsc function MatSchurComplementSetAinvType(petsclib::$UnionPetscLib, S::PetscMat, ainvtype::MatSchurComplementAinvType )

    @chk ccall(
               (:MatSchurComplementSetAinvType, $petsc_library),
               PetscErrorCode,
               (CMat, MatSchurComplementAinvType),
               S, ainvtype,
              )


	return nothing
end 

"""
	ainvtype::MatSchurComplementAinvType = MatSchurComplementGetAinvType(petsclib::PetscLibType,S::PetscMat) 
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
$(_doc_external("Ksp/MatSchurComplementGetAinvType"))
"""
function MatSchurComplementGetAinvType(petsclib::PetscLibType, S::PetscMat) end

@for_petsc function MatSchurComplementGetAinvType(petsclib::$UnionPetscLib, S::PetscMat )
	ainvtype_ = Ref{MatSchurComplementAinvType}()

    @chk ccall(
               (:MatSchurComplementGetAinvType, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{MatSchurComplementAinvType}),
               S, ainvtype_,
              )

	ainvtype = unsafe_string(ainvtype_[])

	return ainvtype
end 

"""
	Sp::PetscMat = MatCreateSchurComplementPmat(petsclib::PetscLibType,A00::PetscMat, A01::PetscMat, A10::PetscMat, A11::PetscMat, ainvtype::MatSchurComplementAinvType, preuse::MatReuse) 
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
$(_doc_external("Ksp/MatCreateSchurComplementPmat"))
"""
function MatCreateSchurComplementPmat(petsclib::PetscLibType, A00::PetscMat, A01::PetscMat, A10::PetscMat, A11::PetscMat, ainvtype::MatSchurComplementAinvType, preuse::MatReuse) end

@for_petsc function MatCreateSchurComplementPmat(petsclib::$UnionPetscLib, A00::PetscMat, A01::PetscMat, A10::PetscMat, A11::PetscMat, ainvtype::MatSchurComplementAinvType, preuse::MatReuse )
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
	MatSchurComplementGetPmat(petsclib::PetscLibType,S::PetscMat, preuse::MatReuse, Sp::PetscMat) 
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
$(_doc_external("Ksp/MatSchurComplementGetPmat"))
"""
function MatSchurComplementGetPmat(petsclib::PetscLibType, S::PetscMat, preuse::MatReuse, Sp::PetscMat) end

@for_petsc function MatSchurComplementGetPmat(petsclib::$UnionPetscLib, S::PetscMat, preuse::MatReuse, Sp::PetscMat )
	Sp_ = Ref(Sp.ptr)

    @chk ccall(
               (:MatSchurComplementGetPmat, $petsc_library),
               PetscErrorCode,
               (CMat, MatReuse, Ptr{CMat}),
               S, preuse, Sp_,
              )

	Sp.ptr = C_NULL

	return nothing
end 

"""
	J::PetscMat = MatCreateSNESMFMore(petsclib::PetscLibType,snes::PetscSNES, x::PetscVec) 
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
$(_doc_external("Snes/MatCreateSNESMFMore"))
"""
function MatCreateSNESMFMore(petsclib::PetscLibType, snes::PetscSNES, x::PetscVec) end

@for_petsc function MatCreateSNESMFMore(petsclib::$UnionPetscLib, snes::PetscSNES, x::PetscVec )
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
	MatSNESMFMoreSetParameters(petsclib::PetscLibType,mat::PetscMat, error::PetscReal, umin::PetscReal, h::PetscReal) 
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
$(_doc_external("Snes/MatSNESMFMoreSetParameters"))
"""
function MatSNESMFMoreSetParameters(petsclib::PetscLibType, mat::PetscMat, error::PetscReal, umin::PetscReal, h::PetscReal) end

@for_petsc function MatSNESMFMoreSetParameters(petsclib::$UnionPetscLib, mat::PetscMat, error::$PetscReal, umin::$PetscReal, h::$PetscReal )

    @chk ccall(
               (:MatSNESMFMoreSetParameters, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscReal, $PetscReal, $PetscReal),
               mat, error, umin, h,
              )


	return nothing
end 

"""
	MatSNESMFGetSNES(petsclib::PetscLibType,J::PetscMat, snes::PetscSNES) 
returns the `SNES` associated with a matrix created with `MatCreateSNESMF()`

Not Collective

Input Parameter:
- `J` - the matrix

Output Parameter:
- `snes` - the `SNES` object

Level: advanced

-seealso: [](ch_snes), `Mat`, `SNES`, `MatCreateSNESMF()`

# External Links
$(_doc_external("Snes/MatSNESMFGetSNES"))
"""
function MatSNESMFGetSNES(petsclib::PetscLibType, J::PetscMat, snes::PetscSNES) end

@for_petsc function MatSNESMFGetSNES(petsclib::$UnionPetscLib, J::PetscMat, snes::PetscSNES )
	snes_ = Ref(snes.ptr)

    @chk ccall(
               (:MatSNESMFGetSNES, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CSNES}),
               J, snes_,
              )

	snes.ptr = C_NULL

	return nothing
end 

"""
	MatSNESMFSetReuseBase(petsclib::PetscLibType,J::PetscMat, use::PetscBool) 
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
$(_doc_external("Snes/MatSNESMFSetReuseBase"))
"""
function MatSNESMFSetReuseBase(petsclib::PetscLibType, J::PetscMat, use::PetscBool) end

@for_petsc function MatSNESMFSetReuseBase(petsclib::$UnionPetscLib, J::PetscMat, use::PetscBool )

    @chk ccall(
               (:MatSNESMFSetReuseBase, $petsc_library),
               PetscErrorCode,
               (CMat, PetscBool),
               J, use,
              )


	return nothing
end 

"""
	use::PetscBool = MatSNESMFGetReuseBase(petsclib::PetscLibType,J::PetscMat) 
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
$(_doc_external("Snes/MatSNESMFGetReuseBase"))
"""
function MatSNESMFGetReuseBase(petsclib::PetscLibType, J::PetscMat) end

@for_petsc function MatSNESMFGetReuseBase(petsclib::$UnionPetscLib, J::PetscMat )
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
	J::PetscMat = MatCreateSNESMF(petsclib::PetscLibType,snes::PetscSNES) 
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
$(_doc_external("Snes/MatCreateSNESMF"))
"""
function MatCreateSNESMF(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function MatCreateSNESMF(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	dm::PetscDM = MatGetDM(petsclib::PetscLibType,A::PetscMat) 
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
function MatGetDM(petsclib::PetscLibType, A::PetscMat) end

@for_petsc function MatGetDM(petsclib::$UnionPetscLib, A::PetscMat )
	dm_ = Ref(CDM)()

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
	MatSetDM(petsclib::PetscLibType,A::PetscMat, dm::PetscDM) 
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
function MatSetDM(petsclib::PetscLibType, A::PetscMat, dm::PetscDM) end

@for_petsc function MatSetDM(petsclib::$UnionPetscLib, A::PetscMat, dm::PetscDM )

    @chk ccall(
               (:MatSetDM, $petsc_library),
               PetscErrorCode,
               (CMat, CDM),
               A, dm,
              )


	return nothing
end 

"""
	MatDFischer(petsclib::PetscLibType,jac::PetscMat, X::PetscVec, Con::PetscVec, XL::PetscVec, XU::PetscVec, T1::PetscVec, T2::PetscVec, Da::PetscVec, Db::PetscVec) 
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
function MatDFischer(petsclib::PetscLibType, jac::PetscMat, X::PetscVec, Con::PetscVec, XL::PetscVec, XU::PetscVec, T1::PetscVec, T2::PetscVec, Da::PetscVec, Db::PetscVec) end

@for_petsc function MatDFischer(petsclib::$UnionPetscLib, jac::PetscMat, X::PetscVec, Con::PetscVec, XL::PetscVec, XU::PetscVec, T1::PetscVec, T2::PetscVec, Da::PetscVec, Db::PetscVec )

    @chk ccall(
               (:MatDFischer, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec, CVec, CVec, CVec, CVec, CVec, CVec),
               jac, X, Con, XL, XU, T1, T2, Da, Db,
              )


	return nothing
end 

"""
	MatDSFischer(petsclib::PetscLibType,jac::PetscMat, X::PetscVec, Con::PetscVec, XL::PetscVec, XU::PetscVec, mu::PetscReal, T1::PetscVec, T2::PetscVec, Da::PetscVec, Db::PetscVec, Dm::PetscVec) 
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
function MatDSFischer(petsclib::PetscLibType, jac::PetscMat, X::PetscVec, Con::PetscVec, XL::PetscVec, XU::PetscVec, mu::PetscReal, T1::PetscVec, T2::PetscVec, Da::PetscVec, Db::PetscVec, Dm::PetscVec) end

@for_petsc function MatDSFischer(petsclib::$UnionPetscLib, jac::PetscMat, X::PetscVec, Con::PetscVec, XL::PetscVec, XU::PetscVec, mu::$PetscReal, T1::PetscVec, T2::PetscVec, Da::PetscVec, Db::PetscVec, Dm::PetscVec )

    @chk ccall(
               (:MatDSFischer, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec, CVec, CVec, $PetscReal, CVec, CVec, CVec, CVec, CVec),
               jac, X, Con, XL, XU, mu, T1, T2, Da, Db, Dm,
              )


	return nothing
end 

"""
	J::PetscMat = MatCreateSubMatrixFree(petsclib::PetscLibType,mat::PetscMat, Rows::IS, Cols::IS) 
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
function MatCreateSubMatrixFree(petsclib::PetscLibType, mat::PetscMat, Rows::IS, Cols::IS) end

@for_petsc function MatCreateSubMatrixFree(petsclib::$UnionPetscLib, mat::PetscMat, Rows::IS, Cols::IS )
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
	MatSetValue(petsclib::PetscLibType,mat::PetscMat, i::PetscInt, j::PetscInt, va::PetscScalar, mode::InsertMode) 

# External Links
$(_doc_external("Mat/MatSetValue"))
"""
function MatSetValue(petsclib::PetscLibType, mat::PetscMat, i::PetscInt, j::PetscInt, va::PetscScalar, mode::InsertMode) end

@for_petsc function MatSetValue(petsclib::$UnionPetscLib, mat::PetscMat, i::$PetscInt, j::$PetscInt, va::$PetscScalar, mode::InsertMode )

    @chk ccall(
               (:MatSetValue, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, $PetscScalar, InsertMode),
               mat, i, j, va, mode,
              )


	return nothing
end 

"""
	va::PetscScalar = MatGetValue(petsclib::PetscLibType,mat::PetscMat, row::PetscInt, col::PetscInt) 

# External Links
$(_doc_external("Mat/MatGetValue"))
"""
function MatGetValue(petsclib::PetscLibType, mat::PetscMat, row::PetscInt, col::PetscInt) end

@for_petsc function MatGetValue(petsclib::$UnionPetscLib, mat::PetscMat, row::$PetscInt, col::$PetscInt )
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
	MatSetValueLocal(petsclib::PetscLibType,mat::PetscMat, i::PetscInt, j::PetscInt, va::PetscScalar, mode::InsertMode) 

# External Links
$(_doc_external("Mat/MatSetValueLocal"))
"""
function MatSetValueLocal(petsclib::PetscLibType, mat::PetscMat, i::PetscInt, j::PetscInt, va::PetscScalar, mode::InsertMode) end

@for_petsc function MatSetValueLocal(petsclib::$UnionPetscLib, mat::PetscMat, i::$PetscInt, j::$PetscInt, va::$PetscScalar, mode::InsertMode )

    @chk ccall(
               (:MatSetValueLocal, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, $PetscInt, $PetscScalar, InsertMode),
               mat, i, j, va, mode,
              )


	return nothing
end 

# NOTE:
# there are a few functions that are not listed by the python routine, but that are relevant
# Below, they are wrapped manually



"""
    MatShellSetOperation(petsclib::PetscLibType, mat::PetscMat, op::MatOperation, g::Ptr)

Allows user to set a matrix operation for a `MATSHELL` shell matrix.

Logically Collective

Input Parameters:
`mat` - the `MATSHELL` shell matrix
`op`  - the name of the operation
`g`   - a pointer to the function that provides the operation created with `@cfunction`

Level: advanced

-seealso: `Mat`, `MATSHELL`, `MatCreateShell()`, `MatShellGetContext()`, `MatShellGetOperation()`, `MatShellSetContext()`, `MatSetOperation()`, `MatShellSetManageScalingShifts()`, `MatShellSetMatProductOperation()`

# External Links
$(_doc_external("Mat/MatShellSetOperation"))
"""
function MatShellSetOperation(petsclib::PetscLibType, mat::AbstractPetscMat, op::MatOperation, g::Ptr) end

@for_petsc function MatShellSetOperation(petsclib::$UnionPetscLib, mat::AbstractPetscMat, op::MatOperation, g::Ptr)

    @chk ccall(
               (:MatShellSetOperation, $petsc_library),
               PetscErrorCode,
               (CMat, MatOperation, Ptr{Cvoid}),
               mat, op, g,
              )

	return nothing
end