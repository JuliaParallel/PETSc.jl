# autodefined type arguments for class ------
mutable struct _n_PetscRegressor end
const PetscRegressor = Ptr{_n_PetscRegressor}
# -------------------------------------------------------
"""
	PetscRegressorRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Adds a method to the `PetscRegressor` package.

Not collective

Input Parameters:
- `sname`    - name of a new user-defined regressor
- `function` - routine to create method context

-seealso: `PetscRegressorRegisterAll()`

# External Links
$(_doc_external("Ml/PetscRegressorRegister"))
"""
function PetscRegressorRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function PetscRegressorRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:PetscRegressorRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	newregressor::PetscRegressor = PetscRegressorCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Creates a `PetscRegressor` object.

Collective

Input Parameter:
- `comm` - the MPI communicator that will share the `PetscRegressor` object

Output Parameter:
- `newregressor` - the new `PetscRegressor` object

Level: beginner

-seealso: `PetscRegressorFit()`, `PetscRegressorPredict()`, `PetscRegressor`

# External Links
$(_doc_external("Ml/PetscRegressorCreate"))
"""
function PetscRegressorCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function PetscRegressorCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	newregressor_ = Ref{PetscRegressor}()

    @chk ccall(
               (:PetscRegressorCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscRegressor}),
               comm, newregressor_,
              )

	newregressor = newregressor_[]

	return newregressor
end 

"""
	PetscRegressorView(petsclib::PetscLibType,regressor::PetscRegressor, viewer::PetscViewer) 
Prints information about the `PetscRegressor` object

Collective

Input Parameters:
- `regressor` - the `PetscRegressor` context
- `viewer`    - a `PetscViewer` context

Options Database Key:
- `-regressor_view` - Calls `PetscRegressorView()` at the end of `PetscRegressorFit()`

Level: beginner

-seealso: [](ch_regressor), `PetscRegressor`, `PetscViewerASCIIOpen()`

# External Links
$(_doc_external("Ml/PetscRegressorView"))
"""
function PetscRegressorView(petsclib::PetscLibType, regressor::PetscRegressor, viewer::PetscViewer) end

@for_petsc function PetscRegressorView(petsclib::$UnionPetscLib, regressor::PetscRegressor, viewer::PetscViewer )

    @chk ccall(
               (:PetscRegressorView, $petsc_library),
               PetscErrorCode,
               (PetscRegressor, PetscViewer),
               regressor, viewer,
              )


	return nothing
end 

"""
	PetscRegressorViewFromOptions(petsclib::PetscLibType,A::PetscRegressor, obj::PetscObject, name::String) 
View a `PetscRegressor` object based on values in the options database

Collective

Input Parameters:
- `A`    - the  `PetscRegressor` context
- `obj`  - Optional object that provides the prefix for the options database
- `name` - command line option

Level: intermediate

-seealso: [](ch_regressor), `PetscRegressor`, `PetscRegressorView`, `PetscObjectViewFromOptions()`, `PetscRegressorCreate()`

# External Links
$(_doc_external("Ml/PetscRegressorViewFromOptions"))
"""
function PetscRegressorViewFromOptions(petsclib::PetscLibType, A::PetscRegressor, obj::PetscObject, name::String) end

@for_petsc function PetscRegressorViewFromOptions(petsclib::$UnionPetscLib, A::PetscRegressor, obj::PetscObject, name::String )

    @chk ccall(
               (:PetscRegressorViewFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscRegressor, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

"""
	PetscRegressorSetFromOptions(petsclib::PetscLibType,regressor::PetscRegressor) 
Sets `PetscRegressor` options from the options database.

Collective

Input Parameter:
- `regressor` - the `PetscRegressor` context

Options Database Keys:
- `-regressor_type <type>` - the particular type of regressor to be used; see `PetscRegressorType` for complete list

Level: beginner

-seealso: `PetscRegressor`, `PetscRegressorCreate()`

# External Links
$(_doc_external("Ml/PetscRegressorSetFromOptions"))
"""
function PetscRegressorSetFromOptions(petsclib::PetscLibType, regressor::PetscRegressor) end

@for_petsc function PetscRegressorSetFromOptions(petsclib::$UnionPetscLib, regressor::PetscRegressor )

    @chk ccall(
               (:PetscRegressorSetFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscRegressor,),
               regressor,
              )


	return nothing
end 

"""
	PetscRegressorSetUp(petsclib::PetscLibType,regressor::PetscRegressor) 
Sets up the internal data structures for the later use of a regressor.

Collective

Input Parameter:
- `regressor` - the `PetscRegressor` context

-seealso: `PetscRegressorCreate()`, `PetscRegressorFit()`, `PetscRegressorDestroy()`

# External Links
$(_doc_external("Ml/PetscRegressorSetUp"))
"""
function PetscRegressorSetUp(petsclib::PetscLibType, regressor::PetscRegressor) end

@for_petsc function PetscRegressorSetUp(petsclib::$UnionPetscLib, regressor::PetscRegressor )

    @chk ccall(
               (:PetscRegressorSetUp, $petsc_library),
               PetscErrorCode,
               (PetscRegressor,),
               regressor,
              )


	return nothing
end 

"""
	PetscRegressorFit(petsclib::PetscLibType,regressor::PetscRegressor, X::PetscMat, y::PetscVec) 
Fit, or train, a regressor from a training dataset

Collective

Input Parameters:
- `regressor` - the `PetscRegressor` context
- `X`         - matrix of training data (of dimension [number of samples] x [number of features])
- `y`         - vector of target values from the training dataset

Level: beginner

-seealso: `PetscRegressorCreate()`, `PetscRegressorSetUp()`, `PetscRegressorDestroy()`, `PetscRegressorPredict()`

# External Links
$(_doc_external("Ml/PetscRegressorFit"))
"""
function PetscRegressorFit(petsclib::PetscLibType, regressor::PetscRegressor, X::PetscMat, y::PetscVec) end

@for_petsc function PetscRegressorFit(petsclib::$UnionPetscLib, regressor::PetscRegressor, X::PetscMat, y::PetscVec )

    @chk ccall(
               (:PetscRegressorFit, $petsc_library),
               PetscErrorCode,
               (PetscRegressor, CMat, CVec),
               regressor, X, y,
              )


	return nothing
end 

"""
	PetscRegressorPredict(petsclib::PetscLibType,regressor::PetscRegressor, X::PetscMat, y::PetscVec) 
Compute predictions (that is, perform inference) using a fitted regression model.

Collective

Input Parameters:
- `regressor` - the `PetscRegressor` context (for which `PetscRegressorFit()` must have been called)
- `X`         - data matrix of unlabeled observations

Output Parameter:
- `y` - vector of predicted labels

Level: beginner

-seealso: `PetscRegressorFit()`, `PetscRegressorDestroy()`

# External Links
$(_doc_external("Ml/PetscRegressorPredict"))
"""
function PetscRegressorPredict(petsclib::PetscLibType, regressor::PetscRegressor, X::PetscMat, y::PetscVec) end

@for_petsc function PetscRegressorPredict(petsclib::$UnionPetscLib, regressor::PetscRegressor, X::PetscMat, y::PetscVec )

    @chk ccall(
               (:PetscRegressorPredict, $petsc_library),
               PetscErrorCode,
               (PetscRegressor, CMat, CVec),
               regressor, X, y,
              )


	return nothing
end 

"""
	PetscRegressorReset(petsclib::PetscLibType,regressor::PetscRegressor) 
Resets a `PetscRegressor` context by removing any allocated `Vec` and `Mat`. Any options set in the object remain.

Collective

Input Parameter:
- `regressor` - context obtained from `PetscRegressorCreate()`

Level: intermediate

-seealso: `PetscRegressorCreate()`, `PetscRegressorSetUp()`, `PetscRegressorFit()`, `PetscRegressorPredict()`, `PetscRegressorDestroy()`

# External Links
$(_doc_external("Ml/PetscRegressorReset"))
"""
function PetscRegressorReset(petsclib::PetscLibType, regressor::PetscRegressor) end

@for_petsc function PetscRegressorReset(petsclib::$UnionPetscLib, regressor::PetscRegressor )

    @chk ccall(
               (:PetscRegressorReset, $petsc_library),
               PetscErrorCode,
               (PetscRegressor,),
               regressor,
              )


	return nothing
end 

"""
	PetscRegressorDestroy(petsclib::PetscLibType,regressor::PetscRegressor) 
Destroys the regressor context that was created with `PetscRegressorCreate()`.

Collective

Input Parameter:
- `regressor` - the `PetscRegressor` context

Level: beginner

-seealso: `PetscRegressorCreate()`, `PetscRegressorSetUp()`, `PetscRegressorReset()`, `PetscRegressor`

# External Links
$(_doc_external("Ml/PetscRegressorDestroy"))
"""
function PetscRegressorDestroy(petsclib::PetscLibType, regressor::PetscRegressor) end

@for_petsc function PetscRegressorDestroy(petsclib::$UnionPetscLib, regressor::PetscRegressor )

    @chk ccall(
               (:PetscRegressorDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscRegressor},),
               regressor,
              )


	return nothing
end 

"""
	PetscRegressorSetType(petsclib::PetscLibType,regressor::PetscRegressor, type::PetscRegressorType) 
Sets the type for the regressor.

Collective

Input Parameters:
- `regressor` - the `PetscRegressor` context
- `type`      - a known regression method

Options Database Key:
- `-regressor_type <type>` - Sets the type of regressor; use -help for a list of available types

Level: intermediate

-seealso: `PetscRegressorType`

# External Links
$(_doc_external("Ml/PetscRegressorSetType"))
"""
function PetscRegressorSetType(petsclib::PetscLibType, regressor::PetscRegressor, type::PetscRegressorType) end

@for_petsc function PetscRegressorSetType(petsclib::$UnionPetscLib, regressor::PetscRegressor, type::PetscRegressorType )

    @chk ccall(
               (:PetscRegressorSetType, $petsc_library),
               PetscErrorCode,
               (PetscRegressor, PetscRegressorType),
               regressor, type,
              )


	return nothing
end 

"""
	type::PetscRegressorType = PetscRegressorGetType(petsclib::PetscLibType,regressor::PetscRegressor) 
Gets the current `PetscRegressorType` being used in the `PetscRegressor` object

Not Collective

Input Parameter:
- `regressor` - the `PetscRegressor` solver context

Output Parameter:
- `type` - the `PetscRegressorType`

Level: intermediate

-seealso: [](ch_regressor), `PetscRegressor`, `PetscRegressorType`, `PetscRegressorSetType()`

# External Links
$(_doc_external("Ml/PetscRegressorGetType"))
"""
function PetscRegressorGetType(petsclib::PetscLibType, regressor::PetscRegressor) end

@for_petsc function PetscRegressorGetType(petsclib::$UnionPetscLib, regressor::PetscRegressor )
	type_ = Ref{PetscRegressorType}()

    @chk ccall(
               (:PetscRegressorGetType, $petsc_library),
               PetscErrorCode,
               (PetscRegressor, Ptr{PetscRegressorType}),
               regressor, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	PetscRegressorSetRegularizerWeight(petsclib::PetscLibType,regressor::PetscRegressor, weight::PetscReal) 
Sets the weight to be used for the regularizer for a `PetscRegressor` context

Logically Collective

Input Parameters:
- `regressor` - the `PetscRegressor` context
- `weight`    - the regularizer weight

Options Database Key:
- `regressor_regularizer_weight <weight>` - sets the regularizer's weight

Level: beginner

-seealso: `PetscRegressorSetType`

# External Links
$(_doc_external("Ml/PetscRegressorSetRegularizerWeight"))
"""
function PetscRegressorSetRegularizerWeight(petsclib::PetscLibType, regressor::PetscRegressor, weight::PetscReal) end

@for_petsc function PetscRegressorSetRegularizerWeight(petsclib::$UnionPetscLib, regressor::PetscRegressor, weight::$PetscReal )

    @chk ccall(
               (:PetscRegressorSetRegularizerWeight, $petsc_library),
               PetscErrorCode,
               (PetscRegressor, $PetscReal),
               regressor, weight,
              )


	return nothing
end 

"""
	PetscRegressorGetTao(petsclib::PetscLibType,regressor::PetscRegressor, tao::Tao) 
Returns the `Tao` context for a `PetscRegressor` object.

Not Collective, but if the `PetscRegressor` is parallel, then the `Tao` object is parallel

Input Parameter:
- `regressor` - the regressor context

Output Parameter:
- `tao` - the `Tao` context

Level: beginner

-seealso: `PetscRegressorLinearGetKSP()`

# External Links
$(_doc_external("Ml/PetscRegressorGetTao"))
"""
function PetscRegressorGetTao(petsclib::PetscLibType, regressor::PetscRegressor, tao::Tao) end

@for_petsc function PetscRegressorGetTao(petsclib::$UnionPetscLib, regressor::PetscRegressor, tao::Tao )

    @chk ccall(
               (:PetscRegressorGetTao, $petsc_library),
               PetscErrorCode,
               (PetscRegressor, Ptr{Tao}),
               regressor, tao,
              )


	return nothing
end 

"""
	PetscRegressorSetOptionsPrefix(petsclib::PetscLibType,regressor::PetscRegressor, p::String) 
Sets the prefix used for searching for all
PetscRegressor options in the database.

Logically Collective

Input Parameters:
- `regressor` - the `PetscRegressor` context
- `p`         - the prefix string to prepend to all PetscRegressor option requests

Level: advanced

-seealso: [](ch_regressor), `PetscRegressor`, `PetscRegressorSetFromOptions()`, `PetscRegressorAppendOptionsPrefix()`, `PetscRegressorGetOptionsPrefix()`

# External Links
$(_doc_external("Ml/PetscRegressorSetOptionsPrefix"))
"""
function PetscRegressorSetOptionsPrefix(petsclib::PetscLibType, regressor::PetscRegressor, p::String) end

@for_petsc function PetscRegressorSetOptionsPrefix(petsclib::$UnionPetscLib, regressor::PetscRegressor, p::String )

    @chk ccall(
               (:PetscRegressorSetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (PetscRegressor, Ptr{Cchar}),
               regressor, p,
              )


	return nothing
end 

"""
	PetscRegressorAppendOptionsPrefix(petsclib::PetscLibType,regressor::PetscRegressor, p::String) 
Appends to the prefix used for searching for all PetscRegressor options in the database.

Logically Collective

Input Parameters:
- `regressor` - the `PetscRegressor` solver context
- `p`         - the prefix string to prepend to all `PetscRegressor` option requests

Level: advanced

-seealso: [](ch_regressor), `PetscRegressor`, `PetscRegressorSetFromOptions()`, `PetscRegressorSetOptionsPrefix()`, `PetscRegressorGetOptionsPrefix()`

# External Links
$(_doc_external("Ml/PetscRegressorAppendOptionsPrefix"))
"""
function PetscRegressorAppendOptionsPrefix(petsclib::PetscLibType, regressor::PetscRegressor, p::String) end

@for_petsc function PetscRegressorAppendOptionsPrefix(petsclib::$UnionPetscLib, regressor::PetscRegressor, p::String )

    @chk ccall(
               (:PetscRegressorAppendOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (PetscRegressor, Ptr{Cchar}),
               regressor, p,
              )


	return nothing
end 

"""
	PetscRegressorGetOptionsPrefix(petsclib::PetscLibType,regressor::PetscRegressor, p::String) 
Gets the prefix used for searching for all
PetscRegressor options in the database

Not Collective

Input Parameter:
- `regressor` - the `PetscRegressor` context

Output Parameter:
- `p` - pointer to the prefix string used is returned

-seealso: [](ch_regressor), `PetscRegressor`, `PetscRegressorSetFromOptions()`, `PetscRegressorSetOptionsPrefix()`, `PetscRegressorAppendOptionsPrefix()`

# External Links
$(_doc_external("Ml/PetscRegressorGetOptionsPrefix"))
"""
function PetscRegressorGetOptionsPrefix(petsclib::PetscLibType, regressor::PetscRegressor, p::String) end

@for_petsc function PetscRegressorGetOptionsPrefix(petsclib::$UnionPetscLib, regressor::PetscRegressor, p::String )
	p_ = Ref(pointer(p))

    @chk ccall(
               (:PetscRegressorGetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (PetscRegressor, Ptr{Ptr{Cchar}}),
               regressor, p_,
              )


	return nothing
end 

"""
	PetscRegressorInitializePackage(petsclib::PetscLibType) 
Initialize `PetscRegressor` package

Logically Collective

Level: developer

-seealso: `PetscRegressorFinalizePackage()`

# External Links
$(_doc_external("Ml/PetscRegressorInitializePackage"))
"""
function PetscRegressorInitializePackage(petsclib::PetscLibType) end

@for_petsc function PetscRegressorInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscRegressorInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscRegressorFinalizePackage(petsclib::PetscLibType) 
Finalize `PetscRegressor` package; it is called from `PetscFinalize()`

Logically Collective

Level: developer

-seealso: `PetscRegressorInitializePackage()`

# External Links
$(_doc_external("Ml/PetscRegressorFinalizePackage"))
"""
function PetscRegressorFinalizePackage(petsclib::PetscLibType) end

@for_petsc function PetscRegressorFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscRegressorFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscRegressorLinearSetFitIntercept(petsclib::PetscLibType,regressor::PetscRegressor, flg::PetscBool) 
Set a flag to indicate that the intercept (also known as the "bias" or "offset") should
be calculated; data are assumed to be mean-centered if false.

Logically Collective

Input Parameters:
- `regressor` - the `PetscRegressor` context
- `flg`       - `PETSC_TRUE` to calculate the intercept, `PETSC_FALSE` to assume mean-centered data (default is `PETSC_TRUE`)

Level: intermediate

Options Database Key:
- `regressor_linear_fit_intercept <true,false>` - fit the intercept

-seealso: `PetscRegressor`, `PetscRegressorFit()`

# External Links
$(_doc_external("Ml/PetscRegressorLinearSetFitIntercept"))
"""
function PetscRegressorLinearSetFitIntercept(petsclib::PetscLibType, regressor::PetscRegressor, flg::PetscBool) end

@for_petsc function PetscRegressorLinearSetFitIntercept(petsclib::$UnionPetscLib, regressor::PetscRegressor, flg::PetscBool )

    @chk ccall(
               (:PetscRegressorLinearSetFitIntercept, $petsc_library),
               PetscErrorCode,
               (PetscRegressor, PetscBool),
               regressor, flg,
              )


	return nothing
end 

"""
	PetscRegressorLinearSetUseKSP(petsclib::PetscLibType,regressor::PetscRegressor, flg::PetscBool) 
Set a flag to indicate that a `KSP` object, instead of a `Tao` one, should be used
to fit the linear regressor

Logically Collective

Input Parameters:
- `regressor` - the `PetscRegressor` context
- `flg`       - `PETSC_TRUE` to use a `KSP`, `PETSC_FALSE` to use a `Tao` object (default is false)

Options Database Key:
- `regressor_linear_use_ksp <true,false>` - use `KSP`

Level: intermediate

-seealso: `PetscRegressor`, `PetscRegressorLinearGetKSP()`, `KSPLSQR`, `PCQR`, `MATSOLVERSPQR`, `MatSolverType`, `MATSEQDENSE`, `PCSVD`

# External Links
$(_doc_external("Ml/PetscRegressorLinearSetUseKSP"))
"""
function PetscRegressorLinearSetUseKSP(petsclib::PetscLibType, regressor::PetscRegressor, flg::PetscBool) end

@for_petsc function PetscRegressorLinearSetUseKSP(petsclib::$UnionPetscLib, regressor::PetscRegressor, flg::PetscBool )

    @chk ccall(
               (:PetscRegressorLinearSetUseKSP, $petsc_library),
               PetscErrorCode,
               (PetscRegressor, PetscBool),
               regressor, flg,
              )


	return nothing
end 

"""
	PetscRegressorLinearGetKSP(petsclib::PetscLibType,regressor::PetscRegressor, ksp::PetscKSP) 
Returns the `KSP` context for a `PETSCREGRESSORLINEAR` object.

Not Collective, but if the `PetscRegressor` is parallel, then the `KSP` object is parallel

Input Parameter:
- `regressor` - the `PetscRegressor` context

Output Parameter:
- `ksp` - the `KSP` context

Level: beginner

-seealso: `PetscRegressorGetTao()`

# External Links
$(_doc_external("Ml/PetscRegressorLinearGetKSP"))
"""
function PetscRegressorLinearGetKSP(petsclib::PetscLibType, regressor::PetscRegressor, ksp::PetscKSP) end

@for_petsc function PetscRegressorLinearGetKSP(petsclib::$UnionPetscLib, regressor::PetscRegressor, ksp::PetscKSP )
	ksp_ = Ref(ksp.ptr)

    @chk ccall(
               (:PetscRegressorLinearGetKSP, $petsc_library),
               PetscErrorCode,
               (PetscRegressor, Ptr{CKSP}),
               regressor, ksp_,
              )

	ksp.ptr = C_NULL

	return nothing
end 

"""
	PetscRegressorLinearSetType(petsclib::PetscLibType,regressor::PetscRegressor, type::PetscRegressorLinearType) 
Sets the type of linear regression to be performed

Logically Collective

Input Parameters:
- `regressor` - the `PetscRegressor` context (should be of type `PETSCREGRESSORLINEAR`)
- `type`      - a known linear regression method

Options Database Key:
- `-regressor_linear_type` - Sets the linear regression method; use -help for a list of available methods
(for instance "-regressor_linear_type ols" or "-regressor_linear_type lasso")

Level: intermediate

-seealso: `PetscRegressorLinearGetType()`, `PetscRegressorLinearType`, `PetscRegressorSetType()`, `REGRESSOR_LINEAR_OLS`,
`REGRESSOR_LINEAR_LASSO`, `REGRESSOR_LINEAR_RIDGE`

# External Links
$(_doc_external("Ml/PetscRegressorLinearSetType"))
"""
function PetscRegressorLinearSetType(petsclib::PetscLibType, regressor::PetscRegressor, type::PetscRegressorLinearType) end

@for_petsc function PetscRegressorLinearSetType(petsclib::$UnionPetscLib, regressor::PetscRegressor, type::PetscRegressorLinearType )

    @chk ccall(
               (:PetscRegressorLinearSetType, $petsc_library),
               PetscErrorCode,
               (PetscRegressor, PetscRegressorLinearType),
               regressor, type,
              )


	return nothing
end 

"""
	type::PetscRegressorLinearType = PetscRegressorLinearGetType(petsclib::PetscLibType,regressor::PetscRegressor) 
Return the type for the `PETSCREGRESSORLINEAR` solver

Input Parameter:
- `regressor` - the `PetscRegressor` solver context

Output Parameter:
- `type` - `PETSCREGRESSORLINEAR` type

Level: advanced

-seealso: `PetscRegressor`, `PETSCREGRESSORLINEAR`, `PetscRegressorLinearSetType()`, `PetscRegressorLinearType`

# External Links
$(_doc_external("Ml/PetscRegressorLinearGetType"))
"""
function PetscRegressorLinearGetType(petsclib::PetscLibType, regressor::PetscRegressor) end

@for_petsc function PetscRegressorLinearGetType(petsclib::$UnionPetscLib, regressor::PetscRegressor )
	type_ = Ref{PetscRegressorLinearType}()

    @chk ccall(
               (:PetscRegressorLinearGetType, $petsc_library),
               PetscErrorCode,
               (PetscRegressor, Ptr{PetscRegressorLinearType}),
               regressor, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

