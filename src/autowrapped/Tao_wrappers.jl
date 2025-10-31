# autodefined type arguments for class ------
mutable struct _n_TaoLineSearch end
const TaoLineSearch = Ptr{_n_TaoLineSearch}

# -------------------------------------------------------
"""
	TaoFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the PETSc/Tao
interface to the Tao package. It is called from `PetscFinalize()`.

Level: developer

-seealso: `TaoInitializePackage()`, `PetscFinalize()`, `TaoRegister()`, `TaoRegisterAll()`

# External Links
$(_doc_external("Tao/TaoFinalizePackage"))
"""
function TaoFinalizePackage(petsclib::PetscLibType) end

@for_petsc function TaoFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TaoFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TaoInitializePackage(petsclib::PetscLibType) 
This function sets up PETSc to use the Tao
package.  When using static or shared libraries, this function is called from the
first entry to `TaoCreate()`; when using shared or static libraries, it is called
from PetscDLLibraryRegister_tao()

Level: developer

-seealso: `TaoCreate()`, `TaoFinalizePackage()`, `TaoRegister()`, `TaoRegisterAll()`

# External Links
$(_doc_external("Tao/TaoInitializePackage"))
"""
function TaoInitializePackage(petsclib::PetscLibType) end

@for_petsc function TaoInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TaoInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TaoSetSolution(petsclib::PetscLibType,tao::Tao, x0::PetscVec) 
Sets the vector holding the initial guess for the solve

Logically Collective

Input Parameters:
- `tao` - the `Tao` context
- `x0`  - the initial guess

Level: beginner

-seealso: [](ch_tao), `Tao`, `TaoCreate()`, `TaoSolve()`, `TaoGetSolution()`

# External Links
$(_doc_external("Tao/TaoSetSolution"))
"""
function TaoSetSolution(petsclib::PetscLibType, tao::Tao, x0::PetscVec) end

@for_petsc function TaoSetSolution(petsclib::$UnionPetscLib, tao::Tao, x0::PetscVec )

    @chk ccall(
               (:TaoSetSolution, $petsc_library),
               PetscErrorCode,
               (CTao, CVec),
               tao, x0,
              )


	return nothing
end 

"""
	TaoTestGradient(petsclib::PetscLibType,tao::Tao, x::PetscVec, g1::PetscVec) 

# External Links
$(_doc_external("Tao/TaoTestGradient"))
"""
function TaoTestGradient(petsclib::PetscLibType, tao::Tao, x::PetscVec, g1::PetscVec) end

@for_petsc function TaoTestGradient(petsclib::$UnionPetscLib, tao::Tao, x::PetscVec, g1::PetscVec )

    @chk ccall(
               (:TaoTestGradient, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, CVec),
               tao, x, g1,
              )


	return nothing
end 

"""
	TaoComputeGradient(petsclib::PetscLibType,tao::Tao, X::PetscVec, G::PetscVec) 
Computes the gradient of the objective function

Collective

Input Parameters:
- `tao` - the `Tao` context
- `X`   - input vector

Output Parameter:
- `G` - gradient vector

Options Database Keys:
- `-tao_test_gradient`      - compare the user provided gradient with one compute via finite differences to check for errors
- `-tao_test_gradient_view` - display the user provided gradient, the finite difference gradient and the difference between them to help users detect the location of errors in the user provided gradient

Level: developer

-seealso: [](ch_tao), `TaoComputeObjective()`, `TaoComputeObjectiveAndGradient()`, `TaoSetGradient()`

# External Links
$(_doc_external("Tao/TaoComputeGradient"))
"""
function TaoComputeGradient(petsclib::PetscLibType, tao::Tao, X::PetscVec, G::PetscVec) end

@for_petsc function TaoComputeGradient(petsclib::$UnionPetscLib, tao::Tao, X::PetscVec, G::PetscVec )

    @chk ccall(
               (:TaoComputeGradient, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, CVec),
               tao, X, G,
              )


	return nothing
end 

"""
	f::PetscReal = TaoComputeObjective(petsclib::PetscLibType,tao::Tao, X::PetscVec) 
Computes the objective function value at a given point

Collective

Input Parameters:
- `tao` - the `Tao` context
- `X`   - input vector

Output Parameter:
- `f` - Objective value at X

Level: developer

-seealso: [](ch_tao), `Tao`, `TaoComputeGradient()`, `TaoComputeObjectiveAndGradient()`, `TaoSetObjective()`

# External Links
$(_doc_external("Tao/TaoComputeObjective"))
"""
function TaoComputeObjective(petsclib::PetscLibType, tao::Tao, X::PetscVec) end

@for_petsc function TaoComputeObjective(petsclib::$UnionPetscLib, tao::Tao, X::PetscVec )
	f_ = Ref{$PetscReal}()

    @chk ccall(
               (:TaoComputeObjective, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, Ptr{$PetscReal}),
               tao, X, f_,
              )

	f = f_[]

	return f
end 

"""
	f::PetscReal = TaoComputeObjectiveAndGradient(petsclib::PetscLibType,tao::Tao, X::PetscVec, G::PetscVec) 
Computes the objective function value at a given point

Collective

Input Parameters:
- `tao` - the `Tao` context
- `X`   - input vector

Output Parameters:
- `f` - Objective value at `X`
- `G` - Gradient vector at `X`

Level: developer

-seealso: [](ch_tao), `TaoComputeGradient()`, `TaoSetObjective()`

# External Links
$(_doc_external("Tao/TaoComputeObjectiveAndGradient"))
"""
function TaoComputeObjectiveAndGradient(petsclib::PetscLibType, tao::Tao, X::PetscVec, G::PetscVec) end

@for_petsc function TaoComputeObjectiveAndGradient(petsclib::$UnionPetscLib, tao::Tao, X::PetscVec, G::PetscVec )
	f_ = Ref{$PetscReal}()

    @chk ccall(
               (:TaoComputeObjectiveAndGradient, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, Ptr{$PetscReal}, CVec),
               tao, X, f_, G,
              )

	f = f_[]

	return f
end 

"""
	TaoSetObjective(petsclib::PetscLibType,tao::Tao, func::external, ctx::Cvoid) 
Sets the function evaluation routine for minimization

Logically Collective

Input Parameters:
- `tao`  - the `Tao` context
- `func` - the objective function
- `ctx`  - [optional] user-defined context for private data for the function evaluation
routine (may be `NULL`)

Calling sequence of `func`:
- `tao` - the optimizer
- `x`   - input vector
- `f`   - function value
- `ctx` - [optional] user-defined function context

Level: beginner

-seealso: [](ch_tao), `TaoSetGradient()`, `TaoSetHessian()`, `TaoSetObjectiveAndGradient()`, `TaoGetObjective()`

# External Links
$(_doc_external("Tao/TaoSetObjective"))
"""
function TaoSetObjective(petsclib::PetscLibType, tao::Tao, func::external, ctx::Cvoid) end

@for_petsc function TaoSetObjective(petsclib::$UnionPetscLib, tao::Tao, func::external, ctx::Cvoid )

    @chk ccall(
               (:TaoSetObjective, $petsc_library),
               PetscErrorCode,
               (CTao, external, Ptr{Cvoid}),
               tao, func, ctx,
              )


	return nothing
end 

"""
	TaoSetResidualRoutine(petsclib::PetscLibType,tao::Tao, res::PetscVec, func::external, ctx::Cvoid) 
Sets the residual evaluation routine for least

Logically Collective

Input Parameters:
- `tao`  - the `Tao` context
- `res`  - the residual vector
- `func` - the residual evaluation routine
- `ctx`  - [optional] user-defined context for private data for the function evaluation
routine (may be `NULL`)

Calling sequence of `func`:
- `tao` - the optimizer
- `x`   - input vector
- `res` - function value vector
- `ctx` - [optional] user-defined function context

Level: beginner

-seealso: [](ch_tao), `Tao`, `TaoSetObjective()`, `TaoSetJacobianRoutine()`

# External Links
$(_doc_external("Tao/TaoSetResidualRoutine"))
"""
function TaoSetResidualRoutine(petsclib::PetscLibType, tao::Tao, res::PetscVec, func::external, ctx::Cvoid) end

@for_petsc function TaoSetResidualRoutine(petsclib::$UnionPetscLib, tao::Tao, res::PetscVec, func::external, ctx::Cvoid )

    @chk ccall(
               (:TaoSetResidualRoutine, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, external, Ptr{Cvoid}),
               tao, res, func, ctx,
              )


	return nothing
end 

"""
	TaoSetResidualWeights(petsclib::PetscLibType,tao::Tao, sigma_v::PetscVec, n::PetscInt, rows::PetscInt, cols::PetscInt, vals::PetscReal) 
Give weights for the residual values. A vector can be used if only diagonal terms are used, otherwise a matrix can be give.

Collective

Input Parameters:
- `tao`     - the `Tao` context
- `sigma_v` - vector of weights (diagonal terms only)
- `n`       - the number of weights (if using off-diagonal)
- `rows`    - index list of rows for `sigma_v`
- `cols`    - index list of columns for `sigma_v`
- `vals`    - array of weights

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoSetResidualRoutine()`

# External Links
$(_doc_external("Tao/TaoSetResidualWeights"))
"""
function TaoSetResidualWeights(petsclib::PetscLibType, tao::Tao, sigma_v::PetscVec, n::PetscInt, rows::PetscInt, cols::PetscInt, vals::PetscReal) end

@for_petsc function TaoSetResidualWeights(petsclib::$UnionPetscLib, tao::Tao, sigma_v::PetscVec, n::$PetscInt, rows::$PetscInt, cols::$PetscInt, vals::$PetscReal )

    @chk ccall(
               (:TaoSetResidualWeights, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscReal}),
               tao, sigma_v, n, rows, cols, vals,
              )


	return nothing
end 

"""
	TaoComputeResidual(petsclib::PetscLibType,tao::Tao, X::PetscVec, F::PetscVec) 
Computes a least

Collective

Input Parameters:
- `tao` - the `Tao` context
- `X`   - input vector

Output Parameter:
- `F` - Objective vector at `X`

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoSetResidualRoutine()`

# External Links
$(_doc_external("Tao/TaoComputeResidual"))
"""
function TaoComputeResidual(petsclib::PetscLibType, tao::Tao, X::PetscVec, F::PetscVec) end

@for_petsc function TaoComputeResidual(petsclib::$UnionPetscLib, tao::Tao, X::PetscVec, F::PetscVec )

    @chk ccall(
               (:TaoComputeResidual, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, CVec),
               tao, X, F,
              )


	return nothing
end 

"""
	TaoSetGradient(petsclib::PetscLibType,tao::Tao, g::PetscVec, func::external, ctx::Cvoid) 
Sets the gradient evaluation routine for the function to be optimized

Logically Collective

Input Parameters:
- `tao`  - the `Tao` context
- `g`    - [optional] the vector to internally hold the gradient computation
- `func` - the gradient function
- `ctx`  - [optional] user-defined context for private data for the gradient evaluation
routine (may be `NULL`)

Calling sequence of `func`:
- `tao` - the optimization solver
- `x`   - input vector
- `g`   - gradient value (output)
- `ctx` - [optional] user-defined function context

Level: beginner

-seealso: [](ch_tao), `Tao`, `TaoSolve()`, `TaoSetObjective()`, `TaoSetHessian()`, `TaoSetObjectiveAndGradient()`, `TaoGetGradient()`

# External Links
$(_doc_external("Tao/TaoSetGradient"))
"""
function TaoSetGradient(petsclib::PetscLibType, tao::Tao, g::PetscVec, func::external, ctx::Cvoid) end

@for_petsc function TaoSetGradient(petsclib::$UnionPetscLib, tao::Tao, g::PetscVec, func::external, ctx::Cvoid )

    @chk ccall(
               (:TaoSetGradient, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, external, Ptr{Cvoid}),
               tao, g, func, ctx,
              )


	return nothing
end 

"""
	TaoSetObjectiveAndGradient(petsclib::PetscLibType,tao::Tao, g::PetscVec, func::external, ctx::Cvoid) 
Sets a combined objective function and gradient evaluation routine for the function to be optimized

Logically Collective

Input Parameters:
- `tao`  - the `Tao` context
- `g`    - [optional] the vector to internally hold the gradient computation
- `func` - the gradient function
- `ctx`  - [optional] user-defined context for private data for the gradient evaluation
routine (may be `NULL`)

Calling sequence of `func`:
- `tao` - the optimization object
- `x`   - input vector
- `f`   - objective value (output)
- `g`   - gradient value (output)
- `ctx` - [optional] user-defined function context

Level: beginner

-seealso: [](ch_tao), `Tao`, `TaoSolve()`, `TaoSetObjective()`, `TaoSetHessian()`, `TaoSetGradient()`, `TaoGetObjectiveAndGradient()`

# External Links
$(_doc_external("Tao/TaoSetObjectiveAndGradient"))
"""
function TaoSetObjectiveAndGradient(petsclib::PetscLibType, tao::Tao, g::PetscVec, func::external, ctx::Cvoid) end

@for_petsc function TaoSetObjectiveAndGradient(petsclib::$UnionPetscLib, tao::Tao, g::PetscVec, func::external, ctx::Cvoid )

    @chk ccall(
               (:TaoSetObjectiveAndGradient, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, external, Ptr{Cvoid}),
               tao, g, func, ctx,
              )


	return nothing
end 

"""
	flg::PetscBool = TaoIsObjectiveDefined(petsclib::PetscLibType,tao::Tao) 
Checks to see if the user has
declared an objective-only routine.  Useful for determining when
it is appropriate to call `TaoComputeObjective()` or
`TaoComputeObjectiveAndGradient()`

Not Collective

Input Parameter:
- `tao` - the `Tao` context

Output Parameter:
- `flg` - `PETSC_TRUE` if function routine is set by user, `PETSC_FALSE` otherwise

Level: developer

-seealso: [](ch_tao), `Tao`, `TaoSetObjective()`, `TaoIsGradientDefined()`, `TaoIsObjectiveAndGradientDefined()`

# External Links
$(_doc_external("Tao/TaoIsObjectiveDefined"))
"""
function TaoIsObjectiveDefined(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoIsObjectiveDefined(petsclib::$UnionPetscLib, tao::Tao )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:TaoIsObjectiveDefined, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{PetscBool}),
               tao, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	flg::PetscBool = TaoIsGradientDefined(petsclib::PetscLibType,tao::Tao) 
Checks to see if the user has
declared an objective-only routine.  Useful for determining when
it is appropriate to call `TaoComputeGradient()` or
`TaoComputeGradientAndGradient()`

Not Collective

Input Parameter:
- `tao` - the `Tao` context

Output Parameter:
- `flg` - `PETSC_TRUE` if function routine is set by user, `PETSC_FALSE` otherwise

Level: developer

-seealso: [](ch_tao), `TaoSetGradient()`, `TaoIsObjectiveDefined()`, `TaoIsObjectiveAndGradientDefined()`

# External Links
$(_doc_external("Tao/TaoIsGradientDefined"))
"""
function TaoIsGradientDefined(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoIsGradientDefined(petsclib::$UnionPetscLib, tao::Tao )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:TaoIsGradientDefined, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{PetscBool}),
               tao, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	flg::PetscBool = TaoIsObjectiveAndGradientDefined(petsclib::PetscLibType,tao::Tao) 
Checks to see if the user has
declared a joint objective/gradient routine.  Useful for determining when
it is appropriate to call `TaoComputeObjective()` or
`TaoComputeObjectiveAndGradient()`

Not Collective

Input Parameter:
- `tao` - the `Tao` context

Output Parameter:
- `flg` - `PETSC_TRUE` if function routine is set by user, `PETSC_FALSE` otherwise

Level: developer

-seealso: [](ch_tao), `TaoSetObjectiveAndGradient()`, `TaoIsObjectiveDefined()`, `TaoIsGradientDefined()`

# External Links
$(_doc_external("Tao/TaoIsObjectiveAndGradientDefined"))
"""
function TaoIsObjectiveAndGradientDefined(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoIsObjectiveAndGradientDefined(petsclib::$UnionPetscLib, tao::Tao )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:TaoIsObjectiveAndGradientDefined, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{PetscBool}),
               tao, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	TaoSetHessian(petsclib::PetscLibType,tao::Tao, H::PetscMat, Hpre::PetscMat, func::external, ctx::Cvoid) 
Sets the function to compute the Hessian as well as the location to store the matrix.

Logically Collective

Input Parameters:
- `tao`  - the `Tao` context
- `H`    - Matrix used for the hessian
- `Hpre` - Matrix that will be used to construct the preconditioner, can be same as `H`
- `func` - Hessian evaluation routine
- `ctx`  - [optional] user-defined context for private data for the
Hessian evaluation routine (may be `NULL`)

Calling sequence of `func`:
- `tao`  - the `Tao`  context
- `x`    - input vector
- `H`    - Hessian matrix
- `Hpre` - matrix used to construct the preconditioner, usually the same as `H`
- `ctx`  - [optional] user-defined Hessian context

Level: beginner

-seealso: [](ch_tao), `Tao`, `TaoTypes`, `TaoSetObjective()`, `TaoSetGradient()`, `TaoSetObjectiveAndGradient()`, `TaoGetHessian()`

# External Links
$(_doc_external("Tao/TaoSetHessian"))
"""
function TaoSetHessian(petsclib::PetscLibType, tao::Tao, H::PetscMat, Hpre::PetscMat, func::external, ctx::Cvoid) end

@for_petsc function TaoSetHessian(petsclib::$UnionPetscLib, tao::Tao, H::PetscMat, Hpre::PetscMat, func::external, ctx::Cvoid )

    @chk ccall(
               (:TaoSetHessian, $petsc_library),
               PetscErrorCode,
               (CTao, CMat, CMat, external, Ptr{Cvoid}),
               tao, H, Hpre, func, ctx,
              )


	return nothing
end 

"""
	TaoTestHessian(petsclib::PetscLibType,tao::Tao) 

# External Links
$(_doc_external("Tao/TaoTestHessian"))
"""
function TaoTestHessian(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoTestHessian(petsclib::$UnionPetscLib, tao::Tao )

    @chk ccall(
               (:TaoTestHessian, $petsc_library),
               PetscErrorCode,
               (CTao,),
               tao,
              )


	return nothing
end 

"""
	TaoComputeHessian(petsclib::PetscLibType,tao::Tao, X::PetscVec, H::PetscMat, Hpre::PetscMat) 
Computes the Hessian matrix that has been
set with `TaoSetHessian()`.

Collective

Input Parameters:
- `tao` - the Tao solver context
- `X`   - input vector

Output Parameters:
- `H`    - Hessian matrix
- `Hpre` - matrix used to construct the preconditioner, usually the same as `H`

Options Database Keys:
- `-tao_test_hessian`                   - compare the user provided Hessian with one compute via finite differences to check for errors
- `-tao_test_hessian <numerical value>` - display entries in the difference between the user provided Hessian and finite difference Hessian that are greater than a certain value to help users detect errors
- `-tao_test_hessian_view`              - display the user provided Hessian, the finite difference Hessian and the difference between them to help users detect the location of errors in the user provided Hessian

Level: developer

-seealso: [](ch_tao), `Tao`, `TaoComputeObjective()`, `TaoComputeObjectiveAndGradient()`, `TaoSetHessian()`

# External Links
$(_doc_external("Tao/TaoComputeHessian"))
"""
function TaoComputeHessian(petsclib::PetscLibType, tao::Tao, X::PetscVec, H::PetscMat, Hpre::PetscMat) end

@for_petsc function TaoComputeHessian(petsclib::$UnionPetscLib, tao::Tao, X::PetscVec, H::PetscMat, Hpre::PetscMat )

    @chk ccall(
               (:TaoComputeHessian, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, CMat, CMat),
               tao, X, H, Hpre,
              )


	return nothing
end 

"""
	TaoComputeJacobian(petsclib::PetscLibType,tao::Tao, X::PetscVec, J::PetscMat, Jpre::PetscMat) 
Computes the Jacobian matrix that has been
set with TaoSetJacobianRoutine().

Collective

Input Parameters:
- `tao` - the Tao solver context
- `X`   - input vector

Output Parameters:
- `J`    - Jacobian matrix
- `Jpre` - matrix used to compute the preconditioner, often the same as `J`

Level: developer

-seealso: [](ch_tao), `TaoComputeObjective()`, `TaoComputeObjectiveAndGradient()`, `TaoSetJacobianRoutine()`

# External Links
$(_doc_external("Tao/TaoComputeJacobian"))
"""
function TaoComputeJacobian(petsclib::PetscLibType, tao::Tao, X::PetscVec, J::PetscMat, Jpre::PetscMat) end

@for_petsc function TaoComputeJacobian(petsclib::$UnionPetscLib, tao::Tao, X::PetscVec, J::PetscMat, Jpre::PetscMat )

    @chk ccall(
               (:TaoComputeJacobian, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, CMat, CMat),
               tao, X, J, Jpre,
              )


	return nothing
end 

"""
	TaoComputeResidualJacobian(petsclib::PetscLibType,tao::Tao, X::PetscVec, J::PetscMat, Jpre::PetscMat) 
Computes the least
set with `TaoSetJacobianResidual()`.

Collective

Input Parameters:
- `tao` - the Tao solver context
- `X`   - input vector

Output Parameters:
- `J`    - Jacobian matrix
- `Jpre` - matrix used to compute the preconditioner, often the same as `J`

Level: developer

-seealso: [](ch_tao), `Tao`, `TaoComputeResidual()`, `TaoSetJacobianResidual()`

# External Links
$(_doc_external("Tao/TaoComputeResidualJacobian"))
"""
function TaoComputeResidualJacobian(petsclib::PetscLibType, tao::Tao, X::PetscVec, J::PetscMat, Jpre::PetscMat) end

@for_petsc function TaoComputeResidualJacobian(petsclib::$UnionPetscLib, tao::Tao, X::PetscVec, J::PetscMat, Jpre::PetscMat )

    @chk ccall(
               (:TaoComputeResidualJacobian, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, CMat, CMat),
               tao, X, J, Jpre,
              )


	return nothing
end 

"""
	TaoComputeJacobianState(petsclib::PetscLibType,tao::Tao, X::PetscVec, J::PetscMat, Jpre::PetscMat, Jinv::PetscMat) 
Computes the Jacobian matrix that has been
set with `TaoSetJacobianStateRoutine()`.

Collective

Input Parameters:
- `tao` - the `Tao` solver context
- `X`   - input vector

Output Parameters:
- `J`    - Jacobian matrix
- `Jpre` - matrix used to construct the preconditioner, often the same as `J`
- `Jinv` - unknown

Level: developer

-seealso: [](ch_tao), `Tao`, `TaoComputeObjective()`, `TaoComputeObjectiveAndGradient()`, `TaoSetJacobianStateRoutine()`, `TaoComputeJacobianDesign()`, `TaoSetStateDesignIS()`

# External Links
$(_doc_external("Tao/TaoComputeJacobianState"))
"""
function TaoComputeJacobianState(petsclib::PetscLibType, tao::Tao, X::PetscVec, J::PetscMat, Jpre::PetscMat, Jinv::PetscMat) end

@for_petsc function TaoComputeJacobianState(petsclib::$UnionPetscLib, tao::Tao, X::PetscVec, J::PetscMat, Jpre::PetscMat, Jinv::PetscMat )

    @chk ccall(
               (:TaoComputeJacobianState, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, CMat, CMat, CMat),
               tao, X, J, Jpre, Jinv,
              )


	return nothing
end 

"""
	TaoComputeJacobianDesign(petsclib::PetscLibType,tao::Tao, X::PetscVec, J::PetscMat) 
Computes the Jacobian matrix that has been
set with `TaoSetJacobianDesignRoutine()`.

Collective

Input Parameters:
- `tao` - the Tao solver context
- `X`   - input vector

Output Parameter:
- `J` - Jacobian matrix

Level: developer

-seealso: [](ch_tao), `Tao`, `TaoComputeObjective()`, `TaoComputeObjectiveAndGradient()`, `TaoSetJacobianDesignRoutine()`, `TaoSetStateDesignIS()`

# External Links
$(_doc_external("Tao/TaoComputeJacobianDesign"))
"""
function TaoComputeJacobianDesign(petsclib::PetscLibType, tao::Tao, X::PetscVec, J::PetscMat) end

@for_petsc function TaoComputeJacobianDesign(petsclib::$UnionPetscLib, tao::Tao, X::PetscVec, J::PetscMat )

    @chk ccall(
               (:TaoComputeJacobianDesign, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, CMat),
               tao, X, J,
              )


	return nothing
end 

"""
	TaoSetJacobianRoutine(petsclib::PetscLibType,tao::Tao, J::PetscMat, Jpre::PetscMat, func::external, ctx::Cvoid) 
Sets the function to compute the Jacobian as well as the location to store the matrix.

Logically Collective

Input Parameters:
- `tao`  - the `Tao` context
- `J`    - Matrix used for the Jacobian
- `Jpre` - Matrix that will be used to construct the preconditioner, can be same as `J`
- `func` - Jacobian evaluation routine
- `ctx`  - [optional] user-defined context for private data for the
Jacobian evaluation routine (may be `NULL`)

Calling sequence of `func`:
- `tao`  - the `Tao` context
- `x`    - input vector
- `J`    - Jacobian matrix
- `Jpre` - matrix used to construct the preconditioner, usually the same as `J`
- `ctx`  - [optional] user-defined Jacobian context

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoSetGradient()`, `TaoSetObjective()`

# External Links
$(_doc_external("Tao/TaoSetJacobianRoutine"))
"""
function TaoSetJacobianRoutine(petsclib::PetscLibType, tao::Tao, J::PetscMat, Jpre::PetscMat, func::external, ctx::Cvoid) end

@for_petsc function TaoSetJacobianRoutine(petsclib::$UnionPetscLib, tao::Tao, J::PetscMat, Jpre::PetscMat, func::external, ctx::Cvoid )

    @chk ccall(
               (:TaoSetJacobianRoutine, $petsc_library),
               PetscErrorCode,
               (CTao, CMat, CMat, external, Ptr{Cvoid}),
               tao, J, Jpre, func, ctx,
              )


	return nothing
end 

"""
	TaoSetJacobianResidualRoutine(petsclib::PetscLibType,tao::Tao, J::PetscMat, Jpre::PetscMat, func::external, ctx::Cvoid) 
Sets the function to compute the least
location to store the matrix.

Logically Collective

Input Parameters:
- `tao`  - the `Tao` context
- `J`    - Matrix used for the jacobian
- `Jpre` - Matrix that will be used to construct the preconditioner, can be same as `J`
- `func` - Jacobian evaluation routine
- `ctx`  - [optional] user-defined context for private data for the
Jacobian evaluation routine (may be `NULL`)

Calling sequence of `func`:
- `tao`  - the `Tao`  context
- `x`    - input vector
- `J`    - Jacobian matrix
- `Jpre` - matrix used to construct the preconditioner, usually the same as `J`
- `ctx`  - [optional] user-defined Jacobian context

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoSetGradient()`, `TaoSetObjective()`

# External Links
$(_doc_external("Tao/TaoSetJacobianResidualRoutine"))
"""
function TaoSetJacobianResidualRoutine(petsclib::PetscLibType, tao::Tao, J::PetscMat, Jpre::PetscMat, func::external, ctx::Cvoid) end

@for_petsc function TaoSetJacobianResidualRoutine(petsclib::$UnionPetscLib, tao::Tao, J::PetscMat, Jpre::PetscMat, func::external, ctx::Cvoid )

    @chk ccall(
               (:TaoSetJacobianResidualRoutine, $petsc_library),
               PetscErrorCode,
               (CTao, CMat, CMat, external, Ptr{Cvoid}),
               tao, J, Jpre, func, ctx,
              )


	return nothing
end 

"""
	TaoSetJacobianStateRoutine(petsclib::PetscLibType,tao::Tao, J::PetscMat, Jpre::PetscMat, Jinv::PetscMat, func::external, ctx::Cvoid) 
Sets the function to compute the Jacobian
(and its inverse) of the constraint function with respect to the state variables.
Used only for PDE-constrained optimization.

Logically Collective

Input Parameters:
- `tao`  - the `Tao` context
- `J`    - Matrix used for the Jacobian
- `Jpre` - Matrix that will be used to construct the preconditioner, can be same as `J`.  Only used if `Jinv` is `NULL`
- `Jinv` - [optional] Matrix used to apply the inverse of the state Jacobian. Use `NULL` to default to PETSc `KSP` solvers to apply the inverse.
- `func` - Jacobian evaluation routine
- `ctx`  - [optional] user-defined context for private data for the
Jacobian evaluation routine (may be `NULL`)

Calling sequence of `func`:
- `tao`  - the `Tao` context
- `x`    - input vector
- `J`    - Jacobian matrix
- `Jpre` - matrix used to construct the preconditioner, usually the same as `J`
- `Jinv` - inverse of `J`
- `ctx`  - [optional] user-defined Jacobian context

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoComputeJacobianState()`, `TaoSetJacobianDesignRoutine()`, `TaoSetStateDesignIS()`

# External Links
$(_doc_external("Tao/TaoSetJacobianStateRoutine"))
"""
function TaoSetJacobianStateRoutine(petsclib::PetscLibType, tao::Tao, J::PetscMat, Jpre::PetscMat, Jinv::PetscMat, func::external, ctx::Cvoid) end

@for_petsc function TaoSetJacobianStateRoutine(petsclib::$UnionPetscLib, tao::Tao, J::PetscMat, Jpre::PetscMat, Jinv::PetscMat, func::external, ctx::Cvoid )

    @chk ccall(
               (:TaoSetJacobianStateRoutine, $petsc_library),
               PetscErrorCode,
               (CTao, CMat, CMat, CMat, external, Ptr{Cvoid}),
               tao, J, Jpre, Jinv, func, ctx,
              )


	return nothing
end 

"""
	TaoSetJacobianDesignRoutine(petsclib::PetscLibType,tao::Tao, J::PetscMat, func::external, ctx::Cvoid) 
Sets the function to compute the Jacobian of
the constraint function with respect to the design variables.  Used only for
PDE-constrained optimization.

Logically Collective

Input Parameters:
- `tao`  - the `Tao` context
- `J`    - Matrix used for the Jacobian
- `func` - Jacobian evaluation routine
- `ctx`  - [optional] user-defined context for private data for the
Jacobian evaluation routine (may be `NULL`)

Calling sequence of `func`:
- `tao` - the `Tao` context
- `x`   - input vector
- `J`   - Jacobian matrix
- `ctx` - [optional] user-defined Jacobian context

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoComputeJacobianDesign()`, `TaoSetJacobianStateRoutine()`, `TaoSetStateDesignIS()`

# External Links
$(_doc_external("Tao/TaoSetJacobianDesignRoutine"))
"""
function TaoSetJacobianDesignRoutine(petsclib::PetscLibType, tao::Tao, J::PetscMat, func::external, ctx::Cvoid) end

@for_petsc function TaoSetJacobianDesignRoutine(petsclib::$UnionPetscLib, tao::Tao, J::PetscMat, func::external, ctx::Cvoid )

    @chk ccall(
               (:TaoSetJacobianDesignRoutine, $petsc_library),
               PetscErrorCode,
               (CTao, CMat, external, Ptr{Cvoid}),
               tao, J, func, ctx,
              )


	return nothing
end 

"""
	TaoSetStateDesignIS(petsclib::PetscLibType,tao::Tao, s_is::IS, d_is::IS) 
Indicate to the `Tao` object which variables in the
solution vector are state variables and which are design.  Only applies to
PDE-constrained optimization.

Logically Collective

Input Parameters:
- `tao`  - The `Tao` context
- `s_is` - the index set corresponding to the state variables
- `d_is` - the index set corresponding to the design variables

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoSetJacobianStateRoutine()`, `TaoSetJacobianDesignRoutine()`

# External Links
$(_doc_external("Tao/TaoSetStateDesignIS"))
"""
function TaoSetStateDesignIS(petsclib::PetscLibType, tao::Tao, s_is::IS, d_is::IS) end

@for_petsc function TaoSetStateDesignIS(petsclib::$UnionPetscLib, tao::Tao, s_is::IS, d_is::IS )

    @chk ccall(
               (:TaoSetStateDesignIS, $petsc_library),
               PetscErrorCode,
               (CTao, CIS, CIS),
               tao, s_is, d_is,
              )


	return nothing
end 

"""
	TaoComputeJacobianEquality(petsclib::PetscLibType,tao::Tao, X::PetscVec, J::PetscMat, Jpre::PetscMat) 
Computes the Jacobian matrix that has been
set with `TaoSetJacobianEqualityRoutine()`.

Collective

Input Parameters:
- `tao` - the `Tao` solver context
- `X`   - input vector

Output Parameters:
- `J`    - Jacobian matrix
- `Jpre` - matrix used to construct the preconditioner, often the same as `J`

Level: developer

-seealso: [](ch_tao), `TaoComputeObjective()`, `TaoComputeObjectiveAndGradient()`, `TaoSetJacobianStateRoutine()`, `TaoComputeJacobianDesign()`, `TaoSetStateDesignIS()`

# External Links
$(_doc_external("Tao/TaoComputeJacobianEquality"))
"""
function TaoComputeJacobianEquality(petsclib::PetscLibType, tao::Tao, X::PetscVec, J::PetscMat, Jpre::PetscMat) end

@for_petsc function TaoComputeJacobianEquality(petsclib::$UnionPetscLib, tao::Tao, X::PetscVec, J::PetscMat, Jpre::PetscMat )

    @chk ccall(
               (:TaoComputeJacobianEquality, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, CMat, CMat),
               tao, X, J, Jpre,
              )


	return nothing
end 

"""
	TaoComputeJacobianInequality(petsclib::PetscLibType,tao::Tao, X::PetscVec, J::PetscMat, Jpre::PetscMat) 
Computes the Jacobian matrix that has been
set with `TaoSetJacobianInequalityRoutine()`.

Collective

Input Parameters:
- `tao` - the `Tao` solver context
- `X`   - input vector

Output Parameters:
- `J`    - Jacobian matrix
- `Jpre` - matrix used to construct the preconditioner

Level: developer

-seealso: [](ch_tao), `Tao`, `TaoComputeObjective()`, `TaoComputeObjectiveAndGradient()`, `TaoSetJacobianStateRoutine()`, `TaoComputeJacobianDesign()`, `TaoSetStateDesignIS()`

# External Links
$(_doc_external("Tao/TaoComputeJacobianInequality"))
"""
function TaoComputeJacobianInequality(petsclib::PetscLibType, tao::Tao, X::PetscVec, J::PetscMat, Jpre::PetscMat) end

@for_petsc function TaoComputeJacobianInequality(petsclib::$UnionPetscLib, tao::Tao, X::PetscVec, J::PetscMat, Jpre::PetscMat )

    @chk ccall(
               (:TaoComputeJacobianInequality, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, CMat, CMat),
               tao, X, J, Jpre,
              )


	return nothing
end 

"""
	TaoSetJacobianEqualityRoutine(petsclib::PetscLibType,tao::Tao, J::PetscMat, Jpre::PetscMat, func::external, ctx::Cvoid) 
Sets the function to compute the Jacobian
(and its inverse) of the constraint function with respect to the equality variables.
Used only for PDE-constrained optimization.

Logically Collective

Input Parameters:
- `tao`  - the `Tao` context
- `J`    - Matrix used for the Jacobian
- `Jpre` - Matrix that will be used to construct the preconditioner, can be same as `J`.
- `func` - Jacobian evaluation routine
- `ctx`  - [optional] user-defined context for private data for the
Jacobian evaluation routine (may be `NULL`)

Calling sequence of `func`:
- `tao`  - the `Tao` context
- `x`    - input vector
- `J`    - Jacobian matrix
- `Jpre` - matrix used to construct the preconditioner, usually the same as `J`
- `ctx`  - [optional] user-defined Jacobian context

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoComputeJacobianEquality()`, `TaoSetJacobianDesignRoutine()`, `TaoSetEqualityDesignIS()`

# External Links
$(_doc_external("Tao/TaoSetJacobianEqualityRoutine"))
"""
function TaoSetJacobianEqualityRoutine(petsclib::PetscLibType, tao::Tao, J::PetscMat, Jpre::PetscMat, func::external, ctx::Cvoid) end

@for_petsc function TaoSetJacobianEqualityRoutine(petsclib::$UnionPetscLib, tao::Tao, J::PetscMat, Jpre::PetscMat, func::external, ctx::Cvoid )

    @chk ccall(
               (:TaoSetJacobianEqualityRoutine, $petsc_library),
               PetscErrorCode,
               (CTao, CMat, CMat, external, Ptr{Cvoid}),
               tao, J, Jpre, func, ctx,
              )


	return nothing
end 

"""
	TaoSetJacobianInequalityRoutine(petsclib::PetscLibType,tao::Tao, J::PetscMat, Jpre::PetscMat, func::external, ctx::Cvoid) 
Sets the function to compute the Jacobian
(and its inverse) of the constraint function with respect to the inequality variables.
Used only for PDE-constrained optimization.

Logically Collective

Input Parameters:
- `tao`  - the `Tao` context
- `J`    - Matrix used for the Jacobian
- `Jpre` - Matrix that will be used to construct the preconditioner, can be same as `J`.
- `func` - Jacobian evaluation routine
- `ctx`  - [optional] user-defined context for private data for the
Jacobian evaluation routine (may be `NULL`)

Calling sequence of `func`:
- `tao`  - the `Tao` context
- `x`    - input vector
- `J`    - Jacobian matrix
- `Jpre` - matrix used to construct the preconditioner, usually the same as `J`
- `ctx`  - [optional] user-defined Jacobian context

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoComputeJacobianInequality()`, `TaoSetJacobianDesignRoutine()`, `TaoSetInequalityDesignIS()`

# External Links
$(_doc_external("Tao/TaoSetJacobianInequalityRoutine"))
"""
function TaoSetJacobianInequalityRoutine(petsclib::PetscLibType, tao::Tao, J::PetscMat, Jpre::PetscMat, func::external, ctx::Cvoid) end

@for_petsc function TaoSetJacobianInequalityRoutine(petsclib::$UnionPetscLib, tao::Tao, J::PetscMat, Jpre::PetscMat, func::external, ctx::Cvoid )

    @chk ccall(
               (:TaoSetJacobianInequalityRoutine, $petsc_library),
               PetscErrorCode,
               (CTao, CMat, CMat, external, Ptr{Cvoid}),
               tao, J, Jpre, func, ctx,
              )


	return nothing
end 

"""
	TaoSetVariableBounds(petsclib::PetscLibType,tao::Tao, XL::PetscVec, XU::PetscVec) 
Sets the upper and lower bounds for the optimization problem

Logically Collective

Input Parameters:
- `tao` - the `Tao` context
- `XL`  - vector of lower bounds
- `XU`  - vector of upper bounds

Level: beginner

-seealso: [](ch_tao), `Tao`, `TaoSetObjective()`, `TaoSetHessian()`, `TaoSetObjectiveAndGradient()`, `TaoGetVariableBounds()`

# External Links
$(_doc_external("Tao/TaoSetVariableBounds"))
"""
function TaoSetVariableBounds(petsclib::PetscLibType, tao::Tao, XL::PetscVec, XU::PetscVec) end

@for_petsc function TaoSetVariableBounds(petsclib::$UnionPetscLib, tao::Tao, XL::PetscVec, XU::PetscVec )

    @chk ccall(
               (:TaoSetVariableBounds, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, CVec),
               tao, XL, XU,
              )


	return nothing
end 

"""
	TaoSetVariableBoundsRoutine(petsclib::PetscLibType,tao::Tao, func::external, ctx::Cvoid) 
Sets a function to be used to compute lower and upper variable bounds for the optimization

Logically Collective

Input Parameters:
- `tao`  - the `Tao` context
- `func` - the bounds computation routine
- `ctx`  - [optional] user-defined context for private data for the bounds computation (may be `NULL`)

Calling sequence of `func`:
- `tao` - the `Tao` solver
- `xl`  - vector of lower bounds
- `xu`  - vector of upper bounds
- `ctx` - the (optional) user-defined function context

Level: beginner

-seealso: [](ch_tao), `Tao`, `TaoSetObjective()`, `TaoSetHessian()`, `TaoSetObjectiveAndGradient()`, `TaoSetVariableBounds()`

# External Links
$(_doc_external("Tao/TaoSetVariableBoundsRoutine"))
"""
function TaoSetVariableBoundsRoutine(petsclib::PetscLibType, tao::Tao, func::external, ctx::Cvoid) end

@for_petsc function TaoSetVariableBoundsRoutine(petsclib::$UnionPetscLib, tao::Tao, func::external, ctx::Cvoid )

    @chk ccall(
               (:TaoSetVariableBoundsRoutine, $petsc_library),
               PetscErrorCode,
               (CTao, external, Ptr{Cvoid}),
               tao, func, ctx,
              )


	return nothing
end 

"""
	TaoGetVariableBounds(petsclib::PetscLibType,tao::Tao, XL::PetscVec, XU::PetscVec) 
Gets the upper and lower bounds vectors set with `TaoSetVariableBounds()`

Not Collective

Input Parameter:
- `tao` - the `Tao` context

Output Parameters:
- `XL` - vector of lower bounds
- `XU` - vector of upper bounds

Level: beginner

-seealso: [](ch_tao), `Tao`, `TaoSetObjective()`, `TaoSetHessian()`, `TaoSetObjectiveAndGradient()`, `TaoSetVariableBounds()`

# External Links
$(_doc_external("Tao/TaoGetVariableBounds"))
"""
function TaoGetVariableBounds(petsclib::PetscLibType, tao::Tao, XL::PetscVec, XU::PetscVec) end

@for_petsc function TaoGetVariableBounds(petsclib::$UnionPetscLib, tao::Tao, XL::PetscVec, XU::PetscVec )
	XL_ = Ref(XL.ptr)
	XU_ = Ref(XU.ptr)

    @chk ccall(
               (:TaoGetVariableBounds, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{CVec}, Ptr{CVec}),
               tao, XL_, XU_,
              )

	XL.ptr = C_NULL
	XU.ptr = C_NULL

	return nothing
end 

"""
	TaoComputeVariableBounds(petsclib::PetscLibType,tao::Tao) 
Compute the variable bounds using the
routine set by `TaoSetVariableBoundsRoutine()`.

Collective

Input Parameter:
- `tao` - the `Tao` context

Level: developer

-seealso: [](ch_tao), `Tao`, `TaoSetVariableBoundsRoutine()`, `TaoSetVariableBounds()`

# External Links
$(_doc_external("Tao/TaoComputeVariableBounds"))
"""
function TaoComputeVariableBounds(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoComputeVariableBounds(petsclib::$UnionPetscLib, tao::Tao )

    @chk ccall(
               (:TaoComputeVariableBounds, $petsc_library),
               PetscErrorCode,
               (CTao,),
               tao,
              )


	return nothing
end 

"""
	TaoSetInequalityBounds(petsclib::PetscLibType,tao::Tao, IL::PetscVec, IU::PetscVec) 
Sets the upper and lower bounds

Logically Collective

Input Parameters:
- `tao` - the `Tao` context
- `IL`  - vector of lower bounds
- `IU`  - vector of upper bounds

Level: beginner

-seealso: [](ch_tao), `Tao`, `TaoSetObjective()`, `TaoSetHessian()`, `TaoSetObjectiveAndGradient()`, `TaoGetInequalityBounds()`

# External Links
$(_doc_external("Tao/TaoSetInequalityBounds"))
"""
function TaoSetInequalityBounds(petsclib::PetscLibType, tao::Tao, IL::PetscVec, IU::PetscVec) end

@for_petsc function TaoSetInequalityBounds(petsclib::$UnionPetscLib, tao::Tao, IL::PetscVec, IU::PetscVec )

    @chk ccall(
               (:TaoSetInequalityBounds, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, CVec),
               tao, IL, IU,
              )


	return nothing
end 

"""
	TaoGetInequalityBounds(petsclib::PetscLibType,tao::Tao, IL::PetscVec, IU::PetscVec) 
Gets the upper and lower bounds set via `TaoSetInequalityBounds()`

Logically Collective

Input Parameter:
- `tao` - the `Tao` context

Output Parameters:
- `IL` - vector of lower bounds
- `IU` - vector of upper bounds

Level: beginner

-seealso: [](ch_tao), `TaoSetObjective()`, `TaoSetHessian()`, `TaoSetObjectiveAndGradient()`, `TaoSetInequalityBounds()`

# External Links
$(_doc_external("Tao/TaoGetInequalityBounds"))
"""
function TaoGetInequalityBounds(petsclib::PetscLibType, tao::Tao, IL::PetscVec, IU::PetscVec) end

@for_petsc function TaoGetInequalityBounds(petsclib::$UnionPetscLib, tao::Tao, IL::PetscVec, IU::PetscVec )
	IL_ = Ref(IL.ptr)
	IU_ = Ref(IU.ptr)

    @chk ccall(
               (:TaoGetInequalityBounds, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{CVec}, Ptr{CVec}),
               tao, IL_, IU_,
              )

	IL.ptr = C_NULL
	IU.ptr = C_NULL

	return nothing
end 

"""
	TaoComputeConstraints(petsclib::PetscLibType,tao::Tao, X::PetscVec, C::PetscVec) 
Compute the variable bounds using the
routine set by `TaoSetConstraintsRoutine()`.

Collective

Input Parameters:
- `tao` - the `Tao` context
- `X`   - location to evaluate the constraints

Output Parameter:
- `C` - the constraints

Level: developer

-seealso: [](ch_tao), `Tao`, `TaoSetConstraintsRoutine()`, `TaoComputeJacobian()`

# External Links
$(_doc_external("Tao/TaoComputeConstraints"))
"""
function TaoComputeConstraints(petsclib::PetscLibType, tao::Tao, X::PetscVec, C::PetscVec) end

@for_petsc function TaoComputeConstraints(petsclib::$UnionPetscLib, tao::Tao, X::PetscVec, C::PetscVec )

    @chk ccall(
               (:TaoComputeConstraints, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, CVec),
               tao, X, C,
              )


	return nothing
end 

"""
	TaoSetConstraintsRoutine(petsclib::PetscLibType,tao::Tao, c::PetscVec, func::external, ctx::Cvoid) 
Sets a function to be used to compute constraints.  Tao only handles constraints under certain conditions, see [](ch_tao) for details

Logically Collective

Input Parameters:
- `tao`  - the `Tao` context
- `c`    - A vector that will be used to store constraint evaluation
- `func` - the bounds computation routine
- `ctx`  - [optional] user-defined context for private data for the constraints computation (may be `NULL`)

Calling sequence of `func`:
- `tao` - the `Tao` solver
- `x`   - point to evaluate constraints
- `c`   - vector constraints evaluated at `x`
- `ctx` - the (optional) user-defined function context

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoSetObjective()`, `TaoSetHessian()`, `TaoSetObjectiveAndGradient()`, `TaoSetVariablevBounds()`

# External Links
$(_doc_external("Tao/TaoSetConstraintsRoutine"))
"""
function TaoSetConstraintsRoutine(petsclib::PetscLibType, tao::Tao, c::PetscVec, func::external, ctx::Cvoid) end

@for_petsc function TaoSetConstraintsRoutine(petsclib::$UnionPetscLib, tao::Tao, c::PetscVec, func::external, ctx::Cvoid )

    @chk ccall(
               (:TaoSetConstraintsRoutine, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, external, Ptr{Cvoid}),
               tao, c, func, ctx,
              )


	return nothing
end 

"""
	TaoComputeDualVariables(petsclib::PetscLibType,tao::Tao, DL::PetscVec, DU::PetscVec) 
Computes the dual vectors corresponding to the bounds
of the variables

Collective

Input Parameter:
- `tao` - the `Tao` context

Output Parameters:
- `DL` - dual variable vector for the lower bounds
- `DU` - dual variable vector for the upper bounds

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoComputeObjective()`, `TaoSetVariableBounds()`

# External Links
$(_doc_external("Tao/TaoComputeDualVariables"))
"""
function TaoComputeDualVariables(petsclib::PetscLibType, tao::Tao, DL::PetscVec, DU::PetscVec) end

@for_petsc function TaoComputeDualVariables(petsclib::$UnionPetscLib, tao::Tao, DL::PetscVec, DU::PetscVec )

    @chk ccall(
               (:TaoComputeDualVariables, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, CVec),
               tao, DL, DU,
              )


	return nothing
end 

"""
	TaoGetDualVariables(petsclib::PetscLibType,tao::Tao, DE::PetscVec, DI::PetscVec) 
Gets the dual vectors

Collective

Input Parameter:
- `tao` - the `Tao` context

Output Parameters:
- `DE` - dual variable vector for the lower bounds
- `DI` - dual variable vector for the upper bounds

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoComputeDualVariables()`

# External Links
$(_doc_external("Tao/TaoGetDualVariables"))
"""
function TaoGetDualVariables(petsclib::PetscLibType, tao::Tao, DE::PetscVec, DI::PetscVec) end

@for_petsc function TaoGetDualVariables(petsclib::$UnionPetscLib, tao::Tao, DE::PetscVec, DI::PetscVec )
	DE_ = Ref(DE.ptr)
	DI_ = Ref(DI.ptr)

    @chk ccall(
               (:TaoGetDualVariables, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{CVec}, Ptr{CVec}),
               tao, DE_, DI_,
              )

	DE.ptr = C_NULL
	DI.ptr = C_NULL

	return nothing
end 

"""
	TaoSetEqualityConstraintsRoutine(petsclib::PetscLibType,tao::Tao, ce::PetscVec, func::external, ctx::Cvoid) 
Sets a function to be used to compute constraints.  Tao only handles constraints under certain conditions, see [](ch_tao) for details

Logically Collective

Input Parameters:
- `tao`  - the `Tao` context
- `ce`   - A vector that will be used to store equality constraint evaluation
- `func` - the bounds computation routine
- `ctx`  - [optional] user-defined context for private data for the equality constraints computation (may be `NULL`)

Calling sequence of `func`:
- `tao` - the `Tao` solver
- `x`   - point to evaluate equality constraints
- `ce`  - vector of equality constraints evaluated at x
- `ctx` - the (optional) user-defined function context

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoSetObjective()`, `TaoSetHessian()`, `TaoSetObjectiveAndGradient()`, `TaoSetVariableBounds()`

# External Links
$(_doc_external("Tao/TaoSetEqualityConstraintsRoutine"))
"""
function TaoSetEqualityConstraintsRoutine(petsclib::PetscLibType, tao::Tao, ce::PetscVec, func::external, ctx::Cvoid) end

@for_petsc function TaoSetEqualityConstraintsRoutine(petsclib::$UnionPetscLib, tao::Tao, ce::PetscVec, func::external, ctx::Cvoid )

    @chk ccall(
               (:TaoSetEqualityConstraintsRoutine, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, external, Ptr{Cvoid}),
               tao, ce, func, ctx,
              )


	return nothing
end 

"""
	TaoSetInequalityConstraintsRoutine(petsclib::PetscLibType,tao::Tao, ci::PetscVec, func::external, ctx::Cvoid) 
Sets a function to be used to compute constraints.  Tao only handles constraints under certain conditions, see [](ch_tao) for details

Logically Collective

Input Parameters:
- `tao`  - the `Tao` context
- `ci`   - A vector that will be used to store inequality constraint evaluation
- `func` - the bounds computation routine
- `ctx`  - [optional] user-defined context for private data for the inequality constraints computation (may be `NULL`)

Calling sequence of `func`:
- `tao` - the `Tao` solver
- `x`   - point to evaluate inequality constraints
- `ci`  - vector of inequality constraints evaluated at x
- `ctx` - the (optional) user-defined function context

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoSetObjective()`, `TaoSetHessian()`, `TaoSetObjectiveAndGradient()`, `TaoSetVariableBounds()`

# External Links
$(_doc_external("Tao/TaoSetInequalityConstraintsRoutine"))
"""
function TaoSetInequalityConstraintsRoutine(petsclib::PetscLibType, tao::Tao, ci::PetscVec, func::external, ctx::Cvoid) end

@for_petsc function TaoSetInequalityConstraintsRoutine(petsclib::$UnionPetscLib, tao::Tao, ci::PetscVec, func::external, ctx::Cvoid )

    @chk ccall(
               (:TaoSetInequalityConstraintsRoutine, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, external, Ptr{Cvoid}),
               tao, ci, func, ctx,
              )


	return nothing
end 

"""
	TaoComputeEqualityConstraints(petsclib::PetscLibType,tao::Tao, X::PetscVec, CE::PetscVec) 
Compute the variable bounds using the
routine set by `TaoSetEqualityConstraintsRoutine()`.

Collective

Input Parameter:
- `tao` - the `Tao` context

Output Parameters:
- `X`  - point the equality constraints were evaluated on
- `CE` - vector of equality constraints evaluated at X

Level: developer

-seealso: [](ch_tao), `Tao`, `TaoSetEqualityConstraintsRoutine()`, `TaoComputeJacobianEquality()`, `TaoComputeInequalityConstraints()`

# External Links
$(_doc_external("Tao/TaoComputeEqualityConstraints"))
"""
function TaoComputeEqualityConstraints(petsclib::PetscLibType, tao::Tao, X::PetscVec, CE::PetscVec) end

@for_petsc function TaoComputeEqualityConstraints(petsclib::$UnionPetscLib, tao::Tao, X::PetscVec, CE::PetscVec )

    @chk ccall(
               (:TaoComputeEqualityConstraints, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, CVec),
               tao, X, CE,
              )


	return nothing
end 

"""
	TaoComputeInequalityConstraints(petsclib::PetscLibType,tao::Tao, X::PetscVec, CI::PetscVec) 
Compute the variable bounds using the
routine set by `TaoSetInequalityConstraintsRoutine()`.

Collective

Input Parameter:
- `tao` - the `Tao` context

Output Parameters:
- `X`  - point the inequality constraints were evaluated on
- `CI` - vector of inequality constraints evaluated at X

Level: developer

-seealso: [](ch_tao), `Tao`, `TaoSetInequalityConstraintsRoutine()`, `TaoComputeJacobianInequality()`, `TaoComputeEqualityConstraints()`

# External Links
$(_doc_external("Tao/TaoComputeInequalityConstraints"))
"""
function TaoComputeInequalityConstraints(petsclib::PetscLibType, tao::Tao, X::PetscVec, CI::PetscVec) end

@for_petsc function TaoComputeInequalityConstraints(petsclib::$UnionPetscLib, tao::Tao, X::PetscVec, CI::PetscVec )

    @chk ccall(
               (:TaoComputeInequalityConstraints, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, CVec),
               tao, X, CI,
              )


	return nothing
end 

"""
	TaoDefaultComputeGradient(petsclib::PetscLibType,tao::Tao, Xin::PetscVec, G::PetscVec, dummy::Cvoid) 
computes the gradient using finite differences.

Collective

Input Parameters:
- `tao`   - the Tao context
- `Xin`   - compute gradient at this point
- `dummy` - not used

Output Parameter:
- `G` - Gradient Vector

Options Database Key:
- `-tao_fd_gradient`      - activates TaoDefaultComputeGradient()
- `-tao_fd_delta <delta>` - change in X used to calculate finite differences

Level: advanced

-seealso: `Tao`, `TaoSetGradient()`

# External Links
$(_doc_external("Tao/TaoDefaultComputeGradient"))
"""
function TaoDefaultComputeGradient(petsclib::PetscLibType, tao::Tao, Xin::PetscVec, G::PetscVec, dummy::Cvoid) end

@for_petsc function TaoDefaultComputeGradient(petsclib::$UnionPetscLib, tao::Tao, Xin::PetscVec, G::PetscVec, dummy::Cvoid )

    @chk ccall(
               (:TaoDefaultComputeGradient, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, CVec, Ptr{Cvoid}),
               tao, Xin, G, dummy,
              )


	return nothing
end 

"""
	TaoDefaultComputeHessian(petsclib::PetscLibType,tao::Tao, V::PetscVec, H::PetscMat, B::PetscMat, dummy::Cvoid) 
Computes the Hessian using finite differences.

Collective

Input Parameters:
- `tao`   - the Tao context
- `V`     - compute Hessian at this point
- `dummy` - not used

Output Parameters:
- `H` - Hessian matrix (not altered in this routine)
- `B` - newly computed Hessian matrix to use with preconditioner (generally the same as H)

Options Database Key:
- `-tao_fd_hessian` - activates TaoDefaultComputeHessian()

Level: advanced

-seealso: `Tao`, `TaoSetHessian()`, `TaoDefaultComputeHessianColor()`, `SNESComputeJacobianDefault()`, `TaoSetGradient()`, `TaoDefaultComputeGradient()`

# External Links
$(_doc_external("Tao/TaoDefaultComputeHessian"))
"""
function TaoDefaultComputeHessian(petsclib::PetscLibType, tao::Tao, V::PetscVec, H::PetscMat, B::PetscMat, dummy::Cvoid) end

@for_petsc function TaoDefaultComputeHessian(petsclib::$UnionPetscLib, tao::Tao, V::PetscVec, H::PetscMat, B::PetscMat, dummy::Cvoid )

    @chk ccall(
               (:TaoDefaultComputeHessian, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, CMat, CMat, Ptr{Cvoid}),
               tao, V, H, B, dummy,
              )


	return nothing
end 

"""
	TaoDefaultComputeHessianColor(petsclib::PetscLibType,tao::Tao, V::PetscVec, H::PetscMat, B::PetscMat, ctx::Cvoid) 
Computes the Hessian using colored finite differences.

Collective

Input Parameters:
- `tao` - the Tao context
- `V`   - compute Hessian at this point
- `ctx` - the color object of type `MatFDColoring`

Output Parameters:
- `H` - Hessian matrix (not altered in this routine)
- `B` - newly computed Hessian matrix to use with preconditioner (generally the same as H)

Level: advanced

-seealso: `Tao`, `MatColoring`, `TaoSetHessian()`, `TaoDefaultComputeHessian()`, `SNESComputeJacobianDefaultColor()`, `TaoSetGradient()`

# External Links
$(_doc_external("Tao/TaoDefaultComputeHessianColor"))
"""
function TaoDefaultComputeHessianColor(petsclib::PetscLibType, tao::Tao, V::PetscVec, H::PetscMat, B::PetscMat, ctx::Cvoid) end

@for_petsc function TaoDefaultComputeHessianColor(petsclib::$UnionPetscLib, tao::Tao, V::PetscVec, H::PetscMat, B::PetscMat, ctx::Cvoid )

    @chk ccall(
               (:TaoDefaultComputeHessianColor, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, CMat, CMat, Ptr{Cvoid}),
               tao, V, H, B, ctx,
              )


	return nothing
end 

"""
	TaoDefaultComputeHessianMFFD(petsclib::PetscLibType,tao::Tao, X::PetscVec, H::PetscMat, B::PetscMat, ctx::Cvoid) 

# External Links
$(_doc_external("Tao/TaoDefaultComputeHessianMFFD"))
"""
function TaoDefaultComputeHessianMFFD(petsclib::PetscLibType, tao::Tao, X::PetscVec, H::PetscMat, B::PetscMat, ctx::Cvoid) end

@for_petsc function TaoDefaultComputeHessianMFFD(petsclib::$UnionPetscLib, tao::Tao, X::PetscVec, H::PetscMat, B::PetscMat, ctx::Cvoid )

    @chk ccall(
               (:TaoDefaultComputeHessianMFFD, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, CMat, CMat, Ptr{Cvoid}),
               tao, X, H, B, ctx,
              )


	return nothing
end 

"""
	TaoParametersInitialize(petsclib::PetscLibType,tao::Tao) 
Sets all the parameters in `tao` to their default value (when `TaoCreate()` was called) if they
currently contain default values. Default values are the parameter values when the object's type is set.

Collective

Input Parameter:
- `tao` - the `Tao` object

Level: developer

-seealso: [](ch_snes), `Tao`, `TaoSolve()`, `TaoDestroy()`,
`PetscObjectParameterSetDefault()`

# External Links
$(_doc_external("Tao/TaoParametersInitialize"))
"""
function TaoParametersInitialize(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoParametersInitialize(petsclib::$UnionPetscLib, tao::Tao )

    @chk ccall(
               (:TaoParametersInitialize, $petsc_library),
               PetscErrorCode,
               (CTao,),
               tao,
              )


	return nothing
end 

"""
	newtao::Tao = TaoCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Creates a Tao solver

Collective

Input Parameter:
- `comm` - MPI communicator

Output Parameter:
- `newtao` - the new `Tao` context

Options Database Key:
- `-tao_type` - select which method Tao should use

Level: beginner

-seealso: [](ch_tao), `Tao`, `TaoSolve()`, `TaoDestroy()`, `TaoSetFromOptions()`, `TaoSetType()`

# External Links
$(_doc_external("Tao/TaoCreate"))
"""
function TaoCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function TaoCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	newtao_ = Ref{CTao}()

    @chk ccall(
               (:TaoCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{CTao}),
               comm, newtao_,
              )

	newtao = Tao(newtao_[], petsclib)

	return newtao
end 

"""
	TaoSolve(petsclib::PetscLibType,tao::Tao) 
Solves an optimization problem min F(x) s.t. l <= x <= u

Collective

Input Parameter:
- `tao` - the `Tao` context

Level: beginner

-seealso: [](ch_tao), `Tao`, `TaoCreate()`, `TaoSetObjective()`, `TaoSetGradient()`, `TaoSetHessian()`, `TaoGetConvergedReason()`, `TaoSetUp()`

# External Links
$(_doc_external("Tao/TaoSolve"))
"""
function TaoSolve(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoSolve(petsclib::$UnionPetscLib, tao::Tao )

    @chk ccall(
               (:TaoSolve, $petsc_library),
               PetscErrorCode,
               (CTao,),
               tao,
              )


	return nothing
end 

"""
	TaoSetUp(petsclib::PetscLibType,tao::Tao) 
Sets up the internal data structures for the later use
of a Tao solver

Collective

Input Parameter:
- `tao` - the `Tao` context

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoCreate()`, `TaoSolve()`

# External Links
$(_doc_external("Tao/TaoSetUp"))
"""
function TaoSetUp(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoSetUp(petsclib::$UnionPetscLib, tao::Tao )

    @chk ccall(
               (:TaoSetUp, $petsc_library),
               PetscErrorCode,
               (CTao,),
               tao,
              )


	return nothing
end 

"""
	TaoDestroy(petsclib::PetscLibType,tao::Tao) 
Destroys the `Tao` context that was created with `TaoCreate()`

Collective

Input Parameter:
- `tao` - the `Tao` context

Level: beginner

-seealso: [](ch_tao), `Tao`, `TaoCreate()`, `TaoSolve()`

# External Links
$(_doc_external("Tao/TaoDestroy"))
"""
function TaoDestroy(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoDestroy(petsclib::$UnionPetscLib, tao::Tao )
	tao_ = Ref(tao.ptr)

    @chk ccall(
               (:TaoDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{CTao},),
               tao_,
              )

	tao.ptr = C_NULL

	return nothing
end 

"""
	TaoKSPSetUseEW(petsclib::PetscLibType,tao::Tao, flag::PetscBool) 
Sets `SNES` to use Eisenstat

Logically Collective

Input Parameters:
- `tao`  - Tao context
- `flag` - `PETSC_TRUE` or `PETSC_FALSE`

Level: advanced

-seealso: [](ch_tao), `Tao`, `SNESKSPSetUseEW()`

# External Links
$(_doc_external("Tao/TaoKSPSetUseEW"))
"""
function TaoKSPSetUseEW(petsclib::PetscLibType, tao::Tao, flag::PetscBool) end

@for_petsc function TaoKSPSetUseEW(petsclib::$UnionPetscLib, tao::Tao, flag::PetscBool )

    @chk ccall(
               (:TaoKSPSetUseEW, $petsc_library),
               PetscErrorCode,
               (CTao, PetscBool),
               tao, flag,
              )


	return nothing
end 

"""
	TaoSetFromOptions(petsclib::PetscLibType,tao::Tao) 
Sets various Tao parameters from the options database

Collective

Input Parameter:
- `tao` - the `Tao` solver context

Options Database Keys:
- `-tao_type <type>`             - The algorithm that Tao uses (lmvm, nls, etc.)
- `-tao_gatol <gatol>`           - absolute error tolerance for ||gradient||
- `-tao_grtol <grtol>`           - relative error tolerance for ||gradient||
- `-tao_gttol <gttol>`           - reduction of ||gradient|| relative to initial gradient
- `-tao_max_it <max>`            - sets maximum number of iterations
- `-tao_max_funcs <max>`         - sets maximum number of function evaluations
- `-tao_fmin <fmin>`             - stop if function value reaches fmin
- `-tao_steptol <tol>`           - stop if trust region radius less than <tol>
- `-tao_trust0 <t>`              - initial trust region radius
- `-tao_view_solution`           - view the solution at the end of the optimization process
- `-tao_monitor`                 - prints function value and residual norm at each iteration
- `-tao_monitor_short`           - same as `-tao_monitor`, but truncates very small values
- `-tao_monitor_constraint_norm` - prints objective value, gradient, and constraint norm at each iteration
- `-tao_monitor_globalization`   - prints information about the globalization at each iteration
- `-tao_monitor_solution`        - prints solution vector at each iteration
- `-tao_monitor_ls_residual`     - prints least-squares residual vector at each iteration
- `-tao_monitor_step`            - prints step vector at each iteration
- `-tao_monitor_gradient`        - prints gradient vector at each iteration
- `-tao_monitor_solution_draw`   - graphically view solution vector at each iteration
- `-tao_monitor_step_draw`       - graphically view step vector at each iteration
- `-tao_monitor_gradient_draw`   - graphically view gradient at each iteration
- `-tao_monitor_cancel`          - cancels all monitors (except those set with command line)
- `-tao_fd_gradient`             - use gradient computed with finite differences
- `-tao_fd_hessian`              - use hessian computed with finite differences
- `-tao_mf_hessian`              - use matrix-free Hessian computed with finite differences
- `-tao_view`                    - prints information about the Tao after solving
- `-tao_converged_reason`        - prints the reason Tao stopped iterating

Level: beginner

-seealso: [](ch_tao), `Tao`, `TaoCreate()`, `TaoSolve()`

# External Links
$(_doc_external("Tao/TaoSetFromOptions"))
"""
function TaoSetFromOptions(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoSetFromOptions(petsclib::$UnionPetscLib, tao::Tao )

    @chk ccall(
               (:TaoSetFromOptions, $petsc_library),
               PetscErrorCode,
               (CTao,),
               tao,
              )


	return nothing
end 

"""
	TaoViewFromOptions(petsclib::PetscLibType,A::Tao, obj::PetscObject, name::String) 
View a `Tao` object based on values in the options database

Collective

Input Parameters:
- `A`    - the  `Tao` context
- `obj`  - Optional object that provides the prefix for the options database
- `name` - command line option

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoView`, `PetscObjectViewFromOptions()`, `TaoCreate()`

# External Links
$(_doc_external("Tao/TaoViewFromOptions"))
"""
function TaoViewFromOptions(petsclib::PetscLibType, A::Tao, obj::PetscObject, name::String) end

@for_petsc function TaoViewFromOptions(petsclib::$UnionPetscLib, A::Tao, obj::PetscObject, name::String )

    @chk ccall(
               (:TaoViewFromOptions, $petsc_library),
               PetscErrorCode,
               (CTao, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

"""
	TaoView(petsclib::PetscLibType,tao::Tao, viewer::PetscViewer) 
Prints information about the `Tao` object

Collective

Input Parameters:
- `tao`    - the `Tao` context
- `viewer` - visualization context

Options Database Key:
- `-tao_view` - Calls `TaoView()` at the end of `TaoSolve()`

Level: beginner

-seealso: [](ch_tao), `Tao`, `PetscViewerASCIIOpen()`

# External Links
$(_doc_external("Tao/TaoView"))
"""
function TaoView(petsclib::PetscLibType, tao::Tao, viewer::PetscViewer) end

@for_petsc function TaoView(petsclib::$UnionPetscLib, tao::Tao, viewer::PetscViewer )

    @chk ccall(
               (:TaoView, $petsc_library),
               PetscErrorCode,
               (CTao, PetscViewer),
               tao, viewer,
              )


	return nothing
end 

"""
	TaoSetRecycleHistory(petsclib::PetscLibType,tao::Tao, recycle::PetscBool) 
Sets the boolean flag to enable/disable re
iterate information from the previous `TaoSolve()`. This feature is disabled by
default.

Logically Collective

Input Parameters:
- `tao`     - the `Tao` context
- `recycle` - boolean flag

Options Database Key:
- `-tao_recycle_history <true,false>` - reuse the history

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoGetRecycleHistory()`, `TAOBNCG`, `TAOBQNLS`, `TAOBQNKLS`, `TAOBQNKTR`, `TAOBQNKTL`

# External Links
$(_doc_external("Tao/TaoSetRecycleHistory"))
"""
function TaoSetRecycleHistory(petsclib::PetscLibType, tao::Tao, recycle::PetscBool) end

@for_petsc function TaoSetRecycleHistory(petsclib::$UnionPetscLib, tao::Tao, recycle::PetscBool )

    @chk ccall(
               (:TaoSetRecycleHistory, $petsc_library),
               PetscErrorCode,
               (CTao, PetscBool),
               tao, recycle,
              )


	return nothing
end 

"""
	recycle::PetscBool = TaoGetRecycleHistory(petsclib::PetscLibType,tao::Tao) 
Retrieve the boolean flag for re
from the previous `TaoSolve()`. This feature is disabled by default.

Logically Collective

Input Parameter:
- `tao` - the `Tao` context

Output Parameter:
- `recycle` - boolean flag

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoSetRecycleHistory()`, `TAOBNCG`, `TAOBQNLS`, `TAOBQNKLS`, `TAOBQNKTR`, `TAOBQNKTL`

# External Links
$(_doc_external("Tao/TaoGetRecycleHistory"))
"""
function TaoGetRecycleHistory(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoGetRecycleHistory(petsclib::$UnionPetscLib, tao::Tao )
	recycle_ = Ref{PetscBool}()

    @chk ccall(
               (:TaoGetRecycleHistory, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{PetscBool}),
               tao, recycle_,
              )

	recycle = recycle_[]

	return recycle
end 

"""
	TaoSetTolerances(petsclib::PetscLibType,tao::Tao, gatol::PetscReal, grtol::PetscReal, gttol::PetscReal) 
Sets parameters used in `TaoSolve()` convergence tests

Logically Collective

Input Parameters:
- `tao`   - the `Tao` context
- `gatol` - stop if norm of gradient is less than this
- `grtol` - stop if relative norm of gradient is less than this
- `gttol` - stop if norm of gradient is reduced by this factor

Options Database Keys:
- `-tao_gatol <gatol>` - Sets gatol
- `-tao_grtol <grtol>` - Sets grtol
- `-tao_gttol <gttol>` - Sets gttol

Stopping Criteria:
-seealso: [](ch_tao), `Tao`, `TaoConvergedReason`, `TaoGetTolerances()`

# External Links
$(_doc_external("Tao/TaoSetTolerances"))
"""
function TaoSetTolerances(petsclib::PetscLibType, tao::Tao, gatol::PetscReal, grtol::PetscReal, gttol::PetscReal) end

@for_petsc function TaoSetTolerances(petsclib::$UnionPetscLib, tao::Tao, gatol::$PetscReal, grtol::$PetscReal, gttol::$PetscReal )

    @chk ccall(
               (:TaoSetTolerances, $petsc_library),
               PetscErrorCode,
               (CTao, $PetscReal, $PetscReal, $PetscReal),
               tao, gatol, grtol, gttol,
              )


	return nothing
end 

"""
	TaoSetConstraintTolerances(petsclib::PetscLibType,tao::Tao, catol::PetscReal, crtol::PetscReal) 
Sets constraint tolerance parameters used in `TaoSolve()` convergence tests

Logically Collective

Input Parameters:
- `tao`   - the `Tao` context
- `catol` - absolute constraint tolerance, constraint norm must be less than `catol` for used for `gatol` convergence criteria
- `crtol` - relative constraint tolerance, constraint norm must be less than `crtol` for used for `gatol`, `gttol` convergence criteria

Options Database Keys:
- `-tao_catol <catol>` - Sets catol
- `-tao_crtol <crtol>` - Sets crtol

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoConvergedReason`, `TaoGetTolerances()`, `TaoGetConstraintTolerances()`, `TaoSetTolerances()`

# External Links
$(_doc_external("Tao/TaoSetConstraintTolerances"))
"""
function TaoSetConstraintTolerances(petsclib::PetscLibType, tao::Tao, catol::PetscReal, crtol::PetscReal) end

@for_petsc function TaoSetConstraintTolerances(petsclib::$UnionPetscLib, tao::Tao, catol::$PetscReal, crtol::$PetscReal )

    @chk ccall(
               (:TaoSetConstraintTolerances, $petsc_library),
               PetscErrorCode,
               (CTao, $PetscReal, $PetscReal),
               tao, catol, crtol,
              )


	return nothing
end 

"""
	catol::PetscReal,crtol::PetscReal = TaoGetConstraintTolerances(petsclib::PetscLibType,tao::Tao) 
Gets constraint tolerance parameters used in `TaoSolve()` convergence tests

Not Collective

Input Parameter:
- `tao` - the `Tao` context

Output Parameters:
- `catol` - absolute constraint tolerance, constraint norm must be less than `catol` for used for `gatol` convergence criteria
- `crtol` - relative constraint tolerance, constraint norm must be less than `crtol` for used for `gatol`, `gttol` convergence criteria

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoConvergedReasons`,`TaoGetTolerances()`, `TaoSetTolerances()`, `TaoSetConstraintTolerances()`

# External Links
$(_doc_external("Tao/TaoGetConstraintTolerances"))
"""
function TaoGetConstraintTolerances(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoGetConstraintTolerances(petsclib::$UnionPetscLib, tao::Tao )
	catol_ = Ref{$PetscReal}()
	crtol_ = Ref{$PetscReal}()

    @chk ccall(
               (:TaoGetConstraintTolerances, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{$PetscReal}, Ptr{$PetscReal}),
               tao, catol_, crtol_,
              )

	catol = catol_[]
	crtol = crtol_[]

	return catol,crtol
end 

"""
	TaoSetFunctionLowerBound(petsclib::PetscLibType,tao::Tao, fmin::PetscReal) 
Sets a bound on the solution objective value.
When an approximate solution with an objective value below this number
has been found, the solver will terminate.

Logically Collective

Input Parameters:
- `tao`  - the Tao solver context
- `fmin` - the tolerance

Options Database Key:
- `-tao_fmin <fmin>` - sets the minimum function value

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoConvergedReason`, `TaoSetTolerances()`

# External Links
$(_doc_external("Tao/TaoSetFunctionLowerBound"))
"""
function TaoSetFunctionLowerBound(petsclib::PetscLibType, tao::Tao, fmin::PetscReal) end

@for_petsc function TaoSetFunctionLowerBound(petsclib::$UnionPetscLib, tao::Tao, fmin::$PetscReal )

    @chk ccall(
               (:TaoSetFunctionLowerBound, $petsc_library),
               PetscErrorCode,
               (CTao, $PetscReal),
               tao, fmin,
              )


	return nothing
end 

"""
	fmin::PetscReal = TaoGetFunctionLowerBound(petsclib::PetscLibType,tao::Tao) 
Gets the bound on the solution objective value.
When an approximate solution with an objective value below this number
has been found, the solver will terminate.

Not Collective

Input Parameter:
- `tao` - the `Tao` solver context

Output Parameter:
- `fmin` - the minimum function value

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoConvergedReason`, `TaoSetFunctionLowerBound()`

# External Links
$(_doc_external("Tao/TaoGetFunctionLowerBound"))
"""
function TaoGetFunctionLowerBound(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoGetFunctionLowerBound(petsclib::$UnionPetscLib, tao::Tao )
	fmin_ = Ref{$PetscReal}()

    @chk ccall(
               (:TaoGetFunctionLowerBound, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{$PetscReal}),
               tao, fmin_,
              )

	fmin = fmin_[]

	return fmin
end 

"""
	TaoSetMaximumFunctionEvaluations(petsclib::PetscLibType,tao::Tao, nfcn::PetscInt) 
Sets a maximum number of function evaluations allowed for a `TaoSolve()`.

Logically Collective

Input Parameters:
- `tao`  - the `Tao` solver context
- `nfcn` - the maximum number of function evaluations (>=0), use `PETSC_UNLIMITED` to have no bound

Options Database Key:
- `-tao_max_funcs <nfcn>` - sets the maximum number of function evaluations

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoSetTolerances()`, `TaoSetMaximumIterations()`

# External Links
$(_doc_external("Tao/TaoSetMaximumFunctionEvaluations"))
"""
function TaoSetMaximumFunctionEvaluations(petsclib::PetscLibType, tao::Tao, nfcn::PetscInt) end

@for_petsc function TaoSetMaximumFunctionEvaluations(petsclib::$UnionPetscLib, tao::Tao, nfcn::$PetscInt )

    @chk ccall(
               (:TaoSetMaximumFunctionEvaluations, $petsc_library),
               PetscErrorCode,
               (CTao, $PetscInt),
               tao, nfcn,
              )


	return nothing
end 

"""
	nfcn::PetscInt = TaoGetMaximumFunctionEvaluations(petsclib::PetscLibType,tao::Tao) 
Gets a maximum number of function evaluations allowed for a `TaoSolve()`

Logically Collective

Input Parameter:
- `tao` - the `Tao` solver context

Output Parameter:
- `nfcn` - the maximum number of function evaluations

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoSetMaximumFunctionEvaluations()`, `TaoGetMaximumIterations()`

# External Links
$(_doc_external("Tao/TaoGetMaximumFunctionEvaluations"))
"""
function TaoGetMaximumFunctionEvaluations(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoGetMaximumFunctionEvaluations(petsclib::$UnionPetscLib, tao::Tao )
	nfcn_ = Ref{$PetscInt}()

    @chk ccall(
               (:TaoGetMaximumFunctionEvaluations, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{$PetscInt}),
               tao, nfcn_,
              )

	nfcn = nfcn_[]

	return nfcn
end 

"""
	nfuncs::PetscInt = TaoGetCurrentFunctionEvaluations(petsclib::PetscLibType,tao::Tao) 
Get current number of function evaluations used by a `Tao` object

Not Collective

Input Parameter:
- `tao` - the `Tao` solver context

Output Parameter:
- `nfuncs` - the current number of function evaluations (maximum between gradient and function evaluations)

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoSetMaximumFunctionEvaluations()`, `TaoGetMaximumFunctionEvaluations()`, `TaoGetMaximumIterations()`

# External Links
$(_doc_external("Tao/TaoGetCurrentFunctionEvaluations"))
"""
function TaoGetCurrentFunctionEvaluations(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoGetCurrentFunctionEvaluations(petsclib::$UnionPetscLib, tao::Tao )
	nfuncs_ = Ref{$PetscInt}()

    @chk ccall(
               (:TaoGetCurrentFunctionEvaluations, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{$PetscInt}),
               tao, nfuncs_,
              )

	nfuncs = nfuncs_[]

	return nfuncs
end 

"""
	TaoSetMaximumIterations(petsclib::PetscLibType,tao::Tao, maxits::PetscInt) 
Sets a maximum number of iterates to be used in `TaoSolve()`

Logically Collective

Input Parameters:
- `tao`    - the `Tao` solver context
- `maxits` - the maximum number of iterates (>=0), use `PETSC_UNLIMITED` to have no bound

Options Database Key:
- `-tao_max_it <its>` - sets the maximum number of iterations

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoSetTolerances()`, `TaoSetMaximumFunctionEvaluations()`

# External Links
$(_doc_external("Tao/TaoSetMaximumIterations"))
"""
function TaoSetMaximumIterations(petsclib::PetscLibType, tao::Tao, maxits::PetscInt) end

@for_petsc function TaoSetMaximumIterations(petsclib::$UnionPetscLib, tao::Tao, maxits::$PetscInt )

    @chk ccall(
               (:TaoSetMaximumIterations, $petsc_library),
               PetscErrorCode,
               (CTao, $PetscInt),
               tao, maxits,
              )


	return nothing
end 

"""
	maxits::PetscInt = TaoGetMaximumIterations(petsclib::PetscLibType,tao::Tao) 
Gets a maximum number of iterates that will be used

Not Collective

Input Parameter:
- `tao` - the `Tao` solver context

Output Parameter:
- `maxits` - the maximum number of iterates

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoSetMaximumIterations()`, `TaoGetMaximumFunctionEvaluations()`

# External Links
$(_doc_external("Tao/TaoGetMaximumIterations"))
"""
function TaoGetMaximumIterations(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoGetMaximumIterations(petsclib::$UnionPetscLib, tao::Tao )
	maxits_ = Ref{$PetscInt}()

    @chk ccall(
               (:TaoGetMaximumIterations, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{$PetscInt}),
               tao, maxits_,
              )

	maxits = maxits_[]

	return maxits
end 

"""
	TaoSetInitialTrustRegionRadius(petsclib::PetscLibType,tao::Tao, radius::PetscReal) 
Sets the initial trust region radius.

Logically Collective

Input Parameters:
- `tao`    - a `Tao` optimization solver
- `radius` - the trust region radius

Options Database Key:
- `-tao_trust0 <t0>` - sets initial trust region radius

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoGetTrustRegionRadius()`, `TaoSetTrustRegionTolerance()`, `TAONTR`

# External Links
$(_doc_external("Tao/TaoSetInitialTrustRegionRadius"))
"""
function TaoSetInitialTrustRegionRadius(petsclib::PetscLibType, tao::Tao, radius::PetscReal) end

@for_petsc function TaoSetInitialTrustRegionRadius(petsclib::$UnionPetscLib, tao::Tao, radius::$PetscReal )

    @chk ccall(
               (:TaoSetInitialTrustRegionRadius, $petsc_library),
               PetscErrorCode,
               (CTao, $PetscReal),
               tao, radius,
              )


	return nothing
end 

"""
	radius::PetscReal = TaoGetInitialTrustRegionRadius(petsclib::PetscLibType,tao::Tao) 
Gets the initial trust region radius.

Not Collective

Input Parameter:
- `tao` - a `Tao` optimization solver

Output Parameter:
- `radius` - the trust region radius

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoSetInitialTrustRegionRadius()`, `TaoGetCurrentTrustRegionRadius()`, `TAONTR`

# External Links
$(_doc_external("Tao/TaoGetInitialTrustRegionRadius"))
"""
function TaoGetInitialTrustRegionRadius(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoGetInitialTrustRegionRadius(petsclib::$UnionPetscLib, tao::Tao )
	radius_ = Ref{$PetscReal}()

    @chk ccall(
               (:TaoGetInitialTrustRegionRadius, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{$PetscReal}),
               tao, radius_,
              )

	radius = radius_[]

	return radius
end 

"""
	radius::PetscReal = TaoGetCurrentTrustRegionRadius(petsclib::PetscLibType,tao::Tao) 
Gets the current trust region radius.

Not Collective

Input Parameter:
- `tao` - a `Tao` optimization solver

Output Parameter:
- `radius` - the trust region radius

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoSetInitialTrustRegionRadius()`, `TaoGetInitialTrustRegionRadius()`, `TAONTR`

# External Links
$(_doc_external("Tao/TaoGetCurrentTrustRegionRadius"))
"""
function TaoGetCurrentTrustRegionRadius(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoGetCurrentTrustRegionRadius(petsclib::$UnionPetscLib, tao::Tao )
	radius_ = Ref{$PetscReal}()

    @chk ccall(
               (:TaoGetCurrentTrustRegionRadius, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{$PetscReal}),
               tao, radius_,
              )

	radius = radius_[]

	return radius
end 

"""
	gatol::PetscReal,grtol::PetscReal,gttol::PetscReal = TaoGetTolerances(petsclib::PetscLibType,tao::Tao) 
gets the current values of some tolerances used for the convergence testing of `TaoSolve()`

Not Collective

Input Parameter:
- `tao` - the `Tao` context

Output Parameters:
- `gatol` - stop if norm of gradient is less than this
- `grtol` - stop if relative norm of gradient is less than this
- `gttol` - stop if norm of gradient is reduced by a this factor

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoSetTolerances()`

# External Links
$(_doc_external("Tao/TaoGetTolerances"))
"""
function TaoGetTolerances(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoGetTolerances(petsclib::$UnionPetscLib, tao::Tao )
	gatol_ = Ref{$PetscReal}()
	grtol_ = Ref{$PetscReal}()
	gttol_ = Ref{$PetscReal}()

    @chk ccall(
               (:TaoGetTolerances, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               tao, gatol_, grtol_, gttol_,
              )

	gatol = gatol_[]
	grtol = grtol_[]
	gttol = gttol_[]

	return gatol,grtol,gttol
end 

"""
	TaoGetKSP(petsclib::PetscLibType,tao::Tao, ksp::PetscKSP) 
Gets the linear solver used by the optimization solver.

Not Collective

Input Parameter:
- `tao` - the `Tao` solver

Output Parameter:
- `ksp` - the `KSP` linear solver used in the optimization solver

Level: intermediate

-seealso: [](ch_tao), `Tao`, `KSP`

# External Links
$(_doc_external("Tao/TaoGetKSP"))
"""
function TaoGetKSP(petsclib::PetscLibType, tao::Tao, ksp::PetscKSP) end

@for_petsc function TaoGetKSP(petsclib::$UnionPetscLib, tao::Tao, ksp::PetscKSP )
	ksp_ = Ref(ksp.ptr)

    @chk ccall(
               (:TaoGetKSP, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{CKSP}),
               tao, ksp_,
              )

	ksp.ptr = C_NULL

	return nothing
end 

"""
	lits::PetscInt = TaoGetLinearSolveIterations(petsclib::PetscLibType,tao::Tao) 
Gets the total number of linear iterations
used by the `Tao` solver

Not Collective

Input Parameter:
- `tao` - the `Tao` context

Output Parameter:
- `lits` - number of linear iterations

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoGetKSP()`

# External Links
$(_doc_external("Tao/TaoGetLinearSolveIterations"))
"""
function TaoGetLinearSolveIterations(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoGetLinearSolveIterations(petsclib::$UnionPetscLib, tao::Tao )
	lits_ = Ref{$PetscInt}()

    @chk ccall(
               (:TaoGetLinearSolveIterations, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{$PetscInt}),
               tao, lits_,
              )

	lits = lits_[]

	return lits
end 

"""
	TaoGetLineSearch(petsclib::PetscLibType,tao::Tao, ls::TaoLineSearch) 
Gets the line search used by the optimization solver.

Not Collective

Input Parameter:
- `tao` - the `Tao` solver

Output Parameter:
- `ls` - the line search used in the optimization solver

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoLineSearch`, `TaoLineSearchType`

# External Links
$(_doc_external("Tao/TaoGetLineSearch"))
"""
function TaoGetLineSearch(petsclib::PetscLibType, tao::Tao, ls::TaoLineSearch) end

@for_petsc function TaoGetLineSearch(petsclib::$UnionPetscLib, tao::Tao, ls::TaoLineSearch )

    @chk ccall(
               (:TaoGetLineSearch, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{TaoLineSearch}),
               tao, ls,
              )


	return nothing
end 

"""
	TaoAddLineSearchCounts(petsclib::PetscLibType,tao::Tao) 
Adds the number of function evaluations spent
in the line search to the running total.

Input Parameters:
- `tao` - the `Tao` solver

Level: developer

-seealso: [](ch_tao), `Tao`, `TaoGetLineSearch()`, `TaoLineSearchApply()`

# External Links
$(_doc_external("Tao/TaoAddLineSearchCounts"))
"""
function TaoAddLineSearchCounts(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoAddLineSearchCounts(petsclib::$UnionPetscLib, tao::Tao )

    @chk ccall(
               (:TaoAddLineSearchCounts, $petsc_library),
               PetscErrorCode,
               (CTao,),
               tao,
              )


	return nothing
end 

"""
	TaoGetSolution(petsclib::PetscLibType,tao::Tao, X::PetscVec) 
Returns the vector with the current solution from the `Tao` object

Not Collective

Input Parameter:
- `tao` - the `Tao` context

Output Parameter:
- `X` - the current solution

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoSetSolution()`, `TaoSolve()`

# External Links
$(_doc_external("Tao/TaoGetSolution"))
"""
function TaoGetSolution(petsclib::PetscLibType, tao::Tao, X::PetscVec) end

@for_petsc function TaoGetSolution(petsclib::$UnionPetscLib, tao::Tao, X::PetscVec )
	X_ = Ref(X.ptr)

    @chk ccall(
               (:TaoGetSolution, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{CVec}),
               tao, X_,
              )

	X.ptr = C_NULL

	return nothing
end 

"""
	TaoResetStatistics(petsclib::PetscLibType,tao::Tao) 
Initialize the statistics collected by the `Tao` object.
These statistics include the iteration number, residual norms, and convergence status.
This routine gets called before solving each optimization problem.

Collective

Input Parameter:
- `tao` - the `Tao` context

Level: developer

-seealso: [](ch_tao), `Tao`, `TaoCreate()`, `TaoSolve()`

# External Links
$(_doc_external("Tao/TaoResetStatistics"))
"""
function TaoResetStatistics(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoResetStatistics(petsclib::$UnionPetscLib, tao::Tao )

    @chk ccall(
               (:TaoResetStatistics, $petsc_library),
               PetscErrorCode,
               (CTao,),
               tao,
              )


	return nothing
end 

"""
	TaoSetUpdate(petsclib::PetscLibType,tao::Tao, func::external, ctx::Cvoid) 
Sets the general
at the beginning of every iteration of the optimization algorithm. Called after the new solution and the gradient
is determined, but before the Hessian is computed (if applicable).

Logically Collective

Input Parameters:
- `tao`  - The `Tao` solver
- `func` - The function
- `ctx`  - The update function context

Calling sequence of `func`:
- `tao` - The optimizer context
- `it`  - The current iteration index
- `ctx` - The update context

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoSolve()`

# External Links
$(_doc_external("Tao/TaoSetUpdate"))
"""
function TaoSetUpdate(petsclib::PetscLibType, tao::Tao, func::external, ctx::Cvoid) end

@for_petsc function TaoSetUpdate(petsclib::$UnionPetscLib, tao::Tao, func::external, ctx::Cvoid )

    @chk ccall(
               (:TaoSetUpdate, $petsc_library),
               PetscErrorCode,
               (CTao, external, Ptr{Cvoid}),
               tao, func, ctx,
              )


	return nothing
end 

"""
	TaoSetConvergenceTest(petsclib::PetscLibType,tao::Tao, conv::external, ctx::Cvoid) 
Sets the function that is to be used to test
for convergence of the iterative minimization solution.  The new convergence
testing routine will replace Tao's default convergence test.

Logically Collective

Input Parameters:
- `tao`  - the `Tao` object
- `conv` - the routine to test for convergence
- `ctx`  - [optional] context for private data for the convergence routine
(may be `NULL`)

Calling sequence of `conv`:
- `tao` - the `Tao` object
- `ctx` - [optional] convergence context

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoSolve()`, `TaoSetConvergedReason()`, `TaoGetSolutionStatus()`, `TaoGetTolerances()`, `TaoMonitorSet()`

# External Links
$(_doc_external("Tao/TaoSetConvergenceTest"))
"""
function TaoSetConvergenceTest(petsclib::PetscLibType, tao::Tao, conv::external, ctx::Cvoid) end

@for_petsc function TaoSetConvergenceTest(petsclib::$UnionPetscLib, tao::Tao, conv::external, ctx::Cvoid )

    @chk ccall(
               (:TaoSetConvergenceTest, $petsc_library),
               PetscErrorCode,
               (CTao, external, Ptr{Cvoid}),
               tao, conv, ctx,
              )


	return nothing
end 

"""
	TaoMonitorSet(petsclib::PetscLibType,tao::Tao, func::external, ctx::Cvoid, dest::PetscCtxDestroyFn) 
Sets an additional function that is to be used at every
iteration of the solver to display the iteration's
progress.

Logically Collective

Input Parameters:
- `tao`  - the `Tao` solver context
- `func` - monitoring routine
- `ctx`  - [optional] user-defined context for private data for the monitor routine (may be `NULL`)
- `dest` - [optional] function to destroy the context when the `Tao` is destroyed, see `PetscCtxDestroyFn` for the calling sequence

Calling sequence of `func`:
- `tao` - the `Tao` solver context
- `ctx` - [optional] monitoring context

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoSolve()`, `TaoMonitorDefault()`, `TaoMonitorCancel()`, `TaoSetDestroyRoutine()`, `TaoView()`, `PetscCtxDestroyFn`

# External Links
$(_doc_external("Tao/TaoMonitorSet"))
"""
function TaoMonitorSet(petsclib::PetscLibType, tao::Tao, func::external, ctx::Cvoid, dest::PetscCtxDestroyFn) end

@for_petsc function TaoMonitorSet(petsclib::$UnionPetscLib, tao::Tao, func::external, ctx::Cvoid, dest::PetscCtxDestroyFn )

    @chk ccall(
               (:TaoMonitorSet, $petsc_library),
               PetscErrorCode,
               (CTao, external, Ptr{Cvoid}, Ptr{PetscCtxDestroyFn}),
               tao, func, ctx, dest,
              )


	return nothing
end 

"""
	TaoMonitorCancel(petsclib::PetscLibType,tao::Tao) 
Clears all the monitor functions for a `Tao` object.

Logically Collective

Input Parameter:
- `tao` - the `Tao` solver context

Options Database Key:
- `-tao_monitor_cancel` - cancels all monitors that have been hardwired
into a code by calls to `TaoMonitorSet()`, but does not cancel those
set via the options database

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoMonitorDefault()`, `TaoMonitorSet()`

# External Links
$(_doc_external("Tao/TaoMonitorCancel"))
"""
function TaoMonitorCancel(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoMonitorCancel(petsclib::$UnionPetscLib, tao::Tao )

    @chk ccall(
               (:TaoMonitorCancel, $petsc_library),
               PetscErrorCode,
               (CTao,),
               tao,
              )


	return nothing
end 

"""
	TaoMonitorDefault(petsclib::PetscLibType,tao::Tao, ctx::Cvoid) 
Default routine for monitoring progress of `TaoSolve()`

Collective

Input Parameters:
- `tao` - the `Tao` context
- `ctx` - `PetscViewer` context or `NULL`

Options Database Key:
- `-tao_monitor` - turn on default monitoring

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoMonitorDefaultShort()`, `TaoMonitorSet()`

# External Links
$(_doc_external("Tao/TaoMonitorDefault"))
"""
function TaoMonitorDefault(petsclib::PetscLibType, tao::Tao, ctx::Cvoid) end

@for_petsc function TaoMonitorDefault(petsclib::$UnionPetscLib, tao::Tao, ctx::Cvoid )

    @chk ccall(
               (:TaoMonitorDefault, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{Cvoid}),
               tao, ctx,
              )


	return nothing
end 

"""
	TaoMonitorGlobalization(petsclib::PetscLibType,tao::Tao, ctx::Cvoid) 
Default routine for monitoring progress of `TaoSolve()` with extra detail on the globalization method.

Collective

Input Parameters:
- `tao` - the `Tao` context
- `ctx` - `PetscViewer` context or `NULL`

Options Database Key:
- `-tao_monitor_globalization` - turn on monitoring with globalization information

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoMonitorDefaultShort()`, `TaoMonitorSet()`

# External Links
$(_doc_external("Tao/TaoMonitorGlobalization"))
"""
function TaoMonitorGlobalization(petsclib::PetscLibType, tao::Tao, ctx::Cvoid) end

@for_petsc function TaoMonitorGlobalization(petsclib::$UnionPetscLib, tao::Tao, ctx::Cvoid )

    @chk ccall(
               (:TaoMonitorGlobalization, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{Cvoid}),
               tao, ctx,
              )


	return nothing
end 

"""
	TaoMonitorDefaultShort(petsclib::PetscLibType,tao::Tao, ctx::Cvoid) 
Routine for monitoring progress of `TaoSolve()` that displays fewer digits than `TaoMonitorDefault()`

Collective

Input Parameters:
- `tao` - the `Tao` context
- `ctx` - `PetscViewer` context of type `PETSCVIEWERASCII`

Options Database Key:
- `-tao_monitor_short` - turn on default short monitoring

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoMonitorDefault()`, `TaoMonitorSet()`

# External Links
$(_doc_external("Tao/TaoMonitorDefaultShort"))
"""
function TaoMonitorDefaultShort(petsclib::PetscLibType, tao::Tao, ctx::Cvoid) end

@for_petsc function TaoMonitorDefaultShort(petsclib::$UnionPetscLib, tao::Tao, ctx::Cvoid )

    @chk ccall(
               (:TaoMonitorDefaultShort, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{Cvoid}),
               tao, ctx,
              )


	return nothing
end 

"""
	TaoMonitorConstraintNorm(petsclib::PetscLibType,tao::Tao, ctx::Cvoid) 
same as `TaoMonitorDefault()` except
it prints the norm of the constraint function.

Collective

Input Parameters:
- `tao` - the `Tao` context
- `ctx` - `PetscViewer` context or `NULL`

Options Database Key:
- `-tao_monitor_constraint_norm` - monitor the constraints

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoMonitorDefault()`, `TaoMonitorSet()`

# External Links
$(_doc_external("Tao/TaoMonitorConstraintNorm"))
"""
function TaoMonitorConstraintNorm(petsclib::PetscLibType, tao::Tao, ctx::Cvoid) end

@for_petsc function TaoMonitorConstraintNorm(petsclib::$UnionPetscLib, tao::Tao, ctx::Cvoid )

    @chk ccall(
               (:TaoMonitorConstraintNorm, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{Cvoid}),
               tao, ctx,
              )


	return nothing
end 

"""
	TaoMonitorSolution(petsclib::PetscLibType,tao::Tao, ctx::Cvoid) 
Views the solution at each iteration of `TaoSolve()`

Collective

Input Parameters:
- `tao` - the `Tao` context
- `ctx` - `PetscViewer` context or `NULL`

Options Database Key:
- `-tao_monitor_solution` - view the solution

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoMonitorDefaultShort()`, `TaoMonitorSet()`

# External Links
$(_doc_external("Tao/TaoMonitorSolution"))
"""
function TaoMonitorSolution(petsclib::PetscLibType, tao::Tao, ctx::Cvoid) end

@for_petsc function TaoMonitorSolution(petsclib::$UnionPetscLib, tao::Tao, ctx::Cvoid )

    @chk ccall(
               (:TaoMonitorSolution, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{Cvoid}),
               tao, ctx,
              )


	return nothing
end 

"""
	TaoMonitorGradient(petsclib::PetscLibType,tao::Tao, ctx::Cvoid) 
Views the gradient at each iteration of `TaoSolve()`

Collective

Input Parameters:
- `tao` - the `Tao` context
- `ctx` - `PetscViewer` context or `NULL`

Options Database Key:
- `-tao_monitor_gradient` - view the gradient at each iteration

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoMonitorDefaultShort()`, `TaoMonitorSet()`

# External Links
$(_doc_external("Tao/TaoMonitorGradient"))
"""
function TaoMonitorGradient(petsclib::PetscLibType, tao::Tao, ctx::Cvoid) end

@for_petsc function TaoMonitorGradient(petsclib::$UnionPetscLib, tao::Tao, ctx::Cvoid )

    @chk ccall(
               (:TaoMonitorGradient, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{Cvoid}),
               tao, ctx,
              )


	return nothing
end 

"""
	TaoMonitorStep(petsclib::PetscLibType,tao::Tao, ctx::Cvoid) 
Views the step

Collective

Input Parameters:
- `tao` - the `Tao` context
- `ctx` - `PetscViewer` context or `NULL`

Options Database Key:
- `-tao_monitor_step` - view the step vector at each iteration

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoMonitorDefaultShort()`, `TaoMonitorSet()`

# External Links
$(_doc_external("Tao/TaoMonitorStep"))
"""
function TaoMonitorStep(petsclib::PetscLibType, tao::Tao, ctx::Cvoid) end

@for_petsc function TaoMonitorStep(petsclib::$UnionPetscLib, tao::Tao, ctx::Cvoid )

    @chk ccall(
               (:TaoMonitorStep, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{Cvoid}),
               tao, ctx,
              )


	return nothing
end 

"""
	TaoMonitorSolutionDraw(petsclib::PetscLibType,tao::Tao, ctx::Cvoid) 
Plots the solution at each iteration of `TaoSolve()`

Collective

Input Parameters:
- `tao` - the `Tao` context
- `ctx` - `TaoMonitorDraw` context

Options Database Key:
- `-tao_monitor_solution_draw` - draw the solution at each iteration

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoMonitorSolution()`, `TaoMonitorSet()`, `TaoMonitorGradientDraw()`, `TaoMonitorDrawCtxCreate()`,
`TaoMonitorDrawCtxDestroy()`

# External Links
$(_doc_external("Tao/TaoMonitorSolutionDraw"))
"""
function TaoMonitorSolutionDraw(petsclib::PetscLibType, tao::Tao, ctx::Cvoid) end

@for_petsc function TaoMonitorSolutionDraw(petsclib::$UnionPetscLib, tao::Tao, ctx::Cvoid )

    @chk ccall(
               (:TaoMonitorSolutionDraw, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{Cvoid}),
               tao, ctx,
              )


	return nothing
end 

"""
	TaoMonitorGradientDraw(petsclib::PetscLibType,tao::Tao, ctx::Cvoid) 
Plots the gradient at each iteration of `TaoSolve()`

Collective

Input Parameters:
- `tao` - the `Tao` context
- `ctx` - `PetscViewer` context

Options Database Key:
- `-tao_monitor_gradient_draw` - draw the gradient at each iteration

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoMonitorGradient()`, `TaoMonitorSet()`, `TaoMonitorSolutionDraw()`

# External Links
$(_doc_external("Tao/TaoMonitorGradientDraw"))
"""
function TaoMonitorGradientDraw(petsclib::PetscLibType, tao::Tao, ctx::Cvoid) end

@for_petsc function TaoMonitorGradientDraw(petsclib::$UnionPetscLib, tao::Tao, ctx::Cvoid )

    @chk ccall(
               (:TaoMonitorGradientDraw, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{Cvoid}),
               tao, ctx,
              )


	return nothing
end 

"""
	TaoMonitorStepDraw(petsclib::PetscLibType,tao::Tao, ctx::Cvoid) 
Plots the step direction at each iteration of `TaoSolve()`

Collective

Input Parameters:
- `tao` - the `Tao` context
- `ctx` - the `PetscViewer` context

Options Database Key:
- `-tao_monitor_step_draw` - draw the step direction at each iteration

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoMonitorSet()`, `TaoMonitorSolutionDraw`

# External Links
$(_doc_external("Tao/TaoMonitorStepDraw"))
"""
function TaoMonitorStepDraw(petsclib::PetscLibType, tao::Tao, ctx::Cvoid) end

@for_petsc function TaoMonitorStepDraw(petsclib::$UnionPetscLib, tao::Tao, ctx::Cvoid )

    @chk ccall(
               (:TaoMonitorStepDraw, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{Cvoid}),
               tao, ctx,
              )


	return nothing
end 

"""
	TaoMonitorResidual(petsclib::PetscLibType,tao::Tao, ctx::Cvoid) 
Views the least

Collective

Input Parameters:
- `tao` - the `Tao` context
- `ctx` - the `PetscViewer` context or `NULL`

Options Database Key:
- `-tao_monitor_ls_residual` - view the residual at each iteration

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoMonitorDefaultShort()`, `TaoMonitorSet()`

# External Links
$(_doc_external("Tao/TaoMonitorResidual"))
"""
function TaoMonitorResidual(petsclib::PetscLibType, tao::Tao, ctx::Cvoid) end

@for_petsc function TaoMonitorResidual(petsclib::$UnionPetscLib, tao::Tao, ctx::Cvoid )

    @chk ccall(
               (:TaoMonitorResidual, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{Cvoid}),
               tao, ctx,
              )


	return nothing
end 

"""
	TaoDefaultConvergenceTest(petsclib::PetscLibType,tao::Tao, dummy::Cvoid) 
Determines whether the solver should continue iterating
or terminate.

Collective

Input Parameters:
- `tao`   - the `Tao` context
- `dummy` - unused dummy context

Level: developer

-seealso: [](ch_tao), `Tao`, `TaoSetTolerances()`, `TaoGetConvergedReason()`, `TaoSetConvergedReason()`

# External Links
$(_doc_external("Tao/TaoDefaultConvergenceTest"))
"""
function TaoDefaultConvergenceTest(petsclib::PetscLibType, tao::Tao, dummy::Cvoid) end

@for_petsc function TaoDefaultConvergenceTest(petsclib::$UnionPetscLib, tao::Tao, dummy::Cvoid )

    @chk ccall(
               (:TaoDefaultConvergenceTest, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{Cvoid}),
               tao, dummy,
              )


	return nothing
end 

"""
	TaoSetOptionsPrefix(petsclib::PetscLibType,tao::Tao, p::String) 
Sets the prefix used for searching for all
Tao options in the database.

Logically Collective

Input Parameters:
- `tao` - the `Tao` context
- `p`   - the prefix string to prepend to all Tao option requests

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoSetFromOptions()`, `TaoAppendOptionsPrefix()`, `TaoGetOptionsPrefix()`

# External Links
$(_doc_external("Tao/TaoSetOptionsPrefix"))
"""
function TaoSetOptionsPrefix(petsclib::PetscLibType, tao::Tao, p::String) end

@for_petsc function TaoSetOptionsPrefix(petsclib::$UnionPetscLib, tao::Tao, p::String )

    @chk ccall(
               (:TaoSetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{Cchar}),
               tao, p,
              )


	return nothing
end 

"""
	TaoAppendOptionsPrefix(petsclib::PetscLibType,tao::Tao, p::String) 
Appends to the prefix used for searching for all Tao options in the database.

Logically Collective

Input Parameters:
- `tao` - the `Tao` solver context
- `p`   - the prefix string to prepend to all `Tao` option requests

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoSetFromOptions()`, `TaoSetOptionsPrefix()`, `TaoGetOptionsPrefix()`

# External Links
$(_doc_external("Tao/TaoAppendOptionsPrefix"))
"""
function TaoAppendOptionsPrefix(petsclib::PetscLibType, tao::Tao, p::String) end

@for_petsc function TaoAppendOptionsPrefix(petsclib::$UnionPetscLib, tao::Tao, p::String )

    @chk ccall(
               (:TaoAppendOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{Cchar}),
               tao, p,
              )


	return nothing
end 

"""
	TaoGetOptionsPrefix(petsclib::PetscLibType,tao::Tao, p::String) 
Gets the prefix used for searching for all
Tao options in the database

Not Collective

Input Parameter:
- `tao` - the `Tao` context

Output Parameter:
- `p` - pointer to the prefix string used is returned

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoSetFromOptions()`, `TaoSetOptionsPrefix()`, `TaoAppendOptionsPrefix()`

# External Links
$(_doc_external("Tao/TaoGetOptionsPrefix"))
"""
function TaoGetOptionsPrefix(petsclib::PetscLibType, tao::Tao, p::String) end

@for_petsc function TaoGetOptionsPrefix(petsclib::$UnionPetscLib, tao::Tao, p::String )
	p_ = Ref(pointer(p))

    @chk ccall(
               (:TaoGetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{Ptr{Cchar}}),
               tao, p_,
              )


	return nothing
end 

"""
	TaoSetType(petsclib::PetscLibType,tao::Tao, type::TaoType) 
Sets the `TaoType` for the minimization solver.

Collective

Input Parameters:
- `tao`  - the `Tao` solver context
- `type` - a known method

Options Database Key:
- `-tao_type <type>` - Sets the method; use -help for a list
of available methods (for instance, "-tao_type lmvm" or "-tao_type tron")

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoCreate()`, `TaoGetType()`, `TaoType`

# External Links
$(_doc_external("Tao/TaoSetType"))
"""
function TaoSetType(petsclib::PetscLibType, tao::Tao, type::TaoType) end

@for_petsc function TaoSetType(petsclib::$UnionPetscLib, tao::Tao, type::TaoType )

    @chk ccall(
               (:TaoSetType, $petsc_library),
               PetscErrorCode,
               (CTao, TaoType),
               tao, type,
              )


	return nothing
end 

"""
	TaoRegister(petsclib::PetscLibType,sname::String, func::external) 
Adds a method to the Tao package for minimization.

Not Collective, No Fortran Support

Input Parameters:
- `sname` - name of a new user-defined solver
- `func`  - routine to Create method context

-seealso: [](ch_tao), `Tao`, `TaoSetType()`, `TaoRegisterAll()`, `TaoRegisterDestroy()`

# External Links
$(_doc_external("Tao/TaoRegister"))
"""
function TaoRegister(petsclib::PetscLibType, sname::String, func::external) end

@for_petsc function TaoRegister(petsclib::$UnionPetscLib, sname::String, func::external )

    @chk ccall(
               (:TaoRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, func,
              )


	return nothing
end 

"""
	TaoRegisterDestroy(petsclib::PetscLibType) 
Frees the list of minimization solvers that were
registered by `TaoRegister()`.

Not Collective

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoRegisterAll()`, `TaoRegister()`

# External Links
$(_doc_external("Tao/TaoRegisterDestroy"))
"""
function TaoRegisterDestroy(petsclib::PetscLibType) end

@for_petsc function TaoRegisterDestroy(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TaoRegisterDestroy, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	iter::PetscInt = TaoGetIterationNumber(petsclib::PetscLibType,tao::Tao) 
Gets the number of `TaoSolve()` iterations completed
at this time.

Not Collective

Input Parameter:
- `tao` - the `Tao` context

Output Parameter:
- `iter` - iteration number

-seealso: [](ch_tao), `Tao`, `TaoGetLinearSolveIterations()`, `TaoGetResidualNorm()`, `TaoGetObjective()`

# External Links
$(_doc_external("Tao/TaoGetIterationNumber"))
"""
function TaoGetIterationNumber(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoGetIterationNumber(petsclib::$UnionPetscLib, tao::Tao )
	iter_ = Ref{$PetscInt}()

    @chk ccall(
               (:TaoGetIterationNumber, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{$PetscInt}),
               tao, iter_,
              )

	iter = iter_[]

	return iter
end 

"""
	value::PetscReal = TaoGetResidualNorm(petsclib::PetscLibType,tao::Tao) 
Gets the current value of the norm of the residual (gradient)
at this time.

Not Collective

Input Parameter:
- `tao` - the `Tao` context

Output Parameter:
- `value` - the current value

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoGetLinearSolveIterations()`, `TaoGetIterationNumber()`, `TaoGetObjective()`

# External Links
$(_doc_external("Tao/TaoGetResidualNorm"))
"""
function TaoGetResidualNorm(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoGetResidualNorm(petsclib::$UnionPetscLib, tao::Tao )
	value_ = Ref{$PetscReal}()

    @chk ccall(
               (:TaoGetResidualNorm, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{$PetscReal}),
               tao, value_,
              )

	value = value_[]

	return value
end 

"""
	TaoSetIterationNumber(petsclib::PetscLibType,tao::Tao, iter::PetscInt) 
Sets the current iteration number.

Logically Collective

Input Parameters:
- `tao`  - the `Tao` context
- `iter` - iteration number

Level: developer

-seealso: [](ch_tao), `Tao`, `TaoGetLinearSolveIterations()`

# External Links
$(_doc_external("Tao/TaoSetIterationNumber"))
"""
function TaoSetIterationNumber(petsclib::PetscLibType, tao::Tao, iter::PetscInt) end

@for_petsc function TaoSetIterationNumber(petsclib::$UnionPetscLib, tao::Tao, iter::$PetscInt )

    @chk ccall(
               (:TaoSetIterationNumber, $petsc_library),
               PetscErrorCode,
               (CTao, $PetscInt),
               tao, iter,
              )


	return nothing
end 

"""
	iter::PetscInt = TaoGetTotalIterationNumber(petsclib::PetscLibType,tao::Tao) 
Gets the total number of `TaoSolve()` iterations
completed. This number keeps accumulating if multiple solves
are called with the `Tao` object.

Not Collective

Input Parameter:
- `tao` - the `Tao` context

Output Parameter:
- `iter` - number of iterations

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoGetLinearSolveIterations()`

# External Links
$(_doc_external("Tao/TaoGetTotalIterationNumber"))
"""
function TaoGetTotalIterationNumber(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoGetTotalIterationNumber(petsclib::$UnionPetscLib, tao::Tao )
	iter_ = Ref{$PetscInt}()

    @chk ccall(
               (:TaoGetTotalIterationNumber, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{$PetscInt}),
               tao, iter_,
              )

	iter = iter_[]

	return iter
end 

"""
	TaoSetTotalIterationNumber(petsclib::PetscLibType,tao::Tao, iter::PetscInt) 
Sets the current total iteration number.

Logically Collective

Input Parameters:
- `tao`  - the `Tao` context
- `iter` - the iteration number

Level: developer

-seealso: [](ch_tao), `Tao`, `TaoGetLinearSolveIterations()`

# External Links
$(_doc_external("Tao/TaoSetTotalIterationNumber"))
"""
function TaoSetTotalIterationNumber(petsclib::PetscLibType, tao::Tao, iter::PetscInt) end

@for_petsc function TaoSetTotalIterationNumber(petsclib::$UnionPetscLib, tao::Tao, iter::$PetscInt )

    @chk ccall(
               (:TaoSetTotalIterationNumber, $petsc_library),
               PetscErrorCode,
               (CTao, $PetscInt),
               tao, iter,
              )


	return nothing
end 

"""
	TaoSetConvergedReason(petsclib::PetscLibType,tao::Tao, reason::TaoConvergedReason) 
Sets the termination flag on a `Tao` object

Logically Collective

Input Parameters:
- `tao`    - the `Tao` context
- `reason` - the `TaoConvergedReason`

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoConvergedReason`

# External Links
$(_doc_external("Tao/TaoSetConvergedReason"))
"""
function TaoSetConvergedReason(petsclib::PetscLibType, tao::Tao, reason::TaoConvergedReason) end

@for_petsc function TaoSetConvergedReason(petsclib::$UnionPetscLib, tao::Tao, reason::TaoConvergedReason )

    @chk ccall(
               (:TaoSetConvergedReason, $petsc_library),
               PetscErrorCode,
               (CTao, TaoConvergedReason),
               tao, reason,
              )


	return nothing
end 

"""
	TaoGetConvergedReason(petsclib::PetscLibType,tao::Tao, reason::TaoConvergedReason) 
Gets the reason the `TaoSolve()` was stopped.

Not Collective

Input Parameter:
- `tao` - the `Tao` solver context

Output Parameter:
- `reason` - value of `TaoConvergedReason`

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoConvergedReason`, `TaoSetConvergenceTest()`, `TaoSetTolerances()`

# External Links
$(_doc_external("Tao/TaoGetConvergedReason"))
"""
function TaoGetConvergedReason(petsclib::PetscLibType, tao::Tao, reason::TaoConvergedReason) end

@for_petsc function TaoGetConvergedReason(petsclib::$UnionPetscLib, tao::Tao, reason::TaoConvergedReason )

    @chk ccall(
               (:TaoGetConvergedReason, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{TaoConvergedReason}),
               tao, reason,
              )


	return nothing
end 

"""
	its::PetscInt,f::PetscReal,gnorm::PetscReal,cnorm::PetscReal,xdiff::PetscReal = TaoGetSolutionStatus(petsclib::PetscLibType,tao::Tao, reason::TaoConvergedReason) 
Get the current iterate, objective value,
residual, infeasibility, and termination from a `Tao` object

Not Collective

Input Parameter:
- `tao` - the `Tao` context

Output Parameters:
- `its`    - the current iterate number (>=0)
- `f`      - the current function value
- `gnorm`  - the square of the gradient norm, duality gap, or other measure indicating distance from optimality.
- `cnorm`  - the infeasibility of the current solution with regard to the constraints.
- `xdiff`  - the step length or trust region radius of the most recent iterate.
- `reason` - The termination reason, which can equal `TAO_CONTINUE_ITERATING`

Level: intermediate

-seealso: [](ch_tao), `TaoMonitor()`, `TaoGetConvergedReason()`

# External Links
$(_doc_external("Tao/TaoGetSolutionStatus"))
"""
function TaoGetSolutionStatus(petsclib::PetscLibType, tao::Tao, reason::TaoConvergedReason) end

@for_petsc function TaoGetSolutionStatus(petsclib::$UnionPetscLib, tao::Tao, reason::TaoConvergedReason )
	its_ = Ref{$PetscInt}()
	f_ = Ref{$PetscReal}()
	gnorm_ = Ref{$PetscReal}()
	cnorm_ = Ref{$PetscReal}()
	xdiff_ = Ref{$PetscReal}()

    @chk ccall(
               (:TaoGetSolutionStatus, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{$PetscInt}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{TaoConvergedReason}),
               tao, its_, f_, gnorm_, cnorm_, xdiff_, reason,
              )

	its = its_[]
	f = f_[]
	gnorm = gnorm_[]
	cnorm = cnorm_[]
	xdiff = xdiff_[]

	return its,f,gnorm,cnorm,xdiff
end 

"""
	type::TaoType = TaoGetType(petsclib::PetscLibType,tao::Tao) 
Gets the current `TaoType` being used in the `Tao` object

Not Collective

Input Parameter:
- `tao` - the `Tao` solver context

Output Parameter:
- `type` - the `TaoType`

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoType`, `TaoSetType()`

# External Links
$(_doc_external("Tao/TaoGetType"))
"""
function TaoGetType(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoGetType(petsclib::$UnionPetscLib, tao::Tao )
	type_ = Ref{TaoType}()

    @chk ccall(
               (:TaoGetType, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{TaoType}),
               tao, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	TaoMonitor(petsclib::PetscLibType,tao::Tao, its::PetscInt, f::PetscReal, res::PetscReal, cnorm::PetscReal, steplength::PetscReal) 
Monitor the solver and the current solution.  This
routine will record the iteration number and residual statistics,
and call any monitors specified by the user.

Input Parameters:
- `tao`        - the `Tao` context
- `its`        - the current iterate number (>=0)
- `f`          - the current objective function value
- `res`        - the gradient norm, square root of the duality gap, or other measure indicating distance from optimality.  This measure will be recorded and
used for some termination tests.
- `cnorm`      - the infeasibility of the current solution with regard to the constraints.
- `steplength` - multiple of the step direction added to the previous iterate.

Options Database Key:
- `-tao_monitor` - Use the default monitor, which prints statistics to standard output

Level: developer

-seealso: [](ch_tao), `Tao`, `TaoGetConvergedReason()`, `TaoMonitorDefault()`, `TaoMonitorSet()`

# External Links
$(_doc_external("Tao/TaoMonitor"))
"""
function TaoMonitor(petsclib::PetscLibType, tao::Tao, its::PetscInt, f::PetscReal, res::PetscReal, cnorm::PetscReal, steplength::PetscReal) end

@for_petsc function TaoMonitor(petsclib::$UnionPetscLib, tao::Tao, its::$PetscInt, f::$PetscReal, res::$PetscReal, cnorm::$PetscReal, steplength::$PetscReal )

    @chk ccall(
               (:TaoMonitor, $petsc_library),
               PetscErrorCode,
               (CTao, $PetscInt, $PetscReal, $PetscReal, $PetscReal, $PetscReal),
               tao, its, f, res, cnorm, steplength,
              )


	return nothing
end 

"""
	TaoSetConvergenceHistory(petsclib::PetscLibType,tao::Tao, obj::Vector{PetscReal}, resid::Vector{PetscReal}, cnorm::Vector{PetscReal}, lits::Vector{PetscInt}, na::PetscInt, reset::PetscBool) 
Sets the array used to hold the convergence history.

Logically Collective

Input Parameters:
- `tao`   - the `Tao` solver context
- `obj`   - array to hold objective value history
- `resid` - array to hold residual history
- `cnorm` - array to hold constraint violation history
- `lits`  - integer array holds the number of linear iterations for each Tao iteration
- `na`    - size of `obj`, `resid`, and `cnorm`
- `reset` - `PETSC_TRUE` indicates each new minimization resets the history counter to zero,
else it continues storing new values for new minimizations after the old ones

Level: intermediate

-seealso: [](ch_tao), `TaoGetConvergenceHistory()`

# External Links
$(_doc_external("Tao/TaoSetConvergenceHistory"))
"""
function TaoSetConvergenceHistory(petsclib::PetscLibType, tao::Tao, obj::Vector{PetscReal}, resid::Vector{PetscReal}, cnorm::Vector{PetscReal}, lits::Vector{PetscInt}, na::PetscInt, reset::PetscBool) end

@for_petsc function TaoSetConvergenceHistory(petsclib::$UnionPetscLib, tao::Tao, obj::Vector{$PetscReal}, resid::Vector{$PetscReal}, cnorm::Vector{$PetscReal}, lits::Vector{$PetscInt}, na::$PetscInt, reset::PetscBool )

    @chk ccall(
               (:TaoSetConvergenceHistory, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscInt}, $PetscInt, PetscBool),
               tao, obj, resid, cnorm, lits, na, reset,
              )


	return nothing
end 

"""
	obj::PetscReal,resid::PetscReal,cnorm::PetscReal,lits::PetscInt,nhist::PetscInt = TaoGetConvergenceHistory(petsclib::PetscLibType,tao::Tao) 
Gets the arrays used that hold the convergence history.

Collective

Input Parameter:
- `tao` - the `Tao` context

Output Parameters:
- `obj`   - array used to hold objective value history
- `resid` - array used to hold residual history
- `cnorm` - array used to hold constraint violation history
- `lits`  - integer array used to hold linear solver iteration count
- `nhist` - size of `obj`, `resid`, `cnorm`, and `lits`

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoSolve()`, `TaoSetConvergenceHistory()`

# External Links
$(_doc_external("Tao/TaoGetConvergenceHistory"))
"""
function TaoGetConvergenceHistory(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoGetConvergenceHistory(petsclib::$UnionPetscLib, tao::Tao )
	obj_ = Ref{$PetscReal}()
	resid_ = Ref{$PetscReal}()
	cnorm_ = Ref{$PetscReal}()
	lits_ = Ref{$PetscInt}()
	nhist_ = Ref{$PetscInt}()

    @chk ccall(
               (:TaoGetConvergenceHistory, $petsc_library),
               PetscErrorCode,
               (CTao, $PetscReal, $PetscReal, $PetscReal, $PetscInt, Ptr{$PetscInt}),
               tao, obj_, resid_, cnorm_, lits_, nhist_,
              )

	obj = obj_[]
	resid = resid_[]
	cnorm = cnorm_[]
	lits = lits_[]
	nhist = nhist_[]

	return obj,resid,cnorm,lits,nhist
end 

"""
	TaoSetApplicationContext(petsclib::PetscLibType,tao::Tao, ctx::Cvoid) 
Sets the optional user
`Tao` callback functions with `TaoGetApplicationContext()`

Logically Collective

Input Parameters:
- `tao` - the `Tao` context
- `ctx` - the user context

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoGetApplicationContext()`

# External Links
$(_doc_external("Tao/TaoSetApplicationContext"))
"""
function TaoSetApplicationContext(petsclib::PetscLibType, tao::Tao, ctx::Cvoid) end

@for_petsc function TaoSetApplicationContext(petsclib::$UnionPetscLib, tao::Tao, ctx::Cvoid )

    @chk ccall(
               (:TaoSetApplicationContext, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{Cvoid}),
               tao, ctx,
              )


	return nothing
end 

"""
	TaoGetApplicationContext(petsclib::PetscLibType,tao::Tao, ctx::PeCtx) 
Gets the user

Not Collective

Input Parameter:
- `tao` - the `Tao` context

Output Parameter:
- `ctx` - a pointer to the user context

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoSetApplicationContext()`

# External Links
$(_doc_external("Tao/TaoGetApplicationContext"))
"""
function TaoGetApplicationContext(petsclib::PetscLibType, tao::Tao, ctx::PeCtx) end

@for_petsc function TaoGetApplicationContext(petsclib::$UnionPetscLib, tao::Tao, ctx::PeCtx )

    @chk ccall(
               (:TaoGetApplicationContext, $petsc_library),
               PetscErrorCode,
               (CTao, PeCtx),
               tao, ctx,
              )


	return nothing
end 

"""
	TaoSetGradientNorm(petsclib::PetscLibType,tao::Tao, M::PetscMat) 
Sets the matrix used to define the norm that measures the size of the gradient in some of the `Tao` algorithms

Collective

Input Parameters:
- `tao` - the `Tao` context
- `M`   - matrix that defines the norm

Level: beginner

-seealso: [](ch_tao), `Tao`, `TaoGetGradientNorm()`, `TaoGradientNorm()`

# External Links
$(_doc_external("Tao/TaoSetGradientNorm"))
"""
function TaoSetGradientNorm(petsclib::PetscLibType, tao::Tao, M::PetscMat) end

@for_petsc function TaoSetGradientNorm(petsclib::$UnionPetscLib, tao::Tao, M::PetscMat )

    @chk ccall(
               (:TaoSetGradientNorm, $petsc_library),
               PetscErrorCode,
               (CTao, CMat),
               tao, M,
              )


	return nothing
end 

"""
	TaoGetGradientNorm(petsclib::PetscLibType,tao::Tao, M::PetscMat) 
Returns the matrix used to define the norm used for measuring the size of the gradient in some of the `Tao` algorithms

Not Collective

Input Parameter:
- `tao` - the `Tao` context

Output Parameter:
- `M` - gradient norm

Level: beginner

-seealso: [](ch_tao), `Tao`, `TaoSetGradientNorm()`, `TaoGradientNorm()`

# External Links
$(_doc_external("Tao/TaoGetGradientNorm"))
"""
function TaoGetGradientNorm(petsclib::PetscLibType, tao::Tao, M::PetscMat) end

@for_petsc function TaoGetGradientNorm(petsclib::$UnionPetscLib, tao::Tao, M::PetscMat )
	M_ = Ref(M.ptr)

    @chk ccall(
               (:TaoGetGradientNorm, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{CMat}),
               tao, M_,
              )

	M.ptr = C_NULL

	return nothing
end 

"""
	gnorm::PetscReal = TaoGradientNorm(petsclib::PetscLibType,tao::Tao, gradient::PetscVec, type::NormType) 
Compute the norm using the `NormType`, the user has selected

Collective

Input Parameters:
- `tao`      - the `Tao` context
- `gradient` - the gradient
- `type`     - the norm type

Output Parameter:
- `gnorm` - the gradient norm

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoSetGradientNorm()`, `TaoGetGradientNorm()`

# External Links
$(_doc_external("Tao/TaoGradientNorm"))
"""
function TaoGradientNorm(petsclib::PetscLibType, tao::Tao, gradient::PetscVec, type::NormType) end

@for_petsc function TaoGradientNorm(petsclib::$UnionPetscLib, tao::Tao, gradient::PetscVec, type::NormType )
	gnorm_ = Ref{$PetscReal}()

    @chk ccall(
               (:TaoGradientNorm, $petsc_library),
               PetscErrorCode,
               (CTao, CVec, NormType, Ptr{$PetscReal}),
               tao, gradient, type, gnorm_,
              )

	gnorm = gnorm_[]

	return gnorm
end 

"""
	TaoSoftThreshold(petsclib::PetscLibType,in::PetscVec, lb::PetscReal, ub::PetscReal, out::PetscVec) 
Calculates soft thresholding routine with input vector
and given lower and upper bound and returns it to output vector.

Input Parameters:
- `in` - input vector to be thresholded
- `lb` - lower bound
- `ub` - upper bound

Output Parameter:
- `out` - Soft thresholded output vector

-seealso: `Tao`, `Vec`

# External Links
$(_doc_external("Tao/TaoSoftThreshold"))
"""
function TaoSoftThreshold(petsclib::PetscLibType, in::PetscVec, lb::PetscReal, ub::PetscReal, out::PetscVec) end

@for_petsc function TaoSoftThreshold(petsclib::$UnionPetscLib, in::PetscVec, lb::$PetscReal, ub::$PetscReal, out::PetscVec )

    @chk ccall(
               (:TaoSoftThreshold, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscReal, $PetscReal, CVec),
               in, lb, ub, out,
              )


	return nothing
end 

"""
	TaoPythonSetType(petsclib::PetscLibType,tao::Tao, pyname::String) 
Initialize a `Tao` object implemented in Python.

Collective

Input Parameters:
- `tao`  - the optimization solver (`Tao`) context.
- `pyname`  - full dotted Python name [package].module[.{class|function}]

Options Database Key:
- `-tao_python_type <pyname>`  - python class

Level: intermediate

-seealso: `TaoCreate()`, `TaoSetType()`, `TAOPYTHON`, `PetscPythonInitialize()`

# External Links
$(_doc_external("Tao/TaoPythonSetType"))
"""
function TaoPythonSetType(petsclib::PetscLibType, tao::Tao, pyname::String) end

@for_petsc function TaoPythonSetType(petsclib::$UnionPetscLib, tao::Tao, pyname::String )

    @chk ccall(
               (:TaoPythonSetType, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{Cchar}),
               tao, pyname,
              )


	return nothing
end 

"""
	pyname::String = TaoPythonGetType(petsclib::PetscLibType,tao::Tao) 
Get the type of a `Tao` object implemented in Python.

Not Collective

Input Parameter:
- `tao`  - the optimization solver (`Tao`) context.

Output Parameter:
- `pyname`  - full dotted Python name [package].module[.{class|function}]

Level: intermediate

-seealso: `TaoCreate()`, `TaoSetType()`, `TaoPYTHON`, `PetscPythonInitialize()`, `TaoPythonSetType()`

# External Links
$(_doc_external("Tao/TaoPythonGetType"))
"""
function TaoPythonGetType(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoPythonGetType(petsclib::$UnionPetscLib, tao::Tao )
	pyname_ = Ref{Ptr{Cchar}}()

    @chk ccall(
               (:TaoPythonGetType, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{Ptr{Cchar}}),
               tao, pyname_,
              )

	pyname = unsafe_wrap(Array, pyname_[], VecGetLocalSize(petsclib, x); own = false)

	return pyname
end 

"""
	TaoShellSetSolve(petsclib::PetscLibType,tao::Tao, solve::external) 
Sets routine to apply as solver

Logically Collective

Input Parameters:
- `tao`   - the nonlinear solver context
- `solve` - the application-provided solver routine

Calling sequence of `solve`:
- `tao` - the optimizer, get the application context with `TaoShellGetContext()`

Level: advanced

-seealso: `Tao`, `TAOSHELL`, `TaoShellSetContext()`, `TaoShellGetContext()`

# External Links
$(_doc_external("Tao/TaoShellSetSolve"))
"""
function TaoShellSetSolve(petsclib::PetscLibType, tao::Tao, solve::external) end

@for_petsc function TaoShellSetSolve(petsclib::$UnionPetscLib, tao::Tao, solve::external )

    @chk ccall(
               (:TaoShellSetSolve, $petsc_library),
               PetscErrorCode,
               (CTao, external),
               tao, solve,
              )


	return nothing
end 

"""
	TaoShellGetContext(petsclib::PetscLibType,tao::Tao, ctx::Cvoid) 
Returns the user

Not Collective

Input Parameter:
- `tao` - should have been created with `TaoSetType`(tao,`TAOSHELL`);

Output Parameter:
- `ctx` - the user provided context

Level: advanced

-seealso: `Tao`, `TAOSHELL`, `TaoShellSetContext()`

# External Links
$(_doc_external("Tao/TaoShellGetContext"))
"""
function TaoShellGetContext(petsclib::PetscLibType, tao::Tao, ctx::Cvoid) end

@for_petsc function TaoShellGetContext(petsclib::$UnionPetscLib, tao::Tao, ctx::Cvoid )

    @chk ccall(
               (:TaoShellGetContext, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{Cvoid}),
               tao, ctx,
              )


	return nothing
end 

"""
	TaoShellSetContext(petsclib::PetscLibType,tao::Tao, ctx::Cvoid) 
sets the context for a `TAOSHELL`

Logically Collective

Input Parameters:
- `tao` - the shell Tao
- `ctx` - the context

Level: advanced

-seealso: `Tao`, `TAOSHELL`, `TaoShellGetContext()`

# External Links
$(_doc_external("Tao/TaoShellSetContext"))
"""
function TaoShellSetContext(petsclib::PetscLibType, tao::Tao, ctx::Cvoid) end

@for_petsc function TaoShellSetContext(petsclib::$UnionPetscLib, tao::Tao, ctx::Cvoid )

    @chk ccall(
               (:TaoShellSetContext, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{Cvoid}),
               tao, ctx,
              )


	return nothing
end 

"""
	TaoBRGNGetDampingVector(petsclib::PetscLibType,tao::Tao, d::PetscVec) 
Get the damping vector {diag}(J^T J) from a `TAOBRGN` with `TAOBRGN_REGULARIZATION_LM` regularization

Collective

Input Parameter:
- `tao` - a `Tao` of type `TAOBRGN` with `TAOBRGN_REGULARIZATION_LM` regularization

Output Parameter:
- `d` - the damping vector

Level: developer

-seealso: [](ch_tao), `Tao`, `TAOBRGN`, `TaoBRGNRegularzationTypes`

# External Links
$(_doc_external("Tao/TaoBRGNGetDampingVector"))
"""
function TaoBRGNGetDampingVector(petsclib::PetscLibType, tao::Tao, d::PetscVec) end

@for_petsc function TaoBRGNGetDampingVector(petsclib::$UnionPetscLib, tao::Tao, d::PetscVec )
	d_ = Ref(d.ptr)

    @chk ccall(
               (:TaoBRGNGetDampingVector, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{CVec}),
               tao, d_,
              )

	d.ptr = C_NULL

	return nothing
end 

"""
	type::TaoBRGNRegularizationType = TaoBRGNGetRegularizationType(petsclib::PetscLibType,tao::Tao) 
Get the `TaoBRGNRegularizationType` of a `TAOBRGN`

Not collective

Input Parameter:
- `tao` - a `Tao` of type `TAOBRGN`

Output Parameter:
- `type` - the `TaoBRGNRegularizationType`

Level: advanced

-seealso: [](ch_tao), `Tao`, `TAOBRGN`, `TaoBRGNRegularizationType`, `TaoBRGNSetRegularizationType()`

# External Links
$(_doc_external("Tao/TaoBRGNGetRegularizationType"))
"""
function TaoBRGNGetRegularizationType(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoBRGNGetRegularizationType(petsclib::$UnionPetscLib, tao::Tao )
	type_ = Ref{TaoBRGNRegularizationType}()

    @chk ccall(
               (:TaoBRGNGetRegularizationType, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{TaoBRGNRegularizationType}),
               tao, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	TaoBRGNSetRegularizationType(petsclib::PetscLibType,tao::Tao, type::TaoBRGNRegularizationType) 
Set the `TaoBRGNRegularizationType` of a `TAOBRGN`

Logically collective

Input Parameters:
- `tao`  - a `Tao` of type `TAOBRGN`
- `type` - the `TaoBRGNRegularizationType`

Level: advanced

-seealso: [](ch_tao), `Tao`, `TAOBRGN`, `TaoBRGNRegularizationType`, `TaoBRGNGetRegularizationType`

# External Links
$(_doc_external("Tao/TaoBRGNSetRegularizationType"))
"""
function TaoBRGNSetRegularizationType(petsclib::PetscLibType, tao::Tao, type::TaoBRGNRegularizationType) end

@for_petsc function TaoBRGNSetRegularizationType(petsclib::$UnionPetscLib, tao::Tao, type::TaoBRGNRegularizationType )

    @chk ccall(
               (:TaoBRGNSetRegularizationType, $petsc_library),
               PetscErrorCode,
               (CTao, TaoBRGNRegularizationType),
               tao, type,
              )


	return nothing
end 

"""
	TaoBRGNGetSubsolver(petsclib::PetscLibType,tao::Tao, subsolver::Tao) 
Get the pointer to the subsolver inside a `TAOBRGN`

Collective

Input Parameters:
- `tao`       - the Tao solver context
- `subsolver` - the `Tao` sub-solver context

Level: advanced

-seealso: `Tao`, `Mat`, `TAOBRGN`

# External Links
$(_doc_external("Tao/TaoBRGNGetSubsolver"))
"""
function TaoBRGNGetSubsolver(petsclib::PetscLibType, tao::Tao, subsolver::Tao) end

@for_petsc function TaoBRGNGetSubsolver(petsclib::$UnionPetscLib, tao::Tao, subsolver::Tao )
	subsolver_ = Ref(subsolver.ptr)

    @chk ccall(
               (:TaoBRGNGetSubsolver, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{CTao}),
               tao, subsolver_,
              )

	subsolver.ptr = C_NULL

	return nothing
end 

"""
	TaoBRGNSetRegularizerWeight(petsclib::PetscLibType,tao::Tao, lambda::PetscReal) 
Set the regularizer weight for the Gauss

Collective

Input Parameters:
- `tao`    - the `Tao` solver context
- `lambda` - L1-norm regularizer weight

Level: beginner

-seealso: `Tao`, `Mat`, `TAOBRGN`

# External Links
$(_doc_external("Tao/TaoBRGNSetRegularizerWeight"))
"""
function TaoBRGNSetRegularizerWeight(petsclib::PetscLibType, tao::Tao, lambda::PetscReal) end

@for_petsc function TaoBRGNSetRegularizerWeight(petsclib::$UnionPetscLib, tao::Tao, lambda::$PetscReal )

    @chk ccall(
               (:TaoBRGNSetRegularizerWeight, $petsc_library),
               PetscErrorCode,
               (CTao, $PetscReal),
               tao, lambda,
              )


	return nothing
end 

"""
	TaoBRGNSetL1SmoothEpsilon(petsclib::PetscLibType,tao::Tao, epsilon::PetscReal) 
Set the L1

Collective

Input Parameters:
- `tao`     - the `Tao` solver context
- `epsilon` - L1-norm smooth approximation parameter

Level: advanced

-seealso: `Tao`, `Mat`, `TAOBRGN`

# External Links
$(_doc_external("Tao/TaoBRGNSetL1SmoothEpsilon"))
"""
function TaoBRGNSetL1SmoothEpsilon(petsclib::PetscLibType, tao::Tao, epsilon::PetscReal) end

@for_petsc function TaoBRGNSetL1SmoothEpsilon(petsclib::$UnionPetscLib, tao::Tao, epsilon::$PetscReal )

    @chk ccall(
               (:TaoBRGNSetL1SmoothEpsilon, $petsc_library),
               PetscErrorCode,
               (CTao, $PetscReal),
               tao, epsilon,
              )


	return nothing
end 

"""
	TaoBRGNSetDictionaryMatrix(petsclib::PetscLibType,tao::Tao, dict::PetscMat) 
bind the dictionary matrix from user application context to gn

Input Parameters:
- `tao`  - the `Tao` context
- `dict` - the user specified dictionary matrix.  We allow to set a `NULL` dictionary, which means identity matrix by default

Level: advanced

-seealso: `Tao`, `Mat`, `TAOBRGN`

# External Links
$(_doc_external("Tao/TaoBRGNSetDictionaryMatrix"))
"""
function TaoBRGNSetDictionaryMatrix(petsclib::PetscLibType, tao::Tao, dict::PetscMat) end

@for_petsc function TaoBRGNSetDictionaryMatrix(petsclib::$UnionPetscLib, tao::Tao, dict::PetscMat )

    @chk ccall(
               (:TaoBRGNSetDictionaryMatrix, $petsc_library),
               PetscErrorCode,
               (CTao, CMat),
               tao, dict,
              )


	return nothing
end 

"""
	TaoBRGNSetRegularizerObjectiveAndGradientRoutine(petsclib::PetscLibType,tao::Tao, func::external, ctx::Cvoid) 
Sets the user
function into the algorithm.

Input Parameters:
- `tao`  - the Tao context
- `func` - function pointer for the regularizer value and gradient evaluation
- `ctx`  - user context for the regularizer

Calling sequence:
- `tao` - the `Tao` context
- `u`   - the location at which to compute the objective and gradient
- `val` - location to store objective function value
- `g`   - location to store gradient
- `ctx` - user context for the regularizer Hessian

Level: advanced

-seealso: `Tao`, `Mat`, `TAOBRGN`

# External Links
$(_doc_external("Tao/TaoBRGNSetRegularizerObjectiveAndGradientRoutine"))
"""
function TaoBRGNSetRegularizerObjectiveAndGradientRoutine(petsclib::PetscLibType, tao::Tao, func::external, ctx::Cvoid) end

@for_petsc function TaoBRGNSetRegularizerObjectiveAndGradientRoutine(petsclib::$UnionPetscLib, tao::Tao, func::external, ctx::Cvoid )

    @chk ccall(
               (:TaoBRGNSetRegularizerObjectiveAndGradientRoutine, $petsc_library),
               PetscErrorCode,
               (CTao, external, Ptr{Cvoid}),
               tao, func, ctx,
              )


	return nothing
end 

"""
	TaoBRGNSetRegularizerHessianRoutine(petsclib::PetscLibType,tao::Tao, Hreg::PetscMat, func::external, ctx::Cvoid) 
Sets the user
function into the algorithm.

Input Parameters:
- `tao`  - the `Tao` context
- `Hreg` - user-created matrix for the Hessian of the regularization term
- `func` - function pointer for the regularizer Hessian evaluation
- `ctx`  - user context for the regularizer Hessian

Calling sequence:
- `tao`  - the `Tao` context
- `u`    - the location at which to compute the Hessian
- `Hreg` - user-created matrix for the Hessian of the regularization term
- `ctx`  - user context for the regularizer Hessian

Level: advanced

-seealso: `Tao`, `Mat`, `TAOBRGN`

# External Links
$(_doc_external("Tao/TaoBRGNSetRegularizerHessianRoutine"))
"""
function TaoBRGNSetRegularizerHessianRoutine(petsclib::PetscLibType, tao::Tao, Hreg::PetscMat, func::external, ctx::Cvoid) end

@for_petsc function TaoBRGNSetRegularizerHessianRoutine(petsclib::$UnionPetscLib, tao::Tao, Hreg::PetscMat, func::external, ctx::Cvoid )

    @chk ccall(
               (:TaoBRGNSetRegularizerHessianRoutine, $petsc_library),
               PetscErrorCode,
               (CTao, CMat, external, Ptr{Cvoid}),
               tao, Hreg, func, ctx,
              )


	return nothing
end 

"""
	TaoADMMSetMisfitHessianChangeStatus(petsclib::PetscLibType,tao::Tao, b::PetscBool) 
Set boolean that determines  whether Hessian matrix of misfit subsolver changes with respect to input vector.

Collective

Input Parameters:
- `tao` - the Tao solver context.
- `b`   - the Hessian matrix change status boolean, `PETSC_FALSE`  when the Hessian matrix does not change, `PETSC_TRUE` otherwise.

Level: advanced

-seealso: `TAOADMM`

# External Links
$(_doc_external("Tao/TaoADMMSetMisfitHessianChangeStatus"))
"""
function TaoADMMSetMisfitHessianChangeStatus(petsclib::PetscLibType, tao::Tao, b::PetscBool) end

@for_petsc function TaoADMMSetMisfitHessianChangeStatus(petsclib::$UnionPetscLib, tao::Tao, b::PetscBool )

    @chk ccall(
               (:TaoADMMSetMisfitHessianChangeStatus, $petsc_library),
               PetscErrorCode,
               (CTao, PetscBool),
               tao, b,
              )


	return nothing
end 

"""
	TaoADMMSetRegHessianChangeStatus(petsclib::PetscLibType,tao::Tao, b::PetscBool) 
Set boolean that determines whether Hessian matrix of regularization subsolver changes with respect to input vector.

Collective

Input Parameters:
- `tao` - the `Tao` solver context
- `b`   - the Hessian matrix change status boolean, `PETSC_FALSE` when the Hessian matrix does not change, `PETSC_TRUE` otherwise.

Level: advanced

-seealso: `TAOADMM`

# External Links
$(_doc_external("Tao/TaoADMMSetRegHessianChangeStatus"))
"""
function TaoADMMSetRegHessianChangeStatus(petsclib::PetscLibType, tao::Tao, b::PetscBool) end

@for_petsc function TaoADMMSetRegHessianChangeStatus(petsclib::$UnionPetscLib, tao::Tao, b::PetscBool )

    @chk ccall(
               (:TaoADMMSetRegHessianChangeStatus, $petsc_library),
               PetscErrorCode,
               (CTao, PetscBool),
               tao, b,
              )


	return nothing
end 

"""
	TaoADMMSetSpectralPenalty(petsclib::PetscLibType,tao::Tao, mu::PetscReal) 
Set the spectral penalty (mu) value

Collective

Input Parameters:
- `tao` - the `Tao` solver context
- `mu`  - spectral penalty

Level: advanced

-seealso: `TaoADMMSetMinimumSpectralPenalty()`, `TAOADMM`

# External Links
$(_doc_external("Tao/TaoADMMSetSpectralPenalty"))
"""
function TaoADMMSetSpectralPenalty(petsclib::PetscLibType, tao::Tao, mu::PetscReal) end

@for_petsc function TaoADMMSetSpectralPenalty(petsclib::$UnionPetscLib, tao::Tao, mu::$PetscReal )

    @chk ccall(
               (:TaoADMMSetSpectralPenalty, $petsc_library),
               PetscErrorCode,
               (CTao, $PetscReal),
               tao, mu,
              )


	return nothing
end 

"""
	mu::PetscReal = TaoADMMGetSpectralPenalty(petsclib::PetscLibType,tao::Tao) 
Get the spectral penalty (mu) value

Collective

Input Parameter:
- `tao` - the `Tao` solver context

Output Parameter:
- `mu` - spectral penalty

Level: advanced

-seealso: `TaoADMMSetMinimumSpectralPenalty()`, `TaoADMMSetSpectralPenalty()`, `TAOADMM`

# External Links
$(_doc_external("Tao/TaoADMMGetSpectralPenalty"))
"""
function TaoADMMGetSpectralPenalty(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoADMMGetSpectralPenalty(petsclib::$UnionPetscLib, tao::Tao )
	mu_ = Ref{$PetscReal}()

    @chk ccall(
               (:TaoADMMGetSpectralPenalty, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{$PetscReal}),
               tao, mu_,
              )

	mu = mu_[]

	return mu
end 

"""
	TaoADMMGetMisfitSubsolver(petsclib::PetscLibType,tao::Tao, misfit::Tao) 
Get the pointer to the misfit subsolver inside `TAOADMM`

Collective

Input Parameter:
- `tao` - the `Tao` solver context

Output Parameter:
- `misfit` - the `Tao` subsolver context

Level: advanced

-seealso: `TAOADMM`, `Tao`

# External Links
$(_doc_external("Tao/TaoADMMGetMisfitSubsolver"))
"""
function TaoADMMGetMisfitSubsolver(petsclib::PetscLibType, tao::Tao, misfit::Tao) end

@for_petsc function TaoADMMGetMisfitSubsolver(petsclib::$UnionPetscLib, tao::Tao, misfit::Tao )
	misfit_ = Ref(misfit.ptr)

    @chk ccall(
               (:TaoADMMGetMisfitSubsolver, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{CTao}),
               tao, misfit_,
              )

	misfit.ptr = C_NULL

	return nothing
end 

"""
	TaoADMMGetRegularizationSubsolver(petsclib::PetscLibType,tao::Tao, reg::Tao) 
Get the pointer to the regularization subsolver inside `TAOADMM`

Collective

Input Parameter:
- `tao` - the `Tao` solver context

Output Parameter:
- `reg` - the `Tao` subsolver context

Level: advanced

-seealso: `TAOADMM`, `Tao`

# External Links
$(_doc_external("Tao/TaoADMMGetRegularizationSubsolver"))
"""
function TaoADMMGetRegularizationSubsolver(petsclib::PetscLibType, tao::Tao, reg::Tao) end

@for_petsc function TaoADMMGetRegularizationSubsolver(petsclib::$UnionPetscLib, tao::Tao, reg::Tao )
	reg_ = Ref(reg.ptr)

    @chk ccall(
               (:TaoADMMGetRegularizationSubsolver, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{CTao}),
               tao, reg_,
              )

	reg.ptr = C_NULL

	return nothing
end 

"""
	TaoADMMSetConstraintVectorRHS(petsclib::PetscLibType,tao::Tao, c::PetscVec) 
Set the RHS constraint vector for `TAOADMM`

Collective

Input Parameters:
- `tao` - the `Tao` solver context
- `c`   - RHS vector

Level: advanced

-seealso: `TAOADMM`

# External Links
$(_doc_external("Tao/TaoADMMSetConstraintVectorRHS"))
"""
function TaoADMMSetConstraintVectorRHS(petsclib::PetscLibType, tao::Tao, c::PetscVec) end

@for_petsc function TaoADMMSetConstraintVectorRHS(petsclib::$UnionPetscLib, tao::Tao, c::PetscVec )

    @chk ccall(
               (:TaoADMMSetConstraintVectorRHS, $petsc_library),
               PetscErrorCode,
               (CTao, CVec),
               tao, c,
              )


	return nothing
end 

"""
	TaoADMMSetMinimumSpectralPenalty(petsclib::PetscLibType,tao::Tao, mu::PetscReal) 
Set the minimum value for the spectral penalty

Collective

Input Parameters:
- `tao` - the `Tao` solver context
- `mu`  - minimum spectral penalty value

Level: advanced

-seealso: `TaoADMMGetSpectralPenalty()`, `TAOADMM`

# External Links
$(_doc_external("Tao/TaoADMMSetMinimumSpectralPenalty"))
"""
function TaoADMMSetMinimumSpectralPenalty(petsclib::PetscLibType, tao::Tao, mu::PetscReal) end

@for_petsc function TaoADMMSetMinimumSpectralPenalty(petsclib::$UnionPetscLib, tao::Tao, mu::$PetscReal )

    @chk ccall(
               (:TaoADMMSetMinimumSpectralPenalty, $petsc_library),
               PetscErrorCode,
               (CTao, $PetscReal),
               tao, mu,
              )


	return nothing
end 

"""
	TaoADMMSetRegularizerCoefficient(petsclib::PetscLibType,tao::Tao, lambda::PetscReal) 
Set the regularization coefficient lambda for L1 norm regularization case

Collective

Input Parameters:
- `tao`    - the `Tao` solver context
- `lambda` - L1-norm regularizer coefficient

Level: advanced

-seealso: `TaoADMMSetMisfitConstraintJacobian()`, `TaoADMMSetRegularizerConstraintJacobian()`, `TAOADMM`

# External Links
$(_doc_external("Tao/TaoADMMSetRegularizerCoefficient"))
"""
function TaoADMMSetRegularizerCoefficient(petsclib::PetscLibType, tao::Tao, lambda::PetscReal) end

@for_petsc function TaoADMMSetRegularizerCoefficient(petsclib::$UnionPetscLib, tao::Tao, lambda::$PetscReal )

    @chk ccall(
               (:TaoADMMSetRegularizerCoefficient, $petsc_library),
               PetscErrorCode,
               (CTao, $PetscReal),
               tao, lambda,
              )


	return nothing
end 

"""
	lambda::PetscReal = TaoADMMGetRegularizerCoefficient(petsclib::PetscLibType,tao::Tao) 
Get the regularization coefficient lambda for L1 norm regularization case

Collective

Input Parameter:
- `tao` - the `Tao` solver context

Output Parameter:
- `lambda` - L1-norm regularizer coefficient

Level: advanced

-seealso: `TaoADMMSetMisfitConstraintJacobian()`, `TaoADMMSetRegularizerConstraintJacobian()`, `TAOADMM`

# External Links
$(_doc_external("Tao/TaoADMMGetRegularizerCoefficient"))
"""
function TaoADMMGetRegularizerCoefficient(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoADMMGetRegularizerCoefficient(petsclib::$UnionPetscLib, tao::Tao )
	lambda_ = Ref{$PetscReal}()

    @chk ccall(
               (:TaoADMMGetRegularizerCoefficient, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{$PetscReal}),
               tao, lambda_,
              )

	lambda = lambda_[]

	return lambda
end 

"""
	TaoADMMSetMisfitConstraintJacobian(petsclib::PetscLibType,tao::Tao, J::PetscMat, Jpre::PetscMat, func::external, ctx::Cvoid) 
Set the constraint matrix B for the `TAOADMM` algorithm. Matrix B constrains the z variable.

Collective

Input Parameters:
- `tao`  - the Tao solver context
- `J`    - user-created regularizer constraint Jacobian matrix
- `Jpre` - user-created regularizer Jacobian constraint matrix for constructing the preconditioner, often this is `J`
- `func` - function pointer for the regularizer constraint Jacobian update function
- `ctx`  - user context for the regularizer Hessian

Level: advanced

-seealso: `TaoADMMSetRegularizerCoefficient()`, `TaoADMMSetRegularizerConstraintJacobian()`, `TAOADMM`

# External Links
$(_doc_external("Tao/TaoADMMSetMisfitConstraintJacobian"))
"""
function TaoADMMSetMisfitConstraintJacobian(petsclib::PetscLibType, tao::Tao, J::PetscMat, Jpre::PetscMat, func::external, ctx::Cvoid) end

@for_petsc function TaoADMMSetMisfitConstraintJacobian(petsclib::$UnionPetscLib, tao::Tao, J::PetscMat, Jpre::PetscMat, func::external, ctx::Cvoid )

    @chk ccall(
               (:TaoADMMSetMisfitConstraintJacobian, $petsc_library),
               PetscErrorCode,
               (CTao, CMat, CMat, external, Ptr{Cvoid}),
               tao, J, Jpre, func, ctx,
              )


	return nothing
end 

"""
	TaoADMMSetRegularizerConstraintJacobian(petsclib::PetscLibType,tao::Tao, J::PetscMat, Jpre::PetscMat, func::external, ctx::Cvoid) 
Set the constraint matrix B for `TAOADMM` algorithm. Matrix B constraints z variable.

Collective

Input Parameters:
- `tao`  - the `Tao` solver context
- `J`    - user-created regularizer constraint Jacobian matrix
- `Jpre` - user-created regularizer Jacobian constraint matrix for constructing the preconditioner, often this is `J`
- `func` - function pointer for the regularizer constraint Jacobian update function
- `ctx`  - user context for the regularizer Hessian

Level: advanced

-seealso: `TaoADMMSetRegularizerCoefficient()`, `TaoADMMSetMisfitConstraintJacobian()`, `TAOADMM`

# External Links
$(_doc_external("Tao/TaoADMMSetRegularizerConstraintJacobian"))
"""
function TaoADMMSetRegularizerConstraintJacobian(petsclib::PetscLibType, tao::Tao, J::PetscMat, Jpre::PetscMat, func::external, ctx::Cvoid) end

@for_petsc function TaoADMMSetRegularizerConstraintJacobian(petsclib::$UnionPetscLib, tao::Tao, J::PetscMat, Jpre::PetscMat, func::external, ctx::Cvoid )

    @chk ccall(
               (:TaoADMMSetRegularizerConstraintJacobian, $petsc_library),
               PetscErrorCode,
               (CTao, CMat, CMat, external, Ptr{Cvoid}),
               tao, J, Jpre, func, ctx,
              )


	return nothing
end 

"""
	TaoADMMSetMisfitObjectiveAndGradientRoutine(petsclib::PetscLibType,tao::Tao, func::external, ctx::Cvoid) 
Sets the user

Collective

Input Parameters:
- `tao`  - the `Tao` context
- `func` - function pointer for the misfit value and gradient evaluation
- `ctx`  - user context for the misfit

Level: advanced

-seealso: `TAOADMM`

# External Links
$(_doc_external("Tao/TaoADMMSetMisfitObjectiveAndGradientRoutine"))
"""
function TaoADMMSetMisfitObjectiveAndGradientRoutine(petsclib::PetscLibType, tao::Tao, func::external, ctx::Cvoid) end

@for_petsc function TaoADMMSetMisfitObjectiveAndGradientRoutine(petsclib::$UnionPetscLib, tao::Tao, func::external, ctx::Cvoid )

    @chk ccall(
               (:TaoADMMSetMisfitObjectiveAndGradientRoutine, $petsc_library),
               PetscErrorCode,
               (CTao, external, Ptr{Cvoid}),
               tao, func, ctx,
              )


	return nothing
end 

"""
	TaoADMMSetMisfitHessianRoutine(petsclib::PetscLibType,tao::Tao, H::PetscMat, Hpre::PetscMat, func::external, ctx::Cvoid) 
Sets the user
function into the algorithm, to be used for subsolverX.

Collective

Input Parameters:
- `tao`  - the `Tao` context
- `H`    - user-created matrix for the Hessian of the misfit term
- `Hpre` - user-created matrix for the preconditioner of Hessian of the misfit term
- `func` - function pointer for the misfit Hessian evaluation
- `ctx`  - user context for the misfit Hessian

Level: advanced

-seealso: `TAOADMM`

# External Links
$(_doc_external("Tao/TaoADMMSetMisfitHessianRoutine"))
"""
function TaoADMMSetMisfitHessianRoutine(petsclib::PetscLibType, tao::Tao, H::PetscMat, Hpre::PetscMat, func::external, ctx::Cvoid) end

@for_petsc function TaoADMMSetMisfitHessianRoutine(petsclib::$UnionPetscLib, tao::Tao, H::PetscMat, Hpre::PetscMat, func::external, ctx::Cvoid )

    @chk ccall(
               (:TaoADMMSetMisfitHessianRoutine, $petsc_library),
               PetscErrorCode,
               (CTao, CMat, CMat, external, Ptr{Cvoid}),
               tao, H, Hpre, func, ctx,
              )


	return nothing
end 

"""
	TaoADMMSetRegularizerObjectiveAndGradientRoutine(petsclib::PetscLibType,tao::Tao, func::external, ctx::Cvoid) 
Sets the user

Collective

Input Parameters:
- `tao`  - the Tao context
- `func` - function pointer for the regularizer value and gradient evaluation
- `ctx`  - user context for the regularizer

Level: advanced

-seealso: `TAOADMM`

# External Links
$(_doc_external("Tao/TaoADMMSetRegularizerObjectiveAndGradientRoutine"))
"""
function TaoADMMSetRegularizerObjectiveAndGradientRoutine(petsclib::PetscLibType, tao::Tao, func::external, ctx::Cvoid) end

@for_petsc function TaoADMMSetRegularizerObjectiveAndGradientRoutine(petsclib::$UnionPetscLib, tao::Tao, func::external, ctx::Cvoid )

    @chk ccall(
               (:TaoADMMSetRegularizerObjectiveAndGradientRoutine, $petsc_library),
               PetscErrorCode,
               (CTao, external, Ptr{Cvoid}),
               tao, func, ctx,
              )


	return nothing
end 

"""
	TaoADMMSetRegularizerHessianRoutine(petsclib::PetscLibType,tao::Tao, H::PetscMat, Hpre::PetscMat, func::external, ctx::Cvoid) 
Sets the user
function, to be used for subsolverZ.

Collective

Input Parameters:
- `tao`  - the `Tao` context
- `H`    - user-created matrix for the Hessian of the regularization term
- `Hpre` - user-created matrix for the preconditioner of Hessian of the regularization term
- `func` - function pointer for the regularizer Hessian evaluation
- `ctx`  - user context for the regularizer Hessian

Level: advanced

-seealso: `TAOADMM`

# External Links
$(_doc_external("Tao/TaoADMMSetRegularizerHessianRoutine"))
"""
function TaoADMMSetRegularizerHessianRoutine(petsclib::PetscLibType, tao::Tao, H::PetscMat, Hpre::PetscMat, func::external, ctx::Cvoid) end

@for_petsc function TaoADMMSetRegularizerHessianRoutine(petsclib::$UnionPetscLib, tao::Tao, H::PetscMat, Hpre::PetscMat, func::external, ctx::Cvoid )

    @chk ccall(
               (:TaoADMMSetRegularizerHessianRoutine, $petsc_library),
               PetscErrorCode,
               (CTao, CMat, CMat, external, Ptr{Cvoid}),
               tao, H, Hpre, func, ctx,
              )


	return nothing
end 

"""
	TaoGetADMMParentTao(petsclib::PetscLibType,tao::Tao, admm_tao::Tao) 
Gets pointer to parent `TAOADMM`, used by inner subsolver.

Collective

Input Parameter:
- `tao` - the `Tao` context

Output Parameter:
- `admm_tao` - the parent `Tao` context

Level: advanced

-seealso: `TAOADMM`

# External Links
$(_doc_external("Tao/TaoGetADMMParentTao"))
"""
function TaoGetADMMParentTao(petsclib::PetscLibType, tao::Tao, admm_tao::Tao) end

@for_petsc function TaoGetADMMParentTao(petsclib::$UnionPetscLib, tao::Tao, admm_tao::Tao )
	admm_tao_ = Ref(admm_tao.ptr)

    @chk ccall(
               (:TaoGetADMMParentTao, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{CTao}),
               tao, admm_tao_,
              )

	admm_tao.ptr = C_NULL

	return nothing
end 

"""
	TaoADMMGetDualVector(petsclib::PetscLibType,tao::Tao, Y::PetscVec) 
Returns the dual vector associated with the current `TAOADMM` state

Not Collective

Input Parameter:
- `tao` - the `Tao` context

Output Parameter:
- `Y` - the current solution

Level: intermediate

-seealso: `TAOADMM`

# External Links
$(_doc_external("Tao/TaoADMMGetDualVector"))
"""
function TaoADMMGetDualVector(petsclib::PetscLibType, tao::Tao, Y::PetscVec) end

@for_petsc function TaoADMMGetDualVector(petsclib::$UnionPetscLib, tao::Tao, Y::PetscVec )
	Y_ = Ref(Y.ptr)

    @chk ccall(
               (:TaoADMMGetDualVector, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{CVec}),
               tao, Y_,
              )

	Y.ptr = C_NULL

	return nothing
end 

"""
	TaoADMMSetRegularizerType(petsclib::PetscLibType,tao::Tao, type::TaoADMMRegularizerType) 
Set regularizer type for `TAOADMM` routine

Not Collective

Input Parameters:
- `tao`  - the `Tao` context
- `type` - regularizer type

Options Database Key:
- `-tao_admm_regularizer_type <admm_regularizer_user,admm_regularizer_soft_thresh>` - select the regularizer

Level: intermediate

-seealso: `TaoADMMGetRegularizerType()`, `TaoADMMRegularizerType`, `TAOADMM`

# External Links
$(_doc_external("Tao/TaoADMMSetRegularizerType"))
"""
function TaoADMMSetRegularizerType(petsclib::PetscLibType, tao::Tao, type::TaoADMMRegularizerType) end

@for_petsc function TaoADMMSetRegularizerType(petsclib::$UnionPetscLib, tao::Tao, type::TaoADMMRegularizerType )

    @chk ccall(
               (:TaoADMMSetRegularizerType, $petsc_library),
               PetscErrorCode,
               (CTao, TaoADMMRegularizerType),
               tao, type,
              )


	return nothing
end 

"""
	type::TaoADMMRegularizerType = TaoADMMGetRegularizerType(petsclib::PetscLibType,tao::Tao) 
Gets the type of regularizer routine for `TAOADMM`

Not Collective

Input Parameter:
- `tao` - the `Tao` context

Output Parameter:
- `type` - the type of regularizer

Level: intermediate

-seealso: `TaoADMMSetRegularizerType()`, `TaoADMMRegularizerType`, `TAOADMM`

# External Links
$(_doc_external("Tao/TaoADMMGetRegularizerType"))
"""
function TaoADMMGetRegularizerType(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoADMMGetRegularizerType(petsclib::$UnionPetscLib, tao::Tao )
	type_ = Ref{TaoADMMRegularizerType}()

    @chk ccall(
               (:TaoADMMGetRegularizerType, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{TaoADMMRegularizerType}),
               tao, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	TaoADMMSetUpdateType(petsclib::PetscLibType,tao::Tao, type::TaoADMMUpdateType) 
Set update routine for `TAOADMM` routine

Not Collective

Input Parameters:
- `tao`  - the `Tao` context
- `type` - spectral parameter update type

Level: intermediate

-seealso: `TaoADMMGetUpdateType()`, `TaoADMMUpdateType`, `TAOADMM`

# External Links
$(_doc_external("Tao/TaoADMMSetUpdateType"))
"""
function TaoADMMSetUpdateType(petsclib::PetscLibType, tao::Tao, type::TaoADMMUpdateType) end

@for_petsc function TaoADMMSetUpdateType(petsclib::$UnionPetscLib, tao::Tao, type::TaoADMMUpdateType )

    @chk ccall(
               (:TaoADMMSetUpdateType, $petsc_library),
               PetscErrorCode,
               (CTao, TaoADMMUpdateType),
               tao, type,
              )


	return nothing
end 

"""
	type::TaoADMMUpdateType = TaoADMMGetUpdateType(petsclib::PetscLibType,tao::Tao) 
Gets the type of spectral penalty update routine for `TAOADMM`

Not Collective

Input Parameter:
- `tao` - the `Tao` context

Output Parameter:
- `type` - the type of spectral penalty update routine

Level: intermediate

-seealso: `TaoADMMSetUpdateType()`, `TaoADMMUpdateType`, `TAOADMM`

# External Links
$(_doc_external("Tao/TaoADMMGetUpdateType"))
"""
function TaoADMMGetUpdateType(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoADMMGetUpdateType(petsclib::$UnionPetscLib, tao::Tao )
	type_ = Ref{TaoADMMUpdateType}()

    @chk ccall(
               (:TaoADMMGetUpdateType, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{TaoADMMUpdateType}),
               tao, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	type::TaoALMMType = TaoALMMGetType(petsclib::PetscLibType,tao::Tao) 
Retrieve the augmented Lagrangian formulation type for the subproblem.

Input Parameter:
- `tao` - the `Tao` context for the `TAOALMM` solver

Output Parameter:
- `type` - augmented Lagragrangian type

Level: advanced

-seealso: `Tao`, `TAOALMM`, `TaoALMMSetType()`, `TaoALMMType`

# External Links
$(_doc_external("Tao/TaoALMMGetType"))
"""
function TaoALMMGetType(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoALMMGetType(petsclib::$UnionPetscLib, tao::Tao )
	type_ = Ref{TaoALMMType}()

    @chk ccall(
               (:TaoALMMGetType, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{TaoALMMType}),
               tao, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	TaoALMMSetType(petsclib::PetscLibType,tao::Tao, type::TaoALMMType) 
Determine the augmented Lagrangian formulation type for the subproblem.

Input Parameters:
- `tao`  - the `Tao` context for the `TAOALMM` solver
- `type` - augmented Lagragrangian type

Level: advanced

-seealso: `Tao`, `TAOALMM`, `TaoALMMGetType()`, `TaoALMMType`

# External Links
$(_doc_external("Tao/TaoALMMSetType"))
"""
function TaoALMMSetType(petsclib::PetscLibType, tao::Tao, type::TaoALMMType) end

@for_petsc function TaoALMMSetType(petsclib::$UnionPetscLib, tao::Tao, type::TaoALMMType )

    @chk ccall(
               (:TaoALMMSetType, $petsc_library),
               PetscErrorCode,
               (CTao, TaoALMMType),
               tao, type,
              )


	return nothing
end 

"""
	TaoALMMGetSubsolver(petsclib::PetscLibType,tao::Tao, subsolver::Tao) 
Retrieve the subsolver being used by `TAOALMM`.

Input Parameter:
- `tao` - the `Tao` context for the `TAOALMM` solver

Output Parameter:
- `subsolver` - the `Tao` context for the subsolver

Level: advanced

-seealso: `Tao`, `TAOALMM`, `TaoALMMSetSubsolver()`

# External Links
$(_doc_external("Tao/TaoALMMGetSubsolver"))
"""
function TaoALMMGetSubsolver(petsclib::PetscLibType, tao::Tao, subsolver::Tao) end

@for_petsc function TaoALMMGetSubsolver(petsclib::$UnionPetscLib, tao::Tao, subsolver::Tao )
	subsolver_ = Ref(subsolver.ptr)

    @chk ccall(
               (:TaoALMMGetSubsolver, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{CTao}),
               tao, subsolver_,
              )

	subsolver.ptr = C_NULL

	return nothing
end 

"""
	TaoALMMSetSubsolver(petsclib::PetscLibType,tao::Tao, subsolver::Tao) 
Changes the subsolver inside `TAOALMM` with the user provided one.

Input Parameters:
- `tao`       - the `Tao` context for the `TAOALMM` solver
- `subsolver` - the Tao context for the subsolver

Level: advanced

-seealso: `Tao`, `TAOALMM`, `TaoALMMGetSubsolver()`

# External Links
$(_doc_external("Tao/TaoALMMSetSubsolver"))
"""
function TaoALMMSetSubsolver(petsclib::PetscLibType, tao::Tao, subsolver::Tao) end

@for_petsc function TaoALMMSetSubsolver(petsclib::$UnionPetscLib, tao::Tao, subsolver::Tao )

    @chk ccall(
               (:TaoALMMSetSubsolver, $petsc_library),
               PetscErrorCode,
               (CTao, CTao),
               tao, subsolver,
              )


	return nothing
end 

"""
	TaoALMMGetMultipliers(petsclib::PetscLibType,tao::Tao, Y::PetscVec) 
Retrieve a pointer to the Lagrange multipliers.

Input Parameter:
- `tao` - the `Tao` context for the `TAOALMM` solver

Output Parameter:
- `Y` - vector of Lagrange multipliers

Level: advanced

-seealso: `TAOALMM`, `Tao`, `TaoALMMSetMultipliers()`, `TaoALMMGetDualIS()`

# External Links
$(_doc_external("Tao/TaoALMMGetMultipliers"))
"""
function TaoALMMGetMultipliers(petsclib::PetscLibType, tao::Tao, Y::PetscVec) end

@for_petsc function TaoALMMGetMultipliers(petsclib::$UnionPetscLib, tao::Tao, Y::PetscVec )
	Y_ = Ref(Y.ptr)

    @chk ccall(
               (:TaoALMMGetMultipliers, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{CVec}),
               tao, Y_,
              )

	Y.ptr = C_NULL

	return nothing
end 

"""
	TaoALMMSetMultipliers(petsclib::PetscLibType,tao::Tao, Y::PetscVec) 
Set user

Input Parameters:
- `tao` - the `Tao` context for the `TAOALMM` solver
- `Y`   - vector of Lagrange multipliers

Level: advanced

-seealso: `TAOALMM`, `Tao`, `TaoALMMGetMultipliers()`

# External Links
$(_doc_external("Tao/TaoALMMSetMultipliers"))
"""
function TaoALMMSetMultipliers(petsclib::PetscLibType, tao::Tao, Y::PetscVec) end

@for_petsc function TaoALMMSetMultipliers(petsclib::$UnionPetscLib, tao::Tao, Y::PetscVec )

    @chk ccall(
               (:TaoALMMSetMultipliers, $petsc_library),
               PetscErrorCode,
               (CTao, CVec),
               tao, Y,
              )


	return nothing
end 

"""
	TaoALMMGetPrimalIS(petsclib::PetscLibType,tao::Tao, opt_is::IS, slack_is::IS) 
Retrieve the index set that identifies optimization
and slack variable components of the subsolver's solution vector.

Input Parameter:
- `tao` - the `Tao` context for the `TAOALMM` solver

Output Parameters:
- `opt_is`   - index set associated with the optimization variables (`NULL` if not needed)
- `slack_is` - index set associated with the slack variables (`NULL` if not needed)

Level: advanced

-seealso: `TAOALMM`, `Tao`, `IS`, `TaoALMMGetPrimalVector()`

# External Links
$(_doc_external("Tao/TaoALMMGetPrimalIS"))
"""
function TaoALMMGetPrimalIS(petsclib::PetscLibType, tao::Tao, opt_is::IS, slack_is::IS) end

@for_petsc function TaoALMMGetPrimalIS(petsclib::$UnionPetscLib, tao::Tao, opt_is::IS, slack_is::IS )
	opt_is_ = Ref(opt_is.ptr)
	slack_is_ = Ref(slack_is.ptr)

    @chk ccall(
               (:TaoALMMGetPrimalIS, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{CIS}, Ptr{CIS}),
               tao, opt_is_, slack_is_,
              )

	opt_is.ptr = C_NULL
	slack_is.ptr = C_NULL

	return nothing
end 

"""
	TaoALMMGetDualIS(petsclib::PetscLibType,tao::Tao, eq_is::IS, ineq_is::IS) 
Retrieve the index set that identifies equality
and inequality constraint components of the dual vector returned
by `TaoALMMGetMultipliers()`.

Input Parameter:
- `tao` - the Tao context for the `TAOALMM` solver

Output Parameters:
- `eq_is`   - index set associated with the equality constraints (`NULL` if not needed)
- `ineq_is` - index set associated with the inequality constraints (`NULL` if not needed)

Level: advanced

-seealso: `TAOALMM`, `Tao`, `TaoALMMGetMultipliers()`

# External Links
$(_doc_external("Tao/TaoALMMGetDualIS"))
"""
function TaoALMMGetDualIS(petsclib::PetscLibType, tao::Tao, eq_is::IS, ineq_is::IS) end

@for_petsc function TaoALMMGetDualIS(petsclib::$UnionPetscLib, tao::Tao, eq_is::IS, ineq_is::IS )
	eq_is_ = Ref(eq_is.ptr)
	ineq_is_ = Ref(ineq_is.ptr)

    @chk ccall(
               (:TaoALMMGetDualIS, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{CIS}, Ptr{CIS}),
               tao, eq_is_, ineq_is_,
              )

	eq_is.ptr = C_NULL
	ineq_is.ptr = C_NULL

	return nothing
end 

"""
	TaoLMVMRecycle(petsclib::PetscLibType,tao::Tao, flg::PetscBool) 
Enable/disable recycling of the QN history between subsequent `TaoSolve()` calls.

Input Parameters:
- `tao` - the `Tao` solver context
- `flg` - Boolean flag for recycling (`PETSC_TRUE` or `PETSC_FALSE`)

Level: intermediate

-seealso: `Tao`, `TAOLMVM`, `TAOBLMVM`

# External Links
$(_doc_external("Tao/TaoLMVMRecycle"))
"""
function TaoLMVMRecycle(petsclib::PetscLibType, tao::Tao, flg::PetscBool) end

@for_petsc function TaoLMVMRecycle(petsclib::$UnionPetscLib, tao::Tao, flg::PetscBool )

    @chk ccall(
               (:TaoLMVMRecycle, $petsc_library),
               PetscErrorCode,
               (CTao, PetscBool),
               tao, flg,
              )


	return nothing
end 

"""
	TaoLMVMSetH0(petsclib::PetscLibType,tao::Tao, H0::PetscMat) 
Set the initial Hessian for the QN approximation

Input Parameters:
- `tao` - the `Tao` solver context
- `H0`  - `Mat` object for the initial Hessian

Level: advanced

-seealso: `Tao`, `TAOLMVM`, `TAOBLMVM`, `TaoLMVMGetH0()`, `TaoLMVMGetH0KSP()`

# External Links
$(_doc_external("Tao/TaoLMVMSetH0"))
"""
function TaoLMVMSetH0(petsclib::PetscLibType, tao::Tao, H0::PetscMat) end

@for_petsc function TaoLMVMSetH0(petsclib::$UnionPetscLib, tao::Tao, H0::PetscMat )

    @chk ccall(
               (:TaoLMVMSetH0, $petsc_library),
               PetscErrorCode,
               (CTao, CMat),
               tao, H0,
              )


	return nothing
end 

"""
	TaoLMVMGetH0(petsclib::PetscLibType,tao::Tao, H0::PetscMat) 
Get the matrix object for the QN initial Hessian

Input Parameter:
- `tao` - the `Tao` solver context

Output Parameter:
- `H0` - `Mat` object for the initial Hessian

Level: advanced

-seealso: `Tao`, `TAOLMVM`, `TAOBLMVM`, `TaoLMVMSetH0()`, `TaoLMVMGetH0KSP()`

# External Links
$(_doc_external("Tao/TaoLMVMGetH0"))
"""
function TaoLMVMGetH0(petsclib::PetscLibType, tao::Tao, H0::PetscMat) end

@for_petsc function TaoLMVMGetH0(petsclib::$UnionPetscLib, tao::Tao, H0::PetscMat )
	H0_ = Ref(H0.ptr)

    @chk ccall(
               (:TaoLMVMGetH0, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{CMat}),
               tao, H0_,
              )

	H0.ptr = C_NULL

	return nothing
end 

"""
	TaoLMVMGetH0KSP(petsclib::PetscLibType,tao::Tao, ksp::PetscKSP) 
Get the iterative solver for applying the inverse of the QN initial Hessian

Input Parameter:
- `tao` - the `Tao` solver context

Output Parameter:
- `ksp` - `KSP` solver context for the initial Hessian

Level: advanced

-seealso: `Tao`, `TAOLMVM`, `TAOBLMVM`, `TaoLMVMGetH0()`

# External Links
$(_doc_external("Tao/TaoLMVMGetH0KSP"))
"""
function TaoLMVMGetH0KSP(petsclib::PetscLibType, tao::Tao, ksp::PetscKSP) end

@for_petsc function TaoLMVMGetH0KSP(petsclib::$UnionPetscLib, tao::Tao, ksp::PetscKSP )
	ksp_ = Ref(ksp.ptr)

    @chk ccall(
               (:TaoLMVMGetH0KSP, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{CKSP}),
               tao, ksp_,
              )

	ksp.ptr = C_NULL

	return nothing
end 

"""
	type::TaoBNCGType = TaoBNCGGetType(petsclib::PetscLibType,tao::Tao) 
Return the type for the `TAOBNCG` solver

Input Parameter:
- `tao` - the `Tao` solver context

Output Parameter:
- `type` - `TAOBNCG` type

Level: advanced

-seealso: `Tao`, `TAOBNCG`, `TaoBNCGSetType()`, `TaoBNCGType`

# External Links
$(_doc_external("Tao/TaoBNCGGetType"))
"""
function TaoBNCGGetType(petsclib::PetscLibType, tao::Tao) end

@for_petsc function TaoBNCGGetType(petsclib::$UnionPetscLib, tao::Tao )
	type_ = Ref{TaoBNCGType}()

    @chk ccall(
               (:TaoBNCGGetType, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{TaoBNCGType}),
               tao, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	TaoBNCGSetType(petsclib::PetscLibType,tao::Tao, type::TaoBNCGType) 
Set the type for the `TAOBNCG` solver

Input Parameters:
- `tao`  - the `Tao` solver context
- `type` - `TAOBNCG` type

Level: advanced

-seealso: `Tao`, `TAOBNCG`, `TaoBNCGGetType()`, `TaoBNCGType`

# External Links
$(_doc_external("Tao/TaoBNCGSetType"))
"""
function TaoBNCGSetType(petsclib::PetscLibType, tao::Tao, type::TaoBNCGType) end

@for_petsc function TaoBNCGSetType(petsclib::$UnionPetscLib, tao::Tao, type::TaoBNCGType )

    @chk ccall(
               (:TaoBNCGSetType, $petsc_library),
               PetscErrorCode,
               (CTao, TaoBNCGType),
               tao, type,
              )


	return nothing
end 

"""
	TaoGetLMVMMatrix(petsclib::PetscLibType,tao::Tao, B::PetscMat) 
Returns a pointer to the internal LMVM matrix. Valid
only for quasi-Newton family of methods.

Input Parameter:
- `tao` - `Tao` solver context

Output Parameter:
- `B` - LMVM matrix

Level: advanced

-seealso: `TAOBQNLS`, `TAOBQNKLS`, `TAOBQNKTL`, `TAOBQNKTR`, `MATLMVM`, `TaoSetLMVMMatrix()`

# External Links
$(_doc_external("Tao/TaoGetLMVMMatrix"))
"""
function TaoGetLMVMMatrix(petsclib::PetscLibType, tao::Tao, B::PetscMat) end

@for_petsc function TaoGetLMVMMatrix(petsclib::$UnionPetscLib, tao::Tao, B::PetscMat )
	B_ = Ref(B.ptr)

    @chk ccall(
               (:TaoGetLMVMMatrix, $petsc_library),
               PetscErrorCode,
               (CTao, Ptr{CMat}),
               tao, B_,
              )

	B.ptr = C_NULL

	return nothing
end 

"""
	TaoSetLMVMMatrix(petsclib::PetscLibType,tao::Tao, B::PetscMat) 
Sets an external LMVM matrix into the Tao solver. Valid
only for quasi-Newton family of methods.

QN family of methods create their own LMVM matrices and users who wish to
manipulate this matrix should use TaoGetLMVMMatrix() instead.

Input Parameters:
- `tao` - Tao solver context
- `B`   - LMVM matrix

Level: advanced

-seealso: `TAOBQNLS`, `TAOBQNKLS`, `TAOBQNKTL`, `TAOBQNKTR`, `MATLMVM`, `TaoGetLMVMMatrix()`

# External Links
$(_doc_external("Tao/TaoSetLMVMMatrix"))
"""
function TaoSetLMVMMatrix(petsclib::PetscLibType, tao::Tao, B::PetscMat) end

@for_petsc function TaoSetLMVMMatrix(petsclib::$UnionPetscLib, tao::Tao, B::PetscMat )

    @chk ccall(
               (:TaoSetLMVMMatrix, $petsc_library),
               PetscErrorCode,
               (CTao, CMat),
               tao, B,
              )


	return nothing
end 

"""
	TaoVecGetSubVec(petsclib::PetscLibType,vfull::PetscVec, is::IS, reduced_type::TaoSubsetType, maskvalue::PetscReal, vreduced::PetscVec) 
Gets a subvector using the `IS`

Input Parameters:
- `vfull`        - the full matrix
- `is`           - the index set for the subvector
- `reduced_type` - the method `Tao` is using for subsetting
- `maskvalue`    - the value to set the unused vector elements to (for `TAO_SUBSET_MASK` or `TAO_SUBSET_MATRIXFREE`)

Output Parameter:
- `vreduced` - the subvector

Level: developer

-seealso: `TaoMatGetSubMat()`, `TaoSubsetType`

# External Links
$(_doc_external("Tao/TaoVecGetSubVec"))
"""
function TaoVecGetSubVec(petsclib::PetscLibType, vfull::PetscVec, is::IS, reduced_type::TaoSubsetType, maskvalue::PetscReal, vreduced::PetscVec) end

@for_petsc function TaoVecGetSubVec(petsclib::$UnionPetscLib, vfull::PetscVec, is::IS, reduced_type::TaoSubsetType, maskvalue::$PetscReal, vreduced::PetscVec )
	vreduced_ = Ref(vreduced.ptr)

    @chk ccall(
               (:TaoVecGetSubVec, $petsc_library),
               PetscErrorCode,
               (CVec, CIS, TaoSubsetType, $PetscReal, Ptr{CVec}),
               vfull, is, reduced_type, maskvalue, vreduced_,
              )

	vreduced.ptr = C_NULL

	return nothing
end 

"""
	TaoMatGetSubMat(petsclib::PetscLibType,M::PetscMat, is::IS, v1::PetscVec, subset_type::TaoSubsetType, Msub::PetscMat) 
Gets a submatrix using the `IS`

Input Parameters:
- `M`           - the full matrix (`n x n`)
- `is`          - the index set for the submatrix (both row and column index sets need to be the same)
- `v1`          - work vector of dimension n, needed for `TAO_SUBSET_MASK` option
- `subset_type` - the method `Tao` is using for subsetting

Output Parameter:
- `Msub` - the submatrix

Level: developer

-seealso: `TaoVecGetSubVec()`, `TaoSubsetType`

# External Links
$(_doc_external("Tao/TaoMatGetSubMat"))
"""
function TaoMatGetSubMat(petsclib::PetscLibType, M::PetscMat, is::IS, v1::PetscVec, subset_type::TaoSubsetType, Msub::PetscMat) end

@for_petsc function TaoMatGetSubMat(petsclib::$UnionPetscLib, M::PetscMat, is::IS, v1::PetscVec, subset_type::TaoSubsetType, Msub::PetscMat )
	Msub_ = Ref(Msub.ptr)

    @chk ccall(
               (:TaoMatGetSubMat, $petsc_library),
               PetscErrorCode,
               (CMat, CIS, CVec, TaoSubsetType, Ptr{CMat}),
               M, is, v1, subset_type, Msub_,
              )

	Msub.ptr = C_NULL

	return nothing
end 

"""
	bound_tol::PetscReal = TaoEstimateActiveBounds(petsclib::PetscLibType,X::PetscVec, XL::PetscVec, XU::PetscVec, G::PetscVec, S::PetscVec, W::PetscVec, steplen::PetscReal, active_lower::IS, active_upper::IS, active_fixed::IS, active::IS, inactive::IS) 
Generates index sets for variables at the lower and upper
bounds, as well as fixed variables where lower and upper bounds equal each other.

Input Parameters:
- `X`       - solution vector
- `XL`      - lower bound vector
- `XU`      - upper bound vector
- `G`       - unprojected gradient
- `S`       - step direction with which the active bounds will be estimated
- `W`       - work vector of type and size of `X`
- `steplen` - the step length at which the active bounds will be estimated (needs to be conservative)

Output Parameters:
- `bound_tol`    - tolerance for the bound estimation
- `active_lower` - index set for active variables at the lower bound
- `active_upper` - index set for active variables at the upper bound
- `active_fixed` - index set for fixed variables
- `active`       - index set for all active variables
- `inactive`     - complementary index set for inactive variables

Level: developer

-seealso: `TAOBNCG`, `TAOBNTL`, `TAOBNTR`, `TaoBoundSolution()`

# External Links
$(_doc_external("Tao/TaoEstimateActiveBounds"))
"""
function TaoEstimateActiveBounds(petsclib::PetscLibType, X::PetscVec, XL::PetscVec, XU::PetscVec, G::PetscVec, S::PetscVec, W::PetscVec, steplen::PetscReal, active_lower::IS, active_upper::IS, active_fixed::IS, active::IS, inactive::IS) end

@for_petsc function TaoEstimateActiveBounds(petsclib::$UnionPetscLib, X::PetscVec, XL::PetscVec, XU::PetscVec, G::PetscVec, S::PetscVec, W::PetscVec, steplen::$PetscReal, active_lower::IS, active_upper::IS, active_fixed::IS, active::IS, inactive::IS )
	bound_tol_ = Ref{$PetscReal}()
	active_lower_ = Ref(active_lower.ptr)
	active_upper_ = Ref(active_upper.ptr)
	active_fixed_ = Ref(active_fixed.ptr)
	active_ = Ref(active.ptr)
	inactive_ = Ref(inactive.ptr)

    @chk ccall(
               (:TaoEstimateActiveBounds, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec, CVec, CVec, CVec, $PetscReal, Ptr{$PetscReal}, Ptr{CIS}, Ptr{CIS}, Ptr{CIS}, Ptr{CIS}, Ptr{CIS}),
               X, XL, XU, G, S, W, steplen, bound_tol_, active_lower_, active_upper_, active_fixed_, active_, inactive_,
              )

	bound_tol = bound_tol_[]
	active_lower.ptr = C_NULL
	active_upper.ptr = C_NULL
	active_fixed.ptr = C_NULL
	active.ptr = C_NULL
	inactive.ptr = C_NULL

	return bound_tol
end 

"""
	TaoBoundStep(petsclib::PetscLibType,X::PetscVec, XL::PetscVec, XU::PetscVec, active_lower::IS, active_upper::IS, active_fixed::IS, scale::PetscReal, S::PetscVec) 
Ensures the correct zero or adjusted step direction values for active
variables.

Input Parameters:
- `X`            - solution vector
- `XL`           - lower bound vector
- `XU`           - upper bound vector
- `active_lower` - index set for lower bounded active variables
- `active_upper` - index set for lower bounded active variables
- `active_fixed` - index set for fixed active variables
- `scale`        - amplification factor for the step that needs to be taken on actively bounded variables

Output Parameter:
- `S` - step direction to be modified

Level: developer

-seealso: `TAOBNCG`, `TAOBNTL`, `TAOBNTR`, `TaoBoundSolution()`

# External Links
$(_doc_external("Tao/TaoBoundStep"))
"""
function TaoBoundStep(petsclib::PetscLibType, X::PetscVec, XL::PetscVec, XU::PetscVec, active_lower::IS, active_upper::IS, active_fixed::IS, scale::PetscReal, S::PetscVec) end

@for_petsc function TaoBoundStep(petsclib::$UnionPetscLib, X::PetscVec, XL::PetscVec, XU::PetscVec, active_lower::IS, active_upper::IS, active_fixed::IS, scale::$PetscReal, S::PetscVec )

    @chk ccall(
               (:TaoBoundStep, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec, CIS, CIS, CIS, $PetscReal, CVec),
               X, XL, XU, active_lower, active_upper, active_fixed, scale, S,
              )


	return nothing
end 

"""
	nDiff::PetscInt = TaoBoundSolution(petsclib::PetscLibType,X::PetscVec, XL::PetscVec, XU::PetscVec, bound_tol::PetscReal, Xout::PetscVec) 
Ensures that the solution vector is snapped into the bounds within a given tolerance.

Collective

Input Parameters:
- `X`         - solution vector
- `XL`        - lower bound vector
- `XU`        - upper bound vector
- `bound_tol` - absolute tolerance in enforcing the bound

Output Parameters:
- `nDiff` - total number of vector entries that have been bounded
- `Xout`  - modified solution vector satisfying bounds to `bound_tol`

Level: developer

-seealso: `TAOBNCG`, `TAOBNTL`, `TAOBNTR`, `TaoBoundStep()`

# External Links
$(_doc_external("Tao/TaoBoundSolution"))
"""
function TaoBoundSolution(petsclib::PetscLibType, X::PetscVec, XL::PetscVec, XU::PetscVec, bound_tol::PetscReal, Xout::PetscVec) end

@for_petsc function TaoBoundSolution(petsclib::$UnionPetscLib, X::PetscVec, XL::PetscVec, XU::PetscVec, bound_tol::$PetscReal, Xout::PetscVec )
	nDiff_ = Ref{$PetscInt}()

    @chk ccall(
               (:TaoBoundSolution, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec, $PetscReal, Ptr{$PetscInt}, CVec),
               X, XL, XU, bound_tol, nDiff_, Xout,
              )

	nDiff = nDiff_[]

	return nDiff
end 

