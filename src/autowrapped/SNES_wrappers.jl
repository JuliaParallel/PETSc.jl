"""
	SNESAddOptionsChecker(petsclib::PetscLibType, snescheck::external) 
Adds an additional function to check for `SNES` options.

Not Collective

Input Parameter:
- `snescheck` - function that checks for options

Calling sequence of `snescheck`:
- `snes` - the `SNES` object for which it is checking options

Level: developer

See also: `SNES`, `SNESSetFromOptions()`

# External Links
$(_doc_external("SNES/SNESAddOptionsChecker"))
"""
function SNESAddOptionsChecker(petsclib::PetscLibType, snescheck::external)
    error("SNESAddOptionsChecker: no generated method for these argument types")
end

@for_petsc function SNESAddOptionsChecker(petsclib::$UnionPetscLib, snescheck::external )

    @chk ccall(
               (:SNESAddOptionsChecker, $petsc_library),
               PetscErrorCode,
               (external,),
               snescheck,
              )


	return nothing
end 

"""
	SNESAppendOptionsPrefix(petsclib::PetscLibType, snes::AbstractSNES, prefix::String) 
Appends to the prefix used for searching for all
`SNES` options in the database.

Logically Collective

Input Parameters:
- `snes`   - the `SNES` context
- `prefix` - the prefix to prepend to all option names

Level: advanced

See also: `SNESGetOptionsPrefix()`, `SNESSetOptionsPrefix()`

# External Links
$(_doc_external("SNES/SNESAppendOptionsPrefix"))
"""
function SNESAppendOptionsPrefix(petsclib::PetscLibType, snes::AbstractSNES, prefix::String)
    error("SNESAppendOptionsPrefix: no generated method for these argument types")
end

@for_petsc function SNESAppendOptionsPrefix(petsclib::$UnionPetscLib, snes::AbstractSNES, prefix::String )

    @chk ccall(
               (:SNESAppendOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Cchar}),
               snes, prefix,
              )


	return nothing
end 

"""
	SNESApplyNPC(petsclib::PetscLibType, snes::AbstractSNES, x::AbstractPetscVec, f::AbstractPetscVec, y::AbstractPetscVec) 
Calls `SNESSolve()` on the preconditioner for the `SNES`

Collective

Input Parameters:
- `snes` - the `SNES` context
- `x`    - input vector
- `f`    - optional; the function evaluation on `x`

Output Parameter:
- `y` - function vector, as set by `SNESSetFunction()`

Level: developer

See also: `SNES`, `SNESGetNPC()`, `SNESSetNPC()`, `SNESComputeFunction()`

# External Links
$(_doc_external("SNES/SNESApplyNPC"))
"""
function SNESApplyNPC(petsclib::PetscLibType, snes::AbstractSNES, x::AbstractPetscVec, f::AbstractPetscVec, y::AbstractPetscVec)
    error("SNESApplyNPC: no generated method for these argument types")
end

@for_petsc function SNESApplyNPC(petsclib::$UnionPetscLib, snes::AbstractSNES, x::AbstractPetscVec, f::AbstractPetscVec, y::AbstractPetscVec )

    @chk ccall(
               (:SNESApplyNPC, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec, CVec),
               snes, x, f, y,
              )


	return nothing
end 

"""
	SNESCompositeAddSNES(petsclib::PetscLibType, snes::AbstractSNES, type::SNESType) 
Adds another `SNES` to the `SNESCOMPOSITE`

Collective

Input Parameters:
- `snes` - the `SNES` context of type `SNESCOMPOSITE`
- `type` - the `SNESType` of the new solver

Level: developer

See also: `SNES`, `SNESCOMPOSITE`, `SNESCompositeGetSNES()`

# External Links
$(_doc_external("SNES/SNESCompositeAddSNES"))
"""
function SNESCompositeAddSNES(petsclib::PetscLibType, snes::AbstractSNES, type::SNESType)
    error("SNESCompositeAddSNES: no generated method for these argument types")
end

@for_petsc function SNESCompositeAddSNES(petsclib::$UnionPetscLib, snes::AbstractSNES, type::SNESType )

    @chk ccall(
               (:SNESCompositeAddSNES, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESType),
               snes, type,
              )


	return nothing
end 

"""
	n::PetscInt = SNESCompositeGetNumber(petsclib::PetscLibType, snes::AbstractSNES) 
Get the number of subsolvers in the `SNESCOMPOSITE`

Logically Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `n` - the number of subsolvers

Level: developer

See also: `SNES`, `SNESCOMPOSITE`, `SNESCompositeAddSNES()`, `SNESCompositeGetSNES()`

# External Links
$(_doc_external("SNES/SNESCompositeGetNumber"))
"""
function SNESCompositeGetNumber(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESCompositeGetNumber: no generated method for these argument types")
end

@for_petsc function SNESCompositeGetNumber(petsclib::$UnionPetscLib, snes::AbstractSNES )
	n_ = Ref{$PetscInt}()

    @chk ccall(
               (:SNESCompositeGetNumber, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscInt}),
               snes, n_,
              )

	n = n_[]

	return n
end 

"""
	subsnes::SNES = SNESCompositeGetSNES(petsclib::PetscLibType, snes::AbstractSNES, n::PetscInt) 
Gets one of the `SNES` objects in the `SNES` of `SNESType` `SNESCOMPOSITE`

Not Collective

Input Parameters:
- `snes` - the `SNES` context
- `n`    - the number of the composed `SNES` requested

Output Parameter:
- `subsnes` - the `SNES` requested

Level: developer

See also: `SNES`, `SNESCOMPOSITE`, `SNESCompositeAddSNES()`, `SNESCompositeGetNumber()`

# External Links
$(_doc_external("SNES/SNESCompositeGetSNES"))
"""
function SNESCompositeGetSNES(petsclib::PetscLibType, snes::AbstractSNES, n::Integer)
    error("SNESCompositeGetSNES: no generated method for these argument types")
end

@for_petsc function SNESCompositeGetSNES(petsclib::$UnionPetscLib, snes::AbstractSNES, n::$PetscInt )
	subsnes_ = Ref{CSNES}()

    @chk ccall(
               (:SNESCompositeGetSNES, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, Ptr{CSNES}),
               snes, n, subsnes_,
              )

	subsnes = SNES(subsnes_[], petsclib)

	return subsnes
end 

"""
	SNESCompositeSetDamping(petsclib::PetscLibType, snes::AbstractSNES, n::PetscInt, dmp::PetscReal) 
Sets the damping of a subsolver when using `SNES_COMPOSITE_ADDITIVE` with a `SNES` of `SNESType` `SNESCOMPOSITE`

Not Collective

Input Parameters:
- `snes` - the `SNES` context
- `n`    - the number of the sub-`SNES` object requested
- `dmp`  - the damping

Level: intermediate

See also: `SNES`, `SNESCOMPOSITE`, `SNESCompositeAddSNES()`, `SNESCompositeGetSNES()`,
`SNES_COMPOSITE_ADDITIVE`, `SNES_COMPOSITE_MULTIPLICATIVE`, `SNESCompositeType`, `SNESCompositeSetType()`

# External Links
$(_doc_external("SNES/SNESCompositeSetDamping"))
"""
function SNESCompositeSetDamping(petsclib::PetscLibType, snes::AbstractSNES, n::Integer, dmp::Real)
    error("SNESCompositeSetDamping: no generated method for these argument types")
end

@for_petsc function SNESCompositeSetDamping(petsclib::$UnionPetscLib, snes::AbstractSNES, n::$PetscInt, dmp::$PetscReal )

    @chk ccall(
               (:SNESCompositeSetDamping, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal),
               snes, n, dmp,
              )


	return nothing
end 

"""
	SNESCompositeSetType(petsclib::PetscLibType, snes::AbstractSNES, type::SNESCompositeType) 
Sets the type of composite preconditioner.

Logically Collective

Input Parameters:
- `snes` - the preconditioner context
- `type` - `SNES_COMPOSITE_ADDITIVE` (default), `SNES_COMPOSITE_MULTIPLICATIVE`, or `SNES_COMPOSITE_ADDITIVEOPTIMAL`

Options Database Key:
- `-snes_composite_type (multiplicative|additive|additive_optimal)` - Sets composite preconditioner type

Level: developer

See also: `SNES_COMPOSITE_ADDITIVE`, `SNES_COMPOSITE_MULTIPLICATIVE`, `SNESCompositeType`, `SNESCOMPOSITE`, `SNES_COMPOSITE_ADDITIVEOPTIMAL`,
`PCCompositeType`

# External Links
$(_doc_external("SNES/SNESCompositeSetType"))
"""
function SNESCompositeSetType(petsclib::PetscLibType, snes::AbstractSNES, type::SNESCompositeType)
    error("SNESCompositeSetType: no generated method for these argument types")
end

@for_petsc function SNESCompositeSetType(petsclib::$UnionPetscLib, snes::AbstractSNES, type::SNESCompositeType )

    @chk ccall(
               (:SNESCompositeSetType, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESCompositeType),
               snes, type,
              )


	return nothing
end 

"""
	SNESComputeFunction(petsclib::PetscLibType, snes::AbstractSNES, x::AbstractPetscVec, f::AbstractPetscVec) 
Calls the function that has been set with `SNESSetFunction()`.

Collective

Input Parameters:
- `snes` - the `SNES` context
- `x`    - input vector

Output Parameter:
- `f` - function vector, as set by `SNESSetFunction()`

Level: developer

See also: `SNES`, `SNESSetFunction()`, `SNESGetFunction()`, `SNESComputeMFFunction()`, `SNESSetFunctionDomainError()`

# External Links
$(_doc_external("SNES/SNESComputeFunction"))
"""
function SNESComputeFunction(petsclib::PetscLibType, snes::AbstractSNES, x::AbstractPetscVec, f::AbstractPetscVec)
    error("SNESComputeFunction: no generated method for these argument types")
end

@for_petsc function SNESComputeFunction(petsclib::$UnionPetscLib, snes::AbstractSNES, x::AbstractPetscVec, f::AbstractPetscVec )

    @chk ccall(
               (:SNESComputeFunction, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec),
               snes, x, f,
              )


	return nothing
end 

"""
	SNESComputeFunctionDefaultNPC(petsclib::PetscLibType, snes::AbstractSNES, X::AbstractPetscVec, F::AbstractPetscVec) 
Compute the residual by applying the attached nonlinear preconditioner when one is present, otherwise defer to `SNESComputeFunction()`

Collective

Input Parameters:
- `snes` - the `SNES` context
- `X`    - the current iterate

Output Parameter:
- `F` - the residual vector produced by the nonlinear preconditioner (or the standard function evaluation)

Level: developer

See also: `SNES`, `SNESSetNPC()`, `SNESApplyNPC()`, `SNESComputeFunction()`, `SNESGetNPCFunction()`

# External Links
$(_doc_external("SNES/SNESComputeFunctionDefaultNPC"))
"""
function SNESComputeFunctionDefaultNPC(petsclib::PetscLibType, snes::AbstractSNES, X::AbstractPetscVec, F::AbstractPetscVec)
    error("SNESComputeFunctionDefaultNPC: no generated method for these argument types")
end

@for_petsc function SNESComputeFunctionDefaultNPC(petsclib::$UnionPetscLib, snes::AbstractSNES, X::AbstractPetscVec, F::AbstractPetscVec )

    @chk ccall(
               (:SNESComputeFunctionDefaultNPC, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec),
               snes, X, F,
              )


	return nothing
end 

"""
	SNESComputeJacobian(petsclib::PetscLibType, snes::AbstractSNES, X::AbstractPetscVec, A::AbstractPetscMat, B::AbstractPetscMat) 
Computes the Jacobian matrix that has been set with `SNESSetJacobian()`.

Collective

Input Parameters:
- `snes` - the `SNES` context
- `X`    - input vector

Output Parameters:
- `A` - Jacobian matrix
- `B` - optional matrix for building the preconditioner, usually the same as `A`

Options Database Keys:
- `-snes_lag_preconditioner lag`          - how often to rebuild preconditioner
- `-snes_lag_jacobian lag`                - how often to rebuild Jacobian
- `-snes_test_jacobian [threshold]`       - compare the user provided Jacobian with one compute via finite differences to check for errors.
If a threshold is given, display only those entries whose difference is greater than the threshold.
- `-snes_test_jacobian_view [viewer]`     - display the user provided Jacobian, the finite difference Jacobian and the difference between them to help users detect the location of errors in the user provided Jacobian
- `-snes_compare_explicit`                - Compare the computed Jacobian to the finite difference Jacobian and output the differences
- `-snes_compare_explicit_draw`           - Compare the computed Jacobian to the finite difference Jacobian and draw the result
- `-snes_compare_explicit_contour`        - Compare the computed Jacobian to the finite difference Jacobian and draw a contour plot with the result
- `-snes_compare_operator`                - Make the comparison options above use the operator instead of the matrix used to construct the preconditioner
- `-snes_compare_coloring`                - Compute the finite difference Jacobian using coloring and display norms of difference
- `-snes_compare_coloring_display`        - Compute the finite difference Jacobian using coloring and display verbose differences
- `-snes_compare_coloring_threshold`      - Display only those matrix entries that differ by more than a given threshold
- `-snes_compare_coloring_threshold_atol` - Absolute tolerance for difference in matrix entries to be displayed by `-snes_compare_coloring_threshold`
- `-snes_compare_coloring_threshold_rtol` - Relative tolerance for difference in matrix entries to be displayed by `-snes_compare_coloring_threshold`
- `-snes_compare_coloring_draw`           - Compute the finite difference Jacobian using coloring and draw differences
- `-snes_compare_coloring_draw_contour`   - Compute the finite difference Jacobian using coloring and show contours of matrices and differences

Level: developer

See also: `SNESSetJacobian()`, `KSPSetOperators()`, `MatStructure`, `SNESSetLagPreconditioner()`, `SNESSetLagJacobian()`,
`SNESSetJacobianDomainError()`, `SNESCheckJacobianDomainError()`, `SNESSetCheckJacobianDomainError()`

# External Links
$(_doc_external("SNES/SNESComputeJacobian"))
"""
function SNESComputeJacobian(petsclib::PetscLibType, snes::AbstractSNES, X::AbstractPetscVec, A::AbstractPetscMat, B::AbstractPetscMat)
    error("SNESComputeJacobian: no generated method for these argument types")
end

@for_petsc function SNESComputeJacobian(petsclib::$UnionPetscLib, snes::AbstractSNES, X::AbstractPetscVec, A::AbstractPetscMat, B::AbstractPetscMat )

    @chk ccall(
               (:SNESComputeJacobian, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CMat, CMat),
               snes, X, A, B,
              )


	return nothing
end 

"""
	SNESComputeJacobianDefault(petsclib::PetscLibType, snes::AbstractSNES, x1::AbstractPetscVec, J::AbstractPetscMat, B::AbstractPetscMat, ctx::Ptr{Cvoid}) 
Computes the Jacobian using finite differences.

Collective

Input Parameters:
- `snes` - the `SNES` context
- `x1`   - compute Jacobian at this point
- `ctx`  - application's function context, as set with `SNESSetFunction()`

Output Parameters:
- `J` - Jacobian matrix (not altered in this routine)
- `B` - newly computed Jacobian matrix to use with preconditioner (generally the same as `J`)

Options Database Keys:
- `-snes_fd`             - Activates `SNESComputeJacobianDefault()`
- `-snes_fd_coloring`    - Activates a faster computation that uses a graph coloring of the matrix
- `-snes_test_err etol`  - Square root of function error tolerance, default square root of machine
epsilon (1.e-8 in double, 3.e-4 in single)
- `-mat_fd_type (wp|ds)` - See `MATMFFD_WP` and `MATMFFD_DS`

Level: intermediate

See also: `SNES`, `SNESSetJacobian()`, `SNESComputeJacobianDefaultColor()`, `MatCreateSNESMF()`

# External Links
$(_doc_external("SNES/SNESComputeJacobianDefault"))
"""
function SNESComputeJacobianDefault(petsclib::PetscLibType, snes::AbstractSNES, x1::AbstractPetscVec, J::AbstractPetscMat, B::AbstractPetscMat, ctx::Ptr{Cvoid})
    error("SNESComputeJacobianDefault: no generated method for these argument types")
end

@for_petsc function SNESComputeJacobianDefault(petsclib::$UnionPetscLib, snes::AbstractSNES, x1::AbstractPetscVec, J::AbstractPetscMat, B::AbstractPetscMat, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESComputeJacobianDefault, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CMat, CMat, Ptr{Cvoid}),
               snes, x1, J, B, ctx,
              )


	return nothing
end 

"""
	SNESComputeJacobianDefaultColor(petsclib::PetscLibType, snes::AbstractSNES, x1::AbstractPetscVec, J::AbstractPetscMat, B::AbstractPetscMat, ctx::Ptr{Cvoid}) 
Computes the Jacobian using
finite differences and coloring to exploit matrix sparsity.

Collective

Input Parameters:
- `snes` - nonlinear solver object
- `x1`   - location at which to evaluate Jacobian
- `ctx`  - `MatFDColoring` context or `NULL`

Output Parameters:
- `J` - Jacobian matrix (not altered in this routine)
- `B` - newly computed Jacobian matrix to use with preconditioner (generally the same as `J`)

Options Database Keys:
- `-snes_fd_color_use_mat`     - use a matrix coloring from the explicit matrix nonzero pattern instead of from the `DM` providing the matrix
- `-snes_fd_color`             - Activates `SNESComputeJacobianDefaultColor()` in `SNESSetFromOptions()`
- `-mat_fd_coloring_err err`   - Sets err (square root of relative error in the function)
- `-mat_fd_coloring_umin umin` - Sets umin, the minimum allowable u-value magnitude
- `-mat_fd_type`               - Either wp or ds (see `MATMFFD_WP` or `MATMFFD_DS`)
- `-snes_mf_operator`          - Use matrix-free application of Jacobian
- `-snes_mf`                   - Use matrix-free Jacobian with no explicit Jacobian representation

See also: `SNES`, `SNESSetJacobian()`, `SNESTestJacobian()`, `SNESComputeJacobianDefault()`, `SNESSetUseMatrixFree()`,
`MatFDColoringCreate()`, `MatFDColoringSetFunction()`

# External Links
$(_doc_external("SNES/SNESComputeJacobianDefaultColor"))
"""
function SNESComputeJacobianDefaultColor(petsclib::PetscLibType, snes::AbstractSNES, x1::AbstractPetscVec, J::AbstractPetscMat, B::AbstractPetscMat, ctx::Ptr{Cvoid})
    error("SNESComputeJacobianDefaultColor: no generated method for these argument types")
end

@for_petsc function SNESComputeJacobianDefaultColor(petsclib::$UnionPetscLib, snes::AbstractSNES, x1::AbstractPetscVec, J::AbstractPetscMat, B::AbstractPetscMat, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESComputeJacobianDefaultColor, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CMat, CMat, Ptr{Cvoid}),
               snes, x1, J, B, ctx,
              )


	return nothing
end 

"""
	SNESComputeMFFunction(petsclib::PetscLibType, snes::AbstractSNES, x::AbstractPetscVec, y::AbstractPetscVec) 
Calls the function that has been set with `DMSNESSetMFFunction()`.

Collective

Input Parameters:
- `snes` - the `SNES` context
- `x`    - input vector

Output Parameter:
- `y` - output vector

Level: developer

See also: `SNES`, `SNESSetFunction()`, `SNESGetFunction()`, `SNESComputeFunction()`, `MatCreateSNESMF()`, `DMSNESSetMFFunction()`

# External Links
$(_doc_external("SNES/SNESComputeMFFunction"))
"""
function SNESComputeMFFunction(petsclib::PetscLibType, snes::AbstractSNES, x::AbstractPetscVec, y::AbstractPetscVec)
    error("SNESComputeMFFunction: no generated method for these argument types")
end

@for_petsc function SNESComputeMFFunction(petsclib::$UnionPetscLib, snes::AbstractSNES, x::AbstractPetscVec, y::AbstractPetscVec )

    @chk ccall(
               (:SNESComputeMFFunction, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec),
               snes, x, y,
              )


	return nothing
end 

"""
	SNESComputeNGS(petsclib::PetscLibType, snes::AbstractSNES, b::AbstractPetscVec, x::AbstractPetscVec) 
Calls the Gauss-Seidel function that has been set with `SNESSetNGS()`.

Collective

Input Parameters:
- `snes` - the `SNES` context
- `x`    - input vector
- `b`    - rhs vector

Output Parameter:
- `x` - new solution vector

Level: developer

See also: `SNESNGSFn`, `SNESSetNGS()`, `SNESComputeFunction()`, `SNESNGS`

# External Links
$(_doc_external("SNES/SNESComputeNGS"))
"""
function SNESComputeNGS(petsclib::PetscLibType, snes::AbstractSNES, b::AbstractPetscVec, x::AbstractPetscVec)
    error("SNESComputeNGS: no generated method for these argument types")
end

@for_petsc function SNESComputeNGS(petsclib::$UnionPetscLib, snes::AbstractSNES, b::AbstractPetscVec, x::AbstractPetscVec )

    @chk ccall(
               (:SNESComputeNGS, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec),
               snes, b, x,
              )


	return nothing
end 

"""
	ob::PetscReal = SNESComputeObjective(petsclib::PetscLibType, snes::AbstractSNES, X::AbstractPetscVec) 
Computes the objective function that has been provided by `SNESSetObjective()`

Collective

Input Parameters:
- `snes` - the `SNES` context
- `X`    - the state vector

Output Parameter:
- `ob` - the objective value

Level: developer

See also: `SNESLineSearch`, `SNES`, `SNESSetObjective()`, `SNESGetSolution()`

# External Links
$(_doc_external("SNES/SNESComputeObjective"))
"""
function SNESComputeObjective(petsclib::PetscLibType, snes::AbstractSNES, X::AbstractPetscVec)
    error("SNESComputeObjective: no generated method for these argument types")
end

@for_petsc function SNESComputeObjective(petsclib::$UnionPetscLib, snes::AbstractSNES, X::AbstractPetscVec )
	ob_ = Ref{$PetscReal}()

    @chk ccall(
               (:SNESComputeObjective, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, Ptr{$PetscReal}),
               snes, X, ob_,
              )

	ob = ob_[]

	return ob
end 

"""
	SNESConverged(petsclib::PetscLibType, snes::AbstractSNES, it::PetscInt, xnorm::PetscReal, snorm::PetscReal, fnorm::PetscReal) 
Run the convergence test and update the `SNESConvergedReason`.

Collective

Input Parameters:
- `snes`  - the `SNES` context
- `it`    - current iteration
- `xnorm` - 2-norm of current iterate
- `snorm` - 2-norm of current step
- `fnorm` - 2-norm of function

Level: developer

See also: `SNES`, `SNESSolve`, `SNESSetConvergenceTest()`

# External Links
$(_doc_external("SNES/SNESConverged"))
"""
function SNESConverged(petsclib::PetscLibType, snes::AbstractSNES, it::Integer, xnorm::Real, snorm::Real, fnorm::Real)
    error("SNESConverged: no generated method for these argument types")
end

@for_petsc function SNESConverged(petsclib::$UnionPetscLib, snes::AbstractSNES, it::$PetscInt, xnorm::$PetscReal, snorm::$PetscReal, fnorm::$PetscReal )

    @chk ccall(
               (:SNESConverged, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, $PetscReal, $PetscReal),
               snes, it, xnorm, snorm, fnorm,
              )


	return nothing
end 

"""
	reason::SNESConvergedReason = SNESConvergedCorrectPressure(petsclib::PetscLibType, snes::AbstractSNES, it::PetscInt, xnorm::PetscReal, gnorm::PetscReal, f::PetscReal, ctx::Ptr{Cvoid}) 
The regular `SNES` convergence test that, up on convergence, adds a vector in the nullspace
to make the continuum integral of the pressure field equal to zero.

Logically Collective

Input Parameters:
- `snes`  - the `SNES` context
- `it`    - the iteration (0 indicates before any Newton steps)
- `xnorm` - 2-norm of current iterate
- `gnorm` - 2-norm of current step
- `f`     - 2-norm of function at current iterate
- `ctx`   - Optional application context

Output Parameter:
- `reason` - `SNES_CONVERGED_ITERATING`, `SNES_CONVERGED_ITS`, or `SNES_DIVERGED_FUNCTION_NANORINF`

Options Database Key:
- `-snes_convergence_test correct_pressure` - see `SNESSetFromOptions()`

Level: advanced

See also: `SNES`, `DM`, `SNESConvergedDefault()`, `SNESSetConvergenceTest()`, `DMSetNullSpaceConstructor()`

# External Links
$(_doc_external("SNES/SNESConvergedCorrectPressure"))
"""
function SNESConvergedCorrectPressure(petsclib::PetscLibType, snes::AbstractSNES, it::Integer, xnorm::Real, gnorm::Real, f::Real, ctx::Ptr{Cvoid})
    error("SNESConvergedCorrectPressure: no generated method for these argument types")
end

@for_petsc function SNESConvergedCorrectPressure(petsclib::$UnionPetscLib, snes::AbstractSNES, it::$PetscInt, xnorm::$PetscReal, gnorm::$PetscReal, f::$PetscReal, ctx::Ptr{Cvoid} )
	reason_ = Ref{SNESConvergedReason}()

    @chk ccall(
               (:SNESConvergedCorrectPressure, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, $PetscReal, $PetscReal, Ptr{SNESConvergedReason}, Ptr{Cvoid}),
               snes, it, xnorm, gnorm, f, reason_, ctx,
              )

	reason = reason_[]

	return reason
end 

"""
	reason::SNESConvergedReason = SNESConvergedDefault(petsclib::PetscLibType, snes::AbstractSNES, it::PetscInt, xnorm::PetscReal, snorm::PetscReal, fnorm::PetscReal, ctx::Ptr{Cvoid}) 
Default convergence test for `SNESSolve()`.

Collective

Input Parameters:
- `snes`  - the `SNES` context
- `it`    - the iteration (0 indicates before any Newton steps)
- `xnorm` - 2-norm of current iterate
- `snorm` - 2-norm of current step
- `fnorm` - 2-norm of function at current iterate
- `ctx`   - unused context

Output Parameter:
- `reason` - converged reason, see `SNESConvergedReason`

Options Database Keys:
- `-snes_convergence_test default`    - see `SNESSetFromOptions()`
- `-snes_stol`                        - convergence tolerance in terms of the norm of the change in the solution between steps
- `-snes_atol abstol`                 - absolute tolerance of residual norm
- `-snes_rtol rtol`                   - relative decrease in tolerance norm from the initial 2-norm of the solution
- `-snes_divergence_tolerance divtol` - if the residual goes above divtol*rnorm0, exit with divergence
- `-snes_max_funcs max_funcs`         - maximum number of function evaluations, use `unlimited` for no maximum
- `-snes_max_fail max_fail`           - maximum number of line search failures allowed before stopping, default is none
- `-snes_max_linear_solve_fail`       - number of linear solver failures before `SNESSolve()` stops

Level: developer

See also: `SNES`, `SNESSolve()`, `SNESSetConvergenceTest()`, `SNESConvergedSkip()`, `SNESSetTolerances()`, `SNESSetDivergenceTolerance()`,
`SNESConvergedReason`

# External Links
$(_doc_external("SNES/SNESConvergedDefault"))
"""
function SNESConvergedDefault(petsclib::PetscLibType, snes::AbstractSNES, it::Integer, xnorm::Real, snorm::Real, fnorm::Real, ctx::Ptr{Cvoid})
    error("SNESConvergedDefault: no generated method for these argument types")
end

@for_petsc function SNESConvergedDefault(petsclib::$UnionPetscLib, snes::AbstractSNES, it::$PetscInt, xnorm::$PetscReal, snorm::$PetscReal, fnorm::$PetscReal, ctx::Ptr{Cvoid} )
	reason_ = Ref{SNESConvergedReason}()

    @chk ccall(
               (:SNESConvergedDefault, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, $PetscReal, $PetscReal, Ptr{SNESConvergedReason}, Ptr{Cvoid}),
               snes, it, xnorm, snorm, fnorm, reason_, ctx,
              )

	reason = reason_[]

	return reason
end 

"""
	SNESConvergedReasonView(petsclib::PetscLibType, snes::AbstractSNES, viewer::PetscViewer) 
Displays the reason a `SNES` solve converged or diverged to a viewer

Collective

Input Parameters:
- `snes`   - iterative context obtained from `SNESCreate()`
- `viewer` - the viewer to display the reason

Options Database Keys:
- `-snes_converged_reason`          - print reason for converged or diverged, also prints number of iterations
- `-snes_converged_reason ::failed` - only print reason and number of iterations when diverged

Level: beginner

See also: `SNESConvergedReason`, `PetscViewer`, `SNES`,
`SNESCreate()`, `SNESSetUp()`, `SNESDestroy()`, `SNESSetTolerances()`, `SNESConvergedDefault()`, `SNESGetConvergedReason()`,
`SNESConvergedReasonViewFromOptions()`,
`PetscViewerPushFormat()`, `PetscViewerPopFormat()`

# External Links
$(_doc_external("SNES/SNESConvergedReasonView"))
"""
function SNESConvergedReasonView(petsclib::PetscLibType, snes::AbstractSNES, viewer::PetscViewer)
    error("SNESConvergedReasonView: no generated method for these argument types")
end

@for_petsc function SNESConvergedReasonView(petsclib::$UnionPetscLib, snes::AbstractSNES, viewer::PetscViewer )

    @chk ccall(
               (:SNESConvergedReasonView, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscViewer),
               snes, viewer,
              )


	return nothing
end 

"""
	SNESConvergedReasonViewCancel(petsclib::PetscLibType, snes::AbstractSNES) 
Clears all the reason view functions for a `SNES` object provided with `SNESConvergedReasonViewSet()` also
removes the default viewer.

Collective

Input Parameter:
- `snes` - the nonlinear iterative solver context obtained from `SNESCreate()`

Level: intermediate

See also: `SNES`, `SNESCreate()`, `SNESDestroy()`, `SNESReset()`, `SNESConvergedReasonViewSet()`

# External Links
$(_doc_external("SNES/SNESConvergedReasonViewCancel"))
"""
function SNESConvergedReasonViewCancel(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESConvergedReasonViewCancel: no generated method for these argument types")
end

@for_petsc function SNESConvergedReasonViewCancel(petsclib::$UnionPetscLib, snes::AbstractSNES )

    @chk ccall(
               (:SNESConvergedReasonViewCancel, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	SNESConvergedReasonViewFromOptions(petsclib::PetscLibType, snes::AbstractSNES) 
Processes command line options to determine if/how a `SNESConvergedReason` is to be viewed at the end of `SNESSolve()`
All the user-provided viewer routines set with `SNESConvergedReasonViewSet()` will be called, if they exist.

Collective

Input Parameter:
- `snes` - the `SNES` object

Level: advanced

See also: `SNES`, `SNESConvergedReason`, `SNESConvergedReasonViewSet()`, `SNESCreate()`, `SNESSetUp()`, `SNESDestroy()`,
`SNESSetTolerances()`, `SNESConvergedDefault()`, `SNESGetConvergedReason()`, `SNESConvergedReasonView()`

# External Links
$(_doc_external("SNES/SNESConvergedReasonViewFromOptions"))
"""
function SNESConvergedReasonViewFromOptions(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESConvergedReasonViewFromOptions: no generated method for these argument types")
end

@for_petsc function SNESConvergedReasonViewFromOptions(petsclib::$UnionPetscLib, snes::AbstractSNES )

    @chk ccall(
               (:SNESConvergedReasonViewFromOptions, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	SNESConvergedReasonViewSet(petsclib::PetscLibType, snes::AbstractSNES, f::external, vctx::Ptr{Cvoid}, reasonviewdestroy::Ptr{Cvoid}) 
Sets an ADDITIONAL function that is to be used at the
end of the nonlinear solver to display the convergence reason of the nonlinear solver.

Logically Collective

Input Parameters:
- `snes`              - the `SNES` context
- `f`                 - the `SNESConvergedReason` view function
- `vctx`              - [optional] user-defined context for private data for the `SNESConvergedReason` view function (use `NULL` if no context is desired)
- `reasonviewdestroy` - [optional] routine that frees the context (may be `NULL`), see `PetscCtxDestroyFn` for the calling sequence

Calling sequence of `f`:
- `snes` - the `SNES` context
- `vctx` - [optional] context for private data for the function

Options Database Keys:
- `-snes_converged_reason`             - sets a default `SNESConvergedReasonView()`
- `-snes_converged_reason_view_cancel` - cancels all converged reason viewers that have been hardwired into a code by
calls to `SNESConvergedReasonViewSet()`, but does not cancel those set via the options database.

Level: intermediate

See also: `SNES`, `SNESSolve()`, `SNESConvergedReason`, `SNESGetConvergedReason()`, `SNESConvergedReasonView()`, `SNESConvergedReasonViewCancel()`,
`PetscCtxDestroyFn`

# External Links
$(_doc_external("SNES/SNESConvergedReasonViewSet"))
"""
function SNESConvergedReasonViewSet(petsclib::PetscLibType, snes::AbstractSNES, f::external, vctx::Ptr{Cvoid}, reasonviewdestroy::Ptr{Cvoid})
    error("SNESConvergedReasonViewSet: no generated method for these argument types")
end

@for_petsc function SNESConvergedReasonViewSet(petsclib::$UnionPetscLib, snes::AbstractSNES, f::external, vctx::Ptr{Cvoid}, reasonviewdestroy::Ptr{Cvoid} )

    @chk ccall(
               (:SNESConvergedReasonViewSet, $petsc_library),
               PetscErrorCode,
               (CSNES, external, Ptr{Cvoid}, Ptr{Cvoid}),
               snes, f, vctx, reasonviewdestroy,
              )


	return nothing
end 

"""
	reason::SNESConvergedReason = SNESConvergedSkip(petsclib::PetscLibType, snes::AbstractSNES, it::PetscInt, xnorm::PetscReal, snorm::PetscReal, fnorm::PetscReal, ctx::Ptr{Cvoid}) 
Convergence test for `SNES` that NEVER returns as
converged, UNLESS the maximum number of iteration have been reached.

Logically Collective

Input Parameters:
- `snes`  - the `SNES` context
- `it`    - the iteration (0 indicates before any Newton steps)
- `xnorm` - 2-norm of current iterate
- `snorm` - 2-norm of current step
- `fnorm` - 2-norm of function at current iterate
- `ctx`   - unused context

Output Parameter:
- `reason` - `SNES_CONVERGED_ITERATING`, `SNES_CONVERGED_ITS`, or `SNES_DIVERGED_FUNCTION_NANORINF`

Options Database Key:
- `-snes_convergence_test skip` - see `SNESSetFromOptions()`

Level: advanced

See also: `SNES`, `SNESSolve()`, `SNESConvergedDefault()`, `SNESSetConvergenceTest()`, `SNESConvergedReason`

# External Links
$(_doc_external("SNES/SNESConvergedSkip"))
"""
function SNESConvergedSkip(petsclib::PetscLibType, snes::AbstractSNES, it::Integer, xnorm::Real, snorm::Real, fnorm::Real, ctx::Ptr{Cvoid})
    error("SNESConvergedSkip: no generated method for these argument types")
end

@for_petsc function SNESConvergedSkip(petsclib::$UnionPetscLib, snes::AbstractSNES, it::$PetscInt, xnorm::$PetscReal, snorm::$PetscReal, fnorm::$PetscReal, ctx::Ptr{Cvoid} )
	reason_ = Ref{SNESConvergedReason}()

    @chk ccall(
               (:SNESConvergedSkip, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, $PetscReal, $PetscReal, Ptr{SNESConvergedReason}, Ptr{Cvoid}),
               snes, it, xnorm, snorm, fnorm, reason_, ctx,
              )

	reason = reason_[]

	return reason
end 

"""
	outsnes::SNES = SNESCreate(petsclib::PetscLibType, comm::MPI_Comm) 
Creates a nonlinear solver context used to manage a set of nonlinear solves

Collective

Input Parameter:
- `comm` - MPI communicator

Output Parameter:
- `outsnes` - the new `SNES` context

Options Database Keys:
- `-snes_mf`          - Activates default matrix-free Jacobian-vector products, and no matrix to construct a preconditioner
- `-snes_mf_operator` - Activates default matrix-free Jacobian-vector products, and a user-provided matrix as set by `SNESSetJacobian()`
- `-snes_fd_coloring` - uses a relative fast computation of the Jacobian using finite differences and a graph coloring
- `-snes_fd`          - Uses (slow!) finite differences to compute Jacobian

Level: beginner

See also: `SNES`, `SNESSolve()`, `SNESDestroy()`, `SNESSetLagPreconditioner()`, `SNESSetLagJacobian()`

# External Links
$(_doc_external("SNES/SNESCreate"))
"""
function SNESCreate(petsclib::PetscLibType, comm::MPI_Comm)
    error("SNESCreate: no generated method for these argument types")
end

@for_petsc function SNESCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	outsnes_ = Ref{CSNES}()

    @chk ccall(
               (:SNESCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{CSNES}),
               comm, outsnes_,
              )

	outsnes = SNES(outsnes_[], petsclib)

	return outsnes
end 

"""
	SNESDestroy(petsclib::PetscLibType, snes::AbstractSNES) 
Destroys the nonlinear solver context that was created
with `SNESCreate()`.

Collective

Input Parameter:
- `snes` - the `SNES` context

Level: beginner

See also: `SNES`, `SNESCreate()`, `SNESSolve()`

# External Links
$(_doc_external("SNES/SNESDestroy"))
"""
function SNESDestroy(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESDestroy: no generated method for these argument types")
end

@for_petsc function SNESDestroy(petsclib::$UnionPetscLib, snes::AbstractSNES )
	snes_ = Ref(snes.ptr)

    @chk ccall(
               (:SNESDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{CSNES},),
               snes_,
              )

	snes.ptr = C_NULL

	return nothing
end 

"""
	Xcoarse::PetscVec = SNESFASCreateCoarseVec(petsclib::PetscLibType, snes::AbstractSNES) 
create a `Vec` corresponding to a state vector on one level coarser than the current level

Collective

Input Parameter:
- `snes` - `SNESFAS` object

Output Parameter:
- `Xcoarse` - vector on level one coarser than the current level

Level: developer

See also: `SNESFASSetRestriction()`, `SNESFASRestrict()`, `SNESFAS`

# External Links
$(_doc_external("SNESFAS/SNESFASCreateCoarseVec"))
"""
function SNESFASCreateCoarseVec(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESFASCreateCoarseVec: no generated method for these argument types")
end

@for_petsc function SNESFASCreateCoarseVec(petsclib::$UnionPetscLib, snes::AbstractSNES )
	Xcoarse_ = Ref{CVec}()

    @chk ccall(
               (:SNESFASCreateCoarseVec, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CVec}),
               snes, Xcoarse_,
              )

	Xcoarse = PetscVec(Xcoarse_[], petsclib)

	return Xcoarse
end 

"""
	correction::SNES = SNESFASCycleGetCorrection(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the coarse correction `SNESFAS` context for this level

Logically Collective

Input Parameter:
- `snes` - the `SNESFAS` obtained with `SNESFASGetCycleSNES()`

Output Parameter:
- `correction` - the coarse correction solve on this level

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASCycleGetSmootherUp()`, `SNESFASCycleGetSmoother()`

# External Links
$(_doc_external("SNESFAS/SNESFASCycleGetCorrection"))
"""
function SNESFASCycleGetCorrection(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESFASCycleGetCorrection: no generated method for these argument types")
end

@for_petsc function SNESFASCycleGetCorrection(petsclib::$UnionPetscLib, snes::AbstractSNES )
	correction_ = Ref{CSNES}()

    @chk ccall(
               (:SNESFASCycleGetCorrection, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CSNES}),
               snes, correction_,
              )

	correction = SNES(correction_[], petsclib)

	return correction
end 

"""
	mat::PetscMat = SNESFASCycleGetInjection(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the injection on a level

Logically Collective

Input Parameter:
- `snes` - the `SNESFAS` obtained with `SNESFASGetCycleSNES()`

Output Parameter:
- `mat` - the restriction operator on this level

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASGetInjection()`, `SNESFASCycleGetRestriction()`

# External Links
$(_doc_external("SNESFAS/SNESFASCycleGetInjection"))
"""
function SNESFASCycleGetInjection(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESFASCycleGetInjection: no generated method for these argument types")
end

@for_petsc function SNESFASCycleGetInjection(petsclib::$UnionPetscLib, snes::AbstractSNES )
	mat_ = Ref{CMat}()

    @chk ccall(
               (:SNESFASCycleGetInjection, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CMat}),
               snes, mat_,
              )

	mat = PetscMat(mat_[], petsclib)

	return mat
end 

"""
	mat::PetscMat = SNESFASCycleGetInterpolation(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the interpolation on a level

Logically Collective

Input Parameter:
- `snes` - the `SNESFAS` obtained with `SNESFASGetCycleSNES()`

Output Parameter:
- `mat` - the interpolation operator on this level

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASCycleGetSmootherUp()`, `SNESFASCycleGetSmoother()`

# External Links
$(_doc_external("SNESFAS/SNESFASCycleGetInterpolation"))
"""
function SNESFASCycleGetInterpolation(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESFASCycleGetInterpolation: no generated method for these argument types")
end

@for_petsc function SNESFASCycleGetInterpolation(petsclib::$UnionPetscLib, snes::AbstractSNES )
	mat_ = Ref{CMat}()

    @chk ccall(
               (:SNESFASCycleGetInterpolation, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CMat}),
               snes, mat_,
              )

	mat = PetscMat(mat_[], petsclib)

	return mat
end 

"""
	vec::PetscVec = SNESFASCycleGetRScale(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the injection scale-factor on a level

Logically Collective

Input Parameter:
- `snes` - the  `SNESFAS` obtained with `SNESFASGetCycleSNES()`

Output Parameter:
- `vec` - the restriction operator on this level

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASCycleGetRestriction()`, `SNESFASGetRScale()`

# External Links
$(_doc_external("SNESFAS/SNESFASCycleGetRScale"))
"""
function SNESFASCycleGetRScale(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESFASCycleGetRScale: no generated method for these argument types")
end

@for_petsc function SNESFASCycleGetRScale(petsclib::$UnionPetscLib, snes::AbstractSNES )
	vec_ = Ref{CVec}()

    @chk ccall(
               (:SNESFASCycleGetRScale, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CVec}),
               snes, vec_,
              )

	vec = PetscVec(vec_[], petsclib)

	return vec
end 

"""
	mat::PetscMat = SNESFASCycleGetRestriction(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the restriction on a level

Logically Collective

Input Parameter:
- `snes` - the `SNESFAS` obtained with `SNESFASGetCycleSNES()`

Output Parameter:
- `mat` - the restriction operator on this level

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASGetRestriction()`, `SNESFASCycleGetInterpolation()`

# External Links
$(_doc_external("SNESFAS/SNESFASCycleGetRestriction"))
"""
function SNESFASCycleGetRestriction(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESFASCycleGetRestriction: no generated method for these argument types")
end

@for_petsc function SNESFASCycleGetRestriction(petsclib::$UnionPetscLib, snes::AbstractSNES )
	mat_ = Ref{CMat}()

    @chk ccall(
               (:SNESFASCycleGetRestriction, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CMat}),
               snes, mat_,
              )

	mat = PetscMat(mat_[], petsclib)

	return mat
end 

"""
	smooth::SNES = SNESFASCycleGetSmoother(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the smoother on a particular cycle level.

Logically Collective

Input Parameter:
- `snes` - the `SNESFAS` obtained with `SNESFASGetCycleSNES()`

Output Parameter:
- `smooth` - the smoother

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASCycleGetSmootherUp()`, `SNESFASCycleGetSmootherDown()`, `SNESFASGetCycleSNES()`

# External Links
$(_doc_external("SNESFAS/SNESFASCycleGetSmoother"))
"""
function SNESFASCycleGetSmoother(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESFASCycleGetSmoother: no generated method for these argument types")
end

@for_petsc function SNESFASCycleGetSmoother(petsclib::$UnionPetscLib, snes::AbstractSNES )
	smooth_ = Ref{CSNES}()

    @chk ccall(
               (:SNESFASCycleGetSmoother, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CSNES}),
               snes, smooth_,
              )

	smooth = SNES(smooth_[], petsclib)

	return smooth
end 

"""
	smoothd::SNES = SNESFASCycleGetSmootherDown(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the down smoother on a particular cycle level.

Logically Collective

Input Parameter:
- `snes` - `SNESFAS` obtained with `SNESFASGetCycleSNES()`

Output Parameter:
- `smoothd` - the smoother

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASCycleGetSmootherUp()`, `SNESFASCycleGetSmoother()`, `SNESFASGetCycleSNES()`

# External Links
$(_doc_external("SNESFAS/SNESFASCycleGetSmootherDown"))
"""
function SNESFASCycleGetSmootherDown(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESFASCycleGetSmootherDown: no generated method for these argument types")
end

@for_petsc function SNESFASCycleGetSmootherDown(petsclib::$UnionPetscLib, snes::AbstractSNES )
	smoothd_ = Ref{CSNES}()

    @chk ccall(
               (:SNESFASCycleGetSmootherDown, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CSNES}),
               snes, smoothd_,
              )

	smoothd = SNES(smoothd_[], petsclib)

	return smoothd
end 

"""
	smoothu::SNES = SNESFASCycleGetSmootherUp(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the up smoother on a particular cycle level.

Logically Collective

Input Parameter:
- `snes` - the `SNESFAS` obtained with `SNESFASGetCycleSNES()`

Output Parameter:
- `smoothu` - the smoother

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASCycleGetSmoother()`, `SNESFASCycleGetSmootherDown()`, `SNESFASGetCycleSNES()`

# External Links
$(_doc_external("SNESFAS/SNESFASCycleGetSmootherUp"))
"""
function SNESFASCycleGetSmootherUp(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESFASCycleGetSmootherUp: no generated method for these argument types")
end

@for_petsc function SNESFASCycleGetSmootherUp(petsclib::$UnionPetscLib, snes::AbstractSNES )
	smoothu_ = Ref{CSNES}()

    @chk ccall(
               (:SNESFASCycleGetSmootherUp, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CSNES}),
               snes, smoothu_,
              )

	smoothu = SNES(smoothu_[], petsclib)

	return smoothu
end 

"""
	flg::PetscBool = SNESFASCycleIsFine(petsclib::PetscLibType, snes::AbstractSNES) 
Determines if a given `SNES` is the finest level in a `SNESFAS`

Logically Collective

Input Parameter:
- `snes` - the `SNESFAS` context obtained with `SNESFASGetCycleSNES()`

Output Parameter:
- `flg` - indicates if this is the fine level or not

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASSetLevels()`

# External Links
$(_doc_external("SNESFAS/SNESFASCycleIsFine"))
"""
function SNESFASCycleIsFine(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESFASCycleIsFine: no generated method for these argument types")
end

@for_petsc function SNESFASCycleIsFine(petsclib::$UnionPetscLib, snes::AbstractSNES )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:SNESFASCycleIsFine, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{PetscBool}),
               snes, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	SNESFASCycleSetCycles(petsclib::PetscLibType, snes::AbstractSNES, cycles::PetscInt) 
Sets the number of cycles for all levels in a `SNESFAS`

Logically Collective

Input Parameters:
- `snes`   - the `SNESFAS` nonlinear multigrid context
- `cycles` - the number of cycles -- 1 for V-cycle, 2 for W-cycle

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASSetCycles()`

# External Links
$(_doc_external("SNESFAS/SNESFASCycleSetCycles"))
"""
function SNESFASCycleSetCycles(petsclib::PetscLibType, snes::AbstractSNES, cycles::Integer)
    error("SNESFASCycleSetCycles: no generated method for these argument types")
end

@for_petsc function SNESFASCycleSetCycles(petsclib::$UnionPetscLib, snes::AbstractSNES, cycles::$PetscInt )

    @chk ccall(
               (:SNESFASCycleSetCycles, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt),
               snes, cycles,
              )


	return nothing
end 

"""
	total::PetscBool = SNESFASFullGetTotal(petsclib::PetscLibType, snes::AbstractSNES) 
Use total residual restriction and total interpolation on the initial down and up sweep of full FAS cycles

Logically Collective

Input Parameter:
- `snes` - the `SNESFAS` nonlinear multigrid context

Output Parameter:
- `total` - whether to use total restriction / interpolatiaon or not (the alternative is defect restriction and correction interpolation)

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASSetNumberSmoothUp()`, `DMInterpolateSolution()`, `SNESFullSetTotal()`

# External Links
$(_doc_external("SNESFAS/SNESFASFullGetTotal"))
"""
function SNESFASFullGetTotal(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESFASFullGetTotal: no generated method for these argument types")
end

@for_petsc function SNESFASFullGetTotal(petsclib::$UnionPetscLib, snes::AbstractSNES )
	total_ = Ref{PetscBool}()

    @chk ccall(
               (:SNESFASFullGetTotal, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{PetscBool}),
               snes, total_,
              )

	total = total_[]

	return total
end 

"""
	SNESFASFullSetDownSweep(petsclib::PetscLibType, snes::AbstractSNES, swp::PetscBool) 
Smooth during the initial downsweep for `SNESFAS`

Logically Collective

Input Parameters:
- `snes` - the `SNESFAS` nonlinear multigrid context
- `swp`  - whether to downsweep or not

Options Database Key:
- `-snes_fas_full_downsweep` - Sets whether to smooth on the initial downsweep

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASSetNumberSmoothUp()`

# External Links
$(_doc_external("SNESFAS/SNESFASFullSetDownSweep"))
"""
function SNESFASFullSetDownSweep(petsclib::PetscLibType, snes::AbstractSNES, swp::PetscBool)
    error("SNESFASFullSetDownSweep: no generated method for these argument types")
end

@for_petsc function SNESFASFullSetDownSweep(petsclib::$UnionPetscLib, snes::AbstractSNES, swp::PetscBool )

    @chk ccall(
               (:SNESFASFullSetDownSweep, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, swp,
              )


	return nothing
end 

"""
	SNESFASFullSetTotal(petsclib::PetscLibType, snes::AbstractSNES, total::PetscBool) 
Use total residual restriction and total interpolation on the initial down and up sweep of full `SNESFAS` cycles

Logically Collective

Input Parameters:
- `snes`  - the `SNESFAS`  nonlinear multigrid context
- `total` - whether to use total restriction / interpolatiaon or not (the alternative is defect restriction and correction interpolation)

Options Database Key:
- `-snes_fas_full_total` - Use total restriction and interpolation on the initial down and up sweeps for the full `SNESFAS` cycle

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASSetNumberSmoothUp()`, `DMInterpolateSolution()`

# External Links
$(_doc_external("SNESFAS/SNESFASFullSetTotal"))
"""
function SNESFASFullSetTotal(petsclib::PetscLibType, snes::AbstractSNES, total::PetscBool)
    error("SNESFASFullSetTotal: no generated method for these argument types")
end

@for_petsc function SNESFASFullSetTotal(petsclib::$UnionPetscLib, snes::AbstractSNES, total::PetscBool )

    @chk ccall(
               (:SNESFASFullSetTotal, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, total,
              )


	return nothing
end 

"""
	SNESFASGalerkinFunctionDefault(petsclib::PetscLibType, snes::AbstractSNES, X::AbstractPetscVec, F::AbstractPetscVec, ctx::Ptr{Cvoid}) 
Computes the Galerkin FAS function

Collective

Input Parameters:
- `snes` - the `SNESFAS` nonlinear solver context
- `X`    - input vector
- `ctx`  - the application context

Output Parameter:
- `F` - output vector

Level: developer

See also: `SNES`, `SNESFAS`, `SNESFASGetGalerkin()`, `SNESFASSetGalerkin()`

# External Links
$(_doc_external("SNESFAS/SNESFASGalerkinFunctionDefault"))
"""
function SNESFASGalerkinFunctionDefault(petsclib::PetscLibType, snes::AbstractSNES, X::AbstractPetscVec, F::AbstractPetscVec, ctx::Ptr{Cvoid})
    error("SNESFASGalerkinFunctionDefault: no generated method for these argument types")
end

@for_petsc function SNESFASGalerkinFunctionDefault(petsclib::$UnionPetscLib, snes::AbstractSNES, X::AbstractPetscVec, F::AbstractPetscVec, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESFASGalerkinFunctionDefault, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec, Ptr{Cvoid}),
               snes, X, F, ctx,
              )


	return nothing
end 

"""
	coarse::SNES = SNESFASGetCoarseSolve(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the coarsest level solver.

Input Parameter:
- `snes` - the `SNESFAS` nonlinear multigrid context

Output Parameter:
- `coarse` - the coarse-level solver

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASSetInjection()`, `SNESFASSetRestriction()`

# External Links
$(_doc_external("SNESFAS/SNESFASGetCoarseSolve"))
"""
function SNESFASGetCoarseSolve(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESFASGetCoarseSolve: no generated method for these argument types")
end

@for_petsc function SNESFASGetCoarseSolve(petsclib::$UnionPetscLib, snes::AbstractSNES )
	coarse_ = Ref{CSNES}()

    @chk ccall(
               (:SNESFASGetCoarseSolve, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CSNES}),
               snes, coarse_,
              )

	coarse = SNES(coarse_[], petsclib)

	return coarse
end 

"""
	lsnes::SNES = SNESFASGetCycleSNES(petsclib::PetscLibType, snes::AbstractSNES, level::PetscInt) 
Gets the `SNES` corresponding to a particular level of the `SNESFAS` hierarchy

Input Parameters:
- `snes`  - the `SNES` nonlinear multigrid context
- `level` - the level to get

Output Parameter:
- `lsnes` - the `SNES` for the requested level

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASSetLevels()`, `SNESFASGetLevels()`

# External Links
$(_doc_external("SNESFAS/SNESFASGetCycleSNES"))
"""
function SNESFASGetCycleSNES(petsclib::PetscLibType, snes::AbstractSNES, level::Integer)
    error("SNESFASGetCycleSNES: no generated method for these argument types")
end

@for_petsc function SNESFASGetCycleSNES(petsclib::$UnionPetscLib, snes::AbstractSNES, level::$PetscInt )
	lsnes_ = Ref{CSNES}()

    @chk ccall(
               (:SNESFASGetCycleSNES, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, Ptr{CSNES}),
               snes, level, lsnes_,
              )

	lsnes = SNES(lsnes_[], petsclib)

	return lsnes
end 

"""
	flg::PetscBool = SNESFASGetGalerkin(petsclib::PetscLibType, snes::AbstractSNES) 
Gets if the coarse problems are formed by projection to the fine problem

Not Collective but the result would be the same on all MPI processes

Input Parameter:
- `snes` - the `SNESFAS` nonlinear solver context

Output Parameter:
- `flg` - `PETSC_TRUE` if the coarse problem is formed by projection

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASSetLevels()`, `SNESFASSetGalerkin()`

# External Links
$(_doc_external("SNESFAS/SNESFASGetGalerkin"))
"""
function SNESFASGetGalerkin(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESFASGetGalerkin: no generated method for these argument types")
end

@for_petsc function SNESFASGetGalerkin(petsclib::$UnionPetscLib, snes::AbstractSNES )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:SNESFASGetGalerkin, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{PetscBool}),
               snes, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	mat::PetscMat = SNESFASGetInjection(petsclib::PetscLibType, snes::AbstractSNES, level::PetscInt) 
Gets the matrix used to calculate the
injection from l-1 to the lth level

Input Parameters:
- `snes`  - the `SNESFAS` nonlinear multigrid context
- `level` - the level (0 is coarsest) to supply [do not supply 0]

Output Parameter:
- `mat` - the injection operator

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASSetInjection()`, `SNESFASGetRestriction()`, `SNESFASGetInterpolation()`, `SNESFASGetRScale()`

# External Links
$(_doc_external("SNESFAS/SNESFASGetInjection"))
"""
function SNESFASGetInjection(petsclib::PetscLibType, snes::AbstractSNES, level::Integer)
    error("SNESFASGetInjection: no generated method for these argument types")
end

@for_petsc function SNESFASGetInjection(petsclib::$UnionPetscLib, snes::AbstractSNES, level::$PetscInt )
	mat_ = Ref{CMat}()

    @chk ccall(
               (:SNESFASGetInjection, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, Ptr{CMat}),
               snes, level, mat_,
              )

	mat = PetscMat(mat_[], petsclib)

	return mat
end 

"""
	mat::PetscMat = SNESFASGetInterpolation(petsclib::PetscLibType, snes::AbstractSNES, level::PetscInt) 
Gets the matrix used to calculate the
interpolation from l-1 to the lth level

Input Parameters:
- `snes`  - the `SNESFAS` nonlinear multigrid context
- `level` - the level (0 is coarsest) to supply [do not supply 0]

Output Parameter:
- `mat` - the interpolation operator

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASSetInterpolation()`, `SNESFASGetInjection()`, `SNESFASGetRestriction()`, `SNESFASGetRScale()`

# External Links
$(_doc_external("SNESFAS/SNESFASGetInterpolation"))
"""
function SNESFASGetInterpolation(petsclib::PetscLibType, snes::AbstractSNES, level::Integer)
    error("SNESFASGetInterpolation: no generated method for these argument types")
end

@for_petsc function SNESFASGetInterpolation(petsclib::$UnionPetscLib, snes::AbstractSNES, level::$PetscInt )
	mat_ = Ref{CMat}()

    @chk ccall(
               (:SNESFASGetInterpolation, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, Ptr{CMat}),
               snes, level, mat_,
              )

	mat = PetscMat(mat_[], petsclib)

	return mat
end 

"""
	levels::PetscInt = SNESFASGetLevels(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the number of levels in a `SNESFAS`, including fine and coarse grids

Input Parameter:
- `snes` - the `SNES` nonlinear solver context of `SNESType` `SNESFAS`

Output Parameter:
- `levels` - the number of levels

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASSetLevels()`, `PCMGGetLevels()`

# External Links
$(_doc_external("SNESFAS/SNESFASGetLevels"))
"""
function SNESFASGetLevels(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESFASGetLevels: no generated method for these argument types")
end

@for_petsc function SNESFASGetLevels(petsclib::$UnionPetscLib, snes::AbstractSNES )
	levels_ = Ref{$PetscInt}()

    @chk ccall(
               (:SNESFASGetLevels, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscInt}),
               snes, levels_,
              )

	levels = levels_[]

	return levels
end 

"""
	mat::PetscMat = SNESFASGetRestriction(petsclib::PetscLibType, snes::AbstractSNES, level::PetscInt) 
Gets the matrix used to calculate the
restriction from l to the l-1th level

Input Parameters:
- `snes`  - the `SNESFAS` nonlinear multigrid context
- `level` - the level (0 is coarsest) to supply [do not supply 0]

Output Parameter:
- `mat` - the interpolation operator

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASSetRestriction()`, `SNESFASGetInjection()`, `SNESFASGetInterpolation()`, `SNESFASGetRScale()`

# External Links
$(_doc_external("SNESFAS/SNESFASGetRestriction"))
"""
function SNESFASGetRestriction(petsclib::PetscLibType, snes::AbstractSNES, level::Integer)
    error("SNESFASGetRestriction: no generated method for these argument types")
end

@for_petsc function SNESFASGetRestriction(petsclib::$UnionPetscLib, snes::AbstractSNES, level::$PetscInt )
	mat_ = Ref{CMat}()

    @chk ccall(
               (:SNESFASGetRestriction, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, Ptr{CMat}),
               snes, level, mat_,
              )

	mat = PetscMat(mat_[], petsclib)

	return mat
end 

"""
	smooth::SNES = SNESFASGetSmoother(petsclib::PetscLibType, snes::AbstractSNES, level::PetscInt) 
Gets the default smoother on a level.

Input Parameters:
- `snes`  - the `SNESFAS` nonlinear multigrid context
- `level` - the level (0 is coarsest) to supply

Output Parameter:
- `smooth` - the smoother

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASSetInjection()`, `SNESFASSetRestriction()`

# External Links
$(_doc_external("SNESFAS/SNESFASGetSmoother"))
"""
function SNESFASGetSmoother(petsclib::PetscLibType, snes::AbstractSNES, level::Integer)
    error("SNESFASGetSmoother: no generated method for these argument types")
end

@for_petsc function SNESFASGetSmoother(petsclib::$UnionPetscLib, snes::AbstractSNES, level::$PetscInt )
	smooth_ = Ref{CSNES}()

    @chk ccall(
               (:SNESFASGetSmoother, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, Ptr{CSNES}),
               snes, level, smooth_,
              )

	smooth = SNES(smooth_[], petsclib)

	return smooth
end 

"""
	smooth::SNES = SNESFASGetSmootherDown(petsclib::PetscLibType, snes::AbstractSNES, level::PetscInt) 
Gets the downsmoother on a level.

Input Parameters:
- `snes`  - the `SNESFAS` nonlinear multigrid context
- `level` - the level (0 is coarsest) to supply

Output Parameter:
- `smooth` - the smoother

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASSetInjection()`, `SNESFASSetRestriction()`

# External Links
$(_doc_external("SNESFAS/SNESFASGetSmootherDown"))
"""
function SNESFASGetSmootherDown(petsclib::PetscLibType, snes::AbstractSNES, level::Integer)
    error("SNESFASGetSmootherDown: no generated method for these argument types")
end

@for_petsc function SNESFASGetSmootherDown(petsclib::$UnionPetscLib, snes::AbstractSNES, level::$PetscInt )
	smooth_ = Ref{CSNES}()

    @chk ccall(
               (:SNESFASGetSmootherDown, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, Ptr{CSNES}),
               snes, level, smooth_,
              )

	smooth = SNES(smooth_[], petsclib)

	return smooth
end 

"""
	smooth::SNES = SNESFASGetSmootherUp(petsclib::PetscLibType, snes::AbstractSNES, level::PetscInt) 
Gets the upsmoother on a level.

Input Parameters:
- `snes`  - the `SNESFAS` nonlinear multigrid context
- `level` - the level (0 is coarsest)

Output Parameter:
- `smooth` - the smoother

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASSetInjection()`, `SNESFASSetRestriction()`

# External Links
$(_doc_external("SNESFAS/SNESFASGetSmootherUp"))
"""
function SNESFASGetSmootherUp(petsclib::PetscLibType, snes::AbstractSNES, level::Integer)
    error("SNESFASGetSmootherUp: no generated method for these argument types")
end

@for_petsc function SNESFASGetSmootherUp(petsclib::$UnionPetscLib, snes::AbstractSNES, level::$PetscInt )
	smooth_ = Ref{CSNES}()

    @chk ccall(
               (:SNESFASGetSmootherUp, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, Ptr{CSNES}),
               snes, level, smooth_,
              )

	smooth = SNES(smooth_[], petsclib)

	return smooth
end 

"""
	fastype::SNESFASType = SNESFASGetType(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the update and correction type used for `SNESFAS`.

Logically Collective

Input Parameter:
- `snes` - `SNESFAS` context

Output Parameter:
- `fastype` - `SNES_FAS_ADDITIVE`, `SNES_FAS_MULTIPLICATIVE`, `SNES_FAS_FULL`, or `SNES_FAS_KASKADE`

Level: intermediate

See also: `SNES`, `SNESFAS`, `PCMGSetType()`, `SNESFASSetType()`, `SNES_FAS_ADDITIVE`, `SNES_FAS_MULTIPLICATIVE`, `SNES_FAS_FULL`, `SNES_FAS_KASKADE`

# External Links
$(_doc_external("SNESFAS/SNESFASGetType"))
"""
function SNESFASGetType(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESFASGetType: no generated method for these argument types")
end

@for_petsc function SNESFASGetType(petsclib::$UnionPetscLib, snes::AbstractSNES )
	fastype_ = Ref{SNESFASType}()

    @chk ccall(
               (:SNESFASGetType, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{SNESFASType}),
               snes, fastype_,
              )

	fastype = fastype_[]

	return fastype
end 

"""
	SNESFASRestrict(petsclib::PetscLibType, fine::AbstractSNES, Xfine::AbstractPetscVec, Xcoarse::AbstractPetscVec) 
restrict a `Vec` to the next coarser level

Collective

Input Parameters:
- `fine`  - `SNES` from which to restrict
- `Xfine` - vector to restrict

Output Parameter:
- `Xcoarse` - result of restriction

Level: developer

See also: `SNES`, `SNESFAS`, `SNESFASSetRestriction()`, `SNESFASSetInjection()`, `SNESFASCreateCoarseVec()`

# External Links
$(_doc_external("SNESFAS/SNESFASRestrict"))
"""
function SNESFASRestrict(petsclib::PetscLibType, fine::AbstractSNES, Xfine::AbstractPetscVec, Xcoarse::AbstractPetscVec)
    error("SNESFASRestrict: no generated method for these argument types")
end

@for_petsc function SNESFASRestrict(petsclib::$UnionPetscLib, fine::AbstractSNES, Xfine::AbstractPetscVec, Xcoarse::AbstractPetscVec )

    @chk ccall(
               (:SNESFASRestrict, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec),
               fine, Xfine, Xcoarse,
              )


	return nothing
end 

"""
	SNESFASSetContinuation(petsclib::PetscLibType, snes::AbstractSNES, continuation::PetscBool) 
Sets the `SNESFAS` cycle to default to using exact Newton solves on the upsweep

Logically Collective

Input Parameters:
- `snes`         - the `SNESFAS` nonlinear multigrid context
- `continuation` - whether to use continuation

Options Database Key:
- `-snes_fas_continuation` - sets continuation to true

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASSetNumberSmoothUp()`

# External Links
$(_doc_external("SNESFAS/SNESFASSetContinuation"))
"""
function SNESFASSetContinuation(petsclib::PetscLibType, snes::AbstractSNES, continuation::PetscBool)
    error("SNESFASSetContinuation: no generated method for these argument types")
end

@for_petsc function SNESFASSetContinuation(petsclib::$UnionPetscLib, snes::AbstractSNES, continuation::PetscBool )

    @chk ccall(
               (:SNESFASSetContinuation, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, continuation,
              )


	return nothing
end 

"""
	SNESFASSetCycles(petsclib::PetscLibType, snes::AbstractSNES, cycles::PetscInt) 
Sets the number of `SNESFAS` multigrid cycles to use each time a grid is visited.  Use `SNESFASSetCyclesOnLevel()` for more
complicated cycling.

Logically Collective

Input Parameters:
- `snes`   - the `SNESFAS` nonlinear multigrid context
- `cycles` - the number of cycles -- 1 for V-cycle, 2 for W-cycle

Options Database Key:
- `-snes_fas_cycles (1|2)` - 1 for V-cycle, 2 for W-cycle

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASSetCyclesOnLevel()`

# External Links
$(_doc_external("SNESFAS/SNESFASSetCycles"))
"""
function SNESFASSetCycles(petsclib::PetscLibType, snes::AbstractSNES, cycles::Integer)
    error("SNESFASSetCycles: no generated method for these argument types")
end

@for_petsc function SNESFASSetCycles(petsclib::$UnionPetscLib, snes::AbstractSNES, cycles::$PetscInt )

    @chk ccall(
               (:SNESFASSetCycles, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt),
               snes, cycles,
              )


	return nothing
end 

"""
	SNESFASSetGalerkin(petsclib::PetscLibType, snes::AbstractSNES, flg::PetscBool) 
Sets coarse problems as formed by projection to the fine problem

Logically Collective

Input Parameters:
- `snes` - the `SNESFAS` nonlinear solver context
- `flg`  - `PETSC_TRUE` to use the projection process

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASSetLevels()`, `SNESFASGetGalerkin()`

# External Links
$(_doc_external("SNESFAS/SNESFASSetGalerkin"))
"""
function SNESFASSetGalerkin(petsclib::PetscLibType, snes::AbstractSNES, flg::PetscBool)
    error("SNESFASSetGalerkin: no generated method for these argument types")
end

@for_petsc function SNESFASSetGalerkin(petsclib::$UnionPetscLib, snes::AbstractSNES, flg::PetscBool )

    @chk ccall(
               (:SNESFASSetGalerkin, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, flg,
              )


	return nothing
end 

"""
	SNESFASSetInjection(petsclib::PetscLibType, snes::AbstractSNES, level::PetscInt, mat::AbstractPetscMat) 
Sets the matrix to be used to inject the solution
from `level` to `level-1`.

Input Parameters:
- `snes`  - the `SNESFAS` nonlinear multigrid context
- `mat`   - the injection matrix
- `level` - the level (0 is coarsest) to supply [Do not supply 0]

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASSetInterpolation()`, `SNESFASSetRestriction()`

# External Links
$(_doc_external("SNESFAS/SNESFASSetInjection"))
"""
function SNESFASSetInjection(petsclib::PetscLibType, snes::AbstractSNES, level::Integer, mat::AbstractPetscMat)
    error("SNESFASSetInjection: no generated method for these argument types")
end

@for_petsc function SNESFASSetInjection(petsclib::$UnionPetscLib, snes::AbstractSNES, level::$PetscInt, mat::AbstractPetscMat )

    @chk ccall(
               (:SNESFASSetInjection, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, CMat),
               snes, level, mat,
              )


	return nothing
end 

"""
	SNESFASSetInterpolation(petsclib::PetscLibType, snes::AbstractSNES, level::PetscInt, mat::AbstractPetscMat) 
Sets the `Mat` to be used to apply the
interpolation from l-1 to the lth level

Input Parameters:
- `snes`  - the `SNESFAS` nonlinear multigrid context
- `mat`   - the interpolation operator
- `level` - the level (0 is coarsest) to supply [do not supply 0]

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASSetInjection()`, `SNESFASSetRestriction()`, `SNESFASSetRScale()`

# External Links
$(_doc_external("SNESFAS/SNESFASSetInterpolation"))
"""
function SNESFASSetInterpolation(petsclib::PetscLibType, snes::AbstractSNES, level::Integer, mat::AbstractPetscMat)
    error("SNESFASSetInterpolation: no generated method for these argument types")
end

@for_petsc function SNESFASSetInterpolation(petsclib::$UnionPetscLib, snes::AbstractSNES, level::$PetscInt, mat::AbstractPetscMat )

    @chk ccall(
               (:SNESFASSetInterpolation, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, CMat),
               snes, level, mat,
              )


	return nothing
end 

"""
	comms::MPI_Comm = SNESFASSetLevels(petsclib::PetscLibType, snes::AbstractSNES, levels::PetscInt) 
Sets the number of levels to use with `SNESFAS`.
Must be called before any other `SNESFAS` routine.

Input Parameters:
- `snes`   - the `SNES` context of `SNESType` `SNESFAS`
- `levels` - the number of levels
- `comms`  - optional communicators for each level; this is to allow solving the coarser
problems on smaller sets of processors.

Level: intermediate

See also: `SNES`, `SNESFAS`, `SNESFASGetLevels()`

# External Links
$(_doc_external("SNESFAS/SNESFASSetLevels"))
"""
function SNESFASSetLevels(petsclib::PetscLibType, snes::AbstractSNES, levels::Integer)
    error("SNESFASSetLevels: no generated method for these argument types")
end

@for_petsc function SNESFASSetLevels(petsclib::$UnionPetscLib, snes::AbstractSNES, levels::$PetscInt )
	comms_ = Ref{MPI.MPI_Comm}()

    @chk ccall(
               (:SNESFASSetLevels, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, Ptr{MPI.MPI_Comm}),
               snes, levels, comms_,
              )

	comms = MPI.Comm(comms_[])

	return comms
end 

"""
	SNESFASSetLog(petsclib::PetscLibType, snes::AbstractSNES, flg::PetscBool) 
Sets or unsets time logging for various `SNESFAS` stages on all levels

Logically Collective

Input Parameters:
- `snes` - the `SNESFAS` context
- `flg`  - whether to log or not

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASSetMonitor()`

# External Links
$(_doc_external("SNESFAS/SNESFASSetLog"))
"""
function SNESFASSetLog(petsclib::PetscLibType, snes::AbstractSNES, flg::PetscBool)
    error("SNESFASSetLog: no generated method for these argument types")
end

@for_petsc function SNESFASSetLog(petsclib::$UnionPetscLib, snes::AbstractSNES, flg::PetscBool )

    @chk ccall(
               (:SNESFASSetLog, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, flg,
              )


	return nothing
end 

"""
	SNESFASSetMonitor(petsclib::PetscLibType, snes::AbstractSNES, vf::Vector{PetscViewerAndFormat}, flg::PetscBool) 
Sets the method-specific cycle monitoring

Logically Collective

Input Parameters:
- `snes` - the `SNESFAS` context
- `vf`   - viewer and format structure (may be `NULL` if `flg` is `PETSC_FALSE`)
- `flg`  - monitor or not

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESMonitorSet()`, `SNESFASSetCyclesOnLevel()`

# External Links
$(_doc_external("SNESFAS/SNESFASSetMonitor"))
"""
function SNESFASSetMonitor(petsclib::PetscLibType, snes::AbstractSNES, vf::Vector{PetscViewerAndFormat}, flg::PetscBool)
    error("SNESFASSetMonitor: no generated method for these argument types")
end

@for_petsc function SNESFASSetMonitor(petsclib::$UnionPetscLib, snes::AbstractSNES, vf::Vector{PetscViewerAndFormat}, flg::PetscBool )

    @chk ccall(
               (:SNESFASSetMonitor, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{PetscViewerAndFormat}, PetscBool),
               snes, vf, flg,
              )


	return nothing
end 

"""
	SNESFASSetNumberSmoothDown(petsclib::PetscLibType, snes::AbstractSNES, n::PetscInt) 
Sets the number of pre-smoothing steps to
use on all levels.

Logically Collective

Input Parameters:
- `snes` - the `SNESFAS` nonlinear multigrid context
- `n`    - the number of smoothing steps to use

Options Database Key:
- `-snes_fas_smoothdown n` - Sets number of pre-smoothing steps

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASSetNumberSmoothUp()`

# External Links
$(_doc_external("SNESFAS/SNESFASSetNumberSmoothDown"))
"""
function SNESFASSetNumberSmoothDown(petsclib::PetscLibType, snes::AbstractSNES, n::Integer)
    error("SNESFASSetNumberSmoothDown: no generated method for these argument types")
end

@for_petsc function SNESFASSetNumberSmoothDown(petsclib::$UnionPetscLib, snes::AbstractSNES, n::$PetscInt )

    @chk ccall(
               (:SNESFASSetNumberSmoothDown, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt),
               snes, n,
              )


	return nothing
end 

"""
	SNESFASSetNumberSmoothUp(petsclib::PetscLibType, snes::AbstractSNES, n::PetscInt) 
Sets the number of post-smoothing steps to
use on all levels.

Logically Collective

Input Parameters:
- `snes` - the `SNES` nonlinear multigrid context
- `n`    - the number of smoothing steps to use

Options Database Key:
- `-snes_fas_smoothup n` - Sets number of pre-smoothing steps

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASSetNumberSmoothDown()`

# External Links
$(_doc_external("SNESFAS/SNESFASSetNumberSmoothUp"))
"""
function SNESFASSetNumberSmoothUp(petsclib::PetscLibType, snes::AbstractSNES, n::Integer)
    error("SNESFASSetNumberSmoothUp: no generated method for these argument types")
end

@for_petsc function SNESFASSetNumberSmoothUp(petsclib::$UnionPetscLib, snes::AbstractSNES, n::$PetscInt )

    @chk ccall(
               (:SNESFASSetNumberSmoothUp, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt),
               snes, n,
              )


	return nothing
end 

"""
	SNESFASSetRScale(petsclib::PetscLibType, snes::AbstractSNES, level::PetscInt, rscale::AbstractPetscVec) 
Sets the scaling factor of the restriction
operator from level l to l-1.

Input Parameters:
- `snes`   - the `SNESFAS` nonlinear multigrid context
- `rscale` - the restriction scaling
- `level`  - the level (0 is coarsest) to supply [Do not supply 0]

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASSetInjection()`, `SNESFASSetRestriction()`

# External Links
$(_doc_external("SNESFAS/SNESFASSetRScale"))
"""
function SNESFASSetRScale(petsclib::PetscLibType, snes::AbstractSNES, level::Integer, rscale::AbstractPetscVec)
    error("SNESFASSetRScale: no generated method for these argument types")
end

@for_petsc function SNESFASSetRScale(petsclib::$UnionPetscLib, snes::AbstractSNES, level::$PetscInt, rscale::AbstractPetscVec )

    @chk ccall(
               (:SNESFASSetRScale, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, CVec),
               snes, level, rscale,
              )


	return nothing
end 

"""
	SNESFASSetRestriction(petsclib::PetscLibType, snes::AbstractSNES, level::PetscInt, mat::AbstractPetscMat) 
Sets the matrix to be used to restrict the defect
from level l to l-1.

Input Parameters:
- `snes`  - the `SNESFAS` nonlinear multigrid context
- `mat`   - the restriction matrix
- `level` - the level (0 is coarsest) to supply [Do not supply 0]

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESFASSetInterpolation()`, `SNESFASSetInjection()`

# External Links
$(_doc_external("SNESFAS/SNESFASSetRestriction"))
"""
function SNESFASSetRestriction(petsclib::PetscLibType, snes::AbstractSNES, level::Integer, mat::AbstractPetscMat)
    error("SNESFASSetRestriction: no generated method for these argument types")
end

@for_petsc function SNESFASSetRestriction(petsclib::$UnionPetscLib, snes::AbstractSNES, level::$PetscInt, mat::AbstractPetscMat )

    @chk ccall(
               (:SNESFASSetRestriction, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, CMat),
               snes, level, mat,
              )


	return nothing
end 

"""
	SNESFASSetType(petsclib::PetscLibType, snes::AbstractSNES, fastype::SNESFASType) 
Sets the update and correction type used for `SNESFAS`.

Logically Collective

Input Parameters:
- `snes`    - `SNESFAS` context
- `fastype` - `SNES_FAS_ADDITIVE`, `SNES_FAS_MULTIPLICATIVE`, `SNES_FAS_FULL`, or `SNES_FAS_KASKADE`

Level: intermediate

See also: `SNES`, `SNESFAS`, `PCMGSetType()`, `SNESFASGetType()`, `SNES_FAS_ADDITIVE`, `SNES_FAS_MULTIPLICATIVE`, `SNES_FAS_FULL`, `SNES_FAS_KASKADE`

# External Links
$(_doc_external("SNESFAS/SNESFASSetType"))
"""
function SNESFASSetType(petsclib::PetscLibType, snes::AbstractSNES, fastype::SNESFASType)
    error("SNESFASSetType: no generated method for these argument types")
end

@for_petsc function SNESFASSetType(petsclib::$UnionPetscLib, snes::AbstractSNES, fastype::SNESFASType )

    @chk ccall(
               (:SNESFASSetType, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESFASType),
               snes, fastype,
              )


	return nothing
end 

"""
	SNESFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the PETSc interface to the `SNES` package. It is
called from `PetscFinalize()`.

Level: developer

See also: `SNES`, `PetscFinalize()`

# External Links
$(_doc_external("SNES/SNESFinalizePackage"))
"""
function SNESFinalizePackage(petsclib::PetscLibType)
    error("SNESFinalizePackage: no generated method for these argument types")
end

@for_petsc function SNESFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:SNESFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	flg::PetscBool = SNESGetAlwaysComputesFinalResidual(petsclib::PetscLibType, snes::AbstractSNES) 
checks if the `SNES` always computes the residual at the final solution

Logically Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `flg` - `PETSC_TRUE` if the residual is computed

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESSolve()`, `SNESSetAlwaysComputesFinalResidual()`

# External Links
$(_doc_external("SNES/SNESGetAlwaysComputesFinalResidual"))
"""
function SNESGetAlwaysComputesFinalResidual(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetAlwaysComputesFinalResidual: no generated method for these argument types")
end

@for_petsc function SNESGetAlwaysComputesFinalResidual(petsclib::$UnionPetscLib, snes::AbstractSNES )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:SNESGetAlwaysComputesFinalResidual, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{PetscBool}),
               snes, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	ctx::Ptr{Cvoid} = SNESGetApplicationContext(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the user-defined context for the
nonlinear solvers set with `SNESGetApplicationContext()` or `SNESSetComputeApplicationContext()`

Not Collective

Input Parameter:
- `snes` - `SNES` context

Output Parameter:
- `ctx` - the application context

Level: intermediate

See also: `SNESSetApplicationContext()`, `SNESSetComputeApplicationContext()`

# External Links
$(_doc_external("SNES/SNESGetApplicationContext"))
"""
function SNESGetApplicationContext(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetApplicationContext: no generated method for these argument types")
end

@for_petsc function SNESGetApplicationContext(petsclib::$UnionPetscLib, snes::AbstractSNES )
	ctx_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:SNESGetApplicationContext, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Cvoid}),
               snes, ctx_,
              )

	ctx = ctx_[]

	return ctx
end 

"""
	flg::PetscBool = SNESGetCheckJacobianDomainError(petsclib::PetscLibType, snes::AbstractSNES) 
Get an indicator whether or not `SNES` is checking Jacobian domain errors after each Jacobian evaluation.

Logically Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `flg` - `PETSC_FALSE` indicates that it is not checking Jacobian domain errors after each Jacobian evaluation

Level: advanced

See also: `SNES`, `SNESCreate()`, `SNESSetFunction()`, `SNESFunctionFn`, `SNESSetFunctionDomainError()`, `SNESSetCheckJacobianDomainError()`

# External Links
$(_doc_external("SNES/SNESGetCheckJacobianDomainError"))
"""
function SNESGetCheckJacobianDomainError(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetCheckJacobianDomainError: no generated method for these argument types")
end

@for_petsc function SNESGetCheckJacobianDomainError(petsclib::$UnionPetscLib, snes::AbstractSNES )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:SNESGetCheckJacobianDomainError, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{PetscBool}),
               snes, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	reason::SNESConvergedReason = SNESGetConvergedReason(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the reason the `SNES` iteration was stopped, which may be due to convergence, divergence, or stagnation

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `reason` - negative value indicates diverged, positive value converged, see `SNESConvergedReason` for the individual convergence tests for complete lists

Options Database Key:
- `-snes_converged_reason` - prints the reason to standard out

Level: intermediate

See also: `SNESSolve()`, `SNESSetConvergenceTest()`, `SNESSetConvergedReason()`, `SNESConvergedReason`, `SNESGetConvergedReasonString()`

# External Links
$(_doc_external("SNES/SNESGetConvergedReason"))
"""
function SNESGetConvergedReason(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetConvergedReason: no generated method for these argument types")
end

@for_petsc function SNESGetConvergedReason(petsclib::$UnionPetscLib, snes::AbstractSNES )
	reason_ = Ref{SNESConvergedReason}()

    @chk ccall(
               (:SNESGetConvergedReason, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{SNESConvergedReason}),
               snes, reason_,
              )

	reason = reason_[]

	return reason
end 

"""
	strreason::String = SNESGetConvergedReasonString(petsclib::PetscLibType, snes::AbstractSNES) 
Return a human readable string for `SNESConvergedReason`

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `strreason` - a human readable string that describes `SNES` converged reason

Level: beginner

See also: `SNES`, `SNESGetConvergedReason()`

# External Links
$(_doc_external("SNES/SNESGetConvergedReasonString"))
"""
function SNESGetConvergedReasonString(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetConvergedReasonString: no generated method for these argument types")
end

@for_petsc function SNESGetConvergedReasonString(petsclib::$UnionPetscLib, snes::AbstractSNES )
	strreason_ = Ref{Ptr{Cchar}}()

    @chk ccall(
               (:SNESGetConvergedReasonString, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Ptr{Cchar}}),
               snes, strreason_,
              )

	strreason = unsafe_string(strreason_[])

	return strreason
end 

"""
	a::Ptr{PetscReal},its::Ptr{PetscInt},na::PetscInt = SNESGetConvergenceHistory(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the arrays used to hold the convergence history.

Not Collective

Input Parameter:
- `snes` - iterative context obtained from `SNESCreate()`

Output Parameters:
- `a`   - array to hold history, usually was set with `SNESSetConvergenceHistory()`
- `its` - integer array holds the number of linear iterations (or
negative if not converged) for each solve.
- `na`  - size of `a` and `its`

Level: intermediate

See also: `SNES`, `SNESSolve()`, `SNESSetConvergenceHistory()`

# External Links
$(_doc_external("SNES/SNESGetConvergenceHistory"))
"""
function SNESGetConvergenceHistory(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetConvergenceHistory: no generated method for these argument types")
end

@for_petsc function SNESGetConvergenceHistory(petsclib::$UnionPetscLib, snes::AbstractSNES )
	a_ = Ref{Ptr{$PetscReal}}()
	its_ = Ref{Ptr{$PetscInt}}()
	na_ = Ref{$PetscInt}()

    @chk ccall(
               (:SNESGetConvergenceHistory, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Ptr{$PetscReal}}, Ptr{Ptr{$PetscInt}}, Ptr{$PetscInt}),
               snes, a_, its_, na_,
              )

	a = a_[]
	its = its_[]
	na = na_[]

	return a,its,na
end 

"""
	dm::PetscDM = SNESGetDM(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the `DM` that may be used by some `SNES` nonlinear solvers/preconditioners

Not Collective but `dm` obtained is parallel on `snes`

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `dm` - the `DM`

Level: intermediate

See also: `DM`, `SNES`, `SNESSetDM()`, `KSPSetDM()`, `KSPGetDM()`

# External Links
$(_doc_external("SNES/SNESGetDM"))
"""
function SNESGetDM(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetDM: no generated method for these argument types")
end

@for_petsc function SNESGetDM(petsclib::$UnionPetscLib, snes::AbstractSNES )
	dm_ = Ref{CDM}()

    @chk ccall(
               (:SNESGetDM, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CDM}),
               snes, dm_,
              )

	dm = PetscDM(dm_[], petsclib)

	return dm
end 

"""
	divtol::PetscReal = SNESGetDivergenceTolerance(petsclib::PetscLibType, snes::AbstractSNES) 
Gets divergence tolerance used in divergence test.

Not Collective

Input Parameters:
- `snes`   - the `SNES` context
- `divtol` - divergence tolerance

Level: intermediate

See also: `SNES`, `SNESSetDivergenceTolerance()`

# External Links
$(_doc_external("SNES/SNESGetDivergenceTolerance"))
"""
function SNESGetDivergenceTolerance(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetDivergenceTolerance: no generated method for these argument types")
end

@for_petsc function SNESGetDivergenceTolerance(petsclib::$UnionPetscLib, snes::AbstractSNES )
	divtol_ = Ref{$PetscReal}()

    @chk ccall(
               (:SNESGetDivergenceTolerance, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscReal}),
               snes, divtol_,
              )

	divtol = divtol_[]

	return divtol
end 

"""
	flag::PetscBool = SNESGetErrorIfNotConverged(petsclib::PetscLibType, snes::AbstractSNES) 
Indicates if `SNESSolve()` will generate an error if the solver does not converge?

Not Collective

Input Parameter:
- `snes` - iterative context obtained from `SNESCreate()`

Output Parameter:
- `flag` - `PETSC_TRUE` if it will generate an error, else `PETSC_FALSE`

Level: intermediate

See also: `SNES`, `SNESSolve()`, `SNESSetErrorIfNotConverged()`, `KSPGetErrorIfNotConverged()`, `KSPSetErrorIfNotConverged()`

# External Links
$(_doc_external("SNES/SNESGetErrorIfNotConverged"))
"""
function SNESGetErrorIfNotConverged(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetErrorIfNotConverged: no generated method for these argument types")
end

@for_petsc function SNESGetErrorIfNotConverged(petsclib::$UnionPetscLib, snes::AbstractSNES )
	flag_ = Ref{PetscBool}()

    @chk ccall(
               (:SNESGetErrorIfNotConverged, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{PetscBool}),
               snes, flag_,
              )

	flag = flag_[]

	return flag
end 

"""
	force::PetscBool = SNESGetForceIteration(petsclib::PetscLibType, snes::AbstractSNES) 
Check whether or not `SNESSolve()` take at least one iteration regardless of the initial residual norm

Logically Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `force` - `PETSC_TRUE` requires at least one iteration.

Level: intermediate

See also: `SNES`, `SNESSetForceIteration()`, `SNESSetDivergenceTolerance()`

# External Links
$(_doc_external("SNES/SNESGetForceIteration"))
"""
function SNESGetForceIteration(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetForceIteration: no generated method for these argument types")
end

@for_petsc function SNESGetForceIteration(petsclib::$UnionPetscLib, snes::AbstractSNES )
	force_ = Ref{PetscBool}()

    @chk ccall(
               (:SNESGetForceIteration, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{PetscBool}),
               snes, force_,
              )

	force = force_[]

	return force
end 

"""
	r::PetscVec,f::Ptr{Cvoid},ctx::Ptr{Cvoid} = SNESGetFunction(petsclib::PetscLibType, snes::AbstractSNES) 
Returns the function that defines the nonlinear system set with `SNESSetFunction()`

Not Collective, but `r` is parallel if `snes` is parallel. Collective if `r` is requested, but has not been created yet.

Input Parameter:
- `snes` - the `SNES` context

Output Parameters:
- `r`   - the vector that is used to store residuals (or `NULL` if you don't want it)
- `f`   - the function (or `NULL` if you don't want it);  for calling sequence see `SNESFunctionFn`
- `ctx` - the function context (or `NULL` if you don't want it)

Level: advanced

See also: `SNES`, `SNESSolve()`, `SNESSetFunction()`, `SNESGetSolution()`, `SNESFunctionFn`

# External Links
$(_doc_external("SNES/SNESGetFunction"))
"""
function SNESGetFunction(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetFunction: no generated method for these argument types")
end

@for_petsc function SNESGetFunction(petsclib::$UnionPetscLib, snes::AbstractSNES )
	r_ = Ref{CVec}()
	f_ = Ref{Ptr{Cvoid}}()
	ctx_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:SNESGetFunction, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CVec}, Ptr{Ptr{Cvoid}}, Ptr{Cvoid}),
               snes, r_, f_, ctx_,
              )

	r = PetscVec(r_[], petsclib)
	f = f_[]
	ctx = ctx_[]

	return r,f,ctx
end 

"""
	norm::PetscReal = SNESGetFunctionNorm(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the last computed norm of the residual

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `norm` - the last computed residual norm

Level: developer

See also: `SNES`, `SNESSetNormSchedule()`, `SNESComputeFunction()`, `VecNorm()`, `SNESSetFunction()`, `SNESSetInitialFunction()`, `SNESNormSchedule`

# External Links
$(_doc_external("SNES/SNESGetFunctionNorm"))
"""
function SNESGetFunctionNorm(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetFunctionNorm: no generated method for these argument types")
end

@for_petsc function SNESGetFunctionNorm(petsclib::$UnionPetscLib, snes::AbstractSNES )
	norm_ = Ref{$PetscReal}()

    @chk ccall(
               (:SNESGetFunctionNorm, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscReal}),
               snes, norm_,
              )

	norm = norm_[]

	return norm
end 

"""
	type::SNESFunctionType = SNESGetFunctionType(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the `SNESFunctionType` used in convergence and monitoring set with `SNESSetFunctionType()`
of the SNES method.

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `type` - the type of the function evaluation, see `SNESSetFunctionType()`

Level: advanced

See also: `SNESSetFunctionType()`, `SNESFunctionType`, `SNESSetNormSchedule()`, `SNESComputeFunction()`, `VecNorm()`, `SNESSetFunction()`, `SNESSetInitialFunction()`, `SNESNormSchedule`

# External Links
$(_doc_external("SNES/SNESGetFunctionType"))
"""
function SNESGetFunctionType(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetFunctionType: no generated method for these argument types")
end

@for_petsc function SNESGetFunctionType(petsclib::$UnionPetscLib, snes::AbstractSNES )
	type_ = Ref{SNESFunctionType}()

    @chk ccall(
               (:SNESGetFunctionType, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{SNESFunctionType}),
               snes, type_,
              )

	type = type_[]

	return type
end 

"""
	steps::PetscInt = SNESGetGridSequence(petsclib::PetscLibType, snes::AbstractSNES) 
gets the number of steps of grid sequencing that `SNES` will do

Logically Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `steps` - the number of refinements to do, defaults to 0

Level: intermediate

See also: `SNESGetLagPreconditioner()`, `SNESSetLagJacobian()`, `SNESGetLagJacobian()`, `SNESSetGridSequence()`

# External Links
$(_doc_external("SNES/SNESGetGridSequence"))
"""
function SNESGetGridSequence(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetGridSequence: no generated method for these argument types")
end

@for_petsc function SNESGetGridSequence(petsclib::$UnionPetscLib, snes::AbstractSNES )
	steps_ = Ref{$PetscInt}()

    @chk ccall(
               (:SNESGetGridSequence, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscInt}),
               snes, steps_,
              )

	steps = steps_[]

	return steps
end 

"""
	iter::PetscInt = SNESGetIterationNumber(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the number of nonlinear iterations completed in the current or most recent `SNESSolve()`

Not Collective

Input Parameter:
- `snes` - `SNES` context

Output Parameter:
- `iter` - iteration number

Level: intermediate

See also: `SNES`, `SNESSolve()`, `SNESSetLagJacobian()`, `SNESGetLinearSolveIterations()`, `SNESSetMonitor()`

# External Links
$(_doc_external("SNES/SNESGetIterationNumber"))
"""
function SNESGetIterationNumber(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetIterationNumber: no generated method for these argument types")
end

@for_petsc function SNESGetIterationNumber(petsclib::$UnionPetscLib, snes::AbstractSNES )
	iter_ = Ref{$PetscInt}()

    @chk ccall(
               (:SNESGetIterationNumber, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscInt}),
               snes, iter_,
              )

	iter = iter_[]

	return iter
end 

"""
	Amat::PetscMat,Pmat::PetscMat,J::Ptr{Cvoid},ctx::Ptr{Cvoid} = SNESGetJacobian(petsclib::PetscLibType, snes::AbstractSNES) 
Returns the Jacobian matrix and optionally the user
provided context for evaluating the Jacobian.

Not Collective, but `Mat` object will be parallel if `SNES` is

Input Parameter:
- `snes` - the nonlinear solver context

Output Parameters:
- `Amat` - location to stash (approximate) Jacobian matrix (or `NULL`)
- `Pmat` - location to stash matrix used to compute the preconditioner (or `NULL`)
- `J`    - location to put Jacobian function (or `NULL`), for calling sequence see `SNESJacobianFn`
- `ctx`  - location to stash Jacobian ctx (or `NULL`)

Level: advanced

See also: `SNES`, `Mat`, `SNESSetJacobian()`, `SNESComputeJacobian()`, `SNESJacobianFn`, `SNESGetFunction()`

# External Links
$(_doc_external("SNES/SNESGetJacobian"))
"""
function SNESGetJacobian(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetJacobian: no generated method for these argument types")
end

@for_petsc function SNESGetJacobian(petsclib::$UnionPetscLib, snes::AbstractSNES )
	Amat_ = Ref{CMat}()
	Pmat_ = Ref{CMat}()
	J_ = Ref{Ptr{Cvoid}}()
	ctx_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:SNESGetJacobian, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CMat}, Ptr{CMat}, Ptr{Ptr{Cvoid}}, Ptr{Cvoid}),
               snes, Amat_, Pmat_, J_, ctx_,
              )

	Amat = PetscMat(Amat_[], petsclib)
	Pmat = PetscMat(Pmat_[], petsclib)
	J = J_[]
	ctx = ctx_[]

	return Amat,Pmat,J,ctx
end 

"""
	ksp::KSP = SNESGetKSP(petsclib::PetscLibType, snes::AbstractSNES) 
Returns the `KSP` context for a `SNES` solver.

Not Collective, but if `snes` is parallel, then `ksp` is parallel

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `ksp` - the `KSP` context

Level: beginner

See also: `SNES`, `KSP`, `PC`, `KSPGetPC()`, `SNESCreate()`, `KSPCreate()`, `SNESSetKSP()`

# External Links
$(_doc_external("SNES/SNESGetKSP"))
"""
function SNESGetKSP(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetKSP: no generated method for these argument types")
end

@for_petsc function SNESGetKSP(petsclib::$UnionPetscLib, snes::AbstractSNES )
	ksp_ = Ref{CKSP}()

    @chk ccall(
               (:SNESGetKSP, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CKSP}),
               snes, ksp_,
              )

	ksp = KSP(ksp_[], petsclib)

	return ksp
end 

"""
	lag::PetscInt = SNESGetLagJacobian(petsclib::PetscLibType, snes::AbstractSNES) 
Get how often the Jacobian is rebuilt. See `SNESGetLagPreconditioner()` to determine when the preconditioner is rebuilt

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `lag` - -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time
the Jacobian is built etc.

Level: intermediate

See also: `SNES`, `SNESSetLagJacobian()`, `SNESSetLagPreconditioner()`, `SNESGetLagPreconditioner()`, `SNESSetLagJacobianPersists()`, `SNESSetLagPreconditionerPersists()`

# External Links
$(_doc_external("SNES/SNESGetLagJacobian"))
"""
function SNESGetLagJacobian(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetLagJacobian: no generated method for these argument types")
end

@for_petsc function SNESGetLagJacobian(petsclib::$UnionPetscLib, snes::AbstractSNES )
	lag_ = Ref{$PetscInt}()

    @chk ccall(
               (:SNESGetLagJacobian, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscInt}),
               snes, lag_,
              )

	lag = lag_[]

	return lag
end 

"""
	lag::PetscInt = SNESGetLagPreconditioner(petsclib::PetscLibType, snes::AbstractSNES) 
Return how often the preconditioner is rebuilt

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `lag` - -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time
the Jacobian is built etc. -2 indicates rebuild preconditioner at next chance but then never rebuild after that

Level: intermediate

See also: `SNES`, `SNESSetLagPreconditioner()`, `SNESSetLagJacobianPersists()`, `SNESSetLagPreconditionerPersists()`

# External Links
$(_doc_external("SNES/SNESGetLagPreconditioner"))
"""
function SNESGetLagPreconditioner(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetLagPreconditioner: no generated method for these argument types")
end

@for_petsc function SNESGetLagPreconditioner(petsclib::$UnionPetscLib, snes::AbstractSNES )
	lag_ = Ref{$PetscInt}()

    @chk ccall(
               (:SNESGetLagPreconditioner, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscInt}),
               snes, lag_,
              )

	lag = lag_[]

	return lag
end 

"""
	linesearch::SNESLineSearch = SNESGetLineSearch(petsclib::PetscLibType, snes::AbstractSNES) 
Returns the line search associated with the `SNES`.

Not Collective

Input Parameter:
- `snes` - iterative context obtained from `SNESCreate()`

Output Parameter:
- `linesearch` - linesearch context

Level: beginner

See also: `SNESLineSearch`, `SNESSetLineSearch()`, `SNESLineSearchCreate()`, `SNESLineSearchSetFromOptions()`

# External Links
$(_doc_external("SNES/SNESGetLineSearch"))
"""
function SNESGetLineSearch(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetLineSearch: no generated method for these argument types")
end

@for_petsc function SNESGetLineSearch(petsclib::$UnionPetscLib, snes::AbstractSNES )
	linesearch_ = Ref{SNESLineSearch}()

    @chk ccall(
               (:SNESGetLineSearch, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{SNESLineSearch}),
               snes, linesearch_,
              )

	linesearch = linesearch_[]

	return linesearch
end 

"""
	nfails::PetscInt = SNESGetLinearSolveFailures(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the number of failed (non-converged)
linear solvers in the current or most recent `SNESSolve()`

Not Collective

Input Parameter:
- `snes` - `SNES` context

Output Parameter:
- `nfails` - number of failed solves

Options Database Key:
- `-snes_max_linear_solve_fail num` - The number of failures before the solve is terminated

Level: intermediate

See also: `SNESGetMaxLinearSolveFailures()`, `SNESGetLinearSolveIterations()`, `SNESSetMaxLinearSolveFailures()`

# External Links
$(_doc_external("SNES/SNESGetLinearSolveFailures"))
"""
function SNESGetLinearSolveFailures(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetLinearSolveFailures: no generated method for these argument types")
end

@for_petsc function SNESGetLinearSolveFailures(petsclib::$UnionPetscLib, snes::AbstractSNES )
	nfails_ = Ref{$PetscInt}()

    @chk ccall(
               (:SNESGetLinearSolveFailures, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscInt}),
               snes, nfails_,
              )

	nfails = nfails_[]

	return nfails
end 

"""
	lits::PetscInt = SNESGetLinearSolveIterations(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the total number of linear iterations
used by the nonlinear solver in the most recent `SNESSolve()`

Not Collective

Input Parameter:
- `snes` - `SNES` context

Output Parameter:
- `lits` - number of linear iterations

Level: intermediate

See also: `SNES`, `SNESGetIterationNumber()`, `SNESGetLinearSolveFailures()`, `SNESGetMaxLinearSolveFailures()`, `SNESSetCountersReset()`

# External Links
$(_doc_external("SNES/SNESGetLinearSolveIterations"))
"""
function SNESGetLinearSolveIterations(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetLinearSolveIterations: no generated method for these argument types")
end

@for_petsc function SNESGetLinearSolveIterations(petsclib::$UnionPetscLib, snes::AbstractSNES )
	lits_ = Ref{$PetscInt}()

    @chk ccall(
               (:SNESGetLinearSolveIterations, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscInt}),
               snes, lits_,
              )

	lits = lits_[]

	return lits
end 

"""
	maxFails::PetscInt = SNESGetMaxLinearSolveFailures(petsclib::PetscLibType, snes::AbstractSNES) 
gets the maximum number of linear solve failures that
are allowed before `SNES` returns as unsuccessful

Not Collective

Input Parameter:
- `snes` - `SNES` context

Output Parameter:
- `maxFails` - maximum of unsuccessful solves allowed

Level: intermediate

See also: `SNESSetErrorIfNotConverged()`, `SNESGetLinearSolveFailures()`, `SNESGetLinearSolveIterations()`, `SNESSetMaxLinearSolveFailures()`

# External Links
$(_doc_external("SNES/SNESGetMaxLinearSolveFailures"))
"""
function SNESGetMaxLinearSolveFailures(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetMaxLinearSolveFailures: no generated method for these argument types")
end

@for_petsc function SNESGetMaxLinearSolveFailures(petsclib::$UnionPetscLib, snes::AbstractSNES )
	maxFails_ = Ref{$PetscInt}()

    @chk ccall(
               (:SNESGetMaxLinearSolveFailures, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscInt}),
               snes, maxFails_,
              )

	maxFails = maxFails_[]

	return maxFails
end 

"""
	maxFails::PetscInt = SNESGetMaxNonlinearStepFailures(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the maximum number of unsuccessful steps
attempted by the nonlinear solver before it gives up and returns unconverged or generates an error

Not Collective

Input Parameter:
- `snes` - `SNES` context

Output Parameter:
- `maxFails` - maximum of unsuccessful steps

Level: intermediate

See also: `SNESSetErrorIfNotConverged()`, `SNESGetMaxLinearSolveFailures()`, `SNESGetLinearSolveIterations()`, `SNESSetMaxLinearSolveFailures()`, `SNESGetLinearSolveFailures()`,
`SNESSetMaxNonlinearStepFailures()`, `SNESGetNonlinearStepFailures()`

# External Links
$(_doc_external("SNES/SNESGetMaxNonlinearStepFailures"))
"""
function SNESGetMaxNonlinearStepFailures(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetMaxNonlinearStepFailures: no generated method for these argument types")
end

@for_petsc function SNESGetMaxNonlinearStepFailures(petsclib::$UnionPetscLib, snes::AbstractSNES )
	maxFails_ = Ref{$PetscInt}()

    @chk ccall(
               (:SNESGetMaxNonlinearStepFailures, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscInt}),
               snes, maxFails_,
              )

	maxFails = maxFails_[]

	return maxFails
end 

"""
	f::Ptr{Cvoid},ctx::Ptr{Cvoid} = SNESGetNGS(petsclib::PetscLibType, snes::AbstractSNES) 
Returns the function and context set with `SNESSetNGS()`

Input Parameter:
- `snes` - the `SNES` context

Output Parameters:
- `f`   - the function (or `NULL`) see `SNESNGSFn` for calling sequence
- `ctx` - the function context (or `NULL`)

Level: advanced

See also: `SNESSetNGS()`, `SNESGetFunction()`, `SNESNGSFn`

# External Links
$(_doc_external("SNES/SNESGetNGS"))
"""
function SNESGetNGS(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetNGS: no generated method for these argument types")
end

@for_petsc function SNESGetNGS(petsclib::$UnionPetscLib, snes::AbstractSNES )
	f_ = Ref{Ptr{Cvoid}}()
	ctx_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:SNESGetNGS, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Ptr{Cvoid}}, Ptr{Cvoid}),
               snes, f_, ctx_,
              )

	f = f_[]
	ctx = ctx_[]

	return f,ctx
end 

"""
	pc::SNES = SNESGetNPC(petsclib::PetscLibType, snes::AbstractSNES) 
Gets a nonlinear preconditioning solver SNES` to be used to precondition the original nonlinear solver.

Not Collective; but any changes to the obtained the `pc` object must be applied collectively

Input Parameter:
- `snes` - iterative context obtained from `SNESCreate()`

Output Parameter:
- `pc` - the `SNES` preconditioner context

Options Database Key:
- `-npc_snes_type type` - set the type of the `SNES` to use as the nonlinear preconditioner

Level: advanced

See also: `SNESSetNPC()`, `SNESHasNPC()`, `SNES`, `SNESCreate()`

# External Links
$(_doc_external("SNES/SNESGetNPC"))
"""
function SNESGetNPC(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetNPC: no generated method for these argument types")
end

@for_petsc function SNESGetNPC(petsclib::$UnionPetscLib, snes::AbstractSNES )
	pc_ = Ref{CSNES}()

    @chk ccall(
               (:SNESGetNPC, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CSNES}),
               snes, pc_,
              )

	pc = SNES(pc_[], petsclib)

	return pc
end 

"""
	fnorm::PetscReal = SNESGetNPCFunction(petsclib::PetscLibType, snes::AbstractSNES, F::AbstractPetscVec) 
Gets the current function value (for the callback function provided by `SNESSetFunction()`,
and its norm from a nonlinear preconditioner after `SNESSolve()` has been called on that `SNES`

Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameters:
- `F`     - function vector
- `fnorm` - the norm of `F`

Level: developer

See also: `SNES`, `SNESGetNPC()`, `SNESSetNPC()`, `SNESComputeFunction()`, `SNESApplyNPC()`, `SNESSolve()`

# External Links
$(_doc_external("SNES/SNESGetNPCFunction"))
"""
function SNESGetNPCFunction(petsclib::PetscLibType, snes::AbstractSNES, F::AbstractPetscVec)
    error("SNESGetNPCFunction: no generated method for these argument types")
end

@for_petsc function SNESGetNPCFunction(petsclib::$UnionPetscLib, snes::AbstractSNES, F::AbstractPetscVec )
	fnorm_ = Ref{$PetscReal}()

    @chk ccall(
               (:SNESGetNPCFunction, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, Ptr{$PetscReal}),
               snes, F, fnorm_,
              )

	fnorm = fnorm_[]

	return fnorm
end 

"""
	side::PCSide = SNESGetNPCSide(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the preconditioning side used by the nonlinear preconditioner inside `SNES`.

Not Collective

Input Parameter:
- `snes` - iterative context obtained from `SNESCreate()`

Output Parameter:
- `side` - the preconditioning side, where side is one of
``
`PC_LEFT` - left preconditioning
`PC_RIGHT` - right preconditioning (default for most nonlinear solvers)
``

Level: intermediate

See also: `SNES`, `SNESGetNPC()`, `SNESSetNPCSide()`, `KSPGetPCSide()`, `PC_LEFT`, `PC_RIGHT`, `PCSide`

# External Links
$(_doc_external("SNES/SNESGetNPCSide"))
"""
function SNESGetNPCSide(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetNPCSide: no generated method for these argument types")
end

@for_petsc function SNESGetNPCSide(petsclib::$UnionPetscLib, snes::AbstractSNES )
	side_ = Ref{PCSide}()

    @chk ccall(
               (:SNESGetNPCSide, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{PCSide}),
               snes, side_,
              )

	side = side_[]

	return side
end 

"""
	nfails::PetscInt = SNESGetNonlinearStepFailures(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the number of unsuccessful steps
taken by the nonlinear solver in the current or most recent `SNESSolve()` .

Not Collective

Input Parameter:
- `snes` - `SNES` context

Output Parameter:
- `nfails` - number of unsuccessful steps attempted

Level: intermediate

See also: `SNES`, `SNESGetMaxLinearSolveFailures()`, `SNESGetLinearSolveIterations()`, `SNESSetMaxLinearSolveFailures()`, `SNESGetLinearSolveFailures()`,
`SNESSetMaxNonlinearStepFailures()`, `SNESGetMaxNonlinearStepFailures()`

# External Links
$(_doc_external("SNES/SNESGetNonlinearStepFailures"))
"""
function SNESGetNonlinearStepFailures(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetNonlinearStepFailures: no generated method for these argument types")
end

@for_petsc function SNESGetNonlinearStepFailures(petsclib::$UnionPetscLib, snes::AbstractSNES )
	nfails_ = Ref{$PetscInt}()

    @chk ccall(
               (:SNESGetNonlinearStepFailures, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscInt}),
               snes, nfails_,
              )

	nfails = nfails_[]

	return nfails
end 

"""
	normschedule::SNESNormSchedule = SNESGetNormSchedule(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the `SNESNormSchedule` used in convergence and monitoring
of the `SNES` method.

Logically Collective

Input Parameters:
- `snes`         - the `SNES` context
- `normschedule` - the type of the norm used

Level: advanced

See also: `SNES`, `SNESSetNormSchedule()`, `SNESComputeFunction()`, `VecNorm()`, `SNESSetFunction()`, `SNESSetInitialFunction()`, `SNESNormSchedule`

# External Links
$(_doc_external("SNES/SNESGetNormSchedule"))
"""
function SNESGetNormSchedule(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetNormSchedule: no generated method for these argument types")
end

@for_petsc function SNESGetNormSchedule(petsclib::$UnionPetscLib, snes::AbstractSNES )
	normschedule_ = Ref{SNESNormSchedule}()

    @chk ccall(
               (:SNESGetNormSchedule, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{SNESNormSchedule}),
               snes, normschedule_,
              )

	normschedule = normschedule_[]

	return normschedule
end 

"""
	nfuncs::PetscInt = SNESGetNumberFunctionEvals(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the number of user provided function evaluations
done by the `SNES` object in the current or most recent `SNESSolve()`

Not Collective

Input Parameter:
- `snes` - `SNES` context

Output Parameter:
- `nfuncs` - number of evaluations

Level: intermediate

See also: `SNES`, `SNESGetMaxLinearSolveFailures()`, `SNESGetLinearSolveIterations()`, `SNESSetMaxLinearSolveFailures()`, `SNESGetLinearSolveFailures()`, `SNESSetCountersReset()`

# External Links
$(_doc_external("SNES/SNESGetNumberFunctionEvals"))
"""
function SNESGetNumberFunctionEvals(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetNumberFunctionEvals: no generated method for these argument types")
end

@for_petsc function SNESGetNumberFunctionEvals(petsclib::$UnionPetscLib, snes::AbstractSNES )
	nfuncs_ = Ref{$PetscInt}()

    @chk ccall(
               (:SNESGetNumberFunctionEvals, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscInt}),
               snes, nfuncs_,
              )

	nfuncs = nfuncs_[]

	return nfuncs
end 

"""
	obj::Ptr{Cvoid},ctx::Ptr{Cvoid} = SNESGetObjective(petsclib::PetscLibType, snes::AbstractSNES) 
Returns the objective function set with `SNESSetObjective()`

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameters:
- `obj` - objective evaluation routine (or `NULL`); see `SNESObjectiveFn` for the calling sequence
- `ctx` - the function context (or `NULL`)

Level: advanced

See also: `SNES`, `SNESSetObjective()`, `SNESGetSolution()`, `SNESObjectiveFn`

# External Links
$(_doc_external("SNES/SNESGetObjective"))
"""
function SNESGetObjective(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetObjective: no generated method for these argument types")
end

@for_petsc function SNESGetObjective(petsclib::$UnionPetscLib, snes::AbstractSNES )
	obj_ = Ref{Ptr{Cvoid}}()
	ctx_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:SNESGetObjective, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Ptr{Cvoid}}, Ptr{Cvoid}),
               snes, obj_, ctx_,
              )

	obj = obj_[]
	ctx = ctx_[]

	return obj,ctx
end 

"""
	prefix::Ptr{Cchar} = SNESGetOptionsPrefix(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the prefix used for searching for all
`SNES` options in the database.

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `prefix` - pointer to the prefix string used

Level: advanced

See also: `SNES`, `SNESSetOptionsPrefix()`, `SNESAppendOptionsPrefix()`

# External Links
$(_doc_external("SNES/SNESGetOptionsPrefix"))
"""
function SNESGetOptionsPrefix(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetOptionsPrefix: no generated method for these argument types")
end

@for_petsc function SNESGetOptionsPrefix(petsclib::$UnionPetscLib, snes::AbstractSNES )
	prefix_ = Ref{Ptr{Cchar}}()

    @chk ccall(
               (:SNESGetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Ptr{Cchar}}),
               snes, prefix_,
              )

	prefix = prefix_[]

	return prefix
end 

"""
	r::PetscVec,f::Ptr{Cvoid},Amat::PetscMat,Pmat::PetscMat,J::Ptr{Cvoid},ctx::Ptr{Cvoid} = SNESGetPicard(petsclib::PetscLibType, snes::AbstractSNES) 
Returns the context for the Picard iteration

Not Collective, but `Vec` is parallel if `SNES` is parallel. Collective if `Vec` is requested, but has not been created yet.

Input Parameter:
- `snes` - the `SNES` context

Output Parameters:
- `r`    - the function (or `NULL`)
- `f`    - the function (or `NULL`);  for calling sequence see `SNESFunctionFn`
- `Amat` - the matrix used to defined the operation A(x) x - b(x) (or `NULL`)
- `Pmat` - the matrix from which the preconditioner will be constructed (or `NULL`)
- `J`    - the function for matrix evaluation (or `NULL`);  for calling sequence see `SNESJacobianFn`
- `ctx`  - the function context (or `NULL`)

Level: advanced

See also: `SNESSetFunction()`, `SNESSetPicard()`, `SNESGetFunction()`, `SNESGetJacobian()`, `SNESGetDM()`, `SNESFunctionFn`, `SNESJacobianFn`

# External Links
$(_doc_external("SNES/SNESGetPicard"))
"""
function SNESGetPicard(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetPicard: no generated method for these argument types")
end

@for_petsc function SNESGetPicard(petsclib::$UnionPetscLib, snes::AbstractSNES )
	r_ = Ref{CVec}()
	f_ = Ref{Ptr{Cvoid}}()
	Amat_ = Ref{CMat}()
	Pmat_ = Ref{CMat}()
	J_ = Ref{Ptr{Cvoid}}()
	ctx_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:SNESGetPicard, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CVec}, Ptr{Ptr{Cvoid}}, Ptr{CMat}, Ptr{CMat}, Ptr{Ptr{Cvoid}}, Ptr{Cvoid}),
               snes, r_, f_, Amat_, Pmat_, J_, ctx_,
              )

	r = PetscVec(r_[], petsclib)
	f = f_[]
	Amat = PetscMat(Amat_[], petsclib)
	Pmat = PetscMat(Pmat_[], petsclib)
	J = J_[]
	ctx = ctx_[]

	return r,f,Amat,Pmat,J,ctx
end 

"""
	rhs::PetscVec = SNESGetRhs(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the vector for solving F(x) = `rhs`. If `rhs` is not set
it assumes a zero right-hand side.

Logically Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `rhs` - the right-hand side vector or `NULL` if there is no right-hand side vector

Level: intermediate

See also: `SNES`, `SNESGetSolution()`, `SNESGetFunction()`, `SNESComputeFunction()`, `SNESSetJacobian()`, `SNESSetFunction()`

# External Links
$(_doc_external("SNES/SNESGetRhs"))
"""
function SNESGetRhs(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetRhs: no generated method for these argument types")
end

@for_petsc function SNESGetRhs(petsclib::$UnionPetscLib, snes::AbstractSNES )
	rhs_ = Ref{CVec}()

    @chk ccall(
               (:SNESGetRhs, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CVec}),
               snes, rhs_,
              )

	rhs = PetscVec(rhs_[], petsclib)

	return rhs
end 

"""
	x::PetscVec = SNESGetSolution(petsclib::PetscLibType, snes::AbstractSNES) 
Returns the vector where the approximate solution is
stored. This is the fine grid solution when using `SNESSetGridSequence()`.

Not Collective, but `x` is parallel if `snes` is parallel

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `x` - the solution

Level: intermediate

See also: `SNESSetSolution()`, `SNESSolve()`, `SNES`, `SNESGetSolutionUpdate()`, `SNESGetFunction()`

# External Links
$(_doc_external("SNES/SNESGetSolution"))
"""
function SNESGetSolution(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetSolution: no generated method for these argument types")
end

@for_petsc function SNESGetSolution(petsclib::$UnionPetscLib, snes::AbstractSNES )
	x_ = Ref{CVec}()

    @chk ccall(
               (:SNESGetSolution, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CVec}),
               snes, x_,
              )

	x = PetscVec(x_[], petsclib)

	return x
end 

"""
	xnorm::PetscReal = SNESGetSolutionNorm(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the last computed norm of the solution

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `xnorm` - the last computed solution norm

Level: developer

See also: `SNES`, `SNESSetNormSchedule()`, `SNESComputeFunction()`, `SNESGetFunctionNorm()`, `SNESGetUpdateNorm()`

# External Links
$(_doc_external("SNES/SNESGetSolutionNorm"))
"""
function SNESGetSolutionNorm(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetSolutionNorm: no generated method for these argument types")
end

@for_petsc function SNESGetSolutionNorm(petsclib::$UnionPetscLib, snes::AbstractSNES )
	xnorm_ = Ref{$PetscReal}()

    @chk ccall(
               (:SNESGetSolutionNorm, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscReal}),
               snes, xnorm_,
              )

	xnorm = xnorm_[]

	return xnorm
end 

"""
	x::PetscVec = SNESGetSolutionUpdate(petsclib::PetscLibType, snes::AbstractSNES) 
Returns the vector where the solution update is
stored.

Not Collective, but `x` is parallel if `snes` is parallel

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `x` - the solution update

Level: advanced

See also: `SNES`, `SNESGetSolution()`, `SNESGetFunction()`

# External Links
$(_doc_external("SNES/SNESGetSolutionUpdate"))
"""
function SNESGetSolutionUpdate(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetSolutionUpdate: no generated method for these argument types")
end

@for_petsc function SNESGetSolutionUpdate(petsclib::$UnionPetscLib, snes::AbstractSNES )
	x_ = Ref{CVec}()

    @chk ccall(
               (:SNESGetSolutionUpdate, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CVec}),
               snes, x_,
              )

	x = PetscVec(x_[], petsclib)

	return x
end 

"""
	atol::PetscReal,rtol::PetscReal,stol::PetscReal,maxit::PetscInt,maxf::PetscInt = SNESGetTolerances(petsclib::PetscLibType, snes::AbstractSNES) 
Gets various parameters used in `SNES` convergence tests.

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameters:
- `atol`  - the absolute convergence tolerance
- `rtol`  - the relative convergence tolerance
- `stol`  - convergence tolerance in terms of the norm of the change in the solution between steps
- `maxit` - the maximum number of iterations allowed
- `maxf`  - the maximum number of function evaluations allowed, `PETSC_UNLIMITED` indicates no bound

Level: intermediate

See also: `SNES`, `SNESSetTolerances()`

# External Links
$(_doc_external("SNES/SNESGetTolerances"))
"""
function SNESGetTolerances(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetTolerances: no generated method for these argument types")
end

@for_petsc function SNESGetTolerances(petsclib::$UnionPetscLib, snes::AbstractSNES )
	atol_ = Ref{$PetscReal}()
	rtol_ = Ref{$PetscReal}()
	stol_ = Ref{$PetscReal}()
	maxit_ = Ref{$PetscInt}()
	maxf_ = Ref{$PetscInt}()

    @chk ccall(
               (:SNESGetTolerances, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscInt}, Ptr{$PetscInt}),
               snes, atol_, rtol_, stol_, maxit_, maxf_,
              )

	atol = atol_[]
	rtol = rtol_[]
	stol = stol_[]
	maxit = maxit_[]
	maxf = maxf_[]

	return atol,rtol,stol,maxit,maxf
end 

"""
	type::SNESType = SNESGetType(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the `SNES` method type and name (as a string).

Not Collective

Input Parameter:
- `snes` - nonlinear solver context

Output Parameter:
- `type` - `SNES` method (a character string)

Level: intermediate

See also: `SNESSetType()`, `SNESType`, `SNESSetFromOptions()`, `SNES`, `PetscObjectTypeCompare()`, `PetscObjectTypeCompareAny()`

# External Links
$(_doc_external("SNES/SNESGetType"))
"""
function SNESGetType(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetType: no generated method for these argument types")
end

@for_petsc function SNESGetType(petsclib::$UnionPetscLib, snes::AbstractSNES )
	type_ = Ref{SNESType}()

    @chk ccall(
               (:SNESGetType, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{SNESType}),
               snes, type_,
              )

	type = type_[] == C_NULL ? "" : unsafe_string(type_[])

	return type
end 

"""
	ynorm::PetscReal = SNESGetUpdateNorm(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the last computed norm of the solution update

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `ynorm` - the last computed update norm

Level: developer

See also: `SNES`, `SNESSetNormSchedule()`, `SNESComputeFunction()`, `SNESGetFunctionNorm()`

# External Links
$(_doc_external("SNES/SNESGetUpdateNorm"))
"""
function SNESGetUpdateNorm(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetUpdateNorm: no generated method for these argument types")
end

@for_petsc function SNESGetUpdateNorm(petsclib::$UnionPetscLib, snes::AbstractSNES )
	ynorm_ = Ref{$PetscReal}()

    @chk ccall(
               (:SNESGetUpdateNorm, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscReal}),
               snes, ynorm_,
              )

	ynorm = ynorm_[]

	return ynorm
end 

"""
	mf_operator::PetscBool,mf::PetscBool = SNESGetUseMatrixFree(petsclib::PetscLibType, snes::AbstractSNES) 
indicates if the `SNES` uses matrix-free finite difference matrix vector products to apply the Jacobian.

Not Collective, but the resulting flags will be the same on all MPI processes

Input Parameter:
- `snes` - `SNES` context

Output Parameters:
- `mf_operator` - use matrix-free only for the Amat used by `SNESSetJacobian()`, this means the user provided Pmat will continue to be used
- `mf`          - use matrix-free for both the Amat and Pmat used by `SNESSetJacobian()`, both the Amat and Pmat set in `SNESSetJacobian()` will be ignored

Level: intermediate

See also: `SNES`, `SNESSetUseMatrixFree()`, `MatCreateSNESMF()`

# External Links
$(_doc_external("SNES/SNESGetUseMatrixFree"))
"""
function SNESGetUseMatrixFree(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESGetUseMatrixFree: no generated method for these argument types")
end

@for_petsc function SNESGetUseMatrixFree(petsclib::$UnionPetscLib, snes::AbstractSNES )
	mf_operator_ = Ref{PetscBool}()
	mf_ = Ref{PetscBool}()

    @chk ccall(
               (:SNESGetUseMatrixFree, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{PetscBool}, Ptr{PetscBool}),
               snes, mf_operator_, mf_,
              )

	mf_operator = mf_operator_[]
	mf = mf_[]

	return mf_operator,mf
end 

"""
	has_npc::PetscBool = SNESHasNPC(petsclib::PetscLibType, snes::AbstractSNES) 
Returns whether a nonlinear preconditioner is associated with the given `SNES`

Not Collective

Input Parameter:
- `snes` - iterative context obtained from `SNESCreate()`

Output Parameter:
- `has_npc` - whether the `SNES` has a nonlinear preconditioner or not

Level: developer

See also: `SNESSetNPC()`, `SNESGetNPC()`

# External Links
$(_doc_external("SNES/SNESHasNPC"))
"""
function SNESHasNPC(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESHasNPC: no generated method for these argument types")
end

@for_petsc function SNESHasNPC(petsclib::$UnionPetscLib, snes::AbstractSNES )
	has_npc_ = Ref{PetscBool}()

    @chk ccall(
               (:SNESHasNPC, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{PetscBool}),
               snes, has_npc_,
              )

	has_npc = has_npc_[]

	return has_npc
end 

"""
	SNESInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `SNES` package. It is called
from PetscDLLibraryRegister_petscsnes() when using dynamic libraries, and on the first call to `SNESCreate()`
when using shared or static libraries.

Level: developer

See also: `SNES`, `PetscInitialize()`

# External Links
$(_doc_external("SNES/SNESInitializePackage"))
"""
function SNESInitializePackage(petsclib::PetscLibType)
    error("SNESInitializePackage: no generated method for these argument types")
end

@for_petsc function SNESInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:SNESInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	version::PetscInt,rtol_0::PetscReal,rtol_max::PetscReal,gamma::PetscReal,alpha::PetscReal,alpha2::PetscReal,threshold::PetscReal = SNESKSPGetParametersEW(petsclib::PetscLibType, snes::AbstractSNES) 
Gets parameters for Eisenstat-Walker
convergence criteria for the linear solvers within an inexact
Newton method.

Not Collective

Input Parameter:
- `snes` - `SNES` context

Output Parameters:
- `version`   - version 1, 2 (default is 2), 3 or 4
- `rtol_0`    - initial relative tolerance (0 <= rtol_0 < 1)
- `rtol_max`  - maximum relative tolerance (0 <= rtol_max < 1)
- `gamma`     - multiplicative factor for version 2 rtol computation (0 <= gamma2 <= 1)
- `alpha`     - power for version 2 rtol computation (1 < alpha <= 2)
- `alpha2`    - power for safeguard
- `threshold` - threshold for imposing safeguard (0 < threshold < 1)

Level: advanced

See also: `SNES`, `SNESKSPSetUseEW()`, `SNESKSPGetUseEW()`, `SNESKSPSetParametersEW()`

# External Links
$(_doc_external("SNES/SNESKSPGetParametersEW"))
"""
function SNESKSPGetParametersEW(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESKSPGetParametersEW: no generated method for these argument types")
end

@for_petsc function SNESKSPGetParametersEW(petsclib::$UnionPetscLib, snes::AbstractSNES )
	version_ = Ref{$PetscInt}()
	rtol_0_ = Ref{$PetscReal}()
	rtol_max_ = Ref{$PetscReal}()
	gamma_ = Ref{$PetscReal}()
	alpha_ = Ref{$PetscReal}()
	alpha2_ = Ref{$PetscReal}()
	threshold_ = Ref{$PetscReal}()

    @chk ccall(
               (:SNESKSPGetParametersEW, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscInt}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               snes, version_, rtol_0_, rtol_max_, gamma_, alpha_, alpha2_, threshold_,
              )

	version = version_[]
	rtol_0 = rtol_0_[]
	rtol_max = rtol_max_[]
	gamma = gamma_[]
	alpha = alpha_[]
	alpha2 = alpha2_[]
	threshold = threshold_[]

	return version,rtol_0,rtol_max,gamma,alpha,alpha2,threshold
end 

"""
	flag::PetscBool = SNESKSPGetUseEW(petsclib::PetscLibType, snes::AbstractSNES) 
Gets if `SNES` is using Eisenstat-Walker method
for computing relative tolerance for linear solvers within an
inexact Newton method.

Not Collective

Input Parameter:
- `snes` - `SNES` context

Output Parameter:
- `flag` - `PETSC_TRUE` or `PETSC_FALSE`

Level: advanced

See also: `SNESKSPSetUseEW()`, `SNESKSPGetParametersEW()`, `SNESKSPSetParametersEW()`

# External Links
$(_doc_external("SNES/SNESKSPGetUseEW"))
"""
function SNESKSPGetUseEW(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESKSPGetUseEW: no generated method for these argument types")
end

@for_petsc function SNESKSPGetUseEW(petsclib::$UnionPetscLib, snes::AbstractSNES )
	flag_ = Ref{PetscBool}()

    @chk ccall(
               (:SNESKSPGetUseEW, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{PetscBool}),
               snes, flag_,
              )

	flag = flag_[]

	return flag
end 

"""
	SNESKSPSetParametersEW(petsclib::PetscLibType, snes::AbstractSNES, version::PetscInt, rtol_0::PetscReal, rtol_max::PetscReal, gamma::PetscReal, alpha::PetscReal, alpha2::PetscReal, threshold::PetscReal) 
Sets parameters for Eisenstat-Walker
convergence criteria for the linear solvers within an inexact
Newton method.

Logically Collective

Input Parameters:
- `snes`      - `SNES` context
- `version`   - version 1, 2 (default is 2), 3 or 4
- `rtol_0`    - initial relative tolerance (0 <= rtol_0 < 1)
- `rtol_max`  - maximum relative tolerance (0 <= rtol_max < 1)
- `gamma`     - multiplicative factor for version 2 rtol computation
(0 <= gamma2 <= 1)
- `alpha`     - power for version 2 rtol computation (1 < alpha <= 2)
- `alpha2`    - power for safeguard
- `threshold` - threshold for imposing safeguard (0 < threshold < 1)

Level: advanced

See also: `SNES`, `SNESKSPSetUseEW()`, `SNESKSPGetUseEW()`, `SNESKSPGetParametersEW()`

# External Links
$(_doc_external("SNES/SNESKSPSetParametersEW"))
"""
function SNESKSPSetParametersEW(petsclib::PetscLibType, snes::AbstractSNES, version::Integer, rtol_0::Real, rtol_max::Real, gamma::Real, alpha::Real, alpha2::Real, threshold::Real)
    error("SNESKSPSetParametersEW: no generated method for these argument types")
end

@for_petsc function SNESKSPSetParametersEW(petsclib::$UnionPetscLib, snes::AbstractSNES, version::$PetscInt, rtol_0::$PetscReal, rtol_max::$PetscReal, gamma::$PetscReal, alpha::$PetscReal, alpha2::$PetscReal, threshold::$PetscReal )

    @chk ccall(
               (:SNESKSPSetParametersEW, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, $PetscReal, $PetscReal, $PetscReal, $PetscReal, $PetscReal),
               snes, version, rtol_0, rtol_max, gamma, alpha, alpha2, threshold,
              )


	return nothing
end 

"""
	SNESKSPSetUseEW(petsclib::PetscLibType, snes::AbstractSNES, flag::PetscBool) 
Sets `SNES` to the use Eisenstat-Walker method for
computing relative tolerance for linear solvers within an inexact
Newton method.

Logically Collective

Input Parameters:
- `snes` - `SNES` context
- `flag` - `PETSC_TRUE` or `PETSC_FALSE`

Options Database Keys:
- `-snes_ksp_ew`                     - use Eisenstat-Walker method for determining linear system convergence
- `-snes_ksp_ew_version ver`         - version of  Eisenstat-Walker method
- `-snes_ksp_ew_rtol0 rtol0`         - Sets rtol0
- `-snes_ksp_ew_rtolmax rtolmax`     - Sets rtolmax
- `-snes_ksp_ew_gamma gamma`         - Sets gamma
- `-snes_ksp_ew_alpha alpha`         - Sets alpha
- `-snes_ksp_ew_alpha2 alpha2`       - Sets alpha2
- `-snes_ksp_ew_threshold threshold` - Sets threshold

Level: advanced

See also: `KSP`, `SNES`, `SNESKSPGetUseEW()`, `SNESKSPGetParametersEW()`, `SNESKSPSetParametersEW()`

# External Links
$(_doc_external("SNES/SNESKSPSetUseEW"))
"""
function SNESKSPSetUseEW(petsclib::PetscLibType, snes::AbstractSNES, flag::PetscBool)
    error("SNESKSPSetUseEW: no generated method for these argument types")
end

@for_petsc function SNESKSPSetUseEW(petsclib::$UnionPetscLib, snes::AbstractSNES, flag::PetscBool )

    @chk ccall(
               (:SNESKSPSetUseEW, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, flag,
              )


	return nothing
end 

"""
	SNESLoad(petsclib::PetscLibType, snes::AbstractSNES, viewer::PetscViewer) 
Loads a `SNES` that has been stored in `PETSCVIEWERBINARY` with `SNESView()`.

Collective

Input Parameters:
- `snes`   - the newly loaded `SNES`, this needs to have been created with `SNESCreate()` or
some related function before a call to `SNESLoad()`.
- `viewer` - binary file viewer, obtained from `PetscViewerBinaryOpen()`

Level: intermediate

See also: `SNES`, `PetscViewer`, `SNESCreate()`, `SNESType`, `PetscViewerBinaryOpen()`, `SNESView()`, `MatLoad()`, `VecLoad()`

# External Links
$(_doc_external("SNES/SNESLoad"))
"""
function SNESLoad(petsclib::PetscLibType, snes::AbstractSNES, viewer::PetscViewer)
    error("SNESLoad: no generated method for these argument types")
end

@for_petsc function SNESLoad(petsclib::$UnionPetscLib, snes::AbstractSNES, viewer::PetscViewer )

    @chk ccall(
               (:SNESLoad, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscViewer),
               snes, viewer,
              )


	return nothing
end 

"""
	SNESMSFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the `SNESMS` package. It is
called from `PetscFinalize()`.

Level: developer

See also: `SNES`, `SNESMS`, `SNESMSRegister()`, `SNESMSRegisterAll()`, `SNESMSInitializePackage()`, `PetscFinalize()`

# External Links
$(_doc_external("SNES/SNESMSFinalizePackage"))
"""
function SNESMSFinalizePackage(petsclib::PetscLibType)
    error("SNESMSFinalizePackage: no generated method for these argument types")
end

@for_petsc function SNESMSFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:SNESMSFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	damping::PetscReal = SNESMSGetDamping(petsclib::PetscLibType, snes::AbstractSNES) 
Get the damping parameter of `SNESMS` multistage scheme

Not Collective

Input Parameter:
- `snes` - nonlinear solver context

Output Parameter:
- `damping` - damping parameter

Level: advanced

See also: `SNESMSSetDamping()`, `SNESMS`

# External Links
$(_doc_external("SNES/SNESMSGetDamping"))
"""
function SNESMSGetDamping(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESMSGetDamping: no generated method for these argument types")
end

@for_petsc function SNESMSGetDamping(petsclib::$UnionPetscLib, snes::AbstractSNES )
	damping_ = Ref{$PetscReal}()

    @chk ccall(
               (:SNESMSGetDamping, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscReal}),
               snes, damping_,
              )

	damping = damping_[]

	return damping
end 

"""
	mstype::SNESMSType = SNESMSGetType(petsclib::PetscLibType, snes::AbstractSNES) 
Get the type of multistage smoother `SNESMS`

Not Collective

Input Parameter:
- `snes` - nonlinear solver context

Output Parameter:
- `mstype` - type of multistage method

Level: advanced

See also: `SNESMS`, `SNESMSSetType()`, `SNESMSType`

# External Links
$(_doc_external("SNES/SNESMSGetType"))
"""
function SNESMSGetType(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESMSGetType: no generated method for these argument types")
end

@for_petsc function SNESMSGetType(petsclib::$UnionPetscLib, snes::AbstractSNES )
	mstype_ = Ref{SNESMSType}()

    @chk ccall(
               (:SNESMSGetType, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{SNESMSType}),
               snes, mstype_,
              )

	mstype = mstype_[] == C_NULL ? "" : unsafe_string(mstype_[])

	return mstype
end 

"""
	SNESMSInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `SNESMS` package. It is called
from `SNESInitializePackage()`.

Level: developer

See also: `SNES`, `SNESMS`, `SNESMSRegister()`, `SNESMSRegisterAll()`, `PetscInitialize()`

# External Links
$(_doc_external("SNES/SNESMSInitializePackage"))
"""
function SNESMSInitializePackage(petsclib::PetscLibType)
    error("SNESMSInitializePackage: no generated method for these argument types")
end

@for_petsc function SNESMSInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:SNESMSInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	SNESMSRegister(petsclib::PetscLibType, name::SNESMSType, nstages::PetscInt, nregisters::PetscInt, stability::PetscReal, gamma::Vector{PetscReal}, delta::Vector{PetscReal}, betasub::Vector{PetscReal}) 
register a multistage scheme for `SNESMS`

Logically Collective, No Fortran Support

Input Parameters:
- `name`       - identifier for method
- `nstages`    - number of stages
- `nregisters` - number of registers used by low-storage implementation
- `stability`  - scaled stability region
- `gamma`      - coefficients, see Ketcheson's paper {cite}`ketcheson2010runge`
- `delta`      - coefficients, see Ketcheson's paper {cite}`ketcheson2010runge`
- `betasub`    - subdiagonal of Shu-Osher form

Level: advanced

See also: `SNES`, `SNESMS`

# External Links
$(_doc_external("SNES/SNESMSRegister"))
"""
function SNESMSRegister(petsclib::PetscLibType, name::SNESMSType, nstages::Integer, nregisters::Integer, stability::Real, gamma::AbstractVector{<:Number}, delta::AbstractVector{<:Number}, betasub::AbstractVector{<:Number})
    error("SNESMSRegister: no generated method for these argument types")
end

@for_petsc function SNESMSRegister(petsclib::$UnionPetscLib, name::SNESMSType, nstages::$PetscInt, nregisters::$PetscInt, stability::$PetscReal, gamma::Vector{$PetscReal}, delta::Vector{$PetscReal}, betasub::Vector{$PetscReal} )

    @chk ccall(
               (:SNESMSRegister, $petsc_library),
               PetscErrorCode,
               (SNESMSType, $PetscInt, $PetscInt, $PetscReal, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               name, nstages, nregisters, stability, gamma, delta, betasub,
              )


	return nothing
end 

"""
	SNESMSRegisterAll(petsclib::PetscLibType) 
Registers all of the multi-stage methods in `SNESMS`

Logically Collective

Level: developer

See also: `SNES`, `SNESMS`, `SNESMSRegisterDestroy()`

# External Links
$(_doc_external("SNES/SNESMSRegisterAll"))
"""
function SNESMSRegisterAll(petsclib::PetscLibType)
    error("SNESMSRegisterAll: no generated method for these argument types")
end

@for_petsc function SNESMSRegisterAll(petsclib::$UnionPetscLib)

    @chk ccall(
               (:SNESMSRegisterAll, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	SNESMSRegisterDestroy(petsclib::PetscLibType) 
Frees the list of schemes that were registered by `SNESMSRegister()`.

Logically Collective

Level: developer

See also: `SNES`, `SNESMS`, `SNESMSRegister()`, `SNESMSRegisterAll()`

# External Links
$(_doc_external("SNES/SNESMSRegisterDestroy"))
"""
function SNESMSRegisterDestroy(petsclib::PetscLibType)
    error("SNESMSRegisterDestroy: no generated method for these argument types")
end

@for_petsc function SNESMSRegisterDestroy(petsclib::$UnionPetscLib)

    @chk ccall(
               (:SNESMSRegisterDestroy, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	SNESMSSetDamping(petsclib::PetscLibType, snes::AbstractSNES, damping::PetscReal) 
Set the damping parameter for a `SNESMS` multistage scheme

Logically Collective

Input Parameters:
- `snes`    - nonlinear solver context
- `damping` - damping parameter

Level: advanced

See also: `SNESMSGetDamping()`, `SNESMS`

# External Links
$(_doc_external("SNES/SNESMSSetDamping"))
"""
function SNESMSSetDamping(petsclib::PetscLibType, snes::AbstractSNES, damping::Real)
    error("SNESMSSetDamping: no generated method for these argument types")
end

@for_petsc function SNESMSSetDamping(petsclib::$UnionPetscLib, snes::AbstractSNES, damping::$PetscReal )

    @chk ccall(
               (:SNESMSSetDamping, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscReal),
               snes, damping,
              )


	return nothing
end 

"""
	SNESMSSetType(petsclib::PetscLibType, snes::AbstractSNES, mstype::SNESMSType) 
Set the type of multistage smoother `SNESMS`

Logically Collective

Input Parameters:
- `snes`   - nonlinear solver context
- `mstype` - type of multistage method

Level: advanced

See also: `SNESMS`, `SNESMSGetType()`, `SNESMSType`

# External Links
$(_doc_external("SNES/SNESMSSetType"))
"""
function SNESMSSetType(petsclib::PetscLibType, snes::AbstractSNES, mstype::SNESMSType)
    error("SNESMSSetType: no generated method for these argument types")
end

@for_petsc function SNESMSSetType(petsclib::$UnionPetscLib, snes::AbstractSNES, mstype::SNESMSType )

    @chk ccall(
               (:SNESMSSetType, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESMSType),
               snes, mstype,
              )


	return nothing
end 

"""
	SNESMonitor(petsclib::PetscLibType, snes::AbstractSNES, iter::PetscInt, rnorm::PetscReal) 
runs any `SNES` monitor routines provided with `SNESMonitor()` or the options database

Collective

Input Parameters:
- `snes`  - nonlinear solver context obtained from `SNESCreate()`
- `iter`  - current iteration number
- `rnorm` - current relative norm of the residual

Level: developer

See also: `SNES`, `SNESMonitorSet()`

# External Links
$(_doc_external("SNES/SNESMonitor"))
"""
function SNESMonitor(petsclib::PetscLibType, snes::AbstractSNES, iter::Integer, rnorm::Real)
    error("SNESMonitor: no generated method for these argument types")
end

@for_petsc function SNESMonitor(petsclib::$UnionPetscLib, snes::AbstractSNES, iter::$PetscInt, rnorm::$PetscReal )

    @chk ccall(
               (:SNESMonitor, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal),
               snes, iter, rnorm,
              )


	return nothing
end 

"""
	SNESMonitorCancel(petsclib::PetscLibType, snes::AbstractSNES) 
Clears all the monitor functions for a `SNES` object.

Logically Collective

Input Parameter:
- `snes` - the `SNES` context

Options Database Key:
- `-snes_monitor_cancel` - cancels all monitors that have been hardwired
into a code by calls to `SNESMonitorSet()`, but does not cancel those
set via the options database

Level: intermediate

See also: `SNES`, `SNESMonitorDefault()`, `SNESMonitorSet()`

# External Links
$(_doc_external("SNES/SNESMonitorCancel"))
"""
function SNESMonitorCancel(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESMonitorCancel: no generated method for these argument types")
end

@for_petsc function SNESMonitorCancel(petsclib::$UnionPetscLib, snes::AbstractSNES )

    @chk ccall(
               (:SNESMonitorCancel, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	SNESMonitorDefault(petsclib::PetscLibType, snes::AbstractSNES, its::PetscInt, fgnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Monitors progress of a `SNESSolve()` (default).

Collective

Input Parameters:
- `snes`   - the `SNES` context
- `its`    - iteration number
- `fgnorm` - 2-norm of residual
- `vf`     - viewer and format structure

Options Database Key:
- `-snes_monitor` - use this function to monitor the convergence of the nonlinear solver

Level: intermediate

See also: `SNESMonitorSet()`, `SNESMonitorSolution()`, `SNESMonitorFunction()`, `SNESMonitorResidual()`,
`SNESMonitorSolutionUpdate()`, `SNESMonitorScaling()`, `SNESMonitorRange()`, `SNESMonitorRatio()`,
`SNESMonitorDefaultField()`, `PetscViewerFormat`, `PetscViewerAndFormat`

# External Links
$(_doc_external("SNES/SNESMonitorDefault"))
"""
function SNESMonitorDefault(petsclib::PetscLibType, snes::AbstractSNES, its::Integer, fgnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("SNESMonitorDefault: no generated method for these argument types")
end

@for_petsc function SNESMonitorDefault(petsclib::$UnionPetscLib, snes::AbstractSNES, its::$PetscInt, fgnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:SNESMonitorDefault, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               snes, its, fgnorm, vf,
              )


	return nothing
end 

"""
	SNESMonitorDefaultField(petsclib::PetscLibType, snes::AbstractSNES, its::PetscInt, fgnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Monitors progress of a `SNESSolve()`, separated into fields.

Collective

Input Parameters:
- `snes`   - the `SNES` context
- `its`    - iteration number
- `fgnorm` - 2-norm of residual
- `vf`     - the PetscViewer

Options Database Key:
- `-snes_monitor_field` - activate this monitor

Level: intermediate

See also: `SNESMonitorSet()`, `SNESMonitorSolution()`, `SNESMonitorDefault()`, `PetscViewerFormat`, `PetscViewerAndFormat`

# External Links
$(_doc_external("SNES/SNESMonitorDefaultField"))
"""
function SNESMonitorDefaultField(petsclib::PetscLibType, snes::AbstractSNES, its::Integer, fgnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("SNESMonitorDefaultField: no generated method for these argument types")
end

@for_petsc function SNESMonitorDefaultField(petsclib::$UnionPetscLib, snes::AbstractSNES, its::$PetscInt, fgnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:SNESMonitorDefaultField, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               snes, its, fgnorm, vf,
              )


	return nothing
end 

"""
	SNESMonitorDefaultSetUp(petsclib::PetscLibType, snes::AbstractSNES, vf::Vector{PetscViewerAndFormat}) 
Prepare the `PetscViewerAndFormat` associated with `SNESMonitorDefault()`, in particular by initializing the underlying `PetscDrawLG` when the viewer format is `PETSC_VIEWER_DRAW_LG`

Collective

Input Parameters:
- `snes` - the `SNES` context
- `vf`   - the viewer/format pair passed to `SNESMonitorSet()` along with `SNESMonitorDefault()`

Level: developer

See also: `SNES`, `SNESMonitorSet()`, `SNESMonitorDefault()`, `PetscViewerAndFormat`, `PetscViewerMonitorLGSetUp()`

# External Links
$(_doc_external("SNES/SNESMonitorDefaultSetUp"))
"""
function SNESMonitorDefaultSetUp(petsclib::PetscLibType, snes::AbstractSNES, vf::Vector{PetscViewerAndFormat})
    error("SNESMonitorDefaultSetUp: no generated method for these argument types")
end

@for_petsc function SNESMonitorDefaultSetUp(petsclib::$UnionPetscLib, snes::AbstractSNES, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:SNESMonitorDefaultSetUp, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{PetscViewerAndFormat}),
               snes, vf,
              )


	return nothing
end 

"""
	SNESMonitorDefaultShort(petsclib::PetscLibType, snes::AbstractSNES, its::PetscInt, fgnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 

# External Links
$(_doc_external("SNES/SNESMonitorDefaultShort"))
"""
function SNESMonitorDefaultShort(petsclib::PetscLibType, snes::AbstractSNES, its::Integer, fgnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("SNESMonitorDefaultShort: no generated method for these argument types")
end

@for_petsc function SNESMonitorDefaultShort(petsclib::$UnionPetscLib, snes::AbstractSNES, its::$PetscInt, fgnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:SNESMonitorDefaultShort, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               snes, its, fgnorm, vf,
              )


	return nothing
end 

"""
	SNESMonitorFields(petsclib::PetscLibType, snes::AbstractSNES, its::PetscInt, fgnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Monitors the residual for each field separately

Collective

Input Parameters:
- `snes`   - the `SNES` context, must have an attached `DM`
- `its`    - iteration number
- `fgnorm` - 2-norm of residual
- `vf`     - `PetscViewerAndFormat` of `PetscViewerType` `PETSCVIEWERASCII`

Level: intermediate

See also: `SNES`, `SNESMonitorSet()`, `SNESMonitorDefault()`

# External Links
$(_doc_external("SNES/SNESMonitorFields"))
"""
function SNESMonitorFields(petsclib::PetscLibType, snes::AbstractSNES, its::Integer, fgnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("SNESMonitorFields: no generated method for these argument types")
end

@for_petsc function SNESMonitorFields(petsclib::$UnionPetscLib, snes::AbstractSNES, its::$PetscInt, fgnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:SNESMonitorFields, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               snes, its, fgnorm, vf,
              )


	return nothing
end 

"""
	SNESMonitorJacUpdateSpectrum(petsclib::PetscLibType, snes::AbstractSNES, it::PetscInt, fnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Monitors the spectrun of the change in the Jacobian from the last Jacobian evaluation of a `SNESSolve()`

Collective

Input Parameters:
- `snes`  - the `SNES` context
- `it`    - iteration number
- `fnorm` - 2-norm of residual
- `vf`    - viewer and format structure

Options Database Key:
- `-snes_monitor_jacupdate_spectrum` - activates this monitor

Level: intermediate

See also: `SNESMonitorSet()`, `SNESMonitorSolution()`, `SNESMonitorRange()`, `PetscViewerFormat`, `PetscViewerAndFormat`

# External Links
$(_doc_external("SNES/SNESMonitorJacUpdateSpectrum"))
"""
function SNESMonitorJacUpdateSpectrum(petsclib::PetscLibType, snes::AbstractSNES, it::Integer, fnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("SNESMonitorJacUpdateSpectrum: no generated method for these argument types")
end

@for_petsc function SNESMonitorJacUpdateSpectrum(petsclib::$UnionPetscLib, snes::AbstractSNES, it::$PetscInt, fnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:SNESMonitorJacUpdateSpectrum, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               snes, it, fnorm, vf,
              )


	return nothing
end 

"""
	SNESMonitorLGRange(petsclib::PetscLibType, snes::AbstractSNES, n::PetscInt, rnorm::PetscReal, monctx::Ptr{Cvoid}) 
Line-graph monitor that plots the residual norm together with residual-range statistics for a `SNESSolve()`

Collective

Input Parameters:
- `snes`   - the `SNES` context
- `n`      - the iteration number
- `rnorm`  - the 2-norm of the residual
- `monctx` - a `PetscViewer` of type `PETSCVIEWERDRAW` set up with `PetscViewerMonitorLGSetUp()`

Level: intermediate

See also: `SNES`, `SNESMonitorSet()`, `SNESMonitorDefault()`, `PetscViewerDrawGetDrawLG()`, `PetscDrawLG`

# External Links
$(_doc_external("SNES/SNESMonitorLGRange"))
"""
function SNESMonitorLGRange(petsclib::PetscLibType, snes::AbstractSNES, n::Integer, rnorm::Real, monctx::Ptr{Cvoid})
    error("SNESMonitorLGRange: no generated method for these argument types")
end

@for_petsc function SNESMonitorLGRange(petsclib::$UnionPetscLib, snes::AbstractSNES, n::$PetscInt, rnorm::$PetscReal, monctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESMonitorLGRange, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, Ptr{Cvoid}),
               snes, n, rnorm, monctx,
              )


	return nothing
end 

"""
	SNESMonitorRange(petsclib::PetscLibType, snes::AbstractSNES, it::PetscInt, rnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Prints the percentage of residual elements that are more than 10 percent of the maximum entry in the residual in each iteration of a `SNESSolve()`

Collective

Input Parameters:
- `snes`  - `SNES` iterative context
- `it`    - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - unused monitor context

Options Database Key:
- `-snes_monitor_range` - Activates `SNESMonitorRange()`

Level: intermediate

See also: `SNESMonitorSet()`, `SNESMonitorDefault()`, `SNESMonitorLGCreate()`, `SNESMonitorScaling()`, `PetscViewerFormat`, `PetscViewerAndFormat`

# External Links
$(_doc_external("SNES/SNESMonitorRange"))
"""
function SNESMonitorRange(petsclib::PetscLibType, snes::AbstractSNES, it::Integer, rnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("SNESMonitorRange: no generated method for these argument types")
end

@for_petsc function SNESMonitorRange(petsclib::$UnionPetscLib, snes::AbstractSNES, it::$PetscInt, rnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:SNESMonitorRange, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               snes, it, rnorm, vf,
              )


	return nothing
end 

"""
	SNESMonitorRatio(petsclib::PetscLibType, snes::AbstractSNES, its::PetscInt, fgnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Monitors progress of a `SNESSolve()` by printing the ratio of residual norm at each iteration to the previous.

Collective

Input Parameters:
- `snes`   - the `SNES` context
- `its`    - iteration number
- `fgnorm` - 2-norm of residual (or gradient)
- `vf`     - context of monitor

Options Database Key:
- `-snes_monitor_ratio` - activate this monitor

Level: intermediate

See also: `SNESMonitorRationSetUp()`, `SNESMonitorSet()`, `SNESMonitorSolution()`, `SNESMonitorDefault()`, `PetscViewerFormat`, `PetscViewerAndFormat`

# External Links
$(_doc_external("SNES/SNESMonitorRatio"))
"""
function SNESMonitorRatio(petsclib::PetscLibType, snes::AbstractSNES, its::Integer, fgnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("SNESMonitorRatio: no generated method for these argument types")
end

@for_petsc function SNESMonitorRatio(petsclib::$UnionPetscLib, snes::AbstractSNES, its::$PetscInt, fgnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:SNESMonitorRatio, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               snes, its, fgnorm, vf,
              )


	return nothing
end 

"""
	SNESMonitorRatioSetUp(petsclib::PetscLibType, snes::AbstractSNES, vf::Vector{PetscViewerAndFormat}) 
Insures the `SNES` object is saving its history since this monitor needs access to it

Collective

Input Parameters:
- `snes` - the `SNES` context
- `vf`   - `PetscViewerAndFormat` (ignored)

Level: intermediate

See also: `SNESMonitorSet()`, `SNESMonitorSolution()`, `SNESMonitorDefault()`, `SNESMonitorRatio()`, `PetscViewerFormat`, `PetscViewerAndFormat`

# External Links
$(_doc_external("SNES/SNESMonitorRatioSetUp"))
"""
function SNESMonitorRatioSetUp(petsclib::PetscLibType, snes::AbstractSNES, vf::Vector{PetscViewerAndFormat})
    error("SNESMonitorRatioSetUp: no generated method for these argument types")
end

@for_petsc function SNESMonitorRatioSetUp(petsclib::$UnionPetscLib, snes::AbstractSNES, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:SNESMonitorRatioSetUp, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{PetscViewerAndFormat}),
               snes, vf,
              )


	return nothing
end 

"""
	SNESMonitorResidual(petsclib::PetscLibType, snes::AbstractSNES, its::PetscInt, fgnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Monitors progress of a `SNESSolve()` by calling
`VecView()` for the residual at each iteration.

Collective

Input Parameters:
- `snes`   - the `SNES` context
- `its`    - iteration number
- `fgnorm` - 2-norm of residual
- `vf`     - a viewer

Options Database Key:
- `-snes_monitor_residual [ascii binary draw][:filename][:viewer format]` - plots residual (not its norm) at each iteration

Level: intermediate

See also: `SNES`, `SNESMonitorSet()`, `SNESMonitorDefault()`, `VecView()`, `SNESMonitor()`

# External Links
$(_doc_external("SNES/SNESMonitorResidual"))
"""
function SNESMonitorResidual(petsclib::PetscLibType, snes::AbstractSNES, its::Integer, fgnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("SNESMonitorResidual: no generated method for these argument types")
end

@for_petsc function SNESMonitorResidual(petsclib::$UnionPetscLib, snes::AbstractSNES, its::$PetscInt, fgnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:SNESMonitorResidual, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               snes, its, fgnorm, vf,
              )


	return nothing
end 

"""
	SNESMonitorSAWs(petsclib::PetscLibType, snes::AbstractSNES, n::PetscInt, rnorm::PetscReal, ctx::Ptr{Cvoid}) 
monitor solution process of `SNES` using SAWs

Collective

Input Parameters:
- `snes`  - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `ctx`   - `PetscViewer` of type `PETSCVIEWERSAWS`

Level: advanced

See also: `PetscViewerSAWsOpen()`, `SNESMonitorSAWsDestroy()`, `SNESMonitorSAWsCreate()`

# External Links
$(_doc_external("SNES/SNESMonitorSAWs"))
"""
function SNESMonitorSAWs(petsclib::PetscLibType, snes::AbstractSNES, n::Integer, rnorm::Real, ctx::Ptr{Cvoid})
    error("SNESMonitorSAWs: no generated method for these argument types")
end

@for_petsc function SNESMonitorSAWs(petsclib::$UnionPetscLib, snes::AbstractSNES, n::$PetscInt, rnorm::$PetscReal, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESMonitorSAWs, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, Ptr{Cvoid}),
               snes, n, rnorm, ctx,
              )


	return nothing
end 

"""
	ctx::Ptr{Cvoid} = SNESMonitorSAWsCreate(petsclib::PetscLibType, snes::AbstractSNES) 
create an SAWs monitor context for `SNES`

Collective

Input Parameter:
- `snes` - `SNES` to monitor

Output Parameter:
- `ctx` - context for monitor

Level: developer

See also: `SNESMonitorSet()`, `SNES`, `SNESMonitorSAWs()`, `SNESMonitorSAWsDestroy()`

# External Links
$(_doc_external("SNES/SNESMonitorSAWsCreate"))
"""
function SNESMonitorSAWsCreate(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESMonitorSAWsCreate: no generated method for these argument types")
end

@for_petsc function SNESMonitorSAWsCreate(petsclib::$UnionPetscLib, snes::AbstractSNES )
	ctx_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:SNESMonitorSAWsCreate, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Ptr{Cvoid}}),
               snes, ctx_,
              )

	ctx = ctx_[]

	return ctx
end 

"""
	SNESMonitorSAWsDestroy(petsclib::PetscLibType, ctx::Ptr{Cvoid}) 
destroy a monitor context created with `SNESMonitorSAWsCreate()`

Collective

Input Parameter:
- `ctx` - monitor context

Level: developer

See also: `SNESMonitorSAWsCreate()`

# External Links
$(_doc_external("SNES/SNESMonitorSAWsDestroy"))
"""
function SNESMonitorSAWsDestroy(petsclib::PetscLibType, ctx::Ptr{Cvoid})
    error("SNESMonitorSAWsDestroy: no generated method for these argument types")
end

@for_petsc function SNESMonitorSAWsDestroy(petsclib::$UnionPetscLib, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESMonitorSAWsDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{Cvoid},),
               ctx,
              )


	return nothing
end 

"""
	SNESMonitorScaling(petsclib::PetscLibType, snes::AbstractSNES, its::PetscInt, fgnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Monitors the largest value in each row of the Jacobian of a `SNESSolve()`

Collective

Input Parameters:
- `snes`   - the `SNES` context
- `its`    - iteration number
- `fgnorm` - 2-norm of residual
- `vf`     - viewer and format structure

Level: intermediate

See also: `SNESMonitorSet()`, `SNESMonitorSolution()`, `SNESMonitorRange()`, `SNESMonitorJacUpdateSpectrum()`,
`PetscViewerFormat`, `PetscViewerAndFormat`

# External Links
$(_doc_external("SNES/SNESMonitorScaling"))
"""
function SNESMonitorScaling(petsclib::PetscLibType, snes::AbstractSNES, its::Integer, fgnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("SNESMonitorScaling: no generated method for these argument types")
end

@for_petsc function SNESMonitorScaling(petsclib::$UnionPetscLib, snes::AbstractSNES, its::$PetscInt, fgnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:SNESMonitorScaling, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               snes, its, fgnorm, vf,
              )


	return nothing
end 

"""
	SNESMonitorSet(petsclib::PetscLibType, snes::AbstractSNES, f::external, mctx::Ptr{Cvoid}, monitordestroy::Ptr{Cvoid}) 
Sets an ADDITIONAL function that is to be used at every
iteration of the `SNES` nonlinear solver to display the iteration's
progress.

Logically Collective

Input Parameters:
- `snes`           - the `SNES` context
- `f`              - the monitor function,  for the calling sequence see `SNESMonitorFunction`
- `mctx`           - [optional] user-defined context for private data for the monitor routine (use `NULL` if no context is desired)
- `monitordestroy` - [optional] routine that frees monitor context (may be `NULL`), see `PetscCtxDestroyFn` for the calling sequence

Calling sequence of f:
- `snes`  - the `SNES` object
- `it`    - the current iteration
- `rnorm` - norm of the residual
- `mctx`  - the optional monitor context

Options Database Keys:
- `-snes_monitor`               - sets `SNESMonitorDefault()`
- `-snes_monitor draw::draw_lg` - sets line graph monitor
- `-snes_monitor_cancel`        - cancels all monitors that have been hardwired into a code by calls to `SNESMonitorSet()`, but does not cancel those set via
the options database.

Level: intermediate

See also: `SNES`, `SNESSolve()`, `SNESMonitorDefault()`, `SNESMonitorCancel()`, `SNESMonitorFunction`, `PetscCtxDestroyFn`

# External Links
$(_doc_external("SNES/SNESMonitorSet"))
"""
function SNESMonitorSet(petsclib::PetscLibType, snes::AbstractSNES, f::external, mctx::Ptr{Cvoid}, monitordestroy::Ptr{Cvoid})
    error("SNESMonitorSet: no generated method for these argument types")
end

@for_petsc function SNESMonitorSet(petsclib::$UnionPetscLib, snes::AbstractSNES, f::external, mctx::Ptr{Cvoid}, monitordestroy::Ptr{Cvoid} )

    @chk ccall(
               (:SNESMonitorSet, $petsc_library),
               PetscErrorCode,
               (CSNES, external, Ptr{Cvoid}, Ptr{Cvoid}),
               snes, f, mctx, monitordestroy,
              )


	return nothing
end 

"""
	SNESMonitorSetFromOptions(petsclib::PetscLibType, snes::AbstractSNES, name::String, help::String, manual::String, monitor::external, monitorsetup::external) 
Sets a monitor function and viewer appropriate for the type indicated by the user

Collective

Input Parameters:
- `snes`         - `SNES` object you wish to monitor
- `name`         - the monitor type one is seeking
- `help`         - message indicating what monitoring is done
- `manual`       - manual page for the monitor
- `monitor`      - the monitor function, this must use a `PetscViewerFormat` as its context
- `monitorsetup` - a function that is called once ONLY if the user selected this monitor that may set additional features of the `SNES` or `PetscViewer` objects

Calling sequence of `monitor`:
- `snes` - the nonlinear solver context
- `it`   - the current iteration
- `r`    - the current function norm
- `vf`   - a `PetscViewerAndFormat` struct that contains the `PetscViewer` and `PetscViewerFormat` to use

Calling sequence of `monitorsetup`:
- `snes` - the nonlinear solver context
- `vf`   - a `PetscViewerAndFormat` struct that contains the `PetscViewer` and `PetscViewerFormat` to use

Options Database Key:
- `-name` - trigger the use of this monitor in `SNESSetFromOptions()`

Level: advanced

See also: `PetscOptionsCreateViewer()`, `PetscOptionsGetReal()`, `PetscOptionsHasName()`, `PetscOptionsGetString()`,
`PetscOptionsGetIntArray()`, `PetscOptionsGetRealArray()`, `PetscOptionsBool()`,
`PetscOptionsInt()`, `PetscOptionsString()`, `PetscOptionsReal()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`

# External Links
$(_doc_external("SNES/SNESMonitorSetFromOptions"))
"""
function SNESMonitorSetFromOptions(petsclib::PetscLibType, snes::AbstractSNES, name::String, help::String, manual::String, monitor::external, monitorsetup::external)
    error("SNESMonitorSetFromOptions: no generated method for these argument types")
end

@for_petsc function SNESMonitorSetFromOptions(petsclib::$UnionPetscLib, snes::AbstractSNES, name::String, help::String, manual::String, monitor::external, monitorsetup::external )

    @chk ccall(
               (:SNESMonitorSetFromOptions, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, external, external),
               snes, name, help, manual, monitor, monitorsetup,
              )


	return nothing
end 

"""
	SNESMonitorSolution(petsclib::PetscLibType, snes::AbstractSNES, its::PetscInt, fgnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Monitors progress of a `SNES` `SNESSolve()` by calling
`VecView()` for the approximate solution at each iteration.

Collective

Input Parameters:
- `snes`   - the `SNES` context
- `its`    - iteration number
- `fgnorm` - 2-norm of residual
- `vf`     - a viewer

Options Database Key:
- `-snes_monitor_solution [ascii binary draw][:filename][:viewer format]` - plots solution at each iteration

Level: intermediate

See also: `SNES`, `SNESMonitorSet()`, `SNESMonitorDefault()`, `VecView()`

# External Links
$(_doc_external("SNES/SNESMonitorSolution"))
"""
function SNESMonitorSolution(petsclib::PetscLibType, snes::AbstractSNES, its::Integer, fgnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("SNESMonitorSolution: no generated method for these argument types")
end

@for_petsc function SNESMonitorSolution(petsclib::$UnionPetscLib, snes::AbstractSNES, its::$PetscInt, fgnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:SNESMonitorSolution, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               snes, its, fgnorm, vf,
              )


	return nothing
end 

"""
	SNESMonitorSolutionUpdate(petsclib::PetscLibType, snes::AbstractSNES, its::PetscInt, fgnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Monitors progress of a `SNESSolve()` by calling
`VecView()` for the UPDATE to the solution at each iteration.

Collective

Input Parameters:
- `snes`   - the `SNES` context
- `its`    - iteration number
- `fgnorm` - 2-norm of residual
- `vf`     - a viewer

Options Database Key:
- `-snes_monitor_solution_update [ascii binary draw][:filename][:viewer format]` - plots update to solution at each iteration

Level: intermediate

See also: `SNESMonitorSet()`, `SNESMonitorDefault()`, `VecView()`, `SNESMonitor()`

# External Links
$(_doc_external("SNES/SNESMonitorSolutionUpdate"))
"""
function SNESMonitorSolutionUpdate(petsclib::PetscLibType, snes::AbstractSNES, its::Integer, fgnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("SNESMonitorSolutionUpdate: no generated method for these argument types")
end

@for_petsc function SNESMonitorSolutionUpdate(petsclib::$UnionPetscLib, snes::AbstractSNES, its::$PetscInt, fgnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:SNESMonitorSolutionUpdate, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               snes, its, fgnorm, vf,
              )


	return nothing
end 

"""
	n::PetscInt,subsnes::Ptr{SNES} = SNESMultiblockGetSubSNES(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the `SNES` contexts for all blocks in a `SNESMULTIBLOCK` solver.

Not Collective but each `SNES` obtained is parallel

Input Parameter:
- `snes` - the solver context

Output Parameters:
- `n`       - the number of blocks
- `subsnes` - the array of `SNES` contexts

Level: advanced

See also: `SNES`, `SNESMULTIBLOCK`, `SNESMultiblockSetIS()`, `SNESMultiblockSetFields()`

# External Links
$(_doc_external("SNES/SNESMultiblockGetSubSNES"))
"""
function SNESMultiblockGetSubSNES(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESMultiblockGetSubSNES: no generated method for these argument types")
end

@for_petsc function SNESMultiblockGetSubSNES(petsclib::$UnionPetscLib, snes::AbstractSNES )
	n_ = Ref{$PetscInt}()
	subsnes_ = Ref{Ptr{SNES}}()

    @chk ccall(
               (:SNESMultiblockGetSubSNES, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscInt}, Ptr{Ptr{CSNES}}),
               snes, n_, subsnes_,
              )

	n = n_[]
	subsnes = subsnes_[]

	return n,subsnes
end 

"""
	SNESMultiblockSetBlockSize(petsclib::PetscLibType, snes::AbstractSNES, bs::PetscInt) 
Sets the block size for structured block division in a `SNESMULTIBLOCK` solver. If not set the matrix block size is used.

Logically Collective

Input Parameters:
- `snes` - the solver context
- `bs`   - the block size

Level: intermediate

See also: `SNES`, `SNESMULTIBLOCK`, `SNESMultiblockGetSubSNES()`, `SNESMultiblockSetFields()`

# External Links
$(_doc_external("SNES/SNESMultiblockSetBlockSize"))
"""
function SNESMultiblockSetBlockSize(petsclib::PetscLibType, snes::AbstractSNES, bs::Integer)
    error("SNESMultiblockSetBlockSize: no generated method for these argument types")
end

@for_petsc function SNESMultiblockSetBlockSize(petsclib::$UnionPetscLib, snes::AbstractSNES, bs::$PetscInt )

    @chk ccall(
               (:SNESMultiblockSetBlockSize, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt),
               snes, bs,
              )


	return nothing
end 

"""
	SNESMultiblockSetFields(petsclib::PetscLibType, snes::AbstractSNES, name::String, n::PetscInt, fields::Vector{PetscInt}) 
Sets the fields for one particular block in a `SNESMULTIBLOCK` solver

Logically Collective

Input Parameters:
- `snes`   - the solver
- `name`   - name of this block, if `NULL` the number of the block is used
- `n`      - the number of fields in this block
- `fields` - the fields in this block

Level: intermediate

See also: `SNES`, `SNESMULTIBLOCK`, `SNESMultiblockGetSubSNES()`, `SNESMultiblockSetBlockSize()`, `SNESMultiblockSetIS()`

# External Links
$(_doc_external("SNES/SNESMultiblockSetFields"))
"""
function SNESMultiblockSetFields(petsclib::PetscLibType, snes::AbstractSNES, name::String, n::Integer, fields::AbstractVector{<:Number})
    error("SNESMultiblockSetFields: no generated method for these argument types")
end

@for_petsc function SNESMultiblockSetFields(petsclib::$UnionPetscLib, snes::AbstractSNES, name::String, n::$PetscInt, fields::Vector{$PetscInt} )

    @chk ccall(
               (:SNESMultiblockSetFields, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Cchar}, $PetscInt, Ptr{$PetscInt}),
               snes, name, n, fields,
              )


	return nothing
end 

"""
	SNESMultiblockSetIS(petsclib::PetscLibType, snes::AbstractSNES, name::String, is::AbstractIS) 
Sets the global row indices for one particular block in a `SNESMULTIBLOCK` solver

Logically Collective

Input Parameters:
- `snes` - the solver context
- `name` - name of this block, if `NULL` the number of the block is used
- `is`   - the index set that defines the global row indices in this block

Level: intermediate

See also: `SNES`, `SNESMULTIBLOCK`, `SNESMultiblockGetSubSNES()`, `SNESMultiblockSetBlockSize()`, `SNESMultiblockSetFields()`

# External Links
$(_doc_external("SNES/SNESMultiblockSetIS"))
"""
function SNESMultiblockSetIS(petsclib::PetscLibType, snes::AbstractSNES, name::String, is::AbstractIS)
    error("SNESMultiblockSetIS: no generated method for these argument types")
end

@for_petsc function SNESMultiblockSetIS(petsclib::$UnionPetscLib, snes::AbstractSNES, name::String, is::AbstractIS )

    @chk ccall(
               (:SNESMultiblockSetIS, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Cchar}, CIS),
               snes, name, is,
              )


	return nothing
end 

"""
	SNESMultiblockSetType(petsclib::PetscLibType, snes::AbstractSNES, type::PCCompositeType) 
Sets the type of block combination used for a `SNESMULTIBLOCK` solver

Logically Collective

Input Parameters:
- `snes` - the solver context
- `type` - `PC_COMPOSITE_ADDITIVE`, `PC_COMPOSITE_MULTIPLICATIVE` (default), `PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE`

Options Database Key:
- `-snes_multiblock_type (multiplicative|additive|symmetric_multiplicative)` - Sets block combination type

Level: advanced

See also: `SNES`, `SNESMULTIBLOCK`, `PCCompositeSetType()`, `PC_COMPOSITE_ADDITIVE`, `PC_COMPOSITE_MULTIPLICATIVE`, `PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE`,
`PCCompositeType`, `SNESCOMPOSITE`, `SNESCompositeSetType()`

# External Links
$(_doc_external("SNES/SNESMultiblockSetType"))
"""
function SNESMultiblockSetType(petsclib::PetscLibType, snes::AbstractSNES, type::PCCompositeType)
    error("SNESMultiblockSetType: no generated method for these argument types")
end

@for_petsc function SNESMultiblockSetType(petsclib::$UnionPetscLib, snes::AbstractSNES, type::PCCompositeType )

    @chk ccall(
               (:SNESMultiblockSetType, $petsc_library),
               PetscErrorCode,
               (CSNES, PCCompositeType),
               snes, type,
              )


	return nothing
end 

"""
	dmp::PetscReal = SNESNASMGetDamping(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the update damping for `SNESNASM` the nonlinear additive Schwarz solver

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `dmp` - damping

Level: intermediate

See also: `SNES`, `SNESNASM`, `SNESNASMSetDamping()`

# External Links
$(_doc_external("SNES/SNESNASMGetDamping"))
"""
function SNESNASMGetDamping(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESNASMGetDamping: no generated method for these argument types")
end

@for_petsc function SNESNASMGetDamping(petsclib::$UnionPetscLib, snes::AbstractSNES )
	dmp_ = Ref{$PetscReal}()

    @chk ccall(
               (:SNESNASMGetDamping, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscReal}),
               snes, dmp_,
              )

	dmp = dmp_[]

	return dmp
end 

"""
	n::PetscInt = SNESNASMGetNumber(petsclib::PetscLibType, snes::AbstractSNES) 
Gets number of subsolvers

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `n` - the number of subsolvers

Level: intermediate

See also: `SNESNASM`, `SNESNASMGetSNES()`

# External Links
$(_doc_external("SNES/SNESNASMGetNumber"))
"""
function SNESNASMGetNumber(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESNASMGetNumber: no generated method for these argument types")
end

@for_petsc function SNESNASMGetNumber(petsclib::$UnionPetscLib, snes::AbstractSNES )
	n_ = Ref{$PetscInt}()

    @chk ccall(
               (:SNESNASMGetNumber, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscInt}),
               snes, n_,
              )

	n = n_[]

	return n
end 

"""
	subsnes::SNES = SNESNASMGetSNES(petsclib::PetscLibType, snes::AbstractSNES, i::PetscInt) 
Gets a subsolver

Not Collective

Input Parameters:
- `snes` - the `SNES` context
- `i`    - the number of the subsnes to get

Output Parameter:
- `subsnes` - the subsolver context

Level: intermediate

See also: `SNESNASM`, `SNESNASMGetNumber()`

# External Links
$(_doc_external("SNES/SNESNASMGetSNES"))
"""
function SNESNASMGetSNES(petsclib::PetscLibType, snes::AbstractSNES, i::Integer)
    error("SNESNASMGetSNES: no generated method for these argument types")
end

@for_petsc function SNESNASMGetSNES(petsclib::$UnionPetscLib, snes::AbstractSNES, i::$PetscInt )
	subsnes_ = Ref{CSNES}()

    @chk ccall(
               (:SNESNASMGetSNES, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, Ptr{CSNES}),
               snes, i, subsnes_,
              )

	subsnes = SNES(subsnes_[], petsclib)

	return subsnes
end 

"""
	n::PetscInt,x::Ptr{PetscVec},y::Ptr{PetscVec},b::Ptr{PetscVec},xl::Ptr{PetscVec} = SNESNASMGetSubdomainVecs(petsclib::PetscLibType, snes::AbstractSNES) 
Get the processor-local subdomain vectors for the nonlinear additive Schwarz solver

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameters:
- `n`  - the number of local subdomains
- `x`  - The subdomain solution vector
- `y`  - The subdomain step vector
- `b`  - The subdomain RHS vector
- `xl` - The subdomain local vectors (ghosted)

Level: developer

See also: `SNES`, `SNESNASM`, `SNESNASMGetSubdomains()`

# External Links
$(_doc_external("SNES/SNESNASMGetSubdomainVecs"))
"""
function SNESNASMGetSubdomainVecs(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESNASMGetSubdomainVecs: no generated method for these argument types")
end

@for_petsc function SNESNASMGetSubdomainVecs(petsclib::$UnionPetscLib, snes::AbstractSNES )
	n_ = Ref{$PetscInt}()
	x_ = Ref{Ptr{PetscVec}}()
	y_ = Ref{Ptr{PetscVec}}()
	b_ = Ref{Ptr{PetscVec}}()
	xl_ = Ref{Ptr{PetscVec}}()

    @chk ccall(
               (:SNESNASMGetSubdomainVecs, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscInt}, Ptr{Ptr{CVec}}, Ptr{Ptr{CVec}}, Ptr{Ptr{CVec}}, Ptr{Ptr{CVec}}),
               snes, n_, x_, y_, b_, xl_,
              )

	n = n_[]
	x = x_[]
	y = y_[]
	b = b_[]
	xl = xl_[]

	return n,x,y,b,xl
end 

"""
	n::PetscInt,subsnes::Ptr{SNES},iscatter::Ptr{VecScatter},oscatter::Ptr{VecScatter},gscatter::Ptr{VecScatter} = SNESNASMGetSubdomains(petsclib::PetscLibType, snes::AbstractSNES) 
Get the local subdomain contexts for the nonlinear additive Schwarz solver

Not Collective but some of the objects returned will be parallel

Input Parameter:
- `snes` - the `SNES` context

Output Parameters:
- `n`        - the number of local subdomains
- `subsnes`  - solvers defined on the local subdomains
- `iscatter` - scatters into the nonoverlapping portions of the local subdomains
- `oscatter` - scatters into the overlapping portions of the local subdomains
- `gscatter` - scatters into the (ghosted) local vector of the local subdomain

Level: intermediate

See also: `SNES`, `SNESNASM`, `SNESNASMSetSubdomains()`

# External Links
$(_doc_external("SNES/SNESNASMGetSubdomains"))
"""
function SNESNASMGetSubdomains(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESNASMGetSubdomains: no generated method for these argument types")
end

@for_petsc function SNESNASMGetSubdomains(petsclib::$UnionPetscLib, snes::AbstractSNES )
	n_ = Ref{$PetscInt}()
	subsnes_ = Ref{Ptr{SNES}}()
	iscatter_ = Ref{Ptr{VecScatter}}()
	oscatter_ = Ref{Ptr{VecScatter}}()
	gscatter_ = Ref{Ptr{VecScatter}}()

    @chk ccall(
               (:SNESNASMGetSubdomains, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscInt}, Ptr{Ptr{CSNES}}, Ptr{Ptr{VecScatter}}, Ptr{Ptr{VecScatter}}, Ptr{Ptr{VecScatter}}),
               snes, n_, subsnes_, iscatter_, oscatter_, gscatter_,
              )

	n = n_[]
	subsnes = subsnes_[]
	iscatter = iscatter_[]
	oscatter = oscatter_[]
	gscatter = gscatter_[]

	return n,subsnes,iscatter,oscatter,gscatter
end 

"""
	type::PCASMType = SNESNASMGetType(petsclib::PetscLibType, snes::AbstractSNES) 
Get the type of subdomain update used for the nonlinear additive Schwarz solver `SNESNASM`

Logically Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `type` - the type of update

Level: intermediate

See also: `SNES`, `SNESNASM`, `SNESNASMSetType()`, `PCASMGetType()`, `PC_ASM_BASIC`, `PC_ASM_RESTRICT`, `PCASMType`

# External Links
$(_doc_external("SNES/SNESNASMGetType"))
"""
function SNESNASMGetType(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESNASMGetType: no generated method for these argument types")
end

@for_petsc function SNESNASMGetType(petsclib::$UnionPetscLib, snes::AbstractSNES )
	type_ = Ref{PCASMType}()

    @chk ccall(
               (:SNESNASMGetType, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{PCASMType}),
               snes, type_,
              )

	type = type_[]

	return type
end 

"""
	SNESNASMSetComputeFinalJacobian(petsclib::PetscLibType, snes::AbstractSNES, flg::PetscBool) 
Schedules the computation of the global and subdomain Jacobians upon convergence for the
nonlinear additive Schwarz solver

Collective

Input Parameters:
- `snes` - the SNES context
- `flg`  - `PETSC_TRUE` to compute the Jacobians

Level: developer

See also: `SNES`, `SNESNASM`, `SNESNASMGetSubdomains()`

# External Links
$(_doc_external("SNES/SNESNASMSetComputeFinalJacobian"))
"""
function SNESNASMSetComputeFinalJacobian(petsclib::PetscLibType, snes::AbstractSNES, flg::PetscBool)
    error("SNESNASMSetComputeFinalJacobian: no generated method for these argument types")
end

@for_petsc function SNESNASMSetComputeFinalJacobian(petsclib::$UnionPetscLib, snes::AbstractSNES, flg::PetscBool )

    @chk ccall(
               (:SNESNASMSetComputeFinalJacobian, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, flg,
              )


	return nothing
end 

"""
	SNESNASMSetDamping(petsclib::PetscLibType, snes::AbstractSNES, dmp::PetscReal) 
Sets the update damping for `SNESNASM` the nonlinear additive Schwarz solver

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `dmp`  - damping

Options Database Key:
- `-snes_nasm_damping dmp` - the new solution is obtained as old solution plus `dmp` times (sum of the solutions on the subdomains)

Level: intermediate

See also: `SNES`, `SNESNASM`, `SNESNASMGetDamping()`

# External Links
$(_doc_external("SNES/SNESNASMSetDamping"))
"""
function SNESNASMSetDamping(petsclib::PetscLibType, snes::AbstractSNES, dmp::Real)
    error("SNESNASMSetDamping: no generated method for these argument types")
end

@for_petsc function SNESNASMSetDamping(petsclib::$UnionPetscLib, snes::AbstractSNES, dmp::$PetscReal )

    @chk ccall(
               (:SNESNASMSetDamping, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscReal),
               snes, dmp,
              )


	return nothing
end 

"""
	SNESNASMSetSubdomains(petsclib::PetscLibType, snes::AbstractSNES, n::PetscInt, subsnes::Vector{<:AbstractSNES}, iscatter::Vector{VecScatter}, oscatter::Vector{VecScatter}, gscatter::Vector{VecScatter}) 
Manually Set the context required to restrict and solve subdomain problems in the nonlinear additive Schwarz solver

Logically Collective

Input Parameters:
- `snes`     - the `SNES` context
- `n`        - the number of local subdomains
- `subsnes`  - solvers defined on the local subdomains
- `iscatter` - scatters into the nonoverlapping portions of the local subdomains
- `oscatter` - scatters into the overlapping portions of the local subdomains
- `gscatter` - scatters into the (ghosted) local vector of the local subdomain

Level: intermediate

See also: `SNES`, `SNESNASM`, `SNESNASMGetSubdomains()`

# External Links
$(_doc_external("SNES/SNESNASMSetSubdomains"))
"""
function SNESNASMSetSubdomains(petsclib::PetscLibType, snes::AbstractSNES, n::Integer, subsnes::Vector{<:AbstractSNES}, iscatter::Vector{VecScatter}, oscatter::Vector{VecScatter}, gscatter::Vector{VecScatter})
    error("SNESNASMSetSubdomains: no generated method for these argument types")
end

@for_petsc function SNESNASMSetSubdomains(petsclib::$UnionPetscLib, snes::AbstractSNES, n::$PetscInt, subsnes::Vector{<:AbstractSNES}, iscatter::Vector{VecScatter}, oscatter::Vector{VecScatter}, gscatter::Vector{VecScatter} )

    @chk ccall(
               (:SNESNASMSetSubdomains, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, Ptr{CSNES}, Ptr{VecScatter}, Ptr{VecScatter}, Ptr{VecScatter}),
               snes, n, subsnes, iscatter, oscatter, gscatter,
              )


	return nothing
end 

"""
	SNESNASMSetType(petsclib::PetscLibType, snes::AbstractSNES, type::PCASMType) 
Set the type of subdomain update used for the nonlinear additive Schwarz solver `SNESNASM`

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `type` - the type of update, `PC_ASM_BASIC` or `PC_ASM_RESTRICT`

Options Database Key:
- `-snes_nasm_type (basic|restrict)` - type of subdomain update used

Level: intermediate

See also: `SNES`, `SNESNASM`, `SNESNASMGetType()`, `PCASMSetType()`, `PC_ASM_BASIC`, `PC_ASM_RESTRICT`, `PCASMType`

# External Links
$(_doc_external("SNES/SNESNASMSetType"))
"""
function SNESNASMSetType(petsclib::PetscLibType, snes::AbstractSNES, type::PCASMType)
    error("SNESNASMSetType: no generated method for these argument types")
end

@for_petsc function SNESNASMSetType(petsclib::$UnionPetscLib, snes::AbstractSNES, type::PCASMType )

    @chk ccall(
               (:SNESNASMSetType, $petsc_library),
               PetscErrorCode,
               (CSNES, PCASMType),
               snes, type,
              )


	return nothing
end 

"""
	SNESNASMSetWeight(petsclib::PetscLibType, snes::AbstractSNES, weight::AbstractPetscVec) 
Sets weight to use when adding overlapping updates

Collective

Input Parameters:
- `snes`   - the `SNES` context
- `weight` - the weights to use (typically 1/N for each dof, where N is the number of patches it appears in)

Level: intermediate

See also: `SNESNASM`

# External Links
$(_doc_external("SNES/SNESNASMSetWeight"))
"""
function SNESNASMSetWeight(petsclib::PetscLibType, snes::AbstractSNES, weight::AbstractPetscVec)
    error("SNESNASMSetWeight: no generated method for these argument types")
end

@for_petsc function SNESNASMSetWeight(petsclib::$UnionPetscLib, snes::AbstractSNES, weight::AbstractPetscVec )

    @chk ccall(
               (:SNESNASMSetWeight, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec),
               snes, weight,
              )


	return nothing
end 

"""
	SNESNCGSetType(petsclib::PetscLibType, snes::AbstractSNES, btype::SNESNCGType) 
Sets the conjugate update type for nonlinear CG `SNESNCG`.

Logically Collective

Input Parameters:
- `snes`  - the iterative context
- `btype` - update type, see `SNESNCGType`

Options Database Key:
- `-snes_ncg_type (prp|fr|hs|dy|cd)` - strategy for selecting algorithm for computing beta

Level: intermediate

See also: `SNES`, `SNESNCG`, `SNESNCGType`, `SNES_NCG_FR`, `SNES_NCG_PRP`, `SNES_NCG_HS`, `SNES_NCG_DY`, `SNES_NCG_CD`

# External Links
$(_doc_external("SNES/SNESNCGSetType"))
"""
function SNESNCGSetType(petsclib::PetscLibType, snes::AbstractSNES, btype::SNESNCGType)
    error("SNESNCGSetType: no generated method for these argument types")
end

@for_petsc function SNESNCGSetType(petsclib::$UnionPetscLib, snes::AbstractSNES, btype::SNESNCGType )

    @chk ccall(
               (:SNESNCGSetType, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESNCGType),
               snes, btype,
              )


	return nothing
end 

"""
	flg::PetscBool = SNESNGMRESGetRestartFmRise(petsclib::PetscLibType, snes::AbstractSNES) 
Get whether `SNESNGMRES` increases the restart count when a step x_M increases the residual F_M

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `flg` - `PETSC_TRUE` if the option is enabled

Level: advanced

See also: `SNES`, `SNESNGMRES`, `SNESNGMRESSetRestartFmRise()`, `SNESNGMRESSetRestartType()`

# External Links
$(_doc_external("SNES/SNESNGMRESGetRestartFmRise"))
"""
function SNESNGMRESGetRestartFmRise(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESNGMRESGetRestartFmRise: no generated method for these argument types")
end

@for_petsc function SNESNGMRESGetRestartFmRise(petsclib::$UnionPetscLib, snes::AbstractSNES )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:SNESNGMRESGetRestartFmRise, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{PetscBool}),
               snes, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	SNESNGMRESSetRestartFmRise(petsclib::PetscLibType, snes::AbstractSNES, flg::PetscBool) 
Increase the restart count if the step x_M increases the residual F_M inside a `SNESNGMRES` solve

Input Parameters:
- `snes` - the `SNES` context.
- `flg`  - boolean value deciding whether to use the option or not, default is `PETSC_FALSE`

Options Database Key:
- `-snes_ngmres_restart_fm_rise (true|false)` - Increase the restart count if the step x_M increases the residual F_M

Level: advanced

See also: `SNES`, `SNES_NGMRES_RESTART_DIFFERENCE`, `SNESNGMRES`, `SNESNGMRESRestartType`, `SNESNGMRESSetRestartType()`

# External Links
$(_doc_external("SNES/SNESNGMRESSetRestartFmRise"))
"""
function SNESNGMRESSetRestartFmRise(petsclib::PetscLibType, snes::AbstractSNES, flg::PetscBool)
    error("SNESNGMRESSetRestartFmRise: no generated method for these argument types")
end

@for_petsc function SNESNGMRESSetRestartFmRise(petsclib::$UnionPetscLib, snes::AbstractSNES, flg::PetscBool )

    @chk ccall(
               (:SNESNGMRESSetRestartFmRise, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, flg,
              )


	return nothing
end 

"""
	SNESNGMRESSetRestartType(petsclib::PetscLibType, snes::AbstractSNES, rtype::SNESNGMRESRestartType) 
Sets the restart type for `SNESNGMRES`.

Logically Collective

Input Parameters:
- `snes`  - the iterative context
- `rtype` - restart type, see `SNESNGMRESRestartType`

Options Database Keys:
- `-snes_ngmres_restart_type (difference|periodic|none)` - set the restart type
- `-snes_ngmres_restart restart`                         - sets the number of iterations before restart for periodic

Level: intermediate

See also: `SNES`, `SNES_NGMRES_RESTART_DIFFERENCE`, `SNESNGMRES`, `SNESNGMRESRestartType`, `SNESNGMRESSetRestartFmRise()`,
`SNESNGMRESSetSelectType()`

# External Links
$(_doc_external("SNES/SNESNGMRESSetRestartType"))
"""
function SNESNGMRESSetRestartType(petsclib::PetscLibType, snes::AbstractSNES, rtype::SNESNGMRESRestartType)
    error("SNESNGMRESSetRestartType: no generated method for these argument types")
end

@for_petsc function SNESNGMRESSetRestartType(petsclib::$UnionPetscLib, snes::AbstractSNES, rtype::SNESNGMRESRestartType )

    @chk ccall(
               (:SNESNGMRESSetRestartType, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESNGMRESRestartType),
               snes, rtype,
              )


	return nothing
end 

"""
	SNESNGMRESSetSelectType(petsclib::PetscLibType, snes::AbstractSNES, stype::SNESNGMRESSelectType) 
Sets the selection type for `SNESNGMRES`.  This determines how the candidate solution and
combined solution are used to create the next iterate.

Logically Collective

Input Parameters:
- `snes`  - the iterative context
- `stype` - selection type, see `SNESNGMRESSelectType`

Options Database Key:
- `-snes_ngmres_select_type (difference|none|linesearch)` - select type

Level: intermediate

See also: `SNES`, `SNESNGMRES`, `SNESNGMRESSelectType`, `SNES_NGMRES_SELECT_NONE`, `SNES_NGMRES_SELECT_DIFFERENCE`, `SNES_NGMRES_SELECT_LINESEARCH`,
`SNESNGMRESSetRestartType()`

# External Links
$(_doc_external("SNES/SNESNGMRESSetSelectType"))
"""
function SNESNGMRESSetSelectType(petsclib::PetscLibType, snes::AbstractSNES, stype::SNESNGMRESSelectType)
    error("SNESNGMRESSetSelectType: no generated method for these argument types")
end

@for_petsc function SNESNGMRESSetSelectType(petsclib::$UnionPetscLib, snes::AbstractSNES, stype::SNESNGMRESSelectType )

    @chk ccall(
               (:SNESNGMRESSetSelectType, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESNGMRESSelectType),
               snes, stype,
              )


	return nothing
end 

"""
	sweeps::PetscInt = SNESNGSGetSweeps(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the number of sweeps nonlinear GS will use in `SNESNCG`

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `sweeps` - the number of sweeps of nonlinear GS to perform.

Level: intermediate

See also: `SNES`, `SNESNCG`, `SNESSetNGS()`, `SNESGetNGS()`, `SNESSetNPC()`, `SNESNGSSetSweeps()`

# External Links
$(_doc_external("SNES/SNESNGSGetSweeps"))
"""
function SNESNGSGetSweeps(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESNGSGetSweeps: no generated method for these argument types")
end

@for_petsc function SNESNGSGetSweeps(petsclib::$UnionPetscLib, snes::AbstractSNES )
	sweeps_ = Ref{$PetscInt}()

    @chk ccall(
               (:SNESNGSGetSweeps, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscInt}),
               snes, sweeps_,
              )

	sweeps = sweeps_[]

	return sweeps
end 

"""
	atol::PetscReal,rtol::PetscReal,stol::PetscReal,maxit::PetscInt = SNESNGSGetTolerances(petsclib::PetscLibType, snes::AbstractSNES) 
Gets various parameters used in convergence tests for nonlinear Gauss-Seidel `SNESNCG`

Not Collective

Input Parameters:
- `snes`  - the `SNES` context
- `atol`  - absolute convergence tolerance
- `rtol`  - relative convergence tolerance
- `stol`  - convergence tolerance in terms of the norm
of the change in the solution between steps
- `maxit` - maximum number of iterations

Level: intermediate

See also: `SNES`, `SNESNCG`, `SNESSetTolerances()`

# External Links
$(_doc_external("SNES/SNESNGSGetTolerances"))
"""
function SNESNGSGetTolerances(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESNGSGetTolerances: no generated method for these argument types")
end

@for_petsc function SNESNGSGetTolerances(petsclib::$UnionPetscLib, snes::AbstractSNES )
	atol_ = Ref{$PetscReal}()
	rtol_ = Ref{$PetscReal}()
	stol_ = Ref{$PetscReal}()
	maxit_ = Ref{$PetscInt}()

    @chk ccall(
               (:SNESNGSGetTolerances, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscInt}),
               snes, atol_, rtol_, stol_, maxit_,
              )

	atol = atol_[]
	rtol = rtol_[]
	stol = stol_[]
	maxit = maxit_[]

	return atol,rtol,stol,maxit
end 

"""
	SNESNGSSetSweeps(petsclib::PetscLibType, snes::AbstractSNES, sweeps::PetscInt) 
Sets the number of sweeps of nonlinear GS to use in `SNESNCG`

Logically Collective

Input Parameters:
- `snes`   - the `SNES` context
- `sweeps` - the number of sweeps of nonlinear GS to perform.

Options Database Key:
- `-snes_ngs_sweeps n` - Number of sweeps of nonlinear GS to apply

Level: intermediate

See also: `SNES`, `SNESNCG`, `SNESSetNGS()`, `SNESGetNGS()`, `SNESSetNPC()`, `SNESNGSGetSweeps()`

# External Links
$(_doc_external("SNES/SNESNGSSetSweeps"))
"""
function SNESNGSSetSweeps(petsclib::PetscLibType, snes::AbstractSNES, sweeps::Integer)
    error("SNESNGSSetSweeps: no generated method for these argument types")
end

@for_petsc function SNESNGSSetSweeps(petsclib::$UnionPetscLib, snes::AbstractSNES, sweeps::$PetscInt )

    @chk ccall(
               (:SNESNGSSetSweeps, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt),
               snes, sweeps,
              )


	return nothing
end 

"""
	SNESNGSSetTolerances(petsclib::PetscLibType, snes::AbstractSNES, abstol::PetscReal, rtol::PetscReal, stol::PetscReal, maxit::PetscInt) 
Sets various parameters used in convergence tests for nonlinear Gauss-Seidel `SNESNCG`

Logically Collective

Input Parameters:
- `snes`   - the `SNES` context
- `abstol` - absolute convergence tolerance
- `rtol`   - relative convergence tolerance
- `stol`   - convergence tolerance in terms of the norm of the change in the solution between steps,  || delta x || < stol*|| x ||
- `maxit`  - maximum number of iterations

Options Database Keys:
- `-snes_ngs_atol abstol` - Sets abstol
- `-snes_ngs_rtol rtol`   - Sets rtol
- `-snes_ngs_stol stol`   - Sets stol
- `-snes_max_it maxit`    - Sets maxit

Level: intermediate

See also: `SNES`, `SNESNCG`

# External Links
$(_doc_external("SNES/SNESNGSSetTolerances"))
"""
function SNESNGSSetTolerances(petsclib::PetscLibType, snes::AbstractSNES, abstol::Real, rtol::Real, stol::Real, maxit::Integer)
    error("SNESNGSSetTolerances: no generated method for these argument types")
end

@for_petsc function SNESNGSSetTolerances(petsclib::$UnionPetscLib, snes::AbstractSNES, abstol::$PetscReal, rtol::$PetscReal, stol::$PetscReal, maxit::$PetscInt )

    @chk ccall(
               (:SNESNGSSetTolerances, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscReal, $PetscReal, $PetscReal, $PetscInt),
               snes, abstol, rtol, stol, maxit,
              )


	return nothing
end 

"""
	SNESNewtonALComputeFunction(petsclib::PetscLibType, snes::AbstractSNES, X::AbstractPetscVec, Q::AbstractPetscVec) 
Calls the function that has been set with `SNESNewtonALSetFunction()`.

Collective

Input Parameters:
- `snes` - the `SNES` context
- `X`    - input vector

Output Parameter:
- `Q` - tangent load vector, as set by `SNESNewtonALSetFunction()`

Level: developer

See also: `SNES`, `SNESNewtonALSetFunction()`, `SNESNewtonALGetFunction()`

# External Links
$(_doc_external("SNES/SNESNewtonALComputeFunction"))
"""
function SNESNewtonALComputeFunction(petsclib::PetscLibType, snes::AbstractSNES, X::AbstractPetscVec, Q::AbstractPetscVec)
    error("SNESNewtonALComputeFunction: no generated method for these argument types")
end

@for_petsc function SNESNewtonALComputeFunction(petsclib::$UnionPetscLib, snes::AbstractSNES, X::AbstractPetscVec, Q::AbstractPetscVec )

    @chk ccall(
               (:SNESNewtonALComputeFunction, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec),
               snes, X, Q,
              )


	return nothing
end 

"""
	SNESNewtonALGetFunction(petsclib::PetscLibType, snes::AbstractSNES, func::Ptr{Ptr{Cvoid}}, ctx::Ptr{Cvoid}) 
Get the user function and context set with `SNESNewtonALSetFunction`

Logically Collective

Input Parameters:
- `snes` - the nonlinear solver object
- `func` - [optional] tangent load function evaluation routine, see `SNESNewtonALSetFunction()` for the call sequence
- `ctx`  - [optional] user-defined context for private data for the function evaluation routine (may be `NULL`)

Level: intermediate

See also: `SNES`, `SNESNEWTONAL`, `SNESNewtonALSetFunction()`

# External Links
$(_doc_external("SNES/SNESNewtonALGetFunction"))
"""
function SNESNewtonALGetFunction(petsclib::PetscLibType, snes::AbstractSNES, func::Ptr{Ptr{Cvoid}}, ctx::Ptr{Cvoid})
    error("SNESNewtonALGetFunction: no generated method for these argument types")
end

@for_petsc function SNESNewtonALGetFunction(petsclib::$UnionPetscLib, snes::AbstractSNES, func::Ptr{Ptr{Cvoid}}, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESNewtonALGetFunction, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Ptr{Cvoid}}, Ptr{Cvoid}),
               snes, func, ctx,
              )


	return nothing
end 

"""
	lambda::PetscReal = SNESNewtonALGetLoadParameter(petsclib::PetscLibType, snes::AbstractSNES) 
Get the value of the load parameter `lambda` for the arc-length continuation method.

Logically Collective

Input Parameter:
- `snes` - the nonlinear solver object

Output Parameter:
- `lambda` - the arc-length parameter

Level: intermediate

See also: `SNES`, `SNESNEWTONAL`, `SNESNewtonALSetFunction()`

# External Links
$(_doc_external("SNES/SNESNewtonALGetLoadParameter"))
"""
function SNESNewtonALGetLoadParameter(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESNewtonALGetLoadParameter: no generated method for these argument types")
end

@for_petsc function SNESNewtonALGetLoadParameter(petsclib::$UnionPetscLib, snes::AbstractSNES )
	lambda_ = Ref{$PetscReal}()

    @chk ccall(
               (:SNESNewtonALGetLoadParameter, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscReal}),
               snes, lambda_,
              )

	lambda = lambda_[]

	return lambda
end 

"""
	SNESNewtonALSetCorrectionType(petsclib::PetscLibType, snes::AbstractSNES, ctype::SNESNewtonALCorrectionType) 
Set the type of correction to use in the arc-length continuation method.

Logically Collective

Input Parameters:
- `snes`  - the nonlinear solver object
- `ctype` - the type of correction to use

Options Database Key:
- `-snes_newtonal_correction_type type` - Set the type of correction to use; use -help for a list of available types

Level: intermediate

See also: `SNES`, `SNESNEWTONAL`, `SNESNewtonALCorrectionType`

# External Links
$(_doc_external("SNES/SNESNewtonALSetCorrectionType"))
"""
function SNESNewtonALSetCorrectionType(petsclib::PetscLibType, snes::AbstractSNES, ctype::SNESNewtonALCorrectionType)
    error("SNESNewtonALSetCorrectionType: no generated method for these argument types")
end

@for_petsc function SNESNewtonALSetCorrectionType(petsclib::$UnionPetscLib, snes::AbstractSNES, ctype::SNESNewtonALCorrectionType )

    @chk ccall(
               (:SNESNewtonALSetCorrectionType, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESNewtonALCorrectionType),
               snes, ctype,
              )


	return nothing
end 

"""
	SNESNewtonALSetDiagonalScaling(petsclib::PetscLibType, snes::AbstractSNES, v::AbstractPetscVec) 
Set the global vector used to rescale DoFs for computation of arc length.

Logically Collective

Input Parameters:
- `snes` - the nonlinear solver object
- `v`    - the `Vec` containing diagonal scaling for each DoF, must be the same size as the solution vector (may be `NULL`)

See also: `SNES`, `SNESNEWTONAL`, `SNESNewtonALSetFunction()`, `SNESNewtonALGetLoadParameter()`

# External Links
$(_doc_external("SNES/SNESNewtonALSetDiagonalScaling"))
"""
function SNESNewtonALSetDiagonalScaling(petsclib::PetscLibType, snes::AbstractSNES, v::AbstractPetscVec)
    error("SNESNewtonALSetDiagonalScaling: no generated method for these argument types")
end

@for_petsc function SNESNewtonALSetDiagonalScaling(petsclib::$UnionPetscLib, snes::AbstractSNES, v::AbstractPetscVec )

    @chk ccall(
               (:SNESNewtonALSetDiagonalScaling, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec),
               snes, v,
              )


	return nothing
end 

"""
	SNESNewtonALSetFunction(petsclib::PetscLibType, snes::AbstractSNES, func::Ptr{Cvoid}, ctx::Ptr{Cvoid}) 
Sets a user function that is called at each function evaluation to
compute the tangent load vector for the arc-length continuation method.

Logically Collective

Input Parameters:
- `snes` - the nonlinear solver object
- `func` - [optional] tangent load function evaluation routine, see `SNESFunctionFn` for the calling sequence. `U` is the current solution vector, `Q` is the output tangent load vector
- `ctx`  - [optional] user-defined context for private data for the function evaluation routine (may be `NULL`)

Level: intermediate

See also: `SNES`, `SNESNEWTONAL`, `SNESNewtonALGetFunction()`, `SNESNewtonALGetLoadParameter()`

# External Links
$(_doc_external("SNES/SNESNewtonALSetFunction"))
"""
function SNESNewtonALSetFunction(petsclib::PetscLibType, snes::AbstractSNES, func::Ptr{Cvoid}, ctx::Ptr{Cvoid})
    error("SNESNewtonALSetFunction: no generated method for these argument types")
end

@for_petsc function SNESNewtonALSetFunction(petsclib::$UnionPetscLib, snes::AbstractSNES, func::Ptr{Cvoid}, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESNewtonALSetFunction, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Cvoid}, Ptr{Cvoid}),
               snes, func, ctx,
              )


	return nothing
end 

"""
	SNESNewtonTRDCGetPostCheck(petsclib::PetscLibType, snes::AbstractSNES, noname::Ptr{Cvoid}) 
Gets the post-check function optionally set with `SNESNewtonTRDCSetPostCheck()`

Not Collective

Input Parameter:
- `snes` - the nonlinear solver context

Output Parameters:
- `func` - [optional] function evaluation routine, for the calling sequence see `SNESNewtonTRDCPostCheck()`
- `ctx`  - [optional] context for private data for the function evaluation routine (may be `NULL`)

Calling sequence of `func`:
- `snes`      - the nonlinear solver object
- `X`         - the current solution value
- `Y`         - the tentative update step
- `W`         - the tentative new solution value
- `changed_y` - output, flag indicated `Y` has been changed by the post-check
- `changed_w` - output, flag indicated `W` has been changed by the post-check
- `ctx`       - the optional application context

Level: intermediate

See also: `SNES`, `SNESNEWTONTRDC`, `SNESNewtonTRDCSetPostCheck()`, `SNESNewtonTRDCPostCheck()`, `SNESNewtonTRDCSetPreCheck()`, `SNESNewtonTRDCGetPreCheck()`

# External Links
$(_doc_external("SNES/SNESNewtonTRDCGetPostCheck"))
"""
function SNESNewtonTRDCGetPostCheck(petsclib::PetscLibType, snes::AbstractSNES, noname::Ptr{Cvoid})
    error("SNESNewtonTRDCGetPostCheck: no generated method for these argument types")
end

@for_petsc function SNESNewtonTRDCGetPostCheck(petsclib::$UnionPetscLib, snes::AbstractSNES, noname::Ptr{Cvoid} )

    @chk ccall(
               (:SNESNewtonTRDCGetPostCheck, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Cvoid}),
               snes, noname,
              )


	return nothing
end 

"""
	SNESNewtonTRDCGetPreCheck(petsclib::PetscLibType, snes::AbstractSNES, noname::Ptr{Cvoid}) 
Gets the pre-check function optionally set with `SNESNewtonTRDCSetPreCheck()`

Not Collective

Input Parameter:
- `snes` - the nonlinear solver context

Output Parameters:
- `func` - [optional] function evaluation routine, for the calling sequence see `SNESNewtonTRDCPreCheck()`
- `ctx`  - [optional] context for private data for the function evaluation routine (may be `NULL`)

Calling sequence of `func`:
- `snes`    - the nonlinear solver object
- `X`       - the current solution value
- `Y`       - the tentative update step
- `changed` - output, flag indicating `Y` has been changed by the pre-check
- `ctx`     - the optional application context

Level: intermediate

See also: `SNES`, `SNESNEWTONTRDC`, `SNESNewtonTRDCSetPreCheck()`, `SNESNewtonTRDCPreCheck()`

# External Links
$(_doc_external("SNES/SNESNewtonTRDCGetPreCheck"))
"""
function SNESNewtonTRDCGetPreCheck(petsclib::PetscLibType, snes::AbstractSNES, noname::Ptr{Cvoid})
    error("SNESNewtonTRDCGetPreCheck: no generated method for these argument types")
end

@for_petsc function SNESNewtonTRDCGetPreCheck(petsclib::$UnionPetscLib, snes::AbstractSNES, noname::Ptr{Cvoid} )

    @chk ccall(
               (:SNESNewtonTRDCGetPreCheck, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Cvoid}),
               snes, noname,
              )


	return nothing
end 

"""
	rho_flag::PetscBool = SNESNewtonTRDCGetRhoFlag(petsclib::PetscLibType, snes::AbstractSNES) 
Get whether the current solution update is within the trust-region.

Logically Collective

Input Parameter:
- `snes` - the nonlinear solver object

Output Parameter:
- `rho_flag` - `PETSC_FALSE` or `PETSC_TRUE`

Level: developer

See also: `SNES`, `SNESNEWTONTRDC`, `SNESNewtonTRDCPreCheck()`, `SNESNewtonTRDCGetPreCheck()`, `SNESNewtonTRDCSetPreCheck()`,
`SNESNewtonTRDCSetPostCheck()`, `SNESNewtonTRDCGetPostCheck()`

# External Links
$(_doc_external("SNES/SNESNewtonTRDCGetRhoFlag"))
"""
function SNESNewtonTRDCGetRhoFlag(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESNewtonTRDCGetRhoFlag: no generated method for these argument types")
end

@for_petsc function SNESNewtonTRDCGetRhoFlag(petsclib::$UnionPetscLib, snes::AbstractSNES )
	rho_flag_ = Ref{PetscBool}()

    @chk ccall(
               (:SNESNewtonTRDCGetRhoFlag, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{PetscBool}),
               snes, rho_flag_,
              )

	rho_flag = rho_flag_[]

	return rho_flag
end 

"""
	SNESNewtonTRDCSetPostCheck(petsclib::PetscLibType, snes::AbstractSNES, func::external, ctx::Ptr{Cvoid}) 
Sets a user function that is called after the search step has been determined but before the next
function evaluation. Allows the user a chance to change or override the decision of the line search routine

Logically Collective

Input Parameters:
- `snes` - the nonlinear solver object
- `func` - [optional] function evaluation routine, for the calling sequence see `SNESNewtonTRDCPostCheck()`
- `ctx`  - [optional] context for private data for the function evaluation routine (may be `NULL`)

Calling sequence of `func`:
- `snes`      - the nonlinear solver object
- `X`         - the current solution value
- `Y`         - the tentative update step
- `W`         - the tentative new solution value
- `changed_y` - output, flag indicated `Y` has been changed by the post-check
- `changed_w` - output, flag indicated `W` has been changed by the post-check
- `ctx`       - the optional application context

Level: intermediate

See also: `SNES`, `SNESNEWTONTRDC`, `SNESNewtonTRDCPostCheck()`, `SNESNewtonTRDCGetPostCheck()`, `SNESNewtonTRDCSetPreCheck()`, `SNESNewtonTRDCGetPreCheck()`

# External Links
$(_doc_external("SNES/SNESNewtonTRDCSetPostCheck"))
"""
function SNESNewtonTRDCSetPostCheck(petsclib::PetscLibType, snes::AbstractSNES, func::external, ctx::Ptr{Cvoid})
    error("SNESNewtonTRDCSetPostCheck: no generated method for these argument types")
end

@for_petsc function SNESNewtonTRDCSetPostCheck(petsclib::$UnionPetscLib, snes::AbstractSNES, func::external, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESNewtonTRDCSetPostCheck, $petsc_library),
               PetscErrorCode,
               (CSNES, external, Ptr{Cvoid}),
               snes, func, ctx,
              )


	return nothing
end 

"""
	SNESNewtonTRDCSetPreCheck(petsclib::PetscLibType, snes::AbstractSNES, func::external, ctx::Ptr{Cvoid}) 
Sets a user function that is called before the search step has been determined.
Allows the user a chance to change or override the trust region decision.

Logically Collective

Input Parameters:
- `snes` - the nonlinear solver object
- `func` - [optional] function evaluation routine, for the calling sequence see `SNESNewtonTRDCPreCheck()`
- `ctx`  - [optional] application context for private data for the function evaluation routine (may be `NULL`)

Calling sequence of `func`:
- `snes`    - the nonlinear solver object
- `X`       - the current solution value
- `Y`       - the tentative update step
- `changed` - output, flag indicating `Y` has been changed by the pre-check
- `ctx`     - the optional application context

Level: intermediate

See also: `SNES`, `SNESNEWTONTRDC`, `SNESNewtonTRDCPreCheck()`, `SNESNewtonTRDCGetPreCheck()`, `SNESNewtonTRDCSetPostCheck()`, `SNESNewtonTRDCGetPostCheck()`,
`SNESNewtonTRDCGetRhoFlag()`

# External Links
$(_doc_external("SNES/SNESNewtonTRDCSetPreCheck"))
"""
function SNESNewtonTRDCSetPreCheck(petsclib::PetscLibType, snes::AbstractSNES, func::external, ctx::Ptr{Cvoid})
    error("SNESNewtonTRDCSetPreCheck: no generated method for these argument types")
end

@for_petsc function SNESNewtonTRDCSetPreCheck(petsclib::$UnionPetscLib, snes::AbstractSNES, func::external, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESNewtonTRDCSetPreCheck, $petsc_library),
               PetscErrorCode,
               (CSNES, external, Ptr{Cvoid}),
               snes, func, ctx,
              )


	return nothing
end 

"""
	SNESNewtonTRGetPostCheck(petsclib::PetscLibType, snes::AbstractSNES, noname::Ptr{Cvoid}) 
Gets the post-check function

Not Collective

Input Parameter:
- `snes` - the nonlinear solver context

Output Parameters:
- `func` - [optional] function evaluation routine, for the calling sequence see `SNESNewtonTRPostCheck()`
- `ctx`  - [optional] user-defined context for private data for the function evaluation routine (may be `NULL`)

Calling sequence of `func`:
- `snes`      - the nonlinear solver object
- `X`         - the current solution value
- `Y`         - the tentative update step
- `W`         - the tentative new solution value
- `changed_Y` - output, flag indicated `Y` has been changed by the post-check
- `changed_W` - output, flag indicated `W` has been changed by the post-check
- `ctx`       - the optional application context

Level: intermediate

See also: `SNESNEWTONTR`, `SNESNewtonTRSetPostCheck()`, `SNESNewtonTRPostCheck()`

# External Links
$(_doc_external("SNES/SNESNewtonTRGetPostCheck"))
"""
function SNESNewtonTRGetPostCheck(petsclib::PetscLibType, snes::AbstractSNES, noname::Ptr{Cvoid})
    error("SNESNewtonTRGetPostCheck: no generated method for these argument types")
end

@for_petsc function SNESNewtonTRGetPostCheck(petsclib::$UnionPetscLib, snes::AbstractSNES, noname::Ptr{Cvoid} )

    @chk ccall(
               (:SNESNewtonTRGetPostCheck, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Cvoid}),
               snes, noname,
              )


	return nothing
end 

"""
	SNESNewtonTRGetPreCheck(petsclib::PetscLibType, snes::AbstractSNES, noname::Ptr{Cvoid}) 
Gets the pre-check function

Not Collective

Input Parameter:
- `snes` - the nonlinear solver context

Output Parameters:
- `func` - [optional] function evaluation routine, for the calling sequence see `SNESNewtonTRPreCheck()`
- `ctx`  - [optional] user-defined context for private data for the function evaluation routine (may be `NULL`)

Calling sequence of `func`:
- `snes`    - the nonlinear solver object
- `X`       - the current solution value
- `Y`       - the tentative update step
- `changed` - output, flag indicating `Y` has been changed by the pre-check
- `ctx`     - the optional application context

Level: intermediate

See also: `SNESNEWTONTR`, `SNESNewtonTRSetPreCheck()`, `SNESNewtonTRPreCheck()`

# External Links
$(_doc_external("SNES/SNESNewtonTRGetPreCheck"))
"""
function SNESNewtonTRGetPreCheck(petsclib::PetscLibType, snes::AbstractSNES, noname::Ptr{Cvoid})
    error("SNESNewtonTRGetPreCheck: no generated method for these argument types")
end

@for_petsc function SNESNewtonTRGetPreCheck(petsclib::$UnionPetscLib, snes::AbstractSNES, noname::Ptr{Cvoid} )

    @chk ccall(
               (:SNESNewtonTRGetPreCheck, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Cvoid}),
               snes, noname,
              )


	return nothing
end 

"""
	delta_min::PetscReal,delta_max::PetscReal,delta_0::PetscReal = SNESNewtonTRGetTolerances(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the trust region parameter tolerances.

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameters:
- `delta_min` - minimum allowed trust region size or `NULL`
- `delta_max` - maximum allowed trust region size or `NULL`
- `delta_0`   - initial trust region size or `NULL`

Level: intermediate

See also: `SNES`, `SNESNEWTONTR`, `SNESNewtonTRSetTolerances()`

# External Links
$(_doc_external("SNES/SNESNewtonTRGetTolerances"))
"""
function SNESNewtonTRGetTolerances(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESNewtonTRGetTolerances: no generated method for these argument types")
end

@for_petsc function SNESNewtonTRGetTolerances(petsclib::$UnionPetscLib, snes::AbstractSNES )
	delta_min_ = Ref{$PetscReal}()
	delta_max_ = Ref{$PetscReal}()
	delta_0_ = Ref{$PetscReal}()

    @chk ccall(
               (:SNESNewtonTRGetTolerances, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               snes, delta_min_, delta_max_, delta_0_,
              )

	delta_min = delta_min_[]
	delta_max = delta_max_[]
	delta_0 = delta_0_[]

	return delta_min,delta_max,delta_0
end 

"""
	eta1::PetscReal,eta2::PetscReal,eta3::PetscReal,t1::PetscReal,t2::PetscReal = SNESNewtonTRGetUpdateParameters(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the trust region update parameters.

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameters:
- `eta1` - acceptance tolerance
- `eta2` - shrinking tolerance
- `eta3` - enlarging tolerance
- `t1`   - shrink factor
- `t2`   - enlarge factor

Level: intermediate

See also: `SNES`, `SNESNEWTONTR`, `SNESNewtonTRSetUpdateParameters()`

# External Links
$(_doc_external("SNES/SNESNewtonTRGetUpdateParameters"))
"""
function SNESNewtonTRGetUpdateParameters(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESNewtonTRGetUpdateParameters: no generated method for these argument types")
end

@for_petsc function SNESNewtonTRGetUpdateParameters(petsclib::$UnionPetscLib, snes::AbstractSNES )
	eta1_ = Ref{$PetscReal}()
	eta2_ = Ref{$PetscReal}()
	eta3_ = Ref{$PetscReal}()
	t1_ = Ref{$PetscReal}()
	t2_ = Ref{$PetscReal}()

    @chk ccall(
               (:SNESNewtonTRGetUpdateParameters, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               snes, eta1_, eta2_, eta3_, t1_, t2_,
              )

	eta1 = eta1_[]
	eta2 = eta2_[]
	eta3 = eta3_[]
	t1 = t1_[]
	t2 = t2_[]

	return eta1,eta2,eta3,t1,t2
end 

"""
	changed_Y::PetscBool,changed_W::PetscBool = SNESNewtonTRPostCheck(petsclib::PetscLibType, snes::AbstractSNES, X::AbstractPetscVec, Y::AbstractPetscVec, W::AbstractPetscVec) 
Runs the postcheck routine

Logically Collective

Input Parameters:
- `snes` - the solver
- `X`    - The last solution
- `Y`    - The full step direction
- `W`    - The updated solution, W = X - Y

Output Parameters:
- `changed_Y` - indicator if step has been changed
- `changed_W` - Indicator if the new candidate solution W has been changed.

See also: `SNESNEWTONTR`, `SNESNewtonTRSetPostCheck()`, `SNESNewtonTRGetPostCheck()`, `SNESNewtonTRPreCheck()`

# External Links
$(_doc_external("SNES/SNESNewtonTRPostCheck"))
"""
function SNESNewtonTRPostCheck(petsclib::PetscLibType, snes::AbstractSNES, X::AbstractPetscVec, Y::AbstractPetscVec, W::AbstractPetscVec)
    error("SNESNewtonTRPostCheck: no generated method for these argument types")
end

@for_petsc function SNESNewtonTRPostCheck(petsclib::$UnionPetscLib, snes::AbstractSNES, X::AbstractPetscVec, Y::AbstractPetscVec, W::AbstractPetscVec )
	changed_Y_ = Ref{PetscBool}()
	changed_W_ = Ref{PetscBool}()

    @chk ccall(
               (:SNESNewtonTRPostCheck, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec, CVec, Ptr{PetscBool}, Ptr{PetscBool}),
               snes, X, Y, W, changed_Y_, changed_W_,
              )

	changed_Y = changed_Y_[]
	changed_W = changed_W_[]

	return changed_Y,changed_W
end 

"""
	changed_Y::PetscBool = SNESNewtonTRPreCheck(petsclib::PetscLibType, snes::AbstractSNES, X::AbstractPetscVec, Y::AbstractPetscVec) 
Runs the precheck routine

Logically Collective

Input Parameters:
- `snes` - the solver
- `X`    - The last solution
- `Y`    - The step direction

Output Parameter:
- `changed_Y` - Indicator that the step direction `Y` has been changed.

Level: intermediate

See also: `SNESNEWTONTR`, `SNESNewtonTRSetPreCheck()`, `SNESNewtonTRGetPreCheck()`, `SNESNewtonTRPostCheck()`

# External Links
$(_doc_external("SNES/SNESNewtonTRPreCheck"))
"""
function SNESNewtonTRPreCheck(petsclib::PetscLibType, snes::AbstractSNES, X::AbstractPetscVec, Y::AbstractPetscVec)
    error("SNESNewtonTRPreCheck: no generated method for these argument types")
end

@for_petsc function SNESNewtonTRPreCheck(petsclib::$UnionPetscLib, snes::AbstractSNES, X::AbstractPetscVec, Y::AbstractPetscVec )
	changed_Y_ = Ref{PetscBool}()

    @chk ccall(
               (:SNESNewtonTRPreCheck, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec, Ptr{PetscBool}),
               snes, X, Y, changed_Y_,
              )

	changed_Y = changed_Y_[]

	return changed_Y
end 

"""
	SNESNewtonTRSetFallbackType(petsclib::PetscLibType, snes::AbstractSNES, ftype::SNESNewtonTRFallbackType) 
Set the type of fallback to use if the solution of the trust region subproblem is outside the radius

Input Parameters:
- `snes`  - the nonlinear solver object
- `ftype` - the fallback type, see `SNESNewtonTRFallbackType`

Level: intermediate

See also: `SNESNEWTONTR`, `SNESNewtonTRPreCheck()`, `SNESNewtonTRGetPreCheck()`, `SNESNewtonTRSetPreCheck()`,
`SNESNewtonTRSetPostCheck()`, `SNESNewtonTRGetPostCheck()`

# External Links
$(_doc_external("SNES/SNESNewtonTRSetFallbackType"))
"""
function SNESNewtonTRSetFallbackType(petsclib::PetscLibType, snes::AbstractSNES, ftype::SNESNewtonTRFallbackType)
    error("SNESNewtonTRSetFallbackType: no generated method for these argument types")
end

@for_petsc function SNESNewtonTRSetFallbackType(petsclib::$UnionPetscLib, snes::AbstractSNES, ftype::SNESNewtonTRFallbackType )

    @chk ccall(
               (:SNESNewtonTRSetFallbackType, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESNewtonTRFallbackType),
               snes, ftype,
              )


	return nothing
end 

"""
	SNESNewtonTRSetNormType(petsclib::PetscLibType, snes::AbstractSNES, norm::NormType) 
Specify the type of norm to use for the computation of the trust region.

Input Parameters:
- `snes` - the nonlinear solver object
- `norm` - the norm type

Level: intermediate

See also: `SNESNEWTONTR`, `NormType`

# External Links
$(_doc_external("SNES/SNESNewtonTRSetNormType"))
"""
function SNESNewtonTRSetNormType(petsclib::PetscLibType, snes::AbstractSNES, norm::NormType)
    error("SNESNewtonTRSetNormType: no generated method for these argument types")
end

@for_petsc function SNESNewtonTRSetNormType(petsclib::$UnionPetscLib, snes::AbstractSNES, norm::NormType )

    @chk ccall(
               (:SNESNewtonTRSetNormType, $petsc_library),
               PetscErrorCode,
               (CSNES, NormType),
               snes, norm,
              )


	return nothing
end 

"""
	SNESNewtonTRSetPostCheck(petsclib::PetscLibType, snes::AbstractSNES, func::external, ctx::Ptr{Cvoid}) 
Sets a user function that is called after the search step has been determined but before the next
function evaluation. Allows the user a chance to change or override the internal decision of the solver

Logically Collective

Input Parameters:
- `snes` - the nonlinear solver object
- `func` - [optional] function evaluation routine, for the calling sequence see `SNESNewtonTRPostCheck()`
- `ctx`  - [optional] user-defined context for private data for the function evaluation routine (may be `NULL`)

Calling sequence of `func`:
- `snes`      - the nonlinear solver object
- `X`         - the current solution value
- `Y`         - the tentative update step
- `W`         - the tentative new solution value
- `changed_Y` - output, flag indicated `Y` has been changed by the post-check
- `changed_W` - output, flag indicated `W` has been changed by the post-check
- `ctx`       - the optional application context

Level: intermediate

See also: `SNESNEWTONTR`, `SNESNewtonTRPostCheck()`, `SNESNewtonTRGetPostCheck()`, `SNESNewtonTRSetPreCheck()`, `SNESNewtonTRGetPreCheck()`

# External Links
$(_doc_external("SNES/SNESNewtonTRSetPostCheck"))
"""
function SNESNewtonTRSetPostCheck(petsclib::PetscLibType, snes::AbstractSNES, func::external, ctx::Ptr{Cvoid})
    error("SNESNewtonTRSetPostCheck: no generated method for these argument types")
end

@for_petsc function SNESNewtonTRSetPostCheck(petsclib::$UnionPetscLib, snes::AbstractSNES, func::external, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESNewtonTRSetPostCheck, $petsc_library),
               PetscErrorCode,
               (CSNES, external, Ptr{Cvoid}),
               snes, func, ctx,
              )


	return nothing
end 

"""
	SNESNewtonTRSetPreCheck(petsclib::PetscLibType, snes::AbstractSNES, func::external, ctx::Ptr{Cvoid}) 
Sets a user function that is called before the search step has been determined.
Allows the user a chance to change or override the trust region decision.

Logically Collective

Input Parameters:
- `snes` - the nonlinear solver object
- `func` - [optional] function evaluation routine, for the calling sequence see `SNESNewtonTRPreCheck()`
- `ctx`  - [optional] user-defined context for private data for the function evaluation routine (may be `NULL`)

Calling sequence of `func`:
- `snes`    - the nonlinear solver object
- `X`       - the current solution value
- `Y`       - the tentative update step
- `changed` - output, flag indicating `Y` has been changed by the pre-check
- `ctx`     - the optional application context

Level: intermediate

See also: `SNESNEWTONTR`, `SNESNewtonTRPreCheck()`, `SNESNewtonTRGetPreCheck()`, `SNESNewtonTRSetPostCheck()`, `SNESNewtonTRGetPostCheck()`

# External Links
$(_doc_external("SNES/SNESNewtonTRSetPreCheck"))
"""
function SNESNewtonTRSetPreCheck(petsclib::PetscLibType, snes::AbstractSNES, func::external, ctx::Ptr{Cvoid})
    error("SNESNewtonTRSetPreCheck: no generated method for these argument types")
end

@for_petsc function SNESNewtonTRSetPreCheck(petsclib::$UnionPetscLib, snes::AbstractSNES, func::external, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESNewtonTRSetPreCheck, $petsc_library),
               PetscErrorCode,
               (CSNES, external, Ptr{Cvoid}),
               snes, func, ctx,
              )


	return nothing
end 

"""
	SNESNewtonTRSetQNType(petsclib::PetscLibType, snes::AbstractSNES, use::SNESNewtonTRQNType) 
Specify to use a quasi-Newton model.

Input Parameters:
- `snes` - the nonlinear solver object
- `use`  - the type of approximations to be used

Level: intermediate

See also: `SNESNEWTONTR`, `SNESNewtonTRQNType`, `MATLMVM`

# External Links
$(_doc_external("SNES/SNESNewtonTRSetQNType"))
"""
function SNESNewtonTRSetQNType(petsclib::PetscLibType, snes::AbstractSNES, use::SNESNewtonTRQNType)
    error("SNESNewtonTRSetQNType: no generated method for these argument types")
end

@for_petsc function SNESNewtonTRSetQNType(petsclib::$UnionPetscLib, snes::AbstractSNES, use::SNESNewtonTRQNType )

    @chk ccall(
               (:SNESNewtonTRSetQNType, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESNewtonTRQNType),
               snes, use,
              )


	return nothing
end 

"""
	SNESNewtonTRSetTolerances(petsclib::PetscLibType, snes::AbstractSNES, delta_min::PetscReal, delta_max::PetscReal, delta_0::PetscReal) 
Sets the trust region parameter tolerances.

Logically Collective

Input Parameters:
- `snes`      - the `SNES` context
- `delta_min` - minimum allowed trust region size
- `delta_max` - maximum allowed trust region size
- `delta_0`   - initial trust region size

Options Database Key:
- `-snes_tr_deltamin tol` - Set minimum size
- `-snes_tr_deltamax tol` - Set maximum size
- `-snes_tr_delta0   tol` - Set initial size

See also: `SNES`, `SNESNEWTONTR`, `SNESNewtonTRGetTolerances()`

# External Links
$(_doc_external("SNES/SNESNewtonTRSetTolerances"))
"""
function SNESNewtonTRSetTolerances(petsclib::PetscLibType, snes::AbstractSNES, delta_min::Real, delta_max::Real, delta_0::Real)
    error("SNESNewtonTRSetTolerances: no generated method for these argument types")
end

@for_petsc function SNESNewtonTRSetTolerances(petsclib::$UnionPetscLib, snes::AbstractSNES, delta_min::$PetscReal, delta_max::$PetscReal, delta_0::$PetscReal )

    @chk ccall(
               (:SNESNewtonTRSetTolerances, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscReal, $PetscReal, $PetscReal),
               snes, delta_min, delta_max, delta_0,
              )


	return nothing
end 

"""
	SNESNewtonTRSetUpdateParameters(petsclib::PetscLibType, snes::AbstractSNES, eta1::PetscReal, eta2::PetscReal, eta3::PetscReal, t1::PetscReal, t2::PetscReal) 
Sets the trust region update parameters.

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `eta1` - acceptance tolerance
- `eta2` - shrinking tolerance
- `eta3` - enlarging tolerance
- `t1`   - shrink factor
- `t2`   - enlarge factor

Options Database Key:
- `-snes_tr_eta1 tol` - Set `eta1`
- `-snes_tr_eta2 tol` - Set `eta2`
- `-snes_tr_eta3 tol` - Set `eta3`
- `-snes_tr_t1   tol` - Set `t1`
- `-snes_tr_t2   tol` - Set `t2`

See also: `SNES`, `SNESNEWTONTR`, `SNESSetObjective()`, `SNESNewtonTRGetUpdateParameters()`

# External Links
$(_doc_external("SNES/SNESNewtonTRSetUpdateParameters"))
"""
function SNESNewtonTRSetUpdateParameters(petsclib::PetscLibType, snes::AbstractSNES, eta1::Real, eta2::Real, eta3::Real, t1::Real, t2::Real)
    error("SNESNewtonTRSetUpdateParameters: no generated method for these argument types")
end

@for_petsc function SNESNewtonTRSetUpdateParameters(petsclib::$UnionPetscLib, snes::AbstractSNES, eta1::$PetscReal, eta2::$PetscReal, eta3::$PetscReal, t1::$PetscReal, t2::$PetscReal )

    @chk ccall(
               (:SNESNewtonTRSetUpdateParameters, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscReal, $PetscReal, $PetscReal, $PetscReal, $PetscReal),
               snes, eta1, eta2, eta3, t1, t2,
              )


	return nothing
end 

"""
	SNESObjectiveComputeFunctionDefaultFD(petsclib::PetscLibType, snes::AbstractSNES, X::AbstractPetscVec, F::AbstractPetscVec, ctx::Ptr{Cvoid}) 
Computes the gradient of a user provided objective function

Collective

Input Parameters:
- `snes` - the `SNES` context
- `X`    - the state vector
- `ctx`  - the (ignored) function context

Output Parameter:
- `F` - the function value

Options Database Keys:
- `-snes_fd_function_eps` - Tolerance for including non-zero entries into the gradient, default is 1.e-6
- `-snes_fd_function`     - Computes function from user provided objective function (set with `SNESSetObjective()`) with finite difference

Level: advanced

See also: `SNESSetObjective()`, `SNESSetFunction()`, `SNESComputeObjective()`, `SNESComputeJacobianDefault()`, `SNESObjectiveFn`

# External Links
$(_doc_external("SNES/SNESObjectiveComputeFunctionDefaultFD"))
"""
function SNESObjectiveComputeFunctionDefaultFD(petsclib::PetscLibType, snes::AbstractSNES, X::AbstractPetscVec, F::AbstractPetscVec, ctx::Ptr{Cvoid})
    error("SNESObjectiveComputeFunctionDefaultFD: no generated method for these argument types")
end

@for_petsc function SNESObjectiveComputeFunctionDefaultFD(petsclib::$UnionPetscLib, snes::AbstractSNES, X::AbstractPetscVec, F::AbstractPetscVec, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESObjectiveComputeFunctionDefaultFD, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec, Ptr{Cvoid}),
               snes, X, F, ctx,
              )


	return nothing
end 

"""
	SNESParametersInitialize(petsclib::PetscLibType, snes::AbstractSNES) 
Sets all the parameters in `snes` to their default value (when `SNESCreate()` was called) if they
currently contain default values

Collective

Input Parameter:
- `snes` - the `SNES` object

Level: developer

See also: `SNES`, `SNESSolve()`, `SNESDestroy()`, `SNESSetLagPreconditioner()`, `SNESSetLagJacobian()`,
`PetscObjectParameterSetDefault()`

# External Links
$(_doc_external("SNES/SNESParametersInitialize"))
"""
function SNESParametersInitialize(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESParametersInitialize: no generated method for these argument types")
end

@for_petsc function SNESParametersInitialize(petsclib::$UnionPetscLib, snes::AbstractSNES )

    @chk ccall(
               (:SNESParametersInitialize, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	SNESPatchSetCellNumbering(petsclib::PetscLibType, snes::AbstractSNES, cellNumbering::PetscSection) 
Set the `PetscSection` that provides a numbering of the cells used to define patches in a `SNESPATCH` solver

Logically Collective

Input Parameters:
- `snes`          - the `SNESPATCH` solver
- `cellNumbering` - the `PetscSection` giving the cell numbering; forwarded to the underlying `PCPATCH` via `PCPatchSetCellNumbering()`

Level: advanced

See also: `SNESPATCH`, `PCPATCH`, `PCPatchSetCellNumbering()`, `PetscSection`

# External Links
$(_doc_external("SNES/SNESPatchSetCellNumbering"))
"""
function SNESPatchSetCellNumbering(petsclib::PetscLibType, snes::AbstractSNES, cellNumbering::PetscSection)
    error("SNESPatchSetCellNumbering: no generated method for these argument types")
end

@for_petsc function SNESPatchSetCellNumbering(petsclib::$UnionPetscLib, snes::AbstractSNES, cellNumbering::PetscSection )

    @chk ccall(
               (:SNESPatchSetCellNumbering, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscSection),
               snes, cellNumbering,
              )


	return nothing
end 

"""
	SNESPatchSetComputeFunction(petsclib::PetscLibType, snes::AbstractSNES, func::external, ctx::Ptr{Cvoid}) 
Set the callback used to compute the per-patch nonlinear residual for a `SNESPATCH` solver

Logically Collective

Input Parameters:
- `snes` - the `SNESPATCH` solver
- `func` - callback that computes the patch residual; forwarded to the underlying `PCPATCH` via `PCPatchSetComputeFunction()`
- `ctx`  - optional application context passed to `func`

Calling sequence of `func`:
- `pc`               - the `PC` associated with the `SNESPATCH` solver
- `point`            - the point
- `x`                - the input solution (not used in linear problems)
- `f`                - the patch residual vector
- `cellIS`           - an array of the cell numbers
- `n`                - the size of `dofsArray`
- `dofsArray`        - the dofmap for the dofs to be solved for
- `dofsArrayWithAll` - the dofmap for all dofs on the patch
- `ctx`              - the application context

Level: advanced

See also: `SNESPATCH`, `PCPATCH`, `PCPatchSetComputeFunction()`, `SNESPatchSetComputeOperator()`

# External Links
$(_doc_external("SNES/SNESPatchSetComputeFunction"))
"""
function SNESPatchSetComputeFunction(petsclib::PetscLibType, snes::AbstractSNES, func::external, ctx::Ptr{Cvoid})
    error("SNESPatchSetComputeFunction: no generated method for these argument types")
end

@for_petsc function SNESPatchSetComputeFunction(petsclib::$UnionPetscLib, snes::AbstractSNES, func::external, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESPatchSetComputeFunction, $petsc_library),
               PetscErrorCode,
               (CSNES, external, Ptr{Cvoid}),
               snes, func, ctx,
              )


	return nothing
end 

"""
	SNESPatchSetComputeOperator(petsclib::PetscLibType, snes::AbstractSNES, func::external, ctx::Ptr{Cvoid}) 
Set the callback used to assemble the per-patch Jacobian for a `SNESPATCH` solver

Logically Collective

Input Parameters:
- `snes` - the `SNESPATCH` solver
- `func` - callback that assembles the patch Jacobian; forwarded to the underlying `PCPATCH` via `PCPatchSetComputeOperator()`
- `ctx`  - optional application context passed to `func`

Calling sequence of `func`:
- `pc`               - the `PC` associated with the `SNESPATCH` solver
- `point`            - the point
- `x`                - the input solution (not used in linear problems)
- `mat`              - the patch matrix
- `cellIS`           - an array of the cell numbers
- `n`                - the size of `dofsArray`
- `dofsArray`        - the dofmap for the dofs to be solved for
- `dofsArrayWithAll` - the dofmap for all dofs on the patch
- `ctx`              - the application context

Level: advanced

See also: `SNESPATCH`, `PCPATCH`, `PCPatchSetComputeOperator()`, `SNESPatchSetComputeFunction()`

# External Links
$(_doc_external("SNES/SNESPatchSetComputeOperator"))
"""
function SNESPatchSetComputeOperator(petsclib::PetscLibType, snes::AbstractSNES, func::external, ctx::Ptr{Cvoid})
    error("SNESPatchSetComputeOperator: no generated method for these argument types")
end

@for_petsc function SNESPatchSetComputeOperator(petsclib::$UnionPetscLib, snes::AbstractSNES, func::external, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESPatchSetComputeOperator, $petsc_library),
               PetscErrorCode,
               (CSNES, external, Ptr{Cvoid}),
               snes, func, ctx,
              )


	return nothing
end 

"""
	SNESPatchSetConstructType(petsclib::PetscLibType, snes::AbstractSNES, ctype::PCPatchConstructType, func::external, ctx::Ptr{Cvoid}) 
Set the way patches are constructed for a `SNESPATCH` solver

Logically Collective

Input Parameters:
- `snes`  - the `SNESPATCH` solver
- `ctype` - the `PCPatchConstructType` selecting the patch construction strategy
- `func`  - user callback that builds the patches, used only when `ctype` is `PC_PATCH_USER` or `PC_PATCH_PYTHON`; may be `NULL` otherwise
- `ctx`   - optional application context passed to `func`

Calling sequence of `func`:
- `pc`                - the `PC` associated with the `SNESPATCH` solver
- `npatch`            - number of patches
- `patches`           - the `IS` that define each patch
- `patchIterationSet` - how the patches are iterated over
- `ctx`               - optional application context

Level: advanced

See also: `SNESPATCH`, `PCPATCH`, `PCPatchSetConstructType()`, `PCPatchConstructType`

# External Links
$(_doc_external("SNES/SNESPatchSetConstructType"))
"""
function SNESPatchSetConstructType(petsclib::PetscLibType, snes::AbstractSNES, ctype::PCPatchConstructType, func::external, ctx::Ptr{Cvoid})
    error("SNESPatchSetConstructType: no generated method for these argument types")
end

@for_petsc function SNESPatchSetConstructType(petsclib::$UnionPetscLib, snes::AbstractSNES, ctype::PCPatchConstructType, func::external, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESPatchSetConstructType, $petsc_library),
               PetscErrorCode,
               (CSNES, PCPatchConstructType, external, Ptr{Cvoid}),
               snes, ctype, func, ctx,
              )


	return nothing
end 

"""
	SNESPatchSetDiscretisationInfo(petsclib::PetscLibType, snes::AbstractSNES, nsubspaces::PetscInt, dms::Vector{<:AbstractPetscDM}, bs::Vector{PetscInt}, nodesPerCell::Vector{PetscInt}, cellNodeMap::PetscInt, subspaceOffsets::Vector{PetscInt}, numGhostBcs::PetscInt, ghostBcNodes::Vector{PetscInt}, numGlobalBcs::PetscInt, globalBcNodes::Vector{PetscInt}) 
Provide the per-subspace discretisation information required by a `SNESPATCH` to build patch problems

Logically Collective

Input Parameters:
- `snes`            - the `SNESPATCH` solver
- `nsubspaces`      - the number of discretisation subspaces (e.g. fields)
- `dms`             - array of length `nsubspaces` of `DM`s, one per subspace
- `bs`              - array of length `nsubspaces` giving the block size of each subspace
- `nodesPerCell`    - array of length `nsubspaces` giving the number of nodes per cell for each subspace
- `cellNodeMap`     - array of length `nsubspaces`; entry `i` is a cell-to-node map for subspace `i`
- `subspaceOffsets` - array of length `nsubspaces + 1` giving the starting global dof offset of each subspace
- `numGhostBcs`     - number of ghost (off-process) boundary-condition dofs
- `ghostBcNodes`    - array of length `numGhostBcs` of the ghost boundary-condition dof indices
- `numGlobalBcs`    - number of global boundary-condition dofs
- `globalBcNodes`   - array of length `numGlobalBcs` of the global boundary-condition dof indices

Level: advanced

See also: `SNESPATCH`, `PCPATCH`, `PCPatchSetDiscretisationInfo()`, `SNESPatchSetComputeOperator()`, `SNESPatchSetComputeFunction()`

# External Links
$(_doc_external("SNES/SNESPatchSetDiscretisationInfo"))
"""
function SNESPatchSetDiscretisationInfo(petsclib::PetscLibType, snes::AbstractSNES, nsubspaces::Integer, dms::Vector{<:AbstractPetscDM}, bs::AbstractVector{<:Number}, nodesPerCell::AbstractVector{<:Number}, cellNodeMap::Integer, subspaceOffsets::AbstractVector{<:Number}, numGhostBcs::Integer, ghostBcNodes::AbstractVector{<:Number}, numGlobalBcs::Integer, globalBcNodes::AbstractVector{<:Number})
    error("SNESPatchSetDiscretisationInfo: no generated method for these argument types")
end

@for_petsc function SNESPatchSetDiscretisationInfo(petsclib::$UnionPetscLib, snes::AbstractSNES, nsubspaces::$PetscInt, dms::Vector{<:AbstractPetscDM}, bs::Vector{$PetscInt}, nodesPerCell::Vector{$PetscInt}, cellNodeMap::$PetscInt, subspaceOffsets::Vector{$PetscInt}, numGhostBcs::$PetscInt, ghostBcNodes::Vector{$PetscInt}, numGlobalBcs::$PetscInt, globalBcNodes::Vector{$PetscInt} )

    @chk ccall(
               (:SNESPatchSetDiscretisationInfo, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, Ptr{CDM}, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}),
               snes, nsubspaces, dms, bs, nodesPerCell, cellNodeMap, subspaceOffsets, numGhostBcs, ghostBcNodes, numGlobalBcs, globalBcNodes,
              )


	return nothing
end 

"""
	SNESPicardComputeFunction(petsclib::PetscLibType, snes::AbstractSNES, x::AbstractPetscVec, f::AbstractPetscVec, ctx::Ptr{Cvoid}) 
Compute the residual A(x) x - b(x) using the callbacks registered by `SNESSetPicard()`

Collective

Input Parameters:
- `snes` - the `SNES` context
- `x`    - the current iterate
- `ctx`  - unused application context; the Picard callbacks are retrieved from the attached `DMSNES`

Output Parameter:
- `f` - the residual vector

Level: developer

See also: `SNES`, `SNESSetPicard()`, `SNESPicardComputeMFFunction()`, `SNESPicardComputeJacobian()`

# External Links
$(_doc_external("SNES/SNESPicardComputeFunction"))
"""
function SNESPicardComputeFunction(petsclib::PetscLibType, snes::AbstractSNES, x::AbstractPetscVec, f::AbstractPetscVec, ctx::Ptr{Cvoid})
    error("SNESPicardComputeFunction: no generated method for these argument types")
end

@for_petsc function SNESPicardComputeFunction(petsclib::$UnionPetscLib, snes::AbstractSNES, x::AbstractPetscVec, f::AbstractPetscVec, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESPicardComputeFunction, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec, Ptr{Cvoid}),
               snes, x, f, ctx,
              )


	return nothing
end 

"""
	SNESPicardComputeJacobian(petsclib::PetscLibType, snes::AbstractSNES, x1::AbstractPetscVec, J::AbstractPetscMat, B::AbstractPetscMat, ctx::Ptr{Cvoid}) 
Trivial Jacobian assembly callback used by `SNESSetPicard()`; the Picard operator is filled in by `SNESPicardComputeFunction()`

Collective

Input Parameters:
- `snes` - the `SNES` context
- `x1`   - the current iterate (unused)
- `J`    - the Jacobian matrix to assemble
- `B`    - the preconditioning matrix (unused)
- `ctx`  - unused application context

Level: developer

See also: `SNES`, `SNESSetPicard()`, `SNESPicardComputeFunction()`, `SNESPicardComputeMFFunction()`

# External Links
$(_doc_external("SNES/SNESPicardComputeJacobian"))
"""
function SNESPicardComputeJacobian(petsclib::PetscLibType, snes::AbstractSNES, x1::AbstractPetscVec, J::AbstractPetscMat, B::AbstractPetscMat, ctx::Ptr{Cvoid})
    error("SNESPicardComputeJacobian: no generated method for these argument types")
end

@for_petsc function SNESPicardComputeJacobian(petsclib::$UnionPetscLib, snes::AbstractSNES, x1::AbstractPetscVec, J::AbstractPetscMat, B::AbstractPetscMat, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESPicardComputeJacobian, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CMat, CMat, Ptr{Cvoid}),
               snes, x1, J, B, ctx,
              )


	return nothing
end 

"""
	SNESPicardComputeMFFunction(petsclib::PetscLibType, snes::AbstractSNES, x::AbstractPetscVec, f::AbstractPetscVec, ctx::Ptr{Cvoid}) 
Matrix-free residual A(x) x - b(x) used by `SNESSetPicard()` when the operator is applied through `-snes_mf_operator`

Collective

Input Parameters:
- `snes` - the `SNES` context
- `x`    - the current iterate
- `ctx`  - unused application context; the Picard callbacks are retrieved from the attached `DMSNES`

Output Parameter:
- `f` - the residual vector

Level: developer

See also: `SNES`, `SNESSetPicard()`, `SNESPicardComputeFunction()`, `SNESPicardComputeJacobian()`

# External Links
$(_doc_external("SNES/SNESPicardComputeMFFunction"))
"""
function SNESPicardComputeMFFunction(petsclib::PetscLibType, snes::AbstractSNES, x::AbstractPetscVec, f::AbstractPetscVec, ctx::Ptr{Cvoid})
    error("SNESPicardComputeMFFunction: no generated method for these argument types")
end

@for_petsc function SNESPicardComputeMFFunction(petsclib::$UnionPetscLib, snes::AbstractSNES, x::AbstractPetscVec, f::AbstractPetscVec, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESPicardComputeMFFunction, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec, Ptr{Cvoid}),
               snes, x, f, ctx,
              )


	return nothing
end 

"""
	SNESPruneJacobianColor(petsclib::PetscLibType, snes::AbstractSNES, J::AbstractPetscMat, B::AbstractPetscMat) 
Remove nondiagonal zeros in the Jacobian matrix and update the `MatMFFD` coloring information based on the new nonzero structure

Collective

Input Parameters:
- `snes` - the `SNES` context
- `J`    - Jacobian matrix (not altered in this routine)
- `B`    - newly computed Jacobian matrix to use with preconditioner (generally the same as `J`)

Level: intermediate

See also: `SNESComputeJacobianDefaultColor()`, `MatEliminateZeros()`, `MatFDColoringCreate()`, `MatFDColoringSetFunction()`

# External Links
$(_doc_external("SNES/SNESPruneJacobianColor"))
"""
function SNESPruneJacobianColor(petsclib::PetscLibType, snes::AbstractSNES, J::AbstractPetscMat, B::AbstractPetscMat)
    error("SNESPruneJacobianColor: no generated method for these argument types")
end

@for_petsc function SNESPruneJacobianColor(petsclib::$UnionPetscLib, snes::AbstractSNES, J::AbstractPetscMat, B::AbstractPetscMat )

    @chk ccall(
               (:SNESPruneJacobianColor, $petsc_library),
               PetscErrorCode,
               (CSNES, CMat, CMat),
               snes, J, B,
              )


	return nothing
end 

"""
	pyname::Ptr{Cchar} = SNESPythonGetType(petsclib::PetscLibType, snes::AbstractSNES) 
Get the type of a `SNES` object implemented in Python set with `SNESPythonSetType()`

Not Collective

Input Parameter:
- `snes`  - the nonlinear solver (`SNES`) context.

Output Parameter:
- `pyname`  - full dotted Python name [package].module[.{class|function}]

Level: intermediate

See also: `SNES`, `SNESCreate()`, `SNESSetType()`, `SNESPYTHON`, `PetscPythonInitialize()`, `SNESPythonSetType()`

# External Links
$(_doc_external("SNES/SNESPythonGetType"))
"""
function SNESPythonGetType(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESPythonGetType: no generated method for these argument types")
end

@for_petsc function SNESPythonGetType(petsclib::$UnionPetscLib, snes::AbstractSNES )
	pyname_ = Ref{Ptr{Cchar}}()

    @chk ccall(
               (:SNESPythonGetType, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Ptr{Cchar}}),
               snes, pyname_,
              )

	pyname = pyname_[]

	return pyname
end 

"""
	SNESPythonSetType(petsclib::PetscLibType, snes::AbstractSNES, pyname::String) 
Initialize a `SNES` object implemented in Python.

Collective

Input Parameters:
- `snes`  - the nonlinear solver (`SNES`) context.
- `pyname`  - full dotted Python name [package].module[.{class|function}]

Options Database Key:
- `-snes_python_type pyname`  - python class

Level: intermediate

See also: `SNES`, `SNESCreate()`, `SNESSetType()`, `SNESPYTHON`, `PetscPythonInitialize()`, `SNESPythonGetType()`

# External Links
$(_doc_external("SNES/SNESPythonSetType"))
"""
function SNESPythonSetType(petsclib::PetscLibType, snes::AbstractSNES, pyname::String)
    error("SNESPythonSetType: no generated method for these argument types")
end

@for_petsc function SNESPythonSetType(petsclib::$UnionPetscLib, snes::AbstractSNES, pyname::String )

    @chk ccall(
               (:SNESPythonSetType, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Cchar}),
               snes, pyname,
              )


	return nothing
end 

"""
	SNESQNSetRestartType(petsclib::PetscLibType, snes::AbstractSNES, rtype::SNESQNRestartType) 
Sets the restart type for `SNESQN`.

Logically Collective

Input Parameters:
- `snes`  - the iterative context
- `rtype` - restart type, see `SNESQNRestartType`

Options Database Keys:
- `-snes_qn_restart_type (powell|periodic|none)` - set the restart type
- `-snes_qn_m m`                                 - sets the number of stored updates and the restart period for periodic

Level: intermediate

See also: `SNES`, `SNESQN`, `SNESQNRestartType`, `SNES_QN_RESTART_NONE`, `SNES_QN_RESTART_POWELL`, `SNES_QN_RESTART_PERIODIC`,
`SNESQNType`, `SNESQNScaleType`

# External Links
$(_doc_external("SNES/SNESQNSetRestartType"))
"""
function SNESQNSetRestartType(petsclib::PetscLibType, snes::AbstractSNES, rtype::SNESQNRestartType)
    error("SNESQNSetRestartType: no generated method for these argument types")
end

@for_petsc function SNESQNSetRestartType(petsclib::$UnionPetscLib, snes::AbstractSNES, rtype::SNESQNRestartType )

    @chk ccall(
               (:SNESQNSetRestartType, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESQNRestartType),
               snes, rtype,
              )


	return nothing
end 

"""
	SNESQNSetScaleType(petsclib::PetscLibType, snes::AbstractSNES, stype::SNESQNScaleType) 
Sets the scaling type for the inner inverse Jacobian in `SNESQN`.

Logically Collective

Input Parameters:
- `snes`  - the nonlinear solver context
- `stype` - scale type, see `SNESQNScaleType`

Options Database Key:
- `-snes_qn_scale_type (diagonal|none|scalar|jacobian)` - Scaling type

Level: intermediate

See also: `SNES`, `SNESQN`, `SNESLineSearch`, `SNESQNScaleType`, `SNESSetJacobian()`, `SNESQNType`, `SNESQNRestartType`

# External Links
$(_doc_external("SNES/SNESQNSetScaleType"))
"""
function SNESQNSetScaleType(petsclib::PetscLibType, snes::AbstractSNES, stype::SNESQNScaleType)
    error("SNESQNSetScaleType: no generated method for these argument types")
end

@for_petsc function SNESQNSetScaleType(petsclib::$UnionPetscLib, snes::AbstractSNES, stype::SNESQNScaleType )

    @chk ccall(
               (:SNESQNSetScaleType, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESQNScaleType),
               snes, stype,
              )


	return nothing
end 

"""
	SNESQNSetType(petsclib::PetscLibType, snes::AbstractSNES, qtype::SNESQNType) 
Sets the quasi-Newton variant to be used in `SNESQN`.

Logically Collective

Input Parameters:
- `snes`  - the iterative context
- `qtype` - variant type, see `SNESQNType`

Options Database Key:
- `-snes_qn_type (lbfgs|broyden|badbroyden)` - quasi-Newton type

Level: intermediate

See also: `SNESQN`, `SNES_QN_LBFGS`, `SNES_QN_BROYDEN`, `SNES_QN_BADBROYDEN`, `SNESQNType`, `SNESQNScaleType`, `TAOLMVM`, `TAOBLMVM`

# External Links
$(_doc_external("SNES/SNESQNSetType"))
"""
function SNESQNSetType(petsclib::PetscLibType, snes::AbstractSNES, qtype::SNESQNType)
    error("SNESQNSetType: no generated method for these argument types")
end

@for_petsc function SNESQNSetType(petsclib::$UnionPetscLib, snes::AbstractSNES, qtype::SNESQNType )

    @chk ccall(
               (:SNESQNSetType, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESQNType),
               snes, qtype,
              )


	return nothing
end 

"""
	SNESRegister(petsclib::PetscLibType, sname::String, fnc::external) 
Adds a method to the nonlinear solver package.

Not Collective

Input Parameters:
- `sname`    - name of a new user-defined solver
- `function` - routine to create method context

Level: advanced

See also: `SNESRegisterAll()`, `SNESRegisterDestroy()`

# External Links
$(_doc_external("SNES/SNESRegister"))
"""
function SNESRegister(petsclib::PetscLibType, sname::String, fnc::external)
    error("SNESRegister: no generated method for these argument types")
end

@for_petsc function SNESRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:SNESRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	SNESReset(petsclib::PetscLibType, snes::AbstractSNES) 
Resets a `SNES` context to the state it was in before `SNESSetUp()` was called and removes any allocated `Vec` and `Mat` from its data structures

Collective

Input Parameter:
- `snes` - the nonlinear iterative solver context obtained from `SNESCreate()`

Level: intermediate

See also: `SNES`, `SNESDestroy()`, `SNESCreate()`, `SNESSetUp()`, `SNESSolve()`

# External Links
$(_doc_external("SNES/SNESReset"))
"""
function SNESReset(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESReset: no generated method for these argument types")
end

@for_petsc function SNESReset(petsclib::$UnionPetscLib, snes::AbstractSNES )

    @chk ccall(
               (:SNESReset, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	SNESResetCounters(petsclib::PetscLibType, snes::AbstractSNES) 
Reset counters for linear iterations and function evaluations.

Logically Collective

Input Parameters:
- `snes` - `SNES` context

Level: developer

See also: `SNESGetNumberFunctionEvals()`, `SNESGetLinearSolveIterations()`, `SNESGetNPC()`

# External Links
$(_doc_external("SNES/SNESResetCounters"))
"""
function SNESResetCounters(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESResetCounters: no generated method for these argument types")
end

@for_petsc function SNESResetCounters(petsclib::$UnionPetscLib, snes::AbstractSNES )

    @chk ccall(
               (:SNESResetCounters, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	SNESResetFromOptions(petsclib::PetscLibType, snes::AbstractSNES) 
Sets various `SNES` and `KSP` parameters from user options ONLY if the `SNESSetFromOptions()` was previously called

Collective

Input Parameter:
- `snes` - the `SNES` context

Level: advanced

See also: `SNES`, `SNESSetFromOptions()`, `SNESSetOptionsPrefix()`

# External Links
$(_doc_external("SNES/SNESResetFromOptions"))
"""
function SNESResetFromOptions(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESResetFromOptions: no generated method for these argument types")
end

@for_petsc function SNESResetFromOptions(petsclib::$UnionPetscLib, snes::AbstractSNES )

    @chk ccall(
               (:SNESResetFromOptions, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	SNESSetAlwaysComputesFinalResidual(petsclib::PetscLibType, snes::AbstractSNES, flg::PetscBool) 
tells the `SNES` to always compute the residual (nonlinear function value) at the final solution

Logically Collective

Input Parameters:
- `snes` - the shell `SNES`
- `flg`  - `PETSC_TRUE` to always compute the residual

Level: advanced

See also: `SNES`, `SNESFAS`, `SNESSolve()`, `SNESGetAlwaysComputesFinalResidual()`

# External Links
$(_doc_external("SNES/SNESSetAlwaysComputesFinalResidual"))
"""
function SNESSetAlwaysComputesFinalResidual(petsclib::PetscLibType, snes::AbstractSNES, flg::PetscBool)
    error("SNESSetAlwaysComputesFinalResidual: no generated method for these argument types")
end

@for_petsc function SNESSetAlwaysComputesFinalResidual(petsclib::$UnionPetscLib, snes::AbstractSNES, flg::PetscBool )

    @chk ccall(
               (:SNESSetAlwaysComputesFinalResidual, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, flg,
              )


	return nothing
end 

"""
	SNESSetApplicationContext(petsclib::PetscLibType, snes::AbstractSNES, ctx::Ptr{Cvoid}) 
Sets the optional user-defined context for the nonlinear solvers.

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `ctx`  - the application context

Level: intermediate

See also: `SNES`, `SNESSetComputeApplicationContext()`, `SNESGetApplicationContext()`

# External Links
$(_doc_external("SNES/SNESSetApplicationContext"))
"""
function SNESSetApplicationContext(petsclib::PetscLibType, snes::AbstractSNES, ctx::Ptr{Cvoid})
    error("SNESSetApplicationContext: no generated method for these argument types")
end

@for_petsc function SNESSetApplicationContext(petsclib::$UnionPetscLib, snes::AbstractSNES, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESSetApplicationContext, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Cvoid}),
               snes, ctx,
              )


	return nothing
end 

"""
	SNESSetCheckJacobianDomainError(petsclib::PetscLibType, snes::AbstractSNES, flg::PetscBool) 
tells `SNESSolve()` whether to check if the user called `SNESSetJacobianDomainError()` to indicate a Jacobian domain error after
each Jacobian evaluation.

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `flg`  - indicates if or not to check Jacobian domain error after each Jacobian evaluation

Level: advanced

See also: `SNES`, `SNESConvergedReason`, `SNESCreate()`, `SNESSetFunction()`, `SNESFunctionFn`, `SNESSetFunctionDomainError()`, `SNESGetCheckJacobianDomainError()`

# External Links
$(_doc_external("SNES/SNESSetCheckJacobianDomainError"))
"""
function SNESSetCheckJacobianDomainError(petsclib::PetscLibType, snes::AbstractSNES, flg::PetscBool)
    error("SNESSetCheckJacobianDomainError: no generated method for these argument types")
end

@for_petsc function SNESSetCheckJacobianDomainError(petsclib::$UnionPetscLib, snes::AbstractSNES, flg::PetscBool )

    @chk ccall(
               (:SNESSetCheckJacobianDomainError, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, flg,
              )


	return nothing
end 

"""
	SNESSetComputeApplicationContext(petsclib::PetscLibType, snes::AbstractSNES, compute::external, destroy::Ptr{Cvoid}) 
Sets an optional function to compute a user-defined context for
the nonlinear solvers.

Logically Collective; No Fortran Support

Input Parameters:
- `snes`    - the `SNES` context
- `compute` - function to compute the context
- `destroy` - function to destroy the context, see `PetscCtxDestroyFn` for the calling sequence

Calling sequence of `compute`:
- `snes` - the `SNES` context
- `ctx`  - context to be computed

Level: intermediate

See also: `SNESGetApplicationContext()`, `SNESSetApplicationContext()`, `PetscCtxDestroyFn`

# External Links
$(_doc_external("SNES/SNESSetComputeApplicationContext"))
"""
function SNESSetComputeApplicationContext(petsclib::PetscLibType, snes::AbstractSNES, compute::external, destroy::Ptr{Cvoid})
    error("SNESSetComputeApplicationContext: no generated method for these argument types")
end

@for_petsc function SNESSetComputeApplicationContext(petsclib::$UnionPetscLib, snes::AbstractSNES, compute::external, destroy::Ptr{Cvoid} )

    @chk ccall(
               (:SNESSetComputeApplicationContext, $petsc_library),
               PetscErrorCode,
               (CSNES, external, Ptr{Cvoid}),
               snes, compute, destroy,
              )


	return nothing
end 

"""
	SNESSetComputeInitialGuess(petsclib::PetscLibType, snes::AbstractSNES, func::Ptr{Cvoid}, ctx::Ptr{Cvoid}) 
Sets a routine used to compute an initial guess for the nonlinear problem

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `func` - function evaluation routine, see `SNESInitialGuessFn` for the calling sequence
- `ctx`  - [optional] user-defined context for private data for the
function evaluation routine (may be `NULL`)

Level: intermediate

See also: `SNES`, `SNESSolve()`, `SNESSetFunction()`, `SNESGetFunction()`, `SNESComputeFunction()`, `SNESSetJacobian()`, `SNESInitialGuessFn`

# External Links
$(_doc_external("SNES/SNESSetComputeInitialGuess"))
"""
function SNESSetComputeInitialGuess(petsclib::PetscLibType, snes::AbstractSNES, func::Ptr{Cvoid}, ctx::Ptr{Cvoid})
    error("SNESSetComputeInitialGuess: no generated method for these argument types")
end

@for_petsc function SNESSetComputeInitialGuess(petsclib::$UnionPetscLib, snes::AbstractSNES, func::Ptr{Cvoid}, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESSetComputeInitialGuess, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Cvoid}, Ptr{Cvoid}),
               snes, func, ctx,
              )


	return nothing
end 

"""
	SNESSetConvergedReason(petsclib::PetscLibType, snes::AbstractSNES, reason::SNESConvergedReason) 
Sets the reason the `SNES` iteration was stopped.

Not Collective

Input Parameters:
- `snes`   - the `SNES` context
- `reason` - negative value indicates diverged, positive value converged, see `SNESConvergedReason` or the
manual pages for the individual convergence tests for complete lists

Level: developer

See also: `SNESGetConvergedReason()`, `SNESSetConvergenceTest()`, `SNESConvergedReason`

# External Links
$(_doc_external("SNES/SNESSetConvergedReason"))
"""
function SNESSetConvergedReason(petsclib::PetscLibType, snes::AbstractSNES, reason::SNESConvergedReason)
    error("SNESSetConvergedReason: no generated method for these argument types")
end

@for_petsc function SNESSetConvergedReason(petsclib::$UnionPetscLib, snes::AbstractSNES, reason::SNESConvergedReason )

    @chk ccall(
               (:SNESSetConvergedReason, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESConvergedReason),
               snes, reason,
              )


	return nothing
end 

"""
	SNESSetConvergenceHistory(petsclib::PetscLibType, snes::AbstractSNES, a::Vector{PetscReal}, its::Vector{PetscInt}, na::PetscInt, reset::PetscBool) 
Sets the arrays used to hold the convergence history.

Logically Collective

Input Parameters:
- `snes`  - iterative context obtained from `SNESCreate()`
- `a`     - array to hold history, this array will contain the function norms computed at each step
- `its`   - integer array holds the number of linear iterations for each solve.
- `na`    - size of `a` and `its`
- `reset` - `PETSC_TRUE` indicates each new nonlinear solve resets the history counter to zero,
else it continues storing new values for new nonlinear solves after the old ones

Level: intermediate

See also: `SNES`, `SNESSolve()`, `SNESGetConvergenceHistory()`

# External Links
$(_doc_external("SNES/SNESSetConvergenceHistory"))
"""
function SNESSetConvergenceHistory(petsclib::PetscLibType, snes::AbstractSNES, a::AbstractVector{<:Number}, its::AbstractVector{<:Number}, na::Integer, reset::PetscBool)
    error("SNESSetConvergenceHistory: no generated method for these argument types")
end

@for_petsc function SNESSetConvergenceHistory(petsclib::$UnionPetscLib, snes::AbstractSNES, a::Vector{$PetscReal}, its::Vector{$PetscInt}, na::$PetscInt, reset::PetscBool )

    @chk ccall(
               (:SNESSetConvergenceHistory, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscReal}, Ptr{$PetscInt}, $PetscInt, PetscBool),
               snes, a, its, na, reset,
              )


	return nothing
end 

"""
	SNESSetConvergenceTest(petsclib::PetscLibType, snes::AbstractSNES, func::external, ctx::Ptr{Cvoid}, destroy::Ptr{Cvoid}) 
Sets the function that is to be used
to test for convergence of the nonlinear iterative solution.

Logically Collective

Input Parameters:
- `snes`    - the `SNES` context
- `func`    - routine to test for convergence
- `ctx`     - [optional] context for private data for the convergence routine  (may be `NULL`)
- `destroy` - [optional] destructor for the context (may be `NULL`; `PETSC_NULL_FUNCTION` in Fortran)

Calling sequence of func:
- `snes`   - the `SNES` context
- `it`     - the current iteration number
- `xnorm`  - the norm of the new solution
- `snorm`  - the norm of the step
- `fnorm`  - the norm of the function value
- `reason` - output, the reason convergence or divergence as declared
- `ctx`    - the optional convergence test context

Level: advanced

See also: `SNES`, `SNESConvergedDefault()`, `SNESConvergedSkip()`

# External Links
$(_doc_external("SNES/SNESSetConvergenceTest"))
"""
function SNESSetConvergenceTest(petsclib::PetscLibType, snes::AbstractSNES, func::external, ctx::Ptr{Cvoid}, destroy::Ptr{Cvoid})
    error("SNESSetConvergenceTest: no generated method for these argument types")
end

@for_petsc function SNESSetConvergenceTest(petsclib::$UnionPetscLib, snes::AbstractSNES, func::external, ctx::Ptr{Cvoid}, destroy::Ptr{Cvoid} )

    @chk ccall(
               (:SNESSetConvergenceTest, $petsc_library),
               PetscErrorCode,
               (CSNES, external, Ptr{Cvoid}, Ptr{Cvoid}),
               snes, func, ctx, destroy,
              )


	return nothing
end 

"""
	SNESSetCountersReset(petsclib::PetscLibType, snes::AbstractSNES, reset::PetscBool) 
Sets whether or not the counters for linear iterations and function evaluations
are reset every time `SNESSolve()` is called.

Logically Collective

Input Parameters:
- `snes`  - `SNES` context
- `reset` - whether to reset the counters or not, defaults to `PETSC_TRUE`

Level: developer

See also: `SNESGetNumberFunctionEvals()`, `SNESGetLinearSolveIterations()`, `SNESGetNPC()`

# External Links
$(_doc_external("SNES/SNESSetCountersReset"))
"""
function SNESSetCountersReset(petsclib::PetscLibType, snes::AbstractSNES, reset::PetscBool)
    error("SNESSetCountersReset: no generated method for these argument types")
end

@for_petsc function SNESSetCountersReset(petsclib::$UnionPetscLib, snes::AbstractSNES, reset::PetscBool )

    @chk ccall(
               (:SNESSetCountersReset, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, reset,
              )


	return nothing
end 

"""
	SNESSetDM(petsclib::PetscLibType, snes::AbstractSNES, dm::AbstractPetscDM) 
Sets the `DM` that may be used by some `SNES` nonlinear solvers or their underlying preconditioners

Logically Collective

Input Parameters:
- `snes` - the nonlinear solver context
- `dm`   - the `DM`, cannot be `NULL`

Level: intermediate

See also: `DM`, `SNES`, `SNESGetDM()`, `KSPSetDM()`, `KSPGetDM()`

# External Links
$(_doc_external("SNES/SNESSetDM"))
"""
function SNESSetDM(petsclib::PetscLibType, snes::AbstractSNES, dm::AbstractPetscDM)
    error("SNESSetDM: no generated method for these argument types")
end

@for_petsc function SNESSetDM(petsclib::$UnionPetscLib, snes::AbstractSNES, dm::AbstractPetscDM )

    @chk ccall(
               (:SNESSetDM, $petsc_library),
               PetscErrorCode,
               (CSNES, CDM),
               snes, dm,
              )


	return nothing
end 

"""
	SNESSetDivergenceTolerance(petsclib::PetscLibType, snes::AbstractSNES, divtol::PetscReal) 
Sets the divergence tolerance used for the `SNES` divergence test.

Logically Collective

Input Parameters:
- `snes`   - the `SNES` context
- `divtol` - the divergence tolerance. Use `PETSC_UNLIMITED` to deactivate the test. If the residual norm  F(x^n) \\ge divtol * F(x^0)  the solver
is stopped due to divergence.

Options Database Key:
- `-snes_divergence_tolerance divtol` - Sets `divtol`

Level: intermediate

See also: `SNES`, `SNESSolve()`, `SNESSetTolerances()`, `SNESGetDivergenceTolerance()`

# External Links
$(_doc_external("SNES/SNESSetDivergenceTolerance"))
"""
function SNESSetDivergenceTolerance(petsclib::PetscLibType, snes::AbstractSNES, divtol::Real)
    error("SNESSetDivergenceTolerance: no generated method for these argument types")
end

@for_petsc function SNESSetDivergenceTolerance(petsclib::$UnionPetscLib, snes::AbstractSNES, divtol::$PetscReal )

    @chk ccall(
               (:SNESSetDivergenceTolerance, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscReal),
               snes, divtol,
              )


	return nothing
end 

"""
	SNESSetErrorIfNotConverged(petsclib::PetscLibType, snes::AbstractSNES, flg::PetscBool) 
Causes `SNESSolve()` to generate an error immediately if the solver has not converged.

Logically Collective

Input Parameters:
- `snes` - iterative context obtained from `SNESCreate()`
- `flg`  - `PETSC_TRUE` indicates you want the error generated

Options Database Key:
- `-snes_error_if_not_converged (true|false)` - cause an immediate error condition and stop the program if the solver does not converge

Level: intermediate

See also: `SNES`, `SNESGetErrorIfNotConverged()`, `KSPGetErrorIfNotConverged()`, `KSPSetErrorIfNotConverged()`

# External Links
$(_doc_external("SNES/SNESSetErrorIfNotConverged"))
"""
function SNESSetErrorIfNotConverged(petsclib::PetscLibType, snes::AbstractSNES, flg::PetscBool)
    error("SNESSetErrorIfNotConverged: no generated method for these argument types")
end

@for_petsc function SNESSetErrorIfNotConverged(petsclib::$UnionPetscLib, snes::AbstractSNES, flg::PetscBool )

    @chk ccall(
               (:SNESSetErrorIfNotConverged, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, flg,
              )


	return nothing
end 

"""
	SNESSetForceIteration(petsclib::PetscLibType, snes::AbstractSNES, force::PetscBool) 
force `SNESSolve()` to take at least one iteration regardless of the initial residual norm

Logically Collective

Input Parameters:
- `snes`  - the `SNES` context
- `force` - `PETSC_TRUE` require at least one iteration

Options Database Key:
- `-snes_force_iteration force` - Sets forcing an iteration

Level: intermediate

See also: `SNES`, `TS`, `SNESSetDivergenceTolerance()`

# External Links
$(_doc_external("SNES/SNESSetForceIteration"))
"""
function SNESSetForceIteration(petsclib::PetscLibType, snes::AbstractSNES, force::PetscBool)
    error("SNESSetForceIteration: no generated method for these argument types")
end

@for_petsc function SNESSetForceIteration(petsclib::$UnionPetscLib, snes::AbstractSNES, force::PetscBool )

    @chk ccall(
               (:SNESSetForceIteration, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, force,
              )


	return nothing
end 

"""
	SNESSetFromOptions(petsclib::PetscLibType, snes::AbstractSNES) 
Sets various `SNES` and `KSP` parameters from user options.

Collective

Input Parameter:
- `snes` - the `SNES` context

Options Database Keys:
- `-snes_type type`                                                              - newtonls, newtontr, ngmres, ncg, nrichardson, qn, vi, fas, `SNESType` for complete list
- `-snes_rtol rtol`                                                              - relative decrease in tolerance norm from initial
- `-snes_atol abstol`                                                            - absolute tolerance of residual norm
- `-snes_stol stol`                                                              - convergence tolerance in terms of the norm of the change in the solution between steps
- `-snes_divergence_tolerance divtol`                                            - if the residual goes above divtol*rnorm0, exit with divergence
- `-snes_max_it max_it`                                                          - maximum number of iterations
- `-snes_max_funcs max_funcs`                                                    - maximum number of function evaluations
- `-snes_force_iteration force`                                                  - force `SNESSolve()` to take at least one iteration
- `-snes_max_fail max_fail`                                                      - maximum number of line search failures allowed before stopping, default is none
- `-snes_max_linear_solve_fail`                                                  - number of linear solver failures before SNESSolve() stops
- `-snes_lag_preconditioner lag`                                                 - how often preconditioner is rebuilt (use -1 to never rebuild)
- `-snes_lag_preconditioner_persists (true|false)`                               - retains the -snes_lag_preconditioner information across multiple SNESSolve()
- `-snes_lag_jacobian lag`                                                       - how often Jacobian is rebuilt (use -1 to never rebuild)
- `-snes_lag_jacobian_persists (true|false)`                                     - retains the -snes_lag_jacobian information across multiple SNESSolve()
- `-snes_convergence_test (default|skip|correct_pressure)`                       - convergence test in nonlinear solver. default `SNESConvergedDefault()`. skip `SNESConvergedSkip()` means continue
iterating until max_it or some other criterion is reached, saving expense of convergence test. correct_pressure
`SNESConvergedCorrectPressure()` has special handling of a pressure null space.
- `-snes_monitor [ascii][:filename][:viewer format]`                             - prints residual norm at each iteration. if no filename given prints to stdout
- `-snes_monitor_solution [ascii binary draw][:filename][:viewer format]`        - plots solution at each iteration
- `-snes_monitor_residual [ascii binary draw][:filename][:viewer format]`        - plots residual (not its norm) at each iteration
- `-snes_monitor_solution_update [ascii binary draw][:filename][:viewer format]` - plots update to solution at each iteration
- `-snes_monitor draw::draw_lg`                                                  - plots residual norm at each iteration
- `-snes_monitor_lg_range`                                                       - plots function range at each iteration
- `-snes_monitor_pause_final`                                                    - Pauses all monitor drawing after the solver ends
- `-snes_fd`                                                                     - use finite differences to compute Jacobian; very slow, only for testing
- `-snes_fd_color`                                                               - use finite differences with coloring to compute Jacobian
- `-snes_mf_ksp_monitor`                                                         - if using matrix-free multiply then print h at each `KSP` iteration
- `-snes_converged_reason`                                                       - print the reason for convergence/divergence after each solve
- `-npc_snes_type type`                                                          - the `SNES` type to use as a nonlinear preconditioner
- `-snes_test_jacobian [threshold]`                                              - compare the user provided Jacobian with one computed via finite differences to check for errors.
If a threshold is given, display only those entries whose difference is greater than the threshold.
- `-snes_test_jacobian_view`                                                     - display the user provided Jacobian, the finite difference Jacobian and the difference between them
to help users detect the location of errors in the user provided Jacobian.

Options Database Keys for Eisenstat-Walker method:
- `-snes_ksp_ew`                     - use Eisenstat-Walker method for determining linear system convergence
- `-snes_ksp_ew_version ver`         - version of  Eisenstat-Walker method
- `-snes_ksp_ew_rtol0 rtol0`         - Sets rtol0
- `-snes_ksp_ew_rtolmax rtolmax`     - Sets rtolmax
- `-snes_ksp_ew_gamma gamma`         - Sets gamma
- `-snes_ksp_ew_alpha alpha`         - Sets alpha
- `-snes_ksp_ew_alpha2 alpha2`       - Sets alpha2
- `-snes_ksp_ew_threshold threshold` - Sets threshold

Level: beginner

See also: `SNESType`, `SNESSetOptionsPrefix()`, `SNESResetFromOptions()`, `SNES`, `SNESCreate()`, `MatCreateSNESMF()`, `MatFDColoring`

# External Links
$(_doc_external("SNES/SNESSetFromOptions"))
"""
function SNESSetFromOptions(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESSetFromOptions: no generated method for these argument types")
end

@for_petsc function SNESSetFromOptions(petsclib::$UnionPetscLib, snes::AbstractSNES )

    @chk ccall(
               (:SNESSetFromOptions, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	SNESSetFunction(petsclib::PetscLibType, snes::AbstractSNES, r::AbstractPetscVec, f::Ptr{Cvoid}, ctx::Ptr{Cvoid}) 
Sets the function evaluation routine and function
vector for use by the `SNES` routines in solving systems of nonlinear
equations.

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `r`    - vector to store function values, may be `NULL`
- `f`    - function evaluation routine;  for calling sequence see `SNESFunctionFn`
- `ctx`  - [optional] user-defined context for private data for the
function evaluation routine (may be `NULL`)

Level: beginner

See also: `SNES`, `SNESGetFunction()`, `SNESComputeFunction()`, `SNESSetJacobian()`, `SNESSetPicard()`, `SNESFunctionFn`

# External Links
$(_doc_external("SNES/SNESSetFunction"))
"""
function SNESSetFunction(petsclib::PetscLibType, snes::AbstractSNES, r::AbstractPetscVec, f::Ptr{Cvoid}, ctx::Ptr{Cvoid})
    error("SNESSetFunction: no generated method for these argument types")
end

@for_petsc function SNESSetFunction(petsclib::$UnionPetscLib, snes::AbstractSNES, r::AbstractPetscVec, f::Ptr{Cvoid}, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESSetFunction, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, Ptr{Cvoid}, Ptr{Cvoid}),
               snes, r, f, ctx,
              )


	return nothing
end 

"""
	SNESSetFunctionDomainError(petsclib::PetscLibType, snes::AbstractSNES) 
tells `SNES` that the input vector, a proposed new solution, to your function you provided to `SNESSetFunction()` is not
in the function's domain. For example, a step with negative pressure.

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Level: advanced

See also: `SNESCreate()`, `SNESSetFunction()`, `SNESFunctionFn`, `SNESSetJacobianDomainError()`, `SNESVISetVariableBounds()`,
`SNESVISetComputeVariableBounds()`, `SNESLineSearchSetPreCheck()`, `SNESLineSearchSetPostCheck()`, `SNESConvergedReason`, `SNESGetConvergedReason()`,
`SNES_DIVERGED_FUNCTION_DOMAIN`, `SNESSetObjectiveDomainError()`, `SNES_DIVERGED_OBJECTIVE_DOMAIN`

# External Links
$(_doc_external("SNES/SNESSetFunctionDomainError"))
"""
function SNESSetFunctionDomainError(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESSetFunctionDomainError: no generated method for these argument types")
end

@for_petsc function SNESSetFunctionDomainError(petsclib::$UnionPetscLib, snes::AbstractSNES )

    @chk ccall(
               (:SNESSetFunctionDomainError, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	SNESSetFunctionNorm(petsclib::PetscLibType, snes::AbstractSNES, norm::PetscReal) 
Sets the last computed residual norm.

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `norm` - the value of the norm

Level: developer

See also: `SNES`, `SNESGetNormSchedule()`, `SNESComputeFunction()`, `VecNorm()`, `SNESSetFunction()`, `SNESSetInitialFunction()`, `SNESNormSchedule`

# External Links
$(_doc_external("SNES/SNESSetFunctionNorm"))
"""
function SNESSetFunctionNorm(petsclib::PetscLibType, snes::AbstractSNES, norm::Real)
    error("SNESSetFunctionNorm: no generated method for these argument types")
end

@for_petsc function SNESSetFunctionNorm(petsclib::$UnionPetscLib, snes::AbstractSNES, norm::$PetscReal )

    @chk ccall(
               (:SNESSetFunctionNorm, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscReal),
               snes, norm,
              )


	return nothing
end 

"""
	SNESSetFunctionType(petsclib::PetscLibType, snes::AbstractSNES, type::SNESFunctionType) 
Sets the `SNESFunctionType`
of the `SNES` method.

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `type` - the function type

Level: developer

Values of the function type:
- `SNES_FUNCTION_DEFAULT`           - the default for the given `SNESType`
- `SNES_FUNCTION_UNPRECONDITIONED`  - an unpreconditioned function evaluation (this is the function provided with `SNESSetFunction()`
- `SNES_FUNCTION_PRECONDITIONED`    - a transformation of the function provided with `SNESSetFunction()`

See also: `SNES`, `SNESFunctionType`, `SNESGetNormSchedule()`, `SNESComputeFunction()`, `VecNorm()`, `SNESSetFunction()`, `SNESSetInitialFunction()`, `SNESNormSchedule`

# External Links
$(_doc_external("SNES/SNESSetFunctionType"))
"""
function SNESSetFunctionType(petsclib::PetscLibType, snes::AbstractSNES, type::SNESFunctionType)
    error("SNESSetFunctionType: no generated method for these argument types")
end

@for_petsc function SNESSetFunctionType(petsclib::$UnionPetscLib, snes::AbstractSNES, type::SNESFunctionType )

    @chk ccall(
               (:SNESSetFunctionType, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESFunctionType),
               snes, type,
              )


	return nothing
end 

"""
	SNESSetGridSequence(petsclib::PetscLibType, snes::AbstractSNES, steps::PetscInt) 
sets the number of steps of grid sequencing that `SNES` will do

Logically Collective

Input Parameters:
- `snes`  - the `SNES` context
- `steps` - the number of refinements to do, defaults to 0

Options Database Key:
- `-snes_grid_sequence steps` - Use grid sequencing to generate initial guess

Level: intermediate

See also: `SNES`, `SNESGetLagPreconditioner()`, `SNESSetLagJacobian()`, `SNESGetLagJacobian()`, `SNESGetGridSequence()`,
`SNESSetDM()`, `SNESSolve()`

# External Links
$(_doc_external("SNES/SNESSetGridSequence"))
"""
function SNESSetGridSequence(petsclib::PetscLibType, snes::AbstractSNES, steps::Integer)
    error("SNESSetGridSequence: no generated method for these argument types")
end

@for_petsc function SNESSetGridSequence(petsclib::$UnionPetscLib, snes::AbstractSNES, steps::$PetscInt )

    @chk ccall(
               (:SNESSetGridSequence, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt),
               snes, steps,
              )


	return nothing
end 

"""
	SNESSetInitialFunction(petsclib::PetscLibType, snes::AbstractSNES, f::AbstractPetscVec) 
Set an already computed function evaluation at the initial guess to be reused by `SNESSolve()`.

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `f`    - vector to store function value

Level: developer

See also: `SNES`, `SNESFAS`, `SNESSetFunction()`, `SNESComputeFunction()`, `SNESSetInitialFunctionNorm()`

# External Links
$(_doc_external("SNES/SNESSetInitialFunction"))
"""
function SNESSetInitialFunction(petsclib::PetscLibType, snes::AbstractSNES, f::AbstractPetscVec)
    error("SNESSetInitialFunction: no generated method for these argument types")
end

@for_petsc function SNESSetInitialFunction(petsclib::$UnionPetscLib, snes::AbstractSNES, f::AbstractPetscVec )

    @chk ccall(
               (:SNESSetInitialFunction, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec),
               snes, f,
              )


	return nothing
end 

"""
	SNESSetIterationNumber(petsclib::PetscLibType, snes::AbstractSNES, iter::PetscInt) 
Sets the current iteration number.

Not Collective

Input Parameters:
- `snes` - `SNES` context
- `iter` - iteration number

Level: developer

See also: `SNESGetLinearSolveIterations()`

# External Links
$(_doc_external("SNES/SNESSetIterationNumber"))
"""
function SNESSetIterationNumber(petsclib::PetscLibType, snes::AbstractSNES, iter::Integer)
    error("SNESSetIterationNumber: no generated method for these argument types")
end

@for_petsc function SNESSetIterationNumber(petsclib::$UnionPetscLib, snes::AbstractSNES, iter::$PetscInt )

    @chk ccall(
               (:SNESSetIterationNumber, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt),
               snes, iter,
              )


	return nothing
end 

"""
	SNESSetJacobian(petsclib::PetscLibType, snes::AbstractSNES, Amat::AbstractPetscMat, Pmat::AbstractPetscMat, J::Ptr{Cvoid}, ctx::Ptr{Cvoid}) 
Sets the function to compute Jacobian as well as the
location to store the matrix.

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `Amat` - the matrix that defines the (approximate) Jacobian
- `Pmat` - the matrix to be used in constructing the preconditioner, usually the same as `Amat`.
- `J`    - Jacobian evaluation routine (if `NULL` then `SNES` retains any previously set value), see `SNESJacobianFn` for details
- `ctx`  - [optional] user-defined context for private data for the
Jacobian evaluation routine (may be `NULL`) (if `NULL` then `SNES` retains any previously set value)

Level: beginner

See also: `SNES`, `KSPSetOperators()`, `SNESSetFunction()`, `MatMFFDComputeJacobian()`, `SNESComputeJacobianDefaultColor()`, `MatStructure`,
`SNESSetPicard()`, `SNESJacobianFn`, `SNESFunctionFn`

# External Links
$(_doc_external("SNES/SNESSetJacobian"))
"""
function SNESSetJacobian(petsclib::PetscLibType, snes::AbstractSNES, Amat::AbstractPetscMat, Pmat::AbstractPetscMat, J::Ptr{Cvoid}, ctx::Ptr{Cvoid})
    error("SNESSetJacobian: no generated method for these argument types")
end

@for_petsc function SNESSetJacobian(petsclib::$UnionPetscLib, snes::AbstractSNES, Amat::AbstractPetscMat, Pmat::AbstractPetscMat, J::Ptr{Cvoid}, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESSetJacobian, $petsc_library),
               PetscErrorCode,
               (CSNES, CMat, CMat, Ptr{Cvoid}, Ptr{Cvoid}),
               snes, Amat, Pmat, J, ctx,
              )


	return nothing
end 

"""
	SNESSetJacobianDomainError(petsclib::PetscLibType, snes::AbstractSNES) 
tells `SNES` that the function you provided to `SNESSetJacobian()` at the proposed step. For example there is a negative element transformation.

Logically Collective

Input Parameter:
- `snes` - the `SNES` context

Level: advanced

See also: `SNESCreate()`, `SNESSetFunction()`, `SNESFunctionFn`, `SNESSetFunctionDomainError()`, `SNESVISetVariableBounds()`,
`SNESVISetComputeVariableBounds()`, `SNESLineSearchSetPreCheck()`, `SNESLineSearchSetPostCheck()`, `SNESConvergedReason`, `SNESGetConvergedReason()`

# External Links
$(_doc_external("SNES/SNESSetJacobianDomainError"))
"""
function SNESSetJacobianDomainError(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESSetJacobianDomainError: no generated method for these argument types")
end

@for_petsc function SNESSetJacobianDomainError(petsclib::$UnionPetscLib, snes::AbstractSNES )

    @chk ccall(
               (:SNESSetJacobianDomainError, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	SNESSetKSP(petsclib::PetscLibType, snes::AbstractSNES, ksp::AbstractKSP) 
Sets a `KSP` context for the `SNES` object to use

Not Collective, but the `SNES` and `KSP` objects must live on the same `MPI_Comm`

Input Parameters:
- `snes` - the `SNES` context
- `ksp`  - the `KSP` context

Level: developer

See also: `SNES`, `KSP`, `KSPGetPC()`, `SNESCreate()`, `KSPCreate()`

# External Links
$(_doc_external("SNES/SNESSetKSP"))
"""
function SNESSetKSP(petsclib::PetscLibType, snes::AbstractSNES, ksp::AbstractKSP)
    error("SNESSetKSP: no generated method for these argument types")
end

@for_petsc function SNESSetKSP(petsclib::$UnionPetscLib, snes::AbstractSNES, ksp::AbstractKSP )

    @chk ccall(
               (:SNESSetKSP, $petsc_library),
               PetscErrorCode,
               (CSNES, CKSP),
               snes, ksp,
              )


	return nothing
end 

"""
	SNESSetLagJacobian(petsclib::PetscLibType, snes::AbstractSNES, lag::PetscInt) 
Set when the Jacobian is rebuilt in the nonlinear solve. See `SNESSetLagPreconditioner()` for determining how
often the preconditioner is rebuilt.

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `lag`  - -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time
the Jacobian is built etc. -2 means rebuild at next chance but then never again

Options Database Keys:
- `-snes_lag_jacobian_persists (true|false)`       - sets the persistence through multiple SNES solves
- `-snes_lag_jacobian (-2|1|2|...)`                - sets the lag
- `-snes_lag_preconditioner_persists (true|false)` - sets the persistence through multiple SNES solves
- `-snes_lag_preconditioner (-2|1|2|...)`          - sets the lag.

Level: intermediate

See also: `SNES`, `SNESGetLagPreconditioner()`, `SNESSetLagPreconditioner()`, `SNESGetLagJacobianPersists()`, `SNESSetLagPreconditionerPersists()`

# External Links
$(_doc_external("SNES/SNESSetLagJacobian"))
"""
function SNESSetLagJacobian(petsclib::PetscLibType, snes::AbstractSNES, lag::Integer)
    error("SNESSetLagJacobian: no generated method for these argument types")
end

@for_petsc function SNESSetLagJacobian(petsclib::$UnionPetscLib, snes::AbstractSNES, lag::$PetscInt )

    @chk ccall(
               (:SNESSetLagJacobian, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt),
               snes, lag,
              )


	return nothing
end 

"""
	SNESSetLagJacobianPersists(petsclib::PetscLibType, snes::AbstractSNES, flg::PetscBool) 
Set whether or not the Jacobian lagging persists through multiple nonlinear solves

Logically collective

Input Parameters:
- `snes` - the `SNES` context
- `flg`  - jacobian lagging persists if true

Options Database Keys:
- `-snes_lag_jacobian_persists (true|false)`       - sets the persistence through multiple SNES solves
- `-snes_lag_jacobian (-2|1|2|...)`                - sets the lag
- `-snes_lag_preconditioner_persists (true|false)` - sets the persistence through multiple SNES solves
- `-snes_lag_preconditioner (-2|1|2|...)`          - sets the lag

Level: advanced

See also: `SNES`, `SNESSetLagPreconditionerPersists()`, `SNESSetLagJacobian()`, `SNESGetLagJacobian()`, `SNESGetNPC()`

# External Links
$(_doc_external("SNES/SNESSetLagJacobianPersists"))
"""
function SNESSetLagJacobianPersists(petsclib::PetscLibType, snes::AbstractSNES, flg::PetscBool)
    error("SNESSetLagJacobianPersists: no generated method for these argument types")
end

@for_petsc function SNESSetLagJacobianPersists(petsclib::$UnionPetscLib, snes::AbstractSNES, flg::PetscBool )

    @chk ccall(
               (:SNESSetLagJacobianPersists, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, flg,
              )


	return nothing
end 

"""
	SNESSetLagPreconditioner(petsclib::PetscLibType, snes::AbstractSNES, lag::PetscInt) 
Sets when the preconditioner is rebuilt in the nonlinear solve `SNESSolve()`.

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `lag`  - 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time
the Jacobian is built etc. -2 indicates rebuild preconditioner at next chance but then never rebuild after that

Options Database Keys:
- `-snes_lag_jacobian_persists (true|false)`       - sets the persistence through multiple `SNESSolve()`
- `-snes_lag_jacobian (-2|1|2|...)`                - sets the lag
- `-snes_lag_preconditioner_persists (true|false)` - sets the persistence through multiple `SNESSolve()`
- `-snes_lag_preconditioner (-2|1|2|...)`          - sets the lag

Level: intermediate

See also: `SNESGetLagPreconditioner()`, `SNESSetLagJacobian()`, `SNESGetLagJacobian()`, `SNESSetLagPreconditionerPersists()`,
`SNESSetLagJacobianPersists()`, `SNES`, `SNESSolve()`

# External Links
$(_doc_external("SNES/SNESSetLagPreconditioner"))
"""
function SNESSetLagPreconditioner(petsclib::PetscLibType, snes::AbstractSNES, lag::Integer)
    error("SNESSetLagPreconditioner: no generated method for these argument types")
end

@for_petsc function SNESSetLagPreconditioner(petsclib::$UnionPetscLib, snes::AbstractSNES, lag::$PetscInt )

    @chk ccall(
               (:SNESSetLagPreconditioner, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt),
               snes, lag,
              )


	return nothing
end 

"""
	SNESSetLagPreconditionerPersists(petsclib::PetscLibType, snes::AbstractSNES, flg::PetscBool) 
Set whether or not the preconditioner lagging persists through multiple nonlinear solves

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `flg`  - preconditioner lagging persists if true

Options Database Keys:
- `-snes_lag_jacobian_persists (true|false)`       - sets the persistence through multiple SNES solves
- `-snes_lag_jacobian (-2|1|2|...)`                - sets the lag
- `-snes_lag_preconditioner_persists (true|false)` - sets the persistence through multiple SNES solves
- `-snes_lag_preconditioner (-2|1|2|...)`          - sets the lag

Level: developer

See also: `SNES`, `SNESSetLagJacobianPersists()`, `SNESSetLagJacobian()`, `SNESGetLagJacobian()`, `SNESGetNPC()`, `SNESSetLagPreconditioner()`

# External Links
$(_doc_external("SNES/SNESSetLagPreconditionerPersists"))
"""
function SNESSetLagPreconditionerPersists(petsclib::PetscLibType, snes::AbstractSNES, flg::PetscBool)
    error("SNESSetLagPreconditionerPersists: no generated method for these argument types")
end

@for_petsc function SNESSetLagPreconditionerPersists(petsclib::$UnionPetscLib, snes::AbstractSNES, flg::PetscBool )

    @chk ccall(
               (:SNESSetLagPreconditionerPersists, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, flg,
              )


	return nothing
end 

"""
	SNESSetLineSearch(petsclib::PetscLibType, snes::AbstractSNES, linesearch::SNESLineSearch) 
Sets the `SNESLineSearch` to be used for a given `SNES`

Collective

Input Parameters:
- `snes`       - iterative context obtained from `SNESCreate()`
- `linesearch` - the linesearch object

Level: developer

See also: `SNES`, `SNESLineSearch`, `SNESGetLineSearch()`

# External Links
$(_doc_external("SNES/SNESSetLineSearch"))
"""
function SNESSetLineSearch(petsclib::PetscLibType, snes::AbstractSNES, linesearch::SNESLineSearch)
    error("SNESSetLineSearch: no generated method for these argument types")
end

@for_petsc function SNESSetLineSearch(petsclib::$UnionPetscLib, snes::AbstractSNES, linesearch::SNESLineSearch )

    @chk ccall(
               (:SNESSetLineSearch, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESLineSearch),
               snes, linesearch,
              )


	return nothing
end 

"""
	SNESSetMaxLinearSolveFailures(petsclib::PetscLibType, snes::AbstractSNES, maxFails::PetscInt) 
the number of failed linear solve attempts
allowed before `SNES` returns with a diverged reason of `SNES_DIVERGED_LINEAR_SOLVE`

Logically Collective

Input Parameters:
- `snes`     - `SNES` context
- `maxFails` - maximum allowed linear solve failures, use `PETSC_UNLIMITED` to have no limit on the number of failures

Options Database Key:
- `-snes_max_linear_solve_fail num` - The number of failures before the solve is terminated

Level: intermediate

See also: `SNESSetErrorIfNotConverged()`, `SNESGetLinearSolveFailures()`, `SNESGetMaxLinearSolveFailures()`, `SNESGetLinearSolveIterations()`

# External Links
$(_doc_external("SNES/SNESSetMaxLinearSolveFailures"))
"""
function SNESSetMaxLinearSolveFailures(petsclib::PetscLibType, snes::AbstractSNES, maxFails::Integer)
    error("SNESSetMaxLinearSolveFailures: no generated method for these argument types")
end

@for_petsc function SNESSetMaxLinearSolveFailures(petsclib::$UnionPetscLib, snes::AbstractSNES, maxFails::$PetscInt )

    @chk ccall(
               (:SNESSetMaxLinearSolveFailures, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt),
               snes, maxFails,
              )


	return nothing
end 

"""
	SNESSetMaxNonlinearStepFailures(petsclib::PetscLibType, snes::AbstractSNES, maxFails::PetscInt) 
Sets the maximum number of unsuccessful steps
attempted by the nonlinear solver before it gives up and returns unconverged or generates an error

Not Collective

Input Parameters:
- `snes`     - `SNES` context
- `maxFails` - maximum of unsuccessful steps allowed, use `PETSC_UNLIMITED` to have no limit on the number of failures

Options Database Key:
- `-snes_max_fail n` - maximum number of unsuccessful steps allowed

Level: intermediate

See also: `SNESSetErrorIfNotConverged()`, `SNESGetMaxLinearSolveFailures()`, `SNESGetLinearSolveIterations()`, `SNESSetMaxLinearSolveFailures()`,
`SNESGetLinearSolveFailures()`, `SNESGetMaxNonlinearStepFailures()`, `SNESGetNonlinearStepFailures()`, `SNESCheckLineSearchFailure()`

# External Links
$(_doc_external("SNES/SNESSetMaxNonlinearStepFailures"))
"""
function SNESSetMaxNonlinearStepFailures(petsclib::PetscLibType, snes::AbstractSNES, maxFails::Integer)
    error("SNESSetMaxNonlinearStepFailures: no generated method for these argument types")
end

@for_petsc function SNESSetMaxNonlinearStepFailures(petsclib::$UnionPetscLib, snes::AbstractSNES, maxFails::$PetscInt )

    @chk ccall(
               (:SNESSetMaxNonlinearStepFailures, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt),
               snes, maxFails,
              )


	return nothing
end 

"""
	SNESSetNGS(petsclib::PetscLibType, snes::AbstractSNES, f::Ptr{Cvoid}, ctx::Ptr{Cvoid}) 
Sets the user nonlinear Gauss-Seidel routine for
use with composed nonlinear solvers.

Input Parameters:
- `snes` - the `SNES` context, usually of the `SNESType` `SNESNGS`
- `f`    - function evaluation routine to apply Gauss-Seidel, see `SNESNGSFn` for calling sequence
- `ctx`  - [optional] user-defined context for private data for the smoother evaluation routine (may be `NULL`)

Level: intermediate

See also: `SNESNGS`, `SNESGetNGS()`, `SNESNCG`, `SNESGetFunction()`, `SNESComputeNGS()`, `SNESNGSFn`

# External Links
$(_doc_external("SNES/SNESSetNGS"))
"""
function SNESSetNGS(petsclib::PetscLibType, snes::AbstractSNES, f::Ptr{Cvoid}, ctx::Ptr{Cvoid})
    error("SNESSetNGS: no generated method for these argument types")
end

@for_petsc function SNESSetNGS(petsclib::$UnionPetscLib, snes::AbstractSNES, f::Ptr{Cvoid}, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESSetNGS, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Cvoid}, Ptr{Cvoid}),
               snes, f, ctx,
              )


	return nothing
end 

"""
	SNESSetNPC(petsclib::PetscLibType, snes::AbstractSNES, npc::AbstractSNES) 
Sets the nonlinear preconditioner to be used.

Collective

Input Parameters:
- `snes` - iterative context obtained from `SNESCreate()`
- `npc`  - the `SNES` nonlinear preconditioner object

Options Database Key:
- `-npc_snes_type type` - set the type of the `SNES` to use as the nonlinear preconditioner

Level: developer

See also: `SNES`, `SNESNGS`, `SNESFAS`, `SNESGetNPC()`, `SNESHasNPC()`

# External Links
$(_doc_external("SNES/SNESSetNPC"))
"""
function SNESSetNPC(petsclib::PetscLibType, snes::AbstractSNES, npc::AbstractSNES)
    error("SNESSetNPC: no generated method for these argument types")
end

@for_petsc function SNESSetNPC(petsclib::$UnionPetscLib, snes::AbstractSNES, npc::AbstractSNES )

    @chk ccall(
               (:SNESSetNPC, $petsc_library),
               PetscErrorCode,
               (CSNES, CSNES),
               snes, npc,
              )


	return nothing
end 

"""
	SNESSetNPCSide(petsclib::PetscLibType, snes::AbstractSNES, side::PCSide) 
Sets the nonlinear preconditioning side used by the nonlinear preconditioner inside `SNES`.

Logically Collective

Input Parameter:
- `snes` - iterative context obtained from `SNESCreate()`

Output Parameter:
- `side` - the preconditioning side, where side is one of
``
PC_LEFT  - left preconditioning
PC_RIGHT - right preconditioning (default for most nonlinear solvers)
``

Options Database Key:
- `-snes_npc_side (right|left)` - nonlinear preconditioner side

Level: intermediate

See also: `SNES`, `SNESGetNPC()`, `SNESNRICHARDSON`, `SNESNCG`, `SNESType`, `SNESGetNPCSide()`, `KSPSetPCSide()`, `PC_LEFT`, `PC_RIGHT`, `PCSide`

# External Links
$(_doc_external("SNES/SNESSetNPCSide"))
"""
function SNESSetNPCSide(petsclib::PetscLibType, snes::AbstractSNES, side::PCSide)
    error("SNESSetNPCSide: no generated method for these argument types")
end

@for_petsc function SNESSetNPCSide(petsclib::$UnionPetscLib, snes::AbstractSNES, side::PCSide )

    @chk ccall(
               (:SNESSetNPCSide, $petsc_library),
               PetscErrorCode,
               (CSNES, PCSide),
               snes, side,
              )


	return nothing
end 

"""
	SNESSetNormSchedule(petsclib::PetscLibType, snes::AbstractSNES, normschedule::SNESNormSchedule) 
Sets the `SNESNormSchedule` used in convergence and monitoring
of the `SNES` method, when norms are computed in the solving process

Logically Collective

Input Parameters:
- `snes`         - the `SNES` context
- `normschedule` - the frequency of norm computation

Options Database Key:
- `-snes_norm_schedule (none|always|initialonly|finalonly|initialfinalonly)` - set the schedule

Level: advanced

See also: `SNESNormSchedule`, `SNESGetNormSchedule()`, `SNESComputeFunction()`, `VecNorm()`, `SNESSetFunction()`, `SNESSetInitialFunction()`

# External Links
$(_doc_external("SNES/SNESSetNormSchedule"))
"""
function SNESSetNormSchedule(petsclib::PetscLibType, snes::AbstractSNES, normschedule::SNESNormSchedule)
    error("SNESSetNormSchedule: no generated method for these argument types")
end

@for_petsc function SNESSetNormSchedule(petsclib::$UnionPetscLib, snes::AbstractSNES, normschedule::SNESNormSchedule )

    @chk ccall(
               (:SNESSetNormSchedule, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESNormSchedule),
               snes, normschedule,
              )


	return nothing
end 

"""
	SNESSetObjective(petsclib::PetscLibType, snes::AbstractSNES, obj::Ptr{Cvoid}, ctx::Ptr{Cvoid}) 
Sets the objective function minimized by some of the `SNES` linesearch methods, used instead of the 2-norm of the residual in the line search

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `obj`  - objective evaluation routine; see `SNESObjectiveFn` for the calling sequence
- `ctx`  - [optional] user-defined context for private data for the objective evaluation routine (may be `NULL`)

Level: intermediate

See also: `SNES`, `SNESLineSearch()`, `SNESGetObjective()`, `SNESComputeObjective()`, `SNESSetFunction()`, `SNESSetJacobian()`,
`SNESObjectiveFn`, `SNESSetObjectiveDomainError()`

# External Links
$(_doc_external("SNES/SNESSetObjective"))
"""
function SNESSetObjective(petsclib::PetscLibType, snes::AbstractSNES, obj::Ptr{Cvoid}, ctx::Ptr{Cvoid})
    error("SNESSetObjective: no generated method for these argument types")
end

@for_petsc function SNESSetObjective(petsclib::$UnionPetscLib, snes::AbstractSNES, obj::Ptr{Cvoid}, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESSetObjective, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Cvoid}, Ptr{Cvoid}),
               snes, obj, ctx,
              )


	return nothing
end 

"""
	SNESSetObjectiveDomainError(petsclib::PetscLibType, snes::AbstractSNES) 
tells `SNES` that the input vector, a proposed new solution, to your function you provided to `SNESSetObjective()` is not
in the function's domain. For example, a step with negative pressure.

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Level: advanced

See also: `SNESCreate()`, `SNESSetFunction()`, `SNESFunctionFn`, `SNESSetJacobianDomainError()`, `SNESVISetVariableBounds()`,
`SNESVISetComputeVariableBounds()`, `SNESLineSearchSetPreCheck()`, `SNESLineSearchSetPostCheck()`, `SNESConvergedReason`, `SNESGetConvergedReason()`,
`SNES_DIVERGED_OBJECTIVE_DOMAIN`, `SNESSetFunctionDomainError()`, `SNES_DIVERGED_FUNCTION_DOMAIN`

# External Links
$(_doc_external("SNES/SNESSetObjectiveDomainError"))
"""
function SNESSetObjectiveDomainError(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESSetObjectiveDomainError: no generated method for these argument types")
end

@for_petsc function SNESSetObjectiveDomainError(petsclib::$UnionPetscLib, snes::AbstractSNES )

    @chk ccall(
               (:SNESSetObjectiveDomainError, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	SNESSetOptionsPrefix(petsclib::PetscLibType, snes::AbstractSNES, prefix::String) 
Sets the prefix used for searching for all
`SNES` options in the database.

Logically Collective

Input Parameters:
- `snes`   - the `SNES` context
- `prefix` - the prefix to prepend to all option names

Level: advanced

See also: `SNES`, `SNESSetFromOptions()`, `SNESAppendOptionsPrefix()`

# External Links
$(_doc_external("SNES/SNESSetOptionsPrefix"))
"""
function SNESSetOptionsPrefix(petsclib::PetscLibType, snes::AbstractSNES, prefix::String)
    error("SNESSetOptionsPrefix: no generated method for these argument types")
end

@for_petsc function SNESSetOptionsPrefix(petsclib::$UnionPetscLib, snes::AbstractSNES, prefix::String )

    @chk ccall(
               (:SNESSetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Cchar}),
               snes, prefix,
              )


	return nothing
end 

"""
	SNESSetPicard(petsclib::PetscLibType, snes::AbstractSNES, r::AbstractPetscVec, bp::Ptr{Cvoid}, Amat::AbstractPetscMat, Pmat::AbstractPetscMat, J::Ptr{Cvoid}, ctx::Ptr{Cvoid}) 
Use `SNES` to solve the system A(x) x = bp(x) + b  via a Picard type iteration (Picard linearization)

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `r`    - vector to store function values, may be `NULL`
- `bp`   - function evaluation routine, may be `NULL`, for the calling sequence see `SNESFunctionFn`
- `Amat` - matrix with which A(x) x - bp(x) - b is to be computed
- `Pmat` - matrix from which preconditioner is computed (usually the same as `Amat`)
- `J`    - function to compute matrix values, for the calling sequence see `SNESJacobianFn`
- `ctx`  - [optional] user-defined context for private data for the function evaluation routine (may be `NULL`)

Level: intermediate

See also: `SNES`, `SNESGetFunction()`, `SNESSetFunction()`, `SNESComputeFunction()`, `SNESSetJacobian()`, `SNESGetPicard()`, `SNESLineSearchPreCheckPicard()`,
`SNESFunctionFn`, `SNESJacobianFn`

# External Links
$(_doc_external("SNES/SNESSetPicard"))
"""
function SNESSetPicard(petsclib::PetscLibType, snes::AbstractSNES, r::AbstractPetscVec, bp::Ptr{Cvoid}, Amat::AbstractPetscMat, Pmat::AbstractPetscMat, J::Ptr{Cvoid}, ctx::Ptr{Cvoid})
    error("SNESSetPicard: no generated method for these argument types")
end

@for_petsc function SNESSetPicard(petsclib::$UnionPetscLib, snes::AbstractSNES, r::AbstractPetscVec, bp::Ptr{Cvoid}, Amat::AbstractPetscMat, Pmat::AbstractPetscMat, J::Ptr{Cvoid}, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESSetPicard, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, Ptr{Cvoid}, CMat, CMat, Ptr{Cvoid}, Ptr{Cvoid}),
               snes, r, bp, Amat, Pmat, J, ctx,
              )


	return nothing
end 

"""
	SNESSetSolution(petsclib::PetscLibType, snes::AbstractSNES, u::AbstractPetscVec) 
Sets the solution vector for use by the `SNES` routines.

Logically Collective

Input Parameters:
- `snes` - the `SNES` context obtained from `SNESCreate()`
- `u`    - the solution vector

Level: beginner

See also: `SNES`, `SNESSolve()`, `SNESGetSolution()`, `Vec`

# External Links
$(_doc_external("SNES/SNESSetSolution"))
"""
function SNESSetSolution(petsclib::PetscLibType, snes::AbstractSNES, u::AbstractPetscVec)
    error("SNESSetSolution: no generated method for these argument types")
end

@for_petsc function SNESSetSolution(petsclib::$UnionPetscLib, snes::AbstractSNES, u::AbstractPetscVec )

    @chk ccall(
               (:SNESSetSolution, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec),
               snes, u,
              )


	return nothing
end 

"""
	SNESSetTolerances(petsclib::PetscLibType, snes::AbstractSNES, abstol::PetscReal, rtol::PetscReal, stol::PetscReal, maxit::PetscInt, maxf::PetscInt) 
Sets various parameters used in `SNES` convergence tests.

Logically Collective

Input Parameters:
- `snes`   - the `SNES` context
- `abstol` - the absolute convergence tolerance,  F(x^n) \\le abstol 
- `rtol`   - the relative convergence tolerance,  F(x^n) \\le reltol * F(x^0) 
- `stol`   - convergence tolerance in terms of the norm of the change in the solution between steps,  || delta x || < stol*|| x ||
- `maxit`  - the maximum number of iterations allowed in the solver, default 50.
- `maxf`   - the maximum number of function evaluations allowed in the solver (use `PETSC_UNLIMITED` indicates no limit), default 10,000

Options Database Keys:
- `-snes_atol abstol`    - Sets `abstol`
- `-snes_rtol rtol`      - Sets `rtol`
- `-snes_stol stol`      - Sets `stol`
- `-snes_max_it maxit`   - Sets `maxit`
- `-snes_max_funcs maxf` - Sets `maxf` (use `unlimited` to have no maximum)

Level: intermediate

See also: `SNESSolve()`, `SNES`, `SNESSetDivergenceTolerance()`, `SNESSetForceIteration()`

# External Links
$(_doc_external("SNES/SNESSetTolerances"))
"""
function SNESSetTolerances(petsclib::PetscLibType, snes::AbstractSNES, abstol::Real, rtol::Real, stol::Real, maxit::Integer, maxf::Integer)
    error("SNESSetTolerances: no generated method for these argument types")
end

@for_petsc function SNESSetTolerances(petsclib::$UnionPetscLib, snes::AbstractSNES, abstol::$PetscReal, rtol::$PetscReal, stol::$PetscReal, maxit::$PetscInt, maxf::$PetscInt )

    @chk ccall(
               (:SNESSetTolerances, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscReal, $PetscReal, $PetscReal, $PetscInt, $PetscInt),
               snes, abstol, rtol, stol, maxit, maxf,
              )


	return nothing
end 

"""
	SNESSetType(petsclib::PetscLibType, snes::AbstractSNES, type::SNESType) 
Sets the algorithm/method to be used to solve the nonlinear system with the given `SNES`

Collective

Input Parameters:
- `snes` - the `SNES` context
- `type` - a known method

Options Database Key:
- `-snes_type type` - Sets the method; see `SNESType`

Level: intermediate

See also: `SNES`, `SNESSolve()`, `SNESType`, `SNESCreate()`, `SNESDestroy()`, `SNESGetType()`, `SNESSetFromOptions()`

# External Links
$(_doc_external("SNES/SNESSetType"))
"""
function SNESSetType(petsclib::PetscLibType, snes::AbstractSNES, type::SNESType)
    error("SNESSetType: no generated method for these argument types")
end

@for_petsc function SNESSetType(petsclib::$UnionPetscLib, snes::AbstractSNES, type::SNESType )

    @chk ccall(
               (:SNESSetType, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESType),
               snes, type,
              )


	return nothing
end 

"""
	SNESSetUp(petsclib::PetscLibType, snes::AbstractSNES) 
Sets up the internal data structures for the later use
of a nonlinear solver `SNESSolve()`.

Collective

Input Parameter:
- `snes` - the `SNES` context

Level: advanced

See also: `SNES`, `SNESCreate()`, `SNESSolve()`, `SNESDestroy()`, `SNESSetFromOptions()`

# External Links
$(_doc_external("SNES/SNESSetUp"))
"""
function SNESSetUp(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESSetUp: no generated method for these argument types")
end

@for_petsc function SNESSetUp(petsclib::$UnionPetscLib, snes::AbstractSNES )

    @chk ccall(
               (:SNESSetUp, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	SNESSetUpMatrices(petsclib::PetscLibType, snes::AbstractSNES) 
ensures that matrices are available for `SNES` Newton-like methods, this is called by `SNESSetUp_XXX()`

Collective

Input Parameter:
- `snes` - `SNES` object to configure

Level: developer

See also: `SNES`, `SNESSetUp()`

# External Links
$(_doc_external("SNES/SNESSetUpMatrices"))
"""
function SNESSetUpMatrices(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESSetUpMatrices: no generated method for these argument types")
end

@for_petsc function SNESSetUpMatrices(petsclib::$UnionPetscLib, snes::AbstractSNES )

    @chk ccall(
               (:SNESSetUpMatrices, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	SNESSetUpdate(petsclib::PetscLibType, snes::AbstractSNES, func::Ptr{Cvoid}) 
Sets the general-purpose update function called
at the beginning of every iteration of the nonlinear solve. Specifically
it is called just before the Jacobian is "evaluated" and after the function
evaluation.

Logically Collective

Input Parameters:
- `snes` - The nonlinear solver context
- `func` - The update function; for calling sequence see `SNESUpdateFn`

Level: advanced

See also: `SNES`, `SNESSolve()`, `SNESSetJacobian()`, `SNESLineSearchSetPreCheck()`, `SNESLineSearchSetPostCheck()`, `SNESNewtonTRSetPreCheck()`, `SNESNewtonTRSetPostCheck()`,
`SNESMonitorSet()`

# External Links
$(_doc_external("SNES/SNESSetUpdate"))
"""
function SNESSetUpdate(petsclib::PetscLibType, snes::AbstractSNES, func::Ptr{Cvoid})
    error("SNESSetUpdate: no generated method for these argument types")
end

@for_petsc function SNESSetUpdate(petsclib::$UnionPetscLib, snes::AbstractSNES, func::Ptr{Cvoid} )

    @chk ccall(
               (:SNESSetUpdate, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Cvoid}),
               snes, func,
              )


	return nothing
end 

"""
	SNESSetUseMatrixFree(petsclib::PetscLibType, snes::AbstractSNES, mf_operator::PetscBool, mf::PetscBool) 
indicates that `SNES` should use matrix-free finite difference matrix-vector products to apply the Jacobian.

Logically Collective

Input Parameters:
- `snes`        - `SNES` context
- `mf_operator` - use matrix-free only for the Amat used by `SNESSetJacobian()`, this means the user provided Pmat will continue to be used
- `mf`          - use matrix-free for both the Amat and Pmat used by `SNESSetJacobian()`, both the Amat and Pmat set in `SNESSetJacobian()` will be ignored. With
this option no matrix-element based preconditioners can be used in the linear solve since the matrix won't be explicitly available

Options Database Keys:
- `-snes_mf_operator` - use matrix-free only for the mat operator
- `-snes_mf`          - use matrix-free for both the mat and pmat operator
- `-snes_fd_color`    - compute the Jacobian via coloring and finite differences.
- `-snes_fd`          - compute the Jacobian via finite differences (slow)

Level: intermediate

See also: `SNES`, `SNESGetUseMatrixFree()`, `MatCreateSNESMF()`, `SNESComputeJacobianDefaultColor()`, `MatFDColoring`

# External Links
$(_doc_external("SNES/SNESSetUseMatrixFree"))
"""
function SNESSetUseMatrixFree(petsclib::PetscLibType, snes::AbstractSNES, mf_operator::PetscBool, mf::PetscBool)
    error("SNESSetUseMatrixFree: no generated method for these argument types")
end

@for_petsc function SNESSetUseMatrixFree(petsclib::$UnionPetscLib, snes::AbstractSNES, mf_operator::PetscBool, mf::PetscBool )

    @chk ccall(
               (:SNESSetUseMatrixFree, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool, PetscBool),
               snes, mf_operator, mf,
              )


	return nothing
end 

"""
	SNESSetWorkVecs(petsclib::PetscLibType, snes::AbstractSNES, nw::PetscInt) 
Allocates a number of work vectors to be used internally by the `SNES` solver

Input Parameters:
- `snes` - the `SNES` context
- `nw`   - number of work vectors to allocate

Level: developer

See also: `SNES`

# External Links
$(_doc_external("SNES/SNESSetWorkVecs"))
"""
function SNESSetWorkVecs(petsclib::PetscLibType, snes::AbstractSNES, nw::Integer)
    error("SNESSetWorkVecs: no generated method for these argument types")
end

@for_petsc function SNESSetWorkVecs(petsclib::$UnionPetscLib, snes::AbstractSNES, nw::$PetscInt )

    @chk ccall(
               (:SNESSetWorkVecs, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt),
               snes, nw,
              )


	return nothing
end 

"""
	ctx::Ptr{Cvoid} = SNESShellGetContext(petsclib::PetscLibType, snes::AbstractSNES) 
Returns the user-provided context associated with a `SNESSHELL`

Not Collective

Input Parameter:
- `snes` - should have been created with `SNESSetType`(snes,`SNESSHELL`);

Output Parameter:
- `ctx` - the user provided context

Level: advanced

See also: `SNES`, `SNESSHELL`, `SNESCreateShell()`, `SNESShellSetContext()`

# External Links
$(_doc_external("SNES/SNESShellGetContext"))
"""
function SNESShellGetContext(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESShellGetContext: no generated method for these argument types")
end

@for_petsc function SNESShellGetContext(petsclib::$UnionPetscLib, snes::AbstractSNES )
	ctx_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:SNESShellGetContext, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Cvoid}),
               snes, ctx_,
              )

	ctx = ctx_[]

	return ctx
end 

"""
	SNESShellSetContext(petsclib::PetscLibType, snes::AbstractSNES, ctx::Ptr{Cvoid}) 
sets the context for a `SNESSHELL`

Logically Collective

Input Parameters:
- `snes` - the `SNESSHELL`
- `ctx`  - the context

Level: advanced

See also: `SNES`, `SNESSHELL`, `SNESCreateShell()`, `SNESShellGetContext()`

# External Links
$(_doc_external("SNES/SNESShellSetContext"))
"""
function SNESShellSetContext(petsclib::PetscLibType, snes::AbstractSNES, ctx::Ptr{Cvoid})
    error("SNESShellSetContext: no generated method for these argument types")
end

@for_petsc function SNESShellSetContext(petsclib::$UnionPetscLib, snes::AbstractSNES, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESShellSetContext, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Cvoid}),
               snes, ctx,
              )


	return nothing
end 

"""
	SNESShellSetSolve(petsclib::PetscLibType, snes::AbstractSNES, solve::external) 
Sets routine to apply as solver to a `SNESSHELL` `SNES` object

Logically Collective

Input Parameters:
- `snes`  - the `SNES` nonlinear solver context
- `solve` - the application-provided solver routine

Calling sequence of `apply`:
- `snes` - the preconditioner, get the application context with `SNESShellGetContext()` provided with `SNESShellSetContext()`
- `xout` - solution vector

Level: advanced

See also: `SNES`, `SNESSHELL`, `SNESShellSetContext()`, `SNESShellGetContext()`

# External Links
$(_doc_external("SNES/SNESShellSetSolve"))
"""
function SNESShellSetSolve(petsclib::PetscLibType, snes::AbstractSNES, solve::external)
    error("SNESShellSetSolve: no generated method for these argument types")
end

@for_petsc function SNESShellSetSolve(petsclib::$UnionPetscLib, snes::AbstractSNES, solve::external )

    @chk ccall(
               (:SNESShellSetSolve, $petsc_library),
               PetscErrorCode,
               (CSNES, external),
               snes, solve,
              )


	return nothing
end 

"""
	SNESSolve(petsclib::PetscLibType, snes::AbstractSNES, b::Union{Ptr, AbstractPetscVec}, x::AbstractPetscVec) 
Solves a nonlinear system F(x) = b  associated with a `SNES` object

Collective

Input Parameters:
- `snes` - the `SNES` context
- `b`    - the constant part of the equation F(x) = b, or `NULL` to use zero.
- `x`    - the solution vector.

Level: beginner

See also: `SNES`, `SNESCreate()`, `SNESDestroy()`, `SNESSetFunction()`, `SNESSetJacobian()`, `SNESSetGridSequence()`, `SNESGetSolution()`,
`SNESNewtonTRSetPreCheck()`, `SNESNewtonTRGetPreCheck()`, `SNESNewtonTRSetPostCheck()`, `SNESNewtonTRGetPostCheck()`,
`SNESLineSearchSetPostCheck()`, `SNESLineSearchGetPostCheck()`, `SNESLineSearchSetPreCheck()`, `SNESLineSearchGetPreCheck()`

# External Links
$(_doc_external("SNES/SNESSolve"))
"""
function SNESSolve(petsclib::PetscLibType, snes::AbstractSNES, b::Union{Ptr, AbstractPetscVec}, x::AbstractPetscVec)
    error("SNESSolve: no generated method for these argument types")
end

@for_petsc function SNESSolve(petsclib::$UnionPetscLib, snes::AbstractSNES, b::Union{Ptr, AbstractPetscVec}, x::AbstractPetscVec )

    @chk ccall(
               (:SNESSolve, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec),
               snes, b, x,
              )


	return nothing
end 

"""
	SNESTSFormFunction(petsclib::PetscLibType, snes::AbstractSNES, U::AbstractPetscVec, F::AbstractPetscVec, ctx::Ptr{Cvoid}) 
Function to evaluate nonlinear residual defined by an ODE solver algorithm implemented within `TS`

Logically Collective

Input Parameters:
- `snes` - nonlinear solver
- `U`    - the current state at which to evaluate the residual
- `ctx`  - application context, must be a `TS`

Output Parameter:
- `F` - the nonlinear residual

Level: developer

See also: `SNESSetFunction()`, `MatFDColoringSetFunction()`

# External Links
$(_doc_external("TS/SNESTSFormFunction"))
"""
function SNESTSFormFunction(petsclib::PetscLibType, snes::AbstractSNES, U::AbstractPetscVec, F::AbstractPetscVec, ctx::Ptr{Cvoid})
    error("SNESTSFormFunction: no generated method for these argument types")
end

@for_petsc function SNESTSFormFunction(petsclib::$UnionPetscLib, snes::AbstractSNES, U::AbstractPetscVec, F::AbstractPetscVec, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESTSFormFunction, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec, Ptr{Cvoid}),
               snes, U, F, ctx,
              )


	return nothing
end 

"""
	SNESTSFormJacobian(petsclib::PetscLibType, snes::AbstractSNES, U::AbstractPetscVec, A::AbstractPetscMat, B::AbstractPetscMat, ctx::Ptr{Cvoid}) 
Function to evaluate the Jacobian defined by an ODE solver algorithm implemented within `TS`

Collective

Input Parameters:
- `snes` - nonlinear solver
- `U`    - the current state at which to evaluate the residual
- `ctx`  - application context, must be a `TS`

Output Parameters:
- `A` - the Jacobian
- `B` - the matrix used to construct the preconditioner (often the same as `A`)

Level: developer

See also: `SNESSetJacobian()`

# External Links
$(_doc_external("TS/SNESTSFormJacobian"))
"""
function SNESTSFormJacobian(petsclib::PetscLibType, snes::AbstractSNES, U::AbstractPetscVec, A::AbstractPetscMat, B::AbstractPetscMat, ctx::Ptr{Cvoid})
    error("SNESTSFormJacobian: no generated method for these argument types")
end

@for_petsc function SNESTSFormJacobian(petsclib::$UnionPetscLib, snes::AbstractSNES, U::AbstractPetscVec, A::AbstractPetscMat, B::AbstractPetscMat, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESTSFormJacobian, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CMat, CMat, Ptr{Cvoid}),
               snes, U, A, B, ctx,
              )


	return nothing
end 

"""
	SNESTestFunction(petsclib::PetscLibType, snes::AbstractSNES) 
Computes the difference between the computed and finite-difference functions

Collective

Input Parameter:
- `snes` - the `SNES` context

Options Database Keys:
- `-snes_test_function`      - compare the user provided function with one compute via finite differences to check for errors.
- `-snes_test_function_view` - display the user provided function, the finite difference function and the difference

Level: developer

See also: `SNESTestJacobian()`, `SNESSetFunction()`, `SNESComputeFunction()`

# External Links
$(_doc_external("SNES/SNESTestFunction"))
"""
function SNESTestFunction(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESTestFunction: no generated method for these argument types")
end

@for_petsc function SNESTestFunction(petsclib::$UnionPetscLib, snes::AbstractSNES )

    @chk ccall(
               (:SNESTestFunction, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	Jnorm::PetscReal,diffNorm::PetscReal = SNESTestJacobian(petsclib::PetscLibType, snes::AbstractSNES) 
Computes the difference between the computed and finite-difference Jacobians

Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameters:
- `Jnorm`    - the Frobenius norm of the computed Jacobian, or `NULL`
- `diffNorm` - the Frobenius norm of the difference of the computed and finite-difference Jacobians, or `NULL`

Options Database Keys:
- `-snes_test_jacobian [threshold]` - compare the user provided Jacobian with one compute via finite differences to check for errors.  If a threshold is given, display only those entries whose difference is greater than the threshold.
- `-snes_test_jacobian_view`        - display the user provided Jacobian, the finite difference Jacobian and the difference

Level: developer

See also: `SNESTestFunction()`, `SNESSetJacobian()`, `SNESComputeJacobian()`

# External Links
$(_doc_external("SNES/SNESTestJacobian"))
"""
function SNESTestJacobian(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESTestJacobian: no generated method for these argument types")
end

@for_petsc function SNESTestJacobian(petsclib::$UnionPetscLib, snes::AbstractSNES )
	Jnorm_ = Ref{$PetscReal}()
	diffNorm_ = Ref{$PetscReal}()

    @chk ccall(
               (:SNESTestJacobian, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscReal}, Ptr{$PetscReal}),
               snes, Jnorm_, diffNorm_,
              )

	Jnorm = Jnorm_[]
	diffNorm = diffNorm_[]

	return Jnorm,diffNorm
end 

"""
	SNESTestLocalMin(petsclib::PetscLibType, snes::AbstractSNES) 
Diagnostic that probes each entry of the current `SNES` solution to check whether the residual norm has a local minimum along the coordinate directions

Collective

Input Parameter:
- `snes` - the `SNES` context

Level: developer

See also: `SNES`, `SNESSolve()`, `SNESComputeFunction()`

# External Links
$(_doc_external("SNES/SNESTestLocalMin"))
"""
function SNESTestLocalMin(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESTestLocalMin: no generated method for these argument types")
end

@for_petsc function SNESTestLocalMin(petsclib::$UnionPetscLib, snes::AbstractSNES )

    @chk ccall(
               (:SNESTestLocalMin, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	SNESVIComputeFunction(petsclib::PetscLibType, snes::AbstractSNES, X::AbstractPetscVec, phi::AbstractPetscVec, functx::Ptr{Cvoid}) 
Provides the function that reformulates a system of nonlinear equations in mixed complementarity form to a system of nonlinear
equations in semismooth form.

Input Parameters:
- `snes`   - the `SNES` context
- `X`      - current iterate
- `functx` - user defined function context

Output Parameter:
- `phi` - the evaluation of semismooth function at `X`

Level: developer

See also: `SNES`, `SNESVINEWTONSSLS`, `SNESVIComputeMeritFunction()`

# External Links
$(_doc_external("SNES/SNESVIComputeFunction"))
"""
function SNESVIComputeFunction(petsclib::PetscLibType, snes::AbstractSNES, X::AbstractPetscVec, phi::AbstractPetscVec, functx::Ptr{Cvoid})
    error("SNESVIComputeFunction: no generated method for these argument types")
end

@for_petsc function SNESVIComputeFunction(petsclib::$UnionPetscLib, snes::AbstractSNES, X::AbstractPetscVec, phi::AbstractPetscVec, functx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESVIComputeFunction, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec, Ptr{Cvoid}),
               snes, X, phi, functx,
              )


	return nothing
end 

"""
	fnorm::PetscReal = SNESVIComputeInactiveSetFnorm(petsclib::PetscLibType, snes::AbstractSNES, F::AbstractPetscVec, X::AbstractPetscVec) 
Computes the function norm for variational inequalities on the inactive set

Input Parameters:
- `snes` - the `SNES` context
- `F`    - the nonlinear function vector
- `X`    - the `SNES` solution vector

Output Parameter:
- `fnorm` - the function norm

Level: developer

See also: `SNES`, `SNESVINEWTONRSLS`, `SNESVINEWTONSSLS`, `SNESLineSearchSetVIFunctions()`

# External Links
$(_doc_external("SNES/SNESVIComputeInactiveSetFnorm"))
"""
function SNESVIComputeInactiveSetFnorm(petsclib::PetscLibType, snes::AbstractSNES, F::AbstractPetscVec, X::AbstractPetscVec)
    error("SNESVIComputeInactiveSetFnorm: no generated method for these argument types")
end

@for_petsc function SNESVIComputeInactiveSetFnorm(petsclib::$UnionPetscLib, snes::AbstractSNES, F::AbstractPetscVec, X::AbstractPetscVec )
	fnorm_ = Ref{$PetscReal}()

    @chk ccall(
               (:SNESVIComputeInactiveSetFnorm, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec, Ptr{$PetscReal}),
               snes, F, X, fnorm_,
              )

	fnorm = fnorm_[]

	return fnorm
end 

"""
	fty::PetscScalar = SNESVIComputeInactiveSetFtY(petsclib::PetscLibType, snes::AbstractSNES, F::AbstractPetscVec, X::AbstractPetscVec, Y::AbstractPetscVec) 
Computes the directional derivative for variational inequalities on the inactive set,
assuming that there exists some G(x) for which the `SNESFunctionFn` F(x) = grad G(x) (relevant for some line search algorithms)

Input Parameters:
- `snes` - the `SNES` context
- `F`    - the nonlinear function vector
- `X`    - the `SNES` solution vector
- `Y`    - the direction vector

Output Parameter:
- `fty` - the directional derivative

Level: developer

See also: `SNES`, `SNESVINEWTONRSLS`, `SNESVINEWTONSSLS`

# External Links
$(_doc_external("SNES/SNESVIComputeInactiveSetFtY"))
"""
function SNESVIComputeInactiveSetFtY(petsclib::PetscLibType, snes::AbstractSNES, F::AbstractPetscVec, X::AbstractPetscVec, Y::AbstractPetscVec)
    error("SNESVIComputeInactiveSetFtY: no generated method for these argument types")
end

@for_petsc function SNESVIComputeInactiveSetFtY(petsclib::$UnionPetscLib, snes::AbstractSNES, F::AbstractPetscVec, X::AbstractPetscVec, Y::AbstractPetscVec )
	fty_ = Ref{$PetscScalar}()

    @chk ccall(
               (:SNESVIComputeInactiveSetFtY, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec, CVec, Ptr{$PetscScalar}),
               snes, F, X, Y, fty_,
              )

	fty = fty_[]

	return fty
end 

"""
	merit::PetscReal,phinorm::PetscReal = SNESVIComputeMeritFunction(petsclib::PetscLibType, phi::AbstractPetscVec) 
Evaluates the merit function for the mixed complementarity problem.

Input Parameter:
- `phi` - the `Vec` holding the evaluation of the semismooth function

Output Parameters:
- `merit`   - the merit function 1/2 ||phi||^2
- `phinorm` - the two-norm of the vector, ||phi||

Level: developer

See also: `SNES`, `SNESVINEWTONSSLS`, `SNESVIComputeFunction()`

# External Links
$(_doc_external("SNES/SNESVIComputeMeritFunction"))
"""
function SNESVIComputeMeritFunction(petsclib::PetscLibType, phi::AbstractPetscVec)
    error("SNESVIComputeMeritFunction: no generated method for these argument types")
end

@for_petsc function SNESVIComputeMeritFunction(petsclib::$UnionPetscLib, phi::AbstractPetscVec )
	merit_ = Ref{$PetscReal}()
	phinorm_ = Ref{$PetscReal}()

    @chk ccall(
               (:SNESVIComputeMeritFunction, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscReal}, Ptr{$PetscReal}),
               phi, merit_, phinorm_,
              )

	merit = merit_[]
	phinorm = phinorm_[]

	return merit,phinorm
end 

"""
	ISact::IS = SNESVIGetActiveSetIS(petsclib::PetscLibType, snes::AbstractSNES, X::AbstractPetscVec, F::AbstractPetscVec) 
Gets the global indices for the active set variables

Input Parameters:
- `snes` - the `SNES` context
- `X`    - the `snes` solution vector
- `F`    - the nonlinear function vector

Output Parameter:
- `ISact` - active set index set

Level: developer

See also: `SNES`, `SNESVINEWTONRSLS`, `SNESVINEWTONSSLS`

# External Links
$(_doc_external("SNES/SNESVIGetActiveSetIS"))
"""
function SNESVIGetActiveSetIS(petsclib::PetscLibType, snes::AbstractSNES, X::AbstractPetscVec, F::AbstractPetscVec)
    error("SNESVIGetActiveSetIS: no generated method for these argument types")
end

@for_petsc function SNESVIGetActiveSetIS(petsclib::$UnionPetscLib, snes::AbstractSNES, X::AbstractPetscVec, F::AbstractPetscVec )
	ISact_ = Ref{CIS}()

    @chk ccall(
               (:SNESVIGetActiveSetIS, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec, Ptr{CIS}),
               snes, X, F, ISact_,
              )

	ISact = IS(ISact_[], petsclib)

	return ISact
end 

"""
	inact::IS = SNESVIGetInactiveSet(petsclib::PetscLibType, snes::AbstractSNES) 
Gets the global indices for the inactive set variables (these correspond to the degrees of freedom the linear
system is solved on)

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `inact` - inactive set index set

Level: advanced

See also: `SNES`, `SNESVINEWTONRSLS`

# External Links
$(_doc_external("SNES/SNESVIGetInactiveSet"))
"""
function SNESVIGetInactiveSet(petsclib::PetscLibType, snes::AbstractSNES)
    error("SNESVIGetInactiveSet: no generated method for these argument types")
end

@for_petsc function SNESVIGetInactiveSet(petsclib::$UnionPetscLib, snes::AbstractSNES )
	inact_ = Ref{CIS}()

    @chk ccall(
               (:SNESVIGetInactiveSet, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CIS}),
               snes, inact_,
              )

	inact = IS(inact_[], petsclib)

	return inact
end 

"""
	SNESVIGetVariableBounds(petsclib::PetscLibType, snes::AbstractSNES, xl::AbstractPetscVec, xu::AbstractPetscVec) 
Gets the lower and upper bounds for the solution vector. `xl` <= x <= `xu`. These are used in solving
(differential) variable inequalities.

Input Parameters:
- `snes` - the `SNES` context.
- `xl`   - lower bound (may be `NULL`)
- `xu`   - upper bound (may be `NULL`)

Level: advanced

See also: [](sec_vi), `SNES`, `SNESVISetVariableBounds()`, `SNESVISetComputeVariableBounds()`, `SNESSetFunctionDomainError()`, `SNESSetJacobianDomainError()`, `SNESVINEWTONRSLS`, `SNESVINEWTONSSLS`, `SNESSetType()`, `PETSC_NINFINITY`, `PETSC_INFINITY`

# External Links
$(_doc_external("SNES/SNESVIGetVariableBounds"))
"""
function SNESVIGetVariableBounds(petsclib::PetscLibType, snes::AbstractSNES, xl::AbstractPetscVec, xu::AbstractPetscVec)
    error("SNESVIGetVariableBounds: no generated method for these argument types")
end

@for_petsc function SNESVIGetVariableBounds(petsclib::$UnionPetscLib, snes::AbstractSNES, xl::AbstractPetscVec, xu::AbstractPetscVec )
	xl_ = Ref(xl.ptr)
	xu_ = Ref(xu.ptr)

    @chk ccall(
               (:SNESVIGetVariableBounds, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CVec}, Ptr{CVec}),
               snes, xl_, xu_,
              )

	xl.ptr = xl_[]
	xu.ptr = xu_[]

	return nothing
end 

"""
	SNESVISetComputeVariableBounds(petsclib::PetscLibType, snes::AbstractSNES, compute::external) 
Sets a function that is called to compute the bounds on variable for
(differential) variable inequalities.

Input Parameters:
- `snes`    - the `SNES` context
- `compute` - function that computes the bounds

Calling sequence of `compute`:
- `snes`   - the `SNES` context
- `lower`  - vector to hold lower bounds
- `higher` - vector to hold upper bounds

Level: advanced

See also: [](sec_vi), `SNES`, `SNESVISetVariableBounds()`, `DMSetVariableBounds()`, `SNESSetFunctionDomainError()`, `SNESSetJacobianDomainError()`, `SNESVINEWTONRSLS`, `SNESVINEWTONSSLS`,
`SNESSetType()`, `PETSC_NINFINITY`, `PETSC_INFINITY`

# External Links
$(_doc_external("SNES/SNESVISetComputeVariableBounds"))
"""
function SNESVISetComputeVariableBounds(petsclib::PetscLibType, snes::AbstractSNES, compute::external)
    error("SNESVISetComputeVariableBounds: no generated method for these argument types")
end

@for_petsc function SNESVISetComputeVariableBounds(petsclib::$UnionPetscLib, snes::AbstractSNES, compute::external )

    @chk ccall(
               (:SNESVISetComputeVariableBounds, $petsc_library),
               PetscErrorCode,
               (CSNES, external),
               snes, compute,
              )


	return nothing
end 

"""
	SNESVISetRedundancyCheck(petsclib::PetscLibType, snes::AbstractSNES, func::external, ctx::Ptr{Cvoid}) 
Provide a function to check for any redundancy in the VI active set

Logically Collective

Input Parameters:
- `snes` - the `SNESVINEWTONRSLS` context
- `func` - the function to check of redundancies
- `ctx`  - optional context used by the function

Calling sequence of func:
- `snes`      - the `SNES` context
- `is_act`    - the set of points in the active sets
- `is_redact` - output, the set of points in the non-redundant active set
- `ctx`       - optional context

Level: advanced

See also: `SNES`, `SNESVINEWTONRSLS`, `SNESVIGetInactiveSet()`, `DMSetVI()`

# External Links
$(_doc_external("SNES/SNESVISetRedundancyCheck"))
"""
function SNESVISetRedundancyCheck(petsclib::PetscLibType, snes::AbstractSNES, func::external, ctx::Ptr{Cvoid})
    error("SNESVISetRedundancyCheck: no generated method for these argument types")
end

@for_petsc function SNESVISetRedundancyCheck(petsclib::$UnionPetscLib, snes::AbstractSNES, func::external, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:SNESVISetRedundancyCheck, $petsc_library),
               PetscErrorCode,
               (CSNES, external, Ptr{Cvoid}),
               snes, func, ctx,
              )


	return nothing
end 

"""
	SNESVISetVariableBounds(petsclib::PetscLibType, snes::AbstractSNES, xl::AbstractPetscVec, xu::AbstractPetscVec) 
Sets the lower and upper bounds for the solution vector. `xl` <= x <= `xu`. This allows solving
(differential) variable inequalities.

Input Parameters:
- `snes` - the `SNES` context.
- `xl`   - lower bound.
- `xu`   - upper bound.

Level: advanced

See also: [](sec_vi), `SNES`, `SNESVIGetVariableBounds()`, `SNESVISetComputeVariableBounds()`, `SNESSetFunctionDomainError()`, `SNESSetJacobianDomainError()`, `SNESVINEWTONRSLS`, `SNESVINEWTONSSLS`, `SNESSetType()`, `PETSC_NINFINITY`, `PETSC_INFINITY`

# External Links
$(_doc_external("SNES/SNESVISetVariableBounds"))
"""
function SNESVISetVariableBounds(petsclib::PetscLibType, snes::AbstractSNES, xl::AbstractPetscVec, xu::AbstractPetscVec)
    error("SNESVISetVariableBounds: no generated method for these argument types")
end

@for_petsc function SNESVISetVariableBounds(petsclib::$UnionPetscLib, snes::AbstractSNES, xl::AbstractPetscVec, xu::AbstractPetscVec )

    @chk ccall(
               (:SNESVISetVariableBounds, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec),
               snes, xl, xu,
              )


	return nothing
end 

"""
	SNESView(petsclib::PetscLibType, snes::AbstractSNES, viewer::PetscViewer) 
Prints or visualizes the `SNES` data structure.

Collective

Input Parameters:
- `snes`   - the `SNES` context
- `viewer` - the `PetscViewer`

Options Database Key:
- `-snes_view` - Calls `SNESView()` at end of `SNESSolve()`

Level: beginner

See also: `SNES`, `SNESLoad()`, `SNESCreate()`, `PetscViewerASCIIOpen()`

# External Links
$(_doc_external("SNES/SNESView"))
"""
function SNESView(petsclib::PetscLibType, snes::AbstractSNES, viewer::PetscViewer)
    error("SNESView: no generated method for these argument types")
end

@for_petsc function SNESView(petsclib::$UnionPetscLib, snes::AbstractSNES, viewer::PetscViewer )

    @chk ccall(
               (:SNESView, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscViewer),
               snes, viewer,
              )


	return nothing
end 

"""
	SNESViewFromOptions(petsclib::PetscLibType, A::AbstractSNES, obj, name::String) 
View a `SNES` based on values in the options database

Collective

Input Parameters:
- `A`    - the `SNES` context
- `obj`  - Optional object that provides the options prefix for the checks
- `name` - command line option

Options Database Key:
- `-name [viewertype][:...]` - option name and values. See `PetscObjectViewFromOptions()` for the possible arguments

Level: intermediate

See also: `SNES`, `SNESView`, `PetscObjectViewFromOptions()`, `SNESCreate()`

# External Links
$(_doc_external("SNES/SNESViewFromOptions"))
"""
function SNESViewFromOptions(petsclib::PetscLibType, A::AbstractSNES, obj, name::String)
    error("SNESViewFromOptions: no generated method for these argument types")
end

@for_petsc function SNESViewFromOptions(petsclib::$UnionPetscLib, A::AbstractSNES, obj, name::String )

    @chk ccall(
               (:SNESViewFromOptions, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

