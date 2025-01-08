"""
	its = KSPGetIterationNumber(ksp::AbstractKSP{PetscLib})

Gets the current iteration number; if the `KSPSolve()` is complete, returns the number of iterations used.

Not Collective

Input Parameter:
===
- `ksp` - the iterative context

Output Parameter:
===
- `its` - number of iterations

Level: intermediate

Note:
During the ith iteration this returns i-1

-seealso: [](ch_ksp), `KSP`, `KSPGetResidualNorm()`, `KSPBuildResidual()`, `KSPGetTotalIterations()`

# External Links
$(_doc_external("KSP/KSPGetIterationNumber"))
"""
function KSPGetIterationNumber(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	its = [PetscInt(1)]

	LibPETSc.KSPGetIterationNumber(
		PetscLib,
		ksp,
		Ref(its,1),
	)

	return its[1]
end
 
 
"""
	rnorm = KSPGetResidualNorm(ksp::AbstractKSP{PetscLib})

Gets the last (possibly approximate and/or preconditioned) residual norm that has been computed.

Not Collective

Input Parameter:
===
- `ksp` - the iterative context

Output Parameter:
===
- `rnorm` - residual norm

Level: intermediate

Notes:
For some methods, such as `KSPGMRES`, the norm is not computed directly from the residual.

The type of norm used by the method can be controlled with `KSPSetNormType()`

Certain solvers, under certain conditions, may not compute the final residual norm in an iteration, in that case the previous norm is returned.

-seealso: [](ch_ksp), `KSP`, `KSPSetNormType()`, `KSPBuildResidual()`, `KSPNormType`

# External Links
$(_doc_external("KSP/KSPGetResidualNorm"))
"""
function KSPGetResidualNorm(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	PetscReal = PetscLib.PetscReal
	rnorm = [PetscReal(1)]

	LibPETSc.KSPGetResidualNorm(
		PetscLib,
		ksp,
		rnorm,
	)

	return rnorm[1]
end
 
 
"""
	type = KSPGetType(ksp::AbstractKSP{PetscLib})

Gets the `KSP` type as a string from the `KSP` object.

Not Collective

Input Parameter:
===
- `ksp` - Krylov context

Output Parameter:
===
- `type` - name of the `KSP` method

Level: intermediate

-seealso: [](ch_ksp), `KSPType`, `KSP`, `KSPSetType()`

# External Links
$(_doc_external("KSP/KSPGetType"))
"""
function KSPGetType(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	r_type = Ref{PETSc.CKSPType}()

	LibPETSc.KSPGetType(
		PetscLib,
		ksp,
		r_type,
	)


	type = unsafe_string(r_type[])
	return type
end

"""
	 KSPSetTolerances(ksp::AbstractKSP{PetscLib},rtol::AbstractFloat,abstol::AbstractFloat,dtol::AbstractFloat,maxits::Int)

Sets the relative, absolute, divergence, and maximum
iteration tolerances used by the default `KSP` convergence testers.

Logically Collective

Input Parameters:
===
- `ksp`    - the Krylov subspace context
- `rtol`   - the relative convergence tolerance, relative decrease in the (possibly preconditioned) residual norm
- `abstol` - the absolute convergence tolerance   absolute size of the (possibly preconditioned) residual norm
- `dtol`   - the divergence tolerance,   amount (possibly preconditioned) residual norm can increase before `KSPConvergedDefault()` concludes that the method is diverging
- `maxits` - maximum number of iterations to use

Options Database Keys:
===
- `-ksp_atol <abstol>`   - Sets `abstol`
- `-ksp_rtol <rtol>`     - Sets `rtol`
- `-ksp_divtol <dtol>`   - Sets `dtol`
- `-ksp_max_it <maxits>` - Sets `maxits`

Level: intermediate

Notes:
All parameters must be non-negative.

Use `PETSC_CURRENT` to retain the current value of any of the parameters. The deprecated `PETSC_DEFAULT` also retains the current value (though the name is confusing).

Use `PETSC_DETERMINE` to use the default value for the given `KSP`. The default value is the value when the object's type is set.

For `dtol` and `maxits` use `PETSC_UMLIMITED` to indicate there is no upper bound on these values

See `KSPConvergedDefault()` for details how these parameters are used in the default convergence test.  See also `KSPSetConvergenceTest()`
for setting user-defined stopping criteria.

Fortran Note:
Use `PETSC_CURRENT_INTEGER`, `PETSC_CURRENT_REAL`, `PETSC_DETERMINE_INTEGER`, or `PETSC_DETERMINE_REAL`

-seealso: [](ch_ksp), `KSPGetTolerances()`, `KSPConvergedDefault()`, `KSPSetConvergenceTest()`, `KSP`, `KSPSetMinimumIterations()`

# External Links
$(_doc_external("KSP/KSPSetTolerances"))
"""
function KSPSetTolerances(ksp::AbstractKSP{PetscLib},rtol::AbstractFloat,abstol::AbstractFloat,dtol::AbstractFloat,maxits::Int) where {PetscLib}

	LibPETSc.KSPSetTolerances(
		PetscLib,
		ksp,
		rtol,
		abstol,
		dtol,
		maxits,
	)

	return nothing
end
 
 

"""
	rtol,abstol,dtol,maxits = KSPGetTolerances(ksp::AbstractKSP{PetscLib})

Gets the relative, absolute, divergence, and maximum
iteration tolerances used by the default `KSP` convergence tests.

Not Collective

Input Parameter:
===
- `ksp` - the Krylov subspace context

Output Parameters:
===
- `rtol`   - the relative convergence tolerance
- `abstol` - the absolute convergence tolerance
- `dtol`   - the divergence tolerance
- `maxits` - maximum number of iterations

Level: intermediate

Note:
The user can specify `NULL` for any parameter that is not needed.

-seealso: [](ch_ksp), `KSPSetTolerances()`, `KSP`, `KSPSetMinimumIterations()`, `KSPGetMinimumIterations()`

# External Links
$(_doc_external("KSP/KSPGetTolerances"))
"""
function KSPGetTolerances(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	PetscReal = PetscLib.PetscReal
	rtol = [PetscReal(1)]
	PetscReal = PetscLib.PetscReal
	abstol = [PetscReal(1)]
	PetscReal = PetscLib.PetscReal
	dtol = [PetscReal(1)]
	maxits = [PetscInt(1)]

	LibPETSc.KSPGetTolerances(
		PetscLib,
		ksp,
		rtol,
		abstol,
		dtol,
		Ref(maxits,1),
	)

	return rtol[1],abstol[1],dtol[1],maxits[1]
end
 
"""
	its = KSPGetTotalIterations(ksp::AbstractKSP{PetscLib})

Gets the total number of iterations this `KSP` object has performed since was created, counted over all linear solves

Not Collective

Input Parameter:
===
- `ksp` - the iterative context

Output Parameter:
===
- `its` - total number of iterations

Level: intermediate

Note:
Use `KSPGetIterationNumber()` to get the count for the most recent solve only
If this is called within a `KSPSolve()` (such as in a `KSPMonitor` routine) then it does not include iterations within that current solve

-seealso: [](ch_ksp), `KSP`, `KSPBuildResidual()`, `KSPGetResidualNorm()`, `KSPGetIterationNumber()`

# External Links
$(_doc_external("KSP/KSPGetTotalIterations"))
"""
function KSPGetTotalIterations(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	its = [PetscInt(1)]

	LibPETSc.KSPGetTotalIterations(
		PetscLib,
		ksp,
		Ref(its,1),
	)

	return its[1]
end
 
 
"""
	 KSPSetOperators(ksp::AbstractKSP{PetscLib},Amat::AbstractMatrix,Pmat::AbstractMatrix)

Sets the matrix associated with the linear system
and a (possibly) different one from which the preconditioner will be built

Collective

Input Parameters:
===
- `ksp`  - the `KSP` context
- `Amat` - the matrix that defines the linear system
- `Pmat` - the matrix to be used in constructing the preconditioner, usually the same as Amat.

Level: beginner

Notes:
If you know the operator Amat has a null space you can use `MatSetNullSpace()` and `MatSetTransposeNullSpace()` to supply the null
space to `Amat` and the `KSP` solvers will automatically use that null space as needed during the solution process.

All future calls to `KSPSetOperators()` must use the same size matrices!

Passing a `NULL` for `Amat` or `Pmat` removes the matrix that is currently used.

If you wish to replace either `Amat` or `Pmat` but leave the other one untouched then
first call `KSPGetOperators()` to get the one you wish to keep, call `PetscObjectReference()`
on it and then pass it back in your call to `KSPSetOperators()`.

Developer Notes:
If the operators have NOT been set with `KSPSetOperators()` then the operators
are created in the `PC` and returned to the user. In this case, if both operators
mat and pmat are requested, two DIFFERENT operators will be returned. If
only one is requested both operators in the `PC` will be the same (i.e. as
if one had called `KSPSetOperators()` with the same argument for both `Mat`s).
The user must set the sizes of the returned matrices and their type etc just
as if the user created them with `MatCreate()`. For example,

-vb
KSPGetOperators(ksp/pc,&mat,NULL); is equivalent to
set size, type, etc of mat

MatCreate(comm,&mat);
KSP/PCSetOperators(ksp/pc,mat,mat);
PetscObjectDereference((PetscObject)mat);
set size, type, etc of mat

and

KSP/PCGetOperators(ksp/pc,&mat,&pmat); is equivalent to
set size, type, etc of mat and pmat

MatCreate(comm,&mat);
MatCreate(comm,&pmat);
KSP/PCSetOperators(ksp/pc,mat,pmat);
PetscObjectDereference((PetscObject)mat);
PetscObjectDereference((PetscObject)pmat);
set size, type, etc of mat and pmat
-ve

The rationale for this support is so that when creating a `TS`, `SNES`, or `KSP` the hierarchy
of underlying objects (i.e. `SNES`, `KSP`, `PC`, `Mat`) and their lifespans can be completely
managed by the top most level object (i.e. the `TS`, `SNES`, or `KSP`). Another way to look
at this is when you create a `SNES` you do not NEED to create a `KSP` and attach it to
the `SNES` object (the `SNES` object manages it for you). Similarly when you create a `KSP`
you do not need to attach a `PC` to it (the `KSP` object manages the `PC` object for you).
Thus, why should YOU have to create the `Mat` and attach it to the `SNES`/`KSP`/`PC`, when
it can be created for you?

-seealso: [](ch_ksp), `KSP`, `Mat`, `KSPSolve()`, `KSPGetPC()`, `PCGetOperators()`, `PCSetOperators()`, `KSPGetOperators()`, `KSPSetComputeOperators()`, `KSPSetComputeInitialGuess()`, `KSPSetComputeRHS()`

# External Links
$(_doc_external("KSP/KSPSetOperators"))
"""
function KSPSetOperators(ksp::AbstractKSP{PetscLib},Amat::AbstractMatrix,Pmat::AbstractMatrix) where {PetscLib}

	LibPETSc.KSPSetOperators(
		PetscLib,
		ksp,
		Amat,
		Pmat,
	)

	return nothing
end

"""
	v = KSPGetSolution(ksp::AbstractKSP{PetscLib})

Gets the location of the solution for the
linear system to be solved.

Not Collective

Input Parameter:
===
- `ksp` - iterative context obtained from `KSPCreate()`

Output Parameter:
===
- `v` - solution vector

Level: developer

Note:
If this is called during a `KSPSolve()` the vector's values may not represent the solution
to the linear system.

-seealso: [](ch_ksp), `KSPGetRhs()`, `KSPBuildSolution()`, `KSPSolve()`, `KSP`

# External Links
$(_doc_external("KSP/KSPGetSolution"))
"""
function KSPGetSolution(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	r_v = Ref{CVec}()

	LibPETSc.KSPGetSolution(
		PetscLib,
		ksp,
		r_v,
	)


	v = VecPtr(PetscLib, r_v[], false)
	return v
end
 

"""
	 KSPSetDM(ksp::AbstractKSP{PetscLib},dm::AbstractDM{PetscLib})

Sets the `DM` that may be used by some preconditioners and that may be used to construct the linear system

Logically Collective

Input Parameters:
===
- `ksp` - the `KSP`
- `dm`  - the `DM`, cannot be `NULL` to remove a previously set `DM`

Level: intermediate

Notes:
If this is used then the `KSP` will attempt to use the `DM` to create the matrix and use the routine set with
`DMKSPSetComputeOperators()`. Use `KSPSetDMActive`(ksp,`PETSC_FALSE`) to instead use the matrix you've provided with
`KSPSetOperators()`.

A `DM` can only be used for solving one problem at a time because information about the problem is stored on the `DM`,
even when not using interfaces like `DMKSPSetComputeOperators()`.  Use `DMClone()` to get a distinct `DM` when solving
different problems using the same function space.

-seealso: [](ch_ksp), `KSP`, `DM`, `KSPGetDM()`, `KSPSetDMActive()`, `KSPSetComputeOperators()`, `KSPSetComputeRHS()`, `KSPSetComputeInitialGuess()`, `DMKSPSetComputeOperators()`, `DMKSPSetComputeRHS()`, `DMKSPSetComputeInitialGuess()`

# External Links
$(_doc_external("DM/KSPSetDM"))
"""
function KSPSetDM(ksp::AbstractKSP{PetscLib},dm::AbstractDM{PetscLib}) where {PetscLib}

	LibPETSc.KSPSetDM(
		PetscLib,
		ksp,
		dm,
	)

	return nothing
end
 

"""
	 UNTESTED !!!
	dm = KSPGetDM(ksp::AbstractKSP{PetscLib})

Gets the `DM` that may be used by some preconditioners and that may be used to construct the linear system

Not Collective

Input Parameter:
===
- `ksp` - the `KSP`

Output Parameter:
===
- `dm` - the `DM`

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `DM`, `KSPSetDM()`, `KSPSetDMActive()`

# External Links
$(_doc_external("DM/KSPGetDM"))
"""
function KSPGetDM(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	petsclib = getlib(PetscLib)
	opts = Options(petsclib)
	dm = DM{PetscLib}(C_NULL, opts, petsclib.age)

	LibPETSc.KSPGetDM(
		PetscLib,
		ksp,
		dm,
	)

	return dm
end
