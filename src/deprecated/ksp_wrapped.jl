"""
	 UNTESTED !!!
	inksp = KSPCreate(comm::MPI_Comm)

Creates the `KSP` context.

Collective

Input Parameter:
===
- `comm` - MPI communicator

Output Parameter:
===
- `inksp` - location to put the `KSP` context

Level: beginner

Note:
The default `KSPType` is `KSPGMRES` with a restart of 30, using modified Gram-Schmidt orthogonalization.

-seealso: [](ch_ksp), `KSPSetUp()`, `KSPSolve()`, `KSPDestroy()`, `KSP`, `KSPGMRES`, `KSPType`

# External Links
$(_doc_external("DM/KSPCreate"))
"""
function KSPCreate(comm::MPI_Comm)

	LibPETSc.KSPCreate(
		PetscLib,
		comm,
		inksp,
	)

	return inksp
end
 
 
"""
	 UNTESTED !!!
	 KSPSetType(ksp::AbstractKSP{PetscLib},type::KSPType)

Builds the `KSP` data structure for a particular `KSPType`

Logically Collective

Input Parameters:
===
- `ksp`  - the Krylov space context
- `type` - a known method

Options Database Key:
===
- `-ksp_type  <method>` - Sets the method; use `-help` for a list  of available methods (for instance, cg or gmres)

Level: intermediate

Notes:
See "petsc/include/petscksp.h" for available methods (for instance, `KSPCG` or `KSPGMRES`).

Normally, it is best to use the `KSPSetFromOptions()` command and
then set the `KSP` type from the options database rather than by using
this routine.  Using the options database provides the user with
maximum flexibility in evaluating the many different Krylov methods.
The `KSPSetType()` routine is provided for those situations where it
is necessary to set the iterative solver independently of the command
line or options database.  This might be the case, for example, when
the choice of iterative solver changes during the execution of the
program, and the user's application is taking responsibility for
choosing the appropriate method.  In other words, this routine is
not for beginners.

Developer Note:
`KSPRegister()` is used to add Krylov types to `KSPList` from which they are accessed by `KSPSetType()`.

-seealso: [](ch_ksp), `PCSetType()`, `KSPType`, `KSPRegister()`, `KSPCreate()`, `KSP`

# External Links
$(_doc_external("DM/KSPSetType"))
"""
function KSPSetType(ksp::AbstractKSP{PetscLib},type::KSPType) where {PetscLib}

	LibPETSc.KSPSetType(
		PetscLib,
		ksp,
		type,
	)

	return nothing
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
$(_doc_external("DM/KSPGetType"))
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
	 UNTESTED !!!
	 KSPSetUp(ksp::AbstractKSP{PetscLib})

Sets up the internal data structures for the
later use of an iterative solver.

Collective

Input Parameter:
===
- `ksp` - iterative context obtained from `KSPCreate()`

Level: developer

-seealso: [](ch_ksp), `KSPCreate()`, `KSPSolve()`, `KSPDestroy()`, `KSP`

# External Links
$(_doc_external("DM/KSPSetUp"))
"""
function KSPSetUp(ksp::AbstractKSP{PetscLib}) where {PetscLib}

	LibPETSc.KSPSetUp(
		PetscLib,
		ksp,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 KSPSetUpOnBlocks(ksp::AbstractKSP{PetscLib})

Sets up the preconditioner for each block in
the block Jacobi, overlapping Schwarz, and fieldsplit methods.

Collective

Input Parameter:
===
- `ksp` - the `KSP` context

Level: advanced

Notes:
`KSPSetUpOnBlocks()` is a routine that the user can optionally call for
more precise profiling (via -log_view) of the setup phase for these
block preconditioners.  If the user does not call `KSPSetUpOnBlocks()`,
it will automatically be called from within `KSPSolve()`.

Calling `KSPSetUpOnBlocks()` is the same as calling `PCSetUpOnBlocks()`
on the PC context within the `KSP` context.

-seealso: [](ch_ksp), `PCSetUpOnBlocks()`, `KSPSetUp()`, `PCSetUp()`, `KSP`

# External Links
$(_doc_external("DM/KSPSetUpOnBlocks"))
"""
function KSPSetUpOnBlocks(ksp::AbstractKSP{PetscLib}) where {PetscLib}

	LibPETSc.KSPSetUpOnBlocks(
		PetscLib,
		ksp,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 KSPSolve(ksp::AbstractKSP{PetscLib},b::AbstractVector,x::AbstractVector)

Solves linear system.

Collective

Input Parameters:
===
- `ksp` - iterative context obtained from `KSPCreate()`
- `b`   - the right-hand side vector
- `x`   - the solution (this may be the same vector as `b`, then `b` will be overwritten with answer)

Options Database Keys:
===
- `-ksp_view_eigenvalues`                      - compute preconditioned operators eigenvalues
- `-ksp_view_eigenvalues_explicit`             - compute the eigenvalues by forming the dense operator and using LAPACK
- `-ksp_view_mat binary`                       - save matrix to the default binary viewer
- `-ksp_view_pmat binary`                      - save matrix used to build preconditioner to the default binary viewer
- `-ksp_view_rhs binary`                       - save right-hand side vector to the default binary viewer
- `-ksp_view_solution binary`                  - save computed solution vector to the default binary viewer
(can be read later with src/ksp/tutorials/ex10.c for testing solvers)
- `-ksp_view_mat_explicit`                     - for matrix-free operators, computes the matrix entries and views them
- `-ksp_view_preconditioned_operator_explicit` - computes the product of the preconditioner and matrix as an explicit matrix and views it
- `-ksp_converged_reason`                      - print reason for converged or diverged, also prints number of iterations
- `-ksp_view_final_residual`                   - print 2-norm of true linear system residual at the end of the solution process
- `-ksp_error_if_not_converged`                - stop the program as soon as an error is detected in a `KSPSolve()`
- `-ksp_view_pre`                              - print the ksp data structure before the system solution
- `-ksp_view`                                  - print the ksp data structure at the end of the system solution

Level: beginner

Notes:
If one uses `KSPSetDM()` then `x` or `b` need not be passed. Use `KSPGetSolution()` to access the solution in this case.

The operator is specified with `KSPSetOperators()`.

`KSPSolve()` will normally return without generating an error regardless of whether the linear system was solved or if constructing the preconditioner failed.
Call `KSPGetConvergedReason()` to determine if the solver converged or failed and why. The option -ksp_error_if_not_converged or function `KSPSetErrorIfNotConverged()`
will cause `KSPSolve()` to error as soon as an error occurs in the linear solver.  In inner `KSPSolve()` `KSP_DIVERGED_ITS` is not treated as an error because when using nested solvers
it may be fine that inner solvers in the preconditioner do not converge during the solution process.

The number of iterations can be obtained from `KSPGetIterationNumber()`.

If you provide a matrix that has a `MatSetNullSpace()` and `MatSetTransposeNullSpace()` this will use that information to solve singular systems
in the least squares sense with a norm minimizing solution.

A x = b   where b = b_p + b_t where b_t is not in the range of A (and hence by the fundamental theorem of linear algebra is in the nullspace(A'), see `MatSetNullSpace()`).

`KSP` first removes b_t producing the linear system  A x = b_p (which has multiple solutions) and solves this to find the ||x|| minimizing solution (and hence
it finds the solution x orthogonal to the nullspace(A). The algorithm is simply in each iteration of the Krylov method we remove the nullspace(A) from the search
direction thus the solution which is a linear combination of the search directions has no component in the nullspace(A).

We recommend always using `KSPGMRES` for such singular systems.
If nullspace(A) = nullspace(A') (note symmetric matrices always satisfy this property) then both left and right preconditioning will work
If nullspace(A) != nullspace(A') then left preconditioning will work but right preconditioning may not work (or it may).

Developer Notes:
The reason we cannot always solve  nullspace(A) != nullspace(A') systems with right preconditioning is because we need to remove at each iteration
the nullspace(AB) from the search direction. While we know the nullspace(A) the nullspace(AB) equals B^-1 times the nullspace(A) but except for trivial preconditioners
such as diagonal scaling we cannot apply the inverse of the preconditioner to a vector and thus cannot compute the nullspace(AB).

If using a direct method (e.g., via the `KSP` solver
`KSPPREONLY` and a preconditioner such as `PCLU` or `PCILU`,
then its=1.  See `KSPSetTolerances()` and `KSPConvergedDefault()`
for more details.

Understanding Convergence:
The routines `KSPMonitorSet()`, `KSPComputeEigenvalues()`, and
`KSPComputeEigenvaluesExplicitly()` provide information on additional
options to monitor convergence and print eigenvalue information.

-seealso: [](ch_ksp), `KSPCreate()`, `KSPSetUp()`, `KSPDestroy()`, `KSPSetTolerances()`, `KSPConvergedDefault()`,
`KSPSolveTranspose()`, `KSPGetIterationNumber()`, `MatNullSpaceCreate()`, `MatSetNullSpace()`, `MatSetTransposeNullSpace()`, `KSP`,
`KSPConvergedReasonView()`, `KSPCheckSolve()`, `KSPSetErrorIfNotConverged()`

# External Links
$(_doc_external("DM/KSPSolve"))
"""
function KSPSolve(ksp::AbstractKSP{PetscLib},b::AbstractVector,x::AbstractVector) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.KSPSolve(
		PetscLib,
		ksp,
		b,
		x,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 KSPSolveTranspose(ksp::AbstractKSP{PetscLib},b::AbstractVector,x::AbstractVector)

Solves a linear system with the transpose of the matrix,  A^T x = b.

Collective

Input Parameters:
===
- `ksp` - iterative context obtained from `KSPCreate()`
- `b`   - right-hand side vector
- `x`   - solution vector

Level: developer

Note:
For complex numbers this solve the non-Hermitian transpose system.

Developer Note:
We need to implement a `KSPSolveHermitianTranspose()`

-seealso: [](ch_ksp), `KSPCreate()`, `KSPSetUp()`, `KSPDestroy()`, `KSPSetTolerances()`, `KSPConvergedDefault()`,
`KSPSolve()`, `KSP`

# External Links
$(_doc_external("DM/KSPSolveTranspose"))
"""
function KSPSolveTranspose(ksp::AbstractKSP{PetscLib},b::AbstractVector,x::AbstractVector) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.KSPSolveTranspose(
		PetscLib,
		ksp,
		b,
		x,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 KSPSetUseExplicitTranspose(ksp::AbstractKSP{PetscLib},flg::PetscBool)

Determines the explicit transpose of the operator is formed in `KSPSolveTranspose()`. In some configurations (like GPUs) it may
be explicitly formed since the solve is much more efficient.

Logically Collective

Input Parameter:
===
- `ksp` - the `KSP` context

Output Parameter:
===
- `flg` - `PETSC_TRUE` to transpose the system in `KSPSolveTranspose()`, `PETSC_FALSE` to not transpose (default)

Level: advanced

-seealso: [](ch_ksp), `KSPSolveTranspose()`, `KSP`

# External Links
$(_doc_external("DM/KSPSetUseExplicitTranspose"))
"""
function KSPSetUseExplicitTranspose(ksp::AbstractKSP{PetscLib},flg::PetscBool) where {PetscLib}

	LibPETSc.KSPSetUseExplicitTranspose(
		PetscLib,
		ksp,
		flg,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 KSPMatSolve(ksp::AbstractKSP{PetscLib},B::AbstractMatrix,X::AbstractMatrix)

Solves a linear system with multiple right

Input Parameters:
===
- `ksp` - iterative context
- `B`   - block of right-hand sides

Output Parameter:
===
- `X` - block of solutions

Level: intermediate

Note:
This is a stripped-down version of `KSPSolve()`, which only handles `-ksp_view`, `-ksp_converged_reason`, `-ksp_converged_rate`, and `-ksp_view_final_residual`.

-seealso: [](ch_ksp), `KSPSolve()`, `MatMatSolve()`, `KSPMatSolveTranspose()`, `MATDENSE`, `KSPHPDDM`, `PCBJACOBI`, `PCASM`

# External Links
$(_doc_external("DM/KSPMatSolve"))
"""
function KSPMatSolve(ksp::AbstractKSP{PetscLib},B::AbstractMatrix,X::AbstractMatrix) where {PetscLib}

	LibPETSc.KSPMatSolve(
		PetscLib,
		ksp,
		B,
		X,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 KSPMatSolveTranspose(ksp::AbstractKSP{PetscLib},B::AbstractMatrix,X::AbstractMatrix)

Solves a linear system with the transposed matrix with multiple right
`B` and `X` must be different matrices and the transposed matrix cannot be assembled explicitly for the user.

Input Parameters:
===
- `ksp` - iterative context
- `B`   - block of right-hand sides

Output Parameter:
===
- `X` - block of solutions

Level: intermediate

Note:
This is a stripped-down version of `KSPSolveTranspose()`, which only handles `-ksp_view`, `-ksp_converged_reason`, `-ksp_converged_rate`, and `-ksp_view_final_residual`.

-seealso: [](ch_ksp), `KSPSolveTranspose()`, `MatMatTransposeSolve()`, `KSPMatSolve()`, `MATDENSE`, `KSPHPDDM`, `PCBJACOBI`, `PCASM`

# External Links
$(_doc_external("DM/KSPMatSolveTranspose"))
"""
function KSPMatSolveTranspose(ksp::AbstractKSP{PetscLib},B::AbstractMatrix,X::AbstractMatrix) where {PetscLib}

	LibPETSc.KSPMatSolveTranspose(
		PetscLib,
		ksp,
		B,
		X,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 KSPReset(ksp::AbstractKSP{PetscLib})

Resets a `KSP` context to the kspsetupcalled = 0 state and removes any allocated Vecs and Mats

Collective

Input Parameter:
===
- `ksp` - iterative context obtained from `KSPCreate()`

Level: beginner

-seealso: [](ch_ksp), `KSPCreate()`, `KSPSetUp()`, `KSPSolve()`, `KSP`

# External Links
$(_doc_external("DM/KSPReset"))
"""
function KSPReset(ksp::AbstractKSP{PetscLib}) where {PetscLib}

	LibPETSc.KSPReset(
		PetscLib,
		ksp,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 KSPResetViewers(ksp::AbstractKSP{PetscLib})

Resets all the viewers set from the options database during `KSPSetFromOptions()`

Collective

Input Parameter:
===
- `ksp` - iterative context obtained from `KSPCreate()`

Level: beginner

-seealso: [](ch_ksp), `KSPCreate()`, `KSPSetUp()`, `KSPSolve()`, `KSPSetFromOptions()`, `KSP`

# External Links
$(_doc_external("DM/KSPResetViewers"))
"""
function KSPResetViewers(ksp::AbstractKSP{PetscLib}) where {PetscLib}

	LibPETSc.KSPResetViewers(
		PetscLib,
		ksp,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	ksp = KSPDestroy()

Destroys a `KSP` context.

Collective

Input Parameter:
===
- `ksp` - iterative context obtained from `KSPCreate()`

Level: beginner

-seealso: [](ch_ksp), `KSPCreate()`, `KSPSetUp()`, `KSPSolve()`, `KSP`

# External Links
$(_doc_external("DM/KSPDestroy"))
"""
function KSPDestroy() 

	LibPETSc.KSPDestroy(
		PetscLib,
		ksp,
	)

	return ksp
end
 
 
"""
	 UNTESTED !!!
	 KSPSetReusePreconditioner(ksp::AbstractKSP{PetscLib},flag::PetscBool)

reuse the current preconditioner, do not construct a new one even if the operator changes

Collective

Input Parameters:
===
- `ksp`  - iterative context obtained from `KSPCreate()`
- `flag` - `PETSC_TRUE` to reuse the current preconditioner

Level: intermediate

-seealso: [](ch_ksp), `KSPCreate()`, `KSPSolve()`, `KSPDestroy()`, `PCSetReusePreconditioner()`, `KSP`

# External Links
$(_doc_external("DM/KSPSetReusePreconditioner"))
"""
function KSPSetReusePreconditioner(ksp::AbstractKSP{PetscLib},flag::PetscBool) where {PetscLib}

	LibPETSc.KSPSetReusePreconditioner(
		PetscLib,
		ksp,
		flag,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	flag = KSPGetReusePreconditioner(ksp::AbstractKSP{PetscLib})

Determines if the `KSP` reuses the current preconditioner even if the operator in the preconditioner has changed.

Collective

Input Parameter:
===
- `ksp` - iterative context obtained from `KSPCreate()`

Output Parameter:
===
- `flag` - the boolean flag

Level: intermediate

-seealso: [](ch_ksp), `KSPCreate()`, `KSPSolve()`, `KSPDestroy()`, `KSPSetReusePreconditioner()`, `KSP`

# External Links
$(_doc_external("DM/KSPGetReusePreconditioner"))
"""
function KSPGetReusePreconditioner(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	flag = Ref{PetscBool}()

	LibPETSc.KSPGetReusePreconditioner(
		PetscLib,
		ksp,
		flag,
	)

	return flag[] == PETSC_TRUE
end
 
 
"""
	 UNTESTED !!!
	 KSPSetSkipPCSetFromOptions(ksp::AbstractKSP{PetscLib},flag::PetscBool)

prevents `KSPSetFromOptions()` from calling `PCSetFromOptions()`. This is used if the same `PC` is shared by more than one `KSP` so its options are not resettable for each `KSP`

Collective

Input Parameters:
===
- `ksp`  - iterative context obtained from `KSPCreate()`
- `flag` - `PETSC_TRUE` to skip calling the `PCSetFromOptions()`

Level: intermediate

-seealso: [](ch_ksp), `KSPCreate()`, `KSPSolve()`, `KSPDestroy()`, `PCSetReusePreconditioner()`, `KSP`

# External Links
$(_doc_external("DM/KSPSetSkipPCSetFromOptions"))
"""
function KSPSetSkipPCSetFromOptions(ksp::AbstractKSP{PetscLib},flag::PetscBool) where {PetscLib}

	LibPETSc.KSPSetSkipPCSetFromOptions(
		PetscLib,
		ksp,
		flag,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	side = KSPGetPCSide(ksp::AbstractKSP{PetscLib})

Gets the preconditioning side.

Not Collective

Input Parameter:
===
- `ksp` - iterative context obtained from `KSPCreate()`

Output Parameter:
===
- `side` - the preconditioning side, where side is one of
-vb
PC_LEFT      - left preconditioning (default)
PC_RIGHT     - right preconditioning
PC_SYMMETRIC - symmetric preconditioning
-ve

Level: intermediate

-seealso: [](ch_ksp), `KSPSetPCSide()`, `KSP`

# External Links
$(_doc_external("DM/KSPGetPCSide"))
"""
function KSPGetPCSide(ksp::AbstractKSP{PetscLib}) where {PetscLib}

	LibPETSc.KSPGetPCSide(
		PetscLib,
		ksp,
		side,
	)

	return side
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
$(_doc_external("DM/KSPSetTolerances"))
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
$(_doc_external("DM/KSPGetTolerances"))
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
	 UNTESTED !!!
	 KSPSetMinimumIterations(ksp::AbstractKSP{PetscLib},minit::Int)

Sets the minimum number of iterations to use, regardless of the tolerances

Logically Collective

Input Parameters:
===
- `ksp`   - the Krylov subspace context
- `minit` - minimum number of iterations to use

Options Database Key:
===
- `-ksp_min_it <minits>` - Sets `minit`

Level: intermediate

Notes:
Use `KSPSetTolerances()` to set a variety of other tolerances

See `KSPConvergedDefault()` for details on how these parameters are used in the default convergence test. See also `KSPSetConvergenceTest()`
for setting user-defined stopping criteria.

-seealso: [](ch_ksp), `KSPGetTolerances()`, `KSPConvergedDefault()`, `KSPSetConvergenceTest()`, `KSP`, `KSPSetTolerances()`, `KSPGetMinimumIterations()`

# External Links
$(_doc_external("DM/KSPSetMinimumIterations"))
"""
function KSPSetMinimumIterations(ksp::AbstractKSP{PetscLib},minit::Int) where {PetscLib}

	LibPETSc.KSPSetMinimumIterations(
		PetscLib,
		ksp,
		minit,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	minit = KSPGetMinimumIterations(ksp::AbstractKSP{PetscLib})

Gets the minimum number of iterations to use, regardless of the tolerances, that was set with `KSPSetMinimumIterations()` or `

Not Collective

Input Parameter:
===
- `ksp` - the Krylov subspace context

Output Parameter:
===
- `minit` - minimum number of iterations to use

Level: intermediate

-seealso: [](ch_ksp), `KSPGetTolerances()`, `KSPConvergedDefault()`, `KSPSetConvergenceTest()`, `KSP`, `KSPSetTolerances()`, `KSPSetMinimumIterations()`

# External Links
$(_doc_external("DM/KSPGetMinimumIterations"))
"""
function KSPGetMinimumIterations(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	minit = [PetscInt(1)]

	LibPETSc.KSPGetMinimumIterations(
		PetscLib,
		ksp,
		Ref(minit,1),
	)

	return minit[1]
end
 
 
"""
	 UNTESTED !!!
	 KSPSetInitialGuessNonzero(ksp::AbstractKSP{PetscLib},flg::PetscBool)

Tells the iterative solver that the
initial guess is nonzero; otherwise `KSP` assumes the initial guess
is to be zero (and thus zeros it out before solving).

Logically Collective

Input Parameters:
===
- `ksp` - iterative context obtained from `KSPCreate()`
- `flg` - ``PETSC_TRUE`` indicates the guess is non-zero, `PETSC_FALSE` indicates the guess is zero

Options Database Key:
===
- `-ksp_initial_guess_nonzero <true,false>` - use nonzero initial guess

Level: beginner

Note:
If this is not called the X vector is zeroed in the call to `KSPSolve()`.

-seealso: [](ch_ksp), `KSPGetInitialGuessNonzero()`, `KSPGuessSetType()`, `KSPGuessType`, `KSP`

# External Links
$(_doc_external("DM/KSPSetInitialGuessNonzero"))
"""
function KSPSetInitialGuessNonzero(ksp::AbstractKSP{PetscLib},flg::PetscBool) where {PetscLib}

	LibPETSc.KSPSetInitialGuessNonzero(
		PetscLib,
		ksp,
		flg,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	flag = KSPGetInitialGuessNonzero(ksp::AbstractKSP{PetscLib})

Determines whether the `KSP` solver is using
a zero initial guess.

Not Collective

Input Parameter:
===
- `ksp` - iterative context obtained from `KSPCreate()`

Output Parameter:
===
- `flag` - `PETSC_TRUE` if guess is nonzero, else `PETSC_FALSE`

Level: intermediate

-seealso: [](ch_ksp), `KSPSetInitialGuessNonzero()`, `KSP`

# External Links
$(_doc_external("DM/KSPGetInitialGuessNonzero"))
"""
function KSPGetInitialGuessNonzero(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	flag = Ref{PetscBool}()

	LibPETSc.KSPGetInitialGuessNonzero(
		PetscLib,
		ksp,
		flag,
	)

	return flag[] == PETSC_TRUE
end
 
 
"""
	 UNTESTED !!!
	 KSPSetErrorIfNotConverged(ksp::AbstractKSP{PetscLib},flg::PetscBool)

Causes `KSPSolve()` to generate an error if the solver has not converged as soon as the error is detected.

Logically Collective

Input Parameters:
===
- `ksp` - iterative context obtained from `KSPCreate()`
- `flg` - `PETSC_TRUE` indicates you want the error generated

Options Database Key:
===
- `-ksp_error_if_not_converged <true,false>` - generate an error and stop the program

Level: intermediate

Notes:
Normally PETSc continues if a linear solver fails to converge, you can call `KSPGetConvergedReason()` after a `KSPSolve()`
to determine if it has converged.

A `KSP_DIVERGED_ITS` will not generate an error in a `KSPSolve()` inside a nested linear solver

-seealso: [](ch_ksp), `KSPGetErrorIfNotConverged()`, `KSP`

# External Links
$(_doc_external("DM/KSPSetErrorIfNotConverged"))
"""
function KSPSetErrorIfNotConverged(ksp::AbstractKSP{PetscLib},flg::PetscBool) where {PetscLib}

	LibPETSc.KSPSetErrorIfNotConverged(
		PetscLib,
		ksp,
		flg,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	flag = KSPGetErrorIfNotConverged(ksp::AbstractKSP{PetscLib})

Will `KSPSolve()` generate an error if the solver does not converge?

Not Collective

Input Parameter:
===
- `ksp` - iterative context obtained from KSPCreate()

Output Parameter:
===
- `flag` - `PETSC_TRUE` if it will generate an error, else `PETSC_FALSE`

Level: intermediate

-seealso: [](ch_ksp), `KSPSetErrorIfNotConverged()`, `KSP`

# External Links
$(_doc_external("DM/KSPGetErrorIfNotConverged"))
"""
function KSPGetErrorIfNotConverged(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	flag = Ref{PetscBool}()

	LibPETSc.KSPGetErrorIfNotConverged(
		PetscLib,
		ksp,
		flag,
	)

	return flag[] == PETSC_TRUE
end
 
 
"""
	 UNTESTED !!!
	 KSPSetComputeEigenvalues(ksp::AbstractKSP{PetscLib},flg::PetscBool)

Sets a flag so that the extreme eigenvalues
values will be calculated via a Lanczos or Arnoldi process as the linear
system is solved.

Logically Collective

Input Parameters:
===
- `ksp` - iterative context obtained from `KSPCreate()`
- `flg` - `PETSC_TRUE` or `PETSC_FALSE`

Level: advanced

Note:
Currently this option is not valid for all iterative methods.

-seealso: [](ch_ksp), `KSPComputeEigenvalues()`, `KSPComputeEigenvaluesExplicitly()`, `KSP`, `KSPSetComputeRitz()`

# External Links
$(_doc_external("DM/KSPSetComputeEigenvalues"))
"""
function KSPSetComputeEigenvalues(ksp::AbstractKSP{PetscLib},flg::PetscBool) where {PetscLib}

	LibPETSc.KSPSetComputeEigenvalues(
		PetscLib,
		ksp,
		flg,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 KSPSetComputeRitz(ksp::AbstractKSP{PetscLib},flg::PetscBool)

Sets a flag so that the Ritz or harmonic Ritz pairs
will be calculated via a Lanczos or Arnoldi process as the linear
system is solved.

Logically Collective

Input Parameters:
===
- `ksp` - iterative context obtained from `KSPCreate()`
- `flg` - `PETSC_TRUE` or `PETSC_FALSE`

Level: advanced

Note:
Currently this option is only valid for the `KSPGMRES` method.

-seealso: [](ch_ksp), `KSPComputeRitz()`, `KSP`

# External Links
$(_doc_external("DM/KSPSetComputeRitz"))
"""
function KSPSetComputeRitz(ksp::AbstractKSP{PetscLib},flg::PetscBool) where {PetscLib}

	LibPETSc.KSPSetComputeRitz(
		PetscLib,
		ksp,
		flg,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	flg = KSPGetComputeEigenvalues(ksp::AbstractKSP{PetscLib})

Gets the flag indicating that the extreme eigenvalues
values will be calculated via a Lanczos or Arnoldi process as the linear
system is solved.

Not Collective

Input Parameter:
===
- `ksp` - iterative context obtained from `KSPCreate()`

Output Parameter:
===
- `flg` - `PETSC_TRUE` or `PETSC_FALSE`

Level: advanced

Note:
Currently this option is not valid for all iterative methods.

-seealso: [](ch_ksp), `KSPComputeEigenvalues()`, `KSPComputeEigenvaluesExplicitly()`, `KSP`, `KSPSetComputeRitz()`

# External Links
$(_doc_external("DM/KSPGetComputeEigenvalues"))
"""
function KSPGetComputeEigenvalues(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	flg = Ref{PetscBool}()

	LibPETSc.KSPGetComputeEigenvalues(
		PetscLib,
		ksp,
		flg,
	)

	return flg[] == PETSC_TRUE
end
 
 
"""
	 UNTESTED !!!
	 KSPSetComputeSingularValues(ksp::AbstractKSP{PetscLib},flg::PetscBool)

Sets a flag so that the extreme singular
values will be calculated via a Lanczos or Arnoldi process as the linear
system is solved.

Logically Collective

Input Parameters:
===
- `ksp` - iterative context obtained from `KSPCreate()`
- `flg` - `PETSC_TRUE` or `PETSC_FALSE`

Options Database Key:
===
- `-ksp_monitor_singular_value` - Activates `KSPSetComputeSingularValues()`

Level: advanced

Notes:
Currently this option is not valid for all iterative methods.

Many users may just want to use the monitoring routine
`KSPMonitorSingularValue()` (which can be set with option -ksp_monitor_singular_value)
to print the singular values at each iteration of the linear solve.

-seealso: [](ch_ksp), `KSPComputeExtremeSingularValues()`, `KSPMonitorSingularValue()`, `KSP`, `KSPSetComputeRitz()`

# External Links
$(_doc_external("DM/KSPSetComputeSingularValues"))
"""
function KSPSetComputeSingularValues(ksp::AbstractKSP{PetscLib},flg::PetscBool) where {PetscLib}

	LibPETSc.KSPSetComputeSingularValues(
		PetscLib,
		ksp,
		flg,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	flg = KSPGetComputeSingularValues(ksp::AbstractKSP{PetscLib})

Gets the flag indicating whether the extreme singular
values will be calculated via a Lanczos or Arnoldi process as the linear
system is solved.

Not Collective

Input Parameter:
===
- `ksp` - iterative context obtained from `KSPCreate()`

Output Parameter:
===
- `flg` - `PETSC_TRUE` or `PETSC_FALSE`

Options Database Key:
===
- `-ksp_monitor_singular_value` - Activates `KSPSetComputeSingularValues()`

Level: advanced

Notes:
Currently this option is not valid for all iterative methods.

Many users may just want to use the monitoring routine
`KSPMonitorSingularValue()` (which can be set with option -ksp_monitor_singular_value)
to print the singular values at each iteration of the linear solve.

-seealso: [](ch_ksp), `KSPComputeExtremeSingularValues()`, `KSPMonitorSingularValue()`, `KSP`

# External Links
$(_doc_external("DM/KSPGetComputeSingularValues"))
"""
function KSPGetComputeSingularValues(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	flg = Ref{PetscBool}()

	LibPETSc.KSPGetComputeSingularValues(
		PetscLib,
		ksp,
		flg,
	)

	return flg[] == PETSC_TRUE
end
 
 
"""
	 UNTESTED !!!
	r = KSPGetRhs(ksp::AbstractKSP{PetscLib})

Gets the right
be solved.

Not Collective

Input Parameter:
===
- `ksp` - iterative context obtained from `KSPCreate()`

Output Parameter:
===
- `r` - right-hand-side vector

Level: developer

-seealso: [](ch_ksp), `KSPGetSolution()`, `KSPSolve()`, `KSP`

# External Links
$(_doc_external("DM/KSPGetRhs"))
"""
function KSPGetRhs(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	r_r = Ref{CVec}()

	LibPETSc.KSPGetRhs(
		PetscLib,
		ksp,
		r_r,
	)


	r = VecPtr(PetscLib, r_r[], false)
	return r
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
$(_doc_external("DM/KSPGetSolution"))
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
$(_doc_external("DM/KSPGetResidualNorm"))
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
$(_doc_external("DM/KSPGetIterationNumber"))
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
$(_doc_external("DM/KSPGetTotalIterations"))
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
	 UNTESTED !!!
	right,left = KSPCreateVecs(ksp::AbstractKSP{PetscLib},rightn::Int,leftn::Int)

Gets a number of work vectors suitably sized for the operator in the `KSP`

Collective

Input Parameters:
===
- `ksp`    - iterative context
- `rightn` - number of right work vectors to allocate
- `leftn`  - number of left work vectors to allocate

Output Parameters:
===
- `right` - the array of vectors created
- `left`  - the array of left vectors

Level: advanced

Notes:
The right vector has as many elements as the matrix has columns. The left
vector has as many elements as the matrix has rows, see `MatSetSizes()` for details on the layout of the vectors.

The vectors are new vectors that are not owned by the `KSP`, they should be destroyed with calls to `VecDestroyVecs()` when no longer needed.

Developer Note:
First tries to duplicate the rhs and solution vectors of the `KSP`, if they do not exist tries to get them from the matrix with `MatCreateVecs()`, if
that does not exist tries to get them from the `DM` (if it is provided) with `DMCreateGlobalVectors()`.

-seealso: [](ch_ksp), `MatCreateVecs()`, `VecDestroyVecs()`, `KSPSetWorkVecs()`

# External Links
$(_doc_external("DM/KSPCreateVecs"))
"""
function KSPCreateVecs(ksp::AbstractKSP{PetscLib},rightn::Int,leftn::Int) where {PetscLib}

	LibPETSc.KSPCreateVecs(
		PetscLib,
		ksp,
		rightn,
		leftn,
		right,
		left,
	)

	return right,left
end
 
 
"""
	 UNTESTED !!!
	pc = KSPGetPC(ksp::AbstractKSP{PetscLib})

Returns a pointer to the preconditioner context with the `KSP`

Not Collective

Input Parameter:
===
- `ksp` - iterative context obtained from `KSPCreate()`

Output Parameter:
===
- `pc` - preconditioner context

Level: developer

Note:
The `PC` is created if it does not already exist.

Developer Note:
Calls `KSPCheckPCMPI()` to check if the `KSP` is effected by `-mpi_linear_solver_server`

-seealso: [](ch_ksp), `KSPSetPC()`, `KSP`, `PC`

# External Links
$(_doc_external("DM/KSPGetPC"))
"""
function KSPGetPC(ksp::AbstractKSP{PetscLib}) where {PetscLib}

	LibPETSc.KSPGetPC(
		PetscLib,
		ksp,
		pc,
	)

	return pc
end
 
 
"""
	 UNTESTED !!!
	 KSPSetNestLevel(ksp::AbstractKSP{PetscLib},level::Int)

sets the amount of nesting the `KSP` has

Collective

Input Parameters:
===
- `ksp`   - the `KSP`
- `level` - the nest level

Level: developer

-seealso: [](ch_ksp), `KSPSetUp()`, `KSPSolve()`, `KSPDestroy()`, `KSP`, `KSPGMRES`, `KSPType`, `KSPGetNestLevel()`, `PCSetKSPNestLevel()`, `PCGetKSPNestLevel()`

# External Links
$(_doc_external("DM/KSPSetNestLevel"))
"""
function KSPSetNestLevel(ksp::AbstractKSP{PetscLib},level::Int) where {PetscLib}

	LibPETSc.KSPSetNestLevel(
		PetscLib,
		ksp,
		level,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	level = KSPGetNestLevel(ksp::AbstractKSP{PetscLib})

gets the amount of nesting the `KSP` has

Not Collective

Input Parameter:
===
- `ksp` - the `KSP`

Output Parameter:
===
- `level` - the nest level

Level: developer

-seealso: [](ch_ksp), `KSPSetUp()`, `KSPSolve()`, `KSPDestroy()`, `KSP`, `KSPGMRES`, `KSPType`, `KSPSetNestLevel()`, `PCSetKSPNestLevel()`, `PCGetKSPNestLevel()`

# External Links
$(_doc_external("DM/KSPGetNestLevel"))
"""
function KSPGetNestLevel(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	PetscInt = PetscLib.PetscInt
	level = [PetscInt(1)]

	LibPETSc.KSPGetNestLevel(
		PetscLib,
		ksp,
		Ref(level,1),
	)

	return level[1]
end
 
 
"""
	 UNTESTED !!!
	 KSPMonitorCancel(ksp::AbstractKSP{PetscLib})

Clears all monitors for a `KSP` object.

Logically Collective

Input Parameter:
===
- `ksp` - iterative context obtained from `KSPCreate()`

Options Database Key:
===
- `-ksp_monitor_cancel` - Cancels all monitors that have been hardwired into a code by calls to `KSPMonitorSet()`, but does not cancel those set via the options database.

Level: intermediate

-seealso: [](ch_ksp), `KSPMonitorResidual()`, `KSPMonitorSet()`, `KSP`

# External Links
$(_doc_external("DM/KSPMonitorCancel"))
"""
function KSPMonitorCancel(ksp::AbstractKSP{PetscLib}) where {PetscLib}

	LibPETSc.KSPMonitorCancel(
		PetscLib,
		ksp,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	ctx = KSPGetMonitorContext(ksp::AbstractKSP{PetscLib})

Gets the monitoring context, as set by `KSPMonitorSet()` for the FIRST monitor only.

Not Collective

Input Parameter:
===
- `ksp` - iterative context obtained from `KSPCreate()`

Output Parameter:
===
- `ctx` - monitoring context

Level: intermediate

-seealso: [](ch_ksp), `KSPMonitorResidual()`, `KSP`

# External Links
$(_doc_external("DM/KSPGetMonitorContext"))
"""
function KSPGetMonitorContext(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	PetscScalar = PetscLib.PetscScalar
	#TODO: your output is a vector; ensure that the size is correct!
	#It may involve: dims = DMStagGetGhostCorners(dm)[4:6]
	#dims = DMStagGetGhostCorners(dm)[4:6]    # dimensions including ghost values; set to 0 if not 2D/3D
	#dims = getindex(dims,findall(dims.>0))   # retrieve non-zero values
	#dmE  = DMStagGetEntriesPerElement(dm)    # dof per element
	dims = (X,)
	r_ctx = PETSc_RefPtr(dims, PetscScalar)

	LibPETSc.KSPGetMonitorContext(
		PetscLib,
		ksp,
		r_ctx,
	)

	ctx = PETSc_unsafe_wrap(r_ctx, dims; own=false)

	return ctx
end
 
 
"""
	 UNTESTED !!!
	na = KSPGetResidualHistory(ksp::AbstractKSP{PetscLib},a::Vector{PetscReal})

Gets the array used to hold the residual history and the number of residuals it contains.

Not Collective

Input Parameter:
===
- `ksp` - iterative context obtained from `KSPCreate()`

Output Parameters:
===
- `a`  - pointer to array to hold history (or `NULL`)
- `na` - number of used entries in a (or `NULL`). Note this has different meanings depending on the `reset` argument to `KSPSetResidualHistory()`

Level: advanced

Note:
This array is borrowed and should not be freed by the caller.

Can only be called after a `KSPSetResidualHistory()` otherwise `a` and `na` are set to `NULL` and zero

When `reset` was `PETSC_TRUE` since a residual is computed before the first iteration, the value of `na` is generally one more than the value
returned with `KSPGetIterationNumber()`.

Some Krylov methods may not compute the final residual norm when convergence is declared because the maximum number of iterations allowed has been reached.
In this situation, when `reset` was `PETSC_TRUE`, `na` will then equal the number of iterations reported with `KSPGetIterationNumber()`

Some Krylov methods (such as `KSPSTCG`), under certain circumstances, do not compute the final residual norm. In this situation, when `reset` was `PETSC_TRUE`,
`na` will then equal the number of iterations reported with `KSPGetIterationNumber()`

`KSPBCGSL` does not record the residual norms for the "subiterations" hence the results from `KSPGetResidualHistory()` and `KSPGetIterationNumber()` will be different

Fortran Note:
The Fortran version of this routine has a calling sequence
-vb
call KSPGetResidualHistory(KSP ksp, integer na, integer ierr)
-ve
note that you have passed a Fortran array into `KSPSetResidualHistory()` and you need
to access the residual values from this Fortran array you provided. Only the `na` (number of
residual norms currently held) is set.

-seealso: [](ch_ksp), `KSPSetResidualHistory()`, `KSP`, `KSPGetIterationNumber()`, `KSPSTCG`, `KSPBCGSL`

# External Links
$(_doc_external("DM/KSPGetResidualHistory"))
"""
function KSPGetResidualHistory(ksp::AbstractKSP{PetscLib},a::Vector{PetscReal}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	PetscInt = PetscLib.PetscInt
	na = [PetscInt(1)]

	LibPETSc.KSPGetResidualHistory(
		PetscLib,
		ksp,
		a,
		Ref(na,1),
	)

	return na[1]
end
 
 
"""
	 UNTESTED !!!
	 KSPSetResidualHistory(ksp::AbstractKSP{PetscLib},a::Vector{PetscReal},na::PetscCount,reset::PetscBool)

Sets the array used to hold the residual history.
If set, this array will contain the residual norms computed at each
iteration of the solver.

Not Collective

Input Parameters:
===
- `ksp`   - iterative context obtained from `KSPCreate()`
- `a`     - array to hold history
- `na`    - size of `a`
- `reset` - `PETSC_TRUE` indicates the history counter is reset to zero
for each new linear solve

Level: advanced

Notes:
If provided, `a` is NOT freed by PETSc so the user needs to keep track of it and destroy once the `KSP` object is destroyed.
If 'a' is `NULL` then space is allocated for the history. If 'na' `PETSC_DECIDE` or (deprecated) `PETSC_DEFAULT` then a
default array of length 10,000 is allocated.

If the array is not long enough then once the iterations is longer than the array length `KSPSolve()` stops recording the history

-seealso: [](ch_ksp), `KSPGetResidualHistory()`, `KSP`

# External Links
$(_doc_external("DM/KSPSetResidualHistory"))
"""
function KSPSetResidualHistory(ksp::AbstractKSP{PetscLib},a::Vector{PetscReal},na::PetscCount,reset::PetscBool) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.KSPSetResidualHistory(
		PetscLib,
		ksp,
		a,
		na,
		reset,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	na = KSPGetErrorHistory(ksp::AbstractKSP{PetscLib},a::Vector{PetscReal})

Gets the array used to hold the error history and the number of residuals it contains.

Not Collective

Input Parameter:
===
- `ksp` - iterative context obtained from `KSPCreate()`

Output Parameters:
===
- `a`  - pointer to array to hold history (or `NULL`)
- `na` - number of used entries in a (or `NULL`)

Level: advanced

Note:
This array is borrowed and should not be freed by the caller.
Can only be called after a `KSPSetErrorHistory()` otherwise `a` and `na` are set to `NULL` and zero

Fortran Note:
The Fortran version of this routine has a calling sequence
-vb
call KSPGetErrorHistory(KSP ksp, integer na, integer ierr)
-ve
note that you have passed a Fortran array into `KSPSetErrorHistory()` and you need
to access the residual values from this Fortran array you provided. Only the `na` (number of
residual norms currently held) is set.

-seealso: [](ch_ksp), `KSPSetErrorHistory()`, `KSPGetResidualHistory()`, `KSP`

# External Links
$(_doc_external("DM/KSPGetErrorHistory"))
"""
function KSPGetErrorHistory(ksp::AbstractKSP{PetscLib},a::Vector{PetscReal}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	PetscInt = PetscLib.PetscInt
	na = [PetscInt(1)]

	LibPETSc.KSPGetErrorHistory(
		PetscLib,
		ksp,
		a,
		Ref(na,1),
	)

	return na[1]
end
 
 
"""
	 UNTESTED !!!
	 KSPSetErrorHistory(ksp::AbstractKSP{PetscLib},a::Vector{PetscReal},na::PetscCount,reset::PetscBool)

Sets the array used to hold the error history. If set, this array will contain the error norms computed at each iteration of the solver.

Not Collective

Input Parameters:
===
- `ksp`   - iterative context obtained from `KSPCreate()`
- `a`     - array to hold history
- `na`    - size of `a`
- `reset` - `PETSC_TRUE` indicates the history counter is reset to zero for each new linear solve

Level: advanced

Notes:
If provided, `a` is NOT freed by PETSc so the user needs to keep track of it and destroy once the `KSP` object is destroyed.
If 'a' is `NULL` then space is allocated for the history. If 'na' is `PETSC_DECIDE` or (deprecated) `PETSC_DEFAULT` then a default array of length 1,0000 is allocated.

If the array is not long enough then once the iterations is longer than the array length `KSPSolve()` stops recording the history

-seealso: [](ch_ksp), `KSPGetErrorHistory()`, `KSPSetResidualHistory()`, `KSP`

# External Links
$(_doc_external("DM/KSPSetErrorHistory"))
"""
function KSPSetErrorHistory(ksp::AbstractKSP{PetscLib},a::Vector{PetscReal},na::PetscCount,reset::PetscBool) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.KSPSetErrorHistory(
		PetscLib,
		ksp,
		a,
		na,
		reset,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	V = KSPBuildSolutionDefault(ksp::AbstractKSP{PetscLib},v::AbstractVector)


# External Links
$(_doc_external("DM/KSPBuildSolutionDefault"))
"""
function KSPBuildSolutionDefault(ksp::AbstractKSP{PetscLib},v::AbstractVector) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	r_V = Ref{CVec}()

	LibPETSc.KSPBuildSolutionDefault(
		PetscLib,
		ksp,
		v,
		r_V,
	)


	V = VecPtr(PetscLib, r_V[], false)
	return V
end
 
 
"""
	 UNTESTED !!!
	V = KSPBuildResidualDefault(ksp::AbstractKSP{PetscLib},t::AbstractVector,v::AbstractVector)

Default code to compute the residual.

Collecive on ksp

Input Parameters:
===
- `ksp` - iterative context
- `t`   - pointer to temporary vector
- `v`   - pointer to user vector

Output Parameter:
===
- `V` - pointer to a vector containing the residual

Level: advanced

Note:
Some `KSP` methods such as `KSPGMRES` do not compute the explicit residual at each iteration, this routine takes the information
they have computed during the previous iterations and uses it to compute the explicit residual via the formula r = b - A*x.

Developer Note:
This is `PETSC_EXTERN` because it may be used by user written plugin `KSPType` implementations

-seealso: [](ch_ksp), `KSP`, `KSPBuildSolutionDefault()`

# External Links
$(_doc_external("DM/KSPBuildResidualDefault"))
"""
function KSPBuildResidualDefault(ksp::AbstractKSP{PetscLib},t::AbstractVector,v::AbstractVector) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	r_V = Ref{CVec}()

	LibPETSc.KSPBuildResidualDefault(
		PetscLib,
		ksp,
		t,
		v,
		r_V,
	)


	V = VecPtr(PetscLib, r_V[], false)
	return V
end
 
 
"""
	 UNTESTED !!!
	 KSPDestroyDefault(ksp::AbstractKSP{PetscLib})


# External Links
$(_doc_external("DM/KSPDestroyDefault"))
"""
function KSPDestroyDefault(ksp::AbstractKSP{PetscLib}) where {PetscLib}

	LibPETSc.KSPDestroyDefault(
		PetscLib,
		ksp,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 KSPSetWorkVecs(ksp::AbstractKSP{PetscLib},nw::Int)

Sets a number of work vectors into a `KSP` object

Collective

Input Parameters:
===
- `ksp` - iterative context
- `nw`  - number of work vectors to allocate

Level: developer

Developer Note:
This is `PETSC_EXTERN` because it may be used by user written plugin `KSPType` implementations

-seealso: [](ch_ksp), `KSP`, `KSPCreateVecs()`

# External Links
$(_doc_external("DM/KSPSetWorkVecs"))
"""
function KSPSetWorkVecs(ksp::AbstractKSP{PetscLib},nw::Int) where {PetscLib}

	LibPETSc.KSPSetWorkVecs(
		PetscLib,
		ksp,
		nw,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	V = KSPBuildSolution(ksp::AbstractKSP{PetscLib},v::AbstractVector)

Builds the approximate solution in a vector provided.

Collective

Input Parameter:
===
- `ksp` - iterative context obtained from `KSPCreate()`

Output Parameter:
===
Provide exactly one of
- `v` - location to stash solution, optional, otherwise pass `NULL`
- `V` - the solution is returned in this location. This vector is created internally. This vector should NOT be destroyed by the user with `VecDestroy()`.

Level: developer

Notes:
This routine can be used in one of two ways
-vb
KSPBuildSolution(ksp,NULL,&V);
or
KSPBuildSolution(ksp,v,NULL); or KSPBuildSolution(ksp,v,&v);
-ve
In the first case an internal vector is allocated to store the solution
(the user cannot destroy this vector). In the second case the solution
is generated in the vector that the user provides. Note that for certain
methods, such as `KSPCG`, the second case requires a copy of the solution,
while in the first case the call is essentially free since it simply
returns the vector where the solution already is stored. For some methods
like `KSPGMRES` during the solve this is a reasonably expensive operation and should only be
used if truly needed.

-seealso: [](ch_ksp), `KSPGetSolution()`, `KSPBuildResidual()`, `KSP`

# External Links
$(_doc_external("DM/KSPBuildSolution"))
"""
function KSPBuildSolution(ksp::AbstractKSP{PetscLib},v::AbstractVector) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	r_V = Ref{CVec}()

	LibPETSc.KSPBuildSolution(
		PetscLib,
		ksp,
		v,
		r_V,
	)


	V = VecPtr(PetscLib, r_V[], false)
	return V
end
 
 
"""
	 UNTESTED !!!
	V = KSPBuildResidual(ksp::AbstractKSP{PetscLib},t::AbstractVector,v::AbstractVector)

Builds the residual in a vector provided.

Collective

Input Parameter:
===
- `ksp` - iterative context obtained from `KSPCreate()`

Output Parameters:
===
- `v` - optional location to stash residual.  If `v` is not provided, then a location is generated.
- `t` - work vector.  If not provided then one is generated.
- `V` - the residual

Level: advanced

Note:
Regardless of whether or not `v` is provided, the residual is
returned in `V`.

-seealso: [](ch_ksp), `KSP`, `KSPBuildSolution()`

# External Links
$(_doc_external("DM/KSPBuildResidual"))
"""
function KSPBuildResidual(ksp::AbstractKSP{PetscLib},t::AbstractVector,v::AbstractVector) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	r_V = Ref{CVec}()

	LibPETSc.KSPBuildResidual(
		PetscLib,
		ksp,
		t,
		v,
		r_V,
	)


	V = VecPtr(PetscLib, r_V[], false)
	return V
end
 
 
"""
	 UNTESTED !!!
	emax,emin = KSPComputeExtremeSingularValues(ksp::AbstractKSP{PetscLib})

Computes the extreme singular values
for the preconditioned operator. Called after or during `KSPSolve()`.

Not Collective

Input Parameter:
===
- `ksp` - iterative context obtained from `KSPCreate()`

Output Parameters:
===
- `emax` - maximum estimated singular value
- `emin` - minimum estimated singular value

Options Database Key:
===
- `-ksp_view_singularvalues` - compute extreme singular values and print when `KSPSolve()` completes.

Level: advanced

Notes:
One must call `KSPSetComputeSingularValues()` before calling `KSPSetUp()`
(or use the option -ksp_view_eigenvalues) in order for this routine to work correctly.

Many users may just want to use the monitoring routine
`KSPMonitorSingularValue()` (which can be set with option -ksp_monitor_singular_value)
to print the extreme singular values at each iteration of the linear solve.

Estimates of the smallest singular value may be very inaccurate, especially if the Krylov method has not converged.
The largest singular value is usually accurate to within a few percent if the method has converged, but is still not
intended for eigenanalysis. Consider the excellent package `SLEPc` if accurate values are required.

Disable restarts if using KSPGMRES, otherwise this estimate will only be using those iterations after the last
restart. See `KSPGMRESSetRestart()` for more details.

-seealso: [](ch_ksp), `KSPSetComputeSingularValues()`, `KSPMonitorSingularValue()`, `KSPComputeEigenvalues()`, `KSP`, `KSPComputeRitz()`

# External Links
$(_doc_external("DM/KSPComputeExtremeSingularValues"))
"""
function KSPComputeExtremeSingularValues(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	PetscReal = PetscLib.PetscReal
	emax = [PetscReal(1)]
	PetscReal = PetscLib.PetscReal
	emin = [PetscReal(1)]

	LibPETSc.KSPComputeExtremeSingularValues(
		PetscLib,
		ksp,
		emax,
		emin,
	)

	return emax[1],emin[1]
end
 
 
"""
	 UNTESTED !!!
	neig = KSPComputeEigenvalues(ksp::AbstractKSP{PetscLib},n::Int,r::Vector{PetscReal},c::Vector{PetscReal})

Computes the extreme eigenvalues for the
preconditioned operator. Called after or during `KSPSolve()`.

Not Collective

Input Parameters:
===
- `ksp` - iterative context obtained from `KSPCreate()`
- `n`   - size of arrays `r` and `c`. The number of eigenvalues computed `neig` will, in
general, be less than this.

Output Parameters:
===
- `r`    - real part of computed eigenvalues, provided by user with a dimension of at least `n`
- `c`    - complex part of computed eigenvalues, provided by user with a dimension of at least `n`
- `neig` - actual number of eigenvalues computed (will be less than or equal to `n`)

Options Database Key:
===
- `-ksp_view_eigenvalues` - Prints eigenvalues to stdout

Level: advanced

Notes:
The number of eigenvalues estimated depends on the size of the Krylov space
generated during the `KSPSolve()` ; for example, with
`KSPCG` it corresponds to the number of CG iterations, for `KSPGMRES` it is the number
of GMRES iterations SINCE the last restart. Any extra space in `r` and `c`
will be ignored.

`KSPComputeEigenvalues()` does not usually provide accurate estimates; it is
intended only for assistance in understanding the convergence of iterative
methods, not for eigenanalysis. For accurate computation of eigenvalues we recommend using
the excellent package SLEPc.

One must call `KSPSetComputeEigenvalues()` before calling `KSPSetUp()`
in order for this routine to work correctly.

Many users may just want to use the monitoring routine
`KSPMonitorSingularValue()` (which can be set with option -ksp_monitor_singular_value)
to print the singular values at each iteration of the linear solve.

`KSPComputeRitz()` provides estimates for both the eigenvalues and their corresponding eigenvectors.

-seealso: [](ch_ksp), `KSPSetComputeEigenvalues()`, `KSPSetComputeSingularValues()`, `KSPMonitorSingularValue()`, `KSPComputeExtremeSingularValues()`, `KSP`, `KSPComputeRitz()`

# External Links
$(_doc_external("DM/KSPComputeEigenvalues"))
"""
function KSPComputeEigenvalues(ksp::AbstractKSP{PetscLib},n::Int,r::Vector{PetscReal},c::Vector{PetscReal}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	PetscInt = PetscLib.PetscInt
	neig = [PetscInt(1)]

	LibPETSc.KSPComputeEigenvalues(
		PetscLib,
		ksp,
		n,
		r,
		c,
		Ref(neig,1),
	)

	return neig[1]
end
 
 
"""
	 UNTESTED !!!
	 KSPComputeEigenvaluesExplicitly(ksp::AbstractKSP{PetscLib},nmax::Int,r::Vector{PetscReal},c::Vector{PetscReal})

Computes all of the eigenvalues of the
preconditioned operator using LAPACK.

Collective

Input Parameters:
===
- `ksp`  - iterative context obtained from `KSPCreate()`
- `nmax` - size of arrays `r` and `c`

Output Parameters:
===
- `r` - real part of computed eigenvalues, provided by user with a dimension at least of `n`
- `c` - complex part of computed eigenvalues, provided by user with a dimension at least of `n`

Level: advanced

Notes:
This approach is very slow but will generally provide accurate eigenvalue
estimates.  This routine explicitly forms a dense matrix representing
the preconditioned operator, and thus will run only for relatively small
problems, say `n` < 500.

Many users may just want to use the monitoring routine
`KSPMonitorSingularValue()` (which can be set with option -ksp_monitor_singular_value)
to print the singular values at each iteration of the linear solve.

The preconditioner operator, rhs vector, and solution vectors should be
set before this routine is called. i.e use `KSPSetOperators()`, `KSPSolve()`

-seealso: [](ch_ksp), `KSP`, `KSPComputeEigenvalues()`, `KSPMonitorSingularValue()`, `KSPComputeExtremeSingularValues()`, `KSPSetOperators()`, `KSPSolve()`

# External Links
$(_doc_external("DM/KSPComputeEigenvaluesExplicitly"))
"""
function KSPComputeEigenvaluesExplicitly(ksp::AbstractKSP{PetscLib},nmax::Int,r::Vector{PetscReal},c::Vector{PetscReal}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.KSPComputeEigenvaluesExplicitly(
		PetscLib,
		ksp,
		nmax,
		r,
		c,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 KSPSetFromOptions(ksp::AbstractKSP{PetscLib})

Sets `KSP` options from the options database.
This routine must be called before `KSPSetUp()` if the user is to be
allowed to set the Krylov type.

Collective

Input Parameter:
===
- `ksp` - the Krylov space context

Options Database Keys:
===
- `-ksp_rtol rtol`                                                          - relative tolerance used in default determination of convergence, i.e.
if residual norm decreases by this factor than convergence is declared
- `-ksp_atol abstol`                                                        - absolute tolerance used in default convergence test, i.e. if residual
norm is less than this then convergence is declared
- `-ksp_divtol tol`                                                         - if residual norm increases by this factor than divergence is declared
- `-ksp_max_it`                                                             - maximum number of linear iterations
- `-ksp_min_it`                                                             - minimum number of linear iterations to use, defaults to zero

- `-ksp_converged_use_initial_residual_norm`                                - see `KSPConvergedDefaultSetUIRNorm()`
- `-ksp_converged_use_min_initial_residual_norm`                            - see `KSPConvergedDefaultSetUMIRNorm()`
- `-ksp_converged_maxits`                                                   - see `KSPConvergedDefaultSetConvergedMaxits()`
- `-ksp_norm_type <none,preconditioned,unpreconditioned,natural>`           - see `KSPSetNormType()`
- `-ksp_check_norm_iteration it`                                            - do not compute residual norm until iteration number it (does compute at 0th iteration)
works only for `KSPBCGS`, `KSPIBCGS`, and `KSPCG`
- `-ksp_lag_norm`                                                           - compute the norm of the residual for the ith iteration on the i+1 iteration;
this means that one can use the norm of the residual for convergence test WITHOUT
an extra `MPI_Allreduce()` limiting global synchronizations.
This will require 1 more iteration of the solver than usual.
- `-ksp_guess_type`                                                         - Type of initial guess generator for repeated linear solves
- `-ksp_fischer_guess <model,size>`                                         - uses the Fischer initial guess generator for repeated linear solves
- `-ksp_constant_null_space`                                                - assume the operator (matrix) has the constant vector in its null space
- `-ksp_test_null_space`                                                    - tests the null space set with `MatSetNullSpace()` to see if it truly is a null space
- `-ksp_knoll`                                                              - compute initial guess by applying the preconditioner to the right-hand side
- `-ksp_monitor_cancel`                                                     - cancel all previous convergene monitor routines set
- `-ksp_monitor`                                                            - print residual norm at each iteration
- `-ksp_monitor draw::draw_lg`                                              - plot residual norm at each iteration
- `-ksp_monitor_true_residual`                                              - print true residual norm at each iteration
- `-all_ksp_monitor <optional filename>`                                    - print residual norm at each iteration for ALL KSP solves, regardless of their prefix. This is
useful for `PCFIELDSPLIT`, `PCMG`, etc that have inner solvers and
you wish to track the convergence of all the solvers
- `-ksp_monitor_solution [ascii binary or draw][:filename][:format option]` - plot solution at each iteration
- `-ksp_monitor_singular_value`                                             - monitor extreme singular values at each iteration
- `-ksp_converged_reason`                                                   - view the convergence state at the end of the solve
- `-ksp_use_explicittranspose`                                              - transpose the system explicitly in KSPSolveTranspose
- `-ksp_error_if_not_converged`                                             - stop the program as soon as an error is detected in a `KSPSolve()`, `KSP_DIVERGED_ITS`
is not treated as an error on inner solves
- `-ksp_converged_rate`                                                     - view the convergence rate at the end of the solve

Level: beginner

Note:
To see all options, run your program with the `-help` option or consult [](ch_ksp)

-seealso: [](ch_ksp), `KSP`, `KSPSetOptionsPrefix()`, `KSPResetFromOptions()`, `KSPSetUseFischerGuess()`

# External Links
$(_doc_external("DM/KSPSetFromOptions"))
"""
function KSPSetFromOptions(ksp::AbstractKSP{PetscLib}) where {PetscLib}

	LibPETSc.KSPSetFromOptions(
		PetscLib,
		ksp,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 KSPResetFromOptions(ksp::AbstractKSP{PetscLib})

Sets `KSP` parameters from user options ONLY if the `KSP` was previously set from options

Collective

Input Parameter:
===
- `ksp` - the `KSP` context

Level: advanced

-seealso: [](ch_ksp), `KSPSetFromOptions()`, `KSPSetOptionsPrefix()`

# External Links
$(_doc_external("DM/KSPResetFromOptions"))
"""
function KSPResetFromOptions(ksp::AbstractKSP{PetscLib}) where {PetscLib}

	LibPETSc.KSPResetFromOptions(
		PetscLib,
		ksp,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	ctx = KSPMonitorSetFromOptions(ksp::AbstractKSP{PetscLib},opt::Vector{Char},name::Vector{Char})

Sets a monitor function and viewer appropriate for the type indicated by the user in the options database

Collective

Input Parameters:
===
- `ksp`  - `KSP` object you wish to monitor
- `opt`  - the command line option for this monitor
- `name` - the monitor type one is seeking
- `ctx`  - An optional user context for the monitor, or `NULL`

Level: developer

-seealso: [](ch_ksp), `KSPMonitorRegister()`, `KSPMonitorSet()`, `PetscOptionsCreateViewer()`, `PetscOptionsGetReal()`, `PetscOptionsHasName()`, `PetscOptionsGetString()`,
`PetscOptionsGetIntArray()`, `PetscOptionsGetRealArray()`, `PetscOptionsBool()`
`PetscOptionsInt()`, `PetscOptionsString()`, `PetscOptionsReal()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`

# External Links
$(_doc_external("DM/KSPMonitorSetFromOptions"))
"""
function KSPMonitorSetFromOptions(ksp::AbstractKSP{PetscLib},opt::Vector{Char},name::Vector{Char}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 
	PetscScalar = PetscLib.PetscScalar
	#TODO: your output is a vector; ensure that the size is correct!
	#It may involve: dims = DMStagGetGhostCorners(dm)[4:6]
	#dims = DMStagGetGhostCorners(dm)[4:6]    # dimensions including ghost values; set to 0 if not 2D/3D
	#dims = getindex(dims,findall(dims.>0))   # retrieve non-zero values
	#dmE  = DMStagGetEntriesPerElement(dm)    # dof per element
	dims = (X,)
	r_ctx = PETSc_RefPtr(dims, PetscScalar)

	LibPETSc.KSPMonitorSetFromOptions(
		PetscLib,
		ksp,
		opt,
		name,
		r_ctx,
	)

	ctx = PETSc_unsafe_wrap(r_ctx, dims; own=false)

	return ctx
end
 

"""
	 UNTESTED !!!
	vf = KSPMonitorResidual(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat)

Print the (possibly preconditioned) residual norm at each iteration of an iterative solver.

Collective

Input Parameters:
===
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
===
- `-ksp_monitor` - Activates `KSPMonitorResidual()`

Level: intermediate

Note:
For some methods, such as `KSPGMRES`, the norm is not computed directly from the residual.

The type of norm used by the method can be controlled with `KSPSetNormType()`

This is not called directly by users, rather one calls `KSPMonitorSet()`, with this function as an argument, to cause the monitor
to be used during the `KSP` solve.

-seealso: [](ch_ksp), `KSP`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPMonitorResidualDraw()`, `KSPMonitorResidualDrawLG()`,
`KSPMonitorResidualRange()`, `KSPMonitorTrueResidualDraw()`, `KSPMonitorTrueResidualDrawLG()`, `KSPMonitorTrueResidualMax()`,
`KSPMonitorSingularValue()`, `KSPMonitorSolutionDrawLG()`, `KSPMonitorSolutionDraw()`, `KSPMonitorSolution()`,
`KSPMonitorErrorDrawLG()`, `KSPMonitorErrorDraw()`, `KSPMonitorError()`

# External Links
$(_doc_external("DM/KSPMonitorResidual"))
"""
function KSPMonitorResidual(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat) where {PetscLib}

	LibPETSc.KSPMonitorResidual(
		PetscLib,
		ksp,
		n,
		rnorm,
		vf,
	)

	return vf
end
 
 
"""
	 UNTESTED !!!
	vf = KSPMonitorResidualDraw(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat)

Plots the (possibly preconditioned) residual at each iteration of an iterative solver.

Collective

Input Parameters:
===
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
===
- `-ksp_monitor draw` - Activates `KSPMonitorResidualDraw()`

Level: intermediate

Note:
This is not called directly by users, rather one calls `KSPMonitorSet()`, with this function as an argument, to cause the monitor
to be used during the `KSP` solve.

-seealso: [](ch_ksp), `KSP`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPMonitorResidual()`, `KSPMonitorResidualDrawLG()`

# External Links
$(_doc_external("DM/KSPMonitorResidualDraw"))
"""
function KSPMonitorResidualDraw(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat) where {PetscLib}

	LibPETSc.KSPMonitorResidualDraw(
		PetscLib,
		ksp,
		n,
		rnorm,
		vf,
	)

	return vf
end
 
 
"""
	 UNTESTED !!!
	vf = KSPMonitorResidualDrawLG(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat)

Plots the (possibly preconditioned) residual norm at each iteration of an iterative solver.

Collective

Input Parameters:
===
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
===
- `-ksp_monitor draw::draw_lg` - Activates `KSPMonitorResidualDrawLG()`

Level: intermediate

Notes:
This is not called directly by users, rather one calls `KSPMonitorSet()`, with this function as an argument, to cause the monitor
to be used during the `KSP` solve.

Use `KSPMonitorResidualDrawLGCreate()` to create the context used with this monitor

-seealso: [](ch_ksp), `KSP`, `PETSCVIEWERDRAW`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPMonitorResidualDraw()`, `KSPMonitorResidual()`

# External Links
$(_doc_external("DM/KSPMonitorResidualDrawLG"))
"""
function KSPMonitorResidualDrawLG(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat) where {PetscLib}

	LibPETSc.KSPMonitorResidualDrawLG(
		PetscLib,
		ksp,
		n,
		rnorm,
		vf,
	)

	return vf
end
 

 
 
"""
	 UNTESTED !!!
	vf = KSPMonitorResidualShort(ksp::AbstractKSP{PetscLib},its::Int,fnorm::AbstractFloat)


# External Links
$(_doc_external("DM/KSPMonitorResidualShort"))
"""
function KSPMonitorResidualShort(ksp::AbstractKSP{PetscLib},its::Int,fnorm::AbstractFloat) where {PetscLib}

	LibPETSc.KSPMonitorResidualShort(
		PetscLib,
		ksp,
		its,
		fnorm,
		vf,
	)

	return vf
end
 
 
"""
	 UNTESTED !!!
	vf = KSPMonitorResidualRange(ksp::AbstractKSP{PetscLib},it::Int,rnorm::AbstractFloat)

Prints the percentage of residual elements that are more than 10 percent of the maximum value.

Collective

Input Parameters:
===
- `ksp`   - iterative context
- `it`    - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
===
- `-ksp_monitor_range` - Activates `KSPMonitorResidualRange()`

Level: intermediate

Note:
This is not called directly by users, rather one calls `KSPMonitorSet()`, with this function as an argument, to cause the monitor
to be used during the `KSP` solve.

-seealso: [](ch_ksp), `KSP`, `KSPMonitorSet()`, `KSPMonitorResidual()`

# External Links
$(_doc_external("DM/KSPMonitorResidualRange"))
"""
function KSPMonitorResidualRange(ksp::AbstractKSP{PetscLib},it::Int,rnorm::AbstractFloat) where {PetscLib}

	LibPETSc.KSPMonitorResidualRange(
		PetscLib,
		ksp,
		it,
		rnorm,
		vf,
	)

	return vf
end
 
 
"""
	 UNTESTED !!!
	vf = KSPMonitorTrueResidual(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat)

Prints the true residual norm, as well as the (possibly preconditioned) approximate residual norm, at each iteration of an iterative solver.

Collective

Input Parameters:
===
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
===
- `-ksp_monitor_true_residual` - Activates `KSPMonitorTrueResidual()`

Level: intermediate

Notes:
When using right preconditioning, these values are equivalent.

This is not called directly by users, rather one calls `KSPMonitorSet()`, with this function as an argument, to cause the monitor
to be used during the `KSP` solve.

-seealso: [](ch_ksp), `KSP`, `KSPMonitorSet()`, `KSPMonitorResidual()`, `KSPMonitorTrueResidualMaxNorm()`, `PetscViewerAndFormat`

# External Links
$(_doc_external("DM/KSPMonitorTrueResidual"))
"""
function KSPMonitorTrueResidual(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat) where {PetscLib}

	LibPETSc.KSPMonitorTrueResidual(
		PetscLib,
		ksp,
		n,
		rnorm,
		vf,
	)

	return vf
end
 
 
"""
	 UNTESTED !!!
	vf = KSPMonitorTrueResidualDraw(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat)

Plots the true residual at each iteration of an iterative solver.

Collective

Input Parameters:
===
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context of type `PETSCVIEWERDRAW`

Options Database Key:
===
- `-ksp_monitor_true_residual draw` - Activates `KSPMonitorResidualDraw()`

Level: intermediate

Note:
This is not called directly by users, rather one calls `KSPMonitorSet()`, with this function as an argument, to cause the monitor
to be used during the `KSP` solve.

-seealso: [](ch_ksp), `PETSCVIEWERDRAW`, `KSP`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPMonitorResidual()`,
`KSPMonitorTrueResidualDrawLG()`, `PetscViewerAndFormat`

# External Links
$(_doc_external("DM/KSPMonitorTrueResidualDraw"))
"""
function KSPMonitorTrueResidualDraw(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat) where {PetscLib}

	LibPETSc.KSPMonitorTrueResidualDraw(
		PetscLib,
		ksp,
		n,
		rnorm,
		vf,
	)

	return vf
end
 
 
"""
	 UNTESTED !!!
	vf = KSPMonitorTrueResidualDrawLG(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat)

Plots the true residual norm at each iteration of an iterative solver.

Collective

Input Parameters:
===
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
===
- `-ksp_monitor_true_residual draw::draw_lg` - Activates `KSPMonitorTrueResidualDrawLG()`

Level: intermediate

Notes:
This is not called directly by users, rather one calls `KSPMonitorSet()`, with this function as an argument, to cause the monitor
to be used during the `KSP` solve.

Call `KSPMonitorTrueResidualDrawLGCreate()` to create the context needed for this monitor

-seealso: [](ch_ksp), `PETSCVIEWERDRAW`, `KSP`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPMonitorTrueResidualDraw()`, `KSPMonitorResidual`,
`KSPMonitorTrueResidualDrawLGCreate()`

# External Links
$(_doc_external("DM/KSPMonitorTrueResidualDrawLG"))
"""
function KSPMonitorTrueResidualDrawLG(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat) where {PetscLib}

	LibPETSc.KSPMonitorTrueResidualDrawLG(
		PetscLib,
		ksp,
		n,
		rnorm,
		vf,
	)

	return vf
end
 
 

"""
	 UNTESTED !!!
	vf = KSPMonitorTrueResidualMax(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat)

Prints the true residual max norm at each iteration of an iterative solver.

Collective

Input Parameters:
===
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
===
- `-ksp_monitor_true_residual_max` - Activates `KSPMonitorTrueResidualMax()`

Level: intermediate

Note:
This is not called directly by users, rather one calls `KSPMonitorSet()`, with this function as an argument, to cause the monitor
to be used during the `KSP` solve.

-seealso: [](ch_ksp), `KSP`, `KSPMonitorSet()`, `KSPMonitorResidual()`, `KSPMonitorTrueResidualMaxNorm()`

# External Links
$(_doc_external("DM/KSPMonitorTrueResidualMax"))
"""
function KSPMonitorTrueResidualMax(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat) where {PetscLib}

	LibPETSc.KSPMonitorTrueResidualMax(
		PetscLib,
		ksp,
		n,
		rnorm,
		vf,
	)

	return vf
end
 
 
"""
	 UNTESTED !!!
	vf = KSPMonitorError(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat)

Prints the error norm, as well as the (possibly preconditioned) residual norm, at each iteration of an iterative solver.

Collective

Input Parameters:
===
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
===
- `-ksp_monitor_error` - Activates `KSPMonitorError()`

Level: intermediate

Note:
This is not called directly by users, rather one calls `KSPMonitorSet()`, with this function as an argument, to cause the monitor
to be used during the `KSP` solve.

-seealso: [](ch_ksp), `KSP`, `KSPMonitorSet()`, `KSPMonitorResidual()`, `KSPMonitorTrueResidualMaxNorm()`

# External Links
$(_doc_external("DM/KSPMonitorError"))
"""
function KSPMonitorError(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat) where {PetscLib}

	LibPETSc.KSPMonitorError(
		PetscLib,
		ksp,
		n,
		rnorm,
		vf,
	)

	return vf
end
 
 
"""
	 UNTESTED !!!
	vf = KSPMonitorErrorDraw(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat)

Plots the error at each iteration of an iterative solver.

Collective

Input Parameters:
===
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
===
- `-ksp_monitor_error draw` - Activates `KSPMonitorErrorDraw()`

Level: intermediate

Note:
This is not called directly by users, rather one calls `KSPMonitorSet()`, with this function as an argument, to cause the monitor
to be used during the `KSP` solve.

-seealso: [](ch_ksp), `PETSCVIEWERDRAW`, `KSP`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPMonitorErrorDrawLG()`

# External Links
$(_doc_external("DM/KSPMonitorErrorDraw"))
"""
function KSPMonitorErrorDraw(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat) where {PetscLib}

	LibPETSc.KSPMonitorErrorDraw(
		PetscLib,
		ksp,
		n,
		rnorm,
		vf,
	)

	return vf
end
 
 
"""
	 UNTESTED !!!
	vf = KSPMonitorErrorDrawLG(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat)

Plots the error and residual norm at each iteration of an iterative solver.

Collective

Input Parameters:
===
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
===
- `-ksp_monitor_error draw::draw_lg` - Activates `KSPMonitorTrueResidualDrawLG()`

Level: intermediate

Notes:
This is not called directly by users, rather one calls `KSPMonitorSet()`, with this function as an argument, to cause the monitor
to be used during the `KSP` solve.

Call `KSPMonitorErrorDrawLGCreate()` to create the context used with this monitor

-seealso: [](ch_ksp), `PETSCVIEWERDRAW`, `KSP`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPMonitorErrorDraw()`

# External Links
$(_doc_external("DM/KSPMonitorErrorDrawLG"))
"""
function KSPMonitorErrorDrawLG(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat) where {PetscLib}

	LibPETSc.KSPMonitorErrorDrawLG(
		PetscLib,
		ksp,
		n,
		rnorm,
		vf,
	)

	return vf
end
  
 
"""
	 UNTESTED !!!
	vf = KSPMonitorSolution(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat)

Print the solution norm at each iteration of an iterative solver.

Collective

Input Parameters:
===
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
===
- `-ksp_monitor_solution` - Activates `KSPMonitorSolution()`

Level: intermediate

Note:
This is not called directly by users, rather one calls `KSPMonitorSet()`, with this function as an argument, to cause the monitor
to be used during the `KSP` solve.

-seealso: [](ch_ksp), `KSPMonitorSet()`, `KSPMonitorTrueResidual()`

# External Links
$(_doc_external("DM/KSPMonitorSolution"))
"""
function KSPMonitorSolution(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat) where {PetscLib}

	LibPETSc.KSPMonitorSolution(
		PetscLib,
		ksp,
		n,
		rnorm,
		vf,
	)

	return vf
end
 
 
"""
	 UNTESTED !!!
	vf = KSPMonitorSolutionDraw(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat)

Plots the solution at each iteration of an iterative solver.

Collective

Input Parameters:
===
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
===
- `-ksp_monitor_solution draw` - Activates `KSPMonitorSolutionDraw()`

Level: intermediate

Note:
This is not called directly by users, rather one calls `KSPMonitorSet()`, with this function as an argument, to cause the monitor
to be used during the `KSP` solve.

-seealso: [](ch_ksp), `KSPMonitorSet()`, `KSPMonitorTrueResidual()`

# External Links
$(_doc_external("DM/KSPMonitorSolutionDraw"))
"""
function KSPMonitorSolutionDraw(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat) where {PetscLib}

	LibPETSc.KSPMonitorSolutionDraw(
		PetscLib,
		ksp,
		n,
		rnorm,
		vf,
	)

	return vf
end
 
 
"""
	 UNTESTED !!!
	vf = KSPMonitorSolutionDrawLG(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat)

Plots the solution norm at each iteration of an iterative solver.

Collective

Input Parameters:
===
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
===
- `-ksp_monitor_solution draw::draw_lg` - Activates `KSPMonitorSolutionDrawLG()`

Level: intermediate

Notes:
This is not called directly by users, rather one calls `KSPMonitorSet()`, with this function as an argument, to cause the monitor
to be used during the `KSP` solve.

Call `KSPMonitorSolutionDrawLGCreate()` to create the context needed with this monitor

-seealso: [](ch_ksp), `KSP`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPMonitorSolutionDrawLGCreate()`

# External Links
$(_doc_external("DM/KSPMonitorSolutionDrawLG"))
"""
function KSPMonitorSolutionDrawLG(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat) where {PetscLib}

	LibPETSc.KSPMonitorSolutionDrawLG(
		PetscLib,
		ksp,
		n,
		rnorm,
		vf,
	)

	return vf
end
 
 
 
"""
	 UNTESTED !!!
	vf = KSPMonitorSingularValue(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat)

Prints the two norm of the true residual and estimation of the extreme singular values of the preconditioned problem at each iteration.

Logically Collective

Input Parameters:
===
- `ksp`   - the iterative context
- `n`     - the iteration
- `rnorm` - the two norm of the residual
- `vf`    - The viewer context

Options Database Key:
===
- `-ksp_monitor_singular_value` - Activates `KSPMonitorSingularValue()`

Level: intermediate

Notes:
The `KSPCG` solver uses the Lanczos technique for eigenvalue computation,
while `KSPGMRES` uses the Arnoldi technique; other iterative methods do
not currently compute singular values.

This is not called directly by users, rather one calls `KSPMonitorSet()`, with this function as an argument, to cause the monitor
to be used during the `KSP` solve.

Call `KSPMonitorSingularValueCreate()` to create the context needed by this monitor

-seealso: [](ch_ksp), `KSP`, `KSPMonitorSet()`, `KSPComputeExtremeSingularValues()`, `KSPMonitorSingularValueCreate()`

# External Links
$(_doc_external("DM/KSPMonitorSingularValue"))
"""
function KSPMonitorSingularValue(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat) where {PetscLib}

	LibPETSc.KSPMonitorSingularValue(
		PetscLib,
		ksp,
		n,
		rnorm,
		vf,
	)

	return vf
end
 
 
 
"""
	 UNTESTED !!!
	ctx = KSPMonitorDynamicTolerance(ksp::AbstractKSP{PetscLib},its::Int,fnorm::AbstractFloat)

A monitor that changes the inner tolerance of nested preconditioners in every outer iteration in an adaptive way.

Collective

Input Parameters:
===
- `ksp`   - iterative context
- `its`   - iteration number (not used)
- `fnorm` - the current residual norm
- `ctx`   - context used by monitor

Options Database Key:
===
- `-sub_ksp_dynamic_tolerance <coef>` - coefficient of dynamic tolerance for inner solver, default is 1.0

Level: advanced

Notes:
Applies for `PCKSP`, `PCBJACOBI`, and `PCDEFLATION` preconditioners

This may be useful for a flexible preconditioned Krylov method, such as `KSPFGMRES`, [](sec_flexibleksp) to
control the accuracy of the inner solves needed to guarantee convergence of the outer iterations.

This is not called directly by users, rather one calls `KSPMonitorSet()`, with this function as an argument, to cause the monitor
to be used during the `KSP` solve.

Use `KSPMonitorDynamicToleranceCreate()` and `KSPMonitorDynamicToleranceSetCoefficient()` to create the context needed by this
monitor function.

Pass the context and `KSPMonitorDynamicToleranceDestroy()` to `KSPMonitorSet()`

-seealso: [](sec_flexibleksp), `KSP`, `KSPMonitorDynamicToleranceCreate()`, `KSPMonitorDynamicToleranceDestroy()`, `KSPMonitorDynamicToleranceSetCoefficient()`

# External Links
$(_doc_external("DM/KSPMonitorDynamicTolerance"))
"""
function KSPMonitorDynamicTolerance(ksp::AbstractKSP{PetscLib},its::Int,fnorm::AbstractFloat) where {PetscLib}
	PetscScalar = PetscLib.PetscScalar
	#TODO: your output is a vector; ensure that the size is correct!
	#It may involve: dims = DMStagGetGhostCorners(dm)[4:6]
	#dims = DMStagGetGhostCorners(dm)[4:6]    # dimensions including ghost values; set to 0 if not 2D/3D
	#dims = getindex(dims,findall(dims.>0))   # retrieve non-zero values
	#dmE  = DMStagGetEntriesPerElement(dm)    # dof per element
	dims = (X,)
	r_ctx = PETSc_RefPtr(dims, PetscScalar)

	LibPETSc.KSPMonitorDynamicTolerance(
		PetscLib,
		ksp,
		its,
		fnorm,
		r_ctx,
	)

	ctx = PETSc_unsafe_wrap(r_ctx, dims; own=false)

	return ctx
end
 
 

 
"""
	 UNTESTED !!!
	 KSPUnwindPreconditioner(ksp::AbstractKSP{PetscLib},vsoln::AbstractVector,vt1::AbstractVector)

Unwinds the preconditioning in the solution. That is,
takes solution to the preconditioned problem and gets the solution to the
original problem from it.

Collective

Input Parameters:
===
- `ksp`   - iterative context
- `vsoln` - solution vector
- `vt1`   - temporary work vector

Output Parameter:
===
- `vsoln` - contains solution on output

Level: advanced

Note:
If preconditioning either symmetrically or on the right, this routine solves
for the correction to the unpreconditioned problem.  If preconditioning on
the left, nothing is done.

-seealso: [](ch_ksp), `KSP`, `KSPSetPCSide()`

# External Links
$(_doc_external("DM/KSPUnwindPreconditioner"))
"""
function KSPUnwindPreconditioner(ksp::AbstractKSP{PetscLib},vsoln::AbstractVector,vt1::AbstractVector) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.KSPUnwindPreconditioner(
		PetscLib,
		ksp,
		vsoln,
		vt1,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 KSPInitialResidual(ksp::AbstractKSP{PetscLib},vsoln::AbstractVector,vt1::AbstractVector,vt2::AbstractVector,vres::AbstractVector,vb::AbstractVector)

Computes the residual. Either b
preconditioning or C*(b - A*x) with left preconditioning; the latter
residual is often called the "preconditioned residual".

Collective

Input Parameters:
===
- `ksp`   - the `KSP` solver object
- `vsoln` - solution to use in computing residual
- `vt1`   - temporary work vector
- `vt2`   - temporary work vector
- `vb`    - right-hand-side vector

Output Parameter:
===
- `vres` - calculated residual

Level: developer

Note:
This routine assumes that an iterative method, designed for  A x = b 
will be used with a preconditioner, C, such that the actual problem is either
-vb
AC u = b (right preconditioning) or
CA x = Cb (left preconditioning).
-ve
This means that the calculated residual will be scaled and/or preconditioned;
the true residual  b-Ax 
is returned in the `vt2` temporary work vector.

-seealso: [](ch_ksp), `KSP`, `KSPSolve()`, `KSPMonitor()`

# External Links
$(_doc_external("DM/KSPInitialResidual"))
"""
function KSPInitialResidual(ksp::AbstractKSP{PetscLib},vsoln::AbstractVector,vt1::AbstractVector,vt2::AbstractVector,vres::AbstractVector,vb::AbstractVector) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.KSPInitialResidual(
		PetscLib,
		ksp,
		vsoln,
		vt1,
		vt2,
		vres,
		vb,
	)

	return nothing
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
$(_doc_external("DM/KSPSetOperators"))
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
	 UNTESTED !!!
	Amat,Pmat = KSPGetOperators(ksp::AbstractKSP{PetscLib})

Gets the matrix associated with the linear system
and a (possibly) different one used to construct the preconditioner.

Collective

Input Parameter:
===
- `ksp` - the `KSP` context

Output Parameters:
===
- `Amat` - the matrix that defines the linear system
- `Pmat` - the matrix to be used in constructing the preconditioner, usually the same as `Amat`.

Level: intermediate

Note:
DOES NOT increase the reference counts of the matrix, so you should NOT destroy them.

-seealso: [](ch_ksp), `KSP`, `KSPSolve()`, `KSPGetPC()`, `PCGetOperators()`, `PCSetOperators()`, `KSPSetOperators()`, `KSPGetOperatorsSet()`

# External Links
$(_doc_external("DM/KSPGetOperators"))
"""
function KSPGetOperators(ksp::AbstractKSP{PetscLib}) where {PetscLib}

	LibPETSc.KSPGetOperators(
		PetscLib,
		ksp,
		Amat,
		Pmat,
	)

	return Amat,Pmat
end
 
 
"""
	 UNTESTED !!!
	mat,pmat = KSPGetOperatorsSet(ksp::AbstractKSP{PetscLib})

Determines if the matrix associated with the linear system and
possibly a different one associated with the preconditioner have been set in the `KSP`.

Not Collective, though the results on all processes should be the same

Input Parameter:
===
- `ksp` - the `KSP` context

Output Parameters:
===
- `mat`  - the matrix associated with the linear system was set
- `pmat` - matrix associated with the preconditioner was set, usually the same as `mat`

Level: intermediate

Note:
This routine exists because if you call `KSPGetOperators()` on a `KSP` that does not yet have operators they are
automatically created in the call.

-seealso: [](ch_ksp), `KSP`, `PCSetOperators()`, `KSPGetOperators()`, `KSPSetOperators()`, `PCGetOperators()`, `PCGetOperatorsSet()`

# External Links
$(_doc_external("DM/KSPGetOperatorsSet"))
"""
function KSPGetOperatorsSet(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	mat = Ref{PetscBool}()
	pmat = Ref{PetscBool}()

	LibPETSc.KSPGetOperatorsSet(
		PetscLib,
		ksp,
		mat,
		pmat,
	)

	return mat[] == PETSC_TRUE,pmat[] == PETSC_TRUE
end
 
 
"""
	 UNTESTED !!!
	 KSPAppendOptionsPrefix(ksp::AbstractKSP{PetscLib},prefix::Vector{Char})

Appends to the prefix used for searching for all
`KSP` options in the database.

Logically Collective

Input Parameters:
===
- `ksp`    - the Krylov context
- `prefix` - the prefix string to prepend to all `KSP` option requests

Level: intermediate

Note:
A hyphen (-) must NOT be given at the beginning of the prefix name.
The first character of all runtime options is AUTOMATICALLY the hyphen.

-seealso: [](ch_ksp), `KSP`, `KSPSetOptionsPrefix()`, `KSPGetOptionsPrefix()`, `KSPSetFromOptions()`

# External Links
$(_doc_external("DM/KSPAppendOptionsPrefix"))
"""
function KSPAppendOptionsPrefix(ksp::AbstractKSP{PetscLib},prefix::Vector{Char}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.KSPAppendOptionsPrefix(
		PetscLib,
		ksp,
		prefix,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 KSPGetOptionsPrefix(ksp::AbstractKSP{PetscLib},prefix::Vector{Char})

Gets the prefix used for searching for all
`KSP` options in the database.

Not Collective

Input Parameter:
===
- `ksp` - the Krylov context

Output Parameter:
===
- `prefix` - pointer to the prefix string used is returned

Level: advanced

Fortran Note:
Pass in a string 'prefix' of
sufficient length to hold the prefix.

-seealso: [](ch_ksp), `KSP`, `KSPSetFromOptions()`, `KSPSetOptionsPrefix()`, `KSPAppendOptionsPrefix()`

# External Links
$(_doc_external("DM/KSPGetOptionsPrefix"))
"""
function KSPGetOptionsPrefix(ksp::AbstractKSP{PetscLib},prefix::Vector{Char}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.KSPGetOptionsPrefix(
		PetscLib,
		ksp,
		prefix,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 KSPSetDiagonalScale(ksp::AbstractKSP{PetscLib},scale::PetscBool)

Tells `KSP` to symmetrically diagonally scale the system
before solving. This actually CHANGES the matrix (and right-hand side).

Logically Collective

Input Parameters:
===
- `ksp`   - the `KSP` context
- `scale` - `PETSC_TRUE` or `PETSC_FALSE`

Options Database Keys:
===
- `-ksp_diagonal_scale`     - perform a diagonal scaling before the solve
- `-ksp_diagonal_scale_fix` - scale the matrix back AFTER the solve

Level: advanced

Notes:
Scales the matrix by  D^{-1/2}  A  D^{-1/2}  [D^{1/2} x ] = D^{-1/2} b 
where D_{ii} is 1/abs(A_{ii})  unless A_{ii} is zero and then it is 1.

BE CAREFUL with this routine: it actually scales the matrix and right
hand side that define the system. After the system is solved the matrix
and right-hand side remain scaled unless you use `KSPSetDiagonalScaleFix()`

This should NOT be used within the `SNES` solves if you are using a line
search.

If you use this with the `PCType` `PCEISENSTAT` preconditioner than you can
use the `PCEisenstatSetNoDiagonalScaling()` option, or `-pc_eisenstat_no_diagonal_scaling`
to save some unneeded, redundant flops.

-seealso: [](ch_ksp), `KSPGetDiagonalScale()`, `KSPSetDiagonalScaleFix()`, `KSP`

# External Links
$(_doc_external("DM/KSPSetDiagonalScale"))
"""
function KSPSetDiagonalScale(ksp::AbstractKSP{PetscLib},scale::PetscBool) where {PetscLib}

	LibPETSc.KSPSetDiagonalScale(
		PetscLib,
		ksp,
		scale,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	scale = KSPGetDiagonalScale(ksp::AbstractKSP{PetscLib})

Checks if `KSP` solver scales the matrix and right

Not Collective

Input Parameter:
===
- `ksp` - the `KSP` context

Output Parameter:
===
- `scale` - `PETSC_TRUE` or `PETSC_FALSE`

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPSetDiagonalScale()`, `KSPSetDiagonalScaleFix()`

# External Links
$(_doc_external("DM/KSPGetDiagonalScale"))
"""
function KSPGetDiagonalScale(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	scale = Ref{PetscBool}()

	LibPETSc.KSPGetDiagonalScale(
		PetscLib,
		ksp,
		scale,
	)

	return scale[] == PETSC_TRUE
end
 
 
"""
	 UNTESTED !!!
	 KSPSetDiagonalScaleFix(ksp::AbstractKSP{PetscLib},fix::PetscBool)

Tells `KSP` to diagonally scale the system back after solving.

Logically Collective

Input Parameters:
===
- `ksp` - the `KSP` context
- `fix` - `PETSC_TRUE` to scale back after the system solve, `PETSC_FALSE` to not
rescale (default)

Level: intermediate

Notes:
Must be called after `KSPSetDiagonalScale()`

Using this will slow things down, because it rescales the matrix before and
after each linear solve. This is intended mainly for testing to allow one
to easily get back the original system to make sure the solution computed is
accurate enough.

-seealso: [](ch_ksp), `KSPGetDiagonalScale()`, `KSPSetDiagonalScale()`, `KSPGetDiagonalScaleFix()`, `KSP`

# External Links
$(_doc_external("DM/KSPSetDiagonalScaleFix"))
"""
function KSPSetDiagonalScaleFix(ksp::AbstractKSP{PetscLib},fix::PetscBool) where {PetscLib}

	LibPETSc.KSPSetDiagonalScaleFix(
		PetscLib,
		ksp,
		fix,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	fix = KSPGetDiagonalScaleFix(ksp::AbstractKSP{PetscLib})

Determines if `KSP` diagonally scales the system back after solving. That is `KSPSetDiagonalScaleFix()` has been called

Not Collective

Input Parameter:
===
- `ksp` - the `KSP` context

Output Parameter:
===
- `fix` - `PETSC_TRUE` to scale back after the system solve, `PETSC_FALSE` to not
rescale (default)

Level: intermediate

-seealso: [](ch_ksp), `KSPGetDiagonalScale()`, `KSPSetDiagonalScale()`, `KSPSetDiagonalScaleFix()`, `KSP`

# External Links
$(_doc_external("DM/KSPGetDiagonalScaleFix"))
"""
function KSPGetDiagonalScaleFix(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	fix = Ref{PetscBool}()

	LibPETSc.KSPGetDiagonalScaleFix(
		PetscLib,
		ksp,
		fix,
	)

	return fix[] == PETSC_TRUE
end
 
 
"""
	 UNTESTED !!!
	 KSPView(ksp::AbstractKSP{PetscLib},viewer::PetscViewer)

Prints the `KSP` data structure.

Collective

Input Parameters:
===
- `ksp`    - the Krylov space context
- `viewer` - visualization context

Options Database Key:
===
- `-ksp_view` - print the `KSP` data structure at the end of each `KSPSolve()` call

Level: beginner

Notes:
The available visualization contexts include
- ``PETSC_VIEWER_STDOUT_SELF``     - standard output (default)
- ``PETSC_VIEWER_STDOUT_WORLD``     - synchronized standard
output where only the first processor opens
the file.  All other processors send their
data to the first processor to print.

The available formats include
- ``PETSC_VIEWER_DEFAULT``     - standard output (default)
- ``PETSC_VIEWER_ASCII_INFO_DETAIL``     - more verbose output for PCBJACOBI and PCASM

The user can open an alternative visualization context with
`PetscViewerASCIIOpen()` - output to a specified file.

In the debugger you can do call `KSPView(ksp,0)` to display the `KSP`. (The same holds for any PETSc object viewer).

-seealso: [](ch_ksp), `KSP`, `PetscViewer`, `PCView()`, `PetscViewerASCIIOpen()`

# External Links
$(_doc_external("DM/KSPView"))
"""
function KSPView(ksp::AbstractKSP{PetscLib},viewer::PetscViewer) where {PetscLib}

	LibPETSc.KSPView(
		PetscLib,
		ksp,
		viewer,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 KSPLoad(newdm::AbstractKSP{PetscLib},viewer::PetscViewer)

Loads a `KSP` that has been stored in a `PETSCVIEWERBINARY`  with `KSPView()`.

Collective

Input Parameters:
===
- `newdm`  - the newly loaded `KSP`, this needs to have been created with `KSPCreate()` or
some related function before a call to `KSPLoad()`.
- `viewer` - binary file viewer, obtained from `PetscViewerBinaryOpen()`

Level: intermediate

Note:
The type is determined by the data in the file, any type set into the `KSP` before this call is ignored.

-seealso: [](ch_ksp), `KSP`, `PetscViewerBinaryOpen()`, `KSPView()`, `MatLoad()`, `VecLoad()`

# External Links
$(_doc_external("DM/KSPLoad"))
"""
function KSPLoad(newdm::AbstractKSP{PetscLib},viewer::PetscViewer) where {PetscLib}

	LibPETSc.KSPLoad(
		PetscLib,
		newdm,
		viewer,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 KSPViewFromOptions(A::AbstractKSP{PetscLib},obj::PetscObject,name::Vector{Char})

View a `KSP` object based on values in the options database

Collective

Input Parameters:
===
- `A`    - Krylov solver context
- `obj`  - Optional object
- `name` - command line option

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPView`, `PetscObjectViewFromOptions()`, `KSPCreate()`

# External Links
$(_doc_external("DM/KSPViewFromOptions"))
"""
function KSPViewFromOptions(A::AbstractKSP{PetscLib},obj::PetscObject,name::Vector{Char}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.KSPViewFromOptions(
		PetscLib,
		A,
		obj,
		name,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 KSPConvergedReasonViewFromOptions(ksp::AbstractKSP{PetscLib})

Processes command line options to determine if/how a `KSPReason` is to be viewed.

Collective

Input Parameter:
===
- `ksp` - the `KSP` object

Level: intermediate

-seealso: [](ch_ksp), `KSPConvergedReasonView()`, `KSPConvergedReasonViewSet()`

# External Links
$(_doc_external("DM/KSPConvergedReasonViewFromOptions"))
"""
function KSPConvergedReasonViewFromOptions(ksp::AbstractKSP{PetscLib}) where {PetscLib}

	LibPETSc.KSPConvergedReasonViewFromOptions(
		PetscLib,
		ksp,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 KSPConvergedReasonViewCancel(ksp::AbstractKSP{PetscLib})

Clears all the reasonview functions for a `KSP` object set with `KSPConvergedReasonViewSet()`
as well as the default viewer.

Collective

Input Parameter:
===
- `ksp` - iterative context obtained from `KSPCreate()`

Level: intermediate

-seealso: [](ch_ksp), `KSPCreate()`, `KSPDestroy()`, `KSPReset()`, `KSPConvergedReasonViewSet()`

# External Links
$(_doc_external("DM/KSPConvergedReasonViewCancel"))
"""
function KSPConvergedReasonViewCancel(ksp::AbstractKSP{PetscLib}) where {PetscLib}

	LibPETSc.KSPConvergedReasonViewCancel(
		PetscLib,
		ksp,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 KSPConvergedRateView(ksp::AbstractKSP{PetscLib},viewer::PetscViewer)

Displays the convergence rate <https://en.wikipedia.org/wiki/Coefficient_of_determination> of `KSPSolve()` to a viewer

Collective

Input Parameters:
===
- `ksp`    - iterative context obtained from `KSPCreate()`
- `viewer` - the viewer to display the reason

Options Database Key:
===
- `-ksp_converged_rate` - print reason for convergence or divergence and the convergence rate (or 0.0 for divergence)

Level: intermediate

Notes:
To change the format of the output, call `PetscViewerPushFormat`(`viewer`,`format`) before this call.

Suppose that the residual is reduced linearly, r_k = c^k r_0, which means log r_k = log r_0 + k log c. After linear regression,
the slope is log c. The coefficient of determination is given by 1 - frac{sum_i (y_i - f(x_i))^2}{sum_i (y_i - bar y)},

-seealso: [](ch_ksp), `KSPConvergedReasonView()`, `KSPGetConvergedRate()`, `KSPSetTolerances()`, `KSPConvergedDefault()`

# External Links
$(_doc_external("DM/KSPConvergedRateView"))
"""
function KSPConvergedRateView(ksp::AbstractKSP{PetscLib},viewer::PetscViewer) where {PetscLib}

	LibPETSc.KSPConvergedRateView(
		PetscLib,
		ksp,
		viewer,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 KSPSetNormType(ksp::AbstractKSP{PetscLib},normtype::KSPNormType)

Sets the norm that is used for convergence testing.

Logically Collective

Input Parameters:
===
- `ksp`      - Krylov solver context
- `normtype` - one of
-vb
KSP_NORM_NONE             - skips computing the norm, this should generally only be used if you are using
the Krylov method as a smoother with a fixed small number of iterations.
Implicitly sets `KSPConvergedSkip()` as the `KSP` convergence test.
Note that certain algorithms such as `KSPGMRES` ALWAYS require the norm calculation,
for these methods the norms are still computed, they are just not used in
the convergence test.
KSP_NORM_PRECONDITIONED   - the default for left-preconditioned solves, uses the l2 norm
of the preconditioned residual  P^{-1}(b - A x).
KSP_NORM_UNPRECONDITIONED - uses the l2 norm of the true  b - Ax residual.
KSP_NORM_NATURAL          - supported by `KSPCG`, `KSPCR`, `KSPCGNE`, `KSPCGS`
-ve

Options Database Key:
===
- `-ksp_norm_type <none,preconditioned,unpreconditioned,natural>` - set `KSP` norm type

Level: advanced

Note:
Not all combinations of preconditioner side (see `KSPSetPCSide()`) and norm type are supported by all Krylov methods.
If only one is set, PETSc tries to automatically change the other to find a compatible pair.  If no such combination
is supported, PETSc will generate an error.

Developer Note:
Supported combinations of norm and preconditioner side are set using `KSPSetSupportedNorm()`.

-seealso: [](ch_ksp), `KSPSetUp()`, `KSPSolve()`, `KSPDestroy()`, `KSPConvergedSkip()`, `KSPSetCheckNormIteration()`, `KSPSetPCSide()`, `KSPGetPCSide()`, `KSPNormType`

# External Links
$(_doc_external("DM/KSPSetNormType"))
"""
function KSPSetNormType(ksp::AbstractKSP{PetscLib},normtype::KSPNormType) where {PetscLib}

	LibPETSc.KSPSetNormType(
		PetscLib,
		ksp,
		normtype,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	normtype = KSPGetNormType(ksp::AbstractKSP{PetscLib})

Gets the norm that is used for convergence testing.

Not Collective

Input Parameter:
===
- `ksp` - Krylov solver context

Output Parameter:
===
- `normtype` - norm that is used for convergence testing

Level: advanced

-seealso: [](ch_ksp), `KSPNormType`, `KSPSetNormType()`, `KSPConvergedSkip()`

# External Links
$(_doc_external("DM/KSPGetNormType"))
"""
function KSPGetNormType(ksp::AbstractKSP{PetscLib}) where {PetscLib}

	LibPETSc.KSPGetNormType(
		PetscLib,
		ksp,
		normtype,
	)

	return normtype
end
 
 
"""
	 UNTESTED !!!
	 KSPSetCheckNormIteration(ksp::AbstractKSP{PetscLib},it::Int)

Sets the first iteration at which the norm of the residual will be
computed and used in the convergence test.

Logically Collective

Input Parameters:
===
- `ksp` - Krylov solver context
- `it`  - use -1 to check at all iterations

Notes:
Currently only works with `KSPCG`, `KSPBCGS` and `KSPIBCGS`

Use `KSPSetNormType`(ksp,`KSP_NORM_NONE`) to never check the norm

On steps where the norm is not computed, the previous norm is still in the variable, so if you run with, for example,
`-ksp_monitor` the residual norm will appear to be unchanged for several iterations (though it is not really unchanged).

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPSetUp()`, `KSPSolve()`, `KSPDestroy()`, `KSPConvergedSkip()`, `KSPSetNormType()`

# External Links
$(_doc_external("DM/KSPSetCheckNormIteration"))
"""
function KSPSetCheckNormIteration(ksp::AbstractKSP{PetscLib},it::Int) where {PetscLib}

	LibPETSc.KSPSetCheckNormIteration(
		PetscLib,
		ksp,
		it,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 KSPSetLagNorm(ksp::AbstractKSP{PetscLib},flg::PetscBool)

Lags the residual norm calculation so that it is computed as part of the `MPI_Allreduce()` for
computing the inner products for the next iteration.  This can reduce communication costs at the expense of doing
one additional iteration.

Logically Collective

Input Parameters:
===
- `ksp` - Krylov solver context
- `flg` - `PETSC_TRUE` or `PETSC_FALSE`

Options Database Key:
===
- `-ksp_lag_norm` - lag the calculated residual norm

Level: advanced

Notes:
Currently only works with `KSPIBCGS`.

Use `KSPSetNormType`(ksp,`KSP_NORM_NONE`) to never check the norm

If you lag the norm and run with, for example, `-ksp_monitor`, the residual norm reported will be the lagged one.

-seealso: [](ch_ksp), `KSPSetUp()`, `KSPSolve()`, `KSPDestroy()`, `KSPConvergedSkip()`, `KSPSetNormType()`, `KSPSetCheckNormIteration()`

# External Links
$(_doc_external("DM/KSPSetLagNorm"))
"""
function KSPSetLagNorm(ksp::AbstractKSP{PetscLib},flg::PetscBool) where {PetscLib}

	LibPETSc.KSPSetLagNorm(
		PetscLib,
		ksp,
		flg,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	ctx = KSPGetConvergenceContext(ksp::AbstractKSP{PetscLib})

Gets the convergence context set with `KSPSetConvergenceTest()`.

Not Collective

Input Parameter:
===
- `ksp` - iterative context obtained from `KSPCreate()`

Output Parameter:
===
- `ctx` - monitoring context

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPConvergedDefault()`, `KSPSetConvergenceTest()`, `KSPGetConvergenceTest()`

# External Links
$(_doc_external("DM/KSPGetConvergenceContext"))
"""
function KSPGetConvergenceContext(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	PetscScalar = PetscLib.PetscScalar
	#TODO: your output is a vector; ensure that the size is correct!
	#It may involve: dims = DMStagGetGhostCorners(dm)[4:6]
	#dims = DMStagGetGhostCorners(dm)[4:6]    # dimensions including ghost values; set to 0 if not 2D/3D
	#dims = getindex(dims,findall(dims.>0))   # retrieve non-zero values
	#dmE  = DMStagGetEntriesPerElement(dm)    # dof per element
	dims = (X,)
	r_ctx = PETSc_RefPtr(dims, PetscScalar)

	LibPETSc.KSPGetConvergenceContext(
		PetscLib,
		ksp,
		r_ctx,
	)

	ctx = PETSc_unsafe_wrap(r_ctx, dims; own=false)

	return ctx
end
 
 
"""
	 UNTESTED !!!
	reason,ctx = KSPConvergedDefault(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat)

Default code to determine convergence of the linear iterative solvers

Collective

Input Parameters:
===
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - residual norm (may be estimated, depending on the method may be the preconditioned residual norm)
- `ctx`   - convergence context which must be created by `KSPConvergedDefaultCreate()`

Output Parameter:
===
- `reason` - the convergence reason; it is positive if the iteration has converged,
negative if the iteration has diverged, and `KSP_CONVERGED_ITERATING` otherwise

Options Database Keys:
===
- `-ksp_max_it`                                  - maximum number of linear iterations
- `-ksp_min_it`                                  - minimum number of linear iterations, defaults to 0
- `-ksp_rtol rtol`                               - relative tolerance used in default determination of convergence, i.e. if residual norm decreases by this factor than convergence is declared
- `-ksp_atol abstol`                             - absolute tolerance used in default convergence test, i.e. if residual norm is less than this then convergence is declared
- `-ksp_divtol tol`                              - if residual norm increases by this factor than divergence is declared
- `-ksp_converged_use_initial_residual_norm`     - see `KSPConvergedDefaultSetUIRNorm()`
- `-ksp_converged_use_min_initial_residual_norm` - see `KSPConvergedDefaultSetUMIRNorm()`
- `-ksp_converged_maxits`                        - see `KSPConvergedDefaultSetConvergedMaxits()`

Level: advanced

Notes:
`KSPConvergedDefault()` reaches convergence when   rnorm < MAX (rtol * rnorm_0, abstol);
Divergence is detected if rnorm > dtol * rnorm_0, or when failures are detected throughout the iteration.
By default, reaching the maximum number of iterations is considered divergence (i.e. KSP_DIVERGED_ITS).
In order to have PETSc declaring convergence in such a case (i.e. `KSP_CONVERGED_ITS`), users can use `KSPConvergedDefaultSetConvergedMaxits()`

where:
- ``rtol``     - relative tolerance,
- ``abstol``     - absolute tolerance.
- ``dtol``     - divergence tolerance,
- ``rnorm_0``     - the two norm of the right-hand side (or the preconditioned norm, depending on what was set with
`KSPSetNormType()`. When initial guess is non-zero you
can call `KSPConvergedDefaultSetUIRNorm()` to use the norm of (b - A*(initial guess))
as the starting point for relative norm convergence testing, that is as `rnorm_0`.
Call `KSPConvergedDefaultSetUMIRNorm()` to use the minimum of the norm of (b - A*(initial guess)) and the norm of b as the starting point.

Use `KSPSetTolerances()` to alter the defaults for `rtol`, `abstol`, `dtol`.

Use `KSPSetNormType()` (or `-ksp_norm_type <none,preconditioned,unpreconditioned,natural>`) to change the norm used for computing rnorm

The precise values of reason are available in `KSPConvergedReason`

This routine is used by `KSP` by default so the user generally never needs call it directly.

Use `KSPSetConvergenceTest()` to provide your own test instead of using this one.

Call `KSPSetConvergenceTest()` with the `ctx`, as created above and the destruction function `KSPConvergedDefaultDestroy()`

-seealso: [](ch_ksp), `KSP`, `KSPSetConvergenceTest()`, `KSPSetTolerances()`, `KSPConvergedSkip()`, `KSPConvergedReason`, `KSPGetConvergedReason()`, `KSPSetMinimumIterations()`,
`KSPConvergedDefaultSetUIRNorm()`, `KSPConvergedDefaultSetUMIRNorm()`, `KSPConvergedDefaultSetConvergedMaxits()`, `KSPConvergedDefaultCreate()`, `KSPConvergedDefaultDestroy()`

# External Links
$(_doc_external("DM/KSPConvergedDefault"))
"""
function KSPConvergedDefault(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat) where {PetscLib}
	PetscScalar = PetscLib.PetscScalar
	#TODO: your output is a vector; ensure that the size is correct!
	#It may involve: dims = DMStagGetGhostCorners(dm)[4:6]
	#dims = DMStagGetGhostCorners(dm)[4:6]    # dimensions including ghost values; set to 0 if not 2D/3D
	#dims = getindex(dims,findall(dims.>0))   # retrieve non-zero values
	#dmE  = DMStagGetEntriesPerElement(dm)    # dof per element
	dims = (X,)
	r_ctx = PETSc_RefPtr(dims, PetscScalar)

	LibPETSc.KSPConvergedDefault(
		PetscLib,
		ksp,
		n,
		rnorm,
		reason,
		r_ctx,
	)

	ctx = PETSc_unsafe_wrap(r_ctx, dims; own=false)

	return reason,ctx
end
 
 
 
 
 
"""
	 UNTESTED !!!
	 KSPConvergedDefaultSetUIRNorm(ksp::AbstractKSP{PetscLib})

makes the default convergence test use  || B*(b
instead of  || B*b ||. In the case of right preconditioner or if `KSPSetNormType`(ksp,`KSP_NORM_UNPRECONDITIONED`)
is used there is no B in the above formula.

Collective

Input Parameters:
===
- `ksp` - iterative context

Options Database Key:
===
- `-ksp_converged_use_initial_residual_norm <bool>` - Use initial residual norm for computing relative convergence

Level: intermediate

Notes:
UIRNorm is short for Use Initial Residual Norm.

Use `KSPSetTolerances()` to alter the defaults for rtol, abstol, dtol.

The precise values of reason are macros such as `KSP_CONVERGED_RTOL`, which
are defined in petscksp.h.

If the convergence test is not `KSPConvergedDefault()` then this is ignored.

If right preconditioning is being used then B does not appear in the above formula.

-seealso: [](ch_ksp), `KSP`, `KSPSetConvergenceTest()`, `KSPSetTolerances()`, `KSPConvergedSkip()`, `KSPConvergedReason`, `KSPGetConvergedReason()`, `KSPConvergedDefaultSetUMIRNorm()`, `KSPConvergedDefaultSetConvergedMaxits()`

# External Links
$(_doc_external("DM/KSPConvergedDefaultSetUIRNorm"))
"""
function KSPConvergedDefaultSetUIRNorm(ksp::AbstractKSP{PetscLib}) where {PetscLib}

	LibPETSc.KSPConvergedDefaultSetUIRNorm(
		PetscLib,
		ksp,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 KSPConvergedDefaultSetUMIRNorm(ksp::AbstractKSP{PetscLib})

makes the default convergence test use min(|| B*(b
In the case of right preconditioner or if `KSPSetNormType`(ksp,`KSP_NORM_UNPRECONDITIONED`)
is used there is no B in the above formula.

Collective

Input Parameters:
===
- `ksp` - iterative context

Options Database Key:
===
- `-ksp_converged_use_min_initial_residual_norm <bool>` - Use minimum of initial residual norm and b for computing relative convergence

Level: intermediate

Notes:
UMIRNorm is short for Use Minimum Initial Residual Norm.

Use `KSPSetTolerances()` to alter the defaults for rtol, abstol, dtol.

-seealso: [](ch_ksp), `KSP`, `KSPSetConvergenceTest()`, `KSPSetTolerances()`, `KSPConvergedSkip()`, `KSPConvergedReason`, `KSPGetConvergedReason()`, `KSPConvergedDefaultSetUIRNorm()`, `KSPConvergedDefaultSetConvergedMaxits()`

# External Links
$(_doc_external("DM/KSPConvergedDefaultSetUMIRNorm"))
"""
function KSPConvergedDefaultSetUMIRNorm(ksp::AbstractKSP{PetscLib}) where {PetscLib}

	LibPETSc.KSPConvergedDefaultSetUMIRNorm(
		PetscLib,
		ksp,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 KSPConvergedDefaultSetConvergedMaxits(ksp::AbstractKSP{PetscLib},flg::PetscBool)

allows the default convergence test to declare convergence and return `KSP_CONVERGED_ITS` if the maximum number of iterations is reached

Collective

Input Parameters:
===
- `ksp` - iterative context
- `flg` - boolean flag

Options Database Key:
===
- `-ksp_converged_maxits <bool>` - Declare convergence if the maximum number of iterations is reached

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPSetConvergenceTest()`, `KSPSetTolerances()`, `KSPConvergedSkip()`, `KSPConvergedReason`, `KSPGetConvergedReason()`, `KSPConvergedDefaultSetUMIRNorm()`, `KSPConvergedDefaultSetUIRNorm()`

# External Links
$(_doc_external("DM/KSPConvergedDefaultSetConvergedMaxits"))
"""
function KSPConvergedDefaultSetConvergedMaxits(ksp::AbstractKSP{PetscLib},flg::PetscBool) where {PetscLib}

	LibPETSc.KSPConvergedDefaultSetConvergedMaxits(
		PetscLib,
		ksp,
		flg,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	reason,dtx = KSPConvergedSkip(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat)

Convergence test that do not return as converged
until the maximum number of iterations is reached.

Collective

Input Parameters:
===
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm residual value (may be estimated)
- `dtx`   - unused convergence context

Output Parameter:
===
- `reason` - `KSP_CONVERGED_ITERATING` or `KSP_CONVERGED_ITS`

Options Database Key:
===
- `-ksp_convergence_test skip` - skips the test

Level: advanced

Note:
This should be used as the convergence test with the option
`KSPSetNormType`(ksp,`KSP_NORM_NONE`), since norms of the residual are
not computed. Convergence is then declared after the maximum number
of iterations have been reached. Useful when one is using `KSPCG` or
`KSPBCGS`. [](sec_flexibleksp)

-seealso: [](ch_ksp), `KSP`, `KSPCG`, `KSPBCGS`, `KSPSetConvergenceTest()`, `KSPSetTolerances()`, `KSPSetNormType()`, [](sec_flexibleksp),
`KSPConvergedReason`

# External Links
$(_doc_external("DM/KSPConvergedSkip"))
"""
function KSPConvergedSkip(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat) where {PetscLib}
	PetscScalar = PetscLib.PetscScalar
	#TODO: your output is a vector; ensure that the size is correct!
	#It may involve: dims = DMStagGetGhostCorners(dm)[4:6]
	#dims = DMStagGetGhostCorners(dm)[4:6]    # dimensions including ghost values; set to 0 if not 2D/3D
	#dims = getindex(dims,findall(dims.>0))   # retrieve non-zero values
	#dmE  = DMStagGetEntriesPerElement(dm)    # dof per element
	dims = (X,)
	r_dtx = PETSc_RefPtr(dims, PetscScalar)

	LibPETSc.KSPConvergedSkip(
		PetscLib,
		ksp,
		n,
		rnorm,
		reason,
		r_dtx,
	)

	dtx = PETSc_unsafe_wrap(r_dtx, dims; own=false)

	return reason,dtx
end
 
 
"""
	 UNTESTED !!!
	reason = KSPGetConvergedReason(ksp::AbstractKSP{PetscLib})

Gets the reason the `KSP` iteration was stopped.

Not Collective

Input Parameter:
===
- `ksp` - the `KSP` context

Output Parameter:
===
- `reason` - negative value indicates diverged, positive value converged, see `KSPConvergedReason` for the possible values

Options Database Key:
===
- `-ksp_converged_reason` - prints the reason to standard out when the solve ends

Level: intermediate

Note:
If this routine is called before or doing the `KSPSolve()` the value of `KSP_CONVERGED_ITERATING` is returned

-seealso: [](ch_ksp), `KSPConvergedReason`, `KSP`, `KSPSetConvergenceTest()`, `KSPConvergedDefault()`, `KSPSetTolerances()`,
`KSPConvergedReasonView()`, `KSPGetConvergedReasonString()`

# External Links
$(_doc_external("DM/KSPGetConvergedReason"))
"""
function KSPGetConvergedReason(ksp::AbstractKSP{PetscLib}) where {PetscLib}

	LibPETSc.KSPGetConvergedReason(
		PetscLib,
		ksp,
		reason,
	)

	return reason
end
 
 
"""
	 UNTESTED !!!
	 KSPGetConvergedReasonString(ksp::AbstractKSP{PetscLib},strreason::Vector{Char})

Return a human readable string for a `KSPConvergedReason`

Not Collective

Input Parameter:
===
- `ksp` - the `KSP` context

Output Parameter:
===
- `strreason` - a human readable string that describes ksp converged reason

Level: beginner

-seealso: [](ch_ksp), `KSP`, `KSPGetConvergedReason()`

# External Links
$(_doc_external("DM/KSPGetConvergedReasonString"))
"""
function KSPGetConvergedReasonString(ksp::AbstractKSP{PetscLib},strreason::Vector{Char}) where {PetscLib}
	# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: 
	# @assert length() == n 
	# You can likely also write a multiple dispatch version of this function where vector length is determined automatically 

	LibPETSc.KSPGetConvergedReasonString(
		PetscLib,
		ksp,
		strreason,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	cr,rRsq,ce,eRsq = KSPComputeConvergenceRate(ksp::AbstractKSP{PetscLib})

Compute the convergence rate for the iteration <https:/en.wikipedia.org/wiki/Coefficient_of_determination>

Not Collective

Input Parameter:
===
- `ksp` - The `KSP`

Output Parameters:
===
- `cr`   - The residual contraction rate
- `rRsq` - The coefficient of determination, R^2, indicating the linearity of the data
- `ce`   - The error contraction rate
- `eRsq` - The coefficient of determination, R^2, indicating the linearity of the data

Level: advanced

Note:
Suppose that the residual is reduced linearly, r_k = c^k r_0, which means log r_k = log r_0 + k log c. After linear regression,
the slope is log c. The coefficient of determination is given by 1 - frac{sum_i (y_i - f(x_i))^2}{sum_i (y_i - bar y)},

-seealso: [](ch_ksp), `KSP`, `KSPConvergedRateView()`

# External Links
$(_doc_external("DM/KSPComputeConvergenceRate"))
"""
function KSPComputeConvergenceRate(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	PetscReal = PetscLib.PetscReal
	cr = [PetscReal(1)]
	PetscReal = PetscLib.PetscReal
	rRsq = [PetscReal(1)]
	PetscReal = PetscLib.PetscReal
	ce = [PetscReal(1)]
	PetscReal = PetscLib.PetscReal
	eRsq = [PetscReal(1)]

	LibPETSc.KSPComputeConvergenceRate(
		PetscLib,
		ksp,
		cr,
		rRsq,
		ce,
		eRsq,
	)

	return cr[1],rRsq[1],ce[1],eRsq[1]
end
 
 
"""
	 UNTESTED !!!
	 KSPSetConvergedNegativeCurvature(ksp::AbstractKSP{PetscLib},flg::PetscBool)

Allows to declare convergence and return `KSP_CONVERGED_NEG_CURVE` when negative curvature is detected

Collective

Input Parameters:
===
- `ksp` - iterative context
- `flg` - the Boolean value

Options Database Key:
===
- `-ksp_converged_neg_curve <bool>` - Declare convergence if negative curvature is detected

Level: advanced

Note:
This is currently used only by a subset of the Krylov solvers, namely `KSPCG`, `KSPSTCG`, `KSPQCG`, `KSPGLTR`, `KSPNASH`, and `KSPMINRES`.

-seealso: [](ch_ksp), `KSP`, `KSPConvergedReason`, `KSPGetConvergedNegativeCurvature()`

# External Links
$(_doc_external("DM/KSPSetConvergedNegativeCurvature"))
"""
function KSPSetConvergedNegativeCurvature(ksp::AbstractKSP{PetscLib},flg::PetscBool) where {PetscLib}

	LibPETSc.KSPSetConvergedNegativeCurvature(
		PetscLib,
		ksp,
		flg,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	flg = KSPGetConvergedNegativeCurvature(ksp::AbstractKSP{PetscLib})

Get the flag to declare convergence if negative curvature is detected

Collective

Input Parameter:
===
- `ksp` - iterative context

Output Parameter:
===
- `flg` - the Boolean value

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPConvergedReason`, `KSPSetConvergedNegativeCurvature()`

# External Links
$(_doc_external("DM/KSPGetConvergedNegativeCurvature"))
"""
function KSPGetConvergedNegativeCurvature(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	flg = Ref{PetscBool}()

	LibPETSc.KSPGetConvergedNegativeCurvature(
		PetscLib,
		ksp,
		flg,
	)

	return flg[] == PETSC_TRUE
end
 
 
"""
	 UNTESTED !!!
	mat = KSPComputeOperator(ksp::AbstractKSP{PetscLib},mattype::MatType)

Computes the explicit preconditioned operator, including diagonal scaling and null
space removal if applicable.

Collective

Input Parameters:
===
- `ksp`     - the Krylov subspace context
- `mattype` - the matrix type to be used

Output Parameter:
===
- `mat` - the explicit preconditioned operator

Level: advanced

Notes:
This computation is done by applying the operators to columns of the
identity matrix.

Currently, this routine uses a dense matrix format for the output operator if `mattype` is `NULL`.
This routine is costly in general, and is recommended for use only with relatively small systems.

-seealso: [](ch_ksp), `KSP`, `KSPSetOperators()`, `KSPComputeEigenvaluesExplicitly()`, `PCComputeOperator()`, `KSPSetDiagonalScale()`, `KSPSetNullSpace()`, `MatType`

# External Links
$(_doc_external("DM/KSPComputeOperator"))
"""
function KSPComputeOperator(ksp::AbstractKSP{PetscLib},mattype::MatType) where {PetscLib}

	LibPETSc.KSPComputeOperator(
		PetscLib,
		ksp,
		mattype,
		mat,
	)

	return mat
end
 
 
"""
	 UNTESTED !!!
	monctx = KSPMonitorLGRange(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat)


# External Links
$(_doc_external("DM/KSPMonitorLGRange"))
"""
function KSPMonitorLGRange(ksp::AbstractKSP{PetscLib},n::Int,rnorm::AbstractFloat) where {PetscLib}
	PetscScalar = PetscLib.PetscScalar
	#TODO: your output is a vector; ensure that the size is correct!
	#It may involve: dims = DMStagGetGhostCorners(dm)[4:6]
	#dims = DMStagGetGhostCorners(dm)[4:6]    # dimensions including ghost values; set to 0 if not 2D/3D
	#dims = getindex(dims,findall(dims.>0))   # retrieve non-zero values
	#dmE  = DMStagGetEntriesPerElement(dm)    # dof per element
	dims = (X,)
	r_monctx = PETSc_RefPtr(dims, PetscScalar)

	LibPETSc.KSPMonitorLGRange(
		PetscLib,
		ksp,
		n,
		rnorm,
		r_monctx,
	)

	monctx = PETSc_unsafe_wrap(r_monctx, dims; own=false)

	return monctx
end
 
 
"""
	 UNTESTED !!!
	 KSPSetGuess(ksp::AbstractKSP{PetscLib},guess::KSPGuess)

Set the initial guess object

Logically Collective

Input Parameters:
===
- `ksp`   - the Krylov context
- `guess` - the object created with `KSPGuessCreate()`

Level: advanced

Notes:
this allows a single `KSP` to be used with several different initial guess generators (likely for different linear
solvers, see `KSPSetPC()`).

This increases the reference count of the guess object, you must destroy the object with `KSPGuessDestroy()`
before the end of the program.

-seealso: [](ch_ksp), `KSP`, `KSPGuess`, `KSPSetOptionsPrefix()`, `KSPAppendOptionsPrefix()`, `KSPSetUseFischerGuess()`, `KSPGetGuess()`

# External Links
$(_doc_external("DM/KSPSetGuess"))
"""
function KSPSetGuess(ksp::AbstractKSP{PetscLib},guess::KSPGuess) where {PetscLib}

	LibPETSc.KSPSetGuess(
		PetscLib,
		ksp,
		guess,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	guess = KSPGetGuess(ksp::AbstractKSP{PetscLib})

Gets the initial guess generator for the `KSP`.

Not Collective

Input Parameter:
===
- `ksp` - the Krylov context

Output Parameter:
===
- `guess` - the object

Level: developer

-seealso: [](ch_ksp), `KSPGuess`, `KSP`, `KSPSetOptionsPrefix()`, `KSPAppendOptionsPrefix()`, `KSPSetUseFischerGuess()`, `KSPSetGuess()`

# External Links
$(_doc_external("DM/KSPGetGuess"))
"""
function KSPGetGuess(ksp::AbstractKSP{PetscLib}) where {PetscLib}

	LibPETSc.KSPGetGuess(
		PetscLib,
		ksp,
		guess,
	)

	return guess
end
 

"""
	 UNTESTED !!!
	 KSPSetUseFischerGuess(ksp::AbstractKSP{PetscLib},model::Int,size::Int)

Use the Paul Fischer algorithm or its variants to compute initial guesses for a set of solves with related right

Logically Collective

Input Parameters:
===
- `ksp`   - the Krylov context
- `model` - use model 1, model 2, model 3, or any other number to turn it off
- `size`  - size of subspace used to generate initial guess

Options Database Key:
===
- `-ksp_fischer_guess <model,size>` - uses the Fischer initial guess generator for repeated linear solves

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPSetOptionsPrefix()`, `KSPAppendOptionsPrefix()`, `KSPSetGuess()`, `KSPGetGuess()`, `KSPGuess`

# External Links
$(_doc_external("DM/KSPSetUseFischerGuess"))
"""
function KSPSetUseFischerGuess(ksp::AbstractKSP{PetscLib},model::Int,size::Int) where {PetscLib}

	LibPETSc.KSPSetUseFischerGuess(
		PetscLib,
		ksp,
		model,
		size,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	 KSPSetInitialGuessKnoll(ksp::AbstractKSP{PetscLib},flg::PetscBool)

Tells the iterative solver to use `PCApply()` to compute the initial guess (The Knoll trick)

Logically Collective

Input Parameters:
===
- `ksp` - iterative context obtained from `KSPCreate()`
- `flg` - `PETSC_TRUE` or `PETSC_FALSE`

Level: advanced

Developer Note:
The Knoll trick is not currently implemented using the `KSPGuess` class

-seealso: [](ch_ksp), `KSPGetInitialGuessKnoll()`, `KSPSetInitialGuessNonzero()`, `KSPGetInitialGuessNonzero()`, `KSP`

# External Links
$(_doc_external("DM/KSPSetInitialGuessKnoll"))
"""
function KSPSetInitialGuessKnoll(ksp::AbstractKSP{PetscLib},flg::PetscBool) where {PetscLib}

	LibPETSc.KSPSetInitialGuessKnoll(
		PetscLib,
		ksp,
		flg,
	)

	return nothing
end
 
 
"""
	 UNTESTED !!!
	flag = KSPGetInitialGuessKnoll(ksp::AbstractKSP{PetscLib})

Determines whether the `KSP` solver is using the Knoll trick (using PCApply(pc,b,...) to compute
the initial guess

Not Collective

Input Parameter:
===
- `ksp` - iterative context obtained from `KSPCreate()`

Output Parameter:
===
- `flag` - `PETSC_TRUE` if using Knoll trick, else `PETSC_FALSE`

Level: advanced

-seealso: [](ch_ksp), `KSPSetInitialGuessKnoll()`, `KSPSetInitialGuessNonzero()`, `KSPGetInitialGuessNonzero()`, `KSP`

# External Links
$(_doc_external("DM/KSPGetInitialGuessKnoll"))
"""
function KSPGetInitialGuessKnoll(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	flag = Ref{PetscBool}()

	LibPETSc.KSPGetInitialGuessKnoll(
		PetscLib,
		ksp,
		flag,
	)

	return flag[] == PETSC_TRUE
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
	 KSPSetDMActive(ksp::AbstractKSP{PetscLib},flg::PetscBool)

Indicates the `DM` should be used to generate the linear system matrix and right

Logically Collective

Input Parameters:
===
- `ksp` - the `KSP`
- `flg` - use the `DM`

Level: intermediate

Note:
By default `KSPSetDM()` sets the `DM` as active, call `KSPSetDMActive`(ksp,`PETSC_FALSE`); after `KSPSetDM`(ksp,dm) to not have the `KSP` object use the `DM` to generate the matrices.

-seealso: [](ch_ksp), `KSP`, `DM`, `KSPGetDM()`, `KSPSetDM()`, `SNESSetDM()`, `KSPSetComputeOperators()`, `KSPSetComputeRHS()`, `KSPSetComputeInitialGuess()`

# External Links
$(_doc_external("DM/KSPSetDMActive"))
"""
function KSPSetDMActive(ksp::AbstractKSP{PetscLib},flg::PetscBool) where {PetscLib}

	LibPETSc.KSPSetDMActive(
		PetscLib,
		ksp,
		flg,
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
 
 
"""
	 UNTESTED !!!
	ctx = KSPSetApplicationContext(ksp::AbstractKSP{PetscLib})

Sets the optional user

Logically Collective

Input Parameters:
===
- `ksp` - the `KSP` context
- `ctx` - optional user context

Level: intermediate

Notes:
The user context is a way for users to attach any information to the `KSP` that they may need later when interacting with the `KSP`

Use `KSPGetApplicationContext()` to get access to the context at a later time.

Fortran Note:
To use this from Fortran you must write a Fortran interface definition for this
function that tells Fortran the Fortran derived data type that you are passing in as the ctx argument.

-seealso: [](ch_ksp), `KSP`, `KSPGetApplicationContext()`

# External Links
$(_doc_external("DM/KSPSetApplicationContext"))
"""
function KSPSetApplicationContext(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	PetscScalar = PetscLib.PetscScalar
	#TODO: your output is a vector; ensure that the size is correct!
	#It may involve: dims = DMStagGetGhostCorners(dm)[4:6]
	#dims = DMStagGetGhostCorners(dm)[4:6]    # dimensions including ghost values; set to 0 if not 2D/3D
	#dims = getindex(dims,findall(dims.>0))   # retrieve non-zero values
	#dmE  = DMStagGetEntriesPerElement(dm)    # dof per element
	dims = (X,)
	r_ctx = PETSc_RefPtr(dims, PetscScalar)

	LibPETSc.KSPSetApplicationContext(
		PetscLib,
		ksp,
		r_ctx,
	)

	ctx = PETSc_unsafe_wrap(r_ctx, dims; own=false)

	return ctx
end
 
 
"""
	 UNTESTED !!!
	ctx = KSPGetApplicationContext(ksp::AbstractKSP{PetscLib})

Gets the user

Not Collective

Input Parameter:
===
- `ksp` - `KSP` context

Output Parameter:
===
- `ctx` - user context

Level: intermediate

Fortran Note:
You may need to write a Fortran interface definition for this
function that tells Fortran the Fortran derived data type that you are passing in as the ctx argument.

-seealso: [](ch_ksp), `KSP`, `KSPSetApplicationContext()`

# External Links
$(_doc_external("DM/KSPGetApplicationContext"))
"""
function KSPGetApplicationContext(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	PetscScalar = PetscLib.PetscScalar
	#TODO: your output is a vector; ensure that the size is correct!
	#It may involve: dims = DMStagGetGhostCorners(dm)[4:6]
	#dims = DMStagGetGhostCorners(dm)[4:6]    # dimensions including ghost values; set to 0 if not 2D/3D
	#dims = getindex(dims,findall(dims.>0))   # retrieve non-zero values
	#dmE  = DMStagGetEntriesPerElement(dm)    # dof per element
	dims = (X,)
	r_ctx = PETSc_RefPtr(dims, PetscScalar)

	LibPETSc.KSPGetApplicationContext(
		PetscLib,
		ksp,
		r_ctx,
	)

	ctx = PETSc_unsafe_wrap(r_ctx, dims; own=false)

	return ctx
end
 
 
"""
	 UNTESTED !!!
	func,ctx = KSPSetComputeRHS(ksp::AbstractKSP{PetscLib})

set routine to compute the right

Logically Collective

Input Parameters:
===
- `ksp`  - the `KSP` context
- `func` - function to compute the right-hand side, see `KSPComputeRHSFn` for the calling sequence
- `ctx`  - optional context

Level: beginner

Note:
The routine you provide will be called EACH you call `KSPSolve()` to prepare the new right-hand side for that solve

-seealso: [](ch_ksp), `KSP`, `KSPSolve()`, `DMKSPSetComputeRHS()`, `KSPSetComputeOperators()`, `KSPSetOperators()`, `KSPComputeRHSFn`

# External Links
$(_doc_external("DM/KSPSetComputeRHS"))
"""
function KSPSetComputeRHS(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	PetscScalar = PetscLib.PetscScalar
	#TODO: your output is a vector; ensure that the size is correct!
	#It may involve: dims = DMStagGetGhostCorners(dm)[4:6]
	#dims = DMStagGetGhostCorners(dm)[4:6]    # dimensions including ghost values; set to 0 if not 2D/3D
	#dims = getindex(dims,findall(dims.>0))   # retrieve non-zero values
	#dmE  = DMStagGetEntriesPerElement(dm)    # dof per element
	dims = (X,)
	r_ctx = PETSc_RefPtr(dims, PetscScalar)

	LibPETSc.KSPSetComputeRHS(
		PetscLib,
		ksp,
		func,
		r_ctx,
	)

	ctx = PETSc_unsafe_wrap(r_ctx, dims; own=false)

	return func,ctx
end
 
 
"""
	 UNTESTED !!!
	func,ctx = KSPSetComputeOperators(ksp::AbstractKSP{PetscLib})

set routine to compute the linear operators

Logically Collective

Input Parameters:
===
- `ksp`  - the `KSP` context
- `func` - function to compute the operators, see `KSPComputeOperatorsFn` for the calling sequence
- `ctx`  - optional context

Level: beginner

Notes:
`func()` will be called automatically at the very next call to `KSPSolve()`. It will NOT be called at future `KSPSolve()` calls
unless either `KSPSetComputeOperators()` or `KSPSetOperators()` is called before that `KSPSolve()` is called. This allows the same system to be solved several times
with different right-hand side functions but is a confusing API since one might expect it to be called for each `KSPSolve()`

To reuse the same preconditioner for the next `KSPSolve()` and not compute a new one based on the most recently computed matrix call `KSPSetReusePreconditioner()`

Developer Note:
Perhaps this routine and `KSPSetComputeRHS()` could be combined into a new API that makes clear when new matrices are computing without requiring call this
routine to indicate when the new matrix should be computed.

-seealso: [](ch_ksp), `KSP`, `KSPSetOperators()`, `KSPSetComputeRHS()`, `DMKSPSetComputeOperators()`, `KSPSetComputeInitialGuess()`, `KSPComputeOperatorsFn`

# External Links
$(_doc_external("DM/KSPSetComputeOperators"))
"""
function KSPSetComputeOperators(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	PetscScalar = PetscLib.PetscScalar
	#TODO: your output is a vector; ensure that the size is correct!
	#It may involve: dims = DMStagGetGhostCorners(dm)[4:6]
	#dims = DMStagGetGhostCorners(dm)[4:6]    # dimensions including ghost values; set to 0 if not 2D/3D
	#dims = getindex(dims,findall(dims.>0))   # retrieve non-zero values
	#dmE  = DMStagGetEntriesPerElement(dm)    # dof per element
	dims = (X,)
	r_ctx = PETSc_RefPtr(dims, PetscScalar)

	LibPETSc.KSPSetComputeOperators(
		PetscLib,
		ksp,
		func,
		r_ctx,
	)

	ctx = PETSc_unsafe_wrap(r_ctx, dims; own=false)

	return func,ctx
end
 
 
"""
	 UNTESTED !!!
	func,ctx = KSPSetComputeInitialGuess(ksp::AbstractKSP{PetscLib})

set routine to compute the initial guess of the linear system

Logically Collective

Input Parameters:
===
- `ksp`  - the `KSP` context
- `func` - function to compute the initial guess, see `KSPComputeInitialGuessFn` for calling sequence
- `ctx`  - optional context

Level: beginner

Note:
This should only be used in conjunction with `KSPSetComputeRHS()` and `KSPSetComputeOperators()`, otherwise
call `KSPSetInitialGuessNonzero()` and set the initial guess values in the solution vector passed to `KSPSolve()` before calling the solver

-seealso: [](ch_ksp), `KSP`, `KSPSolve()`, `KSPSetComputeRHS()`, `KSPSetComputeOperators()`, `DMKSPSetComputeInitialGuess()`, `KSPSetInitialGuessNonzero()`,
`KSPComputeInitialGuessFn`

# External Links
$(_doc_external("DM/KSPSetComputeInitialGuess"))
"""
function KSPSetComputeInitialGuess(ksp::AbstractKSP{PetscLib}) where {PetscLib}
	PetscScalar = PetscLib.PetscScalar
	#TODO: your output is a vector; ensure that the size is correct!
	#It may involve: dims = DMStagGetGhostCorners(dm)[4:6]
	#dims = DMStagGetGhostCorners(dm)[4:6]    # dimensions including ghost values; set to 0 if not 2D/3D
	#dims = getindex(dims,findall(dims.>0))   # retrieve non-zero values
	#dmE  = DMStagGetEntriesPerElement(dm)    # dof per element
	dims = (X,)
	r_ctx = PETSc_RefPtr(dims, PetscScalar)

	LibPETSc.KSPSetComputeInitialGuess(
		PetscLib,
		ksp,
		func,
		r_ctx,
	)

	ctx = PETSc_unsafe_wrap(r_ctx, dims; own=false)

	return func,ctx
end
 
 
"""
	 UNTESTED !!!
	func,ctx = DMKSPSetComputeOperators(dm::AbstractDM{PetscLib})

set `KSP` matrix evaluation function

Not Collective

Input Parameters:
===
- `dm`   - `DM` to be used with `KSP`
- `func` - matrix evaluation function,  for calling sequence see `KSPComputeOperatorsFn`
- `ctx`  - context for matrix evaluation

Level: developer

Note:
`KSPSetComputeOperators()` is normally used, but it calls this function internally because the user context is actually
associated with the `DM`.  This makes the interface consistent regardless of whether the user interacts with a `DM` or
not.

Developer Note:
If `DM` took a more central role at some later date, this could become the primary method of setting the matrix.

-seealso: [](ch_ksp), `DMKSP`, `DM`, `KSP`, `DMKSPSetContext()`, `DMKSPGetComputeOperators()`, `KSPSetOperators()`, `KSPComputeOperatorsFn`

# External Links
$(_doc_external("DM/DMKSPSetComputeOperators"))
"""
function DMKSPSetComputeOperators(dm::AbstractDM{PetscLib}) where {PetscLib}
	PetscScalar = PetscLib.PetscScalar
	#TODO: your output is a vector; ensure that the size is correct!
	#It may involve: dims = DMStagGetGhostCorners(dm)[4:6]
	#dims = DMStagGetGhostCorners(dm)[4:6]    # dimensions including ghost values; set to 0 if not 2D/3D
	#dims = getindex(dims,findall(dims.>0))   # retrieve non-zero values
	#dmE  = DMStagGetEntriesPerElement(dm)    # dof per element
	dims = (X,)
	r_ctx = PETSc_RefPtr(dims, PetscScalar)

	LibPETSc.DMKSPSetComputeOperators(
		PetscLib,
		dm,
		func,
		r_ctx,
	)

	ctx = PETSc_unsafe_wrap(r_ctx, dims; own=false)

	return func,ctx
end
 
 
"""
	 UNTESTED !!!
	func,ctx = DMKSPGetComputeOperators(dm::AbstractDM{PetscLib})

get `KSP` matrix evaluation function

Not Collective

Input Parameter:
===
- `dm` - `DM` used with a `KSP`

Output Parameters:
===
- `func` - matrix evaluation function,  for calling sequence see `KSPComputeOperatorsFn`
- `ctx`  - context for matrix evaluation

Level: developer

-seealso: [](ch_ksp), `DMKSP`, `DM`, `KSP`, `DMKSPSetContext()`, `KSPSetComputeOperators()`, `DMKSPSetComputeOperators()`, `KSPComputeOperatorsFn`

# External Links
$(_doc_external("DM/DMKSPGetComputeOperators"))
"""
function DMKSPGetComputeOperators(dm::AbstractDM{PetscLib}) where {PetscLib}
	PetscScalar = PetscLib.PetscScalar
	#TODO: your output is a vector; ensure that the size is correct!
	#It may involve: dims = DMStagGetGhostCorners(dm)[4:6]
	#dims = DMStagGetGhostCorners(dm)[4:6]    # dimensions including ghost values; set to 0 if not 2D/3D
	#dims = getindex(dims,findall(dims.>0))   # retrieve non-zero values
	#dmE  = DMStagGetEntriesPerElement(dm)    # dof per element
	dims = (X,)
	r_ctx = PETSc_RefPtr(dims, PetscScalar)

	LibPETSc.DMKSPGetComputeOperators(
		PetscLib,
		dm,
		func,
		r_ctx,
	)

	ctx = PETSc_unsafe_wrap(r_ctx, dims; own=false)

	return func,ctx
end
 
 
"""
	 UNTESTED !!!
	func,ctx = DMKSPSetComputeRHS(dm::AbstractDM{PetscLib})

set `KSP` right

Not Collective

Input Parameters:
===
- `dm`   - `DM` used with a `KSP`
- `func` - right-hand side evaluation function,  for calling sequence see `KSPComputeRHSFn`
- `ctx`  - context for right-hand side evaluation

Level: developer

Note:
`KSPSetComputeRHS()` is normally used, but it calls this function internally because the user context is actually
associated with the `DM`.  This makes the interface consistent regardless of whether the user interacts with a `DM` or
not.

Developer Note:
If `DM` took a more central role at some later date, this could become the primary method of setting the matrix.

-seealso: [](ch_ksp), `DMKSP`, `DM`, `KSP`, `DMKSPSetContext()`, `DMKSPGetComputeRHS()`, `KSPSetRHS()`

# External Links
$(_doc_external("DM/DMKSPSetComputeRHS"))
"""
function DMKSPSetComputeRHS(dm::AbstractDM{PetscLib}) where {PetscLib}
	PetscScalar = PetscLib.PetscScalar
	#TODO: your output is a vector; ensure that the size is correct!
	#It may involve: dims = DMStagGetGhostCorners(dm)[4:6]
	#dims = DMStagGetGhostCorners(dm)[4:6]    # dimensions including ghost values; set to 0 if not 2D/3D
	#dims = getindex(dims,findall(dims.>0))   # retrieve non-zero values
	#dmE  = DMStagGetEntriesPerElement(dm)    # dof per element
	dims = (X,)
	r_ctx = PETSc_RefPtr(dims, PetscScalar)

	LibPETSc.DMKSPSetComputeRHS(
		PetscLib,
		dm,
		func,
		r_ctx,
	)

	ctx = PETSc_unsafe_wrap(r_ctx, dims; own=false)

	return func,ctx
end
 
 
"""
	 UNTESTED !!!
	func,ctx = DMKSPGetComputeRHS(dm::AbstractDM{PetscLib})

get `KSP` right

Not Collective

Input Parameter:
===
- `dm` - `DM` to be used with `KSP`

Output Parameters:
===
- `func` - right-hand side evaluation function,  for calling sequence see `KSPComputeRHSFn`
- `ctx`  - context for right-hand side evaluation

Level: advanced

-seealso: [](ch_ksp), `DMKSP`, `DM`, `KSP`, `DMKSPSetContext()`, `KSPSetComputeRHS()`, `DMKSPSetComputeRHS()`, `KSPComputeRHSFn`

# External Links
$(_doc_external("DM/DMKSPGetComputeRHS"))
"""
function DMKSPGetComputeRHS(dm::AbstractDM{PetscLib}) where {PetscLib}
	PetscScalar = PetscLib.PetscScalar
	#TODO: your output is a vector; ensure that the size is correct!
	#It may involve: dims = DMStagGetGhostCorners(dm)[4:6]
	#dims = DMStagGetGhostCorners(dm)[4:6]    # dimensions including ghost values; set to 0 if not 2D/3D
	#dims = getindex(dims,findall(dims.>0))   # retrieve non-zero values
	#dmE  = DMStagGetEntriesPerElement(dm)    # dof per element
	dims = (X,)
	r_ctx = PETSc_RefPtr(dims, PetscScalar)

	LibPETSc.DMKSPGetComputeRHS(
		PetscLib,
		dm,
		func,
		r_ctx,
	)

	ctx = PETSc_unsafe_wrap(r_ctx, dims; own=false)

	return func,ctx
end
 
 
"""
	 UNTESTED !!!
	func,ctx = DMKSPSetComputeInitialGuess(dm::AbstractDM{PetscLib})

set `KSP` initial guess evaluation function

Not Collective

Input Parameters:
===
- `dm`   - `DM` to be used with `KSP`
- `func` - initial guess evaluation function, for calling sequence see `KSPComputeInitialGuessFn`
- `ctx`  - context for initial guess evaluation

Level: developer

Note:
`KSPSetComputeInitialGuess()` is normally used, but it calls this function internally because the user context is actually
associated with the `DM`.

-seealso: [](ch_ksp), `DMKSP`, `DM`, `KSP`, `DMKSPSetContext()`, `DMKSPGetComputeRHS()`, `KSPSetRHS()`, `KSPComputeInitialGuessFn`

# External Links
$(_doc_external("DM/DMKSPSetComputeInitialGuess"))
"""
function DMKSPSetComputeInitialGuess(dm::AbstractDM{PetscLib}) where {PetscLib}
	PetscScalar = PetscLib.PetscScalar
	#TODO: your output is a vector; ensure that the size is correct!
	#It may involve: dims = DMStagGetGhostCorners(dm)[4:6]
	#dims = DMStagGetGhostCorners(dm)[4:6]    # dimensions including ghost values; set to 0 if not 2D/3D
	#dims = getindex(dims,findall(dims.>0))   # retrieve non-zero values
	#dmE  = DMStagGetEntriesPerElement(dm)    # dof per element
	dims = (X,)
	r_ctx = PETSc_RefPtr(dims, PetscScalar)

	LibPETSc.DMKSPSetComputeInitialGuess(
		PetscLib,
		dm,
		func,
		r_ctx,
	)

	ctx = PETSc_unsafe_wrap(r_ctx, dims; own=false)

	return func,ctx
end
 
 
"""
	 UNTESTED !!!
	func,ctx = DMKSPGetComputeInitialGuess(dm::AbstractDM{PetscLib})

get `KSP` initial guess evaluation function

Not Collective

Input Parameter:
===
- `dm` - `DM` used with a `KSP`

Output Parameters:
===
- `func` - initial guess evaluation function, for calling sequence see `KSPComputeInitialGuessFn`
- `ctx`  - context for right-hand side evaluation

Level: advanced

-seealso: [](ch_ksp), `DMKSP`, `DM`, `KSP`, `DMKSPSetContext()`, `KSPSetComputeRHS()`, `DMKSPSetComputeRHS()`, `KSPComputeInitialGuessFn`

# External Links
$(_doc_external("DM/DMKSPGetComputeInitialGuess"))
"""
function DMKSPGetComputeInitialGuess(dm::AbstractDM{PetscLib}) where {PetscLib}
	PetscScalar = PetscLib.PetscScalar
	#TODO: your output is a vector; ensure that the size is correct!
	#It may involve: dims = DMStagGetGhostCorners(dm)[4:6]
	#dims = DMStagGetGhostCorners(dm)[4:6]    # dimensions including ghost values; set to 0 if not 2D/3D
	#dims = getindex(dims,findall(dims.>0))   # retrieve non-zero values
	#dmE  = DMStagGetEntriesPerElement(dm)    # dof per element
	dims = (X,)
	r_ctx = PETSc_RefPtr(dims, PetscScalar)

	LibPETSc.DMKSPGetComputeInitialGuess(
		PetscLib,
		dm,
		func,
		r_ctx,
	)

	ctx = PETSc_unsafe_wrap(r_ctx, dims; own=false)

	return func,ctx
end
 
 
"""
	 UNTESTED !!!
	 KSPSetPC(ksp::AbstractKSP{PetscLib},pc::AbstractPC{PetscLib})

Sets the preconditioner to be used to calculate the
application of the preconditioner on a vector.

Collective

Input Parameters:
===
- `ksp` - iterative context obtained from `KSPCreate()`
- `pc`  - the preconditioner object (can be `NULL`)

Level: developer

Note:
Use `KSPGetPC()` to retrieve the preconditioner context.

-seealso: [](ch_ksp), `KSPGetPC()`, `KSP`

# External Links
$(_doc_external("DM/KSPSetPC"))
"""
function KSPSetPC(ksp::AbstractKSP{PetscLib},pc::AbstractPC{PetscLib}) where {PetscLib}

	LibPETSc.KSPSetPC(
		PetscLib,
		ksp,
		pc,
	)

	return nothing
end
 
 
