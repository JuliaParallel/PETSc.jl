# autodefined type arguments for class ------
mutable struct KSPConvergedReasonViewFn end

mutable struct KSPMonitorFn end

mutable struct KSPConvergenceTestFn end

mutable struct KSPComputeOperatorsFn end

mutable struct KSPComputeRHSFn end

mutable struct KSPComputeInitialGuessFn end

mutable struct _n_PeCtx end
const PeCtx = Ptr{_n_PeCtx}

mutable struct KSPPSolveFn end

mutable struct KSPMonitorRegisterFn end

mutable struct KSPMonitorRegisterCreateFn end

mutable struct KSPMonitorRegisterDestroyFn end

mutable struct _n_KSPGuess end
const KSPGuess = Ptr{_n_KSPGuess}

mutable struct KSPFlexibleModifyPCFn end

# -------------------------------------------------------
"""
	KSPComputeOperator(petsclib::PetscLibType,ksp::PetscKSP, mattype::MatType, mat::PetscMat) 
Computes the explicit preconditioned operator, including diagonal scaling and null
space removal if applicable.

Collective

Input Parameters:
- `ksp`     - the Krylov subspace context
- `mattype` - the matrix type to be used

Output Parameter:
- `mat` - the explicit preconditioned operator

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPSetOperators()`, `KSPComputeEigenvaluesExplicitly()`, `PCComputeOperator()`, `KSPSetDiagonalScale()`, `KSPSetNullSpace()`, `MatType`

# External Links
$(_doc_external("KSP/KSPComputeOperator"))
"""
function KSPComputeOperator(petsclib::PetscLibType, ksp::PetscKSP, mattype::MatType, mat::PetscMat) end

@for_petsc function KSPComputeOperator(petsclib::$UnionPetscLib, ksp::PetscKSP, mattype::MatType, mat::PetscMat )
	mat_ = Ref(mat.ptr)

    @chk ccall(
               (:KSPComputeOperator, $petsc_library),
               PetscErrorCode,
               (CKSP, MatType, Ptr{CMat}),
               ksp, mattype, mat_,
              )

	mat.ptr = C_NULL

	return nothing
end 

"""
	r::Vector{PetscReal},c::Vector{PetscReal} = KSPComputeEigenvaluesExplicitly(petsclib::PetscLibType,ksp::PetscKSP, nmax::PetscInt) 
Computes all of the eigenvalues of the
preconditioned operator using LAPACK.

Collective

Input Parameters:
- `ksp`  - iterative context obtained from `KSPCreate()`
- `nmax` - size of arrays `r` and `c`

Output Parameters:
- `r` - real part of computed eigenvalues, provided by user with a dimension at least of `n`
- `c` - complex part of computed eigenvalues, provided by user with a dimension at least of `n`

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPComputeEigenvalues()`, `KSPMonitorSingularValue()`, `KSPComputeExtremeSingularValues()`, `KSPSetOperators()`, `KSPSolve()`

# External Links
$(_doc_external("KSP/KSPComputeEigenvaluesExplicitly"))
"""
function KSPComputeEigenvaluesExplicitly(petsclib::PetscLibType, ksp::PetscKSP, nmax::PetscInt) end

@for_petsc function KSPComputeEigenvaluesExplicitly(petsclib::$UnionPetscLib, ksp::PetscKSP, nmax::$PetscInt )
	r = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!
	c = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:KSPComputeEigenvaluesExplicitly, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}),
               ksp, nmax, r, c,
              )


	return r,c
end 

"""
	KSPMonitorLGRange(petsclib::PetscLibType,ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, monctx::Cvoid) 

# External Links
$(_doc_external("KSP/KSPMonitorLGRange"))
"""
function KSPMonitorLGRange(petsclib::PetscLibType, ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, monctx::Cvoid) end

@for_petsc function KSPMonitorLGRange(petsclib::$UnionPetscLib, ksp::PetscKSP, n::$PetscInt, rnorm::$PetscReal, monctx::Cvoid )

    @chk ccall(
               (:KSPMonitorLGRange, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{Cvoid}),
               ksp, n, rnorm, monctx,
              )


	return nothing
end 

"""
	emax::PetscReal,emin::PetscReal = KSPComputeExtremeSingularValues(petsclib::PetscLibType,ksp::PetscKSP) 
Computes the extreme singular values
for the preconditioned operator. Called after or during `KSPSolve()`.

Not Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameters:
- `emax` - maximum estimated singular value
- `emin` - minimum estimated singular value

Options Database Key:
- `-ksp_view_singularvalues` - compute extreme singular values and print when `KSPSolve()` completes.

Level: advanced

-seealso: [](ch_ksp), `KSPSetComputeSingularValues()`, `KSPMonitorSingularValue()`, `KSPComputeEigenvalues()`, `KSP`, `KSPComputeRitz()`

# External Links
$(_doc_external("KSP/KSPComputeExtremeSingularValues"))
"""
function KSPComputeExtremeSingularValues(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPComputeExtremeSingularValues(petsclib::$UnionPetscLib, ksp::PetscKSP )
	emax_ = Ref{$PetscReal}()
	emin_ = Ref{$PetscReal}()

    @chk ccall(
               (:KSPComputeExtremeSingularValues, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscReal}, Ptr{$PetscReal}),
               ksp, emax_, emin_,
              )

	emax = emax_[]
	emin = emin_[]

	return emax,emin
end 

"""
	r::Vector{PetscReal},c::Vector{PetscReal},neig::PetscInt = KSPComputeEigenvalues(petsclib::PetscLibType,ksp::PetscKSP, n::PetscInt) 
Computes the extreme eigenvalues for the
preconditioned operator. Called after or during `KSPSolve()`.

Not Collective

Input Parameters:
- `ksp` - iterative solver obtained from `KSPCreate()`
- `n`   - size of arrays `r` and `c`. The number of eigenvalues computed `neig` will, in general, be less than this.

Output Parameters:
- `r`    - real part of computed eigenvalues, provided by user with a dimension of at least `n`
- `c`    - complex part of computed eigenvalues, provided by user with a dimension of at least `n`
- `neig` - actual number of eigenvalues computed (will be less than or equal to `n`)

Options Database Key:
- `-ksp_view_eigenvalues` - Prints eigenvalues to stdout

Level: advanced

-seealso: [](ch_ksp), `KSPSetComputeEigenvalues()`, `KSPSetComputeSingularValues()`, `KSPMonitorSingularValue()`, `KSPComputeExtremeSingularValues()`, `KSP`, `KSPComputeRitz()`

# External Links
$(_doc_external("KSP/KSPComputeEigenvalues"))
"""
function KSPComputeEigenvalues(petsclib::PetscLibType, ksp::PetscKSP, n::PetscInt) end

@for_petsc function KSPComputeEigenvalues(petsclib::$UnionPetscLib, ksp::PetscKSP, n::$PetscInt )
	r = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!
	c = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!
	neig_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPComputeEigenvalues, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscInt}),
               ksp, n, r, c, neig_,
              )

	neig = neig_[]

	return r,c,neig
end 

"""
	nrit::PetscInt,tetar::Vector{PetscReal},tetai::Vector{PetscReal} = KSPComputeRitz(petsclib::PetscLibType,ksp::PetscKSP, ritz::PetscBool, small::PetscBool, S::Vector{PetscVec}) 
Computes the Ritz or harmonic Ritz pairs associated with the
smallest or largest in modulus, for the preconditioned operator.

Not Collective

Input Parameters:
- `ksp`   - iterative solver obtained from `KSPCreate()`
- `ritz`  - `PETSC_TRUE` or `PETSC_FALSE` for Ritz pairs or harmonic Ritz pairs, respectively
- `small` - `PETSC_TRUE` or `PETSC_FALSE` for smallest or largest (harmonic) Ritz values, respectively

Output Parameters:
- `nrit`  - On input number of (harmonic) Ritz pairs to compute; on output, actual number of computed (harmonic) Ritz pairs
- `S`     - an array of the Ritz vectors, pass in an array of vectors of size `nrit`
- `tetar` - real part of the Ritz values, pass in an array of size `nrit`
- `tetai` - imaginary part of the Ritz values, pass in an array of size `nrit`

Level: advanced

-seealso: [](ch_ksp), `KSPSetComputeRitz()`, `KSP`, `KSPGMRES`, `KSPComputeEigenvalues()`, `KSPSetComputeSingularValues()`, `KSPMonitorSingularValue()`

# External Links
$(_doc_external("KSP/KSPComputeRitz"))
"""
function KSPComputeRitz(petsclib::PetscLibType, ksp::PetscKSP, ritz::PetscBool, small::PetscBool, S::Vector{PetscVec}) end

@for_petsc function KSPComputeRitz(petsclib::$UnionPetscLib, ksp::PetscKSP, ritz::PetscBool, small::PetscBool, S::Vector{PetscVec} )
	nrit_ = Ref{$PetscInt}()
	tetar = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!
	tetai = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:KSPComputeRitz, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool, PetscBool, Ptr{$PetscInt}, Ptr{CVec}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               ksp, ritz, small, nrit_, S, tetar, tetai,
              )

	nrit = nrit_[]

	return nrit,tetar,tetai
end 

"""
	KSPSetUpOnBlocks(petsclib::PetscLibType,ksp::PetscKSP) 
Sets up the preconditioner for each block in
the block Jacobi `PCJACOBI`, overlapping Schwarz `PCASM`, and fieldsplit `PCFIELDSPLIT` preconditioners

Collective

Input Parameter:
- `ksp` - the `KSP` context

Level: advanced

-seealso: [](ch_ksp), `PCSetUpOnBlocks()`, `KSPSetUp()`, `PCSetUp()`, `KSP`

# External Links
$(_doc_external("KSP/KSPSetUpOnBlocks"))
"""
function KSPSetUpOnBlocks(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPSetUpOnBlocks(petsclib::$UnionPetscLib, ksp::PetscKSP )

    @chk ccall(
               (:KSPSetUpOnBlocks, $petsc_library),
               PetscErrorCode,
               (CKSP,),
               ksp,
              )


	return nothing
end 

"""
	KSPSetReusePreconditioner(petsclib::PetscLibType,ksp::PetscKSP, flag::PetscBool) 
reuse the current preconditioner for future `KSPSolve()`, do not construct a new preconditioner even if the `Mat` operator
in the `KSP` has different values

Collective

Input Parameters:
- `ksp`  - iterative solver obtained from `KSPCreate()`
- `flag` - `PETSC_TRUE` to reuse the current preconditioner, or `PETSC_FALSE` to construct a new preconditioner

Options Database Key:
- `-ksp_reuse_preconditioner <true,false>` - reuse the previously computed preconditioner

Level: intermediate

-seealso: [](ch_ksp), `KSPCreate()`, `KSPSolve()`, `KSPDestroy()`, `KSP`, `KSPGetReusePreconditioner()`,
`SNESSetLagPreconditioner()`, `SNES`

# External Links
$(_doc_external("KSP/KSPSetReusePreconditioner"))
"""
function KSPSetReusePreconditioner(petsclib::PetscLibType, ksp::PetscKSP, flag::PetscBool) end

@for_petsc function KSPSetReusePreconditioner(petsclib::$UnionPetscLib, ksp::PetscKSP, flag::PetscBool )

    @chk ccall(
               (:KSPSetReusePreconditioner, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flag,
              )


	return nothing
end 

"""
	flag::PetscBool = KSPGetReusePreconditioner(petsclib::PetscLibType,ksp::PetscKSP) 
Determines if the `KSP` reuses the current preconditioner even if the `Mat` operator in the `KSP` has changed.

Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameter:
- `flag` - the boolean flag indicating if the current preconditioner should be reused

Level: intermediate

-seealso: [](ch_ksp), `KSPCreate()`, `KSPSolve()`, `KSPDestroy()`, `KSPSetReusePreconditioner()`, `KSP`

# External Links
$(_doc_external("KSP/KSPGetReusePreconditioner"))
"""
function KSPGetReusePreconditioner(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGetReusePreconditioner(petsclib::$UnionPetscLib, ksp::PetscKSP )
	flag_ = Ref{PetscBool}()

    @chk ccall(
               (:KSPGetReusePreconditioner, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{PetscBool}),
               ksp, flag_,
              )

	flag = flag_[]

	return flag
end 

"""
	KSPSetSkipPCSetFromOptions(petsclib::PetscLibType,ksp::PetscKSP, flag::PetscBool) 
prevents `KSPSetFromOptions()` from calling `PCSetFromOptions()`.
This is used if the same `PC` is shared by more than one `KSP` so its options are not reset for each `KSP`

Collective

Input Parameters:
- `ksp`  - iterative solver obtained from `KSPCreate()`
- `flag` - `PETSC_TRUE` to skip calling the `PCSetFromOptions()`

Level: developer

-seealso: [](ch_ksp), `KSPCreate()`, `KSPSolve()`, `KSPDestroy()`, `PCSetReusePreconditioner()`, `KSP`

# External Links
$(_doc_external("KSP/KSPSetSkipPCSetFromOptions"))
"""
function KSPSetSkipPCSetFromOptions(petsclib::PetscLibType, ksp::PetscKSP, flag::PetscBool) end

@for_petsc function KSPSetSkipPCSetFromOptions(petsclib::$UnionPetscLib, ksp::PetscKSP, flag::PetscBool )

    @chk ccall(
               (:KSPSetSkipPCSetFromOptions, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flag,
              )


	return nothing
end 

"""
	KSPSetUp(petsclib::PetscLibType,ksp::PetscKSP) 
Sets up the internal data structures for the
later use `KSPSolve()` the `KSP` linear iterative solver.

Collective

Input Parameter:
- `ksp` - iterative solver, `KSP`, obtained from `KSPCreate()`

Level: developer

-seealso: [](ch_ksp), `KSPCreate()`, `KSPSolve()`, `KSPDestroy()`, `KSP`, `KSPSetUpOnBlocks()`

# External Links
$(_doc_external("KSP/KSPSetUp"))
"""
function KSPSetUp(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPSetUp(petsclib::$UnionPetscLib, ksp::PetscKSP )

    @chk ccall(
               (:KSPSetUp, $petsc_library),
               PetscErrorCode,
               (CKSP,),
               ksp,
              )


	return nothing
end 

"""
	KSPConvergedReasonView(petsclib::PetscLibType,ksp::PetscKSP, viewer::PetscViewer) 
Displays the reason a `KSP` solve converged or diverged, `KSPConvergedReason` to a `PetscViewer`

Collective

Input Parameters:
- `ksp`    - iterative solver obtained from `KSPCreate()`
- `viewer` - the `PetscViewer` on which to display the reason

Options Database Keys:
- `-ksp_converged_reason`          - print reason for converged or diverged, also prints number of iterations
- `-ksp_converged_reason ::failed` - only print reason and number of iterations when diverged

Level: beginner

-seealso: [](ch_ksp), `KSPConvergedReasonViewFromOptions()`, `KSPCreate()`, `KSPSetUp()`, `KSPDestroy()`, `KSPSetTolerances()`, `KSPConvergedDefault()`,
`KSPSolveTranspose()`, `KSPGetIterationNumber()`, `KSP`, `KSPGetConvergedReason()`, `PetscViewerPushFormat()`, `PetscViewerPopFormat()`

# External Links
$(_doc_external("KSP/KSPConvergedReasonView"))
"""
function KSPConvergedReasonView(petsclib::PetscLibType, ksp::PetscKSP, viewer::PetscViewer) end

@for_petsc function KSPConvergedReasonView(petsclib::$UnionPetscLib, ksp::PetscKSP, viewer::PetscViewer )

    @chk ccall(
               (:KSPConvergedReasonView, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscViewer),
               ksp, viewer,
              )


	return nothing
end 

"""
	KSPConvergedReasonViewSet(petsclib::PetscLibType,ksp::PetscKSP, f::KSPConvergedReasonViewFn, vctx::Cvoid, reasonviewdestroy::PetscCtxDestroyFn) 
Sets an ADDITIONAL function that is to be used at the
end of the linear solver to display the convergence reason of the linear solver.

Logically Collective

Input Parameters:
- `ksp`               - the `KSP` context
- `f`                 - the `ksp` converged reason view function, see `KSPConvergedReasonViewFn`
- `vctx`              - [optional] user-defined context for private data for the
`KSPConvergedReason` view routine (use `NULL` if no context is desired)
- `reasonviewdestroy` - [optional] routine that frees `vctx` (may be `NULL`), see `PetscCtxDestroyFn` for the calling sequence

Options Database Keys:
- `-ksp_converged_reason`             - sets a default `KSPConvergedReasonView()`
- `-ksp_converged_reason_view_cancel` - cancels all converged reason viewers that have been hardwired into a code by
calls to `KSPConvergedReasonViewSet()`, but does not cancel those set via the options database.

Level: intermediate

-seealso: [](ch_ksp), `KSPConvergedReasonView()`, `KSPConvergedReasonViewFn`, `KSPConvergedReasonViewCancel()`, `PetscCtxDestroyFn`

# External Links
$(_doc_external("KSP/KSPConvergedReasonViewSet"))
"""
function KSPConvergedReasonViewSet(petsclib::PetscLibType, ksp::PetscKSP, f::KSPConvergedReasonViewFn, vctx::Cvoid, reasonviewdestroy::PetscCtxDestroyFn) end

@for_petsc function KSPConvergedReasonViewSet(petsclib::$UnionPetscLib, ksp::PetscKSP, f::KSPConvergedReasonViewFn, vctx::Cvoid, reasonviewdestroy::PetscCtxDestroyFn )

    @chk ccall(
               (:KSPConvergedReasonViewSet, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPConvergedReasonViewFn}, Ptr{Cvoid}, Ptr{PetscCtxDestroyFn}),
               ksp, f, vctx, reasonviewdestroy,
              )


	return nothing
end 

"""
	KSPConvergedReasonViewCancel(petsclib::PetscLibType,ksp::PetscKSP) 
Clears all the `KSPConvergedReason` view functions for a `KSP` object set with `KSPConvergedReasonViewSet()`
as well as the default viewer.

Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Level: intermediate

-seealso: [](ch_ksp), `KSPCreate()`, `KSPDestroy()`, `KSPReset()`, `KSPConvergedReasonViewSet()`

# External Links
$(_doc_external("KSP/KSPConvergedReasonViewCancel"))
"""
function KSPConvergedReasonViewCancel(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPConvergedReasonViewCancel(petsclib::$UnionPetscLib, ksp::PetscKSP )

    @chk ccall(
               (:KSPConvergedReasonViewCancel, $petsc_library),
               PetscErrorCode,
               (CKSP,),
               ksp,
              )


	return nothing
end 

"""
	KSPConvergedReasonViewFromOptions(petsclib::PetscLibType,ksp::PetscKSP) 
Processes command line options to determine if/how a `KSPReason` is to be viewed.

Collective

Input Parameter:
- `ksp` - the `KSP` object

Level: intermediate

-seealso: [](ch_ksp), `KSPConvergedReasonView()`, `KSPConvergedReasonViewSet()`

# External Links
$(_doc_external("KSP/KSPConvergedReasonViewFromOptions"))
"""
function KSPConvergedReasonViewFromOptions(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPConvergedReasonViewFromOptions(petsclib::$UnionPetscLib, ksp::PetscKSP )

    @chk ccall(
               (:KSPConvergedReasonViewFromOptions, $petsc_library),
               PetscErrorCode,
               (CKSP,),
               ksp,
              )


	return nothing
end 

"""
	KSPConvergedRateView(petsclib::PetscLibType,ksp::PetscKSP, viewer::PetscViewer) 
Displays the convergence rate <https://en.wikipedia.org/wiki/Coefficient_of_determination> of `KSPSolve()` to a viewer

Collective

Input Parameters:
- `ksp`    - iterative solver obtained from `KSPCreate()`
- `viewer` - the `PetscViewer` to display the reason

Options Database Key:
- `-ksp_converged_rate` - print reason for convergence or divergence and the convergence rate (or 0.0 for divergence)

Level: intermediate

-seealso: [](ch_ksp), `KSPConvergedReasonView()`, `KSPGetConvergedRate()`, `KSPSetTolerances()`, `KSPConvergedDefault()`

# External Links
$(_doc_external("KSP/KSPConvergedRateView"))
"""
function KSPConvergedRateView(petsclib::PetscLibType, ksp::PetscKSP, viewer::PetscViewer) end

@for_petsc function KSPConvergedRateView(petsclib::$UnionPetscLib, ksp::PetscKSP, viewer::PetscViewer )

    @chk ccall(
               (:KSPConvergedRateView, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscViewer),
               ksp, viewer,
              )


	return nothing
end 

"""
	KSPSolve(petsclib::PetscLibType,ksp::PetscKSP, b::PetscVec, x::PetscVec) 
Solves a linear system associated with `KSP` object

Collective

Input Parameters:
- `ksp` - iterative solver obtained from `KSPCreate()`
- `b`   - the right-hand side vector
- `x`   - the solution (this may be the same vector as `b`, then `b` will be overwritten with the answer)

Options Database Keys:
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

-seealso: [](ch_ksp), `KSPCreate()`, `KSPSetUp()`, `KSPDestroy()`, `KSPSetTolerances()`, `KSPConvergedDefault()`,
`KSPSolveTranspose()`, `KSPGetIterationNumber()`, `MatNullSpaceCreate()`, `MatSetNullSpace()`, `MatSetTransposeNullSpace()`, `KSP`,
`KSPConvergedReasonView()`, `KSPCheckSolve()`, `KSPSetErrorIfNotConverged()`

# External Links
$(_doc_external("KSP/KSPSolve"))
"""
function KSPSolve(petsclib::PetscLibType, ksp::PetscKSP, b::Union{PetscVec,Ptr}, x::Union{PetscVec,Ptr}) end

@for_petsc function KSPSolve(petsclib::$UnionPetscLib, ksp::PetscKSP, b::Union{PetscVec,Ptr}, x::Union{PetscVec,Ptr} )

    @chk ccall(
               (:KSPSolve, $petsc_library),
               PetscErrorCode,
               (CKSP, CVec, CVec),
               ksp, b, x,
              )


	return nothing
end 

"""
	KSPSolveTranspose(petsclib::PetscLibType,ksp::PetscKSP, b::PetscVec, x::PetscVec) 
Solves a linear system with the transpose of the matrix associated with the `KSP` object,  A^T x = b.

Collective

Input Parameters:
- `ksp` - iterative solver obtained from `KSPCreate()`
- `b`   - right-hand side vector
- `x`   - solution vector

Level: developer

-seealso: [](ch_ksp), `KSPCreate()`, `KSPSetUp()`, `KSPDestroy()`, `KSPSetTolerances()`, `KSPConvergedDefault()`,
`KSPSolve()`, `KSP`, `KSPSetOperators()`

# External Links
$(_doc_external("KSP/KSPSolveTranspose"))
"""
function KSPSolveTranspose(petsclib::PetscLibType, ksp::PetscKSP, b::PetscVec, x::PetscVec) end

@for_petsc function KSPSolveTranspose(petsclib::$UnionPetscLib, ksp::PetscKSP, b::PetscVec, x::PetscVec )

    @chk ccall(
               (:KSPSolveTranspose, $petsc_library),
               PetscErrorCode,
               (CKSP, CVec, CVec),
               ksp, b, x,
              )


	return nothing
end 

"""
	KSPMatSolve(petsclib::PetscLibType,ksp::PetscKSP, B::PetscMat, X::PetscMat) 
Solves a linear system with multiple right

Input Parameters:
- `ksp` - iterative solver
- `B`   - block of right-hand sides

Output Parameter:
- `X` - block of solutions

Level: intermediate

-seealso: [](ch_ksp), `KSPSolve()`, `MatMatSolve()`, `KSPMatSolveTranspose()`, `MATDENSE`, `KSPHPDDM`, `PCBJACOBI`, `PCASM`, `KSPSetMatSolveBatchSize()`

# External Links
$(_doc_external("KSP/KSPMatSolve"))
"""
function KSPMatSolve(petsclib::PetscLibType, ksp::PetscKSP, B::PetscMat, X::PetscMat) end

@for_petsc function KSPMatSolve(petsclib::$UnionPetscLib, ksp::PetscKSP, B::PetscMat, X::PetscMat )

    @chk ccall(
               (:KSPMatSolve, $petsc_library),
               PetscErrorCode,
               (CKSP, CMat, CMat),
               ksp, B, X,
              )


	return nothing
end 

"""
	KSPMatSolveTranspose(petsclib::PetscLibType,ksp::PetscKSP, B::PetscMat, X::PetscMat) 
Solves a linear system with the transposed matrix with multiple right

Input Parameters:
- `ksp` - iterative solver
- `B`   - block of right-hand sides

Output Parameter:
- `X` - block of solutions

Level: intermediate

-seealso: [](ch_ksp), `KSPSolveTranspose()`, `MatMatTransposeSolve()`, `KSPMatSolve()`, `MATDENSE`, `KSPHPDDM`, `PCBJACOBI`, `PCASM`

# External Links
$(_doc_external("KSP/KSPMatSolveTranspose"))
"""
function KSPMatSolveTranspose(petsclib::PetscLibType, ksp::PetscKSP, B::PetscMat, X::PetscMat) end

@for_petsc function KSPMatSolveTranspose(petsclib::$UnionPetscLib, ksp::PetscKSP, B::PetscMat, X::PetscMat )

    @chk ccall(
               (:KSPMatSolveTranspose, $petsc_library),
               PetscErrorCode,
               (CKSP, CMat, CMat),
               ksp, B, X,
              )


	return nothing
end 

"""
	KSPSetMatSolveBatchSize(petsclib::PetscLibType,ksp::PetscKSP, bs::PetscInt) 
Sets the maximum number of columns treated simultaneously in `KSPMatSolve()`.

Logically Collective

Input Parameters:
- `ksp` - the `KSP` iterative solver
- `bs`  - batch size

Level: advanced

-seealso: [](ch_ksp), `KSPMatSolve()`, `KSPGetMatSolveBatchSize()`, `-mat_mumps_icntl_27`, `-matmatmult_Bbn`

# External Links
$(_doc_external("KSP/KSPSetMatSolveBatchSize"))
"""
function KSPSetMatSolveBatchSize(petsclib::PetscLibType, ksp::PetscKSP, bs::PetscInt) end

@for_petsc function KSPSetMatSolveBatchSize(petsclib::$UnionPetscLib, ksp::PetscKSP, bs::$PetscInt )

    @chk ccall(
               (:KSPSetMatSolveBatchSize, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, bs,
              )


	return nothing
end 

"""
	bs::PetscInt = KSPGetMatSolveBatchSize(petsclib::PetscLibType,ksp::PetscKSP) 
Gets the maximum number of columns treated simultaneously in `KSPMatSolve()`.

Input Parameter:
- `ksp` - iterative solver context

Output Parameter:
- `bs` - batch size

Level: advanced

-seealso: [](ch_ksp), `KSPMatSolve()`, `KSPSetMatSolveBatchSize()`, `-mat_mumps_icntl_27`, `-matmatmult_Bbn`

# External Links
$(_doc_external("KSP/KSPGetMatSolveBatchSize"))
"""
function KSPGetMatSolveBatchSize(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGetMatSolveBatchSize(petsclib::$UnionPetscLib, ksp::PetscKSP )
	bs_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPGetMatSolveBatchSize, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscInt}),
               ksp, bs_,
              )

	bs = bs_[]

	return bs
end 

"""
	KSPResetViewers(petsclib::PetscLibType,ksp::PetscKSP) 
Resets all the viewers set from the options database during `KSPSetFromOptions()`

Collective

Input Parameter:
- `ksp` - the `KSP` iterative solver context obtained from `KSPCreate()`

Level: beginner

-seealso: [](ch_ksp), `KSPCreate()`, `KSPSetUp()`, `KSPSolve()`, `KSPSetFromOptions()`, `KSP`

# External Links
$(_doc_external("KSP/KSPResetViewers"))
"""
function KSPResetViewers(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPResetViewers(petsclib::$UnionPetscLib, ksp::PetscKSP )

    @chk ccall(
               (:KSPResetViewers, $petsc_library),
               PetscErrorCode,
               (CKSP,),
               ksp,
              )


	return nothing
end 

"""
	KSPReset(petsclib::PetscLibType,ksp::PetscKSP) 
Removes any allocated `Vec` and `Mat` from the `KSP` data structures.

Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Level: intermediate

-seealso: [](ch_ksp), `KSPCreate()`, `KSPSetUp()`, `KSPSolve()`, `KSP`

# External Links
$(_doc_external("KSP/KSPReset"))
"""
function KSPReset(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPReset(petsclib::$UnionPetscLib, ksp::PetscKSP )

    @chk ccall(
               (:KSPReset, $petsc_library),
               PetscErrorCode,
               (CKSP,),
               ksp,
              )


	return nothing
end 

"""
	KSPDestroy(petsclib::PetscLibType,ksp::PetscKSP) 
Destroys a `KSP` context.

Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Level: beginner

-seealso: [](ch_ksp), `KSPCreate()`, `KSPSetUp()`, `KSPSolve()`, `KSP`

# External Links
$(_doc_external("KSP/KSPDestroy"))
"""
function KSPDestroy(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPDestroy(petsclib::$UnionPetscLib, ksp::PetscKSP )
	ksp_ = Ref(ksp.ptr)

    @chk ccall(
               (:KSPDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{CKSP},),
               ksp_,
              )

	ksp.ptr = C_NULL

	return nothing
end 

"""
	KSPSetPCSide(petsclib::PetscLibType,ksp::PetscKSP, side::PCSide) 
Sets the preconditioning side.

Logically Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameter:
- `side` - the preconditioning side, where side is one of
-seealso: [](ch_ksp), `KSPGetPCSide()`, `KSPSetNormType()`, `KSPGetNormType()`, `KSP`, `KSPSetPreSolve()`, `KSPSetPostSolve()`

# External Links
$(_doc_external("KSP/KSPSetPCSide"))
"""
function KSPSetPCSide(petsclib::PetscLibType, ksp::PetscKSP, side::PCSide) end

@for_petsc function KSPSetPCSide(petsclib::$UnionPetscLib, ksp::PetscKSP, side::PCSide )

    @chk ccall(
               (:KSPSetPCSide, $petsc_library),
               PetscErrorCode,
               (CKSP, PCSide),
               ksp, side,
              )


	return nothing
end 

"""
	KSPGetPCSide(petsclib::PetscLibType,ksp::PetscKSP, side::PCSide) 
Gets the preconditioning side.

Not Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameter:
- `side` - the preconditioning side, where side is one of
-seealso: [](ch_ksp), `KSPSetPCSide()`, `KSP`

# External Links
$(_doc_external("KSP/KSPGetPCSide"))
"""
function KSPGetPCSide(petsclib::PetscLibType, ksp::PetscKSP, side::PCSide) end

@for_petsc function KSPGetPCSide(petsclib::$UnionPetscLib, ksp::PetscKSP, side::PCSide )

    @chk ccall(
               (:KSPGetPCSide, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{PCSide}),
               ksp, side,
              )


	return nothing
end 

"""
	rtol::PetscReal,abstol::PetscReal,dtol::PetscReal,maxits::PetscInt = KSPGetTolerances(petsclib::PetscLibType,ksp::PetscKSP) 
Gets the relative, absolute, divergence, and maximum
iteration tolerances used by the default `KSP` convergence tests.

Not Collective

Input Parameter:
- `ksp` - the Krylov subspace context

Output Parameters:
- `rtol`   - the relative convergence tolerance
- `abstol` - the absolute convergence tolerance
- `dtol`   - the divergence tolerance
- `maxits` - maximum number of iterations

Level: intermediate

-seealso: [](ch_ksp), `KSPSetTolerances()`, `KSP`, `KSPSetMinimumIterations()`, `KSPGetMinimumIterations()`

# External Links
$(_doc_external("KSP/KSPGetTolerances"))
"""
function KSPGetTolerances(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGetTolerances(petsclib::$UnionPetscLib, ksp::PetscKSP )
	rtol_ = Ref{$PetscReal}()
	abstol_ = Ref{$PetscReal}()
	dtol_ = Ref{$PetscReal}()
	maxits_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPGetTolerances, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscInt}),
               ksp, rtol_, abstol_, dtol_, maxits_,
              )

	rtol = rtol_[]
	abstol = abstol_[]
	dtol = dtol_[]
	maxits = maxits_[]

	return rtol,abstol,dtol,maxits
end 

"""
	KSPSetTolerances(petsclib::PetscLibType,ksp::PetscKSP, rtol::PetscReal, abstol::PetscReal, dtol::PetscReal, maxits::PetscInt) 
Sets the relative, absolute, divergence, and maximum
iteration tolerances used by the default `KSP` convergence testers.

Logically Collective

Input Parameters:
- `ksp`    - the Krylov subspace context
- `rtol`   - the relative convergence tolerance, relative decrease in the (possibly preconditioned) residual norm
- `abstol` - the absolute convergence tolerance   absolute size of the (possibly preconditioned) residual norm
- `dtol`   - the divergence tolerance,   amount (possibly preconditioned) residual norm can increase before `KSPConvergedDefault()` concludes that the method is diverging
- `maxits` - maximum number of iterations to use

Options Database Keys:
- `-ksp_atol <abstol>`   - Sets `abstol`
- `-ksp_rtol <rtol>`     - Sets `rtol`
- `-ksp_divtol <dtol>`   - Sets `dtol`
- `-ksp_max_it <maxits>` - Sets `maxits`

Level: intermediate

-seealso: [](ch_ksp), `KSPGetTolerances()`, `KSPConvergedDefault()`, `KSPSetConvergenceTest()`, `KSP`, `KSPSetMinimumIterations()`

# External Links
$(_doc_external("KSP/KSPSetTolerances"))
"""
function KSPSetTolerances(petsclib::PetscLibType, ksp::PetscKSP, rtol::PetscReal, abstol::PetscReal, dtol::PetscReal, maxits::PetscInt) end

@for_petsc function KSPSetTolerances(petsclib::$UnionPetscLib, ksp::PetscKSP, rtol::$PetscReal, abstol::$PetscReal, dtol::$PetscReal, maxits::$PetscInt )

    @chk ccall(
               (:KSPSetTolerances, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscReal, $PetscReal, $PetscReal, $PetscInt),
               ksp, rtol, abstol, dtol, maxits,
              )


	return nothing
end 

"""
	KSPSetMinimumIterations(petsclib::PetscLibType,ksp::PetscKSP, minit::PetscInt) 
Sets the minimum number of iterations to use, regardless of the tolerances

Logically Collective

Input Parameters:
- `ksp`   - the Krylov subspace context
- `minit` - minimum number of iterations to use

Options Database Key:
- `-ksp_min_it <minits>` - Sets `minit`

Level: intermediate

-seealso: [](ch_ksp), `KSPGetTolerances()`, `KSPConvergedDefault()`, `KSPSetConvergenceTest()`, `KSP`, `KSPSetTolerances()`, `KSPGetMinimumIterations()`

# External Links
$(_doc_external("KSP/KSPSetMinimumIterations"))
"""
function KSPSetMinimumIterations(petsclib::PetscLibType, ksp::PetscKSP, minit::PetscInt) end

@for_petsc function KSPSetMinimumIterations(petsclib::$UnionPetscLib, ksp::PetscKSP, minit::$PetscInt )

    @chk ccall(
               (:KSPSetMinimumIterations, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, minit,
              )


	return nothing
end 

"""
	minit::PetscInt = KSPGetMinimumIterations(petsclib::PetscLibType,ksp::PetscKSP) 
Gets the minimum number of iterations to use, regardless of the tolerances, that was set with `KSPSetMinimumIterations()` or `

Not Collective

Input Parameter:
- `ksp` - the Krylov subspace context

Output Parameter:
- `minit` - minimum number of iterations to use

Level: intermediate

-seealso: [](ch_ksp), `KSPGetTolerances()`, `KSPConvergedDefault()`, `KSPSetConvergenceTest()`, `KSP`, `KSPSetTolerances()`, `KSPSetMinimumIterations()`

# External Links
$(_doc_external("KSP/KSPGetMinimumIterations"))
"""
function KSPGetMinimumIterations(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGetMinimumIterations(petsclib::$UnionPetscLib, ksp::PetscKSP )
	minit_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPGetMinimumIterations, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscInt}),
               ksp, minit_,
              )

	minit = minit_[]

	return minit
end 

"""
	KSPSetInitialGuessNonzero(petsclib::PetscLibType,ksp::PetscKSP, flg::PetscBool) 
Tells the iterative solver that the
initial guess is nonzero; otherwise `KSP` assumes the initial guess
is to be zero (and thus zeros it out before solving).

Logically Collective

Input Parameters:
- `ksp` - iterative solver obtained from `KSPCreate()`
- `flg` - `PETSC_TRUE` indicates the guess is non-zero, `PETSC_FALSE` indicates the guess is zero

Options Database Key:
- `-ksp_initial_guess_nonzero <true,false>` - use nonzero initial guess

Level: beginner

-seealso: [](ch_ksp), `KSPGetInitialGuessNonzero()`, `KSPGuessSetType()`, `KSPGuessType`, `KSP`

# External Links
$(_doc_external("KSP/KSPSetInitialGuessNonzero"))
"""
function KSPSetInitialGuessNonzero(petsclib::PetscLibType, ksp::PetscKSP, flg::PetscBool) end

@for_petsc function KSPSetInitialGuessNonzero(petsclib::$UnionPetscLib, ksp::PetscKSP, flg::PetscBool )

    @chk ccall(
               (:KSPSetInitialGuessNonzero, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flg,
              )


	return nothing
end 

"""
	flag::PetscBool = KSPGetInitialGuessNonzero(petsclib::PetscLibType,ksp::PetscKSP) 
Determines whether the `KSP` solver is using
a zero initial guess.

Not Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameter:
- `flag` - `PETSC_TRUE` if guess is nonzero, else `PETSC_FALSE`

Level: intermediate

-seealso: [](ch_ksp), `KSPSetInitialGuessNonzero()`, `KSP`

# External Links
$(_doc_external("KSP/KSPGetInitialGuessNonzero"))
"""
function KSPGetInitialGuessNonzero(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGetInitialGuessNonzero(petsclib::$UnionPetscLib, ksp::PetscKSP )
	flag_ = Ref{PetscBool}()

    @chk ccall(
               (:KSPGetInitialGuessNonzero, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{PetscBool}),
               ksp, flag_,
              )

	flag = flag_[]

	return flag
end 

"""
	KSPSetErrorIfNotConverged(petsclib::PetscLibType,ksp::PetscKSP, flg::PetscBool) 
Causes `KSPSolve()` to generate an error if the solver has not converged as soon as the error is detected.

Logically Collective

Input Parameters:
- `ksp` - iterative solver obtained from `KSPCreate()`
- `flg` - `PETSC_TRUE` indicates you want the error generated

Options Database Key:
- `-ksp_error_if_not_converged <true,false>` - generate an error and stop the program

Level: intermediate

-seealso: [](ch_ksp), `KSPGetErrorIfNotConverged()`, `KSP`

# External Links
$(_doc_external("KSP/KSPSetErrorIfNotConverged"))
"""
function KSPSetErrorIfNotConverged(petsclib::PetscLibType, ksp::PetscKSP, flg::PetscBool) end

@for_petsc function KSPSetErrorIfNotConverged(petsclib::$UnionPetscLib, ksp::PetscKSP, flg::PetscBool )

    @chk ccall(
               (:KSPSetErrorIfNotConverged, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flg,
              )


	return nothing
end 

"""
	flag::PetscBool = KSPGetErrorIfNotConverged(petsclib::PetscLibType,ksp::PetscKSP) 
Will `KSPSolve()` generate an error if the solver does not converge?

Not Collective

Input Parameter:
- `ksp` - iterative solver obtained from KSPCreate()

Output Parameter:
- `flag` - `PETSC_TRUE` if it will generate an error, else `PETSC_FALSE`

Level: intermediate

-seealso: [](ch_ksp), `KSPSetErrorIfNotConverged()`, `KSP`

# External Links
$(_doc_external("KSP/KSPGetErrorIfNotConverged"))
"""
function KSPGetErrorIfNotConverged(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGetErrorIfNotConverged(petsclib::$UnionPetscLib, ksp::PetscKSP )
	flag_ = Ref{PetscBool}()

    @chk ccall(
               (:KSPGetErrorIfNotConverged, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{PetscBool}),
               ksp, flag_,
              )

	flag = flag_[]

	return flag
end 

"""
	KSPSetInitialGuessKnoll(petsclib::PetscLibType,ksp::PetscKSP, flg::PetscBool) 
Tells the iterative solver to use `PCApply()` on the right hand side vector to compute the initial guess (The Knoll trick)

Logically Collective

Input Parameters:
- `ksp` - iterative solver obtained from `KSPCreate()`
- `flg` - `PETSC_TRUE` or `PETSC_FALSE`

Level: advanced

-seealso: [](ch_ksp), `KSPGetInitialGuessKnoll()`, `KSPGuess`, `KSPSetInitialGuessNonzero()`, `KSPGetInitialGuessNonzero()`, `KSP`

# External Links
$(_doc_external("KSP/KSPSetInitialGuessKnoll"))
"""
function KSPSetInitialGuessKnoll(petsclib::PetscLibType, ksp::PetscKSP, flg::PetscBool) end

@for_petsc function KSPSetInitialGuessKnoll(petsclib::$UnionPetscLib, ksp::PetscKSP, flg::PetscBool )

    @chk ccall(
               (:KSPSetInitialGuessKnoll, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flg,
              )


	return nothing
end 

"""
	flag::PetscBool = KSPGetInitialGuessKnoll(petsclib::PetscLibType,ksp::PetscKSP) 
Determines whether the `KSP` solver is using the Knoll trick (using PCApply(pc,b,...) to compute
the initial guess

Not Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameter:
- `flag` - `PETSC_TRUE` if using Knoll trick, else `PETSC_FALSE`

Level: advanced

-seealso: [](ch_ksp), `KSPSetInitialGuessKnoll()`, `KSPSetInitialGuessNonzero()`, `KSPGetInitialGuessNonzero()`, `KSP`

# External Links
$(_doc_external("KSP/KSPGetInitialGuessKnoll"))
"""
function KSPGetInitialGuessKnoll(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGetInitialGuessKnoll(petsclib::$UnionPetscLib, ksp::PetscKSP )
	flag_ = Ref{PetscBool}()

    @chk ccall(
               (:KSPGetInitialGuessKnoll, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{PetscBool}),
               ksp, flag_,
              )

	flag = flag_[]

	return flag
end 

"""
	flg::PetscBool = KSPGetComputeSingularValues(petsclib::PetscLibType,ksp::PetscKSP) 
Gets the flag indicating whether the extreme singular
values will be calculated via a Lanczos or Arnoldi process as the linear
system is solved.

Not Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameter:
- `flg` - `PETSC_TRUE` or `PETSC_FALSE`

Options Database Key:
- `-ksp_monitor_singular_value` - Activates `KSPSetComputeSingularValues()`

Level: advanced

-seealso: [](ch_ksp), `KSPComputeExtremeSingularValues()`, `KSPMonitorSingularValue()`, `KSP`

# External Links
$(_doc_external("KSP/KSPGetComputeSingularValues"))
"""
function KSPGetComputeSingularValues(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGetComputeSingularValues(petsclib::$UnionPetscLib, ksp::PetscKSP )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:KSPGetComputeSingularValues, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{PetscBool}),
               ksp, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	KSPSetComputeSingularValues(petsclib::PetscLibType,ksp::PetscKSP, flg::PetscBool) 
Sets a flag so that the extreme singular
values will be calculated via a Lanczos or Arnoldi process as the linear
system is solved.

Logically Collective

Input Parameters:
- `ksp` - iterative solver obtained from `KSPCreate()`
- `flg` - `PETSC_TRUE` or `PETSC_FALSE`

Options Database Key:
- `-ksp_monitor_singular_value` - Activates `KSPSetComputeSingularValues()`

Level: advanced

-seealso: [](ch_ksp), `KSPComputeExtremeSingularValues()`, `KSPMonitorSingularValue()`, `KSP`, `KSPSetComputeRitz()`

# External Links
$(_doc_external("KSP/KSPSetComputeSingularValues"))
"""
function KSPSetComputeSingularValues(petsclib::PetscLibType, ksp::PetscKSP, flg::PetscBool) end

@for_petsc function KSPSetComputeSingularValues(petsclib::$UnionPetscLib, ksp::PetscKSP, flg::PetscBool )

    @chk ccall(
               (:KSPSetComputeSingularValues, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flg,
              )


	return nothing
end 

"""
	flg::PetscBool = KSPGetComputeEigenvalues(petsclib::PetscLibType,ksp::PetscKSP) 
Gets the flag indicating that the extreme eigenvalues
values will be calculated via a Lanczos or Arnoldi process as the linear
system is solved.

Not Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameter:
- `flg` - `PETSC_TRUE` or `PETSC_FALSE`

Level: advanced

-seealso: [](ch_ksp), `KSPComputeEigenvalues()`, `KSPComputeEigenvaluesExplicitly()`, `KSP`, `KSPSetComputeRitz()`

# External Links
$(_doc_external("KSP/KSPGetComputeEigenvalues"))
"""
function KSPGetComputeEigenvalues(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGetComputeEigenvalues(petsclib::$UnionPetscLib, ksp::PetscKSP )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:KSPGetComputeEigenvalues, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{PetscBool}),
               ksp, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	KSPSetComputeEigenvalues(petsclib::PetscLibType,ksp::PetscKSP, flg::PetscBool) 
Sets a flag so that the extreme eigenvalues
values will be calculated via a Lanczos or Arnoldi process as the linear
system is solved.

Logically Collective

Input Parameters:
- `ksp` - iterative solver obtained from `KSPCreate()`
- `flg` - `PETSC_TRUE` or `PETSC_FALSE`

Level: advanced

-seealso: [](ch_ksp), `KSPComputeEigenvalues()`, `KSPComputeEigenvaluesExplicitly()`, `KSP`, `KSPSetComputeRitz()`

# External Links
$(_doc_external("KSP/KSPSetComputeEigenvalues"))
"""
function KSPSetComputeEigenvalues(petsclib::PetscLibType, ksp::PetscKSP, flg::PetscBool) end

@for_petsc function KSPSetComputeEigenvalues(petsclib::$UnionPetscLib, ksp::PetscKSP, flg::PetscBool )

    @chk ccall(
               (:KSPSetComputeEigenvalues, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flg,
              )


	return nothing
end 

"""
	KSPSetComputeRitz(petsclib::PetscLibType,ksp::PetscKSP, flg::PetscBool) 
Sets a flag so that the Ritz or harmonic Ritz pairs
will be calculated via a Lanczos or Arnoldi process as the linear
system is solved.

Logically Collective

Input Parameters:
- `ksp` - iterative solver obtained from `KSPCreate()`
- `flg` - `PETSC_TRUE` or `PETSC_FALSE`

Level: advanced

-seealso: [](ch_ksp), `KSPComputeRitz()`, `KSP`, `KSPComputeEigenvalues()`, `KSPComputeExtremeSingularValues()`

# External Links
$(_doc_external("KSP/KSPSetComputeRitz"))
"""
function KSPSetComputeRitz(petsclib::PetscLibType, ksp::PetscKSP, flg::PetscBool) end

@for_petsc function KSPSetComputeRitz(petsclib::$UnionPetscLib, ksp::PetscKSP, flg::PetscBool )

    @chk ccall(
               (:KSPSetComputeRitz, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flg,
              )


	return nothing
end 

"""
	r::PetscVec = KSPGetRhs(petsclib::PetscLibType,ksp::PetscKSP) 
Gets the right
be solved.

Not Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameter:
- `r` - right-hand-side vector

Level: developer

-seealso: [](ch_ksp), `KSPGetSolution()`, `KSPSolve()`, `KSP`

# External Links
$(_doc_external("KSP/KSPGetRhs"))
"""
function KSPGetRhs(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGetRhs(petsclib::$UnionPetscLib, ksp::PetscKSP)
    r_ = Ref{CVec}()

    @chk ccall(
               (:KSPGetRhs, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{CVec}),
               ksp, r_,
              )

	return PetscVec(r_[], petsclib)
end 

"""
	v::PetscVec = KSPGetSolution(petsclib::PetscLibType,ksp::PetscKSP) 
Gets the location of the solution for the linear system to be solved.

Not Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Return Value:
- `v` - solution vector

Level: developer

-seealso: [](ch_ksp), `KSPGetRhs()`, `KSPBuildSolution()`, `KSPSolve()`, `KSP`

# External Links
$(_doc_external("KSP/KSPGetSolution"))
"""
function KSPGetSolution(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGetSolution(petsclib::$UnionPetscLib, ksp::PetscKSP )
	v_ = Ref{CVec}()

    @chk ccall(
               (:KSPGetSolution, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{CVec}),
               ksp, v_,
              )

	v =  PetscVec(v_[], petsclib)
	return v
end 

"""
	KSPSetPC(petsclib::PetscLibType,ksp::PetscKSP, pc::PC) 
Sets the preconditioner to be used to calculate the
application of the preconditioner on a vector into a `KSP`.

Collective

Input Parameters:
- `ksp` - the `KSP` iterative solver obtained from `KSPCreate()`
- `pc`  - the preconditioner object (if `NULL` it returns the `PC` currently held by the `KSP`)

Level: developer

-seealso: [](ch_ksp), `KSPGetPC()`, `KSP`

# External Links
$(_doc_external("KSP/KSPSetPC"))
"""
function KSPSetPC(petsclib::PetscLibType, ksp::PetscKSP, pc::PC) end

@for_petsc function KSPSetPC(petsclib::$UnionPetscLib, ksp::PetscKSP, pc::PC )

    @chk ccall(
               (:KSPSetPC, $petsc_library),
               PetscErrorCode,
               (CKSP, PC),
               ksp, pc,
              )


	return nothing
end 

"""
	KSPGetPC(petsclib::PetscLibType,ksp::PetscKSP, pc::PC) 
Returns a pointer to the preconditioner context with the `KSP`

Not Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameter:
- `pc` - preconditioner context

Level: beginner

-seealso: [](ch_ksp), `KSPSetPC()`, `KSP`, `PC`

# External Links
$(_doc_external("KSP/KSPGetPC"))
"""
function KSPGetPC(petsclib::PetscLibType, ksp::PetscKSP, pc::PC) end

@for_petsc function KSPGetPC(petsclib::$UnionPetscLib, ksp::PetscKSP, pc::PC )

    @chk ccall(
               (:KSPGetPC, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{PC}),
               ksp, pc,
              )


	return nothing
end 

"""
	KSPMonitor(petsclib::PetscLibType,ksp::PetscKSP, it::PetscInt, rnorm::PetscReal) 
runs the user provided monitor routines, if they exist

Collective

Input Parameters:
- `ksp`   - iterative solver obtained from `KSPCreate()`
- `it`    - iteration number
- `rnorm` - relative norm of the residual

Level: developer

-seealso: [](ch_ksp), `KSPMonitorSet()`

# External Links
$(_doc_external("KSP/KSPMonitor"))
"""
function KSPMonitor(petsclib::PetscLibType, ksp::PetscKSP, it::PetscInt, rnorm::PetscReal) end

@for_petsc function KSPMonitor(petsclib::$UnionPetscLib, ksp::PetscKSP, it::$PetscInt, rnorm::$PetscReal )

    @chk ccall(
               (:KSPMonitor, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal),
               ksp, it, rnorm,
              )


	return nothing
end 

"""
	KSPMonitorSet(petsclib::PetscLibType,ksp::PetscKSP, monitor::KSPMonitorFn, ctx::Cvoid, monitordestroy::PetscCtxDestroyFn) 
Sets an ADDITIONAL function to be called at every iteration to monitor, i.e. display in some way, perhaps by printing in the terminal,
the residual norm computed in a `KSPSolve()`

Logically Collective

Input Parameters:
- `ksp`            - iterative solver obtained from `KSPCreate()`
- `monitor`        - pointer to function (if this is `NULL`, it turns off monitoring, see `KSPMonitorFn`
- `ctx`            - [optional] context for private data for the monitor routine (use `NULL` if no context is needed)
- `monitordestroy` - [optional] routine that frees monitor context (may be `NULL`), see `PetscCtxDestroyFn` for the calling sequence

Options Database Keys:
- `-ksp_monitor`                             - sets `KSPMonitorResidual()`
- `-ksp_monitor hdf5:filename`               - sets `KSPMonitorResidualView()` and saves residual
- `-ksp_monitor draw`                        - sets `KSPMonitorResidualView()` and plots residual
- `-ksp_monitor draw::draw_lg`               - sets `KSPMonitorResidualDrawLG()` and plots residual
- `-ksp_monitor_pause_final`                 - Pauses any graphics when the solve finishes (only works for internal monitors)
- `-ksp_monitor_true_residual`               - sets `KSPMonitorTrueResidual()`
- `-ksp_monitor_true_residual draw::draw_lg` - sets `KSPMonitorTrueResidualDrawLG()` and plots residual
- `-ksp_monitor_max`                         - sets `KSPMonitorTrueResidualMax()`
- `-ksp_monitor_singular_value`              - sets `KSPMonitorSingularValue()`
- `-ksp_monitor_cancel`                      - cancels all monitors that have been hardwired into a code by calls to `KSPMonitorSet()`, but
does not cancel those set via the options database.

Level: beginner

-seealso: [](ch_ksp), `KSPMonitorResidual()`, `KSPMonitorRegister()`, `KSPMonitorCancel()`, `KSP`, `PetscCtxDestroyFn`

# External Links
$(_doc_external("KSP/KSPMonitorSet"))
"""
function KSPMonitorSet(petsclib::PetscLibType, ksp::PetscKSP, monitor::KSPMonitorFn, ctx::Cvoid, monitordestroy::PetscCtxDestroyFn) end

@for_petsc function KSPMonitorSet(petsclib::$UnionPetscLib, ksp::PetscKSP, monitor::KSPMonitorFn, ctx::Cvoid, monitordestroy::PetscCtxDestroyFn )

    @chk ccall(
               (:KSPMonitorSet, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPMonitorFn}, Ptr{Cvoid}, Ptr{PetscCtxDestroyFn}),
               ksp, monitor, ctx, monitordestroy,
              )


	return nothing
end 

"""
	KSPMonitorCancel(petsclib::PetscLibType,ksp::PetscKSP) 
Clears all monitors for a `KSP` object.

Logically Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Options Database Key:
- `-ksp_monitor_cancel` - Cancels all monitors that have been hardwired into a code by calls to `KSPMonitorSet()`, but does not cancel those set via the options database.

Level: intermediate

-seealso: [](ch_ksp), `KSPMonitorResidual()`, `KSPMonitorSet()`, `KSP`

# External Links
$(_doc_external("KSP/KSPMonitorCancel"))
"""
function KSPMonitorCancel(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPMonitorCancel(petsclib::$UnionPetscLib, ksp::PetscKSP )

    @chk ccall(
               (:KSPMonitorCancel, $petsc_library),
               PetscErrorCode,
               (CKSP,),
               ksp,
              )


	return nothing
end 

"""
	KSPGetMonitorContext(petsclib::PetscLibType,ksp::PetscKSP, ctx::Cvoid) 
Gets the monitoring context, as set by `KSPMonitorSet()` for the FIRST monitor only.

Not Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameter:
- `ctx` - monitoring context

Level: intermediate

-seealso: [](ch_ksp), `KSPMonitorResidual()`, `KSP`

# External Links
$(_doc_external("KSP/KSPGetMonitorContext"))
"""
function KSPGetMonitorContext(petsclib::PetscLibType, ksp::PetscKSP, ctx::Cvoid) end

@for_petsc function KSPGetMonitorContext(petsclib::$UnionPetscLib, ksp::PetscKSP, ctx::Cvoid )

    @chk ccall(
               (:KSPGetMonitorContext, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Cvoid}),
               ksp, ctx,
              )


	return nothing
end 

"""
	KSPSetResidualHistory(petsclib::PetscLibType,ksp::PetscKSP, a::Vector{PetscReal}, na::PetscCount, reset::PetscBool) 
Sets the array used to hold the residual history.
If set, this array will contain the residual norms computed at each
iteration of the solver.

Not Collective

Input Parameters:
- `ksp`   - iterative solver obtained from `KSPCreate()`
- `a`     - array to hold history
- `na`    - size of `a`
- `reset` - `PETSC_TRUE` indicates the history counter is reset to zero
for each new linear solve

Level: advanced

-seealso: [](ch_ksp), `KSPGetResidualHistory()`, `KSP`

# External Links
$(_doc_external("KSP/KSPSetResidualHistory"))
"""
function KSPSetResidualHistory(petsclib::PetscLibType, ksp::PetscKSP, a::Vector{PetscReal}, na::PetscCount, reset::PetscBool) end

@for_petsc function KSPSetResidualHistory(petsclib::$UnionPetscLib, ksp::PetscKSP, a::Vector{$PetscReal}, na::PetscCount, reset::PetscBool )

    @chk ccall(
               (:KSPSetResidualHistory, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscReal}, PetscCount, PetscBool),
               ksp, a, na, reset,
              )


	return nothing
end 

"""
	a::Vector{PetscReal},na::PetscInt = KSPGetResidualHistory(petsclib::PetscLibType,ksp::PetscKSP) 
Gets the array used to hold the residual history and the number of residuals it contains.

Not Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameters:
- `a`  - pointer to array to hold history (or `NULL`)
- `na` - number of used entries in a (or `NULL`). Note this has different meanings depending on the `reset` argument to `KSPSetResidualHistory()`

Level: advanced

-seealso: [](ch_ksp), `KSPSetResidualHistory()`, `KSP`, `KSPGetIterationNumber()`, `KSPSTCG`, `KSPBCGSL`

# External Links
$(_doc_external("KSP/KSPGetResidualHistory"))
"""
function KSPGetResidualHistory(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGetResidualHistory(petsclib::$UnionPetscLib, ksp::PetscKSP )
	a_ = Ref{Ptr{$PetscReal}}()
	na_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPGetResidualHistory, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Ptr{$PetscReal}}, Ptr{$PetscInt}),
               ksp, a_, na_,
              )

	a = unsafe_wrap(Array, a_[], VecGetLocalSize(petsclib, x); own = false)
	na = na_[]

	return a,na
end 

"""
	KSPSetErrorHistory(petsclib::PetscLibType,ksp::PetscKSP, a::Vector{PetscReal}, na::PetscCount, reset::PetscBool) 
Sets the array used to hold the error history. If set, this array will contain the error norms computed at each iteration of the solver.

Not Collective

Input Parameters:
- `ksp`   - iterative solver obtained from `KSPCreate()`
- `a`     - array to hold history
- `na`    - size of `a`
- `reset` - `PETSC_TRUE` indicates the history counter is reset to zero for each new linear solve

Level: advanced

-seealso: [](ch_ksp), `KSPGetErrorHistory()`, `KSPSetResidualHistory()`, `KSP`

# External Links
$(_doc_external("KSP/KSPSetErrorHistory"))
"""
function KSPSetErrorHistory(petsclib::PetscLibType, ksp::PetscKSP, a::Vector{PetscReal}, na::PetscCount, reset::PetscBool) end

@for_petsc function KSPSetErrorHistory(petsclib::$UnionPetscLib, ksp::PetscKSP, a::Vector{$PetscReal}, na::PetscCount, reset::PetscBool )

    @chk ccall(
               (:KSPSetErrorHistory, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscReal}, PetscCount, PetscBool),
               ksp, a, na, reset,
              )


	return nothing
end 

"""
	a::Vector{PetscReal},na::PetscInt = KSPGetErrorHistory(petsclib::PetscLibType,ksp::PetscKSP) 
Gets the array used to hold the error history and the number of residuals it contains.

Not Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameters:
- `a`  - pointer to array to hold history (or `NULL`)
- `na` - number of used entries in a (or `NULL`)

Level: advanced

-seealso: [](ch_ksp), `KSPSetErrorHistory()`, `KSPGetResidualHistory()`, `KSP`

# External Links
$(_doc_external("KSP/KSPGetErrorHistory"))
"""
function KSPGetErrorHistory(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGetErrorHistory(petsclib::$UnionPetscLib, ksp::PetscKSP )
	a_ = Ref{Ptr{$PetscReal}}()
	na_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPGetErrorHistory, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Ptr{$PetscReal}}, Ptr{$PetscInt}),
               ksp, a_, na_,
              )

	a = unsafe_wrap(Array, a_[], VecGetLocalSize(petsclib, x); own = false)
	na = na_[]

	return a,na
end 

"""
	cr::PetscReal,rRsq::PetscReal,ce::PetscReal,eRsq::PetscReal = KSPComputeConvergenceRate(petsclib::PetscLibType,ksp::PetscKSP) 
Compute the convergence rate for the iteration <https:/en.wikipedia.org/wiki/Coefficient_of_determination>

Not Collective

Input Parameter:
- `ksp` - The `KSP`

Output Parameters:
- `cr`   - The residual contraction rate
- `rRsq` - The coefficient of determination, R^2, indicating the linearity of the data
- `ce`   - The error contraction rate
- `eRsq` - The coefficient of determination, R^2, indicating the linearity of the data

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPConvergedRateView()`

# External Links
$(_doc_external("KSP/KSPComputeConvergenceRate"))
"""
function KSPComputeConvergenceRate(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPComputeConvergenceRate(petsclib::$UnionPetscLib, ksp::PetscKSP )
	cr_ = Ref{$PetscReal}()
	rRsq_ = Ref{$PetscReal}()
	ce_ = Ref{$PetscReal}()
	eRsq_ = Ref{$PetscReal}()

    @chk ccall(
               (:KSPComputeConvergenceRate, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               ksp, cr_, rRsq_, ce_, eRsq_,
              )

	cr = cr_[]
	rRsq = rRsq_[]
	ce = ce_[]
	eRsq = eRsq_[]

	return cr,rRsq,ce,eRsq
end 

"""
	KSPSetConvergenceTest(petsclib::PetscLibType,ksp::PetscKSP, converge::KSPConvergenceTestFn, ctx::Cvoid, destroy::PetscCtxDestroyFn) 
Sets the function to be used to determine convergence of `KSPSolve()`

Logically Collective

Input Parameters:
- `ksp`      - iterative solver obtained from `KSPCreate()`
- `converge` - pointer to the function, see `KSPConvergenceTestFn`
- `ctx`      - context for private data for the convergence routine (may be `NULL`)
- `destroy`  - a routine for destroying the context (may be `NULL`)

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPConvergenceTestFn`, `KSPConvergedDefault()`, `KSPGetConvergenceContext()`, `KSPSetTolerances()`, `KSPGetConvergenceTest()`, `KSPGetAndClearConvergenceTest()`

# External Links
$(_doc_external("KSP/KSPSetConvergenceTest"))
"""
function KSPSetConvergenceTest(petsclib::PetscLibType, ksp::PetscKSP, converge::KSPConvergenceTestFn, ctx::Cvoid, destroy::PetscCtxDestroyFn) end

@for_petsc function KSPSetConvergenceTest(petsclib::$UnionPetscLib, ksp::PetscKSP, converge::KSPConvergenceTestFn, ctx::Cvoid, destroy::PetscCtxDestroyFn )

    @chk ccall(
               (:KSPSetConvergenceTest, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPConvergenceTestFn}, Ptr{Cvoid}, Ptr{PetscCtxDestroyFn}),
               ksp, converge, ctx, destroy,
              )


	return nothing
end 

"""
	KSPGetConvergenceTest(petsclib::PetscLibType,ksp::PetscKSP, converge::KSPConvergenceTestFn, ctx::Cvoid, destroy::PetscCtxDestroyFn) 
Gets the function to be used to determine convergence.

Logically Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameters:
- `converge` - pointer to convergence test function, see `KSPConvergenceTestFn`
- `ctx`      - context for private data for the convergence routine (may be `NULL`)
- `destroy`  - a routine for destroying the context (may be `NULL`)

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPConvergedDefault()`, `KSPGetConvergenceContext()`, `KSPSetTolerances()`, `KSPSetConvergenceTest()`, `KSPGetAndClearConvergenceTest()`

# External Links
$(_doc_external("KSP/KSPGetConvergenceTest"))
"""
function KSPGetConvergenceTest(petsclib::PetscLibType, ksp::PetscKSP, converge::KSPConvergenceTestFn, ctx::Cvoid, destroy::PetscCtxDestroyFn) end

@for_petsc function KSPGetConvergenceTest(petsclib::$UnionPetscLib, ksp::PetscKSP, converge::KSPConvergenceTestFn, ctx::Cvoid, destroy::PetscCtxDestroyFn )

    @chk ccall(
               (:KSPGetConvergenceTest, $petsc_library),
               PetscErrorCode,
               (CKSP, KSPConvergenceTestFn, Cvoid, PetscCtxDestroyFn),
               ksp, converge, ctx, destroy,
              )


	return nothing
end 

"""
	KSPGetAndClearConvergenceTest(petsclib::PetscLibType,ksp::PetscKSP, converge::KSPConvergenceTestFn, ctx::Cvoid, destroy::PetscCtxDestroyFn) 
Gets the function to be used to determine convergence. Removes the current test without calling destroy on the test context

Logically Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameters:
- `converge` - pointer to convergence test function, see `KSPConvergenceTestFn`
- `ctx`      - context for private data for the convergence routine
- `destroy`  - a routine for destroying the context

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPConvergedDefault()`, `KSPGetConvergenceContext()`, `KSPSetTolerances()`, `KSPSetConvergenceTest()`, `KSPGetConvergenceTest()`

# External Links
$(_doc_external("KSP/KSPGetAndClearConvergenceTest"))
"""
function KSPGetAndClearConvergenceTest(petsclib::PetscLibType, ksp::PetscKSP, converge::KSPConvergenceTestFn, ctx::Cvoid, destroy::PetscCtxDestroyFn) end

@for_petsc function KSPGetAndClearConvergenceTest(petsclib::$UnionPetscLib, ksp::PetscKSP, converge::KSPConvergenceTestFn, ctx::Cvoid, destroy::PetscCtxDestroyFn )

    @chk ccall(
               (:KSPGetAndClearConvergenceTest, $petsc_library),
               PetscErrorCode,
               (CKSP, KSPConvergenceTestFn, Cvoid, PetscCtxDestroyFn),
               ksp, converge, ctx, destroy,
              )


	return nothing
end 

"""
	KSPGetConvergenceContext(petsclib::PetscLibType,ksp::PetscKSP, ctx::Cvoid) 
Gets the convergence context set with `KSPSetConvergenceTest()`.

Not Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameter:
- `ctx` - monitoring context

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPConvergedDefault()`, `KSPSetConvergenceTest()`, `KSPGetConvergenceTest()`

# External Links
$(_doc_external("KSP/KSPGetConvergenceContext"))
"""
function KSPGetConvergenceContext(petsclib::PetscLibType, ksp::PetscKSP, ctx::Cvoid) end

@for_petsc function KSPGetConvergenceContext(petsclib::$UnionPetscLib, ksp::PetscKSP, ctx::Cvoid )

    @chk ccall(
               (:KSPGetConvergenceContext, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Cvoid}),
               ksp, ctx,
              )


	return nothing
end 

"""
	KSPBuildSolution(petsclib::PetscLibType,ksp::PetscKSP, v::PetscVec, V::PetscVec) 
Builds the approximate solution in a vector provided.

Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameter:
Provide exactly one of
- `v` - location to stash solution, optional, otherwise pass `NULL`
- `V` - the solution is returned in this location. This vector is created internally. This vector should NOT be destroyed by the user with `VecDestroy()`.

Level: developer

-seealso: [](ch_ksp), `KSPGetSolution()`, `KSPBuildResidual()`, `KSP`

# External Links
$(_doc_external("KSP/KSPBuildSolution"))
"""
function KSPBuildSolution(petsclib::PetscLibType, ksp::PetscKSP, v::PetscVec, V::PetscVec) end

@for_petsc function KSPBuildSolution(petsclib::$UnionPetscLib, ksp::PetscKSP, v::PetscVec, V::PetscVec )
	V_ = Ref(V.ptr)

    @chk ccall(
               (:KSPBuildSolution, $petsc_library),
               PetscErrorCode,
               (CKSP, CVec, Ptr{CVec}),
               ksp, v, V_,
              )

	V.ptr = C_NULL

	return nothing
end 

"""
	KSPBuildResidual(petsclib::PetscLibType,ksp::PetscKSP, t::PetscVec, v::PetscVec, V::PetscVec) 
Builds the residual in a vector provided.

Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameters:
- `t` - work vector.  If not provided then one is generated.
- `v` - optional location to stash residual.  If `v` is not provided, then a location is generated.
- `V` - the residual

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPBuildSolution()`

# External Links
$(_doc_external("KSP/KSPBuildResidual"))
"""
function KSPBuildResidual(petsclib::PetscLibType, ksp::PetscKSP, t::PetscVec, v::PetscVec, V::PetscVec) end

@for_petsc function KSPBuildResidual(petsclib::$UnionPetscLib, ksp::PetscKSP, t::PetscVec, v::PetscVec, V::PetscVec )
	V_ = Ref(V.ptr)

    @chk ccall(
               (:KSPBuildResidual, $petsc_library),
               PetscErrorCode,
               (CKSP, CVec, CVec, Ptr{CVec}),
               ksp, t, v, V_,
              )

	V.ptr = C_NULL

	return nothing
end 

"""
	KSPSetDiagonalScale(petsclib::PetscLibType,ksp::PetscKSP, scale::PetscBool) 
Tells `KSP` to symmetrically diagonally scale the system
before solving. This actually CHANGES the matrix (and right-hand side).

Logically Collective

Input Parameters:
- `ksp`   - the `KSP` context
- `scale` - `PETSC_TRUE` or `PETSC_FALSE`

Options Database Keys:
- `-ksp_diagonal_scale`     - perform a diagonal scaling before the solve
- `-ksp_diagonal_scale_fix` - scale the matrix back AFTER the solve

Level: advanced

-seealso: [](ch_ksp), `KSPGetDiagonalScale()`, `KSPSetDiagonalScaleFix()`, `KSP`

# External Links
$(_doc_external("KSP/KSPSetDiagonalScale"))
"""
function KSPSetDiagonalScale(petsclib::PetscLibType, ksp::PetscKSP, scale::PetscBool) end

@for_petsc function KSPSetDiagonalScale(petsclib::$UnionPetscLib, ksp::PetscKSP, scale::PetscBool )

    @chk ccall(
               (:KSPSetDiagonalScale, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, scale,
              )


	return nothing
end 

"""
	scale::PetscBool = KSPGetDiagonalScale(petsclib::PetscLibType,ksp::PetscKSP) 
Checks if `KSP` solver scales the matrix and right

Not Collective

Input Parameter:
- `ksp` - the `KSP` context

Output Parameter:
- `scale` - `PETSC_TRUE` or `PETSC_FALSE`

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPSetDiagonalScale()`, `KSPSetDiagonalScaleFix()`

# External Links
$(_doc_external("KSP/KSPGetDiagonalScale"))
"""
function KSPGetDiagonalScale(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGetDiagonalScale(petsclib::$UnionPetscLib, ksp::PetscKSP )
	scale_ = Ref{PetscBool}()

    @chk ccall(
               (:KSPGetDiagonalScale, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{PetscBool}),
               ksp, scale_,
              )

	scale = scale_[]

	return scale
end 

"""
	KSPSetDiagonalScaleFix(petsclib::PetscLibType,ksp::PetscKSP, fix::PetscBool) 
Tells `KSP` to diagonally scale the system back after solving.

Logically Collective

Input Parameters:
- `ksp` - the `KSP` context
- `fix` - `PETSC_TRUE` to scale back after the system solve, `PETSC_FALSE` to not
rescale (default)

Level: intermediate

-seealso: [](ch_ksp), `KSPGetDiagonalScale()`, `KSPSetDiagonalScale()`, `KSPGetDiagonalScaleFix()`, `KSP`

# External Links
$(_doc_external("KSP/KSPSetDiagonalScaleFix"))
"""
function KSPSetDiagonalScaleFix(petsclib::PetscLibType, ksp::PetscKSP, fix::PetscBool) end

@for_petsc function KSPSetDiagonalScaleFix(petsclib::$UnionPetscLib, ksp::PetscKSP, fix::PetscBool )

    @chk ccall(
               (:KSPSetDiagonalScaleFix, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, fix,
              )


	return nothing
end 

"""
	fix::PetscBool = KSPGetDiagonalScaleFix(petsclib::PetscLibType,ksp::PetscKSP) 
Determines if `KSP` diagonally scales the system back after solving. That is `KSPSetDiagonalScaleFix()` has been called

Not Collective

Input Parameter:
- `ksp` - the `KSP` context

Output Parameter:
- `fix` - `PETSC_TRUE` to scale back after the system solve, `PETSC_FALSE` to not
rescale (default)

Level: intermediate

-seealso: [](ch_ksp), `KSPGetDiagonalScale()`, `KSPSetDiagonalScale()`, `KSPSetDiagonalScaleFix()`, `KSP`

# External Links
$(_doc_external("KSP/KSPGetDiagonalScaleFix"))
"""
function KSPGetDiagonalScaleFix(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGetDiagonalScaleFix(petsclib::$UnionPetscLib, ksp::PetscKSP )
	fix_ = Ref{PetscBool}()

    @chk ccall(
               (:KSPGetDiagonalScaleFix, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{PetscBool}),
               ksp, fix_,
              )

	fix = fix_[]

	return fix
end 

"""
	KSPSetComputeOperators(petsclib::PetscLibType,ksp::PetscKSP, func::Union{KSPComputeOperatorsFn, Ptr}, ctx::Union{Cvoid, Ptr}) 
set routine to compute the linear operators

Logically Collective

Input Parameters:
- `ksp`  - the `KSP` context
- `func` - function to compute the operators, see `KSPComputeOperatorsFn` for the calling sequence
- `ctx`  - optional context

Level: beginner

-seealso: [](ch_ksp), `KSP`, `KSPSetOperators()`, `KSPSetComputeRHS()`, `DMKSPSetComputeOperators()`, `KSPSetComputeInitialGuess()`, `KSPComputeOperatorsFn`

# External Links
$(_doc_external("KSP/KSPSetComputeOperators"))
"""
function KSPSetComputeOperators(petsclib::PetscLibType, ksp::PetscKSP, func::Union{KSPComputeOperatorsFn, Ptr}, ctx::Union{Cvoid, Ptr}) end

@for_petsc function KSPSetComputeOperators(petsclib::$UnionPetscLib, ksp::PetscKSP, func::Union{KSPComputeOperatorsFn, Ptr}, ctx::Union{Cvoid, Ptr})

    @chk ccall(
               (:KSPSetComputeOperators, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPComputeOperatorsFn}, Ptr{Cvoid}),
               ksp, func, ctx,
              )


	return nothing
end 

"""
	KSPSetComputeRHS(petsclib::PetscLibType,ksp::PetscKSP, func::Union{KSPComputeRHSFn,Ptr}, ctx::Union{Cvoid, Ptr}) 
set routine to compute the right

Logically Collective

Input Parameters:
- `ksp`  - the `KSP` context
- `func` - function to compute the right-hand side, see `KSPComputeRHSFn` for the calling sequence
- `ctx`  - optional context

Level: beginner

-seealso: [](ch_ksp), `KSP`, `KSPSolve()`, `DMKSPSetComputeRHS()`, `KSPSetComputeOperators()`, `KSPSetOperators()`, `KSPComputeRHSFn`

# External Links
$(_doc_external("KSP/KSPSetComputeRHS"))
"""
function KSPSetComputeRHS(petsclib::PetscLibType, ksp::PetscKSP, func::Union{KSPComputeRHSFn,Ptr}, ctx::Union{Cvoid, Ptr}) end

@for_petsc function KSPSetComputeRHS(petsclib::$UnionPetscLib, ksp::PetscKSP, func::Union{KSPComputeRHSFn,Ptr}, ctx::Union{Cvoid, Ptr} )

    @chk ccall(
               (:KSPSetComputeRHS, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPComputeRHSFn}, Ptr{Cvoid}),
               ksp, func, ctx,
              )


	return nothing
end 

"""
	KSPSetComputeInitialGuess(petsclib::PetscLibType,ksp::PetscKSP, func::Union{KSPComputeInitialGuessFn, Ptr}, ctx::Union{Cvoid, Ptr}) 
set routine to compute the initial guess of the linear system

Logically Collective

Input Parameters:
- `ksp`  - the `KSP` context
- `func` - function to compute the initial guess, see `KSPComputeInitialGuessFn` for calling sequence
- `ctx`  - optional context

Level: beginner

-seealso: [](ch_ksp), `KSP`, `KSPSolve()`, `KSPSetComputeRHS()`, `KSPSetComputeOperators()`, `DMKSPSetComputeInitialGuess()`, `KSPSetInitialGuessNonzero()`,
`KSPComputeInitialGuessFn`

# External Links
$(_doc_external("KSP/KSPSetComputeInitialGuess"))
"""
function KSPSetComputeInitialGuess(petsclib::PetscLibType, ksp::PetscKSP, func::Union{KSPComputeInitialGuessFn, Ptr}, ctx::Union{Cvoid, Ptr}) end

@for_petsc function KSPSetComputeInitialGuess(petsclib::$UnionPetscLib, ksp::PetscKSP, func::Union{KSPComputeInitialGuessFn, Ptr}, ctx::Union{Cvoid, Ptr})

    @chk ccall(
               (:KSPSetComputeInitialGuess, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPComputeInitialGuessFn}, Ptr{Cvoid}),
               ksp, func, ctx,
              )


	return nothing
end 

"""
	flg::PetscBool = KSPSetUseExplicitTranspose(petsclib::PetscLibType,ksp::PetscKSP) 
Determines the explicit transpose of the operator is formed in `KSPSolveTranspose()`. In some configurations (like GPUs) it may
be explicitly formed since the solve is much more efficient.

Logically Collective

Input Parameter:
- `ksp` - the `KSP` context

Output Parameter:
- `flg` - `PETSC_TRUE` to transpose the system in `KSPSolveTranspose()`, `PETSC_FALSE` to not transpose (default)

Level: advanced

-seealso: [](ch_ksp), `KSPSolveTranspose()`, `KSP`

# External Links
$(_doc_external("KSP/KSPSetUseExplicitTranspose"))
"""
function KSPSetUseExplicitTranspose(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPSetUseExplicitTranspose(petsclib::$UnionPetscLib, ksp::PetscKSP )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:KSPSetUseExplicitTranspose, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	rnorm::PetscReal = KSPGetResidualNorm(petsclib::PetscLibType,ksp::PetscKSP) 
Gets the last (possibly approximate and/or preconditioned) residual norm that has been computed.

Not Collective

Input Parameter:
- `ksp` - the iterative context

Output Parameter:
- `rnorm` - residual norm

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPSetNormType()`, `KSPBuildResidual()`, `KSPNormType`

# External Links
$(_doc_external("KSP/KSPGetResidualNorm"))
"""
function KSPGetResidualNorm(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGetResidualNorm(petsclib::$UnionPetscLib, ksp::PetscKSP )
	rnorm_ = Ref{$PetscReal}()

    @chk ccall(
               (:KSPGetResidualNorm, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscReal}),
               ksp, rnorm_,
              )

	rnorm = rnorm_[]

	return rnorm
end 

"""
	its::PetscInt = KSPGetIterationNumber(petsclib::PetscLibType,ksp::PetscKSP) 
Gets the current iteration number; if the `KSPSolve()` is complete, returns the number of iterations used.

Not Collective

Input Parameter:
- `ksp` - the iterative context

Output Parameter:
- `its` - number of iterations

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPGetResidualNorm()`, `KSPBuildResidual()`, `KSPGetTotalIterations()`

# External Links
$(_doc_external("KSP/KSPGetIterationNumber"))
"""
function KSPGetIterationNumber(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGetIterationNumber(petsclib::$UnionPetscLib, ksp::PetscKSP )
	its_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPGetIterationNumber, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscInt}),
               ksp, its_,
              )

	its = its_[]

	return its
end 

"""
	its::PetscInt = KSPGetTotalIterations(petsclib::PetscLibType,ksp::PetscKSP) 
Gets the total number of iterations this `KSP` object has performed since was created, counted over all linear solves

Not Collective

Input Parameter:
- `ksp` - the iterative context

Output Parameter:
- `its` - total number of iterations

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPBuildResidual()`, `KSPGetResidualNorm()`, `KSPGetIterationNumber()`

# External Links
$(_doc_external("KSP/KSPGetTotalIterations"))
"""
function KSPGetTotalIterations(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGetTotalIterations(petsclib::$UnionPetscLib, ksp::PetscKSP )
	its_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPGetTotalIterations, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscInt}),
               ksp, its_,
              )

	its = its_[]

	return its
end 

"""
	KSPMonitorResidual(petsclib::PetscLibType,ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) 
Print the (possibly preconditioned, possibly approximate) residual norm at each iteration of an iterative solver.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - (preconditioned) residual norm value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_monitor` - Activates `KSPMonitorResidual()` to print the norm value at each iteration

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPMonitorResidualView()`, `KSPMonitorResidualDrawLG()`,
`KSPMonitorResidualRange()`, `KSPMonitorTrueResidualDraw()`, `KSPMonitorTrueResidualDrawLG()`, `KSPMonitorTrueResidualMax()`,
`KSPMonitorSingularValue()`, `KSPMonitorSolutionDrawLG()`, `KSPMonitorSolutionDraw()`, `KSPMonitorSolution()`,
`KSPMonitorErrorDrawLG()`, `KSPMonitorErrorDraw()`, `KSPMonitorError()`

# External Links
$(_doc_external("KSP/KSPMonitorResidual"))
"""
function KSPMonitorResidual(petsclib::PetscLibType, ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function KSPMonitorResidual(petsclib::$UnionPetscLib, ksp::PetscKSP, n::$PetscInt, rnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:KSPMonitorResidual, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	KSPMonitorResidualView(petsclib::PetscLibType,ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) 
Plots the (possibly preconditioned) residual at each iteration of an iterative solver.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_monitor viewertype` - Activates `KSPMonitorResidualView()`

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPMonitorResidual()`, `KSPMonitorResidualDrawLG()`

# External Links
$(_doc_external("KSP/KSPMonitorResidualView"))
"""
function KSPMonitorResidualView(petsclib::PetscLibType, ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function KSPMonitorResidualView(petsclib::$UnionPetscLib, ksp::PetscKSP, n::$PetscInt, rnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:KSPMonitorResidualView, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	KSPMonitorResidualDrawLG(petsclib::PetscLibType,ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) 
Plots the (possibly preconditioned) residual norm at each iteration of an iterative solver.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_monitor draw::draw_lg` - Activates `KSPMonitorResidualDrawLG()`

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `PETSCVIEWERDRAW`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPMonitorResidualView()`, `KSPMonitorResidual()`

# External Links
$(_doc_external("KSP/KSPMonitorResidualDrawLG"))
"""
function KSPMonitorResidualDrawLG(petsclib::PetscLibType, ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function KSPMonitorResidualDrawLG(petsclib::$UnionPetscLib, ksp::PetscKSP, n::$PetscInt, rnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:KSPMonitorResidualDrawLG, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	vf::PetscViewerAndFormat = KSPMonitorResidualDrawLGCreate(petsclib::PetscLibType,viewer::PetscViewer, format::PetscViewerFormat, ctx::Cvoid) 
Creates the context for the (possibly preconditioned) residual norm monitor `KSPMonitorResidualDrawLG()`

Collective

Input Parameters:
- `viewer` - The `PetscViewer` of type `PETSCVIEWERDRAW`
- `format` - The viewer format
- `ctx`    - An optional user context

Output Parameter:
- `vf` - The viewer context

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `PETSCVIEWERDRAW`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPMonitorResidualDrawLG()`,
`PetscViewerFormat`, `PetscViewer`, `PetscViewerAndFormat`

# External Links
$(_doc_external("KSP/KSPMonitorResidualDrawLGCreate"))
"""
function KSPMonitorResidualDrawLGCreate(petsclib::PetscLibType, viewer::PetscViewer, format::PetscViewerFormat, ctx::Cvoid) end

@for_petsc function KSPMonitorResidualDrawLGCreate(petsclib::$UnionPetscLib, viewer::PetscViewer, format::PetscViewerFormat, ctx::Cvoid )
	vf_ = Ref{PetscViewerAndFormat}()

    @chk ccall(
               (:KSPMonitorResidualDrawLGCreate, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscViewerFormat, Ptr{Cvoid}, PetscViewerAndFormat),
               viewer, format, ctx, vf_,
              )

	vf = vf_[]

	return vf
end 

"""
	KSPMonitorResidualShort(petsclib::PetscLibType,ksp::PetscKSP, its::PetscInt, fnorm::PetscReal, vf::PetscViewerAndFormat) 

# External Links
$(_doc_external("KSP/KSPMonitorResidualShort"))
"""
function KSPMonitorResidualShort(petsclib::PetscLibType, ksp::PetscKSP, its::PetscInt, fnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function KSPMonitorResidualShort(petsclib::$UnionPetscLib, ksp::PetscKSP, its::$PetscInt, fnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:KSPMonitorResidualShort, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, its, fnorm, vf,
              )


	return nothing
end 

"""
	KSPMonitorResidualRange(petsclib::PetscLibType,ksp::PetscKSP, it::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) 
Prints the percentage of residual elements that are more than 10 percent of the maximum value.

Collective

Input Parameters:
- `ksp`   - iterative context
- `it`    - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_monitor_range` - Activates `KSPMonitorResidualRange()`

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPMonitorSet()`, `KSPMonitorResidual()`

# External Links
$(_doc_external("KSP/KSPMonitorResidualRange"))
"""
function KSPMonitorResidualRange(petsclib::PetscLibType, ksp::PetscKSP, it::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function KSPMonitorResidualRange(petsclib::$UnionPetscLib, ksp::PetscKSP, it::$PetscInt, rnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:KSPMonitorResidualRange, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, it, rnorm, vf,
              )


	return nothing
end 

"""
	KSPMonitorTrueResidual(petsclib::PetscLibType,ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) 
Prints the true residual norm, as well as the (possibly preconditioned, possibly approximate) residual norm,
at each iteration of a `KSPSolve()` iterative solver.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_monitor_true_residual` - Activates `KSPMonitorTrueResidual()` to print both norm values at each iteration

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPMonitorSet()`, `KSPMonitorResidual()`, `KSPMonitorTrueResidualMaxNorm()`, `PetscViewerAndFormat`

# External Links
$(_doc_external("KSP/KSPMonitorTrueResidual"))
"""
function KSPMonitorTrueResidual(petsclib::PetscLibType, ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function KSPMonitorTrueResidual(petsclib::$UnionPetscLib, ksp::PetscKSP, n::$PetscInt, rnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:KSPMonitorTrueResidual, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	KSPMonitorTrueResidualView(petsclib::PetscLibType,ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) 
Plots the true residual at each iteration of an iterative solver.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context of type `PETSCVIEWERDRAW`

Options Database Key:
- `-ksp_monitor_true_residual viewertype` - Activates `KSPMonitorTrueResidualView()`

Level: intermediate

-seealso: [](ch_ksp), `PETSCVIEWERDRAW`, `KSP`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPMonitorResidual()`,
`KSPMonitorTrueResidualDrawLG()`, `PetscViewerAndFormat`

# External Links
$(_doc_external("KSP/KSPMonitorTrueResidualView"))
"""
function KSPMonitorTrueResidualView(petsclib::PetscLibType, ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function KSPMonitorTrueResidualView(petsclib::$UnionPetscLib, ksp::PetscKSP, n::$PetscInt, rnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:KSPMonitorTrueResidualView, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	KSPMonitorTrueResidualDrawLG(petsclib::PetscLibType,ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) 
Plots the true residual norm at each iteration of an iterative solver.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_monitor_true_residual draw::draw_lg` - Activates `KSPMonitorTrueResidualDrawLG()`

Level: intermediate

-seealso: [](ch_ksp), `PETSCVIEWERDRAW`, `KSP`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPMonitorTrueResidualDraw()`, `KSPMonitorResidual`,
`KSPMonitorTrueResidualDrawLGCreate()`

# External Links
$(_doc_external("KSP/KSPMonitorTrueResidualDrawLG"))
"""
function KSPMonitorTrueResidualDrawLG(petsclib::PetscLibType, ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function KSPMonitorTrueResidualDrawLG(petsclib::$UnionPetscLib, ksp::PetscKSP, n::$PetscInt, rnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:KSPMonitorTrueResidualDrawLG, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	vf::PetscViewerAndFormat = KSPMonitorTrueResidualDrawLGCreate(petsclib::PetscLibType,viewer::PetscViewer, format::PetscViewerFormat, ctx::Cvoid) 
Creates the context for the true residual monitor `KSPMonitorTrueResidualDrawLG()`

Collective

Input Parameters:
- `viewer` - The `PetscViewer` of type `PETSCVIEWERDRAW`
- `format` - The viewer format
- `ctx`    - An optional user context

Output Parameter:
- `vf` - The viewer context

Level: intermediate

-seealso: [](ch_ksp), `PETSCVIEWERDRAW`, `KSP`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `PetscViewerAndFormat`

# External Links
$(_doc_external("KSP/KSPMonitorTrueResidualDrawLGCreate"))
"""
function KSPMonitorTrueResidualDrawLGCreate(petsclib::PetscLibType, viewer::PetscViewer, format::PetscViewerFormat, ctx::Cvoid) end

@for_petsc function KSPMonitorTrueResidualDrawLGCreate(petsclib::$UnionPetscLib, viewer::PetscViewer, format::PetscViewerFormat, ctx::Cvoid )
	vf_ = Ref{PetscViewerAndFormat}()

    @chk ccall(
               (:KSPMonitorTrueResidualDrawLGCreate, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscViewerFormat, Ptr{Cvoid}, PetscViewerAndFormat),
               viewer, format, ctx, vf_,
              )

	vf = vf_[]

	return vf
end 

"""
	KSPMonitorTrueResidualMax(petsclib::PetscLibType,ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) 
Prints the true residual max norm at each iteration of an iterative solver.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_monitor_true_residual_max` - Activates `KSPMonitorTrueResidualMax()`

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPMonitorSet()`, `KSPMonitorResidual()`, `KSPMonitorTrueResidualMaxNorm()`

# External Links
$(_doc_external("KSP/KSPMonitorTrueResidualMax"))
"""
function KSPMonitorTrueResidualMax(petsclib::PetscLibType, ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function KSPMonitorTrueResidualMax(petsclib::$UnionPetscLib, ksp::PetscKSP, n::$PetscInt, rnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:KSPMonitorTrueResidualMax, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	KSPMonitorError(petsclib::PetscLibType,ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) 
Prints the error norm, as well as the (possibly preconditioned) residual norm, at each iteration of an iterative solver.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_monitor_error` - Activates `KSPMonitorError()`

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPMonitorSet()`, `KSPMonitorResidual()`, `KSPMonitorTrueResidualMaxNorm()`

# External Links
$(_doc_external("KSP/KSPMonitorError"))
"""
function KSPMonitorError(petsclib::PetscLibType, ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function KSPMonitorError(petsclib::$UnionPetscLib, ksp::PetscKSP, n::$PetscInt, rnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:KSPMonitorError, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	KSPMonitorErrorDraw(petsclib::PetscLibType,ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) 
Plots the error at each iteration of an iterative solver.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_monitor_error draw` - Activates `KSPMonitorErrorDraw()`

Level: intermediate

-seealso: [](ch_ksp), `PETSCVIEWERDRAW`, `KSP`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPMonitorErrorDrawLG()`

# External Links
$(_doc_external("KSP/KSPMonitorErrorDraw"))
"""
function KSPMonitorErrorDraw(petsclib::PetscLibType, ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function KSPMonitorErrorDraw(petsclib::$UnionPetscLib, ksp::PetscKSP, n::$PetscInt, rnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:KSPMonitorErrorDraw, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	KSPMonitorErrorDrawLG(petsclib::PetscLibType,ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) 
Plots the error and residual norm at each iteration of an iterative solver.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_monitor_error draw::draw_lg` - Activates `KSPMonitorTrueResidualDrawLG()`

Level: intermediate

-seealso: [](ch_ksp), `PETSCVIEWERDRAW`, `KSP`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPMonitorErrorDraw()`

# External Links
$(_doc_external("KSP/KSPMonitorErrorDrawLG"))
"""
function KSPMonitorErrorDrawLG(petsclib::PetscLibType, ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function KSPMonitorErrorDrawLG(petsclib::$UnionPetscLib, ksp::PetscKSP, n::$PetscInt, rnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:KSPMonitorErrorDrawLG, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	vf::PetscViewerAndFormat = KSPMonitorErrorDrawLGCreate(petsclib::PetscLibType,viewer::PetscViewer, format::PetscViewerFormat, ctx::Cvoid) 
Creates the context for the error and preconditioned residual plotter `KSPMonitorErrorDrawLG()`

Collective

Input Parameters:
- `viewer` - The `PetscViewer`
- `format` - The viewer format
- `ctx`    - An optional user context

Output Parameter:
- `vf` - The viewer context

Level: intermediate

-seealso: [](ch_ksp), `PETSCVIEWERDRAW`, `KSP`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPMonitorErrorDrawLG()`

# External Links
$(_doc_external("KSP/KSPMonitorErrorDrawLGCreate"))
"""
function KSPMonitorErrorDrawLGCreate(petsclib::PetscLibType, viewer::PetscViewer, format::PetscViewerFormat, ctx::Cvoid) end

@for_petsc function KSPMonitorErrorDrawLGCreate(petsclib::$UnionPetscLib, viewer::PetscViewer, format::PetscViewerFormat, ctx::Cvoid )
	vf_ = Ref{PetscViewerAndFormat}()

    @chk ccall(
               (:KSPMonitorErrorDrawLGCreate, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscViewerFormat, Ptr{Cvoid}, PetscViewerAndFormat),
               viewer, format, ctx, vf_,
              )

	vf = vf_[]

	return vf
end 

"""
	KSPMonitorSolution(petsclib::PetscLibType,ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) 
Print the solution norm at each iteration of an iterative solver.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_monitor_solution` - Activates `KSPMonitorSolution()`

Level: intermediate

-seealso: [](ch_ksp), `KSPMonitorSet()`, `KSPMonitorTrueResidual()`

# External Links
$(_doc_external("KSP/KSPMonitorSolution"))
"""
function KSPMonitorSolution(petsclib::PetscLibType, ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function KSPMonitorSolution(petsclib::$UnionPetscLib, ksp::PetscKSP, n::$PetscInt, rnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:KSPMonitorSolution, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	KSPMonitorSolutionDraw(petsclib::PetscLibType,ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) 
Plots the solution at each iteration of an iterative solver.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_monitor_solution draw` - Activates `KSPMonitorSolutionDraw()`

Level: intermediate

-seealso: [](ch_ksp), `KSPMonitorSet()`, `KSPMonitorTrueResidual()`

# External Links
$(_doc_external("KSP/KSPMonitorSolutionDraw"))
"""
function KSPMonitorSolutionDraw(petsclib::PetscLibType, ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function KSPMonitorSolutionDraw(petsclib::$UnionPetscLib, ksp::PetscKSP, n::$PetscInt, rnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:KSPMonitorSolutionDraw, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	KSPMonitorSolutionDrawLG(petsclib::PetscLibType,ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) 
Plots the solution norm at each iteration of an iterative solver.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_monitor_solution draw::draw_lg` - Activates `KSPMonitorSolutionDrawLG()`

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPMonitorSolutionDrawLGCreate()`

# External Links
$(_doc_external("KSP/KSPMonitorSolutionDrawLG"))
"""
function KSPMonitorSolutionDrawLG(petsclib::PetscLibType, ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function KSPMonitorSolutionDrawLG(petsclib::$UnionPetscLib, ksp::PetscKSP, n::$PetscInt, rnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:KSPMonitorSolutionDrawLG, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	vf::PetscViewerAndFormat = KSPMonitorSolutionDrawLGCreate(petsclib::PetscLibType,viewer::PetscViewer, format::PetscViewerFormat, ctx::Cvoid) 
Creates the context for the `KSP` monitor `KSPMonitorSolutionDrawLG()`

Collective

Input Parameters:
- `viewer` - The `PetscViewer`
- `format` - The viewer format
- `ctx`    - An optional user context

Output Parameter:
- `vf` - The viewer context

Level: intermediate

-seealso: [](ch_ksp), `KSPMonitorSet()`, `KSPMonitorTrueResidual()`

# External Links
$(_doc_external("KSP/KSPMonitorSolutionDrawLGCreate"))
"""
function KSPMonitorSolutionDrawLGCreate(petsclib::PetscLibType, viewer::PetscViewer, format::PetscViewerFormat, ctx::Cvoid) end

@for_petsc function KSPMonitorSolutionDrawLGCreate(petsclib::$UnionPetscLib, viewer::PetscViewer, format::PetscViewerFormat, ctx::Cvoid )
	vf_ = Ref{PetscViewerAndFormat}()

    @chk ccall(
               (:KSPMonitorSolutionDrawLGCreate, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscViewerFormat, Ptr{Cvoid}, PetscViewerAndFormat),
               viewer, format, ctx, vf_,
              )

	vf = vf_[]

	return vf
end 

"""
	KSPMonitorSingularValue(petsclib::PetscLibType,ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) 
Prints the two norm of the true residual and estimation of the extreme singular values of the preconditioned problem at each iteration.

Logically Collective

Input Parameters:
- `ksp`   - the iterative context
- `n`     - the iteration
- `rnorm` - the two norm of the residual
- `vf`    - The viewer context

Options Database Key:
- `-ksp_monitor_singular_value` - Activates `KSPMonitorSingularValue()`

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPMonitorSet()`, `KSPComputeExtremeSingularValues()`, `KSPMonitorSingularValueCreate()`

# External Links
$(_doc_external("KSP/KSPMonitorSingularValue"))
"""
function KSPMonitorSingularValue(petsclib::PetscLibType, ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function KSPMonitorSingularValue(petsclib::$UnionPetscLib, ksp::PetscKSP, n::$PetscInt, rnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:KSPMonitorSingularValue, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	vf::PetscViewerAndFormat = KSPMonitorSingularValueCreate(petsclib::PetscLibType,viewer::PetscViewer, format::PetscViewerFormat, ctx::Cvoid) 
Creates the singular value monitor context needed by `KSPMonitorSingularValue()`

Collective

Input Parameters:
- `viewer` - The PetscViewer
- `format` - The viewer format
- `ctx`    - An optional user context

Output Parameter:
- `vf` - The viewer context

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPMonitorSet()`, `KSPMonitorSingularValue()`, `PetscViewer`

# External Links
$(_doc_external("KSP/KSPMonitorSingularValueCreate"))
"""
function KSPMonitorSingularValueCreate(petsclib::PetscLibType, viewer::PetscViewer, format::PetscViewerFormat, ctx::Cvoid) end

@for_petsc function KSPMonitorSingularValueCreate(petsclib::$UnionPetscLib, viewer::PetscViewer, format::PetscViewerFormat, ctx::Cvoid )
	vf_ = Ref{PetscViewerAndFormat}()

    @chk ccall(
               (:KSPMonitorSingularValueCreate, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscViewerFormat, Ptr{Cvoid}, PetscViewerAndFormat),
               viewer, format, ctx, vf_,
              )

	vf = vf_[]

	return vf
end 

"""
	ctx::Cvoid = KSPMonitorDynamicToleranceCreate(petsclib::PetscLibType) 
Creates the context used by `KSPMonitorDynamicTolerance()`

Logically Collective

Output Parameter:
- `ctx` - a void pointer

Options Database Key:
- `-sub_ksp_dynamic_tolerance <coef>` - coefficient of dynamic tolerance for inner solver, default is 1.0

Level: advanced

-seealso: [](sec_flexibleksp), `KSP`, `KSPMonitorDynamicTolerance()`, `KSPMonitorDynamicToleranceDestroy()`, `KSPMonitorDynamicToleranceSetCoefficient()`

# External Links
$(_doc_external("KSP/KSPMonitorDynamicToleranceCreate"))
"""
function KSPMonitorDynamicToleranceCreate(petsclib::PetscLibType) end

@for_petsc function KSPMonitorDynamicToleranceCreate(petsclib::$UnionPetscLib)
	ctx_ = Ref{Cvoid}()

    @chk ccall(
               (:KSPMonitorDynamicToleranceCreate, $petsc_library),
               PetscErrorCode,
               (Ptr{Cvoid},),
               ctx_,
              )

	ctx = ctx_[]

	return ctx
end 

"""
	coeff::PetscReal = KSPMonitorDynamicToleranceSetCoefficient(petsclib::PetscLibType,ctx::Cvoid) 
Sets the coefficient in the context used by `KSPMonitorDynamicTolerance()`

Logically Collective

Output Parameters:
- `ctx`   - the context for `KSPMonitorDynamicTolerance()`
- `coeff` - the coefficient, default is 1.0

Options Database Key:
- `-sub_ksp_dynamic_tolerance <coef>` - coefficient of dynamic tolerance for inner solver, default is 1.0

Level: advanced

-seealso: [](sec_flexibleksp), `KSP`, `KSPMonitorDynamicTolerance()`, `KSPMonitorDynamicToleranceDestroy()`, `KSPMonitorDynamicToleranceCreate()`

# External Links
$(_doc_external("KSP/KSPMonitorDynamicToleranceSetCoefficient"))
"""
function KSPMonitorDynamicToleranceSetCoefficient(petsclib::PetscLibType, ctx::Cvoid) end

@for_petsc function KSPMonitorDynamicToleranceSetCoefficient(petsclib::$UnionPetscLib, ctx::Cvoid )
	coeff_ = Ref{$PetscReal}()

    @chk ccall(
               (:KSPMonitorDynamicToleranceSetCoefficient, $petsc_library),
               PetscErrorCode,
               (Ptr{Cvoid}, $PetscReal),
               ctx, coeff_,
              )

	coeff = coeff_[]

	return coeff
end 

"""
	KSPMonitorDynamicTolerance(petsclib::PetscLibType,ksp::PetscKSP, its::PetscInt, fnorm::PetscReal, ctx::Cvoid) 
A monitor that changes the inner tolerance of nested preconditioners in every outer iteration in an adaptive way.

Collective

Input Parameters:
- `ksp`   - iterative context
- `its`   - iteration number (not used)
- `fnorm` - the current residual norm
- `ctx`   - context used by monitor

Options Database Key:
- `-sub_ksp_dynamic_tolerance <coef>` - coefficient of dynamic tolerance for inner solver, default is 1.0

Level: advanced

-seealso: [](sec_flexibleksp), `KSP`, `KSPMonitorDynamicToleranceCreate()`, `KSPMonitorDynamicToleranceDestroy()`, `KSPMonitorDynamicToleranceSetCoefficient()`

# External Links
$(_doc_external("KSP/KSPMonitorDynamicTolerance"))
"""
function KSPMonitorDynamicTolerance(petsclib::PetscLibType, ksp::PetscKSP, its::PetscInt, fnorm::PetscReal, ctx::Cvoid) end

@for_petsc function KSPMonitorDynamicTolerance(petsclib::$UnionPetscLib, ksp::PetscKSP, its::$PetscInt, fnorm::$PetscReal, ctx::Cvoid )

    @chk ccall(
               (:KSPMonitorDynamicTolerance, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{Cvoid}),
               ksp, its, fnorm, ctx,
              )


	return nothing
end 

"""
	KSPMonitorDynamicToleranceDestroy(petsclib::PetscLibType,ctx::Cvoid) 
Destroy the monitor context used in `KSPMonitorDynamicTolerance()`

Input Parameter:
- `ctx` - the monitor context

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPMonitorDynamicTolerance()`, `KSPMonitorSet()`, `KSPMonitorDynamicToleranceCreate()`

# External Links
$(_doc_external("KSP/KSPMonitorDynamicToleranceDestroy"))
"""
function KSPMonitorDynamicToleranceDestroy(petsclib::PetscLibType, ctx::Cvoid) end

@for_petsc function KSPMonitorDynamicToleranceDestroy(petsclib::$UnionPetscLib, ctx::Cvoid )

    @chk ccall(
               (:KSPMonitorDynamicToleranceDestroy, $petsc_library),
               PetscErrorCode,
               (Cvoid,),
               ctx,
              )


	return nothing
end 

"""
	KSPConvergedSkip(petsclib::PetscLibType,ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, reason::KSPConvergedReason, dtx::Cvoid) 
Convergence test that do not return as converged
until the maximum number of iterations is reached.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm residual value (may be estimated)
- `dtx`   - unused convergence context

Output Parameter:
- `reason` - `KSP_CONVERGED_ITERATING` or `KSP_CONVERGED_ITS`

Options Database Key:
- `-ksp_convergence_test skip` - skips the test

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPCG`, `KSPBCGS`, `KSPConvergenceTestFn`, `KSPSetConvergenceTest()`, `KSPSetTolerances()`, `KSPSetNormType()`, [](sec_flexibleksp),
`KSPConvergedReason`

# External Links
$(_doc_external("KSP/KSPConvergedSkip"))
"""
function KSPConvergedSkip(petsclib::PetscLibType, ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, reason::KSPConvergedReason, dtx::Cvoid) end

@for_petsc function KSPConvergedSkip(petsclib::$UnionPetscLib, ksp::PetscKSP, n::$PetscInt, rnorm::$PetscReal, reason::KSPConvergedReason, dtx::Cvoid )

    @chk ccall(
               (:KSPConvergedSkip, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{KSPConvergedReason}, Ptr{Cvoid}),
               ksp, n, rnorm, reason, dtx,
              )


	return nothing
end 

"""
	KSPSetConvergedNegativeCurvature(petsclib::PetscLibType,ksp::PetscKSP, flg::PetscBool) 
Allows to declare convergence and return `KSP_CONVERGED_NEG_CURVE` when negative curvature is detected

Collective

Input Parameters:
- `ksp` - iterative context
- `flg` - the Boolean value

Options Database Key:
- `-ksp_converged_neg_curve <bool>` - Declare convergence if negative curvature is detected

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPConvergedReason`, `KSPGetConvergedNegativeCurvature()`

# External Links
$(_doc_external("KSP/KSPSetConvergedNegativeCurvature"))
"""
function KSPSetConvergedNegativeCurvature(petsclib::PetscLibType, ksp::PetscKSP, flg::PetscBool) end

@for_petsc function KSPSetConvergedNegativeCurvature(petsclib::$UnionPetscLib, ksp::PetscKSP, flg::PetscBool )

    @chk ccall(
               (:KSPSetConvergedNegativeCurvature, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flg,
              )


	return nothing
end 

"""
	flg::PetscBool = KSPGetConvergedNegativeCurvature(petsclib::PetscLibType,ksp::PetscKSP) 
Get the flag to declare convergence if negative curvature is detected

Collective

Input Parameter:
- `ksp` - iterative context

Output Parameter:
- `flg` - the Boolean value

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPConvergedReason`, `KSPSetConvergedNegativeCurvature()`

# External Links
$(_doc_external("KSP/KSPGetConvergedNegativeCurvature"))
"""
function KSPGetConvergedNegativeCurvature(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGetConvergedNegativeCurvature(petsclib::$UnionPetscLib, ksp::PetscKSP )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:KSPGetConvergedNegativeCurvature, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{PetscBool}),
               ksp, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	ctx::Cvoid = KSPConvergedDefaultCreate(petsclib::PetscLibType) 
Creates and initializes the context used by the `KSPConvergedDefault()` function

Not Collective

Output Parameter:
- `ctx` - convergence context

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPConvergedDefault()`, `KSPConvergedDefaultDestroy()`, `KSPSetConvergenceTest()`, `KSPSetTolerances()`,
`KSPConvergedSkip()`, `KSPConvergedReason`, `KSPGetConvergedReason()`, `KSPConvergedDefaultSetUIRNorm()`, `KSPConvergedDefaultSetUMIRNorm()`,
`KSPConvergedDefaultSetConvergedMaxits()`

# External Links
$(_doc_external("KSP/KSPConvergedDefaultCreate"))
"""
function KSPConvergedDefaultCreate(petsclib::PetscLibType) end

@for_petsc function KSPConvergedDefaultCreate(petsclib::$UnionPetscLib)
	ctx_ = Ref{Cvoid}()

    @chk ccall(
               (:KSPConvergedDefaultCreate, $petsc_library),
               PetscErrorCode,
               (Cvoid,),
               ctx_,
              )

	ctx = ctx_[]

	return ctx
end 

"""
	KSPConvergedDefaultSetUIRNorm(petsclib::PetscLibType,ksp::PetscKSP) 
makes the default convergence test use  || B*(b
instead of  || B*b ||. In the case of right preconditioner or if `KSPSetNormType`(ksp,`KSP_NORM_UNPRECONDITIONED`)
is used there is no B in the above formula.

Collective

Input Parameters:
- `ksp` - iterative context

Options Database Key:
- `-ksp_converged_use_initial_residual_norm <bool>` - Use initial residual norm for computing relative convergence

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPSetConvergenceTest()`, `KSPSetTolerances()`, `KSPConvergedSkip()`, `KSPConvergedReason`, `KSPGetConvergedReason()`, `KSPConvergedDefaultSetUMIRNorm()`, `KSPConvergedDefaultSetConvergedMaxits()`

# External Links
$(_doc_external("KSP/KSPConvergedDefaultSetUIRNorm"))
"""
function KSPConvergedDefaultSetUIRNorm(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPConvergedDefaultSetUIRNorm(petsclib::$UnionPetscLib, ksp::PetscKSP )

    @chk ccall(
               (:KSPConvergedDefaultSetUIRNorm, $petsc_library),
               PetscErrorCode,
               (CKSP,),
               ksp,
              )


	return nothing
end 

"""
	KSPConvergedDefaultSetUMIRNorm(petsclib::PetscLibType,ksp::PetscKSP) 
makes the default convergence test use min(|| B*(b
In the case of right preconditioner or if `KSPSetNormType`(ksp,`KSP_NORM_UNPRECONDITIONED`)
is used there is no B in the above formula.

Collective

Input Parameters:
- `ksp` - iterative context

Options Database Key:
- `-ksp_converged_use_min_initial_residual_norm <bool>` - Use minimum of initial residual norm and b for computing relative convergence

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPSetConvergenceTest()`, `KSPSetTolerances()`, `KSPConvergedSkip()`, `KSPConvergedReason`, `KSPGetConvergedReason()`, `KSPConvergedDefaultSetUIRNorm()`, `KSPConvergedDefaultSetConvergedMaxits()`

# External Links
$(_doc_external("KSP/KSPConvergedDefaultSetUMIRNorm"))
"""
function KSPConvergedDefaultSetUMIRNorm(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPConvergedDefaultSetUMIRNorm(petsclib::$UnionPetscLib, ksp::PetscKSP )

    @chk ccall(
               (:KSPConvergedDefaultSetUMIRNorm, $petsc_library),
               PetscErrorCode,
               (CKSP,),
               ksp,
              )


	return nothing
end 

"""
	KSPConvergedDefaultSetConvergedMaxits(petsclib::PetscLibType,ksp::PetscKSP, flg::PetscBool) 
allows the default convergence test to declare convergence and return `KSP_CONVERGED_ITS` if the maximum number of iterations is reached

Collective

Input Parameters:
- `ksp` - iterative context
- `flg` - boolean flag

Options Database Key:
- `-ksp_converged_maxits <bool>` - Declare convergence if the maximum number of iterations is reached

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPSetConvergenceTest()`, `KSPSetTolerances()`, `KSPConvergedSkip()`, `KSPConvergedReason`, `KSPGetConvergedReason()`, `KSPConvergedDefaultSetUMIRNorm()`, `KSPConvergedDefaultSetUIRNorm()`

# External Links
$(_doc_external("KSP/KSPConvergedDefaultSetConvergedMaxits"))
"""
function KSPConvergedDefaultSetConvergedMaxits(petsclib::PetscLibType, ksp::PetscKSP, flg::PetscBool) end

@for_petsc function KSPConvergedDefaultSetConvergedMaxits(petsclib::$UnionPetscLib, ksp::PetscKSP, flg::PetscBool )

    @chk ccall(
               (:KSPConvergedDefaultSetConvergedMaxits, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flg,
              )


	return nothing
end 

"""
	KSPConvergedDefault(petsclib::PetscLibType,ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, reason::KSPConvergedReason, ctx::Cvoid) 
Default code to determine convergence of the linear iterative solvers

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - residual norm (may be estimated, depending on the method may be the preconditioned residual norm)
- `ctx`   - convergence context which must be created by `KSPConvergedDefaultCreate()`

Output Parameter:
- `reason` - the convergence reason; it is positive if the iteration has converged,
negative if the iteration has diverged, and `KSP_CONVERGED_ITERATING` otherwise

Options Database Keys:
- `-ksp_max_it`                                  - maximum number of linear iterations
- `-ksp_min_it`                                  - minimum number of linear iterations, defaults to 0
- `-ksp_rtol rtol`                               - relative tolerance used in default determination of convergence, i.e. if residual norm decreases by this factor than convergence is declared
- `-ksp_atol abstol`                             - absolute tolerance used in default convergence test, i.e. if residual norm is less than this then convergence is declared
- `-ksp_divtol tol`                              - if residual norm increases by this factor than divergence is declared
- `-ksp_converged_use_initial_residual_norm`     - see `KSPConvergedDefaultSetUIRNorm()`
- `-ksp_converged_use_min_initial_residual_norm` - see `KSPConvergedDefaultSetUMIRNorm()`
- `-ksp_converged_maxits`                        - see `KSPConvergedDefaultSetConvergedMaxits()`

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPSetConvergenceTest()`, `KSPSetTolerances()`, `KSPConvergedSkip()`, `KSPConvergedReason`, `KSPGetConvergedReason()`,
`KSPSetMinimumIterations()`, `KSPConvergenceTestFn`,
`KSPConvergedDefaultSetUIRNorm()`, `KSPConvergedDefaultSetUMIRNorm()`, `KSPConvergedDefaultSetConvergedMaxits()`, `KSPConvergedDefaultCreate()`, `KSPConvergedDefaultDestroy()`

# External Links
$(_doc_external("KSP/KSPConvergedDefault"))
"""
function KSPConvergedDefault(petsclib::PetscLibType, ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, reason::KSPConvergedReason, ctx::Cvoid) end

@for_petsc function KSPConvergedDefault(petsclib::$UnionPetscLib, ksp::PetscKSP, n::$PetscInt, rnorm::$PetscReal, reason::KSPConvergedReason, ctx::Cvoid )

    @chk ccall(
               (:KSPConvergedDefault, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{KSPConvergedReason}, Ptr{Cvoid}),
               ksp, n, rnorm, reason, ctx,
              )


	return nothing
end 

"""
	KSPConvergedDefaultDestroy(petsclib::PetscLibType,ctx::Cvoid) 
Frees the space used by the `KSPConvergedDefault()` function context

Not Collective

Input Parameter:
- `ctx` - convergence context

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPConvergedDefault()`, `KSPConvergedDefaultCreate()`, `KSPSetConvergenceTest()`, `KSPSetTolerances()`, `KSPConvergedSkip()`,
`KSPConvergedReason`, `KSPGetConvergedReason()`, `KSPConvergedDefaultSetUIRNorm()`, `KSPConvergedDefaultSetUMIRNorm()`

# External Links
$(_doc_external("KSP/KSPConvergedDefaultDestroy"))
"""
function KSPConvergedDefaultDestroy(petsclib::PetscLibType, ctx::Cvoid) end

@for_petsc function KSPConvergedDefaultDestroy(petsclib::$UnionPetscLib, ctx::Cvoid )

    @chk ccall(
               (:KSPConvergedDefaultDestroy, $petsc_library),
               PetscErrorCode,
               (Cvoid,),
               ctx,
              )


	return nothing
end 

"""
	KSPBuildSolutionDefault(petsclib::PetscLibType,ksp::PetscKSP, v::PetscVec, V::PetscVec) 

# External Links
$(_doc_external("KSP/KSPBuildSolutionDefault"))
"""
function KSPBuildSolutionDefault(petsclib::PetscLibType, ksp::PetscKSP, v::PetscVec, V::PetscVec) end

@for_petsc function KSPBuildSolutionDefault(petsclib::$UnionPetscLib, ksp::PetscKSP, v::PetscVec, V::PetscVec )
	V_ = Ref(V.ptr)

    @chk ccall(
               (:KSPBuildSolutionDefault, $petsc_library),
               PetscErrorCode,
               (CKSP, CVec, Ptr{CVec}),
               ksp, v, V_,
              )

	V.ptr = C_NULL

	return nothing
end 

"""
	KSPBuildResidualDefault(petsclib::PetscLibType,ksp::PetscKSP, t::PetscVec, v::PetscVec, V::PetscVec) 
Default code to compute the residual.

Collecive on ksp

Input Parameters:
- `ksp` - iterative context
- `t`   - pointer to temporary vector
- `v`   - pointer to user vector

Output Parameter:
- `V` - pointer to a vector containing the residual

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPBuildSolutionDefault()`

# External Links
$(_doc_external("KSP/KSPBuildResidualDefault"))
"""
function KSPBuildResidualDefault(petsclib::PetscLibType, ksp::PetscKSP, t::PetscVec, v::PetscVec, V::PetscVec) end

@for_petsc function KSPBuildResidualDefault(petsclib::$UnionPetscLib, ksp::PetscKSP, t::PetscVec, v::PetscVec, V::PetscVec )
	V_ = Ref(V.ptr)

    @chk ccall(
               (:KSPBuildResidualDefault, $petsc_library),
               PetscErrorCode,
               (CKSP, CVec, CVec, Ptr{CVec}),
               ksp, t, v, V_,
              )

	V.ptr = C_NULL

	return nothing
end 

"""
	right::Vector{PetscVec},left::Vector{PetscVec} = KSPCreateVecs(petsclib::PetscLibType,ksp::PetscKSP, rightn::PetscInt, leftn::PetscInt) 
Gets a number of work vectors suitably sized for the operator in the `KSP`

Collective

Input Parameters:
- `ksp`    - iterative context
- `rightn` - number of right work vectors to allocate
- `leftn`  - number of left work vectors to allocate

Output Parameters:
- `right` - the array of vectors created
- `left`  - the array of left vectors

Level: advanced

-seealso: [](ch_ksp), `MatCreateVecs()`, `VecDestroyVecs()`, `KSPSetWorkVecs()`

# External Links
$(_doc_external("KSP/KSPCreateVecs"))
"""
function KSPCreateVecs(petsclib::PetscLibType, ksp::PetscKSP, rightn::PetscInt, leftn::PetscInt) end

@for_petsc function KSPCreateVecs(petsclib::$UnionPetscLib, ksp::PetscKSP, rightn::$PetscInt, leftn::$PetscInt )
	right_ = Ref{Ptr{CVec}}()
	left_ = Ref{Ptr{CVec}}()

    @chk ccall(
               (:KSPCreateVecs, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, Ptr{Ptr{CVec}}, $PetscInt, Ptr{Ptr{CVec}}),
               ksp, rightn, right_, leftn, left_,
              )

	a_v = unsafe_wrap(Array, right_[], rightn; own = false)
    if rightn != 0
        v = PetscVec(a_v[1], petsclib)
        right = ntuple(i -> similar(v), rightn)
    else
        right = nothing
    end

    
    a_v = unsafe_wrap(Array, left_[], leftn; own = false)
    if leftn != 0
        v = PetscVec(a_v[1], petsclib)
        left = ntuple(i -> similar(v), leftn)
    else
        left = nothing
    end

	return right,left
end 

"""
	KSPSetWorkVecs(petsclib::PetscLibType,ksp::PetscKSP, nw::PetscInt) 
Sets a number of work vectors into a `KSP` object

Collective

Input Parameters:
- `ksp` - iterative context
- `nw`  - number of work vectors to allocate

Level: developer

-seealso: [](ch_ksp), `KSP`, `KSPCreateVecs()`

# External Links
$(_doc_external("KSP/KSPSetWorkVecs"))
"""
function KSPSetWorkVecs(petsclib::PetscLibType, ksp::PetscKSP, nw::PetscInt) end

@for_petsc function KSPSetWorkVecs(petsclib::$UnionPetscLib, ksp::PetscKSP, nw::$PetscInt )

    @chk ccall(
               (:KSPSetWorkVecs, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, nw,
              )


	return nothing
end 

"""
	KSPDestroyDefault(petsclib::PetscLibType,ksp::PetscKSP) 

# External Links
$(_doc_external("KSP/KSPDestroyDefault"))
"""
function KSPDestroyDefault(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPDestroyDefault(petsclib::$UnionPetscLib, ksp::PetscKSP )

    @chk ccall(
               (:KSPDestroyDefault, $petsc_library),
               PetscErrorCode,
               (CKSP,),
               ksp,
              )


	return nothing
end 

"""
	KSPGetConvergedReason(petsclib::PetscLibType,ksp::PetscKSP, reason::KSPConvergedReason) 
Gets the reason the `KSP` iteration was stopped.

Not Collective

Input Parameter:
- `ksp` - the `KSP` context

Output Parameter:
- `reason` - negative value indicates diverged, positive value converged, see `KSPConvergedReason` for the possible values

Options Database Key:
- `-ksp_converged_reason` - prints the reason to standard out when the solve ends

Level: intermediate

-seealso: [](ch_ksp), `KSPConvergedReason`, `KSP`, `KSPSetConvergenceTest()`, `KSPConvergedDefault()`, `KSPSetTolerances()`,
`KSPConvergedReasonView()`, `KSPGetConvergedReasonString()`

# External Links
$(_doc_external("KSP/KSPGetConvergedReason"))
"""
function KSPGetConvergedReason(petsclib::PetscLibType, ksp::PetscKSP, reason::KSPConvergedReason) end

@for_petsc function KSPGetConvergedReason(petsclib::$UnionPetscLib, ksp::PetscKSP, reason::KSPConvergedReason )

    @chk ccall(
               (:KSPGetConvergedReason, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPConvergedReason}),
               ksp, reason,
              )


	return nothing
end 

"""
	KSPGetConvergedReasonString(petsclib::PetscLibType,ksp::PetscKSP, strreason::String) 
Return a human readable string for a `KSPConvergedReason`

Not Collective

Input Parameter:
- `ksp` - the `KSP` context

Output Parameter:
- `strreason` - a human readable string that describes ksp converged reason

Level: beginner

-seealso: [](ch_ksp), `KSP`, `KSPGetConvergedReason()`

# External Links
$(_doc_external("KSP/KSPGetConvergedReasonString"))
"""
function KSPGetConvergedReasonString(petsclib::PetscLibType, ksp::PetscKSP, strreason::String) end

@for_petsc function KSPGetConvergedReasonString(petsclib::$UnionPetscLib, ksp::PetscKSP, strreason::String )
	strreason_ = Ref(pointer(strreason))

    @chk ccall(
               (:KSPGetConvergedReasonString, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Ptr{Cchar}}),
               ksp, strreason_,
              )


	return nothing
end 

"""
	KSPSetDM(petsclib::PetscLibType,ksp::PetscKSP, dm::PetscDM) 
Sets the `DM` that may be used by some preconditioners and that may be used to construct the linear system

Logically Collective

Input Parameters:
- `ksp` - the `KSP`
- `dm`  - the `DM`, cannot be `NULL` to remove a previously set `DM`

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `DM`, `KSPGetDM()`, `KSPSetDMActive()`, `KSPSetComputeOperators()`, `KSPSetComputeRHS()`, `KSPSetComputeInitialGuess()`, `DMKSPSetComputeOperators()`, `DMKSPSetComputeRHS()`, `DMKSPSetComputeInitialGuess()`

# External Links
$(_doc_external("KSP/KSPSetDM"))
"""
function KSPSetDM(petsclib::PetscLibType, ksp::PetscKSP, dm::PetscDM) end

@for_petsc function KSPSetDM(petsclib::$UnionPetscLib, ksp::PetscKSP, dm::PetscDM )

    @chk ccall(
               (:KSPSetDM, $petsc_library),
               PetscErrorCode,
               (CKSP, CDM),
               ksp, dm,
              )


	return nothing
end 

"""
	KSPSetDMActive(petsclib::PetscLibType,ksp::PetscKSP, flg::PetscBool) 
Indicates the `DM` should be used to generate the linear system matrix and right

Logically Collective

Input Parameters:
- `ksp` - the `KSP`
- `flg` - use the `DM`

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `DM`, `KSPGetDM()`, `KSPSetDM()`, `SNESSetDM()`, `KSPSetComputeOperators()`, `KSPSetComputeRHS()`, `KSPSetComputeInitialGuess()`

# External Links
$(_doc_external("KSP/KSPSetDMActive"))
"""
function KSPSetDMActive(petsclib::PetscLibType, ksp::PetscKSP, flg::PetscBool) end

@for_petsc function KSPSetDMActive(petsclib::$UnionPetscLib, ksp::PetscKSP, flg::PetscBool )

    @chk ccall(
               (:KSPSetDMActive, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flg,
              )


	return nothing
end 

"""
    dm::PetscDM = KSPGetDM(petsclib::PetscLibType,ksp::PetscKSP) 
Gets the `DM` that may be used by some preconditioners and that may be used to construct the linear system

Not Collective

Input Parameter:
- `ksp` - the `KSP`

Output Parameter:
- `dm` - the `DM`

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `DM`, `KSPSetDM()`, `KSPSetDMActive()`

# External Links
$(_doc_external("KSP/KSPGetDM"))
"""
function KSPGetDM(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGetDM(petsclib::$UnionPetscLib, ksp::PetscKSP)
	dm_ = Ref{CDM}()

    @chk ccall(
               (:KSPGetDM, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{CDM}),
               ksp, dm_,
              )

    dm = PetscDM(dm_[], petsclib)

	return dm
end 

"""
	KSPSetApplicationContext(petsclib::PetscLibType,ksp::PetscKSP, ctx::Cvoid) 
Sets the optional user

Logically Collective

Input Parameters:
- `ksp` - the `KSP` context
- `ctx` - user context

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPGetApplicationContext()`

# External Links
$(_doc_external("KSP/KSPSetApplicationContext"))
"""
function KSPSetApplicationContext(petsclib::PetscLibType, ksp::PetscKSP, ctx::Cvoid) end

@for_petsc function KSPSetApplicationContext(petsclib::$UnionPetscLib, ksp::PetscKSP, ctx::Cvoid )

    @chk ccall(
               (:KSPSetApplicationContext, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Cvoid}),
               ksp, ctx,
              )


	return nothing
end 

"""
	KSPGetApplicationContext(petsclib::PetscLibType,ksp::PetscKSP, ctx::PeCtx) 
Gets the user

Not Collective

Input Parameter:
- `ksp` - `KSP` context

Output Parameter:
- `ctx` - a pointer to the user context

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPSetApplicationContext()`

# External Links
$(_doc_external("KSP/KSPGetApplicationContext"))
"""
function KSPGetApplicationContext(petsclib::PetscLibType, ksp::PetscKSP, ctx::PeCtx) end

@for_petsc function KSPGetApplicationContext(petsclib::$UnionPetscLib, ksp::PetscKSP, ctx::PeCtx )

    @chk ccall(
               (:KSPGetApplicationContext, $petsc_library),
               PetscErrorCode,
               (CKSP, PeCtx),
               ksp, ctx,
              )


	return nothing
end 

"""
	KSPCheckSolve(petsclib::PetscLibType,ksp::PetscKSP, pc::PC, vec::PetscVec) 
Checks if the `PCSetUp()` or `KSPSolve()` failed and set the error flag for the outer `PC`. A `KSP_DIVERGED_ITS` is
not considered a failure in this context

Collective

Input Parameters:
- `ksp` - the linear solver `KSP` context.
- `pc`  - the preconditioner context
- `vec` - a vector that will be initialized with Inf to indicate lack of convergence

Level: developer

-seealso: [](ch_ksp), `KSP`, `KSPCreate()`, `KSPSetType()`, `KSPCheckNorm()`, `KSPCheckDot()`

# External Links
$(_doc_external("KSP/KSPCheckSolve"))
"""
function KSPCheckSolve(petsclib::PetscLibType, ksp::PetscKSP, pc::PC, vec::PetscVec) end

@for_petsc function KSPCheckSolve(petsclib::$UnionPetscLib, ksp::PetscKSP, pc::PC, vec::PetscVec )

    @chk ccall(
               (:KSPCheckSolve, $petsc_library),
               PetscErrorCode,
               (CKSP, PC, CVec),
               ksp, pc, vec,
              )


	return nothing
end 

"""
	KSPInitialResidual(petsclib::PetscLibType,ksp::PetscKSP, vsoln::PetscVec, vt1::PetscVec, vt2::PetscVec, vres::PetscVec, vb::PetscVec) 
Computes the residual. Either b
preconditioning or C*(b - A*x) with left preconditioning; the latter
residual is often called the "preconditioned residual".

Collective

Input Parameters:
- `ksp`   - the `KSP` solver object
- `vsoln` - solution to use in computing residual
- `vt1`   - temporary work vector
- `vt2`   - temporary work vector
- `vb`    - right-hand-side vector

Output Parameter:
- `vres` - calculated residual

Level: developer

-seealso: [](ch_ksp), `KSP`, `KSPSolve()`, `KSPMonitor()`

# External Links
$(_doc_external("KSP/KSPInitialResidual"))
"""
function KSPInitialResidual(petsclib::PetscLibType, ksp::PetscKSP, vsoln::PetscVec, vt1::PetscVec, vt2::PetscVec, vres::PetscVec, vb::PetscVec) end

@for_petsc function KSPInitialResidual(petsclib::$UnionPetscLib, ksp::PetscKSP, vsoln::PetscVec, vt1::PetscVec, vt2::PetscVec, vres::PetscVec, vb::PetscVec )

    @chk ccall(
               (:KSPInitialResidual, $petsc_library),
               PetscErrorCode,
               (CKSP, CVec, CVec, CVec, CVec, CVec),
               ksp, vsoln, vt1, vt2, vres, vb,
              )


	return nothing
end 

"""
	KSPUnwindPreconditioner(petsclib::PetscLibType,ksp::PetscKSP, vsoln::PetscVec, vt1::PetscVec) 
Unwinds the preconditioning in the solution. That is,
takes solution to the preconditioned problem and gets the solution to the
original problem from it.

Collective

Input Parameters:
- `ksp`   - iterative context
- `vsoln` - solution vector
- `vt1`   - temporary work vector

Output Parameter:
- `vsoln` - contains solution on output

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPSetPCSide()`

# External Links
$(_doc_external("KSP/KSPUnwindPreconditioner"))
"""
function KSPUnwindPreconditioner(petsclib::PetscLibType, ksp::PetscKSP, vsoln::PetscVec, vt1::PetscVec) end

@for_petsc function KSPUnwindPreconditioner(petsclib::$UnionPetscLib, ksp::PetscKSP, vsoln::PetscVec, vt1::PetscVec )

    @chk ccall(
               (:KSPUnwindPreconditioner, $petsc_library),
               PetscErrorCode,
               (CKSP, CVec, CVec),
               ksp, vsoln, vt1,
              )


	return nothing
end 

"""
	KSPLoad(petsclib::PetscLibType,newdm::PetscKSP, viewer::PetscViewer) 
Loads a `KSP` that has been stored in a `PETSCVIEWERBINARY`  with `KSPView()`.

Collective

Input Parameters:
- `newdm`  - the newly loaded `KSP`, this needs to have been created with `KSPCreate()` or
some related function before a call to `KSPLoad()`.
- `viewer` - binary file viewer, obtained from `PetscViewerBinaryOpen()`

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `PetscViewerBinaryOpen()`, `KSPView()`, `MatLoad()`, `VecLoad()`

# External Links
$(_doc_external("KSP/KSPLoad"))
"""
function KSPLoad(petsclib::PetscLibType, newdm::PetscKSP, viewer::PetscViewer) end

@for_petsc function KSPLoad(petsclib::$UnionPetscLib, newdm::PetscKSP, viewer::PetscViewer )

    @chk ccall(
               (:KSPLoad, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscViewer),
               newdm, viewer,
              )


	return nothing
end 

"""
	KSPView(petsclib::PetscLibType,ksp::PetscKSP, viewer::PetscViewer) 
Prints the various parameters currently set in the `KSP` object. For example, the convergence tolerances and `KSPType`.
Also views the `PC` and `Mat` contained by the `KSP` with `PCView()` and `MatView()`.

Collective

Input Parameters:
- `ksp`    - the Krylov space context
- `viewer` - visualization context

Options Database Key:
- `-ksp_view` - print the `KSP` data structure at the end of each `KSPSolve()` call

Level: beginner

-seealso: [](ch_ksp), `KSP`, `PetscViewer`, `PCView()`, `PetscViewerASCIIOpen()`, `KSPViewFromOptions()`

# External Links
$(_doc_external("KSP/KSPView"))
"""
function KSPView(petsclib::PetscLibType, ksp::PetscKSP, viewer::PetscViewer) end

@for_petsc function KSPView(petsclib::$UnionPetscLib, ksp::PetscKSP, viewer::PetscViewer )

    @chk ccall(
               (:KSPView, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscViewer),
               ksp, viewer,
              )


	return nothing
end 

"""
	KSPViewFromOptions(petsclib::PetscLibType,A::PetscKSP, obj::PetscObject, name::String) 
View (print) a `KSP` object based on values in the options database. Also views the `PC` and `Mat` contained by the `KSP`
with `PCView()` and `MatView()`.

Collective

Input Parameters:
- `A`    - Krylov solver context
- `obj`  - Optional object that provides the options prefix used to query the options database
- `name` - command line option

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPView()`, `PetscObjectViewFromOptions()`, `KSPCreate()`

# External Links
$(_doc_external("KSP/KSPViewFromOptions"))
"""
function KSPViewFromOptions(petsclib::PetscLibType, A::PetscKSP, obj::PetscObject, name::String) end

@for_petsc function KSPViewFromOptions(petsclib::$UnionPetscLib, A::PetscKSP, obj::PetscObject, name::String )

    @chk ccall(
               (:KSPViewFromOptions, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

"""
	KSPSetNormType(petsclib::PetscLibType,ksp::PetscKSP, normtype::KSPNormType) 
Sets the type of residual norm that is used for convergence testing in `KSPSolve()` for the given `KSP` context

Logically Collective

Input Parameters:
- `ksp`      - Krylov solver context
- `normtype` - one of
-seealso: [](ch_ksp), `KSPSetUp()`, `KSPSolve()`, `KSPDestroy()`, `KSPConvergedSkip()`, `KSPSetCheckNormIteration()`, `KSPSetPCSide()`, `KSPGetPCSide()`, `KSPNormType`

# External Links
$(_doc_external("KSP/KSPSetNormType"))
"""
function KSPSetNormType(petsclib::PetscLibType, ksp::PetscKSP, normtype::KSPNormType) end

@for_petsc function KSPSetNormType(petsclib::$UnionPetscLib, ksp::PetscKSP, normtype::KSPNormType )

    @chk ccall(
               (:KSPSetNormType, $petsc_library),
               PetscErrorCode,
               (CKSP, KSPNormType),
               ksp, normtype,
              )


	return nothing
end 

"""
	KSPSetCheckNormIteration(petsclib::PetscLibType,ksp::PetscKSP, it::PetscInt) 
Sets the first iteration at which the norm of the residual will be
computed and used in the convergence test of `KSPSolve()` for the given `KSP` context

Logically Collective

Input Parameters:
- `ksp` - Krylov solver context
- `it`  - use -1 to check at all iterations

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPSetUp()`, `KSPSolve()`, `KSPDestroy()`, `KSPConvergedSkip()`, `KSPSetNormType()`, `KSPSetLagNorm()`

# External Links
$(_doc_external("KSP/KSPSetCheckNormIteration"))
"""
function KSPSetCheckNormIteration(petsclib::PetscLibType, ksp::PetscKSP, it::PetscInt) end

@for_petsc function KSPSetCheckNormIteration(petsclib::$UnionPetscLib, ksp::PetscKSP, it::$PetscInt )

    @chk ccall(
               (:KSPSetCheckNormIteration, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, it,
              )


	return nothing
end 

"""
	KSPSetLagNorm(petsclib::PetscLibType,ksp::PetscKSP, flg::PetscBool) 
Lags the residual norm calculation so that it is computed as part of the `MPI_Allreduce()` used for
computing the inner products needed for the next iteration.

Logically Collective

Input Parameters:
- `ksp` - Krylov solver context
- `flg` - `PETSC_TRUE` or `PETSC_FALSE`

Options Database Key:
- `-ksp_lag_norm` - lag the calculated residual norm

Level: advanced

-seealso: [](ch_ksp), `KSPSetUp()`, `KSPSolve()`, `KSPDestroy()`, `KSPConvergedSkip()`, `KSPSetNormType()`, `KSPSetCheckNormIteration()`

# External Links
$(_doc_external("KSP/KSPSetLagNorm"))
"""
function KSPSetLagNorm(petsclib::PetscLibType, ksp::PetscKSP, flg::PetscBool) end

@for_petsc function KSPSetLagNorm(petsclib::$UnionPetscLib, ksp::PetscKSP, flg::PetscBool )

    @chk ccall(
               (:KSPSetLagNorm, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flg,
              )


	return nothing
end 

"""
	KSPSetSupportedNorm(petsclib::PetscLibType,ksp::PetscKSP, normtype::KSPNormType, pcside::PCSide, priority::PetscInt) 
Sets a norm and preconditioner side supported by a `KSPType`

Logically Collective

Input Parameters:
- `ksp`      - Krylov method
- `normtype` - supported norm type of the type `KSPNormType`
- `pcside`   - preconditioner side, of the type `PCSide` that can be used with this `KSPNormType`
- `priority` - positive integer preference for this combination; larger values have higher priority

Level: developer

-seealso: [](ch_ksp), `KSP`, `KSPNormType`, `PCSide`, `KSPSetNormType()`, `KSPSetPCSide()`

# External Links
$(_doc_external("KSP/KSPSetSupportedNorm"))
"""
function KSPSetSupportedNorm(petsclib::PetscLibType, ksp::PetscKSP, normtype::KSPNormType, pcside::PCSide, priority::PetscInt) end

@for_petsc function KSPSetSupportedNorm(petsclib::$UnionPetscLib, ksp::PetscKSP, normtype::KSPNormType, pcside::PCSide, priority::$PetscInt )

    @chk ccall(
               (:KSPSetSupportedNorm, $petsc_library),
               PetscErrorCode,
               (CKSP, KSPNormType, PCSide, $PetscInt),
               ksp, normtype, pcside, priority,
              )


	return nothing
end 

"""
	normtype::KSPNormType = KSPGetNormType(petsclib::PetscLibType,ksp::PetscKSP) 
Gets the `KSPNormType` that is used for convergence testing during `KSPSolve()` for this `KSP` context

Not Collective

Input Parameter:
- `ksp` - Krylov solver context

Output Parameter:
- `normtype` - the `KSPNormType` that is used for convergence testing

Level: advanced

-seealso: [](ch_ksp), `KSPNormType`, `KSPSetNormType()`, `KSPConvergedSkip()`

# External Links
$(_doc_external("KSP/KSPGetNormType"))
"""
function KSPGetNormType(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGetNormType(petsclib::$UnionPetscLib, ksp::PetscKSP )
	normtype_ = Ref{KSPNormType}()

    @chk ccall(
               (:KSPGetNormType, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPNormType}),
               ksp, normtype_,
              )

	normtype = unsafe_string(normtype_[])

	return normtype
end 

"""
	KSPSetOperators(petsclib::PetscLibType,ksp::PetscKSP, Amat::PetscMat, Pmat::PetscMat) 
Sets the matrix associated with the linear system
and a (possibly) different one from which the preconditioner will be built into the `KSP` context. The matrix will then be used during `KSPSolve()`

Collective

Input Parameters:
- `ksp`  - the `KSP` context
- `Amat` - the matrix that defines the linear system
- `Pmat` - the matrix to be used in constructing the preconditioner, usually the same as `Amat`.

Level: beginner

-seealso: [](ch_ksp), `KSP`, `Mat`, `KSPSolve()`, `KSPGetPC()`, `PCGetOperators()`, `PCSetOperators()`, `KSPGetOperators()`, `KSPSetComputeOperators()`, `KSPSetComputeInitialGuess()`, `KSPSetComputeRHS()`

# External Links
$(_doc_external("KSP/KSPSetOperators"))
"""
function KSPSetOperators(petsclib::PetscLibType, ksp::AbstractPetscKSP, Amat::AbstractPetscMat, Pmat::AbstractPetscMat) end

@for_petsc function KSPSetOperators(petsclib::$UnionPetscLib, ksp::AbstractPetscKSP, Amat::AbstractPetscMat, Pmat::AbstractPetscMat )

    @chk ccall(
               (:KSPSetOperators, $petsc_library),
               PetscErrorCode,
               (CKSP, CMat, CMat),
               ksp, Amat, Pmat,
              )


	return nothing
end 

"""
	Amat, PMat = KSPGetOperators(petsclib::PetscLibType,ksp::PetscKSP) 
Gets the matrix associated with the linear system
and a (possibly) different one used to construct the preconditioner from the `KSP` context

Collective

Input Parameter:
- `ksp` - the `KSP` context

Output Parameters:
- `Amat` - the matrix that defines the linear system
- `Pmat` - the matrix to be used in constructing the preconditioner, usually the same as `Amat`.

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPSolve()`, `KSPGetPC()`, `PCSetOperators()`, `KSPSetOperators()`, `KSPGetOperatorsSet()`

# External Links
$(_doc_external("KSP/KSPGetOperators"))
"""
function KSPGetOperators(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGetOperators(petsclib::$UnionPetscLib, ksp::PetscKSP)
	Amat_ = Ref{CVec}()
	Pmat_ = Ref{CVec}()

    @chk ccall(
               (:KSPGetOperators, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{CMat}, Ptr{CMat}),
               ksp, Amat_, Pmat_,
              )

    Amat =  PetscMat(Amat_[], petsclib)
    Pmat =  PetscMat(Pmat_[], petsclib)
    
	return Amat, Pmat
end 

"""
	mat::PetscBool,pmat::PetscBool = KSPGetOperatorsSet(petsclib::PetscLibType,ksp::PetscKSP) 
Determines if the matrix associated with the linear system and
possibly a different one from which the preconditioner will be built have been set in the `KSP` with `KSPSetOperators()`

Not Collective, though the results on all processes will be the same

Input Parameter:
- `ksp` - the `KSP` context

Output Parameters:
- `mat`  - the matrix associated with the linear system was set
- `pmat` - matrix from which the preconditioner will be built, usually the same as `mat` was set

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `PCSetOperators()`, `KSPGetOperators()`, `KSPSetOperators()`, `PCGetOperators()`, `PCGetOperatorsSet()`

# External Links
$(_doc_external("KSP/KSPGetOperatorsSet"))
"""
function KSPGetOperatorsSet(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGetOperatorsSet(petsclib::$UnionPetscLib, ksp::PetscKSP )
	mat_ = Ref{PetscBool}()
	pmat_ = Ref{PetscBool}()

    @chk ccall(
               (:KSPGetOperatorsSet, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{PetscBool}, Ptr{PetscBool}),
               ksp, mat_, pmat_,
              )

	mat = mat_[]
	pmat = pmat_[]

	return mat,pmat
end 

"""
	KSPSetPreSolve(petsclib::PetscLibType,ksp::PetscKSP, presolve::KSPPSolveFn, ctx::Cvoid) 
Sets a function that is called at the beginning of each `KSPSolve()`. Used in conjunction with `KSPSetPostSolve()`.

Logically Collective

Input Parameters:
- `ksp`      - the solver object
- `presolve` - the function to call before the solve, see` KSPPSolveFn`
- `ctx`      - an optional context needed by the function

Level: developer

-seealso: [](ch_ksp), `KSPPSolveFn`, `KSPSetUp()`, `KSPSolve()`, `KSPDestroy()`, `KSP`, `KSPSetPostSolve()`, `PCEISENSTAT`, `PCPreSolve()`, `PCPostSolve()`

# External Links
$(_doc_external("KSP/KSPSetPreSolve"))
"""
function KSPSetPreSolve(petsclib::PetscLibType, ksp::PetscKSP, presolve::KSPPSolveFn, ctx::Cvoid) end

@for_petsc function KSPSetPreSolve(petsclib::$UnionPetscLib, ksp::PetscKSP, presolve::KSPPSolveFn, ctx::Cvoid )

    @chk ccall(
               (:KSPSetPreSolve, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPPSolveFn}, Ptr{Cvoid}),
               ksp, presolve, ctx,
              )


	return nothing
end 

"""
	KSPSetPostSolve(petsclib::PetscLibType,ksp::PetscKSP, postsolve::KSPPSolveFn, ctx::Cvoid) 
Sets a function that is called at the end of each `KSPSolve()` (whether it converges or not). Used in conjunction with `KSPSetPreSolve()`.

Logically Collective

Input Parameters:
- `ksp`       - the solver object
- `postsolve` - the function to call after the solve, see` KSPPSolveFn`
- `ctx`       - an optional context needed by the function

Level: developer

-seealso: [](ch_ksp), `KSPPSolveFn`, `KSPSetUp()`, `KSPSolve()`, `KSPDestroy()`, `KSP`, `KSPSetPreSolve()`, `PCEISENSTAT`

# External Links
$(_doc_external("KSP/KSPSetPostSolve"))
"""
function KSPSetPostSolve(petsclib::PetscLibType, ksp::PetscKSP, postsolve::KSPPSolveFn, ctx::Cvoid) end

@for_petsc function KSPSetPostSolve(petsclib::$UnionPetscLib, ksp::PetscKSP, postsolve::KSPPSolveFn, ctx::Cvoid )

    @chk ccall(
               (:KSPSetPostSolve, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPPSolveFn}, Ptr{Cvoid}),
               ksp, postsolve, ctx,
              )


	return nothing
end 

"""
	KSPSetNestLevel(petsclib::PetscLibType,ksp::PetscKSP, level::PetscInt) 
sets the amount of nesting the `KSP` has. That is the number of levels of `KSP` above this `KSP` in a linear solve.

Collective

Input Parameters:
- `ksp`   - the `KSP`
- `level` - the nest level

Level: developer

-seealso: [](ch_ksp), `KSPSetUp()`, `KSPSolve()`, `KSPDestroy()`, `KSP`, `KSPGMRES`, `KSPType`, `KSPGetNestLevel()`, `PCSetKSPNestLevel()`, `PCGetKSPNestLevel()`

# External Links
$(_doc_external("KSP/KSPSetNestLevel"))
"""
function KSPSetNestLevel(petsclib::PetscLibType, ksp::PetscKSP, level::PetscInt) end

@for_petsc function KSPSetNestLevel(petsclib::$UnionPetscLib, ksp::PetscKSP, level::$PetscInt )

    @chk ccall(
               (:KSPSetNestLevel, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, level,
              )


	return nothing
end 

"""
	level::PetscInt = KSPGetNestLevel(petsclib::PetscLibType,ksp::PetscKSP) 
gets the amount of nesting the `KSP` has

Not Collective

Input Parameter:
- `ksp` - the `KSP`

Output Parameter:
- `level` - the nest level

Level: developer

-seealso: [](ch_ksp), `KSPSetUp()`, `KSPSolve()`, `KSPDestroy()`, `KSP`, `KSPGMRES`, `KSPType`, `KSPSetNestLevel()`, `PCSetKSPNestLevel()`, `PCGetKSPNestLevel()`

# External Links
$(_doc_external("KSP/KSPGetNestLevel"))
"""
function KSPGetNestLevel(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGetNestLevel(petsclib::$UnionPetscLib, ksp::PetscKSP )
	level_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPGetNestLevel, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscInt}),
               ksp, level_,
              )

	level = level_[]

	return level
end 

"""
	inksp::PetscKSP = KSPCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Creates the `KSP` context. This `KSP` context is used in PETSc to solve linear systems with `KSPSolve()`

Collective

Input Parameter:
- `comm` - MPI communicator

Output Parameter:
- `inksp` - location to put the `KSP` context

Level: beginner

-seealso: [](ch_ksp), `KSPSetUp()`, `KSPSolve()`, `KSPDestroy()`, `KSP`, `KSPGMRES`, `KSPType`, `KSPSetType()`

# External Links
$(_doc_external("KSP/KSPCreate"))
"""
function KSPCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function KSPCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	inksp_ = Ref{CKSP}()

    @chk ccall(
               (:KSPCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{CKSP}),
               comm, inksp_,
              )

	inksp = PetscKSP(inksp_[], petsclib)

	return inksp
end 

"""
	KSPSetType(petsclib::PetscLibType,ksp::PetscKSP, type::KSPType) 
Sets the algorithm/method to be used to solve the linear system with the given `KSP`

Logically Collective

Input Parameters:
- `ksp`  - the Krylov space context
- `type` - a known method

Options Database Key:
- `-ksp_type  <method>` - Sets the method; see `KSPGType` or use `-help` for a list  of available methods (for instance, cg or gmres)

Level: intermediate

-seealso: [](ch_ksp), `PCSetType()`, `KSPType`, `KSPRegister()`, `KSPCreate()`, `KSP`

# External Links
$(_doc_external("KSP/KSPSetType"))
"""
function KSPSetType(petsclib::PetscLibType, ksp::PetscKSP, type::KSPType) end

@for_petsc function KSPSetType(petsclib::$UnionPetscLib, ksp::PetscKSP, type::KSPType )

    @chk ccall(
               (:KSPSetType, $petsc_library),
               PetscErrorCode,
               (CKSP, KSPType),
               ksp, type,
              )


	return nothing
end 

"""
	type::KSPType = KSPGetType(petsclib::PetscLibType,ksp::PetscKSP) 
Gets the `KSP` type as a string from the `KSP` object.

Not Collective

Input Parameter:
- `ksp` - Krylov context

Output Parameter:
- `type` - name of the `KSP` method

Level: intermediate

-seealso: [](ch_ksp), `KSPType`, `KSP`, `KSPSetType()`

# External Links
$(_doc_external("KSP/KSPGetType"))
"""
function KSPGetType(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGetType(petsclib::$UnionPetscLib, ksp::PetscKSP )
	type_ = Ref{KSPType}()

    @chk ccall(
               (:KSPGetType, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPType}),
               ksp, type_,
              )
              
    if type_[] == C_NULL
        type = ""
    else
        type = unsafe_string(type_[])
    end

	return type
end 

"""
	KSPRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Adds a method, `KSPType`, to the Krylov subspace solver package.

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - name of a new user-defined solver
- `function` - routine to create method

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPType`, `KSPSetType`, `KSPRegisterAll()`

# External Links
$(_doc_external("KSP/KSPRegister"))
"""
function KSPRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function KSPRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:KSPRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	KSPMonitorRegister(petsclib::PetscLibType,name::String, vtype::PetscViewerType, format::PetscViewerFormat, monitor::KSPMonitorRegisterFn, create::KSPMonitorRegisterCreateFn, destroy::KSPMonitorRegisterDestroyFn) 
Registers a Krylov subspace solver monitor routine that may be accessed with `KSPMonitorSetFromOptions()`

Not Collective

Input Parameters:
- `name`    - name of a new monitor type
- `vtype`   - A `PetscViewerType` for the output
- `format`  - A `PetscViewerFormat` for the output
- `monitor` - Monitor routine, see `KSPMonitorRegisterFn`
- `create`  - Creation routine, or `NULL`
- `destroy` - Destruction routine, or `NULL`

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPMonitorSet()`, `KSPMonitorRegisterAll()`, `KSPMonitorSetFromOptions()`

# External Links
$(_doc_external("KSP/KSPMonitorRegister"))
"""
function KSPMonitorRegister(petsclib::PetscLibType, name::String, vtype::PetscViewerType, format::PetscViewerFormat, monitor::KSPMonitorRegisterFn, create::KSPMonitorRegisterCreateFn, destroy::KSPMonitorRegisterDestroyFn) end

@for_petsc function KSPMonitorRegister(petsclib::$UnionPetscLib, name::String, vtype::PetscViewerType, format::PetscViewerFormat, monitor::KSPMonitorRegisterFn, create::KSPMonitorRegisterCreateFn, destroy::KSPMonitorRegisterDestroyFn )

    @chk ccall(
               (:KSPMonitorRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, PetscViewerType, PetscViewerFormat, Ptr{KSPMonitorRegisterFn}, Ptr{KSPMonitorRegisterCreateFn}, Ptr{KSPMonitorRegisterDestroyFn}),
               name, vtype, format, monitor, create, destroy,
              )


	return nothing
end 

"""
	KSPSetOptionsPrefix(petsclib::PetscLibType,ksp::PetscKSP, prefix::String) 
Sets the prefix used for searching for all
`KSP` options in the database.

Logically Collective

Input Parameters:
- `ksp`    - the Krylov context
- `prefix` - the prefix string to prepend to all `KSP` option requests

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPAppendOptionsPrefix()`, `KSPGetOptionsPrefix()`, `KSPSetFromOptions()`

# External Links
$(_doc_external("KSP/KSPSetOptionsPrefix"))
"""
function KSPSetOptionsPrefix(petsclib::PetscLibType, ksp::PetscKSP, prefix::String) end

@for_petsc function KSPSetOptionsPrefix(petsclib::$UnionPetscLib, ksp::PetscKSP, prefix::String )

    @chk ccall(
               (:KSPSetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Cchar}),
               ksp, prefix,
              )


	return nothing
end 

"""
	KSPAppendOptionsPrefix(petsclib::PetscLibType,ksp::PetscKSP, prefix::String) 
Appends to the prefix used for searching for all
`KSP` options in the database.

Logically Collective

Input Parameters:
- `ksp`    - the Krylov context
- `prefix` - the prefix string to prepend to all `KSP` option requests

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPSetOptionsPrefix()`, `KSPGetOptionsPrefix()`, `KSPSetFromOptions()`

# External Links
$(_doc_external("KSP/KSPAppendOptionsPrefix"))
"""
function KSPAppendOptionsPrefix(petsclib::PetscLibType, ksp::PetscKSP, prefix::String) end

@for_petsc function KSPAppendOptionsPrefix(petsclib::$UnionPetscLib, ksp::PetscKSP, prefix::String )

    @chk ccall(
               (:KSPAppendOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Cchar}),
               ksp, prefix,
              )


	return nothing
end 

"""
	KSPSetUseFischerGuess(petsclib::PetscLibType,ksp::PetscKSP, model::PetscInt, size::PetscInt) 
Use the Paul Fischer algorithm or its variants to compute initial guesses for a set of solves with related right

Logically Collective

Input Parameters:
- `ksp`   - the Krylov context
- `model` - use model 1, model 2, model 3, or any other number to turn it off
- `size`  - size of subspace used to generate initial guess

Options Database Key:
- `-ksp_fischer_guess <model,size>` - uses the Fischer initial guess generator for repeated linear solves

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPSetOptionsPrefix()`, `KSPAppendOptionsPrefix()`, `KSPSetGuess()`, `KSPGetGuess()`, `KSPGuess`

# External Links
$(_doc_external("KSP/KSPSetUseFischerGuess"))
"""
function KSPSetUseFischerGuess(petsclib::PetscLibType, ksp::PetscKSP, model::PetscInt, size::PetscInt) end

@for_petsc function KSPSetUseFischerGuess(petsclib::$UnionPetscLib, ksp::PetscKSP, model::$PetscInt, size::$PetscInt )

    @chk ccall(
               (:KSPSetUseFischerGuess, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscInt),
               ksp, model, size,
              )


	return nothing
end 

"""
	KSPSetGuess(petsclib::PetscLibType,ksp::PetscKSP, guess::KSPGuess) 
Set the initial guess object `KSPGuess` to be used by the `KSP` object to generate initial guesses

Logically Collective

Input Parameters:
- `ksp`   - the Krylov context
- `guess` - the object created with `KSPGuessCreate()`

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPGuess`, `KSPSetOptionsPrefix()`, `KSPAppendOptionsPrefix()`, `KSPSetUseFischerGuess()`, `KSPGetGuess()`

# External Links
$(_doc_external("KSP/KSPSetGuess"))
"""
function KSPSetGuess(petsclib::PetscLibType, ksp::PetscKSP, guess::KSPGuess) end

@for_petsc function KSPSetGuess(petsclib::$UnionPetscLib, ksp::PetscKSP, guess::KSPGuess )

    @chk ccall(
               (:KSPSetGuess, $petsc_library),
               PetscErrorCode,
               (CKSP, KSPGuess),
               ksp, guess,
              )


	return nothing
end 

"""
	KSPGetGuess(petsclib::PetscLibType,ksp::PetscKSP, guess::KSPGuess) 
Gets the initial guess generator for the `KSP`.

Not Collective

Input Parameter:
- `ksp` - the Krylov context

Output Parameter:
- `guess` - the object

Level: developer

-seealso: [](ch_ksp), `KSPGuess`, `KSP`, `KSPSetOptionsPrefix()`, `KSPAppendOptionsPrefix()`, `KSPSetUseFischerGuess()`, `KSPSetGuess()`

# External Links
$(_doc_external("KSP/KSPGetGuess"))
"""
function KSPGetGuess(petsclib::PetscLibType, ksp::PetscKSP, guess::KSPGuess) end

@for_petsc function KSPGetGuess(petsclib::$UnionPetscLib, ksp::PetscKSP, guess::KSPGuess )

    @chk ccall(
               (:KSPGetGuess, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPGuess}),
               ksp, guess,
              )


	return nothing
end 

"""
	prefix::String = KSPGetOptionsPrefix(petsclib::PetscLibType,ksp::PetscKSP) 
Gets the prefix used for searching for all `KSP` options in the database.

Not Collective

Input Parameter:
- `ksp` - the Krylov context

Output Parameter:
- `prefix` - pointer to the prefix string used is returned

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPSetFromOptions()`, `KSPSetOptionsPrefix()`, `KSPAppendOptionsPrefix()`

# External Links
$(_doc_external("KSP/KSPGetOptionsPrefix"))
"""
function KSPGetOptionsPrefix(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGetOptionsPrefix(petsclib::$UnionPetscLib, ksp::PetscKSP )
    
    prefix_ = Ref{Ptr{Int8}}()

    @chk ccall(
               (:KSPGetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Ptr{Cchar}}),
               ksp, prefix_,
              )

    if prefix_[] == C_NULL
        prefix = ""
    else
	    prefix = unsafe_string(prefix_[])
    end

	return prefix
end 

"""
	KSPMonitorSetFromOptions(petsclib::PetscLibType,ksp::PetscKSP, opt::String, name::String, ctx::Cvoid) 
Sets a monitor function and viewer appropriate for the type indicated by the user in the options database

Collective

Input Parameters:
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
$(_doc_external("KSP/KSPMonitorSetFromOptions"))
"""
function KSPMonitorSetFromOptions(petsclib::PetscLibType, ksp::PetscKSP, opt::String, name::String, ctx::Cvoid) end

@for_petsc function KSPMonitorSetFromOptions(petsclib::$UnionPetscLib, ksp::PetscKSP, opt::String, name::String, ctx::Cvoid )

    @chk ccall(
               (:KSPMonitorSetFromOptions, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cvoid}),
               ksp, opt, name, ctx,
              )


	return nothing
end 

"""
	KSPSetFromOptions(petsclib::PetscLibType,ksp::PetscKSP) 
Sets `KSP` options from the options database.
This routine must be called before `KSPSetUp()` if the user is to be
allowed to set the Krylov type.

Collective

Input Parameter:
- `ksp` - the Krylov space context

Options Database Keys:
- `-ksp_rtol rtol`                                                          - relative tolerance used in default determination of convergence, i.e.
if residual norm decreases by this factor than convergence is declared
- `-ksp_atol abstol`                                                        - absolute tolerance used in default convergence test, i.e. if residual
norm is less than this then convergence is declared
- `-ksp_divtol tol`                                                         - if residual norm increases by this factor than divergence is declared
- `-ksp_max_it`                                                             - maximum number of linear iterations
- `-ksp_min_it`                                                             - minimum number of linear iterations to use, defaults to zero

- `-ksp_reuse_preconditioner <true,false>`                                  - reuse the previously computed preconditioner

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
- `-ksp_monitor draw::draw_lg`                                              - plot residual norm at each iteration, see `KSPMonitorResidual()`
- `-ksp_monitor_true_residual`                                              - print the true l2 residual norm at each iteration, see `KSPMonitorTrueResidual()`
- `-all_ksp_monitor <optional filename>`                                    - print residual norm at each iteration for ALL KSP solves, regardless of their prefix. This is
useful for `PCFIELDSPLIT`, `PCMG`, etc that have inner solvers and
you wish to track the convergence of all the solvers
- `-ksp_monitor_solution [ascii binary or draw][:filename][:format option]` - plot solution at each iteration
- `-ksp_monitor_singular_value`                                             - monitor extreme singular values at each iteration
- `-ksp_converged_reason`                                                   - view the convergence state at the end of the solve
- `-ksp_use_explicittranspose`                                              - transpose the system explicitly in `KSPSolveTranspose()`
- `-ksp_error_if_not_converged`                                             - stop the program as soon as an error is detected in a `KSPSolve()`, `KSP_DIVERGED_ITS`
is not treated as an error on inner solves
- `-ksp_converged_rate`                                                     - view the convergence rate at the end of the solve

Level: beginner

-seealso: [](ch_ksp), `KSP`, `KSPSetOptionsPrefix()`, `KSPResetFromOptions()`, `KSPSetUseFischerGuess()`

# External Links
$(_doc_external("KSP/KSPSetFromOptions"))
"""
function KSPSetFromOptions(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPSetFromOptions(petsclib::$UnionPetscLib, ksp::PetscKSP )

    @chk ccall(
               (:KSPSetFromOptions, $petsc_library),
               PetscErrorCode,
               (CKSP,),
               ksp,
              )


	return nothing
end 

"""
	KSPResetFromOptions(petsclib::PetscLibType,ksp::PetscKSP) 
Sets `KSP` parameters from user options ONLY if the `KSP` was previously set from options

Collective

Input Parameter:
- `ksp` - the `KSP` context

Level: advanced

-seealso: [](ch_ksp), `KSPSetFromOptions()`, `KSPSetOptionsPrefix()`

# External Links
$(_doc_external("KSP/KSPResetFromOptions"))
"""
function KSPResetFromOptions(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPResetFromOptions(petsclib::$UnionPetscLib, ksp::PetscKSP )

    @chk ccall(
               (:KSPResetFromOptions, $petsc_library),
               PetscErrorCode,
               (CKSP,),
               ksp,
              )


	return nothing
end 

"""
	KSPFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the `KSP` package. It is
called from `PetscFinalize()`.

Level: developer

-seealso: [](ch_ksp), `PetscFinalize()`, `KSPInitializePackage()`

# External Links
$(_doc_external("KSP/KSPFinalizePackage"))
"""
function KSPFinalizePackage(petsclib::PetscLibType) end

@for_petsc function KSPFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:KSPFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	KSPInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `KSP` package. It is called
from `PetscDLLibraryRegister_petscksp()` when using dynamic libraries, and on the first call to `KSPCreate()`
when using shared or static libraries.

Level: developer

-seealso: [](ch_ksp), `PetscInitialize()`, `KSPFinalizePackage()`

# External Links
$(_doc_external("KSP/KSPInitializePackage"))
"""
function KSPInitializePackage(petsclib::PetscLibType) end

@for_petsc function KSPInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:KSPInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	ctx::Cvoid = KSPMonitorSAWsCreate(petsclib::PetscLibType,ksp::PetscKSP) 
create an SAWs monitor context for `KSP`

Collective

Input Parameter:
- `ksp` - `KSP` to monitor

Output Parameter:
- `ctx` - context for monitor

Level: developer

-seealso: [](ch_ksp), `KSP`, `KSPMonitorSet()`, `KSPMonitorSAWs()`, `KSPMonitorSAWsDestroy()`

# External Links
$(_doc_external("KSP/KSPMonitorSAWsCreate"))
"""
function KSPMonitorSAWsCreate(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPMonitorSAWsCreate(petsclib::$UnionPetscLib, ksp::PetscKSP )
	ctx_ = Ref{Cvoid}()

    @chk ccall(
               (:KSPMonitorSAWsCreate, $petsc_library),
               PetscErrorCode,
               (CKSP, Cvoid),
               ksp, ctx_,
              )

	ctx = ctx_[]

	return ctx
end 

"""
	KSPMonitorSAWsDestroy(petsclib::PetscLibType,ctx::Cvoid) 
destroy a monitor context created with `KSPMonitorSAWsCreate()`

Collective

Input Parameter:
- `ctx` - monitor context

Level: developer

-seealso: [](ch_ksp), `KSP`, `KSPMonitorSet()`, `KSPMonitorSAWsCreate()`

# External Links
$(_doc_external("KSP/KSPMonitorSAWsDestroy"))
"""
function KSPMonitorSAWsDestroy(petsclib::PetscLibType, ctx::Cvoid) end

@for_petsc function KSPMonitorSAWsDestroy(petsclib::$UnionPetscLib, ctx::Cvoid )

    @chk ccall(
               (:KSPMonitorSAWsDestroy, $petsc_library),
               PetscErrorCode,
               (Cvoid,),
               ctx,
              )


	return nothing
end 

"""
	KSPMonitorSAWs(petsclib::PetscLibType,ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, ctx::Cvoid) 
monitor `KSP` solution using SAWs

Logically Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `ctx`   - created with `KSPMonitorSAWsCreate()`

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPMonitorSet()`, `KSPMonitorSAWsCreate()`, `KSPMonitorSAWsDestroy()`, `KSPMonitorSingularValue()`, `KSPComputeExtremeSingularValues()`, `PetscViewerSAWsOpen()`

# External Links
$(_doc_external("KSP/KSPMonitorSAWs"))
"""
function KSPMonitorSAWs(petsclib::PetscLibType, ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, ctx::Cvoid) end

@for_petsc function KSPMonitorSAWs(petsclib::$UnionPetscLib, ksp::PetscKSP, n::$PetscInt, rnorm::$PetscReal, ctx::Cvoid )

    @chk ccall(
               (:KSPMonitorSAWs, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{Cvoid}),
               ksp, n, rnorm, ctx,
              )


	return nothing
end 

"""
	KSPGMRESSetPreAllocateVectors(petsclib::PetscLibType,ksp::PetscKSP) 
Causes `KSPGMRES` and `KSPFGMRES` to preallocate all its
needed work vectors at initial setup rather than the default, which
is to allocate several at a time when needed.

Logically Collective

Input Parameter:
- `ksp` - iterative context obtained from `KSPCreate()`

Options Database Key:
- `-ksp_gmres_preallocate` - Activates `KSPGmresSetPreAllocateVectors()`

Level: intermediate

-seealso: [](ch_ksp), `KSPGMRESSetRestart()`, `KSPGMRESSetOrthogonalization()`, `KSPGMRESGetOrthogonalization()`,
`VecMDot()`, `VecMAXPY()`

# External Links
$(_doc_external("KSP/KSPGMRESSetPreAllocateVectors"))
"""
function KSPGMRESSetPreAllocateVectors(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGMRESSetPreAllocateVectors(petsclib::$UnionPetscLib, ksp::PetscKSP )

    @chk ccall(
               (:KSPGMRESSetPreAllocateVectors, $petsc_library),
               PetscErrorCode,
               (CKSP,),
               ksp,
              )


	return nothing
end 

"""
	KSPGMRESClassicalGramSchmidtOrthogonalization(petsclib::PetscLibType,ksp::PetscKSP, it::PetscInt) 
This is the basic orthogonalization routine
using classical Gram-Schmidt with possible iterative refinement to improve the stability

Collective, No Fortran Support

Input Parameters:
- `ksp` - `KSP` object, must be associated with `KSPGMRES`, `KSPFGMRES`, or `KSPLGMRES` Krylov method
- `it`  - one less than the current GMRES restart iteration, i.e. the size of the Krylov space

Options Database Keys:
- `-ksp_gmres_classicalgramschmidt`                                             - Activates `KSPGMRESClassicalGramSchmidtOrthogonalization()`
- `-ksp_gmres_cgs_refinement_type <refine_never,refine_ifneeded,refine_always>` - determine if iterative refinement is
used to increase the stability of the classical Gram-Schmidt  orthogonalization.

Level: intermediate

-seealso: [](ch_ksp), `KSPGMRESCGSRefinementType`, `KSPGMRESSetOrthogonalization()`, `KSPGMRESSetCGSRefinementType()`,
`KSPGMRESGetCGSRefinementType()`, `KSPGMRESGetOrthogonalization()`, `KSPGMRESModifiedGramSchmidtOrthogonalization()`

# External Links
$(_doc_external("KSP/KSPGMRESClassicalGramSchmidtOrthogonalization"))
"""
function KSPGMRESClassicalGramSchmidtOrthogonalization(petsclib::PetscLibType, ksp::PetscKSP, it::PetscInt) end

@for_petsc function KSPGMRESClassicalGramSchmidtOrthogonalization(petsclib::$UnionPetscLib, ksp::PetscKSP, it::$PetscInt )

    @chk ccall(
               (:KSPGMRESClassicalGramSchmidtOrthogonalization, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, it,
              )


	return nothing
end 

"""
	KSPGMRESSetOrthogonalization(petsclib::PetscLibType,ksp::PetscKSP, fcn::external) 
Sets the orthogonalization routine used by `KSPGMRES` and `KSPFGMRES`.

Logically Collective

Input Parameters:
- `ksp` - iterative context obtained from `KSPCreate()`
- `fcn` - orthogonalization function

Calling sequence of `fcn`:
- `ksp` - the solver context
- `it`  - the current iteration

Options Database Keys:
- `-ksp_gmres_classicalgramschmidt` - Activates KSPGMRESClassicalGramSchmidtOrthogonalization() (default)
- `-ksp_gmres_modifiedgramschmidt`  - Activates KSPGMRESModifiedGramSchmidtOrthogonalization()

Level: intermediate

-seealso: [](ch_ksp), `KSPGMRESSetRestart()`, `KSPGMRESSetPreAllocateVectors()`,
`KSPGMRESSetCGSRefinementType()`, `KSPGMRESModifiedGramSchmidtOrthogonalization()`,
`KSPGMRESClassicalGramSchmidtOrthogonalization()`, `KSPGMRESGetCGSRefinementType()`

# External Links
$(_doc_external("KSP/KSPGMRESSetOrthogonalization"))
"""
function KSPGMRESSetOrthogonalization(petsclib::PetscLibType, ksp::PetscKSP, fcn::external) end

@for_petsc function KSPGMRESSetOrthogonalization(petsclib::$UnionPetscLib, ksp::PetscKSP, fcn::external )

    @chk ccall(
               (:KSPGMRESSetOrthogonalization, $petsc_library),
               PetscErrorCode,
               (CKSP, external),
               ksp, fcn,
              )


	return nothing
end 

"""
	KSPGMRESMonitorKrylov(petsclib::PetscLibType,ksp::PetscKSP, its::PetscInt, fgnorm::PetscReal, dummy::Cvoid) 
Calls `VecView()` to monitor each new direction in the `KSPGMRES` accumulated Krylov space.

Collective

Input Parameters:
- `ksp`    - the `KSP` context
- `its`    - iteration number
- `fgnorm` - 2-norm of residual (or gradient)
- `dummy`  - a collection of viewers created with `PetscViewersCreate()`

Options Database Key:
- `-ksp_gmres_krylov_monitor <bool>` - Plot the Krylov directions

Level: intermediate

-seealso: [](ch_ksp), `KSPGMRES`, `KSPMonitorSet()`, `KSPMonitorResidual()`, `VecView()`, `PetscViewersCreate()`, `PetscViewersDestroy()`

# External Links
$(_doc_external("KSP/KSPGMRESMonitorKrylov"))
"""
function KSPGMRESMonitorKrylov(petsclib::PetscLibType, ksp::PetscKSP, its::PetscInt, fgnorm::PetscReal, dummy::Cvoid) end

@for_petsc function KSPGMRESMonitorKrylov(petsclib::$UnionPetscLib, ksp::PetscKSP, its::$PetscInt, fgnorm::$PetscReal, dummy::Cvoid )

    @chk ccall(
               (:KSPGMRESMonitorKrylov, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{Cvoid}),
               ksp, its, fgnorm, dummy,
              )


	return nothing
end 

"""
	KSPGMRESSetCGSRefinementType(petsclib::PetscLibType,ksp::PetscKSP, type::KSPGMRESCGSRefinementType) 
Sets the type of iterative refinement to use
in the classical Gram-Schmidt orthogonalization used by `KSPGMRES` and other PETSc GMRES implementations.

Logically Collective

Input Parameters:
- `ksp`  - the Krylov space solver context
- `type` - the type of refinement
-seealso: [](ch_ksp), `KSPGMRES`, `KSPGMRESSetOrthogonalization()`, `KSPGMRESCGSRefinementType`, `KSPGMRESClassicalGramSchmidtOrthogonalization()`, `KSPGMRESGetCGSRefinementType()`,
`KSPGMRESGetOrthogonalization()`

# External Links
$(_doc_external("KSP/KSPGMRESSetCGSRefinementType"))
"""
function KSPGMRESSetCGSRefinementType(petsclib::PetscLibType, ksp::PetscKSP, type::KSPGMRESCGSRefinementType) end

@for_petsc function KSPGMRESSetCGSRefinementType(petsclib::$UnionPetscLib, ksp::PetscKSP, type::KSPGMRESCGSRefinementType )

    @chk ccall(
               (:KSPGMRESSetCGSRefinementType, $petsc_library),
               PetscErrorCode,
               (CKSP, KSPGMRESCGSRefinementType),
               ksp, type,
              )


	return nothing
end 

"""
	type::KSPGMRESCGSRefinementType = KSPGMRESGetCGSRefinementType(petsclib::PetscLibType,ksp::PetscKSP) 
Gets the type of iterative refinement to use
in the classical Gram-Schmidt orthogonalization used by `KSPGMRES` and other PETSc GMRES implementations.

Not Collective

Input Parameter:
- `ksp` - the Krylov space solver context

Output Parameter:
- `type` - the type of refinement

Level: intermediate

-seealso: [](ch_ksp), `KSPGMRES`, `KSPGMRESSetOrthogonalization()`, `KSPGMRESCGSRefinementType`, `KSPGMRESClassicalGramSchmidtOrthogonalization()`, `KSPGMRESSetCGSRefinementType()`,
`KSPGMRESGetOrthogonalization()`

# External Links
$(_doc_external("KSP/KSPGMRESGetCGSRefinementType"))
"""
function KSPGMRESGetCGSRefinementType(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGMRESGetCGSRefinementType(petsclib::$UnionPetscLib, ksp::PetscKSP )
	type_ = Ref{KSPGMRESCGSRefinementType}()

    @chk ccall(
               (:KSPGMRESGetCGSRefinementType, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPGMRESCGSRefinementType}),
               ksp, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	KSPGMRESSetRestart(petsclib::PetscLibType,ksp::PetscKSP, restart::PetscInt) 
Sets number of iterations at which GMRES (`KSPGMRES`, `KSPFGMRES`, `KSPPGMRES`, `KSPAGMRES`, `KSPDGMRES`, `KSPPIPEFGMRES`,
and `KSPLGMRES`) restarts.

Logically Collective

Input Parameters:
- `ksp`     - the Krylov space solver context
- `restart` - integer restart value, this corresponds to the number of iterations of GMRES to perform before restarting

Options Database Key:
- `-ksp_gmres_restart <positive integer>` - integer restart value

Level: intermediate

-seealso: [](ch_ksp), `KSPGMRES`, `KSPSetTolerances()`, `KSPGMRESSetOrthogonalization()`, `KSPGMRESSetPreAllocateVectors()`, `KSPGMRESGetRestart()`,
`KSPFGMRES`, `KSPLGMRES`, `KSPPGMRES`, `KSPAGMRES`, `KSPDGMRES`, `KSPPIPEFGMRES`

# External Links
$(_doc_external("KSP/KSPGMRESSetRestart"))
"""
function KSPGMRESSetRestart(petsclib::PetscLibType, ksp::PetscKSP, restart::PetscInt) end

@for_petsc function KSPGMRESSetRestart(petsclib::$UnionPetscLib, ksp::PetscKSP, restart::$PetscInt )

    @chk ccall(
               (:KSPGMRESSetRestart, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, restart,
              )


	return nothing
end 

"""
	restart::PetscInt = KSPGMRESGetRestart(petsclib::PetscLibType,ksp::PetscKSP) 
Gets number of iterations at which GMRES (`KSPGMRES`, `KSPFGMRES`, `KSPPGMRES`, `KSPAGMRES`, `KSPDGMRES`, `KSPPIPEFGMRES`,
and `KSPLGMRES`) restarts.

Not Collective

Input Parameter:
- `ksp` - the Krylov space solver context

Output Parameter:
- `restart` - integer restart value

Level: intermediate

-seealso: [](ch_ksp), `KSPGMRES`, `KSPSetTolerances()`, `KSPGMRESSetOrthogonalization()`, `KSPGMRESSetPreAllocateVectors()`, `KSPGMRESSetRestart()`,
`KSPFGMRES`, `KSPLGMRES`, `KSPPGMRES`, `KSPAGMRES`, `KSPDGMRES`, `KSPPIPEFGMRES`

# External Links
$(_doc_external("KSP/KSPGMRESGetRestart"))
"""
function KSPGMRESGetRestart(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGMRESGetRestart(petsclib::$UnionPetscLib, ksp::PetscKSP )
	restart_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPGMRESGetRestart, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscInt}),
               ksp, restart_,
              )

	restart = restart_[]

	return restart
end 

"""
	KSPGMRESSetHapTol(petsclib::PetscLibType,ksp::PetscKSP, tol::PetscReal) 
Sets the tolerance for detecting a happy ending in GMRES (`KSPGMRES`, `KSPFGMRES` and `KSPLGMRES` and others)

Logically Collective

Input Parameters:
- `ksp` - the Krylov space solver context
- `tol` - the tolerance for detecting a happy ending

Options Database Key:
- `-ksp_gmres_haptol <positive real value>` - set tolerance for determining happy breakdown

Level: intermediate

-seealso: [](ch_ksp), `KSPGMRES`, `KSPSetTolerances()`

# External Links
$(_doc_external("KSP/KSPGMRESSetHapTol"))
"""
function KSPGMRESSetHapTol(petsclib::PetscLibType, ksp::PetscKSP, tol::PetscReal) end

@for_petsc function KSPGMRESSetHapTol(petsclib::$UnionPetscLib, ksp::PetscKSP, tol::$PetscReal )

    @chk ccall(
               (:KSPGMRESSetHapTol, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscReal),
               ksp, tol,
              )


	return nothing
end 

"""
	KSPGMRESSetBreakdownTolerance(petsclib::PetscLibType,ksp::PetscKSP, tol::PetscReal) 
Sets the tolerance for determining divergence breakdown in `KSPGMRES` at restart.

Logically Collective

Input Parameters:
- `ksp` - the Krylov space solver context
- `tol` - the tolerance

Options Database Key:
- `-ksp_gmres_breakdown_tolerance <positive real value>` - set tolerance for determining divergence breakdown

Level: intermediate

-seealso: [](ch_ksp), `KSPGMRES`, `KSPSetTolerances()`, `KSPGMRESSetHapTol()`, `KSPConvergedReason`

# External Links
$(_doc_external("KSP/KSPGMRESSetBreakdownTolerance"))
"""
function KSPGMRESSetBreakdownTolerance(petsclib::PetscLibType, ksp::PetscKSP, tol::PetscReal) end

@for_petsc function KSPGMRESSetBreakdownTolerance(petsclib::$UnionPetscLib, ksp::PetscKSP, tol::$PetscReal )

    @chk ccall(
               (:KSPGMRESSetBreakdownTolerance, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscReal),
               ksp, tol,
              )


	return nothing
end 

"""
	KSPGMRESModifiedGramSchmidtOrthogonalization(petsclib::PetscLibType,ksp::PetscKSP, it::PetscInt) 
This is the basic orthogonalization routine
using modified Gram-Schmidt.

Collective, No Fortran Support

Input Parameters:
- `ksp` - `KSP` object, must be associated with `KSPGMRES`, `KSPFGMRES`, or `KSPLGMRES` Krylov method
- `it`  - one less than the current GMRES restart iteration, i.e. the size of the Krylov space

Options Database Key:
- `-ksp_gmres_modifiedgramschmidt` - Activates `KSPGMRESModifiedGramSchmidtOrthogonalization()`

Level: intermediate

-seealso: [](ch_ksp), `KSPGMRESSetOrthogonalization()`, `KSPGMRESClassicalGramSchmidtOrthogonalization()`, `KSPGMRESGetOrthogonalization()`

# External Links
$(_doc_external("KSP/KSPGMRESModifiedGramSchmidtOrthogonalization"))
"""
function KSPGMRESModifiedGramSchmidtOrthogonalization(petsclib::PetscLibType, ksp::PetscKSP, it::PetscInt) end

@for_petsc function KSPGMRESModifiedGramSchmidtOrthogonalization(petsclib::$UnionPetscLib, ksp::PetscKSP, it::$PetscInt )

    @chk ccall(
               (:KSPGMRESModifiedGramSchmidtOrthogonalization, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, it,
              )


	return nothing
end 

"""
	KSPPIPEFGMRESSetShift(petsclib::PetscLibType,ksp::PetscKSP, shift::PetscScalar) 
Set the shift parameter for the flexible, pipelined `KSPPIPEFGMRES` solver.

Logically Collective

Input Parameters:
- `ksp`   - the Krylov space context
- `shift` - the shift

Options Database Key:
- `-ksp_pipefgmres_shift <shift>` - set the shift parameter

Level: intermediate

-seealso: [](ch_ksp), `KSPPIPEFGMRES`, `KSPComputeEigenvalues()`

# External Links
$(_doc_external("KSP/KSPPIPEFGMRESSetShift"))
"""
function KSPPIPEFGMRESSetShift(petsclib::PetscLibType, ksp::PetscKSP, shift::PetscScalar) end

@for_petsc function KSPPIPEFGMRESSetShift(petsclib::$UnionPetscLib, ksp::PetscKSP, shift::$PetscScalar )

    @chk ccall(
               (:KSPPIPEFGMRESSetShift, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscScalar),
               ksp, shift,
              )


	return nothing
end 

"""
	KSPFGMRESSetModifyPC(petsclib::PetscLibType,ksp::PetscKSP, fcn::KSPFlexibleModifyPCFn, ctx::Cvoid, destroy::PetscCtxDestroyFn) 
Sets the routine used by `KSPFGMRES` to modify the preconditioner. [](sec_flexibleksp)

Logically Collective

Input Parameters:
- `ksp`     - iterative context obtained from `KSPCreate()`
- `fcn`     - function to modify the `PC`, see `KSPFlexibleModifyPCFn`
- `ctx`     - optional context
- `destroy` - optional context destroy routine

Options Database Keys:
- `-ksp_fgmres_modifypcnochange` - do not change the `PC`
- `-ksp_fgmres_modifypcksp`      - changes the inner KSP solver tolerances

Level: intermediate

-seealso: [](ch_ksp), [](sec_flexibleksp), `KSPFGMRES`, `KSPFlexibleModifyPCFn`, `KSPFlexibleSetModifyPC()`, `KSPFGMRESModifyPCNoChange()`, `KSPFGMRESModifyPCKSP()`

# External Links
$(_doc_external("KSP/KSPFGMRESSetModifyPC"))
"""
function KSPFGMRESSetModifyPC(petsclib::PetscLibType, ksp::PetscKSP, fcn::KSPFlexibleModifyPCFn, ctx::Cvoid, destroy::PetscCtxDestroyFn) end

@for_petsc function KSPFGMRESSetModifyPC(petsclib::$UnionPetscLib, ksp::PetscKSP, fcn::KSPFlexibleModifyPCFn, ctx::Cvoid, destroy::PetscCtxDestroyFn )

    @chk ccall(
               (:KSPFGMRESSetModifyPC, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPFlexibleModifyPCFn}, Ptr{Cvoid}, Ptr{PetscCtxDestroyFn}),
               ksp, fcn, ctx, destroy,
              )


	return nothing
end 

"""
	KSPFlexibleSetModifyPC(petsclib::PetscLibType,ksp::PetscKSP, fcn::KSPFlexibleModifyPCFn, ctx::Cvoid, destroy::PetscCtxDestroyFn) 
Sets the routine used by flexible `KSP` methods to modify the preconditioner. [](sec_flexibleksp)

Logically Collective

Input Parameters:
- `ksp`     - iterative context obtained from `KSPCreate()`
- `fcn`     - function to modify the `PC`, see `KSPFlexibleModifyPCFn`
- `ctx`     - optional context
- `destroy` - optional context destroy routine

Level: intermediate

-seealso: [](ch_ksp), [](sec_flexibleksp), `KSPFGMRES`, `KSPFGMRESModifyPCNoChange()`, `KSPFGMRESModifyPCKSP()`

# External Links
$(_doc_external("KSP/KSPFlexibleSetModifyPC"))
"""
function KSPFlexibleSetModifyPC(petsclib::PetscLibType, ksp::PetscKSP, fcn::KSPFlexibleModifyPCFn, ctx::Cvoid, destroy::PetscCtxDestroyFn) end

@for_petsc function KSPFlexibleSetModifyPC(petsclib::$UnionPetscLib, ksp::PetscKSP, fcn::KSPFlexibleModifyPCFn, ctx::Cvoid, destroy::PetscCtxDestroyFn )

    @chk ccall(
               (:KSPFlexibleSetModifyPC, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPFlexibleModifyPCFn}, Ptr{Cvoid}, Ptr{PetscCtxDestroyFn}),
               ksp, fcn, ctx, destroy,
              )


	return nothing
end 

"""
	KSPFGMRESModifyPCNoChange(petsclib::PetscLibType,ksp::PetscKSP, total_its::PetscInt, loc_its::PetscInt, res_norm::PetscReal, ctx::Cvoid) 
this is the default used by `KSPFMGMRES`

Input Parameters:
- `ksp`       - the ksp context being used.
- `total_its` - the total number of `KSPFGMRES` iterations that have occurred.
- `loc_its`   - the number of `KSPFGMRES` iterations since last restart.
- `res_norm`  - the current residual norm.
- `ctx`       - context variable, unused in this routine

Level: intermediate

-seealso: [](ch_ksp), [](sec_flexibleksp), `KSPFGMRES`, `KSPFlexibleModifyPCFn`, `KSPFGMRESSetModifyPC()`, `KSPFGMRESModifyPCKSP()`

# External Links
$(_doc_external("KSP/KSPFGMRESModifyPCNoChange"))
"""
function KSPFGMRESModifyPCNoChange(petsclib::PetscLibType, ksp::PetscKSP, total_its::PetscInt, loc_its::PetscInt, res_norm::PetscReal, ctx::Cvoid) end

@for_petsc function KSPFGMRESModifyPCNoChange(petsclib::$UnionPetscLib, ksp::PetscKSP, total_its::$PetscInt, loc_its::$PetscInt, res_norm::$PetscReal, ctx::Cvoid )

    @chk ccall(
               (:KSPFGMRESModifyPCNoChange, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscInt, $PetscReal, Ptr{Cvoid}),
               ksp, total_its, loc_its, res_norm, ctx,
              )


	return nothing
end 

"""
	KSPFGMRESModifyPCKSP(petsclib::PetscLibType,ksp::PetscKSP, total_its::PetscInt, loc_its::PetscInt, res_norm::PetscReal, ctx::Cvoid) 
modifies the attributes of the `KSPFGMRES` preconditioner, see [](sec_flexibleksp).

Input Parameters:
- `ksp`       - the ksp context being used.
- `total_its` - the total number of `KSPFGMRES` iterations that have occurred.
- `loc_its`   - the number of `KSPFGMRES` iterations since last restart.
- `res_norm`  - the current residual norm.
- `ctx`       - context, not used in this routine

Level: intermediate

-seealso: [](ch_ksp), [](sec_flexibleksp), `KSPFGMRES`, `KSPFlexibleModifyPCFn`, `KSPFGMRESSetModifyPC()`

# External Links
$(_doc_external("KSP/KSPFGMRESModifyPCKSP"))
"""
function KSPFGMRESModifyPCKSP(petsclib::PetscLibType, ksp::PetscKSP, total_its::PetscInt, loc_its::PetscInt, res_norm::PetscReal, ctx::Cvoid) end

@for_petsc function KSPFGMRESModifyPCKSP(petsclib::$UnionPetscLib, ksp::PetscKSP, total_its::$PetscInt, loc_its::$PetscInt, res_norm::$PetscReal, ctx::Cvoid )

    @chk ccall(
               (:KSPFGMRESModifyPCKSP, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscInt, $PetscReal, Ptr{Cvoid}),
               ksp, total_its, loc_its, res_norm, ctx,
              )


	return nothing
end 

"""
	KSPLGMRESSetAugDim(petsclib::PetscLibType,ksp::PetscKSP, dim::PetscInt) 
Set the number of error approximations to include in the approximation space (default is 2) for `KSPLGMRES`

Collective

Input Parameters:
- `ksp` - the `KSP` context
- `dim` - the number of vectors to use

Options Database Key:
- `-ksp_lgmres_augment dim` - the number of error approximations to include

Level: intermediate

-seealso: [](ch_ksp), `KSPLGMRES`, `KSPLGMRESSetConstant()`

# External Links
$(_doc_external("KSP/KSPLGMRESSetAugDim"))
"""
function KSPLGMRESSetAugDim(petsclib::PetscLibType, ksp::PetscKSP, dim::PetscInt) end

@for_petsc function KSPLGMRESSetAugDim(petsclib::$UnionPetscLib, ksp::PetscKSP, dim::$PetscInt )

    @chk ccall(
               (:KSPLGMRESSetAugDim, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, dim,
              )


	return nothing
end 

"""
	KSPLGMRESSetConstant(petsclib::PetscLibType,ksp::PetscKSP) 
keep the error approximation space a constant size for every restart cycle

Collective

Input Parameters:
- `ksp` - the `KSP` context

Options Database Key:
- `-ksp_lgmres_constant` - set the size to be constant

Level: intermediate

-seealso: [](ch_ksp), `KSPLGMRES`, `KSPLGMRESSetAugDim()`

# External Links
$(_doc_external("KSP/KSPLGMRESSetConstant"))
"""
function KSPLGMRESSetConstant(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPLGMRESSetConstant(petsclib::$UnionPetscLib, ksp::PetscKSP )

    @chk ccall(
               (:KSPLGMRESSetConstant, $petsc_library),
               PetscErrorCode,
               (CKSP,),
               ksp,
              )


	return nothing
end 

"""
	KSPChebyshevSetEigenvalues(petsclib::PetscLibType,ksp::PetscKSP, emax::PetscReal, emin::PetscReal) 
Sets estimates for the extreme eigenvalues of the preconditioned problem.

Logically Collective

Input Parameters:
- `ksp`  - the Krylov space context
- `emax` - the eigenvalue maximum estimate
- `emin` - the eigenvalue minimum estimate

Options Database Key:
- `-ksp_chebyshev_eigenvalues emin,emax` - extreme eigenvalues

Level: intermediate

-seealso: [](ch_ksp), `KSPCHEBYSHEV`, `KSPChebyshevEstEigSet()`,

# External Links
$(_doc_external("KSP/KSPChebyshevSetEigenvalues"))
"""
function KSPChebyshevSetEigenvalues(petsclib::PetscLibType, ksp::PetscKSP, emax::PetscReal, emin::PetscReal) end

@for_petsc function KSPChebyshevSetEigenvalues(petsclib::$UnionPetscLib, ksp::PetscKSP, emax::$PetscReal, emin::$PetscReal )

    @chk ccall(
               (:KSPChebyshevSetEigenvalues, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscReal, $PetscReal),
               ksp, emax, emin,
              )


	return nothing
end 

"""
	KSPChebyshevEstEigSet(petsclib::PetscLibType,ksp::PetscKSP, a::PetscReal, b::PetscReal, c::PetscReal, d::PetscReal) 
Automatically estimate the eigenvalues to use for Chebyshev

Logically Collective

Input Parameters:
- `ksp` - the Krylov space context
- `a`   - multiple of min eigenvalue estimate to use for min Chebyshev bound (or `PETSC_DECIDE`)
- `b`   - multiple of max eigenvalue estimate to use for min Chebyshev bound (or `PETSC_DECIDE`)
- `c`   - multiple of min eigenvalue estimate to use for max Chebyshev bound (or `PETSC_DECIDE`)
- `d`   - multiple of max eigenvalue estimate to use for max Chebyshev bound (or `PETSC_DECIDE`)

Options Database Key:
- `-ksp_chebyshev_esteig a,b,c,d` - estimate eigenvalues using a Krylov method, then use this transform for Chebyshev eigenvalue bounds

-seealso: [](ch_ksp), `KSPCHEBYSHEV`, `KSPChebyshevEstEigSetUseNoisy()`, `KSPChebyshevEstEigGetKSP()`

# External Links
$(_doc_external("KSP/KSPChebyshevEstEigSet"))
"""
function KSPChebyshevEstEigSet(petsclib::PetscLibType, ksp::PetscKSP, a::PetscReal, b::PetscReal, c::PetscReal, d::PetscReal) end

@for_petsc function KSPChebyshevEstEigSet(petsclib::$UnionPetscLib, ksp::PetscKSP, a::$PetscReal, b::$PetscReal, c::$PetscReal, d::$PetscReal )

    @chk ccall(
               (:KSPChebyshevEstEigSet, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscReal, $PetscReal, $PetscReal, $PetscReal),
               ksp, a, b, c, d,
              )


	return nothing
end 

"""
	KSPChebyshevEstEigSetUseNoisy(petsclib::PetscLibType,ksp::PetscKSP, use::PetscBool) 
use a noisy random number generated right

Logically Collective

Input Parameters:
- `ksp` - linear solver context
- `use` - `PETSC_TRUE` to use noisy

Options Database Key:
- `-ksp_chebyshev_esteig_noisy <true,false>` - Use noisy right-hand side for estimate

Level: intermediate

-seealso: [](ch_ksp), `KSPCHEBYSHEV`, `KSPChebyshevEstEigSet()`, `KSPChebyshevEstEigGetKSP()`

# External Links
$(_doc_external("KSP/KSPChebyshevEstEigSetUseNoisy"))
"""
function KSPChebyshevEstEigSetUseNoisy(petsclib::PetscLibType, ksp::PetscKSP, use::PetscBool) end

@for_petsc function KSPChebyshevEstEigSetUseNoisy(petsclib::$UnionPetscLib, ksp::PetscKSP, use::PetscBool )

    @chk ccall(
               (:KSPChebyshevEstEigSetUseNoisy, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, use,
              )


	return nothing
end 

"""
	KSPChebyshevEstEigGetKSP(petsclib::PetscLibType,ksp::PetscKSP, kspest::PetscKSP) 
Get the Krylov method context used to estimate the eigenvalues for the Chebyshev method.

Input Parameter:
- `ksp` - the Krylov space context

Output Parameter:
- `kspest` - the eigenvalue estimation Krylov space context

Level: advanced

-seealso: [](ch_ksp), `KSPCHEBYSHEV`, `KSPChebyshevEstEigSet()`

# External Links
$(_doc_external("KSP/KSPChebyshevEstEigGetKSP"))
"""
function KSPChebyshevEstEigGetKSP(petsclib::PetscLibType, ksp::PetscKSP, kspest::PetscKSP) end

@for_petsc function KSPChebyshevEstEigGetKSP(petsclib::$UnionPetscLib, ksp::PetscKSP, kspest::PetscKSP )
	kspest_ = Ref(kspest.ptr)

    @chk ccall(
               (:KSPChebyshevEstEigGetKSP, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{CKSP}),
               ksp, kspest_,
              )

	kspest.ptr = C_NULL

	return nothing
end 

"""
	KSPChebyshevSetKind(petsclib::PetscLibType,ksp::PetscKSP, kind::KSPChebyshevKind) 
set the kind of Chebyshev polynomial to use

Logically Collective

Input Parameters:
- `ksp`  - Linear solver context
- `kind` - The kind of Chebyshev polynomial to use, see `KSPChebyshevKind`, one of `KSP_CHEBYSHEV_FIRST`, `KSP_CHEBYSHEV_FOURTH`, or `KSP_CHEBYSHEV_OPT_FOURTH`

Options Database Key:
- `-ksp_chebyshev_kind <kind>` - which kind of Chebyshev polynomial to use

Level: intermediate

-seealso: [](ch_ksp), `KSPCHEBYSHEV`, `KSPChebyshevKind`, `KSPChebyshevGetKind()`, `KSP_CHEBYSHEV_FIRST`, `KSP_CHEBYSHEV_FOURTH`, `KSP_CHEBYSHEV_OPT_FOURTH`

# External Links
$(_doc_external("KSP/KSPChebyshevSetKind"))
"""
function KSPChebyshevSetKind(petsclib::PetscLibType, ksp::PetscKSP, kind::KSPChebyshevKind) end

@for_petsc function KSPChebyshevSetKind(petsclib::$UnionPetscLib, ksp::PetscKSP, kind::KSPChebyshevKind )

    @chk ccall(
               (:KSPChebyshevSetKind, $petsc_library),
               PetscErrorCode,
               (CKSP, KSPChebyshevKind),
               ksp, kind,
              )


	return nothing
end 

"""
	KSPChebyshevGetKind(petsclib::PetscLibType,ksp::PetscKSP, kind::KSPChebyshevKind) 
get the kind of Chebyshev polynomial to use

Logically Collective

Input Parameters:
- `ksp`  - Linear solver context
- `kind` - The kind of Chebyshev polynomial used

Level: intermediate

-seealso: [](ch_ksp), `KSPCHEBYSHEV`, `KSPChebyshevKind`, `KSPChebyshevSetKind()`, `KSP_CHEBYSHEV_FIRST`, `KSP_CHEBYSHEV_FOURTH`, `KSP_CHEBYSHEV_OPT_FOURTH`

# External Links
$(_doc_external("KSP/KSPChebyshevGetKind"))
"""
function KSPChebyshevGetKind(petsclib::PetscLibType, ksp::PetscKSP, kind::KSPChebyshevKind) end

@for_petsc function KSPChebyshevGetKind(petsclib::$UnionPetscLib, ksp::PetscKSP, kind::KSPChebyshevKind )

    @chk ccall(
               (:KSPChebyshevGetKind, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPChebyshevKind}),
               ksp, kind,
              )


	return nothing
end 

"""
	KSPFCGSetMmax(petsclib::PetscLibType,ksp::PetscKSP, mmax::PetscInt) 
set the maximum number of previous directions `KSPFCG` will store for orthogonalization

Logically Collective

Input Parameters:
- `ksp`  - the Krylov space context
- `mmax` - the maximum number of previous directions to orthogonalize against

Options Database Key:
- `-ksp_fcg_mmax <N>` - maximum number of search directions

Level: intermediate

-seealso: [](ch_ksp), `KSPFCG`, `KSPFCGGetTruncationType()`, `KSPFCGGetNprealloc()`, `KSPFCGetMmax()`

# External Links
$(_doc_external("KSP/KSPFCGSetMmax"))
"""
function KSPFCGSetMmax(petsclib::PetscLibType, ksp::PetscKSP, mmax::PetscInt) end

@for_petsc function KSPFCGSetMmax(petsclib::$UnionPetscLib, ksp::PetscKSP, mmax::$PetscInt )

    @chk ccall(
               (:KSPFCGSetMmax, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, mmax,
              )


	return nothing
end 

"""
	mmax::PetscInt = KSPFCGGetMmax(petsclib::PetscLibType,ksp::PetscKSP) 
get the maximum number of previous directions `KSPFCG` will store

Not Collective

Input Parameter:
- `ksp` - the Krylov space context

Output Parameter:
- `mmax` - the maximum number of previous directions allowed for orthogonalization

Level: intermediate

-seealso: [](ch_ksp), `KSPFCG`, `KSPFCGGetTruncationType()`, `KSPFCGGetNprealloc()`, `KSPFCGSetMmax()`

# External Links
$(_doc_external("KSP/KSPFCGGetMmax"))
"""
function KSPFCGGetMmax(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPFCGGetMmax(petsclib::$UnionPetscLib, ksp::PetscKSP )
	mmax_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPFCGGetMmax, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscInt}),
               ksp, mmax_,
              )

	mmax = mmax_[]

	return mmax
end 

"""
	KSPFCGSetNprealloc(petsclib::PetscLibType,ksp::PetscKSP, nprealloc::PetscInt) 
set the number of directions to preallocate with `KSPFCG`

Logically Collective

Input Parameters:
- `ksp`       - the Krylov space context
- `nprealloc` - the number of vectors to preallocate

Options Database Key:
- `-ksp_fcg_nprealloc <N>` - number of directions to preallocate

Level: advanced

-seealso: [](ch_ksp), `KSPFCG`, `KSPFCGGetTruncationType()`, `KSPFCGGetNprealloc()`, `KSPFCGSetMmax()`, `KSPFCGGetMmax()`

# External Links
$(_doc_external("KSP/KSPFCGSetNprealloc"))
"""
function KSPFCGSetNprealloc(petsclib::PetscLibType, ksp::PetscKSP, nprealloc::PetscInt) end

@for_petsc function KSPFCGSetNprealloc(petsclib::$UnionPetscLib, ksp::PetscKSP, nprealloc::$PetscInt )

    @chk ccall(
               (:KSPFCGSetNprealloc, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, nprealloc,
              )


	return nothing
end 

"""
	nprealloc::PetscInt = KSPFCGGetNprealloc(petsclib::PetscLibType,ksp::PetscKSP) 
get the number of directions preallocate by `KSPFCG`

Not Collective

Input Parameter:
- `ksp` - the Krylov space context

Output Parameter:
- `nprealloc` - the number of directions preallocated

Level: advanced

-seealso: [](ch_ksp), `KSPFCG`, `KSPFCGGetTruncationType()`, `KSPFCGSetNprealloc()`, `KSPFCGSetMmax()`, `KSPFCGGetMmax()`

# External Links
$(_doc_external("KSP/KSPFCGGetNprealloc"))
"""
function KSPFCGGetNprealloc(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPFCGGetNprealloc(petsclib::$UnionPetscLib, ksp::PetscKSP )
	nprealloc_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPFCGGetNprealloc, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscInt}),
               ksp, nprealloc_,
              )

	nprealloc = nprealloc_[]

	return nprealloc
end 

"""
	KSPFCGSetTruncationType(petsclib::PetscLibType,ksp::PetscKSP, truncstrat::KSPFCDTruncationType) 
specify how many of its stored previous directions `KSPFCG` uses during orthoganalization

Logically Collective

Input Parameters:
- `ksp`        - the Krylov space context
- `truncstrat` - the choice of strategy
-seealso: [](ch_ksp), `KSPFCDTruncationType`, `KSPFCGGetTruncationType()`, `KSPFCGSetNprealloc()`, `KSPFCGSetMmax()`, `KSPFCGGetMmax()`,
`KSP_FCD_TRUNC_TYPE_STANDARD`, `KSP_FCD_TRUNC_TYPE_NOTAY`

# External Links
$(_doc_external("KSP/KSPFCGSetTruncationType"))
"""
function KSPFCGSetTruncationType(petsclib::PetscLibType, ksp::PetscKSP, truncstrat::KSPFCDTruncationType) end

@for_petsc function KSPFCGSetTruncationType(petsclib::$UnionPetscLib, ksp::PetscKSP, truncstrat::KSPFCDTruncationType )

    @chk ccall(
               (:KSPFCGSetTruncationType, $petsc_library),
               PetscErrorCode,
               (CKSP, KSPFCDTruncationType),
               ksp, truncstrat,
              )


	return nothing
end 

"""
	truncstrat::KSPFCDTruncationType = KSPFCGGetTruncationType(petsclib::PetscLibType,ksp::PetscKSP) 
get the truncation strategy employed by `KSPFCG`

Not Collective

Input Parameter:
- `ksp` - the Krylov space context

Output Parameter:
- `truncstrat` - the strategy type

Level: intermediate

-seealso: [](ch_ksp), `KSPFCG`, `KSPFCGSetTruncationType()`, `KSPFCDTruncationType`, `KSP_FCD_TRUNC_TYPE_STANDARD`, `KSP_FCD_TRUNC_TYPE_NOTAY`

# External Links
$(_doc_external("KSP/KSPFCGGetTruncationType"))
"""
function KSPFCGGetTruncationType(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPFCGGetTruncationType(petsclib::$UnionPetscLib, ksp::PetscKSP )
	truncstrat_ = Ref{KSPFCDTruncationType}()

    @chk ccall(
               (:KSPFCGGetTruncationType, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPFCDTruncationType}),
               ksp, truncstrat_,
              )

	truncstrat = unsafe_string(truncstrat_[])

	return truncstrat
end 

"""
	KSPPIPEFCGSetMmax(petsclib::PetscLibType,ksp::PetscKSP, mmax::PetscInt) 
set the maximum number of previous directions `KSPPIPEFCG` will store for orthogonalization

Logically Collective

Input Parameters:
- `ksp`  - the Krylov space context
- `mmax` - the maximum number of previous directions to orthogonalize against

Options Database Key:
- `-ksp_pipefcg_mmax <N>` - maximum number of previous directions

Level: intermediate

-seealso: [](ch_ksp), `KSPPIPEFCG`, `KSPPIPEFCGSetTruncationType()`, `KSPPIPEFCGSetNprealloc()`, `KSPFCGSetMmax()`, `KSPFCGGetMmax()`

# External Links
$(_doc_external("KSP/KSPPIPEFCGSetMmax"))
"""
function KSPPIPEFCGSetMmax(petsclib::PetscLibType, ksp::PetscKSP, mmax::PetscInt) end

@for_petsc function KSPPIPEFCGSetMmax(petsclib::$UnionPetscLib, ksp::PetscKSP, mmax::$PetscInt )

    @chk ccall(
               (:KSPPIPEFCGSetMmax, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, mmax,
              )


	return nothing
end 

"""
	mmax::PetscInt = KSPPIPEFCGGetMmax(petsclib::PetscLibType,ksp::PetscKSP) 
get the maximum number of previous directions `KSPPIPEFCG` will store

Not Collective

Input Parameter:
- `ksp` - the Krylov space context

Output Parameter:
- `mmax` - the maximum number of previous directions allowed for orthogonalization

Level: intermediate

-seealso: [](ch_ksp), `KSPPIPEFCG`, `KSPPIPEFCGGetTruncationType()`, `KSPPIPEFCGGetNprealloc()`, `KSPPIPEFCGSetMmax()`, `KSPFCGGetMmax()`, `KSPFCGSetMmax()`

# External Links
$(_doc_external("KSP/KSPPIPEFCGGetMmax"))
"""
function KSPPIPEFCGGetMmax(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPPIPEFCGGetMmax(petsclib::$UnionPetscLib, ksp::PetscKSP )
	mmax_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPPIPEFCGGetMmax, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscInt}),
               ksp, mmax_,
              )

	mmax = mmax_[]

	return mmax
end 

"""
	KSPPIPEFCGSetNprealloc(petsclib::PetscLibType,ksp::PetscKSP, nprealloc::PetscInt) 
set the number of directions to preallocate with `KSPPIPEFCG`

Logically Collective

Input Parameters:
- `ksp`       - the Krylov space context
- `nprealloc` - the number of vectors to preallocate

Options Database Key:
- `-ksp_pipefcg_nprealloc <N>` - the number of vectors to preallocate

Level: advanced

-seealso: [](ch_ksp), `KSPPIPEFCG`, `KSPPIPEFCGSetTruncationType()`, `KSPPIPEFCGGetNprealloc()`, `KSPPIPEFCGSetMmax()`, `KSPPIPEFCGGetMmax()`

# External Links
$(_doc_external("KSP/KSPPIPEFCGSetNprealloc"))
"""
function KSPPIPEFCGSetNprealloc(petsclib::PetscLibType, ksp::PetscKSP, nprealloc::PetscInt) end

@for_petsc function KSPPIPEFCGSetNprealloc(petsclib::$UnionPetscLib, ksp::PetscKSP, nprealloc::$PetscInt )

    @chk ccall(
               (:KSPPIPEFCGSetNprealloc, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, nprealloc,
              )


	return nothing
end 

"""
	nprealloc::PetscInt = KSPPIPEFCGGetNprealloc(petsclib::PetscLibType,ksp::PetscKSP) 
get the number of directions to preallocate by `KSPPIPEFCG`

Not Collective

Input Parameter:
- `ksp` - the Krylov space context

Output Parameter:
- `nprealloc` - the number of directions preallocated

Level: advanced

-seealso: [](ch_ksp), `KSPPIPEFCG`, `KSPPIPEFCGGetTruncationType()`, `KSPPIPEFCGSetNprealloc()`, `KSPPIPEFCGSetMmax()`, `KSPPIPEFCGGetMmax()`

# External Links
$(_doc_external("KSP/KSPPIPEFCGGetNprealloc"))
"""
function KSPPIPEFCGGetNprealloc(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPPIPEFCGGetNprealloc(petsclib::$UnionPetscLib, ksp::PetscKSP )
	nprealloc_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPPIPEFCGGetNprealloc, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscInt}),
               ksp, nprealloc_,
              )

	nprealloc = nprealloc_[]

	return nprealloc
end 

"""
	KSPPIPEFCGSetTruncationType(petsclib::PetscLibType,ksp::PetscKSP, truncstrat::KSPFCDTruncationType) 
specify how many of its stored previous directions `KSPPIPEFCG` uses during orthoganalization

Logically Collective

Input Parameters:
- `ksp`        - the Krylov space context
- `truncstrat` - the choice of strategy
-seealso: [](ch_ksp), `KSPPIPEFCG`, `KSPPIPEFCGGetTruncationType`, `KSPFCDTruncationType`, `KSP_FCD_TRUNC_TYPE_STANDARD`, `KSP_FCD_TRUNC_TYPE_NOTAY`

# External Links
$(_doc_external("KSP/KSPPIPEFCGSetTruncationType"))
"""
function KSPPIPEFCGSetTruncationType(petsclib::PetscLibType, ksp::PetscKSP, truncstrat::KSPFCDTruncationType) end

@for_petsc function KSPPIPEFCGSetTruncationType(petsclib::$UnionPetscLib, ksp::PetscKSP, truncstrat::KSPFCDTruncationType )

    @chk ccall(
               (:KSPPIPEFCGSetTruncationType, $petsc_library),
               PetscErrorCode,
               (CKSP, KSPFCDTruncationType),
               ksp, truncstrat,
              )


	return nothing
end 

"""
	truncstrat::KSPFCDTruncationType = KSPPIPEFCGGetTruncationType(petsclib::PetscLibType,ksp::PetscKSP) 
get the truncation strategy employed by `KSPPIPEFCG`

Not Collective

Input Parameter:
- `ksp` - the Krylov space context

Output Parameter:
- `truncstrat` - the strategy type

Level: intermediate

-seealso: [](ch_ksp), `KSPPIPEFCG`, `KSPPIPEFCGSetTruncationType`, `KSPFCDTruncationType`, `KSP_FCD_TRUNC_TYPE_STANDARD`, `KSP_FCD_TRUNC_TYPE_NOTAY`

# External Links
$(_doc_external("KSP/KSPPIPEFCGGetTruncationType"))
"""
function KSPPIPEFCGGetTruncationType(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPPIPEFCGGetTruncationType(petsclib::$UnionPetscLib, ksp::PetscKSP )
	truncstrat_ = Ref{KSPFCDTruncationType}()

    @chk ccall(
               (:KSPPIPEFCGGetTruncationType, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPFCDTruncationType}),
               ksp, truncstrat_,
              )

	truncstrat = unsafe_string(truncstrat_[])

	return truncstrat
end 

"""
	KSPPythonSetType(petsclib::PetscLibType,ksp::PetscKSP, pyname::String) 
Initialize a `KSP` object to a type implemented in Python.

Collective

Input Parameters:
- `ksp`  - the linear solver `KSP` context.
- `pyname`  - full dotted Python name [package].module[.{class|function}]

Options Database Key:
- `-ksp_python_type <pyname>`  - python class

Level: intermediate

-seealso: [](ch_ksp), `KSPCreate()`, `KSPSetType()`, `KSPPYTHON`, `PetscPythonInitialize()`

# External Links
$(_doc_external("KSP/KSPPythonSetType"))
"""
function KSPPythonSetType(petsclib::PetscLibType, ksp::PetscKSP, pyname::String) end

@for_petsc function KSPPythonSetType(petsclib::$UnionPetscLib, ksp::PetscKSP, pyname::String )

    @chk ccall(
               (:KSPPythonSetType, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Cchar}),
               ksp, pyname,
              )


	return nothing
end 

"""
	pyname::String = KSPPythonGetType(petsclib::PetscLibType,ksp::PetscKSP) 
Get the type of a `KSP` object implemented in Python.

Not Collective

Input Parameter:
- `ksp`  - the linear solver `KSP` context.

Output Parameter:
- `pyname`  - full dotted Python name [package].module[.{class|function}]

Level: intermediate

-seealso: [](ch_ksp), `KSPCreate()`, `KSPSetType()`, `KSPPYTHON`, `PetscPythonInitialize()`, `KSPPythonSetType()`

# External Links
$(_doc_external("KSP/KSPPythonGetType"))
"""
function KSPPythonGetType(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPPythonGetType(petsclib::$UnionPetscLib, ksp::PetscKSP )
	pyname_ = Ref{Ptr{Cchar}}()

    @chk ccall(
               (:KSPPythonGetType, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Ptr{Cchar}}),
               ksp, pyname_,
              )

	pyname = unsafe_wrap(Array, pyname_[], VecGetLocalSize(petsclib, x); own = false)

	return pyname
end 

"""
	KSPCGSetType(petsclib::PetscLibType,ksp::PetscKSP, type::KSPCGType) 
Sets the variant of the conjugate gradient method to
use for solving a linear system with a complex coefficient matrix.
This option is irrelevant when solving a real system.

Logically Collective

Input Parameters:
- `ksp`  - the iterative context
- `type` - the variant of CG to use, one of
-seealso: [](ch_ksp), `KSP`, `KSPCG`

# External Links
$(_doc_external("KSP/KSPCGSetType"))
"""
function KSPCGSetType(petsclib::PetscLibType, ksp::PetscKSP, type::KSPCGType) end

@for_petsc function KSPCGSetType(petsclib::$UnionPetscLib, ksp::PetscKSP, type::KSPCGType )

    @chk ccall(
               (:KSPCGSetType, $petsc_library),
               PetscErrorCode,
               (CKSP, KSPCGType),
               ksp, type,
              )


	return nothing
end 

"""
	KSPCGUseSingleReduction(petsclib::PetscLibType,ksp::PetscKSP, flg::PetscBool) 
Merge the two inner products needed in `KSPCG` into a single `MPI_Allreduce()` call.

Logically Collective

Input Parameters:
- `ksp` - the iterative context
- `flg` - turn on or off the single reduction

Options Database Key:
- `-ksp_cg_single_reduction <bool>` - Merge inner products into single `MPI_Allreduce()`

Level: intermediate

-seealso: [](ch_ksp), [](sec_pipelineksp), `KSP`, `KSPCG`, `KSPGMRES`, `KSPPIPECG`, `KSPPIPECR`, `and KSPGROPPCG`

# External Links
$(_doc_external("KSP/KSPCGUseSingleReduction"))
"""
function KSPCGUseSingleReduction(petsclib::PetscLibType, ksp::PetscKSP, flg::PetscBool) end

@for_petsc function KSPCGUseSingleReduction(petsclib::$UnionPetscLib, ksp::PetscKSP, flg::PetscBool )

    @chk ccall(
               (:KSPCGUseSingleReduction, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flg,
              )


	return nothing
end 

"""
	KSPCGSetRadius(petsclib::PetscLibType,ksp::PetscKSP, radius::PetscReal) 
Sets the radius of the trust region used by the `KSPCG` when the solver is used inside `SNESNEWTONTR`

Logically Collective

Input Parameters:
- `ksp`    - the iterative context
- `radius` - the trust region radius (0 is the default that disable the use of the radius)

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPCG`, `KSPNASH`, `KSPSTCG`, `KSPGLTR`, `SNESNEWTONTR`

# External Links
$(_doc_external("KSP/KSPCGSetRadius"))
"""
function KSPCGSetRadius(petsclib::PetscLibType, ksp::PetscKSP, radius::PetscReal) end

@for_petsc function KSPCGSetRadius(petsclib::$UnionPetscLib, ksp::PetscKSP, radius::$PetscReal )

    @chk ccall(
               (:KSPCGSetRadius, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscReal),
               ksp, radius,
              )


	return nothing
end 

"""
	KSPCGSetObjectiveTarget(petsclib::PetscLibType,ksp::PetscKSP, obj::PetscReal) 
Sets the target value for the CG quadratic model

Logically Collective

Input Parameters:
- `ksp` - the iterative context
- `obj` - the objective value (0 is the default)

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPCG`, `KSPNASH`, `KSPSTCG`, `KSPGLTR`, `SNESNEWTONTR`

# External Links
$(_doc_external("KSP/KSPCGSetObjectiveTarget"))
"""
function KSPCGSetObjectiveTarget(petsclib::PetscLibType, ksp::PetscKSP, obj::PetscReal) end

@for_petsc function KSPCGSetObjectiveTarget(petsclib::$UnionPetscLib, ksp::PetscKSP, obj::$PetscReal )

    @chk ccall(
               (:KSPCGSetObjectiveTarget, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscReal),
               ksp, obj,
              )


	return nothing
end 

"""
	KSPCGGetNormD(petsclib::PetscLibType,ksp::PetscKSP, norm_d::PetscReal) 
Get norm of the direction when the solver is used inside `SNESNEWTONTR`

Not collective

Input Parameters:
- `ksp`    - the iterative context
- `norm_d` - the norm of the direction

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPCG`, `KSPNASH`, `KSPSTCG`, `KSPGLTR`, `SNESNEWTONTR`

# External Links
$(_doc_external("KSP/KSPCGGetNormD"))
"""
function KSPCGGetNormD(petsclib::PetscLibType, ksp::PetscKSP, norm_d::PetscReal) end

@for_petsc function KSPCGGetNormD(petsclib::$UnionPetscLib, ksp::PetscKSP, norm_d::$PetscReal )

    @chk ccall(
               (:KSPCGGetNormD, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscReal}),
               ksp, norm_d,
              )


	return nothing
end 

"""
	KSPCGGetObjFcn(petsclib::PetscLibType,ksp::PetscKSP, o_fcn::PetscReal) 
Get the conjugate gradient objective function value

Not collective

Input Parameters:
- `ksp`   - the iterative context
- `o_fcn` - the objective function value

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPCG`, `KSPNASH`, `KSPSTCG`, `KSPGLTR`, `KSPMonitorSet`

# External Links
$(_doc_external("KSP/KSPCGGetObjFcn"))
"""
function KSPCGGetObjFcn(petsclib::PetscLibType, ksp::PetscKSP, o_fcn::PetscReal) end

@for_petsc function KSPCGGetObjFcn(petsclib::$UnionPetscLib, ksp::PetscKSP, o_fcn::$PetscReal )

    @chk ccall(
               (:KSPCGGetObjFcn, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscReal}),
               ksp, o_fcn,
              )


	return nothing
end 

"""
	e_min::PetscReal = KSPGLTRGetMinEig(petsclib::PetscLibType,ksp::PetscKSP) 
Get minimum eigenvalue computed by `KSPGLTR`

Collective

Input Parameter:
- `ksp` - the iterative context

Output Parameter:
- `e_min` - the minimum eigenvalue

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPGLTR`, `KSPGLTRGetLambda()`

# External Links
$(_doc_external("KSP/KSPGLTRGetMinEig"))
"""
function KSPGLTRGetMinEig(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGLTRGetMinEig(petsclib::$UnionPetscLib, ksp::PetscKSP )
	e_min_ = Ref{$PetscReal}()

    @chk ccall(
               (:KSPGLTRGetMinEig, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscReal}),
               ksp, e_min_,
              )

	e_min = e_min_[]

	return e_min
end 

"""
	lambda::PetscReal = KSPGLTRGetLambda(petsclib::PetscLibType,ksp::PetscKSP) 
Get the multiplier on the trust

Not Collective

Input Parameter:
- `ksp` - the iterative context

Output Parameter:
- `lambda` - the multiplier

Level: advanced

-seealso: [](ch_ksp), `KSP`, `KSPGLTR`, `KSPGLTRGetMinEig()`

# External Links
$(_doc_external("KSP/KSPGLTRGetLambda"))
"""
function KSPGLTRGetLambda(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGLTRGetLambda(petsclib::$UnionPetscLib, ksp::PetscKSP )
	lambda_ = Ref{$PetscReal}()

    @chk ccall(
               (:KSPGLTRGetLambda, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscReal}),
               ksp, lambda_,
              )

	lambda = lambda_[]

	return lambda
end 

"""
	KSPHPDDMSetDeflationMat(petsclib::PetscLibType,ksp::PetscKSP, U::PetscMat) 

# External Links
$(_doc_external("KSP/KSPHPDDMSetDeflationMat"))
"""
function KSPHPDDMSetDeflationMat(petsclib::PetscLibType, ksp::PetscKSP, U::PetscMat) end

@for_petsc function KSPHPDDMSetDeflationMat(petsclib::$UnionPetscLib, ksp::PetscKSP, U::PetscMat )

    @chk ccall(
               (:KSPHPDDMSetDeflationMat, $petsc_library),
               PetscErrorCode,
               (CKSP, CMat),
               ksp, U,
              )


	return nothing
end 

"""
	KSPHPDDMGetDeflationMat(petsclib::PetscLibType,ksp::PetscKSP, U::PetscMat) 

# External Links
$(_doc_external("KSP/KSPHPDDMGetDeflationMat"))
"""
function KSPHPDDMGetDeflationMat(petsclib::PetscLibType, ksp::PetscKSP, U::PetscMat) end

@for_petsc function KSPHPDDMGetDeflationMat(petsclib::$UnionPetscLib, ksp::PetscKSP, U::PetscMat )
	U_ = Ref(U.ptr)

    @chk ccall(
               (:KSPHPDDMGetDeflationMat, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{CMat}),
               ksp, U_,
              )

	U.ptr = C_NULL

	return nothing
end 

"""
	KSPHPDDMSetType(petsclib::PetscLibType,ksp::PetscKSP, type::KSPHPDDMType) 

# External Links
$(_doc_external("KSP/KSPHPDDMSetType"))
"""
function KSPHPDDMSetType(petsclib::PetscLibType, ksp::PetscKSP, type::KSPHPDDMType) end

@for_petsc function KSPHPDDMSetType(petsclib::$UnionPetscLib, ksp::PetscKSP, type::KSPHPDDMType )

    @chk ccall(
               (:KSPHPDDMSetType, $petsc_library),
               PetscErrorCode,
               (CKSP, KSPHPDDMType),
               ksp, type,
              )


	return nothing
end 

"""
	type::KSPHPDDMType = KSPHPDDMGetType(petsclib::PetscLibType,ksp::PetscKSP) 

# External Links
$(_doc_external("KSP/KSPHPDDMGetType"))
"""
function KSPHPDDMGetType(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPHPDDMGetType(petsclib::$UnionPetscLib, ksp::PetscKSP )
	type_ = Ref{KSPHPDDMType}()

    @chk ccall(
               (:KSPHPDDMGetType, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPHPDDMType}),
               ksp, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	KSPRichardsonSetScale(petsclib::PetscLibType,ksp::PetscKSP, scale::PetscReal) 
Set the damping factor; if this routine is not called, the factor defaults to 1.0.

Logically Collective

Input Parameters:
- `ksp`   - the iterative context
- `scale` - the damping factor

Options Database Key:
- `-ksp_richardson_scale <scale>` - Set the scale factor

Level: intermediate

-seealso: [](ch_ksp), `KSPRICHARDSON`, `KSPRichardsonSetSelfScale()`

# External Links
$(_doc_external("KSP/KSPRichardsonSetScale"))
"""
function KSPRichardsonSetScale(petsclib::PetscLibType, ksp::PetscKSP, scale::PetscReal) end

@for_petsc function KSPRichardsonSetScale(petsclib::$UnionPetscLib, ksp::PetscKSP, scale::$PetscReal )

    @chk ccall(
               (:KSPRichardsonSetScale, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscReal),
               ksp, scale,
              )


	return nothing
end 

"""
	KSPRichardsonSetSelfScale(petsclib::PetscLibType,ksp::PetscKSP, scale::PetscBool) 
Sets Richardson to automatically determine optimal scaling at each iteration to minimize the 2
preconditioned residual

Logically Collective

Input Parameters:
- `ksp`   - the iterative context
- `scale` - `PETSC_TRUE` or the default of `PETSC_FALSE`

Options Database Key:
- `-ksp_richardson_self_scale` - Use self-scaling

Level: intermediate

-seealso: [](ch_ksp), `KSPRICHARDSON`, `KSPRichardsonSetScale()`

# External Links
$(_doc_external("KSP/KSPRichardsonSetSelfScale"))
"""
function KSPRichardsonSetSelfScale(petsclib::PetscLibType, ksp::PetscKSP, scale::PetscBool) end

@for_petsc function KSPRichardsonSetSelfScale(petsclib::$UnionPetscLib, ksp::PetscKSP, scale::PetscBool )

    @chk ccall(
               (:KSPRichardsonSetSelfScale, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, scale,
              )


	return nothing
end 

"""
	KSPFETIDPSetPressureOperator(petsclib::PetscLibType,ksp::PetscKSP, P::PetscMat) 
Sets the operator used to set up the pressure preconditioner for the saddle point `KSPFETIDP` solver,

Collective

Input Parameters:
- `ksp` - the `KSPFETIDP` solver
- `P`   - the linear operator to be preconditioned, usually the mass matrix.

Level: advanced

-seealso: [](ch_ksp), `KSPFETIDP`, `MATIS`, `PCBDDC`, `KSPFETIDPGetInnerBDDC()`, `KSPFETIDPGetInnerKSP()`, `KSPSetOperators()`

# External Links
$(_doc_external("KSP/KSPFETIDPSetPressureOperator"))
"""
function KSPFETIDPSetPressureOperator(petsclib::PetscLibType, ksp::PetscKSP, P::PetscMat) end

@for_petsc function KSPFETIDPSetPressureOperator(petsclib::$UnionPetscLib, ksp::PetscKSP, P::PetscMat )

    @chk ccall(
               (:KSPFETIDPSetPressureOperator, $petsc_library),
               PetscErrorCode,
               (CKSP, CMat),
               ksp, P,
              )


	return nothing
end 

"""
	KSPFETIDPGetInnerKSP(petsclib::PetscLibType,ksp::PetscKSP, innerksp::PetscKSP) 
Gets the `KSP` object for the Lagrange multipliers from inside a `KSPFETIDP`

Input Parameter:
- `ksp` - the `KSPFETIDP`

Output Parameter:
- `innerksp` - the `KSP` for the multipliers

Level: advanced

-seealso: [](ch_ksp), `KSPFETIDP`, `MATIS`, `PCBDDC`, `KSPFETIDPSetInnerBDDC()`, `KSPFETIDPGetInnerBDDC()`

# External Links
$(_doc_external("KSP/KSPFETIDPGetInnerKSP"))
"""
function KSPFETIDPGetInnerKSP(petsclib::PetscLibType, ksp::PetscKSP, innerksp::PetscKSP) end

@for_petsc function KSPFETIDPGetInnerKSP(petsclib::$UnionPetscLib, ksp::PetscKSP, innerksp::PetscKSP )
	innerksp_ = Ref(innerksp.ptr)

    @chk ccall(
               (:KSPFETIDPGetInnerKSP, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{CKSP}),
               ksp, innerksp_,
              )

	innerksp.ptr = C_NULL

	return nothing
end 

"""
	KSPFETIDPGetInnerBDDC(petsclib::PetscLibType,ksp::PetscKSP, pc::PC) 
Gets the `PCBDDC` preconditioner used to set up the `KSPFETIDP` matrix for the Lagrange multipliers

Input Parameter:
- `ksp` - the `KSPFETIDP` Krylov solver

Output Parameter:
- `pc` - the `PCBDDC` preconditioner

Level: advanced

-seealso: [](ch_ksp), `MATIS`, `PCBDDC`, `KSPFETIDP`, `KSPFETIDPSetInnerBDDC()`, `KSPFETIDPGetInnerKSP()`

# External Links
$(_doc_external("KSP/KSPFETIDPGetInnerBDDC"))
"""
function KSPFETIDPGetInnerBDDC(petsclib::PetscLibType, ksp::PetscKSP, pc::PC) end

@for_petsc function KSPFETIDPGetInnerBDDC(petsclib::$UnionPetscLib, ksp::PetscKSP, pc::PC )

    @chk ccall(
               (:KSPFETIDPGetInnerBDDC, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{PC}),
               ksp, pc,
              )


	return nothing
end 

"""
	KSPFETIDPSetInnerBDDC(petsclib::PetscLibType,ksp::PetscKSP, pc::PC) 
Provides the `PCBDDC` preconditioner used to set up the `KSPFETIDP` matrix for the Lagrange multipliers

Collective

Input Parameters:
- `ksp` - the `KSPFETIDP` Krylov solver
- `pc`  - the `PCBDDC` preconditioner

Level: advanced

-seealso: [](ch_ksp), `MATIS`, `PCBDDC`, `KSPFETIDPGetInnerBDDC()`, `KSPFETIDPGetInnerKSP()`

# External Links
$(_doc_external("KSP/KSPFETIDPSetInnerBDDC"))
"""
function KSPFETIDPSetInnerBDDC(petsclib::PetscLibType, ksp::PetscKSP, pc::PC) end

@for_petsc function KSPFETIDPSetInnerBDDC(petsclib::$UnionPetscLib, ksp::PetscKSP, pc::PC )

    @chk ccall(
               (:KSPFETIDPSetInnerBDDC, $petsc_library),
               PetscErrorCode,
               (CKSP, PC),
               ksp, pc,
              )


	return nothing
end 

"""
	KSPBCGSLSetXRes(petsclib::PetscLibType,ksp::PetscKSP, delta::PetscReal) 
Sets the parameter governing when
exact residuals will be used instead of computed residuals for `KSPCBGSL`.

Logically Collective

Input Parameters:
- `ksp`   - iterative context of type `KSPBCGSL`
- `delta` - computed residuals are used alone when delta is not positive

Options Database Key:
- `-ksp_bcgsl_xres delta` - Threshold used to decide when to refresh computed residuals

Level: intermediate

-seealso: [](ch_ksp), `KSPBCGSLSetEll()`, `KSPBCGSLSetPol()`, `KSP`, `KSPCBGSL`, `KSPBCGSLSetUsePseudoinverse()`

# External Links
$(_doc_external("KSP/KSPBCGSLSetXRes"))
"""
function KSPBCGSLSetXRes(petsclib::PetscLibType, ksp::PetscKSP, delta::PetscReal) end

@for_petsc function KSPBCGSLSetXRes(petsclib::$UnionPetscLib, ksp::PetscKSP, delta::$PetscReal )

    @chk ccall(
               (:KSPBCGSLSetXRes, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscReal),
               ksp, delta,
              )


	return nothing
end 

"""
	KSPBCGSLSetUsePseudoinverse(petsclib::PetscLibType,ksp::PetscKSP, use_pinv::PetscBool) 
Use pseudoinverse (via SVD) to solve polynomial part of the update in `KSPCBGSL` solver

Logically Collective

Input Parameters:
- `ksp`      - iterative context of type `KSPCBGSL`
- `use_pinv` - set to `PETSC_TRUE` when using pseudoinverse

Options Database Key:
- `-ksp_bcgsl_pinv <true,false>` - use pseudoinverse

Level: intermediate

-seealso: [](ch_ksp), `KSPBCGSLSetEll()`, `KSP`, `KSPCBGSL`, `KSPBCGSLSetPol()`, `KSPBCGSLSetXRes()`

# External Links
$(_doc_external("KSP/KSPBCGSLSetUsePseudoinverse"))
"""
function KSPBCGSLSetUsePseudoinverse(petsclib::PetscLibType, ksp::PetscKSP, use_pinv::PetscBool) end

@for_petsc function KSPBCGSLSetUsePseudoinverse(petsclib::$UnionPetscLib, ksp::PetscKSP, use_pinv::PetscBool )

    @chk ccall(
               (:KSPBCGSLSetUsePseudoinverse, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, use_pinv,
              )


	return nothing
end 

"""
	KSPBCGSLSetPol(petsclib::PetscLibType,ksp::PetscKSP, uMROR::PetscBool) 
Sets the type of polynomial part that will
be used in the `KSPCBGSL` `KSPSolve()`

Logically Collective

Input Parameters:
- `ksp`   - iterative context of type `KSPCBGSL`
- `uMROR` - set to `PETSC_TRUE` when the polynomial is a convex combination of an MR and an OR step.

Options Database Keys:
- `-ksp_bcgsl_cxpoly` - use enhanced polynomial
- `-ksp_bcgsl_mrpoly` - use standard polynomial

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPBCGSL`, `KSPCreate()`, `KSPSetType()`, `KSPCBGSL`, `KSPBCGSLSetUsePseudoinverse()`, `KSPBCGSLSetEll()`, `KSPBCGSLSetXRes()`

# External Links
$(_doc_external("KSP/KSPBCGSLSetPol"))
"""
function KSPBCGSLSetPol(petsclib::PetscLibType, ksp::PetscKSP, uMROR::PetscBool) end

@for_petsc function KSPBCGSLSetPol(petsclib::$UnionPetscLib, ksp::PetscKSP, uMROR::PetscBool )

    @chk ccall(
               (:KSPBCGSLSetPol, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, uMROR,
              )


	return nothing
end 

"""
	KSPBCGSLSetEll(petsclib::PetscLibType,ksp::PetscKSP, ell::PetscInt) 
Sets the number of search directions to use in the `KSPBCGSL` Krylov solver

Logically Collective

Input Parameters:
- `ksp` - iterative context, `KSP`, of type `KSPBCGSL`
- `ell` - number of search directions to use

Options Database Key:
- `-ksp_bcgsl_ell ell` - Number of Krylov search directions

Level: intermediate

-seealso: [](ch_ksp), `KSPBCGSLSetUsePseudoinverse()`, `KSP`, `KSPBCGSL`, `KSPBCGSLSetPol()`, `KSPBCGSLSetXRes()`

# External Links
$(_doc_external("KSP/KSPBCGSLSetEll"))
"""
function KSPBCGSLSetEll(petsclib::PetscLibType, ksp::PetscKSP, ell::PetscInt) end

@for_petsc function KSPBCGSLSetEll(petsclib::$UnionPetscLib, ksp::PetscKSP, ell::$PetscInt )

    @chk ccall(
               (:KSPBCGSLSetEll, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, ell,
              )


	return nothing
end 

"""
	KSPQCGSetTrustRegionRadius(petsclib::PetscLibType,ksp::PetscKSP, delta::PetscReal) 
Sets the radius of the trust region for `KSPQCG`

Logically Collective

Input Parameters:
- `ksp`   - the iterative context
- `delta` - the trust region radius (Infinity is the default)

Options Database Key:
- `-ksp_qcg_trustregionradius <delta>` - trust region radius

Level: advanced

-seealso: [](ch_ksp), `KSPQCG`, `KSPQCGGetTrialStepNorm()`

# External Links
$(_doc_external("KSP/KSPQCGSetTrustRegionRadius"))
"""
function KSPQCGSetTrustRegionRadius(petsclib::PetscLibType, ksp::PetscKSP, delta::PetscReal) end

@for_petsc function KSPQCGSetTrustRegionRadius(petsclib::$UnionPetscLib, ksp::PetscKSP, delta::$PetscReal )

    @chk ccall(
               (:KSPQCGSetTrustRegionRadius, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscReal),
               ksp, delta,
              )


	return nothing
end 

"""
	tsnorm::PetscReal = KSPQCGGetTrialStepNorm(petsclib::PetscLibType,ksp::PetscKSP) 
Gets the norm of a trial step vector in `KSPQCG`.  The WCG step may be
constrained, so this is not necessarily the length of the ultimate step taken in `KSPQCG`.

Not Collective

Input Parameter:
- `ksp` - the iterative context

Output Parameter:
- `tsnorm` - the norm

Level: advanced

-seealso: [](ch_ksp), `KSPQCG`, `KSPQCGSetTrustRegionRadius()`

# External Links
$(_doc_external("KSP/KSPQCGGetTrialStepNorm"))
"""
function KSPQCGGetTrialStepNorm(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPQCGGetTrialStepNorm(petsclib::$UnionPetscLib, ksp::PetscKSP )
	tsnorm_ = Ref{$PetscReal}()

    @chk ccall(
               (:KSPQCGGetTrialStepNorm, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscReal}),
               ksp, tsnorm_,
              )

	tsnorm = tsnorm_[]

	return tsnorm
end 

"""
	quadratic::PetscReal = KSPQCGGetQuadratic(petsclib::PetscLibType,ksp::PetscKSP) 
Gets the value of the quadratic function, evaluated at the new iterate

Collective

Input Parameter:
- `ksp` - the iterative context

Output Parameter:
- `quadratic` - the quadratic function evaluated at the new iterate

Level: advanced

-seealso: [](ch_ksp), `KSPQCG`

# External Links
$(_doc_external("KSP/KSPQCGGetQuadratic"))
"""
function KSPQCGGetQuadratic(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPQCGGetQuadratic(petsclib::$UnionPetscLib, ksp::PetscKSP )
	quadratic_ = Ref{$PetscReal}()

    @chk ccall(
               (:KSPQCGGetQuadratic, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscReal}),
               ksp, quadratic_,
              )

	quadratic = quadratic_[]

	return quadratic
end 

"""
	KSPGCRSetModifyPC(petsclib::PetscLibType,ksp::PetscKSP, fnc::KSPFlexibleModifyPCFn, ctx::Cvoid, destroy::PetscCtxDestroyFn) 
Sets the routine used by `KSPGCR` to modify the preconditioner for each iteration

Logically Collective

Input Parameters:
- `ksp`      - iterative context obtained from `KSPCreate()`
- `function` - user defined function to modify the preconditioner, see `KSPFlexibleModifyPCFn`
- `ctx`      - user provided context for the modify preconditioner function
- `destroy`  - the function to use to destroy the user provided application context.

Level: intermediate

-seealso: [](ch_ksp), `KSP`, `KSPGCR`, `KSPFlexibleModifyPCFn`, `KSPFGMRESModifyPCFn`, [](sec_flexibleksp)

# External Links
$(_doc_external("KSP/KSPGCRSetModifyPC"))
"""
function KSPGCRSetModifyPC(petsclib::PetscLibType, ksp::PetscKSP, fnc::KSPFlexibleModifyPCFn, ctx::Cvoid, destroy::PetscCtxDestroyFn) end

@for_petsc function KSPGCRSetModifyPC(petsclib::$UnionPetscLib, ksp::PetscKSP, fnc::KSPFlexibleModifyPCFn, ctx::Cvoid, destroy::PetscCtxDestroyFn )

    @chk ccall(
               (:KSPGCRSetModifyPC, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPFlexibleModifyPCFn}, Ptr{Cvoid}, Ptr{PetscCtxDestroyFn}),
               ksp, fnc, ctx, destroy,
              )


	return nothing
end 

"""
	KSPGCRSetRestart(petsclib::PetscLibType,ksp::PetscKSP, restart::PetscInt) 
Sets number of iterations at which `KSPGCR` restarts.

Not Collective

Input Parameters:
- `ksp`     - the Krylov space context
- `restart` - integer restart value

Options Database Key:
- `-ksp_gcr_restart <restart>` - the number of stored vectors to orthogonalize against

Level: intermediate

-seealso: [](ch_ksp), `KSPGCR`, `KSPSetTolerances()`, `KSPGCRGetRestart()`, `KSPGMRESSetRestart()`

# External Links
$(_doc_external("KSP/KSPGCRSetRestart"))
"""
function KSPGCRSetRestart(petsclib::PetscLibType, ksp::PetscKSP, restart::PetscInt) end

@for_petsc function KSPGCRSetRestart(petsclib::$UnionPetscLib, ksp::PetscKSP, restart::$PetscInt )

    @chk ccall(
               (:KSPGCRSetRestart, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, restart,
              )


	return nothing
end 

"""
	restart::PetscInt = KSPGCRGetRestart(petsclib::PetscLibType,ksp::PetscKSP) 
Gets number of iterations at which `KSPGCR` restarts.

Not Collective

Input Parameter:
- `ksp` - the Krylov space context

Output Parameter:
- `restart` - integer restart value

Level: intermediate

-seealso: [](ch_ksp), `KSPGCR`, `KSPSetTolerances()`, `KSPGCRSetRestart()`, `KSPGMRESGetRestart()`

# External Links
$(_doc_external("KSP/KSPGCRGetRestart"))
"""
function KSPGCRGetRestart(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPGCRGetRestart(petsclib::$UnionPetscLib, ksp::PetscKSP )
	restart_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPGCRGetRestart, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscInt}),
               ksp, restart_,
              )

	restart = restart_[]

	return restart
end 

"""
	KSPPIPEGCRSetUnrollW(petsclib::PetscLibType,ksp::PetscKSP, unroll_w::PetscBool) 
Set to `PETSC_TRUE` to use `KSPPIPEGCR` with unrolling of the w vector

Logically Collective

Input Parameters:
- `ksp`      - the Krylov space context
- `unroll_w` - use unrolling

Level: intermediate

Options Database Key:
- `-ksp_pipegcr_unroll_w <bool>` - use unrolling

-seealso: [](ch_ksp), `KSPPIPEGCR`, `KSPPIPEGCRSetTruncationType()`, `KSPPIPEGCRSetNprealloc()`, `KSPPIPEGCRGetUnrollW()`

# External Links
$(_doc_external("KSP/KSPPIPEGCRSetUnrollW"))
"""
function KSPPIPEGCRSetUnrollW(petsclib::PetscLibType, ksp::PetscKSP, unroll_w::PetscBool) end

@for_petsc function KSPPIPEGCRSetUnrollW(petsclib::$UnionPetscLib, ksp::PetscKSP, unroll_w::PetscBool )

    @chk ccall(
               (:KSPPIPEGCRSetUnrollW, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, unroll_w,
              )


	return nothing
end 

"""
	unroll_w::PetscBool = KSPPIPEGCRGetUnrollW(petsclib::PetscLibType,ksp::PetscKSP) 
Get information on `KSPPIPEGCR` if it uses unrolling the w vector

Logically Collective

Input Parameter:
- `ksp` - the Krylov space context

Output Parameter:
- `unroll_w` - `KSPPIPEGCR` uses unrolling (bool)

Level: intermediate

-seealso: [](ch_ksp), `KSPPIPEGCR`, `KSPPIPEGCRGetTruncationType()`, `KSPPIPEGCRGetNprealloc()`, `KSPPIPEGCRSetUnrollW()`

# External Links
$(_doc_external("KSP/KSPPIPEGCRGetUnrollW"))
"""
function KSPPIPEGCRGetUnrollW(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPPIPEGCRGetUnrollW(petsclib::$UnionPetscLib, ksp::PetscKSP )
	unroll_w_ = Ref{PetscBool}()

    @chk ccall(
               (:KSPPIPEGCRGetUnrollW, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{PetscBool}),
               ksp, unroll_w_,
              )

	unroll_w = unroll_w_[]

	return unroll_w
end 

"""
	KSPPIPEGCRSetMmax(petsclib::PetscLibType,ksp::PetscKSP, mmax::PetscInt) 
set the maximum number of previous directions `KSPPIPEGCR` will store for orthogonalization

Logically Collective

Input Parameters:
- `ksp`  - the Krylov space context
- `mmax` - the maximum number of previous directions to orthogonalize against

Options Database Key:
- `-ksp_pipegcr_mmax <mmax>` - maximum number of previous directions

Level: intermediate

-seealso: [](ch_ksp), `KSPPIPEGCR`, `KSPPIPEGCRSetTruncationType()`, `KSPPIPEGCRSetNprealloc()`

# External Links
$(_doc_external("KSP/KSPPIPEGCRSetMmax"))
"""
function KSPPIPEGCRSetMmax(petsclib::PetscLibType, ksp::PetscKSP, mmax::PetscInt) end

@for_petsc function KSPPIPEGCRSetMmax(petsclib::$UnionPetscLib, ksp::PetscKSP, mmax::$PetscInt )

    @chk ccall(
               (:KSPPIPEGCRSetMmax, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, mmax,
              )


	return nothing
end 

"""
	mmax::PetscInt = KSPPIPEGCRGetMmax(petsclib::PetscLibType,ksp::PetscKSP) 
get the maximum number of previous directions `KSPPIPEGCR` will store

Not Collective

Input Parameter:
- `ksp` - the Krylov space context

Output Parameter:
- `mmax` - the maximum number of previous directions allowed for orthogonalization

Level: intermediate

-seealso: [](ch_ksp), `KSPPIPEGCR`, `KSPPIPEGCRGetTruncationType()`, `KSPPIPEGCRGetNprealloc()`, `KSPPIPEGCRSetMmax()`

# External Links
$(_doc_external("KSP/KSPPIPEGCRGetMmax"))
"""
function KSPPIPEGCRGetMmax(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPPIPEGCRGetMmax(petsclib::$UnionPetscLib, ksp::PetscKSP )
	mmax_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPPIPEGCRGetMmax, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscInt}),
               ksp, mmax_,
              )

	mmax = mmax_[]

	return mmax
end 

"""
	KSPPIPEGCRSetNprealloc(petsclib::PetscLibType,ksp::PetscKSP, nprealloc::PetscInt) 
set the number of directions to preallocate with `KSPPIPEGCR`

Logically Collective

Input Parameters:
- `ksp`       - the Krylov space context
- `nprealloc` - the number of vectors to preallocate

Level: advanced

Options Database Key:
- `-ksp_pipegcr_nprealloc <N>` - number of vectors to preallocate

-seealso: [](ch_ksp), `KSPPIPEGCR`, `KSPPIPEGCRGetTruncationType()`, `KSPPIPEGCRGetNprealloc()`

# External Links
$(_doc_external("KSP/KSPPIPEGCRSetNprealloc"))
"""
function KSPPIPEGCRSetNprealloc(petsclib::PetscLibType, ksp::PetscKSP, nprealloc::PetscInt) end

@for_petsc function KSPPIPEGCRSetNprealloc(petsclib::$UnionPetscLib, ksp::PetscKSP, nprealloc::$PetscInt )

    @chk ccall(
               (:KSPPIPEGCRSetNprealloc, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, nprealloc,
              )


	return nothing
end 

"""
	nprealloc::PetscInt = KSPPIPEGCRGetNprealloc(petsclib::PetscLibType,ksp::PetscKSP) 
get the number of directions preallocate by `KSPPIPEGCR`

Not Collective

Input Parameter:
- `ksp` - the Krylov space context

Output Parameter:
- `nprealloc` - the number of directions preallocated

Level: advanced

-seealso: [](ch_ksp), `KSPPIPEGCR`, `KSPPIPEGCRGetTruncationType()`, `KSPPIPEGCRSetNprealloc()`

# External Links
$(_doc_external("KSP/KSPPIPEGCRGetNprealloc"))
"""
function KSPPIPEGCRGetNprealloc(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPPIPEGCRGetNprealloc(petsclib::$UnionPetscLib, ksp::PetscKSP )
	nprealloc_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPPIPEGCRGetNprealloc, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscInt}),
               ksp, nprealloc_,
              )

	nprealloc = nprealloc_[]

	return nprealloc
end 

"""
	KSPPIPEGCRSetTruncationType(petsclib::PetscLibType,ksp::PetscKSP, truncstrat::KSPFCDTruncationType) 
specify how many of its stored previous directions `KSPPIPEGCR` uses during orthogonalization

Logically Collective

Input Parameters:
- `ksp`        - the Krylov space context
- `truncstrat` - the choice of strategy
-seealso: [](ch_ksp), `KSPPIPEGCR`, `KSPPIPEGCRTruncationType`, `KSPFCDTruncationType`, `KSP_FCD_TRUNC_TYPE_STANDARD`, `KSP_FCD_TRUNC_TYPE_NOTAY`

# External Links
$(_doc_external("KSP/KSPPIPEGCRSetTruncationType"))
"""
function KSPPIPEGCRSetTruncationType(petsclib::PetscLibType, ksp::PetscKSP, truncstrat::KSPFCDTruncationType) end

@for_petsc function KSPPIPEGCRSetTruncationType(petsclib::$UnionPetscLib, ksp::PetscKSP, truncstrat::KSPFCDTruncationType )

    @chk ccall(
               (:KSPPIPEGCRSetTruncationType, $petsc_library),
               PetscErrorCode,
               (CKSP, KSPFCDTruncationType),
               ksp, truncstrat,
              )


	return nothing
end 

"""
	truncstrat::KSPFCDTruncationType = KSPPIPEGCRGetTruncationType(petsclib::PetscLibType,ksp::PetscKSP) 
get the truncation strategy employed by `KSPPIPEGCR`

Not Collective

Input Parameter:
- `ksp` - the Krylov space context

Output Parameter:
- `truncstrat` - the strategy type
-seealso: [](ch_ksp), `KSPPIPEGCR`, `KSPPIPEGCRSetTruncationType`, `KSPPIPEGCRTruncationType`, `KSPFCDTruncationType`, `KSP_FCD_TRUNC_TYPE_STANDARD`, `KSP_FCD_TRUNC_TYPE_NOTAY`

# External Links
$(_doc_external("KSP/KSPPIPEGCRGetTruncationType"))
"""
function KSPPIPEGCRGetTruncationType(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPPIPEGCRGetTruncationType(petsclib::$UnionPetscLib, ksp::PetscKSP )
	truncstrat_ = Ref{KSPFCDTruncationType}()

    @chk ccall(
               (:KSPPIPEGCRGetTruncationType, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPFCDTruncationType}),
               ksp, truncstrat_,
              )

	truncstrat = unsafe_string(truncstrat_[])

	return truncstrat
end 

"""
	KSPPIPEGCRSetModifyPC(petsclib::PetscLibType,ksp::PetscKSP, fnc::KSPFlexibleModifyPCFn, ctx::Cvoid, destroy::PetscCtxDestroyFn) 
Sets the routine used by `KSPPIPEGCR` to modify the preconditioner at each iteration

Logically Collective

Input Parameters:
- `ksp`      - iterative context obtained from `KSPCreate()`
- `function` - user defined function to modify the preconditioner, see `KSPFlexibleModifyPCFn`
- `ctx`      - user provided context for the modify preconditioner function
- `destroy`  - the function to use to destroy the user provided application context.

Level: intermediate

-seealso: [](ch_ksp), `KSPFlexibleSetModifyPC()`, `KSPFlexibleModifyPCFn`, `KSPPIPEGCR`

# External Links
$(_doc_external("KSP/KSPPIPEGCRSetModifyPC"))
"""
function KSPPIPEGCRSetModifyPC(petsclib::PetscLibType, ksp::PetscKSP, fnc::KSPFlexibleModifyPCFn, ctx::Cvoid, destroy::PetscCtxDestroyFn) end

@for_petsc function KSPPIPEGCRSetModifyPC(petsclib::$UnionPetscLib, ksp::PetscKSP, fnc::KSPFlexibleModifyPCFn, ctx::Cvoid, destroy::PetscCtxDestroyFn )

    @chk ccall(
               (:KSPPIPEGCRSetModifyPC, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPFlexibleModifyPCFn}, Ptr{Cvoid}, Ptr{PetscCtxDestroyFn}),
               ksp, fnc, ctx, destroy,
              )


	return nothing
end 

"""
	KSPMINRESSetUseQLP(petsclib::PetscLibType,ksp::PetscKSP, qlp::PetscBool) 
Use the QLP variant of `KSPMINRES`

Logically Collective

Input Parameters:
- `ksp` - the iterative context
- `qlp` - a Boolean indicating if the QLP variant should be used

Level: beginner

-seealso: [](ch_ksp), `KSP`, `KSPMINRES`, `KSPMINRESGetUseQLP()`

# External Links
$(_doc_external("KSP/KSPMINRESSetUseQLP"))
"""
function KSPMINRESSetUseQLP(petsclib::PetscLibType, ksp::PetscKSP, qlp::PetscBool) end

@for_petsc function KSPMINRESSetUseQLP(petsclib::$UnionPetscLib, ksp::PetscKSP, qlp::PetscBool )

    @chk ccall(
               (:KSPMINRESSetUseQLP, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, qlp,
              )


	return nothing
end 

"""
	KSPMINRESSetRadius(petsclib::PetscLibType,ksp::PetscKSP, radius::PetscReal) 
Set the maximum solution norm allowed for use with trust region methods

Logically Collective

Input Parameters:
- `ksp`    - the iterative context
- `radius` - the value

Level: beginner

Options Database Key:
- `-ksp_minres_radius <real>` - maximum allowed solution norm

-seealso: [](ch_ksp), `KSP`, `KSPMINRES`, `KSPMINRESSetUseQLP()`

# External Links
$(_doc_external("KSP/KSPMINRESSetRadius"))
"""
function KSPMINRESSetRadius(petsclib::PetscLibType, ksp::PetscKSP, radius::PetscReal) end

@for_petsc function KSPMINRESSetRadius(petsclib::$UnionPetscLib, ksp::PetscKSP, radius::$PetscReal )

    @chk ccall(
               (:KSPMINRESSetRadius, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscReal),
               ksp, radius,
              )


	return nothing
end 

"""
	qlp::PetscBool = KSPMINRESGetUseQLP(petsclib::PetscLibType,ksp::PetscKSP) 
Get the flag that indicates if the QLP variant is being used

Logically Collective

Input Parameter:
- `ksp` - the iterative context

Output Parameter:
- `qlp` - a Boolean indicating if the QLP variant is used

Level: beginner

-seealso: [](ch_ksp), `KSP`, `KSPMINRES`, `KSPMINRESSetUseQLP()`

# External Links
$(_doc_external("KSP/KSPMINRESGetUseQLP"))
"""
function KSPMINRESGetUseQLP(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPMINRESGetUseQLP(petsclib::$UnionPetscLib, ksp::PetscKSP )
	qlp_ = Ref{PetscBool}()

    @chk ccall(
               (:KSPMINRESGetUseQLP, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{PetscBool}),
               ksp, qlp_,
              )

	qlp = qlp_[]

	return qlp
end 

"""
	KSPLSQRSetComputeStandardErrorVec(petsclib::PetscLibType,ksp::PetscKSP, flg::PetscBool) 
Compute a vector of standard error estimates during `KSPSolve()` for  `KSPLSQR`.

Logically Collective

Input Parameters:
- `ksp` - iterative context
- `flg` - compute the vector of standard estimates or not

Level: intermediate

-seealso: [](ch_ksp), `KSPSolve()`, `KSPLSQR`, `KSPLSQRGetStandardErrorVec()`

# External Links
$(_doc_external("KSP/KSPLSQRSetComputeStandardErrorVec"))
"""
function KSPLSQRSetComputeStandardErrorVec(petsclib::PetscLibType, ksp::PetscKSP, flg::PetscBool) end

@for_petsc function KSPLSQRSetComputeStandardErrorVec(petsclib::$UnionPetscLib, ksp::PetscKSP, flg::PetscBool )

    @chk ccall(
               (:KSPLSQRSetComputeStandardErrorVec, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flg,
              )


	return nothing
end 

"""
	KSPLSQRSetExactMatNorm(petsclib::PetscLibType,ksp::PetscKSP, flg::PetscBool) 
Compute exact matrix norm instead of iteratively refined estimate.

Not Collective

Input Parameters:
- `ksp` - iterative context
- `flg` - compute exact matrix norm or not

Level: intermediate

-seealso: [](ch_ksp), `KSPSolve()`, `KSPLSQR`, `KSPLSQRGetNorms()`, `KSPLSQRConvergedDefault()`

# External Links
$(_doc_external("KSP/KSPLSQRSetExactMatNorm"))
"""
function KSPLSQRSetExactMatNorm(petsclib::PetscLibType, ksp::PetscKSP, flg::PetscBool) end

@for_petsc function KSPLSQRSetExactMatNorm(petsclib::$UnionPetscLib, ksp::PetscKSP, flg::PetscBool )

    @chk ccall(
               (:KSPLSQRSetExactMatNorm, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flg,
              )


	return nothing
end 

"""
	KSPLSQRGetStandardErrorVec(petsclib::PetscLibType,ksp::PetscKSP, se::PetscVec) 
Get vector of standard error estimates.
Only available if -ksp_lsqr_set_standard_error was set to true
or `KSPLSQRSetComputeStandardErrorVec`(ksp, `PETSC_TRUE`) was called.
Otherwise returns `NULL`.

Not Collective

Input Parameter:
- `ksp` - iterative context

Output Parameter:
- `se` - vector of standard estimates

Level: intermediate

-seealso: [](ch_ksp), `KSPSolve()`, `KSPLSQR`, `KSPLSQRSetComputeStandardErrorVec()`

# External Links
$(_doc_external("KSP/KSPLSQRGetStandardErrorVec"))
"""
function KSPLSQRGetStandardErrorVec(petsclib::PetscLibType, ksp::PetscKSP, se::PetscVec) end

@for_petsc function KSPLSQRGetStandardErrorVec(petsclib::$UnionPetscLib, ksp::PetscKSP, se::PetscVec )
	se_ = Ref(se.ptr)

    @chk ccall(
               (:KSPLSQRGetStandardErrorVec, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{CVec}),
               ksp, se_,
              )

	se.ptr = C_NULL

	return nothing
end 

"""
	arnorm::PetscReal,anorm::PetscReal = KSPLSQRGetNorms(petsclib::PetscLibType,ksp::PetscKSP) 
Get the norm estimates that `KSPLSQR` computes internally during `KSPSolve()`.

Not Collective

Input Parameter:
- `ksp` - iterative context

Output Parameters:
- `arnorm` - good estimate of (A*Pmat^{-T})*r, where r = A x - b, used in specific stopping criterion
- `anorm`  - poor estimate of A*Pmat^{-T}_{frobenius} used in specific stopping criterion

Level: intermediate

-seealso: [](ch_ksp), `KSPSolve()`, `KSPLSQR`, `KSPLSQRSetExactMatNorm()`

# External Links
$(_doc_external("KSP/KSPLSQRGetNorms"))
"""
function KSPLSQRGetNorms(petsclib::PetscLibType, ksp::PetscKSP) end

@for_petsc function KSPLSQRGetNorms(petsclib::$UnionPetscLib, ksp::PetscKSP )
	arnorm_ = Ref{$PetscReal}()
	anorm_ = Ref{$PetscReal}()

    @chk ccall(
               (:KSPLSQRGetNorms, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscReal}, Ptr{$PetscReal}),
               ksp, arnorm_, anorm_,
              )

	arnorm = arnorm_[]
	anorm = anorm_[]

	return arnorm,anorm
end 

"""
	KSPLSQRMonitorResidual(petsclib::PetscLibType,ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) 
Prints the residual norm, as well as the normal equation residual norm, at each iteration of an iterative solver for the `KSPLSQR` solver

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_lsqr_monitor` - Activates `KSPLSQRMonitorResidual()`

Level: intermediate

-seealso: [](ch_ksp), `KSPLSQR`, `KSPMonitorSet()`, `KSPMonitorResidual()`, `KSPMonitorTrueResidualMaxNorm()`, `KSPLSQRMonitorResidualDrawLG()`

# External Links
$(_doc_external("KSP/KSPLSQRMonitorResidual"))
"""
function KSPLSQRMonitorResidual(petsclib::PetscLibType, ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function KSPLSQRMonitorResidual(petsclib::$UnionPetscLib, ksp::PetscKSP, n::$PetscInt, rnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:KSPLSQRMonitorResidual, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	KSPLSQRMonitorResidualDrawLG(petsclib::PetscLibType,ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) 
Plots the true residual norm at each iteration of an iterative solver for the `KSPLSQR` solver

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_lsqr_monitor draw::draw_lg` - Activates `KSPMonitorTrueResidualDrawLG()`

Level: intermediate

-seealso: [](ch_ksp), `KSPLSQR`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPLSQRMonitorResidual()`, `KSPLSQRMonitorResidualDrawLGCreate()`

# External Links
$(_doc_external("KSP/KSPLSQRMonitorResidualDrawLG"))
"""
function KSPLSQRMonitorResidualDrawLG(petsclib::PetscLibType, ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function KSPLSQRMonitorResidualDrawLG(petsclib::$UnionPetscLib, ksp::PetscKSP, n::$PetscInt, rnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:KSPLSQRMonitorResidualDrawLG, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	vf::PetscViewerAndFormat = KSPLSQRMonitorResidualDrawLGCreate(petsclib::PetscLibType,viewer::PetscViewer, format::PetscViewerFormat, ctx::Cvoid) 
Creates the line graph object for the `KSPLSQR` residual and normal equation residual norm

Collective

Input Parameters:
- `viewer` - The `PetscViewer`
- `format` - The viewer format
- `ctx`    - An optional user context

Output Parameter:
- `vf` - The `PetscViewerAndFormat`

Level: intermediate

-seealso: [](ch_ksp), `KSPLSQR`, `KSPMonitorSet()`, `KSPLSQRMonitorResidual()`, `KSPLSQRMonitorResidualDrawLG()`

# External Links
$(_doc_external("KSP/KSPLSQRMonitorResidualDrawLGCreate"))
"""
function KSPLSQRMonitorResidualDrawLGCreate(petsclib::PetscLibType, viewer::PetscViewer, format::PetscViewerFormat, ctx::Cvoid) end

@for_petsc function KSPLSQRMonitorResidualDrawLGCreate(petsclib::$UnionPetscLib, viewer::PetscViewer, format::PetscViewerFormat, ctx::Cvoid )
	vf_ = Ref{PetscViewerAndFormat}()

    @chk ccall(
               (:KSPLSQRMonitorResidualDrawLGCreate, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscViewerFormat, Ptr{Cvoid}, PetscViewerAndFormat),
               viewer, format, ctx, vf_,
              )

	vf = vf_[]

	return vf
end 

"""
	KSPLSQRConvergedDefault(petsclib::PetscLibType,ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, reason::KSPConvergedReason, ctx::Cvoid) 
Determines convergence of the `KSPLSQR` Krylov method, including a check on the residual norm of the normal equations.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm residual value (may be estimated)
- `ctx`   - convergence context which must have been created by `KSPConvergedDefaultCreate()`

Output Parameter:
- `reason` - the convergence reason

Level: advanced

-seealso: [](ch_ksp), `KSPLSQR`, `KSPSetConvergenceTest()`, `KSPSetTolerances()`, `KSPConvergedSkip()`, `KSPConvergedReason`, `KSPGetConvergedReason()`,
`KSPConvergedDefaultSetUIRNorm()`, `KSPConvergedDefaultSetUMIRNorm()`, `KSPConvergedDefaultCreate()`, `KSPConvergedDefaultDestroy()`,
`KSPConvergedDefault()`, `KSPLSQRGetNorms()`, `KSPLSQRSetExactMatNorm()`

# External Links
$(_doc_external("KSP/KSPLSQRConvergedDefault"))
"""
function KSPLSQRConvergedDefault(petsclib::PetscLibType, ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, reason::KSPConvergedReason, ctx::Cvoid) end

@for_petsc function KSPLSQRConvergedDefault(petsclib::$UnionPetscLib, ksp::PetscKSP, n::$PetscInt, rnorm::$PetscReal, reason::KSPConvergedReason, ctx::Cvoid )

    @chk ccall(
               (:KSPLSQRConvergedDefault, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{KSPConvergedReason}, Ptr{Cvoid}),
               ksp, n, rnorm, reason, ctx,
              )


	return nothing
end 

"""
	KSPMonitorSNESResidual(petsclib::PetscLibType,ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) 
Prints the `SNES` residual norm, as well as the `KSP` residual norm, at each iteration of a `KSPSolve()` called within a `SNESSolve()`.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-snes_monitor_ksp` - Activates `KSPMonitorSNESResidual()`

Level: intermediate

-seealso: [](ch_snes), `SNES`, `KSPMonitorSet()`, `KSPMonitorResidual()`, `KSPMonitorTrueResidualMaxNorm()`, `KSPMonitor()`, `SNESMonitor()`, `PetscViewerAndFormat()`

# External Links
$(_doc_external("Snes/KSPMonitorSNESResidual"))
"""
function KSPMonitorSNESResidual(petsclib::PetscLibType, ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function KSPMonitorSNESResidual(petsclib::$UnionPetscLib, ksp::PetscKSP, n::$PetscInt, rnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:KSPMonitorSNESResidual, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	KSPMonitorSNESResidualDrawLG(petsclib::PetscLibType,ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) 
Plots the linear `KSP` residual norm and the `SNES` residual norm of a `KSPSolve()` called within a `SNESSolve()`.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context, created with `KSPMonitorSNESResidualDrawLGCreate()`

Options Database Key:
- `-snes_monitor_ksp draw::draw_lg` - Activates `KSPMonitorSNESResidualDrawLG()`

Level: intermediate

-seealso: [](ch_snes), `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `SNESMonitor()`, `KSPMonitor()`, `KSPMonitorSNESResidualDrawLGCreate()`

# External Links
$(_doc_external("Snes/KSPMonitorSNESResidualDrawLG"))
"""
function KSPMonitorSNESResidualDrawLG(petsclib::PetscLibType, ksp::PetscKSP, n::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function KSPMonitorSNESResidualDrawLG(petsclib::$UnionPetscLib, ksp::PetscKSP, n::$PetscInt, rnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:KSPMonitorSNESResidualDrawLG, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	vf::PetscViewerAndFormat = KSPMonitorSNESResidualDrawLGCreate(petsclib::PetscLibType,viewer::PetscViewer, format::PetscViewerFormat, ctx::Cvoid) 
Creates the `PetscViewer` used by `KSPMonitorSNESResidualDrawLG()`

Collective

Input Parameters:
- `viewer` - The `PetscViewer`
- `format` - The viewer format
- `ctx`    - An optional user context

Output Parameter:
- `vf` - The viewer context

Level: intermediate

-seealso: [](ch_snes), `KSP`, `SNES`, `PetscViewerFormat`, `PetscViewerAndFormat`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`

# External Links
$(_doc_external("Snes/KSPMonitorSNESResidualDrawLGCreate"))
"""
function KSPMonitorSNESResidualDrawLGCreate(petsclib::PetscLibType, viewer::PetscViewer, format::PetscViewerFormat, ctx::Cvoid) end

@for_petsc function KSPMonitorSNESResidualDrawLGCreate(petsclib::$UnionPetscLib, viewer::PetscViewer, format::PetscViewerFormat, ctx::Cvoid )
	vf_ = Ref{PetscViewerAndFormat}()

    @chk ccall(
               (:KSPMonitorSNESResidualDrawLGCreate, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscViewerFormat, Ptr{Cvoid}, PetscViewerAndFormat),
               viewer, format, ctx, vf_,
              )

	vf = vf_[]

	return vf
end 

