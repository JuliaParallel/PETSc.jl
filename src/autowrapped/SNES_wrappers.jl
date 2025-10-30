# autodefined type arguments for class ------
mutable struct SNESFunctionFn end

mutable struct SNESNGSFn end

mutable struct SNESJacobianFn end

mutable struct SNESInitialGuessFn end

mutable struct SNESUpdateFn end

mutable struct _n_SNESLineSearch end
const SNESLineSearch = Ptr{_n_SNESLineSearch}

mutable struct SNESObjectiveFn end

# -------------------------------------------------------
"""
	SNESMonitorSolution(petsclib::PetscLibType,snes::PetscSNES, its::PetscInt, fgnorm::PetscReal, vf::PetscViewerAndFormat) 
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

-seealso: [](ch_snes), `SNES`, `SNESMonitorSet()`, `SNESMonitorDefault()`, `VecView()`

# External Links
$(_doc_external("Snes/SNESMonitorSolution"))
"""
function SNESMonitorSolution(petsclib::PetscLibType, snes::PetscSNES, its::PetscInt, fgnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function SNESMonitorSolution(petsclib::$UnionPetscLib, snes::PetscSNES, its::$PetscInt, fgnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:SNESMonitorSolution, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               snes, its, fgnorm, vf,
              )


	return nothing
end 

"""
	SNESMonitorResidual(petsclib::PetscLibType,snes::PetscSNES, its::PetscInt, fgnorm::PetscReal, vf::PetscViewerAndFormat) 
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

-seealso: [](ch_snes), `SNES`, `SNESMonitorSet()`, `SNESMonitorDefault()`, `VecView()`, `SNESMonitor()`

# External Links
$(_doc_external("Snes/SNESMonitorResidual"))
"""
function SNESMonitorResidual(petsclib::PetscLibType, snes::PetscSNES, its::PetscInt, fgnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function SNESMonitorResidual(petsclib::$UnionPetscLib, snes::PetscSNES, its::$PetscInt, fgnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:SNESMonitorResidual, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               snes, its, fgnorm, vf,
              )


	return nothing
end 

"""
	SNESMonitorSolutionUpdate(petsclib::PetscLibType,snes::PetscSNES, its::PetscInt, fgnorm::PetscReal, vf::PetscViewerAndFormat) 
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

-seealso: [](ch_snes), `SNESMonitorSet()`, `SNESMonitorDefault()`, `VecView()`, `SNESMonitor()`

# External Links
$(_doc_external("Snes/SNESMonitorSolutionUpdate"))
"""
function SNESMonitorSolutionUpdate(petsclib::PetscLibType, snes::PetscSNES, its::PetscInt, fgnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function SNESMonitorSolutionUpdate(petsclib::$UnionPetscLib, snes::PetscSNES, its::$PetscInt, fgnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:SNESMonitorSolutionUpdate, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               snes, its, fgnorm, vf,
              )


	return nothing
end 

"""
	SNESMonitorDefaultSetUp(petsclib::PetscLibType,snes::PetscSNES, vf::PetscViewerAndFormat) 

# External Links
$(_doc_external("Snes/SNESMonitorDefaultSetUp"))
"""
function SNESMonitorDefaultSetUp(petsclib::PetscLibType, snes::PetscSNES, vf::PetscViewerAndFormat) end

@for_petsc function SNESMonitorDefaultSetUp(petsclib::$UnionPetscLib, snes::PetscSNES, vf::PetscViewerAndFormat )

    @chk ccall(
               (:SNESMonitorDefaultSetUp, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{PetscViewerAndFormat}),
               snes, vf,
              )


	return nothing
end 

"""
	SNESMonitorDefault(petsclib::PetscLibType,snes::PetscSNES, its::PetscInt, fgnorm::PetscReal, vf::PetscViewerAndFormat) 
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

-seealso: [](ch_snes), `SNESMonitorSet()`, `SNESMonitorSolution()`, `SNESMonitorFunction()`, `SNESMonitorResidual()`,
`SNESMonitorSolutionUpdate()`, `SNESMonitorScaling()`, `SNESMonitorRange()`, `SNESMonitorRatio()`,
`SNESMonitorDefaultField()`, `PetscViewerFormat`, `PetscViewerAndFormat`

# External Links
$(_doc_external("Snes/SNESMonitorDefault"))
"""
function SNESMonitorDefault(petsclib::PetscLibType, snes::PetscSNES, its::PetscInt, fgnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function SNESMonitorDefault(petsclib::$UnionPetscLib, snes::PetscSNES, its::$PetscInt, fgnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:SNESMonitorDefault, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               snes, its, fgnorm, vf,
              )


	return nothing
end 

"""
	SNESMonitorScaling(petsclib::PetscLibType,snes::PetscSNES, its::PetscInt, fgnorm::PetscReal, vf::PetscViewerAndFormat) 
Monitors the largest value in each row of the Jacobian of a `SNESSolve()`

Collective

Input Parameters:
- `snes`   - the `SNES` context
- `its`    - iteration number
- `fgnorm` - 2-norm of residual
- `vf`     - viewer and format structure

Level: intermediate

-seealso: [](ch_snes), `SNESMonitorSet()`, `SNESMonitorSolution()`, `SNESMonitorRange()`, `SNESMonitorJacUpdateSpectrum()`,
`PetscViewerFormat`, `PetscViewerAndFormat`

# External Links
$(_doc_external("Snes/SNESMonitorScaling"))
"""
function SNESMonitorScaling(petsclib::PetscLibType, snes::PetscSNES, its::PetscInt, fgnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function SNESMonitorScaling(petsclib::$UnionPetscLib, snes::PetscSNES, its::$PetscInt, fgnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:SNESMonitorScaling, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               snes, its, fgnorm, vf,
              )


	return nothing
end 

"""
	SNESMonitorJacUpdateSpectrum(petsclib::PetscLibType,snes::PetscSNES, it::PetscInt, fnorm::PetscReal, vf::PetscViewerAndFormat) 
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

-seealso: [](ch_snes), `SNESMonitorSet()`, `SNESMonitorSolution()`, `SNESMonitorRange()`, `PetscViewerFormat`, `PetscViewerAndFormat`

# External Links
$(_doc_external("Snes/SNESMonitorJacUpdateSpectrum"))
"""
function SNESMonitorJacUpdateSpectrum(petsclib::PetscLibType, snes::PetscSNES, it::PetscInt, fnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function SNESMonitorJacUpdateSpectrum(petsclib::$UnionPetscLib, snes::PetscSNES, it::$PetscInt, fnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:SNESMonitorJacUpdateSpectrum, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               snes, it, fnorm, vf,
              )


	return nothing
end 

"""
	SNESMonitorRange(petsclib::PetscLibType,snes::PetscSNES, it::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) 
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

-seealso: [](ch_snes), `SNESMonitorSet()`, `SNESMonitorDefault()`, `SNESMonitorLGCreate()`, `SNESMonitorScaling()`, `PetscViewerFormat`, `PetscViewerAndFormat`

# External Links
$(_doc_external("Snes/SNESMonitorRange"))
"""
function SNESMonitorRange(petsclib::PetscLibType, snes::PetscSNES, it::PetscInt, rnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function SNESMonitorRange(petsclib::$UnionPetscLib, snes::PetscSNES, it::$PetscInt, rnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:SNESMonitorRange, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               snes, it, rnorm, vf,
              )


	return nothing
end 

"""
	SNESMonitorRatio(petsclib::PetscLibType,snes::PetscSNES, its::PetscInt, fgnorm::PetscReal, vf::PetscViewerAndFormat) 
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

-seealso: [](ch_snes), `SNESMonitorRationSetUp()`, `SNESMonitorSet()`, `SNESMonitorSolution()`, `SNESMonitorDefault()`, `PetscViewerFormat`, `PetscViewerAndFormat`

# External Links
$(_doc_external("Snes/SNESMonitorRatio"))
"""
function SNESMonitorRatio(petsclib::PetscLibType, snes::PetscSNES, its::PetscInt, fgnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function SNESMonitorRatio(petsclib::$UnionPetscLib, snes::PetscSNES, its::$PetscInt, fgnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:SNESMonitorRatio, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               snes, its, fgnorm, vf,
              )


	return nothing
end 

"""
	SNESMonitorRatioSetUp(petsclib::PetscLibType,snes::PetscSNES, vf::PetscViewerAndFormat) 
Insures the `SNES` object is saving its history since this monitor needs access to it

Collective

Input Parameters:
- `snes` - the `SNES` context
- `vf`   - `PetscViewerAndFormat` (ignored)

Level: intermediate

-seealso: [](ch_snes), `SNESMonitorSet()`, `SNESMonitorSolution()`, `SNESMonitorDefault()`, `SNESMonitorRatio()`, `PetscViewerFormat`, `PetscViewerAndFormat`

# External Links
$(_doc_external("Snes/SNESMonitorRatioSetUp"))
"""
function SNESMonitorRatioSetUp(petsclib::PetscLibType, snes::PetscSNES, vf::PetscViewerAndFormat) end

@for_petsc function SNESMonitorRatioSetUp(petsclib::$UnionPetscLib, snes::PetscSNES, vf::PetscViewerAndFormat )

    @chk ccall(
               (:SNESMonitorRatioSetUp, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{PetscViewerAndFormat}),
               snes, vf,
              )


	return nothing
end 

"""
	SNESMonitorDefaultShort(petsclib::PetscLibType,snes::PetscSNES, its::PetscInt, fgnorm::PetscReal, vf::PetscViewerAndFormat) 

# External Links
$(_doc_external("Snes/SNESMonitorDefaultShort"))
"""
function SNESMonitorDefaultShort(petsclib::PetscLibType, snes::PetscSNES, its::PetscInt, fgnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function SNESMonitorDefaultShort(petsclib::$UnionPetscLib, snes::PetscSNES, its::$PetscInt, fgnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:SNESMonitorDefaultShort, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               snes, its, fgnorm, vf,
              )


	return nothing
end 

"""
	SNESMonitorDefaultField(petsclib::PetscLibType,snes::PetscSNES, its::PetscInt, fgnorm::PetscReal, vf::PetscViewerAndFormat) 
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

-seealso: [](ch_snes), `SNESMonitorSet()`, `SNESMonitorSolution()`, `SNESMonitorDefault()`, `PetscViewerFormat`, `PetscViewerAndFormat`

# External Links
$(_doc_external("Snes/SNESMonitorDefaultField"))
"""
function SNESMonitorDefaultField(petsclib::PetscLibType, snes::PetscSNES, its::PetscInt, fgnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function SNESMonitorDefaultField(petsclib::$UnionPetscLib, snes::PetscSNES, its::$PetscInt, fgnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:SNESMonitorDefaultField, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               snes, its, fgnorm, vf,
              )


	return nothing
end 

"""
	SNESConvergedDefault(petsclib::PetscLibType,snes::PetscSNES, it::PetscInt, xnorm::PetscReal, snorm::PetscReal, fnorm::PetscReal, reason::SNESConvergedReason, dummy::Cvoid) 
Default convergence test for `SNESSolve()`.

Collective

Input Parameters:
- `snes`  - the `SNES` context
- `it`    - the iteration (0 indicates before any Newton steps)
- `xnorm` - 2-norm of current iterate
- `snorm` - 2-norm of current step
- `fnorm` - 2-norm of function at current iterate
- `dummy` - unused context

Output Parameter:
- `reason` - converged reason, see `SNESConvergedReason`

Options Database Keys:
- `-snes_convergence_test default`      - see `SNESSetFromOptions()`
- `-snes_stol`                          - convergence tolerance in terms of the norm of the change in the solution between steps
- `-snes_atol <abstol>`                 - absolute tolerance of residual norm
- `-snes_rtol <rtol>`                   - relative decrease in tolerance norm from the initial 2-norm of the solution
- `-snes_divergence_tolerance <divtol>` - if the residual goes above divtol*rnorm0, exit with divergence
- `-snes_max_funcs <max_funcs>`         - maximum number of function evaluations, use `unlimited` for no maximum
- `-snes_max_fail <max_fail>`           - maximum number of line search failures allowed before stopping, default is none
- `-snes_max_linear_solve_fail`         - number of linear solver failures before `SNESSolve()` stops

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESSolve()`, `SNESSetConvergenceTest()`, `SNESConvergedSkip()`, `SNESSetTolerances()`, `SNESSetDivergenceTolerance()`,
`SNESConvergedReason`

# External Links
$(_doc_external("Snes/SNESConvergedDefault"))
"""
function SNESConvergedDefault(petsclib::PetscLibType, snes::PetscSNES, it::PetscInt, xnorm::PetscReal, snorm::PetscReal, fnorm::PetscReal, reason::SNESConvergedReason, dummy::Cvoid) end

@for_petsc function SNESConvergedDefault(petsclib::$UnionPetscLib, snes::PetscSNES, it::$PetscInt, xnorm::$PetscReal, snorm::$PetscReal, fnorm::$PetscReal, reason::SNESConvergedReason, dummy::Cvoid )

    @chk ccall(
               (:SNESConvergedDefault, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, $PetscReal, $PetscReal, Ptr{SNESConvergedReason}, Ptr{Cvoid}),
               snes, it, xnorm, snorm, fnorm, reason, dummy,
              )


	return nothing
end 

"""
	SNESConvergedSkip(petsclib::PetscLibType,snes::PetscSNES, it::PetscInt, xnorm::PetscReal, snorm::PetscReal, fnorm::PetscReal, reason::SNESConvergedReason, dummy::Cvoid) 
Convergence test for `SNES` that NEVER returns as
converged, UNLESS the maximum number of iteration have been reached.

Logically Collective

Input Parameters:
- `snes`  - the `SNES` context
- `it`    - the iteration (0 indicates before any Newton steps)
- `xnorm` - 2-norm of current iterate
- `snorm` - 2-norm of current step
- `fnorm` - 2-norm of function at current iterate
- `dummy` - unused context

Output Parameter:
- `reason` - `SNES_CONVERGED_ITERATING`, `SNES_CONVERGED_ITS`, or `SNES_DIVERGED_FNORM_NAN`

Options Database Key:
- `-snes_convergence_test skip` - see `SNESSetFromOptions()`

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESSolve()`, `SNESConvergedDefault()`, `SNESSetConvergenceTest()`, `SNESConvergedReason`

# External Links
$(_doc_external("Snes/SNESConvergedSkip"))
"""
function SNESConvergedSkip(petsclib::PetscLibType, snes::PetscSNES, it::PetscInt, xnorm::PetscReal, snorm::PetscReal, fnorm::PetscReal, reason::SNESConvergedReason, dummy::Cvoid) end

@for_petsc function SNESConvergedSkip(petsclib::$UnionPetscLib, snes::PetscSNES, it::$PetscInt, xnorm::$PetscReal, snorm::$PetscReal, fnorm::$PetscReal, reason::SNESConvergedReason, dummy::Cvoid )

    @chk ccall(
               (:SNESConvergedSkip, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, $PetscReal, $PetscReal, Ptr{SNESConvergedReason}, Ptr{Cvoid}),
               snes, it, xnorm, snorm, fnorm, reason, dummy,
              )


	return nothing
end 

"""
	SNESSetWorkVecs(petsclib::PetscLibType,snes::PetscSNES, nw::PetscInt) 
Allocates a number of work vectors to be used internally by the `SNES` solver

Input Parameters:
- `snes` - the `SNES` context
- `nw`   - number of work vectors to allocate

Level: developer

-seealso: [](ch_snes), `SNES`

# External Links
$(_doc_external("Snes/SNESSetWorkVecs"))
"""
function SNESSetWorkVecs(petsclib::PetscLibType, snes::PetscSNES, nw::PetscInt) end

@for_petsc function SNESSetWorkVecs(petsclib::$UnionPetscLib, snes::PetscSNES, nw::$PetscInt )

    @chk ccall(
               (:SNESSetWorkVecs, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt),
               snes, nw,
              )


	return nothing
end 

"""
	SNESFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the PETSc interface to the `SNES` package. It is
called from `PetscFinalize()`.

Level: developer

-seealso: [](ch_snes), `SNES`, `PetscFinalize()`

# External Links
$(_doc_external("Snes/SNESFinalizePackage"))
"""
function SNESFinalizePackage(petsclib::PetscLibType) end

@for_petsc function SNESFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:SNESFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	SNESInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `SNES` package. It is called
from PetscDLLibraryRegister_petscsnes() when using dynamic libraries, and on the first call to `SNESCreate()`
when using shared or static libraries.

Level: developer

-seealso: [](ch_snes), `SNES`, `PetscInitialize()`

# External Links
$(_doc_external("Snes/SNESInitializePackage"))
"""
function SNESInitializePackage(petsclib::PetscLibType) end

@for_petsc function SNESInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:SNESInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	SNESSetErrorIfNotConverged(petsclib::PetscLibType,snes::PetscSNES, flg::PetscBool) 
Causes `SNESSolve()` to generate an error immediately if the solver has not converged.

Logically Collective

Input Parameters:
- `snes` - iterative context obtained from `SNESCreate()`
- `flg`  - `PETSC_TRUE` indicates you want the error generated

Options Database Key:
- `-snes_error_if_not_converged <true,false>` - cause an immediate error condition and stop the program if the solver does not converge

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESGetErrorIfNotConverged()`, `KSPGetErrorIfNotConverged()`, `KSPSetErrorIfNotConverged()`

# External Links
$(_doc_external("Snes/SNESSetErrorIfNotConverged"))
"""
function SNESSetErrorIfNotConverged(petsclib::PetscLibType, snes::PetscSNES, flg::PetscBool) end

@for_petsc function SNESSetErrorIfNotConverged(petsclib::$UnionPetscLib, snes::PetscSNES, flg::PetscBool )

    @chk ccall(
               (:SNESSetErrorIfNotConverged, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, flg,
              )


	return nothing
end 

"""
	flag::PetscBool = SNESGetErrorIfNotConverged(petsclib::PetscLibType,snes::PetscSNES) 
Indicates if `SNESSolve()` will generate an error if the solver does not converge?

Not Collective

Input Parameter:
- `snes` - iterative context obtained from `SNESCreate()`

Output Parameter:
- `flag` - `PETSC_TRUE` if it will generate an error, else `PETSC_FALSE`

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESSolve()`, `SNESSetErrorIfNotConverged()`, `KSPGetErrorIfNotConverged()`, `KSPSetErrorIfNotConverged()`

# External Links
$(_doc_external("Snes/SNESGetErrorIfNotConverged"))
"""
function SNESGetErrorIfNotConverged(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESGetErrorIfNotConverged(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESSetAlwaysComputesFinalResidual(petsclib::PetscLibType,snes::PetscSNES, flg::PetscBool) 
tells the `SNES` to always compute the residual (nonlinear function value) at the final solution

Logically Collective

Input Parameters:
- `snes` - the shell `SNES`
- `flg`  - `PETSC_TRUE` to always compute the residual

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESSolve()`, `SNESGetAlwaysComputesFinalResidual()`

# External Links
$(_doc_external("Snes/SNESSetAlwaysComputesFinalResidual"))
"""
function SNESSetAlwaysComputesFinalResidual(petsclib::PetscLibType, snes::PetscSNES, flg::PetscBool) end

@for_petsc function SNESSetAlwaysComputesFinalResidual(petsclib::$UnionPetscLib, snes::PetscSNES, flg::PetscBool )

    @chk ccall(
               (:SNESSetAlwaysComputesFinalResidual, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, flg,
              )


	return nothing
end 

"""
	flg::PetscBool = SNESGetAlwaysComputesFinalResidual(petsclib::PetscLibType,snes::PetscSNES) 
checks if the `SNES` always computes the residual at the final solution

Logically Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `flg` - `PETSC_TRUE` if the residual is computed

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESSolve()`, `SNESSetAlwaysComputesFinalResidual()`

# External Links
$(_doc_external("Snes/SNESGetAlwaysComputesFinalResidual"))
"""
function SNESGetAlwaysComputesFinalResidual(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESGetAlwaysComputesFinalResidual(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESSetFunctionDomainError(petsclib::PetscLibType,snes::PetscSNES) 
tells `SNES` that the input vector, a proposed new solution, to your function you provided to `SNESSetFunction()` is not
in the functions domain. For example, a step with negative pressure.

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Level: advanced

-seealso: [](ch_snes), `SNESCreate()`, `SNESSetFunction()`, `SNESFunctionFn`, `SNESSetJacobianDomainError()`, `SNESVISetVariableBounds()`,
`SNESVISetComputeVariableBounds()`, `SNESLineSearchSetPreCheck()`, `SNESLineSearchSetPostCheck()`, `SNESConvergedReason`, `SNESGetConvergedReason()`,
`SNES_DIVERGED_FUNCTION_DOMAIN`

# External Links
$(_doc_external("Snes/SNESSetFunctionDomainError"))
"""
function SNESSetFunctionDomainError(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESSetFunctionDomainError(petsclib::$UnionPetscLib, snes::PetscSNES )

    @chk ccall(
               (:SNESSetFunctionDomainError, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	SNESSetJacobianDomainError(petsclib::PetscLibType,snes::PetscSNES) 
tells `SNES` that the function you provided to `SNESSetJacobian()` at the proposed step. For example there is a negative element transformation.

Logically Collective

Input Parameter:
- `snes` - the `SNES` context

Level: advanced

-seealso: [](ch_snes), `SNESCreate()`, `SNESSetFunction()`, `SNESFunctionFn`, `SNESSetFunctionDomainError()`, `SNESVISetVariableBounds()`,
`SNESVISetComputeVariableBounds()`, `SNESLineSearchSetPreCheck()`, `SNESLineSearchSetPostCheck()`, `SNESConvergedReason`, `SNESGetConvergedReason()`

# External Links
$(_doc_external("Snes/SNESSetJacobianDomainError"))
"""
function SNESSetJacobianDomainError(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESSetJacobianDomainError(petsclib::$UnionPetscLib, snes::PetscSNES )

    @chk ccall(
               (:SNESSetJacobianDomainError, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	SNESSetCheckJacobianDomainError(petsclib::PetscLibType,snes::PetscSNES, flg::PetscBool) 
tells `SNESSolve()` whether to check if the user called `SNESSetJacobianDomainError()` Jacobian domain error after
each Jacobian evaluation. By default, it checks for the Jacobian domain error in the debug mode, and does not check it in the optimized mode.

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `flg`  - indicates if or not to check Jacobian domain error after each Jacobian evaluation

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESConvergedReason`, `SNESCreate()`, `SNESSetFunction()`, `SNESFunctionFn`, `SNESSetFunctionDomainError()`, `SNESGetCheckJacobianDomainError()`

# External Links
$(_doc_external("Snes/SNESSetCheckJacobianDomainError"))
"""
function SNESSetCheckJacobianDomainError(petsclib::PetscLibType, snes::PetscSNES, flg::PetscBool) end

@for_petsc function SNESSetCheckJacobianDomainError(petsclib::$UnionPetscLib, snes::PetscSNES, flg::PetscBool )

    @chk ccall(
               (:SNESSetCheckJacobianDomainError, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, flg,
              )


	return nothing
end 

"""
	flg::PetscBool = SNESGetCheckJacobianDomainError(petsclib::PetscLibType,snes::PetscSNES) 
Get an indicator whether or not `SNES` is checking Jacobian domain errors after each Jacobian evaluation.

Logically Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `flg` - `PETSC_FALSE` indicates that it is not checking Jacobian domain errors after each Jacobian evaluation

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESCreate()`, `SNESSetFunction()`, `SNESFunctionFn`, `SNESSetFunctionDomainError()`, `SNESSetCheckJacobianDomainError()`

# External Links
$(_doc_external("Snes/SNESGetCheckJacobianDomainError"))
"""
function SNESGetCheckJacobianDomainError(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESGetCheckJacobianDomainError(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	domainerror::PetscBool = SNESGetFunctionDomainError(petsclib::PetscLibType,snes::PetscSNES) 
Gets the status of the domain error after a call to `SNESComputeFunction()`

Not Collective, different MPI processes may return different values

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `domainerror` - Set to `PETSC_TRUE` if there's a domain error; `PETSC_FALSE` otherwise.

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESSetFunctionDomainError()`, `SNESComputeFunction()`

# External Links
$(_doc_external("Snes/SNESGetFunctionDomainError"))
"""
function SNESGetFunctionDomainError(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESGetFunctionDomainError(petsclib::$UnionPetscLib, snes::PetscSNES )
	domainerror_ = Ref{PetscBool}()

    @chk ccall(
               (:SNESGetFunctionDomainError, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{PetscBool}),
               snes, domainerror_,
              )

	domainerror = domainerror_[]

	return domainerror
end 

"""
	domainerror::PetscBool = SNESGetJacobianDomainError(petsclib::PetscLibType,snes::PetscSNES) 
Gets the status of the Jacobian domain error after a call to `SNESComputeJacobian()`

Not Collective, different MPI processes may return different values

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `domainerror` - Set to `PETSC_TRUE` if there's a Jacobian domain error; `PETSC_FALSE` otherwise.

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESSetFunctionDomainError()`, `SNESComputeFunction()`, `SNESGetFunctionDomainError()`

# External Links
$(_doc_external("Snes/SNESGetJacobianDomainError"))
"""
function SNESGetJacobianDomainError(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESGetJacobianDomainError(petsclib::$UnionPetscLib, snes::PetscSNES )
	domainerror_ = Ref{PetscBool}()

    @chk ccall(
               (:SNESGetJacobianDomainError, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{PetscBool}),
               snes, domainerror_,
              )

	domainerror = domainerror_[]

	return domainerror
end 

"""
	SNESLoad(petsclib::PetscLibType,snes::PetscSNES, viewer::PetscViewer) 
Loads a `SNES` that has been stored in `PETSCVIEWERBINARY` with `SNESView()`.

Collective

Input Parameters:
- `snes`   - the newly loaded `SNES`, this needs to have been created with `SNESCreate()` or
some related function before a call to `SNESLoad()`.
- `viewer` - binary file viewer, obtained from `PetscViewerBinaryOpen()`

Level: intermediate

-seealso: [](ch_snes), `SNES`, `PetscViewer`, `SNESCreate()`, `SNESType`, `PetscViewerBinaryOpen()`, `SNESView()`, `MatLoad()`, `VecLoad()`

# External Links
$(_doc_external("Snes/SNESLoad"))
"""
function SNESLoad(petsclib::PetscLibType, snes::PetscSNES, viewer::PetscViewer) end

@for_petsc function SNESLoad(petsclib::$UnionPetscLib, snes::PetscSNES, viewer::PetscViewer )

    @chk ccall(
               (:SNESLoad, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscViewer),
               snes, viewer,
              )


	return nothing
end 

"""
	SNESViewFromOptions(petsclib::PetscLibType,A::PetscSNES, obj::PetscObject, name::String) 
View a `SNES` based on values in the options database

Collective

Input Parameters:
- `A`    - the `SNES` context
- `obj`  - Optional object that provides the options prefix for the checks
- `name` - command line option

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESView`, `PetscObjectViewFromOptions()`, `SNESCreate()`

# External Links
$(_doc_external("Snes/SNESViewFromOptions"))
"""
function SNESViewFromOptions(petsclib::PetscLibType, A::PetscSNES, obj::PetscObject, name::String) end

@for_petsc function SNESViewFromOptions(petsclib::$UnionPetscLib, A::PetscSNES, obj::PetscObject, name::String )

    @chk ccall(
               (:SNESViewFromOptions, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

"""
	SNESView(petsclib::PetscLibType,snes::PetscSNES, viewer::PetscViewer) 
Prints or visualizes the `SNES` data structure.

Collective

Input Parameters:
- `snes`   - the `SNES` context
- `viewer` - the `PetscViewer`

Options Database Key:
- `-snes_view` - Calls `SNESView()` at end of `SNESSolve()`

Level: beginner

-seealso: [](ch_snes), `SNES`, `SNESLoad()`, `SNESCreate()`, `PetscViewerASCIIOpen()`

# External Links
$(_doc_external("Snes/SNESView"))
"""
function SNESView(petsclib::PetscLibType, snes::PetscSNES, viewer::PetscViewer) end

@for_petsc function SNESView(petsclib::$UnionPetscLib, snes::PetscSNES, viewer::PetscViewer )

    @chk ccall(
               (:SNESView, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscViewer),
               snes, viewer,
              )


	return nothing
end 

"""
	SNESAddOptionsChecker(petsclib::PetscLibType,snescheck::external) 
Adds an additional function to check for `SNES` options.

Not Collective

Input Parameter:
- `snescheck` - function that checks for options

Calling sequence of `snescheck`:
- `snes` - the `SNES` object for which it is checking options

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESSetFromOptions()`

# External Links
$(_doc_external("Snes/SNESAddOptionsChecker"))
"""
function SNESAddOptionsChecker(petsclib::PetscLibType, snescheck::external) end

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
	SNESSetUpMatrices(petsclib::PetscLibType,snes::PetscSNES) 
ensures that matrices are available for `SNES` Newton

Collective

Input Parameter:
- `snes` - `SNES` object to configure

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESSetUp()`

# External Links
$(_doc_external("Snes/SNESSetUpMatrices"))
"""
function SNESSetUpMatrices(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESSetUpMatrices(petsclib::$UnionPetscLib, snes::PetscSNES )

    @chk ccall(
               (:SNESSetUpMatrices, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	SNESMonitorSetFromOptions(petsclib::PetscLibType,snes::PetscSNES, name::String, help::String, manual::String, monitor::external, monitorsetup::external) 
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

-seealso: [](ch_snes), `PetscOptionsCreateViewer()`, `PetscOptionsGetReal()`, `PetscOptionsHasName()`, `PetscOptionsGetString()`,
`PetscOptionsGetIntArray()`, `PetscOptionsGetRealArray()`, `PetscOptionsBool()`
`PetscOptionsInt()`, `PetscOptionsString()`, `PetscOptionsReal()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`

# External Links
$(_doc_external("Snes/SNESMonitorSetFromOptions"))
"""
function SNESMonitorSetFromOptions(petsclib::PetscLibType, snes::PetscSNES, name::String, help::String, manual::String, monitor::external, monitorsetup::external) end

@for_petsc function SNESMonitorSetFromOptions(petsclib::$UnionPetscLib, snes::PetscSNES, name::String, help::String, manual::String, monitor::external, monitorsetup::external )

    @chk ccall(
               (:SNESMonitorSetFromOptions, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, external, external),
               snes, name, help, manual, monitor, monitorsetup,
              )


	return nothing
end 

"""
	SNESSetFromOptions(petsclib::PetscLibType,snes::PetscSNES) 
Sets various `SNES` and `KSP` parameters from user options.

Collective

Input Parameter:
- `snes` - the `SNES` context

Options Database Keys:
- `-snes_type <type>`                                                            - newtonls, newtontr, ngmres, ncg, nrichardson, qn, vi, fas, `SNESType` for complete list
- `-snes_rtol <rtol>`                                                            - relative decrease in tolerance norm from initial
- `-snes_atol <abstol>`                                                          - absolute tolerance of residual norm
- `-snes_stol <stol>`                                                            - convergence tolerance in terms of the norm of the change in the solution between steps
- `-snes_divergence_tolerance <divtol>`                                          - if the residual goes above divtol*rnorm0, exit with divergence
- `-snes_max_it <max_it>`                                                        - maximum number of iterations
- `-snes_max_funcs <max_funcs>`                                                  - maximum number of function evaluations
- `-snes_force_iteration <force>`                                                - force `SNESSolve()` to take at least one iteration
- `-snes_max_fail <max_fail>`                                                    - maximum number of line search failures allowed before stopping, default is none
- `-snes_max_linear_solve_fail`                                                  - number of linear solver failures before SNESSolve() stops
- `-snes_lag_preconditioner <lag>`                                               - how often preconditioner is rebuilt (use -1 to never rebuild)
- `-snes_lag_preconditioner_persists <true,false>`                               - retains the -snes_lag_preconditioner information across multiple SNESSolve()
- `-snes_lag_jacobian <lag>`                                                     - how often Jacobian is rebuilt (use -1 to never rebuild)
- `-snes_lag_jacobian_persists <true,false>`                                     - retains the -snes_lag_jacobian information across multiple SNESSolve()
- `-snes_convergence_test <default,skip,correct_pressure>`                       - convergence test in nonlinear solver. default `SNESConvergedDefault()`. skip `SNESConvergedSkip()` means continue iterating until max_it or some other criterion is reached, saving expense of convergence test. correct_pressure `SNESConvergedCorrectPressure()` has special handling of a pressure null space.
- `-snes_monitor [ascii][:filename][:viewer format]`                             - prints residual norm at each iteration. if no filename given prints to stdout
- `-snes_monitor_solution [ascii binary draw][:filename][:viewer format]`        - plots solution at each iteration
- `-snes_monitor_residual [ascii binary draw][:filename][:viewer format]`        - plots residual (not its norm) at each iteration
- `-snes_monitor_solution_update [ascii binary draw][:filename][:viewer format]` - plots update to solution at each iteration
- `-snes_monitor_lg_residualnorm`                                                - plots residual norm at each iteration
- `-snes_monitor_lg_range`                                                       - plots residual norm at each iteration
- `-snes_monitor_pause_final`                                                    - Pauses all monitor drawing after the solver ends
- `-snes_fd`                                                                     - use finite differences to compute Jacobian; very slow, only for testing
- `-snes_fd_color`                                                               - use finite differences with coloring to compute Jacobian
- `-snes_mf_ksp_monitor`                                                         - if using matrix-free multiply then print h at each `KSP` iteration
- `-snes_converged_reason`                                                       - print the reason for convergence/divergence after each solve
- `-npc_snes_type <type>`                                                        - the `SNES` type to use as a nonlinear preconditioner
- `-snes_test_jacobian <optional threshold>`                                     - compare the user provided Jacobian with one computed via finite differences to check for errors.  If a threshold is given, display only those entries whose difference is greater than the threshold.
- `-snes_test_jacobian_view`                                                     - display the user provided Jacobian, the finite difference Jacobian and the difference between them to help users detect the location of errors in the user provided Jacobian.

Options Database Keys for Eisenstat-Walker method:
- `-snes_ksp_ew`                       - use Eisenstat-Walker method for determining linear system convergence
- `-snes_ksp_ew_version ver`           - version of  Eisenstat-Walker method
- `-snes_ksp_ew_rtol0 <rtol0>`         - Sets rtol0
- `-snes_ksp_ew_rtolmax <rtolmax>`     - Sets rtolmax
- `-snes_ksp_ew_gamma <gamma>`         - Sets gamma
- `-snes_ksp_ew_alpha <alpha>`         - Sets alpha
- `-snes_ksp_ew_alpha2 <alpha2>`       - Sets alpha2
- `-snes_ksp_ew_threshold <threshold>` - Sets threshold

Level: beginner

-seealso: [](ch_snes), `SNESType`, `SNESSetOptionsPrefix()`, `SNESResetFromOptions()`, `SNES`, `SNESCreate()`, `MatCreateSNESMF()`, `MatFDColoring`

# External Links
$(_doc_external("Snes/SNESSetFromOptions"))
"""
function SNESSetFromOptions(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESSetFromOptions(petsclib::$UnionPetscLib, snes::PetscSNES )

    @chk ccall(
               (:SNESSetFromOptions, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	SNESResetFromOptions(petsclib::PetscLibType,snes::PetscSNES) 
Sets various `SNES` and `KSP` parameters from user options ONLY if the `SNESSetFromOptions()` was previously called

Collective

Input Parameter:
- `snes` - the `SNES` context

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESSetFromOptions()`, `SNESSetOptionsPrefix()`

# External Links
$(_doc_external("Snes/SNESResetFromOptions"))
"""
function SNESResetFromOptions(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESResetFromOptions(petsclib::$UnionPetscLib, snes::PetscSNES )

    @chk ccall(
               (:SNESResetFromOptions, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	SNESSetComputeApplicationContext(petsclib::PetscLibType,snes::PetscSNES, compute::external, destroy::PetscCtxDestroyFn) 
Sets an optional function to compute a user
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

-seealso: [](ch_snes), `SNESGetApplicationContext()`, `SNESSetApplicationContext()`, `PetscCtxDestroyFn`

# External Links
$(_doc_external("Snes/SNESSetComputeApplicationContext"))
"""
function SNESSetComputeApplicationContext(petsclib::PetscLibType, snes::PetscSNES, compute::external, destroy::PetscCtxDestroyFn) end

@for_petsc function SNESSetComputeApplicationContext(petsclib::$UnionPetscLib, snes::PetscSNES, compute::external, destroy::PetscCtxDestroyFn )

    @chk ccall(
               (:SNESSetComputeApplicationContext, $petsc_library),
               PetscErrorCode,
               (CSNES, external, Ptr{PetscCtxDestroyFn}),
               snes, compute, destroy,
              )


	return nothing
end 

"""
	SNESSetApplicationContext(petsclib::PetscLibType,snes::PetscSNES, ctx::Cvoid) 
Sets the optional user

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `ctx`  - the user context

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESSetComputeApplicationContext()`, `SNESGetApplicationContext()`

# External Links
$(_doc_external("Snes/SNESSetApplicationContext"))
"""
function SNESSetApplicationContext(petsclib::PetscLibType, snes::PetscSNES, ctx::Cvoid) end

@for_petsc function SNESSetApplicationContext(petsclib::$UnionPetscLib, snes::PetscSNES, ctx::Cvoid )

    @chk ccall(
               (:SNESSetApplicationContext, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Cvoid}),
               snes, ctx,
              )


	return nothing
end 

"""
	SNESGetApplicationContext(petsclib::PetscLibType,snes::PetscSNES, ctx::PeCtx) 
Gets the user
nonlinear solvers set with `SNESGetApplicationContext()` or `SNESSetComputeApplicationContext()`

Not Collective

Input Parameter:
- `snes` - `SNES` context

Output Parameter:
- `ctx` - user context

Level: intermediate

-seealso: [](ch_snes), `SNESSetApplicationContext()`, `SNESSetComputeApplicationContext()`

# External Links
$(_doc_external("Snes/SNESGetApplicationContext"))
"""
function SNESGetApplicationContext(petsclib::PetscLibType, snes::PetscSNES, ctx::PeCtx) end

@for_petsc function SNESGetApplicationContext(petsclib::$UnionPetscLib, snes::PetscSNES, ctx::PeCtx )

    @chk ccall(
               (:SNESGetApplicationContext, $petsc_library),
               PetscErrorCode,
               (CSNES, PeCtx),
               snes, ctx,
              )


	return nothing
end 

"""
	SNESSetUseMatrixFree(petsclib::PetscLibType,snes::PetscSNES, mf_operator::PetscBool, mf::PetscBool) 
indicates that `SNES` should use matrix

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

-seealso: [](ch_snes), `SNES`, `SNESGetUseMatrixFree()`, `MatCreateSNESMF()`, `SNESComputeJacobianDefaultColor()`, `MatFDColoring`

# External Links
$(_doc_external("Snes/SNESSetUseMatrixFree"))
"""
function SNESSetUseMatrixFree(petsclib::PetscLibType, snes::PetscSNES, mf_operator::PetscBool, mf::PetscBool) end

@for_petsc function SNESSetUseMatrixFree(petsclib::$UnionPetscLib, snes::PetscSNES, mf_operator::PetscBool, mf::PetscBool )

    @chk ccall(
               (:SNESSetUseMatrixFree, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool, PetscBool),
               snes, mf_operator, mf,
              )


	return nothing
end 

"""
	mf_operator::PetscBool,mf::PetscBool = SNESGetUseMatrixFree(petsclib::PetscLibType,snes::PetscSNES) 
indicates if the `SNES` uses matrix

Not Collective, but the resulting flags will be the same on all MPI processes

Input Parameter:
- `snes` - `SNES` context

Output Parameters:
- `mf_operator` - use matrix-free only for the Amat used by `SNESSetJacobian()`, this means the user provided Pmat will continue to be used
- `mf`          - use matrix-free for both the Amat and Pmat used by `SNESSetJacobian()`, both the Amat and Pmat set in `SNESSetJacobian()` will be ignored

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESSetUseMatrixFree()`, `MatCreateSNESMF()`

# External Links
$(_doc_external("Snes/SNESGetUseMatrixFree"))
"""
function SNESGetUseMatrixFree(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESGetUseMatrixFree(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	iter::PetscInt = SNESGetIterationNumber(petsclib::PetscLibType,snes::PetscSNES) 
Gets the number of nonlinear iterations completed in the current or most recent `SNESSolve()`

Not Collective

Input Parameter:
- `snes` - `SNES` context

Output Parameter:
- `iter` - iteration number

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESSolve()`, `SNESSetLagJacobian()`, `SNESGetLinearSolveIterations()`, `SNESSetMonitor()`

# External Links
$(_doc_external("Snes/SNESGetIterationNumber"))
"""
function SNESGetIterationNumber(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESGetIterationNumber(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESSetIterationNumber(petsclib::PetscLibType,snes::PetscSNES, iter::PetscInt) 
Sets the current iteration number.

Not Collective

Input Parameters:
- `snes` - `SNES` context
- `iter` - iteration number

Level: developer

-seealso: [](ch_snes), `SNESGetLinearSolveIterations()`

# External Links
$(_doc_external("Snes/SNESSetIterationNumber"))
"""
function SNESSetIterationNumber(petsclib::PetscLibType, snes::PetscSNES, iter::PetscInt) end

@for_petsc function SNESSetIterationNumber(petsclib::$UnionPetscLib, snes::PetscSNES, iter::$PetscInt )

    @chk ccall(
               (:SNESSetIterationNumber, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt),
               snes, iter,
              )


	return nothing
end 

"""
	nfails::PetscInt = SNESGetNonlinearStepFailures(petsclib::PetscLibType,snes::PetscSNES) 
Gets the number of unsuccessful steps
attempted by the nonlinear solver in the current or most recent `SNESSolve()` .

Not Collective

Input Parameter:
- `snes` - `SNES` context

Output Parameter:
- `nfails` - number of unsuccessful steps attempted

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESGetMaxLinearSolveFailures()`, `SNESGetLinearSolveIterations()`, `SNESSetMaxLinearSolveFailures()`, `SNESGetLinearSolveFailures()`,
`SNESSetMaxNonlinearStepFailures()`, `SNESGetMaxNonlinearStepFailures()`

# External Links
$(_doc_external("Snes/SNESGetNonlinearStepFailures"))
"""
function SNESGetNonlinearStepFailures(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESGetNonlinearStepFailures(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESSetMaxNonlinearStepFailures(petsclib::PetscLibType,snes::PetscSNES, maxFails::PetscInt) 
Sets the maximum number of unsuccessful steps
attempted by the nonlinear solver before it gives up and returns unconverged or generates an error

Not Collective

Input Parameters:
- `snes`     - `SNES` context
- `maxFails` - maximum of unsuccessful steps allowed, use `PETSC_UNLIMITED` to have no limit on the number of failures

Options Database Key:
- `-snes_max_fail <n>` - maximum number of unsuccessful steps allowed

Level: intermediate

-seealso: [](ch_snes), `SNESSetErrorIfNotConverged()`, `SNESGetMaxLinearSolveFailures()`, `SNESGetLinearSolveIterations()`, `SNESSetMaxLinearSolveFailures()`, `SNESGetLinearSolveFailures()`,
`SNESGetMaxNonlinearStepFailures()`, `SNESGetNonlinearStepFailures()`

# External Links
$(_doc_external("Snes/SNESSetMaxNonlinearStepFailures"))
"""
function SNESSetMaxNonlinearStepFailures(petsclib::PetscLibType, snes::PetscSNES, maxFails::PetscInt) end

@for_petsc function SNESSetMaxNonlinearStepFailures(petsclib::$UnionPetscLib, snes::PetscSNES, maxFails::$PetscInt )

    @chk ccall(
               (:SNESSetMaxNonlinearStepFailures, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt),
               snes, maxFails,
              )


	return nothing
end 

"""
	maxFails::PetscInt = SNESGetMaxNonlinearStepFailures(petsclib::PetscLibType,snes::PetscSNES) 
Gets the maximum number of unsuccessful steps
attempted by the nonlinear solver before it gives up and returns unconverged or generates an error

Not Collective

Input Parameter:
- `snes` - `SNES` context

Output Parameter:
- `maxFails` - maximum of unsuccessful steps

Level: intermediate

-seealso: [](ch_snes), `SNESSetErrorIfNotConverged()`, `SNESGetMaxLinearSolveFailures()`, `SNESGetLinearSolveIterations()`, `SNESSetMaxLinearSolveFailures()`, `SNESGetLinearSolveFailures()`,
`SNESSetMaxNonlinearStepFailures()`, `SNESGetNonlinearStepFailures()`

# External Links
$(_doc_external("Snes/SNESGetMaxNonlinearStepFailures"))
"""
function SNESGetMaxNonlinearStepFailures(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESGetMaxNonlinearStepFailures(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	nfuncs::PetscInt = SNESGetNumberFunctionEvals(petsclib::PetscLibType,snes::PetscSNES) 
Gets the number of user provided function evaluations
done by the `SNES` object in the current or most recent `SNESSolve()`

Not Collective

Input Parameter:
- `snes` - `SNES` context

Output Parameter:
- `nfuncs` - number of evaluations

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESGetMaxLinearSolveFailures()`, `SNESGetLinearSolveIterations()`, `SNESSetMaxLinearSolveFailures()`, `SNESGetLinearSolveFailures()`, `SNESSetCountersReset()`

# External Links
$(_doc_external("Snes/SNESGetNumberFunctionEvals"))
"""
function SNESGetNumberFunctionEvals(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESGetNumberFunctionEvals(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	nfails::PetscInt = SNESGetLinearSolveFailures(petsclib::PetscLibType,snes::PetscSNES) 
Gets the number of failed (non
linear solvers in the current or most recent `SNESSolve()`

Not Collective

Input Parameter:
- `snes` - `SNES` context

Output Parameter:
- `nfails` - number of failed solves

Options Database Key:
- `-snes_max_linear_solve_fail <num>` - The number of failures before the solve is terminated

Level: intermediate

-seealso: [](ch_snes), `SNESGetMaxLinearSolveFailures()`, `SNESGetLinearSolveIterations()`, `SNESSetMaxLinearSolveFailures()`

# External Links
$(_doc_external("Snes/SNESGetLinearSolveFailures"))
"""
function SNESGetLinearSolveFailures(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESGetLinearSolveFailures(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESSetMaxLinearSolveFailures(petsclib::PetscLibType,snes::PetscSNES, maxFails::PetscInt) 
the number of failed linear solve attempts
allowed before `SNES` returns with a diverged reason of `SNES_DIVERGED_LINEAR_SOLVE`

Logically Collective

Input Parameters:
- `snes`     - `SNES` context
- `maxFails` - maximum allowed linear solve failures, use `PETSC_UNLIMITED` to have no limit on the number of failures

Options Database Key:
- `-snes_max_linear_solve_fail <num>` - The number of failures before the solve is terminated

Level: intermediate

-seealso: [](ch_snes), `SNESSetErrorIfNotConverged()`, `SNESGetLinearSolveFailures()`, `SNESGetMaxLinearSolveFailures()`, `SNESGetLinearSolveIterations()`

# External Links
$(_doc_external("Snes/SNESSetMaxLinearSolveFailures"))
"""
function SNESSetMaxLinearSolveFailures(petsclib::PetscLibType, snes::PetscSNES, maxFails::PetscInt) end

@for_petsc function SNESSetMaxLinearSolveFailures(petsclib::$UnionPetscLib, snes::PetscSNES, maxFails::$PetscInt )

    @chk ccall(
               (:SNESSetMaxLinearSolveFailures, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt),
               snes, maxFails,
              )


	return nothing
end 

"""
	maxFails::PetscInt = SNESGetMaxLinearSolveFailures(petsclib::PetscLibType,snes::PetscSNES) 
gets the maximum number of linear solve failures that
are allowed before `SNES` returns as unsuccessful

Not Collective

Input Parameter:
- `snes` - `SNES` context

Output Parameter:
- `maxFails` - maximum of unsuccessful solves allowed

Level: intermediate

-seealso: [](ch_snes), `SNESSetErrorIfNotConverged()`, `SNESGetLinearSolveFailures()`, `SNESGetLinearSolveIterations()`, `SNESSetMaxLinearSolveFailures()`,

# External Links
$(_doc_external("Snes/SNESGetMaxLinearSolveFailures"))
"""
function SNESGetMaxLinearSolveFailures(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESGetMaxLinearSolveFailures(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	lits::PetscInt = SNESGetLinearSolveIterations(petsclib::PetscLibType,snes::PetscSNES) 
Gets the total number of linear iterations
used by the nonlinear solver in the most recent `SNESSolve()`

Not Collective

Input Parameter:
- `snes` - `SNES` context

Output Parameter:
- `lits` - number of linear iterations

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESGetIterationNumber()`, `SNESGetLinearSolveFailures()`, `SNESGetMaxLinearSolveFailures()`, `SNESSetCountersReset()`

# External Links
$(_doc_external("Snes/SNESGetLinearSolveIterations"))
"""
function SNESGetLinearSolveIterations(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESGetLinearSolveIterations(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESSetCountersReset(petsclib::PetscLibType,snes::PetscSNES, reset::PetscBool) 
Sets whether or not the counters for linear iterations and function evaluations
are reset every time `SNESSolve()` is called.

Logically Collective

Input Parameters:
- `snes`  - `SNES` context
- `reset` - whether to reset the counters or not, defaults to `PETSC_TRUE`

Level: developer

-seealso: [](ch_snes), `SNESGetNumberFunctionEvals()`, `SNESGetLinearSolveIterations()`, `SNESGetNPC()`

# External Links
$(_doc_external("Snes/SNESSetCountersReset"))
"""
function SNESSetCountersReset(petsclib::PetscLibType, snes::PetscSNES, reset::PetscBool) end

@for_petsc function SNESSetCountersReset(petsclib::$UnionPetscLib, snes::PetscSNES, reset::PetscBool )

    @chk ccall(
               (:SNESSetCountersReset, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, reset,
              )


	return nothing
end 

"""
	SNESResetCounters(petsclib::PetscLibType,snes::PetscSNES) 
Reset counters for linear iterations and function evaluations.

Logically Collective

Input Parameters:
- `snes` - `SNES` context

Level: developer

-seealso: [](ch_snes), `SNESGetNumberFunctionEvals()`, `SNESGetLinearSolveIterations()`, `SNESGetNPC()`

# External Links
$(_doc_external("Snes/SNESResetCounters"))
"""
function SNESResetCounters(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESResetCounters(petsclib::$UnionPetscLib, snes::PetscSNES )

    @chk ccall(
               (:SNESResetCounters, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	SNESSetKSP(petsclib::PetscLibType,snes::PetscSNES, ksp::PetscKSP) 
Sets a `KSP` context for the `SNES` object to use

Not Collective, but the `SNES` and `KSP` objects must live on the same `MPI_Comm`

Input Parameters:
- `snes` - the `SNES` context
- `ksp`  - the `KSP` context

Level: developer

-seealso: [](ch_snes), `SNES`, `KSP`, `KSPGetPC()`, `SNESCreate()`, `KSPCreate()`

# External Links
$(_doc_external("Snes/SNESSetKSP"))
"""
function SNESSetKSP(petsclib::PetscLibType, snes::PetscSNES, ksp::PetscKSP) end

@for_petsc function SNESSetKSP(petsclib::$UnionPetscLib, snes::PetscSNES, ksp::PetscKSP )

    @chk ccall(
               (:SNESSetKSP, $petsc_library),
               PetscErrorCode,
               (CSNES, CKSP),
               snes, ksp,
              )


	return nothing
end 

"""
	SNESParametersInitialize(petsclib::PetscLibType,snes::PetscSNES) 
Sets all the parameters in `snes` to their default value (when `SNESCreate()` was called) if they
currently contain default values

Collective

Input Parameter:
- `snes` - the `SNES` object

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESSolve()`, `SNESDestroy()`, `SNESSetLagPreconditioner()`, `SNESSetLagJacobian()`,
`PetscObjectParameterSetDefault()`

# External Links
$(_doc_external("Snes/SNESParametersInitialize"))
"""
function SNESParametersInitialize(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESParametersInitialize(petsclib::$UnionPetscLib, snes::PetscSNES )

    @chk ccall(
               (:SNESParametersInitialize, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	outsnes::PetscSNES = SNESCreate(petsclib::PetscLibType,comm::MPI_Comm) 
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

-seealso: [](ch_snes), `SNES`, `SNESSolve()`, `SNESDestroy()`, `SNESSetLagPreconditioner()`, `SNESSetLagJacobian()`

# External Links
$(_doc_external("Snes/SNESCreate"))
"""
function SNESCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function SNESCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	outsnes_ = Ref{CSNES}()

    @chk ccall(
               (:SNESCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{CSNES}),
               comm, outsnes_,
              )

	outsnes = PetscSNES(outsnes_[], petsclib)

	return outsnes
end 

"""
	SNESSetFunction(petsclib::PetscLibType,snes::PetscSNES, r::PetscVec, f::SNESFunctionFn, ctx::Cvoid) 
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

-seealso: [](ch_snes), `SNES`, `SNESGetFunction()`, `SNESComputeFunction()`, `SNESSetJacobian()`, `SNESSetPicard()`, `SNESFunctionFn`

# External Links
$(_doc_external("Snes/SNESSetFunction"))
"""
function SNESSetFunction(petsclib::PetscLibType, snes::PetscSNES, r::PetscVec, f::SNESFunctionFn, ctx::Cvoid) end

@for_petsc function SNESSetFunction(petsclib::$UnionPetscLib, snes::PetscSNES, r::PetscVec, f::SNESFunctionFn, ctx::Cvoid )

    @chk ccall(
               (:SNESSetFunction, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, Ptr{SNESFunctionFn}, Ptr{Cvoid}),
               snes, r, f, ctx,
              )


	return nothing
end 

"""
	SNESSetInitialFunction(petsclib::PetscLibType,snes::PetscSNES, f::PetscVec) 
Set an already computed function evaluation at the initial guess to be reused by `SNESSolve()`.

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `f`    - vector to store function value

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESSetFunction()`, `SNESComputeFunction()`, `SNESSetInitialFunctionNorm()`

# External Links
$(_doc_external("Snes/SNESSetInitialFunction"))
"""
function SNESSetInitialFunction(petsclib::PetscLibType, snes::PetscSNES, f::PetscVec) end

@for_petsc function SNESSetInitialFunction(petsclib::$UnionPetscLib, snes::PetscSNES, f::PetscVec )

    @chk ccall(
               (:SNESSetInitialFunction, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec),
               snes, f,
              )


	return nothing
end 

"""
	SNESSetNormSchedule(petsclib::PetscLibType,snes::PetscSNES, normschedule::SNESNormSchedule) 
Sets the `SNESNormSchedule` used in convergence and monitoring
of the `SNES` method, when norms are computed in the solving process

Logically Collective

Input Parameters:
- `snes`         - the `SNES` context
- `normschedule` - the frequency of norm computation

Options Database Key:
- `-snes_norm_schedule <none, always, initialonly, finalonly, initialfinalonly>` - set the schedule

Level: advanced

-seealso: [](ch_snes), `SNESNormSchedule`, `SNESGetNormSchedule()`, `SNESComputeFunction()`, `VecNorm()`, `SNESSetFunction()`, `SNESSetInitialFunction()`

# External Links
$(_doc_external("Snes/SNESSetNormSchedule"))
"""
function SNESSetNormSchedule(petsclib::PetscLibType, snes::PetscSNES, normschedule::SNESNormSchedule) end

@for_petsc function SNESSetNormSchedule(petsclib::$UnionPetscLib, snes::PetscSNES, normschedule::SNESNormSchedule )

    @chk ccall(
               (:SNESSetNormSchedule, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESNormSchedule),
               snes, normschedule,
              )


	return nothing
end 

"""
	SNESGetNormSchedule(petsclib::PetscLibType,snes::PetscSNES, normschedule::SNESNormSchedule) 
Gets the `SNESNormSchedule` used in convergence and monitoring
of the `SNES` method.

Logically Collective

Input Parameters:
- `snes`         - the `SNES` context
- `normschedule` - the type of the norm used

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESSetNormSchedule()`, `SNESComputeFunction()`, `VecNorm()`, `SNESSetFunction()`, `SNESSetInitialFunction()`, `SNESNormSchedule`

# External Links
$(_doc_external("Snes/SNESGetNormSchedule"))
"""
function SNESGetNormSchedule(petsclib::PetscLibType, snes::PetscSNES, normschedule::SNESNormSchedule) end

@for_petsc function SNESGetNormSchedule(petsclib::$UnionPetscLib, snes::PetscSNES, normschedule::SNESNormSchedule )

    @chk ccall(
               (:SNESGetNormSchedule, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{SNESNormSchedule}),
               snes, normschedule,
              )


	return nothing
end 

"""
	SNESSetFunctionNorm(petsclib::PetscLibType,snes::PetscSNES, norm::PetscReal) 
Sets the last computed residual norm.

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `norm` - the value of the norm

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESGetNormSchedule()`, `SNESComputeFunction()`, `VecNorm()`, `SNESSetFunction()`, `SNESSetInitialFunction()`, `SNESNormSchedule`

# External Links
$(_doc_external("Snes/SNESSetFunctionNorm"))
"""
function SNESSetFunctionNorm(petsclib::PetscLibType, snes::PetscSNES, norm::PetscReal) end

@for_petsc function SNESSetFunctionNorm(petsclib::$UnionPetscLib, snes::PetscSNES, norm::$PetscReal )

    @chk ccall(
               (:SNESSetFunctionNorm, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscReal),
               snes, norm,
              )


	return nothing
end 

"""
	norm::PetscReal = SNESGetFunctionNorm(petsclib::PetscLibType,snes::PetscSNES) 
Gets the last computed norm of the residual

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `norm` - the last computed residual norm

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESSetNormSchedule()`, `SNESComputeFunction()`, `VecNorm()`, `SNESSetFunction()`, `SNESSetInitialFunction()`, `SNESNormSchedule`

# External Links
$(_doc_external("Snes/SNESGetFunctionNorm"))
"""
function SNESGetFunctionNorm(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESGetFunctionNorm(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	ynorm::PetscReal = SNESGetUpdateNorm(petsclib::PetscLibType,snes::PetscSNES) 
Gets the last computed norm of the solution update

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `ynorm` - the last computed update norm

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESSetNormSchedule()`, `SNESComputeFunction()`, `SNESGetFunctionNorm()`

# External Links
$(_doc_external("Snes/SNESGetUpdateNorm"))
"""
function SNESGetUpdateNorm(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESGetUpdateNorm(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	xnorm::PetscReal = SNESGetSolutionNorm(petsclib::PetscLibType,snes::PetscSNES) 
Gets the last computed norm of the solution

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `xnorm` - the last computed solution norm

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESSetNormSchedule()`, `SNESComputeFunction()`, `SNESGetFunctionNorm()`, `SNESGetUpdateNorm()`

# External Links
$(_doc_external("Snes/SNESGetSolutionNorm"))
"""
function SNESGetSolutionNorm(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESGetSolutionNorm(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESSetFunctionType(petsclib::PetscLibType,snes::PetscSNES, type::SNESFunctionType) 
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

-seealso: [](ch_snes), `SNES`, `SNESFunctionType`, `SNESGetNormSchedule()`, `SNESComputeFunction()`, `VecNorm()`, `SNESSetFunction()`, `SNESSetInitialFunction()`, `SNESNormSchedule`

# External Links
$(_doc_external("Snes/SNESSetFunctionType"))
"""
function SNESSetFunctionType(petsclib::PetscLibType, snes::PetscSNES, type::SNESFunctionType) end

@for_petsc function SNESSetFunctionType(petsclib::$UnionPetscLib, snes::PetscSNES, type::SNESFunctionType )

    @chk ccall(
               (:SNESSetFunctionType, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESFunctionType),
               snes, type,
              )


	return nothing
end 

"""
	SNESGetFunctionType(petsclib::PetscLibType,snes::PetscSNES, type::SNESFunctionType) 
Gets the `SNESFunctionType` used in convergence and monitoring set with `SNESSetFunctionType()`
of the SNES method.

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `type` - the type of the function evaluation, see `SNESSetFunctionType()`

Level: advanced

-seealso: [](ch_snes), `SNESSetFunctionType()`, `SNESFunctionType`, `SNESSetNormSchedule()`, `SNESComputeFunction()`, `VecNorm()`, `SNESSetFunction()`, `SNESSetInitialFunction()`, `SNESNormSchedule`

# External Links
$(_doc_external("Snes/SNESGetFunctionType"))
"""
function SNESGetFunctionType(petsclib::PetscLibType, snes::PetscSNES, type::SNESFunctionType) end

@for_petsc function SNESGetFunctionType(petsclib::$UnionPetscLib, snes::PetscSNES, type::SNESFunctionType )

    @chk ccall(
               (:SNESGetFunctionType, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{SNESFunctionType}),
               snes, type,
              )


	return nothing
end 

"""
	SNESSetNGS(petsclib::PetscLibType,snes::PetscSNES, f::SNESNGSFn, ctx::Cvoid) 
Sets the user nonlinear Gauss
use with composed nonlinear solvers.

Input Parameters:
- `snes` - the `SNES` context, usually of the `SNESType` `SNESNGS`
- `f`    - function evaluation routine to apply Gauss-Seidel, see `SNESNGSFn` for calling sequence
- `ctx`  - [optional] user-defined context for private data for the smoother evaluation routine (may be `NULL`)

Level: intermediate

-seealso: [](ch_snes), `SNESNGS`, `SNESGetNGS()`, `SNESNCG`, `SNESGetFunction()`, `SNESComputeNGS()`, `SNESNGSFn`

# External Links
$(_doc_external("Snes/SNESSetNGS"))
"""
function SNESSetNGS(petsclib::PetscLibType, snes::PetscSNES, f::SNESNGSFn, ctx::Cvoid) end

@for_petsc function SNESSetNGS(petsclib::$UnionPetscLib, snes::PetscSNES, f::SNESNGSFn, ctx::Cvoid )

    @chk ccall(
               (:SNESSetNGS, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{SNESNGSFn}, Ptr{Cvoid}),
               snes, f, ctx,
              )


	return nothing
end 

"""
	SNESPicardComputeMFFunction(petsclib::PetscLibType,snes::PetscSNES, x::PetscVec, f::PetscVec, ctx::Cvoid) 

# External Links
$(_doc_external("Snes/SNESPicardComputeMFFunction"))
"""
function SNESPicardComputeMFFunction(petsclib::PetscLibType, snes::PetscSNES, x::PetscVec, f::PetscVec, ctx::Cvoid) end

@for_petsc function SNESPicardComputeMFFunction(petsclib::$UnionPetscLib, snes::PetscSNES, x::PetscVec, f::PetscVec, ctx::Cvoid )

    @chk ccall(
               (:SNESPicardComputeMFFunction, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec, Ptr{Cvoid}),
               snes, x, f, ctx,
              )


	return nothing
end 

"""
	SNESPicardComputeFunction(petsclib::PetscLibType,snes::PetscSNES, x::PetscVec, f::PetscVec, ctx::Cvoid) 

# External Links
$(_doc_external("Snes/SNESPicardComputeFunction"))
"""
function SNESPicardComputeFunction(petsclib::PetscLibType, snes::PetscSNES, x::PetscVec, f::PetscVec, ctx::Cvoid) end

@for_petsc function SNESPicardComputeFunction(petsclib::$UnionPetscLib, snes::PetscSNES, x::PetscVec, f::PetscVec, ctx::Cvoid )

    @chk ccall(
               (:SNESPicardComputeFunction, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec, Ptr{Cvoid}),
               snes, x, f, ctx,
              )


	return nothing
end 

"""
	SNESPicardComputeJacobian(petsclib::PetscLibType,snes::PetscSNES, x1::PetscVec, J::PetscMat, B::PetscMat, ctx::Cvoid) 

# External Links
$(_doc_external("Snes/SNESPicardComputeJacobian"))
"""
function SNESPicardComputeJacobian(petsclib::PetscLibType, snes::PetscSNES, x1::PetscVec, J::PetscMat, B::PetscMat, ctx::Cvoid) end

@for_petsc function SNESPicardComputeJacobian(petsclib::$UnionPetscLib, snes::PetscSNES, x1::PetscVec, J::PetscMat, B::PetscMat, ctx::Cvoid )

    @chk ccall(
               (:SNESPicardComputeJacobian, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CMat, CMat, Ptr{Cvoid}),
               snes, x1, J, B, ctx,
              )


	return nothing
end 

"""
	SNESSetPicard(petsclib::PetscLibType,snes::PetscSNES, r::PetscVec, bp::SNESFunctionFn, Amat::PetscMat, Pmat::PetscMat, J::SNESJacobianFn, ctx::Cvoid) 
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

-seealso: [](ch_snes), `SNES`, `SNESGetFunction()`, `SNESSetFunction()`, `SNESComputeFunction()`, `SNESSetJacobian()`, `SNESGetPicard()`, `SNESLineSearchPreCheckPicard()`,
`SNESFunctionFn`, `SNESJacobianFn`

# External Links
$(_doc_external("Snes/SNESSetPicard"))
"""
function SNESSetPicard(petsclib::PetscLibType, snes::PetscSNES, r::PetscVec, bp::SNESFunctionFn, Amat::PetscMat, Pmat::PetscMat, J::SNESJacobianFn, ctx::Cvoid) end

@for_petsc function SNESSetPicard(petsclib::$UnionPetscLib, snes::PetscSNES, r::PetscVec, bp::SNESFunctionFn, Amat::PetscMat, Pmat::PetscMat, J::SNESJacobianFn, ctx::Cvoid )

    @chk ccall(
               (:SNESSetPicard, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, Ptr{SNESFunctionFn}, CMat, CMat, Ptr{SNESJacobianFn}, Ptr{Cvoid}),
               snes, r, bp, Amat, Pmat, J, ctx,
              )


	return nothing
end 

"""
	SNESGetPicard(petsclib::PetscLibType,snes::PetscSNES, r::PetscVec, f::SNESFunctionFn, Amat::PetscMat, Pmat::PetscMat, J::SNESJacobianFn, ctx::Cvoid) 
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

-seealso: [](ch_snes), `SNESSetFunction()`, `SNESSetPicard()`, `SNESGetFunction()`, `SNESGetJacobian()`, `SNESGetDM()`, `SNESFunctionFn`, `SNESJacobianFn`

# External Links
$(_doc_external("Snes/SNESGetPicard"))
"""
function SNESGetPicard(petsclib::PetscLibType, snes::PetscSNES, r::PetscVec, f::SNESFunctionFn, Amat::PetscMat, Pmat::PetscMat, J::SNESJacobianFn, ctx::Cvoid) end

@for_petsc function SNESGetPicard(petsclib::$UnionPetscLib, snes::PetscSNES, r::PetscVec, f::SNESFunctionFn, Amat::PetscMat, Pmat::PetscMat, J::SNESJacobianFn, ctx::Cvoid )
	r_ = Ref(r.ptr)
	Amat_ = Ref(Amat.ptr)
	Pmat_ = Ref(Pmat.ptr)

    @chk ccall(
               (:SNESGetPicard, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CVec}, SNESFunctionFn, Ptr{CMat}, Ptr{CMat}, SNESJacobianFn, Cvoid),
               snes, r_, f, Amat_, Pmat_, J, ctx,
              )

	r.ptr = C_NULL
	Amat.ptr = C_NULL
	Pmat.ptr = C_NULL

	return nothing
end 

"""
	SNESSetComputeInitialGuess(petsclib::PetscLibType,snes::PetscSNES, func::SNESInitialGuessFn, ctx::Cvoid) 
Sets a routine used to compute an initial guess for the nonlinear problem

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `func` - function evaluation routine, see `SNESInitialGuessFn` for the calling sequence
- `ctx`  - [optional] user-defined context for private data for the
function evaluation routine (may be `NULL`)

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESSolve()`, `SNESSetFunction()`, `SNESGetFunction()`, `SNESComputeFunction()`, `SNESSetJacobian()`, `SNESInitialGuessFn`

# External Links
$(_doc_external("Snes/SNESSetComputeInitialGuess"))
"""
function SNESSetComputeInitialGuess(petsclib::PetscLibType, snes::PetscSNES, func::SNESInitialGuessFn, ctx::Cvoid) end

@for_petsc function SNESSetComputeInitialGuess(petsclib::$UnionPetscLib, snes::PetscSNES, func::SNESInitialGuessFn, ctx::Cvoid )

    @chk ccall(
               (:SNESSetComputeInitialGuess, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{SNESInitialGuessFn}, Ptr{Cvoid}),
               snes, func, ctx,
              )


	return nothing
end 

"""
	SNESGetRhs(petsclib::PetscLibType,snes::PetscSNES, rhs::PetscVec) 
Gets the vector for solving F(x) = `rhs`. If `rhs` is not set
it assumes a zero right-hand side.

Logically Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `rhs` - the right-hand side vector or `NULL` if there is no right-hand side vector

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESGetSolution()`, `SNESGetFunction()`, `SNESComputeFunction()`, `SNESSetJacobian()`, `SNESSetFunction()`

# External Links
$(_doc_external("Snes/SNESGetRhs"))
"""
function SNESGetRhs(petsclib::PetscLibType, snes::PetscSNES, rhs::PetscVec) end

@for_petsc function SNESGetRhs(petsclib::$UnionPetscLib, snes::PetscSNES, rhs::PetscVec )
	rhs_ = Ref(rhs.ptr)

    @chk ccall(
               (:SNESGetRhs, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CVec}),
               snes, rhs_,
              )

	rhs.ptr = C_NULL

	return nothing
end 

"""
	SNESComputeFunction(petsclib::PetscLibType,snes::PetscSNES, x::PetscVec, y::PetscVec) 
Calls the function that has been set with `SNESSetFunction()`.

Collective

Input Parameters:
- `snes` - the `SNES` context
- `x`    - input vector

Output Parameter:
- `y` - function vector, as set by `SNESSetFunction()`

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESSetFunction()`, `SNESGetFunction()`, `SNESComputeMFFunction()`

# External Links
$(_doc_external("Snes/SNESComputeFunction"))
"""
function SNESComputeFunction(petsclib::PetscLibType, snes::PetscSNES, x::PetscVec, y::PetscVec) end

@for_petsc function SNESComputeFunction(petsclib::$UnionPetscLib, snes::PetscSNES, x::PetscVec, y::PetscVec )

    @chk ccall(
               (:SNESComputeFunction, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec),
               snes, x, y,
              )


	return nothing
end 

"""
	SNESComputeMFFunction(petsclib::PetscLibType,snes::PetscSNES, x::PetscVec, y::PetscVec) 
Calls the function that has been set with `DMSNESSetMFFunction()`.

Collective

Input Parameters:
- `snes` - the `SNES` context
- `x`    - input vector

Output Parameter:
- `y` - output vector

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESSetFunction()`, `SNESGetFunction()`, `SNESComputeFunction()`, `MatCreateSNESMF()`, `DMSNESSetMFFunction()`

# External Links
$(_doc_external("Snes/SNESComputeMFFunction"))
"""
function SNESComputeMFFunction(petsclib::PetscLibType, snes::PetscSNES, x::PetscVec, y::PetscVec) end

@for_petsc function SNESComputeMFFunction(petsclib::$UnionPetscLib, snes::PetscSNES, x::PetscVec, y::PetscVec )

    @chk ccall(
               (:SNESComputeMFFunction, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec),
               snes, x, y,
              )


	return nothing
end 

"""
	SNESComputeNGS(petsclib::PetscLibType,snes::PetscSNES, b::PetscVec, x::PetscVec) 
Calls the Gauss

Collective

Input Parameters:
- `snes` - the `SNES` context
- `x`    - input vector
- `b`    - rhs vector

Output Parameter:
- `x` - new solution vector

Level: developer

-seealso: [](ch_snes), `SNESNGSFn`, `SNESSetNGS()`, `SNESComputeFunction()`, `SNESNGS`

# External Links
$(_doc_external("Snes/SNESComputeNGS"))
"""
function SNESComputeNGS(petsclib::PetscLibType, snes::PetscSNES, b::PetscVec, x::PetscVec) end

@for_petsc function SNESComputeNGS(petsclib::$UnionPetscLib, snes::PetscSNES, b::PetscVec, x::PetscVec )

    @chk ccall(
               (:SNESComputeNGS, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec),
               snes, b, x,
              )


	return nothing
end 

"""
	SNESTestFunction(petsclib::PetscLibType,snes::PetscSNES) 
Computes the difference between the computed and finite

Collective

Input Parameter:
- `snes` - the `SNES` context

Options Database Keys:
- `-snes_test_function`      - compare the user provided function with one compute via finite differences to check for errors.
- `-snes_test_function_view` - display the user provided function, the finite difference function and the difference

Level: developer

-seealso: [](ch_snes), `SNESTestJacobian()`, `SNESSetFunction()`, `SNESComputeFunction()`

# External Links
$(_doc_external("Snes/SNESTestFunction"))
"""
function SNESTestFunction(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESTestFunction(petsclib::$UnionPetscLib, snes::PetscSNES )

    @chk ccall(
               (:SNESTestFunction, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	Jnorm::PetscReal,diffNorm::PetscReal = SNESTestJacobian(petsclib::PetscLibType,snes::PetscSNES) 
Computes the difference between the computed and finite

Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameters:
- `Jnorm`    - the Frobenius norm of the computed Jacobian, or `NULL`
- `diffNorm` - the Frobenius norm of the difference of the computed and finite-difference Jacobians, or `NULL`

Options Database Keys:
- `-snes_test_jacobian <optional threshold>` - compare the user provided Jacobian with one compute via finite differences to check for errors.  If a threshold is given, display only those entries whose difference is greater than the threshold.
- `-snes_test_jacobian_view`                 - display the user provided Jacobian, the finite difference Jacobian and the difference

Level: developer

-seealso: [](ch_snes), `SNESTestFunction()`, `SNESSetJacobian()`, `SNESComputeJacobian()`

# External Links
$(_doc_external("Snes/SNESTestJacobian"))
"""
function SNESTestJacobian(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESTestJacobian(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESComputeJacobian(petsclib::PetscLibType,snes::PetscSNES, X::PetscVec, A::PetscMat, B::PetscMat) 
Computes the Jacobian matrix that has been set with `SNESSetJacobian()`.

Collective

Input Parameters:
- `snes` - the `SNES` context
- `X`    - input vector

Output Parameters:
- `A` - Jacobian matrix
- `B` - optional matrix for building the preconditioner, usually the same as `A`

Options Database Keys:
- `-snes_lag_preconditioner <lag>`           - how often to rebuild preconditioner
- `-snes_lag_jacobian <lag>`                 - how often to rebuild Jacobian
- `-snes_test_jacobian <optional threshold>` - compare the user provided Jacobian with one compute via finite differences to check for errors.  If a threshold is given, display only those entries whose difference is greater than the threshold.
- `-snes_test_jacobian_view`                 - display the user provided Jacobian, the finite difference Jacobian and the difference between them to help users detect the location of errors in the user provided Jacobian
- `-snes_compare_explicit`                   - Compare the computed Jacobian to the finite difference Jacobian and output the differences
- `-snes_compare_explicit_draw`              - Compare the computed Jacobian to the finite difference Jacobian and draw the result
- `-snes_compare_explicit_contour`           - Compare the computed Jacobian to the finite difference Jacobian and draw a contour plot with the result
- `-snes_compare_operator`                   - Make the comparison options above use the operator instead of the matrix used to construct the preconditioner
- `-snes_compare_coloring`                   - Compute the finite difference Jacobian using coloring and display norms of difference
- `-snes_compare_coloring_display`           - Compute the finite difference Jacobian using coloring and display verbose differences
- `-snes_compare_coloring_threshold`         - Display only those matrix entries that differ by more than a given threshold
- `-snes_compare_coloring_threshold_atol`    - Absolute tolerance for difference in matrix entries to be displayed by `-snes_compare_coloring_threshold`
- `-snes_compare_coloring_threshold_rtol`    - Relative tolerance for difference in matrix entries to be displayed by `-snes_compare_coloring_threshold`
- `-snes_compare_coloring_draw`              - Compute the finite difference Jacobian using coloring and draw differences
- `-snes_compare_coloring_draw_contour`      - Compute the finite difference Jacobian using coloring and show contours of matrices and differences

Level: developer

-seealso: [](ch_snes), `SNESSetJacobian()`, `KSPSetOperators()`, `MatStructure`, `SNESSetLagPreconditioner()`, `SNESSetLagJacobian()`

# External Links
$(_doc_external("Snes/SNESComputeJacobian"))
"""
function SNESComputeJacobian(petsclib::PetscLibType, snes::PetscSNES, X::PetscVec, A::PetscMat, B::PetscMat) end

@for_petsc function SNESComputeJacobian(petsclib::$UnionPetscLib, snes::PetscSNES, X::PetscVec, A::PetscMat, B::PetscMat )

    @chk ccall(
               (:SNESComputeJacobian, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CMat, CMat),
               snes, X, A, B,
              )


	return nothing
end 

"""
	SNESSetJacobian(petsclib::PetscLibType,snes::PetscSNES, Amat::PetscMat, Pmat::PetscMat, J::SNESJacobianFn, ctx::Cvoid) 
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

-seealso: [](ch_snes), `SNES`, `KSPSetOperators()`, `SNESSetFunction()`, `MatMFFDComputeJacobian()`, `SNESComputeJacobianDefaultColor()`, `MatStructure`,
`SNESSetPicard()`, `SNESJacobianFn`, `SNESFunctionFn`

# External Links
$(_doc_external("Snes/SNESSetJacobian"))
"""
function SNESSetJacobian(petsclib::PetscLibType, snes::PetscSNES, Amat::PetscMat, Pmat::PetscMat, J::SNESJacobianFn, ctx::Cvoid) end

@for_petsc function SNESSetJacobian(petsclib::$UnionPetscLib, snes::PetscSNES, Amat::PetscMat, Pmat::PetscMat, J::SNESJacobianFn, ctx::Cvoid )

    @chk ccall(
               (:SNESSetJacobian, $petsc_library),
               PetscErrorCode,
               (CSNES, CMat, CMat, Ptr{SNESJacobianFn}, Ptr{Cvoid}),
               snes, Amat, Pmat, J, ctx,
              )


	return nothing
end 

"""
	SNESGetJacobian(petsclib::PetscLibType,snes::PetscSNES, Amat::PetscMat, Pmat::PetscMat, J::SNESJacobianFn, ctx::Cvoid) 
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

-seealso: [](ch_snes), `SNES`, `Mat`, `SNESSetJacobian()`, `SNESComputeJacobian()`, `SNESJacobianFn`, `SNESGetFunction()`

# External Links
$(_doc_external("Snes/SNESGetJacobian"))
"""
function SNESGetJacobian(petsclib::PetscLibType, snes::PetscSNES, Amat::PetscMat, Pmat::PetscMat, J::SNESJacobianFn, ctx::Cvoid) end

@for_petsc function SNESGetJacobian(petsclib::$UnionPetscLib, snes::PetscSNES, Amat::PetscMat, Pmat::PetscMat, J::SNESJacobianFn, ctx::Cvoid )
	Amat_ = Ref(Amat.ptr)
	Pmat_ = Ref(Pmat.ptr)

    @chk ccall(
               (:SNESGetJacobian, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CMat}, Ptr{CMat}, SNESJacobianFn, Cvoid),
               snes, Amat_, Pmat_, J, ctx,
              )

	Amat.ptr = C_NULL
	Pmat.ptr = C_NULL

	return nothing
end 

"""
	SNESSetUp(petsclib::PetscLibType,snes::PetscSNES) 
Sets up the internal data structures for the later use
of a nonlinear solver `SNESSolve()`.

Collective

Input Parameter:
- `snes` - the `SNES` context

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESCreate()`, `SNESSolve()`, `SNESDestroy()`, `SNESSetFromOptions()`

# External Links
$(_doc_external("Snes/SNESSetUp"))
"""
function SNESSetUp(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESSetUp(petsclib::$UnionPetscLib, snes::PetscSNES )

    @chk ccall(
               (:SNESSetUp, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	SNESReset(petsclib::PetscLibType,snes::PetscSNES) 
Resets a `SNES` context to the state it was in before `SNESSetUp()` was called and removes any allocated `Vec` and `Mat` from its data structures

Collective

Input Parameter:
- `snes` - the nonlinear iterative solver context obtained from `SNESCreate()`

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESDestroy()`, `SNESCreate()`, `SNESSetUp()`, `SNESSolve()`

# External Links
$(_doc_external("Snes/SNESReset"))
"""
function SNESReset(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESReset(petsclib::$UnionPetscLib, snes::PetscSNES )

    @chk ccall(
               (:SNESReset, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	SNESConvergedReasonViewCancel(petsclib::PetscLibType,snes::PetscSNES) 
Clears all the reason view functions for a `SNES` object provided with `SNESConvergedReasonViewSet()` also
removes the default viewer.

Collective

Input Parameter:
- `snes` - the nonlinear iterative solver context obtained from `SNESCreate()`

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESCreate()`, `SNESDestroy()`, `SNESReset()`, `SNESConvergedReasonViewSet()`

# External Links
$(_doc_external("Snes/SNESConvergedReasonViewCancel"))
"""
function SNESConvergedReasonViewCancel(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESConvergedReasonViewCancel(petsclib::$UnionPetscLib, snes::PetscSNES )

    @chk ccall(
               (:SNESConvergedReasonViewCancel, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	SNESDestroy(petsclib::PetscLibType,snes::PetscSNES) 
Destroys the nonlinear solver context that was created
with `SNESCreate()`.

Collective

Input Parameter:
- `snes` - the `SNES` context

Level: beginner

-seealso: [](ch_snes), `SNES`, `SNESCreate()`, `SNESSolve()`

# External Links
$(_doc_external("Snes/SNESDestroy"))
"""
function SNESDestroy(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESDestroy(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESSetLagPreconditioner(petsclib::PetscLibType,snes::PetscSNES, lag::PetscInt) 
Sets when the preconditioner is rebuilt in the nonlinear solve `SNESSolve()`.

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `lag`  - 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time
the Jacobian is built etc. -2 indicates rebuild preconditioner at next chance but then never rebuild after that

Options Database Keys:
- `-snes_lag_jacobian_persists <true,false>`       - sets the persistence through multiple `SNESSolve()`
- `-snes_lag_jacobian <-2,1,2,...>`                - sets the lag
- `-snes_lag_preconditioner_persists <true,false>` - sets the persistence through multiple `SNESSolve()`
- `-snes_lag_preconditioner <-2,1,2,...>`          - sets the lag

Level: intermediate

-seealso: [](ch_snes), `SNESGetLagPreconditioner()`, `SNESSetLagJacobian()`, `SNESGetLagJacobian()`, `SNESSetLagPreconditionerPersists()`,
`SNESSetLagJacobianPersists()`, `SNES`, `SNESSolve()`

# External Links
$(_doc_external("Snes/SNESSetLagPreconditioner"))
"""
function SNESSetLagPreconditioner(petsclib::PetscLibType, snes::PetscSNES, lag::PetscInt) end

@for_petsc function SNESSetLagPreconditioner(petsclib::$UnionPetscLib, snes::PetscSNES, lag::$PetscInt )

    @chk ccall(
               (:SNESSetLagPreconditioner, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt),
               snes, lag,
              )


	return nothing
end 

"""
	SNESSetGridSequence(petsclib::PetscLibType,snes::PetscSNES, steps::PetscInt) 
sets the number of steps of grid sequencing that `SNES` will do

Logically Collective

Input Parameters:
- `snes`  - the `SNES` context
- `steps` - the number of refinements to do, defaults to 0

Options Database Key:
- `-snes_grid_sequence <steps>` - Use grid sequencing to generate initial guess

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESGetLagPreconditioner()`, `SNESSetLagJacobian()`, `SNESGetLagJacobian()`, `SNESGetGridSequence()`,
`SNESSetDM()`, `SNESSolve()`

# External Links
$(_doc_external("Snes/SNESSetGridSequence"))
"""
function SNESSetGridSequence(petsclib::PetscLibType, snes::PetscSNES, steps::PetscInt) end

@for_petsc function SNESSetGridSequence(petsclib::$UnionPetscLib, snes::PetscSNES, steps::$PetscInt )

    @chk ccall(
               (:SNESSetGridSequence, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt),
               snes, steps,
              )


	return nothing
end 

"""
	steps::PetscInt = SNESGetGridSequence(petsclib::PetscLibType,snes::PetscSNES) 
gets the number of steps of grid sequencing that `SNES` will do

Logically Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `steps` - the number of refinements to do, defaults to 0

Level: intermediate

-seealso: [](ch_snes), `SNESGetLagPreconditioner()`, `SNESSetLagJacobian()`, `SNESGetLagJacobian()`, `SNESSetGridSequence()`

# External Links
$(_doc_external("Snes/SNESGetGridSequence"))
"""
function SNESGetGridSequence(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESGetGridSequence(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	lag::PetscInt = SNESGetLagPreconditioner(petsclib::PetscLibType,snes::PetscSNES) 
Return how often the preconditioner is rebuilt

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `lag` - -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time
the Jacobian is built etc. -2 indicates rebuild preconditioner at next chance but then never rebuild after that

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESSetLagPreconditioner()`, `SNESSetLagJacobianPersists()`, `SNESSetLagPreconditionerPersists()`

# External Links
$(_doc_external("Snes/SNESGetLagPreconditioner"))
"""
function SNESGetLagPreconditioner(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESGetLagPreconditioner(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESSetLagJacobian(petsclib::PetscLibType,snes::PetscSNES, lag::PetscInt) 
Set when the Jacobian is rebuilt in the nonlinear solve. See `SNESSetLagPreconditioner()` for determining how
often the preconditioner is rebuilt.

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `lag`  - -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time
the Jacobian is built etc. -2 means rebuild at next chance but then never again

Options Database Keys:
- `-snes_lag_jacobian_persists <true,false>`       - sets the persistence through multiple SNES solves
- `-snes_lag_jacobian <-2,1,2,...>`                - sets the lag
- `-snes_lag_preconditioner_persists <true,false>` - sets the persistence through multiple SNES solves
- `-snes_lag_preconditioner <-2,1,2,...>`          - sets the lag.

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESGetLagPreconditioner()`, `SNESSetLagPreconditioner()`, `SNESGetLagJacobianPersists()`, `SNESSetLagPreconditionerPersists()`

# External Links
$(_doc_external("Snes/SNESSetLagJacobian"))
"""
function SNESSetLagJacobian(petsclib::PetscLibType, snes::PetscSNES, lag::PetscInt) end

@for_petsc function SNESSetLagJacobian(petsclib::$UnionPetscLib, snes::PetscSNES, lag::$PetscInt )

    @chk ccall(
               (:SNESSetLagJacobian, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt),
               snes, lag,
              )


	return nothing
end 

"""
	lag::PetscInt = SNESGetLagJacobian(petsclib::PetscLibType,snes::PetscSNES) 
Get how often the Jacobian is rebuilt. See `SNESGetLagPreconditioner()` to determine when the preconditioner is rebuilt

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `lag` - -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time
the Jacobian is built etc.

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESSetLagJacobian()`, `SNESSetLagPreconditioner()`, `SNESGetLagPreconditioner()`, `SNESSetLagJacobianPersists()`, `SNESSetLagPreconditionerPersists()`


# External Links
$(_doc_external("Snes/SNESGetLagJacobian"))
"""
function SNESGetLagJacobian(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESGetLagJacobian(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESSetLagJacobianPersists(petsclib::PetscLibType,snes::PetscSNES, flg::PetscBool) 
Set whether or not the Jacobian lagging persists through multiple nonlinear solves

Logically collective

Input Parameters:
- `snes` - the `SNES` context
- `flg`  - jacobian lagging persists if true

Options Database Keys:
- `-snes_lag_jacobian_persists <true,false>`       - sets the persistence through multiple SNES solves
- `-snes_lag_jacobian <-2,1,2,...>`                - sets the lag
- `-snes_lag_preconditioner_persists <true,false>` - sets the persistence through multiple SNES solves
- `-snes_lag_preconditioner <-2,1,2,...>`          - sets the lag

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESSetLagPreconditionerPersists()`, `SNESSetLagJacobian()`, `SNESGetLagJacobian()`, `SNESGetNPC()`

# External Links
$(_doc_external("Snes/SNESSetLagJacobianPersists"))
"""
function SNESSetLagJacobianPersists(petsclib::PetscLibType, snes::PetscSNES, flg::PetscBool) end

@for_petsc function SNESSetLagJacobianPersists(petsclib::$UnionPetscLib, snes::PetscSNES, flg::PetscBool )

    @chk ccall(
               (:SNESSetLagJacobianPersists, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, flg,
              )


	return nothing
end 

"""
	SNESSetLagPreconditionerPersists(petsclib::PetscLibType,snes::PetscSNES, flg::PetscBool) 
Set whether or not the preconditioner lagging persists through multiple nonlinear solves

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `flg`  - preconditioner lagging persists if true

Options Database Keys:
- `-snes_lag_jacobian_persists <true,false>`       - sets the persistence through multiple SNES solves
- `-snes_lag_jacobian <-2,1,2,...>`                - sets the lag
- `-snes_lag_preconditioner_persists <true,false>` - sets the persistence through multiple SNES solves
- `-snes_lag_preconditioner <-2,1,2,...>`          - sets the lag

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESSetLagJacobianPersists()`, `SNESSetLagJacobian()`, `SNESGetLagJacobian()`, `SNESGetNPC()`, `SNESSetLagPreconditioner()`

# External Links
$(_doc_external("Snes/SNESSetLagPreconditionerPersists"))
"""
function SNESSetLagPreconditionerPersists(petsclib::PetscLibType, snes::PetscSNES, flg::PetscBool) end

@for_petsc function SNESSetLagPreconditionerPersists(petsclib::$UnionPetscLib, snes::PetscSNES, flg::PetscBool )

    @chk ccall(
               (:SNESSetLagPreconditionerPersists, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, flg,
              )


	return nothing
end 

"""
	SNESSetForceIteration(petsclib::PetscLibType,snes::PetscSNES, force::PetscBool) 
force `SNESSolve()` to take at least one iteration regardless of the initial residual norm

Logically Collective

Input Parameters:
- `snes`  - the `SNES` context
- `force` - `PETSC_TRUE` require at least one iteration

Options Database Key:
- `-snes_force_iteration <force>` - Sets forcing an iteration

Level: intermediate

-seealso: [](ch_snes), `SNES`, `TS`, `SNESSetDivergenceTolerance()`

# External Links
$(_doc_external("Snes/SNESSetForceIteration"))
"""
function SNESSetForceIteration(petsclib::PetscLibType, snes::PetscSNES, force::PetscBool) end

@for_petsc function SNESSetForceIteration(petsclib::$UnionPetscLib, snes::PetscSNES, force::PetscBool )

    @chk ccall(
               (:SNESSetForceIteration, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, force,
              )


	return nothing
end 

"""
	force::PetscBool = SNESGetForceIteration(petsclib::PetscLibType,snes::PetscSNES) 
Check whether or not `SNESSolve()` take at least one iteration regardless of the initial residual norm

Logically Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `force` - `PETSC_TRUE` requires at least one iteration.

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESSetForceIteration()`, `SNESSetDivergenceTolerance()`

# External Links
$(_doc_external("Snes/SNESGetForceIteration"))
"""
function SNESGetForceIteration(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESGetForceIteration(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESSetTolerances(petsclib::PetscLibType,snes::PetscSNES, abstol::PetscReal, rtol::PetscReal, stol::PetscReal, maxit::PetscInt, maxf::PetscInt) 
Sets various parameters used in `SNES` convergence tests.

Logically Collective

Input Parameters:
- `snes`   - the `SNES` context
- `abstol` - the absolute convergence tolerance,  F(x^n) ≤ abstol 
- `rtol`   - the relative convergence tolerance,  F(x^n) ≤ reltol * F(x^0) 
- `stol`   - convergence tolerance in terms of the norm of the change in the solution between steps,  || delta x || < stol*|| x ||
- `maxit`  - the maximum number of iterations allowed in the solver, default 50.
- `maxf`   - the maximum number of function evaluations allowed in the solver (use `PETSC_UNLIMITED` indicates no limit), default 10,000

Options Database Keys:
- `-snes_atol <abstol>`    - Sets `abstol`
- `-snes_rtol <rtol>`      - Sets `rtol`
- `-snes_stol <stol>`      - Sets `stol`
- `-snes_max_it <maxit>`   - Sets `maxit`
- `-snes_max_funcs <maxf>` - Sets `maxf` (use `unlimited` to have no maximum)

Level: intermediate

-seealso: [](ch_snes), `SNESSolve()`, `SNES`, `SNESSetDivergenceTolerance()`, `SNESSetForceIteration()`

# External Links
$(_doc_external("Snes/SNESSetTolerances"))
"""
function SNESSetTolerances(petsclib::PetscLibType, snes::PetscSNES, abstol::PetscReal, rtol::PetscReal, stol::PetscReal, maxit::PetscInt, maxf::PetscInt) end

@for_petsc function SNESSetTolerances(petsclib::$UnionPetscLib, snes::PetscSNES, abstol::$PetscReal, rtol::$PetscReal, stol::$PetscReal, maxit::$PetscInt, maxf::$PetscInt )

    @chk ccall(
               (:SNESSetTolerances, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscReal, $PetscReal, $PetscReal, $PetscInt, $PetscInt),
               snes, abstol, rtol, stol, maxit, maxf,
              )


	return nothing
end 

"""
	SNESSetDivergenceTolerance(petsclib::PetscLibType,snes::PetscSNES, divtol::PetscReal) 
Sets the divergence tolerance used for the `SNES` divergence test.

Logically Collective

Input Parameters:
- `snes`   - the `SNES` context
- `divtol` - the divergence tolerance. Use `PETSC_UNLIMITED` to deactivate the test. If the residual norm  F(x^n) ≥ divtol * F(x^0)  the solver
is stopped due to divergence.

Options Database Key:
- `-snes_divergence_tolerance <divtol>` - Sets `divtol`

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESSolve()`, `SNESSetTolerances()`, `SNESGetDivergenceTolerance()`

# External Links
$(_doc_external("Snes/SNESSetDivergenceTolerance"))
"""
function SNESSetDivergenceTolerance(petsclib::PetscLibType, snes::PetscSNES, divtol::PetscReal) end

@for_petsc function SNESSetDivergenceTolerance(petsclib::$UnionPetscLib, snes::PetscSNES, divtol::$PetscReal )

    @chk ccall(
               (:SNESSetDivergenceTolerance, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscReal),
               snes, divtol,
              )


	return nothing
end 

"""
	atol::PetscReal,rtol::PetscReal,stol::PetscReal,maxit::PetscInt,maxf::PetscInt = SNESGetTolerances(petsclib::PetscLibType,snes::PetscSNES) 
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

-seealso: [](ch_snes), `SNES`, `SNESSetTolerances()`

# External Links
$(_doc_external("Snes/SNESGetTolerances"))
"""
function SNESGetTolerances(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESGetTolerances(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESGetDivergenceTolerance(petsclib::PetscLibType,snes::PetscSNES, divtol::PetscReal) 
Gets divergence tolerance used in divergence test.

Not Collective

Input Parameters:
- `snes`   - the `SNES` context
- `divtol` - divergence tolerance

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESSetDivergenceTolerance()`

# External Links
$(_doc_external("Snes/SNESGetDivergenceTolerance"))
"""
function SNESGetDivergenceTolerance(petsclib::PetscLibType, snes::PetscSNES, divtol::PetscReal) end

@for_petsc function SNESGetDivergenceTolerance(petsclib::$UnionPetscLib, snes::PetscSNES, divtol::$PetscReal )

    @chk ccall(
               (:SNESGetDivergenceTolerance, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscReal}),
               snes, divtol,
              )


	return nothing
end 

"""
	SNESMonitorLGRange(petsclib::PetscLibType,snes::PetscSNES, n::PetscInt, rnorm::PetscReal, monctx::Cvoid) 

# External Links
$(_doc_external("Snes/SNESMonitorLGRange"))
"""
function SNESMonitorLGRange(petsclib::PetscLibType, snes::PetscSNES, n::PetscInt, rnorm::PetscReal, monctx::Cvoid) end

@for_petsc function SNESMonitorLGRange(petsclib::$UnionPetscLib, snes::PetscSNES, n::$PetscInt, rnorm::$PetscReal, monctx::Cvoid )

    @chk ccall(
               (:SNESMonitorLGRange, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, Ptr{Cvoid}),
               snes, n, rnorm, monctx,
              )


	return nothing
end 

"""
	SNESConverged(petsclib::PetscLibType,snes::PetscSNES, it::PetscInt, xnorm::PetscReal, snorm::PetscReal, fnorm::PetscReal) 
Run the convergence test and update the `SNESConvergedReason`.

Collective

Input Parameters:
- `snes`  - the `SNES` context
- `it`    - current iteration
- `xnorm` - 2-norm of current iterate
- `snorm` - 2-norm of current step
- `fnorm` - 2-norm of function

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESSolve`, `SNESSetConvergenceTest()`

# External Links
$(_doc_external("Snes/SNESConverged"))
"""
function SNESConverged(petsclib::PetscLibType, snes::PetscSNES, it::PetscInt, xnorm::PetscReal, snorm::PetscReal, fnorm::PetscReal) end

@for_petsc function SNESConverged(petsclib::$UnionPetscLib, snes::PetscSNES, it::$PetscInt, xnorm::$PetscReal, snorm::$PetscReal, fnorm::$PetscReal )

    @chk ccall(
               (:SNESConverged, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, $PetscReal, $PetscReal),
               snes, it, xnorm, snorm, fnorm,
              )


	return nothing
end 

"""
	SNESMonitor(petsclib::PetscLibType,snes::PetscSNES, iter::PetscInt, rnorm::PetscReal) 
runs any `SNES` monitor routines provided with `SNESMonitor()` or the options database

Collective

Input Parameters:
- `snes`  - nonlinear solver context obtained from `SNESCreate()`
- `iter`  - current iteration number
- `rnorm` - current relative norm of the residual

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESMonitorSet()`

# External Links
$(_doc_external("Snes/SNESMonitor"))
"""
function SNESMonitor(petsclib::PetscLibType, snes::PetscSNES, iter::PetscInt, rnorm::PetscReal) end

@for_petsc function SNESMonitor(petsclib::$UnionPetscLib, snes::PetscSNES, iter::$PetscInt, rnorm::$PetscReal )

    @chk ccall(
               (:SNESMonitor, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal),
               snes, iter, rnorm,
              )


	return nothing
end 

"""
	SNESMonitorSet(petsclib::PetscLibType,snes::PetscSNES, f::external, mctx::Cvoid, monitordestroy::PetscCtxDestroyFn) 
Sets an ADDITIONAL function that is to be used at every
iteration of the `SNES` nonlinear solver to display the iteration's
progress.

Logically Collective

Input Parameters:
- `snes`           - the `SNES` context
- `f`              - the monitor function,  for the calling sequence see `SNESMonitorFunction`
- `mctx`           - [optional] user-defined context for private data for the monitor routine (use `NULL` if no context is desired)
- `monitordestroy` - [optional] routine that frees monitor context (may be `NULL`), see `PetscCtxDestroyFn` for the calling sequence

Options Database Keys:
- `-snes_monitor`               - sets `SNESMonitorDefault()`
- `-snes_monitor draw::draw_lg` - sets line graph monitor,
- `-snes_monitor_cancel`        - cancels all monitors that have been hardwired into a code by calls to `SNESMonitorSet()`, but does not cancel those set via
the options database.

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESSolve()`, `SNESMonitorDefault()`, `SNESMonitorCancel()`, `SNESMonitorFunction`, `PetscCtxDestroyFn`

# External Links
$(_doc_external("Snes/SNESMonitorSet"))
"""
function SNESMonitorSet(petsclib::PetscLibType, snes::PetscSNES, f::external, mctx::Cvoid, monitordestroy::PetscCtxDestroyFn) end

@for_petsc function SNESMonitorSet(petsclib::$UnionPetscLib, snes::PetscSNES, f::external, mctx::Cvoid, monitordestroy::PetscCtxDestroyFn )

    @chk ccall(
               (:SNESMonitorSet, $petsc_library),
               PetscErrorCode,
               (CSNES, external, Ptr{Cvoid}, Ptr{PetscCtxDestroyFn}),
               snes, f, mctx, monitordestroy,
              )


	return nothing
end 

"""
	SNESMonitorCancel(petsclib::PetscLibType,snes::PetscSNES) 
Clears all the monitor functions for a `SNES` object.

Logically Collective

Input Parameter:
- `snes` - the `SNES` context

Options Database Key:
- `-snes_monitor_cancel` - cancels all monitors that have been hardwired
into a code by calls to `SNESMonitorSet()`, but does not cancel those
set via the options database

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESMonitorDefault()`, `SNESMonitorSet()`

# External Links
$(_doc_external("Snes/SNESMonitorCancel"))
"""
function SNESMonitorCancel(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESMonitorCancel(petsclib::$UnionPetscLib, snes::PetscSNES )

    @chk ccall(
               (:SNESMonitorCancel, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	SNESSetConvergenceTest(petsclib::PetscLibType,snes::PetscSNES, SNESConvergenceTestFunction::external, cctx::Cvoid, destroy::external) 
Sets the function that is to be used
to test for convergence of the nonlinear iterative solution.

Logically Collective

Input Parameters:
- `snes`                        - the `SNES` context
- `SNESConvergenceTestFunction` - routine to test for convergence
- `cctx`                        - [optional] context for private data for the convergence routine  (may be `NULL`)
- `destroy`                     - [optional] destructor for the context (may be `NULL`; `PETSC_NULL_FUNCTION` in Fortran)

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESConvergedDefault()`, `SNESConvergedSkip()`, `SNESConvergenceTestFunction`

# External Links
$(_doc_external("Snes/SNESSetConvergenceTest"))
"""
function SNESSetConvergenceTest(petsclib::PetscLibType, snes::PetscSNES, SNESConvergenceTestFunction::external, cctx::Cvoid, destroy::external) end

@for_petsc function SNESSetConvergenceTest(petsclib::$UnionPetscLib, snes::PetscSNES, SNESConvergenceTestFunction::external, cctx::Cvoid, destroy::external )

    @chk ccall(
               (:SNESSetConvergenceTest, $petsc_library),
               PetscErrorCode,
               (CSNES, external, Ptr{Cvoid}, external),
               snes, SNESConvergenceTestFunction, cctx, destroy,
              )


	return nothing
end 

"""
	SNESGetConvergedReason(petsclib::PetscLibType,snes::PetscSNES, reason::SNESConvergedReason) 
Gets the reason the `SNES` iteration was stopped, which may be due to convergence, divergence, or stagnation

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `reason` - negative value indicates diverged, positive value converged, see `SNESConvergedReason` for the individual convergence tests for complete lists

Options Database Key:
- `-snes_converged_reason` - prints the reason to standard out

Level: intermediate

-seealso: [](ch_snes), `SNESSolve()`, `SNESSetConvergenceTest()`, `SNESSetConvergedReason()`, `SNESConvergedReason`, `SNESGetConvergedReasonString()`

# External Links
$(_doc_external("Snes/SNESGetConvergedReason"))
"""
function SNESGetConvergedReason(petsclib::PetscLibType, snes::PetscSNES, reason::SNESConvergedReason) end

@for_petsc function SNESGetConvergedReason(petsclib::$UnionPetscLib, snes::PetscSNES, reason::SNESConvergedReason )

    @chk ccall(
               (:SNESGetConvergedReason, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{SNESConvergedReason}),
               snes, reason,
              )


	return nothing
end 

"""
	SNESGetConvergedReasonString(petsclib::PetscLibType,snes::PetscSNES, strreason::Cchar) 
Return a human readable string for `SNESConvergedReason`

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `strreason` - a human readable string that describes `SNES` converged reason

Level: beginner

-seealso: [](ch_snes), `SNES`, `SNESGetConvergedReason()`

# External Links
$(_doc_external("Snes/SNESGetConvergedReasonString"))
"""
function SNESGetConvergedReasonString(petsclib::PetscLibType, snes::PetscSNES, strreason::Cchar) end

@for_petsc function SNESGetConvergedReasonString(petsclib::$UnionPetscLib, snes::PetscSNES, strreason::Cchar )

    @chk ccall(
               (:SNESGetConvergedReasonString, $petsc_library),
               PetscErrorCode,
               (CSNES, Cchar),
               snes, strreason,
              )


	return nothing
end 

"""
	SNESSetConvergedReason(petsclib::PetscLibType,snes::PetscSNES, reason::SNESConvergedReason) 
Sets the reason the `SNES` iteration was stopped.

Not Collective

Input Parameters:
- `snes`   - the `SNES` context
- `reason` - negative value indicates diverged, positive value converged, see `SNESConvergedReason` or the
manual pages for the individual convergence tests for complete lists

Level: developer

-seealso: [](ch_snes), `SNESGetConvergedReason()`, `SNESSetConvergenceTest()`, `SNESConvergedReason`

# External Links
$(_doc_external("Snes/SNESSetConvergedReason"))
"""
function SNESSetConvergedReason(petsclib::PetscLibType, snes::PetscSNES, reason::SNESConvergedReason) end

@for_petsc function SNESSetConvergedReason(petsclib::$UnionPetscLib, snes::PetscSNES, reason::SNESConvergedReason )

    @chk ccall(
               (:SNESSetConvergedReason, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESConvergedReason),
               snes, reason,
              )


	return nothing
end 

"""
	SNESSetConvergenceHistory(petsclib::PetscLibType,snes::PetscSNES, a::Vector{PetscReal}, its::Vector{PetscInt}, na::PetscInt, reset::PetscBool) 
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

-seealso: [](ch_snes), `SNES`, `SNESSolve()`, `SNESGetConvergenceHistory()`

# External Links
$(_doc_external("Snes/SNESSetConvergenceHistory"))
"""
function SNESSetConvergenceHistory(petsclib::PetscLibType, snes::PetscSNES, a::Vector{PetscReal}, its::Vector{PetscInt}, na::PetscInt, reset::PetscBool) end

@for_petsc function SNESSetConvergenceHistory(petsclib::$UnionPetscLib, snes::PetscSNES, a::Vector{$PetscReal}, its::Vector{$PetscInt}, na::$PetscInt, reset::PetscBool )

    @chk ccall(
               (:SNESSetConvergenceHistory, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscReal}, Ptr{$PetscInt}, $PetscInt, PetscBool),
               snes, a, its, na, reset,
              )


	return nothing
end 

"""
	a::Vector{PetscReal},its::Vector{PetscInt},na::PetscInt = SNESGetConvergenceHistory(petsclib::PetscLibType,snes::PetscSNES) 
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

-seealso: [](ch_snes), `SNES`, `SNESSolve()`, `SNESSetConvergenceHistory()`

# External Links
$(_doc_external("Snes/SNESGetConvergenceHistory"))
"""
function SNESGetConvergenceHistory(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESGetConvergenceHistory(petsclib::$UnionPetscLib, snes::PetscSNES )
	a_ = Ref{Ptr{$PetscReal}}()
	its_ = Ref{Ptr{$PetscInt}}()
	na_ = Ref{$PetscInt}()

    @chk ccall(
               (:SNESGetConvergenceHistory, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Ptr{$PetscReal}}, Ptr{Ptr{$PetscInt}}, Ptr{$PetscInt}),
               snes, a_, its_, na_,
              )

	a = unsafe_wrap(Array, a_[], VecGetLocalSize(petsclib, x); own = false)
	its = unsafe_wrap(Array, its_[], VecGetLocalSize(petsclib, x); own = false)
	na = na_[]

	return a,its,na
end 

"""
	SNESSetUpdate(petsclib::PetscLibType,snes::PetscSNES, func::SNESUpdateFn) 
Sets the general
at the beginning of every iteration of the nonlinear solve. Specifically
it is called just before the Jacobian is "evaluated" and after the function
evaluation.

Logically Collective

Input Parameters:
- `snes` - The nonlinear solver context
- `func` - The update function; for calling sequence see `SNESUpdateFn`

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESSolve()`, `SNESSetJacobian()`, `SNESLineSearchSetPreCheck()`, `SNESLineSearchSetPostCheck()`, `SNESNewtonTRSetPreCheck()`, `SNESNewtonTRSetPostCheck()`,
`SNESMonitorSet()`

# External Links
$(_doc_external("Snes/SNESSetUpdate"))
"""
function SNESSetUpdate(petsclib::PetscLibType, snes::PetscSNES, func::SNESUpdateFn) end

@for_petsc function SNESSetUpdate(petsclib::$UnionPetscLib, snes::PetscSNES, func::SNESUpdateFn )

    @chk ccall(
               (:SNESSetUpdate, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{SNESUpdateFn}),
               snes, func,
              )


	return nothing
end 

"""
	SNESConvergedReasonView(petsclib::PetscLibType,snes::PetscSNES, viewer::PetscViewer) 
Displays the reason a `SNES` solve converged or diverged to a viewer

Collective

Input Parameters:
- `snes`   - iterative context obtained from `SNESCreate()`
- `viewer` - the viewer to display the reason

Options Database Keys:
- `-snes_converged_reason`          - print reason for converged or diverged, also prints number of iterations
- `-snes_converged_reason ::failed` - only print reason and number of iterations when diverged

Level: beginner

-seealso: [](ch_snes), `SNESConvergedReason`, `PetscViewer`, `SNES`,
`SNESCreate()`, `SNESSetUp()`, `SNESDestroy()`, `SNESSetTolerances()`, `SNESConvergedDefault()`, `SNESGetConvergedReason()`,
`SNESConvergedReasonViewFromOptions()`,
`PetscViewerPushFormat()`, `PetscViewerPopFormat()`

# External Links
$(_doc_external("Snes/SNESConvergedReasonView"))
"""
function SNESConvergedReasonView(petsclib::PetscLibType, snes::PetscSNES, viewer::PetscViewer) end

@for_petsc function SNESConvergedReasonView(petsclib::$UnionPetscLib, snes::PetscSNES, viewer::PetscViewer )

    @chk ccall(
               (:SNESConvergedReasonView, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscViewer),
               snes, viewer,
              )


	return nothing
end 

"""
	SNESConvergedReasonViewSet(petsclib::PetscLibType,snes::PetscSNES, f::external, vctx::Cvoid, reasonviewdestroy::PetscCtxDestroyFn) 
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

-seealso: [](ch_snes), `SNES`, `SNESSolve()`, `SNESConvergedReason`, `SNESGetConvergedReason()`, `SNESConvergedReasonView()`, `SNESConvergedReasonViewCancel()`,
`PetscCtxDestroyFn`

# External Links
$(_doc_external("Snes/SNESConvergedReasonViewSet"))
"""
function SNESConvergedReasonViewSet(petsclib::PetscLibType, snes::PetscSNES, f::external, vctx::Cvoid, reasonviewdestroy::PetscCtxDestroyFn) end

@for_petsc function SNESConvergedReasonViewSet(petsclib::$UnionPetscLib, snes::PetscSNES, f::external, vctx::Cvoid, reasonviewdestroy::PetscCtxDestroyFn )

    @chk ccall(
               (:SNESConvergedReasonViewSet, $petsc_library),
               PetscErrorCode,
               (CSNES, external, Ptr{Cvoid}, Ptr{PetscCtxDestroyFn}),
               snes, f, vctx, reasonviewdestroy,
              )


	return nothing
end 

"""
	SNESConvergedReasonViewFromOptions(petsclib::PetscLibType,snes::PetscSNES) 
Processes command line options to determine if/how a `SNESConvergedReason` is to be viewed at the end of `SNESSolve()`
All the user-provided viewer routines set with `SNESConvergedReasonViewSet()` will be called, if they exist.

Collective

Input Parameter:
- `snes` - the `SNES` object

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESConvergedReason`, `SNESConvergedReasonViewSet()`, `SNESCreate()`, `SNESSetUp()`, `SNESDestroy()`,
`SNESSetTolerances()`, `SNESConvergedDefault()`, `SNESGetConvergedReason()`, `SNESConvergedReasonView()`

# External Links
$(_doc_external("Snes/SNESConvergedReasonViewFromOptions"))
"""
function SNESConvergedReasonViewFromOptions(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESConvergedReasonViewFromOptions(petsclib::$UnionPetscLib, snes::PetscSNES )

    @chk ccall(
               (:SNESConvergedReasonViewFromOptions, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	SNESSolve(petsclib::PetscLibType,snes::PetscSNES, b::PetscVec, x::PetscVec) 
Solves a nonlinear system F(x) = b  associated with a `SNES` object

Collective

Input Parameters:
- `snes` - the `SNES` context
- `b`    - the constant part of the equation F(x) = b, or `NULL` to use zero.
- `x`    - the solution vector.

Level: beginner

-seealso: [](ch_snes), `SNES`, `SNESCreate()`, `SNESDestroy()`, `SNESSetFunction()`, `SNESSetJacobian()`, `SNESSetGridSequence()`, `SNESGetSolution()`,
`SNESNewtonTRSetPreCheck()`, `SNESNewtonTRGetPreCheck()`, `SNESNewtonTRSetPostCheck()`, `SNESNewtonTRGetPostCheck()`,
`SNESLineSearchSetPostCheck()`, `SNESLineSearchGetPostCheck()`, `SNESLineSearchSetPreCheck()`, `SNESLineSearchGetPreCheck()`

# External Links
$(_doc_external("Snes/SNESSolve"))
"""
function SNESSolve(petsclib::PetscLibType, snes::PetscSNES, b::PetscVec, x::PetscVec) end

@for_petsc function SNESSolve(petsclib::$UnionPetscLib, snes::PetscSNES, b::PetscVec, x::PetscVec )

    @chk ccall(
               (:SNESSolve, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec),
               snes, b, x,
              )


	return nothing
end 

"""
	SNESSetType(petsclib::PetscLibType,snes::PetscSNES, type::SNESType) 
Sets the algorithm/method to be used to solve the nonlinear system with the given `SNES`

Collective

Input Parameters:
- `snes` - the `SNES` context
- `type` - a known method

Options Database Key:
- `-snes_type <type>` - Sets the method; use -help for a list
of available methods (for instance, newtonls or newtontr)

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESSolve()`, `SNESType`, `SNESCreate()`, `SNESDestroy()`, `SNESGetType()`, `SNESSetFromOptions()`

# External Links
$(_doc_external("Snes/SNESSetType"))
"""
function SNESSetType(petsclib::PetscLibType, snes::PetscSNES, type::SNESType) end

@for_petsc function SNESSetType(petsclib::$UnionPetscLib, snes::PetscSNES, type::SNESType )

    @chk ccall(
               (:SNESSetType, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESType),
               snes, type,
              )


	return nothing
end 

"""
	type::SNESType = SNESGetType(petsclib::PetscLibType,snes::PetscSNES) 
Gets the `SNES` method type and name (as a string).

Not Collective

Input Parameter:
- `snes` - nonlinear solver context

Output Parameter:
- `type` - `SNES` method (a character string)

Level: intermediate

-seealso: [](ch_snes), `SNESSetType()`, `SNESType`, `SNESSetFromOptions()`, `SNES`

# External Links
$(_doc_external("Snes/SNESGetType"))
"""
function SNESGetType(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESGetType(petsclib::$UnionPetscLib, snes::PetscSNES )
	type_ = Ref{SNESType}()

    @chk ccall(
               (:SNESGetType, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{SNESType}),
               snes, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	SNESSetSolution(petsclib::PetscLibType,snes::PetscSNES, u::PetscVec) 
Sets the solution vector for use by the `SNES` routines.

Logically Collective

Input Parameters:
- `snes` - the `SNES` context obtained from `SNESCreate()`
- `u`    - the solution vector

Level: beginner

-seealso: [](ch_snes), `SNES`, `SNESSolve()`, `SNESGetSolution()`, `Vec`

# External Links
$(_doc_external("Snes/SNESSetSolution"))
"""
function SNESSetSolution(petsclib::PetscLibType, snes::PetscSNES, u::PetscVec) end

@for_petsc function SNESSetSolution(petsclib::$UnionPetscLib, snes::PetscSNES, u::PetscVec )

    @chk ccall(
               (:SNESSetSolution, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec),
               snes, u,
              )


	return nothing
end 

"""
	SNESGetSolution(petsclib::PetscLibType,snes::PetscSNES, x::PetscVec) 
Returns the vector where the approximate solution is
stored. This is the fine grid solution when using `SNESSetGridSequence()`.

Not Collective, but `x` is parallel if `snes` is parallel

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `x` - the solution

Level: intermediate

-seealso: [](ch_snes), `SNESSetSolution()`, `SNESSolve()`, `SNES`, `SNESGetSolutionUpdate()`, `SNESGetFunction()`

# External Links
$(_doc_external("Snes/SNESGetSolution"))
"""
function SNESGetSolution(petsclib::PetscLibType, snes::PetscSNES, x::PetscVec) end

@for_petsc function SNESGetSolution(petsclib::$UnionPetscLib, snes::PetscSNES, x::PetscVec )
	x_ = Ref(x.ptr)

    @chk ccall(
               (:SNESGetSolution, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CVec}),
               snes, x_,
              )

	x.ptr = C_NULL

	return nothing
end 

"""
	SNESGetSolutionUpdate(petsclib::PetscLibType,snes::PetscSNES, x::PetscVec) 
Returns the vector where the solution update is
stored.

Not Collective, but `x` is parallel if `snes` is parallel

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `x` - the solution update

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESGetSolution()`, `SNESGetFunction()`

# External Links
$(_doc_external("Snes/SNESGetSolutionUpdate"))
"""
function SNESGetSolutionUpdate(petsclib::PetscLibType, snes::PetscSNES, x::PetscVec) end

@for_petsc function SNESGetSolutionUpdate(petsclib::$UnionPetscLib, snes::PetscSNES, x::PetscVec )
	x_ = Ref(x.ptr)

    @chk ccall(
               (:SNESGetSolutionUpdate, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CVec}),
               snes, x_,
              )

	x.ptr = C_NULL

	return nothing
end 

"""
	SNESGetFunction(petsclib::PetscLibType,snes::PetscSNES, r::PetscVec, f::SNESFunctionFn, ctx::Cvoid) 
Returns the function that defines the nonlinear system set with `SNESSetFunction()`

Not Collective, but `r` is parallel if `snes` is parallel. Collective if `r` is requested, but has not been created yet.

Input Parameter:
- `snes` - the `SNES` context

Output Parameters:
- `r`   - the vector that is used to store residuals (or `NULL` if you don't want it)
- `f`   - the function (or `NULL` if you don't want it);  for calling sequence see `SNESFunctionFn`
- `ctx` - the function context (or `NULL` if you don't want it)

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESSolve()`, `SNESSetFunction()`, `SNESGetSolution()`, `SNESFunctionFn`

# External Links
$(_doc_external("Snes/SNESGetFunction"))
"""
function SNESGetFunction(petsclib::PetscLibType, snes::PetscSNES, r::PetscVec, f::SNESFunctionFn, ctx::Cvoid) end

@for_petsc function SNESGetFunction(petsclib::$UnionPetscLib, snes::PetscSNES, r::PetscVec, f::SNESFunctionFn, ctx::Cvoid )
	r_ = Ref(r.ptr)

    @chk ccall(
               (:SNESGetFunction, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CVec}, SNESFunctionFn, Cvoid),
               snes, r_, f, ctx,
              )

	r.ptr = C_NULL

	return nothing
end 

"""
	SNESGetNGS(petsclib::PetscLibType,snes::PetscSNES, f::SNESNGSFn, ctx::Cvoid) 
Returns the function and context set with `SNESSetNGS()`

Input Parameter:
- `snes` - the `SNES` context

Output Parameters:
- `f`   - the function (or `NULL`) see `SNESNGSFn` for calling sequence
- `ctx` - the function context (or `NULL`)

Level: advanced

-seealso: [](ch_snes), `SNESSetNGS()`, `SNESGetFunction()`, `SNESNGSFn`

# External Links
$(_doc_external("Snes/SNESGetNGS"))
"""
function SNESGetNGS(petsclib::PetscLibType, snes::PetscSNES, f::SNESNGSFn, ctx::Cvoid) end

@for_petsc function SNESGetNGS(petsclib::$UnionPetscLib, snes::PetscSNES, f::SNESNGSFn, ctx::Cvoid )

    @chk ccall(
               (:SNESGetNGS, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESNGSFn, Cvoid),
               snes, f, ctx,
              )


	return nothing
end 

"""
	SNESSetOptionsPrefix(petsclib::PetscLibType,snes::PetscSNES, prefix::String) 
Sets the prefix used for searching for all
`SNES` options in the database.

Logically Collective

Input Parameters:
- `snes`   - the `SNES` context
- `prefix` - the prefix to prepend to all option names

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESSetFromOptions()`, `SNESAppendOptionsPrefix()`

# External Links
$(_doc_external("Snes/SNESSetOptionsPrefix"))
"""
function SNESSetOptionsPrefix(petsclib::PetscLibType, snes::PetscSNES, prefix::String) end

@for_petsc function SNESSetOptionsPrefix(petsclib::$UnionPetscLib, snes::PetscSNES, prefix::String )

    @chk ccall(
               (:SNESSetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Cchar}),
               snes, prefix,
              )


	return nothing
end 

"""
	SNESAppendOptionsPrefix(petsclib::PetscLibType,snes::PetscSNES, prefix::String) 
Appends to the prefix used for searching for all
`SNES` options in the database.

Logically Collective

Input Parameters:
- `snes`   - the `SNES` context
- `prefix` - the prefix to prepend to all option names

Level: advanced

-seealso: [](ch_snes), `SNESGetOptionsPrefix()`, `SNESSetOptionsPrefix()`

# External Links
$(_doc_external("Snes/SNESAppendOptionsPrefix"))
"""
function SNESAppendOptionsPrefix(petsclib::PetscLibType, snes::PetscSNES, prefix::String) end

@for_petsc function SNESAppendOptionsPrefix(petsclib::$UnionPetscLib, snes::PetscSNES, prefix::String )

    @chk ccall(
               (:SNESAppendOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Cchar}),
               snes, prefix,
              )


	return nothing
end 

"""
	SNESGetOptionsPrefix(petsclib::PetscLibType,snes::PetscSNES, prefix::String) 
Gets the prefix used for searching for all
`SNES` options in the database.

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `prefix` - pointer to the prefix string used

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESSetOptionsPrefix()`, `SNESAppendOptionsPrefix()`

# External Links
$(_doc_external("Snes/SNESGetOptionsPrefix"))
"""
function SNESGetOptionsPrefix(petsclib::PetscLibType, snes::PetscSNES, prefix::String) end

@for_petsc function SNESGetOptionsPrefix(petsclib::$UnionPetscLib, snes::PetscSNES, prefix::String )
	prefix_ = Ref(pointer(prefix))

    @chk ccall(
               (:SNESGetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Ptr{Cchar}}),
               snes, prefix_,
              )


	return nothing
end 

"""
	SNESRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Adds a method to the nonlinear solver package.

Not Collective

Input Parameters:
- `sname`    - name of a new user-defined solver
- `function` - routine to create method context

Level: advanced

-seealso: [](ch_snes), `SNESRegisterAll()`, `SNESRegisterDestroy()`

# External Links
$(_doc_external("Snes/SNESRegister"))
"""
function SNESRegister(petsclib::PetscLibType, sname::String, fnc::external) end

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
	SNESTestLocalMin(petsclib::PetscLibType,snes::PetscSNES) 

# External Links
$(_doc_external("Snes/SNESTestLocalMin"))
"""
function SNESTestLocalMin(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESTestLocalMin(petsclib::$UnionPetscLib, snes::PetscSNES )

    @chk ccall(
               (:SNESTestLocalMin, $petsc_library),
               PetscErrorCode,
               (CSNES,),
               snes,
              )


	return nothing
end 

"""
	SNESKSPSetUseEW(petsclib::PetscLibType,snes::PetscSNES, flag::PetscBool) 
Sets `SNES` to the use Eisenstat
computing relative tolerance for linear solvers within an inexact
Newton method.

Logically Collective

Input Parameters:
- `snes` - `SNES` context
- `flag` - `PETSC_TRUE` or `PETSC_FALSE`

Options Database Keys:
- `-snes_ksp_ew`                       - use Eisenstat-Walker method for determining linear system convergence
- `-snes_ksp_ew_version ver`           - version of  Eisenstat-Walker method
- `-snes_ksp_ew_rtol0 <rtol0>`         - Sets rtol0
- `-snes_ksp_ew_rtolmax <rtolmax>`     - Sets rtolmax
- `-snes_ksp_ew_gamma <gamma>`         - Sets gamma
- `-snes_ksp_ew_alpha <alpha>`         - Sets alpha
- `-snes_ksp_ew_alpha2 <alpha2>`       - Sets alpha2
- `-snes_ksp_ew_threshold <threshold>` - Sets threshold

Level: advanced

-seealso: [](ch_snes), `KSP`, `SNES`, `SNESKSPGetUseEW()`, `SNESKSPGetParametersEW()`, `SNESKSPSetParametersEW()`

# External Links
$(_doc_external("Snes/SNESKSPSetUseEW"))
"""
function SNESKSPSetUseEW(petsclib::PetscLibType, snes::PetscSNES, flag::PetscBool) end

@for_petsc function SNESKSPSetUseEW(petsclib::$UnionPetscLib, snes::PetscSNES, flag::PetscBool )

    @chk ccall(
               (:SNESKSPSetUseEW, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, flag,
              )


	return nothing
end 

"""
	flag::PetscBool = SNESKSPGetUseEW(petsclib::PetscLibType,snes::PetscSNES) 
Gets if `SNES` is using Eisenstat
for computing relative tolerance for linear solvers within an
inexact Newton method.

Not Collective

Input Parameter:
- `snes` - `SNES` context

Output Parameter:
- `flag` - `PETSC_TRUE` or `PETSC_FALSE`

Level: advanced

-seealso: [](ch_snes), `SNESKSPSetUseEW()`, `SNESKSPGetParametersEW()`, `SNESKSPSetParametersEW()`

# External Links
$(_doc_external("Snes/SNESKSPGetUseEW"))
"""
function SNESKSPGetUseEW(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESKSPGetUseEW(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESKSPSetParametersEW(petsclib::PetscLibType,snes::PetscSNES, version::PetscInt, rtol_0::PetscReal, rtol_max::PetscReal, gamma::PetscReal, alpha::PetscReal, alpha2::PetscReal, threshold::PetscReal) 
Sets parameters for Eisenstat
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

-seealso: [](ch_snes), `SNES`, `SNESKSPSetUseEW()`, `SNESKSPGetUseEW()`, `SNESKSPGetParametersEW()`

# External Links
$(_doc_external("Snes/SNESKSPSetParametersEW"))
"""
function SNESKSPSetParametersEW(petsclib::PetscLibType, snes::PetscSNES, version::PetscInt, rtol_0::PetscReal, rtol_max::PetscReal, gamma::PetscReal, alpha::PetscReal, alpha2::PetscReal, threshold::PetscReal) end

@for_petsc function SNESKSPSetParametersEW(petsclib::$UnionPetscLib, snes::PetscSNES, version::$PetscInt, rtol_0::$PetscReal, rtol_max::$PetscReal, gamma::$PetscReal, alpha::$PetscReal, alpha2::$PetscReal, threshold::$PetscReal )

    @chk ccall(
               (:SNESKSPSetParametersEW, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, $PetscReal, $PetscReal, $PetscReal, $PetscReal, $PetscReal),
               snes, version, rtol_0, rtol_max, gamma, alpha, alpha2, threshold,
              )


	return nothing
end 

"""
	version::PetscInt,rtol_0::PetscReal,rtol_max::PetscReal,gamma::PetscReal,alpha::PetscReal,alpha2::PetscReal,threshold::PetscReal = SNESKSPGetParametersEW(petsclib::PetscLibType,snes::PetscSNES) 
Gets parameters for Eisenstat
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

-seealso: [](ch_snes), `SNES`, `SNESKSPSetUseEW()`, `SNESKSPGetUseEW()`, `SNESKSPSetParametersEW()`

# External Links
$(_doc_external("Snes/SNESKSPGetParametersEW"))
"""
function SNESKSPGetParametersEW(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESKSPGetParametersEW(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESGetKSP(petsclib::PetscLibType,snes::PetscSNES, ksp::PetscKSP) 
Returns the `KSP` context for a `SNES` solver.

Not Collective, but if `snes` is parallel, then `ksp` is parallel

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `ksp` - the `KSP` context

Level: beginner

-seealso: [](ch_snes), `SNES`, `KSP`, `PC`, `KSPGetPC()`, `SNESCreate()`, `KSPCreate()`, `SNESSetKSP()`

# External Links
$(_doc_external("Snes/SNESGetKSP"))
"""
function SNESGetKSP(petsclib::PetscLibType, snes::PetscSNES, ksp::PetscKSP) end

@for_petsc function SNESGetKSP(petsclib::$UnionPetscLib, snes::PetscSNES, ksp::PetscKSP )
	ksp_ = Ref(ksp.ptr)

    @chk ccall(
               (:SNESGetKSP, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CKSP}),
               snes, ksp_,
              )

	ksp.ptr = C_NULL

	return nothing
end 

"""
	SNESSetDM(petsclib::PetscLibType,snes::PetscSNES, dm::PetscDM) 
Sets the `DM` that may be used by some `SNES` nonlinear solvers or their underlying preconditioners

Logically Collective

Input Parameters:
- `snes` - the nonlinear solver context
- `dm`   - the `DM`, cannot be `NULL`

Level: intermediate

-seealso: [](ch_snes), `DM`, `SNES`, `SNESGetDM()`, `KSPSetDM()`, `KSPGetDM()`

# External Links
$(_doc_external("Snes/SNESSetDM"))
"""
function SNESSetDM(petsclib::PetscLibType, snes::PetscSNES, dm::PetscDM) end

@for_petsc function SNESSetDM(petsclib::$UnionPetscLib, snes::PetscSNES, dm::PetscDM )

    @chk ccall(
               (:SNESSetDM, $petsc_library),
               PetscErrorCode,
               (CSNES, CDM),
               snes, dm,
              )


	return nothing
end 

"""
	SNESGetDM(petsclib::PetscLibType,snes::PetscSNES, dm::PetscDM) 
Gets the `DM` that may be used by some `SNES` nonlinear solvers/preconditioners

Not Collective but `dm` obtained is parallel on `snes`

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `dm` - the `DM`

Level: intermediate

-seealso: [](ch_snes), `DM`, `SNES`, `SNESSetDM()`, `KSPSetDM()`, `KSPGetDM()`

# External Links
$(_doc_external("Snes/SNESGetDM"))
"""
function SNESGetDM(petsclib::PetscLibType, snes::PetscSNES, dm::PetscDM) end

@for_petsc function SNESGetDM(petsclib::$UnionPetscLib, snes::PetscSNES, dm::PetscDM )
	dm_ = Ref(dm.ptr)

    @chk ccall(
               (:SNESGetDM, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CDM}),
               snes, dm_,
              )

	dm.ptr = C_NULL

	return nothing
end 

"""
	SNESSetNPC(petsclib::PetscLibType,snes::PetscSNES, npc::PetscSNES) 
Sets the nonlinear preconditioner to be used.

Collective

Input Parameters:
- `snes` - iterative context obtained from `SNESCreate()`
- `npc`  - the `SNES` nonlinear preconditioner object

Options Database Key:
- `-npc_snes_type <type>` - set the type of the `SNES` to use as the nonlinear preconditioner

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESNGS`, `SNESFAS`, `SNESGetNPC()`, `SNESHasNPC()`

# External Links
$(_doc_external("Snes/SNESSetNPC"))
"""
function SNESSetNPC(petsclib::PetscLibType, snes::PetscSNES, npc::PetscSNES) end

@for_petsc function SNESSetNPC(petsclib::$UnionPetscLib, snes::PetscSNES, npc::PetscSNES )

    @chk ccall(
               (:SNESSetNPC, $petsc_library),
               PetscErrorCode,
               (CSNES, CSNES),
               snes, npc,
              )


	return nothing
end 

"""
	SNESGetNPC(petsclib::PetscLibType,snes::PetscSNES, pc::PetscSNES) 
Gets a nonlinear preconditioning solver SNES` to be used to precondition the original nonlinear solver.

Not Collective; but any changes to the obtained the `pc` object must be applied collectively

Input Parameter:
- `snes` - iterative context obtained from `SNESCreate()`

Output Parameter:
- `pc` - the `SNES` preconditioner context

Options Database Key:
- `-npc_snes_type <type>` - set the type of the `SNES` to use as the nonlinear preconditioner

Level: advanced

-seealso: [](ch_snes), `SNESSetNPC()`, `SNESHasNPC()`, `SNES`, `SNESCreate()`

# External Links
$(_doc_external("Snes/SNESGetNPC"))
"""
function SNESGetNPC(petsclib::PetscLibType, snes::PetscSNES, pc::PetscSNES) end

@for_petsc function SNESGetNPC(petsclib::$UnionPetscLib, snes::PetscSNES, pc::PetscSNES )
	pc_ = Ref(pc.ptr)

    @chk ccall(
               (:SNESGetNPC, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CSNES}),
               snes, pc_,
              )

	pc.ptr = C_NULL

	return nothing
end 

"""
	has_npc::PetscBool = SNESHasNPC(petsclib::PetscLibType,snes::PetscSNES) 
Returns whether a nonlinear preconditioner is associated with the given `SNES`

Not Collective

Input Parameter:
- `snes` - iterative context obtained from `SNESCreate()`

Output Parameter:
- `has_npc` - whether the `SNES` has a nonlinear preconditioner or not

Level: developer

-seealso: [](ch_snes), `SNESSetNPC()`, `SNESGetNPC()`

# External Links
$(_doc_external("Snes/SNESHasNPC"))
"""
function SNESHasNPC(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESHasNPC(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESSetNPCSide(petsclib::PetscLibType,snes::PetscSNES, side::PCSide) 
Sets the nonlinear preconditioning side used by the nonlinear preconditioner inside `SNES`.

Logically Collective

Input Parameter:
- `snes` - iterative context obtained from `SNESCreate()`

Output Parameter:
- `side` - the preconditioning side, where side is one of
-seealso: [](ch_snes), `SNES`, `SNESGetNPC()`, `SNESNRICHARDSON`, `SNESNCG`, `SNESType`, `SNESGetNPCSide()`, `KSPSetPCSide()`, `PC_LEFT`, `PC_RIGHT`, `PCSide`

# External Links
$(_doc_external("Snes/SNESSetNPCSide"))
"""
function SNESSetNPCSide(petsclib::PetscLibType, snes::PetscSNES, side::PCSide) end

@for_petsc function SNESSetNPCSide(petsclib::$UnionPetscLib, snes::PetscSNES, side::PCSide )

    @chk ccall(
               (:SNESSetNPCSide, $petsc_library),
               PetscErrorCode,
               (CSNES, PCSide),
               snes, side,
              )


	return nothing
end 

"""
	SNESGetNPCSide(petsclib::PetscLibType,snes::PetscSNES, side::PCSide) 
Gets the preconditioning side used by the nonlinear preconditioner inside `SNES`.

Not Collective

Input Parameter:
- `snes` - iterative context obtained from `SNESCreate()`

Output Parameter:
- `side` - the preconditioning side, where side is one of
-seealso: [](ch_snes), `SNES`, `SNESGetNPC()`, `SNESSetNPCSide()`, `KSPGetPCSide()`, `PC_LEFT`, `PC_RIGHT`, `PCSide`

# External Links
$(_doc_external("Snes/SNESGetNPCSide"))
"""
function SNESGetNPCSide(petsclib::PetscLibType, snes::PetscSNES, side::PCSide) end

@for_petsc function SNESGetNPCSide(petsclib::$UnionPetscLib, snes::PetscSNES, side::PCSide )

    @chk ccall(
               (:SNESGetNPCSide, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{PCSide}),
               snes, side,
              )


	return nothing
end 

"""
	SNESSetLineSearch(petsclib::PetscLibType,snes::PetscSNES, linesearch::SNESLineSearch) 
Sets the `SNESLineSearch` to be used for a given `SNES`

Collective

Input Parameters:
- `snes`       - iterative context obtained from `SNESCreate()`
- `linesearch` - the linesearch object

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESGetLineSearch()`

# External Links
$(_doc_external("Snes/SNESSetLineSearch"))
"""
function SNESSetLineSearch(petsclib::PetscLibType, snes::PetscSNES, linesearch::SNESLineSearch) end

@for_petsc function SNESSetLineSearch(petsclib::$UnionPetscLib, snes::PetscSNES, linesearch::SNESLineSearch )

    @chk ccall(
               (:SNESSetLineSearch, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESLineSearch),
               snes, linesearch,
              )


	return nothing
end 

"""
	SNESGetLineSearch(petsclib::PetscLibType,snes::PetscSNES, linesearch::SNESLineSearch) 
Returns the line search associated with the `SNES`.

Not Collective

Input Parameter:
- `snes` - iterative context obtained from `SNESCreate()`

Output Parameter:
- `linesearch` - linesearch context

Level: beginner

-seealso: [](ch_snes), `SNESLineSearch`, `SNESSetLineSearch()`, `SNESLineSearchCreate()`, `SNESLineSearchSetFromOptions()`

# External Links
$(_doc_external("Snes/SNESGetLineSearch"))
"""
function SNESGetLineSearch(petsclib::PetscLibType, snes::PetscSNES, linesearch::SNESLineSearch) end

@for_petsc function SNESGetLineSearch(petsclib::$UnionPetscLib, snes::PetscSNES, linesearch::SNESLineSearch )

    @chk ccall(
               (:SNESGetLineSearch, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{SNESLineSearch}),
               snes, linesearch,
              )


	return nothing
end 

"""
	SNESComputeJacobianDefault(petsclib::PetscLibType,snes::PetscSNES, x1::PetscVec, J::PetscMat, B::PetscMat, ctx::Cvoid) 
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
- `-snes_fd`          - Activates `SNESComputeJacobianDefault()`
- `-snes_fd_coloring` - Activates a faster computation that uses a graph coloring of the matrix
- `-snes_test_err`    - Square root of function error tolerance, default square root of machine
epsilon (1.e-8 in double, 3.e-4 in single)
- `-mat_fd_type`      - Either wp or ds (see `MATMFFD_WP` or `MATMFFD_DS`)

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESSetJacobian()`, `SNESComputeJacobianDefaultColor()`, `MatCreateSNESMF()`

# External Links
$(_doc_external("Snes/SNESComputeJacobianDefault"))
"""
function SNESComputeJacobianDefault(petsclib::PetscLibType, snes::PetscSNES, x1::PetscVec, J::PetscMat, B::PetscMat, ctx::Cvoid) end

@for_petsc function SNESComputeJacobianDefault(petsclib::$UnionPetscLib, snes::PetscSNES, x1::PetscVec, J::PetscMat, B::PetscMat, ctx::Cvoid )

    @chk ccall(
               (:SNESComputeJacobianDefault, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CMat, CMat, Ptr{Cvoid}),
               snes, x1, J, B, ctx,
              )


	return nothing
end 

"""
	SNESApplyNPC(petsclib::PetscLibType,snes::PetscSNES, x::PetscVec, f::PetscVec, y::PetscVec) 
Calls `SNESSolve()` on the preconditioner for the `SNES`

Collective

Input Parameters:
- `snes` - the `SNES` context
- `x`    - input vector
- `f`    - optional; the function evaluation on `x`

Output Parameter:
- `y` - function vector, as set by `SNESSetFunction()`

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESGetNPC()`, `SNESSetNPC()`, `SNESComputeFunction()`

# External Links
$(_doc_external("Snes/SNESApplyNPC"))
"""
function SNESApplyNPC(petsclib::PetscLibType, snes::PetscSNES, x::PetscVec, f::PetscVec, y::PetscVec) end

@for_petsc function SNESApplyNPC(petsclib::$UnionPetscLib, snes::PetscSNES, x::PetscVec, f::PetscVec, y::PetscVec )

    @chk ccall(
               (:SNESApplyNPC, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec, CVec),
               snes, x, f, y,
              )


	return nothing
end 

"""
	SNESComputeFunctionDefaultNPC(petsclib::PetscLibType,snes::PetscSNES, X::PetscVec, F::PetscVec) 

# External Links
$(_doc_external("Snes/SNESComputeFunctionDefaultNPC"))
"""
function SNESComputeFunctionDefaultNPC(petsclib::PetscLibType, snes::PetscSNES, X::PetscVec, F::PetscVec) end

@for_petsc function SNESComputeFunctionDefaultNPC(petsclib::$UnionPetscLib, snes::PetscSNES, X::PetscVec, F::PetscVec )

    @chk ccall(
               (:SNESComputeFunctionDefaultNPC, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec),
               snes, X, F,
              )


	return nothing
end 

"""
	fnorm::PetscReal = SNESGetNPCFunction(petsclib::PetscLibType,snes::PetscSNES, F::PetscVec) 
Gets the current function value and its norm from a nonlinear preconditioner after `SNESSolve()` has been called on that `SNES`

Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameters:
- `F`     - function vector
- `fnorm` - the norm of `F`

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESGetNPC()`, `SNESSetNPC()`, `SNESComputeFunction()`, `SNESApplyNPC()`, `SNESSolve()`

# External Links
$(_doc_external("Snes/SNESGetNPCFunction"))
"""
function SNESGetNPCFunction(petsclib::PetscLibType, snes::PetscSNES, F::PetscVec) end

@for_petsc function SNESGetNPCFunction(petsclib::$UnionPetscLib, snes::PetscSNES, F::PetscVec )
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
	SNESSetObjective(petsclib::PetscLibType,snes::PetscSNES, obj::SNESObjectiveFn, ctx::Cvoid) 
Sets the objective function minimized by some of the `SNES` linesearch methods, used instead of the 2

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `obj`  - objective evaluation routine; see `SNESObjectiveFn` for the calling sequence
- `ctx`  - [optional] user-defined context for private data for the objective evaluation routine (may be `NULL`)

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESLineSearch()`, `SNESGetObjective()`, `SNESComputeObjective()`, `SNESSetFunction()`, `SNESSetJacobian()`,
`SNESObjectiveFn`

# External Links
$(_doc_external("Snes/SNESSetObjective"))
"""
function SNESSetObjective(petsclib::PetscLibType, snes::PetscSNES, obj::SNESObjectiveFn, ctx::Cvoid) end

@for_petsc function SNESSetObjective(petsclib::$UnionPetscLib, snes::PetscSNES, obj::SNESObjectiveFn, ctx::Cvoid )

    @chk ccall(
               (:SNESSetObjective, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{SNESObjectiveFn}, Ptr{Cvoid}),
               snes, obj, ctx,
              )


	return nothing
end 

"""
	SNESGetObjective(petsclib::PetscLibType,snes::PetscSNES, obj::SNESObjectiveFn, ctx::Cvoid) 
Returns the objective function set with `SNESSetObjective()`

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameters:
- `obj` - objective evaluation routine (or `NULL`); see `SNESObjectiveFn` for the calling sequence
- `ctx` - the function context (or `NULL`)

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESSetObjective()`, `SNESGetSolution()`, `SNESObjectiveFn`

# External Links
$(_doc_external("Snes/SNESGetObjective"))
"""
function SNESGetObjective(petsclib::PetscLibType, snes::PetscSNES, obj::SNESObjectiveFn, ctx::Cvoid) end

@for_petsc function SNESGetObjective(petsclib::$UnionPetscLib, snes::PetscSNES, obj::SNESObjectiveFn, ctx::Cvoid )

    @chk ccall(
               (:SNESGetObjective, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESObjectiveFn, Cvoid),
               snes, obj, ctx,
              )


	return nothing
end 

"""
	ob::PetscReal = SNESComputeObjective(petsclib::PetscLibType,snes::PetscSNES, X::PetscVec) 
Computes the objective function that has been provided by `SNESSetObjective()`

Collective

Input Parameters:
- `snes` - the `SNES` context
- `X`    - the state vector

Output Parameter:
- `ob` - the objective value

Level: developer

-seealso: [](ch_snes), `SNESLineSearch`, `SNES`, `SNESSetObjective()`, `SNESGetSolution()`

# External Links
$(_doc_external("Snes/SNESComputeObjective"))
"""
function SNESComputeObjective(petsclib::PetscLibType, snes::PetscSNES, X::PetscVec) end

@for_petsc function SNESComputeObjective(petsclib::$UnionPetscLib, snes::PetscSNES, X::PetscVec )
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
	SNESObjectiveComputeFunctionDefaultFD(petsclib::PetscLibType,snes::PetscSNES, X::PetscVec, F::PetscVec, ctx::Cvoid) 
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

-seealso: [](ch_snes), `SNESSetObjective()`, `SNESSetFunction()`, `SNESComputeObjective()`, `SNESComputeJacobianDefault()`, `SNESObjectiveFn`

# External Links
$(_doc_external("Snes/SNESObjectiveComputeFunctionDefaultFD"))
"""
function SNESObjectiveComputeFunctionDefaultFD(petsclib::PetscLibType, snes::PetscSNES, X::PetscVec, F::PetscVec, ctx::Cvoid) end

@for_petsc function SNESObjectiveComputeFunctionDefaultFD(petsclib::$UnionPetscLib, snes::PetscSNES, X::PetscVec, F::PetscVec, ctx::Cvoid )

    @chk ccall(
               (:SNESObjectiveComputeFunctionDefaultFD, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec, Ptr{Cvoid}),
               snes, X, F, ctx,
              )


	return nothing
end 

"""
	SNESComputeJacobianDefaultColor(petsclib::PetscLibType,snes::PetscSNES, x1::PetscVec, J::PetscMat, B::PetscMat, ctx::Cvoid) 
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
- `-snes_fd_color_use_mat`       - use a matrix coloring from the explicit matrix nonzero pattern instead of from the `DM` providing the matrix
- `-snes_fd_color`               - Activates `SNESComputeJacobianDefaultColor()` in `SNESSetFromOptions()`
- `-mat_fd_coloring_err <err>`   - Sets <err> (square root of relative error in the function)
- `-mat_fd_coloring_umin <umin>` - Sets umin, the minimum allowable u-value magnitude
- `-mat_fd_type`                 - Either wp or ds (see `MATMFFD_WP` or `MATMFFD_DS`)
- `-snes_mf_operator`            - Use matrix-free application of Jacobian
- `-snes_mf`                     - Use matrix-free Jacobian with no explicit Jacobian representation

-seealso: [](ch_snes), `SNES`, `SNESSetJacobian()`, `SNESTestJacobian()`, `SNESComputeJacobianDefault()`, `SNESSetUseMatrixFree()`,
`MatFDColoringCreate()`, `MatFDColoringSetFunction()`

# External Links
$(_doc_external("Snes/SNESComputeJacobianDefaultColor"))
"""
function SNESComputeJacobianDefaultColor(petsclib::PetscLibType, snes::PetscSNES, x1::PetscVec, J::PetscMat, B::PetscMat, ctx::Cvoid) end

@for_petsc function SNESComputeJacobianDefaultColor(petsclib::$UnionPetscLib, snes::PetscSNES, x1::PetscVec, J::PetscMat, B::PetscMat, ctx::Cvoid )

    @chk ccall(
               (:SNESComputeJacobianDefaultColor, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CMat, CMat, Ptr{Cvoid}),
               snes, x1, J, B, ctx,
              )


	return nothing
end 

"""
	SNESPruneJacobianColor(petsclib::PetscLibType,snes::PetscSNES, J::PetscMat, B::PetscMat) 
Remove nondiagonal zeros in the Jacobian matrix and update the `MatMFFD` coloring information based on the new nonzero structure

Collective

Input Parameters:
- `snes` - the `SNES` context
- `J`    - Jacobian matrix (not altered in this routine)
- `B`    - newly computed Jacobian matrix to use with preconditioner (generally the same as `J`)

Level: intermediate

-seealso: [](ch_snes), `SNESComputeJacobianDefaultColor()`, `MatEliminateZeros()`, `MatFDColoringCreate()`, `MatFDColoringSetFunction()`

# External Links
$(_doc_external("Snes/SNESPruneJacobianColor"))
"""
function SNESPruneJacobianColor(petsclib::PetscLibType, snes::PetscSNES, J::PetscMat, B::PetscMat) end

@for_petsc function SNESPruneJacobianColor(petsclib::$UnionPetscLib, snes::PetscSNES, J::PetscMat, B::PetscMat )

    @chk ccall(
               (:SNESPruneJacobianColor, $petsc_library),
               PetscErrorCode,
               (CSNES, CMat, CMat),
               snes, J, B,
              )


	return nothing
end 

"""
	ctx::Cvoid = SNESMonitorSAWsCreate(petsclib::PetscLibType,snes::PetscSNES) 
create an SAWs monitor context for `SNES`

Collective

Input Parameter:
- `snes` - `SNES` to monitor

Output Parameter:
- `ctx` - context for monitor

Level: developer

-seealso: [](ch_snes), `SNESMonitorSet()`, `SNES`, `SNESMonitorSAWs()`, `SNESMonitorSAWsDestroy()`

# External Links
$(_doc_external("Snes/SNESMonitorSAWsCreate"))
"""
function SNESMonitorSAWsCreate(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESMonitorSAWsCreate(petsclib::$UnionPetscLib, snes::PetscSNES )
	ctx_ = Ref{Cvoid}()

    @chk ccall(
               (:SNESMonitorSAWsCreate, $petsc_library),
               PetscErrorCode,
               (CSNES, Cvoid),
               snes, ctx_,
              )

	ctx = ctx_[]

	return ctx
end 

"""
	SNESMonitorSAWsDestroy(petsclib::PetscLibType,ctx::Cvoid) 
destroy a monitor context created with `SNESMonitorSAWsCreate()`

Collective

Input Parameter:
- `ctx` - monitor context

Level: developer

-seealso: [](ch_snes), `SNESMonitorSAWsCreate()`

# External Links
$(_doc_external("Snes/SNESMonitorSAWsDestroy"))
"""
function SNESMonitorSAWsDestroy(petsclib::PetscLibType, ctx::Cvoid) end

@for_petsc function SNESMonitorSAWsDestroy(petsclib::$UnionPetscLib, ctx::Cvoid )

    @chk ccall(
               (:SNESMonitorSAWsDestroy, $petsc_library),
               PetscErrorCode,
               (Cvoid,),
               ctx,
              )


	return nothing
end 

"""
	SNESMonitorSAWs(petsclib::PetscLibType,snes::PetscSNES, n::PetscInt, rnorm::PetscReal, ctx::Cvoid) 
monitor solution process of `SNES` using SAWs

Collective

Input Parameters:
- `snes`  - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `ctx`   - `PetscViewer` of type `PETSCVIEWERSAWS`

Level: advanced

-seealso: [](ch_snes), `PetscViewerSAWsOpen()`, `SNESMonitorSAWsDestroy()`, `SNESMonitorSAWsCreate()`

# External Links
$(_doc_external("Snes/SNESMonitorSAWs"))
"""
function SNESMonitorSAWs(petsclib::PetscLibType, snes::PetscSNES, n::PetscInt, rnorm::PetscReal, ctx::Cvoid) end

@for_petsc function SNESMonitorSAWs(petsclib::$UnionPetscLib, snes::PetscSNES, n::$PetscInt, rnorm::$PetscReal, ctx::Cvoid )

    @chk ccall(
               (:SNESMonitorSAWs, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, Ptr{Cvoid}),
               snes, n, rnorm, ctx,
              )


	return nothing
end 

"""
	SNESNGMRESSetRestartFmRise(petsclib::PetscLibType,snes::PetscSNES, flg::PetscBool) 
Increase the restart count if the step x_M increases the residual F_M inside a `SNESNGMRES` solve

Input Parameters:
- `snes` - the `SNES` context.
- `flg`  - boolean value deciding whether to use the option or not, default is `PETSC_FALSE`

Options Database Key:
- `-snes_ngmres_restart_fm_rise` - Increase the restart count if the step x_M increases the residual F_M

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNES_NGMRES_RESTART_DIFFERENCE`, `SNESNGMRES`, `SNESNGMRESRestartType`, `SNESNGMRESSetRestartType()`

# External Links
$(_doc_external("Snes/SNESNGMRESSetRestartFmRise"))
"""
function SNESNGMRESSetRestartFmRise(petsclib::PetscLibType, snes::PetscSNES, flg::PetscBool) end

@for_petsc function SNESNGMRESSetRestartFmRise(petsclib::$UnionPetscLib, snes::PetscSNES, flg::PetscBool )

    @chk ccall(
               (:SNESNGMRESSetRestartFmRise, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, flg,
              )


	return nothing
end 

"""
	flg::PetscBool = SNESNGMRESGetRestartFmRise(petsclib::PetscLibType,snes::PetscSNES) 

# External Links
$(_doc_external("Snes/SNESNGMRESGetRestartFmRise"))
"""
function SNESNGMRESGetRestartFmRise(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESNGMRESGetRestartFmRise(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESNGMRESSetRestartType(petsclib::PetscLibType,snes::PetscSNES, rtype::SNESNGMRESRestartType) 
Sets the restart type for `SNESNGMRES`.

Logically Collective

Input Parameters:
- `snes`  - the iterative context
- `rtype` - restart type, see `SNESNGMRESRestartType`

Options Database Keys:
- `-snes_ngmres_restart_type<difference,periodic,none>` - set the restart type
- `-snes_ngmres_restart <30>`                           - sets the number of iterations before restart for periodic

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNES_NGMRES_RESTART_DIFFERENCE`, `SNESNGMRES`, `SNESNGMRESRestartType`, `SNESNGMRESSetRestartFmRise()`,
`SNESNGMRESSetSelectType()`

# External Links
$(_doc_external("Snes/SNESNGMRESSetRestartType"))
"""
function SNESNGMRESSetRestartType(petsclib::PetscLibType, snes::PetscSNES, rtype::SNESNGMRESRestartType) end

@for_petsc function SNESNGMRESSetRestartType(petsclib::$UnionPetscLib, snes::PetscSNES, rtype::SNESNGMRESRestartType )

    @chk ccall(
               (:SNESNGMRESSetRestartType, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESNGMRESRestartType),
               snes, rtype,
              )


	return nothing
end 

"""
	SNESNGMRESSetSelectType(petsclib::PetscLibType,snes::PetscSNES, stype::SNESNGMRESSelectType) 
Sets the selection type for `SNESNGMRES`.  This determines how the candidate solution and
combined solution are used to create the next iterate.

Logically Collective

Input Parameters:
- `snes`  - the iterative context
- `stype` - selection type, see `SNESNGMRESSelectType`

Options Database Key:
- `-snes_ngmres_select_type<difference,none,linesearch>` - select type

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESNGMRES`, `SNESNGMRESSelectType`, `SNES_NGMRES_SELECT_NONE`, `SNES_NGMRES_SELECT_DIFFERENCE`, `SNES_NGMRES_SELECT_LINESEARCH`,
`SNESNGMRESSetRestartType()`

# External Links
$(_doc_external("Snes/SNESNGMRESSetSelectType"))
"""
function SNESNGMRESSetSelectType(petsclib::PetscLibType, snes::PetscSNES, stype::SNESNGMRESSelectType) end

@for_petsc function SNESNGMRESSetSelectType(petsclib::$UnionPetscLib, snes::PetscSNES, stype::SNESNGMRESSelectType )

    @chk ccall(
               (:SNESNGMRESSetSelectType, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESNGMRESSelectType),
               snes, stype,
              )


	return nothing
end 

"""
	SNESVISetComputeVariableBounds(petsclib::PetscLibType,snes::PetscSNES, compute::external) 
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

-seealso: [](sec_vi), `SNES`, `SNESVISetVariableBounds()`, `DMSetVariableBounds()`, `SNESSetFunctionDomainError()`, `SNESSetJacobianDomainError()`, `SNESVINEWTONRSLS`, `SNESVINEWTONSSLS`,
`SNESSetType()`, `PETSC_NINFINITY`, `PETSC_INFINITY`

# External Links
$(_doc_external("Snes/SNESVISetComputeVariableBounds"))
"""
function SNESVISetComputeVariableBounds(petsclib::PetscLibType, snes::PetscSNES, compute::external) end

@for_petsc function SNESVISetComputeVariableBounds(petsclib::$UnionPetscLib, snes::PetscSNES, compute::external )

    @chk ccall(
               (:SNESVISetComputeVariableBounds, $petsc_library),
               PetscErrorCode,
               (CSNES, external),
               snes, compute,
              )


	return nothing
end 

"""
	SNESVIGetActiveSetIS(petsclib::PetscLibType,snes::PetscSNES, X::PetscVec, F::PetscVec, ISact::IS) 
Gets the global indices for the active set variables

Input Parameters:
- `snes` - the `SNES` context
- `X`    - the `snes` solution vector
- `F`    - the nonlinear function vector

Output Parameter:
- `ISact` - active set index set

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESVINEWTONRSLS`, `SNESVINEWTONSSLS`

# External Links
$(_doc_external("Snes/SNESVIGetActiveSetIS"))
"""
function SNESVIGetActiveSetIS(petsclib::PetscLibType, snes::PetscSNES, X::PetscVec, F::PetscVec, ISact::IS) end

@for_petsc function SNESVIGetActiveSetIS(petsclib::$UnionPetscLib, snes::PetscSNES, X::PetscVec, F::PetscVec, ISact::IS )

    @chk ccall(
               (:SNESVIGetActiveSetIS, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec, Ptr{IS}),
               snes, X, F, ISact,
              )


	return nothing
end 

"""
	fnorm::PetscReal = SNESVIComputeInactiveSetFnorm(petsclib::PetscLibType,snes::PetscSNES, F::PetscVec, X::PetscVec) 
Computes the function norm for variational inequalities on the inactive set

Input Parameters:
- `snes` - the `SNES` context
- `F`    - the nonlinear function vector
- `X`    - the `SNES` solution vector

Output Parameter:
- `fnorm` - the function norm

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESVINEWTONRSLS`, `SNESVINEWTONSSLS`, `SNESLineSearchSetVIFunctions()`

# External Links
$(_doc_external("Snes/SNESVIComputeInactiveSetFnorm"))
"""
function SNESVIComputeInactiveSetFnorm(petsclib::PetscLibType, snes::PetscSNES, F::PetscVec, X::PetscVec) end

@for_petsc function SNESVIComputeInactiveSetFnorm(petsclib::$UnionPetscLib, snes::PetscSNES, F::PetscVec, X::PetscVec )
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
	fty::PetscScalar = SNESVIComputeInactiveSetFtY(petsclib::PetscLibType,snes::PetscSNES, F::PetscVec, X::PetscVec, Y::PetscVec) 
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

-seealso: [](ch_snes), `SNES`, `SNESVINEWTONRSLS`, `SNESVINEWTONSSLS`

# External Links
$(_doc_external("Snes/SNESVIComputeInactiveSetFtY"))
"""
function SNESVIComputeInactiveSetFtY(petsclib::PetscLibType, snes::PetscSNES, F::PetscVec, X::PetscVec, Y::PetscVec) end

@for_petsc function SNESVIComputeInactiveSetFtY(petsclib::$UnionPetscLib, snes::PetscSNES, F::PetscVec, X::PetscVec, Y::PetscVec )
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
	SNESVISetVariableBounds(petsclib::PetscLibType,snes::PetscSNES, xl::PetscVec, xu::PetscVec) 
Sets the lower and upper bounds for the solution vector. `xl` <= x <= `xu`. This allows solving
(differential) variable inequalities.

Input Parameters:
- `snes` - the `SNES` context.
- `xl`   - lower bound.
- `xu`   - upper bound.

Level: advanced

-seealso: [](sec_vi), `SNES`, `SNESVIGetVariableBounds()`, `SNESVISetComputeVariableBounds()`, `SNESSetFunctionDomainError()`, `SNESSetJacobianDomainError()`, `SNESVINEWTONRSLS`, `SNESVINEWTONSSLS`, `SNESSetType()`, `PETSC_NINFINITY`, `PETSC_INFINITY`

# External Links
$(_doc_external("Snes/SNESVISetVariableBounds"))
"""
function SNESVISetVariableBounds(petsclib::PetscLibType, snes::PetscSNES, xl::PetscVec, xu::PetscVec) end

@for_petsc function SNESVISetVariableBounds(petsclib::$UnionPetscLib, snes::PetscSNES, xl::PetscVec, xu::PetscVec )

    @chk ccall(
               (:SNESVISetVariableBounds, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec),
               snes, xl, xu,
              )


	return nothing
end 

"""
	SNESVIGetVariableBounds(petsclib::PetscLibType,snes::PetscSNES, xl::PetscVec, xu::PetscVec) 
Gets the lower and upper bounds for the solution vector. `xl` <= x <= `xu`. These are used in solving
(differential) variable inequalities.

Input Parameters:
- `snes` - the `SNES` context.
- `xl`   - lower bound (may be `NULL`)
- `xu`   - upper bound (may be `NULL`)

Level: advanced

-seealso: [](sec_vi), `SNES`, `SNESVISetVariableBounds()`, `SNESVISetComputeVariableBounds()`, `SNESSetFunctionDomainError()`, `SNESSetJacobianDomainError()`, `SNESVINEWTONRSLS`, `SNESVINEWTONSSLS`, `SNESSetType()`, `PETSC_NINFINITY`, `PETSC_INFINITY`

# External Links
$(_doc_external("Snes/SNESVIGetVariableBounds"))
"""
function SNESVIGetVariableBounds(petsclib::PetscLibType, snes::PetscSNES, xl::PetscVec, xu::PetscVec) end

@for_petsc function SNESVIGetVariableBounds(petsclib::$UnionPetscLib, snes::PetscSNES, xl::PetscVec, xu::PetscVec )
	xl_ = Ref(xl.ptr)
	xu_ = Ref(xu.ptr)

    @chk ccall(
               (:SNESVIGetVariableBounds, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CVec}, Ptr{CVec}),
               snes, xl_, xu_,
              )

	xl.ptr = C_NULL
	xu.ptr = C_NULL

	return nothing
end 

"""
	SNESVIGetInactiveSet(petsclib::PetscLibType,snes::PetscSNES, inact::IS) 
Gets the global indices for the inactive set variables (these correspond to the degrees of freedom the linear
system is solved on)

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `inact` - inactive set index set

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESVINEWTONRSLS`

# External Links
$(_doc_external("Snes/SNESVIGetInactiveSet"))
"""
function SNESVIGetInactiveSet(petsclib::PetscLibType, snes::PetscSNES, inact::IS) end

@for_petsc function SNESVIGetInactiveSet(petsclib::$UnionPetscLib, snes::PetscSNES, inact::IS )

    @chk ccall(
               (:SNESVIGetInactiveSet, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{IS}),
               snes, inact,
              )


	return nothing
end 

"""
	SNESVISetRedundancyCheck(petsclib::PetscLibType,snes::PetscSNES, func::external, ctx::Cvoid) 
Provide a function to check for any redundancy in the VI active set

Logically Collective

Input Parameters:
- `snes` - the `SNESVINEWTONRSLS` context
- `func` - the function to check of redundancies
- `ctx`  - optional context used by the function

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESVINEWTONRSLS`, `SNESVIGetInactiveSet()`, `DMSetVI()`

# External Links
$(_doc_external("Snes/SNESVISetRedundancyCheck"))
"""
function SNESVISetRedundancyCheck(petsclib::PetscLibType, snes::PetscSNES, func::external, ctx::Cvoid) end

@for_petsc function SNESVISetRedundancyCheck(petsclib::$UnionPetscLib, snes::PetscSNES, func::external, ctx::Cvoid )

    @chk ccall(
               (:SNESVISetRedundancyCheck, $petsc_library),
               PetscErrorCode,
               (CSNES, external, Ptr{Cvoid}),
               snes, func, ctx,
              )


	return nothing
end 

"""
	merit::PetscReal,phinorm::PetscReal = SNESVIComputeMeritFunction(petsclib::PetscLibType,phi::PetscVec) 
Evaluates the merit function for the mixed complementarity problem.

Input Parameter:
- `phi` - the `Vec` holding the evaluation of the semismooth function

Output Parameters:
- `merit`   - the merit function 1/2 ||phi||^2
- `phinorm` - the two-norm of the vector, ||phi||

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESVINEWTONSSLS`, `SNESVIComputeFunction()`

# External Links
$(_doc_external("Snes/SNESVIComputeMeritFunction"))
"""
function SNESVIComputeMeritFunction(petsclib::PetscLibType, phi::PetscVec) end

@for_petsc function SNESVIComputeMeritFunction(petsclib::$UnionPetscLib, phi::PetscVec )
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
	SNESVIComputeFunction(petsclib::PetscLibType,snes::PetscSNES, X::PetscVec, phi::PetscVec, functx::Cvoid) 
Provides the function that reformulates a system of nonlinear equations in mixed complementarity form to a system of nonlinear
equations in semismooth form.

Input Parameters:
- `snes`   - the `SNES` context
- `X`      - current iterate
- `functx` - user defined function context

Output Parameter:
- `phi` - the evaluation of semismooth function at `X`

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESVINEWTONSSLS`, `SNESVIComputeMeritFunction()`

# External Links
$(_doc_external("Snes/SNESVIComputeFunction"))
"""
function SNESVIComputeFunction(petsclib::PetscLibType, snes::PetscSNES, X::PetscVec, phi::PetscVec, functx::Cvoid) end

@for_petsc function SNESVIComputeFunction(petsclib::$UnionPetscLib, snes::PetscSNES, X::PetscVec, phi::PetscVec, functx::Cvoid )

    @chk ccall(
               (:SNESVIComputeFunction, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec, Ptr{Cvoid}),
               snes, X, phi, functx,
              )


	return nothing
end 

"""
	SNESMSRegisterAll(petsclib::PetscLibType) 
Registers all of the multi

Logically Collective

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESMS`, `SNESMSRegisterDestroy()`

# External Links
$(_doc_external("Snes/SNESMSRegisterAll"))
"""
function SNESMSRegisterAll(petsclib::PetscLibType) end

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

-seealso: [](ch_snes), `SNES`, `SNESMS`, `SNESMSRegister()`, `SNESMSRegisterAll()`

# External Links
$(_doc_external("Snes/SNESMSRegisterDestroy"))
"""
function SNESMSRegisterDestroy(petsclib::PetscLibType) end

@for_petsc function SNESMSRegisterDestroy(petsclib::$UnionPetscLib)

    @chk ccall(
               (:SNESMSRegisterDestroy, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	SNESMSInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `SNESMS` package. It is called
from `SNESInitializePackage()`.

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESMS`, `SNESMSRegister()`, `SNESMSRegisterAll()`, `PetscInitialize()`

# External Links
$(_doc_external("Snes/SNESMSInitializePackage"))
"""
function SNESMSInitializePackage(petsclib::PetscLibType) end

@for_petsc function SNESMSInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:SNESMSInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	SNESMSFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the `SNESMS` package. It is
called from `PetscFinalize()`.

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESMS`, `SNESMSRegister()`, `SNESMSRegisterAll()`, `SNESMSInitializePackage()`, `PetscFinalize()`

# External Links
$(_doc_external("Snes/SNESMSFinalizePackage"))
"""
function SNESMSFinalizePackage(petsclib::PetscLibType) end

@for_petsc function SNESMSFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:SNESMSFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	SNESMSRegister(petsclib::PetscLibType,name::SNESMSType, nstages::PetscInt, nregisters::PetscInt, stability::PetscReal, gamma::Vector{PetscReal}, delta::Vector{PetscReal}, betasub::Vector{PetscReal}) 
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

-seealso: [](ch_snes), `SNES`, `SNESMS`

# External Links
$(_doc_external("Snes/SNESMSRegister"))
"""
function SNESMSRegister(petsclib::PetscLibType, name::SNESMSType, nstages::PetscInt, nregisters::PetscInt, stability::PetscReal, gamma::Vector{PetscReal}, delta::Vector{PetscReal}, betasub::Vector{PetscReal}) end

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
	mstype::SNESMSType = SNESMSGetType(petsclib::PetscLibType,snes::PetscSNES) 
Get the type of multistage smoother `SNESMS`

Not Collective

Input Parameter:
- `snes` - nonlinear solver context

Output Parameter:
- `mstype` - type of multistage method

Level: advanced

-seealso: [](ch_snes), `SNESMS`, `SNESMSSetType()`, `SNESMSType`

# External Links
$(_doc_external("Snes/SNESMSGetType"))
"""
function SNESMSGetType(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESMSGetType(petsclib::$UnionPetscLib, snes::PetscSNES )
	mstype_ = Ref{SNESMSType}()

    @chk ccall(
               (:SNESMSGetType, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{SNESMSType}),
               snes, mstype_,
              )

	mstype = unsafe_string(mstype_[])

	return mstype
end 

"""
	SNESMSSetType(petsclib::PetscLibType,snes::PetscSNES, mstype::SNESMSType) 
Set the type of multistage smoother `SNESMS`

Logically Collective

Input Parameters:
- `snes`   - nonlinear solver context
- `mstype` - type of multistage method

Level: advanced

-seealso: [](ch_snes), `SNESMS`, `SNESMSGetType()`, `SNESMSType`

# External Links
$(_doc_external("Snes/SNESMSSetType"))
"""
function SNESMSSetType(petsclib::PetscLibType, snes::PetscSNES, mstype::SNESMSType) end

@for_petsc function SNESMSSetType(petsclib::$UnionPetscLib, snes::PetscSNES, mstype::SNESMSType )

    @chk ccall(
               (:SNESMSSetType, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESMSType),
               snes, mstype,
              )


	return nothing
end 

"""
	damping::PetscReal = SNESMSGetDamping(petsclib::PetscLibType,snes::PetscSNES) 
Get the damping parameter of `SNESMS` multistage scheme

Not Collective

Input Parameter:
- `snes` - nonlinear solver context

Output Parameter:
- `damping` - damping parameter

Level: advanced

-seealso: [](ch_snes), `SNESMSSetDamping()`, `SNESMS`

# External Links
$(_doc_external("Snes/SNESMSGetDamping"))
"""
function SNESMSGetDamping(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESMSGetDamping(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESMSSetDamping(petsclib::PetscLibType,snes::PetscSNES, damping::PetscReal) 
Set the damping parameter for a `SNESMS` multistage scheme

Logically Collective

Input Parameters:
- `snes`    - nonlinear solver context
- `damping` - damping parameter

Level: advanced

-seealso: [](ch_snes), `SNESMSGetDamping()`, `SNESMS`

# External Links
$(_doc_external("Snes/SNESMSSetDamping"))
"""
function SNESMSSetDamping(petsclib::PetscLibType, snes::PetscSNES, damping::PetscReal) end

@for_petsc function SNESMSSetDamping(petsclib::$UnionPetscLib, snes::PetscSNES, damping::$PetscReal )

    @chk ccall(
               (:SNESMSSetDamping, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscReal),
               snes, damping,
              )


	return nothing
end 

"""
	SNESNGSSetTolerances(petsclib::PetscLibType,snes::PetscSNES, abstol::PetscReal, rtol::PetscReal, stol::PetscReal, maxit::PetscInt) 
Sets various parameters used in convergence tests for nonlinear Gauss

Logically Collective

Input Parameters:
- `snes`   - the `SNES` context
- `abstol` - absolute convergence tolerance
- `rtol`   - relative convergence tolerance
- `stol`   - convergence tolerance in terms of the norm of the change in the solution between steps,  || delta x || < stol*|| x ||
- `maxit`  - maximum number of iterations

Options Database Keys:
- `-snes_ngs_atol <abstol>` - Sets abstol
- `-snes_ngs_rtol <rtol>`   - Sets rtol
- `-snes_ngs_stol <stol>`   - Sets stol
- `-snes_max_it <maxit>`    - Sets maxit

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESNCG`

# External Links
$(_doc_external("Snes/SNESNGSSetTolerances"))
"""
function SNESNGSSetTolerances(petsclib::PetscLibType, snes::PetscSNES, abstol::PetscReal, rtol::PetscReal, stol::PetscReal, maxit::PetscInt) end

@for_petsc function SNESNGSSetTolerances(petsclib::$UnionPetscLib, snes::PetscSNES, abstol::$PetscReal, rtol::$PetscReal, stol::$PetscReal, maxit::$PetscInt )

    @chk ccall(
               (:SNESNGSSetTolerances, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscReal, $PetscReal, $PetscReal, $PetscInt),
               snes, abstol, rtol, stol, maxit,
              )


	return nothing
end 

"""
	SNESNGSGetTolerances(petsclib::PetscLibType,snes::PetscSNES, atol::PetscReal, rtol::PetscReal, stol::PetscReal, maxit::PetscInt) 
Gets various parameters used in convergence tests for nonlinear Gauss

Not Collective

Input Parameters:
- `snes`  - the `SNES` context
- `atol`  - absolute convergence tolerance
- `rtol`  - relative convergence tolerance
- `stol`  - convergence tolerance in terms of the norm
of the change in the solution between steps
- `maxit` - maximum number of iterations

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESNCG`, `SNESSetTolerances()`

# External Links
$(_doc_external("Snes/SNESNGSGetTolerances"))
"""
function SNESNGSGetTolerances(petsclib::PetscLibType, snes::PetscSNES, atol::PetscReal, rtol::PetscReal, stol::PetscReal, maxit::PetscInt) end

@for_petsc function SNESNGSGetTolerances(petsclib::$UnionPetscLib, snes::PetscSNES, atol::$PetscReal, rtol::$PetscReal, stol::$PetscReal, maxit::$PetscInt )

    @chk ccall(
               (:SNESNGSGetTolerances, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscInt}),
               snes, atol, rtol, stol, maxit,
              )


	return nothing
end 

"""
	SNESNGSSetSweeps(petsclib::PetscLibType,snes::PetscSNES, sweeps::PetscInt) 
Sets the number of sweeps of nonlinear GS to use in `SNESNCG`

Logically Collective

Input Parameters:
- `snes`   - the `SNES` context
- `sweeps` - the number of sweeps of nonlinear GS to perform.

Options Database Key:
- `-snes_ngs_sweeps <n>` - Number of sweeps of nonlinear GS to apply

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESNCG`, `SNESSetNGS()`, `SNESGetNGS()`, `SNESSetNPC()`, `SNESNGSGetSweeps()`

# External Links
$(_doc_external("Snes/SNESNGSSetSweeps"))
"""
function SNESNGSSetSweeps(petsclib::PetscLibType, snes::PetscSNES, sweeps::PetscInt) end

@for_petsc function SNESNGSSetSweeps(petsclib::$UnionPetscLib, snes::PetscSNES, sweeps::$PetscInt )

    @chk ccall(
               (:SNESNGSSetSweeps, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt),
               snes, sweeps,
              )


	return nothing
end 

"""
	sweeps::PetscInt = SNESNGSGetSweeps(petsclib::PetscLibType,snes::PetscSNES) 
Gets the number of sweeps nonlinear GS will use in `SNESNCG`

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `sweeps` - the number of sweeps of nonlinear GS to perform.

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESNCG`, `SNESSetNGS()`, `SNESGetNGS()`, `SNESSetNPC()`, `SNESNGSSetSweeps()`

# External Links
$(_doc_external("Snes/SNESNGSGetSweeps"))
"""
function SNESNGSGetSweeps(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESNGSGetSweeps(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESNewtonALSetCorrectionType(petsclib::PetscLibType,snes::PetscSNES, ctype::SNESNewtonALCorrectionType) 
Set the type of correction to use in the arc

Logically Collective

Input Parameters:
- `snes`  - the nonlinear solver object
- `ctype` - the type of correction to use

Options Database Key:
- `-snes_newtonal_correction_type <type>` - Set the type of correction to use; use -help for a list of available types

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESNEWTONAL`, `SNESNewtonALCorrectionType`

# External Links
$(_doc_external("Snes/SNESNewtonALSetCorrectionType"))
"""
function SNESNewtonALSetCorrectionType(petsclib::PetscLibType, snes::PetscSNES, ctype::SNESNewtonALCorrectionType) end

@for_petsc function SNESNewtonALSetCorrectionType(petsclib::$UnionPetscLib, snes::PetscSNES, ctype::SNESNewtonALCorrectionType )

    @chk ccall(
               (:SNESNewtonALSetCorrectionType, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESNewtonALCorrectionType),
               snes, ctype,
              )


	return nothing
end 

"""
	SNESNewtonALSetFunction(petsclib::PetscLibType,snes::PetscSNES, func::SNESFunctionFn, ctx::Cvoid) 
Sets a user function that is called at each function evaluation to
compute the tangent load vector for the arc-length continuation method.

Logically Collective

Input Parameters:
- `snes` - the nonlinear solver object
- `func` - [optional] tangent load function evaluation routine, see `SNESFunctionFn` for the calling sequence. `U` is the current solution vector, `Q` is the output tangent load vector
- `ctx`  - [optional] user-defined context for private data for the function evaluation routine (may be `NULL`)

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESNEWTONAL`, `SNESNewtonALGetFunction()`, `SNESNewtonALGetLoadParameter()`

# External Links
$(_doc_external("Snes/SNESNewtonALSetFunction"))
"""
function SNESNewtonALSetFunction(petsclib::PetscLibType, snes::PetscSNES, func::SNESFunctionFn, ctx::Cvoid) end

@for_petsc function SNESNewtonALSetFunction(petsclib::$UnionPetscLib, snes::PetscSNES, func::SNESFunctionFn, ctx::Cvoid )

    @chk ccall(
               (:SNESNewtonALSetFunction, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{SNESFunctionFn}, Ptr{Cvoid}),
               snes, func, ctx,
              )


	return nothing
end 

"""
	SNESNewtonALGetFunction(petsclib::PetscLibType,snes::PetscSNES, func::SNESFunctionFn, ctx::Cvoid) 
Get the user function and context set with `SNESNewtonALSetFunction`

Logically Collective

Input Parameters:
- `snes` - the nonlinear solver object
- `func` - [optional] tangent load function evaluation routine, see `SNESNewtonALSetFunction()` for the call sequence
- `ctx`  - [optional] user-defined context for private data for the function evaluation routine (may be `NULL`)

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESNEWTONAL`, `SNESNewtonALSetFunction()`

# External Links
$(_doc_external("Snes/SNESNewtonALGetFunction"))
"""
function SNESNewtonALGetFunction(petsclib::PetscLibType, snes::PetscSNES, func::SNESFunctionFn, ctx::Cvoid) end

@for_petsc function SNESNewtonALGetFunction(petsclib::$UnionPetscLib, snes::PetscSNES, func::SNESFunctionFn, ctx::Cvoid )

    @chk ccall(
               (:SNESNewtonALGetFunction, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESFunctionFn, Cvoid),
               snes, func, ctx,
              )


	return nothing
end 

"""
	lambda::PetscReal = SNESNewtonALGetLoadParameter(petsclib::PetscLibType,snes::PetscSNES) 
Get the value of the load parameter `lambda` for the arc

Logically Collective

Input Parameter:
- `snes` - the nonlinear solver object

Output Parameter:
- `lambda` - the arc-length parameter

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESNEWTONAL`, `SNESNewtonALSetFunction()`

# External Links
$(_doc_external("Snes/SNESNewtonALGetLoadParameter"))
"""
function SNESNewtonALGetLoadParameter(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESNewtonALGetLoadParameter(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESNewtonALComputeFunction(petsclib::PetscLibType,snes::PetscSNES, X::PetscVec, Q::PetscVec) 
Calls the function that has been set with `SNESNewtonALSetFunction()`.

Collective

Input Parameters:
- `snes` - the `SNES` context
- `X`    - input vector

Output Parameter:
- `Q` - tangent load vector, as set by `SNESNewtonALSetFunction()`

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESNewtonALSetFunction()`, `SNESNewtonALGetFunction()`

# External Links
$(_doc_external("Snes/SNESNewtonALComputeFunction"))
"""
function SNESNewtonALComputeFunction(petsclib::PetscLibType, snes::PetscSNES, X::PetscVec, Q::PetscVec) end

@for_petsc function SNESNewtonALComputeFunction(petsclib::$UnionPetscLib, snes::PetscSNES, X::PetscVec, Q::PetscVec )

    @chk ccall(
               (:SNESNewtonALComputeFunction, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec),
               snes, X, Q,
              )


	return nothing
end 

"""
	SNESPythonSetType(petsclib::PetscLibType,snes::PetscSNES, pyname::String) 
Initialize a `SNES` object implemented in Python.

Collective

Input Parameters:
- `snes`  - the nonlinear solver (`SNES`) context.
- `pyname`  - full dotted Python name [package].module[.{class|function}]

Options Database Key:
- `-snes_python_type <pyname>`  - python class

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESCreate()`, `SNESSetType()`, `SNESPYTHON`, `PetscPythonInitialize()`, `SNESPythonGetType()`

# External Links
$(_doc_external("Snes/SNESPythonSetType"))
"""
function SNESPythonSetType(petsclib::PetscLibType, snes::PetscSNES, pyname::String) end

@for_petsc function SNESPythonSetType(petsclib::$UnionPetscLib, snes::PetscSNES, pyname::String )

    @chk ccall(
               (:SNESPythonSetType, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Cchar}),
               snes, pyname,
              )


	return nothing
end 

"""
	pyname::String = SNESPythonGetType(petsclib::PetscLibType,snes::PetscSNES) 
Get the type of a `SNES` object implemented in Python set with `SNESPythonSetType()`

Not Collective

Input Parameter:
- `snes`  - the nonlinear solver (`SNES`) context.

Output Parameter:
- `pyname`  - full dotted Python name [package].module[.{class|function}]

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESCreate()`, `SNESSetType()`, `SNESPYTHON`, `PetscPythonInitialize()`, `SNESPythonSetType()`

# External Links
$(_doc_external("Snes/SNESPythonGetType"))
"""
function SNESPythonGetType(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESPythonGetType(petsclib::$UnionPetscLib, snes::PetscSNES )
	pyname_ = Ref{Ptr{Cchar}}()

    @chk ccall(
               (:SNESPythonGetType, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Ptr{Cchar}}),
               snes, pyname_,
              )

	pyname = unsafe_wrap(Array, pyname_[], VecGetLocalSize(petsclib, x); own = false)

	return pyname
end 

"""
	SNESShellSetSolve(petsclib::PetscLibType,snes::PetscSNES, solve::external) 
Sets routine to apply as solver to a `SNESSHELL` `SNES` object

Logically Collective

Input Parameters:
- `snes`  - the `SNES` nonlinear solver context
- `solve` - the application-provided solver routine

Calling sequence of `apply`:
- `snes` - the preconditioner, get the application context with `SNESShellGetContext()` provided with `SNESShellSetContext()`
- `xout` - solution vector

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESSHELL`, `SNESShellSetContext()`, `SNESShellGetContext()`

# External Links
$(_doc_external("Snes/SNESShellSetSolve"))
"""
function SNESShellSetSolve(petsclib::PetscLibType, snes::PetscSNES, solve::external) end

@for_petsc function SNESShellSetSolve(petsclib::$UnionPetscLib, snes::PetscSNES, solve::external )

    @chk ccall(
               (:SNESShellSetSolve, $petsc_library),
               PetscErrorCode,
               (CSNES, external),
               snes, solve,
              )


	return nothing
end 

"""
	SNESShellGetContext(petsclib::PetscLibType,snes::PetscSNES, ctx::Cvoid) 
Returns the user

Not Collective

Input Parameter:
- `snes` - should have been created with `SNESSetType`(snes,`SNESSHELL`);

Output Parameter:
- `ctx` - the user provided context

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESSHELL`, `SNESCreateShell()`, `SNESShellSetContext()`

# External Links
$(_doc_external("Snes/SNESShellGetContext"))
"""
function SNESShellGetContext(petsclib::PetscLibType, snes::PetscSNES, ctx::Cvoid) end

@for_petsc function SNESShellGetContext(petsclib::$UnionPetscLib, snes::PetscSNES, ctx::Cvoid )

    @chk ccall(
               (:SNESShellGetContext, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Cvoid}),
               snes, ctx,
              )


	return nothing
end 

"""
	SNESShellSetContext(petsclib::PetscLibType,snes::PetscSNES, ctx::Cvoid) 
sets the context for a `SNESSHELL`

Logically Collective

Input Parameters:
- `snes` - the `SNESSHELL`
- `ctx`  - the context

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESSHELL`, `SNESCreateShell()`, `SNESShellGetContext()`

# External Links
$(_doc_external("Snes/SNESShellSetContext"))
"""
function SNESShellSetContext(petsclib::PetscLibType, snes::PetscSNES, ctx::Cvoid) end

@for_petsc function SNESShellSetContext(petsclib::$UnionPetscLib, snes::PetscSNES, ctx::Cvoid )

    @chk ccall(
               (:SNESShellSetContext, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Cvoid}),
               snes, ctx,
              )


	return nothing
end 

"""
	SNESNASMSetType(petsclib::PetscLibType,snes::PetscSNES, type::PCASMType) 
Set the type of subdomain update used for the nonlinear additive Schwarz solver `SNESNASM`

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `type` - the type of update, `PC_ASM_BASIC` or `PC_ASM_RESTRICT`

Options Database Key:
- `-snes_nasm_type <basic,restrict>` - type of subdomain update used

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESNASM`, `SNESNASMGetType()`, `PCASMSetType()`, `PC_ASM_BASIC`, `PC_ASM_RESTRICT`, `PCASMType`

# External Links
$(_doc_external("Snes/SNESNASMSetType"))
"""
function SNESNASMSetType(petsclib::PetscLibType, snes::PetscSNES, type::PCASMType) end

@for_petsc function SNESNASMSetType(petsclib::$UnionPetscLib, snes::PetscSNES, type::PCASMType )

    @chk ccall(
               (:SNESNASMSetType, $petsc_library),
               PetscErrorCode,
               (CSNES, PCASMType),
               snes, type,
              )


	return nothing
end 

"""
	type::PCASMType = SNESNASMGetType(petsclib::PetscLibType,snes::PetscSNES) 
Get the type of subdomain update used for the nonlinear additive Schwarz solver `SNESNASM`

Logically Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `type` - the type of update

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESNASM`, `SNESNASMSetType()`, `PCASMGetType()`, `PC_ASM_BASIC`, `PC_ASM_RESTRICT`, `PCASMType`

# External Links
$(_doc_external("Snes/SNESNASMGetType"))
"""
function SNESNASMGetType(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESNASMGetType(petsclib::$UnionPetscLib, snes::PetscSNES )
	type_ = Ref{PCASMType}()

    @chk ccall(
               (:SNESNASMGetType, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{PCASMType}),
               snes, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	SNESNASMSetSubdomains(petsclib::PetscLibType,snes::PetscSNES, n::PetscInt, subsnes::Vector{PetscSNES}, iscatter::Vector{VecScatter}, oscatter::Vector{VecScatter}, gscatter::Vector{VecScatter}) 
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

-seealso: [](ch_snes), `SNES`, `SNESNASM`, `SNESNASMGetSubdomains()`

# External Links
$(_doc_external("Snes/SNESNASMSetSubdomains"))
"""
function SNESNASMSetSubdomains(petsclib::PetscLibType, snes::PetscSNES, n::PetscInt, subsnes::Vector{PetscSNES}, iscatter::Vector{VecScatter}, oscatter::Vector{VecScatter}, gscatter::Vector{VecScatter}) end

@for_petsc function SNESNASMSetSubdomains(petsclib::$UnionPetscLib, snes::PetscSNES, n::$PetscInt, subsnes::Vector{PetscSNES}, iscatter::Vector{VecScatter}, oscatter::Vector{VecScatter}, gscatter::Vector{VecScatter} )

    @chk ccall(
               (:SNESNASMSetSubdomains, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, Ptr{CSNES}, Ptr{VecScatter}, Ptr{VecScatter}, Ptr{VecScatter}),
               snes, n, subsnes, iscatter, oscatter, gscatter,
              )


	return nothing
end 

"""
	n::PetscInt = SNESNASMGetSubdomains(petsclib::PetscLibType,snes::PetscSNES, subsnes::Vector{PetscSNES}, iscatter::Vector{VecScatter}, oscatter::Vector{VecScatter}, gscatter::Vector{VecScatter}) 
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

-seealso: [](ch_snes), `SNES`, `SNESNASM`, `SNESNASMSetSubdomains()`

# External Links
$(_doc_external("Snes/SNESNASMGetSubdomains"))
"""
function SNESNASMGetSubdomains(petsclib::PetscLibType, snes::PetscSNES, subsnes::Vector{PetscSNES}, iscatter::Vector{VecScatter}, oscatter::Vector{VecScatter}, gscatter::Vector{VecScatter}) end

@for_petsc function SNESNASMGetSubdomains(petsclib::$UnionPetscLib, snes::PetscSNES, subsnes::Vector{PetscSNES}, iscatter::Vector{VecScatter}, oscatter::Vector{VecScatter}, gscatter::Vector{VecScatter} )
	n_ = Ref{$PetscInt}()
	subsnes_ = Ref(pointer(subsnes))
	iscatter_ = Ref(pointer(iscatter))
	oscatter_ = Ref(pointer(oscatter))
	gscatter_ = Ref(pointer(gscatter))

    @chk ccall(
               (:SNESNASMGetSubdomains, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscInt}, Ptr{Ptr{CSNES}}, Ptr{Ptr{VecScatter}}, Ptr{Ptr{VecScatter}}, Ptr{Ptr{VecScatter}}),
               snes, n_, subsnes_, iscatter_, oscatter_, gscatter_,
              )

	n = n_[]

	return n
end 

"""
	n::PetscInt = SNESNASMGetSubdomainVecs(petsclib::PetscLibType,snes::PetscSNES, x::Vector{PetscVec}, y::Vector{PetscVec}, b::Vector{PetscVec}, xl::Vector{PetscVec}) 
Get the processor

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

-seealso: [](ch_snes), `SNES`, `SNESNASM`, `SNESNASMGetSubdomains()`

# External Links
$(_doc_external("Snes/SNESNASMGetSubdomainVecs"))
"""
function SNESNASMGetSubdomainVecs(petsclib::PetscLibType, snes::PetscSNES, x::Vector{PetscVec}, y::Vector{PetscVec}, b::Vector{PetscVec}, xl::Vector{PetscVec}) end

@for_petsc function SNESNASMGetSubdomainVecs(petsclib::$UnionPetscLib, snes::PetscSNES, x::Vector{PetscVec}, y::Vector{PetscVec}, b::Vector{PetscVec}, xl::Vector{PetscVec} )
	n_ = Ref{$PetscInt}()
	x_ = Ref(pointer(x))
	y_ = Ref(pointer(y))
	b_ = Ref(pointer(b))
	xl_ = Ref(pointer(xl))

    @chk ccall(
               (:SNESNASMGetSubdomainVecs, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscInt}, Ptr{Ptr{CVec}}, Ptr{Ptr{CVec}}, Ptr{Ptr{CVec}}, Ptr{Ptr{CVec}}),
               snes, n_, x_, y_, b_, xl_,
              )

	n = n_[]

	return n
end 

"""
	SNESNASMSetComputeFinalJacobian(petsclib::PetscLibType,snes::PetscSNES, flg::PetscBool) 
Schedules the computation of the global and subdomain Jacobians upon convergence for the
nonlinear additive Schwarz solver

Collective

Input Parameters:
- `snes` - the SNES context
- `flg`  - `PETSC_TRUE` to compute the Jacobians

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESNASM`, `SNESNASMGetSubdomains()`

# External Links
$(_doc_external("Snes/SNESNASMSetComputeFinalJacobian"))
"""
function SNESNASMSetComputeFinalJacobian(petsclib::PetscLibType, snes::PetscSNES, flg::PetscBool) end

@for_petsc function SNESNASMSetComputeFinalJacobian(petsclib::$UnionPetscLib, snes::PetscSNES, flg::PetscBool )

    @chk ccall(
               (:SNESNASMSetComputeFinalJacobian, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, flg,
              )


	return nothing
end 

"""
	SNESNASMSetDamping(petsclib::PetscLibType,snes::PetscSNES, dmp::PetscReal) 
Sets the update damping for `SNESNASM` the nonlinear additive Schwarz solver

Logically Collective

Input Parameters:
- `snes` - the `SNES` context
- `dmp`  - damping

Options Database Key:
- `-snes_nasm_damping <dmp>` - the new solution is obtained as old solution plus `dmp` times (sum of the solutions on the subdomains)

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESNASM`, `SNESNASMGetDamping()`

# External Links
$(_doc_external("Snes/SNESNASMSetDamping"))
"""
function SNESNASMSetDamping(petsclib::PetscLibType, snes::PetscSNES, dmp::PetscReal) end

@for_petsc function SNESNASMSetDamping(petsclib::$UnionPetscLib, snes::PetscSNES, dmp::$PetscReal )

    @chk ccall(
               (:SNESNASMSetDamping, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscReal),
               snes, dmp,
              )


	return nothing
end 

"""
	dmp::PetscReal = SNESNASMGetDamping(petsclib::PetscLibType,snes::PetscSNES) 
Gets the update damping for `SNESNASM` the nonlinear additive Schwarz solver

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `dmp` - damping

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESNASM`, `SNESNASMSetDamping()`

# External Links
$(_doc_external("Snes/SNESNASMGetDamping"))
"""
function SNESNASMGetDamping(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESNASMGetDamping(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESNASMGetSNES(petsclib::PetscLibType,snes::PetscSNES, i::PetscInt, subsnes::PetscSNES) 
Gets a subsolver

Not Collective

Input Parameters:
- `snes` - the `SNES` context
- `i`    - the number of the subsnes to get

Output Parameter:
- `subsnes` - the subsolver context

Level: intermediate

-seealso: [](ch_snes), `SNESNASM`, `SNESNASMGetNumber()`

# External Links
$(_doc_external("Snes/SNESNASMGetSNES"))
"""
function SNESNASMGetSNES(petsclib::PetscLibType, snes::PetscSNES, i::PetscInt, subsnes::PetscSNES) end

@for_petsc function SNESNASMGetSNES(petsclib::$UnionPetscLib, snes::PetscSNES, i::$PetscInt, subsnes::PetscSNES )
	subsnes_ = Ref(subsnes.ptr)

    @chk ccall(
               (:SNESNASMGetSNES, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, Ptr{CSNES}),
               snes, i, subsnes_,
              )

	subsnes.ptr = C_NULL

	return nothing
end 

"""
	n::PetscInt = SNESNASMGetNumber(petsclib::PetscLibType,snes::PetscSNES) 
Gets number of subsolvers

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `n` - the number of subsolvers

Level: intermediate

-seealso: [](ch_snes), `SNESNASM`, `SNESNASMGetSNES()`

# External Links
$(_doc_external("Snes/SNESNASMGetNumber"))
"""
function SNESNASMGetNumber(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESNASMGetNumber(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESNASMSetWeight(petsclib::PetscLibType,snes::PetscSNES, weight::PetscVec) 
Sets weight to use when adding overlapping updates

Collective

Input Parameters:
- `snes`   - the `SNES` context
- `weight` - the weights to use (typically 1/N for each dof, where N is the number of patches it appears in)

Level: intermediate

-seealso: [](ch_snes), `SNESNASM`

# External Links
$(_doc_external("Snes/SNESNASMSetWeight"))
"""
function SNESNASMSetWeight(petsclib::PetscLibType, snes::PetscSNES, weight::PetscVec) end

@for_petsc function SNESNASMSetWeight(petsclib::$UnionPetscLib, snes::PetscSNES, weight::PetscVec )

    @chk ccall(
               (:SNESNASMSetWeight, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec),
               snes, weight,
              )


	return nothing
end 

"""
	SNESQNSetRestartType(petsclib::PetscLibType,snes::PetscSNES, rtype::SNESQNRestartType) 
Sets the restart type for `SNESQN`.

Logically Collective

Input Parameters:
- `snes`  - the iterative context
- `rtype` - restart type, see `SNESQNRestartType`

Options Database Keys:
- `-snes_qn_restart_type <powell,periodic,none>` - set the restart type
- `-snes_qn_m <m>`                               - sets the number of stored updates and the restart period for periodic

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESQN`, `SNESQNRestartType`, `SNES_QN_RESTART_NONE`, `SNES_QN_RESTART_POWELL`, `SNES_QN_RESTART_PERIODIC`,
`SNESQNType`, `SNESQNScaleType`

# External Links
$(_doc_external("Snes/SNESQNSetRestartType"))
"""
function SNESQNSetRestartType(petsclib::PetscLibType, snes::PetscSNES, rtype::SNESQNRestartType) end

@for_petsc function SNESQNSetRestartType(petsclib::$UnionPetscLib, snes::PetscSNES, rtype::SNESQNRestartType )

    @chk ccall(
               (:SNESQNSetRestartType, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESQNRestartType),
               snes, rtype,
              )


	return nothing
end 

"""
	SNESQNSetScaleType(petsclib::PetscLibType,snes::PetscSNES, stype::SNESQNScaleType) 
Sets the scaling type for the inner inverse Jacobian in `SNESQN`.

Logically Collective

Input Parameters:
- `snes`  - the nonlinear solver context
- `stype` - scale type, see `SNESQNScaleType`

Options Database Key:
- `-snes_qn_scale_type <diagonal,none,scalar,jacobian>` - Scaling type

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESQN`, `SNESLineSearch`, `SNESQNScaleType`, `SNESSetJacobian()`, `SNESQNType`, `SNESQNRestartType`

# External Links
$(_doc_external("Snes/SNESQNSetScaleType"))
"""
function SNESQNSetScaleType(petsclib::PetscLibType, snes::PetscSNES, stype::SNESQNScaleType) end

@for_petsc function SNESQNSetScaleType(petsclib::$UnionPetscLib, snes::PetscSNES, stype::SNESQNScaleType )

    @chk ccall(
               (:SNESQNSetScaleType, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESQNScaleType),
               snes, stype,
              )


	return nothing
end 

"""
	SNESQNSetType(petsclib::PetscLibType,snes::PetscSNES, qtype::SNESQNType) 
Sets the quasi

Logically Collective

Input Parameters:
- `snes`  - the iterative context
- `qtype` - variant type, see `SNESQNType`

Options Database Key:
- `-snes_qn_type <lbfgs,broyden,badbroyden>` - quasi-Newton type

Level: intermediate

-seealso: [](ch_snes), `SNESQN`, `SNES_QN_LBFGS`, `SNES_QN_BROYDEN`, `SNES_QN_BADBROYDEN`, `SNESQNType`,  `SNESQNScaleType`, `TAOLMVM`, `TAOBLMVM`

# External Links
$(_doc_external("Snes/SNESQNSetType"))
"""
function SNESQNSetType(petsclib::PetscLibType, snes::PetscSNES, qtype::SNESQNType) end

@for_petsc function SNESQNSetType(petsclib::$UnionPetscLib, snes::PetscSNES, qtype::SNESQNType )

    @chk ccall(
               (:SNESQNSetType, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESQNType),
               snes, qtype,
              )


	return nothing
end 

"""
	SNESCompositeSetType(petsclib::PetscLibType,snes::PetscSNES, type::SNESCompositeType) 
Sets the type of composite preconditioner.

Logically Collective

Input Parameters:
- `snes` - the preconditioner context
- `type` - `SNES_COMPOSITE_ADDITIVE` (default), `SNES_COMPOSITE_MULTIPLICATIVE`, or `SNES_COMPOSITE_ADDITIVEOPTIMAL`

Options Database Key:
- `-snes_composite_type <type: one of multiplicative, additive, additiveoptimal>` - Sets composite preconditioner type

Level: developer

-seealso: [](ch_snes), `SNES_COMPOSITE_ADDITIVE`, `SNES_COMPOSITE_MULTIPLICATIVE`, `SNESCompositeType`, `SNESCOMPOSITE`, `SNES_COMPOSITE_ADDITIVEOPTIMAL`,
`PCCompositeType`

# External Links
$(_doc_external("Snes/SNESCompositeSetType"))
"""
function SNESCompositeSetType(petsclib::PetscLibType, snes::PetscSNES, type::SNESCompositeType) end

@for_petsc function SNESCompositeSetType(petsclib::$UnionPetscLib, snes::PetscSNES, type::SNESCompositeType )

    @chk ccall(
               (:SNESCompositeSetType, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESCompositeType),
               snes, type,
              )


	return nothing
end 

"""
	SNESCompositeAddSNES(petsclib::PetscLibType,snes::PetscSNES, type::SNESType) 
Adds another `SNES` to the `SNESCOMPOSITE`

Collective

Input Parameters:
- `snes` - the `SNES` context of type `SNESCOMPOSITE`
- `type` - the `SNESType` of the new solver

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESCOMPOSITE`, `SNESCompositeGetSNES()`

# External Links
$(_doc_external("Snes/SNESCompositeAddSNES"))
"""
function SNESCompositeAddSNES(petsclib::PetscLibType, snes::PetscSNES, type::SNESType) end

@for_petsc function SNESCompositeAddSNES(petsclib::$UnionPetscLib, snes::PetscSNES, type::SNESType )

    @chk ccall(
               (:SNESCompositeAddSNES, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESType),
               snes, type,
              )


	return nothing
end 

"""
	SNESCompositeGetSNES(petsclib::PetscLibType,snes::PetscSNES, n::PetscInt, subsnes::PetscSNES) 
Gets one of the `SNES` objects in the `SNES` of `SNESType` `SNESCOMPOSITE`

Not Collective

Input Parameters:
- `snes` - the `SNES` context
- `n`    - the number of the composed `SNES` requested

Output Parameter:
- `subsnes` - the `SNES` requested

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESCOMPOSITE`, `SNESCompositeAddSNES()`, `SNESCompositeGetNumber()`

# External Links
$(_doc_external("Snes/SNESCompositeGetSNES"))
"""
function SNESCompositeGetSNES(petsclib::PetscLibType, snes::PetscSNES, n::PetscInt, subsnes::PetscSNES) end

@for_petsc function SNESCompositeGetSNES(petsclib::$UnionPetscLib, snes::PetscSNES, n::$PetscInt, subsnes::PetscSNES )
	subsnes_ = Ref(subsnes.ptr)

    @chk ccall(
               (:SNESCompositeGetSNES, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, Ptr{CSNES}),
               snes, n, subsnes_,
              )

	subsnes.ptr = C_NULL

	return nothing
end 

"""
	n::PetscInt = SNESCompositeGetNumber(petsclib::PetscLibType,snes::PetscSNES) 
Get the number of subsolvers in the `SNESCOMPOSITE`

Logically Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameter:
- `n` - the number of subsolvers

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESCOMPOSITE`, `SNESCompositeAddSNES()`, `SNESCompositeGetSNES()`

# External Links
$(_doc_external("Snes/SNESCompositeGetNumber"))
"""
function SNESCompositeGetNumber(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESCompositeGetNumber(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESCompositeSetDamping(petsclib::PetscLibType,snes::PetscSNES, n::PetscInt, dmp::PetscReal) 
Sets the damping of a subsolver when using `SNES_COMPOSITE_ADDITIVE` with a `SNES` of `SNESType` `SNESCOMPOSITE`

Not Collective

Input Parameters:
- `snes` - the `SNES` context
- `n`    - the number of the sub-`SNES` object requested
- `dmp`  - the damping

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESCOMPOSITE`, `SNESCompositeAddSNES()`, `SNESCompositeGetSNES()`,
`SNES_COMPOSITE_ADDITIVE`, `SNES_COMPOSITE_MULTIPLICATIVE`, `SNESCompositeType`, `SNESCompositeSetType()`

# External Links
$(_doc_external("Snes/SNESCompositeSetDamping"))
"""
function SNESCompositeSetDamping(petsclib::PetscLibType, snes::PetscSNES, n::PetscInt, dmp::PetscReal) end

@for_petsc function SNESCompositeSetDamping(petsclib::$UnionPetscLib, snes::PetscSNES, n::$PetscInt, dmp::$PetscReal )

    @chk ccall(
               (:SNESCompositeSetDamping, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal),
               snes, n, dmp,
              )


	return nothing
end 

"""
	rho_flag::PetscBool = SNESNewtonTRDCGetRhoFlag(petsclib::PetscLibType,snes::PetscSNES) 
Get whether the current solution update is within the trust

Logically Collective

Input Parameter:
- `snes` - the nonlinear solver object

Output Parameter:
- `rho_flag` - `PETSC_FALSE` or `PETSC_TRUE`

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESNEWTONTRDC`, `SNESNewtonTRDCPreCheck()`, `SNESNewtonTRDCGetPreCheck()`, `SNESNewtonTRDCSetPreCheck()`,
`SNESNewtonTRDCSetPostCheck()`, `SNESNewtonTRDCGetPostCheck()`

# External Links
$(_doc_external("Snes/SNESNewtonTRDCGetRhoFlag"))
"""
function SNESNewtonTRDCGetRhoFlag(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESNewtonTRDCGetRhoFlag(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESNewtonTRDCSetPreCheck(petsclib::PetscLibType,snes::PetscSNES, func::external, ctx::Cvoid) 
Sets a user function that is called before the search step has been determined.
Allows the user a chance to change or override the trust region decision.

Logically Collective

Input Parameters:
- `snes` - the nonlinear solver object
- `func` - [optional] function evaluation routine, for the calling sequence see `SNESNewtonTRDCPreCheck()`
- `ctx`  - [optional] user-defined context for private data for the function evaluation routine (may be `NULL`)

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESNEWTONTRDC`, `SNESNewtonTRDCPreCheck()`, `SNESNewtonTRDCGetPreCheck()`, `SNESNewtonTRDCSetPostCheck()`, `SNESNewtonTRDCGetPostCheck()`,
`SNESNewtonTRDCGetRhoFlag()`

# External Links
$(_doc_external("Snes/SNESNewtonTRDCSetPreCheck"))
"""
function SNESNewtonTRDCSetPreCheck(petsclib::PetscLibType, snes::PetscSNES, func::external, ctx::Cvoid) end

@for_petsc function SNESNewtonTRDCSetPreCheck(petsclib::$UnionPetscLib, snes::PetscSNES, func::external, ctx::Cvoid )

    @chk ccall(
               (:SNESNewtonTRDCSetPreCheck, $petsc_library),
               PetscErrorCode,
               (CSNES, external, Ptr{Cvoid}),
               snes, func, ctx,
              )


	return nothing
end 

"""
	SNESNewtonTRDCSetPostCheck(petsclib::PetscLibType,snes::PetscSNES, func::external, ctx::Cvoid) 
Sets a user function that is called after the search step has been determined but before the next
function evaluation. Allows the user a chance to change or override the decision of the line search routine

Logically Collective

Input Parameters:
- `snes` - the nonlinear solver object
- `func` - [optional] function evaluation routine, for the calling sequence see `SNESNewtonTRDCPostCheck()`
- `ctx`  - [optional] user-defined context for private data for the function evaluation routine (may be `NULL`)

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESNEWTONTRDC`, `SNESNewtonTRDCPostCheck()`, `SNESNewtonTRDCGetPostCheck()`, `SNESNewtonTRDCSetPreCheck()`, `SNESNewtonTRDCGetPreCheck()`

# External Links
$(_doc_external("Snes/SNESNewtonTRDCSetPostCheck"))
"""
function SNESNewtonTRDCSetPostCheck(petsclib::PetscLibType, snes::PetscSNES, func::external, ctx::Cvoid) end

@for_petsc function SNESNewtonTRDCSetPostCheck(petsclib::$UnionPetscLib, snes::PetscSNES, func::external, ctx::Cvoid )

    @chk ccall(
               (:SNESNewtonTRDCSetPostCheck, $petsc_library),
               PetscErrorCode,
               (CSNES, external, Ptr{Cvoid}),
               snes, func, ctx,
              )


	return nothing
end 

"""
	bs::PetscInt,nodesPerCell::PetscInt,subspaceOffsets::PetscInt,ghostBcNodes::PetscInt,globalBcNodes::PetscInt = SNESPatchSetDiscretisationInfo(petsclib::PetscLibType,snes::PetscSNES, nsubspaces::PetscInt, dms::PetscDM, cellNodeMap::PetscInt, numGhostBcs::PetscInt, numGlobalBcs::PetscInt) 

# External Links
$(_doc_external("Snes/SNESPatchSetDiscretisationInfo"))
"""
function SNESPatchSetDiscretisationInfo(petsclib::PetscLibType, snes::PetscSNES, nsubspaces::PetscInt, dms::PetscDM, cellNodeMap::PetscInt, numGhostBcs::PetscInt, numGlobalBcs::PetscInt) end

@for_petsc function SNESPatchSetDiscretisationInfo(petsclib::$UnionPetscLib, snes::PetscSNES, nsubspaces::$PetscInt, dms::PetscDM, cellNodeMap::$PetscInt, numGhostBcs::$PetscInt, numGlobalBcs::$PetscInt )
	dms_ = Ref(dms.ptr)
	bs_ = Ref{$PetscInt}()
	nodesPerCell_ = Ref{$PetscInt}()
	subspaceOffsets_ = Ref{$PetscInt}()
	ghostBcNodes_ = Ref{$PetscInt}()
	globalBcNodes_ = Ref{$PetscInt}()

    @chk ccall(
               (:SNESPatchSetDiscretisationInfo, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, Ptr{CDM}, Ptr{$PetscInt}, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}),
               snes, nsubspaces, dms_, bs_, nodesPerCell_, cellNodeMap, subspaceOffsets_, numGhostBcs, ghostBcNodes_, numGlobalBcs, globalBcNodes_,
              )

	dms.ptr = C_NULL
	bs = bs_[]
	nodesPerCell = nodesPerCell_[]
	subspaceOffsets = subspaceOffsets_[]
	ghostBcNodes = ghostBcNodes_[]
	globalBcNodes = globalBcNodes_[]

	return bs,nodesPerCell,subspaceOffsets,ghostBcNodes,globalBcNodes
end 

"""
	SNESPatchSetComputeOperator(petsclib::PetscLibType,snes::PetscSNES, func::external, ctx::Cvoid) 

# External Links
$(_doc_external("Snes/SNESPatchSetComputeOperator"))
"""
function SNESPatchSetComputeOperator(petsclib::PetscLibType, snes::PetscSNES, func::external, ctx::Cvoid) end

@for_petsc function SNESPatchSetComputeOperator(petsclib::$UnionPetscLib, snes::PetscSNES, func::external, ctx::Cvoid )

    @chk ccall(
               (:SNESPatchSetComputeOperator, $petsc_library),
               PetscErrorCode,
               (CSNES, external, Ptr{Cvoid}),
               snes, func, ctx,
              )


	return nothing
end 

"""
	SNESPatchSetComputeFunction(petsclib::PetscLibType,snes::PetscSNES, func::external, ctx::Cvoid) 

# External Links
$(_doc_external("Snes/SNESPatchSetComputeFunction"))
"""
function SNESPatchSetComputeFunction(petsclib::PetscLibType, snes::PetscSNES, func::external, ctx::Cvoid) end

@for_petsc function SNESPatchSetComputeFunction(petsclib::$UnionPetscLib, snes::PetscSNES, func::external, ctx::Cvoid )

    @chk ccall(
               (:SNESPatchSetComputeFunction, $petsc_library),
               PetscErrorCode,
               (CSNES, external, Ptr{Cvoid}),
               snes, func, ctx,
              )


	return nothing
end 

"""
	ctx::Cvoid = SNESPatchSetConstructType(petsclib::PetscLibType,snes::PetscSNES, ctype::PCPatchConstructType, func::external) 

# External Links
$(_doc_external("Snes/SNESPatchSetConstructType"))
"""
function SNESPatchSetConstructType(petsclib::PetscLibType, snes::PetscSNES, ctype::PCPatchConstructType, func::external) end

@for_petsc function SNESPatchSetConstructType(petsclib::$UnionPetscLib, snes::PetscSNES, ctype::PCPatchConstructType, func::external )
	ctx_ = Ref{Cvoid}()

    @chk ccall(
               (:SNESPatchSetConstructType, $petsc_library),
               PetscErrorCode,
               (CSNES, PCPatchConstructType, external, Ptr{Cvoid}),
               snes, ctype, func, ctx_,
              )

	ctx = ctx_[]

	return ctx
end 

"""
	SNESPatchSetCellNumbering(petsclib::PetscLibType,snes::PetscSNES, cellNumbering::PetscSection) 

# External Links
$(_doc_external("Snes/SNESPatchSetCellNumbering"))
"""
function SNESPatchSetCellNumbering(petsclib::PetscLibType, snes::PetscSNES, cellNumbering::PetscSection) end

@for_petsc function SNESPatchSetCellNumbering(petsclib::$UnionPetscLib, snes::PetscSNES, cellNumbering::PetscSection )

    @chk ccall(
               (:SNESPatchSetCellNumbering, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscSection),
               snes, cellNumbering,
              )


	return nothing
end 

"""
	SNESMultiblockSetFields(petsclib::PetscLibType,snes::PetscSNES, name::String, n::PetscInt, fields::PetscInt) 
Sets the fields for one particular block in a `SNESMULTIBLOCK` solver

Logically Collective

Input Parameters:
- `snes`   - the solver
- `name`   - name of this block, if `NULL` the number of the block is used
- `n`      - the number of fields in this block
- `fields` - the fields in this block

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESMULTIBLOCK`, `SNESMultiblockGetSubSNES()`, `SNESMultiblockSetBlockSize()`, `SNESMultiblockSetIS()`

# External Links
$(_doc_external("Snes/SNESMultiblockSetFields"))
"""
function SNESMultiblockSetFields(petsclib::PetscLibType, snes::PetscSNES, name::String, n::PetscInt, fields::PetscInt) end

@for_petsc function SNESMultiblockSetFields(petsclib::$UnionPetscLib, snes::PetscSNES, name::String, n::$PetscInt, fields::$PetscInt )

    @chk ccall(
               (:SNESMultiblockSetFields, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Cchar}, $PetscInt, Ptr{$PetscInt}),
               snes, name, n, fields,
              )


	return nothing
end 

"""
	SNESMultiblockSetIS(petsclib::PetscLibType,snes::PetscSNES, name::String, is::IS) 
Sets the global row indices for one particular block in a `SNESMULTIBLOCK` solver

Logically Collective

Input Parameters:
- `snes` - the solver context
- `name` - name of this block, if `NULL` the number of the block is used
- `is`   - the index set that defines the global row indices in this block

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESMULTIBLOCK`, `SNESMultiblockGetSubSNES()`, `SNESMultiblockSetBlockSize()`, `SNESMultiblockSetFields()`

# External Links
$(_doc_external("Snes/SNESMultiblockSetIS"))
"""
function SNESMultiblockSetIS(petsclib::PetscLibType, snes::PetscSNES, name::String, is::IS) end

@for_petsc function SNESMultiblockSetIS(petsclib::$UnionPetscLib, snes::PetscSNES, name::String, is::IS )

    @chk ccall(
               (:SNESMultiblockSetIS, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{Cchar}, IS),
               snes, name, is,
              )


	return nothing
end 

"""
	SNESMultiblockSetType(petsclib::PetscLibType,snes::PetscSNES, type::PCCompositeType) 
Sets the type of block combination used for a `SNESMULTIBLOCK` solver

Logically Collective

Input Parameters:
- `snes` - the solver context
- `type` - `PC_COMPOSITE_ADDITIVE`, `PC_COMPOSITE_MULTIPLICATIVE` (default), `PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE`

Options Database Key:
- `-snes_multiblock_type <type: one of multiplicative, additive, symmetric_multiplicative>` - Sets block combination type

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESMULTIBLOCK`, `PCCompositeSetType()`, `PC_COMPOSITE_ADDITIVE`, `PC_COMPOSITE_MULTIPLICATIVE`, `PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE`,
`PCCompositeType`, `SNESCOMPOSITE`, `SNESCompositeSetType()`

# External Links
$(_doc_external("Snes/SNESMultiblockSetType"))
"""
function SNESMultiblockSetType(petsclib::PetscLibType, snes::PetscSNES, type::PCCompositeType) end

@for_petsc function SNESMultiblockSetType(petsclib::$UnionPetscLib, snes::PetscSNES, type::PCCompositeType )

    @chk ccall(
               (:SNESMultiblockSetType, $petsc_library),
               PetscErrorCode,
               (CSNES, PCCompositeType),
               snes, type,
              )


	return nothing
end 

"""
	SNESMultiblockSetBlockSize(petsclib::PetscLibType,snes::PetscSNES, bs::PetscInt) 
Sets the block size for structured block division in a `SNESMULTIBLOCK` solver. If not set the matrix block size is used.

Logically Collective

Input Parameters:
- `snes` - the solver context
- `bs`   - the block size

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESMULTIBLOCK`, `SNESMultiblockGetSubSNES()`, `SNESMultiblockSetFields()`

# External Links
$(_doc_external("Snes/SNESMultiblockSetBlockSize"))
"""
function SNESMultiblockSetBlockSize(petsclib::PetscLibType, snes::PetscSNES, bs::PetscInt) end

@for_petsc function SNESMultiblockSetBlockSize(petsclib::$UnionPetscLib, snes::PetscSNES, bs::$PetscInt )

    @chk ccall(
               (:SNESMultiblockSetBlockSize, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt),
               snes, bs,
              )


	return nothing
end 

"""
	n::PetscInt = SNESMultiblockGetSubSNES(petsclib::PetscLibType,snes::PetscSNES, subsnes::Vector{PetscSNES}) 
Gets the `SNES` contexts for all blocks in a `SNESMULTIBLOCK` solver.

Not Collective but each `SNES` obtained is parallel

Input Parameter:
- `snes` - the solver context

Output Parameters:
- `n`       - the number of blocks
- `subsnes` - the array of `SNES` contexts

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESMULTIBLOCK`, `SNESMultiblockSetIS()`, `SNESMultiblockSetFields()`

# External Links
$(_doc_external("Snes/SNESMultiblockGetSubSNES"))
"""
function SNESMultiblockGetSubSNES(petsclib::PetscLibType, snes::PetscSNES, subsnes::Vector{PetscSNES}) end

@for_petsc function SNESMultiblockGetSubSNES(petsclib::$UnionPetscLib, snes::PetscSNES, subsnes::Vector{PetscSNES} )
	n_ = Ref{$PetscInt}()
	subsnes_ = Ref(pointer(subsnes))

    @chk ccall(
               (:SNESMultiblockGetSubSNES, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{$PetscInt}, Ptr{Ptr{CSNES}}),
               snes, n_, subsnes_,
              )

	n = n_[]

	return n
end 

"""
	SNESNCGSetType(petsclib::PetscLibType,snes::PetscSNES, btype::SNESNCGType) 
Sets the conjugate update type for nonlinear CG `SNESNCG`.

Logically Collective

Input Parameters:
- `snes`  - the iterative context
- `btype` - update type, see `SNESNCGType`

Options Database Key:
- `-snes_ncg_type <prp,fr,hs,dy,cd>` - strategy for selecting algorithm for computing beta

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESNCG`, `SNESNCGType`, `SNES_NCG_FR`, `SNES_NCG_PRP`, `SNES_NCG_HS`, `SNES_NCG_DY`, `SNES_NCG_CD`

# External Links
$(_doc_external("Snes/SNESNCGSetType"))
"""
function SNESNCGSetType(petsclib::PetscLibType, snes::PetscSNES, btype::SNESNCGType) end

@for_petsc function SNESNCGSetType(petsclib::$UnionPetscLib, snes::PetscSNES, btype::SNESNCGType )

    @chk ccall(
               (:SNESNCGSetType, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESNCGType),
               snes, btype,
              )


	return nothing
end 

"""
	Xcoarse::PetscVec = SNESFASCreateCoarseVec(petsclib::PetscLibType,snes::PetscSNES) 
create a `Vec` corresponding to a state vector on one level coarser than the current level

Collective

Input Parameter:
- `snes` - `SNESFAS` object

Output Parameter:
- `Xcoarse` - vector on level one coarser than the current level

Level: developer

-seealso: [](ch_snes), `SNESFASSetRestriction()`, `SNESFASRestrict()`, `SNESFAS`

# External Links
$(_doc_external("Snes/SNESFASCreateCoarseVec"))
"""
function SNESFASCreateCoarseVec(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESFASCreateCoarseVec(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESFASRestrict(petsclib::PetscLibType,fine::PetscSNES, Xfine::PetscVec, Xcoarse::PetscVec) 
restrict a `Vec` to the next coarser level

Collective

Input Parameters:
- `fine`  - `SNES` from which to restrict
- `Xfine` - vector to restrict

Output Parameter:
- `Xcoarse` - result of restriction

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASSetRestriction()`, `SNESFASSetInjection()`, `SNESFASCreateCoarseVec()`

# External Links
$(_doc_external("Snes/SNESFASRestrict"))
"""
function SNESFASRestrict(petsclib::PetscLibType, fine::PetscSNES, Xfine::PetscVec, Xcoarse::PetscVec) end

@for_petsc function SNESFASRestrict(petsclib::$UnionPetscLib, fine::PetscSNES, Xfine::PetscVec, Xcoarse::PetscVec )

    @chk ccall(
               (:SNESFASRestrict, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec),
               fine, Xfine, Xcoarse,
              )


	return nothing
end 

"""
	flg::PetscBool = SNESFASGetGalerkin(petsclib::PetscLibType,snes::PetscSNES) 
Gets if the coarse problems are formed by projection to the fine problem

Not Collective but the result would be the same on all MPI processes

Input Parameter:
- `snes` - the `SNESFAS` nonlinear solver context

Output Parameter:
- `flg` - `PETSC_TRUE` if the coarse problem is formed by projection

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASSetLevels()`, `SNESFASSetGalerkin()`

# External Links
$(_doc_external("Snes/SNESFASGetGalerkin"))
"""
function SNESFASGetGalerkin(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESFASGetGalerkin(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESFASSetGalerkin(petsclib::PetscLibType,snes::PetscSNES, flg::PetscBool) 
Sets coarse problems as formed by projection to the fine problem

Logically Collective

Input Parameters:
- `snes` - the `SNESFAS` nonlinear solver context
- `flg`  - `PETSC_TRUE` to use the projection process

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASSetLevels()`, `SNESFASGetGalerkin()`

# External Links
$(_doc_external("Snes/SNESFASSetGalerkin"))
"""
function SNESFASSetGalerkin(petsclib::PetscLibType, snes::PetscSNES, flg::PetscBool) end

@for_petsc function SNESFASSetGalerkin(petsclib::$UnionPetscLib, snes::PetscSNES, flg::PetscBool )

    @chk ccall(
               (:SNESFASSetGalerkin, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, flg,
              )


	return nothing
end 

"""
	SNESFASGalerkinFunctionDefault(petsclib::PetscLibType,snes::PetscSNES, X::PetscVec, F::PetscVec, ctx::Cvoid) 
Computes the Galerkin FAS function

Collective

Input Parameters:
- `snes` - the `SNESFAS` nonlinear solver context
- `X`    - input vector
- `ctx`  - the application context

Output Parameter:
- `F` - output vector

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASGetGalerkin()`, `SNESFASSetGalerkin()`

# External Links
$(_doc_external("Snes/SNESFASGalerkinFunctionDefault"))
"""
function SNESFASGalerkinFunctionDefault(petsclib::PetscLibType, snes::PetscSNES, X::PetscVec, F::PetscVec, ctx::Cvoid) end

@for_petsc function SNESFASGalerkinFunctionDefault(petsclib::$UnionPetscLib, snes::PetscSNES, X::PetscVec, F::PetscVec, ctx::Cvoid )

    @chk ccall(
               (:SNESFASGalerkinFunctionDefault, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec, Ptr{Cvoid}),
               snes, X, F, ctx,
              )


	return nothing
end 

"""
	SNESFASSetType(petsclib::PetscLibType,snes::PetscSNES, fastype::SNESFASType) 
Sets the update and correction type used for `SNESFAS`.

Logically Collective

Input Parameters:
- `snes`    - `SNESFAS` context
- `fastype` - `SNES_FAS_ADDITIVE`, `SNES_FAS_MULTIPLICATIVE`, `SNES_FAS_FULL`, or `SNES_FAS_KASKADE`

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `PCMGSetType()`, `SNESFASGetType()`, `SNES_FAS_ADDITIVE`, `SNES_FAS_MULTIPLICATIVE`, `SNES_FAS_FULL`, `SNES_FAS_KASKADE`

# External Links
$(_doc_external("Snes/SNESFASSetType"))
"""
function SNESFASSetType(petsclib::PetscLibType, snes::PetscSNES, fastype::SNESFASType) end

@for_petsc function SNESFASSetType(petsclib::$UnionPetscLib, snes::PetscSNES, fastype::SNESFASType )

    @chk ccall(
               (:SNESFASSetType, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESFASType),
               snes, fastype,
              )


	return nothing
end 

"""
	fastype::SNESFASType = SNESFASGetType(petsclib::PetscLibType,snes::PetscSNES) 
Gets the update and correction type used for `SNESFAS`.

Logically Collective

Input Parameter:
- `snes` - `SNESFAS` context

Output Parameter:
- `fastype` - `SNES_FAS_ADDITIVE`, `SNES_FAS_MULTIPLICATIVE`, `SNES_FAS_FULL`, or `SNES_FAS_KASKADE`

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `PCMGSetType()`, `SNESFASSetType()`, `SNES_FAS_ADDITIVE`, `SNES_FAS_MULTIPLICATIVE`, `SNES_FAS_FULL`, `SNES_FAS_KASKADE`

# External Links
$(_doc_external("Snes/SNESFASGetType"))
"""
function SNESFASGetType(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESFASGetType(petsclib::$UnionPetscLib, snes::PetscSNES )
	fastype_ = Ref{SNESFASType}()

    @chk ccall(
               (:SNESFASGetType, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{SNESFASType}),
               snes, fastype_,
              )

	fastype = unsafe_string(fastype_[])

	return fastype
end 

"""
	SNESFASSetLevels(petsclib::PetscLibType,snes::PetscSNES, levels::PetscInt, comms::MPI_Comm) 
Sets the number of levels to use with `SNESFAS`.
Must be called before any other `SNESFAS` routine.

Input Parameters:
- `snes`   - the `SNES` context of `SNESType` `SNESFAS`
- `levels` - the number of levels
- `comms`  - optional communicators for each level; this is to allow solving the coarser
problems on smaller sets of processors.

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASGetLevels()`

# External Links
$(_doc_external("Snes/SNESFASSetLevels"))
"""
function SNESFASSetLevels(petsclib::PetscLibType, snes::PetscSNES, levels::PetscInt, comms::MPI_Comm) end

@for_petsc function SNESFASSetLevels(petsclib::$UnionPetscLib, snes::PetscSNES, levels::$PetscInt, comms::MPI_Comm )

    @chk ccall(
               (:SNESFASSetLevels, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, Ptr{MPI_Comm}),
               snes, levels, comms,
              )


	return nothing
end 

"""
	levels::PetscInt = SNESFASGetLevels(petsclib::PetscLibType,snes::PetscSNES) 
Gets the number of levels in a `SNESFAS`, including fine and coarse grids

Input Parameter:
- `snes` - the `SNES` nonlinear solver context of `SNESType` `SNESFAS`

Output Parameter:
- `levels` - the number of levels

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASSetLevels()`, `PCMGGetLevels()`

# External Links
$(_doc_external("Snes/SNESFASGetLevels"))
"""
function SNESFASGetLevels(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESFASGetLevels(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESFASGetCycleSNES(petsclib::PetscLibType,snes::PetscSNES, level::PetscInt, lsnes::PetscSNES) 
Gets the `SNES` corresponding to a particular level of the `SNESFAS` hierarchy

Input Parameters:
- `snes`  - the `SNES` nonlinear multigrid context
- `level` - the level to get

Output Parameter:
- `lsnes` - the `SNES` for the requested level

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASSetLevels()`, `SNESFASGetLevels()`

# External Links
$(_doc_external("Snes/SNESFASGetCycleSNES"))
"""
function SNESFASGetCycleSNES(petsclib::PetscLibType, snes::PetscSNES, level::PetscInt, lsnes::PetscSNES) end

@for_petsc function SNESFASGetCycleSNES(petsclib::$UnionPetscLib, snes::PetscSNES, level::$PetscInt, lsnes::PetscSNES )
	lsnes_ = Ref(lsnes.ptr)

    @chk ccall(
               (:SNESFASGetCycleSNES, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, Ptr{CSNES}),
               snes, level, lsnes_,
              )

	lsnes.ptr = C_NULL

	return nothing
end 

"""
	SNESFASSetNumberSmoothUp(petsclib::PetscLibType,snes::PetscSNES, n::PetscInt) 
Sets the number of post
use on all levels.

Logically Collective

Input Parameters:
- `snes` - the `SNES` nonlinear multigrid context
- `n`    - the number of smoothing steps to use

Options Database Key:
- `-snes_fas_smoothup <n>` - Sets number of pre-smoothing steps

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASSetNumberSmoothDown()`

# External Links
$(_doc_external("Snes/SNESFASSetNumberSmoothUp"))
"""
function SNESFASSetNumberSmoothUp(petsclib::PetscLibType, snes::PetscSNES, n::PetscInt) end

@for_petsc function SNESFASSetNumberSmoothUp(petsclib::$UnionPetscLib, snes::PetscSNES, n::$PetscInt )

    @chk ccall(
               (:SNESFASSetNumberSmoothUp, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt),
               snes, n,
              )


	return nothing
end 

"""
	SNESFASSetNumberSmoothDown(petsclib::PetscLibType,snes::PetscSNES, n::PetscInt) 
Sets the number of pre
use on all levels.

Logically Collective

Input Parameters:
- `snes` - the `SNESFAS` nonlinear multigrid context
- `n`    - the number of smoothing steps to use

Options Database Key:
- `-snes_fas_smoothdown <n>` - Sets number of pre-smoothing steps

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASSetNumberSmoothUp()`

# External Links
$(_doc_external("Snes/SNESFASSetNumberSmoothDown"))
"""
function SNESFASSetNumberSmoothDown(petsclib::PetscLibType, snes::PetscSNES, n::PetscInt) end

@for_petsc function SNESFASSetNumberSmoothDown(petsclib::$UnionPetscLib, snes::PetscSNES, n::$PetscInt )

    @chk ccall(
               (:SNESFASSetNumberSmoothDown, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt),
               snes, n,
              )


	return nothing
end 

"""
	SNESFASSetContinuation(petsclib::PetscLibType,snes::PetscSNES, continuation::PetscBool) 
Sets the `SNESFAS` cycle to default to using exact Newton solves on the upsweep

Logically Collective

Input Parameters:
- `snes`         - the `SNESFAS` nonlinear multigrid context
- `continuation` - whether to use continuation

Options Database Key:
- `-snes_fas_continuation` - sets continuation to true

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASSetNumberSmoothUp()`

# External Links
$(_doc_external("Snes/SNESFASSetContinuation"))
"""
function SNESFASSetContinuation(petsclib::PetscLibType, snes::PetscSNES, continuation::PetscBool) end

@for_petsc function SNESFASSetContinuation(petsclib::$UnionPetscLib, snes::PetscSNES, continuation::PetscBool )

    @chk ccall(
               (:SNESFASSetContinuation, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, continuation,
              )


	return nothing
end 

"""
	SNESFASSetCycles(petsclib::PetscLibType,snes::PetscSNES, cycles::PetscInt) 
Sets the number of `SNESFAS` multigrid cycles to use each time a grid is visited.  Use `SNESFASSetCyclesOnLevel()` for more
complicated cycling.

Logically Collective

Input Parameters:
- `snes`   - the `SNESFAS` nonlinear multigrid context
- `cycles` - the number of cycles -- 1 for V-cycle, 2 for W-cycle

Options Database Key:
- `-snes_fas_cycles <1,2>` - 1 for V-cycle, 2 for W-cycle

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASSetCyclesOnLevel()`

# External Links
$(_doc_external("Snes/SNESFASSetCycles"))
"""
function SNESFASSetCycles(petsclib::PetscLibType, snes::PetscSNES, cycles::PetscInt) end

@for_petsc function SNESFASSetCycles(petsclib::$UnionPetscLib, snes::PetscSNES, cycles::$PetscInt )

    @chk ccall(
               (:SNESFASSetCycles, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt),
               snes, cycles,
              )


	return nothing
end 

"""
	SNESFASSetMonitor(petsclib::PetscLibType,snes::PetscSNES, vf::PetscViewerAndFormat, flg::PetscBool) 
Sets the method

Logically Collective

Input Parameters:
- `snes` - the `SNESFAS` context
- `vf`   - viewer and format structure (may be `NULL` if `flg` is `PETSC_FALSE`)
- `flg`  - monitor or not

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESMonitorSet()`, `SNESFASSetCyclesOnLevel()`

# External Links
$(_doc_external("Snes/SNESFASSetMonitor"))
"""
function SNESFASSetMonitor(petsclib::PetscLibType, snes::PetscSNES, vf::PetscViewerAndFormat, flg::PetscBool) end

@for_petsc function SNESFASSetMonitor(petsclib::$UnionPetscLib, snes::PetscSNES, vf::PetscViewerAndFormat, flg::PetscBool )

    @chk ccall(
               (:SNESFASSetMonitor, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{PetscViewerAndFormat}, PetscBool),
               snes, vf, flg,
              )


	return nothing
end 

"""
	SNESFASSetLog(petsclib::PetscLibType,snes::PetscSNES, flg::PetscBool) 
Sets or unsets time logging for various `SNESFAS` stages on all levels

Logically Collective

Input Parameters:
- `snes` - the `SNESFAS` context
- `flg`  - whether to log or not

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASSetMonitor()`

# External Links
$(_doc_external("Snes/SNESFASSetLog"))
"""
function SNESFASSetLog(petsclib::PetscLibType, snes::PetscSNES, flg::PetscBool) end

@for_petsc function SNESFASSetLog(petsclib::$UnionPetscLib, snes::PetscSNES, flg::PetscBool )

    @chk ccall(
               (:SNESFASSetLog, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, flg,
              )


	return nothing
end 

"""
	SNESFASCycleSetCycles(petsclib::PetscLibType,snes::PetscSNES, cycles::PetscInt) 
Sets the number of cycles for all levels in a `SNESFAS`

Logically Collective

Input Parameters:
- `snes`   - the `SNESFAS` nonlinear multigrid context
- `cycles` - the number of cycles -- 1 for V-cycle, 2 for W-cycle

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASSetCycles()`

# External Links
$(_doc_external("Snes/SNESFASCycleSetCycles"))
"""
function SNESFASCycleSetCycles(petsclib::PetscLibType, snes::PetscSNES, cycles::PetscInt) end

@for_petsc function SNESFASCycleSetCycles(petsclib::$UnionPetscLib, snes::PetscSNES, cycles::$PetscInt )

    @chk ccall(
               (:SNESFASCycleSetCycles, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt),
               snes, cycles,
              )


	return nothing
end 

"""
	SNESFASCycleGetSmoother(petsclib::PetscLibType,snes::PetscSNES, smooth::PetscSNES) 
Gets the smoother on a particular cycle level.

Logically Collective

Input Parameter:
- `snes` - the `SNESFAS` obtained with `SNESFASGetCycleSNES()`

Output Parameter:
- `smooth` - the smoother

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASCycleGetSmootherUp()`, `SNESFASCycleGetSmootherDown()`, `SNESFASGetCycleSNES()`

# External Links
$(_doc_external("Snes/SNESFASCycleGetSmoother"))
"""
function SNESFASCycleGetSmoother(petsclib::PetscLibType, snes::PetscSNES, smooth::PetscSNES) end

@for_petsc function SNESFASCycleGetSmoother(petsclib::$UnionPetscLib, snes::PetscSNES, smooth::PetscSNES )
	smooth_ = Ref(smooth.ptr)

    @chk ccall(
               (:SNESFASCycleGetSmoother, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CSNES}),
               snes, smooth_,
              )

	smooth.ptr = C_NULL

	return nothing
end 

"""
	SNESFASCycleGetSmootherUp(petsclib::PetscLibType,snes::PetscSNES, smoothu::PetscSNES) 
Gets the up smoother on a particular cycle level.

Logically Collective

Input Parameter:
- `snes` - the `SNESFAS` obtained with `SNESFASGetCycleSNES()`

Output Parameter:
- `smoothu` - the smoother

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASCycleGetSmoother()`, `SNESFASCycleGetSmootherDown()`, `SNESFASGetCycleSNES()`

# External Links
$(_doc_external("Snes/SNESFASCycleGetSmootherUp"))
"""
function SNESFASCycleGetSmootherUp(petsclib::PetscLibType, snes::PetscSNES, smoothu::PetscSNES) end

@for_petsc function SNESFASCycleGetSmootherUp(petsclib::$UnionPetscLib, snes::PetscSNES, smoothu::PetscSNES )
	smoothu_ = Ref(smoothu.ptr)

    @chk ccall(
               (:SNESFASCycleGetSmootherUp, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CSNES}),
               snes, smoothu_,
              )

	smoothu.ptr = C_NULL

	return nothing
end 

"""
	SNESFASCycleGetSmootherDown(petsclib::PetscLibType,snes::PetscSNES, smoothd::PetscSNES) 
Gets the down smoother on a particular cycle level.

Logically Collective

Input Parameter:
- `snes` - `SNESFAS` obtained with `SNESFASGetCycleSNES()`

Output Parameter:
- `smoothd` - the smoother

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASCycleGetSmootherUp()`, `SNESFASCycleGetSmoother()`, `SNESFASGetCycleSNES()`

# External Links
$(_doc_external("Snes/SNESFASCycleGetSmootherDown"))
"""
function SNESFASCycleGetSmootherDown(petsclib::PetscLibType, snes::PetscSNES, smoothd::PetscSNES) end

@for_petsc function SNESFASCycleGetSmootherDown(petsclib::$UnionPetscLib, snes::PetscSNES, smoothd::PetscSNES )
	smoothd_ = Ref(smoothd.ptr)

    @chk ccall(
               (:SNESFASCycleGetSmootherDown, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CSNES}),
               snes, smoothd_,
              )

	smoothd.ptr = C_NULL

	return nothing
end 

"""
	SNESFASCycleGetCorrection(petsclib::PetscLibType,snes::PetscSNES, correction::PetscSNES) 
Gets the coarse correction `SNESFAS` context for this level

Logically Collective

Input Parameter:
- `snes` - the `SNESFAS` obtained with `SNESFASGetCycleSNES()`

Output Parameter:
- `correction` - the coarse correction solve on this level

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS` `SNESFASCycleGetSmootherUp()`, `SNESFASCycleGetSmoother()`

# External Links
$(_doc_external("Snes/SNESFASCycleGetCorrection"))
"""
function SNESFASCycleGetCorrection(petsclib::PetscLibType, snes::PetscSNES, correction::PetscSNES) end

@for_petsc function SNESFASCycleGetCorrection(petsclib::$UnionPetscLib, snes::PetscSNES, correction::PetscSNES )
	correction_ = Ref(correction.ptr)

    @chk ccall(
               (:SNESFASCycleGetCorrection, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CSNES}),
               snes, correction_,
              )

	correction.ptr = C_NULL

	return nothing
end 

"""
	SNESFASCycleGetInterpolation(petsclib::PetscLibType,snes::PetscSNES, mat::PetscMat) 
Gets the interpolation on a level

Logically Collective

Input Parameter:
- `snes` - the `SNESFAS` obtained with `SNESFASGetCycleSNES()`

Output Parameter:
- `mat` - the interpolation operator on this level

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASCycleGetSmootherUp()`, `SNESFASCycleGetSmoother()`

# External Links
$(_doc_external("Snes/SNESFASCycleGetInterpolation"))
"""
function SNESFASCycleGetInterpolation(petsclib::PetscLibType, snes::PetscSNES, mat::PetscMat) end

@for_petsc function SNESFASCycleGetInterpolation(petsclib::$UnionPetscLib, snes::PetscSNES, mat::PetscMat )
	mat_ = Ref(mat.ptr)

    @chk ccall(
               (:SNESFASCycleGetInterpolation, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CMat}),
               snes, mat_,
              )

	mat.ptr = C_NULL

	return nothing
end 

"""
	SNESFASCycleGetRestriction(petsclib::PetscLibType,snes::PetscSNES, mat::PetscMat) 
Gets the restriction on a level

Logically Collective

Input Parameter:
- `snes` - the `SNESFAS` obtained with `SNESFASGetCycleSNES()`

Output Parameter:
- `mat` - the restriction operator on this level

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASGetRestriction()`, `SNESFASCycleGetInterpolation()`

# External Links
$(_doc_external("Snes/SNESFASCycleGetRestriction"))
"""
function SNESFASCycleGetRestriction(petsclib::PetscLibType, snes::PetscSNES, mat::PetscMat) end

@for_petsc function SNESFASCycleGetRestriction(petsclib::$UnionPetscLib, snes::PetscSNES, mat::PetscMat )
	mat_ = Ref(mat.ptr)

    @chk ccall(
               (:SNESFASCycleGetRestriction, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CMat}),
               snes, mat_,
              )

	mat.ptr = C_NULL

	return nothing
end 

"""
	SNESFASCycleGetInjection(petsclib::PetscLibType,snes::PetscSNES, mat::PetscMat) 
Gets the injection on a level

Logically Collective

Input Parameter:
- `snes` - the `SNESFAS` obtained with `SNESFASGetCycleSNES()`

Output Parameter:
- `mat` - the restriction operator on this level

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASGetInjection()`, `SNESFASCycleGetRestriction()`

# External Links
$(_doc_external("Snes/SNESFASCycleGetInjection"))
"""
function SNESFASCycleGetInjection(petsclib::PetscLibType, snes::PetscSNES, mat::PetscMat) end

@for_petsc function SNESFASCycleGetInjection(petsclib::$UnionPetscLib, snes::PetscSNES, mat::PetscMat )
	mat_ = Ref(mat.ptr)

    @chk ccall(
               (:SNESFASCycleGetInjection, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CMat}),
               snes, mat_,
              )

	mat.ptr = C_NULL

	return nothing
end 

"""
	SNESFASCycleGetRScale(petsclib::PetscLibType,snes::PetscSNES, vec::PetscVec) 
Gets the injection scale

Logically Collective

Input Parameter:
- `snes` - the  `SNESFAS` obtained with `SNESFASGetCycleSNES()`

Output Parameter:
- `vec` - the restriction operator on this level

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASCycleGetRestriction()`, `SNESFASGetRScale()`

# External Links
$(_doc_external("Snes/SNESFASCycleGetRScale"))
"""
function SNESFASCycleGetRScale(petsclib::PetscLibType, snes::PetscSNES, vec::PetscVec) end

@for_petsc function SNESFASCycleGetRScale(petsclib::$UnionPetscLib, snes::PetscSNES, vec::PetscVec )
	vec_ = Ref(vec.ptr)

    @chk ccall(
               (:SNESFASCycleGetRScale, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CVec}),
               snes, vec_,
              )

	vec.ptr = C_NULL

	return nothing
end 

"""
	flg::PetscBool = SNESFASCycleIsFine(petsclib::PetscLibType,snes::PetscSNES) 
Determines if a given `SNES` is the finest level in a `SNESFAS`

Logically Collective

Input Parameter:
- `snes` - the `SNESFAS` context obtained with `SNESFASGetCycleSNES()`

Output Parameter:
- `flg` - indicates if this is the fine level or not

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASSetLevels()`

# External Links
$(_doc_external("Snes/SNESFASCycleIsFine"))
"""
function SNESFASCycleIsFine(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESFASCycleIsFine(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESFASSetInterpolation(petsclib::PetscLibType,snes::PetscSNES, level::PetscInt, mat::PetscMat) 
Sets the `Mat` to be used to apply the
interpolation from l-1 to the lth level

Input Parameters:
- `snes`  - the `SNESFAS` nonlinear multigrid context
- `mat`   - the interpolation operator
- `level` - the level (0 is coarsest) to supply [do not supply 0]

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASSetInjection()`, `SNESFASSetRestriction()`, `SNESFASSetRScale()`

# External Links
$(_doc_external("Snes/SNESFASSetInterpolation"))
"""
function SNESFASSetInterpolation(petsclib::PetscLibType, snes::PetscSNES, level::PetscInt, mat::PetscMat) end

@for_petsc function SNESFASSetInterpolation(petsclib::$UnionPetscLib, snes::PetscSNES, level::$PetscInt, mat::PetscMat )

    @chk ccall(
               (:SNESFASSetInterpolation, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, CMat),
               snes, level, mat,
              )


	return nothing
end 

"""
	SNESFASGetInterpolation(petsclib::PetscLibType,snes::PetscSNES, level::PetscInt, mat::PetscMat) 
Gets the matrix used to calculate the
interpolation from l-1 to the lth level

Input Parameters:
- `snes`  - the `SNESFAS` nonlinear multigrid context
- `level` - the level (0 is coarsest) to supply [do not supply 0]

Output Parameter:
- `mat` - the interpolation operator

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASSetInterpolation()`, `SNESFASGetInjection()`, `SNESFASGetRestriction()`, `SNESFASGetRScale()`

# External Links
$(_doc_external("Snes/SNESFASGetInterpolation"))
"""
function SNESFASGetInterpolation(petsclib::PetscLibType, snes::PetscSNES, level::PetscInt, mat::PetscMat) end

@for_petsc function SNESFASGetInterpolation(petsclib::$UnionPetscLib, snes::PetscSNES, level::$PetscInt, mat::PetscMat )
	mat_ = Ref(mat.ptr)

    @chk ccall(
               (:SNESFASGetInterpolation, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, Ptr{CMat}),
               snes, level, mat_,
              )

	mat.ptr = C_NULL

	return nothing
end 

"""
	SNESFASSetRestriction(petsclib::PetscLibType,snes::PetscSNES, level::PetscInt, mat::PetscMat) 
Sets the matrix to be used to restrict the defect
from level l to l-1.

Input Parameters:
- `snes`  - the `SNESFAS` nonlinear multigrid context
- `mat`   - the restriction matrix
- `level` - the level (0 is coarsest) to supply [Do not supply 0]

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASSetInterpolation()`, `SNESFASSetInjection()`

# External Links
$(_doc_external("Snes/SNESFASSetRestriction"))
"""
function SNESFASSetRestriction(petsclib::PetscLibType, snes::PetscSNES, level::PetscInt, mat::PetscMat) end

@for_petsc function SNESFASSetRestriction(petsclib::$UnionPetscLib, snes::PetscSNES, level::$PetscInt, mat::PetscMat )

    @chk ccall(
               (:SNESFASSetRestriction, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, CMat),
               snes, level, mat,
              )


	return nothing
end 

"""
	SNESFASGetRestriction(petsclib::PetscLibType,snes::PetscSNES, level::PetscInt, mat::PetscMat) 
Gets the matrix used to calculate the
restriction from l to the l-1th level

Input Parameters:
- `snes`  - the `SNESFAS` nonlinear multigrid context
- `level` - the level (0 is coarsest) to supply [do not supply 0]

Output Parameter:
- `mat` - the interpolation operator

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASSetRestriction()`, `SNESFASGetInjection()`, `SNESFASGetInterpolation()`, `SNESFASGetRScale()`

# External Links
$(_doc_external("Snes/SNESFASGetRestriction"))
"""
function SNESFASGetRestriction(petsclib::PetscLibType, snes::PetscSNES, level::PetscInt, mat::PetscMat) end

@for_petsc function SNESFASGetRestriction(petsclib::$UnionPetscLib, snes::PetscSNES, level::$PetscInt, mat::PetscMat )
	mat_ = Ref(mat.ptr)

    @chk ccall(
               (:SNESFASGetRestriction, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, Ptr{CMat}),
               snes, level, mat_,
              )

	mat.ptr = C_NULL

	return nothing
end 

"""
	SNESFASSetInjection(petsclib::PetscLibType,snes::PetscSNES, level::PetscInt, mat::PetscMat) 
Sets the matrix to be used to inject the solution
from `level` to `level-1`.

Input Parameters:
- `snes`  - the `SNESFAS` nonlinear multigrid context
- `mat`   - the injection matrix
- `level` - the level (0 is coarsest) to supply [Do not supply 0]

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASSetInterpolation()`, `SNESFASSetRestriction()`

# External Links
$(_doc_external("Snes/SNESFASSetInjection"))
"""
function SNESFASSetInjection(petsclib::PetscLibType, snes::PetscSNES, level::PetscInt, mat::PetscMat) end

@for_petsc function SNESFASSetInjection(petsclib::$UnionPetscLib, snes::PetscSNES, level::$PetscInt, mat::PetscMat )

    @chk ccall(
               (:SNESFASSetInjection, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, CMat),
               snes, level, mat,
              )


	return nothing
end 

"""
	SNESFASGetInjection(petsclib::PetscLibType,snes::PetscSNES, level::PetscInt, mat::PetscMat) 
Gets the matrix used to calculate the
injection from l-1 to the lth level

Input Parameters:
- `snes`  - the `SNESFAS` nonlinear multigrid context
- `level` - the level (0 is coarsest) to supply [do not supply 0]

Output Parameter:
- `mat` - the injection operator

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASSetInjection()`, `SNESFASGetRestriction()`, `SNESFASGetInterpolation()`, `SNESFASGetRScale()`

# External Links
$(_doc_external("Snes/SNESFASGetInjection"))
"""
function SNESFASGetInjection(petsclib::PetscLibType, snes::PetscSNES, level::PetscInt, mat::PetscMat) end

@for_petsc function SNESFASGetInjection(petsclib::$UnionPetscLib, snes::PetscSNES, level::$PetscInt, mat::PetscMat )
	mat_ = Ref(mat.ptr)

    @chk ccall(
               (:SNESFASGetInjection, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, Ptr{CMat}),
               snes, level, mat_,
              )

	mat.ptr = C_NULL

	return nothing
end 

"""
	SNESFASSetRScale(petsclib::PetscLibType,snes::PetscSNES, level::PetscInt, rscale::PetscVec) 
Sets the scaling factor of the restriction
operator from level l to l-1.

Input Parameters:
- `snes`   - the `SNESFAS` nonlinear multigrid context
- `rscale` - the restriction scaling
- `level`  - the level (0 is coarsest) to supply [Do not supply 0]

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASSetInjection()`, `SNESFASSetRestriction()`

# External Links
$(_doc_external("Snes/SNESFASSetRScale"))
"""
function SNESFASSetRScale(petsclib::PetscLibType, snes::PetscSNES, level::PetscInt, rscale::PetscVec) end

@for_petsc function SNESFASSetRScale(petsclib::$UnionPetscLib, snes::PetscSNES, level::$PetscInt, rscale::PetscVec )

    @chk ccall(
               (:SNESFASSetRScale, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, CVec),
               snes, level, rscale,
              )


	return nothing
end 

"""
	SNESFASGetSmoother(petsclib::PetscLibType,snes::PetscSNES, level::PetscInt, smooth::PetscSNES) 
Gets the default smoother on a level.

Input Parameters:
- `snes`  - the `SNESFAS` nonlinear multigrid context
- `level` - the level (0 is coarsest) to supply

Output Parameter:
- `smooth` - the smoother

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASSetInjection()`, `SNESFASSetRestriction()`

# External Links
$(_doc_external("Snes/SNESFASGetSmoother"))
"""
function SNESFASGetSmoother(petsclib::PetscLibType, snes::PetscSNES, level::PetscInt, smooth::PetscSNES) end

@for_petsc function SNESFASGetSmoother(petsclib::$UnionPetscLib, snes::PetscSNES, level::$PetscInt, smooth::PetscSNES )
	smooth_ = Ref(smooth.ptr)

    @chk ccall(
               (:SNESFASGetSmoother, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, Ptr{CSNES}),
               snes, level, smooth_,
              )

	smooth.ptr = C_NULL

	return nothing
end 

"""
	SNESFASGetSmootherDown(petsclib::PetscLibType,snes::PetscSNES, level::PetscInt, smooth::PetscSNES) 
Gets the downsmoother on a level.

Input Parameters:
- `snes`  - the `SNESFAS` nonlinear multigrid context
- `level` - the level (0 is coarsest) to supply

Output Parameter:
- `smooth` - the smoother

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASSetInjection()`, `SNESFASSetRestriction()`

# External Links
$(_doc_external("Snes/SNESFASGetSmootherDown"))
"""
function SNESFASGetSmootherDown(petsclib::PetscLibType, snes::PetscSNES, level::PetscInt, smooth::PetscSNES) end

@for_petsc function SNESFASGetSmootherDown(petsclib::$UnionPetscLib, snes::PetscSNES, level::$PetscInt, smooth::PetscSNES )
	smooth_ = Ref(smooth.ptr)

    @chk ccall(
               (:SNESFASGetSmootherDown, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, Ptr{CSNES}),
               snes, level, smooth_,
              )

	smooth.ptr = C_NULL

	return nothing
end 

"""
	SNESFASGetSmootherUp(petsclib::PetscLibType,snes::PetscSNES, level::PetscInt, smooth::PetscSNES) 
Gets the upsmoother on a level.

Input Parameters:
- `snes`  - the `SNESFAS` nonlinear multigrid context
- `level` - the level (0 is coarsest)

Output Parameter:
- `smooth` - the smoother

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASSetInjection()`, `SNESFASSetRestriction()`

# External Links
$(_doc_external("Snes/SNESFASGetSmootherUp"))
"""
function SNESFASGetSmootherUp(petsclib::PetscLibType, snes::PetscSNES, level::PetscInt, smooth::PetscSNES) end

@for_petsc function SNESFASGetSmootherUp(petsclib::$UnionPetscLib, snes::PetscSNES, level::$PetscInt, smooth::PetscSNES )
	smooth_ = Ref(smooth.ptr)

    @chk ccall(
               (:SNESFASGetSmootherUp, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, Ptr{CSNES}),
               snes, level, smooth_,
              )

	smooth.ptr = C_NULL

	return nothing
end 

"""
	SNESFASGetCoarseSolve(petsclib::PetscLibType,snes::PetscSNES, coarse::PetscSNES) 
Gets the coarsest level solver.

Input Parameter:
- `snes` - the `SNESFAS` nonlinear multigrid context

Output Parameter:
- `coarse` - the coarse-level solver

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASSetInjection()`, `SNESFASSetRestriction()`

# External Links
$(_doc_external("Snes/SNESFASGetCoarseSolve"))
"""
function SNESFASGetCoarseSolve(petsclib::PetscLibType, snes::PetscSNES, coarse::PetscSNES) end

@for_petsc function SNESFASGetCoarseSolve(petsclib::$UnionPetscLib, snes::PetscSNES, coarse::PetscSNES )
	coarse_ = Ref(coarse.ptr)

    @chk ccall(
               (:SNESFASGetCoarseSolve, $petsc_library),
               PetscErrorCode,
               (CSNES, Ptr{CSNES}),
               snes, coarse_,
              )

	coarse.ptr = C_NULL

	return nothing
end 

"""
	SNESFASFullSetDownSweep(petsclib::PetscLibType,snes::PetscSNES, swp::PetscBool) 
Smooth during the initial downsweep for `SNESFAS`

Logically Collective

Input Parameters:
- `snes` - the `SNESFAS` nonlinear multigrid context
- `swp`  - whether to downsweep or not

Options Database Key:
- `-snes_fas_full_downsweep` - Sets whether to smooth on the initial downsweep

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASSetNumberSmoothUp()`

# External Links
$(_doc_external("Snes/SNESFASFullSetDownSweep"))
"""
function SNESFASFullSetDownSweep(petsclib::PetscLibType, snes::PetscSNES, swp::PetscBool) end

@for_petsc function SNESFASFullSetDownSweep(petsclib::$UnionPetscLib, snes::PetscSNES, swp::PetscBool )

    @chk ccall(
               (:SNESFASFullSetDownSweep, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, swp,
              )


	return nothing
end 

"""
	SNESFASFullSetTotal(petsclib::PetscLibType,snes::PetscSNES, total::PetscBool) 
Use total residual restriction and total interpolation on the initial down and up sweep of full `SNESFAS` cycles

Logically Collective

Input Parameters:
- `snes`  - the `SNESFAS`  nonlinear multigrid context
- `total` - whether to use total restriction / interpolatiaon or not (the alternative is defect restriction and correction interpolation)

Options Database Key:
- `-snes_fas_full_total` - Use total restriction and interpolation on the initial down and up sweeps for the full `SNESFAS` cycle

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASSetNumberSmoothUp()`, `DMInterpolateSolution()`

# External Links
$(_doc_external("Snes/SNESFASFullSetTotal"))
"""
function SNESFASFullSetTotal(petsclib::PetscLibType, snes::PetscSNES, total::PetscBool) end

@for_petsc function SNESFASFullSetTotal(petsclib::$UnionPetscLib, snes::PetscSNES, total::PetscBool )

    @chk ccall(
               (:SNESFASFullSetTotal, $petsc_library),
               PetscErrorCode,
               (CSNES, PetscBool),
               snes, total,
              )


	return nothing
end 

"""
	total::PetscBool = SNESFASFullGetTotal(petsclib::PetscLibType,snes::PetscSNES) 
Use total residual restriction and total interpolation on the initial down and up sweep of full FAS cycles

Logically Collective

Input Parameter:
- `snes` - the `SNESFAS` nonlinear multigrid context

Output Parameter:
- `total` - whether to use total restriction / interpolatiaon or not (the alternative is defect restriction and correction interpolation)

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESFAS`, `SNESFASSetNumberSmoothUp()`, `DMInterpolateSolution()`, `SNESFullSetTotal()`

# External Links
$(_doc_external("Snes/SNESFASFullGetTotal"))
"""
function SNESFASFullGetTotal(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESFASFullGetTotal(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESNewtonTRSetNormType(petsclib::PetscLibType,snes::PetscSNES, norm::NormType) 
Specify the type of norm to use for the computation of the trust region.

Input Parameters:
- `snes` - the nonlinear solver object
- `norm` - the norm type

Level: intermediate

-seealso: `SNESNEWTONTR`, `NormType`

# External Links
$(_doc_external("Snes/SNESNewtonTRSetNormType"))
"""
function SNESNewtonTRSetNormType(petsclib::PetscLibType, snes::PetscSNES, norm::NormType) end

@for_petsc function SNESNewtonTRSetNormType(petsclib::$UnionPetscLib, snes::PetscSNES, norm::NormType )

    @chk ccall(
               (:SNESNewtonTRSetNormType, $petsc_library),
               PetscErrorCode,
               (CSNES, NormType),
               snes, norm,
              )


	return nothing
end 

"""
	SNESNewtonTRSetQNType(petsclib::PetscLibType,snes::PetscSNES, use::SNESNewtonTRQNType) 
Specify to use a quasi

Input Parameters:
- `snes` - the nonlinear solver object
- `use`  - the type of approximations to be used

Level: intermediate

-seealso: `SNESNEWTONTR`, `SNESNewtonTRQNType`, `MATLMVM`

# External Links
$(_doc_external("Snes/SNESNewtonTRSetQNType"))
"""
function SNESNewtonTRSetQNType(petsclib::PetscLibType, snes::PetscSNES, use::SNESNewtonTRQNType) end

@for_petsc function SNESNewtonTRSetQNType(petsclib::$UnionPetscLib, snes::PetscSNES, use::SNESNewtonTRQNType )

    @chk ccall(
               (:SNESNewtonTRSetQNType, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESNewtonTRQNType),
               snes, use,
              )


	return nothing
end 

"""
	SNESNewtonTRSetFallbackType(petsclib::PetscLibType,snes::PetscSNES, ftype::SNESNewtonTRFallbackType) 
Set the type of fallback to use if the solution of the trust region subproblem is outside the radius

Input Parameters:
- `snes`  - the nonlinear solver object
- `ftype` - the fallback type, see `SNESNewtonTRFallbackType`

Level: intermediate

-seealso: [](ch_snes), `SNESNEWTONTR`, `SNESNewtonTRPreCheck()`, `SNESNewtonTRGetPreCheck()`, `SNESNewtonTRSetPreCheck()`,
`SNESNewtonTRSetPostCheck()`, `SNESNewtonTRGetPostCheck()`

# External Links
$(_doc_external("Snes/SNESNewtonTRSetFallbackType"))
"""
function SNESNewtonTRSetFallbackType(petsclib::PetscLibType, snes::PetscSNES, ftype::SNESNewtonTRFallbackType) end

@for_petsc function SNESNewtonTRSetFallbackType(petsclib::$UnionPetscLib, snes::PetscSNES, ftype::SNESNewtonTRFallbackType )

    @chk ccall(
               (:SNESNewtonTRSetFallbackType, $petsc_library),
               PetscErrorCode,
               (CSNES, SNESNewtonTRFallbackType),
               snes, ftype,
              )


	return nothing
end 

"""
	SNESNewtonTRSetPreCheck(petsclib::PetscLibType,snes::PetscSNES, func::external, ctx::Cvoid) 
Sets a user function that is called before the search step has been determined.
Allows the user a chance to change or override the trust region decision.

Logically Collective

Input Parameters:
- `snes` - the nonlinear solver object
- `func` - [optional] function evaluation routine, for the calling sequence see `SNESNewtonTRPreCheck()`
- `ctx`  - [optional] user-defined context for private data for the function evaluation routine (may be `NULL`)

Level: intermediate

-seealso: [](ch_snes), `SNESNEWTONTR`, `SNESNewtonTRPreCheck()`, `SNESNewtonTRGetPreCheck()`, `SNESNewtonTRSetPostCheck()`, `SNESNewtonTRGetPostCheck()`,

# External Links
$(_doc_external("Snes/SNESNewtonTRSetPreCheck"))
"""
function SNESNewtonTRSetPreCheck(petsclib::PetscLibType, snes::PetscSNES, func::external, ctx::Cvoid) end

@for_petsc function SNESNewtonTRSetPreCheck(petsclib::$UnionPetscLib, snes::PetscSNES, func::external, ctx::Cvoid )

    @chk ccall(
               (:SNESNewtonTRSetPreCheck, $petsc_library),
               PetscErrorCode,
               (CSNES, external, Ptr{Cvoid}),
               snes, func, ctx,
              )


	return nothing
end 

"""
	SNESNewtonTRSetPostCheck(petsclib::PetscLibType,snes::PetscSNES, func::external, ctx::Cvoid) 
Sets a user function that is called after the search step has been determined but before the next
function evaluation. Allows the user a chance to change or override the internal decision of the solver

Logically Collective

Input Parameters:
- `snes` - the nonlinear solver object
- `func` - [optional] function evaluation routine, for the calling sequence see `SNESNewtonTRPostCheck()`
- `ctx`  - [optional] user-defined context for private data for the function evaluation routine (may be `NULL`)

Level: intermediate

-seealso: [](ch_snes), `SNESNEWTONTR`, `SNESNewtonTRPostCheck()`, `SNESNewtonTRGetPostCheck()`, `SNESNewtonTRSetPreCheck()`, `SNESNewtonTRGetPreCheck()`

# External Links
$(_doc_external("Snes/SNESNewtonTRSetPostCheck"))
"""
function SNESNewtonTRSetPostCheck(petsclib::PetscLibType, snes::PetscSNES, func::external, ctx::Cvoid) end

@for_petsc function SNESNewtonTRSetPostCheck(petsclib::$UnionPetscLib, snes::PetscSNES, func::external, ctx::Cvoid )

    @chk ccall(
               (:SNESNewtonTRSetPostCheck, $petsc_library),
               PetscErrorCode,
               (CSNES, external, Ptr{Cvoid}),
               snes, func, ctx,
              )


	return nothing
end 

"""
	changed_Y::PetscBool = SNESNewtonTRPreCheck(petsclib::PetscLibType,snes::PetscSNES, X::PetscVec, Y::PetscVec) 
Runs the precheck routine

Logically Collective

Input Parameters:
- `snes` - the solver
- `X`    - The last solution
- `Y`    - The step direction

Output Parameter:
- `changed_Y` - Indicator that the step direction `Y` has been changed.

Level: intermediate

-seealso: [](ch_snes), `SNESNEWTONTR`, `SNESNewtonTRSetPreCheck()`, `SNESNewtonTRGetPreCheck()`, `SNESNewtonTRPostCheck()`

# External Links
$(_doc_external("Snes/SNESNewtonTRPreCheck"))
"""
function SNESNewtonTRPreCheck(petsclib::PetscLibType, snes::PetscSNES, X::PetscVec, Y::PetscVec) end

@for_petsc function SNESNewtonTRPreCheck(petsclib::$UnionPetscLib, snes::PetscSNES, X::PetscVec, Y::PetscVec )
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
	changed_Y::PetscBool,changed_W::PetscBool = SNESNewtonTRPostCheck(petsclib::PetscLibType,snes::PetscSNES, X::PetscVec, Y::PetscVec, W::PetscVec) 
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

-seealso: [](ch_snes), `SNESNEWTONTR`, `SNESNewtonTRSetPostCheck()`, `SNESNewtonTRGetPostCheck()`, `SNESNewtonTRPreCheck()`

# External Links
$(_doc_external("Snes/SNESNewtonTRPostCheck"))
"""
function SNESNewtonTRPostCheck(petsclib::PetscLibType, snes::PetscSNES, X::PetscVec, Y::PetscVec, W::PetscVec) end

@for_petsc function SNESNewtonTRPostCheck(petsclib::$UnionPetscLib, snes::PetscSNES, X::PetscVec, Y::PetscVec, W::PetscVec )
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
	SNESNewtonTRSetTolerances(petsclib::PetscLibType,snes::PetscSNES, delta_min::PetscReal, delta_max::PetscReal, delta_0::PetscReal) 
Sets the trust region parameter tolerances.

Logically Collective

Input Parameters:
- `snes`      - the `SNES` context
- `delta_min` - minimum allowed trust region size
- `delta_max` - maximum allowed trust region size
- `delta_0`   - initial trust region size

Options Database Key:
- `-snes_tr_deltamin <tol>` - Set minimum size
- `-snes_tr_deltamax <tol>` - Set maximum size
- `-snes_tr_delta0   <tol>` - Set initial size

-seealso: [](ch_snes), `SNES`, `SNESNEWTONTR`, `SNESNewtonTRGetTolerances()`

# External Links
$(_doc_external("Snes/SNESNewtonTRSetTolerances"))
"""
function SNESNewtonTRSetTolerances(petsclib::PetscLibType, snes::PetscSNES, delta_min::PetscReal, delta_max::PetscReal, delta_0::PetscReal) end

@for_petsc function SNESNewtonTRSetTolerances(petsclib::$UnionPetscLib, snes::PetscSNES, delta_min::$PetscReal, delta_max::$PetscReal, delta_0::$PetscReal )

    @chk ccall(
               (:SNESNewtonTRSetTolerances, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscReal, $PetscReal, $PetscReal),
               snes, delta_min, delta_max, delta_0,
              )


	return nothing
end 

"""
	delta_min::PetscReal,delta_max::PetscReal,delta_0::PetscReal = SNESNewtonTRGetTolerances(petsclib::PetscLibType,snes::PetscSNES) 
Gets the trust region parameter tolerances.

Not Collective

Input Parameter:
- `snes` - the `SNES` context

Output Parameters:
- `delta_min` - minimum allowed trust region size or `NULL`
- `delta_max` - maximum allowed trust region size or `NULL`
- `delta_0`   - initial trust region size or `NULL`

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESNEWTONTR`, `SNESNewtonTRSetTolerances()`

# External Links
$(_doc_external("Snes/SNESNewtonTRGetTolerances"))
"""
function SNESNewtonTRGetTolerances(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESNewtonTRGetTolerances(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESNewtonTRSetUpdateParameters(petsclib::PetscLibType,snes::PetscSNES, eta1::PetscReal, eta2::PetscReal, eta3::PetscReal, t1::PetscReal, t2::PetscReal) 
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
- `-snes_tr_eta1 <tol>` - Set `eta1`
- `-snes_tr_eta2 <tol>` - Set `eta2`
- `-snes_tr_eta3 <tol>` - Set `eta3`
- `-snes_tr_t1   <tol>` - Set `t1`
- `-snes_tr_t2   <tol>` - Set `t2`

-seealso: [](ch_snes), `SNES`, `SNESNEWTONTR`, `SNESSetObjective()`, `SNESNewtonTRGetUpdateParameters()`

# External Links
$(_doc_external("Snes/SNESNewtonTRSetUpdateParameters"))
"""
function SNESNewtonTRSetUpdateParameters(petsclib::PetscLibType, snes::PetscSNES, eta1::PetscReal, eta2::PetscReal, eta3::PetscReal, t1::PetscReal, t2::PetscReal) end

@for_petsc function SNESNewtonTRSetUpdateParameters(petsclib::$UnionPetscLib, snes::PetscSNES, eta1::$PetscReal, eta2::$PetscReal, eta3::$PetscReal, t1::$PetscReal, t2::$PetscReal )

    @chk ccall(
               (:SNESNewtonTRSetUpdateParameters, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscReal, $PetscReal, $PetscReal, $PetscReal, $PetscReal),
               snes, eta1, eta2, eta3, t1, t2,
              )


	return nothing
end 

"""
	eta1::PetscReal,eta2::PetscReal,eta3::PetscReal,t1::PetscReal,t2::PetscReal = SNESNewtonTRGetUpdateParameters(petsclib::PetscLibType,snes::PetscSNES) 
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

-seealso: [](ch_snes), `SNES`, `SNESNEWTONTR`, `SNESNewtonTRSetUpdateParameters()`

# External Links
$(_doc_external("Snes/SNESNewtonTRGetUpdateParameters"))
"""
function SNESNewtonTRGetUpdateParameters(petsclib::PetscLibType, snes::PetscSNES) end

@for_petsc function SNESNewtonTRGetUpdateParameters(petsclib::$UnionPetscLib, snes::PetscSNES )
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
	SNESConvergedCorrectPressure(petsclib::PetscLibType,snes::PetscSNES, it::PetscInt, xnorm::PetscReal, gnorm::PetscReal, f::PetscReal, reason::SNESConvergedReason, ctx::Cvoid) 
The regular `SNES` convergence test that, up on convergence, adds a vector in the nullspace
to make the continuum integral of the pressure field equal to zero.

Logically Collective

Input Parameters:
- `snes`  - the `SNES` context
- `it`    - the iteration (0 indicates before any Newton steps)
- `xnorm` - 2-norm of current iterate
- `gnorm` - 2-norm of current step
- `f`     - 2-norm of function at current iterate
- `ctx`   - Optional user context

Output Parameter:
- `reason` - `SNES_CONVERGED_ITERATING`, `SNES_CONVERGED_ITS`, or `SNES_DIVERGED_FNORM_NAN`

Options Database Key:
- `-snes_convergence_test correct_pressure` - see `SNESSetFromOptions()`

Level: advanced

-seealso: [](ch_snes), `SNES`, `DM`, `SNESConvergedDefault()`, `SNESSetConvergenceTest()`, `DMSetNullSpaceConstructor()`

# External Links
$(_doc_external("Snes/SNESConvergedCorrectPressure"))
"""
function SNESConvergedCorrectPressure(petsclib::PetscLibType, snes::PetscSNES, it::PetscInt, xnorm::PetscReal, gnorm::PetscReal, f::PetscReal, reason::SNESConvergedReason, ctx::Cvoid) end

@for_petsc function SNESConvergedCorrectPressure(petsclib::$UnionPetscLib, snes::PetscSNES, it::$PetscInt, xnorm::$PetscReal, gnorm::$PetscReal, f::$PetscReal, reason::SNESConvergedReason, ctx::Cvoid )

    @chk ccall(
               (:SNESConvergedCorrectPressure, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, $PetscReal, $PetscReal, Ptr{SNESConvergedReason}, Ptr{Cvoid}),
               snes, it, xnorm, gnorm, f, reason, ctx,
              )


	return nothing
end 

"""
	SNESMonitorFields(petsclib::PetscLibType,snes::PetscSNES, its::PetscInt, fgnorm::PetscReal, vf::PetscViewerAndFormat) 
Monitors the residual for each field separately

Collective

Input Parameters:
- `snes`   - the `SNES` context, must have an attached `DM`
- `its`    - iteration number
- `fgnorm` - 2-norm of residual
- `vf`     - `PetscViewerAndFormat` of `PetscViewerType` `PETSCVIEWERASCII`

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESMonitorSet()`, `SNESMonitorDefault()`

# External Links
$(_doc_external("Snes/SNESMonitorFields"))
"""
function SNESMonitorFields(petsclib::PetscLibType, snes::PetscSNES, its::PetscInt, fgnorm::PetscReal, vf::PetscViewerAndFormat) end

@for_petsc function SNESMonitorFields(petsclib::$UnionPetscLib, snes::PetscSNES, its::$PetscInt, fgnorm::$PetscReal, vf::PetscViewerAndFormat )

    @chk ccall(
               (:SNESMonitorFields, $petsc_library),
               PetscErrorCode,
               (CSNES, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               snes, its, fgnorm, vf,
              )


	return nothing
end 

"""
	SNESTSFormFunction(petsclib::PetscLibType,snes::PetscSNES, U::PetscVec, F::PetscVec, ctx::Cvoid) 
Function to evaluate nonlinear residual defined by an ODE solver algorithm implemented within `TS`

Logically Collective

Input Parameters:
- `snes` - nonlinear solver
- `U`    - the current state at which to evaluate the residual
- `ctx`  - user context, must be a `TS`

Output Parameter:
- `F` - the nonlinear residual

Level: developer

-seealso: [](ch_ts), `SNESSetFunction()`, `MatFDColoringSetFunction()`

# External Links
$(_doc_external("Ts/SNESTSFormFunction"))
"""
function SNESTSFormFunction(petsclib::PetscLibType, snes::PetscSNES, U::PetscVec, F::PetscVec, ctx::Cvoid) end

@for_petsc function SNESTSFormFunction(petsclib::$UnionPetscLib, snes::PetscSNES, U::PetscVec, F::PetscVec, ctx::Cvoid )

    @chk ccall(
               (:SNESTSFormFunction, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CVec, Ptr{Cvoid}),
               snes, U, F, ctx,
              )


	return nothing
end 

"""
	SNESTSFormJacobian(petsclib::PetscLibType,snes::PetscSNES, U::PetscVec, A::PetscMat, B::PetscMat, ctx::Cvoid) 
Function to evaluate the Jacobian defined by an ODE solver algorithm implemented within `TS`

Collective

Input Parameters:
- `snes` - nonlinear solver
- `U`    - the current state at which to evaluate the residual
- `ctx`  - user context, must be a `TS`

Output Parameters:
- `A` - the Jacobian
- `B` - the matrix used to construct the preconditioner (often the same as `A`)

Level: developer

-seealso: [](ch_ts), `SNESSetJacobian()`

# External Links
$(_doc_external("Ts/SNESTSFormJacobian"))
"""
function SNESTSFormJacobian(petsclib::PetscLibType, snes::PetscSNES, U::PetscVec, A::PetscMat, B::PetscMat, ctx::Cvoid) end

@for_petsc function SNESTSFormJacobian(petsclib::$UnionPetscLib, snes::PetscSNES, U::PetscVec, A::PetscMat, B::PetscMat, ctx::Cvoid )

    @chk ccall(
               (:SNESTSFormJacobian, $petsc_library),
               PetscErrorCode,
               (CSNES, CVec, CMat, CMat, Ptr{Cvoid}),
               snes, U, A, B, ctx,
              )


	return nothing
end 

