# autodefined type arguments for class ------
mutable struct SNESLineSearchVIProjectFn end

mutable struct SNESLineSearchVINormFn end

mutable struct SNESLineSearchVIDirDerivFn end

mutable struct SNESLineSearchShellApplyFn end

mutable struct _n_SNESLineSearch end
const SNESLineSearch = Ptr{_n_SNESLineSearch}

# -------------------------------------------------------
"""
	SNESLineSearchMonitorCancel(petsclib::PetscLibType,ls::SNESLineSearch) 
Clears all the monitor functions for a `SNESLineSearch` object.

Logically Collective

Input Parameter:
- `ls` - the `SNESLineSearch` context

Options Database Key:
- `-snes_linesearch_monitor_cancel` - cancels all monitors that have been hardwired
into a code by calls to `SNESLineSearchMonitorSet()`, but does not cancel those
set via the options database

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESGetLineSearch()`, `SNESLineSearchMonitorDefault()`, `SNESLineSearchMonitorSet()`

# External Links
$(_doc_external("Snes/SNESLineSearchMonitorCancel"))
"""
function SNESLineSearchMonitorCancel(petsclib::PetscLibType, ls::SNESLineSearch) end

@for_petsc function SNESLineSearchMonitorCancel(petsclib::$UnionPetscLib, ls::SNESLineSearch )

    @chk ccall(
               (:SNESLineSearchMonitorCancel, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch,),
               ls,
              )


	return nothing
end 

"""
	SNESLineSearchMonitor(petsclib::PetscLibType,ls::SNESLineSearch) 
runs the user provided monitor routines, if they exist

Collective

Input Parameter:
- `ls` - the linesearch object

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESGetLineSearch()`, `SNESLineSearchMonitorSet()`

# External Links
$(_doc_external("Snes/SNESLineSearchMonitor"))
"""
function SNESLineSearchMonitor(petsclib::PetscLibType, ls::SNESLineSearch) end

@for_petsc function SNESLineSearchMonitor(petsclib::$UnionPetscLib, ls::SNESLineSearch )

    @chk ccall(
               (:SNESLineSearchMonitor, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch,),
               ls,
              )


	return nothing
end 

"""
	SNESLineSearchMonitorSet(petsclib::PetscLibType,ls::SNESLineSearch, f::external, mctx::Cvoid, monitordestroy::PetscCtxDestroyFn) 
Sets an ADDITIONAL function that is to be used at every
iteration of the nonlinear solver to display the iteration's
progress.

Logically Collective

Input Parameters:
- `ls`             - the `SNESLineSearch` context
- `f`              - the monitor function
- `mctx`           - [optional] user-defined context for private data for the monitor routine (use `NULL` if no context is desired)
- `monitordestroy` - [optional] routine that frees monitor context (may be `NULL`), see `PetscCtxDestroyFn` for the calling sequence

Calling sequence of `f`:
- `ls`   - the `SNESLineSearch` context
- `mctx` - [optional] user-defined context for private data for the monitor routine

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESGetLineSearch()`, `SNESLineSearchMonitorDefault()`, `SNESLineSearchMonitorCancel()`, `PetscCtxDestroyFn`

# External Links
$(_doc_external("Snes/SNESLineSearchMonitorSet"))
"""
function SNESLineSearchMonitorSet(petsclib::PetscLibType, ls::SNESLineSearch, f::external, mctx::Cvoid, monitordestroy::PetscCtxDestroyFn) end

@for_petsc function SNESLineSearchMonitorSet(petsclib::$UnionPetscLib, ls::SNESLineSearch, f::external, mctx::Cvoid, monitordestroy::PetscCtxDestroyFn )

    @chk ccall(
               (:SNESLineSearchMonitorSet, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, external, Ptr{Cvoid}, Ptr{PetscCtxDestroyFn}),
               ls, f, mctx, monitordestroy,
              )


	return nothing
end 

"""
	SNESLineSearchMonitorSolutionUpdate(petsclib::PetscLibType,ls::SNESLineSearch, vf::PetscViewerAndFormat) 
Monitors each update of the function value the linesearch tries

Collective

Input Parameters:
- `ls` - the `SNESLineSearch` object
- `vf` - the context for the monitor, in this case it is an `PetscViewerAndFormat`

Options Database Key:
- `-snes_linesearch_monitor_solution_update [viewer:filename:format]` - view each update tried by line search routine

Level: developer

This is not normally called directly but is passed to `SNESLineSearchMonitorSet()`

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESGetLineSearch()`, `SNESLineSearchMonitorSet()`, `SNESMonitorSolution()`

# External Links
$(_doc_external("Snes/SNESLineSearchMonitorSolutionUpdate"))
"""
function SNESLineSearchMonitorSolutionUpdate(petsclib::PetscLibType, ls::SNESLineSearch, vf::PetscViewerAndFormat) end

@for_petsc function SNESLineSearchMonitorSolutionUpdate(petsclib::$UnionPetscLib, ls::SNESLineSearch, vf::PetscViewerAndFormat )

    @chk ccall(
               (:SNESLineSearchMonitorSolutionUpdate, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, Ptr{PetscViewerAndFormat}),
               ls, vf,
              )


	return nothing
end 

"""
	outlinesearch::SNESLineSearch = SNESLineSearchCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Creates a `SNESLineSearch` context.

Logically Collective

Input Parameter:
- `comm` - MPI communicator for the line search (typically from the associated `SNES` context).

Output Parameter:
- `outlinesearch` - the new line search context

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `LineSearchDestroy()`, `SNESGetLineSearch()`

# External Links
$(_doc_external("Snes/SNESLineSearchCreate"))
"""
function SNESLineSearchCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function SNESLineSearchCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	outlinesearch_ = Ref{SNESLineSearch}()

    @chk ccall(
               (:SNESLineSearchCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{SNESLineSearch}),
               comm, outlinesearch_,
              )

	outlinesearch = outlinesearch_[]

	return outlinesearch
end 

"""
	SNESLineSearchSetUp(petsclib::PetscLibType,linesearch::SNESLineSearch) 
Prepares the line search for being applied by allocating
any required vectors.

Collective

Input Parameter:
- `linesearch` - The `SNESLineSearch` instance.

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESGetLineSearch()`, `SNESLineSearchReset()`

# External Links
$(_doc_external("Snes/SNESLineSearchSetUp"))
"""
function SNESLineSearchSetUp(petsclib::PetscLibType, linesearch::SNESLineSearch) end

@for_petsc function SNESLineSearchSetUp(petsclib::$UnionPetscLib, linesearch::SNESLineSearch )

    @chk ccall(
               (:SNESLineSearchSetUp, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch,),
               linesearch,
              )


	return nothing
end 

"""
	SNESLineSearchReset(petsclib::PetscLibType,linesearch::SNESLineSearch) 
Undoes the `SNESLineSearchSetUp()` and deletes any `Vec`s or `Mat`s allocated by the line search.

Collective

Input Parameter:
- `linesearch` - The `SNESLineSearch` instance.

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESGetLineSearch()`, `SNESLineSearchSetUp()`

# External Links
$(_doc_external("Snes/SNESLineSearchReset"))
"""
function SNESLineSearchReset(petsclib::PetscLibType, linesearch::SNESLineSearch) end

@for_petsc function SNESLineSearchReset(petsclib::$UnionPetscLib, linesearch::SNESLineSearch )

    @chk ccall(
               (:SNESLineSearchReset, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch,),
               linesearch,
              )


	return nothing
end 

"""
	SNESLineSearchSetFunction(petsclib::PetscLibType,linesearch::SNESLineSearch, func::external) 
Sets the function evaluation used by the `SNES` line search
`

Input Parameters:
- `linesearch` - the `SNESLineSearch` context
- `func`       - function evaluation routine, this is usually the function provided with `SNESSetFunction()`

Calling sequence of `func`:
- `snes` - the `SNES` with which the `SNESLineSearch` context is associated with
- `x`    - the input vector
- `f`    - the computed value of the function

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESGetLineSearch()`, `SNESSetFunction()`

# External Links
$(_doc_external("Snes/SNESLineSearchSetFunction"))
"""
function SNESLineSearchSetFunction(petsclib::PetscLibType, linesearch::SNESLineSearch, func::external) end

@for_petsc function SNESLineSearchSetFunction(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, func::external )

    @chk ccall(
               (:SNESLineSearchSetFunction, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, external),
               linesearch, func,
              )


	return nothing
end 

"""
	SNESLineSearchSetPreCheck(petsclib::PetscLibType,linesearch::SNESLineSearch, func::external, ctx::Cvoid) 
Sets a function that is called after the initial search direction has been computed but
before the line search routine has been applied. Allows adjusting the result of (usually a linear solve) that
determined the search direction.

Logically Collective

Input Parameters:
- `linesearch` - the `SNESLineSearch` context
- `func`       - [optional] function evaluation routine
- `ctx`        - [optional] user-defined context for private data for the function evaluation routine (may be `NULL`)

Calling sequence of `func`:
- `ls`        - the `SNESLineSearch` context
- `x`         - the current solution
- `d`         - the current search direction
- `changed_d` - indicates if the search direction has been changed
- `ctx`       - the context passed to `SNESLineSearchSetPreCheck()`

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESGetLineSearch()`, `SNESLineSearchPreCheck()`, `SNESLineSearchSetPostCheck()`, `SNESLineSearchGetPostCheck()`, `SNESLineSearchGetPreCheck()`,
`SNESVISetVariableBounds()`, `SNESVISetComputeVariableBounds()`, `SNESSetFunctionDomainError()`, `SNESSetJacobianDomainError()`


# External Links
$(_doc_external("Snes/SNESLineSearchSetPreCheck"))
"""
function SNESLineSearchSetPreCheck(petsclib::PetscLibType, linesearch::SNESLineSearch, func::external, ctx::Cvoid) end

@for_petsc function SNESLineSearchSetPreCheck(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, func::external, ctx::Cvoid )

    @chk ccall(
               (:SNESLineSearchSetPreCheck, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, external, Ptr{Cvoid}),
               linesearch, func, ctx,
              )


	return nothing
end 

"""
	SNESLineSearchSetPostCheck(petsclib::PetscLibType,linesearch::SNESLineSearch, func::external, ctx::Cvoid) 
Sets a user function that is called after the line search has been applied to determine the step
direction and length. Allows the user a chance to change or override the decision of the line search routine

Logically Collective

Input Parameters:
- `linesearch` - the `SNESLineSearch` context
- `func`       - [optional] function evaluation routine
- `ctx`        - [optional] user-defined context for private data for the function evaluation routine (may be `NULL`)

Calling sequence of `func`:
- `ls`        - the `SNESLineSearch` context
- `x`         - the current solution
- `d`         - the current search direction
- `w`         -  w = x + lambda*d  for some lambda
- `changed_d` - indicates if the search direction `d` has been changed
- `changed_w` - indicates `w` has been changed
- `ctx`       - the context passed to `SNESLineSearchSetPreCheck()`

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESGetLineSearch()`, `SNESLineSearchPostCheck()`, `SNESLineSearchSetPreCheck()`, `SNESLineSearchGetPreCheck()`, `SNESLineSearchGetPostCheck()`,
`SNESVISetVariableBounds()`, `SNESVISetComputeVariableBounds()`, `SNESSetFunctionDomainError()`, `SNESSetJacobianDomainError()`

# External Links
$(_doc_external("Snes/SNESLineSearchSetPostCheck"))
"""
function SNESLineSearchSetPostCheck(petsclib::PetscLibType, linesearch::SNESLineSearch, func::external, ctx::Cvoid) end

@for_petsc function SNESLineSearchSetPostCheck(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, func::external, ctx::Cvoid )

    @chk ccall(
               (:SNESLineSearchSetPostCheck, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, external, Ptr{Cvoid}),
               linesearch, func, ctx,
              )


	return nothing
end 

"""
	changed::PetscBool = SNESLineSearchPreCheck(petsclib::PetscLibType,linesearch::SNESLineSearch, X::PetscVec, Y::PetscVec) 
Prepares the line search for being applied.

Logically Collective

Input Parameters:
- `linesearch` - The linesearch instance.
- `X`          - The current solution
- `Y`          - The step direction

Output Parameter:
- `changed` - Indicator that the precheck routine has changed `Y`

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESGetLineSearch()`, `SNESLineSearchPostCheck()`, `SNESLineSearchSetPreCheck()`, `SNESLineSearchGetPreCheck()`, `SNESLineSearchSetPostCheck()`,
`SNESLineSearchGetPostCheck()`

# External Links
$(_doc_external("Snes/SNESLineSearchPreCheck"))
"""
function SNESLineSearchPreCheck(petsclib::PetscLibType, linesearch::SNESLineSearch, X::PetscVec, Y::PetscVec) end

@for_petsc function SNESLineSearchPreCheck(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, X::PetscVec, Y::PetscVec )
	changed_ = Ref{PetscBool}()

    @chk ccall(
               (:SNESLineSearchPreCheck, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, CVec, CVec, Ptr{PetscBool}),
               linesearch, X, Y, changed_,
              )

	changed = changed_[]

	return changed
end 

"""
	changed_Y::PetscBool,changed_W::PetscBool = SNESLineSearchPostCheck(petsclib::PetscLibType,linesearch::SNESLineSearch, X::PetscVec, Y::PetscVec, W::PetscVec) 
Hook to modify step direction or updated solution after a successful linesearch

Logically Collective

Input Parameters:
- `linesearch` - The line search context
- `X`          - The last solution
- `Y`          - The step direction
- `W`          - The updated solution, `W = X - lambda * Y` for some lambda

Output Parameters:
- `changed_Y` - Indicator if the direction `Y` has been changed.
- `changed_W` - Indicator if the new candidate solution `W` has been changed.

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESGetLineSearch()`, `SNESLineSearchPreCheck()`, `SNESLineSearchSetPostCheck()`, `SNESLineSearchGetPostCheck()`, `SNESLineSearchSetPrecheck()`, `SNESLineSearchGetPrecheck()`

# External Links
$(_doc_external("Snes/SNESLineSearchPostCheck"))
"""
function SNESLineSearchPostCheck(petsclib::PetscLibType, linesearch::SNESLineSearch, X::PetscVec, Y::PetscVec, W::PetscVec) end

@for_petsc function SNESLineSearchPostCheck(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, X::PetscVec, Y::PetscVec, W::PetscVec )
	changed_Y_ = Ref{PetscBool}()
	changed_W_ = Ref{PetscBool}()

    @chk ccall(
               (:SNESLineSearchPostCheck, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, CVec, CVec, CVec, Ptr{PetscBool}, Ptr{PetscBool}),
               linesearch, X, Y, W, changed_Y_, changed_W_,
              )

	changed_Y = changed_Y_[]
	changed_W = changed_W_[]

	return changed_Y,changed_W
end 

"""
	changed::PetscBool = SNESLineSearchPreCheckPicard(petsclib::PetscLibType,linesearch::SNESLineSearch, X::PetscVec, Y::PetscVec, ctx::Cvoid) 
Implements a correction that is sometimes useful to improve the convergence rate of Picard iteration {cite}`hindmarsh1996time`

Logically Collective

Input Parameters:
- `linesearch` - the line search context
- `X`          - base state for this step
- `ctx`        - context for this function

Input/Output Parameter:
- `Y` - correction, possibly modified

Output Parameter:
- `changed` - flag indicating that `Y` was modified

Options Database Keys:
- `-snes_linesearch_precheck_picard`       - activate this routine
- `-snes_linesearch_precheck_picard_angle` - angle

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESSetPicard()`, `SNESGetLineSearch()`, `SNESLineSearchSetPreCheck()`, `SNESLineSearchSetPostCheck()`

# External Links
$(_doc_external("Snes/SNESLineSearchPreCheckPicard"))
"""
function SNESLineSearchPreCheckPicard(petsclib::PetscLibType, linesearch::SNESLineSearch, X::PetscVec, Y::PetscVec, ctx::Cvoid) end

@for_petsc function SNESLineSearchPreCheckPicard(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, X::PetscVec, Y::PetscVec, ctx::Cvoid )
	changed_ = Ref{PetscBool}()

    @chk ccall(
               (:SNESLineSearchPreCheckPicard, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, CVec, CVec, Ptr{PetscBool}, Ptr{Cvoid}),
               linesearch, X, Y, changed_, ctx,
              )

	changed = changed_[]

	return changed
end 

"""
	fnorm::PetscReal = SNESLineSearchApply(petsclib::PetscLibType,linesearch::SNESLineSearch, X::PetscVec, F::PetscVec, Y::PetscVec) 
Computes the line

Collective

Input Parameter:
- `linesearch` - The line search context

Input/Output Parameters:
- `X`     - The current solution, on output the new solution
- `F`     - The current function value, on output the new function value at the solution value `X`
- `fnorm` - The current norm of `F`, on output the new norm of `F`
- `Y`     - The current search direction, on output the direction determined by the linesearch, i.e. `Xnew = Xold - lambda*Y`

Options Database Keys:
- `-snes_linesearch_type`                - basic (or equivalently none), bt, secant, cp, nleqerr, bisection, shell
- `-snes_linesearch_monitor [:filename]` - Print progress of line searches
- `-snes_linesearch_damping`             - The linesearch damping parameter, default is 1.0 (no damping)
- `-snes_linesearch_norms`               - Turn on/off the linesearch norms computation (SNESLineSearchSetComputeNorms())
- `-snes_linesearch_keeplambda`          - Keep the previous `lambda` as the initial guess
- `-snes_linesearch_max_it`              - The number of iterations for iterative line searches

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESGetLineSearch()`, `SNESLineSearchCreate()`, `SNESLineSearchGetLambda()`, `SNESLineSearchPreCheck()`, `SNESLineSearchPostCheck()`, `SNESSolve()`, `SNESComputeFunction()`, `SNESLineSearchSetComputeNorms()`,
`SNESLineSearchType`, `SNESLineSearchSetType()`

# External Links
$(_doc_external("Snes/SNESLineSearchApply"))
"""
function SNESLineSearchApply(petsclib::PetscLibType, linesearch::SNESLineSearch, X::PetscVec, F::PetscVec, Y::PetscVec) end

@for_petsc function SNESLineSearchApply(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, X::PetscVec, F::PetscVec, Y::PetscVec )
	fnorm_ = Ref{$PetscReal}()

    @chk ccall(
               (:SNESLineSearchApply, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, CVec, CVec, Ptr{$PetscReal}, CVec),
               linesearch, X, F, fnorm_, Y,
              )

	fnorm = fnorm_[]

	return fnorm
end 

"""
	SNESLineSearchDestroy(petsclib::PetscLibType,linesearch::SNESLineSearch) 
Destroys the line search instance.

Collective

Input Parameter:
- `linesearch` - The line search context

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESGetLineSearch()`, `SNESLineSearchCreate()`, `SNESLineSearchReset()`, `SNESDestroy()`

# External Links
$(_doc_external("Snes/SNESLineSearchDestroy"))
"""
function SNESLineSearchDestroy(petsclib::PetscLibType, linesearch::SNESLineSearch) end

@for_petsc function SNESLineSearchDestroy(petsclib::$UnionPetscLib, linesearch::SNESLineSearch )

    @chk ccall(
               (:SNESLineSearchDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{SNESLineSearch},),
               linesearch,
              )


	return nothing
end 

"""
	SNESLineSearchSetDefaultMonitor(petsclib::PetscLibType,linesearch::SNESLineSearch, viewer::PetscViewer) 
Turns on/off printing useful information and debugging output about the line search.

Logically Collective

Input Parameters:
- `linesearch` - the linesearch object
- `viewer`     - an `PETSCVIEWERASCII` `PetscViewer` or `NULL` to turn off monitor

Options Database Key:
- `-snes_linesearch_monitor [:filename]` - enables the monitor

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `PETSCVIEWERASCII`, `SNESGetLineSearch()`, `SNESLineSearchGetDefaultMonitor()`, `PetscViewer`, `SNESLineSearchSetMonitor()`,
`SNESLineSearchMonitorSetFromOptions()`

# External Links
$(_doc_external("Snes/SNESLineSearchSetDefaultMonitor"))
"""
function SNESLineSearchSetDefaultMonitor(petsclib::PetscLibType, linesearch::SNESLineSearch, viewer::PetscViewer) end

@for_petsc function SNESLineSearchSetDefaultMonitor(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, viewer::PetscViewer )

    @chk ccall(
               (:SNESLineSearchSetDefaultMonitor, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, PetscViewer),
               linesearch, viewer,
              )


	return nothing
end 

"""
	SNESLineSearchGetDefaultMonitor(petsclib::PetscLibType,linesearch::SNESLineSearch, monitor::PetscViewer) 
Gets the `PetscViewer` instance for the default line search monitor that is turned on with `SNESLineSearchSetDefaultMonitor()`

Logically Collective

Input Parameter:
- `linesearch` - the line search context

Output Parameter:
- `monitor` - monitor context

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESGetLineSearch()`, `SNESLineSearchSetDefaultMonitor()`, `PetscViewer`

# External Links
$(_doc_external("Snes/SNESLineSearchGetDefaultMonitor"))
"""
function SNESLineSearchGetDefaultMonitor(petsclib::PetscLibType, linesearch::SNESLineSearch, monitor::PetscViewer) end

@for_petsc function SNESLineSearchGetDefaultMonitor(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, monitor::PetscViewer )

    @chk ccall(
               (:SNESLineSearchGetDefaultMonitor, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, Ptr{PetscViewer}),
               linesearch, monitor,
              )


	return nothing
end 

"""
	SNESLineSearchMonitorSetFromOptions(petsclib::PetscLibType,ls::SNESLineSearch, name::String, help::String, manual::String, monitor::external, monitorsetup::external) 
Sets a monitor function and viewer appropriate for the type indicated in the options database

Collective

Input Parameters:
- `ls`           - `SNESLineSearch` object to monitor
- `name`         - the monitor type
- `help`         - message indicating what monitoring is done
- `manual`       - manual page for the monitor
- `monitor`      - the monitor function, must use `PetscViewerAndFormat` as its context
- `monitorsetup` - a function that is called once ONLY if the user selected this monitor that may set additional features of the `SNESLineSearch` or `PetscViewer`

Calling sequence of `monitor`:
- `ls` - `SNESLineSearch` object being monitored
- `vf` - a `PetscViewerAndFormat` struct that provides the `PetscViewer` and `PetscViewerFormat` being used

Calling sequence of `monitorsetup`:
- `ls` - `SNESLineSearch` object being monitored
- `vf` - a `PetscViewerAndFormat` struct that provides the `PetscViewer` and `PetscViewerFormat` being used

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESLineSearchSetMonitor()`, `PetscOptionsCreateViewer()`, `PetscOptionsGetReal()`, `PetscOptionsHasName()`, `PetscOptionsGetString()`,
`PetscOptionsGetIntArray()`, `PetscOptionsGetRealArray()`, `PetscOptionsBool()`
`PetscOptionsInt()`, `PetscOptionsString()`, `PetscOptionsReal()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`

# External Links
$(_doc_external("Snes/SNESLineSearchMonitorSetFromOptions"))
"""
function SNESLineSearchMonitorSetFromOptions(petsclib::PetscLibType, ls::SNESLineSearch, name::String, help::String, manual::String, monitor::external, monitorsetup::external) end

@for_petsc function SNESLineSearchMonitorSetFromOptions(petsclib::$UnionPetscLib, ls::SNESLineSearch, name::String, help::String, manual::String, monitor::external, monitorsetup::external )

    @chk ccall(
               (:SNESLineSearchMonitorSetFromOptions, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, external, external),
               ls, name, help, manual, monitor, monitorsetup,
              )


	return nothing
end 

"""
	SNESLineSearchSetFromOptions(petsclib::PetscLibType,linesearch::SNESLineSearch) 
Sets options for the line search

Logically Collective

Input Parameter:
- `linesearch` - a `SNESLineSearch` line search context

Options Database Keys:
- `-snes_linesearch_type <type>`                                      - basic (or equivalently none), `bt`, `secant`, `cp`, `nleqerr`, `bisection`, `shell`
- `-snes_linesearch_order <order>`                                    - 1, 2, 3.  Most types only support certain orders (`bt` supports 1, 2 or 3)
- `-snes_linesearch_norms`                                            - Turn on/off the linesearch norms for the basic linesearch typem (`SNESLineSearchSetComputeNorms()`)
- `-snes_linesearch_minlambda`                                        - The minimum `lambda`
- `-snes_linesearch_maxlambda`                                        - The maximum `lambda`
- `-snes_linesearch_rtol`                                             - Relative tolerance for iterative line searches
- `-snes_linesearch_atol`                                             - Absolute tolerance for iterative line searches
- `-snes_linesearch_ltol`                                             - Change in `lambda` tolerance for iterative line searches
- `-snes_linesearch_max_it`                                           - The number of iterations for iterative line searches
- `-snes_linesearch_monitor [:filename]`                              - Print progress of line searches
- `-snes_linesearch_monitor_solution_update [viewer:filename:format]` - view each update tried by line search routine
- `-snes_linesearch_damping`                                          - The linesearch damping parameter
- `-snes_linesearch_keeplambda`                                       - Keep the previous `lambda` as the initial guess.
- `-snes_linesearch_precheck_picard`                                  - Use precheck that speeds up convergence of picard method
- `-snes_linesearch_precheck_picard_angle`                            - Angle used in Picard precheck method

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESGetLineSearch()`, `SNESLineSearchCreate()`, `SNESLineSearchSetOrder()`, `SNESLineSearchSetType()`, `SNESLineSearchSetTolerances()`, `SNESLineSearchSetDamping()`, `SNESLineSearchPreCheckPicard()`,
`SNESLineSearchType`, `SNESLineSearchSetComputeNorms()`

# External Links
$(_doc_external("Snes/SNESLineSearchSetFromOptions"))
"""
function SNESLineSearchSetFromOptions(petsclib::PetscLibType, linesearch::SNESLineSearch) end

@for_petsc function SNESLineSearchSetFromOptions(petsclib::$UnionPetscLib, linesearch::SNESLineSearch )

    @chk ccall(
               (:SNESLineSearchSetFromOptions, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch,),
               linesearch,
              )


	return nothing
end 

"""
	SNESLineSearchView(petsclib::PetscLibType,linesearch::SNESLineSearch, viewer::PetscViewer) 
Prints useful information about the line search

Logically Collective

Input Parameters:
- `linesearch` - line search context
- `viewer`     - the `PetscViewer` to display the line search information to

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `PetscViewer`, `SNESLineSearchCreate()`

# External Links
$(_doc_external("Snes/SNESLineSearchView"))
"""
function SNESLineSearchView(petsclib::PetscLibType, linesearch::SNESLineSearch, viewer::PetscViewer) end

@for_petsc function SNESLineSearchView(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, viewer::PetscViewer )

    @chk ccall(
               (:SNESLineSearchView, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, PetscViewer),
               linesearch, viewer,
              )


	return nothing
end 

"""
	type::SNESLineSearchType = SNESLineSearchGetType(petsclib::PetscLibType,linesearch::SNESLineSearch) 
Gets the `SNESLinesearchType` of a `SNESLineSearch`

Logically Collective

Input Parameter:
- `linesearch` - the line search context

Output Parameter:
- `type` - The type of line search, or `NULL` if not set

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESLineSearchType`, `SNESLineSearchCreate()`, `SNESLineSearchSetFromOptions()`, `SNESLineSearchSetType()`

# External Links
$(_doc_external("Snes/SNESLineSearchGetType"))
"""
function SNESLineSearchGetType(petsclib::PetscLibType, linesearch::SNESLineSearch) end

@for_petsc function SNESLineSearchGetType(petsclib::$UnionPetscLib, linesearch::SNESLineSearch )
	type_ = Ref{SNESLineSearchType}()

    @chk ccall(
               (:SNESLineSearchGetType, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, Ptr{SNESLineSearchType}),
               linesearch, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	SNESLineSearchSetType(petsclib::PetscLibType,linesearch::SNESLineSearch, type::SNESLineSearchType) 
Sets the `SNESLinesearchType` of a `SNESLineSearch` object to indicate the line search algorithm that should be used by a given `SNES` solver

Logically Collective

Input Parameters:
- `linesearch` - the line search context
- `type`       - The type of line search to be used, see `SNESLineSearchType`

Options Database Key:
- `-snes_linesearch_type <type>` - basic (or equivalently none), bt, secant, cp, nleqerr, bisection, shell

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESLineSearchType`, `SNESLineSearchCreate()`, `SNESLineSearchSetFromOptions()`, `SNESLineSearchGetType()`,
`SNESGetLineSearch()`

# External Links
$(_doc_external("Snes/SNESLineSearchSetType"))
"""
function SNESLineSearchSetType(petsclib::PetscLibType, linesearch::SNESLineSearch, type::SNESLineSearchType) end

@for_petsc function SNESLineSearchSetType(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, type::SNESLineSearchType )

    @chk ccall(
               (:SNESLineSearchSetType, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, SNESLineSearchType),
               linesearch, type,
              )


	return nothing
end 

"""
	SNESLineSearchSetSNES(petsclib::PetscLibType,linesearch::SNESLineSearch, snes::PetscSNES) 
Sets the `SNES` for the linesearch for function evaluation.

Input Parameters:
- `linesearch` - the line search context
- `snes`       - The `SNES` instance

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESLineSearchGetSNES()`, `SNESLineSearchSetVecs()`

# External Links
$(_doc_external("Snes/SNESLineSearchSetSNES"))
"""
function SNESLineSearchSetSNES(petsclib::PetscLibType, linesearch::SNESLineSearch, snes::PetscSNES) end

@for_petsc function SNESLineSearchSetSNES(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, snes::PetscSNES )

    @chk ccall(
               (:SNESLineSearchSetSNES, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, CSNES),
               linesearch, snes,
              )


	return nothing
end 

"""
	SNESLineSearchGetSNES(petsclib::PetscLibType,linesearch::SNESLineSearch, snes::PetscSNES) 
Gets the `SNES` instance associated with the line search.

Not Collective

Input Parameter:
- `linesearch` - the line search context

Output Parameter:
- `snes` - The `SNES` instance

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESType`, `SNESLineSearchSetVecs()`

# External Links
$(_doc_external("Snes/SNESLineSearchGetSNES"))
"""
function SNESLineSearchGetSNES(petsclib::PetscLibType, linesearch::SNESLineSearch, snes::PetscSNES) end

@for_petsc function SNESLineSearchGetSNES(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, snes::PetscSNES )
	snes_ = Ref(snes.ptr)

    @chk ccall(
               (:SNESLineSearchGetSNES, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, Ptr{CSNES}),
               linesearch, snes_,
              )

	snes.ptr = C_NULL

	return nothing
end 

"""
	lambda::PetscReal = SNESLineSearchGetLambda(petsclib::PetscLibType,linesearch::SNESLineSearch) 
Gets the last line search `lambda` used

Not Collective

Input Parameter:
- `linesearch` - the line search context

Output Parameter:
- `lambda` - The last `lambda` (scaling of the solution update) computed during `SNESLineSearchApply()`

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESLineSearchSetLambda()`, `SNESLineSearchGetDamping()`, `SNESLineSearchApply()`

# External Links
$(_doc_external("Snes/SNESLineSearchGetLambda"))
"""
function SNESLineSearchGetLambda(petsclib::PetscLibType, linesearch::SNESLineSearch) end

@for_petsc function SNESLineSearchGetLambda(petsclib::$UnionPetscLib, linesearch::SNESLineSearch )
	lambda_ = Ref{$PetscReal}()

    @chk ccall(
               (:SNESLineSearchGetLambda, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, Ptr{$PetscReal}),
               linesearch, lambda_,
              )

	lambda = lambda_[]

	return lambda
end 

"""
	SNESLineSearchSetLambda(petsclib::PetscLibType,linesearch::SNESLineSearch, lambda::PetscReal) 
Sets the line search `lambda` (scaling of the solution update)

Input Parameters:
- `linesearch` - line search context
- `lambda`     - The `lambda` to use

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESLineSearchGetLambda()`

# External Links
$(_doc_external("Snes/SNESLineSearchSetLambda"))
"""
function SNESLineSearchSetLambda(petsclib::PetscLibType, linesearch::SNESLineSearch, lambda::PetscReal) end

@for_petsc function SNESLineSearchSetLambda(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, lambda::$PetscReal )

    @chk ccall(
               (:SNESLineSearchSetLambda, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, $PetscReal),
               linesearch, lambda,
              )


	return nothing
end 

"""
	minlambda::PetscReal,maxlambda::PetscReal,rtol::PetscReal,atol::PetscReal,ltol::PetscReal,max_it::PetscInt = SNESLineSearchGetTolerances(petsclib::PetscLibType,linesearch::SNESLineSearch) 
Gets the tolerances for the line search.

Not Collective

Input Parameter:
- `linesearch` - the line search context

Output Parameters:
- `minlambda` - The minimum `lambda` allowed
- `maxlambda` - The maximum `lambda` allowed
- `rtol`      - The relative tolerance for iterative line searches
- `atol`      - The absolute tolerance for iterative line searches
- `ltol`      - The change in `lambda` tolerance for iterative line searches
- `max_it`    - The maximum number of iterations of the line search

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESLineSearchSetTolerances()`

# External Links
$(_doc_external("Snes/SNESLineSearchGetTolerances"))
"""
function SNESLineSearchGetTolerances(petsclib::PetscLibType, linesearch::SNESLineSearch) end

@for_petsc function SNESLineSearchGetTolerances(petsclib::$UnionPetscLib, linesearch::SNESLineSearch )
	minlambda_ = Ref{$PetscReal}()
	maxlambda_ = Ref{$PetscReal}()
	rtol_ = Ref{$PetscReal}()
	atol_ = Ref{$PetscReal}()
	ltol_ = Ref{$PetscReal}()
	max_it_ = Ref{$PetscInt}()

    @chk ccall(
               (:SNESLineSearchGetTolerances, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscInt}),
               linesearch, minlambda_, maxlambda_, rtol_, atol_, ltol_, max_it_,
              )

	minlambda = minlambda_[]
	maxlambda = maxlambda_[]
	rtol = rtol_[]
	atol = atol_[]
	ltol = ltol_[]
	max_it = max_it_[]

	return minlambda,maxlambda,rtol,atol,ltol,max_it
end 

"""
	SNESLineSearchSetTolerances(petsclib::PetscLibType,linesearch::SNESLineSearch, minlambda::PetscReal, maxlambda::PetscReal, rtol::PetscReal, atol::PetscReal, ltol::PetscReal, max_it::PetscInt) 
Sets the tolerances for the linesearch.

Collective

Input Parameters:
- `linesearch` - the line search context
- `minlambda`  - The minimum `lambda` allowed
- `maxlambda`  - The maximum `lambda` allowed
- `rtol`       - The relative tolerance for iterative line searches
- `atol`       - The absolute tolerance for iterative line searches
- `ltol`       - The change in `lambda` tolerance for iterative line searches
- `max_it`     - The maximum number of iterations of the line search

Options Database Keys:
- `-snes_linesearch_minlambda` - The minimum `lambda` allowed
- `-snes_linesearch_maxlambda` - The maximum `lambda` allowed
- `-snes_linesearch_rtol`      - Relative tolerance for iterative line searches
- `-snes_linesearch_atol`      - Absolute tolerance for iterative line searches
- `-snes_linesearch_ltol`      - Change in `lambda` tolerance for iterative line searches
- `-snes_linesearch_max_it`    - The number of iterations for iterative line searches

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESLineSearchGetTolerances()`

# External Links
$(_doc_external("Snes/SNESLineSearchSetTolerances"))
"""
function SNESLineSearchSetTolerances(petsclib::PetscLibType, linesearch::SNESLineSearch, minlambda::PetscReal, maxlambda::PetscReal, rtol::PetscReal, atol::PetscReal, ltol::PetscReal, max_it::PetscInt) end

@for_petsc function SNESLineSearchSetTolerances(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, minlambda::$PetscReal, maxlambda::$PetscReal, rtol::$PetscReal, atol::$PetscReal, ltol::$PetscReal, max_it::$PetscInt )

    @chk ccall(
               (:SNESLineSearchSetTolerances, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, $PetscReal, $PetscReal, $PetscReal, $PetscReal, $PetscReal, $PetscInt),
               linesearch, minlambda, maxlambda, rtol, atol, ltol, max_it,
              )


	return nothing
end 

"""
	damping::PetscReal = SNESLineSearchGetDamping(petsclib::PetscLibType,linesearch::SNESLineSearch) 
Gets the line search damping parameter.

Input Parameter:
- `linesearch` - the line search context

Output Parameter:
- `damping` - The damping parameter

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESLineSearchGetStepTolerance()`, `SNESQN`

# External Links
$(_doc_external("Snes/SNESLineSearchGetDamping"))
"""
function SNESLineSearchGetDamping(petsclib::PetscLibType, linesearch::SNESLineSearch) end

@for_petsc function SNESLineSearchGetDamping(petsclib::$UnionPetscLib, linesearch::SNESLineSearch )
	damping_ = Ref{$PetscReal}()

    @chk ccall(
               (:SNESLineSearchGetDamping, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, Ptr{$PetscReal}),
               linesearch, damping_,
              )

	damping = damping_[]

	return damping
end 

"""
	SNESLineSearchSetDamping(petsclib::PetscLibType,linesearch::SNESLineSearch, damping::PetscReal) 
Sets the line search damping parameter.

Input Parameters:
- `linesearch` - the line search context
- `damping`    - The damping parameter

Options Database Key:
- `-snes_linesearch_damping <damping>` - the damping value

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESLineSearchGetDamping()`

# External Links
$(_doc_external("Snes/SNESLineSearchSetDamping"))
"""
function SNESLineSearchSetDamping(petsclib::PetscLibType, linesearch::SNESLineSearch, damping::PetscReal) end

@for_petsc function SNESLineSearchSetDamping(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, damping::$PetscReal )

    @chk ccall(
               (:SNESLineSearchSetDamping, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, $PetscReal),
               linesearch, damping,
              )


	return nothing
end 

"""
	order::PetscInt = SNESLineSearchGetOrder(petsclib::PetscLibType,linesearch::SNESLineSearch) 
Gets the line search approximation order.

Input Parameter:
- `linesearch` - the line search context

Output Parameter:
- `order` - The order

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESLineSearchSetOrder()`

# External Links
$(_doc_external("Snes/SNESLineSearchGetOrder"))
"""
function SNESLineSearchGetOrder(petsclib::PetscLibType, linesearch::SNESLineSearch) end

@for_petsc function SNESLineSearchGetOrder(petsclib::$UnionPetscLib, linesearch::SNESLineSearch )
	order_ = Ref{$PetscInt}()

    @chk ccall(
               (:SNESLineSearchGetOrder, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, Ptr{$PetscInt}),
               linesearch, order_,
              )

	order = order_[]

	return order
end 

"""
	SNESLineSearchSetOrder(petsclib::PetscLibType,linesearch::SNESLineSearch, order::PetscInt) 
Sets the maximum order of the polynomial fit used in the line search

Input Parameters:
- `linesearch` - the line search context
- `order`      - The order

Level: intermediate

Values for `order`:
- `1 or `SNES_LINESEARCH_ORDER_LINEAR`  - linear order
- `2 or `SNES_LINESEARCH_ORDER_QUADRATIC`  - quadratic order
- `3 or `SNES_LINESEARCH_ORDER_CUBIC`  - cubic order

Options Database Key:
- `-snes_linesearch_order <order>` - 1, 2, 3.  Most types only support certain orders (`SNESLINESEARCHBT` supports 2 or 3)

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESLineSearchGetOrder()`, `SNESLineSearchSetDamping()`

# External Links
$(_doc_external("Snes/SNESLineSearchSetOrder"))
"""
function SNESLineSearchSetOrder(petsclib::PetscLibType, linesearch::SNESLineSearch, order::PetscInt) end

@for_petsc function SNESLineSearchSetOrder(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, order::$PetscInt )

    @chk ccall(
               (:SNESLineSearchSetOrder, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, $PetscInt),
               linesearch, order,
              )


	return nothing
end 

"""
	xnorm::PetscReal,fnorm::PetscReal,ynorm::PetscReal = SNESLineSearchGetNorms(petsclib::PetscLibType,linesearch::SNESLineSearch) 
Gets the norms for the current solution `X`, the current update `Y`, and the current function value `F`.

Not Collective

Input Parameter:
- `linesearch` - the line search context

Output Parameters:
- `xnorm` - The norm of the current solution
- `fnorm` - The norm of the current function, this is the `norm(function(X))` where `X` is the current solution.
- `ynorm` - The norm of the current update (after scaling by the linesearch computed `lambda`)

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESLineSearchSetNorms()`, `SNESLineSearchGetVecs()`

# External Links
$(_doc_external("Snes/SNESLineSearchGetNorms"))
"""
function SNESLineSearchGetNorms(petsclib::PetscLibType, linesearch::SNESLineSearch) end

@for_petsc function SNESLineSearchGetNorms(petsclib::$UnionPetscLib, linesearch::SNESLineSearch )
	xnorm_ = Ref{$PetscReal}()
	fnorm_ = Ref{$PetscReal}()
	ynorm_ = Ref{$PetscReal}()

    @chk ccall(
               (:SNESLineSearchGetNorms, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               linesearch, xnorm_, fnorm_, ynorm_,
              )

	xnorm = xnorm_[]
	fnorm = fnorm_[]
	ynorm = ynorm_[]

	return xnorm,fnorm,ynorm
end 

"""
	SNESLineSearchSetNorms(petsclib::PetscLibType,linesearch::SNESLineSearch, xnorm::PetscReal, fnorm::PetscReal, ynorm::PetscReal) 
Sets the computed norms for the current solution `X`, the current update `Y`, and the current function value `F`.

Collective

Input Parameters:
- `linesearch` - the line search context
- `xnorm`      - The norm of the current solution
- `fnorm`      - The norm of the current function, this is the `norm(function(X))` where `X` is the current solution
- `ynorm`      - The norm of the current update (after scaling by the linesearch computed `lambda`)

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESLineSearchGetNorms()`, `SNESLineSearchSetVecs()`

# External Links
$(_doc_external("Snes/SNESLineSearchSetNorms"))
"""
function SNESLineSearchSetNorms(petsclib::PetscLibType, linesearch::SNESLineSearch, xnorm::PetscReal, fnorm::PetscReal, ynorm::PetscReal) end

@for_petsc function SNESLineSearchSetNorms(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, xnorm::$PetscReal, fnorm::$PetscReal, ynorm::$PetscReal )

    @chk ccall(
               (:SNESLineSearchSetNorms, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, $PetscReal, $PetscReal, $PetscReal),
               linesearch, xnorm, fnorm, ynorm,
              )


	return nothing
end 

"""
	SNESLineSearchComputeNorms(petsclib::PetscLibType,linesearch::SNESLineSearch) 
Explicitly computes the norms of the current solution `X`, the current update `Y`, and the current function value `F`.

Input Parameter:
- `linesearch` - the line search context

Options Database Key:
- `-snes_linesearch_norms` - turn norm computation on or off

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESLineSearchGetNorms`, `SNESLineSearchSetNorms()`, `SNESLineSearchSetComputeNorms()`

# External Links
$(_doc_external("Snes/SNESLineSearchComputeNorms"))
"""
function SNESLineSearchComputeNorms(petsclib::PetscLibType, linesearch::SNESLineSearch) end

@for_petsc function SNESLineSearchComputeNorms(petsclib::$UnionPetscLib, linesearch::SNESLineSearch )

    @chk ccall(
               (:SNESLineSearchComputeNorms, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch,),
               linesearch,
              )


	return nothing
end 

"""
	SNESLineSearchSetComputeNorms(petsclib::PetscLibType,linesearch::SNESLineSearch, flg::PetscBool) 
Turns on or off the computation of final norms in the line search.

Input Parameters:
- `linesearch` - the line search context
- `flg`        - indicates whether or not to compute norms

Options Database Key:
- `-snes_linesearch_norms <true>` - Turns on/off computation of the norms for basic (none) `SNESLINESEARCHBASIC` line search

Level: intermediate

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESLineSearchGetNorms()`, `SNESLineSearchSetNorms()`, `SNESLineSearchComputeNorms()`, `SNESLINESEARCHBASIC`

# External Links
$(_doc_external("Snes/SNESLineSearchSetComputeNorms"))
"""
function SNESLineSearchSetComputeNorms(petsclib::PetscLibType, linesearch::SNESLineSearch, flg::PetscBool) end

@for_petsc function SNESLineSearchSetComputeNorms(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, flg::PetscBool )

    @chk ccall(
               (:SNESLineSearchSetComputeNorms, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, PetscBool),
               linesearch, flg,
              )


	return nothing
end 

"""
	SNESLineSearchGetVecs(petsclib::PetscLibType,linesearch::SNESLineSearch, X::PetscVec, F::PetscVec, Y::PetscVec, W::PetscVec, G::PetscVec) 
Gets the vectors from the `SNESLineSearch` context

Not Collective but the vectors are parallel

Input Parameter:
- `linesearch` - the line search context

Output Parameters:
- `X` - Solution vector
- `F` - Function vector
- `Y` - Search direction vector
- `W` - Solution work vector
- `G` - Function work vector

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESLineSearchGetNorms()`, `SNESLineSearchSetVecs()`

# External Links
$(_doc_external("Snes/SNESLineSearchGetVecs"))
"""
function SNESLineSearchGetVecs(petsclib::PetscLibType, linesearch::SNESLineSearch, X::PetscVec, F::PetscVec, Y::PetscVec, W::PetscVec, G::PetscVec) end

@for_petsc function SNESLineSearchGetVecs(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, X::PetscVec, F::PetscVec, Y::PetscVec, W::PetscVec, G::PetscVec )
	X_ = Ref(X.ptr)
	F_ = Ref(F.ptr)
	Y_ = Ref(Y.ptr)
	W_ = Ref(W.ptr)
	G_ = Ref(G.ptr)

    @chk ccall(
               (:SNESLineSearchGetVecs, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, Ptr{CVec}, Ptr{CVec}, Ptr{CVec}, Ptr{CVec}, Ptr{CVec}),
               linesearch, X_, F_, Y_, W_, G_,
              )

	X.ptr = C_NULL
	F.ptr = C_NULL
	Y.ptr = C_NULL
	W.ptr = C_NULL
	G.ptr = C_NULL

	return nothing
end 

"""
	SNESLineSearchSetVecs(petsclib::PetscLibType,linesearch::SNESLineSearch, X::PetscVec, F::PetscVec, Y::PetscVec, W::PetscVec, G::PetscVec) 
Sets the vectors on the `SNESLineSearch` context

Logically Collective

Input Parameters:
- `linesearch` - the line search context
- `X`          - Solution vector
- `F`          - Function vector
- `Y`          - Search direction vector
- `W`          - Solution work vector
- `G`          - Function work vector

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESLineSearchSetNorms()`, `SNESLineSearchGetVecs()`

# External Links
$(_doc_external("Snes/SNESLineSearchSetVecs"))
"""
function SNESLineSearchSetVecs(petsclib::PetscLibType, linesearch::SNESLineSearch, X::PetscVec, F::PetscVec, Y::PetscVec, W::PetscVec, G::PetscVec) end

@for_petsc function SNESLineSearchSetVecs(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, X::PetscVec, F::PetscVec, Y::PetscVec, W::PetscVec, G::PetscVec )

    @chk ccall(
               (:SNESLineSearchSetVecs, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, CVec, CVec, CVec, CVec, CVec),
               linesearch, X, F, Y, W, G,
              )


	return nothing
end 

"""
	SNESLineSearchAppendOptionsPrefix(petsclib::PetscLibType,linesearch::SNESLineSearch, prefix::String) 
Appends to the prefix used for searching for all
`SNESLineSearch` options in the database.

Logically Collective

Input Parameters:
- `linesearch` - the `SNESLineSearch` context
- `prefix`     - the prefix to prepend to all option names

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESLineSearch()`, `SNESLineSearchSetFromOptions()`, `SNESGetOptionsPrefix()`

# External Links
$(_doc_external("Snes/SNESLineSearchAppendOptionsPrefix"))
"""
function SNESLineSearchAppendOptionsPrefix(petsclib::PetscLibType, linesearch::SNESLineSearch, prefix::String) end

@for_petsc function SNESLineSearchAppendOptionsPrefix(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, prefix::String )

    @chk ccall(
               (:SNESLineSearchAppendOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, Ptr{Cchar}),
               linesearch, prefix,
              )


	return nothing
end 

"""
	SNESLineSearchGetOptionsPrefix(petsclib::PetscLibType,linesearch::SNESLineSearch, prefix::String) 
Gets the prefix used for searching for all
SNESLineSearch options in the database.

Not Collective

Input Parameter:
- `linesearch` - the `SNESLineSearch` context

Output Parameter:
- `prefix` - pointer to the prefix string used

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESAppendOptionsPrefix()`

# External Links
$(_doc_external("Snes/SNESLineSearchGetOptionsPrefix"))
"""
function SNESLineSearchGetOptionsPrefix(petsclib::PetscLibType, linesearch::SNESLineSearch, prefix::String) end

@for_petsc function SNESLineSearchGetOptionsPrefix(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, prefix::String )
	prefix_ = Ref(pointer(prefix))

    @chk ccall(
               (:SNESLineSearchGetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, Ptr{Ptr{Cchar}}),
               linesearch, prefix_,
              )


	return nothing
end 

"""
	SNESLineSearchSetWorkVecs(petsclib::PetscLibType,linesearch::SNESLineSearch, nwork::PetscInt) 
Sets work vectors for the line search.

Input Parameters:
- `linesearch` - the `SNESLineSearch` context
- `nwork`      - the number of work vectors

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESSetWorkVecs()`

# External Links
$(_doc_external("Snes/SNESLineSearchSetWorkVecs"))
"""
function SNESLineSearchSetWorkVecs(petsclib::PetscLibType, linesearch::SNESLineSearch, nwork::PetscInt) end

@for_petsc function SNESLineSearchSetWorkVecs(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, nwork::$PetscInt )

    @chk ccall(
               (:SNESLineSearchSetWorkVecs, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, $PetscInt),
               linesearch, nwork,
              )


	return nothing
end 

"""
	SNESLineSearchGetReason(petsclib::PetscLibType,linesearch::SNESLineSearch, result::SNESLineSearchReason) 
Gets the success/failure status of the last line search application

Input Parameter:
- `linesearch` - the line search context

Output Parameter:
- `result` - The success or failure status

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESLineSearchSetReason()`, `SNESLineSearchReason`

# External Links
$(_doc_external("Snes/SNESLineSearchGetReason"))
"""
function SNESLineSearchGetReason(petsclib::PetscLibType, linesearch::SNESLineSearch, result::SNESLineSearchReason) end

@for_petsc function SNESLineSearchGetReason(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, result::SNESLineSearchReason )

    @chk ccall(
               (:SNESLineSearchGetReason, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, Ptr{SNESLineSearchReason}),
               linesearch, result,
              )


	return nothing
end 

"""
	SNESLineSearchSetReason(petsclib::PetscLibType,linesearch::SNESLineSearch, result::SNESLineSearchReason) 
Sets the success/failure status of the line search application

Logically Collective; No Fortran Support

Input Parameters:
- `linesearch` - the line search context
- `result`     - The success or failure status

Level: developer

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESLineSearchReason`, `SNESLineSearchGetSResult()`

# External Links
$(_doc_external("Snes/SNESLineSearchSetReason"))
"""
function SNESLineSearchSetReason(petsclib::PetscLibType, linesearch::SNESLineSearch, result::SNESLineSearchReason) end

@for_petsc function SNESLineSearchSetReason(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, result::SNESLineSearchReason )

    @chk ccall(
               (:SNESLineSearchSetReason, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, SNESLineSearchReason),
               linesearch, result,
              )


	return nothing
end 

"""
	SNESLineSearchSetVIFunctions(petsclib::PetscLibType,linesearch::SNESLineSearch, projectfunc::SNESLineSearchVIProjectFn, normfunc::SNESLineSearchVINormFn, dirderivfunc::SNESLineSearchVIDirDerivFn) 
Sets VI

Logically Collective

Input Parameters:
- `linesearch`   - the linesearch object
- `projectfunc`  - function for projecting the function to the bounds, see `SNESLineSearchVIProjectFn` for calling sequence
- `normfunc`     - function for computing the norm of an active set, see `SNESLineSearchVINormFn` for calling sequence
- `dirderivfunc` - function for computing the directional derivative of an active set, see `SNESLineSearchVIDirDerivFn` for calling sequence

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESLineSearchGetVIFunctions()`, `SNESLineSearchSetPostCheck()`, `SNESLineSearchSetPreCheck()`,
`SNESLineSearchVIProjectFn`, `SNESLineSearchVINormFn`, `SNESLineSearchVIDirDerivFn`

# External Links
$(_doc_external("Snes/SNESLineSearchSetVIFunctions"))
"""
function SNESLineSearchSetVIFunctions(petsclib::PetscLibType, linesearch::SNESLineSearch, projectfunc::SNESLineSearchVIProjectFn, normfunc::SNESLineSearchVINormFn, dirderivfunc::SNESLineSearchVIDirDerivFn) end

@for_petsc function SNESLineSearchSetVIFunctions(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, projectfunc::SNESLineSearchVIProjectFn, normfunc::SNESLineSearchVINormFn, dirderivfunc::SNESLineSearchVIDirDerivFn )

    @chk ccall(
               (:SNESLineSearchSetVIFunctions, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, Ptr{SNESLineSearchVIProjectFn}, Ptr{SNESLineSearchVINormFn}, Ptr{SNESLineSearchVIDirDerivFn}),
               linesearch, projectfunc, normfunc, dirderivfunc,
              )


	return nothing
end 

"""
	SNESLineSearchGetVIFunctions(petsclib::PetscLibType,linesearch::SNESLineSearch, projectfunc::SNESLineSearchVIProjectFn, normfunc::SNESLineSearchVINormFn, dirderivfunc::SNESLineSearchVIDirDerivFn) 
Sets VI

Not Collective

Input Parameter:
- `linesearch` - the line search context, obtain with `SNESGetLineSearch()`

Output Parameters:
- `projectfunc`  - function for projecting the function to the bounds, see `SNESLineSearchVIProjectFn` for calling sequence
- `normfunc`     - function for computing the norm of an active set, see `SNESLineSearchVINormFn ` for calling sequence
- `dirderivfunc` - function for computing the directional derivative of an active set, see `SNESLineSearchVIDirDerivFn` for calling sequence

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESLineSearchSetVIFunctions()`, `SNESLineSearchGetPostCheck()`, `SNESLineSearchGetPreCheck()`,
`SNESLineSearchVIProjectFn`, `SNESLineSearchVINormFn`

# External Links
$(_doc_external("Snes/SNESLineSearchGetVIFunctions"))
"""
function SNESLineSearchGetVIFunctions(petsclib::PetscLibType, linesearch::SNESLineSearch, projectfunc::SNESLineSearchVIProjectFn, normfunc::SNESLineSearchVINormFn, dirderivfunc::SNESLineSearchVIDirDerivFn) end

@for_petsc function SNESLineSearchGetVIFunctions(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, projectfunc::SNESLineSearchVIProjectFn, normfunc::SNESLineSearchVINormFn, dirderivfunc::SNESLineSearchVIDirDerivFn )

    @chk ccall(
               (:SNESLineSearchGetVIFunctions, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, SNESLineSearchVIProjectFn, SNESLineSearchVINormFn, SNESLineSearchVIDirDerivFn),
               linesearch, projectfunc, normfunc, dirderivfunc,
              )


	return nothing
end 

"""
	SNESLineSearchRegister(petsclib::PetscLibType,sname::String, fnc::external) 
register a line search type `SNESLineSearchType`

Logically Collective, No Fortran Support

Input Parameters:
- `sname`    - name of the `SNESLineSearchType()`
- `function` - the creation function for that type

Calling sequence of `function`:
- `ls` - the line search context

Level: advanced

-seealso: [](ch_snes), `SNES`, `SNESLineSearch`, `SNESLineSearchType`, `SNESLineSearchSetType()`

# External Links
$(_doc_external("Snes/SNESLineSearchRegister"))
"""
function SNESLineSearchRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function SNESLineSearchRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:SNESLineSearchRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	SNESLineSearchShellSetApply(petsclib::PetscLibType,linesearch::SNESLineSearch, func::SNESLineSearchShellApplyFn, ctx::Cvoid) 
Sets the apply function for the `SNESLINESEARCHSHELL` implementation.

Not Collective

Input Parameters:
- `linesearch` - `SNESLineSearch` context
- `func`       - function implementing the linesearch shell, see `SNESLineSearchShellApplyFn` for calling sequence
- `ctx`        - context for func

Usage:
-seealso: [](ch_snes), `SNESLineSearchShellGetApply()`, `SNESLINESEARCHSHELL`, `SNESLineSearchType`, `SNESLineSearch`,
`SNESLineSearchShellApplyFn`

# External Links
$(_doc_external("Snes/SNESLineSearchShellSetApply"))
"""
function SNESLineSearchShellSetApply(petsclib::PetscLibType, linesearch::SNESLineSearch, func::SNESLineSearchShellApplyFn, ctx::Cvoid) end

@for_petsc function SNESLineSearchShellSetApply(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, func::SNESLineSearchShellApplyFn, ctx::Cvoid )

    @chk ccall(
               (:SNESLineSearchShellSetApply, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, Ptr{SNESLineSearchShellApplyFn}, Ptr{Cvoid}),
               linesearch, func, ctx,
              )


	return nothing
end 

"""
	SNESLineSearchShellGetApply(petsclib::PetscLibType,linesearch::SNESLineSearch, func::SNESLineSearchShellApplyFn, ctx::Cvoid) 
Gets the apply function and context for the `SNESLINESEARCHSHELL`

Not Collective

Input Parameter:
- `linesearch` - the line search object

Output Parameters:
- `func` - the user function; can be `NULL` if it is not needed, see `SNESLineSearchShellApplyFn` for calling sequence
- `ctx`  - the user function context; can be `NULL` if it is not needed

Level: advanced

-seealso: [](ch_snes), `SNESLineSearchShellSetApply()`, `SNESLINESEARCHSHELL`, `SNESLineSearchType`, `SNESLineSearch`,
`SNESLineSearchShellApplyFn`

# External Links
$(_doc_external("Snes/SNESLineSearchShellGetApply"))
"""
function SNESLineSearchShellGetApply(petsclib::PetscLibType, linesearch::SNESLineSearch, func::SNESLineSearchShellApplyFn, ctx::Cvoid) end

@for_petsc function SNESLineSearchShellGetApply(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, func::SNESLineSearchShellApplyFn, ctx::Cvoid )

    @chk ccall(
               (:SNESLineSearchShellGetApply, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, SNESLineSearchShellApplyFn, Cvoid),
               linesearch, func, ctx,
              )


	return nothing
end 

"""
	SNESLineSearchBTSetAlpha(petsclib::PetscLibType,linesearch::SNESLineSearch, alpha::PetscReal) 
Sets the descent parameter, `alpha`, in the `SNESLINESEARCHBT` `SNESLineSearch` variant.

Input Parameters:
- `linesearch` - linesearch context
- `alpha`      - The descent parameter

Level: intermediate

-seealso: [](ch_snes), `SNESLineSearch`, `SNESLineSearchSetLambda()`, `SNESLineSearchGetTolerances()`, `SNESLINESEARCHBT`, `SNESLineSearchBTGetAlpha()`

# External Links
$(_doc_external("Snes/SNESLineSearchBTSetAlpha"))
"""
function SNESLineSearchBTSetAlpha(petsclib::PetscLibType, linesearch::SNESLineSearch, alpha::PetscReal) end

@for_petsc function SNESLineSearchBTSetAlpha(petsclib::$UnionPetscLib, linesearch::SNESLineSearch, alpha::$PetscReal )

    @chk ccall(
               (:SNESLineSearchBTSetAlpha, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, $PetscReal),
               linesearch, alpha,
              )


	return nothing
end 

"""
	alpha::PetscReal = SNESLineSearchBTGetAlpha(petsclib::PetscLibType,linesearch::SNESLineSearch) 
Gets the descent parameter, `alpha`, in the `SNESLINESEARCHBT` variant that was set with `SNESLineSearchBTSetAlpha()`

Input Parameter:
- `linesearch` - linesearch context

Output Parameter:
- `alpha` - The descent parameter

Level: intermediate

-seealso: [](ch_snes), `SNESLineSearch`, `SNESLineSearchGetLambda()`, `SNESLineSearchGetTolerances()`, `SNESLINESEARCHBT`, `SNESLineSearchBTSetAlpha()`

# External Links
$(_doc_external("Snes/SNESLineSearchBTGetAlpha"))
"""
function SNESLineSearchBTGetAlpha(petsclib::PetscLibType, linesearch::SNESLineSearch) end

@for_petsc function SNESLineSearchBTGetAlpha(petsclib::$UnionPetscLib, linesearch::SNESLineSearch )
	alpha_ = Ref{$PetscReal}()

    @chk ccall(
               (:SNESLineSearchBTGetAlpha, $petsc_library),
               PetscErrorCode,
               (SNESLineSearch, Ptr{$PetscReal}),
               linesearch, alpha_,
              )

	alpha = alpha_[]

	return alpha
end 

