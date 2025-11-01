# autodefined type arguments for class ------
mutable struct _n_TaoMonitorDrawCtx end
const TaoMonitorDrawCtx = Ptr{_n_TaoMonitorDrawCtx}
# -------------------------------------------------------
"""
	ctx::TaoMonitorDrawCtx = TaoMonitorDrawCtxCreate(petsclib::PetscLibType,comm::MPI_Comm, host::String, label::String, x::Cint, y::Cint, m::Cint, n::Cint, howoften::PetscInt) 
Creates the monitor context for `TaoMonitorSolutionDraw()`

Collective

Input Parameters:
- `comm`     - the communicator to share the context
- `host`     - the name of the X Windows host that will display the monitor
- `label`    - the label to put at the top of the display window
- `x`        - the horizontal coordinate of the lower left corner of the window to open
- `y`        - the vertical coordinate of the lower left corner of the window to open
- `m`        - the width of the window
- `n`        - the height of the window
- `howoften` - how many `Tao` iterations between displaying the monitor information

Output Parameter:
- `ctx` - the monitor context

Options Database Keys:
- `-tao_monitor_solution_draw` - use `TaoMonitorSolutionDraw()` to monitor the solution
- `-tao_draw_solution_initial` - show initial guess as well as current solution

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoMonitorSet()`, `TaoMonitorDefault()`, `VecView()`, `TaoMonitorDrawCtx()`

# External Links
$(_doc_external("Tao/TaoMonitorDrawCtxCreate"))
"""
function TaoMonitorDrawCtxCreate(petsclib::PetscLibType, comm::MPI_Comm, host::String, label::String, x::Cint, y::Cint, m::Cint, n::Cint, howoften::PetscInt) end

@for_petsc function TaoMonitorDrawCtxCreate(petsclib::$UnionPetscLib, comm::MPI_Comm, host::String, label::String, x::Cint, y::Cint, m::Cint, n::Cint, howoften::$PetscInt )
	ctx_ = Ref{TaoMonitorDrawCtx}()

    @chk ccall(
               (:TaoMonitorDrawCtxCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Cint, Cint, Cint, Cint, $PetscInt, Ptr{TaoMonitorDrawCtx}),
               comm, host, label, x, y, m, n, howoften, ctx_,
              )

	ctx = ctx_[]

	return ctx
end 

"""
	TaoMonitorDrawCtxDestroy(petsclib::PetscLibType,ictx::TaoMonitorDrawCtx) 
Destroys the monitor context for `TaoMonitorSolutionDraw()`

Collective

Input Parameter:
- `ictx` - the monitor context

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoMonitorSet()`, `TaoMonitorDefault()`, `VecView()`, `TaoMonitorSolutionDraw()`

# External Links
$(_doc_external("Tao/TaoMonitorDrawCtxDestroy"))
"""
function TaoMonitorDrawCtxDestroy(petsclib::PetscLibType, ictx::TaoMonitorDrawCtx) end

@for_petsc function TaoMonitorDrawCtxDestroy(petsclib::$UnionPetscLib, ictx::TaoMonitorDrawCtx )

    @chk ccall(
               (:TaoMonitorDrawCtxDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{TaoMonitorDrawCtx},),
               ictx,
              )


	return nothing
end 

"""
	TaoLineSearchViewFromOptions(petsclib::PetscLibType,A::TaoLineSearch, obj::PetscObject, name::String) 
View a `TaoLineSearch` object based on values in the options database

Collective

Input Parameters:
- `A`    - the `Tao` context
- `obj`  - Optional object
- `name` - command line option

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoLineSearch`, `TaoLineSearchView()`, `PetscObjectViewFromOptions()`, `TaoLineSearchCreate()`

# External Links
$(_doc_external("Tao/TaoLineSearchViewFromOptions"))
"""
function TaoLineSearchViewFromOptions(petsclib::PetscLibType, A::TaoLineSearch, obj::PetscObject, name::String) end

@for_petsc function TaoLineSearchViewFromOptions(petsclib::$UnionPetscLib, A::TaoLineSearch, obj::PetscObject, name::String )

    @chk ccall(
               (:TaoLineSearchViewFromOptions, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

"""
	TaoLineSearchView(petsclib::PetscLibType,ls::TaoLineSearch, viewer::PetscViewer) 
Prints information about the `TaoLineSearch`

Collective

Input Parameters:
- `ls`     - the `TaoLineSearch` context
- `viewer` - visualization context

Options Database Key:
- `-tao_ls_view` - Calls `TaoLineSearchView()` at the end of each line search

Level: beginner

-seealso: [](ch_tao), `Tao`, `TaoLineSearch`, `PetscViewerASCIIOpen()`, `TaoLineSearchViewFromOptions()`

# External Links
$(_doc_external("Tao/TaoLineSearchView"))
"""
function TaoLineSearchView(petsclib::PetscLibType, ls::TaoLineSearch, viewer::PetscViewer) end

@for_petsc function TaoLineSearchView(petsclib::$UnionPetscLib, ls::TaoLineSearch, viewer::PetscViewer )

    @chk ccall(
               (:TaoLineSearchView, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, PetscViewer),
               ls, viewer,
              )


	return nothing
end 

"""
	newls::TaoLineSearch = TaoLineSearchCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Creates a `TaoLineSearch` object.  Algorithms in `Tao` that use
line-searches will automatically create one so this all is rarely needed

Collective

Input Parameter:
- `comm` - MPI communicator

Output Parameter:
- `newls` - the new `TaoLineSearch` context

Options Database Key:
- `-tao_ls_type` - select which method `Tao` should use

Level: developer

-seealso: [](ch_tao), `Tao`, `TaoLineSearch`, `TaoLineSearchType`, `TaoLineSearchSetType()`, `TaoLineSearchApply()`, `TaoLineSearchDestroy()`

# External Links
$(_doc_external("Tao/TaoLineSearchCreate"))
"""
function TaoLineSearchCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function TaoLineSearchCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	newls_ = Ref{TaoLineSearch}()

    @chk ccall(
               (:TaoLineSearchCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{TaoLineSearch}),
               comm, newls_,
              )

	newls = newls_[]

	return newls
end 

"""
	TaoLineSearchSetUp(petsclib::PetscLibType,ls::TaoLineSearch) 
Sets up the internal data structures for the later use
of a `TaoLineSearch`

Collective

Input Parameter:
- `ls` - the `TaoLineSearch` context

Level: developer

-seealso: [](ch_tao), `Tao`, `TaoLineSearch`, `TaoLineSearchCreate()`, `TaoLineSearchApply()`

# External Links
$(_doc_external("Tao/TaoLineSearchSetUp"))
"""
function TaoLineSearchSetUp(petsclib::PetscLibType, ls::TaoLineSearch) end

@for_petsc function TaoLineSearchSetUp(petsclib::$UnionPetscLib, ls::TaoLineSearch )

    @chk ccall(
               (:TaoLineSearchSetUp, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch,),
               ls,
              )


	return nothing
end 

"""
	TaoLineSearchReset(petsclib::PetscLibType,ls::TaoLineSearch) 
Some line searches may carry state information
from one `TaoLineSearchApply()` to the next.  This function resets this
state information.

Collective

Input Parameter:
- `ls` - the `TaoLineSearch` context

Level: developer

-seealso: [](ch_tao), `Tao`, `TaoLineSearch`, `TaoLineSearchCreate()`, `TaoLineSearchApply()`

# External Links
$(_doc_external("Tao/TaoLineSearchReset"))
"""
function TaoLineSearchReset(petsclib::PetscLibType, ls::TaoLineSearch) end

@for_petsc function TaoLineSearchReset(petsclib::$UnionPetscLib, ls::TaoLineSearch )

    @chk ccall(
               (:TaoLineSearchReset, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch,),
               ls,
              )


	return nothing
end 

"""
	TaoLineSearchDestroy(petsclib::PetscLibType,ls::TaoLineSearch) 
Destroys the `TaoLineSearch` context that was created with
`TaoLineSearchCreate()`

Collective

Input Parameter:
- `ls` - the `TaoLineSearch` context

Level: developer

-seealso: `TaoLineSearchCreate()`, `TaoLineSearchSolve()`

# External Links
$(_doc_external("Tao/TaoLineSearchDestroy"))
"""
function TaoLineSearchDestroy(petsclib::PetscLibType, ls::TaoLineSearch) end

@for_petsc function TaoLineSearchDestroy(petsclib::$UnionPetscLib, ls::TaoLineSearch )

    @chk ccall(
               (:TaoLineSearchDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{TaoLineSearch},),
               ls,
              )


	return nothing
end 

"""
	f::PetscReal,steplength::PetscReal = TaoLineSearchApply(petsclib::PetscLibType,ls::TaoLineSearch, x::PetscVec, g::PetscVec, s::PetscVec, reason::TaoLineSearchConvergedReason) 
Performs a line
Criteria for acceptable step length depends on the line-search algorithm chosen

Collective

Input Parameters:
- `ls` - the `TaoLineSearch` context
- `s`  - search direction

Output Parameters:
- `x`          - On input the current solution, on output `x` contains the new solution determined by the line search
- `f`          - On input the objective function value at current solution, on output contains the objective function value at new solution
- `g`          - On input the gradient evaluated at `x`, on output contains the gradient at new solution
- `steplength` - scalar multiplier of s used ( x = x0 + steplength * x)
- `reason`     - `TaoLineSearchConvergedReason` reason why the line-search stopped

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoLineSearchConvergedReason`, `TaoLineSearch`, `TaoLineSearchCreate()`, `TaoLineSearchSetType()`,
`TaoLineSearchSetInitialStepLength()`, `TaoAddLineSearchCounts()`

# External Links
$(_doc_external("Tao/TaoLineSearchApply"))
"""
function TaoLineSearchApply(petsclib::PetscLibType, ls::TaoLineSearch, x::PetscVec, g::PetscVec, s::PetscVec, reason::TaoLineSearchConvergedReason) end

@for_petsc function TaoLineSearchApply(petsclib::$UnionPetscLib, ls::TaoLineSearch, x::PetscVec, g::PetscVec, s::PetscVec, reason::TaoLineSearchConvergedReason )
	f_ = Ref{$PetscReal}()
	steplength_ = Ref{$PetscReal}()

    @chk ccall(
               (:TaoLineSearchApply, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, CVec, Ptr{$PetscReal}, CVec, CVec, Ptr{$PetscReal}, Ptr{TaoLineSearchConvergedReason}),
               ls, x, f_, g, s, steplength_, reason,
              )

	f = f_[]
	steplength = steplength_[]

	return f,steplength
end 

"""
	TaoLineSearchSetType(petsclib::PetscLibType,ls::TaoLineSearch, type::TaoLineSearchType) 
Sets the algorithm used in a line search

Collective

Input Parameters:
- `ls`   - the `TaoLineSearch` context
- `type` - the `TaoLineSearchType` selection

Options Database Key:
- `-tao_ls_type <type>` - select which method Tao should use at runtime

Level: beginner

-seealso: [](ch_tao), `Tao`, `TaoLineSearch`, `TaoLineSearchType`, `TaoLineSearchCreate()`, `TaoLineSearchGetType()`,
`TaoLineSearchApply()`

# External Links
$(_doc_external("Tao/TaoLineSearchSetType"))
"""
function TaoLineSearchSetType(petsclib::PetscLibType, ls::TaoLineSearch, type::TaoLineSearchType) end

@for_petsc function TaoLineSearchSetType(petsclib::$UnionPetscLib, ls::TaoLineSearch, type::TaoLineSearchType )

    @chk ccall(
               (:TaoLineSearchSetType, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, TaoLineSearchType),
               ls, type,
              )


	return nothing
end 

"""
	TaoLineSearchMonitor(petsclib::PetscLibType,ls::TaoLineSearch, its::PetscInt, f::PetscReal, step::PetscReal) 
Monitor the line search steps. This routine will output the
iteration number, step length, and function value before calling the implementation
specific monitor.

Input Parameters:
- `ls`   - the `TaoLineSearch` context
- `its`  - the current iterate number (>=0)
- `f`    - the current objective function value
- `step` - the step length

Options Database Key:
- `-tao_ls_monitor` - Use the default monitor, which prints statistics to standard output

Level: developer

-seealso: `TaoLineSearch`

# External Links
$(_doc_external("Tao/TaoLineSearchMonitor"))
"""
function TaoLineSearchMonitor(petsclib::PetscLibType, ls::TaoLineSearch, its::PetscInt, f::PetscReal, step::PetscReal) end

@for_petsc function TaoLineSearchMonitor(petsclib::$UnionPetscLib, ls::TaoLineSearch, its::$PetscInt, f::$PetscReal, step::$PetscReal )

    @chk ccall(
               (:TaoLineSearchMonitor, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, $PetscInt, $PetscReal, $PetscReal),
               ls, its, f, step,
              )


	return nothing
end 

"""
	TaoLineSearchSetFromOptions(petsclib::PetscLibType,ls::TaoLineSearch) 
Sets various `TaoLineSearch` parameters from user
options.

Collective

Input Parameter:
- `ls` - the `TaoLineSearch` context

Options Database Keys:
- `-tao_ls_type <type>`     - The algorithm that `TaoLineSearch` uses (more-thuente, gpcg, unit)
- `-tao_ls_ftol <tol>`      - tolerance for sufficient decrease
- `-tao_ls_gtol <tol>`      - tolerance for curvature condition
- `-tao_ls_rtol <tol>`      - relative tolerance for acceptable step
- `-tao_ls_stepinit <step>` - initial steplength allowed
- `-tao_ls_stepmin <step>`  - minimum steplength allowed
- `-tao_ls_stepmax <step>`  - maximum steplength allowed
- `-tao_ls_max_funcs <n>`   - maximum number of function evaluations allowed
- `-tao_ls_view`            - display line-search results to standard output

Level: beginner

-seealso: `TaoLineSearch`

# External Links
$(_doc_external("Tao/TaoLineSearchSetFromOptions"))
"""
function TaoLineSearchSetFromOptions(petsclib::PetscLibType, ls::TaoLineSearch) end

@for_petsc function TaoLineSearchSetFromOptions(petsclib::$UnionPetscLib, ls::TaoLineSearch )

    @chk ccall(
               (:TaoLineSearchSetFromOptions, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch,),
               ls,
              )


	return nothing
end 

"""
	type::TaoLineSearchType = TaoLineSearchGetType(petsclib::PetscLibType,ls::TaoLineSearch) 
Gets the current line search algorithm

Not Collective

Input Parameter:
- `ls` - the `TaoLineSearch` context

Output Parameter:
- `type` - the line search algorithm in effect

Level: developer

-seealso: `TaoLineSearch`

# External Links
$(_doc_external("Tao/TaoLineSearchGetType"))
"""
function TaoLineSearchGetType(petsclib::PetscLibType, ls::TaoLineSearch) end

@for_petsc function TaoLineSearchGetType(petsclib::$UnionPetscLib, ls::TaoLineSearch )
	type_ = Ref{TaoLineSearchType}()

    @chk ccall(
               (:TaoLineSearchGetType, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, Ptr{TaoLineSearchType}),
               ls, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	nfeval::PetscInt,ngeval::PetscInt,nfgeval::PetscInt = TaoLineSearchGetNumberFunctionEvaluations(petsclib::PetscLibType,ls::TaoLineSearch) 
Gets the number of function and gradient evaluation
routines used by the line search in last application (not cumulative).

Not Collective

Input Parameter:
- `ls` - the `TaoLineSearch` context

Output Parameters:
- `nfeval`  - number of function evaluations
- `ngeval`  - number of gradient evaluations
- `nfgeval` - number of function/gradient evaluations

Level: intermediate

-seealso: `TaoLineSearch`

# External Links
$(_doc_external("Tao/TaoLineSearchGetNumberFunctionEvaluations"))
"""
function TaoLineSearchGetNumberFunctionEvaluations(petsclib::PetscLibType, ls::TaoLineSearch) end

@for_petsc function TaoLineSearchGetNumberFunctionEvaluations(petsclib::$UnionPetscLib, ls::TaoLineSearch )
	nfeval_ = Ref{$PetscInt}()
	ngeval_ = Ref{$PetscInt}()
	nfgeval_ = Ref{$PetscInt}()

    @chk ccall(
               (:TaoLineSearchGetNumberFunctionEvaluations, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
               ls, nfeval_, ngeval_, nfgeval_,
              )

	nfeval = nfeval_[]
	ngeval = ngeval_[]
	nfgeval = nfgeval_[]

	return nfeval,ngeval,nfgeval
end 

"""
	flg::PetscBool = TaoLineSearchIsUsingTaoRoutines(petsclib::PetscLibType,ls::TaoLineSearch) 
Checks whether the line search is using
the standard `Tao` evaluation routines.

Not Collective

Input Parameter:
- `ls` - the `TaoLineSearch` context

Output Parameter:
- `flg` - `PETSC_TRUE` if the line search is using `Tao` evaluation routines,
otherwise `PETSC_FALSE`

Level: developer

-seealso: `TaoLineSearch`

# External Links
$(_doc_external("Tao/TaoLineSearchIsUsingTaoRoutines"))
"""
function TaoLineSearchIsUsingTaoRoutines(petsclib::PetscLibType, ls::TaoLineSearch) end

@for_petsc function TaoLineSearchIsUsingTaoRoutines(petsclib::$UnionPetscLib, ls::TaoLineSearch )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:TaoLineSearchIsUsingTaoRoutines, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, Ptr{PetscBool}),
               ls, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	TaoLineSearchSetObjectiveRoutine(petsclib::PetscLibType,ls::TaoLineSearch, func::external, ctx::Cvoid) 
Sets the function evaluation routine for the line search

Logically Collective

Input Parameters:
- `ls`   - the `TaoLineSearch` context
- `func` - the objective function evaluation routine
- `ctx`  - the (optional) user-defined context for private data

Calling sequence of `func`:
- `ls`  - the line search context
- `x`   - input vector
- `f`   - function value
- `ctx` - (optional) user-defined context

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoLineSearch`, `TaoLineSearchCreate()`, `TaoLineSearchSetGradientRoutine()`, `TaoLineSearchSetObjectiveAndGradientRoutine()`, `TaoLineSearchUseTaoRoutines()`

# External Links
$(_doc_external("Tao/TaoLineSearchSetObjectiveRoutine"))
"""
function TaoLineSearchSetObjectiveRoutine(petsclib::PetscLibType, ls::TaoLineSearch, func::external, ctx::Cvoid) end

@for_petsc function TaoLineSearchSetObjectiveRoutine(petsclib::$UnionPetscLib, ls::TaoLineSearch, func::external, ctx::Cvoid )

    @chk ccall(
               (:TaoLineSearchSetObjectiveRoutine, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, external, Ptr{Cvoid}),
               ls, func, ctx,
              )


	return nothing
end 

"""
	TaoLineSearchSetGradientRoutine(petsclib::PetscLibType,ls::TaoLineSearch, func::external, ctx::Cvoid) 
Sets the gradient evaluation routine for the line search

Logically Collective

Input Parameters:
- `ls`   - the `TaoLineSearch` context
- `func` - the gradient evaluation routine
- `ctx`  - the (optional) user-defined context for private data

Calling sequence of `func`:
- `ls`  - the linesearch object
- `x`   - input vector
- `g`   - gradient vector
- `ctx` - (optional) user-defined context

Level: beginner

-seealso: [](ch_tao), `Tao`, `TaoLineSearch`, `TaoLineSearchCreate()`, `TaoLineSearchSetObjectiveRoutine()`, `TaoLineSearchSetObjectiveAndGradientRoutine()`, `TaoLineSearchUseTaoRoutines()`

# External Links
$(_doc_external("Tao/TaoLineSearchSetGradientRoutine"))
"""
function TaoLineSearchSetGradientRoutine(petsclib::PetscLibType, ls::TaoLineSearch, func::external, ctx::Cvoid) end

@for_petsc function TaoLineSearchSetGradientRoutine(petsclib::$UnionPetscLib, ls::TaoLineSearch, func::external, ctx::Cvoid )

    @chk ccall(
               (:TaoLineSearchSetGradientRoutine, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, external, Ptr{Cvoid}),
               ls, func, ctx,
              )


	return nothing
end 

"""
	TaoLineSearchSetObjectiveAndGradientRoutine(petsclib::PetscLibType,ls::TaoLineSearch, func::external, ctx::Cvoid) 
Sets the objective/gradient evaluation routine for the line search

Logically Collective

Input Parameters:
- `ls`   - the `TaoLineSearch` context
- `func` - the objective and gradient evaluation routine
- `ctx`  - the (optional) user-defined context for private data

Calling sequence of `func`:
- `ls`  - the linesearch object
- `x`   - input vector
- `f`   - function value
- `g`   - gradient vector
- `ctx` - (optional) user-defined context

Level: beginner

-seealso: [](ch_tao), `Tao`, `TaoLineSearch`, `TaoLineSearchCreate()`, `TaoLineSearchSetObjectiveRoutine()`, `TaoLineSearchSetGradientRoutine()`, `TaoLineSearchUseTaoRoutines()`

# External Links
$(_doc_external("Tao/TaoLineSearchSetObjectiveAndGradientRoutine"))
"""
function TaoLineSearchSetObjectiveAndGradientRoutine(petsclib::PetscLibType, ls::TaoLineSearch, func::external, ctx::Cvoid) end

@for_petsc function TaoLineSearchSetObjectiveAndGradientRoutine(petsclib::$UnionPetscLib, ls::TaoLineSearch, func::external, ctx::Cvoid )

    @chk ccall(
               (:TaoLineSearchSetObjectiveAndGradientRoutine, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, external, Ptr{Cvoid}),
               ls, func, ctx,
              )


	return nothing
end 

"""
	TaoLineSearchSetObjectiveAndGTSRoutine(petsclib::PetscLibType,ls::TaoLineSearch, func::external, ctx::Cvoid) 
Sets the objective and
(gradient'*stepdirection) evaluation routine for the line search.

Logically Collective

Input Parameters:
- `ls`   - the `TaoLineSearch` context
- `func` - the objective and gradient evaluation routine
- `ctx`  - the (optional) user-defined context for private data

Calling sequence of `func`:
- `ls`  - the linesearch context
- `x`   - input vector
- `s`   - step direction
- `f`   - function value
- `gts` - inner product of gradient and step direction vectors
- `ctx` - (optional) user-defined context

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoLineSearch`, `TaoLineSearchCreate()`, `TaoLineSearchSetObjective()`, `TaoLineSearchSetGradient()`, `TaoLineSearchUseTaoRoutines()`

# External Links
$(_doc_external("Tao/TaoLineSearchSetObjectiveAndGTSRoutine"))
"""
function TaoLineSearchSetObjectiveAndGTSRoutine(petsclib::PetscLibType, ls::TaoLineSearch, func::external, ctx::Cvoid) end

@for_petsc function TaoLineSearchSetObjectiveAndGTSRoutine(petsclib::$UnionPetscLib, ls::TaoLineSearch, func::external, ctx::Cvoid )

    @chk ccall(
               (:TaoLineSearchSetObjectiveAndGTSRoutine, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, external, Ptr{Cvoid}),
               ls, func, ctx,
              )


	return nothing
end 

"""
	TaoLineSearchUseTaoRoutines(petsclib::PetscLibType,ls::TaoLineSearch, ts::Tao) 
Informs the `TaoLineSearch` to use the
objective and gradient evaluation routines from the given `Tao` object. The default.

Logically Collective

Input Parameters:
- `ls` - the `TaoLineSearch` context
- `ts` - the `Tao` context with defined objective/gradient evaluation routines

Level: developer

-seealso: [](ch_tao), `Tao`, `TaoLineSearch`, `TaoLineSearchCreate()`

# External Links
$(_doc_external("Tao/TaoLineSearchUseTaoRoutines"))
"""
function TaoLineSearchUseTaoRoutines(petsclib::PetscLibType, ls::TaoLineSearch, ts::Tao) end

@for_petsc function TaoLineSearchUseTaoRoutines(petsclib::$UnionPetscLib, ls::TaoLineSearch, ts::Tao )

    @chk ccall(
               (:TaoLineSearchUseTaoRoutines, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, CTao),
               ls, ts,
              )


	return nothing
end 

"""
	f::PetscReal = TaoLineSearchComputeObjective(petsclib::PetscLibType,ls::TaoLineSearch, x::PetscVec) 
Computes the objective function value at a given point

Collective

Input Parameters:
- `ls` - the `TaoLineSearch` context
- `x`  - input vector

Output Parameter:
- `f` - Objective value at `x`

Level: developer

-seealso: [](ch_tao), `Tao`, `TaoLineSearch`, `TaoLineSearchComputeGradient()`, `TaoLineSearchComputeObjectiveAndGradient()`, `TaoLineSearchSetObjectiveRoutine()`

# External Links
$(_doc_external("Tao/TaoLineSearchComputeObjective"))
"""
function TaoLineSearchComputeObjective(petsclib::PetscLibType, ls::TaoLineSearch, x::PetscVec) end

@for_petsc function TaoLineSearchComputeObjective(petsclib::$UnionPetscLib, ls::TaoLineSearch, x::PetscVec )
	f_ = Ref{$PetscReal}()

    @chk ccall(
               (:TaoLineSearchComputeObjective, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, CVec, Ptr{$PetscReal}),
               ls, x, f_,
              )

	f = f_[]

	return f
end 

"""
	f::PetscReal = TaoLineSearchComputeObjectiveAndGradient(petsclib::PetscLibType,ls::TaoLineSearch, x::PetscVec, g::PetscVec) 
Computes the objective function value at a given point

Collective

Input Parameters:
- `ls` - the `TaoLineSearch` context
- `x`  - input vector

Output Parameters:
- `f` - Objective value at `x`
- `g` - Gradient vector at `x`

Level: developer

-seealso: [](ch_tao), `Tao`, `TaoLineSearch`, `TaoLineSearchComputeGradient()`, `TaoLineSearchSetObjectiveRoutine()`

# External Links
$(_doc_external("Tao/TaoLineSearchComputeObjectiveAndGradient"))
"""
function TaoLineSearchComputeObjectiveAndGradient(petsclib::PetscLibType, ls::TaoLineSearch, x::PetscVec, g::PetscVec) end

@for_petsc function TaoLineSearchComputeObjectiveAndGradient(petsclib::$UnionPetscLib, ls::TaoLineSearch, x::PetscVec, g::PetscVec )
	f_ = Ref{$PetscReal}()

    @chk ccall(
               (:TaoLineSearchComputeObjectiveAndGradient, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, CVec, Ptr{$PetscReal}, CVec),
               ls, x, f_, g,
              )

	f = f_[]

	return f
end 

"""
	TaoLineSearchComputeGradient(petsclib::PetscLibType,ls::TaoLineSearch, x::PetscVec, g::PetscVec) 
Computes the gradient of the objective function

Collective

Input Parameters:
- `ls` - the `TaoLineSearch` context
- `x`  - input vector

Output Parameter:
- `g` - gradient vector

Level: developer

-seealso: [](ch_tao), `Tao`, `TaoLineSearch`, `TaoLineSearchComputeObjective()`, `TaoLineSearchComputeObjectiveAndGradient()`, `TaoLineSearchSetGradient()`

# External Links
$(_doc_external("Tao/TaoLineSearchComputeGradient"))
"""
function TaoLineSearchComputeGradient(petsclib::PetscLibType, ls::TaoLineSearch, x::PetscVec, g::PetscVec) end

@for_petsc function TaoLineSearchComputeGradient(petsclib::$UnionPetscLib, ls::TaoLineSearch, x::PetscVec, g::PetscVec )

    @chk ccall(
               (:TaoLineSearchComputeGradient, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, CVec, CVec),
               ls, x, g,
              )


	return nothing
end 

"""
	f::PetscReal,gts::PetscReal = TaoLineSearchComputeObjectiveAndGTS(petsclib::PetscLibType,ls::TaoLineSearch, x::PetscVec) 
Computes the objective function value and inner product of gradient and
step direction at a given point

Collective

Input Parameters:
- `ls` - the `TaoLineSearch` context
- `x`  - input vector

Output Parameters:
- `f`   - Objective value at `x`
- `gts` - inner product of gradient and step direction at `x`

Level: developer

-seealso: [](ch_tao), `Tao`, `TaoLineSearch`, `TaoLineSearchComputeGradient()`, `TaoLineSearchComputeObjectiveAndGradient()`, `TaoLineSearchSetObjectiveRoutine()`

# External Links
$(_doc_external("Tao/TaoLineSearchComputeObjectiveAndGTS"))
"""
function TaoLineSearchComputeObjectiveAndGTS(petsclib::PetscLibType, ls::TaoLineSearch, x::PetscVec) end

@for_petsc function TaoLineSearchComputeObjectiveAndGTS(petsclib::$UnionPetscLib, ls::TaoLineSearch, x::PetscVec )
	f_ = Ref{$PetscReal}()
	gts_ = Ref{$PetscReal}()

    @chk ccall(
               (:TaoLineSearchComputeObjectiveAndGTS, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, CVec, Ptr{$PetscReal}, Ptr{$PetscReal}),
               ls, x, f_, gts_,
              )

	f = f_[]
	gts = gts_[]

	return f,gts
end 

"""
	f::PetscReal,steplength::PetscReal = TaoLineSearchGetSolution(petsclib::PetscLibType,ls::TaoLineSearch, x::PetscVec, g::PetscVec, reason::TaoLineSearchConvergedReason) 
Returns the solution to the line search

Collective

Input Parameter:
- `ls` - the `TaoLineSearch` context

Output Parameters:
- `x`          - the new solution
- `f`          - the objective function value at `x`
- `g`          - the gradient at `x`
- `steplength` - the multiple of the step direction taken by the line search
- `reason`     - the reason why the line search terminated

Level: developer

-seealso: `TaoLineSearchGetStartingVector()`, `TaoLineSearchGetStepDirection()`

# External Links
$(_doc_external("Tao/TaoLineSearchGetSolution"))
"""
function TaoLineSearchGetSolution(petsclib::PetscLibType, ls::TaoLineSearch, x::PetscVec, g::PetscVec, reason::TaoLineSearchConvergedReason) end

@for_petsc function TaoLineSearchGetSolution(petsclib::$UnionPetscLib, ls::TaoLineSearch, x::PetscVec, g::PetscVec, reason::TaoLineSearchConvergedReason )
	f_ = Ref{$PetscReal}()
	steplength_ = Ref{$PetscReal}()

    @chk ccall(
               (:TaoLineSearchGetSolution, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, CVec, Ptr{$PetscReal}, CVec, Ptr{$PetscReal}, Ptr{TaoLineSearchConvergedReason}),
               ls, x, f_, g, steplength_, reason,
              )

	f = f_[]
	steplength = steplength_[]

	return f,steplength
end 

"""
	TaoLineSearchGetStartingVector(petsclib::PetscLibType,ls::TaoLineSearch, x::PetscVec) 
Gets a the initial point of the line
search.

Not Collective

Input Parameter:
- `ls` - the `TaoLineSearch` context

Output Parameter:
- `x` - The initial point of the line search

Level: advanced

-seealso: `TaoLineSearchGetSolution()`, `TaoLineSearchGetStepDirection()`

# External Links
$(_doc_external("Tao/TaoLineSearchGetStartingVector"))
"""
function TaoLineSearchGetStartingVector(petsclib::PetscLibType, ls::TaoLineSearch, x::PetscVec) end

@for_petsc function TaoLineSearchGetStartingVector(petsclib::$UnionPetscLib, ls::TaoLineSearch, x::PetscVec )
	x_ = Ref(x.ptr)

    @chk ccall(
               (:TaoLineSearchGetStartingVector, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, Ptr{CVec}),
               ls, x_,
              )

	x.ptr = C_NULL

	return nothing
end 

"""
	TaoLineSearchGetStepDirection(petsclib::PetscLibType,ls::TaoLineSearch, s::PetscVec) 
Gets the step direction of the line
search.

Not Collective

Input Parameter:
- `ls` - the `TaoLineSearch` context

Output Parameter:
- `s` - the step direction of the line search

Level: advanced

-seealso: `TaoLineSearchGetSolution()`, `TaoLineSearchGetStartingVector()`

# External Links
$(_doc_external("Tao/TaoLineSearchGetStepDirection"))
"""
function TaoLineSearchGetStepDirection(petsclib::PetscLibType, ls::TaoLineSearch, s::PetscVec) end

@for_petsc function TaoLineSearchGetStepDirection(petsclib::$UnionPetscLib, ls::TaoLineSearch, s::PetscVec )
	s_ = Ref(s.ptr)

    @chk ccall(
               (:TaoLineSearchGetStepDirection, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, Ptr{CVec}),
               ls, s_,
              )

	s.ptr = C_NULL

	return nothing
end 

"""
	f_fullstep::PetscReal = TaoLineSearchGetFullStepObjective(petsclib::PetscLibType,ls::TaoLineSearch) 
Returns the objective function value at the full step.  Useful for some minimization algorithms.

Not Collective

Input Parameter:
- `ls` - the `TaoLineSearch` context

Output Parameter:
- `f_fullstep` - the objective value at the full step length

Level: developer

-seealso: `TaoLineSearchGetSolution()`, `TaoLineSearchGetStartingVector()`, `TaoLineSearchGetStepDirection()`

# External Links
$(_doc_external("Tao/TaoLineSearchGetFullStepObjective"))
"""
function TaoLineSearchGetFullStepObjective(petsclib::PetscLibType, ls::TaoLineSearch) end

@for_petsc function TaoLineSearchGetFullStepObjective(petsclib::$UnionPetscLib, ls::TaoLineSearch )
	f_fullstep_ = Ref{$PetscReal}()

    @chk ccall(
               (:TaoLineSearchGetFullStepObjective, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, Ptr{$PetscReal}),
               ls, f_fullstep_,
              )

	f_fullstep = f_fullstep_[]

	return f_fullstep
end 

"""
	TaoLineSearchSetVariableBounds(petsclib::PetscLibType,ls::TaoLineSearch, xl::PetscVec, xu::PetscVec) 
Sets the upper and lower bounds for a bounded line search

Logically Collective

Input Parameters:
- `ls` - the `TaoLineSearch` context
- `xl` - vector of lower bounds
- `xu` - vector of upper bounds

Level: beginner

-seealso: [](ch_tao), `Tao`, `TaoLineSearch`, `TaoSetVariableBounds()`, `TaoLineSearchCreate()`

# External Links
$(_doc_external("Tao/TaoLineSearchSetVariableBounds"))
"""
function TaoLineSearchSetVariableBounds(petsclib::PetscLibType, ls::TaoLineSearch, xl::PetscVec, xu::PetscVec) end

@for_petsc function TaoLineSearchSetVariableBounds(petsclib::$UnionPetscLib, ls::TaoLineSearch, xl::PetscVec, xu::PetscVec )

    @chk ccall(
               (:TaoLineSearchSetVariableBounds, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, CVec, CVec),
               ls, xl, xu,
              )


	return nothing
end 

"""
	TaoLineSearchSetInitialStepLength(petsclib::PetscLibType,ls::TaoLineSearch, s::PetscReal) 
Sets the initial step length of a line
search.  If this value is not set then 1.0 is assumed.

Logically Collective

Input Parameters:
- `ls` - the `TaoLineSearch` context
- `s`  - the initial step size

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoLineSearch`, `TaoLineSearchGetStepLength()`, `TaoLineSearchApply()`

# External Links
$(_doc_external("Tao/TaoLineSearchSetInitialStepLength"))
"""
function TaoLineSearchSetInitialStepLength(petsclib::PetscLibType, ls::TaoLineSearch, s::PetscReal) end

@for_petsc function TaoLineSearchSetInitialStepLength(petsclib::$UnionPetscLib, ls::TaoLineSearch, s::$PetscReal )

    @chk ccall(
               (:TaoLineSearchSetInitialStepLength, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, $PetscReal),
               ls, s,
              )


	return nothing
end 

"""
	s::PetscReal = TaoLineSearchGetStepLength(petsclib::PetscLibType,ls::TaoLineSearch) 
Get the current step length

Not Collective

Input Parameter:
- `ls` - the `TaoLineSearch` context

Output Parameter:
- `s` - the current step length

Level: intermediate

-seealso: [](ch_tao), `Tao`, `TaoLineSearch`, `TaoLineSearchSetInitialStepLength()`, `TaoLineSearchApply()`

# External Links
$(_doc_external("Tao/TaoLineSearchGetStepLength"))
"""
function TaoLineSearchGetStepLength(petsclib::PetscLibType, ls::TaoLineSearch) end

@for_petsc function TaoLineSearchGetStepLength(petsclib::$UnionPetscLib, ls::TaoLineSearch )
	s_ = Ref{$PetscReal}()

    @chk ccall(
               (:TaoLineSearchGetStepLength, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, Ptr{$PetscReal}),
               ls, s_,
              )

	s = s_[]

	return s
end 

"""
	TaoLineSearchRegister(petsclib::PetscLibType,sname::String, func::external) 
Adds a line

Not Collective, No Fortran Support

Input Parameters:
- `sname` - name of a new user-defined solver
- `func`  - routine to Create method context

-seealso: [](ch_tao), `Tao`, `TaoLineSearch`

# External Links
$(_doc_external("Tao/TaoLineSearchRegister"))
"""
function TaoLineSearchRegister(petsclib::PetscLibType, sname::String, func::external) end

@for_petsc function TaoLineSearchRegister(petsclib::$UnionPetscLib, sname::String, func::external )

    @chk ccall(
               (:TaoLineSearchRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, func,
              )


	return nothing
end 

"""
	TaoLineSearchAppendOptionsPrefix(petsclib::PetscLibType,ls::TaoLineSearch, p::String) 
Appends to the prefix used for searching
for all `TaoLineSearch` options in the database.

Collective

Input Parameters:
- `ls` - the `TaoLineSearch` solver context
- `p`  - the prefix string to prepend to all line search requests

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoLineSearch`, `TaoLineSearchSetOptionsPrefix()`, `TaoLineSearchGetOptionsPrefix()`

# External Links
$(_doc_external("Tao/TaoLineSearchAppendOptionsPrefix"))
"""
function TaoLineSearchAppendOptionsPrefix(petsclib::PetscLibType, ls::TaoLineSearch, p::String) end

@for_petsc function TaoLineSearchAppendOptionsPrefix(petsclib::$UnionPetscLib, ls::TaoLineSearch, p::String )

    @chk ccall(
               (:TaoLineSearchAppendOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, Ptr{Cchar}),
               ls, p,
              )


	return nothing
end 

"""
	TaoLineSearchGetOptionsPrefix(petsclib::PetscLibType,ls::TaoLineSearch, p::String) 
Gets the prefix used for searching for all
`TaoLineSearch` options in the database

Not Collective

Input Parameter:
- `ls` - the `TaoLineSearch` context

Output Parameter:
- `p` - pointer to the prefix string used is returned

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoLineSearch`, `TaoLineSearchSetOptionsPrefix()`, `TaoLineSearchAppendOptionsPrefix()`

# External Links
$(_doc_external("Tao/TaoLineSearchGetOptionsPrefix"))
"""
function TaoLineSearchGetOptionsPrefix(petsclib::PetscLibType, ls::TaoLineSearch, p::String) end

@for_petsc function TaoLineSearchGetOptionsPrefix(petsclib::$UnionPetscLib, ls::TaoLineSearch, p::String )
	p_ = Ref(pointer(p))

    @chk ccall(
               (:TaoLineSearchGetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, Ptr{Ptr{Cchar}}),
               ls, p_,
              )


	return nothing
end 

"""
	TaoLineSearchSetOptionsPrefix(petsclib::PetscLibType,ls::TaoLineSearch, p::String) 
Sets the prefix used for searching for all
`TaoLineSearch` options in the database.

Logically Collective

Input Parameters:
- `ls` - the `TaoLineSearch` context
- `p`  - the prefix string to prepend to all `ls` option requests

Level: advanced

-seealso: [](ch_tao), `Tao`, `TaoLineSearch`, `TaoLineSearchAppendOptionsPrefix()`, `TaoLineSearchGetOptionsPrefix()`

# External Links
$(_doc_external("Tao/TaoLineSearchSetOptionsPrefix"))
"""
function TaoLineSearchSetOptionsPrefix(petsclib::PetscLibType, ls::TaoLineSearch, p::String) end

@for_petsc function TaoLineSearchSetOptionsPrefix(petsclib::$UnionPetscLib, ls::TaoLineSearch, p::String )

    @chk ccall(
               (:TaoLineSearchSetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, Ptr{Cchar}),
               ls, p,
              )


	return nothing
end 

"""
	TaoLineSearchFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the `TaoLineSearch` package. It is called from `PetscFinalize()`.

Level: developer

-seealso: `Tao`, `TaoLineSearch`

# External Links
$(_doc_external("Tao/TaoLineSearchFinalizePackage"))
"""
function TaoLineSearchFinalizePackage(petsclib::PetscLibType) end

@for_petsc function TaoLineSearchFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TaoLineSearchFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TaoLineSearchInitializePackage(petsclib::PetscLibType) 
This function registers the line
algorithms in `Tao`.  When using shared or static libraries, this function is called from the
first entry to `TaoCreate()`; when using dynamic, it is called
from PetscDLLibraryRegister_tao()

Level: developer

-seealso: `Tao`, `TaoLineSearch`, `TaoLineSearchCreate()`

# External Links
$(_doc_external("Tao/TaoLineSearchInitializePackage"))
"""
function TaoLineSearchInitializePackage(petsclib::PetscLibType) end

@for_petsc function TaoLineSearchInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TaoLineSearchInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

