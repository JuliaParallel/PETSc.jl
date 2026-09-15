"""
	TSAdaptCandidateAdd(petsclib::PetscLibType,adapt::TSAdapt, name::String, order::PetscInt, stageorder::PetscInt, ccfl::PetscReal, cost::PetscReal, inuse::PetscBool) 
add a candidate scheme for the adaptive controller to select from

Logically Collective; No Fortran Support

Input Parameters:
- `adapt`      - time step adaptivity context, obtained with `TSGetAdapt()` or `TSAdaptCreate()`
- `name`       - name of the candidate scheme to add
- `order`      - order of the candidate scheme
- `stageorder` - stage order of the candidate scheme
- `ccfl`       - stability coefficient relative to explicit Euler, used for CFL constraints
- `cost`       - relative measure of the amount of work required for the candidate scheme
- `inuse`      - indicates that this scheme is the one currently in use, this flag can only be set for one scheme

Level: developer

-seealso: [](ch_ts), `TSAdapt`, `TSAdaptCandidatesClear()`, `TSAdaptChoose()`

# External Links
$(_doc_external("TS/TSAdaptCandidateAdd"))
"""
function TSAdaptCandidateAdd(petsclib::PetscLibType, adapt::TSAdapt, name::String, order::PetscInt, stageorder::PetscInt, ccfl::PetscReal, cost::PetscReal, inuse::PetscBool) end

@for_petsc function TSAdaptCandidateAdd(petsclib::$UnionPetscLib, adapt::TSAdapt, name::String, order::$PetscInt, stageorder::$PetscInt, ccfl::$PetscReal, cost::$PetscReal, inuse::PetscBool )

    @chk ccall(
               (:TSAdaptCandidateAdd, $petsc_library),
               PetscErrorCode,
               (TSAdapt, Ptr{Cchar}, $PetscInt, $PetscInt, $PetscReal, $PetscReal, PetscBool),
               adapt, name, order, stageorder, ccfl, cost, inuse,
              )


	return nothing
end 

"""
	TSAdaptCandidatesClear(petsclib::PetscLibType,adapt::TSAdapt) 
clear any previously set candidate schemes

Logically Collective

Input Parameter:
- `adapt` - adaptive controller

Level: developer

-seealso: [](ch_ts), `TSAdapt`, `TSAdaptCreate()`, `TSAdaptCandidateAdd()`, `TSAdaptChoose()`

# External Links
$(_doc_external("TS/TSAdaptCandidatesClear"))
"""
function TSAdaptCandidatesClear(petsclib::PetscLibType, adapt::TSAdapt) end

@for_petsc function TSAdaptCandidatesClear(petsclib::$UnionPetscLib, adapt::TSAdapt )

    @chk ccall(
               (:TSAdaptCandidatesClear, $petsc_library),
               PetscErrorCode,
               (TSAdapt,),
               adapt,
              )


	return nothing
end 

"""
	n::PetscInt,order::Ptr{PetscInt},stageorder::Ptr{PetscInt},ccfl::Ptr{PetscReal},cost::Ptr{PetscReal} = TSAdaptCandidatesGet(petsclib::PetscLibType,adapt::TSAdapt) 
Get the list of candidate orders of accuracy and cost

Not Collective

Input Parameter:
- `adapt` - time step adaptivity context

Output Parameters:
- `n`          - number of candidate schemes, always at least 1
- `order`      - the order of each candidate scheme
- `stageorder` - the stage order of each candidate scheme
- `ccfl`       - the CFL coefficient of each scheme
- `cost`       - the relative cost of each scheme

Level: developer

-seealso: [](ch_ts), `TSAdapt`, `TSAdaptCandidatesClear()`, `TSAdaptCandidateAdd()`, `TSAdaptChoose()`

# External Links
$(_doc_external("TS/TSAdaptCandidatesGet"))
"""
function TSAdaptCandidatesGet(petsclib::PetscLibType, adapt::TSAdapt) end

@for_petsc function TSAdaptCandidatesGet(petsclib::$UnionPetscLib, adapt::TSAdapt )
	n_ = Ref{$PetscInt}()
	order_ = Ref{Ptr{$PetscInt}}()
	stageorder_ = Ref{Ptr{$PetscInt}}()
	ccfl_ = Ref{Ptr{$PetscReal}}()
	cost_ = Ref{Ptr{$PetscReal}}()

    @chk ccall(
               (:TSAdaptCandidatesGet, $petsc_library),
               PetscErrorCode,
               (TSAdapt, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscReal}}, Ptr{Ptr{$PetscReal}}),
               adapt, n_, order_, stageorder_, ccfl_, cost_,
              )

	n = n_[]
	order = order_[]
	stageorder = stageorder_[]
	ccfl = ccfl_[]
	cost = cost_[]

	return n,order,stageorder,ccfl,cost
end 

"""
	accept::PetscBool = TSAdaptCheckStage(petsclib::PetscLibType,adapt::TSAdapt, ts::AbstractTS, t::PetscReal, Y::AbstractPetscVec) 
checks whether to accept a stage, (e.g. reject and change time step size if nonlinear solve fails or solution vector is infeasible)

Collective

Input Parameters:
- `adapt` - adaptive controller context
- `ts`    - time stepper
- `t`     - Current simulation time
- `Y`     - Current solution vector

Output Parameter:
- `accept` - `PETSC_TRUE` to accept the stage, `PETSC_FALSE` to reject

Level: developer

-seealso: [](ch_ts), `TSAdapt`

# External Links
$(_doc_external("TS/TSAdaptCheckStage"))
"""
function TSAdaptCheckStage(petsclib::PetscLibType, adapt::TSAdapt, ts::AbstractTS, t::PetscReal, Y::AbstractPetscVec) end

@for_petsc function TSAdaptCheckStage(petsclib::$UnionPetscLib, adapt::TSAdapt, ts::AbstractTS, t::$PetscReal, Y::AbstractPetscVec )
	accept_ = Ref{PetscBool}()

    @chk ccall(
               (:TSAdaptCheckStage, $petsc_library),
               PetscErrorCode,
               (TSAdapt, CTS, $PetscReal, CVec, Ptr{PetscBool}),
               adapt, ts, t, Y, accept_,
              )

	accept = accept_[]

	return accept
end 

"""
	next_sc::PetscInt,next_h::PetscReal,accept::PetscBool = TSAdaptChoose(petsclib::PetscLibType,adapt::TSAdapt, ts::AbstractTS, h::PetscReal) 
choose which method and step size to use for the next step

Collective

Input Parameters:
- `adapt` - adaptive controller
- `ts`    - time stepper
- `h`     - current step size

Output Parameters:
- `next_sc` - optional, scheme to use for the next step
- `next_h`  - step size to use for the next step
- `accept`  - `PETSC_TRUE` to accept the current step, `PETSC_FALSE` to repeat the current step with the new step size

Level: developer

-seealso: [](ch_ts), `TSAdapt`, `TSAdaptCandidatesClear()`, `TSAdaptCandidateAdd()`

# External Links
$(_doc_external("TS/TSAdaptChoose"))
"""
function TSAdaptChoose(petsclib::PetscLibType, adapt::TSAdapt, ts::AbstractTS, h::PetscReal) end

@for_petsc function TSAdaptChoose(petsclib::$UnionPetscLib, adapt::TSAdapt, ts::AbstractTS, h::$PetscReal )
	next_sc_ = Ref{$PetscInt}()
	next_h_ = Ref{$PetscReal}()
	accept_ = Ref{PetscBool}()

    @chk ccall(
               (:TSAdaptChoose, $petsc_library),
               PetscErrorCode,
               (TSAdapt, CTS, $PetscReal, Ptr{$PetscInt}, Ptr{$PetscReal}, Ptr{PetscBool}),
               adapt, ts, h, next_sc_, next_h_, accept_,
              )

	next_sc = next_sc_[]
	next_h = next_h_[]
	accept = accept_[]

	return next_sc,next_h,accept
end 

"""
	inadapt::TSAdapt = TSAdaptCreate(petsclib::PetscLibType,comm::MPI_Comm) 
create an adaptive controller context for time stepping

Collective

Input Parameter:
- `comm` - The communicator

Output Parameter:
- `inadapt` - new `TSAdapt` object

Level: developer

-seealso: [](ch_ts), `TSAdapt`, `TSGetAdapt()`, `TSAdaptSetType()`, `TSAdaptDestroy()`

# External Links
$(_doc_external("TS/TSAdaptCreate"))
"""
function TSAdaptCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function TSAdaptCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	inadapt_ = Ref{TSAdapt}()

    @chk ccall(
               (:TSAdaptCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{TSAdapt}),
               comm, inadapt_,
              )

	inadapt = inadapt_[]

	return inadapt
end 

"""
	TSAdaptDSPSetFilter(petsclib::PetscLibType,adapt::TSAdapt, name::String) 
Sets internal parameters corresponding to the named filter {cite}`soderlind2006adaptive` {cite}`soderlind2003digital`

Collective

Input Parameters:
- `adapt` - adaptive controller context
- `name`  - filter name

Options Database Key:
- `-ts_adapt_dsp_filter <name>` - Sets predefined controller by name; use -help for a list of available controllers

Filter names:
- `basic`                       - similar to `TSADAPTBASIC` but with different criteria for step rejections.
- `PI30, PI42, PI33, PI34`      - PI controllers.
- `PC11, PC47, PC36`            - predictive controllers.
- `H0211, H211b, H211PI`        - digital filters with orders dynamics=2, adaptivity=1, filter=1.
- `H0312, H312b, H312PID`       - digital filters with orders dynamics=3, adaptivity=1, filter=2.
- `H0321, H321`                 - digital filters with orders dynamics=3, adaptivity=2, filter=1.

Level: intermediate

-seealso: [](ch_ts), `TSADAPTDSP`, `TS`, `TSAdapt`, `TSGetAdapt()`, `TSAdaptDSPSetPID()`

# External Links
$(_doc_external("TS/TSAdaptDSPSetFilter"))
"""
function TSAdaptDSPSetFilter(petsclib::PetscLibType, adapt::TSAdapt, name::String) end

@for_petsc function TSAdaptDSPSetFilter(petsclib::$UnionPetscLib, adapt::TSAdapt, name::String )

    @chk ccall(
               (:TSAdaptDSPSetFilter, $petsc_library),
               PetscErrorCode,
               (TSAdapt, Ptr{Cchar}),
               adapt, name,
              )


	return nothing
end 

"""
	TSAdaptDSPSetPID(petsclib::PetscLibType,adapt::TSAdapt, kkI::PetscReal, kkP::PetscReal, kkD::PetscReal) 
Set the PID controller parameters {cite}`soderlind2006adaptive`  {cite}`soderlind2003digital`

Input Parameters:
- `adapt` - adaptive controller context
- `kkI`   - Integral parameter
- `kkP`   - Proportional parameter
- `kkD`   - Derivative parameter

Options Database Key:
- `-ts_adapt_dsp_pid <kkI,kkP,kkD>` - Sets PID controller parameters

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSAdapt`, `TSGetAdapt()`, `TSAdaptDSPSetFilter()`

# External Links
$(_doc_external("TS/TSAdaptDSPSetPID"))
"""
function TSAdaptDSPSetPID(petsclib::PetscLibType, adapt::TSAdapt, kkI::PetscReal, kkP::PetscReal, kkD::PetscReal) end

@for_petsc function TSAdaptDSPSetPID(petsclib::$UnionPetscLib, adapt::TSAdapt, kkI::$PetscReal, kkP::$PetscReal, kkD::$PetscReal )

    @chk ccall(
               (:TSAdaptDSPSetPID, $petsc_library),
               PetscErrorCode,
               (TSAdapt, $PetscReal, $PetscReal, $PetscReal),
               adapt, kkI, kkP, kkD,
              )


	return nothing
end 

"""
	TSAdaptDestroy(petsclib::PetscLibType,adapt::Union{TSAdapt, Ref{TSAdapt}}) 

# External Links
$(_doc_external("TS/TSAdaptDestroy"))
"""
function TSAdaptDestroy(petsclib::PetscLibType, adapt::Union{TSAdapt, Ref{TSAdapt}}) end

@for_petsc function TSAdaptDestroy(petsclib::$UnionPetscLib, adapt::Union{TSAdapt, Ref{TSAdapt}} )
	adapt_ = adapt isa Base.RefValue ? adapt : Ref{TSAdapt}(adapt)

    @chk ccall(
               (:TSAdaptDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{TSAdapt},),
               adapt_,
              )


	return nothing
end 

"""
	TSAdaptFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the `TS` package. It is
called from `PetscFinalize()`.

Level: developer

-seealso: [](ch_ts), `PetscFinalize()`

# External Links
$(_doc_external("TS/TSAdaptFinalizePackage"))
"""
function TSAdaptFinalizePackage(petsclib::PetscLibType) end

@for_petsc function TSAdaptFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSAdaptFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	low::PetscReal,high::PetscReal = TSAdaptGetClip(petsclib::PetscLibType,adapt::TSAdapt) 
Gets the admissible decrease/increase factor in step size in the time step adapter

Not Collective

Input Parameter:
- `adapt` - adaptive controller context

Output Parameters:
- `low`  - optional, admissible decrease factor
- `high` - optional, admissible increase factor

Level: intermediate

-seealso: [](ch_ts), `TSAdapt`, `TSAdaptChoose()`, `TSAdaptSetClip()`, `TSAdaptSetScaleSolveFailed()`

# External Links
$(_doc_external("TS/TSAdaptGetClip"))
"""
function TSAdaptGetClip(petsclib::PetscLibType, adapt::TSAdapt) end

@for_petsc function TSAdaptGetClip(petsclib::$UnionPetscLib, adapt::TSAdapt )
	low_ = Ref{$PetscReal}()
	high_ = Ref{$PetscReal}()

    @chk ccall(
               (:TSAdaptGetClip, $petsc_library),
               PetscErrorCode,
               (TSAdapt, Ptr{$PetscReal}, Ptr{$PetscReal}),
               adapt, low_, high_,
              )

	low = low_[]
	high = high_[]

	return low,high
end 

"""
	max_ignore::PetscReal = TSAdaptGetMaxIgnore(petsclib::PetscLibType,adapt::TSAdapt) 
Get error estimation threshold. Solution components below this threshold value will not be considered when computing error norms
for time step adaptivity (in absolute value).

Not Collective

Input Parameter:
- `adapt` - adaptive controller context

Output Parameter:
- `max_ignore` - threshold for solution components that are ignored during error estimation

Level: intermediate

-seealso: [](ch_ts), `TSAdapt`, `TSAdaptSetMaxIgnore()`, `TSAdaptChoose()`

# External Links
$(_doc_external("TS/TSAdaptGetMaxIgnore"))
"""
function TSAdaptGetMaxIgnore(petsclib::PetscLibType, adapt::TSAdapt) end

@for_petsc function TSAdaptGetMaxIgnore(petsclib::$UnionPetscLib, adapt::TSAdapt )
	max_ignore_ = Ref{$PetscReal}()

    @chk ccall(
               (:TSAdaptGetMaxIgnore, $petsc_library),
               PetscErrorCode,
               (TSAdapt, Ptr{$PetscReal}),
               adapt, max_ignore_,
              )

	max_ignore = max_ignore_[]

	return max_ignore
end 

"""
	safety::PetscReal,reject_safety::PetscReal = TSAdaptGetSafety(petsclib::PetscLibType,adapt::TSAdapt) 
Get safety factors for time step adapter

Not Collective

Input Parameter:
- `adapt` - adaptive controller context

Output Parameters:
- `safety`        - safety factor relative to target error/stability goal
- `reject_safety` - extra safety factor to apply if the last step was rejected

Level: intermediate

-seealso: [](ch_ts), `TSAdapt`, `TSAdaptSetSafety()`, `TSAdaptChoose()`

# External Links
$(_doc_external("TS/TSAdaptGetSafety"))
"""
function TSAdaptGetSafety(petsclib::PetscLibType, adapt::TSAdapt) end

@for_petsc function TSAdaptGetSafety(petsclib::$UnionPetscLib, adapt::TSAdapt )
	safety_ = Ref{$PetscReal}()
	reject_safety_ = Ref{$PetscReal}()

    @chk ccall(
               (:TSAdaptGetSafety, $petsc_library),
               PetscErrorCode,
               (TSAdapt, Ptr{$PetscReal}, Ptr{$PetscReal}),
               adapt, safety_, reject_safety_,
              )

	safety = safety_[]
	reject_safety = reject_safety_[]

	return safety,reject_safety
end 

"""
	scale::PetscReal = TSAdaptGetScaleSolveFailed(petsclib::PetscLibType,adapt::TSAdapt) 
Gets the admissible decrease/increase factor in step size

Not Collective

Input Parameter:
- `adapt` - adaptive controller context

Output Parameter:
- `scale` - scale factor

Level: intermediate

-seealso: [](ch_ts), `TSAdapt`, `TSAdaptChoose()`, `TSAdaptSetScaleSolveFailed()`, `TSAdaptSetClip()`

# External Links
$(_doc_external("TS/TSAdaptGetScaleSolveFailed"))
"""
function TSAdaptGetScaleSolveFailed(petsclib::PetscLibType, adapt::TSAdapt) end

@for_petsc function TSAdaptGetScaleSolveFailed(petsclib::$UnionPetscLib, adapt::TSAdapt )
	scale_ = Ref{$PetscReal}()

    @chk ccall(
               (:TSAdaptGetScaleSolveFailed, $petsc_library),
               PetscErrorCode,
               (TSAdapt, Ptr{$PetscReal}),
               adapt, scale_,
              )

	scale = scale_[]

	return scale
end 

"""
	hmin::PetscReal,hmax::PetscReal = TSAdaptGetStepLimits(petsclib::PetscLibType,adapt::TSAdapt) 
Get the minimum and maximum step sizes to be considered by the time step controller

Not Collective

Input Parameter:
- `adapt` - time step adaptivity context, usually gotten with `TSGetAdapt()`

Output Parameters:
- `hmin` - minimum time step
- `hmax` - maximum time step

Level: intermediate

-seealso: [](ch_ts), `TSAdapt`, `TSAdaptSetStepLimits()`, `TSAdaptChoose()`

# External Links
$(_doc_external("TS/TSAdaptGetStepLimits"))
"""
function TSAdaptGetStepLimits(petsclib::PetscLibType, adapt::TSAdapt) end

@for_petsc function TSAdaptGetStepLimits(petsclib::$UnionPetscLib, adapt::TSAdapt )
	hmin_ = Ref{$PetscReal}()
	hmax_ = Ref{$PetscReal}()

    @chk ccall(
               (:TSAdaptGetStepLimits, $petsc_library),
               PetscErrorCode,
               (TSAdapt, Ptr{$PetscReal}, Ptr{$PetscReal}),
               adapt, hmin_, hmax_,
              )

	hmin = hmin_[]
	hmax = hmax_[]

	return hmin,hmax
end 

"""
	type::TSAdaptType = TSAdaptGetType(petsclib::PetscLibType,adapt::TSAdapt) 
gets the `TS` adapter method type (as a string).

Not Collective

Input Parameter:
- `adapt` - The `TS` adapter, most likely obtained with `TSGetAdapt()`

Output Parameter:
- `type` - The name of `TS` adapter method

Level: intermediate

-seealso: `TSAdapt`, `TSAdaptType`, `TSAdaptSetType()`

# External Links
$(_doc_external("TS/TSAdaptGetType"))
"""
function TSAdaptGetType(petsclib::PetscLibType, adapt::TSAdapt) end

@for_petsc function TSAdaptGetType(petsclib::$UnionPetscLib, adapt::TSAdapt )
	type_ = Ref{TSAdaptType}()

    @chk ccall(
               (:TSAdaptGetType, $petsc_library),
               PetscErrorCode,
               (TSAdapt, Ptr{TSAdaptType}),
               adapt, type_,
              )

	type = type_[] == C_NULL ? "" : unsafe_string(type_[])

	return type
end 

"""
	t::PetscReal,dt::PetscReal = TSAdaptHistoryGetStep(petsclib::PetscLibType,adapt::TSAdapt, step::PetscInt) 
Gets time and time step for a given step number in the history

Logically Collective

Input Parameters:
- `adapt` - the TSAdapt context
- `step`  - the step number

Output Parameters:
- `t`  - the time corresponding to the requested step (can be `NULL`)
- `dt` - the time step to be taken at the requested step (can be `NULL`)

Level: advanced

-seealso: [](ch_ts), `TS`, `TSGetAdapt()`, `TSAdaptSetType()`, `TSAdaptHistorySetTrajectory()`, `TSADAPTHISTORY`

# External Links
$(_doc_external("TS/TSAdaptHistoryGetStep"))
"""
function TSAdaptHistoryGetStep(petsclib::PetscLibType, adapt::TSAdapt, step::PetscInt) end

@for_petsc function TSAdaptHistoryGetStep(petsclib::$UnionPetscLib, adapt::TSAdapt, step::$PetscInt )
	t_ = Ref{$PetscReal}()
	dt_ = Ref{$PetscReal}()

    @chk ccall(
               (:TSAdaptHistoryGetStep, $petsc_library),
               PetscErrorCode,
               (TSAdapt, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}),
               adapt, step, t_, dt_,
              )

	t = t_[]
	dt = dt_[]

	return t,dt
end 

"""
	TSAdaptHistorySetHistory(petsclib::PetscLibType,adapt::TSAdapt, n::PetscInt, hist::Vector{PetscReal}, backward::PetscBool) 
Sets the time history in the adaptor

Logically Collective

Input Parameters:
- `adapt`    - the `TSAdapt` context
- `n`        - size of the time history
- `hist`     - the time history
- `backward` - if the time history has to be followed backward

Level: advanced

-seealso: [](ch_ts), `TSGetAdapt()`, `TSAdaptSetType()`, `TSAdaptHistorySetTrajectory()`, `TSADAPTHISTORY`

# External Links
$(_doc_external("TS/TSAdaptHistorySetHistory"))
"""
function TSAdaptHistorySetHistory(petsclib::PetscLibType, adapt::TSAdapt, n::PetscInt, hist::Vector{PetscReal}, backward::PetscBool) end

@for_petsc function TSAdaptHistorySetHistory(petsclib::$UnionPetscLib, adapt::TSAdapt, n::$PetscInt, hist::Vector{$PetscReal}, backward::PetscBool )

    @chk ccall(
               (:TSAdaptHistorySetHistory, $petsc_library),
               PetscErrorCode,
               (TSAdapt, $PetscInt, Ptr{$PetscReal}, PetscBool),
               adapt, n, hist, backward,
              )


	return nothing
end 

"""
	TSAdaptHistorySetTrajectory(petsclib::PetscLibType,adapt::TSAdapt, tj::TSTrajectory, backward::PetscBool) 
Sets a time history in the adaptor from a given `TSTrajectory`

Logically Collective

Input Parameters:
- `adapt`    - the `TSAdapt` context
- `tj`       - the `TSTrajectory` context
- `backward` - if the time history has to be followed backward

Level: advanced

-seealso: [](ch_ts), `TSGetAdapt()`, `TSAdaptSetType()`, `TSAdaptHistorySetHistory()`, `TSADAPTHISTORY`, `TSAdapt`

# External Links
$(_doc_external("TS/TSAdaptHistorySetTrajectory"))
"""
function TSAdaptHistorySetTrajectory(petsclib::PetscLibType, adapt::TSAdapt, tj::TSTrajectory, backward::PetscBool) end

@for_petsc function TSAdaptHistorySetTrajectory(petsclib::$UnionPetscLib, adapt::TSAdapt, tj::TSTrajectory, backward::PetscBool )

    @chk ccall(
               (:TSAdaptHistorySetTrajectory, $petsc_library),
               PetscErrorCode,
               (TSAdapt, TSTrajectory, PetscBool),
               adapt, tj, backward,
              )


	return nothing
end 

"""
	TSAdaptInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `TSAdapt` package. It is
called from `TSInitializePackage()`.

Level: developer

-seealso: [](ch_ts), `PetscInitialize()`

# External Links
$(_doc_external("TS/TSAdaptInitializePackage"))
"""
function TSAdaptInitializePackage(petsclib::PetscLibType) end

@for_petsc function TSAdaptInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSAdaptInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSAdaptLoad(petsclib::PetscLibType,adapt::TSAdapt, viewer::PetscViewer) 
Loads a TSAdapt that has been stored in binary with `TSAdaptView()`.

Collective

Input Parameters:
- `adapt`  - the newly loaded `TSAdapt`, this needs to have been created with `TSAdaptCreate()` or
some related function before a call to `TSAdaptLoad()`.
- `viewer` - binary file viewer, obtained from `PetscViewerBinaryOpen()` or
HDF5 file viewer, obtained from `PetscViewerHDF5Open()`

Level: intermediate

-seealso: [](ch_ts), `PetscViewerBinaryOpen()`, `TSAdaptView()`, `MatLoad()`, `VecLoad()`, `TSAdapt`

# External Links
$(_doc_external("TS/TSAdaptLoad"))
"""
function TSAdaptLoad(petsclib::PetscLibType, adapt::TSAdapt, viewer::PetscViewer) end

@for_petsc function TSAdaptLoad(petsclib::$UnionPetscLib, adapt::TSAdapt, viewer::PetscViewer )

    @chk ccall(
               (:TSAdaptLoad, $petsc_library),
               PetscErrorCode,
               (TSAdapt, PetscViewer),
               adapt, viewer,
              )


	return nothing
end 

"""
	TSAdaptRegister(petsclib::PetscLibType,sname::String, fnc::external) 
adds a TSAdapt implementation

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - name of user-defined adaptivity scheme
- `function` - routine to create method context

Level: advanced

-seealso: [](ch_ts), `TSAdaptRegisterAll()`

# External Links
$(_doc_external("TS/TSAdaptRegister"))
"""
function TSAdaptRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function TSAdaptRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:TSAdaptRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	TSAdaptReset(petsclib::PetscLibType,adapt::TSAdapt) 
Resets a `TSAdapt` context to its defaults

Collective

Input Parameter:
- `adapt` - the `TSAdapt` context obtained from `TSGetAdapt()` or `TSAdaptCreate()`

Level: developer

-seealso: [](ch_ts), `TSGetAdapt()`, `TSAdapt`, `TSAdaptCreate()`, `TSAdaptDestroy()`

# External Links
$(_doc_external("TS/TSAdaptReset"))
"""
function TSAdaptReset(petsclib::PetscLibType, adapt::TSAdapt) end

@for_petsc function TSAdaptReset(petsclib::$UnionPetscLib, adapt::TSAdapt )

    @chk ccall(
               (:TSAdaptReset, $petsc_library),
               PetscErrorCode,
               (TSAdapt,),
               adapt,
              )


	return nothing
end 

"""
	TSAdaptSetAlwaysAccept(petsclib::PetscLibType,adapt::TSAdapt, flag::PetscBool) 
Set whether to always accept steps regardless of
any error or stability condition not meeting the prescribed goal.

Logically Collective

Input Parameters:
- `adapt` - time step adaptivity context, usually gotten with `TSGetAdapt()`
- `flag`  - whether to always accept steps

Options Database Key:
- `-ts_adapt_always_accept` - to always accept steps

Level: intermediate

-seealso: [](ch_ts), `TSAdapt`, `TSGetAdapt()`, `TSAdaptChoose()`

# External Links
$(_doc_external("TS/TSAdaptSetAlwaysAccept"))
"""
function TSAdaptSetAlwaysAccept(petsclib::PetscLibType, adapt::TSAdapt, flag::PetscBool) end

@for_petsc function TSAdaptSetAlwaysAccept(petsclib::$UnionPetscLib, adapt::TSAdapt, flag::PetscBool )

    @chk ccall(
               (:TSAdaptSetAlwaysAccept, $petsc_library),
               PetscErrorCode,
               (TSAdapt, PetscBool),
               adapt, flag,
              )


	return nothing
end 

"""
	TSAdaptSetCheckStage(petsclib::PetscLibType,adapt::TSAdapt, func::external) 
Set a callback to check convergence for a stage

Logically Collective

Input Parameters:
- `adapt` - adaptive controller context
- `func`  - stage check function

Calling sequence:
- `adapt`  - adaptive controller context
- `ts`     - time stepping context
- `t`      - current time
- `Y`      - current solution vector
- `accept` - pending choice of whether to accept, can be modified by this routine

Level: advanced

-seealso: [](ch_ts), `TSAdapt`, `TSGetAdapt()`, `TSAdaptChoose()`

# External Links
$(_doc_external("TS/TSAdaptSetCheckStage"))
"""
function TSAdaptSetCheckStage(petsclib::PetscLibType, adapt::TSAdapt, func::external) end

@for_petsc function TSAdaptSetCheckStage(petsclib::$UnionPetscLib, adapt::TSAdapt, func::external )

    @chk ccall(
               (:TSAdaptSetCheckStage, $petsc_library),
               PetscErrorCode,
               (TSAdapt, external),
               adapt, func,
              )


	return nothing
end 

"""
	TSAdaptSetClip(petsclib::PetscLibType,adapt::TSAdapt, low::PetscReal, high::PetscReal) 
Sets the admissible decrease/increase factor in step size in the time step adapter

Logically collective

Input Parameters:
- `adapt` - adaptive controller context
- `low`   - admissible decrease factor
- `high`  - admissible increase factor

Options Database Key:
- `-ts_adapt_clip <low>,<high>` - to set admissible time step decrease and increase factors

Level: intermediate

-seealso: [](ch_ts), `TSAdapt`, `TSAdaptChoose()`, `TSAdaptGetClip()`, `TSAdaptSetScaleSolveFailed()`

# External Links
$(_doc_external("TS/TSAdaptSetClip"))
"""
function TSAdaptSetClip(petsclib::PetscLibType, adapt::TSAdapt, low::PetscReal, high::PetscReal) end

@for_petsc function TSAdaptSetClip(petsclib::$UnionPetscLib, adapt::TSAdapt, low::$PetscReal, high::$PetscReal )

    @chk ccall(
               (:TSAdaptSetClip, $petsc_library),
               PetscErrorCode,
               (TSAdapt, $PetscReal, $PetscReal),
               adapt, low, high,
              )


	return nothing
end 

"""
	TSAdaptSetFromOptions(petsclib::PetscLibType,adapt::TSAdapt, PetscOptionsObject::PetscOptionItems) 
Sets various `TSAdapt` parameters from user options.

Collective

Input Parameters:
- `adapt`              - the `TSAdapt` context
- `PetscOptionsObject` - object created by `PetscOptionsBegin()`

Options Database Keys:
- `-ts_adapt_type <type>`                - algorithm to use for adaptivity
- `-ts_adapt_always_accept`              - always accept steps regardless of error/stability goals
- `-ts_adapt_safety <safety>`            - safety factor relative to target error/stability goal
- `-ts_adapt_reject_safety <safety>`     - extra safety factor to apply if the last step was rejected
- `-ts_adapt_clip <low,high>`            - admissible time step decrease and increase factors
- `-ts_adapt_dt_min <min>`               - minimum timestep to use
- `-ts_adapt_dt_max <max>`               - maximum timestep to use
- `-ts_adapt_scale_solve_failed <scale>` - scale timestep by this factor if a solve fails
- `-ts_adapt_wnormtype <2 or infinity>`  - type of norm for computing error estimates
- `-ts_adapt_time_step_increase_delay`   - number of timesteps to delay increasing the time step after it has been decreased due to failed solver

Level: advanced

-seealso: [](ch_ts), `TSAdapt`, `TSGetAdapt()`, `TSAdaptSetType()`, `TSAdaptSetAlwaysAccept()`, `TSAdaptSetSafety()`,
`TSAdaptSetClip()`, `TSAdaptSetScaleSolveFailed()`, `TSAdaptSetStepLimits()`, `TSAdaptSetMonitor()`

# External Links
$(_doc_external("TS/TSAdaptSetFromOptions"))
"""
function TSAdaptSetFromOptions(petsclib::PetscLibType, adapt::TSAdapt, PetscOptionsObject::PetscOptionItems) end

@for_petsc function TSAdaptSetFromOptions(petsclib::$UnionPetscLib, adapt::TSAdapt, PetscOptionsObject::PetscOptionItems )

    @chk ccall(
               (:TSAdaptSetFromOptions, $petsc_library),
               PetscErrorCode,
               (TSAdapt, PetscOptionItems),
               adapt, PetscOptionsObject,
              )


	return nothing
end 

"""
	TSAdaptSetMaxIgnore(petsclib::PetscLibType,adapt::TSAdapt, max_ignore::PetscReal) 
Set error estimation threshold. Solution components below this threshold value will not be considered when computing error norms
for time step adaptivity (in absolute value). A negative value (default) of the threshold leads to considering all solution components.

Logically Collective

Input Parameters:
- `adapt`      - adaptive controller context
- `max_ignore` - threshold for solution components that are ignored during error estimation

Options Database Key:
- `-ts_adapt_max_ignore <max_ignore>` - to set the threshold

Level: intermediate

-seealso: [](ch_ts), `TSAdapt`, `TSAdaptGetMaxIgnore()`, `TSAdaptChoose()`

# External Links
$(_doc_external("TS/TSAdaptSetMaxIgnore"))
"""
function TSAdaptSetMaxIgnore(petsclib::PetscLibType, adapt::TSAdapt, max_ignore::PetscReal) end

@for_petsc function TSAdaptSetMaxIgnore(petsclib::$UnionPetscLib, adapt::TSAdapt, max_ignore::$PetscReal )

    @chk ccall(
               (:TSAdaptSetMaxIgnore, $petsc_library),
               PetscErrorCode,
               (TSAdapt, $PetscReal),
               adapt, max_ignore,
              )


	return nothing
end 

"""
	TSAdaptSetMonitor(petsclib::PetscLibType,adapt::TSAdapt, flg::PetscBool) 
Monitor the choices made by the adaptive controller

Collective

Input Parameters:
- `adapt` - adaptive controller context
- `flg`   - `PETSC_TRUE` to active a monitor, `PETSC_FALSE` to disable

Options Database Key:
- `-ts_adapt_monitor` - to turn on monitoring

Level: intermediate

-seealso: [](ch_ts), `TSAdapt`, `TSGetAdapt()`, `TSAdaptChoose()`

# External Links
$(_doc_external("TS/TSAdaptSetMonitor"))
"""
function TSAdaptSetMonitor(petsclib::PetscLibType, adapt::TSAdapt, flg::PetscBool) end

@for_petsc function TSAdaptSetMonitor(petsclib::$UnionPetscLib, adapt::TSAdapt, flg::PetscBool )

    @chk ccall(
               (:TSAdaptSetMonitor, $petsc_library),
               PetscErrorCode,
               (TSAdapt, PetscBool),
               adapt, flg,
              )


	return nothing
end 

"""
	TSAdaptSetOptionsPrefix(petsclib::PetscLibType,adapt::TSAdapt, prefix::String) 

# External Links
$(_doc_external("TS/TSAdaptSetOptionsPrefix"))
"""
function TSAdaptSetOptionsPrefix(petsclib::PetscLibType, adapt::TSAdapt, prefix::String) end

@for_petsc function TSAdaptSetOptionsPrefix(petsclib::$UnionPetscLib, adapt::TSAdapt, prefix::String )

    @chk ccall(
               (:TSAdaptSetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (TSAdapt, Ptr{Cchar}),
               adapt, prefix,
              )


	return nothing
end 

"""
	TSAdaptSetSafety(petsclib::PetscLibType,adapt::TSAdapt, safety::PetscReal, reject_safety::PetscReal) 
Set safety factors for time step adaptor

Logically Collective

Input Parameters:
- `adapt`         - adaptive controller context
- `safety`        - safety factor relative to target error/stability goal
- `reject_safety` - extra safety factor to apply if the last step was rejected

Options Database Keys:
- `-ts_adapt_safety <safety>`               - to set safety factor
- `-ts_adapt_reject_safety <reject_safety>` - to set reject safety factor

Level: intermediate

-seealso: [](ch_ts), `TSAdapt`, `TSAdaptGetSafety()`, `TSAdaptChoose()`

# External Links
$(_doc_external("TS/TSAdaptSetSafety"))
"""
function TSAdaptSetSafety(petsclib::PetscLibType, adapt::TSAdapt, safety::PetscReal, reject_safety::PetscReal) end

@for_petsc function TSAdaptSetSafety(petsclib::$UnionPetscLib, adapt::TSAdapt, safety::$PetscReal, reject_safety::$PetscReal )

    @chk ccall(
               (:TSAdaptSetSafety, $petsc_library),
               PetscErrorCode,
               (TSAdapt, $PetscReal, $PetscReal),
               adapt, safety, reject_safety,
              )


	return nothing
end 

"""
	TSAdaptSetScaleSolveFailed(petsclib::PetscLibType,adapt::TSAdapt, scale::PetscReal) 
Scale step size by this factor if solve fails

Logically Collective

Input Parameters:
- `adapt` - adaptive controller context
- `scale` - scale

Options Database Key:
- `-ts_adapt_scale_solve_failed <scale>` - to set scale step by this factor if solve fails

Level: intermediate

-seealso: [](ch_ts), `TSAdapt`, `TSAdaptChoose()`, `TSAdaptGetScaleSolveFailed()`, `TSAdaptGetClip()`

# External Links
$(_doc_external("TS/TSAdaptSetScaleSolveFailed"))
"""
function TSAdaptSetScaleSolveFailed(petsclib::PetscLibType, adapt::TSAdapt, scale::PetscReal) end

@for_petsc function TSAdaptSetScaleSolveFailed(petsclib::$UnionPetscLib, adapt::TSAdapt, scale::$PetscReal )

    @chk ccall(
               (:TSAdaptSetScaleSolveFailed, $petsc_library),
               PetscErrorCode,
               (TSAdapt, $PetscReal),
               adapt, scale,
              )


	return nothing
end 

"""
	TSAdaptSetStepLimits(petsclib::PetscLibType,adapt::TSAdapt, hmin::PetscReal, hmax::PetscReal) 
Set the minimum and maximum step sizes to be considered by the time step controller

Logically Collective

Input Parameters:
- `adapt` - time step adaptivity context, usually gotten with `TSGetAdapt()`
- `hmin`  - minimum time step
- `hmax`  - maximum time step

Options Database Keys:
- `-ts_adapt_dt_min <min>` - to set minimum time step
- `-ts_adapt_dt_max <max>` - to set maximum time step

Level: intermediate

-seealso: [](ch_ts), `TSAdapt`, `TSAdaptGetStepLimits()`, `TSAdaptChoose()`

# External Links
$(_doc_external("TS/TSAdaptSetStepLimits"))
"""
function TSAdaptSetStepLimits(petsclib::PetscLibType, adapt::TSAdapt, hmin::PetscReal, hmax::PetscReal) end

@for_petsc function TSAdaptSetStepLimits(petsclib::$UnionPetscLib, adapt::TSAdapt, hmin::$PetscReal, hmax::$PetscReal )

    @chk ccall(
               (:TSAdaptSetStepLimits, $petsc_library),
               PetscErrorCode,
               (TSAdapt, $PetscReal, $PetscReal),
               adapt, hmin, hmax,
              )


	return nothing
end 

"""
	TSAdaptSetTimeStepIncreaseDelay(petsclib::PetscLibType,adapt::TSAdapt, cnt::PetscInt) 
The number of timesteps to wait after a decrease in the timestep due to failed solver
before increasing the time step.

Logicially Collective

Input Parameters:
- `adapt` - adaptive controller context
- `cnt`   - the number of timesteps

Options Database Key:
- `-ts_adapt_time_step_increase_delay cnt` - number of steps to delay the increase

Level: advanced

-seealso: [](ch_ts), `TSAdapt`

# External Links
$(_doc_external("TS/TSAdaptSetTimeStepIncreaseDelay"))
"""
function TSAdaptSetTimeStepIncreaseDelay(petsclib::PetscLibType, adapt::TSAdapt, cnt::PetscInt) end

@for_petsc function TSAdaptSetTimeStepIncreaseDelay(petsclib::$UnionPetscLib, adapt::TSAdapt, cnt::$PetscInt )

    @chk ccall(
               (:TSAdaptSetTimeStepIncreaseDelay, $petsc_library),
               PetscErrorCode,
               (TSAdapt, $PetscInt),
               adapt, cnt,
              )


	return nothing
end 

"""
	TSAdaptSetType(petsclib::PetscLibType,adapt::TSAdapt, type::TSAdaptType) 
sets the approach used for the error adapter

Logicially Collective

Input Parameters:
- `adapt` - the `TS` adapter, most likely obtained with `TSGetAdapt()`
- `type`  - one of the `TSAdaptType`

Options Database Key:
- `-ts_adapt_type <basic or dsp or none>` - to set the adapter type

Level: intermediate

-seealso: [](ch_ts), `TSGetAdapt()`, `TSAdaptDestroy()`, `TSAdaptType`, `TSAdaptGetType()`

# External Links
$(_doc_external("TS/TSAdaptSetType"))
"""
function TSAdaptSetType(petsclib::PetscLibType, adapt::TSAdapt, type::TSAdaptType) end

@for_petsc function TSAdaptSetType(petsclib::$UnionPetscLib, adapt::TSAdapt, type::TSAdaptType )

    @chk ccall(
               (:TSAdaptSetType, $petsc_library),
               PetscErrorCode,
               (TSAdapt, TSAdaptType),
               adapt, type,
              )


	return nothing
end 

"""
	TSAdaptView(petsclib::PetscLibType,adapt::TSAdapt, viewer::PetscViewer) 

# External Links
$(_doc_external("TS/TSAdaptView"))
"""
function TSAdaptView(petsclib::PetscLibType, adapt::TSAdapt, viewer::PetscViewer) end

@for_petsc function TSAdaptView(petsclib::$UnionPetscLib, adapt::TSAdapt, viewer::PetscViewer )

    @chk ccall(
               (:TSAdaptView, $petsc_library),
               PetscErrorCode,
               (TSAdapt, PetscViewer),
               adapt, viewer,
              )


	return nothing
end 

"""
	next_sc::PetscInt,next_h::PetscReal,finish::PetscBool = TSGLLEAdaptChoose(petsclib::PetscLibType,adapt::TSGLLEAdapt, n::PetscInt, orders::Vector{PetscInt}, errors::Vector{PetscReal}, cost::Vector{PetscReal}, cur::PetscInt, h::PetscReal, tleft::PetscReal) 

# External Links
$(_doc_external("TS/TSGLLEAdaptChoose"))
"""
function TSGLLEAdaptChoose(petsclib::PetscLibType, adapt::TSGLLEAdapt, n::PetscInt, orders::Vector{PetscInt}, errors::Vector{PetscReal}, cost::Vector{PetscReal}, cur::PetscInt, h::PetscReal, tleft::PetscReal) end

@for_petsc function TSGLLEAdaptChoose(petsclib::$UnionPetscLib, adapt::TSGLLEAdapt, n::$PetscInt, orders::Vector{$PetscInt}, errors::Vector{$PetscReal}, cost::Vector{$PetscReal}, cur::$PetscInt, h::$PetscReal, tleft::$PetscReal )
	next_sc_ = Ref{$PetscInt}()
	next_h_ = Ref{$PetscReal}()
	finish_ = Ref{PetscBool}()

    @chk ccall(
               (:TSGLLEAdaptChoose, $petsc_library),
               PetscErrorCode,
               (TSGLLEAdapt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscReal}, Ptr{$PetscReal}, $PetscInt, $PetscReal, $PetscReal, Ptr{$PetscInt}, Ptr{$PetscReal}, Ptr{PetscBool}),
               adapt, n, orders, errors, cost, cur, h, tleft, next_sc_, next_h_, finish_,
              )

	next_sc = next_sc_[]
	next_h = next_h_[]
	finish = finish_[]

	return next_sc,next_h,finish
end 

"""
	inadapt::TSGLLEAdapt = TSGLLEAdaptCreate(petsclib::PetscLibType,comm::MPI_Comm) 

# External Links
$(_doc_external("TS/TSGLLEAdaptCreate"))
"""
function TSGLLEAdaptCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function TSGLLEAdaptCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	inadapt_ = Ref{TSGLLEAdapt}()

    @chk ccall(
               (:TSGLLEAdaptCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{TSGLLEAdapt}),
               comm, inadapt_,
              )

	inadapt = inadapt_[]

	return inadapt
end 

"""
	TSGLLEAdaptDestroy(petsclib::PetscLibType,adapt::Union{TSGLLEAdapt, Ref{TSGLLEAdapt}}) 

# External Links
$(_doc_external("TS/TSGLLEAdaptDestroy"))
"""
function TSGLLEAdaptDestroy(petsclib::PetscLibType, adapt::Union{TSGLLEAdapt, Ref{TSGLLEAdapt}}) end

@for_petsc function TSGLLEAdaptDestroy(petsclib::$UnionPetscLib, adapt::Union{TSGLLEAdapt, Ref{TSGLLEAdapt}} )
	adapt_ = adapt isa Base.RefValue ? adapt : Ref{TSGLLEAdapt}(adapt)

    @chk ccall(
               (:TSGLLEAdaptDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{TSGLLEAdapt},),
               adapt_,
              )


	return nothing
end 

"""
	TSGLLEAdaptFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the `TSGLLE` package. It is
called from `PetscFinalize()`.

Level: developer

-seealso: [](ch_ts), `PetscFinalize()`, `TSGLLEAdapt`, `TSGLLEAdaptInitializePackage()`

# External Links
$(_doc_external("TS/TSGLLEAdaptFinalizePackage"))
"""
function TSGLLEAdaptFinalizePackage(petsclib::PetscLibType) end

@for_petsc function TSGLLEAdaptFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSGLLEAdaptFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSGLLEAdaptInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `TSGLLEAdapt` package. It is
called from `TSInitializePackage()`.

Level: developer

-seealso: [](ch_ts), `PetscInitialize()`, `TSGLLEAdapt`, `TSGLLEAdaptFinalizePackage()`

# External Links
$(_doc_external("TS/TSGLLEAdaptInitializePackage"))
"""
function TSGLLEAdaptInitializePackage(petsclib::PetscLibType) end

@for_petsc function TSGLLEAdaptInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSGLLEAdaptInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSGLLEAdaptRegister(petsclib::PetscLibType,sname::String, fnc::external) 
adds a `TSGLLEAdapt` implementation

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - name of user-defined adaptivity scheme
- `function` - routine to create method context

Level: advanced

-seealso: [](ch_ts), `TSGLLE`, `TSGLLEAdapt`, `TSGLLEAdaptRegisterAll()`

# External Links
$(_doc_external("TS/TSGLLEAdaptRegister"))
"""
function TSGLLEAdaptRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function TSGLLEAdaptRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:TSGLLEAdaptRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	TSGLLEAdaptSetFromOptions(petsclib::PetscLibType,adapt::TSGLLEAdapt, PetscOptionsObject::PetscOptionItems) 

# External Links
$(_doc_external("TS/TSGLLEAdaptSetFromOptions"))
"""
function TSGLLEAdaptSetFromOptions(petsclib::PetscLibType, adapt::TSGLLEAdapt, PetscOptionsObject::PetscOptionItems) end

@for_petsc function TSGLLEAdaptSetFromOptions(petsclib::$UnionPetscLib, adapt::TSGLLEAdapt, PetscOptionsObject::PetscOptionItems )

    @chk ccall(
               (:TSGLLEAdaptSetFromOptions, $petsc_library),
               PetscErrorCode,
               (TSGLLEAdapt, PetscOptionItems),
               adapt, PetscOptionsObject,
              )


	return nothing
end 

"""
	TSGLLEAdaptSetOptionsPrefix(petsclib::PetscLibType,adapt::TSGLLEAdapt, prefix::String) 

# External Links
$(_doc_external("TS/TSGLLEAdaptSetOptionsPrefix"))
"""
function TSGLLEAdaptSetOptionsPrefix(petsclib::PetscLibType, adapt::TSGLLEAdapt, prefix::String) end

@for_petsc function TSGLLEAdaptSetOptionsPrefix(petsclib::$UnionPetscLib, adapt::TSGLLEAdapt, prefix::String )

    @chk ccall(
               (:TSGLLEAdaptSetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (TSGLLEAdapt, Ptr{Cchar}),
               adapt, prefix,
              )


	return nothing
end 

"""
	TSGLLEAdaptSetType(petsclib::PetscLibType,adapt::TSGLLEAdapt, type::TSGLLEAdaptType) 

# External Links
$(_doc_external("TS/TSGLLEAdaptSetType"))
"""
function TSGLLEAdaptSetType(petsclib::PetscLibType, adapt::TSGLLEAdapt, type::TSGLLEAdaptType) end

@for_petsc function TSGLLEAdaptSetType(petsclib::$UnionPetscLib, adapt::TSGLLEAdapt, type::TSGLLEAdaptType )

    @chk ccall(
               (:TSGLLEAdaptSetType, $petsc_library),
               PetscErrorCode,
               (TSGLLEAdapt, TSGLLEAdaptType),
               adapt, type,
              )


	return nothing
end 

"""
	TSGLLEAdaptView(petsclib::PetscLibType,adapt::TSGLLEAdapt, viewer::PetscViewer) 

# External Links
$(_doc_external("TS/TSGLLEAdaptView"))
"""
function TSGLLEAdaptView(petsclib::PetscLibType, adapt::TSGLLEAdapt, viewer::PetscViewer) end

@for_petsc function TSGLLEAdaptView(petsclib::$UnionPetscLib, adapt::TSGLLEAdapt, viewer::PetscViewer )

    @chk ccall(
               (:TSGLLEAdaptView, $petsc_library),
               PetscErrorCode,
               (TSGLLEAdapt, PetscViewer),
               adapt, viewer,
              )


	return nothing
end 

"""
	ctx::TSMonitorDrawCtx = TSMonitorDrawCtxCreate(petsclib::PetscLibType,comm::MPI_Comm, host::String, label::String, x::Cint, y::Cint, m::Cint, n::Cint, howoften::PetscInt) 
Creates the monitor context for `TSMonitorDrawCtx`

Collective

Input Parameters:
- `comm`     - the MPI communicator to use
- `host`     - the X display to open, or `NULL` for the local machine
- `label`    - the title to put in the title bar
- `x`        - the x screen coordinates of the upper left coordinate of the window
- `y`        - the y screen coordinates of the upper left coordinate of the window
- `m`        - the screen width in pixels
- `n`        - the screen height in pixels
- `howoften` - if positive then determines the frequency of the plotting, if -1 then only at the final time

Output Parameter:
- `ctx` - the monitor context

Options Database Keys:
- `-ts_monitor_draw_solution`         - draw the solution at each time-step
- `-ts_monitor_draw_solution_initial` - show initial solution as well as current solution

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSMonitorDrawCtxDestroy()`, `TSMonitorSet()`, `TSMonitorDefault()`, `VecView()`, `TSMonitorDrawCtx`, `PetscMonitorDrawSolution()`

# External Links
$(_doc_external("TS/TSMonitorDrawCtxCreate"))
"""
function TSMonitorDrawCtxCreate(petsclib::PetscLibType, comm::MPI_Comm, host::String, label::String, x::Cint, y::Cint, m::Cint, n::Cint, howoften::PetscInt) end

@for_petsc function TSMonitorDrawCtxCreate(petsclib::$UnionPetscLib, comm::MPI_Comm, host::String, label::String, x::Cint, y::Cint, m::Cint, n::Cint, howoften::$PetscInt )
	ctx_ = Ref{TSMonitorDrawCtx}()

    @chk ccall(
               (:TSMonitorDrawCtxCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Cint, Cint, Cint, Cint, $PetscInt, Ptr{TSMonitorDrawCtx}),
               comm, host, label, x, y, m, n, howoften, ctx_,
              )

	ctx = ctx_[]

	return ctx
end 

"""
	TSMonitorDrawCtxDestroy(petsclib::PetscLibType,ictx::Union{TSMonitorDrawCtx, Ref{TSMonitorDrawCtx}}) 
Destroys the monitor context for `TSMonitorDrawSolution()`

Collective

Input Parameter:
- `ictx` - the monitor context

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSMonitorSet()`, `TSMonitorDefault()`, `VecView()`, `TSMonitorDrawSolution()`, `TSMonitorDrawError()`, `TSMonitorDrawCtx`

# External Links
$(_doc_external("TS/TSMonitorDrawCtxDestroy"))
"""
function TSMonitorDrawCtxDestroy(petsclib::PetscLibType, ictx::Union{TSMonitorDrawCtx, Ref{TSMonitorDrawCtx}}) end

@for_petsc function TSMonitorDrawCtxDestroy(petsclib::$UnionPetscLib, ictx::Union{TSMonitorDrawCtx, Ref{TSMonitorDrawCtx}} )
	ictx_ = ictx isa Base.RefValue ? ictx : Ref{TSMonitorDrawCtx}(ictx)

    @chk ccall(
               (:TSMonitorDrawCtxDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{TSMonitorDrawCtx},),
               ictx_,
              )


	return nothing
end 

"""
	ctx::TSMonitorEnvelopeCtx = TSMonitorEnvelopeCtxCreate(petsclib::PetscLibType,ts::AbstractTS) 
Creates a context for use with `TSMonitorEnvelope()`

Collective

Input Parameter:
- `ts` - the `TS` solver object

Output Parameter:
- `ctx` - the context

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSMonitorLGTimeStep()`, `TSMonitorSet()`, `TSMonitorLGSolution()`, `TSMonitorLGError()`

# External Links
$(_doc_external("TS/TSMonitorEnvelopeCtxCreate"))
"""
function TSMonitorEnvelopeCtxCreate(petsclib::PetscLibType, ts::AbstractTS) end

@for_petsc function TSMonitorEnvelopeCtxCreate(petsclib::$UnionPetscLib, ts::AbstractTS )
	ctx_ = Ref{TSMonitorEnvelopeCtx}()

    @chk ccall(
               (:TSMonitorEnvelopeCtxCreate, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{TSMonitorEnvelopeCtx}),
               ts, ctx_,
              )

	ctx = ctx_[]

	return ctx
end 

"""
	TSMonitorEnvelopeCtxDestroy(petsclib::PetscLibType,ctx::Union{TSMonitorEnvelopeCtx, Ref{TSMonitorEnvelopeCtx}}) 
Destroys a context that was created  with `TSMonitorEnvelopeCtxCreate()`.

Collective

Input Parameter:
- `ctx` - the monitor context

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSMonitorLGCtxCreate()`, `TSMonitorSet()`, `TSMonitorLGTimeStep()`

# External Links
$(_doc_external("TS/TSMonitorEnvelopeCtxDestroy"))
"""
function TSMonitorEnvelopeCtxDestroy(petsclib::PetscLibType, ctx::Union{TSMonitorEnvelopeCtx, Ref{TSMonitorEnvelopeCtx}}) end

@for_petsc function TSMonitorEnvelopeCtxDestroy(petsclib::$UnionPetscLib, ctx::Union{TSMonitorEnvelopeCtx, Ref{TSMonitorEnvelopeCtx}} )
	ctx_ = ctx isa Base.RefValue ? ctx : Ref{TSMonitorEnvelopeCtx}(ctx)

    @chk ccall(
               (:TSMonitorEnvelopeCtxDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{TSMonitorEnvelopeCtx},),
               ctx_,
              )


	return nothing
end 

"""
	ctx::TSMonitorHGCtx = TSMonitorHGCtxCreate(petsclib::PetscLibType,comm::MPI_Comm, host::String, label::String, x::Cint, y::Cint, m::Cint, n::Cint, howoften::PetscInt, Ns::PetscInt, Nb::PetscInt, velocity::PetscBool) 

# External Links
$(_doc_external("TS/TSMonitorHGCtxCreate"))
"""
function TSMonitorHGCtxCreate(petsclib::PetscLibType, comm::MPI_Comm, host::String, label::String, x::Cint, y::Cint, m::Cint, n::Cint, howoften::PetscInt, Ns::PetscInt, Nb::PetscInt, velocity::PetscBool) end

@for_petsc function TSMonitorHGCtxCreate(petsclib::$UnionPetscLib, comm::MPI_Comm, host::String, label::String, x::Cint, y::Cint, m::Cint, n::Cint, howoften::$PetscInt, Ns::$PetscInt, Nb::$PetscInt, velocity::PetscBool )
	ctx_ = Ref{TSMonitorHGCtx}()

    @chk ccall(
               (:TSMonitorHGCtxCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Cint, Cint, Cint, Cint, $PetscInt, $PetscInt, $PetscInt, PetscBool, Ptr{TSMonitorHGCtx}),
               comm, host, label, x, y, m, n, howoften, Ns, Nb, velocity, ctx_,
              )

	ctx = ctx_[]

	return ctx
end 

"""
	TSMonitorHGCtxDestroy(petsclib::PetscLibType,ctx::Union{TSMonitorHGCtx, Ref{TSMonitorHGCtx}}) 

# External Links
$(_doc_external("TS/TSMonitorHGCtxDestroy"))
"""
function TSMonitorHGCtxDestroy(petsclib::PetscLibType, ctx::Union{TSMonitorHGCtx, Ref{TSMonitorHGCtx}}) end

@for_petsc function TSMonitorHGCtxDestroy(petsclib::$UnionPetscLib, ctx::Union{TSMonitorHGCtx, Ref{TSMonitorHGCtx}} )
	ctx_ = ctx isa Base.RefValue ? ctx : Ref{TSMonitorHGCtx}(ctx)

    @chk ccall(
               (:TSMonitorHGCtxDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{TSMonitorHGCtx},),
               ctx_,
              )


	return nothing
end 

"""
	ctx::TSMonitorLGCtx = TSMonitorLGCtxCreate(petsclib::PetscLibType,comm::MPI_Comm, host::String, label::String, x::Cint, y::Cint, m::Cint, n::Cint, howoften::PetscInt) 
Creates a `TSMonitorLGCtx` context for use with
`TS` to monitor the solution process graphically in various ways

Collective

Input Parameters:
- `comm`     - the MPI communicator to use
- `host`     - the X display to open, or `NULL` for the local machine
- `label`    - the title to put in the title bar
- `x`        - the x screen coordinates of the upper left coordinate of the window
- `y`        - the y screen coordinates of the upper left coordinate of the window
- `m`        - the screen width in pixels
- `n`        - the screen height in pixels
- `howoften` - if positive then determines the frequency of the plotting, if -1 then only at the final time

Output Parameter:
- `ctx` - the context

Options Database Keys:
- `-ts_monitor_lg_timestep`        - automatically sets line graph monitor
- `-ts_monitor_lg_timestep_log`    - automatically sets line graph monitor
- `-ts_monitor_lg_solution`        - monitor the solution (or certain values of the solution by calling `TSMonitorLGSetDisplayVariables()` or `TSMonitorLGCtxSetDisplayVariables()`)
- `-ts_monitor_lg_error`           - monitor the error
- `-ts_monitor_lg_ksp_iterations`  - monitor the number of `KSP` iterations needed for each timestep
- `-ts_monitor_lg_snes_iterations` - monitor the number of `SNES` iterations needed for each timestep
- `-lg_use_markers <true,false>`   - mark the data points (at each time step) on the plot; default is true

Level: intermediate

-seealso: [](ch_ts), `TSMonitorLGTimeStep()`, `TSMonitorSet()`, `TSMonitorLGSolution()`, `TSMonitorLGError()`, `TSMonitorDefault()`, `VecView()`,
`TSMonitorLGCtxSetVariableNames()`, `TSMonitorLGCtxGetVariableNames()`,
`TSMonitorLGSetVariableNames()`, `TSMonitorLGGetVariableNames()`, `TSMonitorLGSetDisplayVariables()`, `TSMonitorLGCtxSetDisplayVariables()`,
`TSMonitorLGCtxSetTransform()`, `TSMonitorLGSetTransform()`, `TSMonitorLGSNESIterations()`, `TSMonitorLGKSPIterations()`,
`TSMonitorEnvelopeCtxCreate()`, `TSMonitorEnvelopeGetBounds()`, `TSMonitorEnvelopeCtxDestroy()`, `TSMonitorEnvelop()`

# External Links
$(_doc_external("TS/TSMonitorLGCtxCreate"))
"""
function TSMonitorLGCtxCreate(petsclib::PetscLibType, comm::MPI_Comm, host::String, label::String, x::Cint, y::Cint, m::Cint, n::Cint, howoften::PetscInt) end

@for_petsc function TSMonitorLGCtxCreate(petsclib::$UnionPetscLib, comm::MPI_Comm, host::String, label::String, x::Cint, y::Cint, m::Cint, n::Cint, howoften::$PetscInt )
	ctx_ = Ref{TSMonitorLGCtx}()

    @chk ccall(
               (:TSMonitorLGCtxCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Cint, Cint, Cint, Cint, $PetscInt, Ptr{TSMonitorLGCtx}),
               comm, host, label, x, y, m, n, howoften, ctx_,
              )

	ctx = ctx_[]

	return ctx
end 

"""
	TSMonitorLGCtxDestroy(petsclib::PetscLibType,ctx::Union{TSMonitorLGCtx, Ref{TSMonitorLGCtx}}) 
Destroys a line graph context that was created with `TSMonitorLGCtxCreate()`.

Collective

Input Parameter:
- `ctx` - the monitor context

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSMonitorLGCtxCreate()`, `TSMonitorSet()`, `TSMonitorLGTimeStep()`

# External Links
$(_doc_external("TS/TSMonitorLGCtxDestroy"))
"""
function TSMonitorLGCtxDestroy(petsclib::PetscLibType, ctx::Union{TSMonitorLGCtx, Ref{TSMonitorLGCtx}}) end

@for_petsc function TSMonitorLGCtxDestroy(petsclib::$UnionPetscLib, ctx::Union{TSMonitorLGCtx, Ref{TSMonitorLGCtx}} )
	ctx_ = ctx isa Base.RefValue ? ctx : Ref{TSMonitorLGCtx}(ctx)

    @chk ccall(
               (:TSMonitorLGCtxDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{TSMonitorLGCtx},),
               ctx_,
              )


	return nothing
end 

"""
	ctx::TSMonitorLGCtxNetwork = TSMonitorLGCtxNetworkCreate(petsclib::PetscLibType,ts::AbstractTS, host::String, label::String, x::Cint, y::Cint, m::Cint, n::Cint, howoften::PetscInt) 

# External Links
$(_doc_external("TS/TSMonitorLGCtxNetworkCreate"))
"""
function TSMonitorLGCtxNetworkCreate(petsclib::PetscLibType, ts::AbstractTS, host::String, label::String, x::Cint, y::Cint, m::Cint, n::Cint, howoften::PetscInt) end

@for_petsc function TSMonitorLGCtxNetworkCreate(petsclib::$UnionPetscLib, ts::AbstractTS, host::String, label::String, x::Cint, y::Cint, m::Cint, n::Cint, howoften::$PetscInt )
	ctx_ = Ref{TSMonitorLGCtxNetwork}()

    @chk ccall(
               (:TSMonitorLGCtxNetworkCreate, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{Cchar}, Ptr{Cchar}, Cint, Cint, Cint, Cint, $PetscInt, Ptr{TSMonitorLGCtxNetwork}),
               ts, host, label, x, y, m, n, howoften, ctx_,
              )

	ctx = ctx_[]

	return ctx
end 

"""
	TSMonitorLGCtxNetworkDestroy(petsclib::PetscLibType,ctx::Union{TSMonitorLGCtxNetwork, Ref{TSMonitorLGCtxNetwork}}) 
Destroys  line graph contexts that where created with `TSMonitorLGCtxNetworkCreate()`.

Collective

Input Parameter:
- `ctx` - the monitor context

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSMonitorLGCtxNetworkSolution()`

# External Links
$(_doc_external("TS/TSMonitorLGCtxNetworkDestroy"))
"""
function TSMonitorLGCtxNetworkDestroy(petsclib::PetscLibType, ctx::Union{TSMonitorLGCtxNetwork, Ref{TSMonitorLGCtxNetwork}}) end

@for_petsc function TSMonitorLGCtxNetworkDestroy(petsclib::$UnionPetscLib, ctx::Union{TSMonitorLGCtxNetwork, Ref{TSMonitorLGCtxNetwork}} )
	ctx_ = ctx isa Base.RefValue ? ctx : Ref{TSMonitorLGCtxNetwork}(ctx)

    @chk ccall(
               (:TSMonitorLGCtxNetworkDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{TSMonitorLGCtxNetwork},),
               ctx_,
              )


	return nothing
end 

"""
	TSMonitorLGCtxNetworkSolution(petsclib::PetscLibType,ts::AbstractTS, step::PetscInt, ptime::PetscReal, u::AbstractPetscVec, dctx::Ptr{Cvoid}) 
Monitors progress of the `TS` solvers for a `DMNETWORK` solution with one window for each vertex and each edge

Collective

Input Parameters:
- `ts`    - the `TS` context
- `step`  - current time-step
- `ptime` - current time
- `u`     - current solution
- `dctx`  - the `TSMonitorLGCtxNetwork` object that contains all the options for the monitoring, this is created with `TSMonitorLGCtxCreateNetwork()`

Options Database Key:
- `-ts_monitor_lg_solution_variables` - monitor solution variables

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSMonitorLGCtxNetworkDestroy()`

# External Links
$(_doc_external("TS/TSMonitorLGCtxNetworkSolution"))
"""
function TSMonitorLGCtxNetworkSolution(petsclib::PetscLibType, ts::AbstractTS, step::PetscInt, ptime::PetscReal, u::AbstractPetscVec, dctx::Ptr{Cvoid}) end

@for_petsc function TSMonitorLGCtxNetworkSolution(petsclib::$UnionPetscLib, ts::AbstractTS, step::$PetscInt, ptime::$PetscReal, u::AbstractPetscVec, dctx::Ptr{Cvoid} )

    @chk ccall(
               (:TSMonitorLGCtxNetworkSolution, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, $PetscReal, CVec, Ptr{Cvoid}),
               ts, step, ptime, u, dctx,
              )


	return nothing
end 

"""
	TSMonitorLGCtxSetDisplayVariables(petsclib::PetscLibType,ctx::TSMonitorLGCtx, displaynames::Cchar) 
Sets the variables that are to be display in the monitor

Collective

Input Parameters:
- `ctx`          - the `TSMonitorLG` context
- `displaynames` - the names of the components, final string must be `NULL`

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSMonitorSet()`, `TSMonitorDefault()`, `VecView()`, `TSMonitorLGSetVariableNames()`

# External Links
$(_doc_external("TS/TSMonitorLGCtxSetDisplayVariables"))
"""
function TSMonitorLGCtxSetDisplayVariables(petsclib::PetscLibType, ctx::TSMonitorLGCtx, displaynames::Cchar) end

@for_petsc function TSMonitorLGCtxSetDisplayVariables(petsclib::$UnionPetscLib, ctx::TSMonitorLGCtx, displaynames::Cchar )

    @chk ccall(
               (:TSMonitorLGCtxSetDisplayVariables, $petsc_library),
               PetscErrorCode,
               (TSMonitorLGCtx, Ptr{Ptr{Cchar}}),
               ctx, displaynames,
              )


	return nothing
end 

"""
	TSMonitorLGCtxSetTransform(petsclib::PetscLibType,ctx::TSMonitorLGCtx, transform::external, destroy::Ptr{Cvoid}, tctx::Ptr{Cvoid}) 
Solution vector will be transformed by provided function before being displayed

Collective

Input Parameters:
- `tctx`      - the `TS` context
- `transform` - the transform function
- `destroy`   - function to destroy the optional context, see `PetscCtxDestroyFn` for its calling sequence
- `ctx`       - optional context used by transform function

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSMonitorSet()`, `TSMonitorDefault()`, `VecView()`, `TSMonitorLGSetVariableNames()`, `TSMonitorLGSetTransform()`, `PetscCtxDestroyFn`

# External Links
$(_doc_external("TS/TSMonitorLGCtxSetTransform"))
"""
function TSMonitorLGCtxSetTransform(petsclib::PetscLibType, ctx::TSMonitorLGCtx, transform::external, destroy::Ptr{Cvoid}, tctx::Ptr{Cvoid}) end

@for_petsc function TSMonitorLGCtxSetTransform(petsclib::$UnionPetscLib, ctx::TSMonitorLGCtx, transform::external, destroy::Ptr{Cvoid}, tctx::Ptr{Cvoid} )

    @chk ccall(
               (:TSMonitorLGCtxSetTransform, $petsc_library),
               PetscErrorCode,
               (TSMonitorLGCtx, external, Ptr{Cvoid}, Ptr{Cvoid}),
               ctx, transform, destroy, tctx,
              )


	return nothing
end 

"""
	TSMonitorLGCtxSetVariableNames(petsclib::PetscLibType,ctx::TSMonitorLGCtx, names::Cchar) 
Sets the name of each component in the solution vector so that it may be displayed in the plot

Collective

Input Parameters:
- `ctx`   - the `TS` context
- `names` - the names of the components, final string must be `NULL`

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSMonitorSet()`, `TSMonitorDefault()`, `VecView()`, `TSMonitorLGSetDisplayVariables()`, `TSMonitorLGSetVariableNames()`

# External Links
$(_doc_external("TS/TSMonitorLGCtxSetVariableNames"))
"""
function TSMonitorLGCtxSetVariableNames(petsclib::PetscLibType, ctx::TSMonitorLGCtx, names::Cchar) end

@for_petsc function TSMonitorLGCtxSetVariableNames(petsclib::$UnionPetscLib, ctx::TSMonitorLGCtx, names::Cchar )

    @chk ccall(
               (:TSMonitorLGCtxSetVariableNames, $petsc_library),
               PetscErrorCode,
               (TSMonitorLGCtx, Ptr{Ptr{Cchar}}),
               ctx, names,
              )


	return nothing
end 

"""
	ctx::TSMonitorSPCtx = TSMonitorSPCtxCreate(petsclib::PetscLibType,comm::MPI_Comm, host::String, label::String, x::Cint, y::Cint, m::Cint, n::Cint, howoften::PetscInt, retain::PetscInt, phase::PetscBool, multispecies::PetscBool) 

# External Links
$(_doc_external("TS/TSMonitorSPCtxCreate"))
"""
function TSMonitorSPCtxCreate(petsclib::PetscLibType, comm::MPI_Comm, host::String, label::String, x::Cint, y::Cint, m::Cint, n::Cint, howoften::PetscInt, retain::PetscInt, phase::PetscBool, multispecies::PetscBool) end

@for_petsc function TSMonitorSPCtxCreate(petsclib::$UnionPetscLib, comm::MPI_Comm, host::String, label::String, x::Cint, y::Cint, m::Cint, n::Cint, howoften::$PetscInt, retain::$PetscInt, phase::PetscBool, multispecies::PetscBool )
	ctx_ = Ref{TSMonitorSPCtx}()

    @chk ccall(
               (:TSMonitorSPCtxCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Cint, Cint, Cint, Cint, $PetscInt, $PetscInt, PetscBool, PetscBool, Ptr{TSMonitorSPCtx}),
               comm, host, label, x, y, m, n, howoften, retain, phase, multispecies, ctx_,
              )

	ctx = ctx_[]

	return ctx
end 

"""
	TSMonitorSPCtxDestroy(petsclib::PetscLibType,ctx::Union{TSMonitorSPCtx, Ref{TSMonitorSPCtx}}) 

# External Links
$(_doc_external("TS/TSMonitorSPCtxDestroy"))
"""
function TSMonitorSPCtxDestroy(petsclib::PetscLibType, ctx::Union{TSMonitorSPCtx, Ref{TSMonitorSPCtx}}) end

@for_petsc function TSMonitorSPCtxDestroy(petsclib::$UnionPetscLib, ctx::Union{TSMonitorSPCtx, Ref{TSMonitorSPCtx}} )
	ctx_ = ctx isa Base.RefValue ? ctx : Ref{TSMonitorSPCtx}(ctx)

    @chk ccall(
               (:TSMonitorSPCtxDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{TSMonitorSPCtx},),
               ctx_,
              )


	return nothing
end 

"""
	ctx::TSMonitorSPEigCtx = TSMonitorSPEigCtxCreate(petsclib::PetscLibType,comm::MPI_Comm, host::String, label::String, x::Cint, y::Cint, m::Cint, n::Cint, howoften::PetscInt) 
Creates a context for use with `TS` to monitor the eigenvalues of the linearized operator

Collective

Input Parameters:
- `comm`     - the communicator to share the monitor
- `host`     - the X display to open, or `NULL` for the local machine
- `label`    - the title to put in the title bar
- `x`        - the horizontal screen coordinates of the upper left coordinate of the window
- `y`        - the vertical coordinates of the upper left coordinate of the window
- `m`        - the screen width in pixels
- `n`        - the screen height in pixels
- `howoften` - if positive then determines the frequency of the plotting, if -1 then only at the final time

Output Parameter:
- `ctx` - the context

Options Database Key:
- `-ts_monitor_sp_eig` - plot egienvalues of linearized right-hand side

Level: intermediate

-seealso: [](ch_ts), `TSMonitorSPEigTimeStep()`, `TSMonitorSet()`, `TSMonitorLGSolution()`, `TSMonitorLGError()`

# External Links
$(_doc_external("TS/TSMonitorSPEigCtxCreate"))
"""
function TSMonitorSPEigCtxCreate(petsclib::PetscLibType, comm::MPI_Comm, host::String, label::String, x::Cint, y::Cint, m::Cint, n::Cint, howoften::PetscInt) end

@for_petsc function TSMonitorSPEigCtxCreate(petsclib::$UnionPetscLib, comm::MPI_Comm, host::String, label::String, x::Cint, y::Cint, m::Cint, n::Cint, howoften::$PetscInt )
	ctx_ = Ref{TSMonitorSPEigCtx}()

    @chk ccall(
               (:TSMonitorSPEigCtxCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Cint, Cint, Cint, Cint, $PetscInt, Ptr{TSMonitorSPEigCtx}),
               comm, host, label, x, y, m, n, howoften, ctx_,
              )

	ctx = ctx_[]

	return ctx
end 

"""
	TSMonitorSPEigCtxDestroy(petsclib::PetscLibType,ctx::Union{TSMonitorSPEigCtx, Ref{TSMonitorSPEigCtx}}) 
Destroys a scatter plot context that was created with `TSMonitorSPEigCtxCreate()`.

Collective

Input Parameter:
- `ctx` - the monitor context

Level: intermediate

-seealso: [](ch_ts), `TSMonitorSPEigCtxCreate()`, `TSMonitorSet()`, `TSMonitorSPEig()`

# External Links
$(_doc_external("TS/TSMonitorSPEigCtxDestroy"))
"""
function TSMonitorSPEigCtxDestroy(petsclib::PetscLibType, ctx::Union{TSMonitorSPEigCtx, Ref{TSMonitorSPEigCtx}}) end

@for_petsc function TSMonitorSPEigCtxDestroy(petsclib::$UnionPetscLib, ctx::Union{TSMonitorSPEigCtx, Ref{TSMonitorSPEigCtx}} )
	ctx_ = ctx isa Base.RefValue ? ctx : Ref{TSMonitorSPEigCtx}(ctx)

    @chk ccall(
               (:TSMonitorSPEigCtxDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{TSMonitorSPEigCtx},),
               ctx_,
              )


	return nothing
end 

"""
	tj::TSTrajectory = TSTrajectoryCreate(petsclib::PetscLibType,comm::MPI_Comm) 
This function creates an empty trajectory object used to store the time dependent solution of an ODE/DAE

Collective

Input Parameter:
- `comm` - the communicator

Output Parameter:
- `tj` - the trajectory object

Level: developer

-seealso: [](ch_ts), `TS`, `TSTrajectory`, `TSTrajectorySetUp()`, `TSTrajectoryDestroy()`, `TSTrajectorySetType()`, `TSTrajectorySetVariableNames()`, `TSGetTrajectory()`, `TSTrajectorySetKeepFiles()`

# External Links
$(_doc_external("TS/TSTrajectoryCreate"))
"""
function TSTrajectoryCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function TSTrajectoryCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	tj_ = Ref{TSTrajectory}()

    @chk ccall(
               (:TSTrajectoryCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{TSTrajectory}),
               comm, tj_,
              )

	tj = tj_[]

	return tj
end 

"""
	TSTrajectoryDestroy(petsclib::PetscLibType,tj::Union{TSTrajectory, Ref{TSTrajectory}}) 
Destroys a trajectory context

Collective

Input Parameter:
- `tj` - the `TSTrajectory` context obtained from `TSTrajectoryCreate()`

Level: developer

-seealso: [](ch_ts), `TSTrajectory`, `TSTrajectoryCreate()`, `TSTrajectorySetUp()`

# External Links
$(_doc_external("TS/TSTrajectoryDestroy"))
"""
function TSTrajectoryDestroy(petsclib::PetscLibType, tj::Union{TSTrajectory, Ref{TSTrajectory}}) end

@for_petsc function TSTrajectoryDestroy(petsclib::$UnionPetscLib, tj::Union{TSTrajectory, Ref{TSTrajectory}} )
	tj_ = tj isa Base.RefValue ? tj : Ref{TSTrajectory}(tj)

    @chk ccall(
               (:TSTrajectoryDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{TSTrajectory},),
               tj_,
              )


	return nothing
end 

"""
	time::PetscReal = TSTrajectoryGet(petsclib::PetscLibType,tj::TSTrajectory, ts::AbstractTS, stepnum::PetscInt) 
Updates the solution vector of a time stepper object by querying the `TSTrajectory`

Collective

Input Parameters:
- `tj`      - the trajectory object
- `ts`      - the time stepper object
- `stepnum` - the step number

Output Parameter:
- `time` - the time associated with the step number

Level: developer

-seealso: [](ch_ts), `TS`, `TSSolve()`, `TSTrajectorySetUp()`, `TSTrajectoryDestroy()`, `TSTrajectorySetType()`, `TSTrajectorySetVariableNames()`, `TSGetTrajectory()`, `TSTrajectorySet()`, `TSTrajectoryGetVecs()`, `TSGetSolution()`

# External Links
$(_doc_external("TS/TSTrajectoryGet"))
"""
function TSTrajectoryGet(petsclib::PetscLibType, tj::TSTrajectory, ts::AbstractTS, stepnum::PetscInt) end

@for_petsc function TSTrajectoryGet(petsclib::$UnionPetscLib, tj::TSTrajectory, ts::AbstractTS, stepnum::$PetscInt )
	time_ = Ref{$PetscReal}()

    @chk ccall(
               (:TSTrajectoryGet, $petsc_library),
               PetscErrorCode,
               (TSTrajectory, CTS, $PetscInt, Ptr{$PetscReal}),
               tj, ts, stepnum, time_,
              )

	time = time_[]

	return time
end 

"""
	steps::PetscInt = TSTrajectoryGetNumSteps(petsclib::PetscLibType,tj::TSTrajectory) 
Return the number of steps registered in the `TSTrajectory` via `TSTrajectorySet()`.

Not Collective.

Input Parameter:
- `tj` - the trajectory object

Output Parameter:
- `steps` - the number of steps

Level: developer

-seealso: [](ch_ts), `TS`, `TSTrajectorySet()`

# External Links
$(_doc_external("TS/TSTrajectoryGetNumSteps"))
"""
function TSTrajectoryGetNumSteps(petsclib::PetscLibType, tj::TSTrajectory) end

@for_petsc function TSTrajectoryGetNumSteps(petsclib::$UnionPetscLib, tj::TSTrajectory )
	steps_ = Ref{$PetscInt}()

    @chk ccall(
               (:TSTrajectoryGetNumSteps, $petsc_library),
               PetscErrorCode,
               (TSTrajectory, Ptr{$PetscInt}),
               tj, steps_,
              )

	steps = steps_[]

	return steps
end 

"""
	solution_only::PetscBool = TSTrajectoryGetSolutionOnly(petsclib::PetscLibType,tj::TSTrajectory) 
Gets the value set with `TSTrajectorySetSolutionOnly()`.

Logically Collective

Input Parameter:
- `tj` - the `TSTrajectory` context

Output Parameter:
- `solution_only` - the boolean flag

Level: developer

-seealso: [](ch_ts), `TSTrajectory`, `TSSetSaveTrajectory()`, `TSTrajectoryCreate()`, `TSTrajectoryDestroy()`, `TSTrajectorySetSolutionOnly()`

# External Links
$(_doc_external("TS/TSTrajectoryGetSolutionOnly"))
"""
function TSTrajectoryGetSolutionOnly(petsclib::PetscLibType, tj::TSTrajectory) end

@for_petsc function TSTrajectoryGetSolutionOnly(petsclib::$UnionPetscLib, tj::TSTrajectory )
	solution_only_ = Ref{PetscBool}()

    @chk ccall(
               (:TSTrajectoryGetSolutionOnly, $petsc_library),
               PetscErrorCode,
               (TSTrajectory, Ptr{PetscBool}),
               tj, solution_only_,
              )

	solution_only = solution_only_[]

	return solution_only
end 

"""
	type::TSTrajectoryType = TSTrajectoryGetType(petsclib::PetscLibType,tj::TSTrajectory, ts::AbstractTS) 
Gets the trajectory type

Collective

Input Parameters:
- `tj` - the `TSTrajectory` context
- `ts` - the `TS` context

Output Parameter:
- `type` - a known method

Level: developer

-seealso: [](ch_ts), `TS`, `TSTrajectory`, `TSTrajectoryCreate()`, `TSTrajectorySetFromOptions()`, `TSTrajectoryDestroy()`, `TSTrajectorySetType()`

# External Links
$(_doc_external("TS/TSTrajectoryGetType"))
"""
function TSTrajectoryGetType(petsclib::PetscLibType, tj::TSTrajectory, ts::AbstractTS) end

@for_petsc function TSTrajectoryGetType(petsclib::$UnionPetscLib, tj::TSTrajectory, ts::AbstractTS )
	type_ = Ref{TSTrajectoryType}()

    @chk ccall(
               (:TSTrajectoryGetType, $petsc_library),
               PetscErrorCode,
               (TSTrajectory, CTS, Ptr{TSTrajectoryType}),
               tj, ts, type_,
              )

	type = type_[] == C_NULL ? "" : unsafe_string(type_[])

	return type
end 

"""
	U::PetscVec,Udot::PetscVec = TSTrajectoryGetUpdatedHistoryVecs(petsclib::PetscLibType,tj::TSTrajectory, ts::AbstractTS, time::PetscReal) 
Get updated state and time

Collective

Input Parameters:
- `tj`   - the `TSTrajectory` context
- `ts`   - the `TS` solver context
- `time` - the requested time

Output Parameters:
- `U`    - state vector at given time (can be interpolated)
- `Udot` - time-derivative vector at given time (can be interpolated)

Level: developer

-seealso: [](ch_ts), `TSTrajectory`, `TSSetSaveTrajectory()`, `TSTrajectoryCreate()`, `TSTrajectoryDestroy()`, `TSTrajectoryRestoreUpdatedHistoryVecs()`, `TSTrajectoryGetVecs()`

# External Links
$(_doc_external("TS/TSTrajectoryGetUpdatedHistoryVecs"))
"""
function TSTrajectoryGetUpdatedHistoryVecs(petsclib::PetscLibType, tj::TSTrajectory, ts::AbstractTS, time::PetscReal) end

@for_petsc function TSTrajectoryGetUpdatedHistoryVecs(petsclib::$UnionPetscLib, tj::TSTrajectory, ts::AbstractTS, time::$PetscReal )
	U_ = Ref{CVec}()
	Udot_ = Ref{CVec}()

    @chk ccall(
               (:TSTrajectoryGetUpdatedHistoryVecs, $petsc_library),
               PetscErrorCode,
               (TSTrajectory, CTS, $PetscReal, Ptr{CVec}, Ptr{CVec}),
               tj, ts, time, U_, Udot_,
              )

	U = PetscVec(U_[], petsclib)
	Udot = PetscVec(Udot_[], petsclib)

	return U,Udot
end 

"""
	time::PetscReal = TSTrajectoryGetVecs(petsclib::PetscLibType,tj::TSTrajectory, ts::AbstractTS, stepnum::PetscInt, U::AbstractPetscVec, Udot::AbstractPetscVec) 
Reconstructs the vector of state and its time derivative using information from the `TSTrajectory` and, possibly, from the `TS`

Collective

Input Parameters:
- `tj`      - the trajectory object
- `ts`      - the time stepper object (optional)
- `stepnum` - the requested step number

Output Parameters:
- `time` - On input time for the step if step number is `PETSC_DECIDE`, on output the time associated with the step number
- `U`    - state vector (can be `NULL`)
- `Udot` - time derivative of state vector (can be `NULL`)

Level: developer

-seealso: [](ch_ts), `TS`, `TSTrajectory`, `TSTrajectorySetUp()`, `TSTrajectoryDestroy()`, `TSTrajectorySetType()`, `TSTrajectorySetVariableNames()`, `TSGetTrajectory()`, `TSTrajectorySet()`, `TSTrajectoryGet()`

# External Links
$(_doc_external("TS/TSTrajectoryGetVecs"))
"""
function TSTrajectoryGetVecs(petsclib::PetscLibType, tj::TSTrajectory, ts::AbstractTS, stepnum::PetscInt, U::AbstractPetscVec, Udot::AbstractPetscVec) end

@for_petsc function TSTrajectoryGetVecs(petsclib::$UnionPetscLib, tj::TSTrajectory, ts::AbstractTS, stepnum::$PetscInt, U::AbstractPetscVec, Udot::AbstractPetscVec )
	time_ = Ref{$PetscReal}()

    @chk ccall(
               (:TSTrajectoryGetVecs, $petsc_library),
               PetscErrorCode,
               (TSTrajectory, CTS, $PetscInt, Ptr{$PetscReal}, CVec, CVec),
               tj, ts, stepnum, time_, U, Udot,
              )

	time = time_[]

	return time
end 

"""
	TSTrajectoryMemorySetType(petsclib::PetscLibType,tj::TSTrajectory, tj_memory_type::TSTrajectoryMemoryType) 
sets the software that is used to generate the checkpointing schedule.

Logically Collective

Input Parameters:
- `tj`             - the `TSTrajectory` context
- `tj_memory_type` - Revolve or CAMS

Options Database Key:
- `-ts_trajectory_memory_type <tj_memory_type>` - petsc, revolve, cams

Level: intermediate

-seealso: [](ch_ts), `TSTrajectory`, `TSTrajectorySetMaxUnitsRAM()`, `TSTrajectoryMemoryType`

# External Links
$(_doc_external("TS/TSTrajectoryMemorySetType"))
"""
function TSTrajectoryMemorySetType(petsclib::PetscLibType, tj::TSTrajectory, tj_memory_type::TSTrajectoryMemoryType) end

@for_petsc function TSTrajectoryMemorySetType(petsclib::$UnionPetscLib, tj::TSTrajectory, tj_memory_type::TSTrajectoryMemoryType )

    @chk ccall(
               (:TSTrajectoryMemorySetType, $petsc_library),
               PetscErrorCode,
               (TSTrajectory, TSTrajectoryMemoryType),
               tj, tj_memory_type,
              )


	return nothing
end 

"""
	TSTrajectoryRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Adds a way of storing trajectories to the `TS` package

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - the name of a new user-defined creation routine
- `function` - the creation routine itself

Level: developer

-seealso: [](ch_ts), `TSTrajectoryRegisterAll()`

# External Links
$(_doc_external("TS/TSTrajectoryRegister"))
"""
function TSTrajectoryRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function TSTrajectoryRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:TSTrajectoryRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	TSTrajectoryRegisterAll(petsclib::PetscLibType) 
Registers all of the `TSTrajectory` storage schecmes in the `TS` package.

Not Collective

Level: developer

-seealso: [](ch_ts), `TSTrajectory`, `TSTrajectoryRegister()`

# External Links
$(_doc_external("TS/TSTrajectoryRegisterAll"))
"""
function TSTrajectoryRegisterAll(petsclib::PetscLibType) end

@for_petsc function TSTrajectoryRegisterAll(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSTrajectoryRegisterAll, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSTrajectoryReset(petsclib::PetscLibType,tj::TSTrajectory) 
Resets a trajectory context

Collective

Input Parameter:
- `tj` - the `TSTrajectory` context obtained from `TSGetTrajectory()`

Level: developer

-seealso: [](ch_ts), `TS`, `TSTrajectory`, `TSTrajectoryCreate()`, `TSTrajectorySetUp()`

# External Links
$(_doc_external("TS/TSTrajectoryReset"))
"""
function TSTrajectoryReset(petsclib::PetscLibType, tj::TSTrajectory) end

@for_petsc function TSTrajectoryReset(petsclib::$UnionPetscLib, tj::TSTrajectory )

    @chk ccall(
               (:TSTrajectoryReset, $petsc_library),
               PetscErrorCode,
               (TSTrajectory,),
               tj,
              )


	return nothing
end 

"""
	TSTrajectoryRestoreUpdatedHistoryVecs(petsclib::PetscLibType,tj::TSTrajectory, U::AbstractPetscVec, Udot::AbstractPetscVec) 
Restores updated state and time

Collective

Input Parameters:
- `tj`   - the `TSTrajectory` context
- `U`    - state vector at given time (can be interpolated)
- `Udot` - time-derivative vector at given time (can be interpolated)

Level: developer

-seealso: [](ch_ts), `TSTrajectory`, `TSTrajectoryGetUpdatedHistoryVecs()`

# External Links
$(_doc_external("TS/TSTrajectoryRestoreUpdatedHistoryVecs"))
"""
function TSTrajectoryRestoreUpdatedHistoryVecs(petsclib::PetscLibType, tj::TSTrajectory, U::AbstractPetscVec, Udot::AbstractPetscVec) end

@for_petsc function TSTrajectoryRestoreUpdatedHistoryVecs(petsclib::$UnionPetscLib, tj::TSTrajectory, U::AbstractPetscVec, Udot::AbstractPetscVec )
	U_ = Ref(U.ptr)
	Udot_ = Ref(Udot.ptr)

    @chk ccall(
               (:TSTrajectoryRestoreUpdatedHistoryVecs, $petsc_library),
               PetscErrorCode,
               (TSTrajectory, Ptr{CVec}, Ptr{CVec}),
               tj, U_, Udot_,
              )

	U.ptr = U_[]
	Udot.ptr = Udot_[]

	return nothing
end 

"""
	TSTrajectorySet(petsclib::PetscLibType,tj::TSTrajectory, ts::AbstractTS, stepnum::PetscInt, time::PetscReal, X::AbstractPetscVec) 
Sets a vector of state in the trajectory object

Collective

Input Parameters:
- `tj`      - the trajectory object
- `ts`      - the time stepper object (optional)
- `stepnum` - the step number
- `time`    - the current time
- `X`       - the current solution

Level: developer

-seealso: [](ch_ts), `TSTrajectorySetUp()`, `TSTrajectoryDestroy()`, `TSTrajectorySetType()`, `TSTrajectorySetVariableNames()`, `TSGetTrajectory()`, `TSTrajectoryGet()`, `TSTrajectoryGetVecs()`

# External Links
$(_doc_external("TS/TSTrajectorySet"))
"""
function TSTrajectorySet(petsclib::PetscLibType, tj::TSTrajectory, ts::AbstractTS, stepnum::PetscInt, time::PetscReal, X::AbstractPetscVec) end

@for_petsc function TSTrajectorySet(petsclib::$UnionPetscLib, tj::TSTrajectory, ts::AbstractTS, stepnum::$PetscInt, time::$PetscReal, X::AbstractPetscVec )

    @chk ccall(
               (:TSTrajectorySet, $petsc_library),
               PetscErrorCode,
               (TSTrajectory, CTS, $PetscInt, $PetscReal, CVec),
               tj, ts, stepnum, time, X,
              )


	return nothing
end 

"""
	TSTrajectorySetDirname(petsclib::PetscLibType,tj::TSTrajectory, dirname::String) 
Specify the name of the directory where `TSTrajectory` disk checkpoints are stored.

Collective

Input Parameters:
- `tj`      - the `TSTrajectory` context
- `dirname` - the directory name

Options Database Key:
- `-ts_trajectory_dirname` - set the directory name

Level: developer

-seealso: [](ch_ts), `TSTrajectory`, `TSTrajectorySetFiletemplate()`, `TSTrajectorySetUp()`

# External Links
$(_doc_external("TS/TSTrajectorySetDirname"))
"""
function TSTrajectorySetDirname(petsclib::PetscLibType, tj::TSTrajectory, dirname::String) end

@for_petsc function TSTrajectorySetDirname(petsclib::$UnionPetscLib, tj::TSTrajectory, dirname::String )

    @chk ccall(
               (:TSTrajectorySetDirname, $petsc_library),
               PetscErrorCode,
               (TSTrajectory, Ptr{Cchar}),
               tj, dirname,
              )


	return nothing
end 

"""
	TSTrajectorySetFiletemplate(petsclib::PetscLibType,tj::TSTrajectory, filetemplate::String) 
Specify the name template for the files storing `TSTrajectory` checkpoints.

Collective

Input Parameters:
- `tj`           - the `TSTrajectory` context
- `filetemplate` - the template

Options Database Key:
- `-ts_trajectory_file_template` - set the file name template

Level: developer

-seealso: [](ch_ts), `TSTrajectory`, `TSTrajectorySetDirname()`, `TSTrajectorySetUp()`

# External Links
$(_doc_external("TS/TSTrajectorySetFiletemplate"))
"""
function TSTrajectorySetFiletemplate(petsclib::PetscLibType, tj::TSTrajectory, filetemplate::String) end

@for_petsc function TSTrajectorySetFiletemplate(petsclib::$UnionPetscLib, tj::TSTrajectory, filetemplate::String )

    @chk ccall(
               (:TSTrajectorySetFiletemplate, $petsc_library),
               PetscErrorCode,
               (TSTrajectory, Ptr{Cchar}),
               tj, filetemplate,
              )


	return nothing
end 

"""
	TSTrajectorySetFromOptions(petsclib::PetscLibType,tj::TSTrajectory, ts::AbstractTS) 
Sets various `TSTrajectory` parameters from user options.

Collective

Input Parameters:
- `tj` - the `TSTrajectory` context obtained from `TSGetTrajectory()`
- `ts` - the `TS` context

Options Database Keys:
- `-ts_trajectory_type <type>`             - basic, memory, singlefile, visualization
- `-ts_trajectory_keep_files <true,false>` - keep the files generated by the code after the program ends. This is true by default for singlefile and visualization
- `-ts_trajectory_monitor`                 - print `TSTrajectory` information

Level: developer

-seealso: [](ch_ts), `TSTrajectory`, `TSSetSaveTrajectory()`, `TSTrajectorySetUp()`

# External Links
$(_doc_external("TS/TSTrajectorySetFromOptions"))
"""
function TSTrajectorySetFromOptions(petsclib::PetscLibType, tj::TSTrajectory, ts::AbstractTS) end

@for_petsc function TSTrajectorySetFromOptions(petsclib::$UnionPetscLib, tj::TSTrajectory, ts::AbstractTS )

    @chk ccall(
               (:TSTrajectorySetFromOptions, $petsc_library),
               PetscErrorCode,
               (TSTrajectory, CTS),
               tj, ts,
              )


	return nothing
end 

"""
	TSTrajectorySetKeepFiles(petsclib::PetscLibType,tj::TSTrajectory, flg::PetscBool) 
Keep the files generated by the `TSTrajectory` once the program is done

Collective

Input Parameters:
- `tj`  - the `TSTrajectory` context
- `flg` - `PETSC_TRUE` to save, `PETSC_FALSE` to disable

Options Database Key:
- `-ts_trajectory_keep_files` - have it keep the files

Level: advanced

-seealso: [](ch_ts), `TSTrajectoryCreate()`, `TSTrajectoryDestroy()`, `TSTrajectorySetUp()`, `TSTrajectorySetMonitor()`

# External Links
$(_doc_external("TS/TSTrajectorySetKeepFiles"))
"""
function TSTrajectorySetKeepFiles(petsclib::PetscLibType, tj::TSTrajectory, flg::PetscBool) end

@for_petsc function TSTrajectorySetKeepFiles(petsclib::$UnionPetscLib, tj::TSTrajectory, flg::PetscBool )

    @chk ccall(
               (:TSTrajectorySetKeepFiles, $petsc_library),
               PetscErrorCode,
               (TSTrajectory, PetscBool),
               tj, flg,
              )


	return nothing
end 

"""
	TSTrajectorySetMaxCpsDisk(petsclib::PetscLibType,tj::TSTrajectory, max_cps_disk::PetscInt) 
Set maximum number of checkpoints on disk

Logically Collective

Input Parameter:
- `tj` - tstrajectory context

Output Parameter:
- `max_cps_disk` - maximum number of checkpoints on disk

Level: intermediate

-seealso: [](ch_ts), `TSTrajectory`, `TSTrajectorySetMaxUnitsDisk()`, `TSTrajectorySetMaxUnitsRAM()`

# External Links
$(_doc_external("TS/TSTrajectorySetMaxCpsDisk"))
"""
function TSTrajectorySetMaxCpsDisk(petsclib::PetscLibType, tj::TSTrajectory, max_cps_disk::PetscInt) end

@for_petsc function TSTrajectorySetMaxCpsDisk(petsclib::$UnionPetscLib, tj::TSTrajectory, max_cps_disk::$PetscInt )

    @chk ccall(
               (:TSTrajectorySetMaxCpsDisk, $petsc_library),
               PetscErrorCode,
               (TSTrajectory, $PetscInt),
               tj, max_cps_disk,
              )


	return nothing
end 

"""
	TSTrajectorySetMaxCpsRAM(petsclib::PetscLibType,tj::TSTrajectory, max_cps_ram::PetscInt) 
Set maximum number of checkpoints in RAM

Logically Collective

Input Parameter:
- `tj` - tstrajectory context

Output Parameter:
- `max_cps_ram` - maximum number of checkpoints in RAM

Level: intermediate

-seealso: [](ch_ts), `TSTrajectory`, `TSTrajectorySetMaxUnitsRAM()`

# External Links
$(_doc_external("TS/TSTrajectorySetMaxCpsRAM"))
"""
function TSTrajectorySetMaxCpsRAM(petsclib::PetscLibType, tj::TSTrajectory, max_cps_ram::PetscInt) end

@for_petsc function TSTrajectorySetMaxCpsRAM(petsclib::$UnionPetscLib, tj::TSTrajectory, max_cps_ram::$PetscInt )

    @chk ccall(
               (:TSTrajectorySetMaxCpsRAM, $petsc_library),
               PetscErrorCode,
               (TSTrajectory, $PetscInt),
               tj, max_cps_ram,
              )


	return nothing
end 

"""
	TSTrajectorySetMaxUnitsDisk(petsclib::PetscLibType,tj::TSTrajectory, max_units_disk::PetscInt) 
Set maximum number of checkpointing units on disk

Logically Collective

Input Parameter:
- `tj` - tstrajectory context

Output Parameter:
- `max_units_disk` - maximum number of checkpointing units on disk

Level: intermediate

-seealso: [](ch_ts), `TSTrajectory`, `TSTrajectorySetMaxCpsDisk()`

# External Links
$(_doc_external("TS/TSTrajectorySetMaxUnitsDisk"))
"""
function TSTrajectorySetMaxUnitsDisk(petsclib::PetscLibType, tj::TSTrajectory, max_units_disk::PetscInt) end

@for_petsc function TSTrajectorySetMaxUnitsDisk(petsclib::$UnionPetscLib, tj::TSTrajectory, max_units_disk::$PetscInt )

    @chk ccall(
               (:TSTrajectorySetMaxUnitsDisk, $petsc_library),
               PetscErrorCode,
               (TSTrajectory, $PetscInt),
               tj, max_units_disk,
              )


	return nothing
end 

"""
	TSTrajectorySetMaxUnitsRAM(petsclib::PetscLibType,tj::TSTrajectory, max_units_ram::PetscInt) 
Set maximum number of checkpointing units in RAM

Logically Collective

Input Parameter:
- `tj` - tstrajectory context

Output Parameter:
- `max_units_ram` - maximum number of checkpointing units in RAM

Level: intermediate

-seealso: [](ch_ts), `TSTrajectory`, `TSTrajectorySetMaxCpsRAM()`

# External Links
$(_doc_external("TS/TSTrajectorySetMaxUnitsRAM"))
"""
function TSTrajectorySetMaxUnitsRAM(petsclib::PetscLibType, tj::TSTrajectory, max_units_ram::PetscInt) end

@for_petsc function TSTrajectorySetMaxUnitsRAM(petsclib::$UnionPetscLib, tj::TSTrajectory, max_units_ram::$PetscInt )

    @chk ccall(
               (:TSTrajectorySetMaxUnitsRAM, $petsc_library),
               PetscErrorCode,
               (TSTrajectory, $PetscInt),
               tj, max_units_ram,
              )


	return nothing
end 

"""
	TSTrajectorySetMonitor(petsclib::PetscLibType,tj::TSTrajectory, flg::PetscBool) 
Monitor the schedules generated by the `TSTrajectory` checkpointing controller

Collective

Input Parameters:
- `tj`  - the `TSTrajectory` context
- `flg` - `PETSC_TRUE` to active a monitor, `PETSC_FALSE` to disable

Options Database Key:
- `-ts_trajectory_monitor` - print `TSTrajectory` information

Level: developer

-seealso: [](ch_ts), `TSTrajectory`, `TSTrajectoryCreate()`, `TSTrajectoryDestroy()`, `TSTrajectorySetUp()`

# External Links
$(_doc_external("TS/TSTrajectorySetMonitor"))
"""
function TSTrajectorySetMonitor(petsclib::PetscLibType, tj::TSTrajectory, flg::PetscBool) end

@for_petsc function TSTrajectorySetMonitor(petsclib::$UnionPetscLib, tj::TSTrajectory, flg::PetscBool )

    @chk ccall(
               (:TSTrajectorySetMonitor, $petsc_library),
               PetscErrorCode,
               (TSTrajectory, PetscBool),
               tj, flg,
              )


	return nothing
end 

"""
	TSTrajectorySetSolutionOnly(petsclib::PetscLibType,tj::TSTrajectory, solution_only::PetscBool) 
Tells the trajectory to store just the solution, and not any intermediate stage information

Collective

Input Parameters:
- `tj`            - the `TSTrajectory` context obtained with `TSGetTrajectory()`
- `solution_only` - the boolean flag

Level: developer

-seealso: [](ch_ts), `TSTrajectory`, `TSSetSaveTrajectory()`, `TSTrajectoryCreate()`, `TSTrajectoryDestroy()`, `TSTrajectoryGetSolutionOnly()`

# External Links
$(_doc_external("TS/TSTrajectorySetSolutionOnly"))
"""
function TSTrajectorySetSolutionOnly(petsclib::PetscLibType, tj::TSTrajectory, solution_only::PetscBool) end

@for_petsc function TSTrajectorySetSolutionOnly(petsclib::$UnionPetscLib, tj::TSTrajectory, solution_only::PetscBool )

    @chk ccall(
               (:TSTrajectorySetSolutionOnly, $petsc_library),
               PetscErrorCode,
               (TSTrajectory, PetscBool),
               tj, solution_only,
              )


	return nothing
end 

"""
	TSTrajectorySetTransform(petsclib::PetscLibType,tj::TSTrajectory, transform::external, destroy::external, tctx::Ptr{Cvoid}) 
Solution vector will be transformed by provided function before being saved to disk

Collective

Input Parameters:
- `tj`        - the `TSTrajectory` context
- `transform` - the transform function
- `destroy`   - function to destroy the optional context
- `tctx`      - optional context used by transform function

Level: intermediate

-seealso: [](ch_ts), `TSTrajectorySetVariableNames()`, `TSTrajectory`, `TSMonitorLGSetTransform()`

# External Links
$(_doc_external("TS/TSTrajectorySetTransform"))
"""
function TSTrajectorySetTransform(petsclib::PetscLibType, tj::TSTrajectory, transform::external, destroy::external, tctx::Ptr{Cvoid}) end

@for_petsc function TSTrajectorySetTransform(petsclib::$UnionPetscLib, tj::TSTrajectory, transform::external, destroy::external, tctx::Ptr{Cvoid} )

    @chk ccall(
               (:TSTrajectorySetTransform, $petsc_library),
               PetscErrorCode,
               (TSTrajectory, external, external, Ptr{Cvoid}),
               tj, transform, destroy, tctx,
              )


	return nothing
end 

"""
	TSTrajectorySetType(petsclib::PetscLibType,tj::TSTrajectory, ts::AbstractTS, type::TSTrajectoryType) 
Sets the storage method to be used as in a trajectory

Collective

Input Parameters:
- `tj`   - the `TSTrajectory` context
- `ts`   - the `TS` context
- `type` - a known method

Options Database Key:
- `-ts_trajectory_type <type>` - Sets the method; use -help for a list of available methods (for instance, basic)

Level: developer

-seealso: [](ch_ts), `TSTrajectory`, `TS`, `TSTrajectoryCreate()`, `TSTrajectorySetFromOptions()`, `TSTrajectoryDestroy()`, `TSTrajectoryGetType()`

# External Links
$(_doc_external("TS/TSTrajectorySetType"))
"""
function TSTrajectorySetType(petsclib::PetscLibType, tj::TSTrajectory, ts::AbstractTS, type::TSTrajectoryType) end

@for_petsc function TSTrajectorySetType(petsclib::$UnionPetscLib, tj::TSTrajectory, ts::AbstractTS, type::TSTrajectoryType )

    @chk ccall(
               (:TSTrajectorySetType, $petsc_library),
               PetscErrorCode,
               (TSTrajectory, CTS, TSTrajectoryType),
               tj, ts, type,
              )


	return nothing
end 

"""
	TSTrajectorySetUp(petsclib::PetscLibType,tj::TSTrajectory, ts::AbstractTS) 
Sets up the internal data structures, e.g. stacks, for the later use
of a `TS` `TSTrajectory`.

Collective

Input Parameters:
- `tj` - the `TSTrajectory` context
- `ts` - the TS context obtained from `TSCreate()`

Level: developer

-seealso: [](ch_ts), `TSTrajectory`, `TSSetSaveTrajectory()`, `TSTrajectoryCreate()`, `TSTrajectoryDestroy()`

# External Links
$(_doc_external("TS/TSTrajectorySetUp"))
"""
function TSTrajectorySetUp(petsclib::PetscLibType, tj::TSTrajectory, ts::AbstractTS) end

@for_petsc function TSTrajectorySetUp(petsclib::$UnionPetscLib, tj::TSTrajectory, ts::AbstractTS )

    @chk ccall(
               (:TSTrajectorySetUp, $petsc_library),
               PetscErrorCode,
               (TSTrajectory, CTS),
               tj, ts,
              )


	return nothing
end 

"""
	TSTrajectorySetUseHistory(petsclib::PetscLibType,tj::TSTrajectory, flg::PetscBool) 
Use `TSHistory` in `TSTrajectory`

Collective

Input Parameters:
- `tj`  - the `TSTrajectory` context
- `flg` - `PETSC_TRUE` to save, `PETSC_FALSE` to disable

Options Database Key:
- `-ts_trajectory_use_history` - have it use `TSHistory`

Level: advanced

-seealso: [](ch_ts), `TSTrajectory`, `TSTrajectoryCreate()`, `TSTrajectoryDestroy()`, `TSTrajectorySetUp()`

# External Links
$(_doc_external("TS/TSTrajectorySetUseHistory"))
"""
function TSTrajectorySetUseHistory(petsclib::PetscLibType, tj::TSTrajectory, flg::PetscBool) end

@for_petsc function TSTrajectorySetUseHistory(petsclib::$UnionPetscLib, tj::TSTrajectory, flg::PetscBool )

    @chk ccall(
               (:TSTrajectorySetUseHistory, $petsc_library),
               PetscErrorCode,
               (TSTrajectory, PetscBool),
               tj, flg,
              )


	return nothing
end 

"""
	TSTrajectorySetVariableNames(petsclib::PetscLibType,ctx::TSTrajectory, names::Cchar) 
Sets the name of each component in the solution vector so that it may be saved with the trajectory

Collective

Input Parameters:
- `ctx`   - the trajectory context
- `names` - the names of the components, final string must be `NULL`

Level: intermediate

-seealso: [](ch_ts), `TSTrajectory`, `TSGetTrajectory()`

# External Links
$(_doc_external("TS/TSTrajectorySetVariableNames"))
"""
function TSTrajectorySetVariableNames(petsclib::PetscLibType, ctx::TSTrajectory, names::Cchar) end

@for_petsc function TSTrajectorySetVariableNames(petsclib::$UnionPetscLib, ctx::TSTrajectory, names::Cchar )

    @chk ccall(
               (:TSTrajectorySetVariableNames, $petsc_library),
               PetscErrorCode,
               (TSTrajectory, Ptr{Ptr{Cchar}}),
               ctx, names,
              )


	return nothing
end 

"""
	TSTrajectoryView(petsclib::PetscLibType,tj::TSTrajectory, viewer::PetscViewer) 
Prints information about the trajectory object

Collective

Input Parameters:
- `tj`     - the `TSTrajectory` context obtained from `TSTrajectoryCreate()`
- `viewer` - visualization context

Options Database Key:
- `-ts_trajectory_view` - calls `TSTrajectoryView()` at end of `TSAdjointStep()`

Level: developer

-seealso: [](ch_ts), `TS`, `TSTrajectory`, `PetscViewer`, `PetscViewerASCIIOpen()`

# External Links
$(_doc_external("TS/TSTrajectoryView"))
"""
function TSTrajectoryView(petsclib::PetscLibType, tj::TSTrajectory, viewer::PetscViewer) end

@for_petsc function TSTrajectoryView(petsclib::$UnionPetscLib, tj::TSTrajectory, viewer::PetscViewer )

    @chk ccall(
               (:TSTrajectoryView, $petsc_library),
               PetscErrorCode,
               (TSTrajectory, PetscViewer),
               tj, viewer,
              )


	return nothing
end 

"""
	TSTrajectoryViewFromOptions(petsclib::PetscLibType,A::TSTrajectory, obj, name::String) 
View a `TSTrajectory` based on values in the options database

Collective

Input Parameters:
- `A`    - the `TSTrajectory` context
- `obj`  - Optional object that provides prefix used for option name
- `name` - command line option

Level: intermediate

-seealso: [](ch_ts), `TSTrajectory`, `TSTrajectoryView`, `PetscObjectViewFromOptions()`, `TSTrajectoryCreate()`

# External Links
$(_doc_external("TS/TSTrajectoryViewFromOptions"))
"""
function TSTrajectoryViewFromOptions(petsclib::PetscLibType, A::TSTrajectory, obj, name::String) end

@for_petsc function TSTrajectoryViewFromOptions(petsclib::$UnionPetscLib, A::TSTrajectory, obj, name::String )

    @chk ccall(
               (:TSTrajectoryViewFromOptions, $petsc_library),
               PetscErrorCode,
               (TSTrajectory, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

