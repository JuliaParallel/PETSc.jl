"""
	TaoLineSearchAppendOptionsPrefix(petsclib::PetscLibType, ls::TaoLineSearch, p::String) 
Appends to the prefix used for searching
for all `TaoLineSearch` options in the database.

Collective

Input Parameters:
- `ls` - the `TaoLineSearch` solver context
- `p`  - the prefix string to prepend to all line search requests

Level: advanced

See also: `Tao`, `TaoLineSearch`, `TaoLineSearchSetOptionsPrefix()`, `TaoLineSearchGetOptionsPrefix()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchAppendOptionsPrefix"))
"""
function TaoLineSearchAppendOptionsPrefix(petsclib::PetscLibType, ls::TaoLineSearch, p::String)
    error("TaoLineSearchAppendOptionsPrefix: no generated method for these argument types")
end

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
	f::PetscReal,steplength::PetscReal,reason::TaoLineSearchConvergedReason = TaoLineSearchApply(petsclib::PetscLibType, ls::TaoLineSearch, x::AbstractPetscVec, g::AbstractPetscVec, s::AbstractPetscVec) 
Performs a line-search in a given step direction.
Criteria for acceptable step length depends on the line-search algorithm chosen

Collective

Input Parameters:
- `ls` - the `TaoLineSearch` context
- `s`  - search direction

Output Parameters:
- `x`          - On input the current solution, on output `x` contains the new solution determined by the line search
- `f`          - On input the objective function value at current solution, on output contains the objective function value at new solution
- `g`          - On input the gradient evaluated at `x`, on output contains the gradient at new solution
- `steplength` - scalar multiplier of `s` used ( x = x_0 + steplength * x)
- `reason`     - `TaoLineSearchConvergedReason` reason why the line-search stopped

Level: advanced

See also: `Tao`, `TaoLineSearchConvergedReason`, `TaoLineSearch`, `TaoLineSearchCreate()`, `TaoLineSearchSetType()`,
`TaoLineSearchSetInitialStepLength()`, `TaoAddLineSearchCounts()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchApply"))
"""
function TaoLineSearchApply(petsclib::PetscLibType, ls::TaoLineSearch, x::AbstractPetscVec, g::AbstractPetscVec, s::AbstractPetscVec)
    error("TaoLineSearchApply: no generated method for these argument types")
end

@for_petsc function TaoLineSearchApply(petsclib::$UnionPetscLib, ls::TaoLineSearch, x::AbstractPetscVec, g::AbstractPetscVec, s::AbstractPetscVec )
	f_ = Ref{$PetscReal}()
	steplength_ = Ref{$PetscReal}()
	reason_ = Ref{TaoLineSearchConvergedReason}()

    @chk ccall(
               (:TaoLineSearchApply, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, CVec, Ptr{$PetscReal}, CVec, CVec, Ptr{$PetscReal}, Ptr{TaoLineSearchConvergedReason}),
               ls, x, f_, g, s, steplength_, reason_,
              )

	f = f_[]
	steplength = steplength_[]
	reason = reason_[]

	return f,steplength,reason
end 

"""
	TaoLineSearchComputeGradient(petsclib::PetscLibType, ls::TaoLineSearch, x::AbstractPetscVec, g::AbstractPetscVec) 
Computes the gradient of the objective function

Collective

Input Parameters:
- `ls` - the `TaoLineSearch` context
- `x`  - input vector

Output Parameter:
- `g` - gradient vector

Level: developer

See also: `Tao`, `TaoLineSearch`, `TaoLineSearchComputeObjective()`, `TaoLineSearchComputeObjectiveAndGradient()`, `TaoLineSearchSetGradient()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchComputeGradient"))
"""
function TaoLineSearchComputeGradient(petsclib::PetscLibType, ls::TaoLineSearch, x::AbstractPetscVec, g::AbstractPetscVec)
    error("TaoLineSearchComputeGradient: no generated method for these argument types")
end

@for_petsc function TaoLineSearchComputeGradient(petsclib::$UnionPetscLib, ls::TaoLineSearch, x::AbstractPetscVec, g::AbstractPetscVec )

    @chk ccall(
               (:TaoLineSearchComputeGradient, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, CVec, CVec),
               ls, x, g,
              )


	return nothing
end 

"""
	f::PetscReal = TaoLineSearchComputeObjective(petsclib::PetscLibType, ls::TaoLineSearch, x::AbstractPetscVec) 
Computes the objective function value at a given point

Collective

Input Parameters:
- `ls` - the `TaoLineSearch` context
- `x`  - input vector

Output Parameter:
- `f` - Objective value at `x`

Level: developer

See also: `Tao`, `TaoLineSearch`, `TaoLineSearchComputeGradient()`, `TaoLineSearchComputeObjectiveAndGradient()`, `TaoLineSearchSetObjectiveRoutine()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchComputeObjective"))
"""
function TaoLineSearchComputeObjective(petsclib::PetscLibType, ls::TaoLineSearch, x::AbstractPetscVec)
    error("TaoLineSearchComputeObjective: no generated method for these argument types")
end

@for_petsc function TaoLineSearchComputeObjective(petsclib::$UnionPetscLib, ls::TaoLineSearch, x::AbstractPetscVec )
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
	f::PetscReal,gts::PetscReal = TaoLineSearchComputeObjectiveAndGTS(petsclib::PetscLibType, ls::TaoLineSearch, x::AbstractPetscVec) 
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

See also: `Tao`, `TaoLineSearch`, `TaoLineSearchComputeGradient()`, `TaoLineSearchComputeObjectiveAndGradient()`, `TaoLineSearchSetObjectiveRoutine()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchComputeObjectiveAndGTS"))
"""
function TaoLineSearchComputeObjectiveAndGTS(petsclib::PetscLibType, ls::TaoLineSearch, x::AbstractPetscVec)
    error("TaoLineSearchComputeObjectiveAndGTS: no generated method for these argument types")
end

@for_petsc function TaoLineSearchComputeObjectiveAndGTS(petsclib::$UnionPetscLib, ls::TaoLineSearch, x::AbstractPetscVec )
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
	f::PetscReal = TaoLineSearchComputeObjectiveAndGradient(petsclib::PetscLibType, ls::TaoLineSearch, x::AbstractPetscVec, g::AbstractPetscVec) 
Computes the objective function value at a given point

Collective

Input Parameters:
- `ls` - the `TaoLineSearch` context
- `x`  - input vector

Output Parameters:
- `f` - Objective value at `x`
- `g` - Gradient vector at `x`

Level: developer

See also: `Tao`, `TaoLineSearch`, `TaoLineSearchComputeGradient()`, `TaoLineSearchSetObjectiveRoutine()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchComputeObjectiveAndGradient"))
"""
function TaoLineSearchComputeObjectiveAndGradient(petsclib::PetscLibType, ls::TaoLineSearch, x::AbstractPetscVec, g::AbstractPetscVec)
    error("TaoLineSearchComputeObjectiveAndGradient: no generated method for these argument types")
end

@for_petsc function TaoLineSearchComputeObjectiveAndGradient(petsclib::$UnionPetscLib, ls::TaoLineSearch, x::AbstractPetscVec, g::AbstractPetscVec )
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
	newls::TaoLineSearch = TaoLineSearchCreate(petsclib::PetscLibType, comm::MPI_Comm) 
Creates a `TaoLineSearch` object.  Algorithms in `Tao` that use
line-searches will automatically create one so this all is rarely needed

Collective

Input Parameter:
- `comm` - MPI communicator

Output Parameter:
- `newls` - the new `TaoLineSearch` context

Options Database Key:
- `-tao_ls_type (unit|more-thuente|gpcg|armijo|owarmijo|ipm)` - select which line search `Tao` should use

Level: developer

See also: `Tao`, `TaoLineSearch`, `TaoLineSearchType`, `TaoLineSearchSetType()`, `TaoLineSearchApply()`, `TaoLineSearchDestroy()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchCreate"))
"""
function TaoLineSearchCreate(petsclib::PetscLibType, comm::MPI_Comm)
    error("TaoLineSearchCreate: no generated method for these argument types")
end

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
	TaoLineSearchDestroy(petsclib::PetscLibType, ls::Union{TaoLineSearch, Ref{TaoLineSearch}}) 
Destroys the `TaoLineSearch` context that was created with
`TaoLineSearchCreate()`

Collective

Input Parameter:
- `ls` - the `TaoLineSearch` context

Level: developer

See also: `TaoLineSearch`, `TaoLineSearchCreate()`, `TaoLineSearchApple()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchDestroy"))
"""
function TaoLineSearchDestroy(petsclib::PetscLibType, ls::Union{TaoLineSearch, Ref{TaoLineSearch}})
    error("TaoLineSearchDestroy: no generated method for these argument types")
end

@for_petsc function TaoLineSearchDestroy(petsclib::$UnionPetscLib, ls::Union{TaoLineSearch, Ref{TaoLineSearch}} )
	ls_ = ls isa Base.RefValue ? ls : Ref{TaoLineSearch}(ls)

    @chk ccall(
               (:TaoLineSearchDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{TaoLineSearch},),
               ls_,
              )


	return nothing
end 

"""
	TaoLineSearchFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the `TaoLineSearch` package. It is called from `PetscFinalize()`.

Level: developer

See also: `Tao`, `TaoLineSearch`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchFinalizePackage"))
"""
function TaoLineSearchFinalizePackage(petsclib::PetscLibType)
    error("TaoLineSearchFinalizePackage: no generated method for these argument types")
end

@for_petsc function TaoLineSearchFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TaoLineSearchFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	f_fullstep::PetscReal = TaoLineSearchGetFullStepObjective(petsclib::PetscLibType, ls::TaoLineSearch) 
Returns the objective function value at the full step.  Useful for some minimization algorithms.

Not Collective

Input Parameter:
- `ls` - the `TaoLineSearch` context

Output Parameter:
- `f_fullstep` - the objective value at the full step length

Level: developer

See also: `TaoLineSearchGetSolution()`, `TaoLineSearchGetStartingVector()`, `TaoLineSearchGetStepDirection()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchGetFullStepObjective"))
"""
function TaoLineSearchGetFullStepObjective(petsclib::PetscLibType, ls::TaoLineSearch)
    error("TaoLineSearchGetFullStepObjective: no generated method for these argument types")
end

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
	nfeval::PetscInt,ngeval::PetscInt,nfgeval::PetscInt = TaoLineSearchGetNumberFunctionEvaluations(petsclib::PetscLibType, ls::TaoLineSearch) 
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

See also: `TaoLineSearch`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchGetNumberFunctionEvaluations"))
"""
function TaoLineSearchGetNumberFunctionEvaluations(petsclib::PetscLibType, ls::TaoLineSearch)
    error("TaoLineSearchGetNumberFunctionEvaluations: no generated method for these argument types")
end

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
	p::Ptr{Cchar} = TaoLineSearchGetOptionsPrefix(petsclib::PetscLibType, ls::TaoLineSearch) 
Gets the prefix used for searching for all
`TaoLineSearch` options in the database

Not Collective

Input Parameter:
- `ls` - the `TaoLineSearch` context

Output Parameter:
- `p` - pointer to the prefix string used is returned

Level: advanced

See also: `Tao`, `TaoLineSearch`, `TaoLineSearchSetOptionsPrefix()`, `TaoLineSearchAppendOptionsPrefix()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchGetOptionsPrefix"))
"""
function TaoLineSearchGetOptionsPrefix(petsclib::PetscLibType, ls::TaoLineSearch)
    error("TaoLineSearchGetOptionsPrefix: no generated method for these argument types")
end

@for_petsc function TaoLineSearchGetOptionsPrefix(petsclib::$UnionPetscLib, ls::TaoLineSearch )
	p_ = Ref{Ptr{Cchar}}()

    @chk ccall(
               (:TaoLineSearchGetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, Ptr{Ptr{Cchar}}),
               ls, p_,
              )

	p = p_[]

	return p
end 

"""
	f::PetscReal,steplength::PetscReal,reason::TaoLineSearchConvergedReason = TaoLineSearchGetSolution(petsclib::PetscLibType, ls::TaoLineSearch, x::AbstractPetscVec, g::AbstractPetscVec) 
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

See also: `TaoLineSearchGetStartingVector()`, `TaoLineSearchGetStepDirection()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchGetSolution"))
"""
function TaoLineSearchGetSolution(petsclib::PetscLibType, ls::TaoLineSearch, x::AbstractPetscVec, g::AbstractPetscVec)
    error("TaoLineSearchGetSolution: no generated method for these argument types")
end

@for_petsc function TaoLineSearchGetSolution(petsclib::$UnionPetscLib, ls::TaoLineSearch, x::AbstractPetscVec, g::AbstractPetscVec )
	f_ = Ref{$PetscReal}()
	steplength_ = Ref{$PetscReal}()
	reason_ = Ref{TaoLineSearchConvergedReason}()

    @chk ccall(
               (:TaoLineSearchGetSolution, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, CVec, Ptr{$PetscReal}, CVec, Ptr{$PetscReal}, Ptr{TaoLineSearchConvergedReason}),
               ls, x, f_, g, steplength_, reason_,
              )

	f = f_[]
	steplength = steplength_[]
	reason = reason_[]

	return f,steplength,reason
end 

"""
	x::PetscVec = TaoLineSearchGetStartingVector(petsclib::PetscLibType, ls::TaoLineSearch) 
Gets a the initial point of the line
search.

Not Collective

Input Parameter:
- `ls` - the `TaoLineSearch` context

Output Parameter:
- `x` - The initial point of the line search

Level: advanced

See also: `TaoLineSearchGetSolution()`, `TaoLineSearchGetStepDirection()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchGetStartingVector"))
"""
function TaoLineSearchGetStartingVector(petsclib::PetscLibType, ls::TaoLineSearch)
    error("TaoLineSearchGetStartingVector: no generated method for these argument types")
end

@for_petsc function TaoLineSearchGetStartingVector(petsclib::$UnionPetscLib, ls::TaoLineSearch )
	x_ = Ref{CVec}()

    @chk ccall(
               (:TaoLineSearchGetStartingVector, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, Ptr{CVec}),
               ls, x_,
              )

	x = PetscVec(x_[], petsclib)

	return x
end 

"""
	s::PetscVec = TaoLineSearchGetStepDirection(petsclib::PetscLibType, ls::TaoLineSearch) 
Gets the step direction of the line
search.

Not Collective

Input Parameter:
- `ls` - the `TaoLineSearch` context

Output Parameter:
- `s` - the step direction of the line search

Level: advanced

See also: `TaoLineSearchGetSolution()`, `TaoLineSearchGetStartingVector()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchGetStepDirection"))
"""
function TaoLineSearchGetStepDirection(petsclib::PetscLibType, ls::TaoLineSearch)
    error("TaoLineSearchGetStepDirection: no generated method for these argument types")
end

@for_petsc function TaoLineSearchGetStepDirection(petsclib::$UnionPetscLib, ls::TaoLineSearch )
	s_ = Ref{CVec}()

    @chk ccall(
               (:TaoLineSearchGetStepDirection, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, Ptr{CVec}),
               ls, s_,
              )

	s = PetscVec(s_[], petsclib)

	return s
end 

"""
	s::PetscReal = TaoLineSearchGetStepLength(petsclib::PetscLibType, ls::TaoLineSearch) 
Get the current step length

Not Collective

Input Parameter:
- `ls` - the `TaoLineSearch` context

Output Parameter:
- `s` - the current step length

Level: intermediate

See also: `Tao`, `TaoLineSearch`, `TaoLineSearchSetInitialStepLength()`, `TaoLineSearchApply()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchGetStepLength"))
"""
function TaoLineSearchGetStepLength(petsclib::PetscLibType, ls::TaoLineSearch)
    error("TaoLineSearchGetStepLength: no generated method for these argument types")
end

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
	type::TaoLineSearchType = TaoLineSearchGetType(petsclib::PetscLibType, ls::TaoLineSearch) 
Gets the current line search algorithm

Not Collective

Input Parameter:
- `ls` - the `TaoLineSearch` context

Output Parameter:
- `type` - the line search algorithm in effect

Level: developer

See also: `TaoLineSearch`, `TaoLineSearchSetType()`, `TaoLineSearchType`, `PetscObjectTypeCompare()`, `PetscObjectTypeCompareAny()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchGetType"))
"""
function TaoLineSearchGetType(petsclib::PetscLibType, ls::TaoLineSearch)
    error("TaoLineSearchGetType: no generated method for these argument types")
end

@for_petsc function TaoLineSearchGetType(petsclib::$UnionPetscLib, ls::TaoLineSearch )
	type_ = Ref{TaoLineSearchType}()

    @chk ccall(
               (:TaoLineSearchGetType, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, Ptr{TaoLineSearchType}),
               ls, type_,
              )

	type = type_[] == C_NULL ? "" : unsafe_string(type_[])

	return type
end 

"""
	TaoLineSearchInitializePackage(petsclib::PetscLibType) 
This function registers the line-search
algorithms in `Tao`.  When using shared or static libraries, this function is called from the
first entry to `TaoCreate()`; when using dynamic, it is called
from PetscDLLibraryRegister_tao()

Level: developer

See also: `Tao`, `TaoLineSearch`, `TaoLineSearchCreate()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchInitializePackage"))
"""
function TaoLineSearchInitializePackage(petsclib::PetscLibType)
    error("TaoLineSearchInitializePackage: no generated method for these argument types")
end

@for_petsc function TaoLineSearchInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TaoLineSearchInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	flg::PetscBool = TaoLineSearchIsUsingTaoRoutines(petsclib::PetscLibType, ls::TaoLineSearch) 
Checks whether the line search is using
the standard `Tao` evaluation routines.

Not Collective

Input Parameter:
- `ls` - the `TaoLineSearch` context

Output Parameter:
- `flg` - `PETSC_TRUE` if the line search is using `Tao` evaluation routines,
otherwise `PETSC_FALSE`

Level: developer

See also: `TaoLineSearch`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchIsUsingTaoRoutines"))
"""
function TaoLineSearchIsUsingTaoRoutines(petsclib::PetscLibType, ls::TaoLineSearch)
    error("TaoLineSearchIsUsingTaoRoutines: no generated method for these argument types")
end

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
	TaoLineSearchMonitor(petsclib::PetscLibType, ls::TaoLineSearch, its::PetscInt, f::PetscReal, step::PetscReal) 
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

See also: `TaoLineSearch`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchMonitor"))
"""
function TaoLineSearchMonitor(petsclib::PetscLibType, ls::TaoLineSearch, its::Integer, f::Real, step::Real)
    error("TaoLineSearchMonitor: no generated method for these argument types")
end

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
	TaoLineSearchRegister(petsclib::PetscLibType, sname::String, func::external) 
Adds a line-search algorithm to the registry

Not Collective, No Fortran Support

Input Parameters:
- `sname` - name of a new user-defined solver
- `func`  - routine to Create method context

Calling sequence of `func`:
- `ls` - the `TaoLineSearch` object to set with the `TaoLineSearchType` specific structure

See also: `Tao`, `TaoLineSearch`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchRegister"))
"""
function TaoLineSearchRegister(petsclib::PetscLibType, sname::String, func::external)
    error("TaoLineSearchRegister: no generated method for these argument types")
end

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
	TaoLineSearchReset(petsclib::PetscLibType, ls::TaoLineSearch) 
Some line searches may carry state information
from one `TaoLineSearchApply()` to the next.  This function resets this
state information.

Collective

Input Parameter:
- `ls` - the `TaoLineSearch` context

Level: developer

See also: `Tao`, `TaoLineSearch`, `TaoLineSearchCreate()`, `TaoLineSearchApply()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchReset"))
"""
function TaoLineSearchReset(petsclib::PetscLibType, ls::TaoLineSearch)
    error("TaoLineSearchReset: no generated method for these argument types")
end

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
	TaoLineSearchSetFromOptions(petsclib::PetscLibType, ls::TaoLineSearch) 
Sets various `TaoLineSearch` parameters from user
options.

Collective

Input Parameter:
- `ls` - the `TaoLineSearch` context

Options Database Keys:
- `-tao_ls_type (unit|more-thuente|gpcg|armijo|owarmijo|ipm)` - select which line search `Tao` should use
- `-tao_ls_ftol tol`                                          - tolerance for sufficient decrease
- `-tao_ls_gtol tol`                                          - tolerance for curvature condition
- `-tao_ls_rtol tol`                                          - relative tolerance for acceptable step
- `-tao_ls_stepinit step`                                     - initial steplength allowed
- `-tao_ls_stepmin step`                                      - minimum steplength allowed
- `-tao_ls_stepmax step`                                      - maximum steplength allowed
- `-tao_ls_max_funcs n`                                       - maximum number of function evaluations allowed
- `-tao_ls_view`                                              - display line-search results

Level: beginner

See also: `Tao`, `TaoLineSearch`, `TaoGetLineSearch()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchSetFromOptions"))
"""
function TaoLineSearchSetFromOptions(petsclib::PetscLibType, ls::TaoLineSearch)
    error("TaoLineSearchSetFromOptions: no generated method for these argument types")
end

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
	TaoLineSearchSetGradientRoutine(petsclib::PetscLibType, ls::TaoLineSearch, func::external, ctx::Ptr{Cvoid}) 
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

See also: `Tao`, `TaoLineSearch`, `TaoLineSearchCreate()`, `TaoLineSearchSetObjectiveRoutine()`, `TaoLineSearchSetObjectiveAndGradientRoutine()`, `TaoLineSearchUseTaoRoutines()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchSetGradientRoutine"))
"""
function TaoLineSearchSetGradientRoutine(petsclib::PetscLibType, ls::TaoLineSearch, func::external, ctx::Ptr{Cvoid})
    error("TaoLineSearchSetGradientRoutine: no generated method for these argument types")
end

@for_petsc function TaoLineSearchSetGradientRoutine(petsclib::$UnionPetscLib, ls::TaoLineSearch, func::external, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:TaoLineSearchSetGradientRoutine, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, external, Ptr{Cvoid}),
               ls, func, ctx,
              )


	return nothing
end 

"""
	TaoLineSearchSetInitialStepLength(petsclib::PetscLibType, ls::TaoLineSearch, s::PetscReal) 
Sets the initial step length of a line
search.  If this value is not set then 1.0 is assumed.

Logically Collective

Input Parameters:
- `ls` - the `TaoLineSearch` context
- `s`  - the initial step size

Level: intermediate

See also: `Tao`, `TaoLineSearch`, `TaoLineSearchGetStepLength()`, `TaoLineSearchApply()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchSetInitialStepLength"))
"""
function TaoLineSearchSetInitialStepLength(petsclib::PetscLibType, ls::TaoLineSearch, s::Real)
    error("TaoLineSearchSetInitialStepLength: no generated method for these argument types")
end

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
	TaoLineSearchSetObjectiveAndGTSRoutine(petsclib::PetscLibType, ls::TaoLineSearch, func::external, ctx::Ptr{Cvoid}) 
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

See also: `Tao`, `TaoLineSearch`, `TaoLineSearchCreate()`, `TaoLineSearchSetObjective()`, `TaoLineSearchSetGradient()`, `TaoLineSearchUseTaoRoutines()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchSetObjectiveAndGTSRoutine"))
"""
function TaoLineSearchSetObjectiveAndGTSRoutine(petsclib::PetscLibType, ls::TaoLineSearch, func::external, ctx::Ptr{Cvoid})
    error("TaoLineSearchSetObjectiveAndGTSRoutine: no generated method for these argument types")
end

@for_petsc function TaoLineSearchSetObjectiveAndGTSRoutine(petsclib::$UnionPetscLib, ls::TaoLineSearch, func::external, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:TaoLineSearchSetObjectiveAndGTSRoutine, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, external, Ptr{Cvoid}),
               ls, func, ctx,
              )


	return nothing
end 

"""
	TaoLineSearchSetObjectiveAndGradientRoutine(petsclib::PetscLibType, ls::TaoLineSearch, func::external, ctx::Ptr{Cvoid}) 
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

See also: `Tao`, `TaoLineSearch`, `TaoLineSearchCreate()`, `TaoLineSearchSetObjectiveRoutine()`, `TaoLineSearchSetGradientRoutine()`, `TaoLineSearchUseTaoRoutines()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchSetObjectiveAndGradientRoutine"))
"""
function TaoLineSearchSetObjectiveAndGradientRoutine(petsclib::PetscLibType, ls::TaoLineSearch, func::external, ctx::Ptr{Cvoid})
    error("TaoLineSearchSetObjectiveAndGradientRoutine: no generated method for these argument types")
end

@for_petsc function TaoLineSearchSetObjectiveAndGradientRoutine(petsclib::$UnionPetscLib, ls::TaoLineSearch, func::external, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:TaoLineSearchSetObjectiveAndGradientRoutine, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, external, Ptr{Cvoid}),
               ls, func, ctx,
              )


	return nothing
end 

"""
	TaoLineSearchSetObjectiveRoutine(petsclib::PetscLibType, ls::TaoLineSearch, func::external, ctx::Ptr{Cvoid}) 
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

See also: `Tao`, `TaoLineSearch`, `TaoLineSearchCreate()`, `TaoLineSearchSetGradientRoutine()`, `TaoLineSearchSetObjectiveAndGradientRoutine()`, `TaoLineSearchUseTaoRoutines()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchSetObjectiveRoutine"))
"""
function TaoLineSearchSetObjectiveRoutine(petsclib::PetscLibType, ls::TaoLineSearch, func::external, ctx::Ptr{Cvoid})
    error("TaoLineSearchSetObjectiveRoutine: no generated method for these argument types")
end

@for_petsc function TaoLineSearchSetObjectiveRoutine(petsclib::$UnionPetscLib, ls::TaoLineSearch, func::external, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:TaoLineSearchSetObjectiveRoutine, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, external, Ptr{Cvoid}),
               ls, func, ctx,
              )


	return nothing
end 

"""
	TaoLineSearchSetOptionsPrefix(petsclib::PetscLibType, ls::TaoLineSearch, p::String) 
Sets the prefix used for searching for all
`TaoLineSearch` options in the database.

Logically Collective

Input Parameters:
- `ls` - the `TaoLineSearch` context
- `p`  - the prefix string to prepend to all `ls` option requests

Level: advanced

See also: `Tao`, `TaoLineSearch`, `TaoLineSearchAppendOptionsPrefix()`, `TaoLineSearchGetOptionsPrefix()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchSetOptionsPrefix"))
"""
function TaoLineSearchSetOptionsPrefix(petsclib::PetscLibType, ls::TaoLineSearch, p::String)
    error("TaoLineSearchSetOptionsPrefix: no generated method for these argument types")
end

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
	TaoLineSearchSetType(petsclib::PetscLibType, ls::TaoLineSearch, type::TaoLineSearchType) 
Sets the algorithm used in a line search

Collective

Input Parameters:
- `ls`   - the `TaoLineSearch` context
- `type` - the `TaoLineSearchType` selection

Options Database Key:
- `-tao_ls_type (unit|more-thuente|gpcg|armijo|owarmijo|ipm)` - select which line search `Tao` should use

Level: beginner

See also: `Tao`, `TaoLineSearch`, `TaoLineSearchType`, `TaoLineSearchCreate()`, `TaoLineSearchGetType()`,
`TaoLineSearchApply()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchSetType"))
"""
function TaoLineSearchSetType(petsclib::PetscLibType, ls::TaoLineSearch, type::TaoLineSearchType)
    error("TaoLineSearchSetType: no generated method for these argument types")
end

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
	TaoLineSearchSetUp(petsclib::PetscLibType, ls::TaoLineSearch) 
Sets up the internal data structures for the later use
of a `TaoLineSearch`

Collective

Input Parameter:
- `ls` - the `TaoLineSearch` context

Level: developer

See also: `Tao`, `TaoLineSearch`, `TaoLineSearchCreate()`, `TaoLineSearchApply()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchSetUp"))
"""
function TaoLineSearchSetUp(petsclib::PetscLibType, ls::TaoLineSearch)
    error("TaoLineSearchSetUp: no generated method for these argument types")
end

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
	TaoLineSearchSetVariableBounds(petsclib::PetscLibType, ls::TaoLineSearch, xl::AbstractPetscVec, xu::AbstractPetscVec) 
Sets the upper and lower bounds for a bounded line search

Logically Collective

Input Parameters:
- `ls` - the `TaoLineSearch` context
- `xl` - vector of lower bounds
- `xu` - vector of upper bounds

Level: beginner

See also: `Tao`, `TaoLineSearch`, `TaoSetVariableBounds()`, `TaoLineSearchCreate()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchSetVariableBounds"))
"""
function TaoLineSearchSetVariableBounds(petsclib::PetscLibType, ls::TaoLineSearch, xl::AbstractPetscVec, xu::AbstractPetscVec)
    error("TaoLineSearchSetVariableBounds: no generated method for these argument types")
end

@for_petsc function TaoLineSearchSetVariableBounds(petsclib::$UnionPetscLib, ls::TaoLineSearch, xl::AbstractPetscVec, xu::AbstractPetscVec )

    @chk ccall(
               (:TaoLineSearchSetVariableBounds, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, CVec, CVec),
               ls, xl, xu,
              )


	return nothing
end 

"""
	TaoLineSearchUseTaoRoutines(petsclib::PetscLibType, ls::TaoLineSearch, ts::AbstractTao) 
Informs the `TaoLineSearch` to use the
objective and gradient evaluation routines from the given `Tao` object. The default.

Logically Collective

Input Parameters:
- `ls` - the `TaoLineSearch` context
- `ts` - the `Tao` context with defined objective/gradient evaluation routines

Level: developer

See also: `Tao`, `TaoLineSearch`, `TaoLineSearchCreate()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchUseTaoRoutines"))
"""
function TaoLineSearchUseTaoRoutines(petsclib::PetscLibType, ls::TaoLineSearch, ts::AbstractTao)
    error("TaoLineSearchUseTaoRoutines: no generated method for these argument types")
end

@for_petsc function TaoLineSearchUseTaoRoutines(petsclib::$UnionPetscLib, ls::TaoLineSearch, ts::AbstractTao )

    @chk ccall(
               (:TaoLineSearchUseTaoRoutines, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, CTao),
               ls, ts,
              )


	return nothing
end 

"""
	TaoLineSearchView(petsclib::PetscLibType, ls::TaoLineSearch, viewer::PetscViewer) 
Prints information about the `TaoLineSearch`

Collective

Input Parameters:
- `ls`     - the `TaoLineSearch` context
- `viewer` - visualization context

Options Database Key:
- `-tao_ls_view` - Calls `TaoLineSearchView()` at the end of each line search

Level: beginner

See also: `Tao`, `TaoLineSearch`, `PetscViewerASCIIOpen()`, `TaoLineSearchViewFromOptions()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchView"))
"""
function TaoLineSearchView(petsclib::PetscLibType, ls::TaoLineSearch, viewer::PetscViewer)
    error("TaoLineSearchView: no generated method for these argument types")
end

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
	TaoLineSearchViewFromOptions(petsclib::PetscLibType, A::TaoLineSearch, obj, name::String) 
View a `TaoLineSearch` object based on values in the options database

Collective

Input Parameters:
- `A`    - the `Tao` context
- `obj`  - Optional object
- `name` - command line option

Options Database Key:
- `-name [viewertype][:...]` - option name and values. See `PetscObjectViewFromOptions()` for the possible arguments

Level: intermediate

See also: `Tao`, `TaoLineSearch`, `TaoLineSearchView()`, `PetscObjectViewFromOptions()`, `TaoLineSearchCreate()`

# External Links
$(_doc_external("TaoLineSearch/TaoLineSearchViewFromOptions"))
"""
function TaoLineSearchViewFromOptions(petsclib::PetscLibType, A::TaoLineSearch, obj, name::String)
    error("TaoLineSearchViewFromOptions: no generated method for these argument types")
end

@for_petsc function TaoLineSearchViewFromOptions(petsclib::$UnionPetscLib, A::TaoLineSearch, obj, name::String )

    @chk ccall(
               (:TaoLineSearchViewFromOptions, $petsc_library),
               PetscErrorCode,
               (TaoLineSearch, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

"""
	ctx::TaoMonitorDrawCtx = TaoMonitorDrawCtxCreate(petsclib::PetscLibType, comm::MPI_Comm, host::String, label::String, x::Cint, y::Cint, m::Cint, n::Cint, howoften::PetscInt) 
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

See also: `Tao`, `TaoMonitorSet()`, `TaoMonitorDefault()`, `VecView()`, `TaoMonitorDrawCtx()`

# External Links
$(_doc_external("Tao/TaoMonitorDrawCtxCreate"))
"""
function TaoMonitorDrawCtxCreate(petsclib::PetscLibType, comm::MPI_Comm, host::String, label::String, x::Cint, y::Cint, m::Cint, n::Cint, howoften::Integer)
    error("TaoMonitorDrawCtxCreate: no generated method for these argument types")
end

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
	TaoMonitorDrawCtxDestroy(petsclib::PetscLibType, ictx::Union{TaoMonitorDrawCtx, Ref{TaoMonitorDrawCtx}}) 
Destroys the monitor context for `TaoMonitorSolutionDraw()`

Collective

Input Parameter:
- `ictx` - the monitor context

Level: intermediate

See also: `Tao`, `TaoMonitorSet()`, `TaoMonitorDefault()`, `VecView()`, `TaoMonitorSolutionDraw()`

# External Links
$(_doc_external("Tao/TaoMonitorDrawCtxDestroy"))
"""
function TaoMonitorDrawCtxDestroy(petsclib::PetscLibType, ictx::Union{TaoMonitorDrawCtx, Ref{TaoMonitorDrawCtx}})
    error("TaoMonitorDrawCtxDestroy: no generated method for these argument types")
end

@for_petsc function TaoMonitorDrawCtxDestroy(petsclib::$UnionPetscLib, ictx::Union{TaoMonitorDrawCtx, Ref{TaoMonitorDrawCtx}} )
	ictx_ = ictx isa Base.RefValue ? ictx : Ref{TaoMonitorDrawCtx}(ictx)

    @chk ccall(
               (:TaoMonitorDrawCtxDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{TaoMonitorDrawCtx},),
               ictx_,
              )


	return nothing
end 

"""
	TaoTermComputeGradient(petsclib::PetscLibType, term::TaoTerm, x::AbstractPetscVec, params::AbstractPetscVec, g::AbstractPetscVec) 
Evaluate the gradient of a `TaoTerm` for a given solution vector and parameter vector

Collective

Input Parameters:
- `term`   - a `TaoTerm` representing a parametric function f(x; p)
- `x`      - the solution variable x in f(x; p)
- `params` - the parameters p in f(x; p) (may be NULL if the term is not parametric)

Output Parameter:
- `g` - the value of \\nabla_x f(x; p)

Level: developer

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermComputeObjective()`,
`TaoTermComputeObjectiveAndGradient()`,
`TaoTermComputeHessian()`,
`TaoTermShellSetGradient()`

# External Links
$(_doc_external("TaoTerm/TaoTermComputeGradient"))
"""
function TaoTermComputeGradient(petsclib::PetscLibType, term::TaoTerm, x::AbstractPetscVec, params::AbstractPetscVec, g::AbstractPetscVec)
    error("TaoTermComputeGradient: no generated method for these argument types")
end

@for_petsc function TaoTermComputeGradient(petsclib::$UnionPetscLib, term::TaoTerm, x::AbstractPetscVec, params::AbstractPetscVec, g::AbstractPetscVec )

    @chk ccall(
               (:TaoTermComputeGradient, $petsc_library),
               PetscErrorCode,
               (TaoTerm, CVec, CVec, CVec),
               term, x, params, g,
              )


	return nothing
end 

"""
	TaoTermComputeGradientFD(petsclib::PetscLibType, term::TaoTerm, x::AbstractPetscVec, params::AbstractPetscVec, g::AbstractPetscVec) 
Approximate the gradient of a `TaoTerm` using finite differences

Collective

Input Parameters:
- `term`   - a `TaoTerm`
- `x`      - a solution vector
- `params` - parameters vector (may be `NULL`, see `TaoTermParametersMode`)

Output Parameter:
- `g` - the computed finite difference approximation to the gradient

Options Database Keys:
- `-tao_term_fd_delta <delta>`       - change in `x` used to calculate finite differences
- `-tao_term_gradient_use_fd <bool>` - Use `TaoTermComputeGradientFD()` in `TaoTermComputeGradient()`

Level: advanced

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermGetFDDelta()`,
`TaoTermSetFDDelta()`,
`TaoTermComputeGradientSetUseFD()`,
`TaoTermComputeGradientGetUseFD()`,
`TaoTermComputeHessianFD()`

# External Links
$(_doc_external("TaoTerm/TaoTermComputeGradientFD"))
"""
function TaoTermComputeGradientFD(petsclib::PetscLibType, term::TaoTerm, x::AbstractPetscVec, params::AbstractPetscVec, g::AbstractPetscVec)
    error("TaoTermComputeGradientFD: no generated method for these argument types")
end

@for_petsc function TaoTermComputeGradientFD(petsclib::$UnionPetscLib, term::TaoTerm, x::AbstractPetscVec, params::AbstractPetscVec, g::AbstractPetscVec )

    @chk ccall(
               (:TaoTermComputeGradientFD, $petsc_library),
               PetscErrorCode,
               (TaoTerm, CVec, CVec, CVec),
               term, x, params, g,
              )


	return nothing
end 

"""
	use_fd::PetscBool = TaoTermComputeGradientGetUseFD(petsclib::PetscLibType, term::TaoTerm) 
Get whether finite differences are used in `TaoTermComputeGradient()`.

Not collective

Input Parameter:
- `term` - a `TaoTerm`

Output Parameter:
- `use_fd` - `PETSC_TRUE` if finite differences are used

Level: advanced

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermGetFDDelta()`,
`TaoTermSetFDDelta()`,
`TaoTermComputeGradientFD()`,
`TaoTermComputeGradientSetUseFD()`,
`TaoTermComputeHessianFD()`,
`TaoTermComputeHessianSetUseFD()`,
`TaoTermComputeHessianGetUseFD()`

# External Links
$(_doc_external("TaoTerm/TaoTermComputeGradientGetUseFD"))
"""
function TaoTermComputeGradientGetUseFD(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermComputeGradientGetUseFD: no generated method for these argument types")
end

@for_petsc function TaoTermComputeGradientGetUseFD(petsclib::$UnionPetscLib, term::TaoTerm )
	use_fd_ = Ref{PetscBool}()

    @chk ccall(
               (:TaoTermComputeGradientGetUseFD, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{PetscBool}),
               term, use_fd_,
              )

	use_fd = use_fd_[]

	return use_fd
end 

"""
	TaoTermComputeGradientSetUseFD(petsclib::PetscLibType, term::TaoTerm, use_fd::PetscBool) 
Set whether to use finite differences instead of the user-provided or built-in gradient method in `TaoTermComputeGradient()`.

Logically collective

Input Parameters:
- `term`   - a `TaoTerm`
- `use_fd` - `PETSC_TRUE` to use finite differences, `PETSC_FALSE` to use the user-provided or built-in gradient method

Options Database Keys:
- `-tao_term_gradient_use_fd <bool>` - use finite differences for gradient computation

Level: advanced

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermGetFDDelta()`,
`TaoTermSetFDDelta()`,
`TaoTermComputeGradientFD()`,
`TaoTermComputeGradientGetUseFD()`,
`TaoTermComputeHessianFD()`,
`TaoTermComputeHessianSetUseFD()`,
`TaoTermComputeHessianGetUseFD()`

# External Links
$(_doc_external("TaoTerm/TaoTermComputeGradientSetUseFD"))
"""
function TaoTermComputeGradientSetUseFD(petsclib::PetscLibType, term::TaoTerm, use_fd::PetscBool)
    error("TaoTermComputeGradientSetUseFD: no generated method for these argument types")
end

@for_petsc function TaoTermComputeGradientSetUseFD(petsclib::$UnionPetscLib, term::TaoTerm, use_fd::PetscBool )

    @chk ccall(
               (:TaoTermComputeGradientSetUseFD, $petsc_library),
               PetscErrorCode,
               (TaoTerm, PetscBool),
               term, use_fd,
              )


	return nothing
end 

"""
	TaoTermComputeHessian(petsclib::PetscLibType, term::TaoTerm, x::AbstractPetscVec, params::AbstractPetscVec, H::AbstractPetscMat, Hpre::AbstractPetscMat) 
Evaluate the Hessian of a `TaoTerm`
(with respect to the solution variables) for a given solution vector and parameter vector

Collective

Input Parameters:
- `term`   - a `TaoTerm` representing a parametric function f(x; p)
- `x`      - the solution variable x in f(x; p)
- `params` - the parameters p in f(x; p) (may be `NULL` if the term is not parametric)

Output Parameters:
- `H`    - Hessian matrix \\nabla_x^2 f(x;p)
- `Hpre` - an (approximate) Hessian from which the preconditioner will be constructed, often the same as `H`

Level: developer

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermComputeObjective()`,
`TaoTermComputeGradient()`,
`TaoTermComputeObjectiveAndGradient()`,
`TaoTermShellSetHessian()`

# External Links
$(_doc_external("TaoTerm/TaoTermComputeHessian"))
"""
function TaoTermComputeHessian(petsclib::PetscLibType, term::TaoTerm, x::AbstractPetscVec, params::AbstractPetscVec, H::AbstractPetscMat, Hpre::AbstractPetscMat)
    error("TaoTermComputeHessian: no generated method for these argument types")
end

@for_petsc function TaoTermComputeHessian(petsclib::$UnionPetscLib, term::TaoTerm, x::AbstractPetscVec, params::AbstractPetscVec, H::AbstractPetscMat, Hpre::AbstractPetscMat )

    @chk ccall(
               (:TaoTermComputeHessian, $petsc_library),
               PetscErrorCode,
               (TaoTerm, CVec, CVec, CMat, CMat),
               term, x, params, H, Hpre,
              )


	return nothing
end 

"""
	TaoTermComputeHessianFD(petsclib::PetscLibType, term::TaoTerm, x::AbstractPetscVec, params::AbstractPetscVec, H::AbstractPetscMat, Hpre::AbstractPetscMat) 
Use finite difference to compute Hessian matrix.

Collective

Input Parameters:
- `term`   - a `TaoTerm`
- `x`      - a solution vector
- `params` - parameters vector (may be `NULL`, see `TaoTermParametersMode`)

Output Parameters:
- `H`    - (optional) Hessian matrix
- `Hpre` - (optional) Hessian preconditioning matrix

Options Database Keys:
- `-tao_term_fd_delta <delta>`      - change in X used to calculate finite differences
- `-tao_term_hessian_use_fd <bool>` - Use `TaoTermComputeHessianFD()` in `TaoTermComputeHessian()`

Level: advanced

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermComputeHessian()`,
`TaoTermGetFDDelta()`,
`TaoTermSetFDDelta()`,
`TaoTermComputeHessianSetUseFD()`,
`TaoTermComputeHessianGetUseFD()`

# External Links
$(_doc_external("TaoTerm/TaoTermComputeHessianFD"))
"""
function TaoTermComputeHessianFD(petsclib::PetscLibType, term::TaoTerm, x::AbstractPetscVec, params::AbstractPetscVec, H::AbstractPetscMat, Hpre::AbstractPetscMat)
    error("TaoTermComputeHessianFD: no generated method for these argument types")
end

@for_petsc function TaoTermComputeHessianFD(petsclib::$UnionPetscLib, term::TaoTerm, x::AbstractPetscVec, params::AbstractPetscVec, H::AbstractPetscMat, Hpre::AbstractPetscMat )

    @chk ccall(
               (:TaoTermComputeHessianFD, $petsc_library),
               PetscErrorCode,
               (TaoTerm, CVec, CVec, CMat, CMat),
               term, x, params, H, Hpre,
              )


	return nothing
end 

"""
	use_fd::PetscBool = TaoTermComputeHessianGetUseFD(petsclib::PetscLibType, term::TaoTerm) 
Get whether finite differences are used in `TaoTermComputeHessian()`.

Not collective

Input Parameter:
- `term` - a `TaoTerm`

Output Parameter:
- `use_fd` - `PETSC_TRUE` if finite differences are used

Level: advanced

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermGetFDDelta()`,
`TaoTermSetFDDelta()`,
`TaoTermComputeGradientFD()`,
`TaoTermComputeGradientSetUseFD()`,
`TaoTermComputeGradientGetUseFD()`,
`TaoTermComputeHessianFD()`,
`TaoTermComputeHessianSetUseFD()`

# External Links
$(_doc_external("TaoTerm/TaoTermComputeHessianGetUseFD"))
"""
function TaoTermComputeHessianGetUseFD(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermComputeHessianGetUseFD: no generated method for these argument types")
end

@for_petsc function TaoTermComputeHessianGetUseFD(petsclib::$UnionPetscLib, term::TaoTerm )
	use_fd_ = Ref{PetscBool}()

    @chk ccall(
               (:TaoTermComputeHessianGetUseFD, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{PetscBool}),
               term, use_fd_,
              )

	use_fd = use_fd_[]

	return use_fd
end 

"""
	TaoTermComputeHessianMFFD(petsclib::PetscLibType, term::TaoTerm, x::AbstractPetscVec, params::AbstractPetscVec, H::AbstractPetscMat, B::AbstractPetscMat) 
Update a matrix-free finite-difference `MATMFFD` Hessian created by
`TaoTermCreateHessianMFFD()` to represent the Hessian of a `TaoTerm` at a given point and parameters.

Collective

Input Parameters:
- `term`   - the `TaoTerm`
- `x`      - the point at which the Hessian is to be applied
- `params` - the current parameter vector for `term`, or `NULL`

Output Parameters:
- `H` - the `MATMFFD` Hessian, reinitialized if needed and updated to base point `x`
- `B` - the preconditioning matrix (unused; retained for API symmetry), or `NULL`

Level: advanced

See also: `TaoTerm`, `TaoTermCreateHessianMFFD()`, `TaoTermComputeHessian()`, `MATMFFD`

# External Links
$(_doc_external("TaoTerm/TaoTermComputeHessianMFFD"))
"""
function TaoTermComputeHessianMFFD(petsclib::PetscLibType, term::TaoTerm, x::AbstractPetscVec, params::AbstractPetscVec, H::AbstractPetscMat, B::AbstractPetscMat)
    error("TaoTermComputeHessianMFFD: no generated method for these argument types")
end

@for_petsc function TaoTermComputeHessianMFFD(petsclib::$UnionPetscLib, term::TaoTerm, x::AbstractPetscVec, params::AbstractPetscVec, H::AbstractPetscMat, B::AbstractPetscMat )

    @chk ccall(
               (:TaoTermComputeHessianMFFD, $petsc_library),
               PetscErrorCode,
               (TaoTerm, CVec, CVec, CMat, CMat),
               term, x, params, H, B,
              )


	return nothing
end 

"""
	TaoTermComputeHessianSetUseFD(petsclib::PetscLibType, term::TaoTerm, use_fd::PetscBool) 
Set whether to use finite differences instead of the user-provided or built-in methods in `TaoTermComputeHessian()`.

Logically collective

Input Parameters:
- `term`   - a `TaoTerm`
- `use_fd` - `PETSC_TRUE` to use finite differences, `PETSC_FALSE` to use the user-provided or built-in Hessian method

Options Database Keys:
- `-tao_term_hessian_use_fd <bool>` - use finite differences for Hessian computation

Level: advanced

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermGetFDDelta()`,
`TaoTermSetFDDelta()`,
`TaoTermComputeGradientFD()`,
`TaoTermComputeGradientSetUseFD()`,
`TaoTermComputeGradientGetUseFD()`,
`TaoTermComputeHessianFD()`,
`TaoTermComputeHessianGetUseFD()`

# External Links
$(_doc_external("TaoTerm/TaoTermComputeHessianSetUseFD"))
"""
function TaoTermComputeHessianSetUseFD(petsclib::PetscLibType, term::TaoTerm, use_fd::PetscBool)
    error("TaoTermComputeHessianSetUseFD: no generated method for these argument types")
end

@for_petsc function TaoTermComputeHessianSetUseFD(petsclib::$UnionPetscLib, term::TaoTerm, use_fd::PetscBool )

    @chk ccall(
               (:TaoTermComputeHessianSetUseFD, $petsc_library),
               PetscErrorCode,
               (TaoTerm, PetscBool),
               term, use_fd,
              )


	return nothing
end 

"""
	value::PetscReal = TaoTermComputeObjective(petsclib::PetscLibType, term::TaoTerm, x::AbstractPetscVec, params::AbstractPetscVec) 
Evaluate a `TaoTerm` for a given solution vector and parameter vector

Collective

Input Parameters:
- `term`   - a `TaoTerm` representing a parametric function f(x; p)
- `x`      - the solution variable x in f(x; p)
- `params` - the parameters p in f(x; p) (may be `NULL` if the term is not parametric)

Output Parameter:
- `value` - the value of f(x; p)

Level: developer

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermComputeGradient()`,
`TaoTermComputeObjectiveAndGradient()`,
`TaoTermComputeHessian()`,
`TaoTermShellSetObjective()`

# External Links
$(_doc_external("TaoTerm/TaoTermComputeObjective"))
"""
function TaoTermComputeObjective(petsclib::PetscLibType, term::TaoTerm, x::AbstractPetscVec, params::AbstractPetscVec)
    error("TaoTermComputeObjective: no generated method for these argument types")
end

@for_petsc function TaoTermComputeObjective(petsclib::$UnionPetscLib, term::TaoTerm, x::AbstractPetscVec, params::AbstractPetscVec )
	value_ = Ref{$PetscReal}()

    @chk ccall(
               (:TaoTermComputeObjective, $petsc_library),
               PetscErrorCode,
               (TaoTerm, CVec, CVec, Ptr{$PetscReal}),
               term, x, params, value_,
              )

	value = value_[]

	return value
end 

"""
	value::PetscReal = TaoTermComputeObjectiveAndGradient(petsclib::PetscLibType, term::TaoTerm, x::AbstractPetscVec, params::AbstractPetscVec, g::AbstractPetscVec) 
Evaluate both the value and gradient of
a `TaoTerm` for a given set of solution vector and parameter vector

Collective

Input Parameters:
- `term`   - a `TaoTerm` representing a parametric function f(x; p)
- `x`      - the solution variable x in f(x; p)
- `params` - the parameters p in f(x; p) (may be NULL if the term is not parametric)

Output Parameters:
- `value` - the value of f(x; p)
- `g`     - the value of \\nabla_x f(x; p)

Level: developer

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermComputeObjective()`,
`TaoTermComputeGradient()`,
`TaoTermComputeHessian()`,
`TaoTermShellSetObjectiveAndGradient()`

# External Links
$(_doc_external("TaoTerm/TaoTermComputeObjectiveAndGradient"))
"""
function TaoTermComputeObjectiveAndGradient(petsclib::PetscLibType, term::TaoTerm, x::AbstractPetscVec, params::AbstractPetscVec, g::AbstractPetscVec)
    error("TaoTermComputeObjectiveAndGradient: no generated method for these argument types")
end

@for_petsc function TaoTermComputeObjectiveAndGradient(petsclib::$UnionPetscLib, term::TaoTerm, x::AbstractPetscVec, params::AbstractPetscVec, g::AbstractPetscVec )
	value_ = Ref{$PetscReal}()

    @chk ccall(
               (:TaoTermComputeObjectiveAndGradient, $petsc_library),
               PetscErrorCode,
               (TaoTerm, CVec, CVec, Ptr{$PetscReal}, CVec),
               term, x, params, value_, g,
              )

	value = value_[]

	return value
end 

"""
	term::TaoTerm = TaoTermCreate(petsclib::PetscLibType, comm::MPI_Comm) 
Create a TaoTerm to use in defining the function `Tao` is to optimize

Collective

Input Parameter:
- `comm` - communicator for MPI processes that compute the term

Output Parameter:
- `term` - a new `TaoTerm`

Level: beginner

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermSetType()`,
`TaoAddTerm()`,
`TaoTermSetFromOptions()`,
`TaoTermSetUp()`,
`TaoTermView()`,
`TaoTermDestroy()`

# External Links
$(_doc_external("TaoTerm/TaoTermCreate"))
"""
function TaoTermCreate(petsclib::PetscLibType, comm::MPI_Comm)
    error("TaoTermCreate: no generated method for these argument types")
end

@for_petsc function TaoTermCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	term_ = Ref{TaoTerm}()

    @chk ccall(
               (:TaoTermCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{TaoTerm}),
               comm, term_,
              )

	term = term_[]

	return term
end 

"""
	term::TaoTerm = TaoTermCreateHalfL2Squared(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, M_N::PetscInt) 
Create a `TaoTerm` for the objective term \\tfrac{1}{2}\\|x - p\\|_2^2, for solution x and parameters p.

Collective

Input Parameters:
- `comm` - the MPI communicator where the `TaoTerm` will be computed
- `n`    - the local size of the x and p vectors (or `PETSC_DECIDE`)
- `N`    - the global size of the x and p vectors (or `PETSC_DECIDE`)

Output Parameter:
- `term` - the `TaoTerm`

Level: beginner

See also: [](sec_tao_term),
`TaoTerm`,
`TAOTERMHALFL2SQUARED`,
`TaoTermCreateL1()`,
`TaoTermCreateQuadratic()`

# External Links
$(_doc_external("TaoTerm/TaoTermCreateHalfL2Squared"))
"""
function TaoTermCreateHalfL2Squared(petsclib::PetscLibType, comm::MPI_Comm, n::Integer, M_N::Integer)
    error("TaoTermCreateHalfL2Squared: no generated method for these argument types")
end

@for_petsc function TaoTermCreateHalfL2Squared(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, M_N::$PetscInt )
	term_ = Ref{TaoTerm}()

    @chk ccall(
               (:TaoTermCreateHalfL2Squared, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{TaoTerm}),
               comm, n, M_N, term_,
              )

	term = term_[]

	return term
end 

"""
	mffd::PetscMat = TaoTermCreateHessianMFFD(petsclib::PetscLibType, term::TaoTerm) 
Create a `MATMFFD` for a matrix-free finite-difference approximation of the Hessian of a `TaoTerm`

Collective

Input Parameter:
- `term` - a `TaoTerm`

Output Parameter:
- `mffd` - a `Mat` of type `MATMFFD`

Level: advanced

See also: [](sec_tao_term), `TaoTerm`, `TaoTermComputeHessianFD()`

# External Links
$(_doc_external("TaoTerm/TaoTermCreateHessianMFFD"))
"""
function TaoTermCreateHessianMFFD(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermCreateHessianMFFD: no generated method for these argument types")
end

@for_petsc function TaoTermCreateHessianMFFD(petsclib::$UnionPetscLib, term::TaoTerm )
	mffd_ = Ref{CMat}()

    @chk ccall(
               (:TaoTermCreateHessianMFFD, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{CMat}),
               term, mffd_,
              )

	mffd = PetscMat(mffd_[], petsclib)

	return mffd
end 

"""
	H::PetscMat,Hpre::PetscMat = TaoTermCreateHessianMatrices(petsclib::PetscLibType, term::TaoTerm) 
Create the matrices that can be inputs to `TaoTermComputeHessian()`

Collective

Input Parameter:
- `term` - a `TaoTerm`

Output Parameters:
- `H`    - (optional) a matrix that can store the Hessian computed in `TaoTermComputeHessian()`
- `Hpre` - (optional) a matrix from which a preconditioner can be computed in `TaoTermComputeHessian()`

Level: advanced

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermComputeHessian()`,
`TaoTermShellSetCreateHessianMatrices()`,
`TaoTermCreateSolutionVec()`,
`TaoTermCreateHessianMatricesDefault()`,
`TaoTermGetCreateHessianMode()`,
`TaoTermSetCreateHessianMode()`,
`TaoTermIsCreateHessianMatricesDefined()`

# External Links
$(_doc_external("TaoTerm/TaoTermCreateHessianMatrices"))
"""
function TaoTermCreateHessianMatrices(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermCreateHessianMatrices: no generated method for these argument types")
end

@for_petsc function TaoTermCreateHessianMatrices(petsclib::$UnionPetscLib, term::TaoTerm )
	H_ = Ref{CMat}()
	Hpre_ = Ref{CMat}()

    @chk ccall(
               (:TaoTermCreateHessianMatrices, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{CMat}, Ptr{CMat}),
               term, H_, Hpre_,
              )

	H = PetscMat(H_[], petsclib)
	Hpre = PetscMat(Hpre_[], petsclib)

	return H,Hpre
end 

"""
	H::PetscMat,Hpre::PetscMat = TaoTermCreateHessianMatricesDefault(petsclib::PetscLibType, term::TaoTerm) 
Default routine for creating Hessian matrices that can be used by many `TaoTerm` implementations

Collective

Input Parameter:
- `term` - a `TaoTerm`

Output Parameters:
- `H`    - (optional) a matrix that can store the Hessian computed in `TaoTermComputeHessian()`
- `Hpre` - (optional) a matrix from which a preconditioner can be computed in `TaoTermComputeHessian()`

Level: developer

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermComputeHessian()`,
`TaoTermCreateHessianMatrices()`,
`TaoTermGetCreateHessianMode()`,
`TaoTermSetCreateHessianMode()`

# External Links
$(_doc_external("TaoTerm/TaoTermCreateHessianMatricesDefault"))
"""
function TaoTermCreateHessianMatricesDefault(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermCreateHessianMatricesDefault: no generated method for these argument types")
end

@for_petsc function TaoTermCreateHessianMatricesDefault(petsclib::$UnionPetscLib, term::TaoTerm )
	H_ = Ref{CMat}()
	Hpre_ = Ref{CMat}()

    @chk ccall(
               (:TaoTermCreateHessianMatricesDefault, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{CMat}, Ptr{CMat}),
               term, H_, Hpre_,
              )

	H = PetscMat(H_[], petsclib)
	Hpre = PetscMat(Hpre_[], petsclib)

	return H,Hpre
end 

"""
	term::TaoTerm = TaoTermCreateL1(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, M_N::PetscInt, epsilon::PetscReal) 
Create a `TaoTerm` for the objective function term \\|x - p\\|_1.

Collective

Input Parameters:
- `comm`    - the MPI communicator where the term will be computed
- `n`       - the local size of the x and p vectors (or `PETSC_DECIDE`)
- `N`       - the global size of the x and p vectors (or `PETSC_DECIDE`)
- `epsilon` - a non-negative smoothing parameter (see `TaoTermL1SetEpsilon()`)

Output Parameter:
- `term` - the `TaoTerm`

Level: beginner

See also: [](sec_tao_term),
`TaoTerm`,
`TAOTERML1`,
`TaoTermL1GetEpsilon()`,
`TaoTermL1SetEpsilon()`,
`TaoTermCreateHalfL2Squared()`,
`TaoTermCreateQuadratic()`

# External Links
$(_doc_external("TaoTerm/TaoTermCreateL1"))
"""
function TaoTermCreateL1(petsclib::PetscLibType, comm::MPI_Comm, n::Integer, M_N::Integer, epsilon::Real)
    error("TaoTermCreateL1: no generated method for these argument types")
end

@for_petsc function TaoTermCreateL1(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, M_N::$PetscInt, epsilon::$PetscReal )
	term_ = Ref{TaoTerm}()

    @chk ccall(
               (:TaoTermCreateL1, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscReal, Ptr{TaoTerm}),
               comm, n, M_N, epsilon, term_,
              )

	term = term_[]

	return term
end 

"""
	parameters::PetscVec = TaoTermCreateParametersVec(petsclib::PetscLibType, term::TaoTerm) 
Create a parameter vector for a `TaoTerm`

Collective

Input Parameter:
- `term` - a `TaoTerm`

Output Parameter:
- `parameters` - a compatible parameter vector for `term`

Level: advanced

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermShellSetCreateParametersVec()`,
`TaoTermGetParametersSizes()`,
`TaoTermSetParametersSizes()`,
`TaoTermSetParametersTemplate()`,
`TaoTermGetParametersVecType()`,
`TaoTermSetParametersVecType()`,
`TaoTermGetParametersLayout()`,
`TaoTermSetParametersLayout()`,
`TaoTermCreateHessianMatrices()`

# External Links
$(_doc_external("TaoTerm/TaoTermCreateParametersVec"))
"""
function TaoTermCreateParametersVec(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermCreateParametersVec: no generated method for these argument types")
end

@for_petsc function TaoTermCreateParametersVec(petsclib::$UnionPetscLib, term::TaoTerm )
	parameters_ = Ref{CVec}()

    @chk ccall(
               (:TaoTermCreateParametersVec, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{CVec}),
               term, parameters_,
              )

	parameters = PetscVec(parameters_[], petsclib)

	return parameters
end 

"""
	term::TaoTerm = TaoTermCreateQuadratic(petsclib::PetscLibType, A::AbstractPetscMat) 
Create a `TAOTERMQUADRATIC` for a given matrix

Collective

Input Parameter:
- `A` - a square matrix

Output Parameter:
- `term` - a `TaoTerm` that implements \\tfrac{1}{2}(x - p)^T A (x - p)

Level: beginner

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermCreate()`,
`TAOTERMQUADRATIC`,
`TaoTermCreateHalfL2Squared()`,
`TaoTermCreateL1()`, `TaoTermQuadraticSetMat()`

# External Links
$(_doc_external("TaoTerm/TaoTermCreateQuadratic"))
"""
function TaoTermCreateQuadratic(petsclib::PetscLibType, A::AbstractPetscMat)
    error("TaoTermCreateQuadratic: no generated method for these argument types")
end

@for_petsc function TaoTermCreateQuadratic(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	term_ = Ref{TaoTerm}()

    @chk ccall(
               (:TaoTermCreateQuadratic, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{TaoTerm}),
               A, term_,
              )

	term = term_[]

	return term
end 

"""
	term::TaoTerm = TaoTermCreateShell(petsclib::PetscLibType, comm::MPI_Comm, ctx::Ptr{Cvoid}, destroy::Ptr{Cvoid}) 
Create a `TaoTerm` of type `TAOTERMSHELL` that is ready to accept user-provided callback operations.

Collective

Input Parameters:
- `comm`    - the MPI communicator for computing the term
- `ctx`     - (optional) a context to be used by routines
- `destroy` - (optional) a routine to destroy the context when `term` is destroyed

Output Parameter:
- `term` - a `TaoTerm` of type `TAOTERMSHELL`

Level: intermediate

See also: [](sec_tao_term), `TaoTerm`, `TAOTERMSHELL`

# External Links
$(_doc_external("TaoTerm/TaoTermCreateShell"))
"""
function TaoTermCreateShell(petsclib::PetscLibType, comm::MPI_Comm, ctx::Ptr{Cvoid}, destroy::Ptr{Cvoid})
    error("TaoTermCreateShell: no generated method for these argument types")
end

@for_petsc function TaoTermCreateShell(petsclib::$UnionPetscLib, comm::MPI_Comm, ctx::Ptr{Cvoid}, destroy::Ptr{Cvoid} )
	term_ = Ref{TaoTerm}()

    @chk ccall(
               (:TaoTermCreateShell, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{TaoTerm}),
               comm, ctx, destroy, term_,
              )

	term = term_[]

	return term
end 

"""
	solution::PetscVec = TaoTermCreateSolutionVec(petsclib::PetscLibType, term::TaoTerm) 
Create a solution vector for a `TaoTerm`

Collective

Input Parameter:
- `term` - a `TaoTerm`

Output Parameter:
- `solution` - a compatible solution vector for `term`

Level: advanced

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermShellSetCreateSolutionVec()`,
`TaoTermGetSolutionSizes()`,
`TaoTermSetSolutionSizes()`,
`TaoTermSetSolutionTemplate()`,
`TaoTermGetSolutionVecType()`,
`TaoTermSetSolutionVecType()`,
`TaoTermGetSolutionLayout()`,
`TaoTermSetSolutionLayout()`,
`TaoTermCreateHessianMatrices()`

# External Links
$(_doc_external("TaoTerm/TaoTermCreateSolutionVec"))
"""
function TaoTermCreateSolutionVec(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermCreateSolutionVec: no generated method for these argument types")
end

@for_petsc function TaoTermCreateSolutionVec(petsclib::$UnionPetscLib, term::TaoTerm )
	solution_ = Ref{CVec}()

    @chk ccall(
               (:TaoTermCreateSolutionVec, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{CVec}),
               term, solution_,
              )

	solution = PetscVec(solution_[], petsclib)

	return solution
end 

"""
	TaoTermDestroy(petsclib::PetscLibType, term::Union{TaoTerm, Ref{TaoTerm}}) 
Destroy a `TaoTerm`.

Collective

Input Parameter:
- `term` - a `TaoTerm`

Level: beginner

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermCreate()`,
`TaoTermSetType()`,
`TaoTermSetFromOptions()`,
`TaoTermSetUp()`,
`TaoTermView()`

# External Links
$(_doc_external("TaoTerm/TaoTermDestroy"))
"""
function TaoTermDestroy(petsclib::PetscLibType, term::Union{TaoTerm, Ref{TaoTerm}})
    error("TaoTermDestroy: no generated method for these argument types")
end

@for_petsc function TaoTermDestroy(petsclib::$UnionPetscLib, term::Union{TaoTerm, Ref{TaoTerm}} )
	term_ = term isa Base.RefValue ? term : Ref{TaoTerm}(term)

    @chk ccall(
               (:TaoTermDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{TaoTerm},),
               term_,
              )


	return nothing
end 

"""
	newterm::TaoTerm = TaoTermDuplicate(petsclib::PetscLibType, term::TaoTerm, opt::TaoTermDuplicateOption) 
Duplicate a `TaoTerm`

Collective

Input Parameters:
- `term` - a `TaoTerm`
- `opt`  - `TAOTERM_DUPLICATE_SIZEONLY` or `TAOTERM_DUPLICATE_TYPE`

Output Parameter:
- `newterm` - the duplicate `TaoTerm`

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermDuplicateOption`

# External Links
$(_doc_external("TaoTerm/TaoTermDuplicate"))
"""
function TaoTermDuplicate(petsclib::PetscLibType, term::TaoTerm, opt::TaoTermDuplicateOption)
    error("TaoTermDuplicate: no generated method for these argument types")
end

@for_petsc function TaoTermDuplicate(petsclib::$UnionPetscLib, term::TaoTerm, opt::TaoTermDuplicateOption )
	newterm_ = Ref{TaoTerm}()

    @chk ccall(
               (:TaoTermDuplicate, $petsc_library),
               PetscErrorCode,
               (TaoTerm, TaoTermDuplicateOption, Ptr{TaoTerm}),
               term, opt, newterm_,
              )

	newterm = newterm_[]

	return newterm
end 

"""
	Hpre_is_H::PetscBool,H_mattype::MatType,Hpre_mattype::MatType = TaoTermGetCreateHessianMode(petsclib::PetscLibType, term::TaoTerm) 
Get the behavior of `TaoTermCreateHessianMatricesDefault()`.

Not collective

Input Parameter:
- `term` - a `TaoTerm`

Output Parameters:
- `Hpre_is_H`    - (optional) should `TaoTermCreateHessianMatricesDefault()` make one matrix for `H` and `Hpre`?
- `H_mattype`    - (optional) the `MatType` to create for `H`
- `Hpre_mattype` - (optional) the `MatType` to create for `Hpre`

Level: developer

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermComputeHessian()`,
`TaoTermCreateHessianMatrices()`,
`TaoTermCreateHessianMatricesDefault()`,
`TaoTermSetCreateHessianMode()`

# External Links
$(_doc_external("TaoTerm/TaoTermGetCreateHessianMode"))
"""
function TaoTermGetCreateHessianMode(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermGetCreateHessianMode: no generated method for these argument types")
end

@for_petsc function TaoTermGetCreateHessianMode(petsclib::$UnionPetscLib, term::TaoTerm )
	Hpre_is_H_ = Ref{PetscBool}()
	H_mattype_ = Ref{MatType}()
	Hpre_mattype_ = Ref{MatType}()

    @chk ccall(
               (:TaoTermGetCreateHessianMode, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{PetscBool}, Ptr{MatType}, Ptr{MatType}),
               term, Hpre_is_H_, H_mattype_, Hpre_mattype_,
              )

	Hpre_is_H = Hpre_is_H_[]
	H_mattype = H_mattype_[] == C_NULL ? "" : unsafe_string(H_mattype_[])
	Hpre_mattype = Hpre_mattype_[] == C_NULL ? "" : unsafe_string(Hpre_mattype_[])

	return Hpre_is_H,H_mattype,Hpre_mattype
end 

"""
	delta::PetscReal = TaoTermGetFDDelta(petsclib::PetscLibType, term::TaoTerm) 
Get the increment used for finite difference derivative approximations in methods like `TaoTermComputeGradientFD()`

Not collective

Input Parameter:
- `term` - a `TaoTerm`

Output Parameter:
- `delta` - the finite difference increment

Options Database Key:
- `-tao_term_fd_delta <delta>` - the above increment

Level: advanced

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermSetFDDelta()`,
`TaoTermComputeGradientFD()`,
`TaoTermComputeGradientSetUseFD()`,
`TaoTermComputeGradientGetUseFD()`

# External Links
$(_doc_external("TaoTerm/TaoTermGetFDDelta"))
"""
function TaoTermGetFDDelta(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermGetFDDelta: no generated method for these argument types")
end

@for_petsc function TaoTermGetFDDelta(petsclib::$UnionPetscLib, term::TaoTerm )
	delta_ = Ref{$PetscReal}()

    @chk ccall(
               (:TaoTermGetFDDelta, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{$PetscReal}),
               term, delta_,
              )

	delta = delta_[]

	return delta
end 

"""
	parameters_layout::PetscLayout = TaoTermGetParametersLayout(petsclib::PetscLibType, term::TaoTerm) 
Get the layouts describing the parameter vectors of a `TaoTerm`.

Not collective

Input Parameter:
- `term` - a `TaoTerm`

Output Parameter:
- `parameters_layout` - the `PetscLayout` for the parameter space

Level: intermediate

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermGetParametersVecType()`,
`TaoTermSetParametersVecType()`,
`TaoTermSetParametersLayout()`,
`TaoTermSetSolutionTemplate()`,
`TaoTermSetParametersTemplate()`,
`TaoTermCreateParametersVec()`

# External Links
$(_doc_external("TaoTerm/TaoTermGetParametersLayout"))
"""
function TaoTermGetParametersLayout(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermGetParametersLayout: no generated method for these argument types")
end

@for_petsc function TaoTermGetParametersLayout(petsclib::$UnionPetscLib, term::TaoTerm )
	parameters_layout_ = Ref{PetscLayout}()

    @chk ccall(
               (:TaoTermGetParametersLayout, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{PetscLayout}),
               term, parameters_layout_,
              )

	parameters_layout = parameters_layout_[]

	return parameters_layout
end 

"""
	parameters_mode::TaoTermParametersMode = TaoTermGetParametersMode(petsclib::PetscLibType, term::TaoTerm) 
Gets the way a `TaoTerm` can accept parameters

Not collective

Input Parameter:
- `term` - a `TaoTerm`

Output Parameter:
- `parameters_mode` - `TAOTERM_PARAMETERS_OPTIONAL`, `TAOTERM_PARAMETERS_NONE`, `TAOTERM_PARAMETERS_REQUIRED`

Level: intermediate

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermParametersMode`,
`TaoTermSetParametersMode()`

# External Links
$(_doc_external("TaoTerm/TaoTermGetParametersMode"))
"""
function TaoTermGetParametersMode(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermGetParametersMode: no generated method for these argument types")
end

@for_petsc function TaoTermGetParametersMode(petsclib::$UnionPetscLib, term::TaoTerm )
	parameters_mode_ = Ref{TaoTermParametersMode}()

    @chk ccall(
               (:TaoTermGetParametersMode, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{TaoTermParametersMode}),
               term, parameters_mode_,
              )

	parameters_mode = parameters_mode_[]

	return parameters_mode
end 

"""
	k::PetscInt,M_K::PetscInt,bs::PetscInt = TaoTermGetParametersSizes(petsclib::PetscLibType, term::TaoTerm) 
Get the sizes describing the layout of the parameter vector space of a `TaoTerm`.

Not collective

Input Parameter:
- `term` - a `TaoTerm`

Output Parameters:
- `k`  - (optional) the size of a parameter vector on the current MPI process
- `K`  - (optional) the global size of a parameter vector
- `bs` - (optional) the block size of a parameter vector

Level: beginner

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermSetParametersSizes()`,
`TaoTermSetParametersTemplate()`,
`TaoTermGetParametersVecType()`,
`TaoTermSetParametersVecType()`,
`TaoTermGetParametersLayout()`,
`TaoTermSetParametersLayout()`,
`TaoTermCreateParametersVec()`

# External Links
$(_doc_external("TaoTerm/TaoTermGetParametersSizes"))
"""
function TaoTermGetParametersSizes(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermGetParametersSizes: no generated method for these argument types")
end

@for_petsc function TaoTermGetParametersSizes(petsclib::$UnionPetscLib, term::TaoTerm )
	k_ = Ref{$PetscInt}()
	M_K_ = Ref{$PetscInt}()
	bs_ = Ref{$PetscInt}()

    @chk ccall(
               (:TaoTermGetParametersSizes, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
               term, k_, M_K_, bs_,
              )

	k = k_[]
	M_K = M_K_[]
	bs = bs_[]

	return k,M_K,bs
end 

"""
	parameters_type::VecType = TaoTermGetParametersVecType(petsclib::PetscLibType, term::TaoTerm) 
Get the vector types of the parameter vector of a `TaoTerm`

Not collective

Input Parameter:
- `term` - a `TaoTerm`

Output Parameter:
- `parameters_type` - the `VecType` for the parameter space

Level: advanced

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermSetParametersVecType()`,
`TaoTermGetParametersLayout()`,
`TaoTermSetParametersLayout()`,
`TaoTermSetSolutionTemplate()`,
`TaoTermSetParametersTemplate()`,
`TaoTermCreateParametersVec()`

# External Links
$(_doc_external("TaoTerm/TaoTermGetParametersVecType"))
"""
function TaoTermGetParametersVecType(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermGetParametersVecType: no generated method for these argument types")
end

@for_petsc function TaoTermGetParametersVecType(petsclib::$UnionPetscLib, term::TaoTerm )
	parameters_type_ = Ref{VecType}()

    @chk ccall(
               (:TaoTermGetParametersVecType, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{VecType}),
               term, parameters_type_,
              )

	parameters_type = parameters_type_[] == C_NULL ? "" : unsafe_string(parameters_type_[])

	return parameters_type
end 

"""
	solution_layout::PetscLayout = TaoTermGetSolutionLayout(petsclib::PetscLibType, term::TaoTerm) 
Get the layouts describing the solution vectors of a `TaoTerm`.

Not collective

Input Parameter:
- `term` - a `TaoTerm`

Output Parameter:
- `solution_layout` - the `PetscLayout` for the solution space

Level: intermediate

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermGetSolutionVecType()`,
`TaoTermSetSolutionVecType()`,
`TaoTermSetSolutionLayout()`,
`TaoTermSetSolutionTemplate()`,
`TaoTermSetParametersTemplate()`,
`TaoTermCreateSolutionVec()`

# External Links
$(_doc_external("TaoTerm/TaoTermGetSolutionLayout"))
"""
function TaoTermGetSolutionLayout(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermGetSolutionLayout: no generated method for these argument types")
end

@for_petsc function TaoTermGetSolutionLayout(petsclib::$UnionPetscLib, term::TaoTerm )
	solution_layout_ = Ref{PetscLayout}()

    @chk ccall(
               (:TaoTermGetSolutionLayout, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{PetscLayout}),
               term, solution_layout_,
              )

	solution_layout = solution_layout_[]

	return solution_layout
end 

"""
	n::PetscInt,M_N::PetscInt,bs::PetscInt = TaoTermGetSolutionSizes(petsclib::PetscLibType, term::TaoTerm) 
Get the sizes describing the layout of the solution vector space of a `TaoTerm`.

Not collective

Input Parameter:
- `term` - a `TaoTerm`

Output Parameters:
- `n`  - (optional) the size of a solution vector on the current MPI process
- `N`  - (optional) the global size of a solution vector
- `bs` - (optional) the block size of a solution vector

Level: beginner

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermSetSolutionSizes()`,
`TaoTermSetSolutionTemplate()`,
`TaoTermGetSolutionVecType()`,
`TaoTermSetSolutionVecType()`,
`TaoTermGetSolutionLayout()`,
`TaoTermSetSolutionLayout()`,
`TaoTermCreateSolutionVec()`

# External Links
$(_doc_external("TaoTerm/TaoTermGetSolutionSizes"))
"""
function TaoTermGetSolutionSizes(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermGetSolutionSizes: no generated method for these argument types")
end

@for_petsc function TaoTermGetSolutionSizes(petsclib::$UnionPetscLib, term::TaoTerm )
	n_ = Ref{$PetscInt}()
	M_N_ = Ref{$PetscInt}()
	bs_ = Ref{$PetscInt}()

    @chk ccall(
               (:TaoTermGetSolutionSizes, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}),
               term, n_, M_N_, bs_,
              )

	n = n_[]
	M_N = M_N_[]
	bs = bs_[]

	return n,M_N,bs
end 

"""
	solution_type::VecType = TaoTermGetSolutionVecType(petsclib::PetscLibType, term::TaoTerm) 
Get the vector types of the solution vector of a `TaoTerm`

Not collective

Input Parameter:
- `term` - a `TaoTerm`

Output Parameter:
- `solution_type` - the `VecType` for the solution space

Level: advanced

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermSetSolutionVecType()`,
`TaoTermGetSolutionLayout()`,
`TaoTermSetSolutionLayout()`,
`TaoTermSetSolutionTemplate()`,
`TaoTermSetParametersTemplate()`,
`TaoTermCreateSolutionVec()`

# External Links
$(_doc_external("TaoTerm/TaoTermGetSolutionVecType"))
"""
function TaoTermGetSolutionVecType(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermGetSolutionVecType: no generated method for these argument types")
end

@for_petsc function TaoTermGetSolutionVecType(petsclib::$UnionPetscLib, term::TaoTerm )
	solution_type_ = Ref{VecType}()

    @chk ccall(
               (:TaoTermGetSolutionVecType, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{VecType}),
               term, solution_type_,
              )

	solution_type = solution_type_[] == C_NULL ? "" : unsafe_string(solution_type_[])

	return solution_type
end 

"""
	type::TaoTermType = TaoTermGetType(petsclib::PetscLibType, term::TaoTerm) 
Get the type of a `TaoTerm`

Not collective

Input Parameter:
- `term` - a `TaoTerm`

Output Parameter:
- `type` - the `TaoTermType`

Level: beginner

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermType`,
`TaoTermCreate()`,
`TaoTermSetType()`,
`TaoTermSetFromOptions()`,
`TaoTermSetUp()`,
`TaoTermView()`,
`TaoTermDestroy()`

# External Links
$(_doc_external("TaoTerm/TaoTermGetType"))
"""
function TaoTermGetType(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermGetType: no generated method for these argument types")
end

@for_petsc function TaoTermGetType(petsclib::$UnionPetscLib, term::TaoTerm )
	type_ = Ref{TaoTermType}()

    @chk ccall(
               (:TaoTermGetType, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{TaoTermType}),
               term, type_,
              )

	type = type_[] == C_NULL ? "" : unsafe_string(type_[])

	return type
end 

"""
	is_fdpossible::PetscBool3 = TaoTermIsComputeHessianFDPossible(petsclib::PetscLibType, term::TaoTerm) 
Whether this term can compute Hessian with finite differences
with either `-tao_term_hessian_use_fd`, `TaoTermComputeHessianSetUseFD()`, or `MATMFFD`.

Not collective

Input Parameter:
- `term` - a `TaoTerm`

Output Parameter:
- `is_fdpossible` - whether Hessian computation with finite differences is possible

Level: developer

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermComputeObjective()`,
`TaoTermShellSetObjective()`,
`TaoTermIsGradientDefined()`,
`TaoTermIsObjectiveAndGradientDefined()`,
`TaoTermIsHessianDefined()`

# External Links
$(_doc_external("TaoTerm/TaoTermIsComputeHessianFDPossible"))
"""
function TaoTermIsComputeHessianFDPossible(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermIsComputeHessianFDPossible: no generated method for these argument types")
end

@for_petsc function TaoTermIsComputeHessianFDPossible(petsclib::$UnionPetscLib, term::TaoTerm )
	is_fdpossible_ = Ref{PetscBool3}()

    @chk ccall(
               (:TaoTermIsComputeHessianFDPossible, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{PetscBool3}),
               term, is_fdpossible_,
              )

	is_fdpossible = is_fdpossible_[]

	return is_fdpossible
end 

"""
	is_defined::PetscBool = TaoTermIsCreateHessianMatricesDefined(petsclib::PetscLibType, term::TaoTerm) 
Whether this term can call `TaoTermCreateHessianMatrices()`.

Not collective

Input Parameter:
- `term` - a `TaoTerm`

Output Parameter:
- `is_defined` - whether the term can create new Hessian matrices

Level: developer

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermCreateHessianMatrices()`,
`TaoTermShellSetCreateHessianMatrices()`,
`TaoTermIsObjectiveDefined()`,
`TaoTermIsGradientDefined()`,
`TaoTermIsObjectiveAndGradientDefined()`,
`TaoTermIsHessianDefined()`

# External Links
$(_doc_external("TaoTerm/TaoTermIsCreateHessianMatricesDefined"))
"""
function TaoTermIsCreateHessianMatricesDefined(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermIsCreateHessianMatricesDefined: no generated method for these argument types")
end

@for_petsc function TaoTermIsCreateHessianMatricesDefined(petsclib::$UnionPetscLib, term::TaoTerm )
	is_defined_ = Ref{PetscBool}()

    @chk ccall(
               (:TaoTermIsCreateHessianMatricesDefined, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{PetscBool}),
               term, is_defined_,
              )

	is_defined = is_defined_[]

	return is_defined
end 

"""
	is_defined::PetscBool = TaoTermIsGradientDefined(petsclib::PetscLibType, term::TaoTerm) 
Whether a standalone gradient operation is defined for this `TaoTerm`

Not collective

Input Parameter:
- `term` - a `TaoTerm`

Output Parameter:
- `is_defined` - whether the gradient is defined

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermComputeGradient()`,
`TaoTermShellSetGradient()`,
`TaoTermIsObjectiveDefined()`,
`TaoTermIsObjectiveAndGradientDefined()`,
`TaoTermIsHessianDefined()`

# External Links
$(_doc_external("TaoTerm/TaoTermIsGradientDefined"))
"""
function TaoTermIsGradientDefined(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermIsGradientDefined: no generated method for these argument types")
end

@for_petsc function TaoTermIsGradientDefined(petsclib::$UnionPetscLib, term::TaoTerm )
	is_defined_ = Ref{PetscBool}()

    @chk ccall(
               (:TaoTermIsGradientDefined, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{PetscBool}),
               term, is_defined_,
              )

	is_defined = is_defined_[]

	return is_defined
end 

"""
	is_defined::PetscBool = TaoTermIsHessianDefined(petsclib::PetscLibType, term::TaoTerm) 
Whether a Hessian operation is defined for this `TaoTerm`

Not collective

Input Parameter:
- `term` - a `TaoTerm`

Output Parameter:
- `is_defined` - whether the Hessian is defined

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermComputeHessian()`,
`TaoTermShellSetHessian()`,
`TaoTermIsObjectiveDefined()`,
`TaoTermIsGradientDefined()`,
`TaoTermIsObjectiveAndGradientDefined()`

# External Links
$(_doc_external("TaoTerm/TaoTermIsHessianDefined"))
"""
function TaoTermIsHessianDefined(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermIsHessianDefined: no generated method for these argument types")
end

@for_petsc function TaoTermIsHessianDefined(petsclib::$UnionPetscLib, term::TaoTerm )
	is_defined_ = Ref{PetscBool}()

    @chk ccall(
               (:TaoTermIsHessianDefined, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{PetscBool}),
               term, is_defined_,
              )

	is_defined = is_defined_[]

	return is_defined
end 

"""
	is_defined::PetscBool = TaoTermIsObjectiveAndGradientDefined(petsclib::PetscLibType, term::TaoTerm) 
Whether a combined objective-and-gradient operation is defined for this `TaoTerm`

Not collective

Input Parameter:
- `term` - a `TaoTerm`

Output Parameter:
- `is_defined` - whether the objective/gradient is defined

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermComputeObjectiveAndGradient()`,
`TaoTermShellSetObjectiveAndGradient()`,
`TaoTermIsObjectiveDefined()`,
`TaoTermIsGradientDefined()`,
`TaoTermIsHessianDefined()`

# External Links
$(_doc_external("TaoTerm/TaoTermIsObjectiveAndGradientDefined"))
"""
function TaoTermIsObjectiveAndGradientDefined(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermIsObjectiveAndGradientDefined: no generated method for these argument types")
end

@for_petsc function TaoTermIsObjectiveAndGradientDefined(petsclib::$UnionPetscLib, term::TaoTerm )
	is_defined_ = Ref{PetscBool}()

    @chk ccall(
               (:TaoTermIsObjectiveAndGradientDefined, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{PetscBool}),
               term, is_defined_,
              )

	is_defined = is_defined_[]

	return is_defined
end 

"""
	is_defined::PetscBool = TaoTermIsObjectiveDefined(petsclib::PetscLibType, term::TaoTerm) 
Whether a standalone objective operation is defined for this `TaoTerm`

Not collective

Input Parameter:
- `term` - a `TaoTerm`

Output Parameter:
- `is_defined` - whether the objective is defined

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermComputeObjective()`,
`TaoTermShellSetObjective()`,
`TaoTermIsGradientDefined()`,
`TaoTermIsObjectiveAndGradientDefined()`,
`TaoTermIsHessianDefined()`

# External Links
$(_doc_external("TaoTerm/TaoTermIsObjectiveDefined"))
"""
function TaoTermIsObjectiveDefined(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermIsObjectiveDefined: no generated method for these argument types")
end

@for_petsc function TaoTermIsObjectiveDefined(petsclib::$UnionPetscLib, term::TaoTerm )
	is_defined_ = Ref{PetscBool}()

    @chk ccall(
               (:TaoTermIsObjectiveDefined, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{PetscBool}),
               term, is_defined_,
              )

	is_defined = is_defined_[]

	return is_defined
end 

"""
	epsilon::PetscReal = TaoTermL1GetEpsilon(petsclib::PetscLibType, term::TaoTerm) 
Get the \\epsilon smoothing parameter set by `TaoTermL1SetEpsilon()`.

Not collective

Input Parameter:
- `term` - a `TaoTerm` of type `TAOTERML1`

Output Parameter:
- `epsilon` - the smoothing parameter

Level: advanced

See also: [](sec_tao_term),
`TaoTerm`,
`TAOTERML1`,
`TaoTermL1SetEpsilon()`

# External Links
$(_doc_external("TaoTerm/TaoTermL1GetEpsilon"))
"""
function TaoTermL1GetEpsilon(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermL1GetEpsilon: no generated method for these argument types")
end

@for_petsc function TaoTermL1GetEpsilon(petsclib::$UnionPetscLib, term::TaoTerm )
	epsilon_ = Ref{$PetscReal}()

    @chk ccall(
               (:TaoTermL1GetEpsilon, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{$PetscReal}),
               term, epsilon_,
              )

	epsilon = epsilon_[]

	return epsilon
end 

"""
	TaoTermL1SetEpsilon(petsclib::PetscLibType, term::TaoTerm, epsilon::PetscReal) 
Set an \\epsilon smoothing parameter.

Logically collective

Input Parameters:
- `term`    - a `TaoTerm` of type `TAOTERML1`
- `epsilon` - a real number \\geq 0

Options Database Keys:
- `-tao_term_l1_epsilon <real>` - \\epsilon

Level: advanced

If \\epsilon = 0 (the default), then `term` computes \\|x - p\\|_1, but if \\epsilon > 0, then it computes
\\sum_{i=0}^{n-1} \\left(\\sqrt{(x_i-p_i)^2 + \\epsilon^2} - \\epsilon\\right).

See also: [](sec_tao_term),
`TaoTerm`,
`TAOTERML1`,
`TaoTermL1GetEpsilon()`

# External Links
$(_doc_external("TaoTerm/TaoTermL1SetEpsilon"))
"""
function TaoTermL1SetEpsilon(petsclib::PetscLibType, term::TaoTerm, epsilon::Real)
    error("TaoTermL1SetEpsilon: no generated method for these argument types")
end

@for_petsc function TaoTermL1SetEpsilon(petsclib::$UnionPetscLib, term::TaoTerm, epsilon::$PetscReal )

    @chk ccall(
               (:TaoTermL1SetEpsilon, $petsc_library),
               PetscErrorCode,
               (TaoTerm, $PetscReal),
               term, epsilon,
              )


	return nothing
end 

"""
	A::PetscMat = TaoTermQuadraticGetMat(petsclib::PetscLibType, term::TaoTerm) 
Get the matrix defining a `TaoTerm` of type `TAOTERMQUADRATIC`

Not collective

Input Parameter:
- `term` - a `TaoTerm` of type `TAOTERMQUADRATIC`

Output Parameter:
- `A` - the matrix

Level: intermediate

See also: [](sec_tao_term),
`TaoTerm`,
`TAOTERMQUADRATIC`,
`TaoTermQuadraticSetMat()`

# External Links
$(_doc_external("TaoTerm/TaoTermQuadraticGetMat"))
"""
function TaoTermQuadraticGetMat(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermQuadraticGetMat: no generated method for these argument types")
end

@for_petsc function TaoTermQuadraticGetMat(petsclib::$UnionPetscLib, term::TaoTerm )
	A_ = Ref{CMat}()

    @chk ccall(
               (:TaoTermQuadraticGetMat, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{CMat}),
               term, A_,
              )

	A = PetscMat(A_[], petsclib)

	return A
end 

"""
	TaoTermQuadraticSetMat(petsclib::PetscLibType, term::TaoTerm, A::AbstractPetscMat) 
Set the matrix defining a `TaoTerm` of type `TAOTERMQUADRATIC`

Collective

Input Parameters:
- `term` - a `TaoTerm` of type `TAOTERMQUADRATIC`
- `A`    - the matrix

Level: intermediate

See also: [](sec_tao_term),
`TaoTerm`,
`TAOTERMQUADRATIC`,
`TaoTermQuadraticGetMat()`

# External Links
$(_doc_external("TaoTerm/TaoTermQuadraticSetMat"))
"""
function TaoTermQuadraticSetMat(petsclib::PetscLibType, term::TaoTerm, A::AbstractPetscMat)
    error("TaoTermQuadraticSetMat: no generated method for these argument types")
end

@for_petsc function TaoTermQuadraticSetMat(petsclib::$UnionPetscLib, term::TaoTerm, A::AbstractPetscMat )

    @chk ccall(
               (:TaoTermQuadraticSetMat, $petsc_library),
               PetscErrorCode,
               (TaoTerm, CMat),
               term, A,
              )


	return nothing
end 

"""
	TaoTermRegister(petsclib::PetscLibType, sname::String, func::external) 
Register an implementation of `TaoTerm`

Not Collective, No Fortran Support

Input Parameters:
- `sname` - name of a new user-defined term
- `func`  - routine to create the context for the `TaoTermType`

See also: [](sec_tao_term), `TaoTerm`, `TaoTermSetType()`

# External Links
$(_doc_external("TaoTerm/TaoTermRegister"))
"""
function TaoTermRegister(petsclib::PetscLibType, sname::String, func::external)
    error("TaoTermRegister: no generated method for these argument types")
end

@for_petsc function TaoTermRegister(petsclib::$UnionPetscLib, sname::String, func::external )

    @chk ccall(
               (:TaoTermRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, func,
              )


	return nothing
end 

"""
	TaoTermSetCreateHessianMode(petsclib::PetscLibType, term::TaoTerm, Hpre_is_H::PetscBool, H_mattype::MatType, Hpre_mattype::MatType) 
Determine the behavior of `TaoTermCreateHessianMatricesDefault()`.

Logically collective

Input Parameters:
- `term`         - a `TaoTerm`
- `Hpre_is_H`    - should `TaoTermCreateHessianMatricesDefault()` make one matrix for `H` and `Hpre`?
- `H_mattype`    - the `MatType` to create for `H`
- `Hpre_mattype` - the `MatType` to create for `Hpre`

Options Database Keys:
- `-tao_term_hessian_pre_is_hessian <bool>` - Whether `TaoTermCreateHessianMatrices()` should make a separate matrix for constructing the preconditioner
- `-tao_term_hessian_mat_type <type>`       - `MatType` for Hessian matrix created by `TaoTermCreateHessianMatrices()`
- `-tao_term_hessian_pre_mat_type <type>`   - `MatType` for matrix from which a preconditioner can be created by `TaoTermCreateHessianMatrices()`

Level: developer

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermComputeHessian()`,
`TaoTermCreateHessianMatrices()`,
`TaoTermCreateHessianMatricesDefault()`,
`TaoTermGetCreateHessianMode()`

# External Links
$(_doc_external("TaoTerm/TaoTermSetCreateHessianMode"))
"""
function TaoTermSetCreateHessianMode(petsclib::PetscLibType, term::TaoTerm, Hpre_is_H::PetscBool, H_mattype::MatType, Hpre_mattype::MatType)
    error("TaoTermSetCreateHessianMode: no generated method for these argument types")
end

@for_petsc function TaoTermSetCreateHessianMode(petsclib::$UnionPetscLib, term::TaoTerm, Hpre_is_H::PetscBool, H_mattype::MatType, Hpre_mattype::MatType )

    @chk ccall(
               (:TaoTermSetCreateHessianMode, $petsc_library),
               PetscErrorCode,
               (TaoTerm, PetscBool, MatType, MatType),
               term, Hpre_is_H, H_mattype, Hpre_mattype,
              )


	return nothing
end 

"""
	TaoTermSetFDDelta(petsclib::PetscLibType, term::TaoTerm, delta::PetscReal) 
Set the increment used for finite difference derivative approximations in methods like `TaoTermComputeGradientFD()`

Logically collective

Input Parameters:
- `term`  - a `TaoTerm`
- `delta` - the finite difference increment

Options Database Key:
- `-tao_term_fd_delta <delta>` - the above increment

Level: advanced

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermGetFDDelta()`,
`TaoTermComputeGradientFD()`,
`TaoTermComputeGradientSetUseFD()`,
`TaoTermComputeGradientGetUseFD()`,
`TaoTermComputeHessianFD()`,
`TaoTermComputeHessianSetUseFD()`,
`TaoTermComputeHessianGetUseFD()`

# External Links
$(_doc_external("TaoTerm/TaoTermSetFDDelta"))
"""
function TaoTermSetFDDelta(petsclib::PetscLibType, term::TaoTerm, delta::Real)
    error("TaoTermSetFDDelta: no generated method for these argument types")
end

@for_petsc function TaoTermSetFDDelta(petsclib::$UnionPetscLib, term::TaoTerm, delta::$PetscReal )

    @chk ccall(
               (:TaoTermSetFDDelta, $petsc_library),
               PetscErrorCode,
               (TaoTerm, $PetscReal),
               term, delta,
              )


	return nothing
end 

"""
	TaoTermSetFromOptions(petsclib::PetscLibType, term::TaoTerm) 
Configure a `TaoTerm` from the PETSc options database

Collective

Input Parameter:
- `term` - a `TaoTerm`

Options Database Keys:
- `-tao_term_type <type>`                              - l1, halfl2squared; see `TaoTermType` for a complete list
- `-tao_term_solution_vec_type <type>`                 - the type of vector to use for the solution, see `VecType` for a complete list of vector types
- `-tao_term_parameters_vec_type <type>`               - the type of vector to use for the parameters, see `VecType` for a complete list of vector types
- `-tao_term_parameters_mode <optional,none,required>` - `TAOTERM_PARAMETERS_OPTIONAL`, `TAOTERM_PARAMETERS_NONE`, `TAOTERM_PARAMETERS_REQUIRED`
- `-tao_term_hessian_pre_is_hessian <bool>`            - Whether `TaoTermCreateHessianMatricesDefault()` should make a separate preconditioning matrix
- `-tao_term_hessian_mat_type <type>`                  - `MatType` for Hessian matrix created by `TaoTermCreateHessianMatricesDefault()`
- `-tao_term_hessian_pre_mat_type <type>`              - `MatType` for approximate Hessian matrix used to construct the preconditioner created by `TaoTermCreateHessianMatricesDefault()`
- `-tao_term_fd_delta <real>`                          - Increment for finite difference derivative approximations in `TaoTermComputeGradientFD()`
- `-tao_term_gradient_use_fd <bool>`                   - Use finite differences in `TaoTermComputeGradient()`, overriding other user-provided or built-in routines
- `-tao_term_hessian_use_fd <bool>`                    - Use finite differences in `TaoTermComputeHessian()`, overriding other user-provided or built-in routines

Level: beginner

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermCreate()`,
`TaoTermSetType()`,
`TaoTermSetUp()`,
`TaoTermView()`,
`TaoTermDestroy()`

# External Links
$(_doc_external("TaoTerm/TaoTermSetFromOptions"))
"""
function TaoTermSetFromOptions(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermSetFromOptions: no generated method for these argument types")
end

@for_petsc function TaoTermSetFromOptions(petsclib::$UnionPetscLib, term::TaoTerm )

    @chk ccall(
               (:TaoTermSetFromOptions, $petsc_library),
               PetscErrorCode,
               (TaoTerm,),
               term,
              )


	return nothing
end 

"""
	TaoTermSetParametersLayout(petsclib::PetscLibType, term::TaoTerm, parameters_layout::PetscLayout) 
Set the layout describing the parameter vector of `TaoTerm`.

Collective

Input Parameters:
- `term`              - a `TaoTerm`
- `parameters_layout` - the `PetscLayout` for the parameter space

Level: intermediate

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermGetParametersVecType()`,
`TaoTermSetParametersVecType()`,
`TaoTermGetParametersLayout()`,
`TaoTermSetSolutionTemplate()`,
`TaoTermSetParametersTemplate()`,
`TaoTermCreateParametersVec()`

# External Links
$(_doc_external("TaoTerm/TaoTermSetParametersLayout"))
"""
function TaoTermSetParametersLayout(petsclib::PetscLibType, term::TaoTerm, parameters_layout::PetscLayout)
    error("TaoTermSetParametersLayout: no generated method for these argument types")
end

@for_petsc function TaoTermSetParametersLayout(petsclib::$UnionPetscLib, term::TaoTerm, parameters_layout::PetscLayout )

    @chk ccall(
               (:TaoTermSetParametersLayout, $petsc_library),
               PetscErrorCode,
               (TaoTerm, PetscLayout),
               term, parameters_layout,
              )


	return nothing
end 

"""
	TaoTermSetParametersMode(petsclib::PetscLibType, term::TaoTerm, parameters_mode::TaoTermParametersMode) 
Sets the way a `TaoTerm` can accept parameters

Logically collective

Input Parameters:
- `term`            - a `TaoTerm`
- `parameters_mode` - `TAOTERM_PARAMETERS_OPTIONAL`, `TAOTERM_PARAMETERS_NONE`, `TAOTERM_PARAMETERS_REQUIRED`

Options Database Keys:
- `-tao_term_parameters_mode <optional,none,required>` - `TAOTERM_PARAMETERS_OPTIONAL`, `TAOTERM_PARAMETERS_NONE`, `TAOTERM_PARAMETERS_REQUIRED`

Level: advanced

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermParametersMode`,
`TaoTermGetParametersMode()`

# External Links
$(_doc_external("TaoTerm/TaoTermSetParametersMode"))
"""
function TaoTermSetParametersMode(petsclib::PetscLibType, term::TaoTerm, parameters_mode::TaoTermParametersMode)
    error("TaoTermSetParametersMode: no generated method for these argument types")
end

@for_petsc function TaoTermSetParametersMode(petsclib::$UnionPetscLib, term::TaoTerm, parameters_mode::TaoTermParametersMode )

    @chk ccall(
               (:TaoTermSetParametersMode, $petsc_library),
               PetscErrorCode,
               (TaoTerm, TaoTermParametersMode),
               term, parameters_mode,
              )


	return nothing
end 

"""
	TaoTermSetParametersSizes(petsclib::PetscLibType, term::TaoTerm, k::PetscInt, M_K::PetscInt, bs::PetscInt) 
Set the sizes describing the layout of the parameter vector space of a `TaoTerm`.

Logically collective

Input Parameters:
- `term` - a `TaoTerm`
- `k`    - the size of a parameter vector on the current MPI process (or `PETSC_DECIDE`)
- `K`    - the global size of a parameter vector (or `PETSC_DECIDE`)
- `bs`   - the block size of a parameter vector (must be >= 1)

Level: beginner

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermGetParametersSizes()`,
`TaoTermSetParametersTemplate()`,
`TaoTermGetParametersVecType()`,
`TaoTermSetParametersVecType()`,
`TaoTermGetParametersLayout()`,
`TaoTermSetParametersLayout()`,
`TaoTermCreateParametersVec()`

# External Links
$(_doc_external("TaoTerm/TaoTermSetParametersSizes"))
"""
function TaoTermSetParametersSizes(petsclib::PetscLibType, term::TaoTerm, k::Integer, M_K::Integer, bs::Integer)
    error("TaoTermSetParametersSizes: no generated method for these argument types")
end

@for_petsc function TaoTermSetParametersSizes(petsclib::$UnionPetscLib, term::TaoTerm, k::$PetscInt, M_K::$PetscInt, bs::$PetscInt )

    @chk ccall(
               (:TaoTermSetParametersSizes, $petsc_library),
               PetscErrorCode,
               (TaoTerm, $PetscInt, $PetscInt, $PetscInt),
               term, k, M_K, bs,
              )


	return nothing
end 

"""
	TaoTermSetParametersTemplate(petsclib::PetscLibType, term::TaoTerm, params_template::AbstractPetscVec) 
Set the parameter vector space to match a template vector

Collective

Input Parameters:
- `term`            - a `TaoTerm`
- `params_template` - a vector with the desired size, layout, and `VecType` of parameter vectors for `TaoTerm`

Level: intermediate

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermGetParametersVecType()`,
`TaoTermSetParametersVecType()`,
`TaoTermSetSolutionTemplate()`,
`TaoTermGetParametersLayout()`,
`TaoTermSetParametersLayout()`,
`TaoTermCreateSolutionVec()`

# External Links
$(_doc_external("TaoTerm/TaoTermSetParametersTemplate"))
"""
function TaoTermSetParametersTemplate(petsclib::PetscLibType, term::TaoTerm, params_template::AbstractPetscVec)
    error("TaoTermSetParametersTemplate: no generated method for these argument types")
end

@for_petsc function TaoTermSetParametersTemplate(petsclib::$UnionPetscLib, term::TaoTerm, params_template::AbstractPetscVec )

    @chk ccall(
               (:TaoTermSetParametersTemplate, $petsc_library),
               PetscErrorCode,
               (TaoTerm, CVec),
               term, params_template,
              )


	return nothing
end 

"""
	TaoTermSetParametersVecType(petsclib::PetscLibType, term::TaoTerm, parameters_type::VecType) 
Set the vector types of the parameters vector of a `TaoTerm`

Logically collective

Input Parameters:
- `term`            - a `TaoTerm`
- `parameters_type` - the `VecType` for the parameters space

Options Database Keys:
- `-tao_term_parameters_vec_type <type>` - `VecType` for complete list of vector types

Level: advanced

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermGetParametersVecType()`,
`TaoTermSetParametersLayout()`,
`TaoTermGetParametersLayout()`,
`TaoTermSetParametersTemplate()`,
`TaoTermCreateParametersVec()`

# External Links
$(_doc_external("TaoTerm/TaoTermSetParametersVecType"))
"""
function TaoTermSetParametersVecType(petsclib::PetscLibType, term::TaoTerm, parameters_type::VecType)
    error("TaoTermSetParametersVecType: no generated method for these argument types")
end

@for_petsc function TaoTermSetParametersVecType(petsclib::$UnionPetscLib, term::TaoTerm, parameters_type::VecType )

    @chk ccall(
               (:TaoTermSetParametersVecType, $petsc_library),
               PetscErrorCode,
               (TaoTerm, VecType),
               term, parameters_type,
              )


	return nothing
end 

"""
	TaoTermSetSolutionLayout(petsclib::PetscLibType, term::TaoTerm, solution_layout::PetscLayout) 
Set the layout describing the solution vector of `TaoTerm`.

Collective

Input Parameters:
- `term`            - a `TaoTerm`
- `solution_layout` - the `PetscLayout` for the solution space

Level: intermediate

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermGetSolutionVecType()`,
`TaoTermSetSolutionVecType()`,
`TaoTermGetSolutionLayout()`,
`TaoTermSetSolutionTemplate()`,
`TaoTermSetParametersTemplate()`,
`TaoTermCreateSolutionVec()`

# External Links
$(_doc_external("TaoTerm/TaoTermSetSolutionLayout"))
"""
function TaoTermSetSolutionLayout(petsclib::PetscLibType, term::TaoTerm, solution_layout::PetscLayout)
    error("TaoTermSetSolutionLayout: no generated method for these argument types")
end

@for_petsc function TaoTermSetSolutionLayout(petsclib::$UnionPetscLib, term::TaoTerm, solution_layout::PetscLayout )

    @chk ccall(
               (:TaoTermSetSolutionLayout, $petsc_library),
               PetscErrorCode,
               (TaoTerm, PetscLayout),
               term, solution_layout,
              )


	return nothing
end 

"""
	TaoTermSetSolutionSizes(petsclib::PetscLibType, term::TaoTerm, n::PetscInt, M_N::PetscInt, bs::PetscInt) 
Set the sizes describing the layout of the solution vector space of a `TaoTerm`.

Logically collective

Input Parameters:
- `term` - a `TaoTerm`
- `n`    - the size of a solution vector on the current MPI process (or `PETSC_DECIDE`)
- `N`    - the global size of a solution vector (or `PETSC_DECIDE`)
- `bs`   - the block size of a solution vector (must be >= 1)

Level: beginner

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermGetSolutionSizes()`,
`TaoTermSetSolutionTemplate()`,
`TaoTermGetSolutionVecType()`,
`TaoTermSetSolutionVecType()`,
`TaoTermGetSolutionLayout()`,
`TaoTermSetSolutionLayout()`,
`TaoTermCreateSolutionVec()`

# External Links
$(_doc_external("TaoTerm/TaoTermSetSolutionSizes"))
"""
function TaoTermSetSolutionSizes(petsclib::PetscLibType, term::TaoTerm, n::Integer, M_N::Integer, bs::Integer)
    error("TaoTermSetSolutionSizes: no generated method for these argument types")
end

@for_petsc function TaoTermSetSolutionSizes(petsclib::$UnionPetscLib, term::TaoTerm, n::$PetscInt, M_N::$PetscInt, bs::$PetscInt )

    @chk ccall(
               (:TaoTermSetSolutionSizes, $petsc_library),
               PetscErrorCode,
               (TaoTerm, $PetscInt, $PetscInt, $PetscInt),
               term, n, M_N, bs,
              )


	return nothing
end 

"""
	TaoTermSetSolutionTemplate(petsclib::PetscLibType, term::TaoTerm, sol_template::AbstractPetscVec) 
Set the solution vector space to match a template vector

Collective

Input Parameters:
- `term`         - a `TaoTerm`
- `sol_template` - a vector with the desired size, layout, and `VecType` of solution vectors for `TaoTerm`

Level: intermediate

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermGetSolutionVecType()`,
`TaoTermSetSolutionVecType()`,
`TaoTermSetParametersTemplate()`,
`TaoTermGetSolutionLayout()`,
`TaoTermSetSolutionLayout()`,
`TaoTermCreateSolutionVec()`

# External Links
$(_doc_external("TaoTerm/TaoTermSetSolutionTemplate"))
"""
function TaoTermSetSolutionTemplate(petsclib::PetscLibType, term::TaoTerm, sol_template::AbstractPetscVec)
    error("TaoTermSetSolutionTemplate: no generated method for these argument types")
end

@for_petsc function TaoTermSetSolutionTemplate(petsclib::$UnionPetscLib, term::TaoTerm, sol_template::AbstractPetscVec )

    @chk ccall(
               (:TaoTermSetSolutionTemplate, $petsc_library),
               PetscErrorCode,
               (TaoTerm, CVec),
               term, sol_template,
              )


	return nothing
end 

"""
	TaoTermSetSolutionVecType(petsclib::PetscLibType, term::TaoTerm, solution_type::VecType) 
Set the vector types of the solution vector of a `TaoTerm`

Logically collective

Input Parameters:
- `term`          - a `TaoTerm`
- `solution_type` - the `VecType` for the solution space

Options Database Keys:
- `-tao_term_solution_vec_type <type>` - `VecType` for complete list of vector types

Level: advanced

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermGetSolutionVecType()`,
`TaoTermSetSolutionLayout()`,
`TaoTermGetSolutionLayout()`,
`TaoTermSetSolutionTemplate()`,
`TaoTermSetParametersTemplate()`,
`TaoTermCreateSolutionVec()`

# External Links
$(_doc_external("TaoTerm/TaoTermSetSolutionVecType"))
"""
function TaoTermSetSolutionVecType(petsclib::PetscLibType, term::TaoTerm, solution_type::VecType)
    error("TaoTermSetSolutionVecType: no generated method for these argument types")
end

@for_petsc function TaoTermSetSolutionVecType(petsclib::$UnionPetscLib, term::TaoTerm, solution_type::VecType )

    @chk ccall(
               (:TaoTermSetSolutionVecType, $petsc_library),
               PetscErrorCode,
               (TaoTerm, VecType),
               term, solution_type,
              )


	return nothing
end 

"""
	TaoTermSetType(petsclib::PetscLibType, term::TaoTerm, type::TaoTermType) 
Set the type of a `TaoTerm`

Collective

Input Parameters:
- `term` - a `TaoTerm`
- `type` - a `TaoTermType`

Options Database Keys:
- `-tao_term_type <type>` - l1, halfl2squared, `TaoTermType` for complete list

Level: beginner

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermType`,
`TaoTermCreate()`,
`TaoTermGetType()`,
`TaoTermSetFromOptions()`,
`TaoTermSetUp()`,
`TaoTermView()`,
`TaoTermDestroy()`

# External Links
$(_doc_external("TaoTerm/TaoTermSetType"))
"""
function TaoTermSetType(petsclib::PetscLibType, term::TaoTerm, type::TaoTermType)
    error("TaoTermSetType: no generated method for these argument types")
end

@for_petsc function TaoTermSetType(petsclib::$UnionPetscLib, term::TaoTerm, type::TaoTermType )

    @chk ccall(
               (:TaoTermSetType, $petsc_library),
               PetscErrorCode,
               (TaoTerm, TaoTermType),
               term, type,
              )


	return nothing
end 

"""
	TaoTermSetUp(petsclib::PetscLibType, term::TaoTerm) 
Set up a `TaoTerm`.

Collective

Input Parameter:
- `term` - a `TaoTerm`

Level: intermediate

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermCreate()`,
`TaoTermSetType()`,
`TaoTermSetFromOptions()`,
`TaoTermView()`,
`TaoTermDestroy()`

# External Links
$(_doc_external("TaoTerm/TaoTermSetUp"))
"""
function TaoTermSetUp(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermSetUp: no generated method for these argument types")
end

@for_petsc function TaoTermSetUp(petsclib::$UnionPetscLib, term::TaoTerm )

    @chk ccall(
               (:TaoTermSetUp, $petsc_library),
               PetscErrorCode,
               (TaoTerm,),
               term,
              )


	return nothing
end 

"""
	ctx::Ptr{Cvoid} = TaoTermShellGetContext(petsclib::PetscLibType, term::TaoTerm) 
Get the context for a `TAOTERMSHELL`

Not collective

Input Parameter:
- `term` - a `TaoTerm` of type `TAOTERMSHELL`

Output Parameter:
- `ctx` - a context

Level: intermediate

See also: [](sec_tao_term), `TaoTerm`, `TAOTERMSHELL`, `TaoTermShellSetContext()`, `TaoTermShellSetContextDestroy()`

# External Links
$(_doc_external("TaoTerm/TaoTermShellGetContext"))
"""
function TaoTermShellGetContext(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermShellGetContext: no generated method for these argument types")
end

@for_petsc function TaoTermShellGetContext(petsclib::$UnionPetscLib, term::TaoTerm )
	ctx_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:TaoTermShellGetContext, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{Cvoid}),
               term, ctx_,
              )

	ctx = ctx_[]

	return ctx
end 

"""
	TaoTermShellSetContext(petsclib::PetscLibType, term::TaoTerm, ctx::Ptr{Cvoid}) 
Set a context for a `TAOTERMSHELL`

Logically collective

Input Parameters:
- `term` - a `TaoTerm` of type `TAOTERMSHELL`
- `ctx`  - a context

Level: intermediate

See also: [](sec_tao_term), `TaoTerm`, `TAOTERMSHELL`, `TaoTermShellGetContext()`, `TaoTermShellSetContextDestroy()`

# External Links
$(_doc_external("TaoTerm/TaoTermShellSetContext"))
"""
function TaoTermShellSetContext(petsclib::PetscLibType, term::TaoTerm, ctx::Ptr{Cvoid})
    error("TaoTermShellSetContext: no generated method for these argument types")
end

@for_petsc function TaoTermShellSetContext(petsclib::$UnionPetscLib, term::TaoTerm, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:TaoTermShellSetContext, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{Cvoid}),
               term, ctx,
              )


	return nothing
end 

"""
	TaoTermShellSetContextDestroy(petsclib::PetscLibType, term::TaoTerm, destroy::Ptr{Cvoid}) 
Set a method to destroy the context resources when a `TAOTERMSHELL` is destroyed

Logically collective

Input Parameters:
- `term`    - a `TaoTerm` of type `TAOTERMSHELL`
- `destroy` - the context destroy function

Level: intermediate

See also: [](sec_tao_term), `TaoTerm`, `TAOTERMSHELL`, `TaoTermShellSetContext()`, `TaoTermShellGetContext()`

# External Links
$(_doc_external("TaoTerm/TaoTermShellSetContextDestroy"))
"""
function TaoTermShellSetContextDestroy(petsclib::PetscLibType, term::TaoTerm, destroy::Ptr{Cvoid})
    error("TaoTermShellSetContextDestroy: no generated method for these argument types")
end

@for_petsc function TaoTermShellSetContextDestroy(petsclib::$UnionPetscLib, term::TaoTerm, destroy::Ptr{Cvoid} )

    @chk ccall(
               (:TaoTermShellSetContextDestroy, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{Cvoid}),
               term, destroy,
              )


	return nothing
end 

"""
	TaoTermShellSetCreateHessianMatrices(petsclib::PetscLibType, term::TaoTerm, createmats::external) 
Set the routine that creates Hessian matrices for a `TaoTerm` of type `TAOTERMSHELL`

Logically collective

Input Parameters:
- `term`       - a `TaoTerm` of type `TAOTERMSHELL`
- `createmats` - a function with the same signature as `TaoTermCreateHessianMatrices()`

Calling sequence of `createmats`:
- `f`    - the `TaoTerm`
- `H`    - (optional) a matrix of the appropriate type and size for the Hessian of `term`
- `Hpre` - (optional) a matrix of the appropriate type and size for constructing a preconditioner for the Hessian of `term`

Level: intermediate

See also: [](sec_tao_term), `TaoTerm`, `TAOTERMSHELL`, `TaoTermShellGetContext()`, `TaoTermShellSetContextDestroy()`,
`TaoTermShellSetCreateSolutionVec()`, `TaoTermShellSetCreateParametersVec()`

# External Links
$(_doc_external("TaoTerm/TaoTermShellSetCreateHessianMatrices"))
"""
function TaoTermShellSetCreateHessianMatrices(petsclib::PetscLibType, term::TaoTerm, createmats::external)
    error("TaoTermShellSetCreateHessianMatrices: no generated method for these argument types")
end

@for_petsc function TaoTermShellSetCreateHessianMatrices(petsclib::$UnionPetscLib, term::TaoTerm, createmats::external )

    @chk ccall(
               (:TaoTermShellSetCreateHessianMatrices, $petsc_library),
               PetscErrorCode,
               (TaoTerm, external),
               term, createmats,
              )


	return nothing
end 

"""
	TaoTermShellSetCreateParametersVec(petsclib::PetscLibType, term::TaoTerm, createparametersvec::external) 
Set the routine that creates parameters vector for a `TaoTerm` of type `TAOTERMSHELL`

Logically collective

Input Parameters:
- `term`                - a `TaoTerm` of type `TAOTERMSHELL`
- `createparametersvec` - a function with the same signature as `TaoTermCreateParametersVec()`

Calling sequence of `createparametersvec`:
- `term`       - the `TaoTerm`
- `parameters` - a parameters vector for `term`

Level: intermediate

See also: [](sec_tao_term), `TaoTerm`, `TAOTERMSHELL`, `TaoTermShellGetContext()`, `TaoTermShellSetContextDestroy()`,
`TaoTermShellSetCreateHessianMatrices()`

# External Links
$(_doc_external("TaoTerm/TaoTermShellSetCreateParametersVec"))
"""
function TaoTermShellSetCreateParametersVec(petsclib::PetscLibType, term::TaoTerm, createparametersvec::external)
    error("TaoTermShellSetCreateParametersVec: no generated method for these argument types")
end

@for_petsc function TaoTermShellSetCreateParametersVec(petsclib::$UnionPetscLib, term::TaoTerm, createparametersvec::external )

    @chk ccall(
               (:TaoTermShellSetCreateParametersVec, $petsc_library),
               PetscErrorCode,
               (TaoTerm, external),
               term, createparametersvec,
              )


	return nothing
end 

"""
	TaoTermShellSetCreateSolutionVec(petsclib::PetscLibType, term::TaoTerm, createsolutionvec::external) 
Set the routine that creates solution vector for a `TaoTerm` of type `TAOTERMSHELL`

Logically collective

Input Parameters:
- `term`              - a `TaoTerm` of type `TAOTERMSHELL`
- `createsolutionvec` - a function with the same signature as `TaoTermCreateSolutionVec()`

Calling sequence of `createsolutionvec`:
- `term`     - the `TaoTerm`
- `solution` - a solution vector for `term`

Level: intermediate

See also: [](sec_tao_term), `TaoTerm`, `TAOTERMSHELL`, `TaoTermShellGetContext()`, `TaoTermShellSetContextDestroy()`,
`TaoTermShellSetCreateHessianMatrices()`

# External Links
$(_doc_external("TaoTerm/TaoTermShellSetCreateSolutionVec"))
"""
function TaoTermShellSetCreateSolutionVec(petsclib::PetscLibType, term::TaoTerm, createsolutionvec::external)
    error("TaoTermShellSetCreateSolutionVec: no generated method for these argument types")
end

@for_petsc function TaoTermShellSetCreateSolutionVec(petsclib::$UnionPetscLib, term::TaoTerm, createsolutionvec::external )

    @chk ccall(
               (:TaoTermShellSetCreateSolutionVec, $petsc_library),
               PetscErrorCode,
               (TaoTerm, external),
               term, createsolutionvec,
              )


	return nothing
end 

"""
	TaoTermShellSetGradient(petsclib::PetscLibType, term::TaoTerm, gradient::Ptr{Cvoid}) 
Set the gradient function of a `TAOTERMSHELL`

Logically collective

Input Parameters:
- `term`     - a `TaoTerm` of type `TAOTERMSHELL`
- `gradient` - a `TaoTermGradientFn` function pointer

Level: intermediate

See also: [](sec_tao_term), `TaoTerm`, `TAOTERMSHELL`, `TaoTermShellGetContext()`, `TaoTermShellSetContextDestroy()`,
`TaoTermShellSetObjective()`,
`TaoTermShellSetObjectiveAndGradient()`,
`TaoTermShellSetHessian()`,
`TaoTermShellSetView()`,
`TaoTermGradientFn`

# External Links
$(_doc_external("TaoTerm/TaoTermShellSetGradient"))
"""
function TaoTermShellSetGradient(petsclib::PetscLibType, term::TaoTerm, gradient::Ptr{Cvoid})
    error("TaoTermShellSetGradient: no generated method for these argument types")
end

@for_petsc function TaoTermShellSetGradient(petsclib::$UnionPetscLib, term::TaoTerm, gradient::Ptr{Cvoid} )

    @chk ccall(
               (:TaoTermShellSetGradient, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{Cvoid}),
               term, gradient,
              )


	return nothing
end 

"""
	TaoTermShellSetHessian(petsclib::PetscLibType, term::TaoTerm, hessian::Ptr{Cvoid}) 
Set the Hessian function of a `TAOTERMSHELL`

Logically collective

Input Parameters:
- `term`    - a `TaoTerm` of type `TAOTERMSHELL`
- `hessian` - a `TaoTermHessianFn` function pointer

Level: intermediate

See also: [](sec_tao_term), `TaoTerm`, `TAOTERMSHELL`, `TaoTermShellGetContext()`, `TaoTermShellSetContextDestroy()`,
`TaoTermShellSetObjective()`,
`TaoTermShellSetGradient()`,
`TaoTermShellSetObjectiveAndGradient()`,
`TaoTermShellSetView()`,
`TaoTermHessianFn`

# External Links
$(_doc_external("TaoTerm/TaoTermShellSetHessian"))
"""
function TaoTermShellSetHessian(petsclib::PetscLibType, term::TaoTerm, hessian::Ptr{Cvoid})
    error("TaoTermShellSetHessian: no generated method for these argument types")
end

@for_petsc function TaoTermShellSetHessian(petsclib::$UnionPetscLib, term::TaoTerm, hessian::Ptr{Cvoid} )

    @chk ccall(
               (:TaoTermShellSetHessian, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{Cvoid}),
               term, hessian,
              )


	return nothing
end 

"""
	TaoTermShellSetIsComputeHessianFDPossible(petsclib::PetscLibType, term::TaoTerm, ispossible::PetscBool3) 
Set whether this term can compute Hessian with finite differences for a `TAOTERMSHELL`

Logically collective

Input Parameters:
- `term`       - a `TaoTerm` of type `TAOTERMSHELL`
- `ispossible` - whether Hessian computation with finite differences is possible

Level: intermediate

See also: [](sec_tao_term), `TaoTerm`, `TAOTERMSHELL`, `TaoTermShellGetContext()`, `TaoTermShellSetContextDestroy()`,
`TaoTermShellSetObjective()`,
`TaoTermShellSetGradient()`,
`TaoTermShellSetObjectiveAndGradient()`,
`TaoTermShellSetHessian()`,
`TaoTermIsComputeHessianFDPossible()`

# External Links
$(_doc_external("TaoTerm/TaoTermShellSetIsComputeHessianFDPossible"))
"""
function TaoTermShellSetIsComputeHessianFDPossible(petsclib::PetscLibType, term::TaoTerm, ispossible::PetscBool3)
    error("TaoTermShellSetIsComputeHessianFDPossible: no generated method for these argument types")
end

@for_petsc function TaoTermShellSetIsComputeHessianFDPossible(petsclib::$UnionPetscLib, term::TaoTerm, ispossible::PetscBool3 )

    @chk ccall(
               (:TaoTermShellSetIsComputeHessianFDPossible, $petsc_library),
               PetscErrorCode,
               (TaoTerm, PetscBool3),
               term, ispossible,
              )


	return nothing
end 

"""
	TaoTermShellSetObjective(petsclib::PetscLibType, term::TaoTerm, objective::Ptr{Cvoid}) 
Set the objective function of a `TAOTERMSHELL`

Logically collective

Input Parameters:
- `term`      - a `TaoTerm` of type `TAOTERMSHELL`
- `objective` - a `TaoTermObjectiveFn` function pointer

Level: intermediate

See also: [](sec_tao_term), `TaoTerm`, `TAOTERMSHELL`, `TaoTermShellGetContext()`, `TaoTermShellSetContextDestroy()`,
`TaoTermShellSetGradient()`,
`TaoTermShellSetObjectiveAndGradient()`,
`TaoTermShellSetHessian()`,
`TaoTermShellSetView()`,
`TaoTermObjectiveFn`

# External Links
$(_doc_external("TaoTerm/TaoTermShellSetObjective"))
"""
function TaoTermShellSetObjective(petsclib::PetscLibType, term::TaoTerm, objective::Ptr{Cvoid})
    error("TaoTermShellSetObjective: no generated method for these argument types")
end

@for_petsc function TaoTermShellSetObjective(petsclib::$UnionPetscLib, term::TaoTerm, objective::Ptr{Cvoid} )

    @chk ccall(
               (:TaoTermShellSetObjective, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{Cvoid}),
               term, objective,
              )


	return nothing
end 

"""
	TaoTermShellSetObjectiveAndGradient(petsclib::PetscLibType, term::TaoTerm, objandgrad::Ptr{Cvoid}) 
Set the objective and gradient function of a `TAOTERMSHELL`

Logically collective

Input Parameters:
- `term`       - a `TaoTerm` of type `TAOTERMSHELL`
- `objandgrad` - a `TaoTermObjectiveAndGradientFn` function pointer

Level: intermediate

See also: [](sec_tao_term), `TaoTerm`, `TAOTERMSHELL`, `TaoTermShellGetContext()`, `TaoTermShellSetContextDestroy()`,
`TaoTermShellSetObjective()`,
`TaoTermShellSetGradient()`,
`TaoTermShellSetHessian()`,
`TaoTermShellSetView()`,
`TaoTermObjectiveAndGradientFn`

# External Links
$(_doc_external("TaoTerm/TaoTermShellSetObjectiveAndGradient"))
"""
function TaoTermShellSetObjectiveAndGradient(petsclib::PetscLibType, term::TaoTerm, objandgrad::Ptr{Cvoid})
    error("TaoTermShellSetObjectiveAndGradient: no generated method for these argument types")
end

@for_petsc function TaoTermShellSetObjectiveAndGradient(petsclib::$UnionPetscLib, term::TaoTerm, objandgrad::Ptr{Cvoid} )

    @chk ccall(
               (:TaoTermShellSetObjectiveAndGradient, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{Cvoid}),
               term, objandgrad,
              )


	return nothing
end 

"""
	TaoTermShellSetView(petsclib::PetscLibType, term::TaoTerm, view::external) 
Set the view function of a `TAOTERMSHELL`

Logically collective

Input Parameters:
- `term` - a `TaoTerm` of type `TAOTERMSHELL`
- `view` - a function with the same signature as `TaoTermView()`

Calling sequence of `view`:
- `term`   - the `TaoTerm`
- `viewer` - a `PetscViewer`

Level: intermediate

See also: [](sec_tao_term), `TaoTerm`, `TAOTERMSHELL`, `TaoTermShellGetContext()`, `TaoTermShellSetContextDestroy()`,
`TaoTermShellSetObjective()`,
`TaoTermShellSetGradient()`,
`TaoTermShellSetObjectiveAndGradient()`,
`TaoTermShellSetHessian()`

# External Links
$(_doc_external("TaoTerm/TaoTermShellSetView"))
"""
function TaoTermShellSetView(petsclib::PetscLibType, term::TaoTerm, view::external)
    error("TaoTermShellSetView: no generated method for these argument types")
end

@for_petsc function TaoTermShellSetView(petsclib::$UnionPetscLib, term::TaoTerm, view::external )

    @chk ccall(
               (:TaoTermShellSetView, $petsc_library),
               PetscErrorCode,
               (TaoTerm, external),
               term, view,
              )


	return nothing
end 

"""
	index::PetscInt = TaoTermSumAddTerm(petsclib::PetscLibType, sumterm::TaoTerm, prefix::String, scale::PetscReal, term::TaoTerm, map::AbstractPetscMat) 
Append a term to the terms being summed

Collective

Input Parameters:
- `sumterm` - a `TaoTerm` of type `TAOTERMSUM`
- `prefix`  - (optional) the prefix used for configuring the term (if `NULL`, the index of the term will be used as a prefix, e.g. `term_0_`, `term_1_`, etc.)
- `scale`   - the coefficient scaling the term in the sum
- `term`    - the `TaoTerm` to add
- `map`     - (optional) a map from the `TAOTERMSUM` solution space to the `term` solution space; if `NULL` the map is assumed to be the identity

Output Parameter:
- `index` - (optional) the index of the newly added term

Level: developer

See also: [](sec_tao_term), `TaoTerm`, `TAOTERMSUM`

# External Links
$(_doc_external("TaoTerm/TaoTermSumAddTerm"))
"""
function TaoTermSumAddTerm(petsclib::PetscLibType, sumterm::TaoTerm, prefix::String, scale::Real, term::TaoTerm, map::AbstractPetscMat)
    error("TaoTermSumAddTerm: no generated method for these argument types")
end

@for_petsc function TaoTermSumAddTerm(petsclib::$UnionPetscLib, sumterm::TaoTerm, prefix::String, scale::$PetscReal, term::TaoTerm, map::AbstractPetscMat )
	index_ = Ref{$PetscInt}()

    @chk ccall(
               (:TaoTermSumAddTerm, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{Cchar}, $PetscReal, TaoTerm, CMat, Ptr{$PetscInt}),
               sumterm, prefix, scale, term, map, index_,
              )

	index = index_[]

	return index
end 

"""
	values::Ptr{PetscReal} = TaoTermSumGetLastTermObjectives(petsclib::PetscLibType, term::TaoTerm) 
Get the contributions from each term to the
last evaluation of `TaoTermComputeObjective()` or `TaoTermComputeObjectiveAndGradient()`

Not collective

Input Parameter:
- `term` - a `TaoTerm` of type `TAOTERMSUM`

Output Parameter:
- `values` - an array of the contributions to the last computed objective value

Level: developer

See also: [](sec_tao_term),
`TaoTerm`,
`TAOTERMSUM`

# External Links
$(_doc_external("TaoTerm/TaoTermSumGetLastTermObjectives"))
"""
function TaoTermSumGetLastTermObjectives(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermSumGetLastTermObjectives: no generated method for these argument types")
end

@for_petsc function TaoTermSumGetLastTermObjectives(petsclib::$UnionPetscLib, term::TaoTerm )
	values_ = Ref{Ptr{$PetscReal}}()

    @chk ccall(
               (:TaoTermSumGetLastTermObjectives, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{Ptr{$PetscReal}}),
               term, values_,
              )

	values = values_[]

	return values
end 

"""
	n_terms::PetscInt = TaoTermSumGetNumberTerms(petsclib::PetscLibType, term::TaoTerm) 
Get the number of terms in the sum

Not collective

Input Parameter:
- `term` - a `TaoTerm` of type `TAOTERMSUM`

Output Parameter:
- `n_terms` - the number of terms that will be in the sum

Level: developer

See also: [](sec_tao_term),
`TaoTerm`,
`TAOTERMSUM`,
`TaoTermSumSetNumberTerms()`

# External Links
$(_doc_external("TaoTerm/TaoTermSumGetNumberTerms"))
"""
function TaoTermSumGetNumberTerms(petsclib::PetscLibType, term::TaoTerm)
    error("TaoTermSumGetNumberTerms: no generated method for these argument types")
end

@for_petsc function TaoTermSumGetNumberTerms(petsclib::$UnionPetscLib, term::TaoTerm )
	n_terms_ = Ref{$PetscInt}()

    @chk ccall(
               (:TaoTermSumGetNumberTerms, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{$PetscInt}),
               term, n_terms_,
              )

	n_terms = n_terms_[]

	return n_terms
end 

"""
	prefix::String,scale::PetscReal,term::TaoTerm,map::PetscMat = TaoTermSumGetTerm(petsclib::PetscLibType, sumterm::TaoTerm, index::PetscInt) 
Get the data for a term in a `TAOTERMSUM`

Not collective

Input Parameters:
- `sumterm` - a `TaoTerm` of type `TAOTERMSUM`
- `index`   - a number 0 \\leq i < n, where n is the number of terms in `TaoTermSumGetNumberTerms()`

Output Parameters:
- `prefix` - (optional) the prefix used for configuring the term
- `scale`  - (optional) the coefficient scaling the term in the sum
- `term`   - the `TaoTerm` at given index of `TAOTERMSUM`
- `map`    - (optional) a map from the `TAOTERMSUM` solution space to the `term` solution space; if `NULL` the map is assumed to be the identity

Level: developer

See also: [](sec_tao_term),
`TaoTerm`,
`TAOTERMSUM`,
`TaoTermSumSetTerm()`,
`TaoTermSumAddTerm()`

# External Links
$(_doc_external("TaoTerm/TaoTermSumGetTerm"))
"""
function TaoTermSumGetTerm(petsclib::PetscLibType, sumterm::TaoTerm, index::Integer)
    error("TaoTermSumGetTerm: no generated method for these argument types")
end

@for_petsc function TaoTermSumGetTerm(petsclib::$UnionPetscLib, sumterm::TaoTerm, index::$PetscInt )
	prefix_ = Ref{Ptr{Cchar}}()
	scale_ = Ref{$PetscReal}()
	term_ = Ref{TaoTerm}()
	map_ = Ref{CMat}()

    @chk ccall(
               (:TaoTermSumGetTerm, $petsc_library),
               PetscErrorCode,
               (TaoTerm, $PetscInt, Ptr{Ptr{Cchar}}, Ptr{$PetscReal}, Ptr{TaoTerm}, Ptr{CMat}),
               sumterm, index, prefix_, scale_, term_, map_,
              )

	prefix = unsafe_string(prefix_[])
	scale = scale_[]
	term = term_[]
	map = PetscMat(map_[], petsclib)

	return prefix,scale,term,map
end 

"""
	unmapped_H::PetscMat,unmapped_Hpre::PetscMat,mapped_H::PetscMat,mapped_Hpre::PetscMat = TaoTermSumGetTermHessianMatrices(petsclib::PetscLibType, term::TaoTerm, index::PetscInt) 
Get Hessian matrices set with `TaoTermSumSetTermHessianMatrices()`.

Not collective

Input Parameters:
- `term`  - a `TaoTerm` of type `TAOTERMSUM`
- `index` - the index for the term from `TaoTermSumSetTerm()` or `TaoTermSumAddTerm()`

Output Parameters:
- `unmapped_H`    - (optional) unmapped Hessian matrix
- `unmapped_Hpre` - (optional) unmapped matrix for constructing the preconditioner for `unmapped_H`
- `mapped_H`      - (optional) Hessian matrix
- `mapped_Hpre`   - (optional) matrix for constructing the preconditioner for `mapped_H`

Level: developer

See also: [](sec_tao_term),
`TaoTerm`,
`TAOTERMSUM`,
`TaoTermComputeHessian()`,
`TaoTermSumSetTermHessianMatrices()`

# External Links
$(_doc_external("TaoTerm/TaoTermSumGetTermHessianMatrices"))
"""
function TaoTermSumGetTermHessianMatrices(petsclib::PetscLibType, term::TaoTerm, index::Integer)
    error("TaoTermSumGetTermHessianMatrices: no generated method for these argument types")
end

@for_petsc function TaoTermSumGetTermHessianMatrices(petsclib::$UnionPetscLib, term::TaoTerm, index::$PetscInt )
	unmapped_H_ = Ref{CMat}()
	unmapped_Hpre_ = Ref{CMat}()
	mapped_H_ = Ref{CMat}()
	mapped_Hpre_ = Ref{CMat}()

    @chk ccall(
               (:TaoTermSumGetTermHessianMatrices, $petsc_library),
               PetscErrorCode,
               (TaoTerm, $PetscInt, Ptr{CMat}, Ptr{CMat}, Ptr{CMat}, Ptr{CMat}),
               term, index, unmapped_H_, unmapped_Hpre_, mapped_H_, mapped_Hpre_,
              )

	unmapped_H = PetscMat(unmapped_H_[], petsclib)
	unmapped_Hpre = PetscMat(unmapped_Hpre_[], petsclib)
	mapped_H = PetscMat(mapped_H_[], petsclib)
	mapped_Hpre = PetscMat(mapped_Hpre_[], petsclib)

	return unmapped_H,unmapped_Hpre,mapped_H,mapped_Hpre
end 

"""
	mask::TaoTermMask = TaoTermSumGetTermMask(petsclib::PetscLibType, term::TaoTerm, index::PetscInt) 
Get the `TaoTermMask` of a term in the sum

Not collective

Input Parameters:
- `term`  - a `TaoTerm` of type `TAOTERMSUM`
- `index` - the index for the term from `TaoTermSumSetTerm()` or `TaoTermSumAddTerm()`

Output Parameter:
- `mask` - a bitmask of `TaoTermMask` evaluation methods to mask (e.g. just `TAOTERM_MASK_OBJECTIVE` or a bitwise-or like `TAOTERM_MASK_OBJECTIVE | TAOTERM_MASK_GRADIENT`)

Level: developer

See also: [](sec_tao_term),
`TaoTerm`,
`TAOTERMSUM`,
`TaoTermSumSetTermMask()`

# External Links
$(_doc_external("TaoTerm/TaoTermSumGetTermMask"))
"""
function TaoTermSumGetTermMask(petsclib::PetscLibType, term::TaoTerm, index::Integer)
    error("TaoTermSumGetTermMask: no generated method for these argument types")
end

@for_petsc function TaoTermSumGetTermMask(petsclib::$UnionPetscLib, term::TaoTerm, index::$PetscInt )
	mask_ = Ref{TaoTermMask}()

    @chk ccall(
               (:TaoTermSumGetTermMask, $petsc_library),
               PetscErrorCode,
               (TaoTerm, $PetscInt, Ptr{TaoTermMask}),
               term, index, mask_,
              )

	mask = mask_[]

	return mask
end 

"""
	params::PetscVec = TaoTermSumParametersPack(petsclib::PetscLibType, term::TaoTerm, p_arr::Vector{<:AbstractPetscVec}) 
Concatenate the parameters for terms into a `VECNEST` parameter vector for a `TAOTERMSUM`

Collective

Input Parameters:
- `term`  - a `TaoTerm` of type `TAOTERMSUM`
- `p_arr` - an array of parameters `Vec`s, one for each term in the sum.  An entry can be `NULL` for a term that doesn't take parameters.

Output Parameter:
- `params` - a `Vec` of type `VECNEST` that concatenates all of the parameters

Level: developer

See also: [](sec_tao_term),
`TaoTerm`,
`TAOTERMSUM`,
`TaoTermSumParametersUnpack()`,
`VECNEST`,
`VecNestGetTaoTermSumParameters()`,
`VecCreateNest()`

# External Links
$(_doc_external("TaoTerm/TaoTermSumParametersPack"))
"""
function TaoTermSumParametersPack(petsclib::PetscLibType, term::TaoTerm, p_arr::Vector{<:AbstractPetscVec})
    error("TaoTermSumParametersPack: no generated method for these argument types")
end

@for_petsc function TaoTermSumParametersPack(petsclib::$UnionPetscLib, term::TaoTerm, p_arr::Vector{<:AbstractPetscVec} )
	params_ = Ref{CVec}()

    @chk ccall(
               (:TaoTermSumParametersPack, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{CVec}, Ptr{CVec}),
               term, p_arr, params_,
              )

	params = PetscVec(params_[], petsclib)

	return params
end 

"""
	TaoTermSumParametersUnpack(petsclib::PetscLibType, term::TaoTerm, params::AbstractPetscVec, p_arr::Vector{<:AbstractPetscVec}) 
Unpack the concatenated parameters created by `TaoTermSumParametersPack()` and destroy the `VECNEST`

Collective

Input Parameters:
- `term`   - a `TaoTerm` of type `TAOTERMSUM`
- `params` - a `Vec` created by `TaoTermSumParametersPack()`

Output Parameter:
- `p_arr` - an array of parameters `Vec`s, one for each term in the sum.  An entry will be `NULL` if `NULL` was passed in the same position of `TaoTermSumParametersPack()`

Level: intermediate

See also: [](sec_tao_term),
`TaoTerm`,
`TAOTERMSUM`,
`TaoTermSumParametersPack()`,
`VecNestGetTaoTermSumParameters()`

# External Links
$(_doc_external("TaoTerm/TaoTermSumParametersUnpack"))
"""
function TaoTermSumParametersUnpack(petsclib::PetscLibType, term::TaoTerm, params::AbstractPetscVec, p_arr::Vector{<:AbstractPetscVec})
    error("TaoTermSumParametersUnpack: no generated method for these argument types")
end

@for_petsc function TaoTermSumParametersUnpack(petsclib::$UnionPetscLib, term::TaoTerm, params::AbstractPetscVec, p_arr::Vector{<:AbstractPetscVec} )
	params_ = Ref(params.ptr)

    @chk ccall(
               (:TaoTermSumParametersUnpack, $petsc_library),
               PetscErrorCode,
               (TaoTerm, Ptr{CVec}, Ptr{CVec}),
               term, params_, p_arr,
              )

	params.ptr = params_[]

	return nothing
end 

"""
	TaoTermSumSetNumberTerms(petsclib::PetscLibType, term::TaoTerm, n_terms::PetscInt) 
Set the number of terms in the sum

Collective

Input Parameters:
- `term`    - a `TaoTerm` of type `TAOTERMSUM`
- `n_terms` - the number of terms that will be in the sum

Level: developer

See also: [](sec_tao_term),
`TaoTerm`,
`TAOTERMSUM`,
`TaoTermSumGetNumberTerms()`

# External Links
$(_doc_external("TaoTerm/TaoTermSumSetNumberTerms"))
"""
function TaoTermSumSetNumberTerms(petsclib::PetscLibType, term::TaoTerm, n_terms::Integer)
    error("TaoTermSumSetNumberTerms: no generated method for these argument types")
end

@for_petsc function TaoTermSumSetNumberTerms(petsclib::$UnionPetscLib, term::TaoTerm, n_terms::$PetscInt )

    @chk ccall(
               (:TaoTermSumSetNumberTerms, $petsc_library),
               PetscErrorCode,
               (TaoTerm, $PetscInt),
               term, n_terms,
              )


	return nothing
end 

"""
	TaoTermSumSetTerm(petsclib::PetscLibType, sumterm::TaoTerm, index::PetscInt, prefix::String, scale::PetscReal, term::TaoTerm, map::AbstractPetscMat) 
Set a term in a sum of terms

Collective

Input Parameters:
- `sumterm` - a `TaoTerm` of type `TAOTERMSUM`
- `index`   - a number 0 \\leq i < n, where n is the number of terms in `TaoTermSumSetNumberTerms()`
- `prefix`  - (optional) the prefix used for configuring the term (if `NULL`, `term_x_` will be the prefix, e.g. "term_0_", "term_1_", etc.)
- `scale`   - the coefficient scaling the term in the sum
- `term`    - the `TaoTerm` to be set in `TAOTERMSUM`
- `map`     - (optional) a map from the `TAOTERMSUM` solution space to the `term` solution space; if `NULL` the map is assumed to be the identity

Level: developer

See also: [](sec_tao_term),
`TaoTerm`,
`TAOTERMSUM`,
`TaoTermSumGetTerm()`,
`TaoTermSumAddTerm()`

# External Links
$(_doc_external("TaoTerm/TaoTermSumSetTerm"))
"""
function TaoTermSumSetTerm(petsclib::PetscLibType, sumterm::TaoTerm, index::Integer, prefix::String, scale::Real, term::TaoTerm, map::AbstractPetscMat)
    error("TaoTermSumSetTerm: no generated method for these argument types")
end

@for_petsc function TaoTermSumSetTerm(petsclib::$UnionPetscLib, sumterm::TaoTerm, index::$PetscInt, prefix::String, scale::$PetscReal, term::TaoTerm, map::AbstractPetscMat )

    @chk ccall(
               (:TaoTermSumSetTerm, $petsc_library),
               PetscErrorCode,
               (TaoTerm, $PetscInt, Ptr{Cchar}, $PetscReal, TaoTerm, CMat),
               sumterm, index, prefix, scale, term, map,
              )


	return nothing
end 

"""
	TaoTermSumSetTermHessianMatrices(petsclib::PetscLibType, term::TaoTerm, index::PetscInt, unmapped_H::AbstractPetscMat, unmapped_Hpre::AbstractPetscMat, mapped_H::AbstractPetscMat, mapped_Hpre::AbstractPetscMat) 
Set Hessian matrices that can be used internally by a `TAOTERMSUM`

Logically collective

Input Parameters:
- `term`          - a `TaoTerm` of type `TAOTERMSUM`
- `index`         - the index for the term from `TaoTermSumSetTerm()` or `TaoTermSumAddTerm()`
- `unmapped_H`    - (optional) unmapped Hessian matrix
- `unmapped_Hpre` - (optional) unmapped matrix for constructing the preconditioner of `unmapped_H`
- `mapped_H`      - (optional) Hessian matrix
- `mapped_Hpre`   - (optional) matrix for constructing the preconditioner of `mapped_H`

Level: developer

See also: [](sec_tao_term),
`TaoTerm`,
`TAOTERMSUM`,
`TaoTermComputeHessian()`,
`TaoTermSumGetTermHessianMatrices()`

# External Links
$(_doc_external("TaoTerm/TaoTermSumSetTermHessianMatrices"))
"""
function TaoTermSumSetTermHessianMatrices(petsclib::PetscLibType, term::TaoTerm, index::Integer, unmapped_H::AbstractPetscMat, unmapped_Hpre::AbstractPetscMat, mapped_H::AbstractPetscMat, mapped_Hpre::AbstractPetscMat)
    error("TaoTermSumSetTermHessianMatrices: no generated method for these argument types")
end

@for_petsc function TaoTermSumSetTermHessianMatrices(petsclib::$UnionPetscLib, term::TaoTerm, index::$PetscInt, unmapped_H::AbstractPetscMat, unmapped_Hpre::AbstractPetscMat, mapped_H::AbstractPetscMat, mapped_Hpre::AbstractPetscMat )

    @chk ccall(
               (:TaoTermSumSetTermHessianMatrices, $petsc_library),
               PetscErrorCode,
               (TaoTerm, $PetscInt, CMat, CMat, CMat, CMat),
               term, index, unmapped_H, unmapped_Hpre, mapped_H, mapped_Hpre,
              )


	return nothing
end 

"""
	TaoTermSumSetTermMask(petsclib::PetscLibType, term::TaoTerm, index::PetscInt, mask::TaoTermMask) 
Set a `TaoTermMask` on a term in the sum

Logically collective

Input Parameters:
- `term`  - a `TaoTerm` of type `TAOTERMSUM`
- `index` - the index for the term from `TaoTermSumSetTerm()` or `TaoTermSumAddTerm()`
- `mask`  - a bitmask of `TaoTermMask` evaluation methods to mask (e.g. just `TAOTERM_MASK_OBJECTIVE` or a bitwise-or like `TAOTERM_MASK_OBJECTIVE | TAOTERM_MASK_GRADIENT`)

Options Database Keys:
- `-tao_term_sum_<prefix_>mask` - a list containing any of `none`, `objective`, `gradient`, and `hessian` to indicate which evaluations to mask for a term with a given prefix (see `TaoTermSumSetTerm()`)

Level: developer

See also: [](sec_tao_term),
`TaoTerm`,
`TAOTERMSUM`,
`TaoTermSumGetTermMask()`

# External Links
$(_doc_external("TaoTerm/TaoTermSumSetTermMask"))
"""
function TaoTermSumSetTermMask(petsclib::PetscLibType, term::TaoTerm, index::Integer, mask::TaoTermMask)
    error("TaoTermSumSetTermMask: no generated method for these argument types")
end

@for_petsc function TaoTermSumSetTermMask(petsclib::$UnionPetscLib, term::TaoTerm, index::$PetscInt, mask::TaoTermMask )

    @chk ccall(
               (:TaoTermSumSetTermMask, $petsc_library),
               PetscErrorCode,
               (TaoTerm, $PetscInt, TaoTermMask),
               term, index, mask,
              )


	return nothing
end 

"""
	TaoTermView(petsclib::PetscLibType, term::TaoTerm, viewer::PetscViewer) 
View a description of a `TaoTerm`.

Collective

Input Parameters:
- `term`   - a `TaoTerm`
- `viewer` - a `PetscViewer`

Level: beginner

See also: [](sec_tao_term),
`TaoTerm`,
`TaoTermCreate()`,
`TaoTermSetType()`,
`TaoTermSetFromOptions()`,
`TaoTermSetUp()`,
`TaoTermDestroy()`,
`PetscViewer`

# External Links
$(_doc_external("TaoTerm/TaoTermView"))
"""
function TaoTermView(petsclib::PetscLibType, term::TaoTerm, viewer::PetscViewer)
    error("TaoTermView: no generated method for these argument types")
end

@for_petsc function TaoTermView(petsclib::$UnionPetscLib, term::TaoTerm, viewer::PetscViewer )

    @chk ccall(
               (:TaoTermView, $petsc_library),
               PetscErrorCode,
               (TaoTerm, PetscViewer),
               term, viewer,
              )


	return nothing
end 

