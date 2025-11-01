# autodefined type arguments for class ------
mutable struct _n_PetscConvEst end
const PetscConvEst = Ptr{_n_PetscConvEst}

# -------------------------------------------------------
"""
	PetscConvEstDestroy(petsclib::PetscLibType,ce::PetscConvEst) 
Destroys a PETSc convergence estimator `PetscConvEst` object

Collective

Input Parameter:
- `ce` - The `PetscConvEst` object

Level: beginner

-seealso: `PetscConvEst`, `PetscConvEstCreate()`, `PetscConvEstGetConvRate()`

# External Links
$(_doc_external("Snes/PetscConvEstDestroy"))
"""
function PetscConvEstDestroy(petsclib::PetscLibType, ce::PetscConvEst) end

@for_petsc function PetscConvEstDestroy(petsclib::$UnionPetscLib, ce::PetscConvEst )

    @chk ccall(
               (:PetscConvEstDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscConvEst},),
               ce,
              )


	return nothing
end 

"""
	PetscConvEstSetFromOptions(petsclib::PetscLibType,ce::PetscConvEst) 
Sets a convergence estimator `PetscConvEst` object based on values in the options database

Collective

Input Parameter:
- `ce` - The `PetscConvEst` object

Level: beginner

-seealso: `PetscConvEst`, `PetscConvEstCreate()`, `PetscConvEstGetConvRate()`

# External Links
$(_doc_external("Snes/PetscConvEstSetFromOptions"))
"""
function PetscConvEstSetFromOptions(petsclib::PetscLibType, ce::PetscConvEst) end

@for_petsc function PetscConvEstSetFromOptions(petsclib::$UnionPetscLib, ce::PetscConvEst )

    @chk ccall(
               (:PetscConvEstSetFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscConvEst,),
               ce,
              )


	return nothing
end 

"""
	PetscConvEstView(petsclib::PetscLibType,ce::PetscConvEst, viewer::PetscViewer) 
Views a `PetscConvEst` object

Collective

Input Parameters:
- `ce`     - The `PetscConvEst` object
- `viewer` - The `PetscViewer`

Level: beginner

-seealso: `PetscConvEst`, `PetscViewer`, `PetscConvEstCreate()`, `PetscConvEstGetConvRate()`

# External Links
$(_doc_external("Snes/PetscConvEstView"))
"""
function PetscConvEstView(petsclib::PetscLibType, ce::PetscConvEst, viewer::PetscViewer) end

@for_petsc function PetscConvEstView(petsclib::$UnionPetscLib, ce::PetscConvEst, viewer::PetscViewer )

    @chk ccall(
               (:PetscConvEstView, $petsc_library),
               PetscErrorCode,
               (PetscConvEst, PetscViewer),
               ce, viewer,
              )


	return nothing
end 

"""
	PetscConvEstGetSolver(petsclib::PetscLibType,ce::PetscConvEst, solver::PetscObject) 
Gets the solver used to produce discrete solutions

Not Collective

Input Parameter:
- `ce` - The `PetscConvEst` object

Output Parameter:
- `solver` - The solver

Level: intermediate

-seealso: `PetscConvEst`, `PetscConvEstSetSolver()`, `PetscConvEstCreate()`, `PetscConvEstGetConvRate()`

# External Links
$(_doc_external("Snes/PetscConvEstGetSolver"))
"""
function PetscConvEstGetSolver(petsclib::PetscLibType, ce::PetscConvEst, solver::PetscObject) end

@for_petsc function PetscConvEstGetSolver(petsclib::$UnionPetscLib, ce::PetscConvEst, solver::PetscObject )

    @chk ccall(
               (:PetscConvEstGetSolver, $petsc_library),
               PetscErrorCode,
               (PetscConvEst, Ptr{PetscObject}),
               ce, solver,
              )


	return nothing
end 

"""
	PetscConvEstSetSolver(petsclib::PetscLibType,ce::PetscConvEst, solver::PetscObject) 
Sets the solver used to produce discrete solutions

Not Collective

Input Parameters:
- `ce`     - The `PetscConvEst` object
- `solver` - The solver, must be a `KSP`, `SNES`, or `TS` object with an attached `DM`/`DS`, that can compute an exact solution

Level: intermediate

-seealso: `PetscConvEst`, `PetscConvEstGetSNES()`, `PetscConvEstCreate()`, `PetscConvEstGetConvRate()`

# External Links
$(_doc_external("Snes/PetscConvEstSetSolver"))
"""
function PetscConvEstSetSolver(petsclib::PetscLibType, ce::PetscConvEst, solver::PetscObject) end

@for_petsc function PetscConvEstSetSolver(petsclib::$UnionPetscLib, ce::PetscConvEst, solver::PetscObject )

    @chk ccall(
               (:PetscConvEstSetSolver, $petsc_library),
               PetscErrorCode,
               (PetscConvEst, PetscObject),
               ce, solver,
              )


	return nothing
end 

"""
	PetscConvEstSetUp(petsclib::PetscLibType,ce::PetscConvEst) 
After the solver is specified, create data structures needed for estimating convergence

Collective

Input Parameter:
- `ce` - The `PetscConvEst` object

Level: beginner

-seealso: `PetscConvEst`, `PetscConvEstCreate()`, `PetscConvEstGetConvRate()`

# External Links
$(_doc_external("Snes/PetscConvEstSetUp"))
"""
function PetscConvEstSetUp(petsclib::PetscLibType, ce::PetscConvEst) end

@for_petsc function PetscConvEstSetUp(petsclib::$UnionPetscLib, ce::PetscConvEst )

    @chk ccall(
               (:PetscConvEstSetUp, $petsc_library),
               PetscErrorCode,
               (PetscConvEst,),
               ce,
              )


	return nothing
end 

"""
	PetscConvEstComputeInitialGuess(petsclib::PetscLibType,ce::PetscConvEst, r::PetscInt, dm::PetscDM, u::PetscVec) 

# External Links
$(_doc_external("Snes/PetscConvEstComputeInitialGuess"))
"""
function PetscConvEstComputeInitialGuess(petsclib::PetscLibType, ce::PetscConvEst, r::PetscInt, dm::PetscDM, u::PetscVec) end

@for_petsc function PetscConvEstComputeInitialGuess(petsclib::$UnionPetscLib, ce::PetscConvEst, r::$PetscInt, dm::PetscDM, u::PetscVec )

    @chk ccall(
               (:PetscConvEstComputeInitialGuess, $petsc_library),
               PetscErrorCode,
               (PetscConvEst, $PetscInt, CDM, CVec),
               ce, r, dm, u,
              )


	return nothing
end 

"""
	PetscConvEstComputeError(petsclib::PetscLibType,ce::PetscConvEst, r::PetscInt, dm::PetscDM, u::PetscVec, errors::Vector{PetscReal}) 

# External Links
$(_doc_external("Snes/PetscConvEstComputeError"))
"""
function PetscConvEstComputeError(petsclib::PetscLibType, ce::PetscConvEst, r::PetscInt, dm::PetscDM, u::PetscVec, errors::Vector{PetscReal}) end

@for_petsc function PetscConvEstComputeError(petsclib::$UnionPetscLib, ce::PetscConvEst, r::$PetscInt, dm::PetscDM, u::PetscVec, errors::Vector{$PetscReal} )

    @chk ccall(
               (:PetscConvEstComputeError, $petsc_library),
               PetscErrorCode,
               (PetscConvEst, $PetscInt, CDM, CVec, Ptr{$PetscReal}),
               ce, r, dm, u, errors,
              )


	return nothing
end 

"""
	PetscConvEstMonitorDefault(petsclib::PetscLibType,ce::PetscConvEst, r::PetscInt) 
Monitors the convergence estimation loop

Collective

Input Parameters:
- `ce` - The `PetscConvEst` object
- `r`  - The refinement level

Options Database Key:
- `-convest_monitor` - Activate the monitor

Level: intermediate

-seealso: `PetscConvEst`, `PetscConvEstCreate()`, `PetscConvEstGetConvRate()`, `SNESSolve()`, `TSSolve()`

# External Links
$(_doc_external("Snes/PetscConvEstMonitorDefault"))
"""
function PetscConvEstMonitorDefault(petsclib::PetscLibType, ce::PetscConvEst, r::PetscInt) end

@for_petsc function PetscConvEstMonitorDefault(petsclib::$UnionPetscLib, ce::PetscConvEst, r::$PetscInt )

    @chk ccall(
               (:PetscConvEstMonitorDefault, $petsc_library),
               PetscErrorCode,
               (PetscConvEst, $PetscInt),
               ce, r,
              )


	return nothing
end 

"""
	alpha::Vector{PetscReal} = PetscConvEstGetConvRate(petsclib::PetscLibType,ce::PetscConvEst) 
Returns an estimate of the convergence rate for the discretization

Not Collective

Input Parameter:
- `ce` - The `PetscConvEst` object

Output Parameter:
- `alpha` - The convergence rate for each field

Options Database Keys:
- `-snes_convergence_estimate` - Execute convergence estimation inside `SNESSolve()` and print out the rate
- `-ts_convergence_estimate`   - Execute convergence estimation inside `TSSolve()` and print out the rate

Level: intermediate

-seealso: `PetscConvEstSetSolver()`, `PetscConvEstCreate()`, `SNESSolve()`, `TSSolve()`

# External Links
$(_doc_external("Snes/PetscConvEstGetConvRate"))
"""
function PetscConvEstGetConvRate(petsclib::PetscLibType, ce::PetscConvEst) end

@for_petsc function PetscConvEstGetConvRate(petsclib::$UnionPetscLib, ce::PetscConvEst )
	alpha = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscConvEstGetConvRate, $petsc_library),
               PetscErrorCode,
               (PetscConvEst, Ptr{$PetscReal}),
               ce, alpha,
              )


	return alpha
end 

"""
	PetscConvEstRateView(petsclib::PetscLibType,ce::PetscConvEst, alpha::Vector{PetscReal}, viewer::PetscViewer) 
Displays the convergence rate obtained from `PetscConvEstGetConvRate()` using a `PetscViewer`

Collective

Input Parameters:
- `ce`     - iterative context obtained from `SNESCreate()`
- `alpha`  - the convergence rate for each field
- `viewer` - the viewer to display the reason

Options Database Key:
- `-snes_convergence_estimate` - print the convergence rate

Level: developer

-seealso: `PetscConvEst`, `PetscConvEstGetConvRate()`

# External Links
$(_doc_external("Snes/PetscConvEstRateView"))
"""
function PetscConvEstRateView(petsclib::PetscLibType, ce::PetscConvEst, alpha::Vector{PetscReal}, viewer::PetscViewer) end

@for_petsc function PetscConvEstRateView(petsclib::$UnionPetscLib, ce::PetscConvEst, alpha::Vector{$PetscReal}, viewer::PetscViewer )

    @chk ccall(
               (:PetscConvEstRateView, $petsc_library),
               PetscErrorCode,
               (PetscConvEst, Ptr{$PetscReal}, PetscViewer),
               ce, alpha, viewer,
              )


	return nothing
end 

"""
	ce::PetscConvEst = PetscConvEstCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Create a `PetscConvEst` object. This is used to study the convergence rate of approximations on grids to a continuum solution

Collective

Input Parameter:
- `comm` - The communicator for the `PetscConvEst` object

Output Parameter:
- `ce` - The `PetscConvEst` object

Level: beginner

-seealso: `PetscConvEst`, `PetscConvEstDestroy()`, `PetscConvEstGetConvRate()`, `DMAdaptorCreate()`, `DMAdaptor`

# External Links
$(_doc_external("Snes/PetscConvEstCreate"))
"""
function PetscConvEstCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function PetscConvEstCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	ce_ = Ref{PetscConvEst}()

    @chk ccall(
               (:PetscConvEstCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscConvEst}),
               comm, ce_,
              )

	ce = ce_[]

	return ce
end 

"""
	PetscConvEstUseTS(petsclib::PetscLibType,ce::PetscConvEst, checkTemporal::PetscBool) 

# External Links
$(_doc_external("Ts/PetscConvEstUseTS"))
"""
function PetscConvEstUseTS(petsclib::PetscLibType, ce::PetscConvEst, checkTemporal::PetscBool) end

@for_petsc function PetscConvEstUseTS(petsclib::$UnionPetscLib, ce::PetscConvEst, checkTemporal::PetscBool )

    @chk ccall(
               (:PetscConvEstUseTS, $petsc_library),
               PetscErrorCode,
               (PetscConvEst, PetscBool),
               ce, checkTemporal,
              )


	return nothing
end 

