"""
	PetscConvEstComputeError(petsclib::PetscLibType, ce::PetscConvEst, r::PetscInt, dm::AbstractPetscDM, u::AbstractPetscVec, errors::Vector{PetscReal}) 
Compute per-field discretization errors of `u` on refinement level `r`

Collective

Input Parameters:
- `ce` - the `PetscConvEst` object
- `r`  - the refinement level
- `dm` - the `DM` on which `u` is defined (may be `NULL`)
- `u`  - the computed solution

Output Parameter:
- `errors` - array of length `Nf` (number of fields in the DS) filled with the error in each field

Level: developer

See also: `PetscConvEst`, `PetscConvEstComputeInitialGuess()`, `PetscConvEstGetConvRate()`

# External Links
$(_doc_external("SNES/PetscConvEstComputeError"))
"""
function PetscConvEstComputeError(petsclib::PetscLibType, ce::PetscConvEst, r::Integer, dm::AbstractPetscDM, u::AbstractPetscVec, errors::AbstractVector{<:Number})
    error("PetscConvEstComputeError: no generated method for these argument types")
end

@for_petsc function PetscConvEstComputeError(petsclib::$UnionPetscLib, ce::PetscConvEst, r::$PetscInt, dm::AbstractPetscDM, u::AbstractPetscVec, errors::Vector{$PetscReal} )

    @chk ccall(
               (:PetscConvEstComputeError, $petsc_library),
               PetscErrorCode,
               (PetscConvEst, $PetscInt, CDM, CVec, Ptr{$PetscReal}),
               ce, r, dm, u, errors,
              )


	return nothing
end 

"""
	PetscConvEstComputeInitialGuess(petsclib::PetscLibType, ce::PetscConvEst, r::PetscInt, dm::AbstractPetscDM, u::AbstractPetscVec) 
Fill `u` with the initial guess to use on refinement level `r` of a convergence-estimation run

Collective

Input Parameters:
- `ce` - the `PetscConvEst` object
- `r`  - the refinement level
- `dm` - the `DM` on which `u` is defined (may be `NULL`)

Output Parameter:
- `u` - the initial-guess vector

Level: developer

See also: `PetscConvEst`, `PetscConvEstComputeError()`, `PetscConvEstGetConvRate()`

# External Links
$(_doc_external("SNES/PetscConvEstComputeInitialGuess"))
"""
function PetscConvEstComputeInitialGuess(petsclib::PetscLibType, ce::PetscConvEst, r::Integer, dm::AbstractPetscDM, u::AbstractPetscVec)
    error("PetscConvEstComputeInitialGuess: no generated method for these argument types")
end

@for_petsc function PetscConvEstComputeInitialGuess(petsclib::$UnionPetscLib, ce::PetscConvEst, r::$PetscInt, dm::AbstractPetscDM, u::AbstractPetscVec )

    @chk ccall(
               (:PetscConvEstComputeInitialGuess, $petsc_library),
               PetscErrorCode,
               (PetscConvEst, $PetscInt, CDM, CVec),
               ce, r, dm, u,
              )


	return nothing
end 

"""
	ce::PetscConvEst = PetscConvEstCreate(petsclib::PetscLibType, comm::MPI_Comm) 
Create a `PetscConvEst` object. This is used to study the convergence rate of approximations on grids to a continuum solution

Collective

Input Parameter:
- `comm` - The communicator for the `PetscConvEst` object

Output Parameter:
- `ce` - The `PetscConvEst` object

Level: beginner

See also: `PetscConvEst`, `PetscConvEstDestroy()`, `PetscConvEstGetConvRate()`, `DMAdaptorCreate()`, `DMAdaptor`

# External Links
$(_doc_external("SNES/PetscConvEstCreate"))
"""
function PetscConvEstCreate(petsclib::PetscLibType, comm::MPI_Comm)
    error("PetscConvEstCreate: no generated method for these argument types")
end

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
	PetscConvEstDestroy(petsclib::PetscLibType, ce::Union{PetscConvEst, Ref{PetscConvEst}}) 
Destroys a PETSc convergence estimator `PetscConvEst` object

Collective

Input Parameter:
- `ce` - The `PetscConvEst` object

Level: beginner

See also: `PetscConvEst`, `PetscConvEstCreate()`, `PetscConvEstGetConvRate()`

# External Links
$(_doc_external("SNES/PetscConvEstDestroy"))
"""
function PetscConvEstDestroy(petsclib::PetscLibType, ce::Union{PetscConvEst, Ref{PetscConvEst}})
    error("PetscConvEstDestroy: no generated method for these argument types")
end

@for_petsc function PetscConvEstDestroy(petsclib::$UnionPetscLib, ce::Union{PetscConvEst, Ref{PetscConvEst}} )
	ce_ = ce isa Base.RefValue ? ce : Ref{PetscConvEst}(ce)

    @chk ccall(
               (:PetscConvEstDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscConvEst},),
               ce_,
              )


	return nothing
end 

"""
	PetscConvEstGetConvRate(petsclib::PetscLibType, ce::PetscConvEst, alpha::Vector{PetscReal}) 
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

See also: `PetscConvEstSetSolver()`, `PetscConvEstCreate()`, `SNESSolve()`, `TSSolve()`

# External Links
$(_doc_external("SNES/PetscConvEstGetConvRate"))
"""
function PetscConvEstGetConvRate(petsclib::PetscLibType, ce::PetscConvEst, alpha::AbstractVector{<:Number})
    error("PetscConvEstGetConvRate: no generated method for these argument types")
end

@for_petsc function PetscConvEstGetConvRate(petsclib::$UnionPetscLib, ce::PetscConvEst, alpha::Vector{$PetscReal} )

    @chk ccall(
               (:PetscConvEstGetConvRate, $petsc_library),
               PetscErrorCode,
               (PetscConvEst, Ptr{$PetscReal}),
               ce, alpha,
              )


	return nothing
end 

"""
	solver::PetscObject = PetscConvEstGetSolver(petsclib::PetscLibType, ce::PetscConvEst) 
Gets the solver used to produce discrete solutions

Not Collective

Input Parameter:
- `ce` - The `PetscConvEst` object

Output Parameter:
- `solver` - The solver

Level: intermediate

See also: `PetscConvEst`, `PetscConvEstSetSolver()`, `PetscConvEstCreate()`, `PetscConvEstGetConvRate()`

# External Links
$(_doc_external("SNES/PetscConvEstGetSolver"))
"""
function PetscConvEstGetSolver(petsclib::PetscLibType, ce::PetscConvEst)
    error("PetscConvEstGetSolver: no generated method for these argument types")
end

@for_petsc function PetscConvEstGetSolver(petsclib::$UnionPetscLib, ce::PetscConvEst )
	solver_ = Ref{PetscObject}()

    @chk ccall(
               (:PetscConvEstGetSolver, $petsc_library),
               PetscErrorCode,
               (PetscConvEst, Ptr{PetscObject}),
               ce, solver_,
              )

	solver = solver_[]

	return solver
end 

"""
	PetscConvEstMonitorDefault(petsclib::PetscLibType, ce::PetscConvEst, r::PetscInt) 
Monitors the convergence estimation loop

Collective

Input Parameters:
- `ce` - The `PetscConvEst` object
- `r`  - The refinement level

Options Database Key:
- `-convest_monitor` - Activate the monitor

Level: intermediate

See also: `PetscConvEst`, `PetscConvEstCreate()`, `PetscConvEstGetConvRate()`, `SNESSolve()`, `TSSolve()`

# External Links
$(_doc_external("SNES/PetscConvEstMonitorDefault"))
"""
function PetscConvEstMonitorDefault(petsclib::PetscLibType, ce::PetscConvEst, r::Integer)
    error("PetscConvEstMonitorDefault: no generated method for these argument types")
end

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
	PetscConvEstRateView(petsclib::PetscLibType, ce::PetscConvEst, alpha::Vector{PetscReal}, viewer::PetscViewer) 
Displays the convergence rate obtained from `PetscConvEstGetConvRate()` using a `PetscViewer`

Collective

Input Parameters:
- `ce`     - iterative context obtained from `SNESCreate()`
- `alpha`  - the convergence rate for each field
- `viewer` - the viewer to display the reason

Options Database Key:
- `-snes_convergence_estimate` - print the convergence rate

Level: developer

See also: `PetscConvEst`, `PetscConvEstGetConvRate()`

# External Links
$(_doc_external("SNES/PetscConvEstRateView"))
"""
function PetscConvEstRateView(petsclib::PetscLibType, ce::PetscConvEst, alpha::AbstractVector{<:Number}, viewer::PetscViewer)
    error("PetscConvEstRateView: no generated method for these argument types")
end

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
	PetscConvEstSetFromOptions(petsclib::PetscLibType, ce::PetscConvEst) 
Sets a convergence estimator `PetscConvEst` object based on values in the options database

Collective

Input Parameter:
- `ce` - The `PetscConvEst` object

Level: beginner

See also: `PetscConvEst`, `PetscConvEstCreate()`, `PetscConvEstGetConvRate()`

# External Links
$(_doc_external("SNES/PetscConvEstSetFromOptions"))
"""
function PetscConvEstSetFromOptions(petsclib::PetscLibType, ce::PetscConvEst)
    error("PetscConvEstSetFromOptions: no generated method for these argument types")
end

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
	PetscConvEstSetSolver(petsclib::PetscLibType, ce::PetscConvEst, solver) 
Sets the solver used to produce discrete solutions

Not Collective

Input Parameters:
- `ce`     - The `PetscConvEst` object
- `solver` - The solver, must be a `KSP`, `SNES`, or `TS` object with an attached `DM`/`DS`, that can compute an exact solution

Level: intermediate

See also: `PetscConvEst`, `PetscConvEstGetSNES()`, `PetscConvEstCreate()`, `PetscConvEstGetConvRate()`

# External Links
$(_doc_external("SNES/PetscConvEstSetSolver"))
"""
function PetscConvEstSetSolver(petsclib::PetscLibType, ce::PetscConvEst, solver)
    error("PetscConvEstSetSolver: no generated method for these argument types")
end

@for_petsc function PetscConvEstSetSolver(petsclib::$UnionPetscLib, ce::PetscConvEst, solver )

    @chk ccall(
               (:PetscConvEstSetSolver, $petsc_library),
               PetscErrorCode,
               (PetscConvEst, PetscObject),
               ce, solver,
              )


	return nothing
end 

"""
	PetscConvEstSetUp(petsclib::PetscLibType, ce::PetscConvEst) 
After the solver is specified, create data structures needed for estimating convergence

Collective

Input Parameter:
- `ce` - The `PetscConvEst` object

Level: beginner

See also: `PetscConvEst`, `PetscConvEstCreate()`, `PetscConvEstGetConvRate()`

# External Links
$(_doc_external("SNES/PetscConvEstSetUp"))
"""
function PetscConvEstSetUp(petsclib::PetscLibType, ce::PetscConvEst)
    error("PetscConvEstSetUp: no generated method for these argument types")
end

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
	PetscConvEstUseTS(petsclib::PetscLibType, ce::PetscConvEst, checkTemporal::PetscBool) 
Configure a `PetscConvEst` object to use a `TS` solver for a convergence study

Not Collective

Input Parameters:
- `ce`            - the convergence estimator
- `checkTemporal` - `PETSC_TRUE` to run a temporal convergence study (refining the time step), `PETSC_FALSE` to run a spatial convergence study (refining the mesh)

Level: intermediate

See also: `PetscConvEst`, `TS`, `PetscConvEstGetConvRate()`

# External Links
$(_doc_external("TS/PetscConvEstUseTS"))
"""
function PetscConvEstUseTS(petsclib::PetscLibType, ce::PetscConvEst, checkTemporal::PetscBool)
    error("PetscConvEstUseTS: no generated method for these argument types")
end

@for_petsc function PetscConvEstUseTS(petsclib::$UnionPetscLib, ce::PetscConvEst, checkTemporal::PetscBool )

    @chk ccall(
               (:PetscConvEstUseTS, $petsc_library),
               PetscErrorCode,
               (PetscConvEst, PetscBool),
               ce, checkTemporal,
              )


	return nothing
end 

"""
	PetscConvEstView(petsclib::PetscLibType, ce::PetscConvEst, viewer::PetscViewer) 
Views a `PetscConvEst` object

Collective

Input Parameters:
- `ce`     - The `PetscConvEst` object
- `viewer` - The `PetscViewer`

Level: beginner

See also: `PetscConvEst`, `PetscViewer`, `PetscConvEstCreate()`, `PetscConvEstGetConvRate()`

# External Links
$(_doc_external("SNES/PetscConvEstView"))
"""
function PetscConvEstView(petsclib::PetscLibType, ce::PetscConvEst, viewer::PetscViewer)
    error("PetscConvEstView: no generated method for these argument types")
end

@for_petsc function PetscConvEstView(petsclib::$UnionPetscLib, ce::PetscConvEst, viewer::PetscViewer )

    @chk ccall(
               (:PetscConvEstView, $petsc_library),
               PetscErrorCode,
               (PetscConvEst, PetscViewer),
               ce, viewer,
              )


	return nothing
end 

