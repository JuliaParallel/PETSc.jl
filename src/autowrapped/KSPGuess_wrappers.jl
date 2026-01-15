"""
	KSPGuessRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Registers a method for initial guess computation in Krylov subspace solver package.

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - name of a new user-defined solver
- `function` - routine to create method context

-seealso: [](ch_ksp), `KSPGuess`, `KSPGuessRegisterAll()`

# External Links
$(_doc_external("Ksp/KSPGuessRegister"))
"""
function KSPGuessRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function KSPGuessRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:KSPGuessRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	KSPGuessSetFromOptions(petsclib::PetscLibType,guess::KSPGuess) 
Sets the options for a `KSPGuess` from the options database

Collective

Input Parameter:
- `guess` - `KSPGuess` object

Options Database Keys:
- `-ksp_guess_type <method>`       - Turns on generation of initial guesses and sets the method; use -help for a list of available methods
- `-ksp_guess_view <viewer>`       - view the `KSPGuess` object
- `-ksp_guess_fischer_model <a,b>` - set details for the Fischer models
- `-ksp_guess_fischer_monitor`     - monitor the Fischer models
- `-ksp_guess_fischer_tol <tol>`   - set the tolerance for the Fischer models
- `-ksp_guess_pod_size <size>`     - Number of snapshots
- `-ksp_guess_pod_monitor true`    - monitor the pod initial guess processing
- `-ksp_guess_pod_tol <tol>`       - Tolerance to retain eigenvectors
- `-ksp_guess_pod_Ainner true`     - Use the operator as inner product (must be SPD)

Level: developer

-seealso: [](ch_ksp), `KSPGuess`, `KSPGetGuess()`, `KSPGuessSetType()`, `KSPGuessType`

# External Links
$(_doc_external("Ksp/KSPGuessSetFromOptions"))
"""
function KSPGuessSetFromOptions(petsclib::PetscLibType, guess::KSPGuess) end

@for_petsc function KSPGuessSetFromOptions(petsclib::$UnionPetscLib, guess::KSPGuess )

    @chk ccall(
               (:KSPGuessSetFromOptions, $petsc_library),
               PetscErrorCode,
               (KSPGuess,),
               guess,
              )


	return nothing
end 

"""
	KSPGuessSetTolerance(petsclib::PetscLibType,guess::KSPGuess, tol::PetscReal) 
Sets the relative tolerance used in either eigenvalue (POD) or singular value (Fischer type 3) calculations.

Collective

Input Parameters:
- `guess` - `KSPGuess` object
- `tol`   - the tolerance

Options Database Key:
- `-ksp_guess_fischer_tol <tol>` - set the tolerance for the Fischer models
- `-ksp_guess_pod_tol <tol>`     - set the tolerance for the Pod models

Level: developer

-seealso: [](ch_ksp), `KSPGuess`, `KSPGuessType`, `KSPGuessSetFromOptions()`

# External Links
$(_doc_external("Ksp/KSPGuessSetTolerance"))
"""
function KSPGuessSetTolerance(petsclib::PetscLibType, guess::KSPGuess, tol::PetscReal) end

@for_petsc function KSPGuessSetTolerance(petsclib::$UnionPetscLib, guess::KSPGuess, tol::$PetscReal )

    @chk ccall(
               (:KSPGuessSetTolerance, $petsc_library),
               PetscErrorCode,
               (KSPGuess, $PetscReal),
               guess, tol,
              )


	return nothing
end 

"""
	KSPGuessDestroy(petsclib::PetscLibType,guess::KSPGuess) 
Destroys `KSPGuess` context.

Collective

Input Parameter:
- `guess` - initial guess object

Level: developer

-seealso: [](ch_ksp), `KSPGuessCreate()`, `KSPGuess`, `KSPGuessType`

# External Links
$(_doc_external("Ksp/KSPGuessDestroy"))
"""
function KSPGuessDestroy(petsclib::PetscLibType, guess::KSPGuess) end

@for_petsc function KSPGuessDestroy(petsclib::$UnionPetscLib, guess::KSPGuess )

    @chk ccall(
               (:KSPGuessDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{KSPGuess},),
               guess,
              )


	return nothing
end 

"""
	KSPGuessView(petsclib::PetscLibType,guess::KSPGuess, view::PetscViewer) 
View the `KSPGuess` object

Logically Collective

Input Parameters:
- `guess` - the initial guess object for the Krylov method
- `view`  - the viewer object

Options Database Key:
- `-ksp_guess_view viewer` - view the `KSPGuess` object

Level: developer

-seealso: [](ch_ksp), `KSP`, `KSPGuess`, `KSPGuessType`, `KSPGuessRegister()`, `KSPGuessCreate()`, `PetscViewer`

# External Links
$(_doc_external("Ksp/KSPGuessView"))
"""
function KSPGuessView(petsclib::PetscLibType, guess::KSPGuess, view::PetscViewer) end

@for_petsc function KSPGuessView(petsclib::$UnionPetscLib, guess::KSPGuess, view::PetscViewer )

    @chk ccall(
               (:KSPGuessView, $petsc_library),
               PetscErrorCode,
               (KSPGuess, PetscViewer),
               guess, view,
              )


	return nothing
end 

"""
	guess::KSPGuess = KSPGuessCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Creates a `KSPGuess` context.

Collective

Input Parameter:
- `comm` - MPI communicator

Output Parameter:
- `guess` - location to put the `KSPGuess` context

Options Database Keys:
- `-ksp_guess_type  <method>`      - Turns on generation of initial guesses and sets the method; use -help for a list of available methods
- `-ksp_guess_view <viewer>`       - view the `KSPGuess` object
- `-ksp_guess_fischer_model <a,b>` - set details for the Fischer models
- `-ksp_guess_fischer_monitor`     - monitor the fischer models
- `-ksp_guess_fischer_tol <tol>`   - set the tolerance for the Fischer models
- `-ksp_guess_pod_size <size>`     - Number of snapshots
- `-ksp_guess_pod_monitor true`    - monitor the pod initial guess processing
- `-ksp_guess_pod_tol <tol>`       - Tolerance to retain eigenvectors
- `-ksp_guess_pod_Ainner true`     - Use the operator as inner product (must be SPD)

Level: developer

-seealso: [](ch_ksp), `KSPSolve()`, `KSPGuessDestroy()`, `KSPGuess`, `KSPGuessType`, `KSP`

# External Links
$(_doc_external("Ksp/KSPGuessCreate"))
"""
function KSPGuessCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function KSPGuessCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	guess_ = Ref{KSPGuess}()

    @chk ccall(
               (:KSPGuessCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{KSPGuess}),
               comm, guess_,
              )

	guess = guess_[]

	return guess
end 

"""
	KSPGuessSetType(petsclib::PetscLibType,guess::KSPGuess, type::KSPGuessType) 
Sets the type of a `KSPGuess`. Each `KSPGuessType` provides a different algorithm for computing the initial guess.

Logically Collective

Input Parameters:
- `guess` - the initial guess object for the Krylov method
- `type`  - a known `KSPGuessType`

Options Database Key:
- `-ksp_guess_type  <method>` - Turns on generation of initial guesses and sets the method; see `KSPGuessType` for a list of available types

Level: developer

-seealso: [](ch_ksp), `KSP`, `KSPGuess`, `KSPGuessType`, `KSPGuessRegister()`, `KSPGuessCreate()`, `KSPGUESSFISCHER`, `KSPGUESSPOD`

# External Links
$(_doc_external("Ksp/KSPGuessSetType"))
"""
function KSPGuessSetType(petsclib::PetscLibType, guess::KSPGuess, type::KSPGuessType) end

@for_petsc function KSPGuessSetType(petsclib::$UnionPetscLib, guess::KSPGuess, type::KSPGuessType )

    @chk ccall(
               (:KSPGuessSetType, $petsc_library),
               PetscErrorCode,
               (KSPGuess, KSPGuessType),
               guess, type,
              )


	return nothing
end 

"""
	type::KSPGuessType = KSPGuessGetType(petsclib::PetscLibType,guess::KSPGuess) 
Gets the `KSPGuessType` as a string from the `KSPGuess` object.

Not Collective

Input Parameter:
- `guess` - the initial guess context

Output Parameter:
- `type` - type of `KSPGuess` method

Level: developer

-seealso: [](ch_ksp), `KSPGuess`, `KSPGuessSetType()`

# External Links
$(_doc_external("Ksp/KSPGuessGetType"))
"""
function KSPGuessGetType(petsclib::PetscLibType, guess::KSPGuess) end

@for_petsc function KSPGuessGetType(petsclib::$UnionPetscLib, guess::KSPGuess )
	type_ = Ref{KSPGuessType}()

    @chk ccall(
               (:KSPGuessGetType, $petsc_library),
               PetscErrorCode,
               (KSPGuess, Ptr{KSPGuessType}),
               guess, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	KSPGuessUpdate(petsclib::PetscLibType,guess::KSPGuess, rhs::PetscVec, sol::PetscVec) 
Updates the guess object with the current solution and rhs vector

Collective

Input Parameters:
- `guess` - the initial guess context
- `rhs`   - the corresponding rhs
- `sol`   - the computed solution

Level: developer

-seealso: [](ch_ksp), `KSPGuessCreate()`, `KSPGuess`

# External Links
$(_doc_external("Ksp/KSPGuessUpdate"))
"""
function KSPGuessUpdate(petsclib::PetscLibType, guess::KSPGuess, rhs::PetscVec, sol::PetscVec) end

@for_petsc function KSPGuessUpdate(petsclib::$UnionPetscLib, guess::KSPGuess, rhs::PetscVec, sol::PetscVec )

    @chk ccall(
               (:KSPGuessUpdate, $petsc_library),
               PetscErrorCode,
               (KSPGuess, CVec, CVec),
               guess, rhs, sol,
              )


	return nothing
end 

"""
	KSPGuessFormGuess(petsclib::PetscLibType,guess::KSPGuess, rhs::PetscVec, sol::PetscVec) 
Form the initial guess

Collective

Input Parameters:
- `guess` - the initial guess context
- `rhs`   - the current right-hand side vector
- `sol`   - the initial guess vector

Level: developer

-seealso: [](ch_ksp), `KSPGuessCreate()`, `KSPGuess`

# External Links
$(_doc_external("Ksp/KSPGuessFormGuess"))
"""
function KSPGuessFormGuess(petsclib::PetscLibType, guess::KSPGuess, rhs::PetscVec, sol::PetscVec) end

@for_petsc function KSPGuessFormGuess(petsclib::$UnionPetscLib, guess::KSPGuess, rhs::PetscVec, sol::PetscVec )

    @chk ccall(
               (:KSPGuessFormGuess, $petsc_library),
               PetscErrorCode,
               (KSPGuess, CVec, CVec),
               guess, rhs, sol,
              )


	return nothing
end 

"""
	KSPGuessSetUp(petsclib::PetscLibType,guess::KSPGuess) 
Setup the initial guess object

Collective

Input Parameter:
- `guess` - the initial guess context

Level: developer

-seealso: [](ch_ksp), `KSPGuessCreate()`, `KSPGuess`

# External Links
$(_doc_external("Ksp/KSPGuessSetUp"))
"""
function KSPGuessSetUp(petsclib::PetscLibType, guess::KSPGuess) end

@for_petsc function KSPGuessSetUp(petsclib::$UnionPetscLib, guess::KSPGuess )

    @chk ccall(
               (:KSPGuessSetUp, $petsc_library),
               PetscErrorCode,
               (KSPGuess,),
               guess,
              )


	return nothing
end 

"""
	KSPGuessFischerSetModel(petsclib::PetscLibType,guess::KSPGuess, model::PetscInt, size::PetscInt) 
Set the Paul Fischer algorithm or its variants to compute the initial guess for a `KSPSolve()`

Logically Collective

Input Parameters:
- `guess` - the initial guess context
- `model` - use model 1, model 2, model 3, or any other number to turn it off
- `size`  - size of subspace used to generate initial guess

Options Database Key:
- `-ksp_guess_fischer_model <model,size>` - uses the Fischer initial guess generator for repeated linear solves

Level: advanced

-seealso: [](ch_ksp), `KSPGuess`, `KSPGuessCreate()`, `KSPSetUseFischerGuess()`, `KSPSetGuess()`, `KSPGetGuess()`, `KSP`

# External Links
$(_doc_external("Ksp/KSPGuessFischerSetModel"))
"""
function KSPGuessFischerSetModel(petsclib::PetscLibType, guess::KSPGuess, model::PetscInt, size::PetscInt) end

@for_petsc function KSPGuessFischerSetModel(petsclib::$UnionPetscLib, guess::KSPGuess, model::$PetscInt, size::$PetscInt )

    @chk ccall(
               (:KSPGuessFischerSetModel, $petsc_library),
               PetscErrorCode,
               (KSPGuess, $PetscInt, $PetscInt),
               guess, model, size,
              )


	return nothing
end 

