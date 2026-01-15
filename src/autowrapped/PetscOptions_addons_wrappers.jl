# autodefined type arguments for class ------
mutable struct _n_PetscOptionsHelpPrCinted end
const PetscOptionsHelpPrCinted = Ptr{_n_PetscOptionsHelpPrCinted}

# -------------------------------------------------------
# autodefined type arguments for class ------
# -------------------------------------------------------
# autodefined type arguments for class ------
# -------------------------------------------------------
"""
	PetscOptionsHelpPrintedDestroy(petsclib::PetscLibType,hp::PetscOptionsHelpPrCinted) 

# External Links
$(_doc_external("Sys/PetscOptionsHelpPrintedDestroy"))
"""
function PetscOptionsHelpPrintedDestroy(petsclib::PetscLibType, hp::PetscOptionsHelpPrCinted) end

@for_petsc function PetscOptionsHelpPrintedDestroy(petsclib::$UnionPetscLib, hp::PetscOptionsHelpPrCinted )

    @chk ccall(
               (:PetscOptionsHelpPrintedDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscOptionsHelpPrCinted},),
               hp,
              )


	return nothing
end 

"""
	hp::PetscOptionsHelpPrCinted = PetscOptionsHelpPrintedCreate(petsclib::PetscLibType) 
Creates an object used to manage tracking which help messages have
been printed so they will not be printed again.

Output Parameter:
- `hp` - the created object

Not Collective

Level: developer

-seealso: `PetscOptionsHelpPrintedCheck()`, `PetscOptionsHelpPrintChecked()`

# External Links
$(_doc_external("Sys/PetscOptionsHelpPrintedCreate"))
"""
function PetscOptionsHelpPrintedCreate(petsclib::PetscLibType) end

@for_petsc function PetscOptionsHelpPrintedCreate(petsclib::$UnionPetscLib)
	hp_ = Ref{PetscOptionsHelpPrCinted}()

    @chk ccall(
               (:PetscOptionsHelpPrintedCreate, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscOptionsHelpPrCinted},),
               hp_,
              )

	hp = hp_[]

	return hp
end 

"""
	found::PetscBool = PetscOptionsHelpPrintedCheck(petsclib::PetscLibType,hp::PetscOptionsHelpPrCinted, pre::String, name::String) 
Checks if a particular pre, name pair has previous been entered (meaning the help message was printed)

Not Collective

Input Parameters:
- `hp`   - the object used to manage tracking what help messages have been printed
- `pre`  - the prefix part of the string, many be `NULL`
- `name` - the string to look for (cannot be `NULL`)

Output Parameter:
- `found` - `PETSC_TRUE` if the string was already set

Level: intermediate

-seealso: `PetscOptionsHelpPrintedCreate()`

# External Links
$(_doc_external("Sys/PetscOptionsHelpPrintedCheck"))
"""
function PetscOptionsHelpPrintedCheck(petsclib::PetscLibType, hp::PetscOptionsHelpPrCinted, pre::String, name::String) end

@for_petsc function PetscOptionsHelpPrintedCheck(petsclib::$UnionPetscLib, hp::PetscOptionsHelpPrCinted, pre::String, name::String )
	found_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsHelpPrintedCheck, $petsc_library),
               PetscErrorCode,
               (PetscOptionsHelpPrCinted, Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscBool}),
               hp, pre, name, found_,
              )

	found = found_[]

	return found
end 

