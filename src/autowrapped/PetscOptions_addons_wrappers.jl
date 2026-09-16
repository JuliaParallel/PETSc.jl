"""
	found::PetscBool = PetscOptionsHelpPrintedCheck(petsclib::PetscLibType, hp::PetscOptionsHelpPrCinted, pre::String, name::String) 
Checks if a particular pre, name pair has previous been entered (meaning the help message was printed)

Not Collective

Input Parameters:
- `hp`   - the object used to manage tracking what help messages have been printed
- `pre`  - the prefix part of the string, many be `NULL`
- `name` - the string to look for (cannot be `NULL`)

Output Parameter:
- `found` - `PETSC_TRUE` if the string was already set

Level: intermediate

See also: `PetscOptionsHelpPrintedCreate()`

# External Links
$(_doc_external("Viewer/PetscOptionsHelpPrintedCheck"))
"""
function PetscOptionsHelpPrintedCheck(petsclib::PetscLibType, hp::PetscOptionsHelpPrCinted, pre::String, name::String)
    error("PetscOptionsHelpPrintedCheck: no generated method for these argument types")
end

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

"""
	hp::PetscOptionsHelpPrCinted = PetscOptionsHelpPrintedCreate(petsclib::PetscLibType) 
Creates an object used to manage tracking which help messages have
been printed so they will not be printed again.

Output Parameter:
- `hp` - the created object

Not Collective

Level: developer

See also: `PetscOptionsHelpPrintedCheck()`, `PetscOptionsHelpPrintChecked()`

# External Links
$(_doc_external("Viewer/PetscOptionsHelpPrintedCreate"))
"""
function PetscOptionsHelpPrintedCreate(petsclib::PetscLibType)
    error("PetscOptionsHelpPrintedCreate: no generated method for these argument types")
end

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
	PetscOptionsHelpPrintedDestroy(petsclib::PetscLibType, hp::Union{PetscOptionsHelpPrCinted, Ref{PetscOptionsHelpPrCinted}}) 
Destroys the object used to track which help messages have already been printed

Not Collective

Input Parameter:
- `hp` - pointer to the `PetscOptionsHelpPrinted` object to destroy; set to `NULL` on return

Level: developer

See also: `PetscOptionsHelpPrintedCreate()`, `PetscOptionsHelpPrintedCheck()`

# External Links
$(_doc_external("Viewer/PetscOptionsHelpPrintedDestroy"))
"""
function PetscOptionsHelpPrintedDestroy(petsclib::PetscLibType, hp::Union{PetscOptionsHelpPrCinted, Ref{PetscOptionsHelpPrCinted}})
    error("PetscOptionsHelpPrintedDestroy: no generated method for these argument types")
end

@for_petsc function PetscOptionsHelpPrintedDestroy(petsclib::$UnionPetscLib, hp::Union{PetscOptionsHelpPrCinted, Ref{PetscOptionsHelpPrCinted}} )
	hp_ = hp isa Base.RefValue ? hp : Ref{PetscOptionsHelpPrCinted}(hp)

    @chk ccall(
               (:PetscOptionsHelpPrintedDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscOptionsHelpPrCinted},),
               hp_,
              )


	return nothing
end 

