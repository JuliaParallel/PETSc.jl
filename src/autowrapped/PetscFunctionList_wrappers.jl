"""
	PetscFunctionListDestroy(petsclib::PetscLibType,fl::PetscFunctionList) 
Destroys a list of registered routines.

Input Parameter:
- `fl` - pointer to list

Level: developer

-seealso: `PetscFunctionListAdd()`, `PetscFunctionList`, `PetscFunctionListClear()`

# External Links
$(_doc_external("Sys/PetscFunctionListDestroy"))
"""
function PetscFunctionListDestroy(petsclib::PetscLibType, fl::PetscFunctionList) end

@for_petsc function PetscFunctionListDestroy(petsclib::$UnionPetscLib, fl::PetscFunctionList )

    @chk ccall(
               (:PetscFunctionListDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscFunctionList},),
               fl,
              )


	return nothing
end 

"""
	PetscFunctionListClear(petsclib::PetscLibType,fl::PetscFunctionList) 
Clear a `PetscFunctionList`

Not Collective

Input Parameter:
- `fl` - The `PetscFunctionList` to clear

Level: developer

-seealso: `PetscFunctionList`, `PetscFunctionListDestroy()`, `PetscFunctionListAdd()`

# External Links
$(_doc_external("Sys/PetscFunctionListClear"))
"""
function PetscFunctionListClear(petsclib::PetscLibType, fl::PetscFunctionList) end

@for_petsc function PetscFunctionListClear(petsclib::$UnionPetscLib, fl::PetscFunctionList )

    @chk ccall(
               (:PetscFunctionListClear, $petsc_library),
               PetscErrorCode,
               (PetscFunctionList,),
               fl,
              )


	return nothing
end 

"""
	PetscFunctionListPrintAll(petsclib::PetscLibType) 

# External Links
$(_doc_external("Sys/PetscFunctionListPrintAll"))
"""
function PetscFunctionListPrintAll(petsclib::PetscLibType) end

@for_petsc function PetscFunctionListPrintAll(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscFunctionListPrintAll, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscFunctionListPrintNonEmpty(petsclib::PetscLibType,fl::PetscFunctionList) 
Print composed names for non `NULL` function pointers

Logically Collective, No Fortran Support

Input Parameter:
- `fl` - the function list

Level: developer

-seealso: `PetscFunctionListAdd()`, `PetscFunctionList`, `PetscObjectQueryFunction()`

# External Links
$(_doc_external("Sys/PetscFunctionListPrintNonEmpty"))
"""
function PetscFunctionListPrintNonEmpty(petsclib::PetscLibType, fl::PetscFunctionList) end

@for_petsc function PetscFunctionListPrintNonEmpty(petsclib::$UnionPetscLib, fl::PetscFunctionList )

    @chk ccall(
               (:PetscFunctionListPrintNonEmpty, $petsc_library),
               PetscErrorCode,
               (PetscFunctionList,),
               fl,
              )


	return nothing
end 

"""
	PetscFunctionListView(petsclib::PetscLibType,list::PetscFunctionList, viewer::PetscViewer) 
prints out contents of a `PetscFunctionList`

Collective

Input Parameters:
- `list`   - the list of functions
- `viewer` - the `PetscViewer` used to view the `PetscFunctionList`

Level: developer

-seealso: `PetscFunctionListAdd()`, `PetscFunctionListPrintTypes()`, `PetscFunctionList`

# External Links
$(_doc_external("Sys/PetscFunctionListView"))
"""
function PetscFunctionListView(petsclib::PetscLibType, list::PetscFunctionList, viewer::PetscViewer) end

@for_petsc function PetscFunctionListView(petsclib::$UnionPetscLib, list::PetscFunctionList, viewer::PetscViewer )

    @chk ccall(
               (:PetscFunctionListView, $petsc_library),
               PetscErrorCode,
               (PetscFunctionList, PetscViewer),
               list, viewer,
              )


	return nothing
end 

"""
	PetscFunctionListGet(petsclib::PetscLibType,list::PetscFunctionList, array::Cchar, n::Cint) 
Gets an array the contains the entries in `PetscFunctionList`, this is used
by help etc.

Not Collective, No Fortran Support

Input Parameter:
- `list` - list of types

Output Parameters:
- `array` - array of names
- `n`     - length of `array`

Level: developer

-seealso: `PetscFunctionListAdd()`, `PetscFunctionList`

# External Links
$(_doc_external("Sys/PetscFunctionListGet"))
"""
function PetscFunctionListGet(petsclib::PetscLibType, list::PetscFunctionList, array::Cchar, n::Cint) end

@for_petsc function PetscFunctionListGet(petsclib::$UnionPetscLib, list::PetscFunctionList, array::Cchar, n::Cint )

    @chk ccall(
               (:PetscFunctionListGet, $petsc_library),
               PetscErrorCode,
               (PetscFunctionList, Cchar, Ptr{Cint}),
               list, array, n,
              )


	return nothing
end 

"""
	PetscFunctionListPrintTypes(petsclib::PetscLibType,comm::MPI_Comm, fd::Libc.FILE, prefix::String, name::String, text::String, man::String, list::PetscFunctionList, def::String, newv::String) 
Prints the methods available in a list of functions

Collective, No Fortran Support

Input Parameters:
- `comm`   - the communicator (usually `MPI_COMM_WORLD`)
- `fd`     - file to print to, usually `stdout`
- `prefix` - prefix to prepend to name (optional)
- `name`   - option string (for example, `-ksp_type`)
- `text`   - short description of the object (for example, "Krylov solvers")
- `man`    - name of manual page that discusses the object (for example, `KSPCreate`)
- `list`   - list of types
- `def`    - default (current) value
- `newv`   - new value

Level: developer

-seealso: `PetscFunctionListAdd()`, `PetscFunctionList`

# External Links
$(_doc_external("Sys/PetscFunctionListPrintTypes"))
"""
function PetscFunctionListPrintTypes(petsclib::PetscLibType, comm::MPI_Comm, fd::Libc.FILE, prefix::String, name::String, text::String, man::String, list::PetscFunctionList, def::String, newv::String) end

@for_petsc function PetscFunctionListPrintTypes(petsclib::$UnionPetscLib, comm::MPI_Comm, fd::Libc.FILE, prefix::String, name::String, text::String, man::String, list::PetscFunctionList, def::String, newv::String )

    @chk ccall(
               (:PetscFunctionListPrintTypes, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Libc.FILE}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, PetscFunctionList, Ptr{Cchar}, Ptr{Cchar}),
               comm, fd, prefix, name, text, man, list, def, newv,
              )


	return nothing
end 

"""
	nl::PetscFunctionList = PetscFunctionListDuplicate(petsclib::PetscLibType,fl::PetscFunctionList) 
Creates a new list from a given function list `PetscFunctionList`.

Input Parameter:
- `fl` - pointer to list

Output Parameter:
- `nl` - the new list (should point to `NULL` to start, otherwise appends)

Level: developer

-seealso: `PetscFunctionList`, `PetscFunctionListAdd()`, `PetscFlistDestroy()`

# External Links
$(_doc_external("Sys/PetscFunctionListDuplicate"))
"""
function PetscFunctionListDuplicate(petsclib::PetscLibType, fl::PetscFunctionList) end

@for_petsc function PetscFunctionListDuplicate(petsclib::$UnionPetscLib, fl::PetscFunctionList )
	nl_ = Ref{PetscFunctionList}()

    @chk ccall(
               (:PetscFunctionListDuplicate, $petsc_library),
               PetscErrorCode,
               (PetscFunctionList, Ptr{PetscFunctionList}),
               fl, nl_,
              )

	nl = nl_[]

	return nl
end 

