# autodefined type arguments for class ------
mutable struct _n_PetscDLLibrary end
const PetscDLLibrary = Ptr{_n_PetscDLLibrary}

# -------------------------------------------------------
"""
	PetscDLLibraryPrintPath(petsclib::PetscLibType,libs::PetscDLLibrary) 

# External Links
$(_doc_external("Sys/PetscDLLibraryPrintPath"))
"""
function PetscDLLibraryPrintPath(petsclib::PetscLibType, libs::PetscDLLibrary) end

@for_petsc function PetscDLLibraryPrintPath(petsclib::$UnionPetscLib, libs::PetscDLLibrary )

    @chk ccall(
               (:PetscDLLibraryPrintPath, $petsc_library),
               PetscErrorCode,
               (PetscDLLibrary,),
               libs,
              )


	return nothing
end 

"""
	found::PetscBool = PetscDLLibraryRetrieve(petsclib::PetscLibType,comm::MPI_Comm, libname::String, lname::String, llen::Csize_t) 
Copies a PETSc dynamic library from a remote location
(if it is remote), then indicates if it exits and its local name.

Collective

Input Parameters:
- `comm`    - MPI processes that will be opening the library
- `libname` - name of the library, can be a relative or absolute path and be a URL
- `llen`    - length of the `name` buffer

Output Parameters:
- `lname` - actual name of the file on local filesystem if `found`
- `found` - true if the file exists

Level: developer

-seealso: `PetscFileRetrieve()`

# External Links
$(_doc_external("Sys/PetscDLLibraryRetrieve"))
"""
function PetscDLLibraryRetrieve(petsclib::PetscLibType, comm::MPI_Comm, libname::String, lname::String, llen::Csize_t) end

@for_petsc function PetscDLLibraryRetrieve(petsclib::$UnionPetscLib, comm::MPI_Comm, libname::String, lname::String, llen::Csize_t )
	found_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscDLLibraryRetrieve, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Csize_t, Ptr{PetscBool}),
               comm, libname, lname, llen, found_,
              )

	found = found_[]

	return found
end 

"""
	PetscDLLibraryOpen(petsclib::PetscLibType,comm::MPI_Comm, path::String, entry::PetscDLLibrary) 
Opens a PETSc dynamic link library

Collective, No Fortran Support

Input Parameters:
- `comm` - MPI processes that are opening the library
- `path` - name of the library, can be a relative or absolute path

Output Parameter:
- `entry` - a PETSc dynamic link library entry

Level: developer

-seealso: `PetscDLLibrary`, `PetscLoadDynamicLibrary()`, `PetscDLLibraryAppend()`, `PetscDLLibraryRetrieve()`, `PetscDLLibrarySym()`, `PetscDLLibraryClose()`

# External Links
$(_doc_external("Sys/PetscDLLibraryOpen"))
"""
function PetscDLLibraryOpen(petsclib::PetscLibType, comm::MPI_Comm, path::String, entry::PetscDLLibrary) end

@for_petsc function PetscDLLibraryOpen(petsclib::$UnionPetscLib, comm::MPI_Comm, path::String, entry::PetscDLLibrary )

    @chk ccall(
               (:PetscDLLibraryOpen, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{PetscDLLibrary}),
               comm, path, entry,
              )


	return nothing
end 

"""
	PetscDLLibrarySym(petsclib::PetscLibType,comm::MPI_Comm, outlist::PetscDLLibrary, path::String, insymbol::String, value::Cvoid) 
Load a symbol from a list of dynamic link libraries.

Collective, No Fortran Support

Input Parameters:
- `comm`     - the MPI communicator that will load the symbol
- `outlist`  - list of already open libraries that may contain symbol (can be `NULL` and only the executable is searched for the function)
- `path`     - optional complete library name (if provided it checks here before checking `outlist`)
- `insymbol` - name of symbol

Output Parameter:
- `value` - if symbol not found then this value is set to `NULL`

Level: developer

-seealso: `PetscDLLibrary`, `PetscLoadDynamicLibrary()`, `PetscDLLibraryAppend()`, `PetscDLLibraryRetrieve()`, `PetscDLLibraryOpen()`, `PetscDLLibraryClose()`

# External Links
$(_doc_external("Sys/PetscDLLibrarySym"))
"""
function PetscDLLibrarySym(petsclib::PetscLibType, comm::MPI_Comm, outlist::PetscDLLibrary, path::String, insymbol::String, value::Cvoid) end

@for_petsc function PetscDLLibrarySym(petsclib::$UnionPetscLib, comm::MPI_Comm, outlist::PetscDLLibrary, path::String, insymbol::String, value::Cvoid )

    @chk ccall(
               (:PetscDLLibrarySym, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscDLLibrary}, Ptr{Cchar}, Ptr{Cchar}, Cvoid),
               comm, outlist, path, insymbol, value,
              )


	return nothing
end 

"""
	PetscDLLibraryAppend(petsclib::PetscLibType,comm::MPI_Comm, outlist::PetscDLLibrary, path::String) 
Appends another dynamic link library to the end  of the search list

Collective, No Fortran Support

Input Parameters:
- `comm` - MPI communicator
- `path` - name of the library

Output Parameter:
- `outlist` - list of libraries

Level: developer

-seealso: `PetscDLLibrary`, `PetscDLLibraryOpen()`, `PetscLoadDynamicLibrary()`, `PetscDLLibraryRetrieve()`, `PetscDLLibraryPrepend()`

# External Links
$(_doc_external("Sys/PetscDLLibraryAppend"))
"""
function PetscDLLibraryAppend(petsclib::PetscLibType, comm::MPI_Comm, outlist::PetscDLLibrary, path::String) end

@for_petsc function PetscDLLibraryAppend(petsclib::$UnionPetscLib, comm::MPI_Comm, outlist::PetscDLLibrary, path::String )

    @chk ccall(
               (:PetscDLLibraryAppend, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscDLLibrary}, Ptr{Cchar}),
               comm, outlist, path,
              )


	return nothing
end 

"""
	PetscDLLibraryPrepend(petsclib::PetscLibType,comm::MPI_Comm, outlist::PetscDLLibrary, path::String) 
Add another dynamic library to search for symbols to the beginning of the search list

Collective, No Fortran Support

Input Parameters:
- `comm` - MPI communicator
- `path` - name of the library

Output Parameter:
- `outlist` - list of libraries

Level: developer

-seealso: `PetscDLLibrary`, `PetscDLLibraryOpen()`, `PetscLoadDynamicLibrary()`, `PetscDLLibraryRetrieve()`, `PetscDLLibraryAppend()`

# External Links
$(_doc_external("Sys/PetscDLLibraryPrepend"))
"""
function PetscDLLibraryPrepend(petsclib::PetscLibType, comm::MPI_Comm, outlist::PetscDLLibrary, path::String) end

@for_petsc function PetscDLLibraryPrepend(petsclib::$UnionPetscLib, comm::MPI_Comm, outlist::PetscDLLibrary, path::String )

    @chk ccall(
               (:PetscDLLibraryPrepend, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscDLLibrary}, Ptr{Cchar}),
               comm, outlist, path,
              )


	return nothing
end 

"""
	PetscDLLibraryClose(petsclib::PetscLibType,list::PetscDLLibrary) 
Destroys the search path of dynamic libraries and closes the libraries.

Collective, No Fortran Support

Input Parameter:
- `list` - library list

Level: developer

-seealso: `PetscDLLibrary`, `PetscDLLibraryOpen()`, `PetscLoadDynamicLibrary()`, `PetscDLLibraryRetrieve()`, `PetscDLLibraryAppend()`,
`PetscDLLibraryPrepend()`

# External Links
$(_doc_external("Sys/PetscDLLibraryClose"))
"""
function PetscDLLibraryClose(petsclib::PetscLibType, list::PetscDLLibrary) end

@for_petsc function PetscDLLibraryClose(petsclib::$UnionPetscLib, list::PetscDLLibrary )

    @chk ccall(
               (:PetscDLLibraryClose, $petsc_library),
               PetscErrorCode,
               (PetscDLLibrary,),
               list,
              )


	return nothing
end 

