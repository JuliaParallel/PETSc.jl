"""
	outlist::PetscDLLibrary = PetscDLLibraryAppend(petsclib::PetscLibType, comm::MPI_Comm, path::String) 
Appends another dynamic link library to the end  of the search list

Collective, No Fortran Support

Input Parameters:
- `comm` - MPI communicator
- `path` - name of the library

Output Parameter:
- `outlist` - list of libraries

Level: developer

See also: `PetscDLLibrary`, `PetscDLLibraryOpen()`, `PetscLoadDynamicLibrary()`, `PetscDLLibraryRetrieve()`, `PetscDLLibraryPrepend()`

# External Links
$(_doc_external("Sys/PetscDLLibraryAppend"))
"""
function PetscDLLibraryAppend(petsclib::PetscLibType, comm::MPI_Comm, path::String)
    error("PetscDLLibraryAppend: no generated method for these argument types")
end

@for_petsc function PetscDLLibraryAppend(petsclib::$UnionPetscLib, comm::MPI_Comm, path::String )
	outlist_ = Ref{PetscDLLibrary}()

    @chk ccall(
               (:PetscDLLibraryAppend, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscDLLibrary}, Ptr{Cchar}),
               comm, outlist_, path,
              )

	outlist = outlist_[]

	return outlist
end 

"""
	PetscDLLibraryClose(petsclib::PetscLibType, list::PetscDLLibrary) 
Destroys the search path of dynamic libraries and closes the libraries.

Collective, No Fortran Support

Input Parameter:
- `list` - library list

Level: developer

See also: `PetscDLLibrary`, `PetscDLLibraryOpen()`, `PetscLoadDynamicLibrary()`, `PetscDLLibraryRetrieve()`, `PetscDLLibraryAppend()`,
`PetscDLLibraryPrepend()`

# External Links
$(_doc_external("Sys/PetscDLLibraryClose"))
"""
function PetscDLLibraryClose(petsclib::PetscLibType, list::PetscDLLibrary)
    error("PetscDLLibraryClose: no generated method for these argument types")
end

@for_petsc function PetscDLLibraryClose(petsclib::$UnionPetscLib, list::PetscDLLibrary )

    @chk ccall(
               (:PetscDLLibraryClose, $petsc_library),
               PetscErrorCode,
               (PetscDLLibrary,),
               list,
              )


	return nothing
end 

"""
	entry::PetscDLLibrary = PetscDLLibraryOpen(petsclib::PetscLibType, comm::MPI_Comm, path::String) 
Opens a PETSc dynamic link library

Collective, No Fortran Support

Input Parameters:
- `comm` - MPI processes that are opening the library
- `path` - name of the library, can be a relative or absolute path

Output Parameter:
- `entry` - a PETSc dynamic link library entry

Level: developer

See also: `PetscDLLibrary`, `PetscLoadDynamicLibrary()`, `PetscDLLibraryAppend()`, `PetscDLLibraryRetrieve()`, `PetscDLLibrarySym()`, `PetscDLLibraryClose()`

# External Links
$(_doc_external("Sys/PetscDLLibraryOpen"))
"""
function PetscDLLibraryOpen(petsclib::PetscLibType, comm::MPI_Comm, path::String)
    error("PetscDLLibraryOpen: no generated method for these argument types")
end

@for_petsc function PetscDLLibraryOpen(petsclib::$UnionPetscLib, comm::MPI_Comm, path::String )
	entry_ = Ref{PetscDLLibrary}()

    @chk ccall(
               (:PetscDLLibraryOpen, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{PetscDLLibrary}),
               comm, path, entry_,
              )

	entry = entry_[]

	return entry
end 

"""
	outlist::PetscDLLibrary = PetscDLLibraryPrepend(petsclib::PetscLibType, comm::MPI_Comm, path::String) 
Add another dynamic library to search for symbols to the beginning of the search list

Collective, No Fortran Support

Input Parameters:
- `comm` - MPI communicator
- `path` - name of the library

Output Parameter:
- `outlist` - list of libraries

Level: developer

See also: `PetscDLLibrary`, `PetscDLLibraryOpen()`, `PetscLoadDynamicLibrary()`, `PetscDLLibraryRetrieve()`, `PetscDLLibraryAppend()`

# External Links
$(_doc_external("Sys/PetscDLLibraryPrepend"))
"""
function PetscDLLibraryPrepend(petsclib::PetscLibType, comm::MPI_Comm, path::String)
    error("PetscDLLibraryPrepend: no generated method for these argument types")
end

@for_petsc function PetscDLLibraryPrepend(petsclib::$UnionPetscLib, comm::MPI_Comm, path::String )
	outlist_ = Ref{PetscDLLibrary}()

    @chk ccall(
               (:PetscDLLibraryPrepend, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscDLLibrary}, Ptr{Cchar}),
               comm, outlist_, path,
              )

	outlist = outlist_[]

	return outlist
end 

"""
	PetscDLLibraryPrintPath(petsclib::PetscLibType, libs::PetscDLLibrary) 
Prints the names of all dynamic libraries in a `PetscDLLibrary` list to the PETSc error output stream

Not Collective

Input Parameter:
- `libs` - the linked list of currently opened dynamic libraries

Level: developer

See also: `PetscDLLibrary`, `PetscDLLibraryOpen()`, `PetscDLLibrarySym()`, `PetscDLLibraryAppend()`, `PetscDLLibraryClose()`

# External Links
$(_doc_external("Sys/PetscDLLibraryPrintPath"))
"""
function PetscDLLibraryPrintPath(petsclib::PetscLibType, libs::PetscDLLibrary)
    error("PetscDLLibraryPrintPath: no generated method for these argument types")
end

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
	found::PetscBool = PetscDLLibraryRetrieve(petsclib::PetscLibType, comm::MPI_Comm, libname::String, lname::String, llen::Csize_t) 
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

See also: `PetscFileRetrieve()`

# External Links
$(_doc_external("Sys/PetscDLLibraryRetrieve"))
"""
function PetscDLLibraryRetrieve(petsclib::PetscLibType, comm::MPI_Comm, libname::String, lname::String, llen::Csize_t)
    error("PetscDLLibraryRetrieve: no generated method for these argument types")
end

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
	value::Ptr{Cvoid} = PetscDLLibrarySym(petsclib::PetscLibType, comm::MPI_Comm, outlist::PetscDLLibrary, path::String, insymbol::String) 
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

See also: `PetscDLLibrary`, `PetscLoadDynamicLibrary()`, `PetscDLLibraryAppend()`, `PetscDLLibraryRetrieve()`, `PetscDLLibraryOpen()`, `PetscDLLibraryClose()`

# External Links
$(_doc_external("Sys/PetscDLLibrarySym"))
"""
function PetscDLLibrarySym(petsclib::PetscLibType, comm::MPI_Comm, outlist::PetscDLLibrary, path::String, insymbol::String)
    error("PetscDLLibrarySym: no generated method for these argument types")
end

@for_petsc function PetscDLLibrarySym(petsclib::$UnionPetscLib, comm::MPI_Comm, outlist::PetscDLLibrary, path::String, insymbol::String )
	value_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:PetscDLLibrarySym, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscDLLibrary}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Ptr{Cvoid}}),
               comm, outlist, path, insymbol, value_,
              )

	value = value_[]

	return value
end 

