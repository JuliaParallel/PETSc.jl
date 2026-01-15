# autodefined type arguments for class ------
# -------------------------------------------------------

"""
	PetscOptionsPushCreateViewerOff(petsclib::PetscLibType,flg::PetscBool) 
sets if `PetscOptionsCreateViewer()`, `PetscOptionsViewer()`, and `PetscOptionsCreateViewers()` return viewers.

Logically Collective

Input Parameter:
- `flg` - `PETSC_TRUE` to turn off viewer creation, `PETSC_FALSE` to turn it on.

Level: developer

-seealso: [](sec_viewers), `PetscOptionsCreateViewer()`, `PetscOptionsPopCreateViewerOff()`

# External Links
$(_doc_external("Sys/PetscOptionsPushCreateViewerOff"))
"""
function PetscOptionsPushCreateViewerOff(petsclib::PetscLibType, flg::PetscBool) end

@for_petsc function PetscOptionsPushCreateViewerOff(petsclib::$UnionPetscLib, flg::PetscBool )

    @chk ccall(
               (:PetscOptionsPushCreateViewerOff, $petsc_library),
               PetscErrorCode,
               (PetscBool,),
               flg,
              )


	return nothing
end 

"""
	PetscOptionsPopCreateViewerOff(petsclib::PetscLibType) 
reset whether `PetscOptionsCreateViewer()` returns a viewer.

Logically Collective

Level: developer

-seealso: [](sec_viewers), `PetscOptionsCreateViewer()`, `PetscOptionsPushCreateViewerOff()`

# External Links
$(_doc_external("Sys/PetscOptionsPopCreateViewerOff"))
"""
function PetscOptionsPopCreateViewerOff(petsclib::PetscLibType) end

@for_petsc function PetscOptionsPopCreateViewerOff(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscOptionsPopCreateViewerOff, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	flg::PetscBool = PetscOptionsGetCreateViewerOff(petsclib::PetscLibType) 
do `PetscOptionsCreateViewer()`, `PetscOptionsViewer()`, and `PetscOptionsCreateViewers()` return viewers

Logically Collective

Output Parameter:
- `flg` - whether viewers are returned.

Level: developer

-seealso: [](sec_viewers), `PetscOptionsCreateViewer()`, `PetscOptionsPushCreateViewerOff()`, `PetscOptionsPopCreateViewerOff()`

# External Links
$(_doc_external("Sys/PetscOptionsGetCreateViewerOff"))
"""
function PetscOptionsGetCreateViewerOff(petsclib::PetscLibType) end

@for_petsc function PetscOptionsGetCreateViewerOff(petsclib::$UnionPetscLib)
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsGetCreateViewerOff, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscBool},),
               flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	viewer::PetscViewer,format::PetscViewerFormat,set::PetscBool = PetscOptionsCreateViewer(petsclib::PetscLibType,comm::MPI_Comm, options::Union{Ptr,PetscOptions}, pre::String, name::String) 
Creates a viewer appropriate for the type indicated by the user

Collective

Input Parameters:
- `comm`    - the communicator to own the viewer
- `options` - options database, use `NULL` for default global database
- `pre`     - the string to prepend to the name or `NULL`
- `name`    - the options database name that will be checked for

Output Parameters:
- `viewer` - the viewer, pass `NULL` if not needed
- `format` - the `PetscViewerFormat` requested by the user, pass `NULL` if not needed
- `set`    - `PETSC_TRUE` if found, else `PETSC_FALSE`

Level: intermediate

-seealso: [](sec_viewers), `PetscViewerDestroy()`, `PetscOptionsGetReal()`, `PetscOptionsHasName()`, `PetscOptionsGetString()`,
`PetscOptionsGetIntArray()`, `PetscOptionsGetRealArray()`, `PetscOptionsBool()`
`PetscOptionsInt()`, `PetscOptionsString()`, `PetscOptionsReal()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`, `PetscOptionsPushCreateViewerOff()`, `PetscOptionsPopCreateViewerOff()`,
`PetscOptionsCreateViewerOff()`

# External Links
$(_doc_external("Sys/PetscOptionsCreateViewer"))
"""
function PetscOptionsCreateViewer(petsclib::PetscLibType, comm::MPI_Comm, options::Union{Ptr,PetscOptions}, pre::String, name::String) end

@for_petsc function PetscOptionsCreateViewer(petsclib::$UnionPetscLib, comm::MPI_Comm, options::Union{Ptr,PetscOptions}, pre::String, name::String)
	viewer_ = Ref{PetscViewer}()
	format_ = Ref{PetscViewerFormat}()
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsCreateViewer, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, COptions, Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscViewer}, Ptr{PetscViewerFormat}, Ptr{PetscBool}),
               comm, options, pre, name, viewer_, format_, set_,
              )

	viewer = viewer_[]
	format = format_[]
	set = set_[]

	return viewer,format,set
end 

"""
	viewers::Vector{PetscViewer},formats::Vector{PetscViewerFormat},set::PetscBool = PetscOptionsCreateViewers(petsclib::PetscLibType,comm::MPI_Comm, options::PetscOptions, pre::String, name::String, n_max::PetscInt) 
Create multiple viewers from a comma

Collective

Input Parameters:
- `comm`    - the communicator to own the viewers
- `options` - options database, use `NULL` for default global database
- `pre`     - the string to prepend to the name or `NULL`
- `name`    - the options database name that will be checked for
- `n_max`   - on input: the maximum number of viewers; on output: the number of viewers in the comma-separated list

Output Parameters:
- `viewers` - an array to hold at least `n_max` `PetscViewer`s, or `NULL` if not needed; on output: if not `NULL`, the
first `n_max` entries are initialized `PetscViewer`s
- `formats` - an array to hold at least `n_max` `PetscViewerFormat`s, or `NULL` if not needed; on output: if not `NULL`, the first `n_max` entries are valid `PetscViewerFormat`s
- `set`     - `PETSC_TRUE` if found, else `PETSC_FALSE`

Level: intermediate

-seealso: [](sec_viewers), `PetscOptionsCreateViewer()`

# External Links
$(_doc_external("Sys/PetscOptionsCreateViewers"))
"""
function PetscOptionsCreateViewers(petsclib::PetscLibType, comm::MPI_Comm, options::PetscOptions, pre::String, name::String, n_max::PetscInt) end

@for_petsc function PetscOptionsCreateViewers(petsclib::$UnionPetscLib, comm::MPI_Comm, options::PetscOptions, pre::String, name::String, n_max::$PetscInt )
	viewers = Vector{PetscViewer}(undef, ni);  # CHECK SIZE!!
	formats = Vector{PetscViewerFormat}(undef, ni);  # CHECK SIZE!!
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsCreateViewers, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, COptions, Ptr{Cchar}, Ptr{Cchar}, Ptr{$PetscInt}, Ptr{PetscViewer}, Ptr{PetscViewerFormat}, Ptr{PetscBool}),
               comm, options, pre, name, n_max, viewers, formats, set_,
              )

	set = set_[]

	return viewers,formats,set
end 

"""
	PetscOptionsInsertStringYAML(petsclib::PetscLibType,options::PetscOptions, in_str::String) 
Inserts YAML

Logically Collective

Input Parameters:
- `options` - options database, use `NULL` for default global database
- `in_str`  - YAML-formatted string options

Level: intermediate

-seealso: `PetscOptionsSetValue()`, `PetscOptionsView()`, `PetscOptionsHasName()`, `PetscOptionsGetInt()`,
`PetscOptionsGetReal()`, `PetscOptionsGetString()`, `PetscOptionsGetIntArray()`, `PetscOptionsBool()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`, `PetscOptionsInsertFile()`, `PetscOptionsInsertFileYAML()`

# External Links
$(_doc_external("Sys/PetscOptionsInsertStringYAML"))
"""
function PetscOptionsInsertStringYAML(petsclib::PetscLibType, options::PetscOptions, in_str::String) end

@for_petsc function PetscOptionsInsertStringYAML(petsclib::$UnionPetscLib, options::PetscOptions, in_str::String )

    @chk ccall(
               (:PetscOptionsInsertStringYAML, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Cchar}),
               options, in_str,
              )


	return nothing
end 

"""
	PetscOptionsInsertFileYAML(petsclib::PetscLibType,comm::MPI_Comm, options::PetscOptions, file::String, require::PetscBool) 
Insert a YAML

Collective

Input Parameters:
- `comm`    - the processes that will share the options (usually `PETSC_COMM_WORLD`)
- `options` - options database, use `NULL` for default global database
- `file`    - name of file
- `require` - if `PETSC_TRUE` will generate an error if the file does not exist

Level: intermediate

-seealso: `PetscOptionsSetValue()`, `PetscOptionsView()`, `PetscOptionsHasName()`, `PetscOptionsGetInt()`,
`PetscOptionsGetReal()`, `PetscOptionsGetString()`, `PetscOptionsGetIntArray()`, `PetscOptionsBool()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`, `PetscOptionsInsertFile()`, `PetscOptionsInsertStringYAML()`

# External Links
$(_doc_external("Sys/PetscOptionsInsertFileYAML"))
"""
function PetscOptionsInsertFileYAML(petsclib::PetscLibType, comm::MPI_Comm, options::PetscOptions, file::String, require::PetscBool) end

@for_petsc function PetscOptionsInsertFileYAML(petsclib::$UnionPetscLib, comm::MPI_Comm, options::PetscOptions, file::String, require::PetscBool )

    @chk ccall(
               (:PetscOptionsInsertFileYAML, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, COptions, Ptr{Cchar}, PetscBool),
               comm, options, file, require,
              )


	return nothing
end 

"""
	options::PetscOptions = PetscOptionsCreate(petsclib::PetscLibType) 
Creates an empty options database.

Logically Collective

Output Parameter:
- `options` - Options database object

Level: advanced

-seealso: `PetscOptionsDestroy()`, `PetscOptionsPush()`, `PetscOptionsPop()`, `PetscOptionsInsert()`, `PetscOptionsSetValue()`

# External Links
$(_doc_external("Sys/PetscOptionsCreate"))
"""
function PetscOptionsCreate(petsclib::PetscLibType) end

@for_petsc function PetscOptionsCreate(petsclib::$UnionPetscLib)
	options_ = Ref{COptions}()

    @chk ccall(
               (:PetscOptionsCreate, $petsc_library),
               PetscErrorCode,
               (Ptr{COptions},),
               options_,
              )

	options = PetscOptions(options_[], petsclib)

	return options
end 

"""
	PetscOptionsDestroy(petsclib::PetscLibType,options::PetscOptions) 
Destroys an option database.

Logically Collective on whatever communicator was associated with the call to `PetscOptionsCreate()`

Input Parameter:
- `options` - the `PetscOptions` object

Level: advanced

-seealso: `PetscOptionsInsert()`, `PetscOptionsPush()`, `PetscOptionsPop()`, `PetscOptionsSetValue()`

# External Links
$(_doc_external("Sys/PetscOptionsDestroy"))
"""
function PetscOptionsDestroy(petsclib::PetscLibType, options::PetscOptions) end

@for_petsc function PetscOptionsDestroy(petsclib::$UnionPetscLib, options::PetscOptions )
	options_ = Ref(options.ptr)

    @chk ccall(
               (:PetscOptionsDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{COptions},),
               options_,
              )

	options.ptr = C_NULL

	return nothing
end 

"""
	PetscOptionsCreateDefault(petsclib::PetscLibType) 

# External Links
$(_doc_external("Sys/PetscOptionsCreateDefault"))
"""
function PetscOptionsCreateDefault(petsclib::PetscLibType) end

@for_petsc function PetscOptionsCreateDefault(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscOptionsCreateDefault, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscOptionsPush(petsclib::PetscLibType,opt::PetscOptions) 
Push a new `PetscOptions` object as the default provider of options
Allows using different parts of a code to use different options databases

Logically Collective

Input Parameter:
- `opt` - the options obtained with `PetscOptionsCreate()`

Level: advanced

-seealso: `PetscOptionsPop()`, `PetscOptionsCreate()`, `PetscOptionsInsert()`, `PetscOptionsSetValue()`, `PetscOptionsLeft()`

# External Links
$(_doc_external("Sys/PetscOptionsPush"))
"""
function PetscOptionsPush(petsclib::PetscLibType, opt::PetscOptions) end

@for_petsc function PetscOptionsPush(petsclib::$UnionPetscLib, opt::PetscOptions )

    @chk ccall(
               (:PetscOptionsPush, $petsc_library),
               PetscErrorCode,
               (COptions,),
               opt,
              )


	return nothing
end 

"""
	PetscOptionsPop(petsclib::PetscLibType) 
Pop the most recent `PetscOptionsPush()` to return to the previous default options

Logically Collective on whatever communicator was associated with the call to `PetscOptionsCreate()`

Level: advanced

-seealso: `PetscOptionsCreate()`, `PetscOptionsInsert()`, `PetscOptionsSetValue()`, `PetscOptionsLeft()`

# External Links
$(_doc_external("Sys/PetscOptionsPop"))
"""
function PetscOptionsPop(petsclib::PetscLibType) end

@for_petsc function PetscOptionsPop(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscOptionsPop, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscOptionsDestroyDefault(petsclib::PetscLibType) 

# External Links
$(_doc_external("Sys/PetscOptionsDestroyDefault"))
"""
function PetscOptionsDestroyDefault(petsclib::PetscLibType) end

@for_petsc function PetscOptionsDestroyDefault(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscOptionsDestroyDefault, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	valid::PetscBool = PetscOptionsValidKey(petsclib::PetscLibType,key::String) 
PETSc Options database keys must begin with one or two dashes (

Not Collective

Input Parameter:
- `key` - string to check if valid

Output Parameter:
- `valid` - `PETSC_TRUE` if a valid key

Level: intermediate

-seealso: `PetscOptionsCreate()`, `PetscOptionsInsert()`

# External Links
$(_doc_external("Sys/PetscOptionsValidKey"))
"""
function PetscOptionsValidKey(petsclib::PetscLibType, key::String) end

@for_petsc function PetscOptionsValidKey(petsclib::$UnionPetscLib, key::String )
	valid_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsValidKey, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{PetscBool}),
               key, valid_,
              )

	valid = valid_[]

	return valid
end 

"""
	PetscOptionsInsertString(petsclib::PetscLibType,options::PetscOptions, in_str::String) 
Inserts options into the database from a string

Logically Collective

Input Parameters:
- `options` - options object
- `in_str`  - string that contains options separated by blanks

Level: intermediate

The collectivity of this routine is complex; only the MPI processes that call this routine will
have the affect of these options. If some processes that create objects call this routine and others do
not the code may fail in complicated ways because the same parallel solvers may incorrectly use different options
on different ranks.

Contributed by Boyana Norris

-seealso: `PetscOptionsSetValue()`, `PetscOptionsView()`, `PetscOptionsHasName()`, `PetscOptionsGetInt()`,
`PetscOptionsGetReal()`, `PetscOptionsGetString()`, `PetscOptionsGetIntArray()`, `PetscOptionsBool()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`, `PetscOptionsInsertFile()`

# External Links
$(_doc_external("Sys/PetscOptionsInsertString"))
"""
function PetscOptionsInsertString(petsclib::PetscLibType, options::PetscOptions, in_str::String) end

@for_petsc function PetscOptionsInsertString(petsclib::$UnionPetscLib, options::PetscOptions, in_str::String )

    @chk ccall(
               (:PetscOptionsInsertString, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Cchar}),
               options, in_str,
              )


	return nothing
end 

"""
	PetscOptionsInsertFile(petsclib::PetscLibType,comm::MPI_Comm, options::PetscOptions, file::String, require::PetscBool) 
Inserts options into the database from a file.

Collective

Input Parameters:
- `comm`    - the processes that will share the options (usually `PETSC_COMM_WORLD`)
- `options` - options database, use `NULL` for default global database
- `file`    - name of file,
".yml" and ".yaml" filename extensions are inserted as YAML options,
append ":yaml" to filename to force YAML options.
- `require` - if `PETSC_TRUE` will generate an error if the file does not exist

Level: developer

-seealso: `PetscOptionsSetValue()`, `PetscOptionsView()`, `PetscOptionsHasName()`, `PetscOptionsGetInt()`,
`PetscOptionsGetReal()`, `PetscOptionsGetString()`, `PetscOptionsGetIntArray()`, `PetscOptionsBool()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`

# External Links
$(_doc_external("Sys/PetscOptionsInsertFile"))
"""
function PetscOptionsInsertFile(petsclib::PetscLibType, comm::MPI_Comm, options::PetscOptions, file::String, require::PetscBool) end

@for_petsc function PetscOptionsInsertFile(petsclib::$UnionPetscLib, comm::MPI_Comm, options::PetscOptions, file::String, require::PetscBool )

    @chk ccall(
               (:PetscOptionsInsertFile, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, COptions, Ptr{Cchar}, PetscBool),
               comm, options, file, require,
              )


	return nothing
end 

"""
	PetscOptionsInsertArgs(petsclib::PetscLibType,options::PetscOptions, argc::Cint, args::String) 
Inserts options into the database from a array of strings

Logically Collective

Input Parameters:
- `options` - options object
- `argc`    - the array length
- `args`    - the string array

Level: intermediate

-seealso: `PetscOptions`, `PetscOptionsInsertString()`, `PetscOptionsInsertFile()`

# External Links
$(_doc_external("Sys/PetscOptionsInsertArgs"))
"""
function PetscOptionsInsertArgs(petsclib::PetscLibType, options::PetscOptions, argc::Cint, args::String) end

@for_petsc function PetscOptionsInsertArgs(petsclib::$UnionPetscLib, options::PetscOptions, argc::Cint, args::String )
	args_ = Ref(pointer(args))

    @chk ccall(
               (:PetscOptionsInsertArgs, $petsc_library),
               PetscErrorCode,
               (COptions, Cint, Ptr{Ptr{Cchar}}),
               options, argc, args_,
              )


	return nothing
end 

"""
	PetscOptionsInsert(petsclib::PetscLibType,options::PetscOptions, argc::Cint, args::Cchar, file::String) 
Inserts into the options database from the command line,
the environmental variable and a file.

Collective on `PETSC_COMM_WORLD`

Input Parameters:
- `options` - options database or `NULL` for the default global database
- `argc`    - count of number of command line arguments
- `args`    - the command line arguments
- `file`    - [optional] PETSc database file, append ":yaml" to filename to specify YAML options format.
Use `NULL` or empty string to not check for code specific file.
Also checks ~/.petscrc, .petscrc and petscrc.
Use -skip_petscrc in the code specific file (or command line) to skip ~/.petscrc, .petscrc and petscrc files.

Options Database Keys:
- `-options_file <filename>`      - read options from a file
- `-options_file_yaml <filename>` - read options from a YAML file

Level: advanced

-seealso: `PetscOptionsDestroy()`, `PetscOptionsView()`, `PetscOptionsInsertString()`, `PetscOptionsInsertFile()`,
`PetscInitialize()`

# External Links
$(_doc_external("Sys/PetscOptionsInsert"))
"""
function PetscOptionsInsert(petsclib::PetscLibType, options::PetscOptions, argc::Cint, args::Cchar, file::String) end

@for_petsc function PetscOptionsInsert(petsclib::$UnionPetscLib, options::PetscOptions, argc::Cint, args::Cchar, file::String )

    @chk ccall(
               (:PetscOptionsInsert, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Cint}, Cchar, Ptr{Cchar}),
               options, argc, args, file,
              )


	return nothing
end 

"""
	PetscOptionsView(petsclib::PetscLibType,options::PetscOptions, viewer::PetscViewerFormat) 
Prints the options that have been loaded. This is
useful for debugging purposes.

Logically Collective, No Fortran Support

Input Parameters:
- `options` - options database, use `NULL` for default global database
- `viewer`  - must be an `PETSCVIEWERASCII` viewer, can be `PETSC_VIEWER_DEFAULT`

Options Database Key:
- `-options_view` - Activates `PetscOptionsView()` within `PetscFinalize()`

Level: advanced

-seealso: `PetscOptionsAllUsed()`

# External Links
$(_doc_external("Sys/PetscOptionsView"))
"""
function PetscOptionsView(petsclib::PetscLibType, options::PetscOptions, viewer::PetscViewer) end

@for_petsc function PetscOptionsView(petsclib::$UnionPetscLib, options::PetscOptions, viewer::PetscViewerFormat )

    @chk ccall(
               (:PetscOptionsView, $petsc_library),
               PetscErrorCode,
               (COptions, PetscViewerFormat),
               options, viewer,
              )


	return nothing
end 

"""
	PetscOptionsLeftError(petsclib::PetscLibType) 

# External Links
$(_doc_external("Sys/PetscOptionsLeftError"))
"""
function PetscOptionsLeftError(petsclib::PetscLibType) end

@for_petsc function PetscOptionsLeftError(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscOptionsLeftError, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscOptionsPrefixPush(petsclib::PetscLibType,options::PetscOptions, prefix::String) 
Designate a prefix to be used by all options insertions to follow.

Logically Collective

Input Parameters:
- `options` - options database, or `NULL` for the default global database
- `prefix`  - The string to append to the existing prefix

Options Database Keys:
- `-prefix_push <some_prefix_>` - push the given prefix
- `-prefix_pop`                 - pop the last prefix

Level: advanced

-seealso: `PetscOptionsPrefixPop()`, `PetscOptionsPush()`, `PetscOptionsPop()`, `PetscOptionsCreate()`, `PetscOptionsSetValue()`

# External Links
$(_doc_external("Sys/PetscOptionsPrefixPush"))
"""
function PetscOptionsPrefixPush(petsclib::PetscLibType, options::PetscOptions, prefix::String) end

@for_petsc function PetscOptionsPrefixPush(petsclib::$UnionPetscLib, options::PetscOptions, prefix::String )

    @chk ccall(
               (:PetscOptionsPrefixPush, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Cchar}),
               options, prefix,
              )


	return nothing
end 

"""
	PetscOptionsPrefixPop(petsclib::PetscLibType,options::PetscOptions) 
Remove the latest options prefix, see `PetscOptionsPrefixPush()` for details

Logically Collective on the `MPI_Comm` used when called `PetscOptionsPrefixPush()`

Input Parameter:
- `options` - options database, or `NULL` for the default global database

Level: advanced

-seealso: `PetscOptionsPrefixPush()`, `PetscOptionsPush()`, `PetscOptionsPop()`, `PetscOptionsCreate()`, `PetscOptionsSetValue()`

# External Links
$(_doc_external("Sys/PetscOptionsPrefixPop"))
"""
function PetscOptionsPrefixPop(petsclib::PetscLibType, options::PetscOptions) end

@for_petsc function PetscOptionsPrefixPop(petsclib::$UnionPetscLib, options::PetscOptions )

    @chk ccall(
               (:PetscOptionsPrefixPop, $petsc_library),
               PetscErrorCode,
               (COptions,),
               options,
              )


	return nothing
end 

"""
	PetscOptionsClear(petsclib::PetscLibType,options::PetscOptions) 
Removes all options form the database leaving it empty.

Logically Collective

Input Parameter:
- `options` - options database, use `NULL` for the default global database

Level: developer

-seealso: `PetscOptionsInsert()`

# External Links
$(_doc_external("Sys/PetscOptionsClear"))
"""
function PetscOptionsClear(petsclib::PetscLibType, options::PetscOptions) end

@for_petsc function PetscOptionsClear(petsclib::$UnionPetscLib, options::PetscOptions )

    @chk ccall(
               (:PetscOptionsClear, $petsc_library),
               PetscErrorCode,
               (COptions,),
               options,
              )


	return nothing
end 

"""
	PetscOptionsSetAlias(petsclib::PetscLibType,options::PetscOptions, newname::String, oldname::String) 
Makes a key and alias for another key

Logically Collective

Input Parameters:
- `options` - options database, or `NULL` for default global database
- `newname` - the alias
- `oldname` - the name that alias will refer to

Level: advanced

-seealso: `PetscOptionsGetInt()`, `PetscOptionsGetReal()`, `PetscOptionsHasName()`,
`PetscOptionsGetString()`, `PetscOptionsGetIntArray()`, `PetscOptionsGetRealArray()`, `PetscOptionsBool()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`

# External Links
$(_doc_external("Sys/PetscOptionsSetAlias"))
"""
function PetscOptionsSetAlias(petsclib::PetscLibType, options::PetscOptions, newname::String, oldname::String) end

@for_petsc function PetscOptionsSetAlias(petsclib::$UnionPetscLib, options::PetscOptions, newname::String, oldname::String )

    @chk ccall(
               (:PetscOptionsSetAlias, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Cchar}, Ptr{Cchar}),
               options, newname, oldname,
              )


	return nothing
end 

"""
	PetscOptionsSetValue(petsclib::PetscLibType,options::PetscOptions, name::String, value::Union{Ptr,String}) 
Sets an option name
database, overriding whatever is already present.

Logically Collective

Input Parameters:
- `options` - options database, use `NULL` for the default global database
- `name`    - name of option, this SHOULD have the - prepended
- `value`   - the option value (not used for all options, so can be `NULL`)

Level: intermediate

-seealso: `PetscOptionsInsert()`, `PetscOptionsClearValue()`

# External Links
$(_doc_external("Sys/PetscOptionsSetValue"))
"""
function PetscOptionsSetValue(petsclib::PetscLibType, options::PetscOptions, name::String, value::Union{Ptr,String}) end

@for_petsc function PetscOptionsSetValue(petsclib::$UnionPetscLib, options::PetscOptions, name::String, value::Union{Ptr,String} )

    @chk ccall(
               (:PetscOptionsSetValue, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Cchar}, Ptr{Cchar}),
               options, name, value,
              )


	return nothing
end 

"""
	PetscOptionsClearValue(petsclib::PetscLibType,options::PetscOptions, name::String) 
Clears an option name
database, overriding whatever is already present.

Logically Collective

Input Parameters:
- `options` - options database, use `NULL` for the default global database
- `name`    - name of option, this SHOULD have the - prepended

Level: intermediate

-seealso: `PetscOptionsInsert()`

# External Links
$(_doc_external("Sys/PetscOptionsClearValue"))
"""
function PetscOptionsClearValue(petsclib::PetscLibType, options::PetscOptions, name::String) end

@for_petsc function PetscOptionsClearValue(petsclib::$UnionPetscLib, options::PetscOptions, name::String )

    @chk ccall(
               (:PetscOptionsClearValue, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Cchar}),
               options, name,
              )


	return nothing
end 

"""
	set::PetscBool = PetscOptionsFindPair(petsclib::PetscLibType,options::PetscOptions, pre::String, name::String, value::String) 
Gets an option name

Not Collective

Input Parameters:
- `options` - options database, use `NULL` for the default global database
- `pre`     - the string to prepend to the name or `NULL`, this SHOULD NOT have the "-" prepended
- `name`    - name of option, this SHOULD have the "-" prepended

Output Parameters:
- `value` - the option value (optional, not used for all options)
- `set`   - whether the option is set (optional)

Level: developer

-seealso: `PetscOptionsSetValue()`, `PetscOptionsClearValue()`

# External Links
$(_doc_external("Sys/PetscOptionsFindPair"))
"""
function PetscOptionsFindPair(petsclib::PetscLibType, options::PetscOptions, pre::String, name::String, value::String) end

@for_petsc function PetscOptionsFindPair(petsclib::$UnionPetscLib, options::PetscOptions, pre::String, name::String, value::String )
	value_ = Ref(pointer(value))
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsFindPair, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Cchar}, Ptr{Cchar}, Ptr{Ptr{Cchar}}, Ptr{PetscBool}),
               options, pre, name, value_, set_,
              )

	set = set_[]

	return set
end 

"""
	PetscOptionsReject(petsclib::PetscLibType,options::PetscOptions, pre::String, name::String, mess::String) 
Generates an error if a certain option is given.

Not Collective

Input Parameters:
- `options` - options database, use `NULL` for default global database
- `pre`     - the option prefix (may be `NULL`)
- `name`    - the option name one is seeking
- `mess`    - error message (may be `NULL`)

Level: advanced

-seealso: `PetscOptionsGetInt()`, `PetscOptionsGetReal()`, `PetscOptionsHasName()`,
`PetscOptionsGetString()`, `PetscOptionsGetIntArray()`, `PetscOptionsGetRealArray()`, `PetscOptionsBool()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`

# External Links
$(_doc_external("Sys/PetscOptionsReject"))
"""
function PetscOptionsReject(petsclib::PetscLibType, options::PetscOptions, pre::String, name::String, mess::String) end

@for_petsc function PetscOptionsReject(petsclib::$UnionPetscLib, options::PetscOptions, pre::String, name::String, mess::String )

    @chk ccall(
               (:PetscOptionsReject, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}),
               options, pre, name, mess,
              )


	return nothing
end 

"""
	set::PetscBool = PetscOptionsHasHelp(petsclib::PetscLibType,options::PetscOptions) 
Determines whether the "

Not Collective

Input Parameter:
- `options` - options database, use `NULL` for default global database

Output Parameter:
- `set` - `PETSC_TRUE` if found else `PETSC_FALSE`.

Level: advanced

-seealso: `PetscOptionsHasName()`

# External Links
$(_doc_external("Sys/PetscOptionsHasHelp"))
"""
function PetscOptionsHasHelp(petsclib::PetscLibType, options::PetscOptions) end

@for_petsc function PetscOptionsHasHelp(petsclib::$UnionPetscLib, options::PetscOptions )
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsHasHelp, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{PetscBool}),
               options, set_,
              )

	set = set_[]

	return set
end 

"""
	set::PetscBool = PetscOptionsHasName(petsclib::PetscLibType,options::PetscOptions, pre::String, name::String) 
Determines whether a certain option is given in the database. This returns true whether the option is a number, string or Boolean, even
if its value is set to false.

Not Collective

Input Parameters:
- `options` - options database, use `NULL` for default global database
- `pre`     - string to prepend to the name or `NULL`
- `name`    - the option one is seeking

Output Parameter:
- `set` - `PETSC_TRUE` if found else `PETSC_FALSE`.

Level: beginner

-seealso: `PetscOptionsGetInt()`, `PetscOptionsGetReal()`,
`PetscOptionsGetString()`, `PetscOptionsGetIntArray()`, `PetscOptionsGetRealArray()`, `PetscOptionsBool()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`

# External Links
$(_doc_external("Sys/PetscOptionsHasName"))
"""
function PetscOptionsHasName(petsclib::PetscLibType, options::PetscOptions, pre::String, name::String) end

@for_petsc function PetscOptionsHasName(petsclib::$UnionPetscLib, options::PetscOptions, pre::String, name::String )
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsHasName, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscBool}),
               options, pre, name, set_,
              )

	set = set_[]

	return set
end 

"""
	PetscOptionsGetAll(petsclib::PetscLibType,options::PetscOptions, copts::String) 
Lists all the options the program was run with in a single string.

Not Collective

Input Parameter:
- `options` - the options database, use `NULL` for the default global database

Output Parameter:
- `copts` - pointer where string pointer is stored

Level: advanced

-seealso: `PetscOptionsAllUsed()`, `PetscOptionsView()`, `PetscOptionsPush()`, `PetscOptionsPop()`

# External Links
$(_doc_external("Sys/PetscOptionsGetAll"))
"""
function PetscOptionsGetAll(petsclib::PetscLibType, options::PetscOptions, copts::String) end

@for_petsc function PetscOptionsGetAll(petsclib::$UnionPetscLib, options::PetscOptions, copts::String )
	copts_ = Ref(pointer(copts))

    @chk ccall(
               (:PetscOptionsGetAll, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Ptr{Cchar}}),
               options, copts_,
              )


	return nothing
end 

"""
	used::PetscBool = PetscOptionsUsed(petsclib::PetscLibType,options::PetscOptions, name::String) 
Indicates if PETSc has used a particular option set in the database

Not Collective

Input Parameters:
- `options` - options database, use `NULL` for default global database
- `name`    - string name of option

Output Parameter:
- `used` - `PETSC_TRUE` if the option was used, otherwise false, including if option was not found in options database

Level: advanced

-seealso: `PetscOptionsView()`, `PetscOptionsLeft()`, `PetscOptionsAllUsed()`

# External Links
$(_doc_external("Sys/PetscOptionsUsed"))
"""
function PetscOptionsUsed(petsclib::PetscLibType, options::PetscOptions, name::String) end

@for_petsc function PetscOptionsUsed(petsclib::$UnionPetscLib, options::PetscOptions, name::String )
	used_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsUsed, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Cchar}, Ptr{PetscBool}),
               options, name, used_,
              )

	used = used_[]

	return used
end 

"""
	N::PetscInt = PetscOptionsAllUsed(petsclib::PetscLibType,options::PetscOptions) 
Returns a count of the number of options in the
database that have never been selected.

Not Collective

Input Parameter:
- `options` - options database, use `NULL` for default global database

Output Parameter:
- `N` - count of options not used

Level: advanced

-seealso: `PetscOptionsView()`

# External Links
$(_doc_external("Sys/PetscOptionsAllUsed"))
"""
function PetscOptionsAllUsed(petsclib::PetscLibType, options::PetscOptions) end

@for_petsc function PetscOptionsAllUsed(petsclib::$UnionPetscLib, options::PetscOptions )
	N_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscOptionsAllUsed, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{$PetscInt}),
               options, N_,
              )

	N = N_[]

	return N
end 

"""
	PetscOptionsLeft(petsclib::PetscLibType,options::PetscOptions) 
Prints to screen any options that were set and never used.

Not Collective

Input Parameter:
- `options` - options database; use `NULL` for default global database

Options Database Key:
- `-options_left` - activates `PetscOptionsAllUsed()` within `PetscFinalize()`

Level: advanced

-seealso: `PetscOptionsAllUsed()`

# External Links
$(_doc_external("Sys/PetscOptionsLeft"))
"""
function PetscOptionsLeft(petsclib::PetscLibType, options::PetscOptions) end

@for_petsc function PetscOptionsLeft(petsclib::$UnionPetscLib, options::PetscOptions )

    @chk ccall(
               (:PetscOptionsLeft, $petsc_library),
               PetscErrorCode,
               (COptions,),
               options,
              )


	return nothing
end 

"""
	N::PetscInt = PetscOptionsLeftGet(petsclib::PetscLibType,options::PetscOptions, names::String, values::String) 
Returns all options that were set and never used.

Not Collective

Input Parameter:
- `options` - options database, use `NULL` for default global database

Output Parameters:
- `N`      - count of options not used
- `names`  - names of options not used
- `values` - values of options not used

Level: advanced

-seealso: `PetscOptionsAllUsed()`, `PetscOptionsLeft()`

# External Links
$(_doc_external("Sys/PetscOptionsLeftGet"))
"""
function PetscOptionsLeftGet(petsclib::PetscLibType, options::PetscOptions, names::String, values::String) end

@for_petsc function PetscOptionsLeftGet(petsclib::$UnionPetscLib, options::PetscOptions, names::String, values::String )
	N_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscOptionsLeftGet, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{$PetscInt}, Ptr{Cchar}, Ptr{Cchar}),
               options, N_, names, values,
              )

	N = N_[]

	return N
end 

"""
	PetscOptionsLeftRestore(petsclib::PetscLibType,options::PetscOptions, N::PetscInt, names::String, values::String) 
Free memory for the unused PETSc options obtained using `PetscOptionsLeftGet()`.

Not Collective

Input Parameters:
- `options` - options database, use `NULL` for default global database
- `N`       - count of options not used
- `names`   - names of options not used
- `values`  - values of options not used

Level: advanced

-seealso: `PetscOptionsAllUsed()`, `PetscOptionsLeft()`, `PetscOptionsLeftGet()`

# External Links
$(_doc_external("Sys/PetscOptionsLeftRestore"))
"""
function PetscOptionsLeftRestore(petsclib::PetscLibType, options::PetscOptions, N::PetscInt, names::String, values::String) end

@for_petsc function PetscOptionsLeftRestore(petsclib::$UnionPetscLib, options::PetscOptions, N::$PetscInt, names::String, values::String )

    @chk ccall(
               (:PetscOptionsLeftRestore, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{$PetscInt}, Ptr{Cchar}, Ptr{Cchar}),
               options, N, names, values,
              )


	return nothing
end 

"""
	PetscOptionsMonitorDefault(petsclib::PetscLibType,name::String, value::String, source::PetscOptionSource, ctx::Cvoid) 
Print all options set value events using the supplied `PetscViewer`.

Logically Collective

Input Parameters:
- `name`   - option name string
- `value`  - option value string
- `source` - The source for the option
- `ctx`    - a `PETSCVIEWERASCII` or `NULL`

Level: intermediate

-seealso: `PetscOptionsMonitorSet()`

# External Links
$(_doc_external("Sys/PetscOptionsMonitorDefault"))
"""
function PetscOptionsMonitorDefault(petsclib::PetscLibType, name::String, value::String, source::PetscOptionSource, ctx::Cvoid) end

@for_petsc function PetscOptionsMonitorDefault(petsclib::$UnionPetscLib, name::String, value::String, source::PetscOptionSource, ctx::Cvoid )

    @chk ccall(
               (:PetscOptionsMonitorDefault, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cchar}, PetscOptionSource, Ptr{Cvoid}),
               name, value, source, ctx,
              )


	return nothing
end 

"""
	PetscOptionsMonitorSet(petsclib::PetscLibType,monitor::external, mctx::Cvoid, monitordestroy::PetscCtxDestroyFn) 
Sets an ADDITIONAL function to be called at every method that
modified the PETSc options database.

Not Collective

Input Parameters:
- `monitor`        - pointer to function (if this is `NULL`, it turns off monitoring
- `mctx`           - [optional] context for private data for the monitor routine (use `NULL` if
no context is desired)
- `monitordestroy` - [optional] routine that frees monitor context (may be `NULL`), see `PetscCtxDestroyFn` for its calling sequence

Calling sequence of `monitor`:
- `name`   - option name string
- `value`  - option value string, a value of `NULL` indicates the option is being removed from the database. A value
of "" indicates the option is in the database but has no value.
- `source` - option source
- `mctx`   - optional monitoring context, as set by `PetscOptionsMonitorSet()`

Options Database Keys:
- `-options_monitor <viewer>` - turn on default monitoring of changes to the options database
- `-options_monitor_cancel`   - turn off any option monitors except the default monitor obtained with `-options_monitor`

Level: intermediate

-seealso: `PetscOptionsMonitorDefault()`, `PetscInitialize()`, `PetscCtxDestroyFn`

# External Links
$(_doc_external("Sys/PetscOptionsMonitorSet"))
"""
function PetscOptionsMonitorSet(petsclib::PetscLibType, monitor::external, mctx::Cvoid, monitordestroy::PetscCtxDestroyFn) end

@for_petsc function PetscOptionsMonitorSet(petsclib::$UnionPetscLib, monitor::external, mctx::Cvoid, monitordestroy::PetscCtxDestroyFn )

    @chk ccall(
               (:PetscOptionsMonitorSet, $petsc_library),
               PetscErrorCode,
               (external, Ptr{Cvoid}, Ptr{PetscCtxDestroyFn}),
               monitor, mctx, monitordestroy,
              )


	return nothing
end 

"""
	a::PetscBool = PetscOptionsStringToBool(petsclib::PetscLibType,value::String) 

# External Links
$(_doc_external("Sys/PetscOptionsStringToBool"))
"""
function PetscOptionsStringToBool(petsclib::PetscLibType, value::String) end

@for_petsc function PetscOptionsStringToBool(petsclib::$UnionPetscLib, value::String )
	a_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsStringToBool, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{PetscBool}),
               value, a_,
              )

	a = a_[]

	return a
end 

"""
	a::PetscInt = PetscOptionsStringToInt(petsclib::PetscLibType,name::String) 

# External Links
$(_doc_external("Sys/PetscOptionsStringToInt"))
"""
function PetscOptionsStringToInt(petsclib::PetscLibType, name::String) end

@for_petsc function PetscOptionsStringToInt(petsclib::$UnionPetscLib, name::String )
	a_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscOptionsStringToInt, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{$PetscInt}),
               name, a_,
              )

	a = a_[]

	return a
end 

"""
	a::PetscReal = PetscOptionsStringToReal(petsclib::PetscLibType,name::String) 

# External Links
$(_doc_external("Sys/PetscOptionsStringToReal"))
"""
function PetscOptionsStringToReal(petsclib::PetscLibType, name::String) end

@for_petsc function PetscOptionsStringToReal(petsclib::$UnionPetscLib, name::String )
	a_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscOptionsStringToReal, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{$PetscReal}),
               name, a_,
              )

	a = a_[]

	return a
end 

"""
	a::PetscScalar = PetscOptionsStringToScalar(petsclib::PetscLibType,name::String) 

# External Links
$(_doc_external("Sys/PetscOptionsStringToScalar"))
"""
function PetscOptionsStringToScalar(petsclib::PetscLibType, name::String) end

@for_petsc function PetscOptionsStringToScalar(petsclib::$UnionPetscLib, name::String )
	a_ = Ref{$PetscScalar}()

    @chk ccall(
               (:PetscOptionsStringToScalar, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{$PetscScalar}),
               name, a_,
              )

	a = a_[]

	return a
end 

"""
	ivalue::PetscBool,set::PetscBool = PetscOptionsGetBool(petsclib::PetscLibType,options::PetscOptions, pre::String, name::String) 
Gets the Logical (true or false) value for a particular
option in the database.

Not Collective

Input Parameters:
- `options` - options database, use `NULL` for default global database
- `pre`     - the string to prepend to the name or `NULL`
- `name`    - the option one is seeking

Output Parameters:
- `ivalue` - the logical value to return
- `set`    - `PETSC_TRUE`  if found, else `PETSC_FALSE`

Level: beginner

-seealso: `PetscOptionsGetBool3()`, `PetscOptionsGetReal()`, `PetscOptionsHasName()`, `PetscOptionsGetString()`,
`PetscOptionsGetIntArray()`, `PetscOptionsGetRealArray()`, `PetscOptionsGetInt()`, `PetscOptionsBool()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`

# External Links
$(_doc_external("Sys/PetscOptionsGetBool"))
"""
function PetscOptionsGetBool(petsclib::PetscLibType, options::PetscOptions, pre::String, name::String) end

@for_petsc function PetscOptionsGetBool(petsclib::$UnionPetscLib, options::PetscOptions, pre::String, name::String )
	ivalue_ = Ref{PetscBool}()
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsGetBool, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscBool}, Ptr{PetscBool}),
               options, pre, name, ivalue_, set_,
              )

	ivalue = ivalue_[]
	set = set_[]

	return ivalue,set
end 

"""
	set::PetscBool = PetscOptionsGetBool3(petsclib::PetscLibType,options::PetscOptions, pre::String, name::String, ivalue::PetscBool3) 
Gets the ternary logical (true, false or unknown) value for a particular
option in the database.

Not Collective

Input Parameters:
- `options` - options database, use `NULL` for default global database
- `pre`     - the string to prepend to the name or `NULL`
- `name`    - the option one is seeking

Output Parameters:
- `ivalue` - the ternary logical value to return
- `set`    - `PETSC_TRUE`  if found, else `PETSC_FALSE`

Level: beginner

-seealso: `PetscOptionsGetBool()`, `PetscOptionsGetReal()`, `PetscOptionsHasName()`, `PetscOptionsGetString()`,
`PetscOptionsGetIntArray()`, `PetscOptionsGetRealArray()`, `PetscOptionsGetInt()`, `PetscOptionsBool()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`

# External Links
$(_doc_external("Sys/PetscOptionsGetBool3"))
"""
function PetscOptionsGetBool3(petsclib::PetscLibType, options::PetscOptions, pre::String, name::String, ivalue::PetscBool3) end

@for_petsc function PetscOptionsGetBool3(petsclib::$UnionPetscLib, options::PetscOptions, pre::String, name::String, ivalue::PetscBool3 )
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsGetBool3, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscBool3}, Ptr{PetscBool}),
               options, pre, name, ivalue, set_,
              )

	set = set_[]

	return set
end 

"""
	value::PetscInt,set::PetscBool = PetscOptionsGetEList(petsclib::PetscLibType,options::PetscOptions, pre::String, opt::String, list::String, ntext::PetscInt) 
Puts a list of option values that a single one may be selected from

Not Collective

Input Parameters:
- `options` - options database, use `NULL` for default global database
- `pre`     - the string to prepend to the name or `NULL`
- `opt`     - option name
- `list`    - the possible choices (one of these must be selected, anything else is invalid)
- `ntext`   - number of choices

Output Parameters:
- `value` - the index of the value to return (defaults to zero if the option name is given but no choice is listed)
- `set`   - `PETSC_TRUE` if found, else `PETSC_FALSE`

Level: intermediate

-seealso: `PetscOptionsGetInt()`, `PetscOptionsGetReal()`,
`PetscOptionsHasName()`, `PetscOptionsGetIntArray()`, `PetscOptionsGetRealArray()`, `PetscOptionsBool()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`

# External Links
$(_doc_external("Sys/PetscOptionsGetEList"))
"""
function PetscOptionsGetEList(petsclib::PetscLibType, options::PetscOptions, pre::String, opt::String, list::String, ntext::PetscInt) end

@for_petsc function PetscOptionsGetEList(petsclib::$UnionPetscLib, options::PetscOptions, pre::String, opt::String, list::String, ntext::$PetscInt )
	list_ = Ref(pointer(list))
	value_ = Ref{$PetscInt}()
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsGetEList, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Cchar}, Ptr{Cchar}, Ptr{Ptr{Cchar}}, $PetscInt, Ptr{$PetscInt}, Ptr{PetscBool}),
               options, pre, opt, list_, ntext, value_, set_,
              )

	value = value_[]
	set = set_[]

	return value,set
end 

"""
	set::PetscBool = PetscOptionsGetEnum(petsclib::PetscLibType,options::PetscOptions, pre::String, opt::String, list::String, value::PetscEnum) 
Gets the enum value for a particular option in the database.

Not Collective

Input Parameters:
- `options` - options database, use `NULL` for default global database
- `pre`     - option prefix or `NULL`
- `opt`     - option name
- `list`    - array containing the list of choices, followed by the enum name, followed by the enum prefix, followed by a null

Output Parameters:
- `value` - the value to return
- `set`   - `PETSC_TRUE` if found, else `PETSC_FALSE`

Level: beginner

-seealso: `PetscOptionsGetReal()`, `PetscOptionsHasName()`, `PetscOptionsGetString()`, `PetscOptionsGetInt()`,
`PetscOptionsGetIntArray()`, `PetscOptionsGetRealArray()`, `PetscOptionsBool()`
`PetscOptionsInt()`, `PetscOptionsString()`, `PetscOptionsReal()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`, `PetscOptionsGetEList()`, `PetscOptionsEnum()`

# External Links
$(_doc_external("Sys/PetscOptionsGetEnum"))
"""
function PetscOptionsGetEnum(petsclib::PetscLibType, options::PetscOptions, pre::String, opt::String, list::String, value::PetscEnum) end

@for_petsc function PetscOptionsGetEnum(petsclib::$UnionPetscLib, options::PetscOptions, pre::String, opt::String, list::String, value::PetscEnum )
	list_ = Ref(pointer(list))
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsGetEnum, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Cchar}, Ptr{Cchar}, Ptr{Ptr{Cchar}}, Ptr{PetscEnum}, Ptr{PetscBool}),
               options, pre, opt, list_, value, set_,
              )

	set = set_[]

	return set
end 

"""
	ivalue::PetscInt,set::PetscBool = PetscOptionsGetInt(petsclib::PetscLibType,options::PetscOptions, pre::String, name::String) 
Gets the integer value for a particular option in the database.

Not Collective

Input Parameters:
- `options` - options database, use `NULL` for default global database
- `pre`     - the string to prepend to the name or `NULL`
- `name`    - the option one is seeking

Output Parameters:
- `ivalue` - the integer value to return
- `set`    - `PETSC_TRUE` if found, else `PETSC_FALSE`

Level: beginner

-seealso: `PetscOptionsGetReal()`, `PetscOptionsHasName()`, `PetscOptionsGetString()`,
`PetscOptionsGetIntArray()`, `PetscOptionsGetRealArray()`, `PetscOptionsBool()`
`PetscOptionsInt()`, `PetscOptionsString()`, `PetscOptionsReal()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`

# External Links
$(_doc_external("Sys/PetscOptionsGetInt"))
"""
function PetscOptionsGetInt(petsclib::PetscLibType, options::PetscOptions, pre::String, name::String) end

@for_petsc function PetscOptionsGetInt(petsclib::$UnionPetscLib, options::PetscOptions, pre::String, name::String )
	ivalue_ = Ref{$PetscInt}()
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsGetInt, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Cchar}, Ptr{Cchar}, Ptr{$PetscInt}, Ptr{PetscBool}),
               options, pre, name, ivalue_, set_,
              )

	ivalue = ivalue_[]
	set = set_[]

	return ivalue,set
end 

"""
	set::PetscBool = PetscOptionsGetMPIInt(petsclib::PetscLibType,options::PetscOptions, pre::String, name::String, ivalue::PetscMPIInt) 
Gets the MPI integer value for a particular option in the database.

Not Collective

Input Parameters:
- `options` - options database, use `NULL` for default global database
- `pre`     - the string to prepend to the name or `NULL`
- `name`    - the option one is seeking

Output Parameters:
- `ivalue` - the MPI integer value to return
- `set`    - `PETSC_TRUE` if found, else `PETSC_FALSE`

Level: beginner

-seealso: `PetscOptionsGetReal()`, `PetscOptionsHasName()`, `PetscOptionsGetString()`,
`PetscOptionsGetIntArray()`, `PetscOptionsGetRealArray()`, `PetscOptionsBool()`
`PetscOptionsInt()`, `PetscOptionsString()`, `PetscOptionsReal()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`

# External Links
$(_doc_external("Sys/PetscOptionsGetMPIInt"))
"""
function PetscOptionsGetMPIInt(petsclib::PetscLibType, options::PetscOptions, pre::String, name::String, ivalue::PetscMPIInt) end

@for_petsc function PetscOptionsGetMPIInt(petsclib::$UnionPetscLib, options::PetscOptions, pre::String, name::String, ivalue::PetscMPIInt )
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsGetMPIInt, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscMPIInt}, Ptr{PetscBool}),
               options, pre, name, ivalue, set_,
              )

	set = set_[]

	return set
end 

"""
	dvalue::PetscReal,set::PetscBool = PetscOptionsGetReal(petsclib::PetscLibType,options::PetscOptions, pre::String, name::String) 
Gets the double precision value for a particular
option in the database.

Not Collective

Input Parameters:
- `options` - options database, use `NULL` for default global database
- `pre`     - string to prepend to each name or `NULL`
- `name`    - the option one is seeking

Output Parameters:
- `dvalue` - the double value to return
- `set`    - `PETSC_TRUE` if found, `PETSC_FALSE` if not found

Level: beginner

-seealso: `PetscOptionsGetInt()`, `PetscOptionsHasName()`,
`PetscOptionsGetString()`, `PetscOptionsGetIntArray()`, `PetscOptionsGetRealArray()`, `PetscOptionsBool()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`

# External Links
$(_doc_external("Sys/PetscOptionsGetReal"))
"""
function PetscOptionsGetReal(petsclib::PetscLibType, options::PetscOptions, pre::String, name::String) end

@for_petsc function PetscOptionsGetReal(petsclib::$UnionPetscLib, options::PetscOptions, pre::String, name::String )
	dvalue_ = Ref{$PetscReal}()
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsGetReal, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Cchar}, Ptr{Cchar}, Ptr{$PetscReal}, Ptr{PetscBool}),
               options, pre, name, dvalue_, set_,
              )

	dvalue = dvalue_[]
	set = set_[]

	return dvalue,set
end 

"""
	dvalue::PetscScalar,set::PetscBool = PetscOptionsGetScalar(petsclib::PetscLibType,options::PetscOptions, pre::String, name::String) 
Gets the scalar value for a particular
option in the database.

Not Collective

Input Parameters:
- `options` - options database, use `NULL` for default global database
- `pre`     - string to prepend to each name or `NULL`
- `name`    - the option one is seeking

Output Parameters:
- `dvalue` - the scalar value to return
- `set`    - `PETSC_TRUE` if found, else `PETSC_FALSE`

Level: beginner

-seealso: `PetscOptionsGetInt()`, `PetscOptionsHasName()`,
`PetscOptionsGetString()`, `PetscOptionsGetIntArray()`, `PetscOptionsGetRealArray()`, `PetscOptionsBool()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`

# External Links
$(_doc_external("Sys/PetscOptionsGetScalar"))
"""
function PetscOptionsGetScalar(petsclib::PetscLibType, options::PetscOptions, pre::String, name::String) end

@for_petsc function PetscOptionsGetScalar(petsclib::$UnionPetscLib, options::PetscOptions, pre::String, name::String )
	dvalue_ = Ref{$PetscScalar}()
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsGetScalar, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Cchar}, Ptr{Cchar}, Ptr{$PetscScalar}, Ptr{PetscBool}),
               options, pre, name, dvalue_, set_,
              )

	dvalue = dvalue_[]
	set = set_[]

	return dvalue,set
end 

"""
	string::Union{Bool,String} = PetscOptionsGetString(petsclib::PetscLibType,options::PetscOptions, pre::String, name::String, string::String, len::Csize_t) 
Gets the string value for a particular option in
the database.

Not Collective

Input Parameters:
- `options` - options database, use `NULL` for default global database
- `pre`     - string to prepend to name or `NULL`
- `name`    - the option one is seeking
- `len`     - maximum length of the string including null termination

Output Parameters:
- `string` - returns the value of the parameter ifn set, otherwise `false`

Level: beginner

-seealso: `PetscOptionsGetInt()`, `PetscOptionsGetReal()`,
`PetscOptionsHasName()`, `PetscOptionsGetIntArray()`, `PetscOptionsGetRealArray()`, `PetscOptionsBool()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`

# External Links
$(_doc_external("Sys/PetscOptionsGetString"))
"""
function PetscOptionsGetString(petsclib::PetscLibType, options::PetscOptions, pre::String, name::String) end

@for_petsc function PetscOptionsGetString(petsclib::$UnionPetscLib, options::PetscOptions, pre::Union{Ptr,String}, name::String)
	set_ = Ref{PetscBool}()
    val = Vector{UInt8}(undef, 256)

    @chk ccall(
               (:PetscOptionsGetString, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Csize_t, Ptr{PetscBool}),
               options, pre, name, val, sizeof(val), set_,
              )

	set = set_[]
    if set
        val = GC.@preserve val unsafe_string(pointer(val))
    else
        val = false
    end
  

	return val
end 

"""
	dvalue::Vector{PetscBool},nmax::PetscInt,set::PetscBool = PetscOptionsGetBoolArray(petsclib::PetscLibType,options::PetscOptions, pre::String, name::String) 
Gets an array of Logical (true or false) values for a particular
option in the database.  The values must be separated with commas with no intervening spaces.

Not Collective

Input Parameters:
- `options` - options database, use `NULL` for default global database
- `pre`     - string to prepend to each name or `NULL`
- `name`    - the option one is seeking

Output Parameters:
- `dvalue` - the Boolean values to return
- `nmax`   - On input maximum number of values to retrieve, on output the actual number of values retrieved
- `set`    - `PETSC_TRUE` if found, else `PETSC_FALSE`

Level: beginner

-seealso: `PetscOptionsGetInt()`, `PetscOptionsHasName()`,
`PetscOptionsGetString()`, `PetscOptionsGetRealArray()`, `PetscOptionsBool()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`

# External Links
$(_doc_external("Sys/PetscOptionsGetBoolArray"))
"""
function PetscOptionsGetBoolArray(petsclib::PetscLibType, options::PetscOptions, pre::String, name::String) end

@for_petsc function PetscOptionsGetBoolArray(petsclib::$UnionPetscLib, options::PetscOptions, pre::String, name::String )
	dvalue = Vector{PetscBool}(undef, ni);  # CHECK SIZE!!
	nmax_ = Ref{$PetscInt}()
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsGetBoolArray, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscBool}, Ptr{$PetscInt}, Ptr{PetscBool}),
               options, pre, name, dvalue, nmax_, set_,
              )

	nmax = nmax_[]
	set = set_[]

	return dvalue,nmax,set
end 

"""
	nmax::PetscInt,set::PetscBool = PetscOptionsGetEnumArray(petsclib::PetscLibType,options::PetscOptions, pre::String, name::String, list::String, ivalue::Vector{PetscEnum}) 
Gets an array of enum values for a particular option in the database.

Not Collective

Input Parameters:
- `options` - options database, use `NULL` for default global database
- `pre`     - option prefix or `NULL`
- `name`    - option name
- `list`    - array containing the list of choices, followed by the enum name, followed by the enum prefix, followed by a null

Output Parameters:
- `ivalue` - the  enum values to return
- `nmax`   - On input maximum number of values to retrieve, on output the actual number of values retrieved
- `set`    - `PETSC_TRUE` if found, else `PETSC_FALSE`

Level: beginner

-seealso: `PetscOptionsGetReal()`, `PetscOptionsHasName()`, `PetscOptionsGetString()`, `PetscOptionsGetInt()`,
`PetscOptionsGetEnum()`, `PetscOptionsGetIntArray()`, `PetscOptionsGetRealArray()`, `PetscOptionsBool()`
`PetscOptionsInt()`, `PetscOptionsString()`, `PetscOptionsReal()`, `PetscOptionsName()`,
`PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`, `PetscOptionsStringArray()`, `PetscOptionsRealArray()`,
`PetscOptionsScalar()`, `PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`, `PetscOptionsGetEList()`, `PetscOptionsEnum()`

# External Links
$(_doc_external("Sys/PetscOptionsGetEnumArray"))
"""
function PetscOptionsGetEnumArray(petsclib::PetscLibType, options::PetscOptions, pre::String, name::String, list::String, ivalue::Vector{PetscEnum}) end

@for_petsc function PetscOptionsGetEnumArray(petsclib::$UnionPetscLib, options::PetscOptions, pre::String, name::String, list::String, ivalue::Vector{PetscEnum} )
	list_ = Ref(pointer(list))
	nmax_ = Ref{$PetscInt}()
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsGetEnumArray, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Cchar}, Ptr{Cchar}, Ptr{Ptr{Cchar}}, Ptr{PetscEnum}, Ptr{$PetscInt}, Ptr{PetscBool}),
               options, pre, name, list_, ivalue, nmax_, set_,
              )

	nmax = nmax_[]
	set = set_[]

	return nmax,set
end 

"""
	ivalue::Vector{PetscInt},nmax::PetscInt,set::PetscBool = PetscOptionsGetIntArray(petsclib::PetscLibType,options::PetscOptions, pre::String, name::String) 
Gets an array of integer values for a particular option in the database.

Not Collective

Input Parameters:
- `options` - options database, use `NULL` for default global database
- `pre`     - string to prepend to each name or `NULL`
- `name`    - the option one is seeking

Output Parameters:
- `ivalue` - the integer values to return
- `nmax`   - On input maximum number of values to retrieve, on output the actual number of values retrieved
- `set`    - `PETSC_TRUE` if found, else `PETSC_FALSE`

Level: beginner

-seealso: `PetscOptionsGetInt()`, `PetscOptionsHasName()`,
`PetscOptionsGetString()`, `PetscOptionsGetRealArray()`, `PetscOptionsBool()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`

# External Links
$(_doc_external("Sys/PetscOptionsGetIntArray"))
"""
function PetscOptionsGetIntArray(petsclib::PetscLibType, options::PetscOptions, pre::String, name::String) end

@for_petsc function PetscOptionsGetIntArray(petsclib::$UnionPetscLib, options::PetscOptions, pre::String, name::String )
	ivalue = Vector{$PetscInt}(undef, ni);  # CHECK SIZE!!
	nmax_ = Ref{$PetscInt}()
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsGetIntArray, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Cchar}, Ptr{Cchar}, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{PetscBool}),
               options, pre, name, ivalue, nmax_, set_,
              )

	nmax = nmax_[]
	set = set_[]

	return ivalue,nmax,set
end 

"""
	dvalue::Vector{PetscReal},nmax::PetscInt,set::PetscBool = PetscOptionsGetRealArray(petsclib::PetscLibType,options::PetscOptions, pre::String, name::String) 
Gets an array of double precision values for a
particular option in the database.  The values must be separated with commas with no intervening spaces.

Not Collective

Input Parameters:
- `options` - options database, use `NULL` for default global database
- `pre`     - string to prepend to each name or `NULL`
- `name`    - the option one is seeking

Output Parameters:
- `dvalue` - the double values to return
- `nmax`   - On input maximum number of values to retrieve, on output the actual number of values retrieved
- `set`    - `PETSC_TRUE` if found, else `PETSC_FALSE`

Level: beginner

-seealso: `PetscOptionsGetInt()`, `PetscOptionsHasName()`,
`PetscOptionsGetString()`, `PetscOptionsGetIntArray()`, `PetscOptionsBool()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`

# External Links
$(_doc_external("Sys/PetscOptionsGetRealArray"))
"""
function PetscOptionsGetRealArray(petsclib::PetscLibType, options::PetscOptions, pre::String, name::String) end

@for_petsc function PetscOptionsGetRealArray(petsclib::$UnionPetscLib, options::PetscOptions, pre::String, name::String )
	dvalue = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!
	nmax_ = Ref{$PetscInt}()
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsGetRealArray, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Cchar}, Ptr{Cchar}, Ptr{$PetscReal}, Ptr{$PetscInt}, Ptr{PetscBool}),
               options, pre, name, dvalue, nmax_, set_,
              )

	nmax = nmax_[]
	set = set_[]

	return dvalue,nmax,set
end 

"""
	dvalue::Vector{PetscScalar},nmax::PetscInt,set::PetscBool = PetscOptionsGetScalarArray(petsclib::PetscLibType,options::PetscOptions, pre::String, name::String) 
Gets an array of scalars for a
particular option in the database.  The values must be separated with commas with no intervening spaces.

Not Collective

Input Parameters:
- `options` - options database, use `NULL` for default global database
- `pre`     - string to prepend to each name or `NULL`
- `name`    - the option one is seeking

Output Parameters:
- `dvalue` - the scalar values to return
- `nmax`   - On input maximum number of values to retrieve, on output the actual number of values retrieved
- `set`    - `PETSC_TRUE` if found, else `PETSC_FALSE`

Level: beginner

-seealso: `PetscOptionsGetInt()`, `PetscOptionsHasName()`,
`PetscOptionsGetString()`, `PetscOptionsGetIntArray()`, `PetscOptionsBool()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`

# External Links
$(_doc_external("Sys/PetscOptionsGetScalarArray"))
"""
function PetscOptionsGetScalarArray(petsclib::PetscLibType, options::PetscOptions, pre::String, name::String) end

@for_petsc function PetscOptionsGetScalarArray(petsclib::$UnionPetscLib, options::PetscOptions, pre::String, name::String )
	dvalue = Vector{$PetscScalar}(undef, ni);  # CHECK SIZE!!
	nmax_ = Ref{$PetscInt}()
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsGetScalarArray, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Cchar}, Ptr{Cchar}, Ptr{$PetscScalar}, Ptr{$PetscInt}, Ptr{PetscBool}),
               options, pre, name, dvalue, nmax_, set_,
              )

	nmax = nmax_[]
	set = set_[]

	return dvalue,nmax,set
end 

"""
	nmax::PetscInt,set::PetscBool = PetscOptionsGetStringArray(petsclib::PetscLibType,options::PetscOptions, pre::String, name::String, strings::String) 
Gets an array of string values for a particular
option in the database. The values must be separated with commas with no intervening spaces.

Not Collective; No Fortran Support

Input Parameters:
- `options` - options database, use `NULL` for default global database
- `pre`     - string to prepend to name or `NULL`
- `name`    - the option one is seeking

Output Parameters:
- `strings` - location to copy strings
- `nmax`    - On input maximum number of strings, on output the actual number of strings found
- `set`     - `PETSC_TRUE` if found, else `PETSC_FALSE`

Level: beginner

-seealso: `PetscOptionsGetInt()`, `PetscOptionsGetReal()`,
`PetscOptionsHasName()`, `PetscOptionsGetIntArray()`, `PetscOptionsGetRealArray()`, `PetscOptionsBool()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`

# External Links
$(_doc_external("Sys/PetscOptionsGetStringArray"))
"""
function PetscOptionsGetStringArray(petsclib::PetscLibType, options::PetscOptions, pre::String, name::String, strings::String) end

@for_petsc function PetscOptionsGetStringArray(petsclib::$UnionPetscLib, options::PetscOptions, pre::String, name::String, strings::String )
	strings_ = Ref(pointer(strings))
	nmax_ = Ref{$PetscInt}()
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsGetStringArray, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Cchar}, Ptr{Cchar}, Ptr{Ptr{Cchar}}, Ptr{$PetscInt}, Ptr{PetscBool}),
               options, pre, name, strings_, nmax_, set_,
              )

	nmax = nmax_[]
	set = set_[]

	return nmax,set
end 

"""
	flag::PetscBool = PetscOptionsGetenv(petsclib::PetscLibType,comm::MPI_Comm, name::String, env::String, len::Csize_t) 
Gets an environmental variable, broadcasts to all
processors in communicator from MPI rank zero

Collective

Input Parameters:
- `comm` - communicator to share variable
- `name` - name of environmental variable
- `len`  - amount of space allocated to hold variable

Output Parameters:
- `flag` - if not `NULL` indicates if the variable was found
- `env`  - value of variable

Level: advanced

-seealso: `PetscOptionsHasName()`

# External Links
$(_doc_external("Sys/PetscOptionsGetenv"))
"""
function PetscOptionsGetenv(petsclib::PetscLibType, comm::MPI_Comm, name::String, env::String, len::Csize_t) end

@for_petsc function PetscOptionsGetenv(petsclib::$UnionPetscLib, comm::MPI_Comm, name::String, env::String, len::Csize_t )
	flag_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsGetenv, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Csize_t, Ptr{PetscBool}),
               comm, name, env, len, flag_,
              )

	flag = flag_[]

	return flag
end 

"""
	set::PetscBool = PetscOptionsGetVec(petsclib::PetscLibType,options::PetscOptions, prefix::String, key::String, v::PetscVec) 

# External Links
$(_doc_external("Vec/PetscOptionsGetVec"))
"""
function PetscOptionsGetVec(petsclib::PetscLibType, options::PetscOptions, prefix::String, key::String, v::PetscVec) end

@for_petsc function PetscOptionsGetVec(petsclib::$UnionPetscLib, options::PetscOptions, prefix::String, key::String, v::PetscVec )
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOptionsGetVec, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Cchar}, Ptr{Cchar}, CVec, Ptr{PetscBool}),
               options, prefix, key, v, set_,
              )

	set = set_[]

	return set
end 

