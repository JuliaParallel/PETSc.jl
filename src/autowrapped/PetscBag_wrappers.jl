# autodefined type arguments for class ------
mutable struct _n_PetscBag end
const PetscBag = Ptr{_n_PetscBag}

# -------------------------------------------------------
# autodefined type arguments for class ------
# -------------------------------------------------------
"""
	PetscBagRegisterEnum(petsclib::PetscLibType,bag::PetscBag, addr::Cvoid, list::Cchar, mdefault::PetscEnum, name::String, help::String) 
add an enum value to a `PetscBag`

Logically Collective

Input Parameters:
- `bag`      - the bag of values
- `addr`     - location of enum in struct, for example `&params->dt`
- `list`     - array of strings containing names of enum values followed by enum name followed by enum prefix
- `mdefault` - the initial value, cast with (`PetscEnum`)
- `name`     - the name of the item
- `help`     - longer string with more information about the value

Level: beginner

-seealso: `PetscBag`, `PetscBagSetName()`, `PetscBagView()`, `PetscBagLoad()`, `PetscBagGetData()`
`PetscBagRegisterInt()`, `PetscBagRegisterBool()`, `PetscBagRegisterScalar()`
`PetscBagSetFromOptions()`, `PetscBagCreate()`, `PetscBagGetName()`

# External Links
$(_doc_external("Sys/PetscBagRegisterEnum"))
"""
function PetscBagRegisterEnum(petsclib::PetscLibType, bag::PetscBag, addr::Cvoid, list::Cchar, mdefault::PetscEnum, name::String, help::String) end

@for_petsc function PetscBagRegisterEnum(petsclib::$UnionPetscLib, bag::PetscBag, addr::Cvoid, list::Cchar, mdefault::PetscEnum, name::String, help::String )

    @chk ccall(
               (:PetscBagRegisterEnum, $petsc_library),
               PetscErrorCode,
               (PetscBag, Ptr{Cvoid}, Cchar, PetscEnum, Ptr{Cchar}, Ptr{Cchar}),
               bag, addr, list, mdefault, name, help,
              )


	return nothing
end 

"""
	PetscBagRegisterIntArray(petsclib::PetscLibType,bag::PetscBag, addr::Cvoid, msize::PetscInt, name::String, help::String) 
add a `PetscInt` array to a `PetscBag`

Logically Collective

Input Parameters:
- `bag`   - the bag of values
- `addr`  - location of integer in struct, for example `&params->i`
- `msize` - number of entries in array
- `name`  - name of the array
- `help`  - longer string with more information about the value

Level: beginner

-seealso: `PetscBag`, `PetscBagSetName()`, `PetscBagView()`, `PetscBagLoad()`, `PetscBagGetData()`
`PetscBagRegisterInt()`, `PetscBagRegisterBool()`, `PetscBagRegisterScalar()`
`PetscBagSetFromOptions()`, `PetscBagCreate()`, `PetscBagGetName()`, `PetscBagRegisterEnum()`

# External Links
$(_doc_external("Sys/PetscBagRegisterIntArray"))
"""
function PetscBagRegisterIntArray(petsclib::PetscLibType, bag::PetscBag, addr::Cvoid, msize::PetscInt, name::String, help::String) end

@for_petsc function PetscBagRegisterIntArray(petsclib::$UnionPetscLib, bag::PetscBag, addr::Cvoid, msize::$PetscInt, name::String, help::String )

    @chk ccall(
               (:PetscBagRegisterIntArray, $petsc_library),
               PetscErrorCode,
               (PetscBag, Ptr{Cvoid}, $PetscInt, Ptr{Cchar}, Ptr{Cchar}),
               bag, addr, msize, name, help,
              )


	return nothing
end 

"""
	PetscBagRegisterRealArray(petsclib::PetscLibType,bag::PetscBag, addr::Cvoid, msize::PetscInt, name::String, help::String) 
add a `PetscReal` array to a `PetscBag`

Logically Collective

Input Parameters:
- `bag`   - the bag of values
- `addr`  - location of real array in struct, for example `&params->d`
- `msize` - number of entries in the array
- `name`  - name of the array
- `help`  - longer string with more information about the value

Level: beginner

-seealso: `PetscBag`, `PetscBagSetName()`, `PetscBagView()`, `PetscBagLoad()`, `PetscBagGetData()`
`PetscBagRegisterInt()`, `PetscBagRegisterBool()`, `PetscBagRegisterScalar()`
`PetscBagSetFromOptions()`, `PetscBagCreate()`, `PetscBagGetName()`, `PetscBagRegisterEnum()`

# External Links
$(_doc_external("Sys/PetscBagRegisterRealArray"))
"""
function PetscBagRegisterRealArray(petsclib::PetscLibType, bag::PetscBag, addr::Cvoid, msize::PetscInt, name::String, help::String) end

@for_petsc function PetscBagRegisterRealArray(petsclib::$UnionPetscLib, bag::PetscBag, addr::Cvoid, msize::$PetscInt, name::String, help::String )

    @chk ccall(
               (:PetscBagRegisterRealArray, $petsc_library),
               PetscErrorCode,
               (PetscBag, Ptr{Cvoid}, $PetscInt, Ptr{Cchar}, Ptr{Cchar}),
               bag, addr, msize, name, help,
              )


	return nothing
end 

"""
	PetscBagRegisterInt(petsclib::PetscLibType,bag::PetscBag, addr::Cvoid, mdefault::PetscInt, name::String, help::String) 
add a `PetscInt` value to a `PetscBag`

Logically Collective

Input Parameters:
- `bag`      - the bag of values
- `addr`     - location of integer in struct, for example `&params->i`
- `mdefault` - the initial value
- `name`     - name of the integer
- `help`     - longer string with more information about the value

Level: beginner

-seealso: `PetscBag`, `PetscBagSetName()`, `PetscBagView()`, `PetscBagLoad()`, `PetscBagGetData()`
`PetscBagRegisterInt64()`, `PetscBagRegisterBool()`, `PetscBagRegisterScalar()`
`PetscBagSetFromOptions()`, `PetscBagCreate()`, `PetscBagGetName()`, `PetscBagRegisterEnum()`

# External Links
$(_doc_external("Sys/PetscBagRegisterInt"))
"""
function PetscBagRegisterInt(petsclib::PetscLibType, bag::PetscBag, addr::Cvoid, mdefault::PetscInt, name::String, help::String) end

@for_petsc function PetscBagRegisterInt(petsclib::$UnionPetscLib, bag::PetscBag, addr::Cvoid, mdefault::$PetscInt, name::String, help::String )

    @chk ccall(
               (:PetscBagRegisterInt, $petsc_library),
               PetscErrorCode,
               (PetscBag, Ptr{Cvoid}, $PetscInt, Ptr{Cchar}, Ptr{Cchar}),
               bag, addr, mdefault, name, help,
              )


	return nothing
end 

"""
	PetscBagRegisterInt64(petsclib::PetscLibType,bag::PetscBag, addr::Cvoid, mdefault::PetscInt64, name::String, help::String) 
add a `PetscInt64` value to a `PetscBag`

Logically Collective

Input Parameters:
- `bag`      - the bag of values
- `addr`     - location of integer in struct, for example `&params->i`
- `mdefault` - the initial value
- `name`     - name of the integer
- `help`     - longer string with more information about the value

Level: beginner

-seealso: `PetscBag`, `PetscBagSetName()`, `PetscBagView()`, `PetscBagLoad()`, `PetscBagGetData()`
`PetscBagRegisterInt()`, `PetscBagRegisterBool()`, `PetscBagRegisterScalar()`
`PetscBagSetFromOptions()`, `PetscBagCreate()`, `PetscBagGetName()`, `PetscBagRegisterEnum()`

# External Links
$(_doc_external("Sys/PetscBagRegisterInt64"))
"""
function PetscBagRegisterInt64(petsclib::PetscLibType, bag::PetscBag, addr::Cvoid, mdefault::PetscInt64, name::String, help::String) end

@for_petsc function PetscBagRegisterInt64(petsclib::$UnionPetscLib, bag::PetscBag, addr::Cvoid, mdefault::$PetscInt64, name::String, help::String )

    @chk ccall(
               (:PetscBagRegisterInt64, $petsc_library),
               PetscErrorCode,
               (PetscBag, Ptr{Cvoid}, $PetscInt64, Ptr{Cchar}, Ptr{Cchar}),
               bag, addr, mdefault, name, help,
              )


	return nothing
end 

"""
	PetscBagRegisterBoolArray(petsclib::PetscLibType,bag::PetscBag, addr::Cvoid, msize::PetscInt, name::String, help::String) 
add a n `PetscBool` values to a `PetscBag`

Logically Collective

Input Parameters:
- `bag`   - the bag of values
- `addr`  - location of boolean array in struct, for example `&params->b`
- `msize` - number of entries in array
- `name`  - name of the boolean array
- `help`  - longer string with more information about the value

Level: beginner

-seealso: `PetscBag`, `PetscBagSetName()`, `PetscBagView()`, `PetscBagLoad()`, `PetscBagGetData()`
`PetscBagRegisterInt()`, `PetscBagRegisterBool()`, `PetscBagRegisterScalar()`
`PetscBagSetFromOptions()`, `PetscBagCreate()`, `PetscBagGetName()`, `PetscBagRegisterEnum()`

# External Links
$(_doc_external("Sys/PetscBagRegisterBoolArray"))
"""
function PetscBagRegisterBoolArray(petsclib::PetscLibType, bag::PetscBag, addr::Cvoid, msize::PetscInt, name::String, help::String) end

@for_petsc function PetscBagRegisterBoolArray(petsclib::$UnionPetscLib, bag::PetscBag, addr::Cvoid, msize::$PetscInt, name::String, help::String )

    @chk ccall(
               (:PetscBagRegisterBoolArray, $petsc_library),
               PetscErrorCode,
               (PetscBag, Ptr{Cvoid}, $PetscInt, Ptr{Cchar}, Ptr{Cchar}),
               bag, addr, msize, name, help,
              )


	return nothing
end 

"""
	PetscBagRegisterString(petsclib::PetscLibType,bag::PetscBag, addr::Cvoid, msize::PetscInt, mdefault::String, name::String, help::String) 
add a string value to a `PetscBag`

Logically Collective

Input Parameters:
- `bag`      - the bag of values
- `addr`     - location of start of string in struct, for example `&params->mystring`
- `msize`    - length of the string space in the struct
- `mdefault` - the initial value
- `name`     - name of the string
- `help`     - longer string with more information about the value

Level: beginner

-seealso: `PetscBag`, `PetscBagSetName()`, `PetscBagView()`, `PetscBagLoad()`, `PetscBagGetData()`
`PetscBagRegisterInt()`, `PetscBagRegisterBool()`, `PetscBagRegisterScalar()`
`PetscBagSetFromOptions()`, `PetscBagCreate()`, `PetscBagGetName()`, `PetscBagRegisterEnum()`

# External Links
$(_doc_external("Sys/PetscBagRegisterString"))
"""
function PetscBagRegisterString(petsclib::PetscLibType, bag::PetscBag, addr::Cvoid, msize::PetscInt, mdefault::String, name::String, help::String) end

@for_petsc function PetscBagRegisterString(petsclib::$UnionPetscLib, bag::PetscBag, addr::Cvoid, msize::$PetscInt, mdefault::String, name::String, help::String )

    @chk ccall(
               (:PetscBagRegisterString, $petsc_library),
               PetscErrorCode,
               (PetscBag, Ptr{Cvoid}, $PetscInt, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}),
               bag, addr, msize, mdefault, name, help,
              )


	return nothing
end 

"""
	PetscBagRegisterReal(petsclib::PetscLibType,bag::PetscBag, addr::Cvoid, mdefault::PetscReal, name::String, help::String) 
add a `PetscReal` value to a `PetscBag`

Logically Collective

Input Parameters:
- `bag`      - the bag of values
- `addr`     - location of `PetscReal` in struct, for example `&params->r`
- `mdefault` - the initial value
- `name`     - name of the variable
- `help`     - longer string with more information about the value

Level: beginner

-seealso: `PetscBag`, `PetscBagSetName()`, `PetscBagView()`, `PetscBagLoad()`, `PetscBagGetData()`
`PetscBagRegisterInt()`, `PetscBagRegisterBool()`, `PetscBagRegisterScalar()`
`PetscBagSetFromOptions()`, `PetscBagCreate()`, `PetscBagGetName()`, `PetscBagRegisterEnum()`

# External Links
$(_doc_external("Sys/PetscBagRegisterReal"))
"""
function PetscBagRegisterReal(petsclib::PetscLibType, bag::PetscBag, addr::Cvoid, mdefault::PetscReal, name::String, help::String) end

@for_petsc function PetscBagRegisterReal(petsclib::$UnionPetscLib, bag::PetscBag, addr::Cvoid, mdefault::$PetscReal, name::String, help::String )

    @chk ccall(
               (:PetscBagRegisterReal, $petsc_library),
               PetscErrorCode,
               (PetscBag, Ptr{Cvoid}, $PetscReal, Ptr{Cchar}, Ptr{Cchar}),
               bag, addr, mdefault, name, help,
              )


	return nothing
end 

"""
	PetscBagRegisterScalar(petsclib::PetscLibType,bag::PetscBag, addr::Cvoid, mdefault::PetscScalar, name::String, help::String) 
add a `PetscScalar` value to a `PetscBag`

Logically Collective

Input Parameters:
- `bag`      - the bag of values
- `addr`     - location of `PetscScalar` in struct, for example `&params->c`
- `mdefault` - the initial value
- `name`     - name of the variable
- `help`     - longer string with more information about the value

Level: beginner

-seealso: `PetscBag`, `PetscBagSetName()`, `PetscBagView()`, `PetscBagLoad()`, `PetscBagGetData()`
`PetscBagRegisterInt()`, `PetscBagRegisterBool()`, `PetscBagSetFromOptions()`,
`PetscBagCreate()`, `PetscBagGetName()`, `PetscBagRegisterEnum()`

# External Links
$(_doc_external("Sys/PetscBagRegisterScalar"))
"""
function PetscBagRegisterScalar(petsclib::PetscLibType, bag::PetscBag, addr::Cvoid, mdefault::PetscScalar, name::String, help::String) end

@for_petsc function PetscBagRegisterScalar(petsclib::$UnionPetscLib, bag::PetscBag, addr::Cvoid, mdefault::$PetscScalar, name::String, help::String )

    @chk ccall(
               (:PetscBagRegisterScalar, $petsc_library),
               PetscErrorCode,
               (PetscBag, Ptr{Cvoid}, $PetscScalar, Ptr{Cchar}, Ptr{Cchar}),
               bag, addr, mdefault, name, help,
              )


	return nothing
end 

"""
	PetscBagRegisterBool(petsclib::PetscLibType,bag::PetscBag, addr::Cvoid, mdefault::PetscBool, name::String, help::String) 
add a `PetscBool` to a `PetscBag`

Logically Collective

Input Parameters:
- `bag`      - the bag of values
- `addr`     - location of `PetscBool` in struct, for example `&params->b`
- `mdefault` - the initial value, either `PETSC_FALSE` or `PETSC_TRUE`
- `name`     - name of the variable
- `help`     - longer string with more information about the value

Level: beginner

-seealso: `PetscBag`, `PetscBagSetName()`, `PetscBagView()`, `PetscBagLoad()`, `PetscBagGetData()`
`PetscBagRegisterInt()`, `PetscBagRegisterScalar()`
`PetscBagSetFromOptions()`, `PetscBagCreate()`, `PetscBagGetName()`, `PetscBagRegisterEnum()`

# External Links
$(_doc_external("Sys/PetscBagRegisterBool"))
"""
function PetscBagRegisterBool(petsclib::PetscLibType, bag::PetscBag, addr::Cvoid, mdefault::PetscBool, name::String, help::String) end

@for_petsc function PetscBagRegisterBool(petsclib::$UnionPetscLib, bag::PetscBag, addr::Cvoid, mdefault::PetscBool, name::String, help::String )

    @chk ccall(
               (:PetscBagRegisterBool, $petsc_library),
               PetscErrorCode,
               (PetscBag, Ptr{Cvoid}, PetscBool, Ptr{Cchar}, Ptr{Cchar}),
               bag, addr, mdefault, name, help,
              )


	return nothing
end 

"""
	PetscBagDestroy(petsclib::PetscLibType,bag::PetscBag) 
Destroys a `PetscBag`

Collective

Input Parameter:
- `bag` - the bag of values

Level: beginner

-seealso: `PetscBag`, `PetscBagSetName()`, `PetscBagView()`, `PetscBagLoad()`, `PetscBagGetData()`
`PetscBagRegisterReal()`, `PetscBagRegisterInt()`, `PetscBagRegisterBool()`, `PetscBagRegisterScalar()`
`PetscBagSetFromOptions()`, `PetscBagCreate()`, `PetscBagGetName()`, `PetscBagRegisterEnum()`

# External Links
$(_doc_external("Sys/PetscBagDestroy"))
"""
function PetscBagDestroy(petsclib::PetscLibType, bag::PetscBag) end

@for_petsc function PetscBagDestroy(petsclib::$UnionPetscLib, bag::PetscBag )

    @chk ccall(
               (:PetscBagDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscBag},),
               bag,
              )


	return nothing
end 

"""
	PetscBagSetFromOptions(petsclib::PetscLibType,bag::PetscBag) 
Allows setting entries to a `PetscBag` using the options database

Collective

Input Parameter:
- `bag` - the bag of values

Level: beginner

-seealso: `PetscBag`, `PetscBagSetName()`, `PetscBagDestroy()`, `PetscBagLoad()`, `PetscBagGetData()`
`PetscBagRegisterReal()`, `PetscBagRegisterInt()`, `PetscBagRegisterBool()`, `PetscBagRegisterScalar()`
`PetscBagCreate()`, `PetscBagGetName()`, `PetscBagView()`, `PetscBagRegisterEnum()`

# External Links
$(_doc_external("Sys/PetscBagSetFromOptions"))
"""
function PetscBagSetFromOptions(petsclib::PetscLibType, bag::PetscBag) end

@for_petsc function PetscBagSetFromOptions(petsclib::$UnionPetscLib, bag::PetscBag )

    @chk ccall(
               (:PetscBagSetFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscBag,),
               bag,
              )


	return nothing
end 

"""
	PetscBagView(petsclib::PetscLibType,bag::PetscBag, view::PetscViewer) 
Views a bag of values as either ASCII text or a binary file

Collective

Input Parameters:
- `bag`  - the bag of values
- `view` - location to view the values

Level: beginner

-seealso: `PetscBag`, `PetscBagSetName()`, `PetscBagDestroy()`, `PetscBagLoad()`, `PetscBagGetData()`
`PetscBagRegisterReal()`, `PetscBagRegisterInt()`, `PetscBagRegisterBool()`, `PetscBagRegisterScalar()`, `PetscBagRegisterEnum()`
`PetscBagSetFromOptions()`, `PetscBagCreate()`, `PetscBagGetName()`

# External Links
$(_doc_external("Sys/PetscBagView"))
"""
function PetscBagView(petsclib::PetscLibType, bag::PetscBag, view::PetscViewer) end

@for_petsc function PetscBagView(petsclib::$UnionPetscLib, bag::PetscBag, view::PetscViewer )

    @chk ccall(
               (:PetscBagView, $petsc_library),
               PetscErrorCode,
               (PetscBag, PetscViewer),
               bag, view,
              )


	return nothing
end 

"""
	PetscBagViewFromOptions(petsclib::PetscLibType,bag::PetscBag, bobj::PetscObject, optionname::String) 
Processes command line options to determine if/how a `PetscBag` is to be viewed.

Collective

Input Parameters:
- `bag`        - the object
- `bobj`       - optional other object that provides prefix (if `NULL` then the prefix in obj is used)
- `optionname` - option to activate viewing

Level: intermediate

-seealso: `PetscBagCreate()`, `PetscBag`, `PetscViewer`

# External Links
$(_doc_external("Sys/PetscBagViewFromOptions"))
"""
function PetscBagViewFromOptions(petsclib::PetscLibType, bag::PetscBag, bobj::PetscObject, optionname::String) end

@for_petsc function PetscBagViewFromOptions(petsclib::$UnionPetscLib, bag::PetscBag, bobj::PetscObject, optionname::String )

    @chk ccall(
               (:PetscBagViewFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscBag, PetscObject, Ptr{Cchar}),
               bag, bobj, optionname,
              )


	return nothing
end 

"""
	PetscBagLoad(petsclib::PetscLibType,view::PetscViewer, bag::PetscBag) 
Loads a bag of values from a binary file

Collective

Input Parameters:
- `view` - file to load values from
- `bag`  - the bag of values

Level: beginner

-seealso: `PetscBag`, `PetscBagSetName()`, `PetscBagDestroy()`, `PetscBagView()`, `PetscBagGetData()`
`PetscBagRegisterReal()`, `PetscBagRegisterInt()`, `PetscBagRegisterBool()`, `PetscBagRegisterScalar()`
`PetscBagSetFromOptions()`, `PetscBagCreate()`, `PetscBagGetName()`, `PetscBagRegisterEnum()`

# External Links
$(_doc_external("Sys/PetscBagLoad"))
"""
function PetscBagLoad(petsclib::PetscLibType, view::PetscViewer, bag::PetscBag) end

@for_petsc function PetscBagLoad(petsclib::$UnionPetscLib, view::PetscViewer, bag::PetscBag )

    @chk ccall(
               (:PetscBagLoad, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscBag),
               view, bag,
              )


	return nothing
end 

"""
	bag::PetscBag = PetscBagCreate(petsclib::PetscLibType,comm::MPI_Comm, bagsize::Csize_t) 
Create a bag of values. A `PetscBag` is a representation of a C struct that can be saved to and read from files,
can have values set from the options database

Collective

Input Parameters:
- `comm`    - communicator to share bag
- `bagsize` - size of the C structure holding the values, for example `sizeof(mystruct)`

Output Parameter:
- `bag` - the bag of values

Level: intermediate

-seealso: `PetscBag`, `PetscBagGetName()`, `PetscBagView()`, `PetscBagLoad()`, `PetscBagGetData()`
`PetscBagRegisterReal()`, `PetscBagRegisterInt()`, `PetscBagRegisterBool()`, `PetscBagRegisterScalar()`
`PetscBagSetFromOptions()`, `PetscBagDestroy()`, `PetscBagRegisterEnum()`

# External Links
$(_doc_external("Sys/PetscBagCreate"))
"""
function PetscBagCreate(petsclib::PetscLibType, comm::MPI_Comm, bagsize::Csize_t) end

@for_petsc function PetscBagCreate(petsclib::$UnionPetscLib, comm::MPI_Comm, bagsize::Csize_t )
	bag_ = Ref{PetscBag}()

    @chk ccall(
               (:PetscBagCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Csize_t, Ptr{PetscBag}),
               comm, bagsize, bag_,
              )

	bag = bag_[]

	return bag
end 

"""
	PetscBagSetName(petsclib::PetscLibType,bag::PetscBag, name::String, help::String) 
Sets the name of a bag of values

Not Collective

Level: intermediate

Input Parameters:
- `bag`  - the bag of values
- `name` - the name assigned to the bag
- `help` - help message for bag

-seealso: `PetscBag`, `PetscBagGetName()`, `PetscBagView()`, `PetscBagLoad()`, `PetscBagGetData()`
`PetscBagRegisterReal()`, `PetscBagRegisterInt()`, `PetscBagRegisterBool()`, `PetscBagRegisterScalar()`
`PetscBagSetFromOptions()`, `PetscBagCreate()`, `PetscBagDestroy()`, `PetscBagRegisterEnum()`

# External Links
$(_doc_external("Sys/PetscBagSetName"))
"""
function PetscBagSetName(petsclib::PetscLibType, bag::PetscBag, name::String, help::String) end

@for_petsc function PetscBagSetName(petsclib::$UnionPetscLib, bag::PetscBag, name::String, help::String )

    @chk ccall(
               (:PetscBagSetName, $petsc_library),
               PetscErrorCode,
               (PetscBag, Ptr{Cchar}, Ptr{Cchar}),
               bag, name, help,
              )


	return nothing
end 

"""
	PetscBagGetName(petsclib::PetscLibType,bag::PetscBag, name::Cchar) 
Gets the name of a bag of values

Not Collective

Level: intermediate

Input Parameter:
- `bag` - the bag of values

Output Parameter:
- `name` - the name assigned to the bag

-seealso: `PetscBag`, `PetscBagSetName()`, `PetscBagView()`, `PetscBagLoad()`, `PetscBagGetData()`
`PetscBagRegisterReal()`, `PetscBagRegisterInt()`, `PetscBagRegisterBool()`, `PetscBagRegisterScalar()`
`PetscBagSetFromOptions()`, `PetscBagCreate()`, `PetscBagDestroy()`, `PetscBagRegisterEnum()`

# External Links
$(_doc_external("Sys/PetscBagGetName"))
"""
function PetscBagGetName(petsclib::PetscLibType, bag::PetscBag, name::Cchar) end

@for_petsc function PetscBagGetName(petsclib::$UnionPetscLib, bag::PetscBag, name::Cchar )

    @chk ccall(
               (:PetscBagGetName, $petsc_library),
               PetscErrorCode,
               (PetscBag, Cchar),
               bag, name,
              )


	return nothing
end 

"""
	PetscBagGetData(petsclib::PetscLibType,bag::PetscBag, data::PeCtx) 
Gives back the user
can be used for storing user-data-structure

Not Collective

Input Parameter:
- `bag` - the bag of values

Output Parameter:
- `data` - pointer to memory that will have user-data-structure, this can be cast to a pointer of the type the C struct used in
defining the bag

Level: intermediate

-seealso: `PetscBag`, `PetscBagSetName()`, `PetscBagView()`, `PetscBagLoad()`
`PetscBagRegisterReal()`, `PetscBagRegisterInt()`, `PetscBagRegisterBool()`, `PetscBagRegisterScalar()`
`PetscBagSetFromOptions()`, `PetscBagCreate()`, `PetscBagDestroy()`, `PetscBagRegisterEnum()`

# External Links
$(_doc_external("Sys/PetscBagGetData"))
"""
function PetscBagGetData(petsclib::PetscLibType, bag::PetscBag, data::PeCtx) end

@for_petsc function PetscBagGetData(petsclib::$UnionPetscLib, bag::PetscBag, data::PeCtx )

    @chk ccall(
               (:PetscBagGetData, $petsc_library),
               PetscErrorCode,
               (PetscBag, PeCtx),
               bag, data,
              )


	return nothing
end 

"""
	PetscBagSetOptionsPrefix(petsclib::PetscLibType,bag::PetscBag, pre::String) 
Sets the prefix used for searching for all
`PetscBag` items in the options database.

Logically Collective

Level: intermediate

Input Parameters:
- `bag` - the bag of values
- `pre` - the prefix to prepend all Bag item names with.

-seealso: `PetscBag`, `PetscBagRegisterReal()`, `PetscBagRegisterInt()`, `PetscBagRegisterBool()`, `PetscBagRegisterScalar()`
`PetscBagSetFromOptions()`, `PetscBagCreate()`, `PetscBagDestroy()`, `PetscBagRegisterEnum()`

# External Links
$(_doc_external("Sys/PetscBagSetOptionsPrefix"))
"""
function PetscBagSetOptionsPrefix(petsclib::PetscLibType, bag::PetscBag, pre::String) end

@for_petsc function PetscBagSetOptionsPrefix(petsclib::$UnionPetscLib, bag::PetscBag, pre::String )

    @chk ccall(
               (:PetscBagSetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (PetscBag, Ptr{Cchar}),
               bag, pre,
              )


	return nothing
end 

"""
	PetscBagGetNames(petsclib::PetscLibType,bag::PetscBag, names::String) 
Get the names of all entries in the bag

Not Collective

Input Parameter:
- `bag` - the bag of values

Output Parameter:
- `names` - pass in an array of char pointers to hold the names. The array must be as long as the number of items in the bag.

Level: intermediate

-seealso: `PetscBag`, `PetscBagGetName()`, `PetscBagSetName()`, `PetscBagCreate()`, `PetscBagGetData()`
`PetscBagRegisterReal()`, `PetscBagRegisterInt()`, `PetscBagRegisterBool()`, `PetscBagRegisterScalar()`, `PetscBagRegisterEnum()`

# External Links
$(_doc_external("Sys/PetscBagGetNames"))
"""
function PetscBagGetNames(petsclib::PetscLibType, bag::PetscBag, names::String) end

@for_petsc function PetscBagGetNames(petsclib::$UnionPetscLib, bag::PetscBag, names::String )
	names_ = Ref(pointer(names))

    @chk ccall(
               (:PetscBagGetNames, $petsc_library),
               PetscErrorCode,
               (PetscBag, Ptr{Ptr{Cchar}}),
               bag, names_,
              )


	return nothing
end 

