# autodefined type arguments for class ------
mutable struct _n_PetscBench end
const PetscBench = Ptr{_n_PetscBench}

# -------------------------------------------------------
"""
	PetscBenchInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `PetscBench` package.

Level: developer

-seealso: `PetscInitialize()`, `PetscBenchCreate()`, `PetscBench`, `PetscBenchType`

# External Links
$(_doc_external("Sys/PetscBenchInitializePackage"))
"""
function PetscBenchInitializePackage(petsclib::PetscLibType) end

@for_petsc function PetscBenchInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscBenchInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscBenchRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Adds a benchmark test, `PetscBenchType`, to the `PetscBench` package

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - name of a new benchmark
- `function` - routine to create benchmark

Level: advanced

-seealso: `PetscBenchInitializePackage()`, `PetscBenchCreate()`, `PetscBench`, `PetscBenchType`, `PetscBenchSetType()`, `PetscBenchGetType()`

# External Links
$(_doc_external("Sys/PetscBenchRegister"))
"""
function PetscBenchRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function PetscBenchRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:PetscBenchRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	PetscBenchReset(petsclib::PetscLibType,bm::PetscBench) 
removes all the intermediate data structures in a `PetscBench`

Collective

Input Parameter:
- `bm` - the `PetscBench`

Level: advanced

-seealso: `PetscBench`, `PetscBenchView()`, `PetscBenchSetFromOptions()`, `PetscBenchCreate()`, `PetscBenchDestroy()`, `PetscBenchSetUp()`, `PetscBenchSetType()`

# External Links
$(_doc_external("Sys/PetscBenchReset"))
"""
function PetscBenchReset(petsclib::PetscLibType, bm::PetscBench) end

@for_petsc function PetscBenchReset(petsclib::$UnionPetscLib, bm::PetscBench )

    @chk ccall(
               (:PetscBenchReset, $petsc_library),
               PetscErrorCode,
               (PetscBench,),
               bm,
              )


	return nothing
end 

"""
	PetscBenchDestroy(petsclib::PetscLibType,bm::PetscBench) 
Destroys a `PetscBench`

Collective

Input Parameter:
- `bm` - the `PetscBench`

Level: advanced

-seealso: `PetscBench`, `PetscBenchView()`, `PetscBenchSetFromOptions()`, `PetscBenchCreate()`

# External Links
$(_doc_external("Sys/PetscBenchDestroy"))
"""
function PetscBenchDestroy(petsclib::PetscLibType, bm::PetscBench) end

@for_petsc function PetscBenchDestroy(petsclib::$UnionPetscLib, bm::PetscBench )

    @chk ccall(
               (:PetscBenchDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscBench},),
               bm,
              )


	return nothing
end 

"""
	PetscBenchSetUp(petsclib::PetscLibType,bm::PetscBench) 
sets up the `PetscBench`

Collective

Input Parameter:
- `bm` - the `PetscBench`

Level: advanced

-seealso: `PetscBench`, `PetscBenchView()`, `PetscBenchSetFromOptions()`, `PetscBenchCreate()`, `PetscBenchDestroy()`, `PetscBenchSetType()`,
`PetscBenchRun()`, `PetscBenchSetSize()`, `PetscBenchGetSize()`

# External Links
$(_doc_external("Sys/PetscBenchSetUp"))
"""
function PetscBenchSetUp(petsclib::PetscLibType, bm::PetscBench) end

@for_petsc function PetscBenchSetUp(petsclib::$UnionPetscLib, bm::PetscBench )

    @chk ccall(
               (:PetscBenchSetUp, $petsc_library),
               PetscErrorCode,
               (PetscBench,),
               bm,
              )


	return nothing
end 

"""
	PetscBenchRun(petsclib::PetscLibType,bm::PetscBench) 
runs the `PetscBench`

Collective

Input Parameter:
- `bm` - the `PetscBench`

Level: advanced

-seealso: `PetscBench`, `PetscBenchView()`, `PetscBenchSetFromOptions()`, `PetscBenchCreate()`, `PetscBenchDestroy()`, `PetscBenchSetUp()`, `PetscBenchSetType()`,
`PetscBenchSetSize()`, `PetscBenchGetSize()`

# External Links
$(_doc_external("Sys/PetscBenchRun"))
"""
function PetscBenchRun(petsclib::PetscLibType, bm::PetscBench) end

@for_petsc function PetscBenchRun(petsclib::$UnionPetscLib, bm::PetscBench )

    @chk ccall(
               (:PetscBenchRun, $petsc_library),
               PetscErrorCode,
               (PetscBench,),
               bm,
              )


	return nothing
end 

"""
	PetscBenchSetFromOptions(petsclib::PetscLibType,bm::PetscBench) 
Sets options to a `PetscBench` using the options database

Collective

Input Parameter:
- `bm` - the `PetscBench`

Level: advanced

-seealso: `PetscBench`, `PetscBenchView()`, `PetscBenchRun()`, `PetscBenchCreate()`, `PetscBenchDestroy()`, `PetscBenchSetUp()`, `PetscBenchSetType()`,
`PetscBenchSetSize()`, `PetscBenchGetSize()`

# External Links
$(_doc_external("Sys/PetscBenchSetFromOptions"))
"""
function PetscBenchSetFromOptions(petsclib::PetscLibType, bm::PetscBench) end

@for_petsc function PetscBenchSetFromOptions(petsclib::$UnionPetscLib, bm::PetscBench )

    @chk ccall(
               (:PetscBenchSetFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscBench,),
               bm,
              )


	return nothing
end 

"""
	PetscBenchView(petsclib::PetscLibType,bm::PetscBench, viewer::PetscViewer) 
Views a PETSc benchmark `PetscBench`

Collective

Input Parameters:
- `bm`     - the `PetscBench`
- `viewer` - location to view `bm`

Level: advanced

-seealso: `PetscBench`, `PetscBenchSetFromOptions()`, `PetscBenchRun()`, `PetscBenchCreate()`, `PetscBenchDestroy()`, `PetscBenchSetUp()`, `PetscBenchSetType()`,
`PetscBenchSetSize()`, `PetscBenchGetSize()`, `PetscBenchViewFromOptions()`

# External Links
$(_doc_external("Sys/PetscBenchView"))
"""
function PetscBenchView(petsclib::PetscLibType, bm::PetscBench, viewer::PetscViewer) end

@for_petsc function PetscBenchView(petsclib::$UnionPetscLib, bm::PetscBench, viewer::PetscViewer )

    @chk ccall(
               (:PetscBenchView, $petsc_library),
               PetscErrorCode,
               (PetscBench, PetscViewer),
               bm, viewer,
              )


	return nothing
end 

"""
	PetscBenchViewFromOptions(petsclib::PetscLibType,bm::PetscBench, bobj::PetscObject, optionname::String) 
Processes command line options to determine if/how a `PetscBench` is to be viewed.

Collective

Input Parameters:
- `bm`         - the object
- `bobj`       - optional other object that provides prefix (if `NULL` then the prefix in `bm` is used)
- `optionname` - option to activate viewing

Level: advanced

-seealso: `PetscBench`, `PetscBenchSetFromOptions()`, `PetscBenchRun()`, `PetscBenchCreate()`, `PetscBenchDestroy()`, `PetscBenchSetUp()`, `PetscBenchSetType()`,
`PetscBenchSetSize()`, `PetscBenchGetSize()`

# External Links
$(_doc_external("Sys/PetscBenchViewFromOptions"))
"""
function PetscBenchViewFromOptions(petsclib::PetscLibType, bm::PetscBench, bobj::PetscObject, optionname::String) end

@for_petsc function PetscBenchViewFromOptions(petsclib::$UnionPetscLib, bm::PetscBench, bobj::PetscObject, optionname::String )

    @chk ccall(
               (:PetscBenchViewFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscBench, PetscObject, Ptr{Cchar}),
               bm, bobj, optionname,
              )


	return nothing
end 

"""
	bm::PetscBench = PetscBenchCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Create a PETSc benchmark `PetscBench` object

Collective

Input Parameter:
- `comm` - communicator to share the `PetscBench`

Output Parameter:
- `bm` - the `PetscBench`

Level: advanced

-seealso: `PetscBench`, `PetscBenchSetFromOptions()`, `PetscBenchRun()`, `PetscBenchViewFromOptions()`, `PetscBenchDestroy()`, `PetscBenchSetUp()`, `PetscBenchSetType()`,
`PetscBenchSetSize()`, `PetscBenchGetSize()`

# External Links
$(_doc_external("Sys/PetscBenchCreate"))
"""
function PetscBenchCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function PetscBenchCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	bm_ = Ref{PetscBench}()

    @chk ccall(
               (:PetscBenchCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscBench}),
               comm, bm_,
              )

	bm = bm_[]

	return bm
end 

"""
	PetscBenchSetOptionsPrefix(petsclib::PetscLibType,bm::PetscBench, pre::String) 
Sets the prefix used for searching for all `PetscBench` items in the options database.

Logically Collective

Input Parameters:
- `bm`  - the `PetscBench`
- `pre` - the prefix to prepend all `PetscBench` option names

Level: advanced

-seealso: `PetscBench`, `PetscBenchSetFromOptions()`, `PetscBenchRun()`, `PetscBenchViewFromOptions()`, `PetscBenchDestroy()`, `PetscBenchSetUp()`, `PetscBenchSetType()`,
`PetscBenchSetSize()`, `PetscBenchGetSize()`

# External Links
$(_doc_external("Sys/PetscBenchSetOptionsPrefix"))
"""
function PetscBenchSetOptionsPrefix(petsclib::PetscLibType, bm::PetscBench, pre::String) end

@for_petsc function PetscBenchSetOptionsPrefix(petsclib::$UnionPetscLib, bm::PetscBench, pre::String )

    @chk ccall(
               (:PetscBenchSetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (PetscBench, Ptr{Cchar}),
               bm, pre,
              )


	return nothing
end 

"""
	PetscBenchSetSize(petsclib::PetscLibType,bm::PetscBench, n::PetscInt) 
Sets the size of the `PetscBench` benchmark to run

Logically Collective

Input Parameters:
- `bm` - the `PetscBench`
- `n`  - the size

Level: advanced

-seealso: `PetscBench`, `PetscBenchSetFromOptions()`, `PetscBenchRun()`, `PetscBenchViewFromOptions()`, `PetscBenchDestroy()`, `PetscBenchSetUp()`, `PetscBenchSetType()`,
`PetscBenchSetOptionsPrefix()`, `PetscBenchGetSize()`

# External Links
$(_doc_external("Sys/PetscBenchSetSize"))
"""
function PetscBenchSetSize(petsclib::PetscLibType, bm::PetscBench, n::PetscInt) end

@for_petsc function PetscBenchSetSize(petsclib::$UnionPetscLib, bm::PetscBench, n::$PetscInt )

    @chk ccall(
               (:PetscBenchSetSize, $petsc_library),
               PetscErrorCode,
               (PetscBench, $PetscInt),
               bm, n,
              )


	return nothing
end 

"""
	n::PetscInt = PetscBenchGetSize(petsclib::PetscLibType,bm::PetscBench) 
Gets the size of the `PetscBench` benchmark to run

Logically Collective

Input Parameter:
- `bm` - the `PetscBench`

Output Parameter:
- `n` - the size

Level: advanced

-seealso: `PetscBench`, `PetscBenchSetFromOptions()`, `PetscBenchRun()`, `PetscBenchViewFromOptions()`, `PetscBenchDestroy()`, `PetscBenchSetUp()`, `PetscBenchSetType()`,
`PetscBenchSetOptionsPrefix()`, `PetscBenchSetSize()`

# External Links
$(_doc_external("Sys/PetscBenchGetSize"))
"""
function PetscBenchGetSize(petsclib::PetscLibType, bm::PetscBench) end

@for_petsc function PetscBenchGetSize(petsclib::$UnionPetscLib, bm::PetscBench )
	n_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscBenchGetSize, $petsc_library),
               PetscErrorCode,
               (PetscBench, Ptr{$PetscInt}),
               bm, n_,
              )

	n = n_[]

	return n
end 

"""
	PetscBenchSetType(petsclib::PetscLibType,bm::PetscBench, type::PetscBenchType) 
set the type of `PetscBench` benchmark to run

Collective

Input Parameters:
- `bm`   - the `PetscBench`
- `type` - a known method

Options Database Key:
- `-bm_type <type>` - Sets `PetscBench` type

Level: advanced

-seealso: `PetscBench`, `PetscBenchSetFromOptions()`, `PetscBenchRun()`, `PetscBenchViewFromOptions()`, `PetscBenchDestroy()`, `PetscBenchSetUp()`, `PetscBenchGetSize()`,
`PetscBenchSetOptionsPrefix()`, `PetscBenchSetSize()`, `PetscBenchGetType()`, `PetscBenchCreate()`

# External Links
$(_doc_external("Sys/PetscBenchSetType"))
"""
function PetscBenchSetType(petsclib::PetscLibType, bm::PetscBench, type::PetscBenchType) end

@for_petsc function PetscBenchSetType(petsclib::$UnionPetscLib, bm::PetscBench, type::PetscBenchType )

    @chk ccall(
               (:PetscBenchSetType, $petsc_library),
               PetscErrorCode,
               (PetscBench, PetscBenchType),
               bm, type,
              )


	return nothing
end 

"""
	type::PetscBenchType = PetscBenchGetType(petsclib::PetscLibType,bm::PetscBench) 
Gets the `PetscBenchType` (as a string) from the `PetscBench`
context.

Not Collective

Input Parameter:
- `bm` - the `PetscBench`

Output Parameter:
- `type` - name of benchmark method

Level: intermediate

-seealso: `PetscBench`, `PetscBenchType`, `PetscBenchSetType()`, `PetscBenchCreate()`

# External Links
$(_doc_external("Sys/PetscBenchGetType"))
"""
function PetscBenchGetType(petsclib::PetscLibType, bm::PetscBench) end

@for_petsc function PetscBenchGetType(petsclib::$UnionPetscLib, bm::PetscBench )
	type_ = Ref{PetscBenchType}()

    @chk ccall(
               (:PetscBenchGetType, $petsc_library),
               PetscErrorCode,
               (PetscBench, Ptr{PetscBenchType}),
               bm, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

