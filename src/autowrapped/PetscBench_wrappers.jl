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
$(_doc_external("BM/PetscBenchCreate"))
"""
function PetscBenchCreate(petsclib::PetscLibType, comm::MPI_Comm)
    error("PetscBenchCreate: no generated method for these argument types")
end

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
	PetscBenchDestroy(petsclib::PetscLibType,bm::Union{PetscBench, Ref{PetscBench}}) 
Destroys a `PetscBench`

Collective

Input Parameter:
- `bm` - the `PetscBench`

Level: advanced

-seealso: `PetscBench`, `PetscBenchView()`, `PetscBenchSetFromOptions()`, `PetscBenchCreate()`

# External Links
$(_doc_external("BM/PetscBenchDestroy"))
"""
function PetscBenchDestroy(petsclib::PetscLibType, bm::Union{PetscBench, Ref{PetscBench}})
    error("PetscBenchDestroy: no generated method for these argument types")
end

@for_petsc function PetscBenchDestroy(petsclib::$UnionPetscLib, bm::Union{PetscBench, Ref{PetscBench}} )
	bm_ = bm isa Base.RefValue ? bm : Ref{PetscBench}(bm)

    @chk ccall(
               (:PetscBenchDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscBench},),
               bm_,
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
$(_doc_external("BM/PetscBenchGetSize"))
"""
function PetscBenchGetSize(petsclib::PetscLibType, bm::PetscBench)
    error("PetscBenchGetSize: no generated method for these argument types")
end

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
$(_doc_external("BM/PetscBenchGetType"))
"""
function PetscBenchGetType(petsclib::PetscLibType, bm::PetscBench)
    error("PetscBenchGetType: no generated method for these argument types")
end

@for_petsc function PetscBenchGetType(petsclib::$UnionPetscLib, bm::PetscBench )
	type_ = Ref{PetscBenchType}()

    @chk ccall(
               (:PetscBenchGetType, $petsc_library),
               PetscErrorCode,
               (PetscBench, Ptr{PetscBenchType}),
               bm, type_,
              )

	type = type_[] == C_NULL ? "" : unsafe_string(type_[])

	return type
end 

"""
	PetscBenchInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `PetscBench` package.

Level: developer

-seealso: `PetscInitialize()`, `PetscBenchCreate()`, `PetscBench`, `PetscBenchType`

# External Links
$(_doc_external("BM/PetscBenchInitializePackage"))
"""
function PetscBenchInitializePackage(petsclib::PetscLibType)
    error("PetscBenchInitializePackage: no generated method for these argument types")
end

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

Calling sequence of function:
- `bm` - the `PetscBench` to be created

Level: advanced

-seealso: `PetscBenchInitializePackage()`, `PetscBenchCreate()`, `PetscBench`, `PetscBenchType`, `PetscBenchSetType()`, `PetscBenchGetType()`

# External Links
$(_doc_external("BM/PetscBenchRegister"))
"""
function PetscBenchRegister(petsclib::PetscLibType, sname::String, fnc::external)
    error("PetscBenchRegister: no generated method for these argument types")
end

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
$(_doc_external("BM/PetscBenchReset"))
"""
function PetscBenchReset(petsclib::PetscLibType, bm::PetscBench)
    error("PetscBenchReset: no generated method for these argument types")
end

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
	PetscBenchRun(petsclib::PetscLibType,bm::PetscBench) 
runs the `PetscBench`

Collective

Input Parameter:
- `bm` - the `PetscBench`

Level: advanced

-seealso: `PetscBench`, `PetscBenchView()`, `PetscBenchSetFromOptions()`, `PetscBenchCreate()`, `PetscBenchDestroy()`, `PetscBenchSetUp()`, `PetscBenchSetType()`,
`PetscBenchSetSize()`, `PetscBenchGetSize()`

# External Links
$(_doc_external("BM/PetscBenchRun"))
"""
function PetscBenchRun(petsclib::PetscLibType, bm::PetscBench)
    error("PetscBenchRun: no generated method for these argument types")
end

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
$(_doc_external("BM/PetscBenchSetFromOptions"))
"""
function PetscBenchSetFromOptions(petsclib::PetscLibType, bm::PetscBench)
    error("PetscBenchSetFromOptions: no generated method for these argument types")
end

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
$(_doc_external("BM/PetscBenchSetOptionsPrefix"))
"""
function PetscBenchSetOptionsPrefix(petsclib::PetscLibType, bm::PetscBench, pre::String)
    error("PetscBenchSetOptionsPrefix: no generated method for these argument types")
end

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
$(_doc_external("BM/PetscBenchSetSize"))
"""
function PetscBenchSetSize(petsclib::PetscLibType, bm::PetscBench, n::Integer)
    error("PetscBenchSetSize: no generated method for these argument types")
end

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
	PetscBenchSetType(petsclib::PetscLibType,bm::PetscBench, type::PetscBenchType) 
set the type of `PetscBench` benchmark to run

Collective

Input Parameters:
- `bm`   - the `PetscBench`
- `type` - a known method

Options Database Key:
- `-bm_type type` - Sets `PetscBench` type

Level: advanced

-seealso: `PetscBench`, `PetscBenchSetFromOptions()`, `PetscBenchRun()`, `PetscBenchViewFromOptions()`, `PetscBenchDestroy()`, `PetscBenchSetUp()`, `PetscBenchGetSize()`,
`PetscBenchSetOptionsPrefix()`, `PetscBenchSetSize()`, `PetscBenchGetType()`, `PetscBenchCreate()`

# External Links
$(_doc_external("BM/PetscBenchSetType"))
"""
function PetscBenchSetType(petsclib::PetscLibType, bm::PetscBench, type::PetscBenchType)
    error("PetscBenchSetType: no generated method for these argument types")
end

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
	PetscBenchSetUp(petsclib::PetscLibType,bm::PetscBench) 
sets up the `PetscBench`

Collective

Input Parameter:
- `bm` - the `PetscBench`

Level: advanced

-seealso: `PetscBench`, `PetscBenchView()`, `PetscBenchSetFromOptions()`, `PetscBenchCreate()`, `PetscBenchDestroy()`, `PetscBenchSetType()`,
`PetscBenchRun()`, `PetscBenchSetSize()`, `PetscBenchGetSize()`

# External Links
$(_doc_external("BM/PetscBenchSetUp"))
"""
function PetscBenchSetUp(petsclib::PetscLibType, bm::PetscBench)
    error("PetscBenchSetUp: no generated method for these argument types")
end

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
$(_doc_external("BM/PetscBenchView"))
"""
function PetscBenchView(petsclib::PetscLibType, bm::PetscBench, viewer::PetscViewer)
    error("PetscBenchView: no generated method for these argument types")
end

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
	PetscBenchViewFromOptions(petsclib::PetscLibType,bm::PetscBench, bobj, name::String) 
Processes command line options to determine if/how a `PetscBench` is to be viewed.

Collective

Input Parameters:
- `bm`   - the object
- `bobj` - optional other object that provides prefix (if `NULL` then the prefix in `bm` is used)
- `name` - option to activate viewing

Options Database Key:
- `-name [viewertype][:...]` - option name and values. See `PetscObjectViewFromOptions()` for the possible arguments

Level: advanced

-seealso: `PetscBench`, `PetscBenchSetFromOptions()`, `PetscBenchRun()`, `PetscBenchCreate()`, `PetscBenchDestroy()`, `PetscBenchSetUp()`, `PetscBenchSetType()`,
`PetscBenchSetSize()`, `PetscBenchGetSize()`, `PetscObjectViewFromOptions()`, `PetscViewer`, `PetscBenchView()`

# External Links
$(_doc_external("BM/PetscBenchViewFromOptions"))
"""
function PetscBenchViewFromOptions(petsclib::PetscLibType, bm::PetscBench, bobj, name::String)
    error("PetscBenchViewFromOptions: no generated method for these argument types")
end

@for_petsc function PetscBenchViewFromOptions(petsclib::$UnionPetscLib, bm::PetscBench, bobj, name::String )

    @chk ccall(
               (:PetscBenchViewFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscBench, PetscObject, Ptr{Cchar}),
               bm, bobj, name,
              )


	return nothing
end 

