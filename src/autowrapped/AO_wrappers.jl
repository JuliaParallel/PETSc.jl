# autodefined type arguments for class ------
# -------------------------------------------------------
"""
	AOFinalizePackage(petsclib::PetscLibType) 
This function finalizes everything in the `AO` package. It is called
from `PetscFinalize()`.

Level: developer

-seealso: `AOInitializePackage()`, `PetscInitialize()`

# External Links
$(_doc_external("Vec/AOFinalizePackage"))
"""
function AOFinalizePackage(petsclib::PetscLibType) end

@for_petsc function AOFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:AOFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	AOInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `AO` package. It is called
from `PetscDLLibraryRegister_petscvec()` when using dynamic libraries, and on the first call to `AOCreate()`
when using static or shared libraries.

Level: developer

-seealso: `AOFinalizePackage()`, `PetscInitialize()`

# External Links
$(_doc_external("Vec/AOInitializePackage"))
"""
function AOInitializePackage(petsclib::PetscLibType) end

@for_petsc function AOInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:AOInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	AOSetType(petsclib::PetscLibType,ao::AO, method::AOType) 
Builds an application ordering for a particular `AOType`

Collective

Input Parameters:
- `ao`     - The `AO` object
- `method` - The name of the AO type

Options Database Key:
- `-ao_type <type>` - Sets the `AO` type; use -help for a list of available types

Level: intermediate

-seealso: `AO`, `AOType`, `AOCreateBasic()`, `AOCreateMemoryScalable()`, `AOGetType()`, `AOCreate()`

# External Links
$(_doc_external("Vec/AOSetType"))
"""
function AOSetType(petsclib::PetscLibType, ao::AO, method::AOType) end

@for_petsc function AOSetType(petsclib::$UnionPetscLib, ao::AO, method::AOType )

    @chk ccall(
               (:AOSetType, $petsc_library),
               PetscErrorCode,
               (CAO, AOType),
               ao, method,
              )


	return nothing
end 

"""
	type::AOType = AOGetType(petsclib::PetscLibType,ao::AO) 
Gets the `AO` type name (as a string) from the AO.

Not Collective

Input Parameter:
- `ao` - The vector

Output Parameter:
- `type` - The `AO` type name

Level: intermediate

-seealso: `AO`, `AOType`, `AOSetType()`, `AOCreate()`

# External Links
$(_doc_external("Vec/AOGetType"))
"""
function AOGetType(petsclib::PetscLibType, ao::AO) end

@for_petsc function AOGetType(petsclib::$UnionPetscLib, ao::AO )
	type_ = Ref{AOType}()

    @chk ccall(
               (:AOGetType, $petsc_library),
               PetscErrorCode,
               (CAO, Ptr{AOType}),
               ao, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	AORegister(petsclib::PetscLibType,sname::String, fnc::external) 
Register  an application ordering method

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - the name (`AOType`) of the `AO` scheme
- `function` - the create routine for the application ordering method

Level: advanced

-seealso: `AO`, `AOType`, `AOCreate()`, `AORegisterAll()`, `AOBASIC`, `AOADVANCED`, `AOMAPPING`, `AOMEMORYSCALABLE`

# External Links
$(_doc_external("Vec/AORegister"))
"""
function AORegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function AORegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:AORegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	AORegisterAll(petsclib::PetscLibType) 
Registers all of the application ordering components in the `AO` package.

Not Collective

Level: advanced

-seealso: `AO`, `AOType`, `AORegister()`, `AORegisterDestroy()`

# External Links
$(_doc_external("Vec/AORegisterAll"))
"""
function AORegisterAll(petsclib::PetscLibType) end

@for_petsc function AORegisterAll(petsclib::$UnionPetscLib)

    @chk ccall(
               (:AORegisterAll, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	AOView(petsclib::PetscLibType,ao::AO, viewer::PetscViewer) 
Displays an application ordering.

Collective

Input Parameters:
- `ao`     - the application ordering context
- `viewer` - viewer used for display

Level: intermediate

Options Database Key:
- `-ao_view` - calls `AOView()` at end of `AOCreate()`

-seealso: [](sec_ao), `AO`, `PetscViewerASCIIOpen()`, `AOViewFromOptions()`

# External Links
$(_doc_external("Vec/AOView"))
"""
function AOView(petsclib::PetscLibType, ao::AO, viewer::PetscViewer) end

@for_petsc function AOView(petsclib::$UnionPetscLib, ao::AO, viewer::PetscViewer )

    @chk ccall(
               (:AOView, $petsc_library),
               PetscErrorCode,
               (CAO, PetscViewer),
               ao, viewer,
              )


	return nothing
end 

"""
	AOViewFromOptions(petsclib::PetscLibType,ao::AO, obj::PetscObject, name::String) 
View an `AO` based on values in the options database

Collective

Input Parameters:
- `ao`   - the application ordering context
- `obj`  - optional object that provides the prefix used to search the options database
- `name` - command line option

Level: intermediate

-seealso: [](sec_ao), `AO`, `AOView()`, `PetscObjectViewFromOptions()`, `AOCreate()`

# External Links
$(_doc_external("Vec/AOViewFromOptions"))
"""
function AOViewFromOptions(petsclib::PetscLibType, ao::AO, obj::PetscObject, name::String) end

@for_petsc function AOViewFromOptions(petsclib::$UnionPetscLib, ao::AO, obj::PetscObject, name::String )

    @chk ccall(
               (:AOViewFromOptions, $petsc_library),
               PetscErrorCode,
               (CAO, PetscObject, Ptr{Cchar}),
               ao, obj, name,
              )


	return nothing
end 

"""
	AODestroy(petsclib::PetscLibType,ao::AO) 
Destroys an application ordering.

Collective

Input Parameter:
- `ao` - the application ordering context

Level: beginner

-seealso: [](sec_ao), `AO`, `AOCreate()`

# External Links
$(_doc_external("Vec/AODestroy"))
"""
function AODestroy(petsclib::PetscLibType, ao::AO) end

@for_petsc function AODestroy(petsclib::$UnionPetscLib, ao::AO )
	ao_ = Ref(ao.ptr)

    @chk ccall(
               (:AODestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{CAO},),
               ao_,
              )

	ao.ptr = C_NULL

	return nothing
end 

"""
	AOPetscToApplicationIS(petsclib::PetscLibType,ao::AO, is::IS) 
Maps an index set in the PETSc ordering to
the application-defined ordering.

Collective

Input Parameters:
- `ao` - the application ordering context
- `is` - the index set; this is replaced with its mapped values

Output Parameter:
- `is` - the mapped index set

Level: intermediate

-seealso: [](sec_ao), `AO`, `AOCreateBasic()`, `AOView()`, `AOApplicationToPetsc()`,
`AOApplicationToPetscIS()`, `AOPetscToApplication()`, `ISSTRIDE`, `ISBLOCK`

# External Links
$(_doc_external("Vec/AOPetscToApplicationIS"))
"""
function AOPetscToApplicationIS(petsclib::PetscLibType, ao::AO, is::IS) end

@for_petsc function AOPetscToApplicationIS(petsclib::$UnionPetscLib, ao::AO, is::IS )

    @chk ccall(
               (:AOPetscToApplicationIS, $petsc_library),
               PetscErrorCode,
               (CAO, CIS),
               ao, is,
              )


	return nothing
end 

"""
	AOApplicationToPetscIS(petsclib::PetscLibType,ao::AO, is::IS) 
Maps an index set in the application
ordering to the PETSc ordering.

Collective

Input Parameters:
- `ao` - the application ordering context
- `is` - the index set; this is replaced with its mapped values

Output Parameter:
- `is` - the mapped index set

Level: beginner

-seealso: [](sec_ao), `AO`, `AOCreateBasic()`, `AOView()`, `AOPetscToApplication()`,
`AOPetscToApplicationIS()`, `AOApplicationToPetsc()`, `ISSTRIDE`, `ISBLOCK`

# External Links
$(_doc_external("Vec/AOApplicationToPetscIS"))
"""
function AOApplicationToPetscIS(petsclib::PetscLibType, ao::AO, is::IS) end

@for_petsc function AOApplicationToPetscIS(petsclib::$UnionPetscLib, ao::AO, is::IS )

    @chk ccall(
               (:AOApplicationToPetscIS, $petsc_library),
               PetscErrorCode,
               (CAO, CIS),
               ao, is,
              )


	return nothing
end 

"""
	AOPetscToApplication(petsclib::PetscLibType,ao::AO, n::PetscInt, ia::Vector{PetscInt}) 
Maps a set of integers in the PETSc ordering to
the application-defined ordering.

Collective

Input Parameters:
- `ao` - the application ordering context
- `n`  - the number of integers
- `ia` - the integers; these are replaced with their mapped value

Output Parameter:
- `ia` - the mapped integers

Level: beginner

-seealso: [](sec_ao), `AO`, `AOCreateBasic()`, `AOView()`, `AOApplicationToPetsc()`,
`AOPetscToApplicationIS()`

# External Links
$(_doc_external("Vec/AOPetscToApplication"))
"""
function AOPetscToApplication(petsclib::PetscLibType, ao::AO, n::PetscInt, ia::Vector{PetscInt}) end

@for_petsc function AOPetscToApplication(petsclib::$UnionPetscLib, ao::AO, n::$PetscInt, ia::Vector{$PetscInt} )

    @chk ccall(
               (:AOPetscToApplication, $petsc_library),
               PetscErrorCode,
               (CAO, $PetscInt, Ptr{$PetscInt}),
               ao, n, ia,
              )


	return nothing
end 

"""
	AOApplicationToPetsc(petsclib::PetscLibType,ao::AO, n::PetscInt, ia::Vector{PetscInt}) 
Maps a set of integers in the application
ordering to the PETSc ordering.

Collective

Input Parameters:
- `ao` - the application ordering context
- `n`  - the number of integers
- `ia` - the integers; these are replaced with their mapped value

Output Parameter:
- `ia` - the mapped integers

Level: beginner

-seealso: [](sec_ao), `AOCreateBasic()`, `AOView()`, `AOPetscToApplication()`,
`AOPetscToApplicationIS()`

# External Links
$(_doc_external("Vec/AOApplicationToPetsc"))
"""
function AOApplicationToPetsc(petsclib::PetscLibType, ao::AO, n::PetscInt, ia::Vector{PetscInt}) end

@for_petsc function AOApplicationToPetsc(petsclib::$UnionPetscLib, ao::AO, n::$PetscInt, ia::Vector{$PetscInt} )

    @chk ccall(
               (:AOApplicationToPetsc, $petsc_library),
               PetscErrorCode,
               (CAO, $PetscInt, Ptr{$PetscInt}),
               ao, n, ia,
              )


	return nothing
end 

"""
	AOPetscToApplicationPermuteInt(petsclib::PetscLibType,ao::AO, block::PetscInt, array::Vector{PetscInt}) 
Permutes an array of blocks of integers
in the PETSc ordering to the application-defined ordering.

Collective

Input Parameters:
- `ao`    - The application ordering context
- `block` - The block size
- `array` - The integer array

Output Parameter:
- `array` - The permuted array

Level: beginner

-seealso: [](sec_ao), `AO`, `AOCreateBasic()`, `AOView()`, `AOApplicationToPetsc()`, `AOPetscToApplicationIS()`

# External Links
$(_doc_external("Vec/AOPetscToApplicationPermuteInt"))
"""
function AOPetscToApplicationPermuteInt(petsclib::PetscLibType, ao::AO, block::PetscInt, array::Vector{PetscInt}) end

@for_petsc function AOPetscToApplicationPermuteInt(petsclib::$UnionPetscLib, ao::AO, block::$PetscInt, array::Vector{$PetscInt} )

    @chk ccall(
               (:AOPetscToApplicationPermuteInt, $petsc_library),
               PetscErrorCode,
               (CAO, $PetscInt, Ptr{$PetscInt}),
               ao, block, array,
              )


	return nothing
end 

"""
	AOApplicationToPetscPermuteInt(petsclib::PetscLibType,ao::AO, block::PetscInt, array::Vector{PetscInt}) 
Permutes an array of blocks of integers
in the application-defined ordering to the PETSc ordering.

Collective

Input Parameters:
- `ao`    - The application ordering context
- `block` - The block size
- `array` - The integer array

Output Parameter:
- `array` - The permuted array

Level: beginner

-seealso: [](sec_ao), `AO`, `AOCreateBasic()`, `AOView()`, `AOPetscToApplicationIS()`, `AOApplicationToPetsc()`

# External Links
$(_doc_external("Vec/AOApplicationToPetscPermuteInt"))
"""
function AOApplicationToPetscPermuteInt(petsclib::PetscLibType, ao::AO, block::PetscInt, array::Vector{PetscInt}) end

@for_petsc function AOApplicationToPetscPermuteInt(petsclib::$UnionPetscLib, ao::AO, block::$PetscInt, array::Vector{$PetscInt} )

    @chk ccall(
               (:AOApplicationToPetscPermuteInt, $petsc_library),
               PetscErrorCode,
               (CAO, $PetscInt, Ptr{$PetscInt}),
               ao, block, array,
              )


	return nothing
end 

"""
	AOPetscToApplicationPermuteReal(petsclib::PetscLibType,ao::AO, block::PetscInt, array::Vector{PetscReal}) 
Permutes an array of blocks of reals
in the PETSc ordering to the application-defined ordering.

Collective

Input Parameters:
- `ao`    - The application ordering context
- `block` - The block size
- `array` - The integer array

Output Parameter:
- `array` - The permuted array

Level: beginner

-seealso: [](sec_ao), `AO`, `AOCreateBasic()`, `AOView()`, `AOApplicationToPetsc()`, `AOPetscToApplicationIS()`

# External Links
$(_doc_external("Vec/AOPetscToApplicationPermuteReal"))
"""
function AOPetscToApplicationPermuteReal(petsclib::PetscLibType, ao::AO, block::PetscInt, array::Vector{PetscReal}) end

@for_petsc function AOPetscToApplicationPermuteReal(petsclib::$UnionPetscLib, ao::AO, block::$PetscInt, array::Vector{$PetscReal} )

    @chk ccall(
               (:AOPetscToApplicationPermuteReal, $petsc_library),
               PetscErrorCode,
               (CAO, $PetscInt, Ptr{$PetscReal}),
               ao, block, array,
              )


	return nothing
end 

"""
	AOApplicationToPetscPermuteReal(petsclib::PetscLibType,ao::AO, block::PetscInt, array::Vector{PetscReal}) 
Permutes an array of blocks of reals
in the application-defined ordering to the PETSc ordering.

Collective

Input Parameters:
- `ao`    - The application ordering context
- `block` - The block size
- `array` - The integer array

Output Parameter:
- `array` - The permuted array

Level: beginner

-seealso: [](sec_ao), `AO`, `AOCreateBasic()`, `AOView()`, `AOApplicationToPetsc()`, `AOPetscToApplicationIS()`

# External Links
$(_doc_external("Vec/AOApplicationToPetscPermuteReal"))
"""
function AOApplicationToPetscPermuteReal(petsclib::PetscLibType, ao::AO, block::PetscInt, array::Vector{PetscReal}) end

@for_petsc function AOApplicationToPetscPermuteReal(petsclib::$UnionPetscLib, ao::AO, block::$PetscInt, array::Vector{$PetscReal} )

    @chk ccall(
               (:AOApplicationToPetscPermuteReal, $petsc_library),
               PetscErrorCode,
               (CAO, $PetscInt, Ptr{$PetscReal}),
               ao, block, array,
              )


	return nothing
end 

"""
	AOSetFromOptions(petsclib::PetscLibType,ao::AO) 
Sets `AO` options from the options database.

Collective

Input Parameter:
- `ao` - the application ordering

Options Database Key:
- `-ao_type <basic, memoryscalable>` - sets the type of the `AO`

Level: beginner

-seealso: [](sec_ao), `AO`, `AOCreate()`, `AOSetType()`, `AODestroy()`, `AOPetscToApplication()`, `AOApplicationToPetsc()`

# External Links
$(_doc_external("Vec/AOSetFromOptions"))
"""
function AOSetFromOptions(petsclib::PetscLibType, ao::AO) end

@for_petsc function AOSetFromOptions(petsclib::$UnionPetscLib, ao::AO )

    @chk ccall(
               (:AOSetFromOptions, $petsc_library),
               PetscErrorCode,
               (CAO,),
               ao,
              )


	return nothing
end 

"""
	AOSetIS(petsclib::PetscLibType,ao::AO, isapp::IS, ispetsc::IS) 
Sets the `IS` associated with the application ordering.

Collective

Input Parameters:
- `ao`      - the application ordering
- `isapp`   - index set that defines an ordering
- `ispetsc` - index set that defines another ordering (may be `NULL` to use the natural ordering)

Level: beginner

-seealso: [](sec_ao), [](sec_scatter), `AO`, `AOCreate()`, `AODestroy()`, `AOPetscToApplication()`, `AOApplicationToPetsc()`

# External Links
$(_doc_external("Vec/AOSetIS"))
"""
function AOSetIS(petsclib::PetscLibType, ao::AO, isapp::IS, ispetsc::IS) end

@for_petsc function AOSetIS(petsclib::$UnionPetscLib, ao::AO, isapp::IS, ispetsc::IS )

    @chk ccall(
               (:AOSetIS, $petsc_library),
               PetscErrorCode,
               (CAO, CIS, CIS),
               ao, isapp, ispetsc,
              )


	return nothing
end 

"""
	ao::AO = AOCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Creates an application ordering. That is an object that maps from an application ordering to a PETSc ordering and vice versa

Collective

Input Parameter:
- `comm` - MPI communicator that is to share the `AO`

Output Parameter:
- `ao` - the new application ordering

Options Database Key:
- `-ao_type <aotype>` - create `AO` with particular format
- `-ao_view`          - call `AOView()` at the conclusion of `AOCreate()`

Level: beginner

-seealso: [](sec_ao), `AO`, `AOView()`, `AOSetIS()`, `AODestroy()`, `AOPetscToApplication()`, `AOApplicationToPetsc()`

# External Links
$(_doc_external("Vec/AOCreate"))
"""
function AOCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function AOCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	ao_ = Ref{CAO}()

    @chk ccall(
               (:AOCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{CAO}),
               comm, ao_,
              )

	ao = AO(ao_[], petsclib)

	return ao
end 

"""
	hasIndex::PetscBool = AOMappingHasApplicationIndex(petsclib::PetscLibType,ao::AO, idex::PetscInt) 
Checks if an `AO` has a requested application index.

Not Collective

Input Parameters:
- `ao`   - The `AO`
- `idex` - The application index

Output Parameter:
- `hasIndex` - Flag is `PETSC_TRUE` if the index exists

Level: intermediate

-seealso: [](sec_ao), `AOMappingHasPetscIndex()`, `AOCreateMapping()`, `AO`

# External Links
$(_doc_external("Vec/AOMappingHasApplicationIndex"))
"""
function AOMappingHasApplicationIndex(petsclib::PetscLibType, ao::AO, idex::PetscInt) end

@for_petsc function AOMappingHasApplicationIndex(petsclib::$UnionPetscLib, ao::AO, idex::$PetscInt )
	hasIndex_ = Ref{PetscBool}()

    @chk ccall(
               (:AOMappingHasApplicationIndex, $petsc_library),
               PetscErrorCode,
               (CAO, $PetscInt, Ptr{PetscBool}),
               ao, idex, hasIndex_,
              )

	hasIndex = hasIndex_[]

	return hasIndex
end 

"""
	hasIndex::PetscBool = AOMappingHasPetscIndex(petsclib::PetscLibType,ao::AO, idex::PetscInt) 
checks if an `AO` has a requested PETSc index.

Not Collective

Input Parameters:
- `ao`   - The `AO`
- `idex` - The PETSc index

Output Parameter:
- `hasIndex` - Flag is `PETSC_TRUE` if the index exists

Level: intermediate

-seealso: [](sec_ao), `AOMappingHasApplicationIndex()`, `AOCreateMapping()`

# External Links
$(_doc_external("Vec/AOMappingHasPetscIndex"))
"""
function AOMappingHasPetscIndex(petsclib::PetscLibType, ao::AO, idex::PetscInt) end

@for_petsc function AOMappingHasPetscIndex(petsclib::$UnionPetscLib, ao::AO, idex::$PetscInt )
	hasIndex_ = Ref{PetscBool}()

    @chk ccall(
               (:AOMappingHasPetscIndex, $petsc_library),
               PetscErrorCode,
               (CAO, $PetscInt, Ptr{PetscBool}),
               ao, idex, hasIndex_,
              )

	hasIndex = hasIndex_[]

	return hasIndex
end 

"""
	aoout::AO = AOCreateMapping(petsclib::PetscLibType,comm::MPI_Comm, napp::PetscInt, myapp::Vector{PetscInt}, mypetsc::Vector{PetscInt}) 
Creates an application mapping using two integer arrays.

Input Parameters:
- `comm`    - MPI communicator that is to share the `AO`
- `napp`    - size of integer arrays
- `myapp`   - integer array that defines an ordering
- `mypetsc` - integer array that defines another ordering (may be `NULL` to indicate the identity ordering)

Output Parameter:
- `aoout` - the new application mapping

Options Database Key:
- `-ao_view` - call `AOView()` at the conclusion of `AOCreateMapping()`

Level: beginner

-seealso: [](sec_ao), `AOCreateBasic()`, `AOCreateMappingIS()`, `AODestroy()`

# External Links
$(_doc_external("Vec/AOCreateMapping"))
"""
function AOCreateMapping(petsclib::PetscLibType, comm::MPI_Comm, napp::PetscInt, myapp::Vector{PetscInt}, mypetsc::Vector{PetscInt}) end

@for_petsc function AOCreateMapping(petsclib::$UnionPetscLib, comm::MPI_Comm, napp::$PetscInt, myapp::Vector{$PetscInt}, mypetsc::Vector{$PetscInt} )
	aoout_ = Ref{CAO}()

    @chk ccall(
               (:AOCreateMapping, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{CAO}),
               comm, napp, myapp, mypetsc, aoout_,
              )

	aoout = AO(aoout_[], petsclib)

	return aoout
end 

"""
	aoout::AO = AOCreateMappingIS(petsclib::PetscLibType,isapp::IS, ispetsc::IS) 
Creates an application mapping using two index sets.

Input Parameters:
- `isapp`   - index set that defines an ordering
- `ispetsc` - index set that defines another ordering, maybe `NULL` for identity `IS`

Output Parameter:
- `aoout` - the new application ordering

Options Database Key:
- `-ao_view` - call `AOView()` at the conclusion of `AOCreateMappingIS()`

Level: beginner

-seealso: [](sec_ao), [](sec_scatter), `AOCreateBasic()`, `AOCreateMapping()`, `AODestroy()`

# External Links
$(_doc_external("Vec/AOCreateMappingIS"))
"""
function AOCreateMappingIS(petsclib::PetscLibType, isapp::IS, ispetsc::IS) end

@for_petsc function AOCreateMappingIS(petsclib::$UnionPetscLib, isapp::IS, ispetsc::IS )
	aoout_ = Ref{CAO}()

    @chk ccall(
               (:AOCreateMappingIS, $petsc_library),
               PetscErrorCode,
               (CIS, CIS, Ptr{CAO}),
               isapp, ispetsc, aoout_,
              )

	aoout = AO(aoout_[], petsclib)

	return aoout
end 

"""
	aoout::AO = AOCreateMemoryScalable(petsclib::PetscLibType,comm::MPI_Comm, napp::PetscInt, myapp::Vector{PetscInt}, mypetsc::Vector{PetscInt}) 
Creates a memory scalable application ordering using two integer arrays.

Collective

Input Parameters:
- `comm`    - MPI communicator that is to share the `AO`
- `napp`    - size of `myapp` and `mypetsc`
- `myapp`   - integer array that defines an ordering
- `mypetsc` - integer array that defines another ordering (may be `NULL` to indicate the natural ordering, that is 0,1,2,3,...)

Output Parameter:
- `aoout` - the new application ordering

Level: beginner

-seealso: [](sec_ao), [](sec_scatter), `AO`, `AOCreateMemoryScalableIS()`, `AODestroy()`, `AOPetscToApplication()`, `AOApplicationToPetsc()`

# External Links
$(_doc_external("Vec/AOCreateMemoryScalable"))
"""
function AOCreateMemoryScalable(petsclib::PetscLibType, comm::MPI_Comm, napp::PetscInt, myapp::Vector{PetscInt}, mypetsc::Vector{PetscInt}) end

@for_petsc function AOCreateMemoryScalable(petsclib::$UnionPetscLib, comm::MPI_Comm, napp::$PetscInt, myapp::Vector{$PetscInt}, mypetsc::Vector{$PetscInt} )
	aoout_ = Ref{CAO}()

    @chk ccall(
               (:AOCreateMemoryScalable, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{CAO}),
               comm, napp, myapp, mypetsc, aoout_,
              )

	aoout = AO(aoout_[], petsclib)

	return aoout
end 

"""
	aoout::AO = AOCreateMemoryScalableIS(petsclib::PetscLibType,isapp::IS, ispetsc::IS) 
Creates a memory scalable application ordering using two index sets.

Collective

Input Parameters:
- `isapp`   - index set that defines an ordering
- `ispetsc` - index set that defines another ordering (may be `NULL` to use the natural ordering)

Output Parameter:
- `aoout` - the new application ordering

Level: beginner

-seealso: [](sec_ao), [](sec_scatter), `AO`, `AOCreateBasicIS()`, `AOCreateMemoryScalable()`, `AODestroy()`

# External Links
$(_doc_external("Vec/AOCreateMemoryScalableIS"))
"""
function AOCreateMemoryScalableIS(petsclib::PetscLibType, isapp::IS, ispetsc::IS) end

@for_petsc function AOCreateMemoryScalableIS(petsclib::$UnionPetscLib, isapp::IS, ispetsc::IS )
	aoout_ = Ref{CAO}()

    @chk ccall(
               (:AOCreateMemoryScalableIS, $petsc_library),
               PetscErrorCode,
               (CIS, CIS, Ptr{CAO}),
               isapp, ispetsc, aoout_,
              )

	aoout = AO(aoout_[], petsclib)

	return aoout
end 

"""
	aoout::AO = AOCreateBasic(petsclib::PetscLibType,comm::MPI_Comm, napp::PetscInt, myapp::Vector{PetscInt}, mypetsc::Vector{PetscInt}) 
Creates a basic application ordering using two integer arrays.

Collective

Input Parameters:
- `comm`    - MPI communicator that is to share `AO`
- `napp`    - size of `myapp` and `mypetsc`
- `myapp`   - integer array that defines an ordering
- `mypetsc` - integer array that defines another ordering (may be `NULL` to
indicate the natural ordering, that is 0,1,2,3,...)

Output Parameter:
- `aoout` - the new application ordering

Level: beginner

-seealso: [](sec_ao), [](sec_scatter), `AO`, `AOCreateBasicIS()`, `AODestroy()`, `AOPetscToApplication()`, `AOApplicationToPetsc()`

# External Links
$(_doc_external("Vec/AOCreateBasic"))
"""
function AOCreateBasic(petsclib::PetscLibType, comm::MPI_Comm, napp::PetscInt, myapp::Vector{PetscInt}, mypetsc::Vector{PetscInt}) end

@for_petsc function AOCreateBasic(petsclib::$UnionPetscLib, comm::MPI_Comm, napp::$PetscInt, myapp::Vector{$PetscInt}, mypetsc::Vector{$PetscInt} )
	aoout_ = Ref{CAO}()

    @chk ccall(
               (:AOCreateBasic, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{CAO}),
               comm, napp, myapp, mypetsc, aoout_,
              )

	aoout = AO(aoout_[], petsclib)

	return aoout
end 

"""
	aoout::AO = AOCreateBasicIS(petsclib::PetscLibType,isapp::IS, ispetsc::IS) 
Creates a basic application ordering using two `IS` index sets.

Collective

Input Parameters:
- `isapp`   - index set that defines an ordering
- `ispetsc` - index set that defines another ordering (may be `NULL` to use the natural ordering)

Output Parameter:
- `aoout` - the new application ordering

Level: beginner

-seealso: [](sec_ao), [](sec_scatter), `IS`, `AO`, `AOCreateBasic()`, `AODestroy()`

# External Links
$(_doc_external("Vec/AOCreateBasicIS"))
"""
function AOCreateBasicIS(petsclib::PetscLibType, isapp::IS, ispetsc::IS) end

@for_petsc function AOCreateBasicIS(petsclib::$UnionPetscLib, isapp::IS, ispetsc::IS )
	aoout_ = Ref{CAO}()

    @chk ccall(
               (:AOCreateBasicIS, $petsc_library),
               PetscErrorCode,
               (CIS, CIS, Ptr{CAO}),
               isapp, ispetsc, aoout_,
              )

	aoout = AO(aoout_[], petsclib)

	return aoout
end 

