"""
	AOApplicationToPetsc(petsclib::PetscLibType, ao::AbstractAO, n::PetscInt, ia::Vector{PetscInt}) 
Maps a set of integers in the application-defined
ordering to the PETSc ordering.

Collective

Input Parameters:
- `ao` - the application ordering context
- `n`  - the number of integers
- `ia` - the integers; these are replaced with their mapped value

Output Parameter:
- `ia` - the mapped integers

Level: beginner

See also: `AOCreateBasic()`, `AOView()`, `AOPetscToApplication()`,
`AOPetscToApplicationIS()`

# External Links
$(_doc_external("AO/AOApplicationToPetsc"))
"""
function AOApplicationToPetsc(petsclib::PetscLibType, ao::AbstractAO, n::Integer, ia::AbstractVector{<:Number})
    error("AOApplicationToPetsc: no generated method for these argument types")
end

@for_petsc function AOApplicationToPetsc(petsclib::$UnionPetscLib, ao::AbstractAO, n::$PetscInt, ia::Vector{$PetscInt} )

    @chk ccall(
               (:AOApplicationToPetsc, $petsc_library),
               PetscErrorCode,
               (CAO, $PetscInt, Ptr{$PetscInt}),
               ao, n, ia,
              )


	return nothing
end 

"""
	AOApplicationToPetscIS(petsclib::PetscLibType, ao::AbstractAO, is::AbstractIS) 
Maps an index set in the application-defined
ordering to the PETSc ordering.

Collective

Input Parameters:
- `ao` - the application ordering context
- `is` - the index set; this is replaced with its mapped values

Output Parameter:
- `is` - the mapped index set

Level: beginner

See also: `AO`, `AOCreateBasic()`, `AOView()`, `AOPetscToApplication()`,
`AOPetscToApplicationIS()`, `AOApplicationToPetsc()`, `ISSTRIDE`, `ISBLOCK`

# External Links
$(_doc_external("AO/AOApplicationToPetscIS"))
"""
function AOApplicationToPetscIS(petsclib::PetscLibType, ao::AbstractAO, is::AbstractIS)
    error("AOApplicationToPetscIS: no generated method for these argument types")
end

@for_petsc function AOApplicationToPetscIS(petsclib::$UnionPetscLib, ao::AbstractAO, is::AbstractIS )

    @chk ccall(
               (:AOApplicationToPetscIS, $petsc_library),
               PetscErrorCode,
               (CAO, CIS),
               ao, is,
              )


	return nothing
end 

"""
	AOApplicationToPetscPermuteInt(petsclib::PetscLibType, ao::AbstractAO, block::PetscInt, array::Vector{PetscInt}) 
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

See also: `AO`, `AOCreateBasic()`, `AOView()`, `AOPetscToApplicationIS()`, `AOApplicationToPetsc()`

# External Links
$(_doc_external("AO/AOApplicationToPetscPermuteInt"))
"""
function AOApplicationToPetscPermuteInt(petsclib::PetscLibType, ao::AbstractAO, block::Integer, array::AbstractVector{<:Number})
    error("AOApplicationToPetscPermuteInt: no generated method for these argument types")
end

@for_petsc function AOApplicationToPetscPermuteInt(petsclib::$UnionPetscLib, ao::AbstractAO, block::$PetscInt, array::Vector{$PetscInt} )

    @chk ccall(
               (:AOApplicationToPetscPermuteInt, $petsc_library),
               PetscErrorCode,
               (CAO, $PetscInt, Ptr{$PetscInt}),
               ao, block, array,
              )


	return nothing
end 

"""
	AOApplicationToPetscPermuteReal(petsclib::PetscLibType, ao::AbstractAO, block::PetscInt, array::Vector{PetscReal}) 
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

See also: `AO`, `AOCreateBasic()`, `AOView()`, `AOApplicationToPetsc()`, `AOPetscToApplicationIS()`

# External Links
$(_doc_external("AO/AOApplicationToPetscPermuteReal"))
"""
function AOApplicationToPetscPermuteReal(petsclib::PetscLibType, ao::AbstractAO, block::Integer, array::AbstractVector{<:Number})
    error("AOApplicationToPetscPermuteReal: no generated method for these argument types")
end

@for_petsc function AOApplicationToPetscPermuteReal(petsclib::$UnionPetscLib, ao::AbstractAO, block::$PetscInt, array::Vector{$PetscReal} )

    @chk ccall(
               (:AOApplicationToPetscPermuteReal, $petsc_library),
               PetscErrorCode,
               (CAO, $PetscInt, Ptr{$PetscReal}),
               ao, block, array,
              )


	return nothing
end 

"""
	ao::AO = AOCreate(petsclib::PetscLibType, comm::MPI_Comm) 
Creates an application ordering. That is an object that maps from an application ordering to a PETSc ordering and vice versa

Collective

Input Parameter:
- `comm` - MPI communicator that is to share the `AO`

Output Parameter:
- `ao` - the new application ordering

Options Database Key:
- `-ao_type (basic|advanced|mapping|memoryscalable)` - Sets the `AO` type; see `AOType`
- `-ao_view`                                         - call `AOView()` at the conclusion of `AOCreate()`

Level: beginner

See also: `AO`, `AOView()`, `AOSetIS()`, `AODestroy()`, `AOPetscToApplication()`, `AOApplicationToPetsc()`

# External Links
$(_doc_external("AO/AOCreate"))
"""
function AOCreate(petsclib::PetscLibType, comm::MPI_Comm)
    error("AOCreate: no generated method for these argument types")
end

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
	aoout::AO = AOCreateBasic(petsclib::PetscLibType, comm::MPI_Comm, napp::PetscInt, myapp::Vector{PetscInt}, mypetsc::Vector{PetscInt}) 
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

See also: `AO`, `AOCreateBasicIS()`, `AODestroy()`, `AOPetscToApplication()`, `AOApplicationToPetsc()`

# External Links
$(_doc_external("AO/AOCreateBasic"))
"""
function AOCreateBasic(petsclib::PetscLibType, comm::MPI_Comm, napp::Integer, myapp::AbstractVector{<:Number}, mypetsc::AbstractVector{<:Number})
    error("AOCreateBasic: no generated method for these argument types")
end

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
	aoout::AO = AOCreateBasicIS(petsclib::PetscLibType, isapp::AbstractIS, ispetsc::AbstractIS) 
Creates a basic application ordering using two `IS` index sets.

Collective

Input Parameters:
- `isapp`   - index set that defines an ordering
- `ispetsc` - index set that defines another ordering (may be `NULL` to use the natural ordering)

Output Parameter:
- `aoout` - the new application ordering

Level: beginner

See also: `IS`, `AO`, `AOCreateBasic()`, `AODestroy()`

# External Links
$(_doc_external("AO/AOCreateBasicIS"))
"""
function AOCreateBasicIS(petsclib::PetscLibType, isapp::AbstractIS, ispetsc::AbstractIS)
    error("AOCreateBasicIS: no generated method for these argument types")
end

@for_petsc function AOCreateBasicIS(petsclib::$UnionPetscLib, isapp::AbstractIS, ispetsc::AbstractIS )
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

"""
	aoout::AO = AOCreateMapping(petsclib::PetscLibType, comm::MPI_Comm, napp::PetscInt, myapp::Vector{PetscInt}, mypetsc::Vector{PetscInt}) 
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

See also: `AOCreateBasic()`, `AOCreateMappingIS()`, `AODestroy()`

# External Links
$(_doc_external("AO/AOCreateMapping"))
"""
function AOCreateMapping(petsclib::PetscLibType, comm::MPI_Comm, napp::Integer, myapp::AbstractVector{<:Number}, mypetsc::AbstractVector{<:Number})
    error("AOCreateMapping: no generated method for these argument types")
end

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
	aoout::AO = AOCreateMappingIS(petsclib::PetscLibType, isapp::AbstractIS, ispetsc::AbstractIS) 
Creates an application mapping using two index sets.

Input Parameters:
- `isapp`   - index set that defines an ordering
- `ispetsc` - index set that defines another ordering, maybe `NULL` for identity `IS`

Output Parameter:
- `aoout` - the new application ordering

Options Database Key:
- `-ao_view` - call `AOView()` at the conclusion of `AOCreateMappingIS()`

Level: beginner

See also: `AOCreateBasic()`, `AOCreateMapping()`, `AODestroy()`

# External Links
$(_doc_external("AO/AOCreateMappingIS"))
"""
function AOCreateMappingIS(petsclib::PetscLibType, isapp::AbstractIS, ispetsc::AbstractIS)
    error("AOCreateMappingIS: no generated method for these argument types")
end

@for_petsc function AOCreateMappingIS(petsclib::$UnionPetscLib, isapp::AbstractIS, ispetsc::AbstractIS )
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
	aoout::AO = AOCreateMemoryScalable(petsclib::PetscLibType, comm::MPI_Comm, napp::PetscInt, myapp::Vector{PetscInt}, mypetsc::Vector{PetscInt}) 
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

See also: `AO`, `AOCreateMemoryScalableIS()`, `AODestroy()`, `AOPetscToApplication()`, `AOApplicationToPetsc()`

# External Links
$(_doc_external("AO/AOCreateMemoryScalable"))
"""
function AOCreateMemoryScalable(petsclib::PetscLibType, comm::MPI_Comm, napp::Integer, myapp::AbstractVector{<:Number}, mypetsc::AbstractVector{<:Number})
    error("AOCreateMemoryScalable: no generated method for these argument types")
end

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
	aoout::AO = AOCreateMemoryScalableIS(petsclib::PetscLibType, isapp::AbstractIS, ispetsc::AbstractIS) 
Creates a memory scalable application ordering using two index sets.

Collective

Input Parameters:
- `isapp`   - index set that defines an ordering
- `ispetsc` - index set that defines another ordering (may be `NULL` to use the natural ordering)

Output Parameter:
- `aoout` - the new application ordering

Level: beginner

See also: `AO`, `AOCreateBasicIS()`, `AOCreateMemoryScalable()`, `AODestroy()`

# External Links
$(_doc_external("AO/AOCreateMemoryScalableIS"))
"""
function AOCreateMemoryScalableIS(petsclib::PetscLibType, isapp::AbstractIS, ispetsc::AbstractIS)
    error("AOCreateMemoryScalableIS: no generated method for these argument types")
end

@for_petsc function AOCreateMemoryScalableIS(petsclib::$UnionPetscLib, isapp::AbstractIS, ispetsc::AbstractIS )
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
	AODestroy(petsclib::PetscLibType, ao::AbstractAO) 
Destroys an application ordering.

Collective

Input Parameter:
- `ao` - the application ordering context

Level: beginner

See also: `AO`, `AOCreate()`

# External Links
$(_doc_external("AO/AODestroy"))
"""
function AODestroy(petsclib::PetscLibType, ao::AbstractAO)
    error("AODestroy: no generated method for these argument types")
end

@for_petsc function AODestroy(petsclib::$UnionPetscLib, ao::AbstractAO )
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
	AOFinalizePackage(petsclib::PetscLibType) 
This function finalizes everything in the `AO` package. It is called
from `PetscFinalize()`.

Level: developer

See also: `AOInitializePackage()`, `PetscInitialize()`

# External Links
$(_doc_external("AO/AOFinalizePackage"))
"""
function AOFinalizePackage(petsclib::PetscLibType)
    error("AOFinalizePackage: no generated method for these argument types")
end

@for_petsc function AOFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:AOFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	type::String = AOGetType(petsclib::PetscLibType, ao::AbstractAO) 
Gets the `AO` type name (as a string) from the `AO`.

Not Collective

Input Parameter:
- `ao` - The vector

Output Parameter:
- `type` - The `AO` type name

Level: intermediate

See also: `AO`, `AOType`, `AOSetType()`, `AOCreate()`, `PetscObjectTypeCompare()`, `PetscObjectTypeCompareAny()`

# External Links
$(_doc_external("AO/AOGetType"))
"""
function AOGetType(petsclib::PetscLibType, ao::AbstractAO)
    error("AOGetType: no generated method for these argument types")
end

@for_petsc function AOGetType(petsclib::$UnionPetscLib, ao::AbstractAO )
	type_ = Ref{AOType}()

    @chk ccall(
               (:AOGetType, $petsc_library),
               PetscErrorCode,
               (CAO, Ptr{AOType}),
               ao, type_,
              )

	type = type_[] == C_NULL ? "" : unsafe_string(type_[])

	return type
end 

"""
	AOInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `AO` package. It is called
from `PetscDLLibraryRegister_petscvec()` when using dynamic libraries, and on the first call to `AOCreate()`
when using static or shared libraries.

Level: developer

See also: `AOFinalizePackage()`, `PetscInitialize()`

# External Links
$(_doc_external("Sys/AOInitializePackage"))
"""
function AOInitializePackage(petsclib::PetscLibType)
    error("AOInitializePackage: no generated method for these argument types")
end

@for_petsc function AOInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:AOInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	hasIndex::PetscBool = AOMappingHasApplicationIndex(petsclib::PetscLibType, ao::AbstractAO, idex::PetscInt) 
Checks if an `AO` has a requested application index.

Not Collective

Input Parameters:
- `ao`   - The `AO`
- `idex` - The application index

Output Parameter:
- `hasIndex` - Flag is `PETSC_TRUE` if the index exists

Level: intermediate

See also: `AOMappingHasPetscIndex()`, `AOCreateMapping()`, `AO`

# External Links
$(_doc_external("AO/AOMappingHasApplicationIndex"))
"""
function AOMappingHasApplicationIndex(petsclib::PetscLibType, ao::AbstractAO, idex::Integer)
    error("AOMappingHasApplicationIndex: no generated method for these argument types")
end

@for_petsc function AOMappingHasApplicationIndex(petsclib::$UnionPetscLib, ao::AbstractAO, idex::$PetscInt )
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
	hasIndex::PetscBool = AOMappingHasPetscIndex(petsclib::PetscLibType, ao::AbstractAO, idex::PetscInt) 
checks if an `AO` has a requested PETSc index.

Not Collective

Input Parameters:
- `ao`   - The `AO`
- `idex` - The PETSc index

Output Parameter:
- `hasIndex` - Flag is `PETSC_TRUE` if the index exists

Level: intermediate

See also: `AOMappingHasApplicationIndex()`, `AOCreateMapping()`

# External Links
$(_doc_external("AO/AOMappingHasPetscIndex"))
"""
function AOMappingHasPetscIndex(petsclib::PetscLibType, ao::AbstractAO, idex::Integer)
    error("AOMappingHasPetscIndex: no generated method for these argument types")
end

@for_petsc function AOMappingHasPetscIndex(petsclib::$UnionPetscLib, ao::AbstractAO, idex::$PetscInt )
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
	AOPetscToApplication(petsclib::PetscLibType, ao::AbstractAO, n::PetscInt, ia::Vector{PetscInt}) 
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

See also: `AO`, `AOCreateBasic()`, `AOView()`, `AOApplicationToPetsc()`,
`AOPetscToApplicationIS()`

# External Links
$(_doc_external("AO/AOPetscToApplication"))
"""
function AOPetscToApplication(petsclib::PetscLibType, ao::AbstractAO, n::Integer, ia::AbstractVector{<:Number})
    error("AOPetscToApplication: no generated method for these argument types")
end

@for_petsc function AOPetscToApplication(petsclib::$UnionPetscLib, ao::AbstractAO, n::$PetscInt, ia::Vector{$PetscInt} )

    @chk ccall(
               (:AOPetscToApplication, $petsc_library),
               PetscErrorCode,
               (CAO, $PetscInt, Ptr{$PetscInt}),
               ao, n, ia,
              )


	return nothing
end 

"""
	AOPetscToApplicationIS(petsclib::PetscLibType, ao::AbstractAO, is::AbstractIS) 
Maps an index set in the PETSc ordering to
the application-defined ordering.

Collective

Input Parameters:
- `ao` - the application ordering context
- `is` - the index set; this is replaced with its mapped values

Output Parameter:
- `is` - the mapped index set

Level: intermediate

See also: `AO`, `AOCreateBasic()`, `AOView()`, `AOApplicationToPetsc()`,
`AOApplicationToPetscIS()`, `AOPetscToApplication()`, `ISSTRIDE`, `ISBLOCK`

# External Links
$(_doc_external("AO/AOPetscToApplicationIS"))
"""
function AOPetscToApplicationIS(petsclib::PetscLibType, ao::AbstractAO, is::AbstractIS)
    error("AOPetscToApplicationIS: no generated method for these argument types")
end

@for_petsc function AOPetscToApplicationIS(petsclib::$UnionPetscLib, ao::AbstractAO, is::AbstractIS )

    @chk ccall(
               (:AOPetscToApplicationIS, $petsc_library),
               PetscErrorCode,
               (CAO, CIS),
               ao, is,
              )


	return nothing
end 

"""
	AOPetscToApplicationPermuteInt(petsclib::PetscLibType, ao::AbstractAO, block::PetscInt, array::Vector{PetscInt}) 
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

See also: `AO`, `AOCreateBasic()`, `AOView()`, `AOApplicationToPetsc()`, `AOPetscToApplicationIS()`

# External Links
$(_doc_external("AO/AOPetscToApplicationPermuteInt"))
"""
function AOPetscToApplicationPermuteInt(petsclib::PetscLibType, ao::AbstractAO, block::Integer, array::AbstractVector{<:Number})
    error("AOPetscToApplicationPermuteInt: no generated method for these argument types")
end

@for_petsc function AOPetscToApplicationPermuteInt(petsclib::$UnionPetscLib, ao::AbstractAO, block::$PetscInt, array::Vector{$PetscInt} )

    @chk ccall(
               (:AOPetscToApplicationPermuteInt, $petsc_library),
               PetscErrorCode,
               (CAO, $PetscInt, Ptr{$PetscInt}),
               ao, block, array,
              )


	return nothing
end 

"""
	AOPetscToApplicationPermuteReal(petsclib::PetscLibType, ao::AbstractAO, block::PetscInt, array::Vector{PetscReal}) 
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

See also: `AO`, `AOCreateBasic()`, `AOView()`, `AOApplicationToPetsc()`, `AOPetscToApplicationIS()`

# External Links
$(_doc_external("AO/AOPetscToApplicationPermuteReal"))
"""
function AOPetscToApplicationPermuteReal(petsclib::PetscLibType, ao::AbstractAO, block::Integer, array::AbstractVector{<:Number})
    error("AOPetscToApplicationPermuteReal: no generated method for these argument types")
end

@for_petsc function AOPetscToApplicationPermuteReal(petsclib::$UnionPetscLib, ao::AbstractAO, block::$PetscInt, array::Vector{$PetscReal} )

    @chk ccall(
               (:AOPetscToApplicationPermuteReal, $petsc_library),
               PetscErrorCode,
               (CAO, $PetscInt, Ptr{$PetscReal}),
               ao, block, array,
              )


	return nothing
end 

"""
	AORegister(petsclib::PetscLibType, sname::String, fnc::external) 
Register  an application ordering method

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - the name (`AOType`) of the `AO` scheme
- `function` - the create routine for the application ordering method

Level: advanced

See also: `AO`, `AOType`, `AOCreate()`, `AORegisterAll()`, `AOBASIC`, `AOADVANCED`, `AOMAPPING`, `AOMEMORYSCALABLE`

# External Links
$(_doc_external("AO/AORegister"))
"""
function AORegister(petsclib::PetscLibType, sname::String, fnc::external)
    error("AORegister: no generated method for these argument types")
end

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

See also: `AO`, `AOType`, `AORegister()`, `AORegisterDestroy()`

# External Links
$(_doc_external("AO/AORegisterAll"))
"""
function AORegisterAll(petsclib::PetscLibType)
    error("AORegisterAll: no generated method for these argument types")
end

@for_petsc function AORegisterAll(petsclib::$UnionPetscLib)

    @chk ccall(
               (:AORegisterAll, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	AOSetFromOptions(petsclib::PetscLibType, ao::AbstractAO) 
Sets `AO` options from the options database.

Collective

Input Parameter:
- `ao` - the application ordering

Options Database Key:
- `-ao_type (basic|memoryscalable)` - sets the type of the `AO`

Level: beginner

See also: `AO`, `AOCreate()`, `AOSetType()`, `AODestroy()`, `AOPetscToApplication()`, `AOApplicationToPetsc()`

# External Links
$(_doc_external("AO/AOSetFromOptions"))
"""
function AOSetFromOptions(petsclib::PetscLibType, ao::AbstractAO)
    error("AOSetFromOptions: no generated method for these argument types")
end

@for_petsc function AOSetFromOptions(petsclib::$UnionPetscLib, ao::AbstractAO )

    @chk ccall(
               (:AOSetFromOptions, $petsc_library),
               PetscErrorCode,
               (CAO,),
               ao,
              )


	return nothing
end 

"""
	AOSetIS(petsclib::PetscLibType, ao::AbstractAO, isapp::AbstractIS, ispetsc::AbstractIS) 
Sets the `IS` associated with the application ordering.

Collective

Input Parameters:
- `ao`      - the application ordering
- `isapp`   - index set that defines an ordering
- `ispetsc` - index set that defines another ordering (may be `NULL` to use the natural ordering)

Level: beginner

See also: `AO`, `AOCreate()`, `AODestroy()`, `AOPetscToApplication()`, `AOApplicationToPetsc()`

# External Links
$(_doc_external("AO/AOSetIS"))
"""
function AOSetIS(petsclib::PetscLibType, ao::AbstractAO, isapp::AbstractIS, ispetsc::AbstractIS)
    error("AOSetIS: no generated method for these argument types")
end

@for_petsc function AOSetIS(petsclib::$UnionPetscLib, ao::AbstractAO, isapp::AbstractIS, ispetsc::AbstractIS )

    @chk ccall(
               (:AOSetIS, $petsc_library),
               PetscErrorCode,
               (CAO, CIS, CIS),
               ao, isapp, ispetsc,
              )


	return nothing
end 

"""
	AOSetType(petsclib::PetscLibType, ao::AbstractAO, method::String) 
Builds an application ordering for a particular `AOType`

Collective

Input Parameters:
- `ao`     - The `AO` object
- `method` - The name of the AO type

Options Database Key:
- `-ao_type (basic|advanced|mapping|memoryscalable)` - Sets the `AO` type; see `AOType`

Level: intermediate

See also: `AO`, `AOType`, `AOCreateBasic()`, `AOCreateMemoryScalable()`, `AOGetType()`, `AOCreate()`

# External Links
$(_doc_external("AO/AOSetType"))
"""
function AOSetType(petsclib::PetscLibType, ao::AbstractAO, method::String)
    error("AOSetType: no generated method for these argument types")
end

@for_petsc function AOSetType(petsclib::$UnionPetscLib, ao::AbstractAO, method::String )

    @chk ccall(
               (:AOSetType, $petsc_library),
               PetscErrorCode,
               (CAO, AOType),
               ao, method,
              )


	return nothing
end 

"""
	AOView(petsclib::PetscLibType, ao::AbstractAO, viewer::PetscViewer) 
Displays an application ordering.

Collective

Input Parameters:
- `ao`     - the application ordering context
- `viewer` - viewer used for display

Options Database Key:
- `-ao_view` - calls `AOView()` at end of `AOCreate()`

Level: intermediate

See also: `AO`, `PetscViewer`, `PetscViewerASCIIOpen()`, `AOViewFromOptions()`

# External Links
$(_doc_external("AO/AOView"))
"""
function AOView(petsclib::PetscLibType, ao::AbstractAO, viewer::PetscViewer)
    error("AOView: no generated method for these argument types")
end

@for_petsc function AOView(petsclib::$UnionPetscLib, ao::AbstractAO, viewer::PetscViewer )

    @chk ccall(
               (:AOView, $petsc_library),
               PetscErrorCode,
               (CAO, PetscViewer),
               ao, viewer,
              )


	return nothing
end 

"""
	AOViewFromOptions(petsclib::PetscLibType, ao::AbstractAO, obj, name::String) 
View an `AO` based on values in the options database

Collective

Input Parameters:
- `ao`   - the application ordering context
- `obj`  - optional object that provides the prefix used to search the options database
- `name` - command line option

Options Database Key:
- `-name [viewertype][:...]` - option name and values. See `PetscObjectViewFromOptions()` for the possible arguments

Level: intermediate

See also: `AO`, `AOView()`, `PetscObjectViewFromOptions()`, `AOCreate()`

# External Links
$(_doc_external("AO/AOViewFromOptions"))
"""
function AOViewFromOptions(petsclib::PetscLibType, ao::AbstractAO, obj, name::String)
    error("AOViewFromOptions: no generated method for these argument types")
end

@for_petsc function AOViewFromOptions(petsclib::$UnionPetscLib, ao::AbstractAO, obj, name::String )

    @chk ccall(
               (:AOViewFromOptions, $petsc_library),
               PetscErrorCode,
               (CAO, PetscObject, Ptr{Cchar}),
               ao, obj, name,
              )


	return nothing
end 

