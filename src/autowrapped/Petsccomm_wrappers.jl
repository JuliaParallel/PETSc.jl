# autodefined type arguments for class ------
mutable struct _n_PetscSubcomm end
const PetscSubcomm = Ptr{_n_PetscSubcomm}

# -------------------------------------------------------
# autodefined type arguments for class ------
mutable struct _n_PetscShmComm end
const PetscShmComm = Ptr{_n_PetscShmComm}
# -------------------------------------------------------

"""
	PetscSubcommSetFromOptions(petsclib::PetscLibType,psubcomm::PetscSubcomm) 
Allows setting options for a `PetscSubcomm`

Collective

Input Parameter:
- `psubcomm` - `PetscSubcomm` context

Level: beginner

-seealso: `PetscSubcomm`, `PetscSubcommCreate()`

# External Links
$(_doc_external("Sys/PetscSubcommSetFromOptions"))
"""
function PetscSubcommSetFromOptions(petsclib::PetscLibType, psubcomm::PetscSubcomm) end

@for_petsc function PetscSubcommSetFromOptions(petsclib::$UnionPetscLib, psubcomm::PetscSubcomm )

    @chk ccall(
               (:PetscSubcommSetFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscSubcomm,),
               psubcomm,
              )


	return nothing
end 

"""
	PetscSubcommSetOptionsPrefix(petsclib::PetscLibType,psubcomm::PetscSubcomm, pre::String) 
Sets the prefix used for searching for options in the options database for this object

Logically Collective

Level: intermediate

Input Parameters:
- `psubcomm` - `PetscSubcomm` context
- `pre`      - the prefix to prepend all `PetscSubcomm` item names with.

-seealso: `PetscSubcomm`, `PetscSubcommCreate()`

# External Links
$(_doc_external("Sys/PetscSubcommSetOptionsPrefix"))
"""
function PetscSubcommSetOptionsPrefix(petsclib::PetscLibType, psubcomm::PetscSubcomm, pre::String) end

@for_petsc function PetscSubcommSetOptionsPrefix(petsclib::$UnionPetscLib, psubcomm::PetscSubcomm, pre::String )

    @chk ccall(
               (:PetscSubcommSetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (PetscSubcomm, Ptr{Cchar}),
               psubcomm, pre,
              )


	return nothing
end 

"""
	PetscSubcommView(petsclib::PetscLibType,psubcomm::PetscSubcomm, viewer::PetscViewer) 
Views a `PetscSubcomm`

Collective

Input Parameters:
- `psubcomm` - `PetscSubcomm` context
- `viewer`   - `PetscViewer` to display the information

Level: beginner

-seealso: `PetscSubcomm`, `PetscSubcommCreate()`, `PetscViewer`

# External Links
$(_doc_external("Sys/PetscSubcommView"))
"""
function PetscSubcommView(petsclib::PetscLibType, psubcomm::PetscSubcomm, viewer::PetscViewer) end

@for_petsc function PetscSubcommView(petsclib::$UnionPetscLib, psubcomm::PetscSubcomm, viewer::PetscViewer )

    @chk ccall(
               (:PetscSubcommView, $petsc_library),
               PetscErrorCode,
               (PetscSubcomm, PetscViewer),
               psubcomm, viewer,
              )


	return nothing
end 

"""
	PetscSubcommSetNumber(petsclib::PetscLibType,psubcomm::PetscSubcomm, nsubcomm::PetscInt) 
Set total number of subcommunicators desired in the given `PetscSubcomm`

Collective

Input Parameters:
- `psubcomm` - `PetscSubcomm` context
- `nsubcomm` - the total number of subcommunicators in psubcomm

Level: advanced

-seealso: `PetscSubcomm`, `PetscSubcommCreate()`, `PetscSubcommDestroy()`, `PetscSubcommSetType()`, `PetscSubcommSetTypeGeneral()`

# External Links
$(_doc_external("Sys/PetscSubcommSetNumber"))
"""
function PetscSubcommSetNumber(petsclib::PetscLibType, psubcomm::PetscSubcomm, nsubcomm::PetscInt) end

@for_petsc function PetscSubcommSetNumber(petsclib::$UnionPetscLib, psubcomm::PetscSubcomm, nsubcomm::$PetscInt )

    @chk ccall(
               (:PetscSubcommSetNumber, $petsc_library),
               PetscErrorCode,
               (PetscSubcomm, $PetscInt),
               psubcomm, nsubcomm,
              )


	return nothing
end 

"""
	PetscSubcommSetType(petsclib::PetscLibType,psubcomm::PetscSubcomm, subcommtype::PetscSubcommType) 
Set the way the original MPI communicator is divided up in the `PetscSubcomm`

Collective

Input Parameters:
- `psubcomm`    - `PetscSubcomm` context
- `subcommtype` - `PetscSubcommType` `PETSC_SUBCOMM_CONTIGUOUS` or `PETSC_SUBCOMM_INTERLACED`

Level: advanced

-seealso: `PetscSubcommType`, `PETSC_SUBCOMM_CONTIGUOUS`, `PETSC_SUBCOMM_INTERLACED`,
`PetscSubcommCreate()`, `PetscSubcommDestroy()`, `PetscSubcommSetNumber()`, `PetscSubcommSetTypeGeneral()`

# External Links
$(_doc_external("Sys/PetscSubcommSetType"))
"""
function PetscSubcommSetType(petsclib::PetscLibType, psubcomm::PetscSubcomm, subcommtype::PetscSubcommType) end

@for_petsc function PetscSubcommSetType(petsclib::$UnionPetscLib, psubcomm::PetscSubcomm, subcommtype::PetscSubcommType )

    @chk ccall(
               (:PetscSubcommSetType, $petsc_library),
               PetscErrorCode,
               (PetscSubcomm, PetscSubcommType),
               psubcomm, subcommtype,
              )


	return nothing
end 

"""
	PetscSubcommSetTypeGeneral(petsclib::PetscLibType,psubcomm::PetscSubcomm, color::PetscMPIInt, subrank::PetscMPIInt) 
Divides up a communicator based on a specific user's specification

Collective

Input Parameters:
- `psubcomm` - `PetscSubcomm` context
- `color`    - control of subset assignment (nonnegative integer). Processes with the same color are in the same subcommunicator.
- `subrank`  - rank in the subcommunicator

Level: advanced

-seealso: `PetscSubcommType`, `PETSC_SUBCOMM_CONTIGUOUS`, `PETSC_SUBCOMM_INTERLACED`, `PetscSubcommCreate()`, `PetscSubcommDestroy()`, `PetscSubcommSetNumber()`, `PetscSubcommSetType()`

# External Links
$(_doc_external("Sys/PetscSubcommSetTypeGeneral"))
"""
function PetscSubcommSetTypeGeneral(petsclib::PetscLibType, psubcomm::PetscSubcomm, color::PetscMPIInt, subrank::PetscMPIInt) end

@for_petsc function PetscSubcommSetTypeGeneral(petsclib::$UnionPetscLib, psubcomm::PetscSubcomm, color::PetscMPIInt, subrank::PetscMPIInt )

    @chk ccall(
               (:PetscSubcommSetTypeGeneral, $petsc_library),
               PetscErrorCode,
               (PetscSubcomm, PetscMPIInt, PetscMPIInt),
               psubcomm, color, subrank,
              )


	return nothing
end 

"""
	PetscSubcommDestroy(petsclib::PetscLibType,psubcomm::PetscSubcomm) 
Destroys a `PetscSubcomm` object

Collective

Input Parameter:
- `psubcomm` - the `PetscSubcomm` context

Level: advanced

-seealso: `PetscSubcommCreate()`, `PetscSubcommSetType()`

# External Links
$(_doc_external("Sys/PetscSubcommDestroy"))
"""
function PetscSubcommDestroy(petsclib::PetscLibType, psubcomm::PetscSubcomm) end

@for_petsc function PetscSubcommDestroy(petsclib::$UnionPetscLib, psubcomm::PetscSubcomm )

    @chk ccall(
               (:PetscSubcommDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscSubcomm},),
               psubcomm,
              )


	return nothing
end 

"""
	psubcomm::PetscSubcomm = PetscSubcommCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Create a `PetscSubcomm` context. This object is used to manage the division of a `MPI_Comm` into subcommunicators

Collective

Input Parameter:
- `comm` - MPI communicator

Output Parameter:
- `psubcomm` - location to store the `PetscSubcomm` context

Level: advanced

-seealso: `PetscSubcomm`, `PetscSubcommDestroy()`, `PetscSubcommSetTypeGeneral()`, `PetscSubcommSetFromOptions()`, `PetscSubcommSetType()`,
`PetscSubcommSetNumber()`

# External Links
$(_doc_external("Sys/PetscSubcommCreate"))
"""
function PetscSubcommCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function PetscSubcommCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	psubcomm_ = Ref{PetscSubcomm}()

    @chk ccall(
               (:PetscSubcommCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscSubcomm}),
               comm, psubcomm_,
              )

	psubcomm = psubcomm_[]

	return psubcomm
end 

"""
	PetscSubcommGetParent(petsclib::PetscLibType,scomm::PetscSubcomm, pcomm::MPI_Comm) 
Gets the communicator that was used to create the `PetscSubcomm`

Collective

Input Parameter:
- `scomm` - the `PetscSubcomm`

Output Parameter:
- `pcomm` - location to store the parent communicator

Level: intermediate

-seealso: `PetscSubcommDestroy()`, `PetscSubcommSetTypeGeneral()`, `PetscSubcommSetFromOptions()`, `PetscSubcommSetType()`,
`PetscSubcommSetNumber()`, `PetscSubcommGetChild()`, `PetscSubcommContiguousParent()`

# External Links
$(_doc_external("Sys/PetscSubcommGetParent"))
"""
function PetscSubcommGetParent(petsclib::PetscLibType, scomm::PetscSubcomm, pcomm::MPI_Comm) end

@for_petsc function PetscSubcommGetParent(petsclib::$UnionPetscLib, scomm::PetscSubcomm, pcomm::MPI_Comm )

    @chk ccall(
               (:PetscSubcommGetParent, $petsc_library),
               PetscErrorCode,
               (PetscSubcomm, Ptr{MPI_Comm}),
               scomm, pcomm,
              )


	return nothing
end 

"""
	PetscSubcommGetContiguousParent(petsclib::PetscLibType,scomm::PetscSubcomm, pcomm::MPI_Comm) 
Gets a communicator that is a duplicate of the parent but has the ranks
reordered by the order they are in the children

Collective

Input Parameter:
- `scomm` - the `PetscSubcomm`

Output Parameter:
- `pcomm` - location to store the parent communicator

Level: intermediate

-seealso: `PetscSubcommDestroy()`, `PetscSubcommSetTypeGeneral()`, `PetscSubcommSetFromOptions()`, `PetscSubcommSetType()`,
`PetscSubcommSetNumber()`, `PetscSubcommGetChild()`, `PetscSubcommContiguousParent()`

# External Links
$(_doc_external("Sys/PetscSubcommGetContiguousParent"))
"""
function PetscSubcommGetContiguousParent(petsclib::PetscLibType, scomm::PetscSubcomm, pcomm::MPI_Comm) end

@for_petsc function PetscSubcommGetContiguousParent(petsclib::$UnionPetscLib, scomm::PetscSubcomm, pcomm::MPI_Comm )

    @chk ccall(
               (:PetscSubcommGetContiguousParent, $petsc_library),
               PetscErrorCode,
               (PetscSubcomm, Ptr{MPI_Comm}),
               scomm, pcomm,
              )


	return nothing
end 

"""
	PetscSubcommGetChild(petsclib::PetscLibType,scomm::PetscSubcomm, ccomm::MPI_Comm) 
Gets the communicator created by the `PetscSubcomm`. This is part of one of the subcommunicators created by the `PetscSubcomm`

Collective

Input Parameter:
- `scomm` - the `PetscSubcomm`

Output Parameter:
- `ccomm` - location to store the child communicator

Level: intermediate

-seealso: `PetscSubcommDestroy()`, `PetscSubcommSetTypeGeneral()`, `PetscSubcommSetFromOptions()`, `PetscSubcommSetType()`,
`PetscSubcommSetNumber()`, `PetscSubcommGetParent()`, `PetscSubcommContiguousParent()`

# External Links
$(_doc_external("Sys/PetscSubcommGetChild"))
"""
function PetscSubcommGetChild(petsclib::PetscLibType, scomm::PetscSubcomm, ccomm::MPI_Comm) end

@for_petsc function PetscSubcommGetChild(petsclib::$UnionPetscLib, scomm::PetscSubcomm, ccomm::MPI_Comm )

    @chk ccall(
               (:PetscSubcommGetChild, $petsc_library),
               PetscErrorCode,
               (PetscSubcomm, Ptr{MPI_Comm}),
               scomm, ccomm,
              )


	return nothing
end 

"""
	PetscShmCommGet(petsclib::PetscLibType,globcomm::MPI_Comm, pshmcomm::PetscShmComm) 
Returns a sub

Collective.

Input Parameter:
- `globcomm` - `MPI_Comm`, which can be a user `MPI_Comm` or a PETSc inner `MPI_Comm`

Output Parameter:
- `pshmcomm` - the PETSc shared memory communicator object

Level: developer

-seealso: `PetscShmCommGlobalToLocal()`, `PetscShmCommLocalToGlobal()`, `PetscShmCommGetMpiShmComm()`

# External Links
$(_doc_external("Sys/PetscShmCommGet"))
"""
function PetscShmCommGet(petsclib::PetscLibType, globcomm::MPI_Comm, pshmcomm::PetscShmComm) end

@for_petsc function PetscShmCommGet(petsclib::$UnionPetscLib, globcomm::MPI_Comm, pshmcomm::PetscShmComm )

    @chk ccall(
               (:PetscShmCommGet, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscShmComm}),
               globcomm, pshmcomm,
              )


	return nothing
end 

"""
	PetscShmCommGlobalToLocal(petsclib::PetscLibType,pshmcomm::PetscShmComm, grank::PetscMPIInt, lrank::PetscMPIInt) 
Given a global rank returns the local rank in the shared memory communicator

Input Parameters:
- `pshmcomm` - the shared memory communicator object
- `grank`    - the global rank

Output Parameter:
- `lrank` - the local rank, or `MPI_PROC_NULL` if it does not exist

Level: developer

-seealso: `PetscShmCommGet()`, `PetscShmCommLocalToGlobal()`, `PetscShmCommGetMpiShmComm()`

# External Links
$(_doc_external("Sys/PetscShmCommGlobalToLocal"))
"""
function PetscShmCommGlobalToLocal(petsclib::PetscLibType, pshmcomm::PetscShmComm, grank::PetscMPIInt, lrank::PetscMPIInt) end

@for_petsc function PetscShmCommGlobalToLocal(petsclib::$UnionPetscLib, pshmcomm::PetscShmComm, grank::PetscMPIInt, lrank::PetscMPIInt )

    @chk ccall(
               (:PetscShmCommGlobalToLocal, $petsc_library),
               PetscErrorCode,
               (PetscShmComm, PetscMPIInt, Ptr{PetscMPIInt}),
               pshmcomm, grank, lrank,
              )


	return nothing
end 

"""
	PetscShmCommLocalToGlobal(petsclib::PetscLibType,pshmcomm::PetscShmComm, lrank::PetscMPIInt, grank::PetscMPIInt) 
Given a local rank in the shared memory communicator returns the global rank

Input Parameters:
- `pshmcomm` - the shared memory communicator object
- `lrank`    - the local rank in the shared memory communicator

Output Parameter:
- `grank` - the global rank in the global communicator where the shared memory communicator is built

Level: developer

-seealso: `PetscShmCommGlobalToLocal()`, `PetscShmCommGet()`, `PetscShmCommGetMpiShmComm()`

# External Links
$(_doc_external("Sys/PetscShmCommLocalToGlobal"))
"""
function PetscShmCommLocalToGlobal(petsclib::PetscLibType, pshmcomm::PetscShmComm, lrank::PetscMPIInt, grank::PetscMPIInt) end

@for_petsc function PetscShmCommLocalToGlobal(petsclib::$UnionPetscLib, pshmcomm::PetscShmComm, lrank::PetscMPIInt, grank::PetscMPIInt )

    @chk ccall(
               (:PetscShmCommLocalToGlobal, $petsc_library),
               PetscErrorCode,
               (PetscShmComm, PetscMPIInt, Ptr{PetscMPIInt}),
               pshmcomm, lrank, grank,
              )


	return nothing
end 

"""
	PetscShmCommGetMpiShmComm(petsclib::PetscLibType,pshmcomm::PetscShmComm, comm::MPI_Comm) 
Returns the MPI communicator that represents all processes with common shared memory

Input Parameter:
- `pshmcomm` - PetscShmComm object obtained with PetscShmCommGet()

Output Parameter:
- `comm` - the MPI communicator

Level: developer

-seealso: `PetscShmCommGlobalToLocal()`, `PetscShmCommGet()`, `PetscShmCommLocalToGlobal()`

# External Links
$(_doc_external("Sys/PetscShmCommGetMpiShmComm"))
"""
function PetscShmCommGetMpiShmComm(petsclib::PetscLibType, pshmcomm::PetscShmComm, comm::MPI_Comm) end

@for_petsc function PetscShmCommGetMpiShmComm(petsclib::$UnionPetscLib, pshmcomm::PetscShmComm, comm::MPI_Comm )

    @chk ccall(
               (:PetscShmCommGetMpiShmComm, $petsc_library),
               PetscErrorCode,
               (PetscShmComm, Ptr{MPI_Comm}),
               pshmcomm, comm,
              )


	return nothing
end 

