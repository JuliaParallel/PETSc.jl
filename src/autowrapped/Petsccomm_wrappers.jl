"""
	pshmcomm::PetscShmComm = PetscShmCommGet(petsclib::PetscLibType,globcomm::MPI_Comm) 
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
function PetscShmCommGet(petsclib::PetscLibType, globcomm::MPI_Comm) end

@for_petsc function PetscShmCommGet(petsclib::$UnionPetscLib, globcomm::MPI_Comm )
	pshmcomm_ = Ref{PetscShmComm}()

    @chk ccall(
               (:PetscShmCommGet, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscShmComm}),
               globcomm, pshmcomm_,
              )

	pshmcomm = pshmcomm_[]

	return pshmcomm
end 

"""
	comm::MPI_Comm = PetscShmCommGetMpiShmComm(petsclib::PetscLibType,pshmcomm::PetscShmComm) 
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
function PetscShmCommGetMpiShmComm(petsclib::PetscLibType, pshmcomm::PetscShmComm) end

@for_petsc function PetscShmCommGetMpiShmComm(petsclib::$UnionPetscLib, pshmcomm::PetscShmComm )
	comm_ = Ref{MPI.MPI_Comm}()

    @chk ccall(
               (:PetscShmCommGetMpiShmComm, $petsc_library),
               PetscErrorCode,
               (PetscShmComm, Ptr{MPI.MPI_Comm}),
               pshmcomm, comm_,
              )

	comm = MPI.Comm(comm_[])

	return comm
end 

"""
	lrank::PetscMPIInt = PetscShmCommGlobalToLocal(petsclib::PetscLibType,pshmcomm::PetscShmComm, grank::PetscMPIInt) 
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
function PetscShmCommGlobalToLocal(petsclib::PetscLibType, pshmcomm::PetscShmComm, grank::PetscMPIInt) end

@for_petsc function PetscShmCommGlobalToLocal(petsclib::$UnionPetscLib, pshmcomm::PetscShmComm, grank::PetscMPIInt )
	lrank_ = Ref{PetscMPIInt}()

    @chk ccall(
               (:PetscShmCommGlobalToLocal, $petsc_library),
               PetscErrorCode,
               (PetscShmComm, PetscMPIInt, Ptr{PetscMPIInt}),
               pshmcomm, grank, lrank_,
              )

	lrank = lrank_[]

	return lrank
end 

"""
	grank::PetscMPIInt = PetscShmCommLocalToGlobal(petsclib::PetscLibType,pshmcomm::PetscShmComm, lrank::PetscMPIInt) 
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
function PetscShmCommLocalToGlobal(petsclib::PetscLibType, pshmcomm::PetscShmComm, lrank::PetscMPIInt) end

@for_petsc function PetscShmCommLocalToGlobal(petsclib::$UnionPetscLib, pshmcomm::PetscShmComm, lrank::PetscMPIInt )
	grank_ = Ref{PetscMPIInt}()

    @chk ccall(
               (:PetscShmCommLocalToGlobal, $petsc_library),
               PetscErrorCode,
               (PetscShmComm, PetscMPIInt, Ptr{PetscMPIInt}),
               pshmcomm, lrank, grank_,
              )

	grank = grank_[]

	return grank
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
	PetscSubcommDestroy(petsclib::PetscLibType,psubcomm::Union{PetscSubcomm, Ref{PetscSubcomm}}) 
Destroys a `PetscSubcomm` object

Collective

Input Parameter:
- `psubcomm` - the `PetscSubcomm` context

Level: advanced

-seealso: `PetscSubcommCreate()`, `PetscSubcommSetType()`

# External Links
$(_doc_external("Sys/PetscSubcommDestroy"))
"""
function PetscSubcommDestroy(petsclib::PetscLibType, psubcomm::Union{PetscSubcomm, Ref{PetscSubcomm}}) end

@for_petsc function PetscSubcommDestroy(petsclib::$UnionPetscLib, psubcomm::Union{PetscSubcomm, Ref{PetscSubcomm}} )
	psubcomm_ = psubcomm isa Base.RefValue ? psubcomm : Ref{PetscSubcomm}(psubcomm)

    @chk ccall(
               (:PetscSubcommDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscSubcomm},),
               psubcomm_,
              )


	return nothing
end 

"""
	ccomm::MPI_Comm = PetscSubcommGetChild(petsclib::PetscLibType,scomm::PetscSubcomm) 
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
function PetscSubcommGetChild(petsclib::PetscLibType, scomm::PetscSubcomm) end

@for_petsc function PetscSubcommGetChild(petsclib::$UnionPetscLib, scomm::PetscSubcomm )
	ccomm_ = Ref{MPI.MPI_Comm}()

    @chk ccall(
               (:PetscSubcommGetChild, $petsc_library),
               PetscErrorCode,
               (PetscSubcomm, Ptr{MPI.MPI_Comm}),
               scomm, ccomm_,
              )

	ccomm = MPI.Comm(ccomm_[])

	return ccomm
end 

"""
	pcomm::MPI_Comm = PetscSubcommGetContiguousParent(petsclib::PetscLibType,scomm::PetscSubcomm) 
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
function PetscSubcommGetContiguousParent(petsclib::PetscLibType, scomm::PetscSubcomm) end

@for_petsc function PetscSubcommGetContiguousParent(petsclib::$UnionPetscLib, scomm::PetscSubcomm )
	pcomm_ = Ref{MPI.MPI_Comm}()

    @chk ccall(
               (:PetscSubcommGetContiguousParent, $petsc_library),
               PetscErrorCode,
               (PetscSubcomm, Ptr{MPI.MPI_Comm}),
               scomm, pcomm_,
              )

	pcomm = MPI.Comm(pcomm_[])

	return pcomm
end 

"""
	pcomm::MPI_Comm = PetscSubcommGetParent(petsclib::PetscLibType,scomm::PetscSubcomm) 
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
function PetscSubcommGetParent(petsclib::PetscLibType, scomm::PetscSubcomm) end

@for_petsc function PetscSubcommGetParent(petsclib::$UnionPetscLib, scomm::PetscSubcomm )
	pcomm_ = Ref{MPI.MPI_Comm}()

    @chk ccall(
               (:PetscSubcommGetParent, $petsc_library),
               PetscErrorCode,
               (PetscSubcomm, Ptr{MPI.MPI_Comm}),
               scomm, pcomm_,
              )

	pcomm = MPI.Comm(pcomm_[])

	return pcomm
end 

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

