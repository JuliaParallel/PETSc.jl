# autodefined type arguments for class ------
mutable struct _n_DMPlexStorageVersion end
const DMPlexStorageVersion = Ptr{_n_DMPlexStorageVersion}

mutable struct _n_PetscViewers end
const PetscViewers = Ptr{_n_PetscViewers}
# -------------------------------------------------------

"""
	inviewer::PetscViewer = PetscViewerCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Creates a viewing context. A `PetscViewer` represents a file, a graphical window, a Unix socket or a variety of other ways
of viewing a PETSc object

Collective

Input Parameter:
- `comm` - MPI communicator

Output Parameter:
- `inviewer` - location to put the `PetscViewer` context

Level: advanced

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewerDestroy()`, `PetscViewerSetType()`, `PetscViewerType`

# External Links
$(_doc_external("Sys/PetscViewerCreate"))
"""
function PetscViewerCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function PetscViewerCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	inviewer_ = Ref{PetscViewer}()

    @chk ccall(
               (:PetscViewerCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscViewer}),
               comm, inviewer_,
              )

	inviewer = inviewer_[]

	return inviewer
end 

"""
	PetscViewerRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Adds a viewer to those available for use with `PetscViewerSetType()`

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - name of a new user-defined viewer
- `function` - routine to create method context

Level: developer

-seealso: [](sec_viewers), `PetscViewerRegisterAll()`

# External Links
$(_doc_external("Sys/PetscViewerRegister"))
"""
function PetscViewerRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function PetscViewerRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:PetscViewerRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	mcnt::PetscInt,cnt::PetscInt = PetscViewerFlowControlStart(petsclib::PetscLibType,viewer::PetscViewer) 

# External Links
$(_doc_external("Sys/PetscViewerFlowControlStart"))
"""
function PetscViewerFlowControlStart(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerFlowControlStart(petsclib::$UnionPetscLib, viewer::PetscViewer )
	mcnt_ = Ref{$PetscInt}()
	cnt_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscViewerFlowControlStart, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{$PetscInt}, Ptr{$PetscInt}),
               viewer, mcnt_, cnt_,
              )

	mcnt = mcnt_[]
	cnt = cnt_[]

	return mcnt,cnt
end 

"""
	mcnt::PetscInt = PetscViewerFlowControlStepMain(petsclib::PetscLibType,viewer::PetscViewer, i::PetscInt, cnt::PetscInt) 

# External Links
$(_doc_external("Sys/PetscViewerFlowControlStepMain"))
"""
function PetscViewerFlowControlStepMain(petsclib::PetscLibType, viewer::PetscViewer, i::PetscInt, cnt::PetscInt) end

@for_petsc function PetscViewerFlowControlStepMain(petsclib::$UnionPetscLib, viewer::PetscViewer, i::$PetscInt, cnt::$PetscInt )
	mcnt_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscViewerFlowControlStepMain, $petsc_library),
               PetscErrorCode,
               (PetscViewer, $PetscInt, Ptr{$PetscInt}, $PetscInt),
               viewer, i, mcnt_, cnt,
              )

	mcnt = mcnt_[]

	return mcnt
end 

"""
	mcnt::PetscInt = PetscViewerFlowControlEndMain(petsclib::PetscLibType,viewer::PetscViewer) 

# External Links
$(_doc_external("Sys/PetscViewerFlowControlEndMain"))
"""
function PetscViewerFlowControlEndMain(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerFlowControlEndMain(petsclib::$UnionPetscLib, viewer::PetscViewer )
	mcnt_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscViewerFlowControlEndMain, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{$PetscInt}),
               viewer, mcnt_,
              )

	mcnt = mcnt_[]

	return mcnt
end 

"""
	mcnt::PetscInt = PetscViewerFlowControlStepWorker(petsclib::PetscLibType,viewer::PetscViewer, rank::PetscMPIInt) 

# External Links
$(_doc_external("Sys/PetscViewerFlowControlStepWorker"))
"""
function PetscViewerFlowControlStepWorker(petsclib::PetscLibType, viewer::PetscViewer, rank::PetscMPIInt) end

@for_petsc function PetscViewerFlowControlStepWorker(petsclib::$UnionPetscLib, viewer::PetscViewer, rank::PetscMPIInt )
	mcnt_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscViewerFlowControlStepWorker, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscMPIInt, Ptr{$PetscInt}),
               viewer, rank, mcnt_,
              )

	mcnt = mcnt_[]

	return mcnt
end 

"""
	mcnt::PetscInt = PetscViewerFlowControlEndWorker(petsclib::PetscLibType,viewer::PetscViewer) 

# External Links
$(_doc_external("Sys/PetscViewerFlowControlEndWorker"))
"""
function PetscViewerFlowControlEndWorker(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerFlowControlEndWorker(petsclib::$UnionPetscLib, viewer::PetscViewer )
	mcnt_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscViewerFlowControlEndWorker, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{$PetscInt}),
               viewer, mcnt_,
              )

	mcnt = mcnt_[]

	return mcnt
end 

"""
	PetscViewerFlush(petsclib::PetscLibType,viewer::PetscViewer) 
Flushes a `PetscViewer` (i.e. tries to dump all the
data that has been printed through a `PetscViewer`).

Collective

Input Parameter:
- `viewer` - the `PetscViewer` to be flushed

Level: intermediate

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewerWriteable()`, `PetscViewerSocketOpen()`, `PetscViewerASCIIOpen()`, `PetscViewerDrawOpen()`, `PetscViewerCreate()`, `PetscViewerDestroy()`,
`PetscViewerSetType()`

# External Links
$(_doc_external("Sys/PetscViewerFlush"))
"""
function PetscViewerFlush(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerFlush(petsclib::$UnionPetscLib, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerFlush, $petsc_library),
               PetscErrorCode,
               (PetscViewer,),
               viewer,
              )


	return nothing
end 

"""
	PetscViewerPushFormat(petsclib::PetscLibType,viewer::PetscViewer, format::PetscViewerFormat) 
Sets the format for a `PetscViewer`.

Logically Collective

Input Parameters:
- `viewer` - the `PetscViewer`
- `format` - the format

Level: intermediate

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewerFormat`, `PetscViewerASCIIOpen()`, `PetscViewerBinaryOpen()`, `MatView()`, `VecView()`,
`PetscViewerSetFormat()`, `PetscViewerPopFormat()`

# External Links
$(_doc_external("Sys/PetscViewerPushFormat"))
"""
function PetscViewerPushFormat(petsclib::PetscLibType, viewer::PetscViewer, format::PetscViewerFormat) end

@for_petsc function PetscViewerPushFormat(petsclib::$UnionPetscLib, viewer::PetscViewer, format::PetscViewerFormat )

    @chk ccall(
               (:PetscViewerPushFormat, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscViewerFormat),
               viewer, format,
              )


	return nothing
end 

"""
	PetscViewerPopFormat(petsclib::PetscLibType,viewer::PetscViewer) 
Resets the format for a `PetscViewer` to the value it had before the previous call to `PetscViewerPushFormat()`

Logically Collective

Input Parameter:
- `viewer` - the `PetscViewer`

Level: intermediate

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewerFormat`, `PetscViewerASCIIOpen()`, `PetscViewerBinaryOpen()`, `MatView()`, `VecView()`,
`PetscViewerSetFormat()`, `PetscViewerPushFormat()`

# External Links
$(_doc_external("Sys/PetscViewerPopFormat"))
"""
function PetscViewerPopFormat(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerPopFormat(petsclib::$UnionPetscLib, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerPopFormat, $petsc_library),
               PetscErrorCode,
               (PetscViewer,),
               viewer,
              )


	return nothing
end 

"""
	PetscViewerGetFormat(petsclib::PetscLibType,viewer::PetscViewer, format::PetscViewerFormat) 
Gets the current format for `PetscViewer`.

Not Collective

Input Parameter:
- `viewer` - the `PetscViewer`

Output Parameter:
- `format` - the format

Level: intermediate

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewerASCIIOpen()`, `PetscViewerBinaryOpen()`, `MatView()`, `VecView()`, `PetscViewerType`,
`PetscViewerPushFormat()`, `PetscViewerPopFormat()`, `PetscViewerDrawOpen()`, `PetscViewerSocketOpen()`

# External Links
$(_doc_external("Sys/PetscViewerGetFormat"))
"""
function PetscViewerGetFormat(petsclib::PetscLibType, viewer::PetscViewer, format::PetscViewerFormat) end

@for_petsc function PetscViewerGetFormat(petsclib::$UnionPetscLib, viewer::PetscViewer, format::PetscViewerFormat )

    @chk ccall(
               (:PetscViewerGetFormat, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{PetscViewerFormat}),
               viewer, format,
              )


	return nothing
end 

"""
	PetscViewerGetSubViewer(petsclib::PetscLibType,viewer::PetscViewer, comm::MPI_Comm, outviewer::PetscViewer) 
Creates a new `PetscViewer` (same type as the old)
that lives on a subcommunicator of the original viewer's communicator

Collective

Input Parameters:
- `viewer` - the `PetscViewer` to be reproduced
- `comm`   - the sub communicator to use

Output Parameter:
- `outviewer` - new `PetscViewer`

Level: advanced

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewerSocketOpen()`, `PetscViewerASCIIOpen()`, `PetscViewerDrawOpen()`,
`PetscViewerFlush()`, `PetscViewerRestoreSubViewer()`

# External Links
$(_doc_external("Sys/PetscViewerGetSubViewer"))
"""
function PetscViewerGetSubViewer(petsclib::PetscLibType, viewer::PetscViewer, comm::MPI_Comm, outviewer::PetscViewer) end

@for_petsc function PetscViewerGetSubViewer(petsclib::$UnionPetscLib, viewer::PetscViewer, comm::MPI_Comm, outviewer::PetscViewer )

    @chk ccall(
               (:PetscViewerGetSubViewer, $petsc_library),
               PetscErrorCode,
               (PetscViewer, MPI_Comm, Ptr{PetscViewer}),
               viewer, comm, outviewer,
              )


	return nothing
end 

"""
	PetscViewerRestoreSubViewer(petsclib::PetscLibType,viewer::PetscViewer, comm::MPI_Comm, outviewer::PetscViewer) 
Restores a  `PetscViewer` obtained with `PetscViewerGetSubViewer()`.

Collective

Input Parameters:
- `viewer`    - the `PetscViewer` that was reproduced
- `comm`      - the sub communicator
- `outviewer` - the subviewer to be returned

Level: advanced

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewerSocketOpen()`, `PetscViewerASCIIOpen()`, `PetscViewerDrawOpen()`, `PetscViewerGetSubViewer()`,
`PetscViewerFlush()`

# External Links
$(_doc_external("Sys/PetscViewerRestoreSubViewer"))
"""
function PetscViewerRestoreSubViewer(petsclib::PetscLibType, viewer::PetscViewer, comm::MPI_Comm, outviewer::PetscViewer) end

@for_petsc function PetscViewerRestoreSubViewer(petsclib::$UnionPetscLib, viewer::PetscViewer, comm::MPI_Comm, outviewer::PetscViewer )

    @chk ccall(
               (:PetscViewerRestoreSubViewer, $petsc_library),
               PetscErrorCode,
               (PetscViewer, MPI_Comm, Ptr{PetscViewer}),
               viewer, comm, outviewer,
              )


	return nothing
end 

"""
	PetscViewerFinalizePackage(petsclib::PetscLibType) 
This function destroys any global objects created in PETSc viewers. It is
called from `PetscFinalize()`.

Level: developer

-seealso: [](sec_viewers), `PetscViewer`, `PetscFinalize()`, `PetscViewerInitializePackage()`

# External Links
$(_doc_external("Sys/PetscViewerFinalizePackage"))
"""
function PetscViewerFinalizePackage(petsclib::PetscLibType) end

@for_petsc function PetscViewerFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscViewerFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscViewerInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `PetscViewer` package.

Level: developer

-seealso: [](sec_viewers), `PetscViewer`, `PetscInitialize()`, `PetscViewerFinalizePackage()`

# External Links
$(_doc_external("Sys/PetscViewerInitializePackage"))
"""
function PetscViewerInitializePackage(petsclib::PetscLibType) end

@for_petsc function PetscViewerInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscViewerInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscViewerDestroy(petsclib::PetscLibType,viewer::PetscViewer) 
Destroys a `PetscViewer`.

Collective

Input Parameter:
- `viewer` - the `PetscViewer` to be destroyed.

Level: beginner

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewerCreate()`, `PetscViewerSocketOpen()`, `PetscViewerASCIIOpen()`, `PetscViewerDrawOpen()`

# External Links
$(_doc_external("Sys/PetscViewerDestroy"))
"""
function PetscViewerDestroy(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerDestroy(petsclib::$UnionPetscLib, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscViewer},),
               viewer,
              )


	return nothing
end 

"""
	vf::PetscViewerAndFormat = PetscViewerAndFormatCreate(petsclib::PetscLibType,viewer::PetscViewer, format::PetscViewerFormat) 
Creates a `PetscViewerAndFormat` struct.

Collective

Input Parameters:
- `viewer` - the viewer
- `format` - the format

Output Parameter:
- `vf` - viewer and format object

Level: developer

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewerAndFormat`, `PetscViewerFormat`, `PetscViewerSocketOpen()`, `PetscViewerASCIIOpen()`, `PetscViewerCreate()`,
`PetscViewerDrawOpen()`, `PetscViewerAndFormatDestroy()`

# External Links
$(_doc_external("Sys/PetscViewerAndFormatCreate"))
"""
function PetscViewerAndFormatCreate(petsclib::PetscLibType, viewer::PetscViewer, format::PetscViewerFormat) end

@for_petsc function PetscViewerAndFormatCreate(petsclib::$UnionPetscLib, viewer::PetscViewer, format::PetscViewerFormat )
	vf_ = Ref{PetscViewerAndFormat}()

    @chk ccall(
               (:PetscViewerAndFormatCreate, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscViewerFormat, PetscViewerAndFormat),
               viewer, format, vf_,
              )

	vf = vf_[]

	return vf
end 

"""
	PetscViewerAndFormatDestroy(petsclib::PetscLibType,vf::PetscViewerAndFormat) 
Destroys a `PetscViewerAndFormat` struct created with `PetscViewerAndFormatCreate()`

Collective

Input Parameter:
- `vf` - the `PetscViewerAndFormat` to be destroyed.

Level: developer

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewerAndFormat`, `PetscViewerFormat`, `PetscViewerAndFormatCreate()`, `PetscViewerSocketOpen()`,
`PetscViewerASCIIOpen()`, `PetscViewerCreate()`, `PetscViewerDrawOpen()`

# External Links
$(_doc_external("Sys/PetscViewerAndFormatDestroy"))
"""
function PetscViewerAndFormatDestroy(petsclib::PetscLibType, vf::PetscViewerAndFormat) end

@for_petsc function PetscViewerAndFormatDestroy(petsclib::$UnionPetscLib, vf::PetscViewerAndFormat )

    @chk ccall(
               (:PetscViewerAndFormatDestroy, $petsc_library),
               PetscErrorCode,
               (PetscViewerAndFormat,),
               vf,
              )


	return nothing
end 

"""
	type::PetscViewerType = PetscViewerGetType(petsclib::PetscLibType,viewer::PetscViewer) 
Returns the type of a `PetscViewer`.

Not Collective

Input Parameter:
- `viewer` - the `PetscViewer`

Output Parameter:
- `type` - `PetscViewerType`

Level: intermediate

-seealso: [](sec_viewers), `PetscViewerType`, `PetscViewer`, `PetscViewerCreate()`, `PetscViewerSetType()`

# External Links
$(_doc_external("Sys/PetscViewerGetType"))
"""
function PetscViewerGetType(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerGetType(petsclib::$UnionPetscLib, viewer::PetscViewer )
	type_ = Ref{PetscViewerType}()

    @chk ccall(
               (:PetscViewerGetType, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{PetscViewerType}),
               viewer, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	PetscViewerAppendOptionsPrefix(petsclib::PetscLibType,viewer::PetscViewer, prefix::String) 
Appends to the prefix used for searching for
`PetscViewer` options in the database during `PetscViewerSetFromOptions()`.

Logically Collective

Input Parameters:
- `viewer` - the `PetscViewer` context
- `prefix` - the prefix to prepend to all option names

Level: advanced

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewerGetOptionsPrefix()`, `PetscViewerSetOptionsPrefix()`

# External Links
$(_doc_external("Sys/PetscViewerAppendOptionsPrefix"))
"""
function PetscViewerAppendOptionsPrefix(petsclib::PetscLibType, viewer::PetscViewer, prefix::String) end

@for_petsc function PetscViewerAppendOptionsPrefix(petsclib::$UnionPetscLib, viewer::PetscViewer, prefix::String )

    @chk ccall(
               (:PetscViewerAppendOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cchar}),
               viewer, prefix,
              )


	return nothing
end 

"""
	PetscViewerGetOptionsPrefix(petsclib::PetscLibType,viewer::PetscViewer, prefix::String) 
Gets the prefix used for searching for
`PetscViewer` options in the database during `PetscViewerSetFromOptions()`.

Not Collective

Input Parameter:
- `viewer` - the `PetscViewer` context

Output Parameter:
- `prefix` - pointer to the prefix string used

Level: advanced

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewerAppendOptionsPrefix()`, `PetscViewerSetOptionsPrefix()`

# External Links
$(_doc_external("Sys/PetscViewerGetOptionsPrefix"))
"""
function PetscViewerGetOptionsPrefix(petsclib::PetscLibType, viewer::PetscViewer, prefix::String) end

@for_petsc function PetscViewerGetOptionsPrefix(petsclib::$UnionPetscLib, viewer::PetscViewer, prefix::String )
	prefix_ = Ref(pointer(prefix))

    @chk ccall(
               (:PetscViewerGetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Ptr{Cchar}}),
               viewer, prefix_,
              )


	return nothing
end 

"""
	PetscViewerViewFromOptions(petsclib::PetscLibType,A::PetscViewer, obj::PetscObject, name::String) 
View from the viewer based on options in the options database

Collective

Input Parameters:
- `A`    - the `PetscViewer` context
- `obj`  - Optional object that provides the prefix for the option names
- `name` - command line option

Level: intermediate

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewerView`, `PetscObjectViewFromOptions()`, `PetscViewerCreate()`

# External Links
$(_doc_external("Sys/PetscViewerViewFromOptions"))
"""
function PetscViewerViewFromOptions(petsclib::PetscLibType, A::PetscViewer, obj::PetscObject, name::String) end

@for_petsc function PetscViewerViewFromOptions(petsclib::$UnionPetscLib, A::PetscViewer, obj::PetscObject, name::String )

    @chk ccall(
               (:PetscViewerViewFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

"""
	PetscViewerView(petsclib::PetscLibType,v::PetscViewer, viewer::PetscViewer) 
Visualizes a viewer object.

Collective

Input Parameters:
- `v`      - the viewer to be viewed
- `viewer` - visualization context

Level: beginner

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewerPushFormat()`, `PetscViewerASCIIOpen()`, `PetscViewerDrawOpen()`,
`PetscViewerSocketOpen()`, `PetscViewerBinaryOpen()`, `PetscViewerLoad()`

# External Links
$(_doc_external("Sys/PetscViewerView"))
"""
function PetscViewerView(petsclib::PetscLibType, v::PetscViewer, viewer::PetscViewer) end

@for_petsc function PetscViewerView(petsclib::$UnionPetscLib, v::PetscViewer, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerView, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscViewer),
               v, viewer,
              )


	return nothing
end 

"""
	count::PetscInt = PetscViewerRead(petsclib::PetscLibType,viewer::PetscViewer, data::Cvoid, num::PetscInt, dtype::PetscDataType) 
Reads data from a `PetscViewer`

Collective

Input Parameters:
- `viewer` - The viewer
- `data`   - Location to write the data, treated as an array of the type defined by `datatype`
- `num`    - Number of items of data to read
- `dtype`  - Type of data to read

Output Parameter:
- `count` - number of items of data actually read, or `NULL`

Level: beginner

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewerASCIIOpen()`, `PetscViewerPushFormat()`, `PetscViewerDestroy()`,
`PetscViewerReadable()`, `PetscViewerBinaryGetDescriptor()`,
`PetscViewerBinaryGetInfoPointer()`, `PetscFileMode`

# External Links
$(_doc_external("Sys/PetscViewerRead"))
"""
function PetscViewerRead(petsclib::PetscLibType, viewer::PetscViewer, data::Cvoid, num::PetscInt, dtype::PetscDataType) end

@for_petsc function PetscViewerRead(petsclib::$UnionPetscLib, viewer::PetscViewer, data::Cvoid, num::$PetscInt, dtype::PetscDataType )
	count_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscViewerRead, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cvoid}, $PetscInt, Ptr{$PetscInt}, PetscDataType),
               viewer, data, num, count_, dtype,
              )

	count = count_[]

	return count
end 

"""
	flg::PetscBool = PetscViewerReadable(petsclib::PetscLibType,viewer::PetscViewer) 
Return a flag whether the viewer can be read from with `PetscViewerRead()`

Not Collective

Input Parameter:
- `viewer` - the `PetscViewer` context

Output Parameter:
- `flg` - `PETSC_TRUE` if the viewer is readable, `PETSC_FALSE` otherwise

Level: intermediate

-seealso: [](sec_viewers), `PetscViewerRead()`, `PetscViewer`, `PetscViewerWritable()`, `PetscViewerCheckReadable()`, `PetscViewerCreate()`, `PetscViewerFileSetMode()`, `PetscViewerFileSetType()`

# External Links
$(_doc_external("Sys/PetscViewerReadable"))
"""
function PetscViewerReadable(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerReadable(petsclib::$UnionPetscLib, viewer::PetscViewer )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscViewerReadable, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{PetscBool}),
               viewer, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	flg::PetscBool = PetscViewerWritable(petsclib::PetscLibType,viewer::PetscViewer) 
Return a flag whether the viewer can be written to with `PetscViewerWrite()`

Not Collective

Input Parameter:
- `viewer` - the `PetscViewer` context

Output Parameter:
- `flg` - `PETSC_TRUE` if the viewer is writable, `PETSC_FALSE` otherwise

Level: intermediate

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewerReadable()`, `PetscViewerCheckWritable()`, `PetscViewerCreate()`, `PetscViewerFileSetMode()`, `PetscViewerFileSetType()`

# External Links
$(_doc_external("Sys/PetscViewerWritable"))
"""
function PetscViewerWritable(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerWritable(petsclib::$UnionPetscLib, viewer::PetscViewer )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscViewerWritable, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{PetscBool}),
               viewer, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	PetscViewerCheckReadable(petsclib::PetscLibType,viewer::PetscViewer) 
Check whether the viewer can be read from, generates an error if not

Collective

Input Parameter:
- `viewer` - the `PetscViewer` context

Level: intermediate

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewerReadable()`, `PetscViewerCheckWritable()`, `PetscViewerCreate()`, `PetscViewerFileSetMode()`, `PetscViewerFileSetType()`

# External Links
$(_doc_external("Sys/PetscViewerCheckReadable"))
"""
function PetscViewerCheckReadable(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerCheckReadable(petsclib::$UnionPetscLib, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerCheckReadable, $petsc_library),
               PetscErrorCode,
               (PetscViewer,),
               viewer,
              )


	return nothing
end 

"""
	PetscViewerCheckWritable(petsclib::PetscLibType,viewer::PetscViewer) 
Check whether the viewer can be written to, generates an error if not

Collective

Input Parameter:
- `viewer` - the `PetscViewer` context

Level: intermediate

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewerWritable()`, `PetscViewerCheckReadable()`, `PetscViewerCreate()`, `PetscViewerFileSetMode()`, `PetscViewerFileSetType()`

# External Links
$(_doc_external("Sys/PetscViewerCheckWritable"))
"""
function PetscViewerCheckWritable(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerCheckWritable(petsclib::$UnionPetscLib, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerCheckWritable, $petsc_library),
               PetscErrorCode,
               (PetscViewer,),
               viewer,
              )


	return nothing
end 

"""
	PetscViewerASCIIGetPointer(petsclib::PetscLibType,viewer::PetscViewer, fd::Libc.FILE) 
Extracts the file pointer from an ASCII `PetscViewer`.

Not Collective, depending on the viewer the value may be meaningless except for process 0 of the viewer; No Fortran Support

Input Parameter:
- `viewer` - `PetscViewer` context, obtained from `PetscViewerASCIIOpen()`

Output Parameter:
- `fd` - file pointer

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERASCII`, `PetscViewerASCIIOpen()`, `PetscViewerDestroy()`, `PetscViewerSetType()`,
`PetscViewerCreate()`, `PetscViewerASCIIPrintf()`, `PetscViewerASCIISynchronizedPrintf()`, `PetscViewerFlush()`

# External Links
$(_doc_external("Sys/PetscViewerASCIIGetPointer"))
"""
function PetscViewerASCIIGetPointer(petsclib::PetscLibType, viewer::PetscViewer, fd::Libc.FILE) end

@for_petsc function PetscViewerASCIIGetPointer(petsclib::$UnionPetscLib, viewer::PetscViewer, fd::Libc.FILE )

    @chk ccall(
               (:PetscViewerASCIIGetPointer, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Libc.FILE),
               viewer, fd,
              )


	return nothing
end 

"""
	PetscViewerASCIISetTab(petsclib::PetscLibType,viewer::PetscViewer, tabs::PetscInt) 
Causes `PetscViewer` to tab in a number of times before printing

Not Collective, but only first processor in set has any effect; No Fortran Support

Input Parameters:
- `viewer` - obtained with `PetscViewerASCIIOpen()`
- `tabs`   - number of tabs

Level: developer

-seealso: [](sec_viewers), `PETSCVIEWERASCII`, `PetscPrintf()`, `PetscSynchronizedPrintf()`, `PetscViewerASCIIPrintf()`,
`PetscViewerASCIIGetTab()`,
`PetscViewerASCIIPopTab()`, `PetscViewerASCIISynchronizedPrintf()`, `PetscViewerASCIIOpen()`,
`PetscViewerCreate()`, `PetscViewerDestroy()`, `PetscViewerSetType()`, `PetscViewerASCIIGetPointer()`,
`PetscViewerASCIIPushTab()`

# External Links
$(_doc_external("Sys/PetscViewerASCIISetTab"))
"""
function PetscViewerASCIISetTab(petsclib::PetscLibType, viewer::PetscViewer, tabs::PetscInt) end

@for_petsc function PetscViewerASCIISetTab(petsclib::$UnionPetscLib, viewer::PetscViewer, tabs::$PetscInt )

    @chk ccall(
               (:PetscViewerASCIISetTab, $petsc_library),
               PetscErrorCode,
               (PetscViewer, $PetscInt),
               viewer, tabs,
              )


	return nothing
end 

"""
	tabs::PetscInt = PetscViewerASCIIGetTab(petsclib::PetscLibType,viewer::PetscViewer) 
Return the number of tabs used by `PetscViewer`.

Not Collective, meaningful on first processor only; No Fortran Support

Input Parameter:
- `viewer` - obtained with `PetscViewerASCIIOpen()`

Output Parameter:
- `tabs` - number of tabs

Level: developer

-seealso: [](sec_viewers), `PETSCVIEWERASCII`, `PetscPrintf()`, `PetscSynchronizedPrintf()`, `PetscViewerASCIIPrintf()`,
`PetscViewerASCIISetTab()`,
`PetscViewerASCIIPopTab()`, `PetscViewerASCIISynchronizedPrintf()`, `PetscViewerASCIIOpen()`,
`PetscViewerCreate()`, `PetscViewerDestroy()`, `PetscViewerSetType()`, `PetscViewerASCIIGetPointer()`, `PetscViewerASCIIPushTab()`

# External Links
$(_doc_external("Sys/PetscViewerASCIIGetTab"))
"""
function PetscViewerASCIIGetTab(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerASCIIGetTab(petsclib::$UnionPetscLib, viewer::PetscViewer )
	tabs_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscViewerASCIIGetTab, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{$PetscInt}),
               viewer, tabs_,
              )

	tabs = tabs_[]

	return tabs
end 

"""
	PetscViewerASCIIAddTab(petsclib::PetscLibType,viewer::PetscViewer, tabs::PetscInt) 
Add to the number of times a `PETSCVIEWERASCII` viewer tabs before printing

Not Collective, but only first processor in set has any effect; No Fortran Support

Input Parameters:
- `viewer` - obtained with `PetscViewerASCIIOpen()`
- `tabs`   - number of tabs

Level: developer

-seealso: [](sec_viewers), `PETSCVIEWERASCII`, `PetscPrintf()`, `PetscSynchronizedPrintf()`, `PetscViewerASCIIPrintf()`,
`PetscViewerASCIIPopTab()`, `PetscViewerASCIISynchronizedPrintf()`, `PetscViewerASCIIOpen()`,
`PetscViewerCreate()`, `PetscViewerDestroy()`, `PetscViewerSetType()`, `PetscViewerASCIIGetPointer()`, `PetscViewerASCIIPushTab()`

# External Links
$(_doc_external("Sys/PetscViewerASCIIAddTab"))
"""
function PetscViewerASCIIAddTab(petsclib::PetscLibType, viewer::PetscViewer, tabs::PetscInt) end

@for_petsc function PetscViewerASCIIAddTab(petsclib::$UnionPetscLib, viewer::PetscViewer, tabs::$PetscInt )

    @chk ccall(
               (:PetscViewerASCIIAddTab, $petsc_library),
               PetscErrorCode,
               (PetscViewer, $PetscInt),
               viewer, tabs,
              )


	return nothing
end 

"""
	PetscViewerASCIISubtractTab(petsclib::PetscLibType,viewer::PetscViewer, tabs::PetscInt) 
Subtracts from the number of times a `PETSCVIEWERASCII` viewer tabs before printing

Not Collective, but only first processor in set has any effect; No Fortran Support

Input Parameters:
- `viewer` - obtained with `PetscViewerASCIIOpen()`
- `tabs`   - number of tabs

Level: developer

-seealso: [](sec_viewers), `PETSCVIEWERASCII`, `PetscPrintf()`, `PetscSynchronizedPrintf()`, `PetscViewerASCIIPrintf()`,
`PetscViewerASCIIPopTab()`, `PetscViewerASCIISynchronizedPrintf()`, `PetscViewerASCIIOpen()`,
`PetscViewerCreate()`, `PetscViewerDestroy()`, `PetscViewerSetType()`, `PetscViewerASCIIGetPointer()`,
`PetscViewerASCIIPushTab()`

# External Links
$(_doc_external("Sys/PetscViewerASCIISubtractTab"))
"""
function PetscViewerASCIISubtractTab(petsclib::PetscLibType, viewer::PetscViewer, tabs::PetscInt) end

@for_petsc function PetscViewerASCIISubtractTab(petsclib::$UnionPetscLib, viewer::PetscViewer, tabs::$PetscInt )

    @chk ccall(
               (:PetscViewerASCIISubtractTab, $petsc_library),
               PetscErrorCode,
               (PetscViewer, $PetscInt),
               viewer, tabs,
              )


	return nothing
end 

"""
	PetscViewerASCIIPushSynchronized(petsclib::PetscLibType,viewer::PetscViewer) 
Allows calls to `PetscViewerASCIISynchronizedPrintf()` for this viewer

Collective

Input Parameter:
- `viewer` - obtained with `PetscViewerASCIIOpen()`

Level: intermediate

-seealso: [](sec_viewers), `PetscViewerASCIISynchronizedPrintf()`, `PetscViewerFlush()`, `PetscViewerASCIIPopSynchronized()`,
`PetscSynchronizedPrintf()`, `PetscViewerASCIIPrintf()`, `PetscViewerASCIIOpen()`,
`PetscViewerCreate()`, `PetscViewerDestroy()`, `PetscViewerSetType()`

# External Links
$(_doc_external("Sys/PetscViewerASCIIPushSynchronized"))
"""
function PetscViewerASCIIPushSynchronized(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerASCIIPushSynchronized(petsclib::$UnionPetscLib, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerASCIIPushSynchronized, $petsc_library),
               PetscErrorCode,
               (PetscViewer,),
               viewer,
              )


	return nothing
end 

"""
	PetscViewerASCIIPopSynchronized(petsclib::PetscLibType,viewer::PetscViewer) 
Undoes most recent `PetscViewerASCIIPushSynchronized()` for this viewer

Collective

Input Parameter:
- `viewer` - obtained with `PetscViewerASCIIOpen()`

Level: intermediate

-seealso: [](sec_viewers), `PetscViewerASCIIPushSynchronized()`, `PetscViewerASCIISynchronizedPrintf()`, `PetscViewerFlush()`,
`PetscSynchronizedPrintf()`, `PetscViewerASCIIPrintf()`, `PetscViewerASCIIOpen()`,
`PetscViewerCreate()`, `PetscViewerDestroy()`, `PetscViewerSetType()`

# External Links
$(_doc_external("Sys/PetscViewerASCIIPopSynchronized"))
"""
function PetscViewerASCIIPopSynchronized(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerASCIIPopSynchronized(petsclib::$UnionPetscLib, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerASCIIPopSynchronized, $petsc_library),
               PetscErrorCode,
               (PetscViewer,),
               viewer,
              )


	return nothing
end 

"""
	PetscViewerASCIIPushTab(petsclib::PetscLibType,viewer::PetscViewer) 
Adds one more tab to the amount that `PetscViewerASCIIPrintf()`
lines are tabbed.

Not Collective, but only first MPI rank in the viewer has any effect; No Fortran Support

Input Parameter:
- `viewer` - obtained with `PetscViewerASCIIOpen()`

Level: developer

-seealso: [](sec_viewers), `PetscPrintf()`, `PetscSynchronizedPrintf()`, `PetscViewerASCIIPrintf()`,
`PetscViewerASCIIPopTab()`, `PetscViewerASCIISynchronizedPrintf()`, `PetscViewerASCIIOpen()`,
`PetscViewerCreate()`, `PetscViewerDestroy()`, `PetscViewerSetType()`, `PetscViewerASCIIGetPointer()`

# External Links
$(_doc_external("Sys/PetscViewerASCIIPushTab"))
"""
function PetscViewerASCIIPushTab(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerASCIIPushTab(petsclib::$UnionPetscLib, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerASCIIPushTab, $petsc_library),
               PetscErrorCode,
               (PetscViewer,),
               viewer,
              )


	return nothing
end 

"""
	PetscViewerASCIIPopTab(petsclib::PetscLibType,viewer::PetscViewer) 
Removes one tab from the amount that `PetscViewerASCIIPrintf()` lines are tabbed that was provided by
`PetscViewerASCIIPushTab()`

Not Collective, but only first MPI rank in the viewer has any effect; No Fortran Support

Input Parameter:
- `viewer` - obtained with `PetscViewerASCIIOpen()`

Level: developer

-seealso: [](sec_viewers), `PetscPrintf()`, `PetscSynchronizedPrintf()`, `PetscViewerASCIIPrintf()`,
`PetscViewerASCIIPushTab()`, `PetscViewerASCIISynchronizedPrintf()`, `PetscViewerASCIIOpen()`,
`PetscViewerCreate()`, `PetscViewerDestroy()`, `PetscViewerSetType()`, `PetscViewerASCIIGetPointer()`

# External Links
$(_doc_external("Sys/PetscViewerASCIIPopTab"))
"""
function PetscViewerASCIIPopTab(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerASCIIPopTab(petsclib::$UnionPetscLib, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerASCIIPopTab, $petsc_library),
               PetscErrorCode,
               (PetscViewer,),
               viewer,
              )


	return nothing
end 

"""
	PetscViewerASCIIUseTabs(petsclib::PetscLibType,viewer::PetscViewer, flg::PetscBool) 
Turns on or off the use of tabs with the `PETSCVIEWERASCII` `PetscViewer`

Not Collective, but only first MPI rank in the viewer has any effect; No Fortran Support

Input Parameters:
- `viewer` - obtained with `PetscViewerASCIIOpen()`
- `flg`    - `PETSC_TRUE` or `PETSC_FALSE`

Level: developer

-seealso: [](sec_viewers), `PetscPrintf()`, `PetscSynchronizedPrintf()`, `PetscViewerASCIIPrintf()`,
`PetscViewerASCIIPopTab()`, `PetscViewerASCIISynchronizedPrintf()`, `PetscViewerASCIIPushTab()`, `PetscViewerASCIIOpen()`,
`PetscViewerCreate()`, `PetscViewerDestroy()`, `PetscViewerSetType()`, `PetscViewerASCIIGetPointer()`

# External Links
$(_doc_external("Sys/PetscViewerASCIIUseTabs"))
"""
function PetscViewerASCIIUseTabs(petsclib::PetscLibType, viewer::PetscViewer, flg::PetscBool) end

@for_petsc function PetscViewerASCIIUseTabs(petsclib::$UnionPetscLib, viewer::PetscViewer, flg::PetscBool )

    @chk ccall(
               (:PetscViewerASCIIUseTabs, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscBool),
               viewer, flg,
              )


	return nothing
end 

"""
	PetscViewerASCIIGetStdout(petsclib::PetscLibType,comm::MPI_Comm, viewer::PetscViewer) 
Creates a `PETSCVIEWERASCII` `PetscViewer` shared by all processes
in a communicator that prints to `stdout`. Error returning version of `PETSC_VIEWER_STDOUT_()`

Collective

Input Parameter:
- `comm` - the MPI communicator to share the `PetscViewer`

Output Parameter:
- `viewer` - the viewer

Level: beginner

-seealso: [](sec_viewers), `PetscViewerASCIIGetStderr()`, `PETSC_VIEWER_DRAW_()`, `PetscViewerASCIIOpen()`, `PETSC_VIEWER_STDERR_`, `PETSC_VIEWER_STDOUT_WORLD`,
`PETSC_VIEWER_STDOUT_SELF`

# External Links
$(_doc_external("Sys/PetscViewerASCIIGetStdout"))
"""
function PetscViewerASCIIGetStdout(petsclib::PetscLibType, comm::MPI_Comm, viewer::PetscViewer) end

@for_petsc function PetscViewerASCIIGetStdout(petsclib::$UnionPetscLib, comm::MPI_Comm, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerASCIIGetStdout, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscViewer}),
               comm, viewer,
              )


	return nothing
end 

"""
	PetscViewerFileSetName(petsclib::PetscLibType,viewer::PetscViewer, name::String) 
Sets the name of the file the `PetscViewer` should use.

Collective

Input Parameters:
- `viewer` - the `PetscViewer`; for example, of type `PETSCVIEWERASCII` or `PETSCVIEWERBINARY`
- `name`   - the name of the file it should use

Level: advanced

-seealso: [](sec_viewers), `PetscViewerCreate()`, `PetscViewerSetType()`, `PetscViewerASCIIOpen()`, `PetscViewerBinaryOpen()`, `PetscViewerDestroy()`,
`PetscViewerASCIIGetPointer()`, `PetscViewerASCIIPrintf()`, `PetscViewerASCIISynchronizedPrintf()`

# External Links
$(_doc_external("Sys/PetscViewerFileSetName"))
"""
function PetscViewerFileSetName(petsclib::PetscLibType, viewer::PetscViewer, name::String) end

@for_petsc function PetscViewerFileSetName(petsclib::$UnionPetscLib, viewer::PetscViewer, name::String )

    @chk ccall(
               (:PetscViewerFileSetName, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cchar}),
               viewer, name,
              )


	return nothing
end 

"""
	PetscViewerFileGetName(petsclib::PetscLibType,viewer::PetscViewer, name::String) 
Gets the name of the file the `PetscViewer` is using

Not Collective

Input Parameter:
- `viewer` - the `PetscViewer`

Output Parameter:
- `name` - the name of the file it is using

Level: advanced

-seealso: [](sec_viewers), `PetscViewerCreate()`, `PetscViewerSetType()`, `PetscViewerASCIIOpen()`, `PetscViewerBinaryOpen()`, `PetscViewerFileSetName()`

# External Links
$(_doc_external("Sys/PetscViewerFileGetName"))
"""
function PetscViewerFileGetName(petsclib::PetscLibType, viewer::PetscViewer, name::String) end

@for_petsc function PetscViewerFileGetName(petsclib::$UnionPetscLib, viewer::PetscViewer, name::String )
	name_ = Ref(pointer(name))

    @chk ccall(
               (:PetscViewerFileGetName, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Ptr{Cchar}}),
               viewer, name_,
              )


	return nothing
end 

"""
	count::PetscInt = PetscViewerASCIIRead(petsclib::PetscLibType,viewer::PetscViewer, data::Cvoid, num::PetscInt, dtype::PetscDataType) 
Reads from a `PETSCVIEWERASCII` file

Only MPI rank 0 in the `PetscViewer` may call this

Input Parameters:
- `viewer` - the `PETSCVIEWERASCII` viewer
- `data`   - location to write the data, treated as an array of type indicated by `datatype`
- `num`    - number of items of data to read
- `dtype`  - type of data to read

Output Parameter:
- `count` - number of items of data actually read, or `NULL`

Level: beginner

-seealso: [](sec_viewers), `PetscViewerASCIIOpen()`, `PetscViewerPushFormat()`, `PetscViewerDestroy()`, `PetscViewerCreate()`, `PetscViewerFileSetMode()`, `PetscViewerFileSetName()`
`VecView()`, `MatView()`, `VecLoad()`, `MatLoad()`, `PetscViewerBinaryGetDescriptor()`,
`PetscViewerBinaryGetInfoPointer()`, `PetscFileMode`, `PetscViewer`, `PetscViewerBinaryRead()`

# External Links
$(_doc_external("Sys/PetscViewerASCIIRead"))
"""
function PetscViewerASCIIRead(petsclib::PetscLibType, viewer::PetscViewer, data::Cvoid, num::PetscInt, dtype::PetscDataType) end

@for_petsc function PetscViewerASCIIRead(petsclib::$UnionPetscLib, viewer::PetscViewer, data::Cvoid, num::$PetscInt, dtype::PetscDataType )
	count_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscViewerASCIIRead, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cvoid}, $PetscInt, Ptr{$PetscInt}, PetscDataType),
               viewer, data, num, count_, dtype,
              )

	count = count_[]

	return count
end 

"""
	PetscViewerASCIIGetStderr(petsclib::PetscLibType,comm::MPI_Comm, viewer::PetscViewer) 
Creates a `PETSCVIEWERASCII` `PetscViewer` shared by all MPI processes
in a communicator that prints to `stderr`. Error returning version of `PETSC_VIEWER_STDERR_()`

Collective

Input Parameter:
- `comm` - the MPI communicator to share the `PetscViewer`

Output Parameter:
- `viewer` - the viewer

Level: beginner

-seealso: [](sec_viewers), `PetscViewerASCIIGetStdout()`, `PETSC_VIEWER_DRAW_()`, `PetscViewerASCIIOpen()`, `PETSC_VIEWER_STDERR_`, `PETSC_VIEWER_STDERR_WORLD`,
`PETSC_VIEWER_STDERR_SELF`

# External Links
$(_doc_external("Sys/PetscViewerASCIIGetStderr"))
"""
function PetscViewerASCIIGetStderr(petsclib::PetscLibType, comm::MPI_Comm, viewer::PetscViewer) end

@for_petsc function PetscViewerASCIIGetStderr(petsclib::$UnionPetscLib, comm::MPI_Comm, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerASCIIGetStderr, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscViewer}),
               comm, viewer,
              )


	return nothing
end 

"""
	PetscViewerASCIIOpen(petsclib::PetscLibType,comm::MPI_Comm, name::String, viewer::PetscViewer) 
Opens an ASCII file for writing as a `PETSCVIEWERASCII` `PetscViewer`.

Collective

Input Parameters:
- `comm` - the communicator
- `name` - the file name

Output Parameter:
- `viewer` - the `PetscViewer` to use with the specified file

Level: beginner

-seealso: [](sec_viewers), `MatView()`, `VecView()`, `PetscViewerDestroy()`, `PetscViewerBinaryOpen()`, `PetscViewerASCIIRead()`, `PETSCVIEWERASCII`
`PetscViewerASCIIGetPointer()`, `PetscViewerPushFormat()`, `PETSC_VIEWER_STDOUT_`, `PETSC_VIEWER_STDERR_`,
`PETSC_VIEWER_STDOUT_WORLD`, `PETSC_VIEWER_STDOUT_SELF`, `PetscViewerASCIIGetStdout()`, `PetscViewerASCIIGetStderr()`

# External Links
$(_doc_external("Sys/PetscViewerASCIIOpen"))
"""
function PetscViewerASCIIOpen(petsclib::PetscLibType, comm::MPI_Comm, name::String, viewer::PetscViewer) end

@for_petsc function PetscViewerASCIIOpen(petsclib::$UnionPetscLib, comm::MPI_Comm, name::String, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerASCIIOpen, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{PetscViewer}),
               comm, name, viewer,
              )


	return nothing
end 

"""
	PetscViewerASCIIOpenWithFILE(petsclib::PetscLibType,comm::MPI_Comm, fd::Libc.FILE, viewer::PetscViewer) 
Given an open file creates an `PETSCVIEWERASCII` viewer that prints to it.

Collective

Input Parameters:
- `comm` - the communicator
- `fd`   - the `FILE` pointer

Output Parameter:
- `viewer` - the `PetscViewer` to use with the specified file

Level: beginner

-seealso: [](sec_viewers), `MatView()`, `VecView()`, `PetscViewerDestroy()`, `PetscViewerBinaryOpen()`, `PetscViewerASCIIOpenWithFileUnit()`,
`PetscViewerASCIIGetPointer()`, `PetscViewerPushFormat()`, `PETSC_VIEWER_STDOUT_`, `PETSC_VIEWER_STDERR_`,
`PETSC_VIEWER_STDOUT_WORLD`, `PETSC_VIEWER_STDOUT_SELF`, `PetscViewerASCIIOpen()`, `PetscViewerASCIISetFILE()`, `PETSCVIEWERASCII`

# External Links
$(_doc_external("Sys/PetscViewerASCIIOpenWithFILE"))
"""
function PetscViewerASCIIOpenWithFILE(petsclib::PetscLibType, comm::MPI_Comm, fd::Libc.FILE, viewer::PetscViewer) end

@for_petsc function PetscViewerASCIIOpenWithFILE(petsclib::$UnionPetscLib, comm::MPI_Comm, fd::Libc.FILE, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerASCIIOpenWithFILE, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Libc.FILE}, Ptr{PetscViewer}),
               comm, fd, viewer,
              )


	return nothing
end 

"""
	PetscViewerASCIISetFILE(petsclib::PetscLibType,viewer::PetscViewer, fd::Libc.FILE) 
Given an open file sets the `PETSCVIEWERASCII` viewer to use the file for output

Not Collective

Input Parameters:
- `viewer` - the `PetscViewer` to use with the specified file
- `fd`     - the `FILE` pointer

Level: beginner

-seealso: `MatView()`, `VecView()`, `PetscViewerDestroy()`, `PetscViewerBinaryOpen()`, `PetscViewerASCIISetFileUnit()`,
`PetscViewerASCIIGetPointer()`, `PetscViewerPushFormat()`, `PETSC_VIEWER_STDOUT_`, `PETSC_VIEWER_STDERR_`,
`PETSC_VIEWER_STDOUT_WORLD`, `PETSC_VIEWER_STDOUT_SELF`, `PetscViewerASCIIOpen()`, `PetscViewerASCIIOpenWithFILE()`, `PETSCVIEWERASCII`

# External Links
$(_doc_external("Sys/PetscViewerASCIISetFILE"))
"""
function PetscViewerASCIISetFILE(petsclib::PetscLibType, viewer::PetscViewer, fd::Libc.FILE) end

@for_petsc function PetscViewerASCIISetFILE(petsclib::$UnionPetscLib, viewer::PetscViewer, fd::Libc.FILE )

    @chk ccall(
               (:PetscViewerASCIISetFILE, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Libc.FILE}),
               viewer, fd,
              )


	return nothing
end 

"""
	PetscViewerVUGetPointer(petsclib::PetscLibType,viewer::PetscViewer, fd::Libc.FILE) 
Extracts the file pointer from a `PETSCVIEWERVU` `PetscViewer`.

Not Collective

Input Parameter:
- `viewer` - The `PetscViewer`

Output Parameter:
- `fd` - The file pointer

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERVU`, `PetscViewerASCIIGetPointer()`

# External Links
$(_doc_external("Sys/PetscViewerVUGetPointer"))
"""
function PetscViewerVUGetPointer(petsclib::PetscLibType, viewer::PetscViewer, fd::Libc.FILE) end

@for_petsc function PetscViewerVUGetPointer(petsclib::$UnionPetscLib, viewer::PetscViewer, fd::Libc.FILE )

    @chk ccall(
               (:PetscViewerVUGetPointer, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Libc.FILE),
               viewer, fd,
              )


	return nothing
end 

"""
	PetscViewerVUSetVecSeen(petsclib::PetscLibType,viewer::PetscViewer, vecSeen::PetscBool) 
Sets the flag which indicates whether we have viewed
a vector. This is usually called internally rather than by a user.

Not Collective

Input Parameters:
- `viewer`  - The `PETSCVIEWERVU` `PetscViewer`
- `vecSeen` - The flag which indicates whether we have viewed a vector

Level: developer

-seealso: [](sec_viewers), `PETSCVIEWERVU`, `PetscViewerVUGetVecSeen()`

# External Links
$(_doc_external("Sys/PetscViewerVUSetVecSeen"))
"""
function PetscViewerVUSetVecSeen(petsclib::PetscLibType, viewer::PetscViewer, vecSeen::PetscBool) end

@for_petsc function PetscViewerVUSetVecSeen(petsclib::$UnionPetscLib, viewer::PetscViewer, vecSeen::PetscBool )

    @chk ccall(
               (:PetscViewerVUSetVecSeen, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscBool),
               viewer, vecSeen,
              )


	return nothing
end 

"""
	vecSeen::PetscBool = PetscViewerVUGetVecSeen(petsclib::PetscLibType,viewer::PetscViewer) 
Gets the flag which indicates whether we have viewed
a vector. This is usually called internally rather than by a user.

Not Collective

Input Parameter:
- `viewer` - The `PETSCVIEWERVU` `PetscViewer`

Output Parameter:
- `vecSeen` - The flag which indicates whether we have viewed a vector

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERVU`

# External Links
$(_doc_external("Sys/PetscViewerVUGetVecSeen"))
"""
function PetscViewerVUGetVecSeen(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerVUGetVecSeen(petsclib::$UnionPetscLib, viewer::PetscViewer )
	vecSeen_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscViewerVUGetVecSeen, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{PetscBool}),
               viewer, vecSeen_,
              )

	vecSeen = vecSeen_[]

	return vecSeen
end 

"""
	PetscViewerVUFlushDeferred(petsclib::PetscLibType,viewer::PetscViewer) 
Flushes the deferred write cache to the file.

Not Collective

Input Parameter:
- `viewer` - The `PETSCVIEWERVU` `PetscViewer`

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERVU`, `PetscViewerVUPrintDeferred()`

# External Links
$(_doc_external("Sys/PetscViewerVUFlushDeferred"))
"""
function PetscViewerVUFlushDeferred(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerVUFlushDeferred(petsclib::$UnionPetscLib, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerVUFlushDeferred, $petsc_library),
               PetscErrorCode,
               (PetscViewer,),
               viewer,
              )


	return nothing
end 

"""
	PetscViewerPythonSetType(petsclib::PetscLibType,viewer::PetscViewer, pyname::String) 
Initialize a `PetscViewer` object implemented in Python.

Collective

Input Parameters:
- `viewer`  - the viewer object.
- `pyname`  - full dotted Python name [package].module[.{class|function}]

Options Database Key:
- `-viewer_python_type <pyname>`  - python class

Level: intermediate

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewerType`, `PetscViewerCreate()`, `PetscViewerSetType()`, `PETSCVIEWERPYTHON`, `PetscPythonInitialize()`

# External Links
$(_doc_external("Sys/PetscViewerPythonSetType"))
"""
function PetscViewerPythonSetType(petsclib::PetscLibType, viewer::PetscViewer, pyname::String) end

@for_petsc function PetscViewerPythonSetType(petsclib::$UnionPetscLib, viewer::PetscViewer, pyname::String )

    @chk ccall(
               (:PetscViewerPythonSetType, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cchar}),
               viewer, pyname,
              )


	return nothing
end 

"""
	pyname::String = PetscViewerPythonGetType(petsclib::PetscLibType,viewer::PetscViewer) 
Get the Python name of a `PetscViewer` object implemented in Python.

Not Collective

Input Parameter:
- `viewer`  - the viewer

Output Parameter:
- `pyname`  - full dotted Python name [package].module[.{class|function}]

Level: intermediate

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewerType`, `PetscViewerCreate()`, `PetscViewerSetType()`, `PETSCVIEWERPYTHON`, `PetscPythonInitialize()`, `PetscViewerPythonSetType()`

# External Links
$(_doc_external("Sys/PetscViewerPythonGetType"))
"""
function PetscViewerPythonGetType(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerPythonGetType(petsclib::$UnionPetscLib, viewer::PetscViewer )
	pyname_ = Ref{Ptr{Cchar}}()

    @chk ccall(
               (:PetscViewerPythonGetType, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Ptr{Cchar}}),
               viewer, pyname_,
              )

	pyname = unsafe_wrap(Array, pyname_[], VecGetLocalSize(petsclib, x); own = false)

	return pyname
end 

"""
	PetscViewerPythonViewObject(petsclib::PetscLibType,viewer::PetscViewer, obj::PetscObject) 
View a `PetscObject`.

Collective

Input Parameters:
- `viewer`  - the viewer object.
- `obj`  - the object to be viewed.

Level: intermediate

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewerPythonCreate()`

# External Links
$(_doc_external("Sys/PetscViewerPythonViewObject"))
"""
function PetscViewerPythonViewObject(petsclib::PetscLibType, viewer::PetscViewer, obj::PetscObject) end

@for_petsc function PetscViewerPythonViewObject(petsclib::$UnionPetscLib, viewer::PetscViewer, obj::PetscObject )

    @chk ccall(
               (:PetscViewerPythonViewObject, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscObject),
               viewer, obj,
              )


	return nothing
end 

"""
	viewer::PetscViewer = PetscViewerPythonCreate(petsclib::PetscLibType,comm::MPI_Comm, pyname::String) 
Create a `PetscViewer` object implemented in Python.

Collective

Input Parameters:
- `comm`  - MPI communicator
- `pyname`  - full dotted Python name [package].module[.{class|function}]

Output Parameter:
- `viewer`  - the viewer

Level: intermediate

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewerType`, `PETSCVIEWERPYTHON`, `PetscViewerPythonSetType()`, `PetscPythonInitialize()`, `PetscViewerPythonViewObject()`

# External Links
$(_doc_external("Sys/PetscViewerPythonCreate"))
"""
function PetscViewerPythonCreate(petsclib::PetscLibType, comm::MPI_Comm, pyname::String) end

@for_petsc function PetscViewerPythonCreate(petsclib::$UnionPetscLib, comm::MPI_Comm, pyname::String )
	viewer_ = Ref{PetscViewer}()

    @chk ccall(
               (:PetscViewerPythonCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{PetscViewer}),
               comm, pyname, viewer_,
              )

	viewer = viewer_[]

	return viewer
end 

"""
	PetscViewerHDF5GetGroup(petsclib::PetscLibType,viewer::PetscViewer, path::String, abspath::String) 
Get the current HDF5 group name (full path), set with `PetscViewerHDF5PushGroup()`/`PetscViewerHDF5PopGroup()`.

Not Collective

Input Parameters:
- `viewer` - the `PetscViewer` of type `PETSCVIEWERHDF5`
- `path`   - (Optional) The path relative to the pushed group

Output Parameter:
- `abspath` - The absolute HDF5 path (group)

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerHDF5Open()`, `PetscViewerHDF5PushGroup()`, `PetscViewerHDF5PopGroup()`, `PetscViewerHDF5OpenGroup()`, `PetscViewerHDF5WriteGroup()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5GetGroup"))
"""
function PetscViewerHDF5GetGroup(petsclib::PetscLibType, viewer::PetscViewer, path::String, abspath::String) end

@for_petsc function PetscViewerHDF5GetGroup(petsclib::$UnionPetscLib, viewer::PetscViewer, path::String, abspath::String )
	abspath_ = Ref(pointer(abspath))

    @chk ccall(
               (:PetscViewerHDF5GetGroup, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cchar}, Ptr{Ptr{Cchar}}),
               viewer, path, abspath_,
              )


	return nothing
end 

"""
	PetscViewerHDF5SetBaseDimension2(petsclib::PetscLibType,viewer::PetscViewer, flg::PetscBool) 
Vectors of 1 dimension (i.e. bs/dof is 1) will be saved in the HDF5 file with a
dimension of 2.

Logically Collective

Input Parameters:
- `viewer` - the `PetscViewer`; if it is a `PETSCVIEWERHDF5` then this command is ignored
- `flg`    - if `PETSC_TRUE` the vector will always have at least a dimension of 2 even if that first dimension is of size 1

Options Database Key:
- `-viewer_hdf5_base_dimension2` - turns on (true) or off (false) using a dimension of 2 in the HDF5 file even if the bs/dof of the vector is 1

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerFileSetMode()`, `PetscViewerCreate()`, `PetscViewerSetType()`, `PetscViewerBinaryOpen()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5SetBaseDimension2"))
"""
function PetscViewerHDF5SetBaseDimension2(petsclib::PetscLibType, viewer::PetscViewer, flg::PetscBool) end

@for_petsc function PetscViewerHDF5SetBaseDimension2(petsclib::$UnionPetscLib, viewer::PetscViewer, flg::PetscBool )

    @chk ccall(
               (:PetscViewerHDF5SetBaseDimension2, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscBool),
               viewer, flg,
              )


	return nothing
end 

"""
	flg::PetscBool = PetscViewerHDF5GetBaseDimension2(petsclib::PetscLibType,viewer::PetscViewer) 
Vectors of 1 dimension (i.e. bs/dof is 1) will be saved in the HDF5 file with a
dimension of 2.

Logically Collective

Input Parameter:
- `viewer` - the `PetscViewer`, must be `PETSCVIEWERHDF5`

Output Parameter:
- `flg` - if `PETSC_TRUE` the vector will always have at least a dimension of 2 even if that first dimension is of size 1

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerFileSetMode()`, `PetscViewerCreate()`, `PetscViewerSetType()`, `PetscViewerBinaryOpen()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5GetBaseDimension2"))
"""
function PetscViewerHDF5GetBaseDimension2(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerHDF5GetBaseDimension2(petsclib::$UnionPetscLib, viewer::PetscViewer )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscViewerHDF5GetBaseDimension2, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{PetscBool}),
               viewer, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	PetscViewerHDF5SetSPOutput(petsclib::PetscLibType,viewer::PetscViewer, flg::PetscBool) 
Data is written to disk in single precision even if PETSc is
compiled with double precision `PetscReal`.

Logically Collective

Input Parameters:
- `viewer` - the PetscViewer; if it is a `PETSCVIEWERHDF5` then this command is ignored
- `flg`    - if `PETSC_TRUE` the data will be written to disk with single precision

Options Database Key:
- `-viewer_hdf5_sp_output` - turns on (true) or off (false) output in single precision

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerFileSetMode()`, `PetscViewerCreate()`, `PetscViewerSetType()`, `PetscViewerBinaryOpen()`,
`PetscReal`, `PetscViewerHDF5GetSPOutput()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5SetSPOutput"))
"""
function PetscViewerHDF5SetSPOutput(petsclib::PetscLibType, viewer::PetscViewer, flg::PetscBool) end

@for_petsc function PetscViewerHDF5SetSPOutput(petsclib::$UnionPetscLib, viewer::PetscViewer, flg::PetscBool )

    @chk ccall(
               (:PetscViewerHDF5SetSPOutput, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscBool),
               viewer, flg,
              )


	return nothing
end 

"""
	flg::PetscBool = PetscViewerHDF5GetSPOutput(petsclib::PetscLibType,viewer::PetscViewer) 
Data is written to disk in single precision even if PETSc is
compiled with double precision `PetscReal`.

Logically Collective

Input Parameter:
- `viewer` - the PetscViewer, must be of type `PETSCVIEWERHDF5`

Output Parameter:
- `flg` - if `PETSC_TRUE` the data will be written to disk with single precision

Level: intermediate

-seealso: [](sec_viewers), `PetscViewerFileSetMode()`, `PetscViewerCreate()`, `PetscViewerSetType()`, `PetscViewerBinaryOpen()`,
`PetscReal`, `PetscViewerHDF5SetSPOutput()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5GetSPOutput"))
"""
function PetscViewerHDF5GetSPOutput(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerHDF5GetSPOutput(petsclib::$UnionPetscLib, viewer::PetscViewer )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscViewerHDF5GetSPOutput, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{PetscBool}),
               viewer, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	PetscViewerHDF5SetCollective(petsclib::PetscLibType,viewer::PetscViewer, flg::PetscBool) 
Use collective MPI

Logically Collective; flg must contain common value

Input Parameters:
- `viewer` - the `PetscViewer`; if it is not `PETSCVIEWERHDF5` then this command is ignored
- `flg`    - `PETSC_TRUE` for collective mode; `PETSC_FALSE` for independent mode (default)

Options Database Key:
- `-viewer_hdf5_collective` - turns on (true) or off (false) collective transfers

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerHDF5GetCollective()`, `PetscViewerCreate()`, `PetscViewerSetType()`, `PetscViewerHDF5Open()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5SetCollective"))
"""
function PetscViewerHDF5SetCollective(petsclib::PetscLibType, viewer::PetscViewer, flg::PetscBool) end

@for_petsc function PetscViewerHDF5SetCollective(petsclib::$UnionPetscLib, viewer::PetscViewer, flg::PetscBool )

    @chk ccall(
               (:PetscViewerHDF5SetCollective, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscBool),
               viewer, flg,
              )


	return nothing
end 

"""
	flg::PetscBool = PetscViewerHDF5GetCollective(petsclib::PetscLibType,viewer::PetscViewer) 
Return flag whether collective MPI

Not Collective

Input Parameter:
- `viewer` - the `PETSCVIEWERHDF5` `PetscViewer`

Output Parameter:
- `flg` - the flag

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerHDF5SetCollective()`, `PetscViewerCreate()`, `PetscViewerSetType()`, `PetscViewerHDF5Open()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5GetCollective"))
"""
function PetscViewerHDF5GetCollective(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerHDF5GetCollective(petsclib::$UnionPetscLib, viewer::PetscViewer )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscViewerHDF5GetCollective, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{PetscBool}),
               viewer, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	PetscViewerHDF5SetDefaultTimestepping(petsclib::PetscLibType,viewer::PetscViewer, flg::PetscBool) 
Set the flag for default timestepping

Logically Collective

Input Parameters:
- `viewer` - the `PetscViewer`; if it is not `PETSCVIEWERHDF5` then this command is ignored
- `flg`    - if `PETSC_TRUE` we will assume that timestepping is on

Options Database Key:
- `-viewer_hdf5_default_timestepping` - turns on (true) or off (false) default timestepping

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerHDF5GetDefaultTimestepping()`, `PetscViewerHDF5PushTimestepping()`, `PetscViewerHDF5GetTimestep()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5SetDefaultTimestepping"))
"""
function PetscViewerHDF5SetDefaultTimestepping(petsclib::PetscLibType, viewer::PetscViewer, flg::PetscBool) end

@for_petsc function PetscViewerHDF5SetDefaultTimestepping(petsclib::$UnionPetscLib, viewer::PetscViewer, flg::PetscBool )

    @chk ccall(
               (:PetscViewerHDF5SetDefaultTimestepping, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscBool),
               viewer, flg,
              )


	return nothing
end 

"""
	flg::PetscBool = PetscViewerHDF5GetDefaultTimestepping(petsclib::PetscLibType,viewer::PetscViewer) 
Get the flag for default timestepping

Not Collective

Input Parameter:
- `viewer` - the `PetscViewer` of type `PETSCVIEWERHDF5`

Output Parameter:
- `flg` - if `PETSC_TRUE` we will assume that timestepping is on

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerHDF5SetDefaultTimestepping()`, `PetscViewerHDF5PushTimestepping()`, `PetscViewerHDF5GetTimestep()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5GetDefaultTimestepping"))
"""
function PetscViewerHDF5GetDefaultTimestepping(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerHDF5GetDefaultTimestepping(petsclib::$UnionPetscLib, viewer::PetscViewer )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscViewerHDF5GetDefaultTimestepping, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{PetscBool}),
               viewer, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	PetscViewerHDF5SetCompress(petsclib::PetscLibType,viewer::PetscViewer, flg::PetscBool) 
Set the flag for compression

Logically Collective

Input Parameters:
- `viewer` - the `PetscViewer`; if it is not `PETSCVIEWERHDF5` then this command is ignored
- `flg`    - if `PETSC_TRUE` we will turn on compression

Options Database Key:
- `-viewer_hdf5_compress` - turns on (true) or off (false) compression

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerHDF5GetCompress()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5SetCompress"))
"""
function PetscViewerHDF5SetCompress(petsclib::PetscLibType, viewer::PetscViewer, flg::PetscBool) end

@for_petsc function PetscViewerHDF5SetCompress(petsclib::$UnionPetscLib, viewer::PetscViewer, flg::PetscBool )

    @chk ccall(
               (:PetscViewerHDF5SetCompress, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscBool),
               viewer, flg,
              )


	return nothing
end 

"""
	flg::PetscBool = PetscViewerHDF5GetCompress(petsclib::PetscLibType,viewer::PetscViewer) 
Get the flag for compression

Not Collective

Input Parameter:
- `viewer` - the `PetscViewer` of type `PETSCVIEWERHDF5`

Output Parameter:
- `flg` - if `PETSC_TRUE` we will turn on compression

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerHDF5SetCompress()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5GetCompress"))
"""
function PetscViewerHDF5GetCompress(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerHDF5GetCompress(petsclib::$UnionPetscLib, viewer::PetscViewer )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscViewerHDF5GetCompress, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{PetscBool}),
               viewer, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	PetscViewerHDF5Open(petsclib::PetscLibType,comm::MPI_Comm, name::String, type::PetscFileMode, hdf5v::PetscViewer) 
Opens a file for HDF5 input/output as a `PETSCVIEWERHDF5` `PetscViewer`

Collective

Input Parameters:
- `comm` - MPI communicator
- `name` - name of file
- `type` - type of file

Output Parameter:
- `hdf5v` - `PetscViewer` for HDF5 input/output to use with the specified file

Options Database Keys:
- `-viewer_hdf5_base_dimension2` - turns on (true) or off (false) using a dimension of 2 in the HDF5 file even if the bs/dof of the vector is 1
- `-viewer_hdf5_sp_output`       - forces (if true) the viewer to write data in single precision independent on the precision of PetscReal

Level: beginner

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerASCIIOpen()`, `PetscViewerPushFormat()`, `PetscViewerDestroy()`, `PetscViewerHDF5SetBaseDimension2()`,
`PetscViewerHDF5SetSPOutput()`, `PetscViewerHDF5GetBaseDimension2()`, `VecView()`, `MatView()`, `VecLoad()`,
`MatLoad()`, `PetscFileMode`, `PetscViewer`, `PetscViewerSetType()`, `PetscViewerFileSetMode()`, `PetscViewerFileSetName()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5Open"))
"""
function PetscViewerHDF5Open(petsclib::PetscLibType, comm::MPI_Comm, name::String, type::PetscFileMode, hdf5v::PetscViewer) end

@for_petsc function PetscViewerHDF5Open(petsclib::$UnionPetscLib, comm::MPI_Comm, name::String, type::PetscFileMode, hdf5v::PetscViewer )

    @chk ccall(
               (:PetscViewerHDF5Open, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, PetscFileMode, Ptr{PetscViewer}),
               comm, name, type, hdf5v,
              )


	return nothing
end 

"""
	PetscViewerHDF5GetFileId(petsclib::PetscLibType,viewer::PetscViewer, file_id::hid_t) 
Retrieve the file id, this file ID then can be used in direct HDF5 calls

Not Collective

Input Parameter:
- `viewer` - the `PetscViewer`

Output Parameter:
- `file_id` - The file id

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerHDF5Open()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5GetFileId"))
"""
function PetscViewerHDF5GetFileId(petsclib::PetscLibType, viewer::PetscViewer, file_id::hid_t) end

@for_petsc function PetscViewerHDF5GetFileId(petsclib::$UnionPetscLib, viewer::PetscViewer, file_id::hid_t )

    @chk ccall(
               (:PetscViewerHDF5GetFileId, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{hid_t}),
               viewer, file_id,
              )


	return nothing
end 

"""
	PetscViewerHDF5PushGroup(petsclib::PetscLibType,viewer::PetscViewer, name::String) 
Set the current HDF5 group for output

Not Collective

Input Parameters:
- `viewer` - the `PetscViewer` of type `PETSCVIEWERHDF5`
- `name`   - The group name

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerHDF5Open()`, `PetscViewerHDF5PopGroup()`, `PetscViewerHDF5GetGroup()`, `PetscViewerHDF5OpenGroup()`, `PetscViewerHDF5WriteGroup()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5PushGroup"))
"""
function PetscViewerHDF5PushGroup(petsclib::PetscLibType, viewer::PetscViewer, name::String) end

@for_petsc function PetscViewerHDF5PushGroup(petsclib::$UnionPetscLib, viewer::PetscViewer, name::String )

    @chk ccall(
               (:PetscViewerHDF5PushGroup, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cchar}),
               viewer, name,
              )


	return nothing
end 

"""
	PetscViewerHDF5PopGroup(petsclib::PetscLibType,viewer::PetscViewer) 
Return the current HDF5 group for output to the previous value

Not Collective

Input Parameter:
- `viewer` - the `PetscViewer` of type `PETSCVIEWERHDF5`

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerHDF5Open()`, `PetscViewerHDF5PushGroup()`, `PetscViewerHDF5GetGroup()`, `PetscViewerHDF5OpenGroup()`, `PetscViewerHDF5WriteGroup()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5PopGroup"))
"""
function PetscViewerHDF5PopGroup(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerHDF5PopGroup(petsclib::$UnionPetscLib, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerHDF5PopGroup, $petsc_library),
               PetscErrorCode,
               (PetscViewer,),
               viewer,
              )


	return nothing
end 

"""
	PetscViewerHDF5OpenGroup(petsclib::PetscLibType,viewer::PetscViewer, path::String, fileId::hid_t, groupId::hid_t) 
Open the HDF5 group with the name (full path) returned by `PetscViewerHDF5GetGroup()`,
and return this group's ID and file ID.
If `PetscViewerHDF5GetGroup()` yields NULL, then group ID is file ID.

Not Collective

Input Parameters:
- `viewer` - the `PetscViewer` of type `PETSCVIEWERHDF5`
- `path`   - (Optional) The path relative to the pushed group

Output Parameters:
- `fileId`  - The HDF5 file ID
- `groupId` - The HDF5 group ID

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerHDF5Open()`, `PetscViewerHDF5PushGroup()`, `PetscViewerHDF5PopGroup()`, `PetscViewerHDF5GetGroup()`, `PetscViewerHDF5WriteGroup()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5OpenGroup"))
"""
function PetscViewerHDF5OpenGroup(petsclib::PetscLibType, viewer::PetscViewer, path::String, fileId::hid_t, groupId::hid_t) end

@for_petsc function PetscViewerHDF5OpenGroup(petsclib::$UnionPetscLib, viewer::PetscViewer, path::String, fileId::hid_t, groupId::hid_t )

    @chk ccall(
               (:PetscViewerHDF5OpenGroup, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cchar}, Ptr{hid_t}, Ptr{hid_t}),
               viewer, path, fileId, groupId,
              )


	return nothing
end 

"""
	PetscViewerHDF5WriteGroup(petsclib::PetscLibType,viewer::PetscViewer, path::String) 
Ensure the HDF5 group exists in the HDF5 file

Not Collective

Input Parameters:
- `viewer` - the `PetscViewer` of type `PETSCVIEWERHDF5`
- `path`   - (Optional) The path relative to the pushed group

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerHDF5Open()`, `PetscViewerHDF5PushGroup()`, `PetscViewerHDF5PopGroup()`, `PetscViewerHDF5GetGroup()`, `PetscViewerHDF5OpenGroup()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5WriteGroup"))
"""
function PetscViewerHDF5WriteGroup(petsclib::PetscLibType, viewer::PetscViewer, path::String) end

@for_petsc function PetscViewerHDF5WriteGroup(petsclib::$UnionPetscLib, viewer::PetscViewer, path::String )

    @chk ccall(
               (:PetscViewerHDF5WriteGroup, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cchar}),
               viewer, path,
              )


	return nothing
end 

"""
	PetscViewerHDF5PushTimestepping(petsclib::PetscLibType,viewer::PetscViewer) 
Activate timestepping mode for subsequent HDF5 reading and writing.

Not Collective

Input Parameter:
- `viewer` - the `PetscViewer` of type `PETSCVIEWERHDF5`

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerHDF5Open()`, `PetscViewerHDF5PopTimestepping()`, `PetscViewerHDF5IsTimestepping()`, `PetscViewerHDF5SetTimestep()`, `PetscViewerHDF5IncrementTimestep()`, `PetscViewerHDF5GetTimestep()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5PushTimestepping"))
"""
function PetscViewerHDF5PushTimestepping(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerHDF5PushTimestepping(petsclib::$UnionPetscLib, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerHDF5PushTimestepping, $petsc_library),
               PetscErrorCode,
               (PetscViewer,),
               viewer,
              )


	return nothing
end 

"""
	PetscViewerHDF5PopTimestepping(petsclib::PetscLibType,viewer::PetscViewer) 
Deactivate timestepping mode for subsequent HDF5 reading and writing.

Not Collective

Input Parameter:
- `viewer` - the `PetscViewer` of type `PETSCVIEWERHDF5`

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerHDF5Open()`, `PetscViewerHDF5PushTimestepping()`, `PetscViewerHDF5IsTimestepping()`, `PetscViewerHDF5SetTimestep()`, `PetscViewerHDF5IncrementTimestep()`, `PetscViewerHDF5GetTimestep()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5PopTimestepping"))
"""
function PetscViewerHDF5PopTimestepping(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerHDF5PopTimestepping(petsclib::$UnionPetscLib, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerHDF5PopTimestepping, $petsc_library),
               PetscErrorCode,
               (PetscViewer,),
               viewer,
              )


	return nothing
end 

"""
	flg::PetscBool = PetscViewerHDF5IsTimestepping(petsclib::PetscLibType,viewer::PetscViewer) 
Ask the viewer whether it is in timestepping mode currently.

Not Collective

Input Parameter:
- `viewer` - the `PetscViewer` of type `PETSCVIEWERHDF5`

Output Parameter:
- `flg` - is timestepping active?

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerHDF5Open()`, `PetscViewerHDF5PushTimestepping()`, `PetscViewerHDF5PopTimestepping()`, `PetscViewerHDF5SetTimestep()`, `PetscViewerHDF5IncrementTimestep()`, `PetscViewerHDF5GetTimestep()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5IsTimestepping"))
"""
function PetscViewerHDF5IsTimestepping(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerHDF5IsTimestepping(petsclib::$UnionPetscLib, viewer::PetscViewer )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscViewerHDF5IsTimestepping, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{PetscBool}),
               viewer, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	PetscViewerHDF5IncrementTimestep(petsclib::PetscLibType,viewer::PetscViewer) 
Increments current timestep for the HDF5 output. Fields are stacked in time.

Not Collective

Input Parameter:
- `viewer` - the `PetscViewer` of type `PETSCVIEWERHDF5`

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerHDF5Open()`, `PetscViewerHDF5PushTimestepping()`, `PetscViewerHDF5SetTimestep()`, `PetscViewerHDF5GetTimestep()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5IncrementTimestep"))
"""
function PetscViewerHDF5IncrementTimestep(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerHDF5IncrementTimestep(petsclib::$UnionPetscLib, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerHDF5IncrementTimestep, $petsc_library),
               PetscErrorCode,
               (PetscViewer,),
               viewer,
              )


	return nothing
end 

"""
	PetscViewerHDF5SetTimestep(petsclib::PetscLibType,viewer::PetscViewer, timestep::PetscInt) 
Set the current timestep for the HDF5 output. Fields are stacked in time.

Logically Collective

Input Parameters:
- `viewer`   - the `PetscViewer` of type `PETSCVIEWERHDF5`
- `timestep` - The timestep

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerHDF5Open()`, `PetscViewerHDF5PushTimestepping()`, `PetscViewerHDF5IncrementTimestep()`, `PetscViewerHDF5GetTimestep()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5SetTimestep"))
"""
function PetscViewerHDF5SetTimestep(petsclib::PetscLibType, viewer::PetscViewer, timestep::PetscInt) end

@for_petsc function PetscViewerHDF5SetTimestep(petsclib::$UnionPetscLib, viewer::PetscViewer, timestep::$PetscInt )

    @chk ccall(
               (:PetscViewerHDF5SetTimestep, $petsc_library),
               PetscErrorCode,
               (PetscViewer, $PetscInt),
               viewer, timestep,
              )


	return nothing
end 

"""
	timestep::PetscInt = PetscViewerHDF5GetTimestep(petsclib::PetscLibType,viewer::PetscViewer) 
Get the current timestep for the HDF5 output. Fields are stacked in time.

Not Collective

Input Parameter:
- `viewer` - the `PetscViewer` of type `PETSCVIEWERHDF5`

Output Parameter:
- `timestep` - The timestep

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerHDF5Open()`, `PetscViewerHDF5PushTimestepping()`, `PetscViewerHDF5IncrementTimestep()`, `PetscViewerHDF5SetTimestep()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5GetTimestep"))
"""
function PetscViewerHDF5GetTimestep(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerHDF5GetTimestep(petsclib::$UnionPetscLib, viewer::PetscViewer )
	timestep_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscViewerHDF5GetTimestep, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{$PetscInt}),
               viewer, timestep_,
              )

	timestep = timestep_[]

	return timestep
end 

"""
	PetscViewerHDF5WriteAttribute(petsclib::PetscLibType,viewer::PetscViewer, parent::String, name::String, datatype::PetscDataType, value::Cvoid) 
Write an attribute

Collective

Input Parameters:
- `viewer`   - The `PETSCVIEWERHDF5` viewer
- `parent`   - The parent dataset/group name
- `name`     - The attribute name
- `datatype` - The attribute type
- `value`    - The attribute value

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerHDF5Open()`, `PetscViewerHDF5WriteObjectAttribute()`, `PetscViewerHDF5ReadAttribute()`, `PetscViewerHDF5HasAttribute()`,
`PetscViewerHDF5PushGroup()`, `PetscViewerHDF5PopGroup()`, `PetscViewerHDF5GetGroup()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5WriteAttribute"))
"""
function PetscViewerHDF5WriteAttribute(petsclib::PetscLibType, viewer::PetscViewer, parent::String, name::String, datatype::PetscDataType, value::Cvoid) end

@for_petsc function PetscViewerHDF5WriteAttribute(petsclib::$UnionPetscLib, viewer::PetscViewer, parent::String, name::String, datatype::PetscDataType, value::Cvoid )

    @chk ccall(
               (:PetscViewerHDF5WriteAttribute, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cchar}, Ptr{Cchar}, PetscDataType, Ptr{Cvoid}),
               viewer, parent, name, datatype, value,
              )


	return nothing
end 

"""
	PetscViewerHDF5WriteObjectAttribute(petsclib::PetscLibType,viewer::PetscViewer, obj::PetscObject, name::String, datatype::PetscDataType, value::Cvoid) 
Write an attribute to the dataset matching the given `PetscObject` by name

Collective

Input Parameters:
- `viewer`   - The `PETSCVIEWERHDF5` viewer
- `obj`      - The object whose name is used to lookup the parent dataset, relative to the current group.
- `name`     - The attribute name
- `datatype` - The attribute type
- `value`    - The attribute value

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerHDF5Open()`, `PetscViewerHDF5WriteAttribute()`, `PetscViewerHDF5ReadObjectAttribute()`, `PetscViewerHDF5HasObjectAttribute()`,
`PetscViewerHDF5HasObject()`, `PetscViewerHDF5PushGroup()`, `PetscViewerHDF5PopGroup()`, `PetscViewerHDF5GetGroup()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5WriteObjectAttribute"))
"""
function PetscViewerHDF5WriteObjectAttribute(petsclib::PetscLibType, viewer::PetscViewer, obj::PetscObject, name::String, datatype::PetscDataType, value::Cvoid) end

@for_petsc function PetscViewerHDF5WriteObjectAttribute(petsclib::$UnionPetscLib, viewer::PetscViewer, obj::PetscObject, name::String, datatype::PetscDataType, value::Cvoid )

    @chk ccall(
               (:PetscViewerHDF5WriteObjectAttribute, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscObject, Ptr{Cchar}, PetscDataType, Ptr{Cvoid}),
               viewer, obj, name, datatype, value,
              )


	return nothing
end 

"""
	PetscViewerHDF5ReadAttribute(petsclib::PetscLibType,viewer::PetscViewer, parent::String, name::String, datatype::PetscDataType, defaultValue::Cvoid, value::Cvoid) 
Read an attribute

Collective

Input Parameters:
- `viewer`       - The `PETSCVIEWERHDF5` viewer
- `parent`       - The parent dataset/group name
- `name`         - The attribute name
- `datatype`     - The attribute type
- `defaultValue` - The pointer to the default value

Output Parameter:
- `value` - The pointer to the read HDF5 attribute value

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerHDF5Open()`, `PetscViewerHDF5ReadObjectAttribute()`, `PetscViewerHDF5WriteAttribute()`, `PetscViewerHDF5HasAttribute()`, `PetscViewerHDF5HasObject()`, `PetscViewerHDF5PushGroup()`, `PetscViewerHDF5PopGroup()`, `PetscViewerHDF5GetGroup()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5ReadAttribute"))
"""
function PetscViewerHDF5ReadAttribute(petsclib::PetscLibType, viewer::PetscViewer, parent::String, name::String, datatype::PetscDataType, defaultValue::Cvoid, value::Cvoid) end

@for_petsc function PetscViewerHDF5ReadAttribute(petsclib::$UnionPetscLib, viewer::PetscViewer, parent::String, name::String, datatype::PetscDataType, defaultValue::Cvoid, value::Cvoid )

    @chk ccall(
               (:PetscViewerHDF5ReadAttribute, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cchar}, Ptr{Cchar}, PetscDataType, Ptr{Cvoid}, Ptr{Cvoid}),
               viewer, parent, name, datatype, defaultValue, value,
              )


	return nothing
end 

"""
	PetscViewerHDF5ReadObjectAttribute(petsclib::PetscLibType,viewer::PetscViewer, obj::PetscObject, name::String, datatype::PetscDataType, defaultValue::Cvoid, value::Cvoid) 
Read an attribute from the dataset matching the given `PetscObject` by name

Collective

Input Parameters:
- `viewer`       - The `PETSCVIEWERHDF5` viewer
- `obj`          - The object whose name is used to lookup the parent dataset, relative to the current group.
- `name`         - The attribute name
- `datatype`     - The attribute type
- `defaultValue` - The default attribute value

Output Parameter:
- `value` - The attribute value

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerHDF5Open()`, `PetscViewerHDF5ReadAttribute()` `PetscViewerHDF5WriteObjectAttribute()`, `PetscViewerHDF5HasObjectAttribute()`, `PetscViewerHDF5PushGroup()`, `PetscViewerHDF5PopGroup()`, `PetscViewerHDF5GetGroup()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5ReadObjectAttribute"))
"""
function PetscViewerHDF5ReadObjectAttribute(petsclib::PetscLibType, viewer::PetscViewer, obj::PetscObject, name::String, datatype::PetscDataType, defaultValue::Cvoid, value::Cvoid) end

@for_petsc function PetscViewerHDF5ReadObjectAttribute(petsclib::$UnionPetscLib, viewer::PetscViewer, obj::PetscObject, name::String, datatype::PetscDataType, defaultValue::Cvoid, value::Cvoid )

    @chk ccall(
               (:PetscViewerHDF5ReadObjectAttribute, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscObject, Ptr{Cchar}, PetscDataType, Ptr{Cvoid}, Ptr{Cvoid}),
               viewer, obj, name, datatype, defaultValue, value,
              )


	return nothing
end 

"""
	has::PetscBool = PetscViewerHDF5HasGroup(petsclib::PetscLibType,viewer::PetscViewer, path::String) 
Check whether the current (pushed) group exists in the HDF5 file

Collective

Input Parameters:
- `viewer` - The `PETSCVIEWERHDF5` viewer
- `path`   - (Optional) The path relative to the pushed group

Output Parameter:
- `has` - Flag for group existence

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerHDF5HasAttribute()`, `PetscViewerHDF5HasDataset()`, `PetscViewerHDF5PushGroup()`, `PetscViewerHDF5PopGroup()`, `PetscViewerHDF5GetGroup()`, `PetscViewerHDF5OpenGroup()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5HasGroup"))
"""
function PetscViewerHDF5HasGroup(petsclib::PetscLibType, viewer::PetscViewer, path::String) end

@for_petsc function PetscViewerHDF5HasGroup(petsclib::$UnionPetscLib, viewer::PetscViewer, path::String )
	has_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscViewerHDF5HasGroup, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cchar}, Ptr{PetscBool}),
               viewer, path, has_,
              )

	has = has_[]

	return has
end 

"""
	has::PetscBool = PetscViewerHDF5HasDataset(petsclib::PetscLibType,viewer::PetscViewer, path::String) 
Check whether a given dataset exists in the HDF5 file

Collective

Input Parameters:
- `viewer` - The `PETSCVIEWERHDF5` viewer
- `path`   - The dataset path

Output Parameter:
- `has` - Flag whether dataset exists

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerHDF5HasObject()`, `PetscViewerHDF5HasAttribute()`, `PetscViewerHDF5HasGroup()`, `PetscViewerHDF5PushGroup()`, `PetscViewerHDF5PopGroup()`, `PetscViewerHDF5GetGroup()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5HasDataset"))
"""
function PetscViewerHDF5HasDataset(petsclib::PetscLibType, viewer::PetscViewer, path::String) end

@for_petsc function PetscViewerHDF5HasDataset(petsclib::$UnionPetscLib, viewer::PetscViewer, path::String )
	has_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscViewerHDF5HasDataset, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cchar}, Ptr{PetscBool}),
               viewer, path, has_,
              )

	has = has_[]

	return has
end 

"""
	has::PetscBool = PetscViewerHDF5HasObject(petsclib::PetscLibType,viewer::PetscViewer, obj::PetscObject) 
Check whether a dataset with the same name as given object exists in the HDF5 file under current group

Collective

Input Parameters:
- `viewer` - The `PETSCVIEWERHDF5` viewer
- `obj`    - The named object

Output Parameter:
- `has` - Flag for dataset existence

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerHDF5Open()`, `PetscViewerHDF5HasDataset()`, `PetscViewerHDF5HasAttribute()`, `PetscViewerHDF5PushGroup()`, `PetscViewerHDF5PopGroup()`, `PetscViewerHDF5GetGroup()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5HasObject"))
"""
function PetscViewerHDF5HasObject(petsclib::PetscLibType, viewer::PetscViewer, obj::PetscObject) end

@for_petsc function PetscViewerHDF5HasObject(petsclib::$UnionPetscLib, viewer::PetscViewer, obj::PetscObject )
	has_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscViewerHDF5HasObject, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscObject, Ptr{PetscBool}),
               viewer, obj, has_,
              )

	has = has_[]

	return has
end 

"""
	has::PetscBool = PetscViewerHDF5HasAttribute(petsclib::PetscLibType,viewer::PetscViewer, parent::String, name::String) 
Check whether an attribute exists

Collective

Input Parameters:
- `viewer` - The `PETSCVIEWERHDF5` viewer
- `parent` - The parent dataset/group name
- `name`   - The attribute name

Output Parameter:
- `has` - Flag for attribute existence

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerHDF5Open()`, `PetscViewerHDF5HasObjectAttribute()`, `PetscViewerHDF5WriteAttribute()`, `PetscViewerHDF5ReadAttribute()`, `PetscViewerHDF5PushGroup()`, `PetscViewerHDF5PopGroup()`, `PetscViewerHDF5GetGroup()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5HasAttribute"))
"""
function PetscViewerHDF5HasAttribute(petsclib::PetscLibType, viewer::PetscViewer, parent::String, name::String) end

@for_petsc function PetscViewerHDF5HasAttribute(petsclib::$UnionPetscLib, viewer::PetscViewer, parent::String, name::String )
	has_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscViewerHDF5HasAttribute, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscBool}),
               viewer, parent, name, has_,
              )

	has = has_[]

	return has
end 

"""
	has::PetscBool = PetscViewerHDF5HasObjectAttribute(petsclib::PetscLibType,viewer::PetscViewer, obj::PetscObject, name::String) 
Check whether an attribute is attached to the dataset matching the given `PetscObject` by name

Collective

Input Parameters:
- `viewer` - The `PETSCVIEWERHDF5` viewer
- `obj`    - The object whose name is used to lookup the parent dataset, relative to the current group.
- `name`   - The attribute name

Output Parameter:
- `has` - Flag for attribute existence

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERHDF5`, `PetscViewerHDF5Open()`, `PetscViewerHDF5HasAttribute()`, `PetscViewerHDF5WriteObjectAttribute()`, `PetscViewerHDF5ReadObjectAttribute()`, `PetscViewerHDF5HasObject()`, `PetscViewerHDF5PushGroup()`, `PetscViewerHDF5PopGroup()`, `PetscViewerHDF5GetGroup()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5HasObjectAttribute"))
"""
function PetscViewerHDF5HasObjectAttribute(petsclib::PetscLibType, viewer::PetscViewer, obj::PetscObject, name::String) end

@for_petsc function PetscViewerHDF5HasObjectAttribute(petsclib::$UnionPetscLib, viewer::PetscViewer, obj::PetscObject, name::String )
	has_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscViewerHDF5HasObjectAttribute, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscObject, Ptr{Cchar}, Ptr{PetscBool}),
               viewer, obj, name, has_,
              )

	has = has_[]

	return has
end 

"""
	PetscViewerCGNSOpen(petsclib::PetscLibType,comm::MPI_Comm, name::String, type::PetscFileMode, viewer::PetscViewer) 
Opens a file for CGNS input/output.

Collective

Input Parameters:
- `comm` - MPI communicator
- `name` - name of file
- `type` - type of file
-seealso: `PETSCVIEWERCGNS`, `PetscViewer`, `PetscViewerPushFormat()`, `PetscViewerDestroy()`,
`DMLoad()`, `PetscFileMode`, `PetscViewerSetType()`, `PetscViewerFileSetMode()`, `PetscViewerFileSetName()`

# External Links
$(_doc_external("Sys/PetscViewerCGNSOpen"))
"""
function PetscViewerCGNSOpen(petsclib::PetscLibType, comm::MPI_Comm, name::String, type::PetscFileMode, viewer::PetscViewer) end

@for_petsc function PetscViewerCGNSOpen(petsclib::$UnionPetscLib, comm::MPI_Comm, name::String, type::PetscFileMode, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerCGNSOpen, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, PetscFileMode, Ptr{PetscViewer}),
               comm, name, type, viewer,
              )


	return nothing
end 

"""
	PetscViewerCGNSSetSolutionIndex(petsclib::PetscLibType,viewer::PetscViewer, solution_id::PetscInt) 
Set index of solution

Not Collective

Input Parameters:
- `viewer`      - `PETSCVIEWERCGNS` `PetscViewer` for CGNS input/output to use with the specified file
- `solution_id` - Index of the solution id, or `-1` for the last solution on the file

Level: intermediate

-seealso: `PETSCVIEWERCGNS`, `PetscViewerCGNSGetSolutionIndex()`, `PetscViewerCGNSGetSolutionInfo()`


# External Links
$(_doc_external("Sys/PetscViewerCGNSSetSolutionIndex"))
"""
function PetscViewerCGNSSetSolutionIndex(petsclib::PetscLibType, viewer::PetscViewer, solution_id::PetscInt) end

@for_petsc function PetscViewerCGNSSetSolutionIndex(petsclib::$UnionPetscLib, viewer::PetscViewer, solution_id::$PetscInt )

    @chk ccall(
               (:PetscViewerCGNSSetSolutionIndex, $petsc_library),
               PetscErrorCode,
               (PetscViewer, $PetscInt),
               viewer, solution_id,
              )


	return nothing
end 

"""
	solution_id::PetscInt = PetscViewerCGNSGetSolutionIndex(petsclib::PetscLibType,viewer::PetscViewer) 
Get index of solution

Not Collective

Input Parameter:
- `viewer` - `PETSCVIEWERCGNS` `PetscViewer` for CGNS input/output to use with the specified file

Output Parameter:
- `solution_id` - Index of the solution id

Level: intermediate

-seealso: `PETSCVIEWERCGNS`, `PetscViewerCGNSSetSolutionIndex()`, `PetscViewerCGNSGetSolutionInfo()`


# External Links
$(_doc_external("Sys/PetscViewerCGNSGetSolutionIndex"))
"""
function PetscViewerCGNSGetSolutionIndex(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerCGNSGetSolutionIndex(petsclib::$UnionPetscLib, viewer::PetscViewer )
	solution_id_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscViewerCGNSGetSolutionIndex, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{$PetscInt}),
               viewer, solution_id_,
              )

	solution_id = solution_id_[]

	return solution_id
end 

"""
	time::PetscReal,set::PetscBool = PetscViewerCGNSGetSolutionTime(petsclib::PetscLibType,viewer::PetscViewer) 
Gets the solution time for the FlowSolution of the viewer

Collective

Input Parameter:
- `viewer` - `PETSCVIEWERCGNS` `PetscViewer` for CGNS input/output to use with the specified file

Output Parameters:
- `time` - Solution time of the FlowSolution_t node
- `set`  - Whether the time data is in the file

Level: intermediate

-seealso: `PETSCVIEWERCGNS`, `PetscViewer`, `PetscViewerCGNSGetSolutionIteration()`, `PetscViewerCGNSSetSolutionIndex()`, `PetscViewerCGNSGetSolutionIndex()`, `PetscViewerCGNSGetSolutionName()`

# External Links
$(_doc_external("Sys/PetscViewerCGNSGetSolutionTime"))
"""
function PetscViewerCGNSGetSolutionTime(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerCGNSGetSolutionTime(petsclib::$UnionPetscLib, viewer::PetscViewer )
	time_ = Ref{$PetscReal}()
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscViewerCGNSGetSolutionTime, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{$PetscReal}, Ptr{PetscBool}),
               viewer, time_, set_,
              )

	time = time_[]
	set = set_[]

	return time,set
end 

"""
	iteration::PetscInt,set::PetscBool = PetscViewerCGNSGetSolutionIteration(petsclib::PetscLibType,viewer::PetscViewer) 
Gets the solution iteration for the FlowSolution of the viewer

Collective

Input Parameter:
- `viewer` - `PETSCVIEWERCGNS` `PetscViewer` for CGNS input/output to use with the specified file

Output Parameters:
- `iteration` - Solution iteration of the FlowSolution_t node
- `set`       - Whether the time data is in the file

Level: intermediate

-seealso: `PETSCVIEWERCGNS`, `PetscViewer`, `PetscViewerCGNSGetSolutionTime()`, `PetscViewerCGNSSetSolutionIndex()`, `PetscViewerCGNSGetSolutionIndex()`, `PetscViewerCGNSGetSolutionName()`

# External Links
$(_doc_external("Sys/PetscViewerCGNSGetSolutionIteration"))
"""
function PetscViewerCGNSGetSolutionIteration(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerCGNSGetSolutionIteration(petsclib::$UnionPetscLib, viewer::PetscViewer )
	iteration_ = Ref{$PetscInt}()
	set_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscViewerCGNSGetSolutionIteration, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{$PetscInt}, Ptr{PetscBool}),
               viewer, iteration_, set_,
              )

	iteration = iteration_[]
	set = set_[]

	return iteration,set
end 

"""
	PetscViewerCGNSGetSolutionName(petsclib::PetscLibType,viewer::PetscViewer, name::String) 
Gets name of FlowSolution of the viewer

Collective

Input Parameter:
- `viewer` - `PETSCVIEWERCGNS` `PetscViewer` for CGNS input/output to use with the specified file

Output Parameter:
- `name` - Name of the FlowSolution_t node corresponding to the solution index

Level: intermediate

-seealso: `PETSCVIEWERCGNS`, `PetscViewer`, `PetscViewerCGNSSetSolutionIndex()`, `PetscViewerCGNSGetSolutionIndex()`, `PetscViewerCGNSGetSolutionTime()`

# External Links
$(_doc_external("Sys/PetscViewerCGNSGetSolutionName"))
"""
function PetscViewerCGNSGetSolutionName(petsclib::PetscLibType, viewer::PetscViewer, name::String) end

@for_petsc function PetscViewerCGNSGetSolutionName(petsclib::$UnionPetscLib, viewer::PetscViewer, name::String )
	name_ = Ref(pointer(name))

    @chk ccall(
               (:PetscViewerCGNSGetSolutionName, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Ptr{Cchar}}),
               viewer, name_,
              )


	return nothing
end 

"""
	PetscViewerMathematicaFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the PETSc interface to Mathematica. It is
called from PetscFinalize().

Level: developer

-seealso: `PetscFinalize()`

# External Links
$(_doc_external("Sys/PetscViewerMathematicaFinalizePackage"))
"""
function PetscViewerMathematicaFinalizePackage(petsclib::PetscLibType) end

@for_petsc function PetscViewerMathematicaFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscViewerMathematicaFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscViewerMathematicaInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the PETSc interface to Mathematica. It is
called from `PetscViewerInitializePackage()`.

Level: developer

-seealso: `PetscSysInitializePackage()`, `PetscInitialize()`

# External Links
$(_doc_external("Sys/PetscViewerMathematicaInitializePackage"))
"""
function PetscViewerMathematicaInitializePackage(petsclib::PetscLibType) end

@for_petsc function PetscViewerMathematicaInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscViewerMathematicaInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscViewerMathematicaOpen(petsclib::PetscLibType,comm::MPI_Comm, port::Cint, machine::String, mode::String, v::PetscViewer) 
Communicates with Mathemtica using MathLink.

Collective

Input Parameters:
- `comm`    - The MPI communicator
- `port`    - [optional] The port to connect on, or PETSC_DECIDE
- `machine` - [optional] The machine to run Mathematica on, or NULL
- `mode`    - [optional] The connection mode, or NULL

Output Parameter:
- `v` - The Mathematica viewer

Options Database Keys:
- `-viewer_math_linkhost <machine>` - The host machine for the kernel
- `-viewer_math_linkname <name>`    - The full link name for the connection
- `-viewer_math_linkport <port>`    - The port for the connection
- `-viewer_math_mode <mode>`        - The mode, e.g. Launch, Connect
- `-viewer_math_type <type>`        - The plot type, e.g. Triangulation, Vector
- `-viewer_math_graphics <output>`  - The output type, e.g. Motif, PS, PSFile

Level: intermediate

-seealso: `PETSCVIEWERMATHEMATICA`, `MatView()`, `VecView()`

# External Links
$(_doc_external("Sys/PetscViewerMathematicaOpen"))
"""
function PetscViewerMathematicaOpen(petsclib::PetscLibType, comm::MPI_Comm, port::Cint, machine::String, mode::String, v::PetscViewer) end

@for_petsc function PetscViewerMathematicaOpen(petsclib::$UnionPetscLib, comm::MPI_Comm, port::Cint, machine::String, mode::String, v::PetscViewer )

    @chk ccall(
               (:PetscViewerMathematicaOpen, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Cint, Ptr{Cchar}, Ptr{Cchar}, Ptr{PetscViewer}),
               comm, port, machine, mode, v,
              )


	return nothing
end 

"""
	PetscViewerMathematicaSkipPackets(petsclib::PetscLibType,viewer::PetscViewer, type::Cint) 
Discard packets sent by Mathematica until a certain packet type is received

Input Parameters:
- `viewer` - The Mathematica viewer
- `type`   - The packet type to search for, e.g RETURNPKT

Level: advanced

-seealso: `PetscViewerMathematicaSetName()`, `PetscViewerMathematicaGetVector()`

# External Links
$(_doc_external("Sys/PetscViewerMathematicaSkipPackets"))
"""
function PetscViewerMathematicaSkipPackets(petsclib::PetscLibType, viewer::PetscViewer, type::Cint) end

@for_petsc function PetscViewerMathematicaSkipPackets(petsclib::$UnionPetscLib, viewer::PetscViewer, type::Cint )

    @chk ccall(
               (:PetscViewerMathematicaSkipPackets, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Cint),
               viewer, type,
              )


	return nothing
end 

"""
	PetscViewerMathematicaGetName(petsclib::PetscLibType,viewer::PetscViewer, name::Cchar) 
Retrieve the default name for objects communicated to Mathematica via `PETSCVIEWERMATHEMATICA`

Input Parameter:
- `viewer` - The Mathematica viewer

Output Parameter:
- `name` - The name for new objects created in Mathematica

Level: intermediate

-seealso: `PETSCVIEWERMATHEMATICA`, `PetscViewerMathematicaSetName()`, `PetscViewerMathematicaClearName()`

# External Links
$(_doc_external("Sys/PetscViewerMathematicaGetName"))
"""
function PetscViewerMathematicaGetName(petsclib::PetscLibType, viewer::PetscViewer, name::Cchar) end

@for_petsc function PetscViewerMathematicaGetName(petsclib::$UnionPetscLib, viewer::PetscViewer, name::Cchar )

    @chk ccall(
               (:PetscViewerMathematicaGetName, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Cchar),
               viewer, name,
              )


	return nothing
end 

"""
	PetscViewerMathematicaSetName(petsclib::PetscLibType,viewer::PetscViewer, name::String) 
Override the default name for objects communicated to Mathematica via `PETSCVIEWERMATHEMATICA`

Input Parameters:
- `viewer` - The Mathematica viewer
- `name`   - The name for new objects created in Mathematica

Level: intermediate

-seealso: `PETSCVIEWERMATHEMATICA`, `PetscViewerMathematicaClearName()`

# External Links
$(_doc_external("Sys/PetscViewerMathematicaSetName"))
"""
function PetscViewerMathematicaSetName(petsclib::PetscLibType, viewer::PetscViewer, name::String) end

@for_petsc function PetscViewerMathematicaSetName(petsclib::$UnionPetscLib, viewer::PetscViewer, name::String )

    @chk ccall(
               (:PetscViewerMathematicaSetName, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cchar}),
               viewer, name,
              )


	return nothing
end 

"""
	PetscViewerMathematicaClearName(petsclib::PetscLibType,viewer::PetscViewer) 
Use the default name for objects communicated to Mathematica

Input Parameter:
- `viewer` - The Mathematica viewer

Level: intermediate

-seealso: `PETSCVIEWERMATHEMATICA`, `PetscViewerMathematicaGetName()`, `PetscViewerMathematicaSetName()`

# External Links
$(_doc_external("Sys/PetscViewerMathematicaClearName"))
"""
function PetscViewerMathematicaClearName(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerMathematicaClearName(petsclib::$UnionPetscLib, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerMathematicaClearName, $petsc_library),
               PetscErrorCode,
               (PetscViewer,),
               viewer,
              )


	return nothing
end 

"""
	a::PetscReal = PetscViewerMathematicaPutMatrix(petsclib::PetscLibType,viewer::PetscViewer, m::Cint, n::Cint) 

# External Links
$(_doc_external("Sys/PetscViewerMathematicaPutMatrix"))
"""
function PetscViewerMathematicaPutMatrix(petsclib::PetscLibType, viewer::PetscViewer, m::Cint, n::Cint) end

@for_petsc function PetscViewerMathematicaPutMatrix(petsclib::$UnionPetscLib, viewer::PetscViewer, m::Cint, n::Cint )
	a_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscViewerMathematicaPutMatrix, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Cint, Cint, Ptr{$PetscReal}),
               viewer, m, n, a_,
              )

	a = a_[]

	return a
end 

"""
	a::PetscReal = PetscViewerMathematicaPutCSRMatrix(petsclib::PetscLibType,viewer::PetscViewer, m::Cint, n::Cint, i::Cint, j::Cint) 

# External Links
$(_doc_external("Sys/PetscViewerMathematicaPutCSRMatrix"))
"""
function PetscViewerMathematicaPutCSRMatrix(petsclib::PetscLibType, viewer::PetscViewer, m::Cint, n::Cint, i::Cint, j::Cint) end

@for_petsc function PetscViewerMathematicaPutCSRMatrix(petsclib::$UnionPetscLib, viewer::PetscViewer, m::Cint, n::Cint, i::Cint, j::Cint )
	a_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscViewerMathematicaPutCSRMatrix, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{$PetscReal}),
               viewer, m, n, i, j, a_,
              )

	a = a_[]

	return a
end 

"""
	PetscViewerVTKAddField(petsclib::PetscLibType,viewer::PetscViewer, dm::PetscObject, PetscViewerVTKWriteFunction::external, fieldnum::PetscInt, fieldtype::PetscViewerVTKFieldType, checkdm::PetscBool, vec::PetscObject) 
Add a field to the viewer

Collective

Input Parameters:
- `viewer`                      - `PETSCVIEWERVTK`
- `dm`                          - `DM` on which `Vec` lives
- `PetscViewerVTKWriteFunction` - function to write this `Vec`
- `fieldnum`                    - which field of the `DM` to write (`PETSC_DEFAULT` if the whole vector should be written)
- `fieldtype`                   - Either `PETSC_VTK_POINT_FIELD` or `PETSC_VTK_CELL_FIELD`
- `checkdm`                     - whether to check for identical dm arguments as fields are added
- `vec`                         - `Vec` from which to write

Level: developer

-seealso: [](sec_viewers), `PETSCVIEWERVTK`, `PetscViewerVTKOpen()`, `DMDAVTKWriteAll()`, `PetscViewerVTKWriteFunction`, `PetscViewerVTKGetDM()`

# External Links
$(_doc_external("Sys/PetscViewerVTKAddField"))
"""
function PetscViewerVTKAddField(petsclib::PetscLibType, viewer::PetscViewer, dm::PetscObject, PetscViewerVTKWriteFunction::external, fieldnum::PetscInt, fieldtype::PetscViewerVTKFieldType, checkdm::PetscBool, vec::PetscObject) end

@for_petsc function PetscViewerVTKAddField(petsclib::$UnionPetscLib, viewer::PetscViewer, dm::PetscObject, PetscViewerVTKWriteFunction::external, fieldnum::$PetscInt, fieldtype::PetscViewerVTKFieldType, checkdm::PetscBool, vec::PetscObject )

    @chk ccall(
               (:PetscViewerVTKAddField, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscObject, external, $PetscInt, PetscViewerVTKFieldType, PetscBool, PetscObject),
               viewer, dm, PetscViewerVTKWriteFunction, fieldnum, fieldtype, checkdm, vec,
              )


	return nothing
end 

"""
	PetscViewerVTKGetDM(petsclib::PetscLibType,viewer::PetscViewer, dm::PetscObject) 
get the `DM` associated with the `PETSCVIEWERVTK` viewer

Collective

Input Parameters:
- `viewer` - `PETSCVIEWERVTK` viewer
- `dm`     - `DM` associated with the viewer (as a `PetscObject`)

Level: developer

-seealso: [](sec_viewers), `PETSCVIEWERVTK`, `PetscViewerVTKOpen()`, `DMDAVTKWriteAll()`, `PetscViewerVTKWriteFunction`, `PetscViewerVTKAddField()`

# External Links
$(_doc_external("Sys/PetscViewerVTKGetDM"))
"""
function PetscViewerVTKGetDM(petsclib::PetscLibType, viewer::PetscViewer, dm::PetscObject) end

@for_petsc function PetscViewerVTKGetDM(petsclib::$UnionPetscLib, viewer::PetscViewer, dm::PetscObject )

    @chk ccall(
               (:PetscViewerVTKGetDM, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{PetscObject}),
               viewer, dm,
              )


	return nothing
end 

"""
	PetscViewerVTKOpen(petsclib::PetscLibType,comm::MPI_Comm, name::String, type::PetscFileMode, vtk::PetscViewer) 
Opens a `PETSCVIEWERVTK` viewer file.

Collective

Input Parameters:
- `comm` - MPI communicator
- `name` - name of file
- `type` - type of file
-seealso: [](sec_viewers), `PETSCVIEWERVTK`, `PetscViewerASCIIOpen()`, `PetscViewerPushFormat()`, `PetscViewerDestroy()`,
`VecView()`, `MatView()`, `VecLoad()`, `MatLoad()`,
`PetscFileMode`, `PetscViewer`

# External Links
$(_doc_external("Sys/PetscViewerVTKOpen"))
"""
function PetscViewerVTKOpen(petsclib::PetscLibType, comm::MPI_Comm, name::String, type::PetscFileMode, vtk::PetscViewer) end

@for_petsc function PetscViewerVTKOpen(petsclib::$UnionPetscLib, comm::MPI_Comm, name::String, type::PetscFileMode, vtk::PetscViewer )

    @chk ccall(
               (:PetscViewerVTKOpen, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, PetscFileMode, Ptr{PetscViewer}),
               comm, name, type, vtk,
              )


	return nothing
end 

"""
	PetscViewerMatlabPutArray(petsclib::PetscLibType,mfile::PetscViewer, m::Cint, n::Cint, array::PetscScalar, name::String) 
Puts an array into the `PETSCVIEWERMATLAB` viewer.

Not Collective, only processor zero saves `array`

Input Parameters:
- `mfile` - the viewer
- `m`     - the first dimensions of `array`
- `n`     - the second dimensions of `array`
- `array` - the array (represented in one dimension)
- `name`  - the MATLAB name of `array`

Level: advanced

-seealso: `PETSCVIEWERMATLAB`, `PetscViewerMatlabGetArray()`

# External Links
$(_doc_external("Sys/PetscViewerMatlabPutArray"))
"""
function PetscViewerMatlabPutArray(petsclib::PetscLibType, mfile::PetscViewer, m::Cint, n::Cint, array::PetscScalar, name::String) end

@for_petsc function PetscViewerMatlabPutArray(petsclib::$UnionPetscLib, mfile::PetscViewer, m::Cint, n::Cint, array::$PetscScalar, name::String )

    @chk ccall(
               (:PetscViewerMatlabPutArray, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Cint, Cint, Ptr{$PetscScalar}, Ptr{Cchar}),
               mfile, m, n, array, name,
              )


	return nothing
end 

"""
	PetscViewerMatlabPutVariable(petsclib::PetscLibType,viewer::PetscViewer, name::String, mat::Cvoid) 

# External Links
$(_doc_external("Sys/PetscViewerMatlabPutVariable"))
"""
function PetscViewerMatlabPutVariable(petsclib::PetscLibType, viewer::PetscViewer, name::String, mat::Cvoid) end

@for_petsc function PetscViewerMatlabPutVariable(petsclib::$UnionPetscLib, viewer::PetscViewer, name::String, mat::Cvoid )

    @chk ccall(
               (:PetscViewerMatlabPutVariable, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cchar}, Ptr{Cvoid}),
               viewer, name, mat,
              )


	return nothing
end 

"""
	PetscViewerMatlabGetArray(petsclib::PetscLibType,mfile::PetscViewer, m::Cint, n::Cint, array::Vector{PetscScalar}, name::String) 
Gets a variable from a `PETSCVIEWERMATLAB` viewer into an array

Not Collective; only processor zero reads in the array

Input Parameters:
- `mfile` - the MATLAB file viewer
- `m`     - the first dimensions of `array`
- `n`     - the second dimensions of `array`
- `array` - the array (represented in one dimension), must of be length `m` * `n`
- `name`  - the MATLAB name of `array`

Level: advanced

-seealso: `PETSCVIEWERMATLAB`, `PetscViewerMatlabPutArray()`

# External Links
$(_doc_external("Sys/PetscViewerMatlabGetArray"))
"""
function PetscViewerMatlabGetArray(petsclib::PetscLibType, mfile::PetscViewer, m::Cint, n::Cint, array::Vector{PetscScalar}, name::String) end

@for_petsc function PetscViewerMatlabGetArray(petsclib::$UnionPetscLib, mfile::PetscViewer, m::Cint, n::Cint, array::Vector{$PetscScalar}, name::String )

    @chk ccall(
               (:PetscViewerMatlabGetArray, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Cint, Cint, Ptr{$PetscScalar}, Ptr{Cchar}),
               mfile, m, n, array, name,
              )


	return nothing
end 

"""
	PetscViewerMatlabOpen(petsclib::PetscLibType,comm::MPI_Comm, name::String, type::PetscFileMode, binv::PetscViewer) 
Opens a MATLAB .mat file for output

Collective

Input Parameters:
- `comm` - MPI communicator
- `name` - name of file
- `type` - type of file
-seealso: `PETSCVIEWERMATLAB`, `PetscViewerASCIIOpen()`, `PetscViewerPushFormat()`, `PetscViewerDestroy()`, `PETSCVIEWERBINARY`, `PetscViewerBinaryOpen()`
`VecView()`, `MatView()`, `VecLoad()`, `MatLoad()`

# External Links
$(_doc_external("Sys/PetscViewerMatlabOpen"))
"""
function PetscViewerMatlabOpen(petsclib::PetscLibType, comm::MPI_Comm, name::String, type::PetscFileMode, binv::PetscViewer) end

@for_petsc function PetscViewerMatlabOpen(petsclib::$UnionPetscLib, comm::MPI_Comm, name::String, type::PetscFileMode, binv::PetscViewer )

    @chk ccall(
               (:PetscViewerMatlabOpen, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, PetscFileMode, Ptr{PetscViewer}),
               comm, name, type, binv,
              )


	return nothing
end 

"""
	PetscViewerGLVisSetPrecision(petsclib::PetscLibType,viewer::PetscViewer, prec::PetscInt) 
Set the number of digits for floating point values to be displayed

Not Collective

Input Parameters:
- `viewer` - the `PetscViewer` of type `PETSCVIEWERGLVIS`
- `prec`   - the number of digits required

Level: beginner

-seealso: [](sec_viewers), `PETSCVIEWERGLVIS`, `PetscViewerGLVisOpen()`, `PetscViewerGLVisSetFields()`, `PetscViewerCreate()`, `PetscViewerSetType()`

# External Links
$(_doc_external("Sys/PetscViewerGLVisSetPrecision"))
"""
function PetscViewerGLVisSetPrecision(petsclib::PetscLibType, viewer::PetscViewer, prec::PetscInt) end

@for_petsc function PetscViewerGLVisSetPrecision(petsclib::$UnionPetscLib, viewer::PetscViewer, prec::$PetscInt )

    @chk ccall(
               (:PetscViewerGLVisSetPrecision, $petsc_library),
               PetscErrorCode,
               (PetscViewer, $PetscInt),
               viewer, prec,
              )


	return nothing
end 

"""
	PetscViewerGLVisSetSnapId(petsclib::PetscLibType,viewer::PetscViewer, id::PetscInt) 
Set the snapshot id. Only relevant when the `PetscViewerGLVisType` is `PETSC_VIEWER_GLVIS_DUMP`

Logically Collective

Input Parameters:
- `viewer` - the `PetscViewer` of type `PETSCVIEWERGLVIS`
- `id`     - the current snapshot id in a time-dependent simulation

Level: beginner

-seealso: [](sec_viewers), `PETSCVIEWERGLVIS`, `PetscViewerGLVisOpen()`, `PetscViewerGLVisSetFields()`, `PetscViewerCreate()`, `PetscViewerSetType()`

# External Links
$(_doc_external("Sys/PetscViewerGLVisSetSnapId"))
"""
function PetscViewerGLVisSetSnapId(petsclib::PetscLibType, viewer::PetscViewer, id::PetscInt) end

@for_petsc function PetscViewerGLVisSetSnapId(petsclib::$UnionPetscLib, viewer::PetscViewer, id::$PetscInt )

    @chk ccall(
               (:PetscViewerGLVisSetSnapId, $petsc_library),
               PetscErrorCode,
               (PetscViewer, $PetscInt),
               viewer, id,
              )


	return nothing
end 

"""
	PetscViewerGLVisSetFields(petsclib::PetscLibType,viewer::PetscViewer, nf::PetscInt, fec_type::String, dim::Vector{PetscInt}, g2l::external, Vfield::Vector{PetscObject}, ctx::Cvoid, destroyctx::external) 
Sets the required information to visualize different fields from a vector.

Logically Collective

Input Parameters:
- `viewer`     - the `PetscViewer` of type `PETSCVIEWERGLVIS`
- `nf`         - number of fields to be visualized
- `fec_type`   - the type of finite element to be used to visualize the data (see FiniteElementCollection::Name() in MFEM)
- `dim`        - array of space dimension for field vectors (used to initialize the scene)
- `g2l`        - User routine to compute the local field vectors to be visualized; PetscObject is used in place of Vec on the prototype
- `Vfield`     - array of work vectors, one for each field
- `ctx`        - User context to store the relevant data to apply g2lfields
- `destroyctx` - Destroy function for userctx

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERGLVIS`, `PetscViewerGLVisOpen()`, `PetscViewerCreate()`, `PetscViewerSetType()`, `PetscObjectSetName()`

# External Links
$(_doc_external("Sys/PetscViewerGLVisSetFields"))
"""
function PetscViewerGLVisSetFields(petsclib::PetscLibType, viewer::PetscViewer, nf::PetscInt, fec_type::String, dim::Vector{PetscInt}, g2l::external, Vfield::Vector{PetscObject}, ctx::Cvoid, destroyctx::external) end

@for_petsc function PetscViewerGLVisSetFields(petsclib::$UnionPetscLib, viewer::PetscViewer, nf::$PetscInt, fec_type::String, dim::Vector{$PetscInt}, g2l::external, Vfield::Vector{PetscObject}, ctx::Cvoid, destroyctx::external )
	fec_type_ = Ref(pointer(fec_type))

    @chk ccall(
               (:PetscViewerGLVisSetFields, $petsc_library),
               PetscErrorCode,
               (PetscViewer, $PetscInt, Ptr{Ptr{Cchar}}, Ptr{$PetscInt}, external, Ptr{PetscObject}, Ptr{Cvoid}, external),
               viewer, nf, fec_type_, dim, g2l, Vfield, ctx, destroyctx,
              )


	return nothing
end 

"""
	PetscViewerGLVisOpen(petsclib::PetscLibType,comm::MPI_Comm, type::PetscViewerGLVisType, name::String, port::PetscInt, viewer::PetscViewer) 
Opens a `PETSCVIEWERGLVIS` `PetscViewer`

Collective; No Fortran Support

Input Parameters:
- `comm` - the MPI communicator
- `type` - the viewer type: `PETSC_VIEWER_GLVIS_SOCKET` for real-time visualization or `PETSC_VIEWER_GLVIS_DUMP` for dumping to a file
- `name` - either the hostname where the GLVis server is running or the base filename for dumping the data for subsequent visualizations
- `port` - socket port where the GLVis server is listening. Not referenced when type is `PETSC_VIEWER_GLVIS_DUMP`

Output Parameter:
- `viewer` - the `PetscViewer` object

Options Database Keys:
- `-glvis_precision <precision>` - Sets number of digits for floating point values
- `-glvis_size <width,height>`   - Sets the window size (in pixels)
- `-glvis_pause <pause>`         - Sets time (in seconds) that the program pauses after each visualization
(0 is default, -1 implies every visualization)
- `-glvis_keys`                  - Additional keys to configure visualization
- `-glvis_exec`                  - Additional commands to configure visualization

Level: beginner

-seealso: [](sec_viewers), `PETSCVIEWERGLVIS`, `PetscViewerCreate()`, `PetscViewerSetType()`, `PetscViewerGLVisType`

# External Links
$(_doc_external("Sys/PetscViewerGLVisOpen"))
"""
function PetscViewerGLVisOpen(petsclib::PetscLibType, comm::MPI_Comm, type::PetscViewerGLVisType, name::String, port::PetscInt, viewer::PetscViewer) end

@for_petsc function PetscViewerGLVisOpen(petsclib::$UnionPetscLib, comm::MPI_Comm, type::PetscViewerGLVisType, name::String, port::$PetscInt, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerGLVisOpen, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, PetscViewerGLVisType, Ptr{Cchar}, $PetscInt, Ptr{PetscViewer}),
               comm, type, name, port, viewer,
              )


	return nothing
end 

"""
	PetscViewerBinaryGetMPIIOOffset(petsclib::PetscLibType,viewer::PetscViewer, off::MPI_Offset) 
Gets the current global offset that should be passed to `MPI_File_set_view()` or `MPI_File_{write|read}_at[_all]()`

Not Collective; No Fortran Support

Input Parameter:
- `viewer` - `PetscViewer` context, obtained from `PetscViewerBinaryOpen()`

Output Parameter:
- `off` - the current global offset

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERBINARY`, `PetscViewerBinaryOpen()`, `PetscViewerBinaryGetInfoPointer()`, `PetscViewerBinaryGetUseMPIIO()`, `PetscViewerBinarySetUseMPIIO()`, `PetscViewerBinaryAddMPIIOOffset()`

# External Links
$(_doc_external("Sys/PetscViewerBinaryGetMPIIOOffset"))
"""
function PetscViewerBinaryGetMPIIOOffset(petsclib::PetscLibType, viewer::PetscViewer, off::MPI_Offset) end

@for_petsc function PetscViewerBinaryGetMPIIOOffset(petsclib::$UnionPetscLib, viewer::PetscViewer, off::MPI_Offset )

    @chk ccall(
               (:PetscViewerBinaryGetMPIIOOffset, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{MPI_Offset}),
               viewer, off,
              )


	return nothing
end 

"""
	PetscViewerBinaryAddMPIIOOffset(petsclib::PetscLibType,viewer::PetscViewer, off::MPI_Offset) 
Adds to the current global offset

Logically Collective; No Fortran Support

Input Parameters:
- `viewer` - `PetscViewer` context, obtained from `PetscViewerBinaryOpen()`
- `off`    - the addition to the global offset

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERBINARY`, `PetscViewerBinaryOpen()`, `PetscViewerBinaryGetInfoPointer()`, `PetscViewerBinaryGetUseMPIIO()`, `PetscViewerBinarySetUseMPIIO()`, `PetscViewerBinaryGetMPIIOOffset()`

# External Links
$(_doc_external("Sys/PetscViewerBinaryAddMPIIOOffset"))
"""
function PetscViewerBinaryAddMPIIOOffset(petsclib::PetscLibType, viewer::PetscViewer, off::MPI_Offset) end

@for_petsc function PetscViewerBinaryAddMPIIOOffset(petsclib::$UnionPetscLib, viewer::PetscViewer, off::MPI_Offset )

    @chk ccall(
               (:PetscViewerBinaryAddMPIIOOffset, $petsc_library),
               PetscErrorCode,
               (PetscViewer, MPI_Offset),
               viewer, off,
              )


	return nothing
end 

"""
	PetscViewerBinaryGetMPIIODescriptor(petsclib::PetscLibType,viewer::PetscViewer, fdes::MPI_File) 
Extracts the MPI IO file descriptor from a `PetscViewer`.

Not Collective; No Fortran Support

Input Parameter:
- `viewer` - `PetscViewer` context, obtained from `PetscViewerBinaryOpen()`

Output Parameter:
- `fdes` - file descriptor

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERBINARY`, `PetscViewerBinaryOpen()`, `PetscViewerBinaryGetInfoPointer()`, `PetscViewerBinaryGetUseMPIIO()`, `PetscViewerBinarySetUseMPIIO()`, `PetscViewerBinaryGetMPIIOOffset()`

# External Links
$(_doc_external("Sys/PetscViewerBinaryGetMPIIODescriptor"))
"""
function PetscViewerBinaryGetMPIIODescriptor(petsclib::PetscLibType, viewer::PetscViewer, fdes::MPI_File) end

@for_petsc function PetscViewerBinaryGetMPIIODescriptor(petsclib::$UnionPetscLib, viewer::PetscViewer, fdes::MPI_File )

    @chk ccall(
               (:PetscViewerBinaryGetMPIIODescriptor, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{MPI_File}),
               viewer, fdes,
              )


	return nothing
end 

"""
	PetscViewerBinarySetUseMPIIO(petsclib::PetscLibType,viewer::PetscViewer, use::PetscBool) 
Sets a binary viewer to use MPI
before `PetscViewerFileSetName()`

Logically Collective

Input Parameters:
- `viewer` - the `PetscViewer`; must be a `PETSCVIEWERBINARY`
- `use`    - `PETSC_TRUE` means MPI-IO will be used

Options Database Key:
- `-viewer_binary_mpiio` - <true or false> flag for using MPI-IO

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERBINARY`, `PetscViewerFileSetMode()`, `PetscViewerCreate()`, `PetscViewerSetType()`, `PetscViewerBinaryOpen()`,
`PetscViewerBinaryGetUseMPIIO()`

# External Links
$(_doc_external("Sys/PetscViewerBinarySetUseMPIIO"))
"""
function PetscViewerBinarySetUseMPIIO(petsclib::PetscLibType, viewer::PetscViewer, use::PetscBool) end

@for_petsc function PetscViewerBinarySetUseMPIIO(petsclib::$UnionPetscLib, viewer::PetscViewer, use::PetscBool )

    @chk ccall(
               (:PetscViewerBinarySetUseMPIIO, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscBool),
               viewer, use,
              )


	return nothing
end 

"""
	use::PetscBool = PetscViewerBinaryGetUseMPIIO(petsclib::PetscLibType,viewer::PetscViewer) 
Returns `PETSC_TRUE` if the binary viewer uses MPI

Not Collective

Input Parameter:
- `viewer` - `PetscViewer` context, obtained from `PetscViewerBinaryOpen()`; must be a `PETSCVIEWERBINARY`

Output Parameter:
- `use` - `PETSC_TRUE` if MPI-IO is being used

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERBINARY`, `PetscViewerBinaryOpen()`, `PetscViewerBinaryGetInfoPointer()`, `PetscViewerBinarySetUseMPIIO()`, `PetscViewerBinaryGetMPIIOOffset()`

# External Links
$(_doc_external("Sys/PetscViewerBinaryGetUseMPIIO"))
"""
function PetscViewerBinaryGetUseMPIIO(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerBinaryGetUseMPIIO(petsclib::$UnionPetscLib, viewer::PetscViewer )
	use_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscViewerBinaryGetUseMPIIO, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{PetscBool}),
               viewer, use_,
              )

	use = use_[]

	return use
end 

"""
	PetscViewerBinarySetFlowControl(petsclib::PetscLibType,viewer::PetscViewer, fc::PetscInt) 
Sets how many messages are allowed to be outstanding at the same time during parallel IO reads/writes

Not Collective

Input Parameters:
- `viewer` - `PetscViewer` context, obtained from `PetscViewerBinaryOpen()`
- `fc`     - the number of messages, defaults to 256 if this function was not called

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERBINARY`, `PetscViewerBinaryOpen()`, `PetscViewerBinaryGetInfoPointer()`, `PetscViewerBinaryGetFlowControl()`

# External Links
$(_doc_external("Sys/PetscViewerBinarySetFlowControl"))
"""
function PetscViewerBinarySetFlowControl(petsclib::PetscLibType, viewer::PetscViewer, fc::PetscInt) end

@for_petsc function PetscViewerBinarySetFlowControl(petsclib::$UnionPetscLib, viewer::PetscViewer, fc::$PetscInt )

    @chk ccall(
               (:PetscViewerBinarySetFlowControl, $petsc_library),
               PetscErrorCode,
               (PetscViewer, $PetscInt),
               viewer, fc,
              )


	return nothing
end 

"""
	fc::PetscInt = PetscViewerBinaryGetFlowControl(petsclib::PetscLibType,viewer::PetscViewer) 
Returns how many messages are allowed to be outstanding at the same time during parallel IO reads/writes

Not Collective

Input Parameter:
- `viewer` - `PetscViewer` context, obtained from `PetscViewerBinaryOpen()`

Output Parameter:
- `fc` - the number of messages

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERBINARY`, `PetscViewerBinaryOpen()`, `PetscViewerBinaryGetInfoPointer()`, `PetscViewerBinarySetFlowControl()`

# External Links
$(_doc_external("Sys/PetscViewerBinaryGetFlowControl"))
"""
function PetscViewerBinaryGetFlowControl(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerBinaryGetFlowControl(petsclib::$UnionPetscLib, viewer::PetscViewer )
	fc_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscViewerBinaryGetFlowControl, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{$PetscInt}),
               viewer, fc_,
              )

	fc = fc_[]

	return fc
end 

"""
	PetscViewerBinaryGetDescriptor(petsclib::PetscLibType,viewer::PetscViewer, fdes::Cint) 
Extracts the file descriptor from a `PetscViewer` of `PetscViewerType` `PETSCVIEWERBINARY`.

Collective because it may trigger a `PetscViewerSetUp()` call; No Fortran Support

Input Parameter:
- `viewer` - `PetscViewer` context, obtained from `PetscViewerBinaryOpen()`

Output Parameter:
- `fdes` - file descriptor

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERBINARY`, `PetscViewerBinaryOpen()`, `PetscViewerBinaryGetInfoPointer()`

# External Links
$(_doc_external("Sys/PetscViewerBinaryGetDescriptor"))
"""
function PetscViewerBinaryGetDescriptor(petsclib::PetscLibType, viewer::PetscViewer, fdes::Cint) end

@for_petsc function PetscViewerBinaryGetDescriptor(petsclib::$UnionPetscLib, viewer::PetscViewer, fdes::Cint )

    @chk ccall(
               (:PetscViewerBinaryGetDescriptor, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cint}),
               viewer, fdes,
              )


	return nothing
end 

"""
	PetscViewerBinarySkipInfo(petsclib::PetscLibType,viewer::PetscViewer) 
Binary file will not have `.info` file created with it

Not Collective

Input Parameter:
- `viewer` - `PetscViewer` context, obtained from `PetscViewerCreate()`

Options Database Key:
- `-viewer_binary_skip_info` - true indicates do not generate `.info` file

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERBINARY`, `PetscViewerBinaryOpen()`, `PetscViewerBinaryGetDescriptor()`, `PetscViewerBinarySetSkipOptions()`,
`PetscViewerBinaryGetSkipOptions()`, `PetscViewerBinaryGetSkipInfo()`

# External Links
$(_doc_external("Sys/PetscViewerBinarySkipInfo"))
"""
function PetscViewerBinarySkipInfo(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerBinarySkipInfo(petsclib::$UnionPetscLib, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerBinarySkipInfo, $petsc_library),
               PetscErrorCode,
               (PetscViewer,),
               viewer,
              )


	return nothing
end 

"""
	PetscViewerBinarySetSkipInfo(petsclib::PetscLibType,viewer::PetscViewer, skip::PetscBool) 
Binary file will not have `.info` file created with it

Not Collective

Input Parameters:
- `viewer` - PetscViewer context, obtained from `PetscViewerCreate()`
- `skip`   - `PETSC_TRUE` implies the `.info` file will not be generated

Options Database Key:
- `-viewer_binary_skip_info` - true indicates do not generate `.info` file

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERBINARY`, `PetscViewerBinaryOpen()`, `PetscViewerBinaryGetDescriptor()`, `PetscViewerBinarySetSkipOptions()`,
`PetscViewerBinaryGetSkipOptions()`, `PetscViewerBinaryGetSkipInfo()`, `PetscViewerBinaryGetInfoPointer()`

# External Links
$(_doc_external("Sys/PetscViewerBinarySetSkipInfo"))
"""
function PetscViewerBinarySetSkipInfo(petsclib::PetscLibType, viewer::PetscViewer, skip::PetscBool) end

@for_petsc function PetscViewerBinarySetSkipInfo(petsclib::$UnionPetscLib, viewer::PetscViewer, skip::PetscBool )

    @chk ccall(
               (:PetscViewerBinarySetSkipInfo, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscBool),
               viewer, skip,
              )


	return nothing
end 

"""
	skip::PetscBool = PetscViewerBinaryGetSkipInfo(petsclib::PetscLibType,viewer::PetscViewer) 
check if viewer wrote a `.info` file

Not Collective

Input Parameter:
- `viewer` - `PetscViewer` context, obtained from `PetscViewerBinaryOpen()`

Output Parameter:
- `skip` - `PETSC_TRUE` implies the `.info` file was not generated

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERBINARY`, `PetscViewerBinaryOpen()`, `PetscViewerBinaryGetDescriptor()`,
`PetscViewerBinarySetSkipOptions()`, `PetscViewerBinarySetSkipInfo()`, `PetscViewerBinaryGetInfoPointer()`

# External Links
$(_doc_external("Sys/PetscViewerBinaryGetSkipInfo"))
"""
function PetscViewerBinaryGetSkipInfo(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerBinaryGetSkipInfo(petsclib::$UnionPetscLib, viewer::PetscViewer )
	skip_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscViewerBinaryGetSkipInfo, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{PetscBool}),
               viewer, skip_,
              )

	skip = skip_[]

	return skip
end 

"""
	PetscViewerBinarySetSkipOptions(petsclib::PetscLibType,viewer::PetscViewer, skip::PetscBool) 
do not use values in the PETSc options database when loading objects

Not Collective

Input Parameters:
- `viewer` - `PetscViewer` context, obtained from `PetscViewerBinaryOpen()`
- `skip`   - `PETSC_TRUE` means do not use the options from the options database

Options Database Key:
- `-viewer_binary_skip_options <true or false>` - true means do not use the options from the options database

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERBINARY`, `PetscViewerBinaryOpen()`, `PetscViewerBinaryGetDescriptor()`, `PetscViewerBinarySkipInfo()`,
`PetscViewerBinaryGetSkipOptions()`

# External Links
$(_doc_external("Sys/PetscViewerBinarySetSkipOptions"))
"""
function PetscViewerBinarySetSkipOptions(petsclib::PetscLibType, viewer::PetscViewer, skip::PetscBool) end

@for_petsc function PetscViewerBinarySetSkipOptions(petsclib::$UnionPetscLib, viewer::PetscViewer, skip::PetscBool )

    @chk ccall(
               (:PetscViewerBinarySetSkipOptions, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscBool),
               viewer, skip,
              )


	return nothing
end 

"""
	skip::PetscBool = PetscViewerBinaryGetSkipOptions(petsclib::PetscLibType,viewer::PetscViewer) 
checks if viewer uses the PETSc options database when loading objects

Not Collective

Input Parameter:
- `viewer` - `PetscViewer` context, obtained from `PetscViewerBinaryOpen()`

Output Parameter:
- `skip` - `PETSC_TRUE` means do not use

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERBINARY`, `PetscViewerBinaryOpen()`, `PetscViewerBinaryGetDescriptor()`, `PetscViewerBinarySkipInfo()`,
`PetscViewerBinarySetSkipOptions()`

# External Links
$(_doc_external("Sys/PetscViewerBinaryGetSkipOptions"))
"""
function PetscViewerBinaryGetSkipOptions(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerBinaryGetSkipOptions(petsclib::$UnionPetscLib, viewer::PetscViewer )
	skip_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscViewerBinaryGetSkipOptions, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{PetscBool}),
               viewer, skip_,
              )

	skip = skip_[]

	return skip
end 

"""
	PetscViewerBinarySetSkipHeader(petsclib::PetscLibType,viewer::PetscViewer, skip::PetscBool) 
do not write a header with size information on output, just raw data

Not Collective

Input Parameters:
- `viewer` - `PetscViewer` context, obtained from `PetscViewerBinaryOpen()`
- `skip`   - `PETSC_TRUE` means do not write header

Options Database Key:
- `-viewer_binary_skip_header <true or false>` - true means do not write header

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERBINARY`, `PetscViewerBinaryOpen()`, `PetscViewerBinaryGetDescriptor()`, `PetscViewerBinarySkipInfo()`,
`PetscViewerBinaryGetSkipHeader()`

# External Links
$(_doc_external("Sys/PetscViewerBinarySetSkipHeader"))
"""
function PetscViewerBinarySetSkipHeader(petsclib::PetscLibType, viewer::PetscViewer, skip::PetscBool) end

@for_petsc function PetscViewerBinarySetSkipHeader(petsclib::$UnionPetscLib, viewer::PetscViewer, skip::PetscBool )

    @chk ccall(
               (:PetscViewerBinarySetSkipHeader, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscBool),
               viewer, skip,
              )


	return nothing
end 

"""
	skip::PetscBool = PetscViewerBinaryGetSkipHeader(petsclib::PetscLibType,viewer::PetscViewer) 
checks whether to write a header with size information on output, or just raw data

Not Collective

Input Parameter:
- `viewer` - `PetscViewer` context, obtained from `PetscViewerBinaryOpen()`

Output Parameter:
- `skip` - `PETSC_TRUE` means do not write header

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERBINARY`, `PetscViewerBinaryOpen()`, `PetscViewerBinaryGetDescriptor()`, `PetscViewerBinarySkipInfo()`,
`PetscViewerBinarySetSkipHeader()`

# External Links
$(_doc_external("Sys/PetscViewerBinaryGetSkipHeader"))
"""
function PetscViewerBinaryGetSkipHeader(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerBinaryGetSkipHeader(petsclib::$UnionPetscLib, viewer::PetscViewer )
	skip_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscViewerBinaryGetSkipHeader, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{PetscBool}),
               viewer, skip_,
              )

	skip = skip_[]

	return skip
end 

"""
	PetscViewerBinaryGetInfoPointer(petsclib::PetscLibType,viewer::PetscViewer, file::Libc.FILE) 
Extracts the file pointer for the ASCII
`.info` file associated with a binary file.

Not Collective; No Fortran Support

Input Parameter:
- `viewer` - `PetscViewer` context, obtained from `PetscViewerBinaryOpen()`

Output Parameter:
- `file` - file pointer  Always returns `NULL` if not a binary viewer

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERBINARY`, `PetscViewerBinaryOpen()`, `PetscViewerBinaryGetDescriptor()`, `PetscViewerBinaryGetSkipInfo()`,
`PetscViewerBinarySetSkipInfo()`

# External Links
$(_doc_external("Sys/PetscViewerBinaryGetInfoPointer"))
"""
function PetscViewerBinaryGetInfoPointer(petsclib::PetscLibType, viewer::PetscViewer, file::Libc.FILE) end

@for_petsc function PetscViewerBinaryGetInfoPointer(petsclib::$UnionPetscLib, viewer::PetscViewer, file::Libc.FILE )

    @chk ccall(
               (:PetscViewerBinaryGetInfoPointer, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Libc.FILE),
               viewer, file,
              )


	return nothing
end 

"""
	PetscViewerBinaryOpen(petsclib::PetscLibType,comm::MPI_Comm, name::String, mode::PetscFileMode, viewer::PetscViewer) 
Opens a file for binary input/output.

Collective

Input Parameters:
- `comm` - MPI communicator
- `name` - name of file
- `mode` - open mode of file
-seealso: [](sec_viewers), `PETSCVIEWERBINARY`, `PetscViewerASCIIOpen()`, `PetscViewerPushFormat()`, `PetscViewerDestroy()`,
`VecView()`, `MatView()`, `VecLoad()`, `MatLoad()`, `PetscViewerBinaryGetDescriptor()`,
`PetscViewerBinaryGetInfoPointer()`, `PetscFileMode`, `PetscViewer`, `PetscViewerBinaryRead()`, `PetscViewerBinarySetUseMPIIO()`,
`PetscViewerBinaryGetUseMPIIO()`, `PetscViewerBinaryGetMPIIOOffset()`

# External Links
$(_doc_external("Sys/PetscViewerBinaryOpen"))
"""
function PetscViewerBinaryOpen(petsclib::PetscLibType, comm::MPI_Comm, name::String, mode::PetscFileMode, viewer::PetscViewer) end

@for_petsc function PetscViewerBinaryOpen(petsclib::$UnionPetscLib, comm::MPI_Comm, name::String, mode::PetscFileMode, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerBinaryOpen, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, PetscFileMode, Ptr{PetscViewer}),
               comm, name, mode, viewer,
              )


	return nothing
end 

"""
	count::PetscInt = PetscViewerBinaryRead(petsclib::PetscLibType,viewer::PetscViewer, data::Cvoid, num::PetscInt, dtype::PetscDataType) 
Reads from a binary file, all processors get the same result

Collective; No Fortran Support

Input Parameters:
- `viewer` - the `PETSCVIEWERBINARY` viewer
- `num`    - number of items of data to read
- `dtype`  - type of data to read

Output Parameters:
- `data`  - location of the read data, treated as an array of the type indicated by `dtype`
- `count` - number of items of data actually read, or `NULL`.

Level: beginner

-seealso: [](sec_viewers), `PETSCVIEWERBINARY`, `PetscViewerASCIIOpen()`, `PetscViewerPushFormat()`, `PetscViewerDestroy()`,
`VecView()`, `MatView()`, `VecLoad()`, `MatLoad()`, `PetscViewerBinaryGetDescriptor()`,
`PetscViewerBinaryGetInfoPointer()`, `PetscFileMode`, `PetscViewer`

# External Links
$(_doc_external("Sys/PetscViewerBinaryRead"))
"""
function PetscViewerBinaryRead(petsclib::PetscLibType, viewer::PetscViewer, data::Cvoid, num::PetscInt, dtype::PetscDataType) end

@for_petsc function PetscViewerBinaryRead(petsclib::$UnionPetscLib, viewer::PetscViewer, data::Cvoid, num::$PetscInt, dtype::PetscDataType )
	count_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscViewerBinaryRead, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cvoid}, $PetscInt, Ptr{$PetscInt}, PetscDataType),
               viewer, data, num, count_, dtype,
              )

	count = count_[]

	return count
end 

"""
	PetscViewerBinaryWrite(petsclib::PetscLibType,viewer::PetscViewer, data::Cvoid, count::PetscInt, dtype::PetscDataType) 
writes to a binary file, only from the first MPI rank

Collective; No Fortran Support

Input Parameters:
- `viewer` - the `PETSCVIEWERBINARY` viewer
- `data`   - location of data, treated as an array of the type indicated by `dtype`
- `count`  - number of items of data to write
- `dtype`  - type of data to write

Level: beginner

-seealso: [](sec_viewers), `PETSCVIEWERBINARY`, `PetscViewerASCIIOpen()`, `PetscViewerPushFormat()`, `PetscViewerDestroy()`,
`VecView()`, `MatView()`, `VecLoad()`, `MatLoad()`, `PetscViewerBinaryGetDescriptor()`, `PetscDataType`
`PetscViewerBinaryGetInfoPointer()`, `PetscFileMode`, `PetscViewer`, `PetscViewerBinaryRead()`

# External Links
$(_doc_external("Sys/PetscViewerBinaryWrite"))
"""
function PetscViewerBinaryWrite(petsclib::PetscLibType, viewer::PetscViewer, data::Cvoid, count::PetscInt, dtype::PetscDataType) end

@for_petsc function PetscViewerBinaryWrite(petsclib::$UnionPetscLib, viewer::PetscViewer, data::Cvoid, count::$PetscInt, dtype::PetscDataType )

    @chk ccall(
               (:PetscViewerBinaryWrite, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cvoid}, $PetscInt, PetscDataType),
               viewer, data, count, dtype,
              )


	return nothing
end 

"""
	PetscViewerBinaryReadAll(petsclib::PetscLibType,viewer::PetscViewer, data::Cvoid, count::PetscCount, start::PetscCount, total::PetscCount, dtype::PetscDataType) 
reads from a binary file from all MPI processes, each rank receives its own portion of the data

Collective; No Fortran Support

Input Parameters:
- `viewer` - the `PETSCVIEWERBINARY` viewer
- `count`  - local number of items of data to read
- `start`  - local start, can be `PETSC_DETERMINE`
- `total`  - global number of items of data to read, can be `PETSC_DETERMINE`
- `dtype`  - type of data to read

Output Parameter:
- `data` - location of data, treated as an array of type indicated by `dtype`

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERBINARY`, `PetscViewerBinaryOpen()`, `PetscViewerBinarySetUseMPIIO()`, `PetscViewerBinaryRead()`, `PetscViewerBinaryWriteAll()`

# External Links
$(_doc_external("Sys/PetscViewerBinaryReadAll"))
"""
function PetscViewerBinaryReadAll(petsclib::PetscLibType, viewer::PetscViewer, data::Cvoid, count::PetscCount, start::PetscCount, total::PetscCount, dtype::PetscDataType) end

@for_petsc function PetscViewerBinaryReadAll(petsclib::$UnionPetscLib, viewer::PetscViewer, data::Cvoid, count::PetscCount, start::PetscCount, total::PetscCount, dtype::PetscDataType )

    @chk ccall(
               (:PetscViewerBinaryReadAll, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cvoid}, PetscCount, PetscCount, PetscCount, PetscDataType),
               viewer, data, count, start, total, dtype,
              )


	return nothing
end 

"""
	PetscViewerBinaryWriteAll(petsclib::PetscLibType,viewer::PetscViewer, data::Cvoid, count::PetscCount, start::PetscCount, total::PetscCount, dtype::PetscDataType) 
writes to a binary file from all MPI processes, each rank writes its own portion of the data

Collective; No Fortran Support

Input Parameters:
- `viewer` - the `PETSCVIEWERBINARY` viewer
- `data`   - location of data
- `count`  - local number of items of data to write, treated as an array of type indicated by `dtype`
- `start`  - local start, can be `PETSC_DETERMINE`
- `total`  - global number of items of data to write, can be `PETSC_DETERMINE`
- `dtype`  - type of data to write

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERBINARY`, `PetscViewerBinaryOpen()`, `PetscViewerBinarySetUseMPIIO()`, `PetscViewerBinaryReadAll()`

# External Links
$(_doc_external("Sys/PetscViewerBinaryWriteAll"))
"""
function PetscViewerBinaryWriteAll(petsclib::PetscLibType, viewer::PetscViewer, data::Cvoid, count::PetscCount, start::PetscCount, total::PetscCount, dtype::PetscDataType) end

@for_petsc function PetscViewerBinaryWriteAll(petsclib::$UnionPetscLib, viewer::PetscViewer, data::Cvoid, count::PetscCount, start::PetscCount, total::PetscCount, dtype::PetscDataType )

    @chk ccall(
               (:PetscViewerBinaryWriteAll, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cvoid}, PetscCount, PetscCount, PetscCount, PetscDataType),
               viewer, data, count, start, total, dtype,
              )


	return nothing
end 

"""
	PetscViewerBinaryWriteStringArray(petsclib::PetscLibType,viewer::PetscViewer, data::String) 
writes to a binary file, only from the first MPI rank, an array of strings

Collective; No Fortran Support

Input Parameters:
- `viewer` - the `PETSCVIEWERBINARY` viewer
- `data`   - location of the array of strings

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERBINARY`, `PetscViewerASCIIOpen()`, `PetscViewerPushFormat()`, `PetscViewerDestroy()`,
`VecView()`, `MatView()`, `VecLoad()`, `MatLoad()`, `PetscViewerBinaryGetDescriptor()`,
`PetscViewerBinaryGetInfoPointer()`, `PetscFileMode`, `PetscViewer`, `PetscViewerBinaryRead()`

# External Links
$(_doc_external("Sys/PetscViewerBinaryWriteStringArray"))
"""
function PetscViewerBinaryWriteStringArray(petsclib::PetscLibType, viewer::PetscViewer, data::String) end

@for_petsc function PetscViewerBinaryWriteStringArray(petsclib::$UnionPetscLib, viewer::PetscViewer, data::String )
	data_ = Ref(pointer(data))

    @chk ccall(
               (:PetscViewerBinaryWriteStringArray, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Ptr{Cchar}}),
               viewer, data_,
              )


	return nothing
end 

"""
	PetscViewerBinaryReadStringArray(petsclib::PetscLibType,viewer::PetscViewer, data::Cchar) 
reads a binary file an array of strings to all MPI processes

Collective; No Fortran Support

Input Parameter:
- `viewer` - the `PETSCVIEWERBINARY` viewer

Output Parameter:
- `data` - location of the array of strings

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERBINARY`, `PetscViewerASCIIOpen()`, `PetscViewerPushFormat()`, `PetscViewerDestroy()`,
`VecView()`, `MatView()`, `VecLoad()`, `MatLoad()`, `PetscViewerBinaryGetDescriptor()`,
`PetscViewerBinaryGetInfoPointer()`, `PetscFileMode`, `PetscViewer`, `PetscViewerBinaryRead()`

# External Links
$(_doc_external("Sys/PetscViewerBinaryReadStringArray"))
"""
function PetscViewerBinaryReadStringArray(petsclib::PetscLibType, viewer::PetscViewer, data::Cchar) end

@for_petsc function PetscViewerBinaryReadStringArray(petsclib::$UnionPetscLib, viewer::PetscViewer, data::Cchar )

    @chk ccall(
               (:PetscViewerBinaryReadStringArray, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Cchar),
               viewer, data,
              )


	return nothing
end 

"""
	PetscViewerFileSetMode(petsclib::PetscLibType,viewer::PetscViewer, mode::PetscFileMode) 
Sets the open mode of file

Logically Collective

Input Parameters:
- `viewer` - the `PetscViewer`; must be a `PETSCVIEWERBINARY`, `PETSCVIEWERMATLAB`, `PETSCVIEWERHDF5`, or `PETSCVIEWERASCII`  `PetscViewer`
- `mode`   - open mode of file
-seealso: [](sec_viewers), `PetscViewerCreate()`, `PetscViewerSetType()`, `PetscViewerBinaryOpen()`

# External Links
$(_doc_external("Sys/PetscViewerFileSetMode"))
"""
function PetscViewerFileSetMode(petsclib::PetscLibType, viewer::PetscViewer, mode::PetscFileMode) end

@for_petsc function PetscViewerFileSetMode(petsclib::$UnionPetscLib, viewer::PetscViewer, mode::PetscFileMode )

    @chk ccall(
               (:PetscViewerFileSetMode, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscFileMode),
               viewer, mode,
              )


	return nothing
end 

"""
	PetscViewerFileGetMode(petsclib::PetscLibType,viewer::PetscViewer, mode::PetscFileMode) 
Gets the open mode of a file associated with a `PetscViewer`

Not Collective

Input Parameter:
- `viewer` - the `PetscViewer`; must be a `PETSCVIEWERBINARY`, `PETSCVIEWERMATLAB`, `PETSCVIEWERHDF5`, or `PETSCVIEWERASCII`  `PetscViewer`

Output Parameter:
- `mode` - open mode of file
-seealso: [](sec_viewers), `PetscViewerFileSetMode()`, `PetscViewerCreate()`, `PetscViewerSetType()`, `PetscViewerBinaryOpen()`

# External Links
$(_doc_external("Sys/PetscViewerFileGetMode"))
"""
function PetscViewerFileGetMode(petsclib::PetscLibType, viewer::PetscViewer, mode::PetscFileMode) end

@for_petsc function PetscViewerFileGetMode(petsclib::$UnionPetscLib, viewer::PetscViewer, mode::PetscFileMode )

    @chk ccall(
               (:PetscViewerFileGetMode, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{PetscFileMode}),
               viewer, mode,
              )


	return nothing
end 

"""
	PetscViewerDrawBaseAdd(petsclib::PetscLibType,viewer::PetscViewer, windownumber::PetscInt) 
add to the base integer that is added to the `windownumber` passed to `PetscViewerDrawGetDraw()`

Logically Collective

Input Parameters:
- `viewer`       - the `PetscViewer` (created with `PetscViewerDrawOpen()`)
- `windownumber` - how much to add to the base

Level: developer

-seealso: [](sec_viewers), `PetscViewerDrawGetLG()`, `PetscViewerDrawGetAxis()`, `PetscViewerDrawOpen()`, `PetscViewerDrawGetDraw()`, `PetscViewerDrawBaseSet()`

# External Links
$(_doc_external("Sys/PetscViewerDrawBaseAdd"))
"""
function PetscViewerDrawBaseAdd(petsclib::PetscLibType, viewer::PetscViewer, windownumber::PetscInt) end

@for_petsc function PetscViewerDrawBaseAdd(petsclib::$UnionPetscLib, viewer::PetscViewer, windownumber::$PetscInt )

    @chk ccall(
               (:PetscViewerDrawBaseAdd, $petsc_library),
               PetscErrorCode,
               (PetscViewer, $PetscInt),
               viewer, windownumber,
              )


	return nothing
end 

"""
	PetscViewerDrawBaseSet(petsclib::PetscLibType,viewer::PetscViewer, windownumber::PetscInt) 
sets the base integer that is added to the `windownumber` passed to `PetscViewerDrawGetDraw()`

Logically Collective

Input Parameters:
- `viewer`       - the `PetscViewer` (created with `PetscViewerDrawOpen()`)
- `windownumber` - value to set the base

Level: developer

-seealso: [](sec_viewers), `PetscViewerDrawGetLG()`, `PetscViewerDrawGetAxis()`, `PetscViewerDrawOpen()`, `PetscViewerDrawGetDraw()`, `PetscViewerDrawBaseAdd()`

# External Links
$(_doc_external("Sys/PetscViewerDrawBaseSet"))
"""
function PetscViewerDrawBaseSet(petsclib::PetscLibType, viewer::PetscViewer, windownumber::PetscInt) end

@for_petsc function PetscViewerDrawBaseSet(petsclib::$UnionPetscLib, viewer::PetscViewer, windownumber::$PetscInt )

    @chk ccall(
               (:PetscViewerDrawBaseSet, $petsc_library),
               PetscErrorCode,
               (PetscViewer, $PetscInt),
               viewer, windownumber,
              )


	return nothing
end 

"""
	PetscViewerDrawResize(petsclib::PetscLibType,v::PetscViewer, w::Cint, h::Cint) 

# External Links
$(_doc_external("Sys/PetscViewerDrawResize"))
"""
function PetscViewerDrawResize(petsclib::PetscLibType, v::PetscViewer, w::Cint, h::Cint) end

@for_petsc function PetscViewerDrawResize(petsclib::$UnionPetscLib, v::PetscViewer, w::Cint, h::Cint )

    @chk ccall(
               (:PetscViewerDrawResize, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Cint, Cint),
               v, w, h,
              )


	return nothing
end 

"""
	PetscViewerDrawSetInfo(petsclib::PetscLibType,v::PetscViewer, display::String, title::String, x::Cint, y::Cint, w::Cint, h::Cint) 

# External Links
$(_doc_external("Sys/PetscViewerDrawSetInfo"))
"""
function PetscViewerDrawSetInfo(petsclib::PetscLibType, v::PetscViewer, display::String, title::String, x::Cint, y::Cint, w::Cint, h::Cint) end

@for_petsc function PetscViewerDrawSetInfo(petsclib::$UnionPetscLib, v::PetscViewer, display::String, title::String, x::Cint, y::Cint, w::Cint, h::Cint )

    @chk ccall(
               (:PetscViewerDrawSetInfo, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cchar}, Ptr{Cchar}, Cint, Cint, Cint, Cint),
               v, display, title, x, y, w, h,
              )


	return nothing
end 

"""
	PetscViewerDrawSetTitle(petsclib::PetscLibType,v::PetscViewer, title::String) 

# External Links
$(_doc_external("Sys/PetscViewerDrawSetTitle"))
"""
function PetscViewerDrawSetTitle(petsclib::PetscLibType, v::PetscViewer, title::String) end

@for_petsc function PetscViewerDrawSetTitle(petsclib::$UnionPetscLib, v::PetscViewer, title::String )

    @chk ccall(
               (:PetscViewerDrawSetTitle, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cchar}),
               v, title,
              )


	return nothing
end 

"""
	PetscViewerDrawGetTitle(petsclib::PetscLibType,v::PetscViewer, title::String) 

# External Links
$(_doc_external("Sys/PetscViewerDrawGetTitle"))
"""
function PetscViewerDrawGetTitle(petsclib::PetscLibType, v::PetscViewer, title::String) end

@for_petsc function PetscViewerDrawGetTitle(petsclib::$UnionPetscLib, v::PetscViewer, title::String )
	title_ = Ref(pointer(title))

    @chk ccall(
               (:PetscViewerDrawGetTitle, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Ptr{Cchar}}),
               v, title_,
              )


	return nothing
end 

"""
	PetscViewerDrawOpen(petsclib::PetscLibType,comm::MPI_Comm, display::String, title::String, x::Cint, y::Cint, w::Cint, h::Cint, viewer::PetscViewer) 
Opens a `PetscDraw` window for use as a `PetscViewer` with type
`PETSCVIEWERDRAW`.

Collective

Input Parameters:
- `comm`    - communicator that will share window
- `display` - the X display on which to open, or `NULL` for the local machine
- `title`   - the title to put in the title bar, or `NULL` for no title
- `x`       - horizontal screen coordinate of the upper left corner of window, or use `PETSC_DECIDE`
- `y`       - vertical screen coordinate of the upper left corner of window, or use `PETSC_DECIDE`
- `w`       - window width in pixels, or may use `PETSC_DECIDE` or `PETSC_DRAW_FULL_SIZE`, `PETSC_DRAW_HALF_SIZE`,`PETSC_DRAW_THIRD_SIZE`, `PETSC_DRAW_QUARTER_SIZE`
- `h`       - window height in pixels, or may use `PETSC_DECIDE` or `PETSC_DRAW_FULL_SIZE`, `PETSC_DRAW_HALF_SIZE`,`PETSC_DRAW_THIRD_SIZE`, `PETSC_DRAW_QUARTER_SIZE`

Output Parameter:
- `viewer` - the `PetscViewer`

Options Database Keys:
- `-draw_type`          - use x or null
- `-nox`                - Disables all x-windows output
- `-display <name>`     - Specifies name of machine for the X display
- `-geometry <x,y,w,h>` - allows setting the window location and size
- `-draw_pause <pause>` - Sets time (in seconds) that the
program pauses after `PetscDrawPause()` has been called
(0 is default, -1 implies until user input).

Level: beginner

-seealso: [](sec_viewers), `PETSCVIEWERDRAW`, `PetscDrawCreate()`, `PetscViewerDestroy()`, `PetscViewerDrawGetDraw()`, `PetscViewerCreate()`, `PETSC_VIEWER_DRAW_`,
`PETSC_VIEWER_DRAW_WORLD`, `PETSC_VIEWER_DRAW_SELF`

# External Links
$(_doc_external("Sys/PetscViewerDrawOpen"))
"""
function PetscViewerDrawOpen(petsclib::PetscLibType, comm::MPI_Comm, display::String, title::String, x::Cint, y::Cint, w::Cint, h::Cint, viewer::PetscViewer) end

@for_petsc function PetscViewerDrawOpen(petsclib::$UnionPetscLib, comm::MPI_Comm, display::String, title::String, x::Cint, y::Cint, w::Cint, h::Cint, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerDrawOpen, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Cint, Cint, Cint, Cint, Ptr{PetscViewer}),
               comm, display, title, x, y, w, h, viewer,
              )


	return nothing
end 

"""
	PetscViewerDrawClear(petsclib::PetscLibType,viewer::PetscViewer) 
Clears a `PetscDraw` graphic associated with a `PetscViewer`.

Not Collective

Input Parameter:
- `viewer` - the `PetscViewer`

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERDRAW`, `PetscViewerDrawOpen()`, `PetscViewerDrawGetDraw()`,

# External Links
$(_doc_external("Sys/PetscViewerDrawClear"))
"""
function PetscViewerDrawClear(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerDrawClear(petsclib::$UnionPetscLib, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerDrawClear, $petsc_library),
               PetscErrorCode,
               (PetscViewer,),
               viewer,
              )


	return nothing
end 

"""
	pause::PetscReal = PetscViewerDrawGetPause(petsclib::PetscLibType,viewer::PetscViewer) 
Gets the pause value (how long to pause before an image is changed)  in the `PETSCVIEWERDRAW` `PetscViewer`

Not Collective

Input Parameter:
- `viewer` - the `PetscViewer`

Output Parameter:
- `pause` - the pause value

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERDRAW`, `PetscViewerDrawOpen()`, `PetscViewerDrawGetDraw()`,

# External Links
$(_doc_external("Sys/PetscViewerDrawGetPause"))
"""
function PetscViewerDrawGetPause(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerDrawGetPause(petsclib::$UnionPetscLib, viewer::PetscViewer )
	pause_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscViewerDrawGetPause, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{$PetscReal}),
               viewer, pause_,
              )

	pause = pause_[]

	return pause
end 

"""
	PetscViewerDrawSetPause(petsclib::PetscLibType,viewer::PetscViewer, pause::PetscReal) 
Sets a pause for each `PetscDraw` in the `PETSCVIEWERDRAW` `PetscViewer`

Not Collective

Input Parameters:
- `viewer` - the `PetscViewer`
- `pause`  - the pause value

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERDRAW`, `PetscViewerDrawOpen()`, `PetscViewerDrawGetDraw()`,

# External Links
$(_doc_external("Sys/PetscViewerDrawSetPause"))
"""
function PetscViewerDrawSetPause(petsclib::PetscLibType, viewer::PetscViewer, pause::PetscReal) end

@for_petsc function PetscViewerDrawSetPause(petsclib::$UnionPetscLib, viewer::PetscViewer, pause::$PetscReal )

    @chk ccall(
               (:PetscViewerDrawSetPause, $petsc_library),
               PetscErrorCode,
               (PetscViewer, $PetscReal),
               viewer, pause,
              )


	return nothing
end 

"""
	PetscViewerDrawSetHold(petsclib::PetscLibType,viewer::PetscViewer, hold::PetscBool) 
Holds previous image when drawing new image in a `PETSCVIEWERDRAW`

Not Collective

Input Parameters:
- `viewer` - the `PetscViewer`
- `hold`   - `PETSC_TRUE` indicates to hold the previous image

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERDRAW`, `PetscViewerDrawOpen()`, `PetscViewerDrawGetDraw()`,

# External Links
$(_doc_external("Sys/PetscViewerDrawSetHold"))
"""
function PetscViewerDrawSetHold(petsclib::PetscLibType, viewer::PetscViewer, hold::PetscBool) end

@for_petsc function PetscViewerDrawSetHold(petsclib::$UnionPetscLib, viewer::PetscViewer, hold::PetscBool )

    @chk ccall(
               (:PetscViewerDrawSetHold, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscBool),
               viewer, hold,
              )


	return nothing
end 

"""
	hold::PetscBool = PetscViewerDrawGetHold(petsclib::PetscLibType,viewer::PetscViewer) 
Checks if the `PETSCVIEWERDRAW` `PetscViewer` holds previous image when drawing new image

Not Collective

Input Parameter:
- `viewer` - the `PetscViewer`

Output Parameter:
- `hold` - indicates to hold or not

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERDRAW`, `PetscViewerDrawOpen()`, `PetscViewerDrawGetDraw()`,

# External Links
$(_doc_external("Sys/PetscViewerDrawGetHold"))
"""
function PetscViewerDrawGetHold(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerDrawGetHold(petsclib::$UnionPetscLib, viewer::PetscViewer )
	hold_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscViewerDrawGetHold, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{PetscBool}),
               viewer, hold_,
              )

	hold = hold_[]

	return hold
end 

"""
	PetscViewerDrawSetBounds(petsclib::PetscLibType,viewer::PetscViewer, nbounds::PetscInt, bounds::PetscReal) 
sets the upper and lower bounds to be used in plotting in a `PETSCVIEWERDRAW` `PetscViewer`

Collective

Input Parameters:
- `viewer`  - the `PetscViewer` (created with `PetscViewerDrawOpen()`)
- `nbounds` - number of plots that can be made with this viewer, for example the dof passed to `DMDACreate()`
- `bounds`  - the actual bounds, the size of this is 2*`nbounds`, the values are stored in the order min F_0, max F_0, min F_1, max F_1, .....

Options Database Key:
- `-draw_bounds  minF0,maxF0,minF1,maxF1` - the lower left and upper right bounds

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERDRAW`, `PetscViewerDrawGetLG()`, `PetscViewerDrawGetAxis()`, `PetscViewerDrawOpen()`

# External Links
$(_doc_external("Sys/PetscViewerDrawSetBounds"))
"""
function PetscViewerDrawSetBounds(petsclib::PetscLibType, viewer::PetscViewer, nbounds::PetscInt, bounds::PetscReal) end

@for_petsc function PetscViewerDrawSetBounds(petsclib::$UnionPetscLib, viewer::PetscViewer, nbounds::$PetscInt, bounds::$PetscReal )

    @chk ccall(
               (:PetscViewerDrawSetBounds, $petsc_library),
               PetscErrorCode,
               (PetscViewer, $PetscInt, Ptr{$PetscReal}),
               viewer, nbounds, bounds,
              )


	return nothing
end 

"""
	nbounds::PetscInt,bounds::Vector{PetscReal} = PetscViewerDrawGetBounds(petsclib::PetscLibType,viewer::PetscViewer) 
gets the upper and lower bounds to be used in plotting set with `PetscViewerDrawSetBounds()`

Collective

Input Parameter:
- `viewer` - the `PetscViewer` (created with `PetscViewerDrawOpen()`)

Output Parameters:
- `nbounds` - number of plots that can be made with this viewer, for example the dof passed to `DMDACreate()`
- `bounds`  - the actual bounds, the size of this is 2*`nbounds`, the values are stored in the order min F_0, max F_0, min F_1, max F_1, .....

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERDRAW`, `PetscViewerDrawGetLG()`, `PetscViewerDrawGetAxis()`, `PetscViewerDrawOpen()`, `PetscViewerDrawSetBounds()`

# External Links
$(_doc_external("Sys/PetscViewerDrawGetBounds"))
"""
function PetscViewerDrawGetBounds(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerDrawGetBounds(petsclib::$UnionPetscLib, viewer::PetscViewer )
	nbounds_ = Ref{$PetscInt}()
	bounds_ = Ref{Ptr{$PetscReal}}()

    @chk ccall(
               (:PetscViewerDrawGetBounds, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{$PetscInt}, Ptr{Ptr{$PetscReal}}),
               viewer, nbounds_, bounds_,
              )

	nbounds = nbounds_[]
	bounds = unsafe_wrap(Array, bounds_[], VecGetLocalSize(petsclib, x); own = false)

	return nbounds,bounds
end 

"""
	PetscViewerMonitorLGSetUp(petsclib::PetscLibType,viewer::PetscViewer, host::String, title::String, metric::String, l::PetscInt, names::String, x::Cint, y::Cint, m::Cint, n::Cint) 
sets up a viewer to be used by line graph monitoring routines such as `KSPMonitorResidualDrawLG()`

Collective

Input Parameters:
- `viewer` - the viewer in which to display the line graphs, it not a `PETSCVIEWERDRAW` it will set to that `PetscViewerType`
- `host`   - the host to open the window on, 'NULL' indicates the local host
- `title`  - the title at the top of the window
- `metric` - the label above the graph
- `l`      - the number of curves
- `names`  - the names of each curve to be used in displaying the legend. May be 'NULL'
- `x`      - horizontal screen coordinate of the upper left corner of window, or use `PETSC_DECIDE`
- `y`      - vertical screen coordinate of the upper left corner of window, or use `PETSC_DECIDE`
- `m`      - window width in pixels, or may use `PETSC_DECIDE` or `PETSC_DRAW_FULL_SIZE`, `PETSC_DRAW_HALF_SIZE`,`PETSC_DRAW_THIRD_SIZE`, `PETSC_DRAW_QUARTER_SIZE`
- `n`      - window height in pixels, or may use `PETSC_DECIDE` or `PETSC_DRAW_FULL_SIZE`, `PETSC_DRAW_HALF_SIZE`,`PETSC_DRAW_THIRD_SIZE`, `PETSC_DRAW_QUARTER_SIZE`

Level: developer

-seealso: `PetscViewer()`, `PETSCVIEWERDRAW`, `PetscViewerDrawGetDrawLG()`, `PetscViewerDrawOpen()`, `PetscViewerDrawSetInfo()`

# External Links
$(_doc_external("Sys/PetscViewerMonitorLGSetUp"))
"""
function PetscViewerMonitorLGSetUp(petsclib::PetscLibType, viewer::PetscViewer, host::String, title::String, metric::String, l::PetscInt, names::String, x::Cint, y::Cint, m::Cint, n::Cint) end

@for_petsc function PetscViewerMonitorLGSetUp(petsclib::$UnionPetscLib, viewer::PetscViewer, host::String, title::String, metric::String, l::$PetscInt, names::String, x::Cint, y::Cint, m::Cint, n::Cint )
	names_ = Ref(pointer(names))

    @chk ccall(
               (:PetscViewerMonitorLGSetUp, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, $PetscInt, Ptr{Ptr{Cchar}}, Cint, Cint, Cint, Cint),
               viewer, host, title, metric, l, names_, x, y, m, n,
              )


	return nothing
end 

"""
	PetscViewerDrawGetDraw(petsclib::PetscLibType,viewer::PetscViewer, windownumber::PetscInt, draw::PetscDraw) 
Returns `PetscDraw` object from `PETSCVIEWERDRAW` `PetscViewer` object.
This `PetscDraw` object may then be used to perform graphics using `PetscDraw` commands.

Collective

Input Parameters:
- `viewer`       - the `PetscViewer` (created with `PetscViewerDrawOpen()` of type `PETSCVIEWERDRAW`)
- `windownumber` - indicates which subwindow (usually 0) to obtain

Output Parameter:
- `draw` - the draw object

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERDRAW`, `PetscViewerDrawGetLG()`, `PetscViewerDrawGetAxis()`, `PetscViewerDrawOpen()`

# External Links
$(_doc_external("Sys/PetscViewerDrawGetDraw"))
"""
function PetscViewerDrawGetDraw(petsclib::PetscLibType, viewer::PetscViewer, windownumber::PetscInt, draw::PetscDraw) end

@for_petsc function PetscViewerDrawGetDraw(petsclib::$UnionPetscLib, viewer::PetscViewer, windownumber::$PetscInt, draw::PetscDraw )

    @chk ccall(
               (:PetscViewerDrawGetDraw, $petsc_library),
               PetscErrorCode,
               (PetscViewer, $PetscInt, Ptr{PetscDraw}),
               viewer, windownumber, draw,
              )


	return nothing
end 

"""
	PetscViewerDrawGetDrawLG(petsclib::PetscLibType,viewer::PetscViewer, windownumber::PetscInt, drawlg::PetscDrawLG) 
Returns a `PetscDrawLG` object from `PetscViewer` object of type `PETSCVIEWERDRAW`.
This `PetscDrawLG` object may then be used to perform graphics using `PetscDrawLG` commands.

Collective

Input Parameters:
- `viewer`       - the `PetscViewer` (created with `PetscViewerDrawOpen()`)
- `windownumber` - indicates which subwindow (usually 0)

Output Parameter:
- `drawlg` - the draw line graph object

Level: intermediate

-seealso: [](sec_viewers), `PetscDrawLG`, `PetscViewerDrawGetDraw()`, `PetscViewerDrawGetAxis()`, `PetscViewerDrawOpen()`

# External Links
$(_doc_external("Sys/PetscViewerDrawGetDrawLG"))
"""
function PetscViewerDrawGetDrawLG(petsclib::PetscLibType, viewer::PetscViewer, windownumber::PetscInt, drawlg::PetscDrawLG) end

@for_petsc function PetscViewerDrawGetDrawLG(petsclib::$UnionPetscLib, viewer::PetscViewer, windownumber::$PetscInt, drawlg::PetscDrawLG )

    @chk ccall(
               (:PetscViewerDrawGetDrawLG, $petsc_library),
               PetscErrorCode,
               (PetscViewer, $PetscInt, Ptr{PetscDrawLG}),
               viewer, windownumber, drawlg,
              )


	return nothing
end 

"""
	PetscViewerDrawGetDrawAxis(petsclib::PetscLibType,viewer::PetscViewer, windownumber::PetscInt, drawaxis::PetscDrawAxis) 
Returns a `PetscDrawAxis` object from a `PetscViewer` object of type `PETSCVIEWERDRAW`.
This `PetscDrawAxis` object may then be used to perform graphics using `PetscDrawAxis` commands.

Collective

Input Parameters:
- `viewer`       - the `PetscViewer` (created with `PetscViewerDrawOpen()`)
- `windownumber` - indicates which subwindow (usually 0)

Output Parameter:
- `drawaxis` - the draw axis object

Level: advanced

-seealso: [](sec_viewers), `PetscViewerDrawGetDraw()`, `PetscViewerDrawGetLG()`, `PetscViewerDrawOpen()`

# External Links
$(_doc_external("Sys/PetscViewerDrawGetDrawAxis"))
"""
function PetscViewerDrawGetDrawAxis(petsclib::PetscLibType, viewer::PetscViewer, windownumber::PetscInt, drawaxis::PetscDrawAxis) end

@for_petsc function PetscViewerDrawGetDrawAxis(petsclib::$UnionPetscLib, viewer::PetscViewer, windownumber::$PetscInt, drawaxis::PetscDrawAxis )

    @chk ccall(
               (:PetscViewerDrawGetDrawAxis, $petsc_library),
               PetscErrorCode,
               (PetscViewer, $PetscInt, Ptr{PetscDrawAxis}),
               viewer, windownumber, drawaxis,
              )


	return nothing
end 

"""
	PetscViewerDrawSetDrawType(petsclib::PetscLibType,v::PetscViewer, drawtype::PetscDrawType) 

# External Links
$(_doc_external("Sys/PetscViewerDrawSetDrawType"))
"""
function PetscViewerDrawSetDrawType(petsclib::PetscLibType, v::PetscViewer, drawtype::PetscDrawType) end

@for_petsc function PetscViewerDrawSetDrawType(petsclib::$UnionPetscLib, v::PetscViewer, drawtype::PetscDrawType )

    @chk ccall(
               (:PetscViewerDrawSetDrawType, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscDrawType),
               v, drawtype,
              )


	return nothing
end 

"""
	drawtype::PetscDrawType = PetscViewerDrawGetDrawType(petsclib::PetscLibType,v::PetscViewer) 

# External Links
$(_doc_external("Sys/PetscViewerDrawGetDrawType"))
"""
function PetscViewerDrawGetDrawType(petsclib::PetscLibType, v::PetscViewer) end

@for_petsc function PetscViewerDrawGetDrawType(petsclib::$UnionPetscLib, v::PetscViewer )
	drawtype_ = Ref{PetscDrawType}()

    @chk ccall(
               (:PetscViewerDrawGetDrawType, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{PetscDrawType}),
               v, drawtype_,
              )

	drawtype = unsafe_string(drawtype_[])

	return drawtype
end 

"""
	PetscViewerADIOSOpen(petsclib::PetscLibType,comm::MPI_Comm, name::String, type::PetscFileMode, adiosv::PetscViewer) 
Opens a file for ADIOS input/output.

Collective

Input Parameters:
- `comm` - MPI communicator
- `name` - name of file
- `type` - type of file
-seealso: `PetscViewerASCIIOpen()`, `PetscViewerPushFormat()`, `PetscViewerDestroy()`, `PetscViewerHDF5Open()`,
`VecView()`, `MatView()`, `VecLoad()`, `PetscViewerSetType()`, `PetscViewerFileSetMode()`, `PetscViewerFileSetName()`
`MatLoad()`, `PetscFileMode`, `PetscViewer`

# External Links
$(_doc_external("Sys/PetscViewerADIOSOpen"))
"""
function PetscViewerADIOSOpen(petsclib::PetscLibType, comm::MPI_Comm, name::String, type::PetscFileMode, adiosv::PetscViewer) end

@for_petsc function PetscViewerADIOSOpen(petsclib::$UnionPetscLib, comm::MPI_Comm, name::String, type::PetscFileMode, adiosv::PetscViewer )

    @chk ccall(
               (:PetscViewerADIOSOpen, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, PetscFileMode, Ptr{PetscViewer}),
               comm, name, type, adiosv,
              )


	return nothing
end 

"""
	PetscViewerHDF5Load(petsclib::PetscLibType,viewer::PetscViewer, name::String, map::PetscLayout, datatype::hid_t, newarr::Cvoid) 
Read a raw array from the `PETSCVIEWERHDF5` dataset in parallel

Collective; No Fortran Support

Input Parameters:
- `viewer`   - The `PETSCVIEWERHDF5` viewer
- `name`     - The dataset name
- `datatype` - The HDF5 datatype of the items in the dataset

Input/Output Parameter:
- `map` - The layout which specifies array partitioning, on output the
set up layout (with global size and blocksize according to dataset)

Output Parameter:
- `newarr` - The partitioned array, a memory image of the given dataset

Level: developer

-seealso: `PetscViewer`, `PETSCVIEWERHDF5`, `PetscViewerHDF5Open()`, `PetscViewerHDF5PushGroup()`, `PetscViewerHDF5OpenGroup()`, `PetscViewerHDF5ReadSizes()`,
`VecLoad()`, `ISLoad()`, `PetscLayout`

# External Links
$(_doc_external("Sys/PetscViewerHDF5Load"))
"""
function PetscViewerHDF5Load(petsclib::PetscLibType, viewer::PetscViewer, name::String, map::PetscLayout, datatype::hid_t, newarr::Cvoid) end

@for_petsc function PetscViewerHDF5Load(petsclib::$UnionPetscLib, viewer::PetscViewer, name::String, map::PetscLayout, datatype::hid_t, newarr::Cvoid )

    @chk ccall(
               (:PetscViewerHDF5Load, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cchar}, PetscLayout, hid_t, Cvoid),
               viewer, name, map, datatype, newarr,
              )


	return nothing
end 

"""
	bs::PetscInt,N::PetscInt = PetscViewerHDF5ReadSizes(petsclib::PetscLibType,viewer::PetscViewer, name::String) 
Read block size and global size of a `Vec` or `IS` stored in an HDF5 file.

Input Parameters:
- `viewer` - The `PETSCVIEWERHDF5` viewer
- `name`   - The dataset name

Output Parameters:
- `bs` - block size
- `N`  - global size

Level: advanced

-seealso: `PetscViewer`, `PETSCVIEWERHDF5`, `PetscViewerHDF5Open()`, `VecLoad()`, `ISLoad()`, `VecGetSize()`, `ISGetSize()`, `PetscViewerHDF5SetBaseDimension2()`

# External Links
$(_doc_external("Sys/PetscViewerHDF5ReadSizes"))
"""
function PetscViewerHDF5ReadSizes(petsclib::PetscLibType, viewer::PetscViewer, name::String) end

@for_petsc function PetscViewerHDF5ReadSizes(petsclib::$UnionPetscLib, viewer::PetscViewer, name::String )
	bs_ = Ref{$PetscInt}()
	N_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscViewerHDF5ReadSizes, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cchar}, Ptr{$PetscInt}, Ptr{$PetscInt}),
               viewer, name, bs_, N_,
              )

	bs = bs_[]
	N = N_[]

	return bs,N
end 

"""
	PetscViewerExodusIIGetId(petsclib::PetscLibType,viewer::PetscViewer, exoid::Cint) 
Get the file id of the `PETSCVIEWEREXODUSII` file

Logically Collective

Input Parameter:
- `viewer` - the `PetscViewer`

Output Parameter:
- `exoid` - The ExodusII file id

Level: intermediate

-seealso: `PETSCVIEWEREXODUSII`, `PetscViewer`, `PetscViewerFileSetMode()`, `PetscViewerCreate()`, `PetscViewerSetType()`, `PetscViewerBinaryOpen()`

# External Links
$(_doc_external("Dm/PetscViewerExodusIIGetId"))
"""
function PetscViewerExodusIIGetId(petsclib::PetscLibType, viewer::PetscViewer, exoid::Cint) end

@for_petsc function PetscViewerExodusIIGetId(petsclib::$UnionPetscLib, viewer::PetscViewer, exoid::Cint )

    @chk ccall(
               (:PetscViewerExodusIIGetId, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cint}),
               viewer, exoid,
              )


	return nothing
end 

"""
	PetscViewerExodusIISetOrder(petsclib::PetscLibType,viewer::PetscViewer, order::PetscInt) 
Set the elements order in the ExodusII file.

Collective

Input Parameters:
- `viewer` - the `PETSCVIEWEREXODUSII` viewer
- `order`  - elements order

Output Parameter:

Level: beginner

-seealso: `PETSCVIEWEREXODUSII`, `PetscViewer`, `PetscViewerExodusIIGetId()`, `PetscViewerExodusIIGetOrder()`

# External Links
$(_doc_external("Dm/PetscViewerExodusIISetOrder"))
"""
function PetscViewerExodusIISetOrder(petsclib::PetscLibType, viewer::PetscViewer, order::PetscInt) end

@for_petsc function PetscViewerExodusIISetOrder(petsclib::$UnionPetscLib, viewer::PetscViewer, order::$PetscInt )

    @chk ccall(
               (:PetscViewerExodusIISetOrder, $petsc_library),
               PetscErrorCode,
               (PetscViewer, $PetscInt),
               viewer, order,
              )


	return nothing
end 

"""
	PetscViewerExodusIIGetOrder(petsclib::PetscLibType,viewer::PetscViewer, order::PetscInt) 
Get the elements order in the ExodusII file.

Collective

Input Parameters:
- `viewer` - the `PETSCVIEWEREXODUSII` viewer
- `order`  - elements order

Output Parameter:

Level: beginner

-seealso: `PETSCVIEWEREXODUSII`, `PetscViewer`, `PetscViewerExodusIIGetId()`, `PetscViewerExodusIISetOrder()`

# External Links
$(_doc_external("Dm/PetscViewerExodusIIGetOrder"))
"""
function PetscViewerExodusIIGetOrder(petsclib::PetscLibType, viewer::PetscViewer, order::PetscInt) end

@for_petsc function PetscViewerExodusIIGetOrder(petsclib::$UnionPetscLib, viewer::PetscViewer, order::$PetscInt )

    @chk ccall(
               (:PetscViewerExodusIIGetOrder, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{$PetscInt}),
               viewer, order,
              )


	return nothing
end 

"""
	PetscViewerExodusIIOpen(petsclib::PetscLibType,comm::MPI_Comm, name::String, type::PetscFileMode, exo::PetscViewer) 
Opens a file for ExodusII input/output.

Collective

Input Parameters:
- `comm` - MPI communicator
- `name` - name of file
- `type` - type of file
-seealso: `PETSCVIEWEREXODUSII`, `PetscViewer`, `PetscViewerPushFormat()`, `PetscViewerDestroy()`,
`DMLoad()`, `PetscFileMode`, `PetscViewerSetType()`, `PetscViewerFileSetMode()`, `PetscViewerFileSetName()`

# External Links
$(_doc_external("Dm/PetscViewerExodusIIOpen"))
"""
function PetscViewerExodusIIOpen(petsclib::PetscLibType, comm::MPI_Comm, name::String, type::PetscFileMode, exo::PetscViewer) end

@for_petsc function PetscViewerExodusIIOpen(petsclib::$UnionPetscLib, comm::MPI_Comm, name::String, type::PetscFileMode, exo::PetscViewer )

    @chk ccall(
               (:PetscViewerExodusIIOpen, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, PetscFileMode, Ptr{PetscViewer}),
               comm, name, type, exo,
              )


	return nothing
end 

"""
	PetscViewerHDF5SetDMPlexStorageVersionWriting(petsclib::PetscLibType,viewer::PetscViewer, version::DMPlexStorageVersion) 
Set the storage version for writing

Logically collective

Input Parameters:
- `viewer`  - The `PetscViewer`
- `version` - The storage format version

Level: advanced

Note:
The version has major, minor, and subminor integers. Parallel operations are only available for version 3.0.0.

See also: 
=== 
`DM`, `PetscViewerHDF5GetDMPlexStorageVersionWriting()`, `PetscViewerHDF5GetDMPlexStorageVersionReading()`, `PetscViewerHDF5SetDMPlexStorageVersionReading()`

# External Links
$(_doc_external("Dm/PetscViewerHDF5SetDMPlexStorageVersionWriting"))
"""
function PetscViewerHDF5SetDMPlexStorageVersionWriting(petsclib::PetscLibType, viewer::PetscViewer, version::DMPlexStorageVersion) end

@for_petsc function PetscViewerHDF5SetDMPlexStorageVersionWriting(petsclib::$UnionPetscLib, viewer::PetscViewer, version::DMPlexStorageVersion )

    @chk ccall(
               (:PetscViewerHDF5SetDMPlexStorageVersionWriting, $petsc_library),
               PetscErrorCode,
               (PetscViewer, DMPlexStorageVersion),
               viewer, version,
              )


	return nothing
end 

"""
	PetscViewerHDF5GetDMPlexStorageVersionWriting(petsclib::PetscLibType,viewer::PetscViewer, version::DMPlexStorageVersion) 
Get the storage version for writing

Logically collective

Input Parameter:
- `viewer` - The `PetscViewer`

Output Parameter:
- `version` - The storage format version

Options Database Keys:
- `-dm_plex_view_hdf5_storage_version <num>` - Overrides the storage format version

Level: advanced

Note:
The version has major, minor, and subminor integers. Parallel operations are only available for version 3.0.0.

See also: 
=== 
`DM`, `PetscViewerHDF5SetDMPlexStorageVersionWriting()`, `PetscViewerHDF5GetDMPlexStorageVersionReading()`, `PetscViewerHDF5SetDMPlexStorageVersionReading()`

# External Links
$(_doc_external("Dm/PetscViewerHDF5GetDMPlexStorageVersionWriting"))
"""
function PetscViewerHDF5GetDMPlexStorageVersionWriting(petsclib::PetscLibType, viewer::PetscViewer, version::DMPlexStorageVersion) end

@for_petsc function PetscViewerHDF5GetDMPlexStorageVersionWriting(petsclib::$UnionPetscLib, viewer::PetscViewer, version::DMPlexStorageVersion )

    @chk ccall(
               (:PetscViewerHDF5GetDMPlexStorageVersionWriting, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{DMPlexStorageVersion}),
               viewer, version,
              )


	return nothing
end 

"""
	PetscViewerHDF5SetDMPlexStorageVersionReading(petsclib::PetscLibType,viewer::PetscViewer, version::DMPlexStorageVersion) 
Set the storage version for reading

Logically collective

Input Parameters:
- `viewer`  - The `PetscViewer`
- `version` - The storage format version

Level: advanced

Note:
The version has major, minor, and subminor integers. Parallel operations are only available for version 3.0.0.

See also: 
=== 
`DM`, `PetscViewerHDF5GetDMPlexStorageVersionReading()`, `PetscViewerHDF5GetDMPlexStorageVersionWriting()`, `PetscViewerHDF5SetDMPlexStorageVersionWriting()`

# External Links
$(_doc_external("Dm/PetscViewerHDF5SetDMPlexStorageVersionReading"))
"""
function PetscViewerHDF5SetDMPlexStorageVersionReading(petsclib::PetscLibType, viewer::PetscViewer, version::DMPlexStorageVersion) end

@for_petsc function PetscViewerHDF5SetDMPlexStorageVersionReading(petsclib::$UnionPetscLib, viewer::PetscViewer, version::DMPlexStorageVersion )

    @chk ccall(
               (:PetscViewerHDF5SetDMPlexStorageVersionReading, $petsc_library),
               PetscErrorCode,
               (PetscViewer, DMPlexStorageVersion),
               viewer, version,
              )


	return nothing
end 

"""
	PetscViewerHDF5GetDMPlexStorageVersionReading(petsclib::PetscLibType,viewer::PetscViewer, version::DMPlexStorageVersion) 
Get the storage version for reading

Logically collective

Input Parameter:
- `viewer` - The `PetscViewer`

Output Parameter:
- `version` - The storage format version

Options Database Keys:
- `-dm_plex_view_hdf5_storage_version <num>` - Overrides the storage format version

Level: advanced

Note:
The version has major, minor, and subminor integers. Parallel operations are only available for version 3.0.0.

See also: 
=== 
`DM`, `PetscViewerHDF5SetDMPlexStorageVersionReading()`, `PetscViewerHDF5GetDMPlexStorageVersionWriting()`, `PetscViewerHDF5SetDMPlexStorageVersionWriting()`

# External Links
$(_doc_external("Dm/PetscViewerHDF5GetDMPlexStorageVersionReading"))
"""
function PetscViewerHDF5GetDMPlexStorageVersionReading(petsclib::PetscLibType, viewer::PetscViewer, version::DMPlexStorageVersion) end

@for_petsc function PetscViewerHDF5GetDMPlexStorageVersionReading(petsclib::$UnionPetscLib, viewer::PetscViewer, version::DMPlexStorageVersion )

    @chk ccall(
               (:PetscViewerHDF5GetDMPlexStorageVersionReading, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{DMPlexStorageVersion}),
               viewer, version,
              )


	return nothing
end 

"""
	PetscViewerExodusIISetZonalVariable(petsclib::PetscLibType,viewer::PetscViewer, num::PetscExodusIIInt) 
Sets the number of zonal variables in an ExodusII file

Collective;

Input Parameters:
- `viewer` - a `PetscViewer` of type `PETSCVIEWEREXODUSII`
- `num`    - the number of zonal variables in the ExodusII file

Level: intermediate

-seealso: `PETSCVIEWEREXODUSII`, `PetscViewer`, `PetscViewerCreate()`, `PetscViewerDestroy()`, `PetscViewerExodusIIOpen()`, `PetscViewerSetType()`, `PetscViewerType`, `PetscViewerExodusIIGetZonalVariable()`

# External Links
$(_doc_external("Dm/PetscViewerExodusIISetZonalVariable"))
"""
function PetscViewerExodusIISetZonalVariable(petsclib::PetscLibType, viewer::PetscViewer, num::PetscExodusIIInt) end

@for_petsc function PetscViewerExodusIISetZonalVariable(petsclib::$UnionPetscLib, viewer::PetscViewer, num::PetscExodusIIInt )

    @chk ccall(
               (:PetscViewerExodusIISetZonalVariable, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscExodusIIInt),
               viewer, num,
              )


	return nothing
end 

"""
	PetscViewerExodusIISetNodalVariable(petsclib::PetscLibType,viewer::PetscViewer, num::PetscExodusIIInt) 
Sets the number of nodal variables in an ExodusII file

Collective;

Input Parameters:
- `viewer` - a `PetscViewer` of type `PETSCVIEWEREXODUSII`
- `num`    - the number of nodal variables in the ExodusII file

Level: intermediate

-seealso: `PETSCVIEWEREXODUSII`, `PetscViewer`, `PetscViewerCreate()`, `PetscViewerDestroy()`, `PetscViewerExodusIIOpen()`, `PetscViewerSetType()`, `PetscViewerType`, `PetscViewerExodusIIGetNodalVariable()`

# External Links
$(_doc_external("Dm/PetscViewerExodusIISetNodalVariable"))
"""
function PetscViewerExodusIISetNodalVariable(petsclib::PetscLibType, viewer::PetscViewer, num::PetscExodusIIInt) end

@for_petsc function PetscViewerExodusIISetNodalVariable(petsclib::$UnionPetscLib, viewer::PetscViewer, num::PetscExodusIIInt )

    @chk ccall(
               (:PetscViewerExodusIISetNodalVariable, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscExodusIIInt),
               viewer, num,
              )


	return nothing
end 

"""
	PetscViewerExodusIIGetZonalVariable(petsclib::PetscLibType,viewer::PetscViewer, num::PetscExodusIIInt) 
Gets the number of zonal variables in an ExodusII file

Collective

Input Parameters:
- `viewer` - a `PetscViewer` of type `PETSCVIEWEREXODUSII`

Output Parameter:
- `num` - the number variables in the ExodusII file

Level: intermediate

-seealso: `PETSCVIEWEREXODUSII`, `PetscViewer`, `PetscViewerCreate()`, `PetscViewerDestroy()`, `PetscViewerExodusIIOpen()`, `PetscViewerSetType()`, `PetscViewerType`, `PetscViewerExodusIIsetZonalVariable()`

# External Links
$(_doc_external("Dm/PetscViewerExodusIIGetZonalVariable"))
"""
function PetscViewerExodusIIGetZonalVariable(petsclib::PetscLibType, viewer::PetscViewer, num::PetscExodusIIInt) end

@for_petsc function PetscViewerExodusIIGetZonalVariable(petsclib::$UnionPetscLib, viewer::PetscViewer, num::PetscExodusIIInt )

    @chk ccall(
               (:PetscViewerExodusIIGetZonalVariable, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{PetscExodusIIInt}),
               viewer, num,
              )


	return nothing
end 

"""
	PetscViewerExodusIIGetNodalVariable(petsclib::PetscLibType,viewer::PetscViewer, num::PetscExodusIIInt) 
Gets the number of nodal variables in an ExodusII file

Collective

Input Parameters:
- `viewer` - a `PetscViewer` of type `PETSCVIEWEREXODUSII`

Output Parameter:
- `num` - the number variables in the ExodusII file

Level: intermediate

-seealso: `PETSCVIEWEREXODUSII`, `PetscViewer`, `PetscViewerCreate()`, `PetscViewerDestroy()`, `PetscViewerExodusIIOpen()`, `PetscViewerSetType()`, `PetscViewerType`, `PetscViewerExodusIISetNodalVariable()`

# External Links
$(_doc_external("Dm/PetscViewerExodusIIGetNodalVariable"))
"""
function PetscViewerExodusIIGetNodalVariable(petsclib::PetscLibType, viewer::PetscViewer, num::PetscExodusIIInt) end

@for_petsc function PetscViewerExodusIIGetNodalVariable(petsclib::$UnionPetscLib, viewer::PetscViewer, num::PetscExodusIIInt )

    @chk ccall(
               (:PetscViewerExodusIIGetNodalVariable, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{PetscExodusIIInt}),
               viewer, num,
              )


	return nothing
end 

"""
	PetscViewerExodusIISetZonalVariableName(petsclib::PetscLibType,viewer::PetscViewer, idx::PetscExodusIIInt, name::String) 
Sets the name of a zonal variable.

Collective;

Input Parameters:
- `viewer` - a `PetscViewer` of type `PETSCVIEWEREXODUSII`
- `idx`    - the index for which you want to save the name
- `name`   - string containing the name characters

Level: intermediate

-seealso: `PETSCVIEWEREXODUSII`, `PetscViewer`, `PetscViewerCreate()`, `PetscViewerDestroy()`, `PetscViewerExodusIIOpen()`, `PetscViewerSetType()`, `PetscViewerType`, `PetscViewerExodusIIGetZonalVariableName()`

# External Links
$(_doc_external("Dm/PetscViewerExodusIISetZonalVariableName"))
"""
function PetscViewerExodusIISetZonalVariableName(petsclib::PetscLibType, viewer::PetscViewer, idx::PetscExodusIIInt, name::String) end

@for_petsc function PetscViewerExodusIISetZonalVariableName(petsclib::$UnionPetscLib, viewer::PetscViewer, idx::PetscExodusIIInt, name::String )

    @chk ccall(
               (:PetscViewerExodusIISetZonalVariableName, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscExodusIIInt, Ptr{Cchar}),
               viewer, idx, name,
              )


	return nothing
end 

"""
	PetscViewerExodusIISetNodalVariableName(petsclib::PetscLibType,viewer::PetscViewer, idx::PetscExodusIIInt, name::String) 
Sets the name of a nodal variable.

Collective;

Input Parameters:
- `viewer` - a `PetscViewer` of type `PETSCVIEWEREXODUSII`
- `idx`    - the index for which you want to save the name
- `name`   - string containing the name characters

Level: intermediate

-seealso: `PETSCVIEWEREXODUSII`, `PetscViewer`, `PetscViewerCreate()`, `PetscViewerDestroy()`, `PetscViewerExodusIIOpen()`, `PetscViewerSetType()`, `PetscViewerType`, `PetscViewerExodusIIGetNodalVariableName()`

# External Links
$(_doc_external("Dm/PetscViewerExodusIISetNodalVariableName"))
"""
function PetscViewerExodusIISetNodalVariableName(petsclib::PetscLibType, viewer::PetscViewer, idx::PetscExodusIIInt, name::String) end

@for_petsc function PetscViewerExodusIISetNodalVariableName(petsclib::$UnionPetscLib, viewer::PetscViewer, idx::PetscExodusIIInt, name::String )

    @chk ccall(
               (:PetscViewerExodusIISetNodalVariableName, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscExodusIIInt, Ptr{Cchar}),
               viewer, idx, name,
              )


	return nothing
end 

"""
	PetscViewerExodusIIGetZonalVariableName(petsclib::PetscLibType,viewer::PetscViewer, idx::PetscExodusIIInt, name::String) 
Gets the name of a zonal variable.

Collective;

Input Parameters:
- `viewer` - a `PetscViewer` of type `PETSCVIEWEREXODUSII`
- `idx`    - the index for which you want to get the name

Output Parameter:
- `name` - pointer to the string containing the name characters

Level: intermediate

-seealso: `PETSCVIEWEREXODUSII`, `PetscViewer`, `PetscViewerCreate()`, `PetscViewerDestroy()`, `PetscViewerExodusIIOpen()`, `PetscViewerSetType()`, `PetscViewerType`, `PetscViewerExodusIISetZonalVariableName()`

# External Links
$(_doc_external("Dm/PetscViewerExodusIIGetZonalVariableName"))
"""
function PetscViewerExodusIIGetZonalVariableName(petsclib::PetscLibType, viewer::PetscViewer, idx::PetscExodusIIInt, name::String) end

@for_petsc function PetscViewerExodusIIGetZonalVariableName(petsclib::$UnionPetscLib, viewer::PetscViewer, idx::PetscExodusIIInt, name::String )
	name_ = Ref(pointer(name))

    @chk ccall(
               (:PetscViewerExodusIIGetZonalVariableName, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscExodusIIInt, Ptr{Ptr{Cchar}}),
               viewer, idx, name_,
              )


	return nothing
end 

"""
	PetscViewerExodusIIGetNodalVariableName(petsclib::PetscLibType,viewer::PetscViewer, idx::PetscExodusIIInt, name::String) 
Gets the name of a nodal variable.

Collective;

Input Parameters:
- `viewer` - a `PetscViewer` of type `PETSCVIEWEREXODUSII`
- `idx`    - the index for which you want to save the name

Output Parameter:
- `name` - string array containing name characters

Level: intermediate

-seealso: `PETSCVIEWEREXODUSII`, `PetscViewer`, `PetscViewerCreate()`, `PetscViewerDestroy()`, `PetscViewerExodusIIOpen()`, `PetscViewerSetType()`, `PetscViewerType`, `PetscViewerExodusIISetNodalVariableName()`

# External Links
$(_doc_external("Dm/PetscViewerExodusIIGetNodalVariableName"))
"""
function PetscViewerExodusIIGetNodalVariableName(petsclib::PetscLibType, viewer::PetscViewer, idx::PetscExodusIIInt, name::String) end

@for_petsc function PetscViewerExodusIIGetNodalVariableName(petsclib::$UnionPetscLib, viewer::PetscViewer, idx::PetscExodusIIInt, name::String )
	name_ = Ref(pointer(name))

    @chk ccall(
               (:PetscViewerExodusIIGetNodalVariableName, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscExodusIIInt, Ptr{Ptr{Cchar}}),
               viewer, idx, name_,
              )


	return nothing
end 

"""
	PetscViewerExodusIISetZonalVariableNames(petsclib::PetscLibType,viewer::PetscViewer, names::String) 
Sets the names of all nodal variables

Collective; No Fortran Support

Input Parameters:
- `viewer` - a `PetscViewer` of type `PETSCVIEWEREXODUSII`
- `names`  - an array of string names to be set, the strings are copied into the `PetscViewer`

Level: intermediate

-seealso: `PETSCVIEWEREXODUSII`, `PetscViewer`, `PetscViewerCreate()`, `PetscViewerDestroy()`, `PetscViewerExodusIIOpen()`, `PetscViewerSetType()`, `PetscViewerType`, `PetscViewerExodusIIGetZonalVariableNames()`

# External Links
$(_doc_external("Dm/PetscViewerExodusIISetZonalVariableNames"))
"""
function PetscViewerExodusIISetZonalVariableNames(petsclib::PetscLibType, viewer::PetscViewer, names::String) end

@for_petsc function PetscViewerExodusIISetZonalVariableNames(petsclib::$UnionPetscLib, viewer::PetscViewer, names::String )
	names_ = Ref(pointer(names))

    @chk ccall(
               (:PetscViewerExodusIISetZonalVariableNames, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Ptr{Cchar}}),
               viewer, names_,
              )


	return nothing
end 

"""
	PetscViewerExodusIISetNodalVariableNames(petsclib::PetscLibType,viewer::PetscViewer, names::String) 
Sets the names of all nodal variables.

Collective; No Fortran Support

Input Parameters:
- `viewer` - a `PetscViewer` of type `PETSCVIEWEREXODUSII`
- `names`  - an array of string names to be set, the strings are copied into the `PetscViewer`

Level: intermediate

-seealso: `PETSCVIEWEREXODUSII`, `PetscViewer`, `PetscViewerCreate()`, `PetscViewerDestroy()`, `PetscViewerExodusIIOpen()`, `PetscViewerSetType()`, `PetscViewerType`, `PetscViewerExodusIIGetNodalVariableNames()`

# External Links
$(_doc_external("Dm/PetscViewerExodusIISetNodalVariableNames"))
"""
function PetscViewerExodusIISetNodalVariableNames(petsclib::PetscLibType, viewer::PetscViewer, names::String) end

@for_petsc function PetscViewerExodusIISetNodalVariableNames(petsclib::$UnionPetscLib, viewer::PetscViewer, names::String )
	names_ = Ref(pointer(names))

    @chk ccall(
               (:PetscViewerExodusIISetNodalVariableNames, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Ptr{Cchar}}),
               viewer, names_,
              )


	return nothing
end 

"""
	PetscViewerExodusIIGetZonalVariableNames(petsclib::PetscLibType,viewer::PetscViewer, numVars::PetscExodusIIInt, varNames::String) 
Gets the names of all zonal variables.

Collective; No Fortran Support

Input Parameters:
- `viewer`  - a `PetscViewer` of type `PETSCVIEWEREXODUSII`
- `numVars` - the number of zonal variable names to retrieve

Output Parameter:
- `varNames` - returns an array of char pointers where the zonal variable names are

Level: intermediate

-seealso: `PETSCVIEWEREXODUSII`, `PetscViewer`, `PetscViewerCreate()`, `PetscViewerDestroy()`, `PetscViewerExodusIIOpen()`, `PetscViewerSetType()`, `PetscViewerType`, `PetscViewerExodusIISetZonalVariableNames()`

# External Links
$(_doc_external("Dm/PetscViewerExodusIIGetZonalVariableNames"))
"""
function PetscViewerExodusIIGetZonalVariableNames(petsclib::PetscLibType, viewer::PetscViewer, numVars::PetscExodusIIInt, varNames::String) end

@for_petsc function PetscViewerExodusIIGetZonalVariableNames(petsclib::$UnionPetscLib, viewer::PetscViewer, numVars::PetscExodusIIInt, varNames::String )

    @chk ccall(
               (:PetscViewerExodusIIGetZonalVariableNames, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{PetscExodusIIInt}, Ptr{Cchar}),
               viewer, numVars, varNames,
              )


	return nothing
end 

"""
	PetscViewerExodusIIGetNodalVariableNames(petsclib::PetscLibType,viewer::PetscViewer, numVars::PetscExodusIIInt, varNames::String) 
Gets the names of all nodal variables.

Collective; No Fortran Support

Input Parameters:
- `viewer`  - a `PetscViewer` of type `PETSCVIEWEREXODUSII`
- `numVars` - the number of nodal variable names to retrieve

Output Parameter:
- `varNames` - returns an array of char pointers where the nodal variable names are

Level: intermediate

-seealso: `PETSCVIEWEREXODUSII`, `PetscViewer`, `PetscViewerCreate()`, `PetscViewerDestroy()`, `PetscViewerExodusIIOpen()`, `PetscViewerSetType()`, `PetscViewerType`, `PetscViewerExodusIISetNodalVariableNames()`

# External Links
$(_doc_external("Dm/PetscViewerExodusIIGetNodalVariableNames"))
"""
function PetscViewerExodusIIGetNodalVariableNames(petsclib::PetscLibType, viewer::PetscViewer, numVars::PetscExodusIIInt, varNames::String) end

@for_petsc function PetscViewerExodusIIGetNodalVariableNames(petsclib::$UnionPetscLib, viewer::PetscViewer, numVars::PetscExodusIIInt, varNames::String )

    @chk ccall(
               (:PetscViewerExodusIIGetNodalVariableNames, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{PetscExodusIIInt}, Ptr{Cchar}),
               viewer, numVars, varNames,
              )


	return nothing
end 

"""
	PetscViewerExodusIIGetNodalVariableIndex(petsclib::PetscLibType,viewer::PetscViewer, name::String, varIndex::PetscExodusIIInt) 
return the location of a nodal variable in an ExodusII file given its name

Collective

Input Parameters:
- `viewer` - a `PetscViewer` of type `PETSCVIEWEREXODUSII`
- `name`   - the name of the result

Output Parameter:
- `varIndex` - the location of the variable in the exodus file or -1 if the variable is not found

Level: beginner

-seealso: `PetscViewerExodusIISetNodalVariable()`, `PetscViewerExodusIIGetNodalVariable()`, `PetscViewerExodusIISetNodalVariableName()`, `PetscViewerExodusIISetNodalVariableNames()`, `PetscViewerExodusIIGetNodalVariableName()`, `PetscViewerExodusIIGetNodalVariableNames()`

# External Links
$(_doc_external("Dm/PetscViewerExodusIIGetNodalVariableIndex"))
"""
function PetscViewerExodusIIGetNodalVariableIndex(petsclib::PetscLibType, viewer::PetscViewer, name::String, varIndex::PetscExodusIIInt) end

@for_petsc function PetscViewerExodusIIGetNodalVariableIndex(petsclib::$UnionPetscLib, viewer::PetscViewer, name::String, varIndex::PetscExodusIIInt )

    @chk ccall(
               (:PetscViewerExodusIIGetNodalVariableIndex, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cchar}, Ptr{PetscExodusIIInt}),
               viewer, name, varIndex,
              )


	return nothing
end 

"""
	PetscViewerExodusIIGetZonalVariableIndex(petsclib::PetscLibType,viewer::PetscViewer, name::String, varIndex::Cint) 
return the location of a zonal variable in an ExodusII file given its name

Collective

Input Parameters:
- `viewer` - a `PetscViewer` of type `PETSCVIEWEREXODUSII`
- `name`   - the name of the result

Output Parameter:
- `varIndex` - the location of the variable in the exodus file or -1 if the variable is not found

Level: beginner

-seealso: `PetscViewerExodusIISetNodalVariable()`, `PetscViewerExodusIIGetNodalVariable()`, `PetscViewerExodusIISetNodalVariableName()`, `PetscViewerExodusIISetNodalVariableNames()`, `PetscViewerExodusIIGetNodalVariableName()`, `PetscViewerExodusIIGetNodalVariableNames()`

# External Links
$(_doc_external("Dm/PetscViewerExodusIIGetZonalVariableIndex"))
"""
function PetscViewerExodusIIGetZonalVariableIndex(petsclib::PetscLibType, viewer::PetscViewer, name::String, varIndex::Cint) end

@for_petsc function PetscViewerExodusIIGetZonalVariableIndex(petsclib::$UnionPetscLib, viewer::PetscViewer, name::String, varIndex::Cint )

    @chk ccall(
               (:PetscViewerExodusIIGetZonalVariableIndex, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cchar}, Ptr{Cint}),
               viewer, name, varIndex,
              )


	return nothing
end 

"""
	has::PetscBool = PetscViewerHDF5PathIsRelative(petsclib::PetscLibType,path::String, emptyIsRelative::PetscBool) 

# External Links
$(_doc_external("Viewer/PetscViewerHDF5PathIsRelative"))
"""
function PetscViewerHDF5PathIsRelative(petsclib::PetscLibType, path::String, emptyIsRelative::PetscBool) end

@for_petsc function PetscViewerHDF5PathIsRelative(petsclib::$UnionPetscLib, path::String, emptyIsRelative::PetscBool )
	has_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscViewerHDF5PathIsRelative, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, PetscBool, Ptr{PetscBool}),
               path, emptyIsRelative, has_,
              )

	has = has_[]

	return has
end 

"""
	PetscViewerSetType(petsclib::PetscLibType,viewer::PetscViewer, type::PetscViewerType) 
Builds `PetscViewer` for a particular implementation.

Collective

Input Parameters:
- `viewer` - the `PetscViewer` context obtained with `PetscViewerCreate()`
- `type`   - for example, `PETSCVIEWERASCII`

Options Database Key:
- `-viewer_type  <type>` - Sets the type; use -help for a list of available methods (for instance, ascii)

Level: advanced

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewerCreate()`, `PetscViewerGetType()`, `PetscViewerType`, `PetscViewerPushFormat()`

# External Links
$(_doc_external("Sys/PetscViewerSetType"))
"""
function PetscViewerSetType(petsclib::PetscLibType, viewer::PetscViewer, type::PetscViewerType) end

@for_petsc function PetscViewerSetType(petsclib::$UnionPetscLib, viewer::PetscViewer, type::PetscViewerType )

    @chk ccall(
               (:PetscViewerSetType, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscViewerType),
               viewer, type,
              )


	return nothing
end 

"""
	PetscViewerSetFromOptions(petsclib::PetscLibType,viewer::PetscViewer) 
Sets various options for a viewer based on values in the options database.

Collective

Input Parameter:
- `viewer` - the viewer context

Level: intermediate

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewerCreate()`, `PetscViewerSetType()`, `PetscViewerType`

# External Links
$(_doc_external("Sys/PetscViewerSetFromOptions"))
"""
function PetscViewerSetFromOptions(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerSetFromOptions(petsclib::$UnionPetscLib, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerSetFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscViewer,),
               viewer,
              )


	return nothing
end 

"""
	PetscViewerSetOptionsPrefix(petsclib::PetscLibType,viewer::PetscViewer, prefix::String) 
Sets the prefix used for searching for
`PetscViewer` options in the database during `PetscViewerSetFromOptions()`.

Logically Collective

Input Parameters:
- `viewer` - the `PetscViewer` context
- `prefix` - the prefix to prepend to all option names

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewerSetFromOptions()`, `PetscViewerAppendOptionsPrefix()`

# External Links
$(_doc_external("Sys/PetscViewerSetOptionsPrefix"))
"""
function PetscViewerSetOptionsPrefix(petsclib::PetscLibType, viewer::PetscViewer, prefix::String) end

@for_petsc function PetscViewerSetOptionsPrefix(petsclib::$UnionPetscLib, viewer::PetscViewer, prefix::String )

    @chk ccall(
               (:PetscViewerSetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cchar}),
               viewer, prefix,
              )


	return nothing
end 

"""
	PetscViewerSetUp(petsclib::PetscLibType,viewer::PetscViewer) 
Sets up the internal viewer data structures for the later use.

Collective

Input Parameter:
- `viewer` - the `PetscViewer` context

Level: advanced

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewerCreate()`, `PetscViewerDestroy()`

# External Links
$(_doc_external("Sys/PetscViewerSetUp"))
"""
function PetscViewerSetUp(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerSetUp(petsclib::$UnionPetscLib, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerSetUp, $petsc_library),
               PetscErrorCode,
               (PetscViewer,),
               viewer,
              )


	return nothing
end 

"""
	PetscViewersDestroy(petsclib::PetscLibType,v::PetscViewers) 
Destroys a set of `PetscViewer`s created with `PetscViewersCreate()`.

Collective

Input Parameter:
- `v` - the `PetscViewers` to be destroyed.

Level: intermediate

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewerDestroy()`, `PetscViewers`, `PetscViewerSocketOpen()`, `PetscViewerASCIIOpen()`, `PetscViewerCreate()`, `PetscViewerDrawOpen()`, `PetscViewersCreate()`

# External Links
$(_doc_external("Sys/PetscViewersDestroy"))
"""
function PetscViewersDestroy(petsclib::PetscLibType, v::PetscViewers) end

@for_petsc function PetscViewersDestroy(petsclib::$UnionPetscLib, v::PetscViewers )

    @chk ccall(
               (:PetscViewersDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscViewers},),
               v,
              )


	return nothing
end 

"""
	v::PetscViewers = PetscViewersCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Creates a container to hold a set of `PetscViewer`'s. The container is essentially a sparse, growable in length array of `PetscViewer`s

Collective

Input Parameter:
- `comm` - the MPI communicator

Output Parameter:
- `v` - the collection of `PetscViewers`

Level: intermediate

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewers`, `PetscViewerCreate()`, `PetscViewersDestroy()`

# External Links
$(_doc_external("Sys/PetscViewersCreate"))
"""
function PetscViewersCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function PetscViewersCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	v_ = Ref{PetscViewers}()

    @chk ccall(
               (:PetscViewersCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscViewers}),
               comm, v_,
              )

	v = v_[]

	return v
end 

"""
	PetscViewersGetViewer(petsclib::PetscLibType,viewers::PetscViewers, n::PetscInt, viewer::PetscViewer) 
Gets a `PetscViewer` from a `PetscViewers` collection

Collective if the viewer has not previously be obtained.

Input Parameters:
- `viewers` - object created with `PetscViewersCreate()`
- `n`       - number of `PetscViewer` you want

Output Parameter:
- `viewer` - the `PetscViewer`

Level: intermediate

-seealso: [](sec_viewers), `PetscViewer`, `PetscViewers`, `PetscViewersCreate()`, `PetscViewersDestroy()`

# External Links
$(_doc_external("Sys/PetscViewersGetViewer"))
"""
function PetscViewersGetViewer(petsclib::PetscLibType, viewers::PetscViewers, n::PetscInt, viewer::PetscViewer) end

@for_petsc function PetscViewersGetViewer(petsclib::$UnionPetscLib, viewers::PetscViewers, n::$PetscInt, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewersGetViewer, $petsc_library),
               PetscErrorCode,
               (PetscViewers, $PetscInt, Ptr{PetscViewer}),
               viewers, n, viewer,
              )


	return nothing
end 

"""
	PetscViewerSAWsOpen(petsclib::PetscLibType,comm::MPI_Comm, lab::PetscViewer) 
Opens an SAWs `PetscViewer`.

Collective; No Fortran Support

Input Parameter:
- `comm` - the MPI communicator

Output Parameter:
- `lab` - the `PetscViewer`

Options Database Keys:
- `-saws_port <port number>` - port number where you are running SAWs client
- `-xxx_view saws`           - publish the object xxx
- `-xxx_saws_block`          - blocks the program at the end of a critical point (for `KSP` and `SNES` it is the end of a solve) until
the user unblocks the problem with an external tool that access the object with SAWS

Level: advanced

-seealso: [](sec_viewers), `PetscViewerDestroy()`, `PetscViewerStringSPrintf()`, `PETSC_VIEWER_SAWS_()`, `PetscObjectSAWsBlock()`,
`PetscObjectSAWsViewOff()`, `PetscObjectSAWsTakeAccess()`, `PetscObjectSAWsGrantAccess()`

# External Links
$(_doc_external("Sys/PetscViewerSAWsOpen"))
"""
function PetscViewerSAWsOpen(petsclib::PetscLibType, comm::MPI_Comm, lab::PetscViewer) end

@for_petsc function PetscViewerSAWsOpen(petsclib::$UnionPetscLib, comm::MPI_Comm, lab::PetscViewer )

    @chk ccall(
               (:PetscViewerSAWsOpen, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscViewer}),
               comm, lab,
              )


	return nothing
end 

"""
	PetscViewerStringOpen(petsclib::PetscLibType,comm::MPI_Comm, string::String, len::Csize_t, lab::PetscViewer) 
Opens a string as a `PETSCVIEWERSTRING` `PetscViewer`. This is a very
simple `PetscViewer`; information on the object is simply stored into
the string in a fairly nice way.

Collective; No Fortran Support

Input Parameters:
- `comm`   - the communicator
- `string` - the string to use
- `len`    - the string length

Output Parameter:
- `lab` - the `PetscViewer`

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERSTRING`, `PetscViewerDestroy()`, `PetscViewerStringSPrintf()`, `PetscViewerStringGetStringRead()`, `PetscViewerStringSetString()`

# External Links
$(_doc_external("Sys/PetscViewerStringOpen"))
"""
function PetscViewerStringOpen(petsclib::PetscLibType, comm::MPI_Comm, string::String, len::Csize_t, lab::PetscViewer) end

@for_petsc function PetscViewerStringOpen(petsclib::$UnionPetscLib, comm::MPI_Comm, string::String, len::Csize_t, lab::PetscViewer )

    @chk ccall(
               (:PetscViewerStringOpen, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Csize_t, Ptr{PetscViewer}),
               comm, string, len, lab,
              )


	return nothing
end 

"""
	PetscViewerStringGetStringRead(petsclib::PetscLibType,viewer::PetscViewer, string::String, len::Csize_t) 
Returns the string that a `PETSCVIEWERSTRING` uses

Logically Collective

Input Parameter:
- `viewer` - `PETSCVIEWERSTRING` viewer

Output Parameters:
- `string` - the string, optional use `NULL` if you do not need
- `len`    - the length of the string, optional use `NULL` if you do not need it

Level: advanced

-seealso: [](sec_viewers), `PetscViewerStringOpen()`, `PETSCVIEWERSTRING`, `PetscViewerStringSetString()`, `PetscViewerStringSPrintf()`,
`PetscViewerStringSetOwnString()`

# External Links
$(_doc_external("Sys/PetscViewerStringGetStringRead"))
"""
function PetscViewerStringGetStringRead(petsclib::PetscLibType, viewer::PetscViewer, string::String, len::Csize_t) end

@for_petsc function PetscViewerStringGetStringRead(petsclib::$UnionPetscLib, viewer::PetscViewer, string::String, len::Csize_t )
	string_ = Ref(pointer(string))

    @chk ccall(
               (:PetscViewerStringGetStringRead, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Ptr{Cchar}}, Ptr{Csize_t}),
               viewer, string_, len,
              )


	return nothing
end 

"""
	PetscViewerStringSetString(petsclib::PetscLibType,viewer::PetscViewer, string::String, len::Csize_t) 
sets the string that a string viewer will print to

Logically Collective

Input Parameters:
- `viewer` - string viewer you wish to attach string to
- `string` - the string to print data into
- `len`    - the length of the string

Level: advanced

-seealso: [](sec_viewers), `PetscViewerStringOpen()`, `PETSCVIEWERSTRING`, `PetscViewerStringGetStringRead()`, `PetscViewerStringSPrintf()`,
`PetscViewerStringSetOwnString()`

# External Links
$(_doc_external("Sys/PetscViewerStringSetString"))
"""
function PetscViewerStringSetString(petsclib::PetscLibType, viewer::PetscViewer, string::String, len::Csize_t) end

@for_petsc function PetscViewerStringSetString(petsclib::$UnionPetscLib, viewer::PetscViewer, string::String, len::Csize_t )

    @chk ccall(
               (:PetscViewerStringSetString, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cchar}, Csize_t),
               viewer, string, len,
              )


	return nothing
end 

"""
	PetscViewerStringSetOwnString(petsclib::PetscLibType,viewer::PetscViewer) 
tells the viewer that it now owns the string and is responsible for freeing it with `PetscFree()`

Logically Collective

Input Parameter:
- `viewer` - string viewer

Level: advanced

-seealso: [](sec_viewers), `PetscViewerStringOpen()`, `PETSCVIEWERSTRING`, `PetscViewerStringGetStringRead()`, `PetscViewerStringSPrintf()`,
`PetscViewerStringSetString()`

# External Links
$(_doc_external("Sys/PetscViewerStringSetOwnString"))
"""
function PetscViewerStringSetOwnString(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscViewerStringSetOwnString(petsclib::$UnionPetscLib, viewer::PetscViewer )

    @chk ccall(
               (:PetscViewerStringSetOwnString, $petsc_library),
               PetscErrorCode,
               (PetscViewer,),
               viewer,
              )


	return nothing
end 

"""
	PetscViewerSocketOpen(petsclib::PetscLibType,comm::MPI_Comm, machine::String, port::Cint, lab::PetscViewer) 
Opens a connection to a MATLAB or other socket based server.

Collective

Input Parameters:
- `comm`    - the MPI communicator
- `machine` - the machine the server is running on, use `NULL` for the local machine, use "server" to passively wait for
a connection from elsewhere
- `port`    - the port to connect to, use `PETSC_DEFAULT` for the default

Output Parameter:
- `lab` - a context to use when communicating with the server

Options Database Keys:
For use with  `PETSC_VIEWER_SOCKET_WORLD`, `PETSC_VIEWER_SOCKET_SELF`,
`PETSC_VIEWER_SOCKET_()` or if
`NULL` is passed for machine or PETSC_DEFAULT is passed for port
- `-viewer_socket_machine <machine>` - the machine where the socket is available
- `-viewer_socket_port <port>`       - the socket to connect to

Environmental variables:
- `PETSC_VIEWER_SOCKET_MACHINE`   - machine name
- `PETSC_VIEWER_SOCKET_PORT`   - portnumber

Level: intermediate

-seealso: [](sec_viewers), `PETSCVIEWERBINARY`, `PETSCVIEWERSOCKET`, `MatView()`, `VecView()`, `PetscViewerDestroy()`, `PetscViewerCreate()`, `PetscViewerSetType()`,
`PetscViewerSocketSetConnection()`, `PETSC_VIEWER_SOCKET_`, `PETSC_VIEWER_SOCKET_WORLD`,
`PETSC_VIEWER_SOCKET_SELF`, `PetscViewerBinaryWrite()`, `PetscViewerBinaryRead()`, `PetscViewerBinaryWriteStringArray()`,
`PetscBinaryViewerGetDescriptor()`, `PetscMatlabEngineCreate()`

# External Links
$(_doc_external("Sys/PetscViewerSocketOpen"))
"""
function PetscViewerSocketOpen(petsclib::PetscLibType, comm::MPI_Comm, machine::String, port::Cint, lab::PetscViewer) end

@for_petsc function PetscViewerSocketOpen(petsclib::$UnionPetscLib, comm::MPI_Comm, machine::String, port::Cint, lab::PetscViewer )

    @chk ccall(
               (:PetscViewerSocketOpen, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Cint, Ptr{PetscViewer}),
               comm, machine, port, lab,
              )


	return nothing
end 

"""
	PetscViewerSocketSetConnection(petsclib::PetscLibType,v::PetscViewer, machine::String, port::Cint) 
Sets the machine and port that a PETSc socket
viewer is to use

Logically Collective

Input Parameters:
- `v`       - viewer to connect
- `machine` - host to connect to, use `NULL` for the local machine,use "server" to passively wait for
a connection from elsewhere
- `port`    - the port on the machine one is connecting to, use `PETSC_DEFAULT` for default

Level: advanced

-seealso: [](sec_viewers), `PETSCVIEWERMATLAB`, `PETSCVIEWERSOCKET`, `PetscViewerSocketOpen()`

# External Links
$(_doc_external("Sys/PetscViewerSocketSetConnection"))
"""
function PetscViewerSocketSetConnection(petsclib::PetscLibType, v::PetscViewer, machine::String, port::Cint) end

@for_petsc function PetscViewerSocketSetConnection(petsclib::$UnionPetscLib, v::PetscViewer, machine::String, port::Cint )

    @chk ccall(
               (:PetscViewerSocketSetConnection, $petsc_library),
               PetscErrorCode,
               (PetscViewer, Ptr{Cchar}, Cint),
               v, machine, port,
              )


	return nothing
end 

