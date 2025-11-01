# autodefined type arguments for class ------
mutable struct _n_DMAdaptor end
const DMAdaptor = Ptr{_n_DMAdaptor}

# -------------------------------------------------------
# autodefined type arguments for class ------
mutable struct _n_DMSwarmDataField end
const DMSwarmDataField = Ptr{_n_DMSwarmDataField}

# -------------------------------------------------------
# autodefined type arguments for class ------
mutable struct _n_DMSwarmDataBucket end
const DMSwarmDataBucket = Ptr{_n_DMSwarmDataBucket}

mutable struct _n_DMSwarmCellDM end
const DMSwarmCellDM = Ptr{_n_DMSwarmCellDM}

# -------------------------------------------------------
# autodefined type arguments for class ------
mutable struct _n_DMSwarmSort end
const DMSwarmSort = Ptr{_n_DMSwarmSort}

mutable struct _n_DMPlexPoCintQueue end
const DMPlexPoCintQueue = Ptr{_n_DMPlexPoCintQueue}

mutable struct _n_DMNetworkMonitor end
const DMNetworkMonitor = Ptr{_n_DMNetworkMonitor}

mutable struct _n_DMField end
const DMField = Ptr{_n_DMField}

mutable struct _n_DMPlexTransform end
const DMPlexTransform = Ptr{_n_DMPlexTransform}

mutable struct _n_PetscSimplePoCintFn end
const PetscSimplePoCintFn = Ptr{_n_PetscSimplePoCintFn}


# -------------------------------------------------------
# autodefined type arguments for class ------
# -------------------------------------------------------
"""
	DMAdaptorRegisterAll(petsclib::PetscLibType) 
Registers all of the adaptor components in the `DM` package.

Not Collective

Level: advanced

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `DMAdaptorType`, `DMRegisterAll()`, `DMAdaptorRegisterDestroy()`

# External Links
$(_doc_external("Dm/DMAdaptorRegisterAll"))
"""
function DMAdaptorRegisterAll(petsclib::PetscLibType) end

@for_petsc function DMAdaptorRegisterAll(petsclib::$UnionPetscLib)

    @chk ccall(
               (:DMAdaptorRegisterAll, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	DMAdaptorRegisterDestroy(petsclib::PetscLibType) 
This function destroys the registered `DMAdaptorType`. It is called from `PetscFinalize()`.

Not collective

Level: developer

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `DMAdaptorRegisterAll()`, `DMAdaptorType`, `PetscFinalize()`

# External Links
$(_doc_external("Dm/DMAdaptorRegisterDestroy"))
"""
function DMAdaptorRegisterDestroy(petsclib::PetscLibType) end

@for_petsc function DMAdaptorRegisterDestroy(petsclib::$UnionPetscLib)

    @chk ccall(
               (:DMAdaptorRegisterDestroy, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	DMAdaptorMonitorRegister(petsclib::PetscLibType,name::String, vtype::PetscViewerType, format::PetscViewerFormat, monitor::external, create::external, destroy::external) 
Registers a mesh adaptation monitor routine that may be accessed with `DMAdaptorMonitorSetFromOptions()`

Not Collective

Input Parameters:
- `name`    - name of a new monitor routine
- `vtype`   - A `PetscViewerType` for the output
- `format`  - A `PetscViewerFormat` for the output
- `monitor` - Monitor routine
- `create`  - Creation routine, or `NULL`
- `destroy` - Destruction routine, or `NULL`

Level: advanced

-seealso: [](ch_snes), `DMAdaptor`, `DMAdaptorMonitorSet()`, `DMAdaptorMonitorRegisterAll()`, `DMAdaptorMonitorSetFromOptions()`

# External Links
$(_doc_external("Dm/DMAdaptorMonitorRegister"))
"""
function DMAdaptorMonitorRegister(petsclib::PetscLibType, name::String, vtype::PetscViewerType, format::PetscViewerFormat, monitor::external, create::external, destroy::external) end

@for_petsc function DMAdaptorMonitorRegister(petsclib::$UnionPetscLib, name::String, vtype::PetscViewerType, format::PetscViewerFormat, monitor::external, create::external, destroy::external )

    @chk ccall(
               (:DMAdaptorMonitorRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, PetscViewerType, PetscViewerFormat, external, external, external),
               name, vtype, format, monitor, create, destroy,
              )


	return nothing
end 

"""
	DMAdaptorMonitorRegisterDestroy(petsclib::PetscLibType) 
This function destroys the registered monitors for `DMAdaptor`. It is called from `PetscFinalize()`.

Not collective

Level: developer

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `DMAdaptorMonitorRegisterAll()`, `DMAdaptor`, `PetscFinalize()`

# External Links
$(_doc_external("Dm/DMAdaptorMonitorRegisterDestroy"))
"""
function DMAdaptorMonitorRegisterDestroy(petsclib::PetscLibType) end

@for_petsc function DMAdaptorMonitorRegisterDestroy(petsclib::$UnionPetscLib)

    @chk ccall(
               (:DMAdaptorMonitorRegisterDestroy, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	adaptor::DMAdaptor = DMAdaptorCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Create a `DMAdaptor` object. Its purpose is to construct a adaptation `DMLabel` or metric `Vec` that can be used to modify the `DM`.

Collective

Input Parameter:
- `comm` - The communicator for the `DMAdaptor` object

Output Parameter:
- `adaptor` - The `DMAdaptor` object

Level: beginner

See also: 
=== 
`DM`, `DMAdaptor`, `DMAdaptorDestroy()`, `DMAdaptorAdapt()`, `PetscConvEst`, `PetscConvEstCreate()`

# External Links
$(_doc_external("Dm/DMAdaptorCreate"))
"""
function DMAdaptorCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function DMAdaptorCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	adaptor_ = Ref{DMAdaptor}()

    @chk ccall(
               (:DMAdaptorCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{DMAdaptor}),
               comm, adaptor_,
              )

	adaptor = adaptor_[]

	return adaptor
end 

"""
	DMAdaptorDestroy(petsclib::PetscLibType,adaptor::DMAdaptor) 
Destroys a `DMAdaptor` object

Collective

Input Parameter:
- `adaptor` - The `DMAdaptor` object

Level: beginner

See also: 
=== 
`DM`, `DMAdaptor`, `DMAdaptorCreate()`, `DMAdaptorAdapt()`

# External Links
$(_doc_external("Dm/DMAdaptorDestroy"))
"""
function DMAdaptorDestroy(petsclib::PetscLibType, adaptor::DMAdaptor) end

@for_petsc function DMAdaptorDestroy(petsclib::$UnionPetscLib, adaptor::DMAdaptor )

    @chk ccall(
               (:DMAdaptorDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{DMAdaptor},),
               adaptor,
              )


	return nothing
end 

"""
	DMAdaptorSetType(petsclib::PetscLibType,adaptor::DMAdaptor, method::DMAdaptorType) 
Sets the particular implementation for a adaptor.

Collective

Input Parameters:
- `adaptor` - The `DMAdaptor`
- `method`  - The name of the adaptor type

Options Database Key:
- `-adaptor_type <type>` - Sets the adaptor type; see `DMAdaptorType`

Level: intermediate

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `DMAdaptor`, `DMAdaptorType`, `DMAdaptorGetType()`, `DMAdaptorCreate()`

# External Links
$(_doc_external("Dm/DMAdaptorSetType"))
"""
function DMAdaptorSetType(petsclib::PetscLibType, adaptor::DMAdaptor, method::DMAdaptorType) end

@for_petsc function DMAdaptorSetType(petsclib::$UnionPetscLib, adaptor::DMAdaptor, method::DMAdaptorType )

    @chk ccall(
               (:DMAdaptorSetType, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, DMAdaptorType),
               adaptor, method,
              )


	return nothing
end 

"""
	type::DMAdaptorType = DMAdaptorGetType(petsclib::PetscLibType,adaptor::DMAdaptor) 
Gets the type name (as a string) from the adaptor.

Not Collective

Input Parameter:
- `adaptor` - The `DMAdaptor`

Output Parameter:
- `type` - The `DMAdaptorType` name

Level: intermediate

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `DMAdaptor`, `DMAdaptorType`, `DMAdaptorSetType()`, `DMAdaptorCreate()`

# External Links
$(_doc_external("Dm/DMAdaptorGetType"))
"""
function DMAdaptorGetType(petsclib::PetscLibType, adaptor::DMAdaptor) end

@for_petsc function DMAdaptorGetType(petsclib::$UnionPetscLib, adaptor::DMAdaptor )
	type_ = Ref{DMAdaptorType}()

    @chk ccall(
               (:DMAdaptorGetType, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, Ptr{DMAdaptorType}),
               adaptor, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	DMAdaptorMonitorSet(petsclib::PetscLibType,adaptor::DMAdaptor, monitor::external, ctx::Cvoid, monitordestroy::PetscCtxDestroyFn) 
Sets an ADDITIONAL function to be called at every iteration to monitor
the error etc.

Logically Collective

Input Parameters:
- `adaptor`        - the `DMAdaptor`
- `monitor`        - pointer to function (if this is `NULL`, it turns off monitoring
- `ctx`            - [optional] context for private data for the monitor routine (use `NULL` if no context is needed)
- `monitordestroy` - [optional] routine that frees monitor context (may be `NULL`), see `PetscCtxDestroyFn` for its calling sequence

Calling sequence of `monitor`:
- `adaptor` - the `DMAdaptor`
- `it`      - iteration number
- `odm`     - the original `DM`
- `adm`     - the adapted `DM`
- `Nf`      - number of fields
- `enorms`  - (estimated) 2-norm of the error for each field
- `error`   - `Vec` of cellwise errors
- `ctx`     - optional monitoring context, as set by `DMAdaptorMonitorSet()`

Options Database Keys:
- `-adaptor_monitor_size`                - sets `DMAdaptorMonitorSize()`
- `-adaptor_monitor_error`               - sets `DMAdaptorMonitorError()`
- `-adaptor_monitor_error draw`          - sets `DMAdaptorMonitorErrorDraw()` and plots error
- `-adaptor_monitor_error draw::draw_lg` - sets `DMAdaptorMonitorErrorDrawLG()` and plots error
- `-dm_adaptor_monitor_cancel`           - Cancels all monitors that have been hardwired into a code by calls to `DMAdaptorMonitorSet()`, but does not cancel those set via the options database.

Level: beginner

-seealso: [](ch_snes), `DMAdaptorMonitorError()`, `DMAdaptor`, `PetscCtxDestroyFn`

# External Links
$(_doc_external("Dm/DMAdaptorMonitorSet"))
"""
function DMAdaptorMonitorSet(petsclib::PetscLibType, adaptor::DMAdaptor, monitor::external, ctx::Cvoid, monitordestroy::PetscCtxDestroyFn) end

@for_petsc function DMAdaptorMonitorSet(petsclib::$UnionPetscLib, adaptor::DMAdaptor, monitor::external, ctx::Cvoid, monitordestroy::PetscCtxDestroyFn )

    @chk ccall(
               (:DMAdaptorMonitorSet, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, external, Ptr{Cvoid}, Ptr{PetscCtxDestroyFn}),
               adaptor, monitor, ctx, monitordestroy,
              )


	return nothing
end 

"""
	DMAdaptorMonitorCancel(petsclib::PetscLibType,adaptor::DMAdaptor) 
Clears all monitors for a `DMAdaptor` object.

Logically Collective

Input Parameter:
- `adaptor` - the `DMAdaptor`

Options Database Key:
- `-dm_adaptor_monitor_cancel` - Cancels all monitors that have been hardwired into a code by calls to `DMAdaptorMonitorSet()`, but does not cancel those set via the options database.

Level: intermediate

-seealso: [](ch_snes), `DMAdaptorMonitorError()`, `DMAdaptorMonitorSet()`, `DMAdaptor`

# External Links
$(_doc_external("Dm/DMAdaptorMonitorCancel"))
"""
function DMAdaptorMonitorCancel(petsclib::PetscLibType, adaptor::DMAdaptor) end

@for_petsc function DMAdaptorMonitorCancel(petsclib::$UnionPetscLib, adaptor::DMAdaptor )

    @chk ccall(
               (:DMAdaptorMonitorCancel, $petsc_library),
               PetscErrorCode,
               (DMAdaptor,),
               adaptor,
              )


	return nothing
end 

"""
	DMAdaptorMonitorSetFromOptions(petsclib::PetscLibType,adaptor::DMAdaptor, opt::String, name::String, ctx::Cvoid) 
Sets a monitor function and viewer appropriate for the type indicated by the user in the options database

Collective

Input Parameters:
- `adaptor` - `DMadaptor` object you wish to monitor
- `opt`     - the command line option for this monitor
- `name`    - the monitor type one is seeking
- `ctx`     - An optional user context for the monitor, or `NULL`

Level: developer

-seealso: [](ch_snes), `DMAdaptorMonitorRegister()`, `DMAdaptorMonitorSet()`, `PetscOptionsGetViewer()`

# External Links
$(_doc_external("Dm/DMAdaptorMonitorSetFromOptions"))
"""
function DMAdaptorMonitorSetFromOptions(petsclib::PetscLibType, adaptor::DMAdaptor, opt::String, name::String, ctx::Cvoid) end

@for_petsc function DMAdaptorMonitorSetFromOptions(petsclib::$UnionPetscLib, adaptor::DMAdaptor, opt::String, name::String, ctx::Cvoid )

    @chk ccall(
               (:DMAdaptorMonitorSetFromOptions, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cvoid}),
               adaptor, opt, name, ctx,
              )


	return nothing
end 

"""
	DMAdaptorSetOptionsPrefix(petsclib::PetscLibType,adaptor::DMAdaptor, prefix::String) 
Sets the prefix used for searching for all `DMAdaptor` options in the database.

Logically Collective

Input Parameters:
- `adaptor` - the `DMAdaptor`
- `prefix`  - the prefix to prepend to all option names

Level: advanced

-seealso: [](ch_snes), `DMAdaptor`, `SNESSetOptionsPrefix()`, `DMAdaptorSetFromOptions()`

# External Links
$(_doc_external("Dm/DMAdaptorSetOptionsPrefix"))
"""
function DMAdaptorSetOptionsPrefix(petsclib::PetscLibType, adaptor::DMAdaptor, prefix::String) end

@for_petsc function DMAdaptorSetOptionsPrefix(petsclib::$UnionPetscLib, adaptor::DMAdaptor, prefix::String )

    @chk ccall(
               (:DMAdaptorSetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, Ptr{Cchar}),
               adaptor, prefix,
              )


	return nothing
end 

"""
	DMAdaptorSetFromOptions(petsclib::PetscLibType,adaptor::DMAdaptor) 
Sets properties of a `DMAdaptor` object from values in the options database

Collective

Input Parameter:
- `adaptor` - The `DMAdaptor` object

Options Database Keys:
- `-adaptor_monitor_size`                - Monitor the mesh size
- `-adaptor_monitor_error`               - Monitor the solution error
- `-adaptor_sequence_num <num>`          - Number of adaptations to generate an optimal grid
- `-adaptor_target_num <num>`            - Set the target number of vertices N_adapt, -1 for automatic determination
- `-adaptor_refinement_factor <r>`       - Set r such that N_adapt = r^dim N_orig
- `-adaptor_mixed_setup_function <func>` - Set the function func that sets up the mixed problem

Level: beginner

See also: 
=== 
`DM`, `DMAdaptor`, `DMAdaptorCreate()`, `DMAdaptorAdapt()`

# External Links
$(_doc_external("Dm/DMAdaptorSetFromOptions"))
"""
function DMAdaptorSetFromOptions(petsclib::PetscLibType, adaptor::DMAdaptor) end

@for_petsc function DMAdaptorSetFromOptions(petsclib::$UnionPetscLib, adaptor::DMAdaptor )

    @chk ccall(
               (:DMAdaptorSetFromOptions, $petsc_library),
               PetscErrorCode,
               (DMAdaptor,),
               adaptor,
              )


	return nothing
end 

"""
	DMAdaptorView(petsclib::PetscLibType,adaptor::DMAdaptor, viewer::PetscViewer) 
Views a `DMAdaptor` object

Collective

Input Parameters:
- `adaptor` - The `DMAdaptor` object
- `viewer`  - The `PetscViewer` object

Level: beginner

See also: 
=== 
`DM`, `DMAdaptor`, `DMAdaptorCreate()`, `DMAdaptorAdapt()`

# External Links
$(_doc_external("Dm/DMAdaptorView"))
"""
function DMAdaptorView(petsclib::PetscLibType, adaptor::DMAdaptor, viewer::PetscViewer) end

@for_petsc function DMAdaptorView(petsclib::$UnionPetscLib, adaptor::DMAdaptor, viewer::PetscViewer )

    @chk ccall(
               (:DMAdaptorView, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, PetscViewer),
               adaptor, viewer,
              )


	return nothing
end 

"""
	DMAdaptorGetSolver(petsclib::PetscLibType,adaptor::DMAdaptor, snes::PetscSNES) 
Gets the solver used to produce discrete solutions

Not Collective

Input Parameter:
- `adaptor` - The `DMAdaptor` object

Output Parameter:
- `snes` - The solver

Level: intermediate

See also: 
=== 
`DM`, `DMAdaptor`, `DMAdaptorSetSolver()`, `DMAdaptorCreate()`, `DMAdaptorAdapt()`

# External Links
$(_doc_external("Dm/DMAdaptorGetSolver"))
"""
function DMAdaptorGetSolver(petsclib::PetscLibType, adaptor::DMAdaptor, snes::PetscSNES) end

@for_petsc function DMAdaptorGetSolver(petsclib::$UnionPetscLib, adaptor::DMAdaptor, snes::PetscSNES )
	snes_ = Ref(snes.ptr)

    @chk ccall(
               (:DMAdaptorGetSolver, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, Ptr{CSNES}),
               adaptor, snes_,
              )

	snes.ptr = C_NULL

	return nothing
end 

"""
	DMAdaptorSetSolver(petsclib::PetscLibType,adaptor::DMAdaptor, snes::PetscSNES) 
Sets the solver used to produce discrete solutions

Not Collective

Input Parameters:
- `adaptor` - The `DMAdaptor` object
- `snes`    - The solver, this MUST have an attached `DM`/`PetscDS`, so that the exact solution can be computed

Level: intermediate

See also: 
=== 
`DMAdaptor`, `DMAdaptorGetSolver()`, `DMAdaptorCreate()`, `DMAdaptorAdapt()`

# External Links
$(_doc_external("Dm/DMAdaptorSetSolver"))
"""
function DMAdaptorSetSolver(petsclib::PetscLibType, adaptor::DMAdaptor, snes::PetscSNES) end

@for_petsc function DMAdaptorSetSolver(petsclib::$UnionPetscLib, adaptor::DMAdaptor, snes::PetscSNES )

    @chk ccall(
               (:DMAdaptorSetSolver, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, CSNES),
               adaptor, snes,
              )


	return nothing
end 

"""
	num::PetscInt = DMAdaptorGetSequenceLength(petsclib::PetscLibType,adaptor::DMAdaptor) 
Gets the number of sequential adaptations used by an adapter

Not Collective

Input Parameter:
- `adaptor` - The `DMAdaptor` object

Output Parameter:
- `num` - The number of adaptations

Level: intermediate

See also: 
=== 
`DMAdaptor`, `DMAdaptorSetSequenceLength()`, `DMAdaptorCreate()`, `DMAdaptorAdapt()`

# External Links
$(_doc_external("Dm/DMAdaptorGetSequenceLength"))
"""
function DMAdaptorGetSequenceLength(petsclib::PetscLibType, adaptor::DMAdaptor) end

@for_petsc function DMAdaptorGetSequenceLength(petsclib::$UnionPetscLib, adaptor::DMAdaptor )
	num_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMAdaptorGetSequenceLength, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, Ptr{$PetscInt}),
               adaptor, num_,
              )

	num = num_[]

	return num
end 

"""
	DMAdaptorSetSequenceLength(petsclib::PetscLibType,adaptor::DMAdaptor, num::PetscInt) 
Sets the number of sequential adaptations

Not Collective

Input Parameters:
- `adaptor` - The `DMAdaptor` object
- `num`     - The number of adaptations

Level: intermediate

See also: 
=== 
`DMAdaptorGetSequenceLength()`, `DMAdaptorCreate()`, `DMAdaptorAdapt()`

# External Links
$(_doc_external("Dm/DMAdaptorSetSequenceLength"))
"""
function DMAdaptorSetSequenceLength(petsclib::PetscLibType, adaptor::DMAdaptor, num::PetscInt) end

@for_petsc function DMAdaptorSetSequenceLength(petsclib::$UnionPetscLib, adaptor::DMAdaptor, num::$PetscInt )

    @chk ccall(
               (:DMAdaptorSetSequenceLength, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, $PetscInt),
               adaptor, num,
              )


	return nothing
end 

"""
	DMAdaptorSetUp(petsclib::PetscLibType,adaptor::DMAdaptor) 
After the solver is specified, creates data structures for controlling adaptivity

Collective

Input Parameter:
- `adaptor` - The `DMAdaptor` object

Level: beginner

See also: 
=== 
`DMAdaptor`, `DMAdaptorCreate()`, `DMAdaptorAdapt()`

# External Links
$(_doc_external("Dm/DMAdaptorSetUp"))
"""
function DMAdaptorSetUp(petsclib::PetscLibType, adaptor::DMAdaptor) end

@for_petsc function DMAdaptorSetUp(petsclib::$UnionPetscLib, adaptor::DMAdaptor )

    @chk ccall(
               (:DMAdaptorSetUp, $petsc_library),
               PetscErrorCode,
               (DMAdaptor,),
               adaptor,
              )


	return nothing
end 

"""
	DMAdaptorSetTransferFunction(petsclib::PetscLibType,adaptor::DMAdaptor, tfunc::external) 

# External Links
$(_doc_external("Dm/DMAdaptorSetTransferFunction"))
"""
function DMAdaptorSetTransferFunction(petsclib::PetscLibType, adaptor::DMAdaptor, tfunc::external) end

@for_petsc function DMAdaptorSetTransferFunction(petsclib::$UnionPetscLib, adaptor::DMAdaptor, tfunc::external )

    @chk ccall(
               (:DMAdaptorSetTransferFunction, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, external),
               adaptor, tfunc,
              )


	return nothing
end 

"""
	DMAdaptorMonitor(petsclib::PetscLibType,adaptor::DMAdaptor, it::PetscInt, odm::PetscDM, adm::PetscDM, Nf::PetscInt, enorms::Vector{PetscReal}, error::PetscVec) 
runs the user provided monitor routines, if they exist

Collective

Input Parameters:
- `adaptor` - the `DMAdaptor`
- `it`      - iteration number
- `odm`     - the original `DM`
- `adm`     - the adapted `DM`
- `Nf`      - the number of fields
- `enorms`  - the 2-norm error values for each field
- `error`   - `Vec` of cellwise errors

Level: developer

-seealso: [](ch_snes), `DMAdaptorMonitorSet()`

# External Links
$(_doc_external("Dm/DMAdaptorMonitor"))
"""
function DMAdaptorMonitor(petsclib::PetscLibType, adaptor::DMAdaptor, it::PetscInt, odm::PetscDM, adm::PetscDM, Nf::PetscInt, enorms::Vector{PetscReal}, error::PetscVec) end

@for_petsc function DMAdaptorMonitor(petsclib::$UnionPetscLib, adaptor::DMAdaptor, it::$PetscInt, odm::PetscDM, adm::PetscDM, Nf::$PetscInt, enorms::Vector{$PetscReal}, error::PetscVec )

    @chk ccall(
               (:DMAdaptorMonitor, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, $PetscInt, CDM, CDM, $PetscInt, Ptr{$PetscReal}, CVec),
               adaptor, it, odm, adm, Nf, enorms, error,
              )


	return nothing
end 

"""
	DMAdaptorMonitorSize(petsclib::PetscLibType,adaptor::DMAdaptor, n::PetscInt, odm::PetscDM, adm::PetscDM, Nf::PetscInt, enorms::Vector{PetscReal}, error::PetscVec, vf::PetscViewerAndFormat) 
Prints the mesh sizes at each iteration of an adaptation loop.

Collective

Input Parameters:
- `adaptor` - the `DMAdaptor`
- `n`       - iteration number
- `odm`     - the original `DM`
- `adm`     - the adapted `DM`
- `Nf`      - number of fields
- `enorms`  - 2-norm error values for each field (may be estimated).
- `error`   - `Vec` of cellwise errors
- `vf`      - The viewer context

Options Database Key:
- `-adaptor_monitor_size` - Activates `DMAdaptorMonitorSize()`

Level: intermediate

-seealso: [](ch_snes), `DMAdaptor`, `DMAdaptorMonitorSet()`, `DMAdaptorMonitorError()`, `DMAdaptorMonitorErrorDraw()`, `DMAdaptorMonitorErrorDrawLG()`

# External Links
$(_doc_external("Dm/DMAdaptorMonitorSize"))
"""
function DMAdaptorMonitorSize(petsclib::PetscLibType, adaptor::DMAdaptor, n::PetscInt, odm::PetscDM, adm::PetscDM, Nf::PetscInt, enorms::Vector{PetscReal}, error::PetscVec, vf::PetscViewerAndFormat) end

@for_petsc function DMAdaptorMonitorSize(petsclib::$UnionPetscLib, adaptor::DMAdaptor, n::$PetscInt, odm::PetscDM, adm::PetscDM, Nf::$PetscInt, enorms::Vector{$PetscReal}, error::PetscVec, vf::PetscViewerAndFormat )

    @chk ccall(
               (:DMAdaptorMonitorSize, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, $PetscInt, CDM, CDM, $PetscInt, Ptr{$PetscReal}, CVec, Ptr{PetscViewerAndFormat}),
               adaptor, n, odm, adm, Nf, enorms, error, vf,
              )


	return nothing
end 

"""
	DMAdaptorMonitorError(petsclib::PetscLibType,adaptor::DMAdaptor, n::PetscInt, odm::PetscDM, adm::PetscDM, Nf::PetscInt, enorms::Vector{PetscReal}, error::PetscVec, vf::PetscViewerAndFormat) 
Prints the error norm at each iteration of an adaptation loop.

Collective

Input Parameters:
- `adaptor` - the `DMAdaptor`
- `n`       - iteration number
- `odm`     - the original `DM`
- `adm`     - the adapted `DM`
- `Nf`      - number of fields
- `enorms`  - 2-norm error values for each field (may be estimated).
- `error`   - `Vec` of cellwise errors
- `vf`      - The viewer context

Options Database Key:
- `-adaptor_monitor_error` - Activates `DMAdaptorMonitorError()`

Level: intermediate

-seealso: [](ch_snes), `DMAdaptor`, `DMAdaptorMonitorSet()`, `DMAdaptorMonitorErrorDraw()`, `DMAdaptorMonitorErrorDrawLG()`

# External Links
$(_doc_external("Dm/DMAdaptorMonitorError"))
"""
function DMAdaptorMonitorError(petsclib::PetscLibType, adaptor::DMAdaptor, n::PetscInt, odm::PetscDM, adm::PetscDM, Nf::PetscInt, enorms::Vector{PetscReal}, error::PetscVec, vf::PetscViewerAndFormat) end

@for_petsc function DMAdaptorMonitorError(petsclib::$UnionPetscLib, adaptor::DMAdaptor, n::$PetscInt, odm::PetscDM, adm::PetscDM, Nf::$PetscInt, enorms::Vector{$PetscReal}, error::PetscVec, vf::PetscViewerAndFormat )

    @chk ccall(
               (:DMAdaptorMonitorError, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, $PetscInt, CDM, CDM, $PetscInt, Ptr{$PetscReal}, CVec, Ptr{PetscViewerAndFormat}),
               adaptor, n, odm, adm, Nf, enorms, error, vf,
              )


	return nothing
end 

"""
	DMAdaptorMonitorErrorDraw(petsclib::PetscLibType,adaptor::DMAdaptor, n::PetscInt, odm::PetscDM, adm::PetscDM, Nf::PetscInt, enorms::Vector{PetscReal}, error::PetscVec, vf::PetscViewerAndFormat) 
Plots the error at each iteration of an iterative solver.

Collective

Input Parameters:
- `adaptor` - the `DMAdaptor`
- `n`       - iteration number
- `odm`     - the original `DM`
- `adm`     - the adapted `DM`
- `Nf`      - number of fields
- `enorms`  - 2-norm error values for each field (may be estimated).
- `error`   - `Vec` of cellwise errors
- `vf`      - The viewer context

Options Database Key:
- `-adaptor_monitor_error draw` - Activates `DMAdaptorMonitorErrorDraw()`

Level: intermediate

-seealso: [](ch_snes), `PETSCVIEWERDRAW`, `DMAdaptor`, `DMAdaptorMonitorSet()`, `DMAdaptorMonitorErrorDrawLG()`

# External Links
$(_doc_external("Dm/DMAdaptorMonitorErrorDraw"))
"""
function DMAdaptorMonitorErrorDraw(petsclib::PetscLibType, adaptor::DMAdaptor, n::PetscInt, odm::PetscDM, adm::PetscDM, Nf::PetscInt, enorms::Vector{PetscReal}, error::PetscVec, vf::PetscViewerAndFormat) end

@for_petsc function DMAdaptorMonitorErrorDraw(petsclib::$UnionPetscLib, adaptor::DMAdaptor, n::$PetscInt, odm::PetscDM, adm::PetscDM, Nf::$PetscInt, enorms::Vector{$PetscReal}, error::PetscVec, vf::PetscViewerAndFormat )

    @chk ccall(
               (:DMAdaptorMonitorErrorDraw, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, $PetscInt, CDM, CDM, $PetscInt, Ptr{$PetscReal}, CVec, Ptr{PetscViewerAndFormat}),
               adaptor, n, odm, adm, Nf, enorms, error, vf,
              )


	return nothing
end 

"""
	vf::PetscViewerAndFormat = DMAdaptorMonitorErrorDrawLGCreate(petsclib::PetscLibType,viewer::PetscViewer, format::PetscViewerFormat, ctx::Cvoid) 
Creates the context for the error plotter `DMAdaptorMonitorErrorDrawLG()`

Collective

Input Parameters:
- `viewer` - The `PetscViewer`
- `format` - The viewer format
- `ctx`    - An optional user context

Output Parameter:
- `vf` - The viewer context

Level: intermediate

-seealso: [](ch_snes), `PETSCVIEWERDRAW`, `PetscViewerMonitorGLSetUp()`, `DMAdaptor`, `DMAdaptorMonitorSet()`, `DMAdaptorMonitorErrorDrawLG()`

# External Links
$(_doc_external("Dm/DMAdaptorMonitorErrorDrawLGCreate"))
"""
function DMAdaptorMonitorErrorDrawLGCreate(petsclib::PetscLibType, viewer::PetscViewer, format::PetscViewerFormat, ctx::Cvoid) end

@for_petsc function DMAdaptorMonitorErrorDrawLGCreate(petsclib::$UnionPetscLib, viewer::PetscViewer, format::PetscViewerFormat, ctx::Cvoid )
	vf_ = Ref{PetscViewerAndFormat}()

    @chk ccall(
               (:DMAdaptorMonitorErrorDrawLGCreate, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscViewerFormat, Ptr{Cvoid}, PetscViewerAndFormat),
               viewer, format, ctx, vf_,
              )

	vf = vf_[]

	return vf
end 

"""
	DMAdaptorMonitorErrorDrawLG(petsclib::PetscLibType,adaptor::DMAdaptor, n::PetscInt, odm::PetscDM, adm::PetscDM, Nf::PetscInt, enorms::Vector{PetscReal}, error::PetscVec, vf::PetscViewerAndFormat) 
Plots the error norm at each iteration of an adaptive loop.

Collective

Input Parameters:
- `adaptor` - the `DMAdaptor`
- `n`       - iteration number
- `odm`     - the original `DM`
- `adm`     - the adapted `DM`
- `Nf`      - number of fields
- `enorms`  - 2-norm error values for each field (may be estimated).
- `error`   - `Vec` of cellwise errors
- `vf`      - The viewer context, obtained via `DMAdaptorMonitorErrorDrawLGCreate()`

Options Database Key:
- `-adaptor_error draw::draw_lg` - Activates `DMAdaptorMonitorErrorDrawLG()`

Level: intermediate

-seealso: [](ch_snes), `PETSCVIEWERDRAW`, `DMAdaptor`, `DMAdaptorMonitorSet()`, `DMAdaptorMonitorErrorDraw()`, `DMAdaptorMonitorError()`,
`DMAdaptorMonitorTrueResidualDrawLGCreate()`

# External Links
$(_doc_external("Dm/DMAdaptorMonitorErrorDrawLG"))
"""
function DMAdaptorMonitorErrorDrawLG(petsclib::PetscLibType, adaptor::DMAdaptor, n::PetscInt, odm::PetscDM, adm::PetscDM, Nf::PetscInt, enorms::Vector{PetscReal}, error::PetscVec, vf::PetscViewerAndFormat) end

@for_petsc function DMAdaptorMonitorErrorDrawLG(petsclib::$UnionPetscLib, adaptor::DMAdaptor, n::$PetscInt, odm::PetscDM, adm::PetscDM, Nf::$PetscInt, enorms::Vector{$PetscReal}, error::PetscVec, vf::PetscViewerAndFormat )

    @chk ccall(
               (:DMAdaptorMonitorErrorDrawLG, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, $PetscInt, CDM, CDM, $PetscInt, Ptr{$PetscReal}, CVec, Ptr{PetscViewerAndFormat}),
               adaptor, n, odm, adm, Nf, enorms, error, vf,
              )


	return nothing
end 

"""
	DMAdaptorMonitorRegisterAll(petsclib::PetscLibType) 
Registers all of the mesh adaptation monitors in the `SNES` package.

Not Collective

Level: advanced

-seealso: [](ch_snes), `SNES`, `DM`, `DMAdaptorMonitorRegister()`, `DMAdaptorRegister()`

# External Links
$(_doc_external("Dm/DMAdaptorMonitorRegisterAll"))
"""
function DMAdaptorMonitorRegisterAll(petsclib::PetscLibType) end

@for_petsc function DMAdaptorMonitorRegisterAll(petsclib::$UnionPetscLib)

    @chk ccall(
               (:DMAdaptorMonitorRegisterAll, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	DMAdaptorAdapt(petsclib::PetscLibType,adaptor::DMAdaptor, x::PetscVec, strategy::DMAdaptationStrategy, adm::PetscDM, ax::PetscVec) 
Creates a new `DM` that is adapted to the problem

Not Collective

Input Parameters:
- `adaptor`  - The `DMAdaptor` object
- `x`        - The global approximate solution
- `strategy` - The adaptation strategy, see `DMAdaptationStrategy`

Output Parameters:
- `adm` - The adapted `DM`
- `ax`  - The adapted solution

Options Database Keys:
- `-snes_adapt <strategy>` - initial, sequential, multigrid
- `-adapt_gradient_view`   - View the Clement interpolant of the solution gradient
- `-adapt_hessian_view`    - View the Clement interpolant of the solution Hessian
- `-adapt_metric_view`     - View the metric tensor for adaptive mesh refinement

Level: intermediate

See also: 
=== 
`DMAdaptor`, `DMAdaptationStrategy`, `DMAdaptorSetSolver()`, `DMAdaptorCreate()`

# External Links
$(_doc_external("Dm/DMAdaptorAdapt"))
"""
function DMAdaptorAdapt(petsclib::PetscLibType, adaptor::DMAdaptor, x::PetscVec, strategy::DMAdaptationStrategy, adm::PetscDM, ax::PetscVec) end

@for_petsc function DMAdaptorAdapt(petsclib::$UnionPetscLib, adaptor::DMAdaptor, x::PetscVec, strategy::DMAdaptationStrategy, adm::PetscDM, ax::PetscVec )
	adm_ = Ref(adm.ptr)
	ax_ = Ref(ax.ptr)

    @chk ccall(
               (:DMAdaptorAdapt, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, CVec, DMAdaptationStrategy, Ptr{CDM}, Ptr{CVec}),
               adaptor, x, strategy, adm_, ax_,
              )

	adm.ptr = C_NULL
	ax.ptr = C_NULL

	return nothing
end 

"""
	DMAdaptorSetMixedSetupFunction(petsclib::PetscLibType,adaptor::DMAdaptor, setupFunc::external) 
Set the function setting up the mixed problem

Not Collective

Input Parameters:
- `adaptor`   - the `DMAdaptor`
- `setupFunc` - the function setting up the mixed problem

Level: advanced

-seealso: `DMAdaptor`, `DMAdaptorGetMixedSetupFunction()`, `DMAdaptorAdapt()`

# External Links
$(_doc_external("Dm/DMAdaptorSetMixedSetupFunction"))
"""
function DMAdaptorSetMixedSetupFunction(petsclib::PetscLibType, adaptor::DMAdaptor, setupFunc::external) end

@for_petsc function DMAdaptorSetMixedSetupFunction(petsclib::$UnionPetscLib, adaptor::DMAdaptor, setupFunc::external )

    @chk ccall(
               (:DMAdaptorSetMixedSetupFunction, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, external),
               adaptor, setupFunc,
              )


	return nothing
end 

"""
	DMAdaptorGetCriterion(petsclib::PetscLibType,adaptor::DMAdaptor, criterion::DMAdaptationCriterion) 
Get the adaptation criterion

Not Collective

Input Parameter:
- `adaptor` - the `DMAdaptor`

Output Parameter:
- `criterion` - the criterion for adaptation

Level: advanced

-seealso: `DMAdaptor`, `DMAdaptorSetCriterion()`, `DMAdaptationCriterion`

# External Links
$(_doc_external("Dm/DMAdaptorGetCriterion"))
"""
function DMAdaptorGetCriterion(petsclib::PetscLibType, adaptor::DMAdaptor, criterion::DMAdaptationCriterion) end

@for_petsc function DMAdaptorGetCriterion(petsclib::$UnionPetscLib, adaptor::DMAdaptor, criterion::DMAdaptationCriterion )

    @chk ccall(
               (:DMAdaptorGetCriterion, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, Ptr{DMAdaptationCriterion}),
               adaptor, criterion,
              )


	return nothing
end 

"""
	DMAdaptorSetCriterion(petsclib::PetscLibType,adaptor::DMAdaptor, criterion::DMAdaptationCriterion) 
Set the adaptation criterion

Not Collective

Input Parameters:
- `adaptor`   - the `DMAdaptor`
- `criterion` - the adaptation criterion

Level: advanced

-seealso: `DMAdaptor`, `DMAdaptorGetCriterion()`, `DMAdaptationCriterion`

# External Links
$(_doc_external("Dm/DMAdaptorSetCriterion"))
"""
function DMAdaptorSetCriterion(petsclib::PetscLibType, adaptor::DMAdaptor, criterion::DMAdaptationCriterion) end

@for_petsc function DMAdaptorSetCriterion(petsclib::$UnionPetscLib, adaptor::DMAdaptor, criterion::DMAdaptationCriterion )

    @chk ccall(
               (:DMAdaptorSetCriterion, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, DMAdaptationCriterion),
               adaptor, criterion,
              )


	return nothing
end 

"""
	DMSwarmDataFieldGetEntries(petsclib::PetscLibType,gfield::DMSwarmDataField, data::Cvoid) 

# External Links
$(_doc_external("Dm/DMSwarmDataFieldGetEntries"))
"""
function DMSwarmDataFieldGetEntries(petsclib::PetscLibType, gfield::DMSwarmDataField, data::Cvoid) end

@for_petsc function DMSwarmDataFieldGetEntries(petsclib::$UnionPetscLib, gfield::DMSwarmDataField, data::Cvoid )

    @chk ccall(
               (:DMSwarmDataFieldGetEntries, $petsc_library),
               PetscErrorCode,
               (DMSwarmDataField, Cvoid),
               gfield, data,
              )


	return nothing
end 

"""
	DMSwarmDataFieldRestoreEntries(petsclib::PetscLibType,gfield::DMSwarmDataField, data::Cvoid) 

# External Links
$(_doc_external("Dm/DMSwarmDataFieldRestoreEntries"))
"""
function DMSwarmDataFieldRestoreEntries(petsclib::PetscLibType, gfield::DMSwarmDataField, data::Cvoid) end

@for_petsc function DMSwarmDataFieldRestoreEntries(petsclib::$UnionPetscLib, gfield::DMSwarmDataField, data::Cvoid )

    @chk ccall(
               (:DMSwarmDataFieldRestoreEntries, $petsc_library),
               PetscErrorCode,
               (DMSwarmDataField, Cvoid),
               gfield, data,
              )


	return nothing
end 

"""
	idx::PetscInt = DMSwarmDataBucketGetDMSwarmDataFieldIdByName(petsclib::PetscLibType,db::DMSwarmDataBucket, name::String) 

# External Links
$(_doc_external("Dm/DMSwarmDataBucketGetDMSwarmDataFieldIdByName"))
"""
function DMSwarmDataBucketGetDMSwarmDataFieldIdByName(petsclib::PetscLibType, db::DMSwarmDataBucket, name::String) end

@for_petsc function DMSwarmDataBucketGetDMSwarmDataFieldIdByName(petsclib::$UnionPetscLib, db::DMSwarmDataBucket, name::String )
	idx_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMSwarmDataBucketGetDMSwarmDataFieldIdByName, $petsc_library),
               PetscErrorCode,
               (DMSwarmDataBucket, Ptr{Cchar}, Ptr{$PetscInt}),
               db, name, idx_,
              )

	idx = idx_[]

	return idx
end 

"""
	DMSwarmDataBucketGetDMSwarmDataFieldByName(petsclib::PetscLibType,db::DMSwarmDataBucket, name::String, gfield::DMSwarmDataField) 

# External Links
$(_doc_external("Dm/DMSwarmDataBucketGetDMSwarmDataFieldByName"))
"""
function DMSwarmDataBucketGetDMSwarmDataFieldByName(petsclib::PetscLibType, db::DMSwarmDataBucket, name::String, gfield::DMSwarmDataField) end

@for_petsc function DMSwarmDataBucketGetDMSwarmDataFieldByName(petsclib::$UnionPetscLib, db::DMSwarmDataBucket, name::String, gfield::DMSwarmDataField )

    @chk ccall(
               (:DMSwarmDataBucketGetDMSwarmDataFieldByName, $petsc_library),
               PetscErrorCode,
               (DMSwarmDataBucket, Ptr{Cchar}, Ptr{DMSwarmDataField}),
               db, name, gfield,
              )


	return nothing
end 

"""
	found::PetscBool = DMSwarmDataBucketQueryDMSwarmDataFieldByName(petsclib::PetscLibType,db::DMSwarmDataBucket, name::String) 

# External Links
$(_doc_external("Dm/DMSwarmDataBucketQueryDMSwarmDataFieldByName"))
"""
function DMSwarmDataBucketQueryDMSwarmDataFieldByName(petsclib::PetscLibType, db::DMSwarmDataBucket, name::String) end

@for_petsc function DMSwarmDataBucketQueryDMSwarmDataFieldByName(petsclib::$UnionPetscLib, db::DMSwarmDataBucket, name::String )
	found_ = Ref{PetscBool}()

    @chk ccall(
               (:DMSwarmDataBucketQueryDMSwarmDataFieldByName, $petsc_library),
               PetscErrorCode,
               (DMSwarmDataBucket, Ptr{Cchar}, Ptr{PetscBool}),
               db, name, found_,
              )

	found = found_[]

	return found
end 

"""
	DMSwarmSortDestroy(petsclib::PetscLibType,_ctx::DMSwarmSort) 

# External Links
$(_doc_external("Dm/DMSwarmSortDestroy"))
"""
function DMSwarmSortDestroy(petsclib::PetscLibType, _ctx::DMSwarmSort) end

@for_petsc function DMSwarmSortDestroy(petsclib::$UnionPetscLib, _ctx::DMSwarmSort )

    @chk ccall(
               (:DMSwarmSortDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{DMSwarmSort},),
               _ctx,
              )


	return nothing
end 

"""
	npoints::PetscInt = DMSwarmSortGetNumberOfPointsPerCell(petsclib::PetscLibType,sw::PetscDM, cell::PetscInt) 
Returns the number of points in a cell

Not Collective

Input Parameters:
- `sw`   - a `DMSWARM` objects
- `cell` - the cell number in the cell `DM`

Output Parameter:
- `npoints` - the number of points in the cell

Level: advanced

-seealso: `DMSWARM`, `DMSwarmSetType()`, `DMSwarmSortGetAccess()`, `DMSwarmSortGetPointsPerCell()`

# External Links
$(_doc_external("Dm/DMSwarmSortGetNumberOfPointsPerCell"))
"""
function DMSwarmSortGetNumberOfPointsPerCell(petsclib::PetscLibType, sw::PetscDM, cell::PetscInt) end

@for_petsc function DMSwarmSortGetNumberOfPointsPerCell(petsclib::$UnionPetscLib, sw::PetscDM, cell::$PetscInt )
	npoints_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMSwarmSortGetNumberOfPointsPerCell, $petsc_library),
               PetscErrorCode,
               (CDM, $PetscInt, Ptr{$PetscInt}),
               sw, cell, npoints_,
              )

	npoints = npoints_[]

	return npoints
end 

"""
	DMSwarmSortGetPointsPerCell(petsclib::PetscLibType,sw::PetscDM, cell::PetscInt, npoints::PetscInt, pidlist::PetscInt) 
Creates an array of point indices for all points in a cell

Not Collective

Input Parameters:
- `sw`      - a `DMSWARM` object
- `cell`    - the cell number in the cell `DM`
- `npoints` - the number of points in the cell
- `pidlist` - array of the indices identifying all points in cell e

Level: advanced

-seealso: `DMSWARM`, `DMSwarmSetType()`, `DMSwarmRestorePointsPerCell()`, `DMSwarmSortGetAccess()`, `DMSwarmSortGetNumberOfPointsPerCell()`

# External Links
$(_doc_external("Dm/DMSwarmSortGetPointsPerCell"))
"""
function DMSwarmSortGetPointsPerCell(petsclib::PetscLibType, sw::PetscDM, cell::PetscInt, npoints::PetscInt, pidlist::PetscInt) end

@for_petsc function DMSwarmSortGetPointsPerCell(petsclib::$UnionPetscLib, sw::PetscDM, cell::$PetscInt, npoints::$PetscInt, pidlist::$PetscInt )

    @chk ccall(
               (:DMSwarmSortGetPointsPerCell, $petsc_library),
               PetscErrorCode,
               (CDM, $PetscInt, Ptr{$PetscInt}, $PetscInt),
               sw, cell, npoints, pidlist,
              )


	return nothing
end 

"""
	DMSwarmSortRestorePointsPerCell(petsclib::PetscLibType,dm::PetscDM, e::PetscInt, npoints::PetscInt, pidlist::PetscInt) 
Restores an array of point indices for all points in a cell

Not Collective

Input Parameters:
- `dm`      - a `DMSWARM` object
- `e`       - the index of the cell
- `npoints` - the number of points in the cell
- `pidlist` - array of the indices identifying all points in cell e

Level: advanced

-seealso: `DMSWARM`, `DMSwarmSetType()`, `DMSwarmSortGetPointsPerCell()`, `DMSwarmSortGetAccess()`, `DMSwarmSortGetNumberOfPointsPerCell()`

# External Links
$(_doc_external("Dm/DMSwarmSortRestorePointsPerCell"))
"""
function DMSwarmSortRestorePointsPerCell(petsclib::PetscLibType, dm::PetscDM, e::PetscInt, npoints::PetscInt, pidlist::PetscInt) end

@for_petsc function DMSwarmSortRestorePointsPerCell(petsclib::$UnionPetscLib, dm::PetscDM, e::$PetscInt, npoints::$PetscInt, pidlist::$PetscInt )

    @chk ccall(
               (:DMSwarmSortRestorePointsPerCell, $petsc_library),
               PetscErrorCode,
               (CDM, $PetscInt, Ptr{$PetscInt}, $PetscInt),
               dm, e, npoints, pidlist,
              )


	return nothing
end 

"""
	DMSwarmSortGetAccess(petsclib::PetscLibType,sw::PetscDM) 
Setups up a `DMSWARM` point sort context for efficient traversal of points within a cell

Not Collective

Input Parameter:
- `sw` - a `DMSWARM` object

Level: advanced

-seealso: `DMSWARM`, `DMSwarmSetType()`, `DMSwarmSortRestoreAccess()`

# External Links
$(_doc_external("Dm/DMSwarmSortGetAccess"))
"""
function DMSwarmSortGetAccess(petsclib::PetscLibType, sw::PetscDM) end

@for_petsc function DMSwarmSortGetAccess(petsclib::$UnionPetscLib, sw::PetscDM )

    @chk ccall(
               (:DMSwarmSortGetAccess, $petsc_library),
               PetscErrorCode,
               (CDM,),
               sw,
              )


	return nothing
end 

"""
	DMSwarmSortRestoreAccess(petsclib::PetscLibType,sw::PetscDM) 
Invalidates the `DMSWARM` point sorting context previously computed with `DMSwarmSortGetAccess()`

Not Collective

Input Parameter:
- `sw` - a `DMSWARM` object

Level: advanced

-seealso: `DMSWARM`, `DMSwarmSetType()`, `DMSwarmSortGetAccess()`

# External Links
$(_doc_external("Dm/DMSwarmSortRestoreAccess"))
"""
function DMSwarmSortRestoreAccess(petsclib::PetscLibType, sw::PetscDM) end

@for_petsc function DMSwarmSortRestoreAccess(petsclib::$UnionPetscLib, sw::PetscDM )

    @chk ccall(
               (:DMSwarmSortRestoreAccess, $petsc_library),
               PetscErrorCode,
               (CDM,),
               sw,
              )


	return nothing
end 

"""
	isvalid::PetscBool = DMSwarmSortGetIsValid(petsclib::PetscLibType,sw::PetscDM) 
Gets the isvalid flag associated with a `DMSWARM` point sorting context

Not Collective

Input Parameter:
- `sw` - a `DMSWARM` object

Output Parameter:
- `isvalid` - flag indicating whether the sort context is up-to-date

Level: advanced

-seealso: `DMSWARM`, `DMSwarmSetType()`, `DMSwarmSortGetAccess()`

# External Links
$(_doc_external("Dm/DMSwarmSortGetIsValid"))
"""
function DMSwarmSortGetIsValid(petsclib::PetscLibType, sw::PetscDM) end

@for_petsc function DMSwarmSortGetIsValid(petsclib::$UnionPetscLib, sw::PetscDM )
	isvalid_ = Ref{PetscBool}()

    @chk ccall(
               (:DMSwarmSortGetIsValid, $petsc_library),
               PetscErrorCode,
               (CDM, Ptr{PetscBool}),
               sw, isvalid_,
              )

	isvalid = isvalid_[]

	return isvalid
end 

"""
	ncells::PetscInt,npoints::PetscInt = DMSwarmSortGetSizes(petsclib::PetscLibType,sw::PetscDM) 
Gets the sizes associated with a `DMSWARM` point sorting context

Not Collective

Input Parameter:
- `sw` - a `DMSWARM` object

Output Parameters:
- `ncells`  - number of cells within the sort context (pass `NULL` to ignore)
- `npoints` - number of points used to create the sort context (pass `NULL` to ignore)

Level: advanced

-seealso: `DMSWARM`, `DMSwarmSetType()`, `DMSwarmSortGetAccess()`

# External Links
$(_doc_external("Dm/DMSwarmSortGetSizes"))
"""
function DMSwarmSortGetSizes(petsclib::PetscLibType, sw::PetscDM) end

@for_petsc function DMSwarmSortGetSizes(petsclib::$UnionPetscLib, sw::PetscDM )
	ncells_ = Ref{$PetscInt}()
	npoints_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMSwarmSortGetSizes, $petsc_library),
               PetscErrorCode,
               (CDM, Ptr{$PetscInt}, Ptr{$PetscInt}),
               sw, ncells_, npoints_,
              )

	ncells = ncells_[]
	npoints = npoints_[]

	return ncells,npoints
end 

"""
	DMSwarmCellDMDestroy(petsclib::PetscLibType,celldm::DMSwarmCellDM) 
destroy a `DMSwarmCellDM`

Collective

Input Parameter:
- `celldm` - address of `DMSwarmCellDM`

Level: advanced

-seealso: `DMSwarmCellDM`, `DMSwarmCellDMCreate()`

# External Links
$(_doc_external("Dm/DMSwarmCellDMDestroy"))
"""
function DMSwarmCellDMDestroy(petsclib::PetscLibType, celldm::DMSwarmCellDM) end

@for_petsc function DMSwarmCellDMDestroy(petsclib::$UnionPetscLib, celldm::DMSwarmCellDM )

    @chk ccall(
               (:DMSwarmCellDMDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{DMSwarmCellDM},),
               celldm,
              )


	return nothing
end 

"""
	DMSwarmCellDMView(petsclib::PetscLibType,celldm::DMSwarmCellDM, viewer::PetscViewer) 
view a `DMSwarmCellDM`

Collective

Input Parameters:
- `celldm` - `DMSwarmCellDM`
- `viewer` - viewer to display field, for example `PETSC_VIEWER_STDOUT_WORLD`

Level: advanced

-seealso: `DMSwarmCellDM`, `DMSwarmCellDMCreate()`

# External Links
$(_doc_external("Dm/DMSwarmCellDMView"))
"""
function DMSwarmCellDMView(petsclib::PetscLibType, celldm::DMSwarmCellDM, viewer::PetscViewer) end

@for_petsc function DMSwarmCellDMView(petsclib::$UnionPetscLib, celldm::DMSwarmCellDM, viewer::PetscViewer )

    @chk ccall(
               (:DMSwarmCellDMView, $petsc_library),
               PetscErrorCode,
               (DMSwarmCellDM, PetscViewer),
               celldm, viewer,
              )


	return nothing
end 

"""
	DMSwarmCellDMGetDM(petsclib::PetscLibType,celldm::DMSwarmCellDM, dm::PetscDM) 
Returns the background `DM` for the `DMSwarm`

Not Collective

Input Parameter:
- `celldm` - The `DMSwarmCellDM` object

Output Parameter:
- `dm` - The `DM` object

Level: intermediate

-seealso: `DMSwarmCellDM`, `DM`, `DMSwarmSetCellDM()`

# External Links
$(_doc_external("Dm/DMSwarmCellDMGetDM"))
"""
function DMSwarmCellDMGetDM(petsclib::PetscLibType, celldm::DMSwarmCellDM, dm::PetscDM) end

@for_petsc function DMSwarmCellDMGetDM(petsclib::$UnionPetscLib, celldm::DMSwarmCellDM, dm::PetscDM )
	dm_ = Ref(dm.ptr)

    @chk ccall(
               (:DMSwarmCellDMGetDM, $petsc_library),
               PetscErrorCode,
               (DMSwarmCellDM, Ptr{CDM}),
               celldm, dm_,
              )

	dm.ptr = C_NULL

	return nothing
end 

"""
	Nf::PetscInt = DMSwarmCellDMGetFields(petsclib::PetscLibType,celldm::DMSwarmCellDM, names::String) 
Returns the `DM` fields for the `DMSwarm`

Not Collective

Input Parameter:
- `celldm` - The `DMSwarmCellDM` object

Output Parameters:
- `Nf`    - The number of fields
- `names` - The array of field names in the `DMSWARM`

Level: intermediate

-seealso: `DMSwarmCellDM`, `DM`, `DMSwarmSetCellDM()`

# External Links
$(_doc_external("Dm/DMSwarmCellDMGetFields"))
"""
function DMSwarmCellDMGetFields(petsclib::PetscLibType, celldm::DMSwarmCellDM, names::String) end

@for_petsc function DMSwarmCellDMGetFields(petsclib::$UnionPetscLib, celldm::DMSwarmCellDM, names::String )
	Nf_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMSwarmCellDMGetFields, $petsc_library),
               PetscErrorCode,
               (DMSwarmCellDM, Ptr{$PetscInt}, Ptr{Cchar}),
               celldm, Nf_, names,
              )

	Nf = Nf_[]

	return Nf
end 

"""
	Nfc::PetscInt = DMSwarmCellDMGetCoordinateFields(petsclib::PetscLibType,celldm::DMSwarmCellDM, names::String) 
Returns the `DM` coordinate fields for the `DMSwarm`

Not Collective

Input Parameter:
- `celldm` - The `DMSwarmCellDM` object

Output Parameters:
- `Nfc`   - The number of coordinate fields
- `names` - The array of coordinate field names in the `DMSWARM`

Level: intermediate

-seealso: `DMSwarmCellDM`, `DM`, `DMSwarmSetCellDM()`

# External Links
$(_doc_external("Dm/DMSwarmCellDMGetCoordinateFields"))
"""
function DMSwarmCellDMGetCoordinateFields(petsclib::PetscLibType, celldm::DMSwarmCellDM, names::String) end

@for_petsc function DMSwarmCellDMGetCoordinateFields(petsclib::$UnionPetscLib, celldm::DMSwarmCellDM, names::String )
	Nfc_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMSwarmCellDMGetCoordinateFields, $petsc_library),
               PetscErrorCode,
               (DMSwarmCellDM, Ptr{$PetscInt}, Ptr{Cchar}),
               celldm, Nfc_, names,
              )

	Nfc = Nfc_[]

	return Nfc
end 

"""
	DMSwarmCellDMGetCellID(petsclib::PetscLibType,celldm::DMSwarmCellDM, cellid::String) 
Returns the cell id field name for the `DMSwarm`

Not Collective

Input Parameter:
- `celldm` - The `DMSwarmCellDM` object

Output Parameters:
- `cellid` - The cell id field name in the `DMSWARM`

Level: intermediate

-seealso: `DMSwarmCellDM`, `DM`, `DMSwarmSetCellDM()`

# External Links
$(_doc_external("Dm/DMSwarmCellDMGetCellID"))
"""
function DMSwarmCellDMGetCellID(petsclib::PetscLibType, celldm::DMSwarmCellDM, cellid::String) end

@for_petsc function DMSwarmCellDMGetCellID(petsclib::$UnionPetscLib, celldm::DMSwarmCellDM, cellid::String )
	cellid_ = Ref(pointer(cellid))

    @chk ccall(
               (:DMSwarmCellDMGetCellID, $petsc_library),
               PetscErrorCode,
               (DMSwarmCellDM, Ptr{Ptr{Cchar}}),
               celldm, cellid_,
              )


	return nothing
end 

"""
	DMSwarmCellDMGetSort(petsclib::PetscLibType,celldm::DMSwarmCellDM, sort::DMSwarmSort) 
Returns the sort context over the active `DMSwarmCellDM` for the `DMSwarm`

Not Collective

Input Parameter:
- `celldm` - The `DMSwarmCellDM` object

Output Parameter:
- `sort` - The `DMSwarmSort` object

Level: intermediate

-seealso: `DMSwarmCellDM`, `DM`, `DMSwarmSetCellDM()`

# External Links
$(_doc_external("Dm/DMSwarmCellDMGetSort"))
"""
function DMSwarmCellDMGetSort(petsclib::PetscLibType, celldm::DMSwarmCellDM, sort::DMSwarmSort) end

@for_petsc function DMSwarmCellDMGetSort(petsclib::$UnionPetscLib, celldm::DMSwarmCellDM, sort::DMSwarmSort )

    @chk ccall(
               (:DMSwarmCellDMGetSort, $petsc_library),
               PetscErrorCode,
               (DMSwarmCellDM, Ptr{DMSwarmSort}),
               celldm, sort,
              )


	return nothing
end 

"""
	DMSwarmCellDMSetSort(petsclib::PetscLibType,celldm::DMSwarmCellDM, sort::DMSwarmSort) 
Sets the sort context over the active `DMSwarmCellDM` for the `DMSwarm`

Not Collective

Input Parameters:
- `celldm` - The `DMSwarmCellDM` object
- `sort`   - The `DMSwarmSort` object

Level: intermediate

-seealso: `DMSwarmCellDM`, `DM`, `DMSwarmSetCellDM()`

# External Links
$(_doc_external("Dm/DMSwarmCellDMSetSort"))
"""
function DMSwarmCellDMSetSort(petsclib::PetscLibType, celldm::DMSwarmCellDM, sort::DMSwarmSort) end

@for_petsc function DMSwarmCellDMSetSort(petsclib::$UnionPetscLib, celldm::DMSwarmCellDM, sort::DMSwarmSort )

    @chk ccall(
               (:DMSwarmCellDMSetSort, $petsc_library),
               PetscErrorCode,
               (DMSwarmCellDM, DMSwarmSort),
               celldm, sort,
              )


	return nothing
end 

"""
	bs::PetscInt = DMSwarmCellDMGetBlockSize(petsclib::PetscLibType,celldm::DMSwarmCellDM, sw::PetscDM) 
Returns the total blocksize for the `DM` fields

Not Collective

Input Parameters:
- `celldm` - The `DMSwarmCellDM` object
- `sw`     - The `DMSwarm` object

Output Parameter:
- `bs` - The total block size

Level: intermediate

-seealso: `DMSwarmCellDM`, `DM`, `DMSwarmSetCellDM()`

# External Links
$(_doc_external("Dm/DMSwarmCellDMGetBlockSize"))
"""
function DMSwarmCellDMGetBlockSize(petsclib::PetscLibType, celldm::DMSwarmCellDM, sw::PetscDM) end

@for_petsc function DMSwarmCellDMGetBlockSize(petsclib::$UnionPetscLib, celldm::DMSwarmCellDM, sw::PetscDM )
	bs_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMSwarmCellDMGetBlockSize, $petsc_library),
               PetscErrorCode,
               (DMSwarmCellDM, CDM, Ptr{$PetscInt}),
               celldm, sw, bs_,
              )

	bs = bs_[]

	return bs
end 

"""
	celldm::DMSwarmCellDM = DMSwarmCellDMCreate(petsclib::PetscLibType,dm::PetscDM, Nf::PetscInt, dmFields::String, Nfc::PetscInt, coordFields::String) 
create a `DMSwarmCellDM`

Collective

Input Parameters:
- `dm`          - The background `DM` for the `DMSwarm`
- `Nf`          - The number of swarm fields defined over `dm`
- `dmFields`    - The swarm field names for the `dm` fields
- `Nfc`         - The number of swarm fields to use for coordinates over `dm`
- `coordFields` - The swarm field names for the `dm` coordinate fields

Output Parameter:
- `celldm` - The new `DMSwarmCellDM`

Level: advanced

-seealso: `DMSwarmCellDM`, `DMSWARM`, `DMSetType()`

# External Links
$(_doc_external("Dm/DMSwarmCellDMCreate"))
"""
function DMSwarmCellDMCreate(petsclib::PetscLibType, dm::PetscDM, Nf::PetscInt, dmFields::String, Nfc::PetscInt, coordFields::String) end

@for_petsc function DMSwarmCellDMCreate(petsclib::$UnionPetscLib, dm::PetscDM, Nf::$PetscInt, dmFields::String, Nfc::$PetscInt, coordFields::String )
	dmFields_ = Ref(pointer(dmFields))
	coordFields_ = Ref(pointer(coordFields))
	celldm_ = Ref{DMSwarmCellDM}()

    @chk ccall(
               (:DMSwarmCellDMCreate, $petsc_library),
               PetscErrorCode,
               (CDM, $PetscInt, Ptr{Ptr{Cchar}}, $PetscInt, Ptr{Ptr{Cchar}}, Ptr{DMSwarmCellDM}),
               dm, Nf, dmFields_, Nfc, coordFields_, celldm_,
              )

	celldm = celldm_[]

	return celldm
end 

"""
	queue::DMPlexPoCintQueue = DMPlexPointQueueCreate(petsclib::PetscLibType,size::PetscInt) 

# External Links
$(_doc_external("Dm/DMPlexPointQueueCreate"))
"""
function DMPlexPointQueueCreate(petsclib::PetscLibType, size::PetscInt) end

@for_petsc function DMPlexPointQueueCreate(petsclib::$UnionPetscLib, size::$PetscInt )
	queue_ = Ref{DMPlexPoCintQueue}()

    @chk ccall(
               (:DMPlexPointQueueCreate, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{DMPlexPoCintQueue}),
               size, queue_,
              )

	queue = queue_[]

	return queue
end 

"""
	DMPlexPointQueueDestroy(petsclib::PetscLibType,queue::DMPlexPoCintQueue) 

# External Links
$(_doc_external("Dm/DMPlexPointQueueDestroy"))
"""
function DMPlexPointQueueDestroy(petsclib::PetscLibType, queue::DMPlexPoCintQueue) end

@for_petsc function DMPlexPointQueueDestroy(petsclib::$UnionPetscLib, queue::DMPlexPoCintQueue )

    @chk ccall(
               (:DMPlexPointQueueDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{DMPlexPoCintQueue},),
               queue,
              )


	return nothing
end 

"""
	DMPlexPointQueueEnsureSize(petsclib::PetscLibType,queue::DMPlexPoCintQueue) 

# External Links
$(_doc_external("Dm/DMPlexPointQueueEnsureSize"))
"""
function DMPlexPointQueueEnsureSize(petsclib::PetscLibType, queue::DMPlexPoCintQueue) end

@for_petsc function DMPlexPointQueueEnsureSize(petsclib::$UnionPetscLib, queue::DMPlexPoCintQueue )

    @chk ccall(
               (:DMPlexPointQueueEnsureSize, $petsc_library),
               PetscErrorCode,
               (DMPlexPoCintQueue,),
               queue,
              )


	return nothing
end 

"""
	DMPlexPointQueueEnqueue(petsclib::PetscLibType,queue::DMPlexPoCintQueue, p::PetscInt) 

# External Links
$(_doc_external("Dm/DMPlexPointQueueEnqueue"))
"""
function DMPlexPointQueueEnqueue(petsclib::PetscLibType, queue::DMPlexPoCintQueue, p::PetscInt) end

@for_petsc function DMPlexPointQueueEnqueue(petsclib::$UnionPetscLib, queue::DMPlexPoCintQueue, p::$PetscInt )

    @chk ccall(
               (:DMPlexPointQueueEnqueue, $petsc_library),
               PetscErrorCode,
               (DMPlexPoCintQueue, $PetscInt),
               queue, p,
              )


	return nothing
end 

"""
	p::PetscInt = DMPlexPointQueueDequeue(petsclib::PetscLibType,queue::DMPlexPoCintQueue) 

# External Links
$(_doc_external("Dm/DMPlexPointQueueDequeue"))
"""
function DMPlexPointQueueDequeue(petsclib::PetscLibType, queue::DMPlexPoCintQueue) end

@for_petsc function DMPlexPointQueueDequeue(petsclib::$UnionPetscLib, queue::DMPlexPoCintQueue )
	p_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMPlexPointQueueDequeue, $petsc_library),
               PetscErrorCode,
               (DMPlexPoCintQueue, Ptr{$PetscInt}),
               queue, p_,
              )

	p = p_[]

	return p
end 

"""
	p::PetscInt = DMPlexPointQueueFront(petsclib::PetscLibType,queue::DMPlexPoCintQueue) 

# External Links
$(_doc_external("Dm/DMPlexPointQueueFront"))
"""
function DMPlexPointQueueFront(petsclib::PetscLibType, queue::DMPlexPoCintQueue) end

@for_petsc function DMPlexPointQueueFront(petsclib::$UnionPetscLib, queue::DMPlexPoCintQueue )
	p_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMPlexPointQueueFront, $petsc_library),
               PetscErrorCode,
               (DMPlexPoCintQueue, Ptr{$PetscInt}),
               queue, p_,
              )

	p = p_[]

	return p
end 

"""
	p::PetscInt = DMPlexPointQueueBack(petsclib::PetscLibType,queue::DMPlexPoCintQueue) 

# External Links
$(_doc_external("Dm/DMPlexPointQueueBack"))
"""
function DMPlexPointQueueBack(petsclib::PetscLibType, queue::DMPlexPoCintQueue) end

@for_petsc function DMPlexPointQueueBack(petsclib::$UnionPetscLib, queue::DMPlexPoCintQueue )
	p_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMPlexPointQueueBack, $petsc_library),
               PetscErrorCode,
               (DMPlexPoCintQueue, Ptr{$PetscInt}),
               queue, p_,
              )

	p = p_[]

	return p
end 

"""
	empty::PetscBool = DMPlexPointQueueEmptyCollective(petsclib::PetscLibType,obj::PetscObject, queue::DMPlexPoCintQueue) 

# External Links
$(_doc_external("Dm/DMPlexPointQueueEmptyCollective"))
"""
function DMPlexPointQueueEmptyCollective(petsclib::PetscLibType, obj::PetscObject, queue::DMPlexPoCintQueue) end

@for_petsc function DMPlexPointQueueEmptyCollective(petsclib::$UnionPetscLib, obj::PetscObject, queue::DMPlexPoCintQueue )
	empty_ = Ref{PetscBool}()

    @chk ccall(
               (:DMPlexPointQueueEmptyCollective, $petsc_library),
               PetscErrorCode,
               (PetscObject, DMPlexPoCintQueue, Ptr{PetscBool}),
               obj, queue, empty_,
              )

	empty = empty_[]

	return empty
end 

"""
	DMFieldDestroy(petsclib::PetscLibType,field::DMField) 
destroy a `DMField`

Collective

Input Parameter:
- `field` - address of `DMField`

Level: advanced

-seealso: `DMField`, `DMFieldCreate()`

# External Links
$(_doc_external("Dm/DMFieldDestroy"))
"""
function DMFieldDestroy(petsclib::PetscLibType, field::DMField) end

@for_petsc function DMFieldDestroy(petsclib::$UnionPetscLib, field::DMField )

    @chk ccall(
               (:DMFieldDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{DMField},),
               field,
              )


	return nothing
end 

"""
	DMFieldView(petsclib::PetscLibType,field::DMField, viewer::PetscViewer) 
view a `DMField`

Collective

Input Parameters:
- `field`  - `DMField`
- `viewer` - viewer to display field, for example `PETSC_VIEWER_STDOUT_WORLD`

Level: advanced

-seealso: `DMField`, `DMFieldCreate()`

# External Links
$(_doc_external("Dm/DMFieldView"))
"""
function DMFieldView(petsclib::PetscLibType, field::DMField, viewer::PetscViewer) end

@for_petsc function DMFieldView(petsclib::$UnionPetscLib, field::DMField, viewer::PetscViewer )

    @chk ccall(
               (:DMFieldView, $petsc_library),
               PetscErrorCode,
               (DMField, PetscViewer),
               field, viewer,
              )


	return nothing
end 

"""
	DMFieldSetType(petsclib::PetscLibType,field::DMField, type::DMFieldType) 
set the `DMField` implementation

Collective

Input Parameters:
- `field` - the `DMField` context
- `type`  - a known method

Level: advanced

-seealso: `DMField`, `DMFieldGetType()`, `DMFieldType`,

# External Links
$(_doc_external("Dm/DMFieldSetType"))
"""
function DMFieldSetType(petsclib::PetscLibType, field::DMField, type::DMFieldType) end

@for_petsc function DMFieldSetType(petsclib::$UnionPetscLib, field::DMField, type::DMFieldType )

    @chk ccall(
               (:DMFieldSetType, $petsc_library),
               PetscErrorCode,
               (DMField, DMFieldType),
               field, type,
              )


	return nothing
end 

"""
	type::DMFieldType = DMFieldGetType(petsclib::PetscLibType,field::DMField) 
Gets the `DMFieldType` name (as a string) from the `DMField`.

Not Collective

Input Parameter:
- `field` - The `DMField` context

Output Parameter:
- `type` - The `DMFieldType` name

Level: advanced

-seealso: `DMField`, `DMFieldSetType()`, `DMFieldType`

# External Links
$(_doc_external("Dm/DMFieldGetType"))
"""
function DMFieldGetType(petsclib::PetscLibType, field::DMField) end

@for_petsc function DMFieldGetType(petsclib::$UnionPetscLib, field::DMField )
	type_ = Ref{DMFieldType}()

    @chk ccall(
               (:DMFieldGetType, $petsc_library),
               PetscErrorCode,
               (DMField, Ptr{DMFieldType}),
               field, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	nc::PetscInt = DMFieldGetNumComponents(petsclib::PetscLibType,field::DMField) 
Returns the number of components in the field

Not Collective

Input Parameter:
- `field` - The `DMField` object

Output Parameter:
- `nc` - The number of field components

Level: intermediate

-seealso: `DMField`, `DMFieldEvaluate()`

# External Links
$(_doc_external("Dm/DMFieldGetNumComponents"))
"""
function DMFieldGetNumComponents(petsclib::PetscLibType, field::DMField) end

@for_petsc function DMFieldGetNumComponents(petsclib::$UnionPetscLib, field::DMField )
	nc_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMFieldGetNumComponents, $petsc_library),
               PetscErrorCode,
               (DMField, Ptr{$PetscInt}),
               field, nc_,
              )

	nc = nc_[]

	return nc
end 

"""
	DMFieldGetDM(petsclib::PetscLibType,field::DMField, dm::PetscDM) 
Returns the `DM` for the manifold over which the field is defined.

Not Collective

Input Parameter:
- `field` - The `DMField` object

Output Parameter:
- `dm` - The `DM` object

Level: intermediate

-seealso: `DMField`, `DM`, `DMFieldEvaluate()`

# External Links
$(_doc_external("Dm/DMFieldGetDM"))
"""
function DMFieldGetDM(petsclib::PetscLibType, field::DMField, dm::PetscDM) end

@for_petsc function DMFieldGetDM(petsclib::$UnionPetscLib, field::DMField, dm::PetscDM )
	dm_ = Ref(dm.ptr)

    @chk ccall(
               (:DMFieldGetDM, $petsc_library),
               PetscErrorCode,
               (DMField, Ptr{CDM}),
               field, dm_,
              )

	dm.ptr = C_NULL

	return nothing
end 

"""
	DMFieldEvaluate(petsclib::PetscLibType,field::DMField, points::PetscVec, datatype::PetscDataType, B::Cvoid, D::Cvoid, H::Cvoid) 
Evaluate the field and its derivatives on a set of points

Collective

Input Parameters:
- `field`    - The `DMField` object
- `points`   - The points at which to evaluate the field.  Should have size d x n,
where d is the coordinate dimension of the manifold and n is the number
of points
- `datatype` - The PetscDataType of the output arrays: either `PETSC_REAL` or `PETSC_SCALAR`.
If the field is complex and datatype is `PETSC_REAL`, the real part of the
field is returned.

Output Parameters:
- `B` - pointer to data of size c * n * sizeof(datatype), where c is the number of components in the field.
If B is not `NULL`, the values of the field are written in this array, varying first by component,
then by point.
- `D` - pointer to data of size d * c * n * sizeof(datatype).
If `D` is not `NULL`, the values of the field's spatial derivatives are written in this array,
varying first by the partial derivative component, then by field component, then by point.
- `H` - pointer to data of size d * d * c * n * sizeof(datatype).
If `H` is not `NULL`, the values of the field's second spatial derivatives are written in this array,
varying first by the second partial derivative component, then by field component, then by point.

Level: intermediate

-seealso: `DMField`, `DMFieldGetDM()`, `DMFieldGetNumComponents()`, `DMFieldEvaluateFE()`, `DMFieldEvaluateFV()`, `PetscDataType`

# External Links
$(_doc_external("Dm/DMFieldEvaluate"))
"""
function DMFieldEvaluate(petsclib::PetscLibType, field::DMField, points::PetscVec, datatype::PetscDataType, B::Cvoid, D::Cvoid, H::Cvoid) end

@for_petsc function DMFieldEvaluate(petsclib::$UnionPetscLib, field::DMField, points::PetscVec, datatype::PetscDataType, B::Cvoid, D::Cvoid, H::Cvoid )

    @chk ccall(
               (:DMFieldEvaluate, $petsc_library),
               PetscErrorCode,
               (DMField, CVec, PetscDataType, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
               field, points, datatype, B, D, H,
              )


	return nothing
end 

"""
	DMFieldEvaluateFE(petsclib::PetscLibType,field::DMField, cellIS::IS, points::PetscQuadrature, datatype::PetscDataType, B::Cvoid, D::Cvoid, H::Cvoid) 
Evaluate the field and its derivatives on a set of points mapped from
quadrature points on a reference point.  The derivatives are taken with respect to the
reference coordinates.

Not Collective

Input Parameters:
- `field`    - The `DMField` object
- `cellIS`   - Index set for cells on which to evaluate the field
- `points`   - The quadature containing the points in the reference cell at which to evaluate the field.
- `datatype` - The PetscDataType of the output arrays: either `PETSC_REAL` or `PETSC_SCALAR`.
If the field is complex and datatype is `PETSC_REAL`, the real part of the
field is returned.

Output Parameters:
- `B` - pointer to data of size c * n * sizeof(datatype), where c is the number of components in the field.
If B is not `NULL`, the values of the field are written in this array, varying first by component,
then by point.
- `D` - pointer to data of size d * c * n * sizeof(datatype).
If D is not `NULL`, the values of the field's spatial derivatives are written in this array,
varying first by the partial derivative component, then by field component, then by point.
- `H` - pointer to data of size d * d * c * n * sizeof(datatype).
If H is not `NULL`, the values of the field's second spatial derivatives are written in this array,
varying first by the second partial derivative component, then by field component, then by point.

Level: intermediate

-seealso: `DMField`, `DM`, `DMFieldGetNumComponents()`, `DMFieldEvaluate()`, `DMFieldEvaluateFV()`

# External Links
$(_doc_external("Dm/DMFieldEvaluateFE"))
"""
function DMFieldEvaluateFE(petsclib::PetscLibType, field::DMField, cellIS::IS, points::PetscQuadrature, datatype::PetscDataType, B::Cvoid, D::Cvoid, H::Cvoid) end

@for_petsc function DMFieldEvaluateFE(petsclib::$UnionPetscLib, field::DMField, cellIS::IS, points::PetscQuadrature, datatype::PetscDataType, B::Cvoid, D::Cvoid, H::Cvoid )

    @chk ccall(
               (:DMFieldEvaluateFE, $petsc_library),
               PetscErrorCode,
               (DMField, CIS, PetscQuadrature, PetscDataType, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
               field, cellIS, points, datatype, B, D, H,
              )


	return nothing
end 

"""
	DMFieldEvaluateFV(petsclib::PetscLibType,field::DMField, cellIS::IS, datatype::PetscDataType, B::Cvoid, D::Cvoid, H::Cvoid) 
Evaluate the mean of a field and its finite volume derivatives on a set of points.

Not Collective

Input Parameters:
- `field`    - The `DMField` object
- `cellIS`   - Index set for cells on which to evaluate the field
- `datatype` - The PetscDataType of the output arrays: either `PETSC_REAL` or `PETSC_SCALAR`.
If the field is complex and datatype is `PETSC_REAL`, the real part of the
field is returned.

Output Parameters:
- `B` - pointer to data of size c * n * sizeof(datatype), where c is the number of components in the field.
If B is not `NULL`, the values of the field are written in this array, varying first by component,
then by point.
- `D` - pointer to data of size d * c * n * sizeof(datatype).
If D is not `NULL`, the values of the field's spatial derivatives are written in this array,
varying first by the partial derivative component, then by field component, then by point.
- `H` - pointer to data of size d * d * c * n * sizeof(datatype).
If H is not `NULL`, the values of the field's second spatial derivatives are written in this array,
varying first by the second partial derivative component, then by field component, then by point.

Level: intermediate

-seealso: `DMField`, `IS`, `DMFieldGetNumComponents()`, `DMFieldEvaluate()`, `DMFieldEvaluateFE()`, `PetscDataType`

# External Links
$(_doc_external("Dm/DMFieldEvaluateFV"))
"""
function DMFieldEvaluateFV(petsclib::PetscLibType, field::DMField, cellIS::IS, datatype::PetscDataType, B::Cvoid, D::Cvoid, H::Cvoid) end

@for_petsc function DMFieldEvaluateFV(petsclib::$UnionPetscLib, field::DMField, cellIS::IS, datatype::PetscDataType, B::Cvoid, D::Cvoid, H::Cvoid )

    @chk ccall(
               (:DMFieldEvaluateFV, $petsc_library),
               PetscErrorCode,
               (DMField, CIS, PetscDataType, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
               field, cellIS, datatype, B, D, H,
              )


	return nothing
end 

"""
	minDegree::PetscInt,maxDegree::PetscInt = DMFieldGetDegree(petsclib::PetscLibType,field::DMField, cellIS::IS) 
Get the polynomial degree of a field when pulled back onto the
reference element

Not Collective

Input Parameters:
- `field`  - the `DMField` object
- `cellIS` - the index set of points over which we want know the invariance

Output Parameters:
- `minDegree` - the degree of the largest polynomial space contained in the field on each element
- `maxDegree` - the largest degree of the smallest polynomial space containing the field on any element

Level: intermediate

-seealso: `DMField`, `IS`, `DMFieldEvaluateFE()`

# External Links
$(_doc_external("Dm/DMFieldGetDegree"))
"""
function DMFieldGetDegree(petsclib::PetscLibType, field::DMField, cellIS::IS) end

@for_petsc function DMFieldGetDegree(petsclib::$UnionPetscLib, field::DMField, cellIS::IS )
	minDegree_ = Ref{$PetscInt}()
	maxDegree_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMFieldGetDegree, $petsc_library),
               PetscErrorCode,
               (DMField, CIS, Ptr{$PetscInt}, Ptr{$PetscInt}),
               field, cellIS, minDegree_, maxDegree_,
              )

	minDegree = minDegree_[]
	maxDegree = maxDegree_[]

	return minDegree,maxDegree
end 

"""
	quad::PetscQuadrature = DMFieldCreateDefaultQuadrature(petsclib::PetscLibType,field::DMField, pointIS::IS) 
Creates a quadrature sufficient to integrate the field on the selected
points via pullback onto the reference element

Not Collective

Input Parameters:
- `field`   - the `DMField` object
- `pointIS` - the index set of points over which we wish to integrate the field

Output Parameter:
- `quad` - a `PetscQuadrature` object

Level: developer

-seealso: `DMFieldCreateDefaultFaceQuadrature()`, `DMField`, `PetscQuadrature`, `IS`, `DMFieldEvaluteFE()`, `DMFieldGetDegree()`

# External Links
$(_doc_external("Dm/DMFieldCreateDefaultQuadrature"))
"""
function DMFieldCreateDefaultQuadrature(petsclib::PetscLibType, field::DMField, pointIS::IS) end

@for_petsc function DMFieldCreateDefaultQuadrature(petsclib::$UnionPetscLib, field::DMField, pointIS::IS )
	quad_ = Ref{PetscQuadrature}()

    @chk ccall(
               (:DMFieldCreateDefaultQuadrature, $petsc_library),
               PetscErrorCode,
               (DMField, CIS, Ptr{PetscQuadrature}),
               field, pointIS, quad_,
              )

	quad = quad_[]

	return quad
end 

"""
	quad::PetscQuadrature = DMFieldCreateDefaultFaceQuadrature(petsclib::PetscLibType,field::DMField, pointIS::IS) 
Creates a quadrature sufficient to integrate the field on all faces of the selected cells via pullback onto the reference element

Not Collective

Input Parameters:
- `field`   - the `DMField` object
- `pointIS` - the index set of points over which we wish to integrate the field over faces

Output Parameter:
- `quad` - a `PetscQuadrature` object

Level: developer

-seealso: `DMFieldCreateDefaultQuadrature()`, `DMField`, `PetscQuadrature`, `IS`, `DMFieldEvaluteFE()`, `DMFieldGetDegree()`

# External Links
$(_doc_external("Dm/DMFieldCreateDefaultFaceQuadrature"))
"""
function DMFieldCreateDefaultFaceQuadrature(petsclib::PetscLibType, field::DMField, pointIS::IS) end

@for_petsc function DMFieldCreateDefaultFaceQuadrature(petsclib::$UnionPetscLib, field::DMField, pointIS::IS )
	quad_ = Ref{PetscQuadrature}()

    @chk ccall(
               (:DMFieldCreateDefaultFaceQuadrature, $petsc_library),
               PetscErrorCode,
               (DMField, CIS, Ptr{PetscQuadrature}),
               field, pointIS, quad_,
              )

	quad = quad_[]

	return quad
end 

"""
	geom::PetscFEGeom = DMFieldCreateFEGeom(petsclib::PetscLibType,field::DMField, pointIS::IS, quad::PetscQuadrature, mode::PetscFEGeomMode) 
Compute and create the geometric factors of a coordinate field

Not Collective

Input Parameters:
- `field`   - the `DMField` object
- `pointIS` - the index set of points over which we wish to integrate the field
- `quad`    - the quadrature points at which to evaluate the geometric factors
- `mode`    - Type of geometry data to store

Output Parameter:
- `geom` - the geometric factors

Level: developer

-seealso: `DMField`, `PetscQuadrature`, `IS`, `PetscFEGeom`, `DMFieldEvaluateFE()`, `DMFieldCreateDefaulteQuadrature()`, `DMFieldGetDegree()`

# External Links
$(_doc_external("Dm/DMFieldCreateFEGeom"))
"""
function DMFieldCreateFEGeom(petsclib::PetscLibType, field::DMField, pointIS::IS, quad::PetscQuadrature, mode::PetscFEGeomMode) end

@for_petsc function DMFieldCreateFEGeom(petsclib::$UnionPetscLib, field::DMField, pointIS::IS, quad::PetscQuadrature, mode::PetscFEGeomMode )
	geom_ = Ref{PetscFEGeom}()

    @chk ccall(
               (:DMFieldCreateFEGeom, $petsc_library),
               PetscErrorCode,
               (DMField, CIS, PetscQuadrature, PetscFEGeomMode, PetscFEGeom),
               field, pointIS, quad, mode, geom_,
              )

	geom = geom_[]

	return geom
end 

"""
	DMFieldInitializePackage(petsclib::PetscLibType) 
Initialize `DMField` package

Logically Collective

Level: developer

-seealso: `DMFieldFinalizePackage()`

# External Links
$(_doc_external("Dm/DMFieldInitializePackage"))
"""
function DMFieldInitializePackage(petsclib::PetscLibType) end

@for_petsc function DMFieldInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:DMFieldInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	DMFieldFinalizePackage(petsclib::PetscLibType) 
Finalize `DMField` package, it is called from `PetscFinalize()`

Logically Collective

Level: developer

-seealso: `DMFieldInitializePackage()`

# External Links
$(_doc_external("Dm/DMFieldFinalizePackage"))
"""
function DMFieldFinalizePackage(petsclib::PetscLibType) end

@for_petsc function DMFieldFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:DMFieldFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	DMFieldRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Adds an implementation of the `DMField` object.

Not collective, No Fortran Support

Input Parameters:
- `sname`    - name of a new user-defined implementation
- `function` - routine to create method context

-seealso: `DMField`, `DMFieldRegisterAll()`, `DMFieldRegisterDestroy()`

# External Links
$(_doc_external("Dm/DMFieldRegister"))
"""
function DMFieldRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function DMFieldRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:DMFieldRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	cornerValues::PetscScalar,field::DMField = DMFieldCreateDA(petsclib::PetscLibType,dm::PetscDM, nc::PetscInt) 

# External Links
$(_doc_external("Dm/DMFieldCreateDA"))
"""
function DMFieldCreateDA(petsclib::PetscLibType, dm::PetscDM, nc::PetscInt) end

@for_petsc function DMFieldCreateDA(petsclib::$UnionPetscLib, dm::PetscDM, nc::$PetscInt )
	cornerValues_ = Ref{$PetscScalar}()
	field_ = Ref{DMField}()

    @chk ccall(
               (:DMFieldCreateDA, $petsc_library),
               PetscErrorCode,
               (CDM, $PetscInt, Ptr{$PetscScalar}, Ptr{DMField}),
               dm, nc, cornerValues_, field_,
              )

	cornerValues = cornerValues_[]
	field = field_[]

	return cornerValues,field
end 

"""
	field::DMField = DMFieldCreateDSWithDG(petsclib::PetscLibType,dm::PetscDM, dmDG::PetscDM, fieldNum::PetscInt, vec::PetscVec, vecDG::PetscVec) 

# External Links
$(_doc_external("Dm/DMFieldCreateDSWithDG"))
"""
function DMFieldCreateDSWithDG(petsclib::PetscLibType, dm::PetscDM, dmDG::PetscDM, fieldNum::PetscInt, vec::PetscVec, vecDG::PetscVec) end

@for_petsc function DMFieldCreateDSWithDG(petsclib::$UnionPetscLib, dm::PetscDM, dmDG::PetscDM, fieldNum::$PetscInt, vec::PetscVec, vecDG::PetscVec )
	field_ = Ref{DMField}()

    @chk ccall(
               (:DMFieldCreateDSWithDG, $petsc_library),
               PetscErrorCode,
               (CDM, CDM, $PetscInt, CVec, CVec, Ptr{DMField}),
               dm, dmDG, fieldNum, vec, vecDG, field_,
              )

	field = field_[]

	return field
end 

"""
	field::DMField = DMFieldCreateDS(petsclib::PetscLibType,dm::PetscDM, fieldNum::PetscInt, vec::PetscVec) 

# External Links
$(_doc_external("Dm/DMFieldCreateDS"))
"""
function DMFieldCreateDS(petsclib::PetscLibType, dm::PetscDM, fieldNum::PetscInt, vec::PetscVec) end

@for_petsc function DMFieldCreateDS(petsclib::$UnionPetscLib, dm::PetscDM, fieldNum::$PetscInt, vec::PetscVec )
	field_ = Ref{DMField}()

    @chk ccall(
               (:DMFieldCreateDS, $petsc_library),
               PetscErrorCode,
               (CDM, $PetscInt, CVec, Ptr{DMField}),
               dm, fieldNum, vec, field_,
              )

	field = field_[]

	return field
end 

"""
	DMFieldShellGetContext(petsclib::PetscLibType,field::DMField, ctx::Cvoid) 

# External Links
$(_doc_external("Dm/DMFieldShellGetContext"))
"""
function DMFieldShellGetContext(petsclib::PetscLibType, field::DMField, ctx::Cvoid) end

@for_petsc function DMFieldShellGetContext(petsclib::$UnionPetscLib, field::DMField, ctx::Cvoid )

    @chk ccall(
               (:DMFieldShellGetContext, $petsc_library),
               PetscErrorCode,
               (DMField, Ptr{Cvoid}),
               field, ctx,
              )


	return nothing
end 

"""
	DMFieldShellEvaluateFEDefault(petsclib::PetscLibType,field::DMField, pointIS::IS, quad::PetscQuadrature, type::PetscDataType, B::Cvoid, D::Cvoid, H::Cvoid) 

# External Links
$(_doc_external("Dm/DMFieldShellEvaluateFEDefault"))
"""
function DMFieldShellEvaluateFEDefault(petsclib::PetscLibType, field::DMField, pointIS::IS, quad::PetscQuadrature, type::PetscDataType, B::Cvoid, D::Cvoid, H::Cvoid) end

@for_petsc function DMFieldShellEvaluateFEDefault(petsclib::$UnionPetscLib, field::DMField, pointIS::IS, quad::PetscQuadrature, type::PetscDataType, B::Cvoid, D::Cvoid, H::Cvoid )

    @chk ccall(
               (:DMFieldShellEvaluateFEDefault, $petsc_library),
               PetscErrorCode,
               (DMField, CIS, PetscQuadrature, PetscDataType, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
               field, pointIS, quad, type, B, D, H,
              )


	return nothing
end 

"""
	DMFieldShellEvaluateFVDefault(petsclib::PetscLibType,field::DMField, pointIS::IS, type::PetscDataType, B::Cvoid, D::Cvoid, H::Cvoid) 

# External Links
$(_doc_external("Dm/DMFieldShellEvaluateFVDefault"))
"""
function DMFieldShellEvaluateFVDefault(petsclib::PetscLibType, field::DMField, pointIS::IS, type::PetscDataType, B::Cvoid, D::Cvoid, H::Cvoid) end

@for_petsc function DMFieldShellEvaluateFVDefault(petsclib::$UnionPetscLib, field::DMField, pointIS::IS, type::PetscDataType, B::Cvoid, D::Cvoid, H::Cvoid )

    @chk ccall(
               (:DMFieldShellEvaluateFVDefault, $petsc_library),
               PetscErrorCode,
               (DMField, CIS, PetscDataType, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
               field, pointIS, type, B, D, H,
              )


	return nothing
end 

"""
	DMFieldShellSetDestroy(petsclib::PetscLibType,field::DMField, destroy::external) 

# External Links
$(_doc_external("Dm/DMFieldShellSetDestroy"))
"""
function DMFieldShellSetDestroy(petsclib::PetscLibType, field::DMField, destroy::external) end

@for_petsc function DMFieldShellSetDestroy(petsclib::$UnionPetscLib, field::DMField, destroy::external )

    @chk ccall(
               (:DMFieldShellSetDestroy, $petsc_library),
               PetscErrorCode,
               (DMField, external),
               field, destroy,
              )


	return nothing
end 

"""
	DMFieldShellSetEvaluate(petsclib::PetscLibType,field::DMField, evaluate::external) 

# External Links
$(_doc_external("Dm/DMFieldShellSetEvaluate"))
"""
function DMFieldShellSetEvaluate(petsclib::PetscLibType, field::DMField, evaluate::external) end

@for_petsc function DMFieldShellSetEvaluate(petsclib::$UnionPetscLib, field::DMField, evaluate::external )

    @chk ccall(
               (:DMFieldShellSetEvaluate, $petsc_library),
               PetscErrorCode,
               (DMField, external),
               field, evaluate,
              )


	return nothing
end 

"""
	DMFieldShellSetEvaluateFE(petsclib::PetscLibType,field::DMField, evaluateFE::external) 

# External Links
$(_doc_external("Dm/DMFieldShellSetEvaluateFE"))
"""
function DMFieldShellSetEvaluateFE(petsclib::PetscLibType, field::DMField, evaluateFE::external) end

@for_petsc function DMFieldShellSetEvaluateFE(petsclib::$UnionPetscLib, field::DMField, evaluateFE::external )

    @chk ccall(
               (:DMFieldShellSetEvaluateFE, $petsc_library),
               PetscErrorCode,
               (DMField, external),
               field, evaluateFE,
              )


	return nothing
end 

"""
	DMFieldShellSetEvaluateFV(petsclib::PetscLibType,field::DMField, evaluateFV::external) 

# External Links
$(_doc_external("Dm/DMFieldShellSetEvaluateFV"))
"""
function DMFieldShellSetEvaluateFV(petsclib::PetscLibType, field::DMField, evaluateFV::external) end

@for_petsc function DMFieldShellSetEvaluateFV(petsclib::$UnionPetscLib, field::DMField, evaluateFV::external )

    @chk ccall(
               (:DMFieldShellSetEvaluateFV, $petsc_library),
               PetscErrorCode,
               (DMField, external),
               field, evaluateFV,
              )


	return nothing
end 

"""
	DMFieldShellSetGetDegree(petsclib::PetscLibType,field::DMField, getDegree::external) 

# External Links
$(_doc_external("Dm/DMFieldShellSetGetDegree"))
"""
function DMFieldShellSetGetDegree(petsclib::PetscLibType, field::DMField, getDegree::external) end

@for_petsc function DMFieldShellSetGetDegree(petsclib::$UnionPetscLib, field::DMField, getDegree::external )

    @chk ccall(
               (:DMFieldShellSetGetDegree, $petsc_library),
               PetscErrorCode,
               (DMField, external),
               field, getDegree,
              )


	return nothing
end 

"""
	DMFieldShellSetCreateDefaultQuadrature(petsclib::PetscLibType,field::DMField, createDefaultQuadrature::external) 

# External Links
$(_doc_external("Dm/DMFieldShellSetCreateDefaultQuadrature"))
"""
function DMFieldShellSetCreateDefaultQuadrature(petsclib::PetscLibType, field::DMField, createDefaultQuadrature::external) end

@for_petsc function DMFieldShellSetCreateDefaultQuadrature(petsclib::$UnionPetscLib, field::DMField, createDefaultQuadrature::external )

    @chk ccall(
               (:DMFieldShellSetCreateDefaultQuadrature, $petsc_library),
               PetscErrorCode,
               (DMField, external),
               field, createDefaultQuadrature,
              )


	return nothing
end 

"""
	ctx::Cvoid,field::DMField = DMFieldCreateShell(petsclib::PetscLibType,dm::PetscDM, numComponents::PetscInt, continuity::DMFieldContinuity) 

# External Links
$(_doc_external("Dm/DMFieldCreateShell"))
"""
function DMFieldCreateShell(petsclib::PetscLibType, dm::PetscDM, numComponents::PetscInt, continuity::DMFieldContinuity) end

@for_petsc function DMFieldCreateShell(petsclib::$UnionPetscLib, dm::PetscDM, numComponents::$PetscInt, continuity::DMFieldContinuity )
	ctx_ = Ref{Cvoid}()
	field_ = Ref{DMField}()

    @chk ccall(
               (:DMFieldCreateShell, $petsc_library),
               PetscErrorCode,
               (CDM, $PetscInt, DMFieldContinuity, Ptr{Cvoid}, Ptr{DMField}),
               dm, numComponents, continuity, ctx_, field_,
              )

	ctx = ctx_[]
	field = field_[]

	return ctx,field
end 

"""
	DMPlexTransformRegisterAll(petsclib::PetscLibType) 
Registers all of the transform components in the `DM` package.

Not Collective

Level: advanced

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `DMPlexTransformType`, `DMRegisterAll()`, `DMPlexTransformRegisterDestroy()`

# External Links
$(_doc_external("Dm/DMPlexTransformRegisterAll"))
"""
function DMPlexTransformRegisterAll(petsclib::PetscLibType) end

@for_petsc function DMPlexTransformRegisterAll(petsclib::$UnionPetscLib)

    @chk ccall(
               (:DMPlexTransformRegisterAll, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	DMPlexTransformRegisterDestroy(petsclib::PetscLibType) 
This function destroys the registered `DMPlexTransformType`. It is called from `PetscFinalize()`.

Not collective

Level: developer

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `DMRegisterAll()`, `DMPlexTransformType`, `PetscInitialize()`

# External Links
$(_doc_external("Dm/DMPlexTransformRegisterDestroy"))
"""
function DMPlexTransformRegisterDestroy(petsclib::PetscLibType) end

@for_petsc function DMPlexTransformRegisterDestroy(petsclib::$UnionPetscLib)

    @chk ccall(
               (:DMPlexTransformRegisterDestroy, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	tr::DMPlexTransform = DMPlexTransformCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Creates an empty transform object. The type can then be set with `DMPlexTransformSetType()`.

Collective

Input Parameter:
- `comm` - The communicator for the transform object

Output Parameter:
- `tr` - The transform object

Level: beginner

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformType`, `DMPlexTransformSetType()`, `DMPLEXREFINEREGULAR`, `DMPLEXTRANSFORMFILTER`

# External Links
$(_doc_external("Dm/DMPlexTransformCreate"))
"""
function DMPlexTransformCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function DMPlexTransformCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	tr_ = Ref{DMPlexTransform}()

    @chk ccall(
               (:DMPlexTransformCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{DMPlexTransform}),
               comm, tr_,
              )

	tr = tr_[]

	return tr
end 

"""
	DMPlexTransformSetType(petsclib::PetscLibType,tr::DMPlexTransform, method::DMPlexTransformType) 
Sets the particular implementation for a transform.

Collective

Input Parameters:
- `tr`     - The transform
- `method` - The name of the transform type

Options Database Key:
- `-dm_plex_transform_type <type>` - Sets the transform type; see `DMPlexTransformType`

Level: intermediate

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformType`, `DMPlexTransformGetType()`, `DMPlexTransformCreate()`

# External Links
$(_doc_external("Dm/DMPlexTransformSetType"))
"""
function DMPlexTransformSetType(petsclib::PetscLibType, tr::DMPlexTransform, method::DMPlexTransformType) end

@for_petsc function DMPlexTransformSetType(petsclib::$UnionPetscLib, tr::DMPlexTransform, method::DMPlexTransformType )

    @chk ccall(
               (:DMPlexTransformSetType, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, DMPlexTransformType),
               tr, method,
              )


	return nothing
end 

"""
	type::DMPlexTransformType = DMPlexTransformGetType(petsclib::PetscLibType,tr::DMPlexTransform) 
Gets the type name (as a string) from the transform.

Not Collective

Input Parameter:
- `tr` - The `DMPlexTransform`

Output Parameter:
- `type` - The `DMPlexTransformType` name

Level: intermediate

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformType`, `DMPlexTransformSetType()`, `DMPlexTransformCreate()`

# External Links
$(_doc_external("Dm/DMPlexTransformGetType"))
"""
function DMPlexTransformGetType(petsclib::PetscLibType, tr::DMPlexTransform) end

@for_petsc function DMPlexTransformGetType(petsclib::$UnionPetscLib, tr::DMPlexTransform )
	type_ = Ref{DMPlexTransformType}()

    @chk ccall(
               (:DMPlexTransformGetType, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{DMPlexTransformType}),
               tr, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	DMPlexTransformView(petsclib::PetscLibType,tr::DMPlexTransform, v::PetscViewer) 
Views a `DMPlexTransform`

Collective

Input Parameters:
- `tr` - the `DMPlexTransform` object to view
- `v`  - the viewer

Level: beginner

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformType`, `PetscViewer`, `DMPlexTransformDestroy()`, `DMPlexTransformCreate()`

# External Links
$(_doc_external("Dm/DMPlexTransformView"))
"""
function DMPlexTransformView(petsclib::PetscLibType, tr::DMPlexTransform, v::PetscViewer) end

@for_petsc function DMPlexTransformView(petsclib::$UnionPetscLib, tr::DMPlexTransform, v::PetscViewer )

    @chk ccall(
               (:DMPlexTransformView, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, PetscViewer),
               tr, v,
              )


	return nothing
end 

"""
	DMPlexTransformSetFromOptions(petsclib::PetscLibType,tr::DMPlexTransform) 
Sets parameters in a transform from values in the options database

Collective

Input Parameter:
- `tr` - the `DMPlexTransform` object to set options for

Options Database Keys:
- `-dm_plex_transform_type`                      - Set the transform type, e.g. refine_regular
- `-dm_plex_transform_label_match_strata`        - Only label points of the same stratum as the producing point
- `-dm_plex_transform_label_replica_inc <inc>`   - Increment for the label value to be multiplied by the replica number, so that the new label value is oldValue + r * inc
- `-dm_plex_transform_active <name>`             - Name for active mesh label
- `-dm_plex_transform_active_values <v0,v1,...>` - Values in the active label

Level: intermediate

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformView()`, `DMPlexTransformCreate()`

# External Links
$(_doc_external("Dm/DMPlexTransformSetFromOptions"))
"""
function DMPlexTransformSetFromOptions(petsclib::PetscLibType, tr::DMPlexTransform) end

@for_petsc function DMPlexTransformSetFromOptions(petsclib::$UnionPetscLib, tr::DMPlexTransform )

    @chk ccall(
               (:DMPlexTransformSetFromOptions, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform,),
               tr,
              )


	return nothing
end 

"""
	DMPlexTransformDestroy(petsclib::PetscLibType,tr::DMPlexTransform) 
Destroys a `DMPlexTransform`

Collective

Input Parameter:
- `tr` - the transform object to destroy

Level: beginner

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformView()`, `DMPlexTransformCreate()`

# External Links
$(_doc_external("Dm/DMPlexTransformDestroy"))
"""
function DMPlexTransformDestroy(petsclib::PetscLibType, tr::DMPlexTransform) end

@for_petsc function DMPlexTransformDestroy(petsclib::$UnionPetscLib, tr::DMPlexTransform )

    @chk ccall(
               (:DMPlexTransformDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{DMPlexTransform},),
               tr,
              )


	return nothing
end 

"""
	DMPlexTransformSetUp(petsclib::PetscLibType,tr::DMPlexTransform) 
Create the tables that drive the transform

Input Parameter:
- `tr` - The `DMPlexTransform` object

Level: intermediate

-seealso: [](plex_transform_table), [](ch_unstructured), `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformApply()`, `DMPlexTransformCreate()`

# External Links
$(_doc_external("Dm/DMPlexTransformSetUp"))
"""
function DMPlexTransformSetUp(petsclib::PetscLibType, tr::DMPlexTransform) end

@for_petsc function DMPlexTransformSetUp(petsclib::$UnionPetscLib, tr::DMPlexTransform )

    @chk ccall(
               (:DMPlexTransformSetUp, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform,),
               tr,
              )


	return nothing
end 

"""
	DMPlexTransformGetDM(petsclib::PetscLibType,tr::DMPlexTransform, dm::PetscDM) 
Get the base `DM` for the transform

Input Parameter:
- `tr` - The `DMPlexTransform` object

Output Parameter:
- `dm` - The original `DM` which will be transformed

Level: intermediate

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformSetDM()`, `DMPlexTransformApply()`, `DMPlexTransformCreate()`

# External Links
$(_doc_external("Dm/DMPlexTransformGetDM"))
"""
function DMPlexTransformGetDM(petsclib::PetscLibType, tr::DMPlexTransform, dm::PetscDM) end

@for_petsc function DMPlexTransformGetDM(petsclib::$UnionPetscLib, tr::DMPlexTransform, dm::PetscDM )
	dm_ = Ref(dm.ptr)

    @chk ccall(
               (:DMPlexTransformGetDM, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{CDM}),
               tr, dm_,
              )

	dm.ptr = C_NULL

	return nothing
end 

"""
	DMPlexTransformSetDM(petsclib::PetscLibType,tr::DMPlexTransform, dm::PetscDM) 
Set the base `DM` for the transform

Input Parameters:
- `tr` - The `DMPlexTransform` object
- `dm` - The original `DM` which will be transformed

Level: intermediate

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformGetDM()`, `DMPlexTransformApply()`, `DMPlexTransformCreate()`

# External Links
$(_doc_external("Dm/DMPlexTransformSetDM"))
"""
function DMPlexTransformSetDM(petsclib::PetscLibType, tr::DMPlexTransform, dm::PetscDM) end

@for_petsc function DMPlexTransformSetDM(petsclib::$UnionPetscLib, tr::DMPlexTransform, dm::PetscDM )

    @chk ccall(
               (:DMPlexTransformSetDM, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, CDM),
               tr, dm,
              )


	return nothing
end 

"""
	DMPlexTransformGetActive(petsclib::PetscLibType,tr::DMPlexTransform, active::DMLabel) 
Get the `DMLabel` marking the active points for the transform

Input Parameter:
- `tr` - The `DMPlexTransform` object

Output Parameter:
- `active` - The `DMLabel` indicating which points will be transformed

Level: intermediate

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformSetActive()`, `DMPlexTransformApply()`, `DMPlexTransformCreate()`

# External Links
$(_doc_external("Dm/DMPlexTransformGetActive"))
"""
function DMPlexTransformGetActive(petsclib::PetscLibType, tr::DMPlexTransform, active::DMLabel) end

@for_petsc function DMPlexTransformGetActive(petsclib::$UnionPetscLib, tr::DMPlexTransform, active::DMLabel )

    @chk ccall(
               (:DMPlexTransformGetActive, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{DMLabel}),
               tr, active,
              )


	return nothing
end 

"""
	DMPlexTransformSetActive(petsclib::PetscLibType,tr::DMPlexTransform, active::DMLabel) 
Set the `DMLabel` marking the active points for the transform

Input Parameters:
- `tr`     - The `DMPlexTransform` object
- `active` - The `DMLabel` indicating which points will be transformed

Level: intermediate

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformGetActive()`, `DMPlexTransformApply()`, `DMPlexTransformCreate()`

# External Links
$(_doc_external("Dm/DMPlexTransformSetActive"))
"""
function DMPlexTransformSetActive(petsclib::PetscLibType, tr::DMPlexTransform, active::DMLabel) end

@for_petsc function DMPlexTransformSetActive(petsclib::$UnionPetscLib, tr::DMPlexTransform, active::DMLabel )

    @chk ccall(
               (:DMPlexTransformSetActive, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, DMLabel),
               tr, active,
              )


	return nothing
end 

"""
	trType::DMLabel = DMPlexTransformGetTransformTypes(petsclib::PetscLibType,tr::DMPlexTransform) 
Get the `DMLabel` marking the transform type of each point for the transform

Input Parameter:
- `tr` - The `DMPlexTransform` object

Output Parameter:
- `trType` - The `DMLabel` indicating the transform type for each point

Level: intermediate

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexSetTransformType()`, `DMPlexTransformGetActive()`, `DMPlexTransformApply()`, `DMPlexTransformCreate()`

# External Links
$(_doc_external("Dm/DMPlexTransformGetTransformTypes"))
"""
function DMPlexTransformGetTransformTypes(petsclib::PetscLibType, tr::DMPlexTransform) end

@for_petsc function DMPlexTransformGetTransformTypes(petsclib::$UnionPetscLib, tr::DMPlexTransform )
	trType_ = Ref{DMLabel}()

    @chk ccall(
               (:DMPlexTransformGetTransformTypes, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{DMLabel}),
               tr, trType_,
              )

	trType = trType_[]

	return trType
end 

"""
	DMPlexTransformSetTransformTypes(petsclib::PetscLibType,tr::DMPlexTransform, trType::DMLabel) 
Set the `DMLabel` marking the transform type of each point for the transform

Input Parameters:
- `tr`     - The `DMPlexTransform` object
- `trType` - The original `DM` which will be transformed

Level: intermediate

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformGetTransformTypes()`, `DMPlexTransformGetActive())`, `DMPlexTransformApply()`, `DMPlexTransformCreate()`

# External Links
$(_doc_external("Dm/DMPlexTransformSetTransformTypes"))
"""
function DMPlexTransformSetTransformTypes(petsclib::PetscLibType, tr::DMPlexTransform, trType::DMLabel) end

@for_petsc function DMPlexTransformSetTransformTypes(petsclib::$UnionPetscLib, tr::DMPlexTransform, trType::DMLabel )

    @chk ccall(
               (:DMPlexTransformSetTransformTypes, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, DMLabel),
               tr, trType,
              )


	return nothing
end 

"""
	DMPlexTransformSetDimensions(petsclib::PetscLibType,tr::DMPlexTransform, dm::PetscDM, tdm::PetscDM) 
Set the dimensions for the transformed `DM`

Input Parameters:
- `tr` - The `DMPlexTransform` object
- `dm` - The original `DM`

Output Parameter:
- `tdm` - The transformed `DM`

Level: advanced

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformApply()`, `DMPlexTransformCreate()`

# External Links
$(_doc_external("Dm/DMPlexTransformSetDimensions"))
"""
function DMPlexTransformSetDimensions(petsclib::PetscLibType, tr::DMPlexTransform, dm::PetscDM, tdm::PetscDM) end

@for_petsc function DMPlexTransformSetDimensions(petsclib::$UnionPetscLib, tr::DMPlexTransform, dm::PetscDM, tdm::PetscDM )

    @chk ccall(
               (:DMPlexTransformSetDimensions, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, CDM, CDM),
               tr, dm, tdm,
              )


	return nothing
end 

"""
	pStart::PetscInt,pEnd::PetscInt = DMPlexTransformGetChart(petsclib::PetscLibType,tr::DMPlexTransform) 

# External Links
$(_doc_external("Dm/DMPlexTransformGetChart"))
"""
function DMPlexTransformGetChart(petsclib::PetscLibType, tr::DMPlexTransform) end

@for_petsc function DMPlexTransformGetChart(petsclib::$UnionPetscLib, tr::DMPlexTransform )
	pStart_ = Ref{$PetscInt}()
	pEnd_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMPlexTransformGetChart, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{$PetscInt}, Ptr{$PetscInt}),
               tr, pStart_, pEnd_,
              )

	pStart = pStart_[]
	pEnd = pEnd_[]

	return pStart,pEnd
end 

"""
	celltype::DMPolytopeType = DMPlexTransformGetCellType(petsclib::PetscLibType,tr::DMPlexTransform, cell::PetscInt) 

# External Links
$(_doc_external("Dm/DMPlexTransformGetCellType"))
"""
function DMPlexTransformGetCellType(petsclib::PetscLibType, tr::DMPlexTransform, cell::PetscInt) end

@for_petsc function DMPlexTransformGetCellType(petsclib::$UnionPetscLib, tr::DMPlexTransform, cell::$PetscInt )
	celltype_ = Ref{DMPolytopeType}()

    @chk ccall(
               (:DMPlexTransformGetCellType, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, $PetscInt, Ptr{DMPolytopeType}),
               tr, cell, celltype_,
              )

	celltype = unsafe_string(celltype_[])

	return celltype
end 

"""
	start::PetscInt,end_::PetscInt = DMPlexTransformGetCellTypeStratum(petsclib::PetscLibType,tr::DMPlexTransform, celltype::DMPolytopeType) 

# External Links
$(_doc_external("Dm/DMPlexTransformGetCellTypeStratum"))
"""
function DMPlexTransformGetCellTypeStratum(petsclib::PetscLibType, tr::DMPlexTransform, celltype::DMPolytopeType) end

@for_petsc function DMPlexTransformGetCellTypeStratum(petsclib::$UnionPetscLib, tr::DMPlexTransform, celltype::DMPolytopeType )
	start_ = Ref{$PetscInt}()
	end__ = Ref{$PetscInt}()

    @chk ccall(
               (:DMPlexTransformGetCellTypeStratum, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, DMPolytopeType, Ptr{$PetscInt}, Ptr{$PetscInt}),
               tr, celltype, start_, end__,
              )

	start = start_[]
	end_ = end__[]

	return start,end_
end 

"""
	depth::PetscInt = DMPlexTransformGetDepth(petsclib::PetscLibType,tr::DMPlexTransform) 

# External Links
$(_doc_external("Dm/DMPlexTransformGetDepth"))
"""
function DMPlexTransformGetDepth(petsclib::PetscLibType, tr::DMPlexTransform) end

@for_petsc function DMPlexTransformGetDepth(petsclib::$UnionPetscLib, tr::DMPlexTransform )
	depth_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMPlexTransformGetDepth, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{$PetscInt}),
               tr, depth_,
              )

	depth = depth_[]

	return depth
end 

"""
	start::PetscInt,end_::PetscInt = DMPlexTransformGetDepthStratum(petsclib::PetscLibType,tr::DMPlexTransform, depth::PetscInt) 

# External Links
$(_doc_external("Dm/DMPlexTransformGetDepthStratum"))
"""
function DMPlexTransformGetDepthStratum(petsclib::PetscLibType, tr::DMPlexTransform, depth::PetscInt) end

@for_petsc function DMPlexTransformGetDepthStratum(petsclib::$UnionPetscLib, tr::DMPlexTransform, depth::$PetscInt )
	start_ = Ref{$PetscInt}()
	end__ = Ref{$PetscInt}()

    @chk ccall(
               (:DMPlexTransformGetDepthStratum, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
               tr, depth, start_, end__,
              )

	start = start_[]
	end_ = end__[]

	return start,end_
end 

"""
	match::PetscBool = DMPlexTransformGetMatchStrata(petsclib::PetscLibType,tr::DMPlexTransform) 
Get the flag which determines what points get added to the transformed labels

Not Collective

Input Parameter:
- `tr` - The `DMPlexTransform`

Output Parameter:
- `match` - If `PETSC_TRUE`, only add produced points at the same stratum as the original point to new labels

Level: intermediate

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformSetMatchStrata()`, `DMPlexGetPointDepth()`

# External Links
$(_doc_external("Dm/DMPlexTransformGetMatchStrata"))
"""
function DMPlexTransformGetMatchStrata(petsclib::PetscLibType, tr::DMPlexTransform) end

@for_petsc function DMPlexTransformGetMatchStrata(petsclib::$UnionPetscLib, tr::DMPlexTransform )
	match_ = Ref{PetscBool}()

    @chk ccall(
               (:DMPlexTransformGetMatchStrata, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{PetscBool}),
               tr, match_,
              )

	match = match_[]

	return match
end 

"""
	DMPlexTransformSetMatchStrata(petsclib::PetscLibType,tr::DMPlexTransform, match::PetscBool) 
Set the flag which determines what points get added to the transformed labels

Not Collective

Input Parameters:
- `tr`    - The `DMPlexTransform`
- `match` - If `PETSC_TRUE`, only add produced points at the same stratum as the original point to new labels

Level: intermediate

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformGetMatchStrata()`, `DMPlexGetPointDepth()`

# External Links
$(_doc_external("Dm/DMPlexTransformSetMatchStrata"))
"""
function DMPlexTransformSetMatchStrata(petsclib::PetscLibType, tr::DMPlexTransform, match::PetscBool) end

@for_petsc function DMPlexTransformSetMatchStrata(petsclib::$UnionPetscLib, tr::DMPlexTransform, match::PetscBool )

    @chk ccall(
               (:DMPlexTransformSetMatchStrata, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, PetscBool),
               tr, match,
              )


	return nothing
end 

"""
	pNew::PetscInt = DMPlexTransformGetTargetPoint(petsclib::PetscLibType,tr::DMPlexTransform, ct::DMPolytopeType, ctNew::DMPolytopeType, p::PetscInt, r::PetscInt) 
Get the number of a point in the transformed mesh based on information from the original mesh.

Not Collective

Input Parameters:
- `tr`    - The `DMPlexTransform`
- `ct`    - The type of the original point which produces the new point
- `ctNew` - The type of the new point
- `p`     - The original point which produces the new point
- `r`     - The replica number of the new point, meaning it is the rth point of type `ctNew` produced from `p`

Output Parameter:
- `pNew` - The new point number

Level: developer

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `DMPlexTransform`, `DMPolytopeType`, `DMPlexTransformGetSourcePoint()`, `DMPlexTransformCellTransform()`

# External Links
$(_doc_external("Dm/DMPlexTransformGetTargetPoint"))
"""
function DMPlexTransformGetTargetPoint(petsclib::PetscLibType, tr::DMPlexTransform, ct::DMPolytopeType, ctNew::DMPolytopeType, p::PetscInt, r::PetscInt) end

@for_petsc function DMPlexTransformGetTargetPoint(petsclib::$UnionPetscLib, tr::DMPlexTransform, ct::DMPolytopeType, ctNew::DMPolytopeType, p::$PetscInt, r::$PetscInt )
	pNew_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMPlexTransformGetTargetPoint, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, DMPolytopeType, DMPolytopeType, $PetscInt, $PetscInt, Ptr{$PetscInt}),
               tr, ct, ctNew, p, r, pNew_,
              )

	pNew = pNew_[]

	return pNew
end 

"""
	p::PetscInt,r::PetscInt = DMPlexTransformGetSourcePoint(petsclib::PetscLibType,tr::DMPlexTransform, pNew::PetscInt, ct::DMPolytopeType, ctNew::DMPolytopeType) 
Get the number of a point in the original mesh based on information from the transformed mesh.

Not Collective

Input Parameters:
- `tr`   - The `DMPlexTransform`
- `pNew` - The new point number

Output Parameters:
- `ct`    - The type of the original point which produces the new point
- `ctNew` - The type of the new point
- `p`     - The original point which produces the new point
- `r`     - The replica number of the new point, meaning it is the rth point of type ctNew produced from p

Level: developer

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `DMPlexTransform`, `DMPolytopeType`, `DMPlexTransformGetTargetPoint()`, `DMPlexTransformCellTransform()`

# External Links
$(_doc_external("Dm/DMPlexTransformGetSourcePoint"))
"""
function DMPlexTransformGetSourcePoint(petsclib::PetscLibType, tr::DMPlexTransform, pNew::PetscInt, ct::DMPolytopeType, ctNew::DMPolytopeType) end

@for_petsc function DMPlexTransformGetSourcePoint(petsclib::$UnionPetscLib, tr::DMPlexTransform, pNew::$PetscInt, ct::DMPolytopeType, ctNew::DMPolytopeType )
	p_ = Ref{$PetscInt}()
	r_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMPlexTransformGetSourcePoint, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, $PetscInt, Ptr{DMPolytopeType}, Ptr{DMPolytopeType}, Ptr{$PetscInt}, Ptr{$PetscInt}),
               tr, pNew, ct, ctNew, p_, r_,
              )

	p = p_[]
	r = r_[]

	return p,r
end 

"""
	rt::PetscInt,Nt::PetscInt,size::Vector{PetscInt},cone::Vector{PetscInt},ornt::Vector{PetscInt} = DMPlexTransformCellTransform(petsclib::PetscLibType,tr::DMPlexTransform, source::DMPolytopeType, p::PetscInt, target::Vector{DMPolytopeType}) 
Describes the transform of a given source cell into a set of other target cells. These produced cells become the new mesh.

Input Parameters:
- `tr`     - The `DMPlexTransform` object
- `source` - The source cell type
- `p`      - The source point, which can also determine the refine type

Output Parameters:
- `rt`     - The refine type for this point
- `Nt`     - The number of types produced by this point
- `target` - An array of length `Nt` giving the types produced
- `size`   - An array of length `Nt` giving the number of cells of each type produced
- `cone`   - An array of length `Nt`*size[t]*coneSize[t] giving the cell type for each point in the cone of each produced point
- `ornt`   - An array of length `Nt`*size[t]*coneSize[t] giving the orientation for each point in the cone of each produced point

Level: advanced

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `DMPlexTransform`, `DMPolytopeType`, `DMPlexTransformApply()`, `DMPlexTransformCreate()`

# External Links
$(_doc_external("Dm/DMPlexTransformCellTransform"))
"""
function DMPlexTransformCellTransform(petsclib::PetscLibType, tr::DMPlexTransform, source::DMPolytopeType, p::PetscInt, target::Vector{DMPolytopeType}) end

@for_petsc function DMPlexTransformCellTransform(petsclib::$UnionPetscLib, tr::DMPlexTransform, source::DMPolytopeType, p::$PetscInt, target::Vector{DMPolytopeType} )
	rt_ = Ref{$PetscInt}()
	Nt_ = Ref{$PetscInt}()
	target_ = Ref(pointer(target))
	size_ = Ref{Ptr{$PetscInt}}()
	cone_ = Ref{Ptr{$PetscInt}}()
	ornt_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:DMPlexTransformCellTransform, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, DMPolytopeType, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{DMPolytopeType}}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}),
               tr, source, p, rt_, Nt_, target_, size_, cone_, ornt_,
              )

	rt = rt_[]
	Nt = Nt_[]
	size = unsafe_wrap(Array, size_[], VecGetLocalSize(petsclib, x); own = false)
	cone = unsafe_wrap(Array, cone_[], VecGetLocalSize(petsclib, x); own = false)
	ornt = unsafe_wrap(Array, ornt_[], VecGetLocalSize(petsclib, x); own = false)

	return rt,Nt,size,cone,ornt
end 

"""
	rnew::PetscInt,onew::PetscInt = DMPlexTransformGetSubcellOrientationIdentity(petsclib::PetscLibType,tr::DMPlexTransform, sct::DMPolytopeType, sp::PetscInt, so::PetscInt, tct::DMPolytopeType, r::PetscInt, o::PetscInt) 

# External Links
$(_doc_external("Dm/DMPlexTransformGetSubcellOrientationIdentity"))
"""
function DMPlexTransformGetSubcellOrientationIdentity(petsclib::PetscLibType, tr::DMPlexTransform, sct::DMPolytopeType, sp::PetscInt, so::PetscInt, tct::DMPolytopeType, r::PetscInt, o::PetscInt) end

@for_petsc function DMPlexTransformGetSubcellOrientationIdentity(petsclib::$UnionPetscLib, tr::DMPlexTransform, sct::DMPolytopeType, sp::$PetscInt, so::$PetscInt, tct::DMPolytopeType, r::$PetscInt, o::$PetscInt )
	rnew_ = Ref{$PetscInt}()
	onew_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMPlexTransformGetSubcellOrientationIdentity, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, DMPolytopeType, $PetscInt, $PetscInt, DMPolytopeType, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
               tr, sct, sp, so, tct, r, o, rnew_, onew_,
              )

	rnew = rnew_[]
	onew = onew_[]

	return rnew,onew
end 

"""
	rt::PetscInt,Nt::PetscInt,size::Vector{PetscInt},cone::Vector{PetscInt},ornt::Vector{PetscInt} = DMPlexTransformCellTransformIdentity(petsclib::PetscLibType,tr::DMPlexTransform, source::DMPolytopeType, p::PetscInt, target::Vector{DMPolytopeType}) 

# External Links
$(_doc_external("Dm/DMPlexTransformCellTransformIdentity"))
"""
function DMPlexTransformCellTransformIdentity(petsclib::PetscLibType, tr::DMPlexTransform, source::DMPolytopeType, p::PetscInt, target::Vector{DMPolytopeType}) end

@for_petsc function DMPlexTransformCellTransformIdentity(petsclib::$UnionPetscLib, tr::DMPlexTransform, source::DMPolytopeType, p::$PetscInt, target::Vector{DMPolytopeType} )
	rt_ = Ref{$PetscInt}()
	Nt_ = Ref{$PetscInt}()
	target_ = Ref(pointer(target))
	size_ = Ref{Ptr{$PetscInt}}()
	cone_ = Ref{Ptr{$PetscInt}}()
	ornt_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:DMPlexTransformCellTransformIdentity, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, DMPolytopeType, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{DMPolytopeType}}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}),
               tr, source, p, rt_, Nt_, target_, size_, cone_, ornt_,
              )

	rt = rt_[]
	Nt = Nt_[]
	size = unsafe_wrap(Array, size_[], VecGetLocalSize(petsclib, x); own = false)
	cone = unsafe_wrap(Array, cone_[], VecGetLocalSize(petsclib, x); own = false)
	ornt = unsafe_wrap(Array, ornt_[], VecGetLocalSize(petsclib, x); own = false)

	return rt,Nt,size,cone,ornt
end 

"""
	rnew::PetscInt,onew::PetscInt = DMPlexTransformGetSubcellOrientation(petsclib::PetscLibType,tr::DMPlexTransform, sct::DMPolytopeType, sp::PetscInt, so::PetscInt, tct::DMPolytopeType, r::PetscInt, o::PetscInt) 
Transform the replica number and orientation for a target point according to the group action for the source point

Not Collective

Input Parameters:
- `tr`  - The `DMPlexTransform`
- `sct` - The source point cell type, from whom the new cell is being produced
- `sp`  - The source point
- `so`  - The orientation of the source point in its enclosing parent
- `tct` - The target point cell type
- `r`   - The replica number requested for the produced cell type
- `o`   - The orientation of the replica

Output Parameters:
- `rnew` - The replica number, given the orientation of the parent
- `onew` - The replica orientation, given the orientation of the parent

Level: advanced

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `DMPlexTransform`, `DMPolytopeType`, `DMPlexTransformCellTransform()`, `DMPlexTransformApply()`

# External Links
$(_doc_external("Dm/DMPlexTransformGetSubcellOrientation"))
"""
function DMPlexTransformGetSubcellOrientation(petsclib::PetscLibType, tr::DMPlexTransform, sct::DMPolytopeType, sp::PetscInt, so::PetscInt, tct::DMPolytopeType, r::PetscInt, o::PetscInt) end

@for_petsc function DMPlexTransformGetSubcellOrientation(petsclib::$UnionPetscLib, tr::DMPlexTransform, sct::DMPolytopeType, sp::$PetscInt, so::$PetscInt, tct::DMPolytopeType, r::$PetscInt, o::$PetscInt )
	rnew_ = Ref{$PetscInt}()
	onew_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMPlexTransformGetSubcellOrientation, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, DMPolytopeType, $PetscInt, $PetscInt, DMPolytopeType, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
               tr, sct, sp, so, tct, r, o, rnew_, onew_,
              )

	rnew = rnew_[]
	onew = onew_[]

	return rnew,onew
end 

"""
	coneSize::PetscInt = DMPlexTransformGetConeSize(petsclib::PetscLibType,tr::DMPlexTransform, q::PetscInt) 

# External Links
$(_doc_external("Dm/DMPlexTransformGetConeSize"))
"""
function DMPlexTransformGetConeSize(petsclib::PetscLibType, tr::DMPlexTransform, q::PetscInt) end

@for_petsc function DMPlexTransformGetConeSize(petsclib::$UnionPetscLib, tr::DMPlexTransform, q::$PetscInt )
	coneSize_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMPlexTransformGetConeSize, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, $PetscInt, Ptr{$PetscInt}),
               tr, q, coneSize_,
              )

	coneSize = coneSize_[]

	return coneSize
end 

"""
	cone::Vector{PetscInt},ornt::Vector{PetscInt} = DMPlexTransformGetConeOriented(petsclib::PetscLibType,tr::DMPlexTransform, q::PetscInt, po::PetscInt) 

# External Links
$(_doc_external("Dm/DMPlexTransformGetConeOriented"))
"""
function DMPlexTransformGetConeOriented(petsclib::PetscLibType, tr::DMPlexTransform, q::PetscInt, po::PetscInt) end

@for_petsc function DMPlexTransformGetConeOriented(petsclib::$UnionPetscLib, tr::DMPlexTransform, q::$PetscInt, po::$PetscInt )
	cone_ = Ref{Ptr{$PetscInt}}()
	ornt_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:DMPlexTransformGetConeOriented, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, $PetscInt, $PetscInt, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}),
               tr, q, po, cone_, ornt_,
              )

	cone = unsafe_wrap(Array, cone_[], VecGetLocalSize(petsclib, x); own = false)
	ornt = unsafe_wrap(Array, ornt_[], VecGetLocalSize(petsclib, x); own = false)

	return cone,ornt
end 

"""
	cone::Vector{PetscInt},ornt::Vector{PetscInt} = DMPlexTransformGetCone(petsclib::PetscLibType,tr::DMPlexTransform, q::PetscInt) 

# External Links
$(_doc_external("Dm/DMPlexTransformGetCone"))
"""
function DMPlexTransformGetCone(petsclib::PetscLibType, tr::DMPlexTransform, q::PetscInt) end

@for_petsc function DMPlexTransformGetCone(petsclib::$UnionPetscLib, tr::DMPlexTransform, q::$PetscInt )
	cone_ = Ref{Ptr{$PetscInt}}()
	ornt_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:DMPlexTransformGetCone, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, $PetscInt, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}),
               tr, q, cone_, ornt_,
              )

	cone = unsafe_wrap(Array, cone_[], VecGetLocalSize(petsclib, x); own = false)
	ornt = unsafe_wrap(Array, ornt_[], VecGetLocalSize(petsclib, x); own = false)

	return cone,ornt
end 

"""
	cone::Vector{PetscInt},ornt::Vector{PetscInt} = DMPlexTransformRestoreCone(petsclib::PetscLibType,tr::DMPlexTransform, q::PetscInt) 

# External Links
$(_doc_external("Dm/DMPlexTransformRestoreCone"))
"""
function DMPlexTransformRestoreCone(petsclib::PetscLibType, tr::DMPlexTransform, q::PetscInt) end

@for_petsc function DMPlexTransformRestoreCone(petsclib::$UnionPetscLib, tr::DMPlexTransform, q::$PetscInt )
	cone_ = Ref{Ptr{$PetscInt}}()
	ornt_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:DMPlexTransformRestoreCone, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, $PetscInt, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}),
               tr, q, cone_, ornt_,
              )

	cone = unsafe_wrap(Array, cone_[], VecGetLocalSize(petsclib, x); own = false)
	ornt = unsafe_wrap(Array, ornt_[], VecGetLocalSize(petsclib, x); own = false)

	return cone,ornt
end 

"""
	Nv::PetscInt,trVerts::Vector{PetscScalar} = DMPlexTransformGetCellVertices(petsclib::PetscLibType,tr::DMPlexTransform, ct::DMPolytopeType) 
Get the set of transformed vertices lying in the closure of a reference cell of given type

Input Parameters:
- `tr` - The `DMPlexTransform` object
- `ct` - The cell type

Output Parameters:
- `Nv`      - The number of transformed vertices in the closure of the reference cell of given type
- `trVerts` - The coordinates of these vertices in the reference cell

Level: developer

-seealso: `DMPLEX`, `DMPlexTransform`, `DMPolytopeType`, `DMPlexTransformGetSubcellVertices()`

# External Links
$(_doc_external("Dm/DMPlexTransformGetCellVertices"))
"""
function DMPlexTransformGetCellVertices(petsclib::PetscLibType, tr::DMPlexTransform, ct::DMPolytopeType) end

@for_petsc function DMPlexTransformGetCellVertices(petsclib::$UnionPetscLib, tr::DMPlexTransform, ct::DMPolytopeType )
	Nv_ = Ref{$PetscInt}()
	trVerts_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:DMPlexTransformGetCellVertices, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, DMPolytopeType, Ptr{$PetscInt}, Ptr{Ptr{$PetscScalar}}),
               tr, ct, Nv_, trVerts_,
              )

	Nv = Nv_[]
	trVerts = unsafe_wrap(Array, trVerts_[], VecGetLocalSize(petsclib, x); own = false)

	return Nv,trVerts
end 

"""
	subVerts::Vector{PetscInt} = DMPlexTransformGetSubcellVertices(petsclib::PetscLibType,tr::DMPlexTransform, ct::DMPolytopeType, rct::DMPolytopeType, r::PetscInt) 
Get the set of transformed vertices defining a subcell in the reference cell of given type

Input Parameters:
- `tr`  - The `DMPlexTransform` object
- `ct`  - The cell type
- `rct` - The subcell type
- `r`   - The subcell index

Output Parameter:
- `subVerts` - The indices of these vertices in the set of vertices returned by `DMPlexTransformGetCellVertices()`

Level: developer

-seealso: `DMPLEX`, `DMPlexTransform`, `DMPolytopeType`, `DMPlexTransformGetCellVertices()`

# External Links
$(_doc_external("Dm/DMPlexTransformGetSubcellVertices"))
"""
function DMPlexTransformGetSubcellVertices(petsclib::PetscLibType, tr::DMPlexTransform, ct::DMPolytopeType, rct::DMPolytopeType, r::PetscInt) end

@for_petsc function DMPlexTransformGetSubcellVertices(petsclib::$UnionPetscLib, tr::DMPlexTransform, ct::DMPolytopeType, rct::DMPolytopeType, r::$PetscInt )
	subVerts_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:DMPlexTransformGetSubcellVertices, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, DMPolytopeType, DMPolytopeType, $PetscInt, Ptr{Ptr{$PetscInt}}),
               tr, ct, rct, r, subVerts_,
              )

	subVerts = unsafe_wrap(Array, subVerts_[], VecGetLocalSize(petsclib, x); own = false)

	return subVerts
end 

"""
	out::Vector{PetscScalar} = DMPlexTransformMapCoordinates(petsclib::PetscLibType,tr::DMPlexTransform, pct::DMPolytopeType, ct::DMPolytopeType, p::PetscInt, r::PetscInt, Nv::PetscInt, dE::PetscInt, in::Vector{PetscScalar}) 
Calculate new coordinates for produced points

Not collective

Input Parameters:
- `tr`  - The `DMPlexTransform`
- `pct` - The cell type of the parent, from whom the new cell is being produced
- `ct`  - The type being produced
- `p`   - The original point
- `r`   - The replica number requested for the produced cell type
- `Nv`  - Number of vertices in the closure of the parent cell
- `dE`  - Spatial dimension
- `in`  - array of size Nv*dE, holding coordinates of the vertices in the closure of the parent cell

Output Parameter:
- `out` - The coordinates of the new vertices

Level: intermediate

-seealso: [](ch_unstructured), `DM`, `DMPLEX`, `DMPlexTransform`, `DMPolytopeType`, `DMPlexTransformApply()`

# External Links
$(_doc_external("Dm/DMPlexTransformMapCoordinates"))
"""
function DMPlexTransformMapCoordinates(petsclib::PetscLibType, tr::DMPlexTransform, pct::DMPolytopeType, ct::DMPolytopeType, p::PetscInt, r::PetscInt, Nv::PetscInt, dE::PetscInt, in::Vector{PetscScalar}) end

@for_petsc function DMPlexTransformMapCoordinates(petsclib::$UnionPetscLib, tr::DMPlexTransform, pct::DMPolytopeType, ct::DMPolytopeType, p::$PetscInt, r::$PetscInt, Nv::$PetscInt, dE::$PetscInt, in::Vector{$PetscScalar} )
	out = Vector{$PetscScalar}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:DMPlexTransformMapCoordinates, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, DMPolytopeType, DMPolytopeType, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{$PetscScalar}),
               tr, pct, ct, p, r, Nv, dE, in, out,
              )


	return out
end 

"""
	DMPlexTransformCreateDiscLabels(petsclib::PetscLibType,tr::DMPlexTransform, rdm::PetscDM) 

# External Links
$(_doc_external("Dm/DMPlexTransformCreateDiscLabels"))
"""
function DMPlexTransformCreateDiscLabels(petsclib::PetscLibType, tr::DMPlexTransform, rdm::PetscDM) end

@for_petsc function DMPlexTransformCreateDiscLabels(petsclib::$UnionPetscLib, tr::DMPlexTransform, rdm::PetscDM )

    @chk ccall(
               (:DMPlexTransformCreateDiscLabels, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, CDM),
               tr, rdm,
              )


	return nothing
end 

"""
	DMPlexTransformApply(petsclib::PetscLibType,tr::DMPlexTransform, dm::PetscDM, tdm::PetscDM) 
Execute the transformation, producing another `DM`

Collective

Input Parameters:
- `tr` - The `DMPlexTransform` object
- `dm` - The original `DM`

Output Parameter:
- `tdm` - The transformed `DM`

Level: intermediate

Options Database Keys:
- `-dm_plex_transform_label_match_strata`      - Only label points of the same stratum as the producing point
- `-dm_plex_transform_label_replica_inc <num>` - Increment for the label value to be multiplied by the replica number
- `-dm_plex_transform_active <name>`           - Name for active mesh label

-seealso: [](plex_transform_table), [](ch_unstructured), `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformCreate()`, `DMPlexTransformSetDM()`

# External Links
$(_doc_external("Dm/DMPlexTransformApply"))
"""
function DMPlexTransformApply(petsclib::PetscLibType, tr::DMPlexTransform, dm::PetscDM, tdm::PetscDM) end

@for_petsc function DMPlexTransformApply(petsclib::$UnionPetscLib, tr::DMPlexTransform, dm::PetscDM, tdm::PetscDM )
	tdm_ = Ref(tdm.ptr)

    @chk ccall(
               (:DMPlexTransformApply, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, CDM, Ptr{CDM}),
               tr, dm, tdm_,
              )

	tdm.ptr = C_NULL

	return nothing
end 

"""
	DMPlexTransformAdaptLabel(petsclib::PetscLibType,dm::PetscDM, metric::PetscVec, adaptLabel::DMLabel, rgLabel::DMLabel, rdm::PetscDM) 

# External Links
$(_doc_external("Dm/DMPlexTransformAdaptLabel"))
"""
function DMPlexTransformAdaptLabel(petsclib::PetscLibType, dm::PetscDM, metric::PetscVec, adaptLabel::DMLabel, rgLabel::DMLabel, rdm::PetscDM) end

@for_petsc function DMPlexTransformAdaptLabel(petsclib::$UnionPetscLib, dm::PetscDM, metric::PetscVec, adaptLabel::DMLabel, rgLabel::DMLabel, rdm::PetscDM )
	rdm_ = Ref(rdm.ptr)

    @chk ccall(
               (:DMPlexTransformAdaptLabel, $petsc_library),
               PetscErrorCode,
               (CDM, CVec, DMLabel, DMLabel, Ptr{CDM}),
               dm, metric, adaptLabel, rgLabel, rdm_,
              )

	rdm.ptr = C_NULL

	return nothing
end 

"""
	useTensor::PetscBool = DMPlexTransformCohesiveExtrudeGetTensor(petsclib::PetscLibType,tr::DMPlexTransform) 
Get the flag to use tensor cells

Not Collective

Input Parameter:
- `tr` - The `DMPlexTransform`

Output Parameter:
- `useTensor` - The flag to use tensor cells

-seealso: `DMPlexTransform`, `DMPlexTransformCohesiveExtrudeSetTensor()`, `DMPlexTransformExtrudeGetTensor()`

# External Links
$(_doc_external("Dm/DMPlexTransformCohesiveExtrudeGetTensor"))
"""
function DMPlexTransformCohesiveExtrudeGetTensor(petsclib::PetscLibType, tr::DMPlexTransform) end

@for_petsc function DMPlexTransformCohesiveExtrudeGetTensor(petsclib::$UnionPetscLib, tr::DMPlexTransform )
	useTensor_ = Ref{PetscBool}()

    @chk ccall(
               (:DMPlexTransformCohesiveExtrudeGetTensor, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{PetscBool}),
               tr, useTensor_,
              )

	useTensor = useTensor_[]

	return useTensor
end 

"""
	DMPlexTransformCohesiveExtrudeSetTensor(petsclib::PetscLibType,tr::DMPlexTransform, useTensor::PetscBool) 
Set the flag to use tensor cells

Not Collective

Input Parameters:
- `tr`        - The `DMPlexTransform`
- `useTensor` - The flag for tensor cells

-seealso: `DMPlexTransform`, `DMPlexTransformCohesiveExtrudeGetTensor()`, `DMPlexTransformExtrudeSetTensor()`

# External Links
$(_doc_external("Dm/DMPlexTransformCohesiveExtrudeSetTensor"))
"""
function DMPlexTransformCohesiveExtrudeSetTensor(petsclib::PetscLibType, tr::DMPlexTransform, useTensor::PetscBool) end

@for_petsc function DMPlexTransformCohesiveExtrudeSetTensor(petsclib::$UnionPetscLib, tr::DMPlexTransform, useTensor::PetscBool )

    @chk ccall(
               (:DMPlexTransformCohesiveExtrudeSetTensor, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, PetscBool),
               tr, useTensor,
              )


	return nothing
end 

"""
	width::PetscReal = DMPlexTransformCohesiveExtrudeGetWidth(petsclib::PetscLibType,tr::DMPlexTransform) 
Get the width of extruded cells

Not Collective

Input Parameter:
- `tr` - The `DMPlexTransform`

Output Parameter:
- `width` - The width of extruded cells, or 0.

Level: intermediate

-seealso: `DMPlexTransform`, `DMPlexTransformCohesiveExtrudeSetWidth()`

# External Links
$(_doc_external("Dm/DMPlexTransformCohesiveExtrudeGetWidth"))
"""
function DMPlexTransformCohesiveExtrudeGetWidth(petsclib::PetscLibType, tr::DMPlexTransform) end

@for_petsc function DMPlexTransformCohesiveExtrudeGetWidth(petsclib::$UnionPetscLib, tr::DMPlexTransform )
	width_ = Ref{$PetscReal}()

    @chk ccall(
               (:DMPlexTransformCohesiveExtrudeGetWidth, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{$PetscReal}),
               tr, width_,
              )

	width = width_[]

	return width
end 

"""
	DMPlexTransformCohesiveExtrudeSetWidth(petsclib::PetscLibType,tr::DMPlexTransform, width::PetscReal) 
Set the width of extruded cells

Not Collective

Input Parameters:
- `tr`    - The `DMPlexTransform`
- `width` - The width of the extruded cells, or 0.

Level: intermediate

-seealso: `DMPlexTransform`, `DMPlexTransformCohesiveExtrudeGetWidth()`

# External Links
$(_doc_external("Dm/DMPlexTransformCohesiveExtrudeSetWidth"))
"""
function DMPlexTransformCohesiveExtrudeSetWidth(petsclib::PetscLibType, tr::DMPlexTransform, width::PetscReal) end

@for_petsc function DMPlexTransformCohesiveExtrudeSetWidth(petsclib::$UnionPetscLib, tr::DMPlexTransform, width::$PetscReal )

    @chk ccall(
               (:DMPlexTransformCohesiveExtrudeSetWidth, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, $PetscReal),
               tr, width,
              )


	return nothing
end 

"""
	DMPlexTransformCohesiveExtrudeGetUnsplit(petsclib::PetscLibType,tr::DMPlexTransform, unsplit::DMLabel) 
Get a new label marking the unsplit points in the transformed mesh

Not Collective

Input Parameter:
- `tr` - The `DMPlexTransform`

Output Parameter:
- `unsplit` - A new `DMLabel` marking the unsplit points in the transformed mesh

Level: intermediate

-seealso: `DMPlexTransform`, `DMPlexTransformGetTransformTypes()`

# External Links
$(_doc_external("Dm/DMPlexTransformCohesiveExtrudeGetUnsplit"))
"""
function DMPlexTransformCohesiveExtrudeGetUnsplit(petsclib::PetscLibType, tr::DMPlexTransform, unsplit::DMLabel) end

@for_petsc function DMPlexTransformCohesiveExtrudeGetUnsplit(petsclib::$UnionPetscLib, tr::DMPlexTransform, unsplit::DMLabel )

    @chk ccall(
               (:DMPlexTransformCohesiveExtrudeGetUnsplit, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{DMLabel}),
               tr, unsplit,
              )


	return nothing
end 

"""
	layers::PetscInt = DMPlexTransformExtrudeGetLayers(petsclib::PetscLibType,tr::DMPlexTransform) 
Get the number of extruded layers.

Not Collective

Input Parameter:
- `tr` - The `DMPlexTransform`

Output Parameter:
- `layers` - The number of layers

Level: intermediate

-seealso: `DMPlexTransform`, `DMPlexTransformExtrudeSetLayers()`

# External Links
$(_doc_external("Dm/DMPlexTransformExtrudeGetLayers"))
"""
function DMPlexTransformExtrudeGetLayers(petsclib::PetscLibType, tr::DMPlexTransform) end

@for_petsc function DMPlexTransformExtrudeGetLayers(petsclib::$UnionPetscLib, tr::DMPlexTransform )
	layers_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMPlexTransformExtrudeGetLayers, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{$PetscInt}),
               tr, layers_,
              )

	layers = layers_[]

	return layers
end 

"""
	DMPlexTransformExtrudeSetLayers(petsclib::PetscLibType,tr::DMPlexTransform, layers::PetscInt) 
Set the number of extruded layers.

Not Collective

Input Parameters:
- `tr`     - The `DMPlexTransform`
- `layers` - The number of layers

Level: intermediate

-seealso: `DMPlexTransform`, `DMPlexTransformExtrudeGetLayers()`

# External Links
$(_doc_external("Dm/DMPlexTransformExtrudeSetLayers"))
"""
function DMPlexTransformExtrudeSetLayers(petsclib::PetscLibType, tr::DMPlexTransform, layers::PetscInt) end

@for_petsc function DMPlexTransformExtrudeSetLayers(petsclib::$UnionPetscLib, tr::DMPlexTransform, layers::$PetscInt )

    @chk ccall(
               (:DMPlexTransformExtrudeSetLayers, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, $PetscInt),
               tr, layers,
              )


	return nothing
end 

"""
	thickness::PetscReal = DMPlexTransformExtrudeGetThickness(petsclib::PetscLibType,tr::DMPlexTransform) 
Get the total thickness of the layers

Not Collective

Input Parameter:
- `tr` - The `DMPlexTransform`

Output Parameter:
- `thickness` - The total thickness of the layers

Level: intermediate

-seealso: `DMPlexTransform`, `DMPlexTransformExtrudeSetThickness()`

# External Links
$(_doc_external("Dm/DMPlexTransformExtrudeGetThickness"))
"""
function DMPlexTransformExtrudeGetThickness(petsclib::PetscLibType, tr::DMPlexTransform) end

@for_petsc function DMPlexTransformExtrudeGetThickness(petsclib::$UnionPetscLib, tr::DMPlexTransform )
	thickness_ = Ref{$PetscReal}()

    @chk ccall(
               (:DMPlexTransformExtrudeGetThickness, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{$PetscReal}),
               tr, thickness_,
              )

	thickness = thickness_[]

	return thickness
end 

"""
	DMPlexTransformExtrudeSetThickness(petsclib::PetscLibType,tr::DMPlexTransform, thickness::PetscReal) 
Set the total thickness of the layers

Not Collective

Input Parameters:
- `tr`        - The `DMPlexTransform`
- `thickness` - The total thickness of the layers

Level: intermediate

-seealso: `DMPlexTransform`, `DMPlexTransformExtrudeGetThickness()`

# External Links
$(_doc_external("Dm/DMPlexTransformExtrudeSetThickness"))
"""
function DMPlexTransformExtrudeSetThickness(petsclib::PetscLibType, tr::DMPlexTransform, thickness::PetscReal) end

@for_petsc function DMPlexTransformExtrudeSetThickness(petsclib::$UnionPetscLib, tr::DMPlexTransform, thickness::$PetscReal )

    @chk ccall(
               (:DMPlexTransformExtrudeSetThickness, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, $PetscReal),
               tr, thickness,
              )


	return nothing
end 

"""
	useTensor::PetscBool = DMPlexTransformExtrudeGetTensor(petsclib::PetscLibType,tr::DMPlexTransform) 
Get the flag to use tensor cells

Not Collective

Input Parameter:
- `tr` - The `DMPlexTransform`

Output Parameter:
- `useTensor` - The flag to use tensor cells

-seealso: `DMPlexTransform`, `DMPlexTransformExtrudeSetTensor()`

# External Links
$(_doc_external("Dm/DMPlexTransformExtrudeGetTensor"))
"""
function DMPlexTransformExtrudeGetTensor(petsclib::PetscLibType, tr::DMPlexTransform) end

@for_petsc function DMPlexTransformExtrudeGetTensor(petsclib::$UnionPetscLib, tr::DMPlexTransform )
	useTensor_ = Ref{PetscBool}()

    @chk ccall(
               (:DMPlexTransformExtrudeGetTensor, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{PetscBool}),
               tr, useTensor_,
              )

	useTensor = useTensor_[]

	return useTensor
end 

"""
	DMPlexTransformExtrudeSetTensor(petsclib::PetscLibType,tr::DMPlexTransform, useTensor::PetscBool) 
Set the flag to use tensor cells

Not Collective

Input Parameters:
- `tr`        - The `DMPlexTransform`
- `useTensor` - The flag for tensor cells

-seealso: `DMPlexTransform`, `DMPlexTransformExtrudeGetTensor()`

# External Links
$(_doc_external("Dm/DMPlexTransformExtrudeSetTensor"))
"""
function DMPlexTransformExtrudeSetTensor(petsclib::PetscLibType, tr::DMPlexTransform, useTensor::PetscBool) end

@for_petsc function DMPlexTransformExtrudeSetTensor(petsclib::$UnionPetscLib, tr::DMPlexTransform, useTensor::PetscBool )

    @chk ccall(
               (:DMPlexTransformExtrudeSetTensor, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, PetscBool),
               tr, useTensor,
              )


	return nothing
end 

"""
	symmetric::PetscBool = DMPlexTransformExtrudeGetSymmetric(petsclib::PetscLibType,tr::DMPlexTransform) 
Get the flag to extrude symmetrically from the initial surface

Not Collective

Input Parameter:
- `tr` - The `DMPlexTransform`

Output Parameter:
- `symmetric` - The flag to extrude symmetrically

Level: intermediate

-seealso: `DMPlexTransform`, `DMPlexTransformExtrudeSetSymmetric()`, `DMPlexTransformExtrudeGetPeriodic()`

# External Links
$(_doc_external("Dm/DMPlexTransformExtrudeGetSymmetric"))
"""
function DMPlexTransformExtrudeGetSymmetric(petsclib::PetscLibType, tr::DMPlexTransform) end

@for_petsc function DMPlexTransformExtrudeGetSymmetric(petsclib::$UnionPetscLib, tr::DMPlexTransform )
	symmetric_ = Ref{PetscBool}()

    @chk ccall(
               (:DMPlexTransformExtrudeGetSymmetric, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{PetscBool}),
               tr, symmetric_,
              )

	symmetric = symmetric_[]

	return symmetric
end 

"""
	DMPlexTransformExtrudeSetSymmetric(petsclib::PetscLibType,tr::DMPlexTransform, symmetric::PetscBool) 
Set the flag to extrude symmetrically from the initial surface

Not Collective

Input Parameters:
- `tr`        - The `DMPlexTransform`
- `symmetric` - The flag to extrude symmetrically

Level: intermediate

-seealso: `DMPlexTransform`, `DMPlexTransformExtrudeGetSymmetric()`, `DMPlexTransformExtrudeSetPeriodic()`

# External Links
$(_doc_external("Dm/DMPlexTransformExtrudeSetSymmetric"))
"""
function DMPlexTransformExtrudeSetSymmetric(petsclib::PetscLibType, tr::DMPlexTransform, symmetric::PetscBool) end

@for_petsc function DMPlexTransformExtrudeSetSymmetric(petsclib::$UnionPetscLib, tr::DMPlexTransform, symmetric::PetscBool )

    @chk ccall(
               (:DMPlexTransformExtrudeSetSymmetric, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, PetscBool),
               tr, symmetric,
              )


	return nothing
end 

"""
	periodic::PetscBool = DMPlexTransformExtrudeGetPeriodic(petsclib::PetscLibType,tr::DMPlexTransform) 
Get the flag to extrude periodically from the initial surface

Not Collective

Input Parameter:
- `tr` - The `DMPlexTransform`

Output Parameter:
- `periodic` - The flag to extrude periodically

Level: intermediate

-seealso: `DMPlexTransform`, `DMPlexTransformExtrudeSetPeriodic()`, `DMPlexTransformExtrudeGetSymmetric()`

# External Links
$(_doc_external("Dm/DMPlexTransformExtrudeGetPeriodic"))
"""
function DMPlexTransformExtrudeGetPeriodic(petsclib::PetscLibType, tr::DMPlexTransform) end

@for_petsc function DMPlexTransformExtrudeGetPeriodic(petsclib::$UnionPetscLib, tr::DMPlexTransform )
	periodic_ = Ref{PetscBool}()

    @chk ccall(
               (:DMPlexTransformExtrudeGetPeriodic, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{PetscBool}),
               tr, periodic_,
              )

	periodic = periodic_[]

	return periodic
end 

"""
	DMPlexTransformExtrudeSetPeriodic(petsclib::PetscLibType,tr::DMPlexTransform, periodic::PetscBool) 
Set the flag to extrude periodically from the initial surface

Not Collective

Input Parameters:
- `tr`       - The `DMPlexTransform`
- `periodic` - The flag to extrude periodically

Level: intermediate

-seealso: `DMPlexTransform`, `DMPlexTransformExtrudeGetPeriodic()`, `DMPlexTransformExtrudeSetSymmetric()`

# External Links
$(_doc_external("Dm/DMPlexTransformExtrudeSetPeriodic"))
"""
function DMPlexTransformExtrudeSetPeriodic(petsclib::PetscLibType, tr::DMPlexTransform, periodic::PetscBool) end

@for_petsc function DMPlexTransformExtrudeSetPeriodic(petsclib::$UnionPetscLib, tr::DMPlexTransform, periodic::PetscBool )

    @chk ccall(
               (:DMPlexTransformExtrudeSetPeriodic, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, PetscBool),
               tr, periodic,
              )


	return nothing
end 

"""
	normal::Vector{PetscReal} = DMPlexTransformExtrudeGetNormal(petsclib::PetscLibType,tr::DMPlexTransform) 
Get the extrusion normal vector

Not Collective

Input Parameter:
- `tr` - The `DMPlexTransform`

Output Parameter:
- `normal` - The extrusion direction

-seealso: `DMPlexTransform`, `DMPlexTransformExtrudeSetNormal()`

# External Links
$(_doc_external("Dm/DMPlexTransformExtrudeGetNormal"))
"""
function DMPlexTransformExtrudeGetNormal(petsclib::PetscLibType, tr::DMPlexTransform) end

@for_petsc function DMPlexTransformExtrudeGetNormal(petsclib::$UnionPetscLib, tr::DMPlexTransform )
	normal = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:DMPlexTransformExtrudeGetNormal, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{$PetscReal}),
               tr, normal,
              )


	return normal
end 

"""
	DMPlexTransformExtrudeSetNormal(petsclib::PetscLibType,tr::DMPlexTransform, normal::Vector{PetscReal}) 
Set the extrusion normal

Not Collective

Input Parameters:
- `tr`     - The `DMPlexTransform`
- `normal` - The extrusion direction

Level: intermediate

-seealso: `DMPlexTransform`, `DMPlexTransformExtrudeGetNormal()`

# External Links
$(_doc_external("Dm/DMPlexTransformExtrudeSetNormal"))
"""
function DMPlexTransformExtrudeSetNormal(petsclib::PetscLibType, tr::DMPlexTransform, normal::Vector{PetscReal}) end

@for_petsc function DMPlexTransformExtrudeSetNormal(petsclib::$UnionPetscLib, tr::DMPlexTransform, normal::Vector{$PetscReal} )

    @chk ccall(
               (:DMPlexTransformExtrudeSetNormal, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{$PetscReal}),
               tr, normal,
              )


	return nothing
end 

"""
	DMPlexTransformExtrudeSetNormalFunction(petsclib::PetscLibType,tr::DMPlexTransform, normalFunc::PetscSimplePoCintFn) 
Set a function to determine the extrusion normal

Not Collective

Input Parameters:
- `tr`         - The `DMPlexTransform`
- `normalFunc` - A function determining the extrusion direction, see `PetscSimplePointFn` for the calling sequence

Level: intermediate

-seealso: `DMPlexTransform`, `DMPlexTransformExtrudeGetNormal()`, `PetscSimplePointFn`

# External Links
$(_doc_external("Dm/DMPlexTransformExtrudeSetNormalFunction"))
"""
function DMPlexTransformExtrudeSetNormalFunction(petsclib::PetscLibType, tr::DMPlexTransform, normalFunc::PetscSimplePoCintFn) end

@for_petsc function DMPlexTransformExtrudeSetNormalFunction(petsclib::$UnionPetscLib, tr::DMPlexTransform, normalFunc::PetscSimplePoCintFn )

    @chk ccall(
               (:DMPlexTransformExtrudeSetNormalFunction, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{PetscSimplePoCintFn}),
               tr, normalFunc,
              )


	return nothing
end 

"""
	DMPlexTransformExtrudeSetThicknesses(petsclib::PetscLibType,tr::DMPlexTransform, Nth::PetscInt, thicknesses::Vector{PetscReal}) 
Set the thickness of each layer

Not Collective

Input Parameters:
- `tr`          - The `DMPlexTransform`
- `Nth`         - The number of thicknesses
- `thicknesses` - The array of thicknesses

Level: intermediate

-seealso: `DMPlexTransform`, `DMPlexTransformExtrudeSetThickness()`, `DMPlexTransformExtrudeGetThickness()`

# External Links
$(_doc_external("Dm/DMPlexTransformExtrudeSetThicknesses"))
"""
function DMPlexTransformExtrudeSetThicknesses(petsclib::PetscLibType, tr::DMPlexTransform, Nth::PetscInt, thicknesses::Vector{PetscReal}) end

@for_petsc function DMPlexTransformExtrudeSetThicknesses(petsclib::$UnionPetscLib, tr::DMPlexTransform, Nth::$PetscInt, thicknesses::Vector{$PetscReal} )

    @chk ccall(
               (:DMPlexTransformExtrudeSetThicknesses, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, $PetscInt, Ptr{$PetscReal}),
               tr, Nth, thicknesses,
              )


	return nothing
end 

"""
	monitorptr::DMNetworkMonitor = DMNetworkMonitorCreate(petsclib::PetscLibType,network::PetscDM) 
Creates a network monitor context

Collective

Input Parameter:
- `network` - network to monitor

Output Parameter:
- `monitorptr` - the `DMNetworkMonitor` object

Level: intermediate

-seealso: `DM`, `DMNETWORK`, `DMNetworkMonitor`, `DMNetworkMonitorDestroy()`, `DMNetworkMonitorAdd()`

# External Links
$(_doc_external("Dm/DMNetworkMonitorCreate"))
"""
function DMNetworkMonitorCreate(petsclib::PetscLibType, network::PetscDM) end

@for_petsc function DMNetworkMonitorCreate(petsclib::$UnionPetscLib, network::PetscDM )
	monitorptr_ = Ref{DMNetworkMonitor}()

    @chk ccall(
               (:DMNetworkMonitorCreate, $petsc_library),
               PetscErrorCode,
               (CDM, Ptr{DMNetworkMonitor}),
               network, monitorptr_,
              )

	monitorptr = monitorptr_[]

	return monitorptr
end 

"""
	DMNetworkMonitorDestroy(petsclib::PetscLibType,monitor::DMNetworkMonitor) 
Destroys a network monitor and all associated viewers

Collective

Input Parameter:
- `monitor` - monitor to destroy

Level: intermediate

-seealso: `DM`, `DMNETWORK`, `DMNetworkMonitor`, `DMNetworkMonitorCreate()`, `DMNetworkMonitorAdd()`

# External Links
$(_doc_external("Dm/DMNetworkMonitorDestroy"))
"""
function DMNetworkMonitorDestroy(petsclib::PetscLibType, monitor::DMNetworkMonitor) end

@for_petsc function DMNetworkMonitorDestroy(petsclib::$UnionPetscLib, monitor::DMNetworkMonitor )

    @chk ccall(
               (:DMNetworkMonitorDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{DMNetworkMonitor},),
               monitor,
              )


	return nothing
end 

"""
	DMNetworkMonitorPop(petsclib::PetscLibType,monitor::DMNetworkMonitor) 
Removes the most recently added viewer to a `DMNetworkMonitor`

Collective

Input Parameter:
- `monitor` - the monitor

Level: intermediate

-seealso: `DM`, `DMNETWORK`, `DMNetworkMonitor`, `DMNetworkMonitorCreate()`, `DMNetworkMonitorDestroy()`

# External Links
$(_doc_external("Dm/DMNetworkMonitorPop"))
"""
function DMNetworkMonitorPop(petsclib::PetscLibType, monitor::DMNetworkMonitor) end

@for_petsc function DMNetworkMonitorPop(petsclib::$UnionPetscLib, monitor::DMNetworkMonitor )

    @chk ccall(
               (:DMNetworkMonitorPop, $petsc_library),
               PetscErrorCode,
               (DMNetworkMonitor,),
               monitor,
              )


	return nothing
end 

"""
	DMNetworkMonitorAdd(petsclib::PetscLibType,monitor::DMNetworkMonitor, name::String, element::PetscInt, nodes::PetscInt, start::PetscInt, blocksize::PetscInt, xmin::PetscReal, xmax::PetscReal, ymin::PetscReal, ymax::PetscReal, hold::PetscBool) 
Adds a new viewer to a `DMNetworkMonitor`

Collective

Input Parameters:
- `monitor`   - the monitor
- `name`      - name of viewer
- `element`   - vertex / edge number
- `nodes`     - number of nodes
- `start`     - variable starting offset
- `blocksize` - variable blocksize
- `xmin`      - xmin (or `PETSC_DECIDE`) for viewer
- `xmax`      - xmax (or `PETSC_DECIDE`) for viewer
- `ymin`      - ymin for viewer
- `ymax`      - ymax for viewer
- `hold`      - determines if plot limits should be held

Level: intermediate

-seealso: `DM`, `DMNETWORK`, `DMNetworkMonitor`, `DMNetworkMonitorCreate()`, `DMNetworkMonitorDestroy()`

# External Links
$(_doc_external("Dm/DMNetworkMonitorAdd"))
"""
function DMNetworkMonitorAdd(petsclib::PetscLibType, monitor::DMNetworkMonitor, name::String, element::PetscInt, nodes::PetscInt, start::PetscInt, blocksize::PetscInt, xmin::PetscReal, xmax::PetscReal, ymin::PetscReal, ymax::PetscReal, hold::PetscBool) end

@for_petsc function DMNetworkMonitorAdd(petsclib::$UnionPetscLib, monitor::DMNetworkMonitor, name::String, element::$PetscInt, nodes::$PetscInt, start::$PetscInt, blocksize::$PetscInt, xmin::$PetscReal, xmax::$PetscReal, ymin::$PetscReal, ymax::$PetscReal, hold::PetscBool )

    @chk ccall(
               (:DMNetworkMonitorAdd, $petsc_library),
               PetscErrorCode,
               (DMNetworkMonitor, Ptr{Cchar}, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscReal, $PetscReal, $PetscReal, $PetscReal, PetscBool),
               monitor, name, element, nodes, start, blocksize, xmin, xmax, ymin, ymax, hold,
              )


	return nothing
end 

"""
	DMNetworkMonitorView(petsclib::PetscLibType,monitor::DMNetworkMonitor, x::PetscVec) 
A `DMNETWORK` specific monitor function for `TSMonitorSet()`

Collective, No Fortran support

Input Parameters:
- `monitor` - `DMNetworkMonitor` object
- `x`       - `TS` solution vector

Level: intermediate

-seealso: `DM`, `DMNETWORK`, `DMNetworkMonitor`, `DMNetworkMonitorCreate()`, `DMNetworkMonitorDestroy()`, `DMNetworkMonitorAdd()`

# External Links
$(_doc_external("Dm/DMNetworkMonitorView"))
"""
function DMNetworkMonitorView(petsclib::PetscLibType, monitor::DMNetworkMonitor, x::PetscVec) end

@for_petsc function DMNetworkMonitorView(petsclib::$UnionPetscLib, monitor::DMNetworkMonitor, x::PetscVec )

    @chk ccall(
               (:DMNetworkMonitorView, $petsc_library),
               PetscErrorCode,
               (DMNetworkMonitor, CVec),
               monitor, x,
              )


	return nothing
end 

"""
	label::DMLabel = DMLabelCreate(petsclib::PetscLibType,comm::MPI_Comm, name::String) 
Create a `DMLabel` object, which is a multimap

Collective

Input Parameters:
- `comm` - The communicator, usually `PETSC_COMM_SELF`
- `name` - The label name

Output Parameter:
- `label` - The `DMLabel`

Level: beginner

-seealso: `DMLabel`, `DM`, `DMLabelDestroy()`

# External Links
$(_doc_external("Dm/DMLabelCreate"))
"""
function DMLabelCreate(petsclib::PetscLibType, comm::MPI_Comm, name::String) end

@for_petsc function DMLabelCreate(petsclib::$UnionPetscLib, comm::MPI_Comm, name::String )
	label_ = Ref{DMLabel}()

    @chk ccall(
               (:DMLabelCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{DMLabel}),
               comm, name, label_,
              )

	label = label_[]

	return label
end 

"""
	DMLabelSetUp(petsclib::PetscLibType,label::DMLabel) 
SetUp a `DMLabel` object

Collective

Input Parameters:
- `label` - The `DMLabel`

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelDestroy()`

# External Links
$(_doc_external("Dm/DMLabelSetUp"))
"""
function DMLabelSetUp(petsclib::PetscLibType, label::DMLabel) end

@for_petsc function DMLabelSetUp(petsclib::$UnionPetscLib, label::DMLabel )

    @chk ccall(
               (:DMLabelSetUp, $petsc_library),
               PetscErrorCode,
               (DMLabel,),
               label,
              )


	return nothing
end 

"""
	DMLabelAddStratum(petsclib::PetscLibType,label::DMLabel, value::PetscInt) 
Adds a new stratum value in a `DMLabel`

Input Parameters:
- `label` - The `DMLabel`
- `value` - The stratum value

Level: beginner

-seealso: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelDestroy()`

# External Links
$(_doc_external("Dm/DMLabelAddStratum"))
"""
function DMLabelAddStratum(petsclib::PetscLibType, label::DMLabel, value::PetscInt) end

@for_petsc function DMLabelAddStratum(petsclib::$UnionPetscLib, label::DMLabel, value::$PetscInt )

    @chk ccall(
               (:DMLabelAddStratum, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt),
               label, value,
              )


	return nothing
end 

"""
	DMLabelAddStrata(petsclib::PetscLibType,label::DMLabel, numStrata::PetscInt, stratumValues::Vector{PetscInt}) 
Adds new stratum values in a `DMLabel`

Not Collective

Input Parameters:
- `label`         - The `DMLabel`
- `numStrata`     - The number of stratum values
- `stratumValues` - The stratum values

Level: beginner

-seealso: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelDestroy()`

# External Links
$(_doc_external("Dm/DMLabelAddStrata"))
"""
function DMLabelAddStrata(petsclib::PetscLibType, label::DMLabel, numStrata::PetscInt, stratumValues::Vector{PetscInt}) end

@for_petsc function DMLabelAddStrata(petsclib::$UnionPetscLib, label::DMLabel, numStrata::$PetscInt, stratumValues::Vector{$PetscInt} )

    @chk ccall(
               (:DMLabelAddStrata, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, Ptr{$PetscInt}),
               label, numStrata, stratumValues,
              )


	return nothing
end 

"""
	DMLabelAddStrataIS(petsclib::PetscLibType,label::DMLabel, valueIS::IS) 
Adds new stratum values in a `DMLabel`

Not Collective

Input Parameters:
- `label`   - The `DMLabel`
- `valueIS` - Index set with stratum values

Level: beginner

-seealso: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelDestroy()`

# External Links
$(_doc_external("Dm/DMLabelAddStrataIS"))
"""
function DMLabelAddStrataIS(petsclib::PetscLibType, label::DMLabel, valueIS::IS) end

@for_petsc function DMLabelAddStrataIS(petsclib::$UnionPetscLib, label::DMLabel, valueIS::IS )

    @chk ccall(
               (:DMLabelAddStrataIS, $petsc_library),
               PetscErrorCode,
               (DMLabel, CIS),
               label, valueIS,
              )


	return nothing
end 

"""
	DMLabelView(petsclib::PetscLibType,label::DMLabel, viewer::PetscViewer) 
View the label

Collective

Input Parameters:
- `label`  - The `DMLabel`
- `viewer` - The `PetscViewer`

Level: intermediate

-seealso: `DMLabel`, `PetscViewer`, `DM`, `DMLabelCreate()`, `DMLabelDestroy()`

# External Links
$(_doc_external("Dm/DMLabelView"))
"""
function DMLabelView(petsclib::PetscLibType, label::DMLabel, viewer::PetscViewer) end

@for_petsc function DMLabelView(petsclib::$UnionPetscLib, label::DMLabel, viewer::PetscViewer )

    @chk ccall(
               (:DMLabelView, $petsc_library),
               PetscErrorCode,
               (DMLabel, PetscViewer),
               label, viewer,
              )


	return nothing
end 

"""
	DMLabelReset(petsclib::PetscLibType,label::DMLabel) 
Destroys internal data structures in a `DMLabel`

Not Collective

Input Parameter:
- `label` - The `DMLabel`

Level: beginner

-seealso: `DMLabel`, `DM`, `DMLabelDestroy()`, `DMLabelCreate()`

# External Links
$(_doc_external("Dm/DMLabelReset"))
"""
function DMLabelReset(petsclib::PetscLibType, label::DMLabel) end

@for_petsc function DMLabelReset(petsclib::$UnionPetscLib, label::DMLabel )

    @chk ccall(
               (:DMLabelReset, $petsc_library),
               PetscErrorCode,
               (DMLabel,),
               label,
              )


	return nothing
end 

"""
	DMLabelDestroy(petsclib::PetscLibType,label::DMLabel) 
Destroys a `DMLabel`

Collective

Input Parameter:
- `label` - The `DMLabel`

Level: beginner

-seealso: `DMLabel`, `DM`, `DMLabelReset()`, `DMLabelCreate()`

# External Links
$(_doc_external("Dm/DMLabelDestroy"))
"""
function DMLabelDestroy(petsclib::PetscLibType, label::DMLabel) end

@for_petsc function DMLabelDestroy(petsclib::$UnionPetscLib, label::DMLabel )

    @chk ccall(
               (:DMLabelDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{DMLabel},),
               label,
              )


	return nothing
end 

"""
	labelnew::DMLabel = DMLabelDuplicate(petsclib::PetscLibType,label::DMLabel) 
Duplicates a `DMLabel`

Collective

Input Parameter:
- `label` - The `DMLabel`

Output Parameter:
- `labelnew` - new label

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelDestroy()`

# External Links
$(_doc_external("Dm/DMLabelDuplicate"))
"""
function DMLabelDuplicate(petsclib::PetscLibType, label::DMLabel) end

@for_petsc function DMLabelDuplicate(petsclib::$UnionPetscLib, label::DMLabel )
	labelnew_ = Ref{DMLabel}()

    @chk ccall(
               (:DMLabelDuplicate, $petsc_library),
               PetscErrorCode,
               (DMLabel, Ptr{DMLabel}),
               label, labelnew_,
              )

	labelnew = labelnew_[]

	return labelnew
end 

"""
	equal::PetscBool = DMLabelCompare(petsclib::PetscLibType,comm::MPI_Comm, l0::DMLabel, l1::DMLabel, message::Cchar) 
Compare two `DMLabel` objects

Collective; No Fortran Support

Input Parameters:
- `comm` - Comm over which to compare labels
- `l0`   - First `DMLabel`
- `l1`   - Second `DMLabel`

Output Parameters:
- `equal`   - (Optional) Flag whether the two labels are equal
- `message` - (Optional) Message describing the difference

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMCompareLabels()`, `DMLabelGetNumValues()`, `DMLabelGetDefaultValue()`, `DMLabelGetNonEmptyStratumValuesIS()`, `DMLabelGetStratumIS()`

# External Links
$(_doc_external("Dm/DMLabelCompare"))
"""
function DMLabelCompare(petsclib::PetscLibType, comm::MPI_Comm, l0::DMLabel, l1::DMLabel, message::Cchar) end

@for_petsc function DMLabelCompare(petsclib::$UnionPetscLib, comm::MPI_Comm, l0::DMLabel, l1::DMLabel, message::Cchar )
	equal_ = Ref{PetscBool}()

    @chk ccall(
               (:DMLabelCompare, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, DMLabel, DMLabel, Ptr{PetscBool}, Cchar),
               comm, l0, l1, equal_, message,
              )

	equal = equal_[]

	return equal
end 

"""
	DMLabelComputeIndex(petsclib::PetscLibType,label::DMLabel) 
Create an index structure for membership determination, automatically determining the bounds

Not Collective

Input Parameter:
- `label` - The `DMLabel`

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelHasPoint()`, `DMLabelCreateIndex()`, `DMLabelDestroyIndex()`, `DMLabelGetValue()`, `DMLabelSetValue()`

# External Links
$(_doc_external("Dm/DMLabelComputeIndex"))
"""
function DMLabelComputeIndex(petsclib::PetscLibType, label::DMLabel) end

@for_petsc function DMLabelComputeIndex(petsclib::$UnionPetscLib, label::DMLabel )

    @chk ccall(
               (:DMLabelComputeIndex, $petsc_library),
               PetscErrorCode,
               (DMLabel,),
               label,
              )


	return nothing
end 

"""
	DMLabelCreateIndex(petsclib::PetscLibType,label::DMLabel, pStart::PetscInt, pEnd::PetscInt) 
Create an index structure for membership determination

Not Collective

Input Parameters:
- `label`  - The `DMLabel`
- `pStart` - The smallest point
- `pEnd`   - The largest point + 1

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelHasPoint()`, `DMLabelComputeIndex()`, `DMLabelDestroyIndex()`, `DMLabelGetValue()`, `DMLabelSetValue()`

# External Links
$(_doc_external("Dm/DMLabelCreateIndex"))
"""
function DMLabelCreateIndex(petsclib::PetscLibType, label::DMLabel, pStart::PetscInt, pEnd::PetscInt) end

@for_petsc function DMLabelCreateIndex(petsclib::$UnionPetscLib, label::DMLabel, pStart::$PetscInt, pEnd::$PetscInt )

    @chk ccall(
               (:DMLabelCreateIndex, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, $PetscInt),
               label, pStart, pEnd,
              )


	return nothing
end 

"""
	DMLabelDestroyIndex(petsclib::PetscLibType,label::DMLabel) 
Destroy the index structure

Not Collective

Input Parameter:
- `label` - the `DMLabel`

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelHasPoint()`, `DMLabelCreateIndex()`, `DMLabelGetValue()`, `DMLabelSetValue()`

# External Links
$(_doc_external("Dm/DMLabelDestroyIndex"))
"""
function DMLabelDestroyIndex(petsclib::PetscLibType, label::DMLabel) end

@for_petsc function DMLabelDestroyIndex(petsclib::$UnionPetscLib, label::DMLabel )

    @chk ccall(
               (:DMLabelDestroyIndex, $petsc_library),
               PetscErrorCode,
               (DMLabel,),
               label,
              )


	return nothing
end 

"""
	pStart::PetscInt,pEnd::PetscInt = DMLabelGetBounds(petsclib::PetscLibType,label::DMLabel) 
Return the smallest and largest point in the label

Not Collective

Input Parameter:
- `label` - the `DMLabel`

Output Parameters:
- `pStart` - The smallest point
- `pEnd`   - The largest point + 1

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelHasPoint()`, `DMLabelCreateIndex()`, `DMLabelGetValue()`, `DMLabelSetValue()`

# External Links
$(_doc_external("Dm/DMLabelGetBounds"))
"""
function DMLabelGetBounds(petsclib::PetscLibType, label::DMLabel) end

@for_petsc function DMLabelGetBounds(petsclib::$UnionPetscLib, label::DMLabel )
	pStart_ = Ref{$PetscInt}()
	pEnd_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMLabelGetBounds, $petsc_library),
               PetscErrorCode,
               (DMLabel, Ptr{$PetscInt}, Ptr{$PetscInt}),
               label, pStart_, pEnd_,
              )

	pStart = pStart_[]
	pEnd = pEnd_[]

	return pStart,pEnd
end 

"""
	contains::PetscBool = DMLabelHasValue(petsclib::PetscLibType,label::DMLabel, value::PetscInt) 
Determine whether a label assigns the value to any point

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `value` - the value

Output Parameter:
- `contains` - Flag indicating whether the label maps this value to any point

Level: developer

-seealso: `DMLabel`, `DM`, `DMLabelHasPoint()`, `DMLabelGetValue()`, `DMLabelSetValue()`

# External Links
$(_doc_external("Dm/DMLabelHasValue"))
"""
function DMLabelHasValue(petsclib::PetscLibType, label::DMLabel, value::PetscInt) end

@for_petsc function DMLabelHasValue(petsclib::$UnionPetscLib, label::DMLabel, value::$PetscInt )
	contains_ = Ref{PetscBool}()

    @chk ccall(
               (:DMLabelHasValue, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, Ptr{PetscBool}),
               label, value, contains_,
              )

	contains = contains_[]

	return contains
end 

"""
	contains::PetscBool = DMLabelHasPoint(petsclib::PetscLibType,label::DMLabel, point::PetscInt) 
Determine whether a label assigns a value to a point

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `point` - the point

Output Parameter:
- `contains` - Flag indicating whether the label maps this point to a value

Level: developer

-seealso: `DMLabel`, `DM`, `DMLabelCreateIndex()`, `DMLabelGetValue()`, `DMLabelSetValue()`

# External Links
$(_doc_external("Dm/DMLabelHasPoint"))
"""
function DMLabelHasPoint(petsclib::PetscLibType, label::DMLabel, point::PetscInt) end

@for_petsc function DMLabelHasPoint(petsclib::$UnionPetscLib, label::DMLabel, point::$PetscInt )
	contains_ = Ref{PetscBool}()

    @chk ccall(
               (:DMLabelHasPoint, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, Ptr{PetscBool}),
               label, point, contains_,
              )

	contains = contains_[]

	return contains
end 

"""
	contains::PetscBool = DMLabelStratumHasPoint(petsclib::PetscLibType,label::DMLabel, value::PetscInt, point::PetscInt) 
Return true if the stratum contains a point

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `value` - the stratum value
- `point` - the point

Output Parameter:
- `contains` - true if the stratum contains the point

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("Dm/DMLabelStratumHasPoint"))
"""
function DMLabelStratumHasPoint(petsclib::PetscLibType, label::DMLabel, value::PetscInt, point::PetscInt) end

@for_petsc function DMLabelStratumHasPoint(petsclib::$UnionPetscLib, label::DMLabel, value::$PetscInt, point::$PetscInt )
	contains_ = Ref{PetscBool}()

    @chk ccall(
               (:DMLabelStratumHasPoint, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, $PetscInt, Ptr{PetscBool}),
               label, value, point, contains_,
              )

	contains = contains_[]

	return contains
end 

"""
	defaultValue::PetscInt = DMLabelGetDefaultValue(petsclib::PetscLibType,label::DMLabel) 
Get the default value returned by `DMLabelGetValue()` if a point has not been explicitly given a value.
When a label is created, it is initialized to -1.

Not Collective

Input Parameter:
- `label` - a `DMLabel` object

Output Parameter:
- `defaultValue` - the default value

Level: beginner

-seealso: `DMLabel`, `DM`, `DMLabelSetDefaultValue()`, `DMLabelGetValue()`, `DMLabelSetValue()`

# External Links
$(_doc_external("Dm/DMLabelGetDefaultValue"))
"""
function DMLabelGetDefaultValue(petsclib::PetscLibType, label::DMLabel) end

@for_petsc function DMLabelGetDefaultValue(petsclib::$UnionPetscLib, label::DMLabel )
	defaultValue_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMLabelGetDefaultValue, $petsc_library),
               PetscErrorCode,
               (DMLabel, Ptr{$PetscInt}),
               label, defaultValue_,
              )

	defaultValue = defaultValue_[]

	return defaultValue
end 

"""
	defaultValue::PetscInt = DMLabelSetDefaultValue(petsclib::PetscLibType,label::DMLabel) 
Set the default value returned by `DMLabelGetValue()` if a point has not been explicitly given a value.
When a label is created, it is initialized to -1.

Not Collective

Input Parameter:
- `label` - a `DMLabel` object

Output Parameter:
- `defaultValue` - the default value

Level: beginner

-seealso: `DMLabel`, `DM`, `DMLabelGetDefaultValue()`, `DMLabelGetValue()`, `DMLabelSetValue()`

# External Links
$(_doc_external("Dm/DMLabelSetDefaultValue"))
"""
function DMLabelSetDefaultValue(petsclib::PetscLibType, label::DMLabel) end

@for_petsc function DMLabelSetDefaultValue(petsclib::$UnionPetscLib, label::DMLabel )
	defaultValue_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMLabelSetDefaultValue, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt),
               label, defaultValue_,
              )

	defaultValue = defaultValue_[]

	return defaultValue
end 

"""
	value::PetscInt = DMLabelGetValue(petsclib::PetscLibType,label::DMLabel, point::PetscInt) 
Return the value a label assigns to a point, or the label's default value (which is initially
`DMLabelSetDefaultValue()`)

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `point` - the point

Output Parameter:
- `value` - The point value, or the default value (-1 by default)

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelSetValue()`, `DMLabelClearValue()`, `DMLabelGetDefaultValue()`, `DMLabelSetDefaultValue()`

# External Links
$(_doc_external("Dm/DMLabelGetValue"))
"""
function DMLabelGetValue(petsclib::PetscLibType, label::DMLabel, point::PetscInt) end

@for_petsc function DMLabelGetValue(petsclib::$UnionPetscLib, label::DMLabel, point::$PetscInt )
	value_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMLabelGetValue, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, Ptr{$PetscInt}),
               label, point, value_,
              )

	value = value_[]

	return value
end 

"""
	DMLabelSetValue(petsclib::PetscLibType,label::DMLabel, point::PetscInt, value::PetscInt) 
Set the value a label assigns to a point.  If the value is the same as the label's default value (which is initially
be changed with `DMLabelSetDefaultValue()` to something different), then this function will do nothing.

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `point` - the point
- `value` - The point value

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelClearValue()`, `DMLabelGetDefaultValue()`, `DMLabelSetDefaultValue()`

# External Links
$(_doc_external("Dm/DMLabelSetValue"))
"""
function DMLabelSetValue(petsclib::PetscLibType, label::DMLabel, point::PetscInt, value::PetscInt) end

@for_petsc function DMLabelSetValue(petsclib::$UnionPetscLib, label::DMLabel, point::$PetscInt, value::$PetscInt )

    @chk ccall(
               (:DMLabelSetValue, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, $PetscInt),
               label, point, value,
              )


	return nothing
end 

"""
	DMLabelClearValue(petsclib::PetscLibType,label::DMLabel, point::PetscInt, value::PetscInt) 
Clear the value a label assigns to a point

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `point` - the point
- `value` - The point value

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`

# External Links
$(_doc_external("Dm/DMLabelClearValue"))
"""
function DMLabelClearValue(petsclib::PetscLibType, label::DMLabel, point::PetscInt, value::PetscInt) end

@for_petsc function DMLabelClearValue(petsclib::$UnionPetscLib, label::DMLabel, point::$PetscInt, value::$PetscInt )

    @chk ccall(
               (:DMLabelClearValue, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, $PetscInt),
               label, point, value,
              )


	return nothing
end 

"""
	DMLabelInsertIS(petsclib::PetscLibType,label::DMLabel, is::IS, value::PetscInt) 
Set all points in the `IS` to a value

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `is`    - the point `IS`
- `value` - The point value

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("Dm/DMLabelInsertIS"))
"""
function DMLabelInsertIS(petsclib::PetscLibType, label::DMLabel, is::IS, value::PetscInt) end

@for_petsc function DMLabelInsertIS(petsclib::$UnionPetscLib, label::DMLabel, is::IS, value::$PetscInt )

    @chk ccall(
               (:DMLabelInsertIS, $petsc_library),
               PetscErrorCode,
               (DMLabel, CIS, $PetscInt),
               label, is, value,
              )


	return nothing
end 

"""
	numValues::PetscInt = DMLabelGetNumValues(petsclib::PetscLibType,label::DMLabel) 
Get the number of values that the `DMLabel` takes

Not Collective

Input Parameter:
- `label` - the `DMLabel`

Output Parameter:
- `numValues` - the number of values

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("Dm/DMLabelGetNumValues"))
"""
function DMLabelGetNumValues(petsclib::PetscLibType, label::DMLabel) end

@for_petsc function DMLabelGetNumValues(petsclib::$UnionPetscLib, label::DMLabel )
	numValues_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMLabelGetNumValues, $petsc_library),
               PetscErrorCode,
               (DMLabel, Ptr{$PetscInt}),
               label, numValues_,
              )

	numValues = numValues_[]

	return numValues
end 

"""
	DMLabelGetValueIS(petsclib::PetscLibType,label::DMLabel, values::IS) 
Get an `IS` of all values that the `DMlabel` takes

Not Collective

Input Parameter:
- `label` - the `DMLabel`

Output Parameter:
- `values` - the value `IS`

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelGetNonEmptyStratumValuesIS()`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("Dm/DMLabelGetValueIS"))
"""
function DMLabelGetValueIS(petsclib::PetscLibType, label::DMLabel, values::IS) end

@for_petsc function DMLabelGetValueIS(petsclib::$UnionPetscLib, label::DMLabel, values::IS )
	values_ = Ref(values.ptr)

    @chk ccall(
               (:DMLabelGetValueIS, $petsc_library),
               PetscErrorCode,
               (DMLabel, Ptr{CIS}),
               label, values_,
              )

	values.ptr = C_NULL

	return nothing
end 

"""
	minValue::PetscInt,maxValue::PetscInt = DMLabelGetValueBounds(petsclib::PetscLibType,label::DMLabel) 
Return the smallest and largest value in the label

Not Collective

Input Parameter:
- `label` - the `DMLabel`

Output Parameters:
- `minValue` - The smallest value
- `maxValue` - The largest value

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelGetBounds()`, `DMLabelGetValue()`, `DMLabelSetValue()`

# External Links
$(_doc_external("Dm/DMLabelGetValueBounds"))
"""
function DMLabelGetValueBounds(petsclib::PetscLibType, label::DMLabel) end

@for_petsc function DMLabelGetValueBounds(petsclib::$UnionPetscLib, label::DMLabel )
	minValue_ = Ref{$PetscInt}()
	maxValue_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMLabelGetValueBounds, $petsc_library),
               PetscErrorCode,
               (DMLabel, Ptr{$PetscInt}, Ptr{$PetscInt}),
               label, minValue_, maxValue_,
              )

	minValue = minValue_[]
	maxValue = maxValue_[]

	return minValue,maxValue
end 

"""
	DMLabelGetNonEmptyStratumValuesIS(petsclib::PetscLibType,label::DMLabel, values::IS) 
Get an `IS` of all values that the `DMlabel` takes

Not Collective

Input Parameter:
- `label` - the `DMLabel`

Output Parameter:
- `values` - the value `IS`

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelGetValueIS()`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("Dm/DMLabelGetNonEmptyStratumValuesIS"))
"""
function DMLabelGetNonEmptyStratumValuesIS(petsclib::PetscLibType, label::DMLabel, values::IS) end

@for_petsc function DMLabelGetNonEmptyStratumValuesIS(petsclib::$UnionPetscLib, label::DMLabel, values::IS )
	values_ = Ref(values.ptr)

    @chk ccall(
               (:DMLabelGetNonEmptyStratumValuesIS, $petsc_library),
               PetscErrorCode,
               (DMLabel, Ptr{CIS}),
               label, values_,
              )

	values.ptr = C_NULL

	return nothing
end 

"""
	index::PetscInt = DMLabelGetValueIndex(petsclib::PetscLibType,label::DMLabel, value::PetscInt) 
Get the index of a given value in the list of values for the `DMlabel`, or

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `value` - the value

Output Parameter:
- `index` - the index of value in the list of values

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelGetValueIS()`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("Dm/DMLabelGetValueIndex"))
"""
function DMLabelGetValueIndex(petsclib::PetscLibType, label::DMLabel, value::PetscInt) end

@for_petsc function DMLabelGetValueIndex(petsclib::$UnionPetscLib, label::DMLabel, value::$PetscInt )
	index_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMLabelGetValueIndex, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, Ptr{$PetscInt}),
               label, value, index_,
              )

	index = index_[]

	return index
end 

"""
	exists::PetscBool = DMLabelHasStratum(petsclib::PetscLibType,label::DMLabel, value::PetscInt) 
Determine whether points exist with the given value

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `value` - the stratum value

Output Parameter:
- `exists` - Flag saying whether points exist

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("Dm/DMLabelHasStratum"))
"""
function DMLabelHasStratum(petsclib::PetscLibType, label::DMLabel, value::PetscInt) end

@for_petsc function DMLabelHasStratum(petsclib::$UnionPetscLib, label::DMLabel, value::$PetscInt )
	exists_ = Ref{PetscBool}()

    @chk ccall(
               (:DMLabelHasStratum, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, Ptr{PetscBool}),
               label, value, exists_,
              )

	exists = exists_[]

	return exists
end 

"""
	size::PetscInt = DMLabelGetStratumSize(petsclib::PetscLibType,label::DMLabel, value::PetscInt) 
Get the size of a stratum

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `value` - the stratum value

Output Parameter:
- `size` - The number of points in the stratum

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("Dm/DMLabelGetStratumSize"))
"""
function DMLabelGetStratumSize(petsclib::PetscLibType, label::DMLabel, value::PetscInt) end

@for_petsc function DMLabelGetStratumSize(petsclib::$UnionPetscLib, label::DMLabel, value::$PetscInt )
	size_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMLabelGetStratumSize, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, Ptr{$PetscInt}),
               label, value, size_,
              )

	size = size_[]

	return size
end 

"""
	start::PetscInt,end_::PetscInt = DMLabelGetStratumBounds(petsclib::PetscLibType,label::DMLabel, value::PetscInt) 
Get the largest and smallest point of a stratum

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `value` - the stratum value

Output Parameters:
- `start` - the smallest point in the stratum
- `end`   - the largest point in the stratum

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("Dm/DMLabelGetStratumBounds"))
"""
function DMLabelGetStratumBounds(petsclib::PetscLibType, label::DMLabel, value::PetscInt) end

@for_petsc function DMLabelGetStratumBounds(petsclib::$UnionPetscLib, label::DMLabel, value::$PetscInt )
	start_ = Ref{$PetscInt}()
	end__ = Ref{$PetscInt}()

    @chk ccall(
               (:DMLabelGetStratumBounds, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
               label, value, start_, end__,
              )

	start = start_[]
	end_ = end__[]

	return start,end_
end 

"""
	DMLabelGetStratumIS(petsclib::PetscLibType,label::DMLabel, value::PetscInt, points::IS) 
Get an `IS` with the stratum points

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `value` - the stratum value

Output Parameter:
- `points` - The stratum points

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("Dm/DMLabelGetStratumIS"))
"""
function DMLabelGetStratumIS(petsclib::PetscLibType, label::DMLabel, value::PetscInt, points::IS) end

@for_petsc function DMLabelGetStratumIS(petsclib::$UnionPetscLib, label::DMLabel, value::$PetscInt, points::IS )
	points_ = Ref(points.ptr)

    @chk ccall(
               (:DMLabelGetStratumIS, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, Ptr{CIS}),
               label, value, points_,
              )

	points.ptr = C_NULL

	return nothing
end 

"""
	DMLabelSetStratumIS(petsclib::PetscLibType,label::DMLabel, value::PetscInt, is::IS) 
Set the stratum points using an `IS`

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `value` - the stratum value
- `is`    - The stratum points

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("Dm/DMLabelSetStratumIS"))
"""
function DMLabelSetStratumIS(petsclib::PetscLibType, label::DMLabel, value::PetscInt, is::IS) end

@for_petsc function DMLabelSetStratumIS(petsclib::$UnionPetscLib, label::DMLabel, value::$PetscInt, is::IS )

    @chk ccall(
               (:DMLabelSetStratumIS, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, CIS),
               label, value, is,
              )


	return nothing
end 

"""
	DMLabelClearStratum(petsclib::PetscLibType,label::DMLabel, value::PetscInt) 
Remove a stratum

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `value` - the stratum value

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("Dm/DMLabelClearStratum"))
"""
function DMLabelClearStratum(petsclib::PetscLibType, label::DMLabel, value::PetscInt) end

@for_petsc function DMLabelClearStratum(petsclib::$UnionPetscLib, label::DMLabel, value::$PetscInt )

    @chk ccall(
               (:DMLabelClearStratum, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt),
               label, value,
              )


	return nothing
end 

"""
	DMLabelSetStratumBounds(petsclib::PetscLibType,label::DMLabel, value::PetscInt, pStart::PetscInt, pEnd::PetscInt) 
Efficiently give a contiguous set of points a given label value

Not Collective

Input Parameters:
- `label`  - The `DMLabel`
- `value`  - The label value for all points
- `pStart` - The first point
- `pEnd`   - A point beyond all marked points

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelSetStratumIS()`, `DMLabelGetStratumIS()`

# External Links
$(_doc_external("Dm/DMLabelSetStratumBounds"))
"""
function DMLabelSetStratumBounds(petsclib::PetscLibType, label::DMLabel, value::PetscInt, pStart::PetscInt, pEnd::PetscInt) end

@for_petsc function DMLabelSetStratumBounds(petsclib::$UnionPetscLib, label::DMLabel, value::$PetscInt, pStart::$PetscInt, pEnd::$PetscInt )

    @chk ccall(
               (:DMLabelSetStratumBounds, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, $PetscInt, $PetscInt),
               label, value, pStart, pEnd,
              )


	return nothing
end 

"""
	index::PetscInt = DMLabelGetStratumPointIndex(petsclib::PetscLibType,label::DMLabel, value::PetscInt, p::PetscInt) 
Get the index of a point in a given stratum

Not Collective

Input Parameters:
- `label` - The `DMLabel`
- `value` - The label value
- `p`     - A point with this value

Output Parameter:
- `index` - The index of this point in the stratum, or -1 if the point is not in the stratum or the stratum does not exist

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelGetValueIndex()`, `DMLabelGetStratumIS()`, `DMLabelCreate()`

# External Links
$(_doc_external("Dm/DMLabelGetStratumPointIndex"))
"""
function DMLabelGetStratumPointIndex(petsclib::PetscLibType, label::DMLabel, value::PetscInt, p::PetscInt) end

@for_petsc function DMLabelGetStratumPointIndex(petsclib::$UnionPetscLib, label::DMLabel, value::$PetscInt, p::$PetscInt )
	index_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMLabelGetStratumPointIndex, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, $PetscInt, Ptr{$PetscInt}),
               label, value, p, index_,
              )

	index = index_[]

	return index
end 

"""
	DMLabelFilter(petsclib::PetscLibType,label::DMLabel, start::PetscInt, end_::PetscInt) 
Remove all points outside of [`start`, `end`)

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `start` - the first point kept
- `end`   - one more than the last point kept

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("Dm/DMLabelFilter"))
"""
function DMLabelFilter(petsclib::PetscLibType, label::DMLabel, start::PetscInt, end_::PetscInt) end

@for_petsc function DMLabelFilter(petsclib::$UnionPetscLib, label::DMLabel, start::$PetscInt, end_::$PetscInt )

    @chk ccall(
               (:DMLabelFilter, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, $PetscInt),
               label, start, end_,
              )


	return nothing
end 

"""
	DMLabelPermute(petsclib::PetscLibType,label::DMLabel, permutation::IS, labelNew::DMLabel) 
Create a new label with permuted points

Not Collective

Input Parameters:
- `label`       - the `DMLabel`
- `permutation` - the point permutation

Output Parameter:
- `labelNew` - the new label containing the permuted points

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("Dm/DMLabelPermute"))
"""
function DMLabelPermute(petsclib::PetscLibType, label::DMLabel, permutation::IS, labelNew::DMLabel) end

@for_petsc function DMLabelPermute(petsclib::$UnionPetscLib, label::DMLabel, permutation::IS, labelNew::DMLabel )

    @chk ccall(
               (:DMLabelPermute, $petsc_library),
               PetscErrorCode,
               (DMLabel, CIS, Ptr{DMLabel}),
               label, permutation, labelNew,
              )


	return nothing
end 

"""
	DMLabelPermuteValues(petsclib::PetscLibType,label::DMLabel, permutation::IS) 
Permute the values in a label

Not collective

Input Parameters:
- `label`       - the `DMLabel`
- `permutation` - the value permutation, permutation[old value] = new value

Output Parameter:
- `label` - the `DMLabel` now with permuted values

-seealso: `DMLabelRewriteValues()`, `DMLabel`, `DM`, `DMLabelPermute()`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("Dm/DMLabelPermuteValues"))
"""
function DMLabelPermuteValues(petsclib::PetscLibType, label::DMLabel, permutation::IS) end

@for_petsc function DMLabelPermuteValues(petsclib::$UnionPetscLib, label::DMLabel, permutation::IS )

    @chk ccall(
               (:DMLabelPermuteValues, $petsc_library),
               PetscErrorCode,
               (DMLabel, CIS),
               label, permutation,
              )


	return nothing
end 

"""
	DMLabelRewriteValues(petsclib::PetscLibType,label::DMLabel, permutation::IS) 
Permute the values in a label, but some may be omitted

Not collective

Input Parameters:
- `label`       - the `DMLabel`
- `permutation` - the value permutation, permutation[old value] = new value, but some maybe omitted

Output Parameter:
- `label` - the `DMLabel` now with permuted values

-seealso: `DMLabelPermuteValues()`, `DMLabel`, `DM`, `DMLabelPermute()`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("Dm/DMLabelRewriteValues"))
"""
function DMLabelRewriteValues(petsclib::PetscLibType, label::DMLabel, permutation::IS) end

@for_petsc function DMLabelRewriteValues(petsclib::$UnionPetscLib, label::DMLabel, permutation::IS )

    @chk ccall(
               (:DMLabelRewriteValues, $petsc_library),
               PetscErrorCode,
               (DMLabel, CIS),
               label, permutation,
              )


	return nothing
end 

"""
	DMLabelDistribute(petsclib::PetscLibType,label::DMLabel, sf::PetscSF, labelNew::DMLabel) 
Create a new label pushed forward over the `PetscSF`

Collective

Input Parameters:
- `label` - the `DMLabel`
- `sf`    - the map from old to new distribution

Output Parameter:
- `labelNew` - the new redistributed label

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("Dm/DMLabelDistribute"))
"""
function DMLabelDistribute(petsclib::PetscLibType, label::DMLabel, sf::PetscSF, labelNew::DMLabel) end

@for_petsc function DMLabelDistribute(petsclib::$UnionPetscLib, label::DMLabel, sf::PetscSF, labelNew::DMLabel )

    @chk ccall(
               (:DMLabelDistribute, $petsc_library),
               PetscErrorCode,
               (DMLabel, PetscSF, Ptr{DMLabel}),
               label, sf, labelNew,
              )


	return nothing
end 

"""
	DMLabelGather(petsclib::PetscLibType,label::DMLabel, sf::PetscSF, labelNew::DMLabel) 
Gather all label values from leafs into roots

Collective

Input Parameters:
- `label` - the `DMLabel`
- `sf`    - the `PetscSF` communication map

Output Parameter:
- `labelNew` - the new `DMLabel` with localised leaf values

Level: developer

-seealso: `DMLabel`, `DM`, `DMLabelDistribute()`

# External Links
$(_doc_external("Dm/DMLabelGather"))
"""
function DMLabelGather(petsclib::PetscLibType, label::DMLabel, sf::PetscSF, labelNew::DMLabel) end

@for_petsc function DMLabelGather(petsclib::$UnionPetscLib, label::DMLabel, sf::PetscSF, labelNew::DMLabel )

    @chk ccall(
               (:DMLabelGather, $petsc_library),
               PetscErrorCode,
               (DMLabel, PetscSF, Ptr{DMLabel}),
               label, sf, labelNew,
              )


	return nothing
end 

"""
	DMLabelPropagateBegin(petsclib::PetscLibType,label::DMLabel, sf::PetscSF) 
Setup a cycle of label propagation

Collective

Input Parameters:
- `label` - The `DMLabel` to propagate across processes
- `sf`    - The `PetscSF` describing parallel layout of the label points

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelPropagateEnd()`, `DMLabelPropagatePush()`

# External Links
$(_doc_external("Dm/DMLabelPropagateBegin"))
"""
function DMLabelPropagateBegin(petsclib::PetscLibType, label::DMLabel, sf::PetscSF) end

@for_petsc function DMLabelPropagateBegin(petsclib::$UnionPetscLib, label::DMLabel, sf::PetscSF )

    @chk ccall(
               (:DMLabelPropagateBegin, $petsc_library),
               PetscErrorCode,
               (DMLabel, PetscSF),
               label, sf,
              )


	return nothing
end 

"""
	DMLabelPropagateEnd(petsclib::PetscLibType,label::DMLabel, pointSF::PetscSF) 
Tear down a cycle of label propagation

Collective

Input Parameters:
- `label`   - The `DMLabel` to propagate across processes
- `pointSF` - The `PetscSF` describing parallel layout of the label points

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelPropagateBegin()`, `DMLabelPropagatePush()`

# External Links
$(_doc_external("Dm/DMLabelPropagateEnd"))
"""
function DMLabelPropagateEnd(petsclib::PetscLibType, label::DMLabel, pointSF::PetscSF) end

@for_petsc function DMLabelPropagateEnd(petsclib::$UnionPetscLib, label::DMLabel, pointSF::PetscSF )

    @chk ccall(
               (:DMLabelPropagateEnd, $petsc_library),
               PetscErrorCode,
               (DMLabel, PetscSF),
               label, pointSF,
              )


	return nothing
end 

"""
	DMLabelPropagatePush(petsclib::PetscLibType,label::DMLabel, pointSF::PetscSF, markPoint::external, ctx::Cvoid) 
Tear down a cycle of label propagation

Collective

Input Parameters:
- `label`     - The `DMLabel` to propagate across processes
- `pointSF`   - The `PetscSF` describing parallel layout of the label points
- `markPoint` - An optional callback that is called when a point is marked, or `NULL`
- `ctx`       - An optional user context for the callback, or `NULL`

Calling sequence of `markPoint`:
- `label` - The `DMLabel`
- `p`     - The point being marked
- `val`   - The label value for `p`
- `ctx`   - An optional user context

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelPropagateBegin()`, `DMLabelPropagateEnd()`

# External Links
$(_doc_external("Dm/DMLabelPropagatePush"))
"""
function DMLabelPropagatePush(petsclib::PetscLibType, label::DMLabel, pointSF::PetscSF, markPoint::external, ctx::Cvoid) end

@for_petsc function DMLabelPropagatePush(petsclib::$UnionPetscLib, label::DMLabel, pointSF::PetscSF, markPoint::external, ctx::Cvoid )

    @chk ccall(
               (:DMLabelPropagatePush, $petsc_library),
               PetscErrorCode,
               (DMLabel, PetscSF, external, Ptr{Cvoid}),
               label, pointSF, markPoint, ctx,
              )


	return nothing
end 

"""
	DMLabelConvertToSection(petsclib::PetscLibType,label::DMLabel, section::PetscSection, is::IS) 
Make a `PetscSection`/`IS` pair that encodes the label

Not Collective

Input Parameter:
- `label` - the `DMLabel`

Output Parameters:
- `section` - the section giving offsets for each stratum
- `is`      - An `IS` containing all the label points

Level: developer

-seealso: `DMLabel`, `DM`, `DMLabelDistribute()`

# External Links
$(_doc_external("Dm/DMLabelConvertToSection"))
"""
function DMLabelConvertToSection(petsclib::PetscLibType, label::DMLabel, section::PetscSection, is::IS) end

@for_petsc function DMLabelConvertToSection(petsclib::$UnionPetscLib, label::DMLabel, section::PetscSection, is::IS )
	is_ = Ref(is.ptr)

    @chk ccall(
               (:DMLabelConvertToSection, $petsc_library),
               PetscErrorCode,
               (DMLabel, Ptr{PetscSection}, Ptr{CIS}),
               label, section, is_,
              )

	is.ptr = C_NULL

	return nothing
end 

"""
	DMLabelRegisterAll(petsclib::PetscLibType) 
Registers all of the `DMLabel` implementations in the `DM` package.

Not Collective

Level: advanced

-seealso: `DMLabel`, `DM`, `DMRegisterAll()`, `DMLabelRegisterDestroy()`

# External Links
$(_doc_external("Dm/DMLabelRegisterAll"))
"""
function DMLabelRegisterAll(petsclib::PetscLibType) end

@for_petsc function DMLabelRegisterAll(petsclib::$UnionPetscLib)

    @chk ccall(
               (:DMLabelRegisterAll, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	DMLabelRegisterDestroy(petsclib::PetscLibType) 
This function destroys the `DMLabel` registry. It is called from `PetscFinalize()`.

Level: developer

-seealso: `DMLabel`, `DM`, `PetscInitialize()`

# External Links
$(_doc_external("Dm/DMLabelRegisterDestroy"))
"""
function DMLabelRegisterDestroy(petsclib::PetscLibType) end

@for_petsc function DMLabelRegisterDestroy(petsclib::$UnionPetscLib)

    @chk ccall(
               (:DMLabelRegisterDestroy, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	DMLabelSetType(petsclib::PetscLibType,label::DMLabel, method::DMLabelType) 
Sets the particular implementation for a label.

Collective

Input Parameters:
- `label`  - The label
- `method` - The name of the label type

Options Database Key:
- `-dm_label_type <type>` - Sets the label type; use -help for a list of available types or see `DMLabelType`

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelGetType()`, `DMLabelCreate()`

# External Links
$(_doc_external("Dm/DMLabelSetType"))
"""
function DMLabelSetType(petsclib::PetscLibType, label::DMLabel, method::DMLabelType) end

@for_petsc function DMLabelSetType(petsclib::$UnionPetscLib, label::DMLabel, method::DMLabelType )

    @chk ccall(
               (:DMLabelSetType, $petsc_library),
               PetscErrorCode,
               (DMLabel, DMLabelType),
               label, method,
              )


	return nothing
end 

"""
	type::DMLabelType = DMLabelGetType(petsclib::PetscLibType,label::DMLabel) 
Gets the type name (as a string) from the label.

Not Collective

Input Parameter:
- `label` - The `DMLabel`

Output Parameter:
- `type` - The `DMLabel` type name

Level: intermediate

-seealso: `DMLabel`, `DM`, `DMLabelSetType()`, `DMLabelCreate()`

# External Links
$(_doc_external("Dm/DMLabelGetType"))
"""
function DMLabelGetType(petsclib::PetscLibType, label::DMLabel) end

@for_petsc function DMLabelGetType(petsclib::$UnionPetscLib, label::DMLabel )
	type_ = Ref{DMLabelType}()

    @chk ccall(
               (:DMLabelGetType, $petsc_library),
               PetscErrorCode,
               (DMLabel, Ptr{DMLabelType}),
               label, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	DMLabelEphemeralGetLabel(petsclib::PetscLibType,label::DMLabel, olabel::DMLabel) 
Get the base label for this ephemeral label

Not Collective

Input Parameter:
- `label` - the `DMLabel`

Output Parameter:
- `olabel` - the base label for this ephemeral label

Level: intermediate

-seealso: `DMLabelEphemeralSetLabel()`, `DMLabelEphemeralGetTransform()`, `DMLabelSetType()`

# External Links
$(_doc_external("Dm/DMLabelEphemeralGetLabel"))
"""
function DMLabelEphemeralGetLabel(petsclib::PetscLibType, label::DMLabel, olabel::DMLabel) end

@for_petsc function DMLabelEphemeralGetLabel(petsclib::$UnionPetscLib, label::DMLabel, olabel::DMLabel )

    @chk ccall(
               (:DMLabelEphemeralGetLabel, $petsc_library),
               PetscErrorCode,
               (DMLabel, Ptr{DMLabel}),
               label, olabel,
              )


	return nothing
end 

"""
	DMLabelEphemeralSetLabel(petsclib::PetscLibType,label::DMLabel, olabel::DMLabel) 
Set the base label for this ephemeral label

Not Collective

Input Parameters:
- `label`  - the `DMLabel`
- `olabel` - the base label for this ephemeral label

Level: intermediate

-seealso: `DMLabelEphemeralGetLabel()`, `DMLabelEphemeralSetTransform()`, `DMLabelSetType()`

# External Links
$(_doc_external("Dm/DMLabelEphemeralSetLabel"))
"""
function DMLabelEphemeralSetLabel(petsclib::PetscLibType, label::DMLabel, olabel::DMLabel) end

@for_petsc function DMLabelEphemeralSetLabel(petsclib::$UnionPetscLib, label::DMLabel, olabel::DMLabel )

    @chk ccall(
               (:DMLabelEphemeralSetLabel, $petsc_library),
               PetscErrorCode,
               (DMLabel, DMLabel),
               label, olabel,
              )


	return nothing
end 

"""
	DMLabelEphemeralGetTransform(petsclib::PetscLibType,label::DMLabel, tr::DMPlexTransform) 
Get the transform for this ephemeral label

Not Collective

Input Parameter:
- `label` - the `DMLabel`

Output Parameter:
- `tr` - the transform for this ephemeral label

Level: intermediate

-seealso: `DMLabelEphemeralSetTransform()`, `DMLabelEphemeralGetLabel()`, `DMLabelSetType()`

# External Links
$(_doc_external("Dm/DMLabelEphemeralGetTransform"))
"""
function DMLabelEphemeralGetTransform(petsclib::PetscLibType, label::DMLabel, tr::DMPlexTransform) end

@for_petsc function DMLabelEphemeralGetTransform(petsclib::$UnionPetscLib, label::DMLabel, tr::DMPlexTransform )

    @chk ccall(
               (:DMLabelEphemeralGetTransform, $petsc_library),
               PetscErrorCode,
               (DMLabel, Ptr{DMPlexTransform}),
               label, tr,
              )


	return nothing
end 

"""
	DMLabelEphemeralSetTransform(petsclib::PetscLibType,label::DMLabel, tr::DMPlexTransform) 
Set the transform for this ephemeral label

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `tr`    - the transform for this ephemeral label

Level: intermediate

-seealso: `DMLabelEphemeralGetTransform()`, `DMLabelEphemeralSetLabel()`, `DMLabelSetType()`

# External Links
$(_doc_external("Dm/DMLabelEphemeralSetTransform"))
"""
function DMLabelEphemeralSetTransform(petsclib::PetscLibType, label::DMLabel, tr::DMPlexTransform) end

@for_petsc function DMLabelEphemeralSetTransform(petsclib::$UnionPetscLib, label::DMLabel, tr::DMPlexTransform )

    @chk ccall(
               (:DMLabelEphemeralSetTransform, $petsc_library),
               PetscErrorCode,
               (DMLabel, DMPlexTransform),
               label, tr,
              )


	return nothing
end 

