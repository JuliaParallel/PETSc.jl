"""
	handler::PetscLogHandler = PetscLogHandlerCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Create a log handler for profiling events and stages.  PETSc
provides several implementations of `PetscLogHandler` that interface to different ways to
summarize or visualize profiling data: see `PetscLogHandlerType` for a list.

Collective

Input Parameter:
- `comm` - the communicator for synchronizing and viewing events with this handler

Output Parameter:
- `handler` - the `PetscLogHandler`

Level: developer

-seealso: [](ch_profiling), `PetscLogHandler`, `PetscLogHandlerSetType()`, `PetscLogHandlerStart()`, `PetscLogHandlerStop()`

# External Links
$(_doc_external("Log/PetscLogHandlerCreate"))
"""
function PetscLogHandlerCreate(petsclib::PetscLibType, comm::MPI_Comm)
    error("PetscLogHandlerCreate: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	handler_ = Ref{PetscLogHandler}()

    @chk ccall(
               (:PetscLogHandlerCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscLogHandler}),
               comm, handler_,
              )

	handler = handler_[]

	return handler
end 

"""
	handler::PetscLogHandler = PetscLogHandlerCreateLegacy(petsclib::PetscLibType,comm::MPI_Comm, PetscLogPLB::external, PetscLogPLE::external, PetscLogPHC::external, PetscLogPHD::external) 
Create a `PetscLogHandler` from callbacks matching PETSc's legacy log handler callbacks

Collective

Input Parameters:
- `comm`        - an MPI communicator
- `PetscLogPLB` - a function to call during `PetscLogHandlerEventBegin()` (or `NULL`)
- `PetscLogPLE` - a function to call during `PetscLogHandlerEventEnd()` (or `NULL`)
- `PetscLogPHC` - a function to call during `PetscLogHandlerObjectCreate()` (or `NULL`)
- `PetscLogPHD` - a function to call during `PetscLogHandlerObjectDestroy()` (or `NULL`)

Output Parameter:
- `handler` - a `PetscLogHandler`

Calling sequence of `PetscLogPLB`:
- `e`  - a `PetscLogEvent` that is beginning
- `_i` - deprecated, unused
- `o1` - a `PetscObject` associated with `e` (or `NULL`)
- `o2` - a `PetscObject` associated with `e` (or `NULL`)
- `o3` - a `PetscObject` associated with `e` (or `NULL`)
- `o4` - a `PetscObject` associated with `e` (or `NULL`)

Calling sequence of `PetscLogPLE`:
- `e`  - a `PetscLogEvent` that is beginning
- `_i` - deprecated, unused
- `o1` - a `PetscObject` associated with `e` (or `NULL`)
- `o2` - a `PetscObject` associated with `e` (or `NULL`)
- `o3` - a `PetscObject` associated with `e` (or `NULL`)
- `o4` - a `PetscObject` associated with `e` (or `NULL`)

Calling sequence of `PetscLogPHC`:
- `o` - a `PetscObject` that has just been created

Calling sequence of `PetscLogPHD`:
- `o` - a `PetscObject` that is about to be destroyed

Level: developer

-seealso: [](ch_profiling)

# External Links
$(_doc_external("Log/PetscLogHandlerCreateLegacy"))
"""
function PetscLogHandlerCreateLegacy(petsclib::PetscLibType, comm::MPI_Comm, PetscLogPLB::external, PetscLogPLE::external, PetscLogPHC::external, PetscLogPHD::external)
    error("PetscLogHandlerCreateLegacy: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerCreateLegacy(petsclib::$UnionPetscLib, comm::MPI_Comm, PetscLogPLB::external, PetscLogPLE::external, PetscLogPHC::external, PetscLogPHD::external )
	handler_ = Ref{PetscLogHandler}()

    @chk ccall(
               (:PetscLogHandlerCreateLegacy, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, external, external, external, external, Ptr{PetscLogHandler}),
               comm, PetscLogPLB, PetscLogPLE, PetscLogPHC, PetscLogPHD, handler_,
              )

	handler = handler_[]

	return handler
end 

"""
	file::Libc.FILE,handler::PetscLogHandler = PetscLogHandlerCreateTrace(petsclib::PetscLibType,comm::MPI_Comm) 
Create a logger that traces events and stages to a given file descriptor

Collective, No Fortran Support

Input Parameters:
- `comm` - an MPI communicator
- `file` - a file descriptor

Output Parameters:
- `handler` - a `PetscLogHandler of type `PETSCLOGHANDLERTRACE`

Level: developer

-seealso: [](ch_profiling), `PetscLogHandler`, `PetscLogHandlerTraceBegin()`

# External Links
$(_doc_external("Log/PetscLogHandlerCreateTrace"))
"""
function PetscLogHandlerCreateTrace(petsclib::PetscLibType, comm::MPI_Comm)
    error("PetscLogHandlerCreateTrace: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerCreateTrace(petsclib::$UnionPetscLib, comm::MPI_Comm )
	file_ = Ref{Libc.FILE}()
	handler_ = Ref{PetscLogHandler}()

    @chk ccall(
               (:PetscLogHandlerCreateTrace, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Libc.FILE}, Ptr{PetscLogHandler}),
               comm, file_, handler_,
              )

	file = file_[]
	handler = handler_[]

	return file,handler
end 

"""
	PetscLogHandlerDestroy(petsclib::PetscLibType,handler::Union{PetscLogHandler, Ref{PetscLogHandler}}) 
Destroy a `PetscLogHandler`

Logically collective

Input Parameter:
- `handler` - handler to be destroyed

Level: developer

-seealso: [](ch_profiling), `PetscLogHandler`, `PetscLogHandlerCreate()`

# External Links
$(_doc_external("Log/PetscLogHandlerDestroy"))
"""
function PetscLogHandlerDestroy(petsclib::PetscLibType, handler::Union{PetscLogHandler, Ref{PetscLogHandler}})
    error("PetscLogHandlerDestroy: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerDestroy(petsclib::$UnionPetscLib, handler::Union{PetscLogHandler, Ref{PetscLogHandler}} )
	handler_ = handler isa Base.RefValue ? handler : Ref{PetscLogHandler}(handler)

    @chk ccall(
               (:PetscLogHandlerDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscLogHandler},),
               handler_,
              )


	return nothing
end 

"""
	PetscLogHandlerDump(petsclib::PetscLibType,handler::PetscLogHandler, sname::String) 
Dump the records of a log handler to file

Not collective

Input Parameters:
- `handler` - a `PetscLogHandler`
- `sname`   - the name of the file to dump log data to

Level: developer

-seealso: [](ch_profiling)

# External Links
$(_doc_external("Log/PetscLogHandlerDump"))
"""
function PetscLogHandlerDump(petsclib::PetscLibType, handler::PetscLogHandler, sname::String)
    error("PetscLogHandlerDump: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerDump(petsclib::$UnionPetscLib, handler::PetscLogHandler, sname::String )

    @chk ccall(
               (:PetscLogHandlerDump, $petsc_library),
               PetscErrorCode,
               (PetscLogHandler, Ptr{Cchar}),
               handler, sname,
              )


	return nothing
end 

"""
	PetscLogHandlerEventBegin(petsclib::PetscLibType,h::PetscLogHandler, e::PetscLogEvent, o1, o2, o3, o4) 
Record the beginning of an event in a log handler

Not collective

Input Parameters:
- `h`  - the `PetscLogHandler`
- `e`  - a registered `PetscLogEvent`
- `o1` - `PetscObject` associated with the event (may be `NULL`)
- `o2` - `PetscObject` associated with the event (may be `NULL`)
- `o3` - `PetscObject` associated with the event (may be `NULL`)
- `o4` - `PetscObject` associated with the event (may be `NULL`)

Level: developer

-seealso: [](ch_profiling), `PetscLogHandler`, `PetscLogEventBegin()`, `PetscLogEventEnd()`, `PetscLogEventSync()`, `PetscLogHandlerEventEnd()`, `PetscLogHandlerEventSync()`

# External Links
$(_doc_external("Log/PetscLogHandlerEventBegin"))
"""
function PetscLogHandlerEventBegin(petsclib::PetscLibType, h::PetscLogHandler, e::PetscLogEvent, o1, o2, o3, o4)
    error("PetscLogHandlerEventBegin: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerEventBegin(petsclib::$UnionPetscLib, h::PetscLogHandler, e::PetscLogEvent, o1, o2, o3, o4 )

    @chk ccall(
               (:PetscLogHandlerEventBegin, $petsc_library),
               PetscErrorCode,
               (PetscLogHandler, PetscLogEvent, PetscObject, PetscObject, PetscObject, PetscObject),
               h, e, o1, o2, o3, o4,
              )


	return nothing
end 

"""
	PetscLogHandlerEventDeactivatePop(petsclib::PetscLibType,handler::PetscLogHandler, stage::PetscLogStage, event::PetscLogEvent) 
Undo temporary deactivation a logging event for a log handler

Not collective

Input Parameters:
- `handler` - a `PetscLogHandler`
- `stage`   - a `PetscLogStage` (or `PETSC_DEFAULT` for the current stage)
- `event`   - a `PetscLogEvent`

Level: developer

-seealso: [](ch_profiling), `PetscLogHandlerEventDeactivatePush()`

# External Links
$(_doc_external("Log/PetscLogHandlerEventDeactivatePop"))
"""
function PetscLogHandlerEventDeactivatePop(petsclib::PetscLibType, handler::PetscLogHandler, stage::PetscLogStage, event::PetscLogEvent)
    error("PetscLogHandlerEventDeactivatePop: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerEventDeactivatePop(petsclib::$UnionPetscLib, handler::PetscLogHandler, stage::PetscLogStage, event::PetscLogEvent )

    @chk ccall(
               (:PetscLogHandlerEventDeactivatePop, $petsc_library),
               PetscErrorCode,
               (PetscLogHandler, PetscLogStage, PetscLogEvent),
               handler, stage, event,
              )


	return nothing
end 

"""
	PetscLogHandlerEventDeactivatePush(petsclib::PetscLibType,handler::PetscLogHandler, stage::PetscLogStage, event::PetscLogEvent) 
Temporarily deactivate a logging event for a log handler

Not collective

Input Parameters:
- `handler` - a `PetscLogHandler`
- `stage`   - a `PetscLogStage` (or `PETSC_DEFAULT` for the current stage)
- `event`   - a `PetscLogEvent`

Level: developer

-seealso: [](ch_profiling), `PetscLogHandlerEventDeactivatePop()`

# External Links
$(_doc_external("Log/PetscLogHandlerEventDeactivatePush"))
"""
function PetscLogHandlerEventDeactivatePush(petsclib::PetscLibType, handler::PetscLogHandler, stage::PetscLogStage, event::PetscLogEvent)
    error("PetscLogHandlerEventDeactivatePush: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerEventDeactivatePush(petsclib::$UnionPetscLib, handler::PetscLogHandler, stage::PetscLogStage, event::PetscLogEvent )

    @chk ccall(
               (:PetscLogHandlerEventDeactivatePush, $petsc_library),
               PetscErrorCode,
               (PetscLogHandler, PetscLogStage, PetscLogEvent),
               handler, stage, event,
              )


	return nothing
end 

"""
	PetscLogHandlerEventEnd(petsclib::PetscLibType,h::PetscLogHandler, e::PetscLogEvent, o1, o2, o3, o4) 
Record the end of an event in a log handler

Not collective

Input Parameters:
- `h`  - the `PetscLogHandler`
- `e`  - a registered `PetscLogEvent`
- `o1` - `PetscObject` associated with the event (may be `NULL`)
- `o2` - `PetscObject` associated with the event (may be `NULL`)
- `o3` - `PetscObject` associated with the event (may be `NULL`)
- `o4` - `PetscObject` associated with the event (may be `NULL`)

Level: developer

-seealso: [](ch_profiling), `PetscLogHandler`, `PetscLogEventBegin()`, `PetscLogEventEnd()`, `PetscLogEventSync()`, `PetscLogHandlerEventBegin()`, `PetscLogHandlerEventSync()`

# External Links
$(_doc_external("Log/PetscLogHandlerEventEnd"))
"""
function PetscLogHandlerEventEnd(petsclib::PetscLibType, h::PetscLogHandler, e::PetscLogEvent, o1, o2, o3, o4)
    error("PetscLogHandlerEventEnd: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerEventEnd(petsclib::$UnionPetscLib, h::PetscLogHandler, e::PetscLogEvent, o1, o2, o3, o4 )

    @chk ccall(
               (:PetscLogHandlerEventEnd, $petsc_library),
               PetscErrorCode,
               (PetscLogHandler, PetscLogEvent, PetscObject, PetscObject, PetscObject, PetscObject),
               h, e, o1, o2, o3, o4,
              )


	return nothing
end 

"""
	PetscLogHandlerEventSync(petsclib::PetscLibType,h::PetscLogHandler, e::PetscLogEvent, comm::MPI_Comm) 
Synchronize a logging event

Collective

Input Parameters:
- `h`    - the `PetscLogHandler`
- `e`    - a registered `PetscLogEvent`
- `comm` - the communicator over which to synchronize `e`

Level: developer

-seealso: [](ch_profiling), `PetscLogHandler`, `PetscLogEventBegin()`, `PetscLogEventEnd()`, `PetscLogEventSync()`, `PetscLogHandlerEventBegin()`, `PetscLogHandlerEventEnd()`

# External Links
$(_doc_external("Log/PetscLogHandlerEventSync"))
"""
function PetscLogHandlerEventSync(petsclib::PetscLibType, h::PetscLogHandler, e::PetscLogEvent, comm::MPI_Comm)
    error("PetscLogHandlerEventSync: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerEventSync(petsclib::$UnionPetscLib, h::PetscLogHandler, e::PetscLogEvent, comm::MPI_Comm )

    @chk ccall(
               (:PetscLogHandlerEventSync, $petsc_library),
               PetscErrorCode,
               (PetscLogHandler, PetscLogEvent, MPI_Comm),
               h, e, comm,
              )


	return nothing
end 

"""
	PetscLogHandlerEventsPause(petsclib::PetscLibType,handler::PetscLogHandler) 
Put event logging into "paused" mode (see `PetscLogEventsPause()` for details.) for a log handler

Not collective

Input Parameter:
- `handler` - a `PetscLogHandler`

Level: developer

-seealso: [](ch_profiling), `PetscLogHandlerEventsResume()`

# External Links
$(_doc_external("Log/PetscLogHandlerEventsPause"))
"""
function PetscLogHandlerEventsPause(petsclib::PetscLibType, handler::PetscLogHandler)
    error("PetscLogHandlerEventsPause: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerEventsPause(petsclib::$UnionPetscLib, handler::PetscLogHandler )

    @chk ccall(
               (:PetscLogHandlerEventsPause, $petsc_library),
               PetscErrorCode,
               (PetscLogHandler,),
               handler,
              )


	return nothing
end 

"""
	PetscLogHandlerEventsResume(petsclib::PetscLibType,handler::PetscLogHandler) 
Resume event logging that had been put into "paused" mode (see `PetscLogEventsPause()` for details.) for a log handler

Not collective

Input Parameter:
- `handler` - a `PetscLogHandler`

Level: developer

-seealso: [](ch_profiling), `PetscLogHandlerEventsPause()`

# External Links
$(_doc_external("Log/PetscLogHandlerEventsResume"))
"""
function PetscLogHandlerEventsResume(petsclib::PetscLibType, handler::PetscLogHandler)
    error("PetscLogHandlerEventsResume: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerEventsResume(petsclib::$UnionPetscLib, handler::PetscLogHandler )

    @chk ccall(
               (:PetscLogHandlerEventsResume, $petsc_library),
               PetscErrorCode,
               (PetscLogHandler,),
               handler,
              )


	return nothing
end 

"""
	event_info::Ptr{PetscEventPerfInfo} = PetscLogHandlerGetEventPerfInfo(petsclib::PetscLibType,handler::PetscLogHandler, stage::PetscLogStage, event::PetscLogEvent) 
Get a direct reference to the `PetscEventPerfInfo` of a stage and event

Not collective, No Fortran Support

Input Parameters:
- `handler` - a `PetscLogHandler`
- `stage`   - a `PetscLogStage` (or `PETSC_DEFAULT` for the current stage)
- `event`   - a `PetscLogEvent`

Output Parameter:
- `event_info` - a pointer to a performance log for `event` during `stage` (or `NULL` if this handler does not use
`PetscEventPerfInfo` to record performance data); writing to `event_info` will change the record in
`handler`

Level: developer

-seealso: [](ch_profiling), `PetscLogEventGetPerfInfo()`, `PETSCLOGHANDLERDEFAULT`

# External Links
$(_doc_external("Log/PetscLogHandlerGetEventPerfInfo"))
"""
function PetscLogHandlerGetEventPerfInfo(petsclib::PetscLibType, handler::PetscLogHandler, stage::PetscLogStage, event::PetscLogEvent)
    error("PetscLogHandlerGetEventPerfInfo: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerGetEventPerfInfo(petsclib::$UnionPetscLib, handler::PetscLogHandler, stage::PetscLogStage, event::PetscLogEvent )
	event_info_ = Ref{Ptr{PetscEventPerfInfo}}()

    @chk ccall(
               (:PetscLogHandlerGetEventPerfInfo, $petsc_library),
               PetscErrorCode,
               (PetscLogHandler, PetscLogStage, PetscLogEvent, Ptr{Ptr{PetscEventPerfInfo}}),
               handler, stage, event, event_info_,
              )

	event_info = event_info_[]

	return event_info
end 

"""
	num_objects::PetscInt = PetscLogHandlerGetNumObjects(petsclib::PetscLibType,handler::PetscLogHandler) 
Get the number of objects that were logged with a log handler

Not Collective

Input Parameter:
- `handler` - a `PetscLogHandler`

Output Parameter:
- `num_objects` - the number of objects whose creations and destructions were logged with `handler`
(`PetscLogHandlerObjectCreate()` / `PetscLogHandlerObjectDestroy()`), or -1
if the handler does not keep track of this number.

Level: developer

-seealso: [](ch_profiling)

# External Links
$(_doc_external("Log/PetscLogHandlerGetNumObjects"))
"""
function PetscLogHandlerGetNumObjects(petsclib::PetscLibType, handler::PetscLogHandler)
    error("PetscLogHandlerGetNumObjects: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerGetNumObjects(petsclib::$UnionPetscLib, handler::PetscLogHandler )
	num_objects_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscLogHandlerGetNumObjects, $petsc_library),
               PetscErrorCode,
               (PetscLogHandler, Ptr{$PetscInt}),
               handler, num_objects_,
              )

	num_objects = num_objects_[]

	return num_objects
end 

"""
	stage_info::Ptr{PetscEventPerfInfo} = PetscLogHandlerGetStagePerfInfo(petsclib::PetscLibType,handler::PetscLogHandler, stage::PetscLogStage) 
Get a direct reference to the `PetscEventPerfInfo` of a stage

Not collective, No Fortran Support

Input Parameters:
- `handler` - a `PetscLogHandler`
- `stage`   - a `PetscLogStage` (or `PETSC_DEFAULT` for the current stage)

Output Parameter:
- `stage_info` - a pointer to a performance log for `stage` (or `NULL` if this handler does not use `PetscEventPerfInfo`
to record performance data); writing to `stage_info` will change the record in `handler`

Level: developer

-seealso: [](ch_profiling), `PetscLogEventGetPerfInfo()`, `PETSCLOGHANDLERDEFAULT`

# External Links
$(_doc_external("Log/PetscLogHandlerGetStagePerfInfo"))
"""
function PetscLogHandlerGetStagePerfInfo(petsclib::PetscLibType, handler::PetscLogHandler, stage::PetscLogStage)
    error("PetscLogHandlerGetStagePerfInfo: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerGetStagePerfInfo(petsclib::$UnionPetscLib, handler::PetscLogHandler, stage::PetscLogStage )
	stage_info_ = Ref{Ptr{PetscEventPerfInfo}}()

    @chk ccall(
               (:PetscLogHandlerGetStagePerfInfo, $petsc_library),
               PetscErrorCode,
               (PetscLogHandler, PetscLogStage, Ptr{Ptr{PetscEventPerfInfo}}),
               handler, stage, stage_info_,
              )

	stage_info = stage_info_[]

	return stage_info
end 

"""
	state::PetscLogState = PetscLogHandlerGetState(petsclib::PetscLibType,h::PetscLogHandler) 
Get the logging state that provides the stream of events and stages for a log handler.

Logically collective

Input Parameter:
- `h` - the `PetscLogHandler`

Output Parameter:
- `state` - the `PetscLogState`

Level: developer

-seealso: [](ch_profiling), `PetscLogHandler`, `PetscLogState`, `PetscLogEventBegin()`, `PetscLogHandlerStart()`

# External Links
$(_doc_external("Log/PetscLogHandlerGetState"))
"""
function PetscLogHandlerGetState(petsclib::PetscLibType, h::PetscLogHandler)
    error("PetscLogHandlerGetState: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerGetState(petsclib::$UnionPetscLib, h::PetscLogHandler )
	state_ = Ref{PetscLogState}()

    @chk ccall(
               (:PetscLogHandlerGetState, $petsc_library),
               PetscErrorCode,
               (PetscLogHandler, Ptr{PetscLogState}),
               h, state_,
              )

	state = state_[]

	return state
end 

"""
	name::PetscLogHandlerType = PetscLogHandlerGetType(petsclib::PetscLibType,handler::PetscLogHandler) 
Gets the `PetscLoagHandlerType` (as a string) from the `PetscLogHandler` object.

Not collective

Input Parameter:
- `handler` - the `PetscLogHandler`

Output Parameter:
- `name` - The `PetscLogHandlerType` name

Level: developer

-seealso: [](ch_profiling), `PetscLogHandler`, `PetscLogHandlerCreate()`, `PetscLogHandlerRegister()`, `PetscLogHandlerSetType()`

# External Links
$(_doc_external("Log/PetscLogHandlerGetType"))
"""
function PetscLogHandlerGetType(petsclib::PetscLibType, handler::PetscLogHandler)
    error("PetscLogHandlerGetType: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerGetType(petsclib::$UnionPetscLib, handler::PetscLogHandler )
	name_ = Ref{PetscLogHandlerType}()

    @chk ccall(
               (:PetscLogHandlerGetType, $petsc_library),
               PetscErrorCode,
               (PetscLogHandler, Ptr{PetscLogHandlerType}),
               handler, name_,
              )

	name = name_[] == C_NULL ? "" : unsafe_string(name_[])

	return name
end 

"""
	PetscLogHandlerObjectCreate(petsclib::PetscLibType,h::PetscLogHandler, obj) 
Record the creation of an object in a log handler.

Not collective

Input Parameters:
- `h`   - the `PetscLogHandler`
- `obj` - a newly created `PetscObject`

Level: developer

-seealso: [](ch_profiling), `PetscLogHandler`, `PetscLogObjectCreate()`, `PetscLogObjectDestroy()`, `PetscLogHandlerObjectDestroy()`

# External Links
$(_doc_external("Log/PetscLogHandlerObjectCreate"))
"""
function PetscLogHandlerObjectCreate(petsclib::PetscLibType, h::PetscLogHandler, obj)
    error("PetscLogHandlerObjectCreate: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerObjectCreate(petsclib::$UnionPetscLib, h::PetscLogHandler, obj )

    @chk ccall(
               (:PetscLogHandlerObjectCreate, $petsc_library),
               PetscErrorCode,
               (PetscLogHandler, PetscObject),
               h, obj,
              )


	return nothing
end 

"""
	PetscLogHandlerObjectDestroy(petsclib::PetscLibType,h::PetscLogHandler, obj) 
Record the destruction of an object in a log handler.

Not collective

Input Parameters:
- `h`   - the `PetscLogHandler`
- `obj` - a newly created `PetscObject`

Level: developer

-seealso: [](ch_profiling), `PetscLogHandler`, `PetscLogObjectCreate()`, `PetscLogObjectDestroy()`, `PetscLogHandlerObjectCreate()`

# External Links
$(_doc_external("Log/PetscLogHandlerObjectDestroy"))
"""
function PetscLogHandlerObjectDestroy(petsclib::PetscLibType, h::PetscLogHandler, obj)
    error("PetscLogHandlerObjectDestroy: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerObjectDestroy(petsclib::$UnionPetscLib, h::PetscLogHandler, obj )

    @chk ccall(
               (:PetscLogHandlerObjectDestroy, $petsc_library),
               PetscErrorCode,
               (PetscLogHandler, PetscObject),
               h, obj,
              )


	return nothing
end 

"""
	PetscLogHandlerRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Register a new `PetscLogHandler`

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - The name of a new user-defined creation routine
- `function` - The creation routine

-seealso: [](ch_profiling), `PetscLogHandler`, `PetscLogHandlerCreate()`, `PetscLogHandlerSetType()`, `PetscLogHandlerGetType()`

# External Links
$(_doc_external("Log/PetscLogHandlerRegister"))
"""
function PetscLogHandlerRegister(petsclib::PetscLibType, sname::String, fnc::external)
    error("PetscLogHandlerRegister: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:PetscLogHandlerRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	PetscLogHandlerSetLogActions(petsclib::PetscLibType,handler::PetscLogHandler, flag::PetscBool) 
Determines whether actions are logged for a log handler.

Not Collective

Input Parameters:
- `handler` - a `PetscLogHandler`
- `flag`    - `PETSC_TRUE` if actions are to be logged (ignored if `handler` does not log actions)

Level: developer

-seealso: [](ch_profiling), `PetscLogSetLogActions()`, `PetscLogStagePush()`, `PetscLogStagePop()`, `PetscLogGetDefaultHandler()`

# External Links
$(_doc_external("Log/PetscLogHandlerSetLogActions"))
"""
function PetscLogHandlerSetLogActions(petsclib::PetscLibType, handler::PetscLogHandler, flag::PetscBool)
    error("PetscLogHandlerSetLogActions: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerSetLogActions(petsclib::$UnionPetscLib, handler::PetscLogHandler, flag::PetscBool )

    @chk ccall(
               (:PetscLogHandlerSetLogActions, $petsc_library),
               PetscErrorCode,
               (PetscLogHandler, PetscBool),
               handler, flag,
              )


	return nothing
end 

"""
	PetscLogHandlerSetLogObjects(petsclib::PetscLibType,handler::PetscLogHandler, flag::PetscBool) 
Determines whether objects are logged for a log handler.

Not Collective

Input Parameters:
- `handler` - a `PetscLogHandler`
- `flag`    - `PETSC_TRUE` if objects are to be logged (ignored if `handler` does not log objects)

Level: developer

-seealso: [](ch_profiling), `PetscLogSetLogObjects()`, `PetscLogStagePush()`, `PetscLogStagePop()`, `PetscLogGetDefaultHandler()`

# External Links
$(_doc_external("Log/PetscLogHandlerSetLogObjects"))
"""
function PetscLogHandlerSetLogObjects(petsclib::PetscLibType, handler::PetscLogHandler, flag::PetscBool)
    error("PetscLogHandlerSetLogObjects: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerSetLogObjects(petsclib::$UnionPetscLib, handler::PetscLogHandler, flag::PetscBool )

    @chk ccall(
               (:PetscLogHandlerSetLogObjects, $petsc_library),
               PetscErrorCode,
               (PetscLogHandler, PetscBool),
               handler, flag,
              )


	return nothing
end 

"""
	PetscLogHandlerSetState(petsclib::PetscLibType,h::PetscLogHandler, state::PetscLogState) 
Set the logging state that provides the stream of events and stages for a log handler.

Logically collective

Input Parameters:
- `h`     - the `PetscLogHandler`
- `state` - the `PetscLogState`

Level: developer

-seealso: [](ch_profiling), `PetscLogHandler`, `PetscLogState`, `PetscLogEventBegin()`, `PetscLogHandlerStart()`

# External Links
$(_doc_external("Log/PetscLogHandlerSetState"))
"""
function PetscLogHandlerSetState(petsclib::PetscLibType, h::PetscLogHandler, state::PetscLogState)
    error("PetscLogHandlerSetState: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerSetState(petsclib::$UnionPetscLib, h::PetscLogHandler, state::PetscLogState )

    @chk ccall(
               (:PetscLogHandlerSetState, $petsc_library),
               PetscErrorCode,
               (PetscLogHandler, PetscLogState),
               h, state,
              )


	return nothing
end 

"""
	PetscLogHandlerSetType(petsclib::PetscLibType,handler::PetscLogHandler, name::PetscLogHandlerType) 
Set the type of a `PetscLogHandler`

Input Parameters:
- `handler` - the `PetscLogHandler`
- `name`    - The kind of log handler

Level: developer

-seealso: [](ch_profiling), `PetscLogHandler`, `PetscLogHandlerCreate()`, `PetscLogHandlerRegister()`, `PetscLogHandlerGetType()`

# External Links
$(_doc_external("Log/PetscLogHandlerSetType"))
"""
function PetscLogHandlerSetType(petsclib::PetscLibType, handler::PetscLogHandler, name::PetscLogHandlerType)
    error("PetscLogHandlerSetType: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerSetType(petsclib::$UnionPetscLib, handler::PetscLogHandler, name::PetscLogHandlerType )

    @chk ccall(
               (:PetscLogHandlerSetType, $petsc_library),
               PetscErrorCode,
               (PetscLogHandler, PetscLogHandlerType),
               handler, name,
              )


	return nothing
end 

"""
	isVisible::PetscBool = PetscLogHandlerStageGetVisible(petsclib::PetscLibType,handler::PetscLogHandler, stage::PetscLogStage) 
Get the visibility of logging stage in `PetscLogHandlerView()` for a log handler

Not collective

Input Parameters:
- `handler` - a `PetscLogHandler`
- `stage`   - a `PetscLogStage`

Output Parameter:
- `isVisible` - the visibility flag, `PETSC_TRUE` to print, else `PETSC_FALSE` (defaults to `PETSC_TRUE`)

Level: developer

-seealso: [](ch_profiling), `PetscLogHandlerStageSetVisible()`

# External Links
$(_doc_external("Log/PetscLogHandlerStageGetVisible"))
"""
function PetscLogHandlerStageGetVisible(petsclib::PetscLibType, handler::PetscLogHandler, stage::PetscLogStage)
    error("PetscLogHandlerStageGetVisible: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerStageGetVisible(petsclib::$UnionPetscLib, handler::PetscLogHandler, stage::PetscLogStage )
	isVisible_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscLogHandlerStageGetVisible, $petsc_library),
               PetscErrorCode,
               (PetscLogHandler, PetscLogStage, Ptr{PetscBool}),
               handler, stage, isVisible_,
              )

	isVisible = isVisible_[]

	return isVisible
end 

"""
	PetscLogHandlerStagePop(petsclib::PetscLibType,h::PetscLogHandler, stage::PetscLogStage) 
End the current logging stage in a log handler.

Not collective

Input Parameters:
- `h`     - the `PetscLogHandler`
- `stage` - a registered `PetscLogStage`

Level: developer

-seealso: [](ch_profiling), `PetscLogHandler`, `PetscLogStagePush()`, `PetscLogStagePop()`, `PetscLogHandlerStagePush()`

# External Links
$(_doc_external("Log/PetscLogHandlerStagePop"))
"""
function PetscLogHandlerStagePop(petsclib::PetscLibType, h::PetscLogHandler, stage::PetscLogStage)
    error("PetscLogHandlerStagePop: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerStagePop(petsclib::$UnionPetscLib, h::PetscLogHandler, stage::PetscLogStage )

    @chk ccall(
               (:PetscLogHandlerStagePop, $petsc_library),
               PetscErrorCode,
               (PetscLogHandler, PetscLogStage),
               h, stage,
              )


	return nothing
end 

"""
	PetscLogHandlerStagePush(petsclib::PetscLibType,h::PetscLogHandler, stage::PetscLogStage) 
Begin a new logging stage in a log handler.

Not collective

Input Parameters:
- `h`     - the `PetscLogHandler`
- `stage` - a registered `PetscLogStage`

Level: developer

-seealso: [](ch_profiling), `PetscLogHandler`, `PetscLogStagePush()`, `PetscLogStagePop()`, `PetscLogHandlerStagePop()`

# External Links
$(_doc_external("Log/PetscLogHandlerStagePush"))
"""
function PetscLogHandlerStagePush(petsclib::PetscLibType, h::PetscLogHandler, stage::PetscLogStage)
    error("PetscLogHandlerStagePush: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerStagePush(petsclib::$UnionPetscLib, h::PetscLogHandler, stage::PetscLogStage )

    @chk ccall(
               (:PetscLogHandlerStagePush, $petsc_library),
               PetscErrorCode,
               (PetscLogHandler, PetscLogStage),
               h, stage,
              )


	return nothing
end 

"""
	PetscLogHandlerStageSetVisible(petsclib::PetscLibType,handler::PetscLogHandler, stage::PetscLogStage, isVisible::PetscBool) 
Set the visibility of logging stage in `PetscLogHandlerView()` for a log handler

Not collective

Input Parameters:
- `handler`   - a `PetscLogHandler`
- `stage`     - a `PetscLogStage`
- `isVisible` - the visibility flag, `PETSC_TRUE` to print, else `PETSC_FALSE` (defaults to `PETSC_TRUE`)

Level: developer

-seealso: [](ch_profiling), `PetscLogHandlerStageGetVisible()`

# External Links
$(_doc_external("Log/PetscLogHandlerStageSetVisible"))
"""
function PetscLogHandlerStageSetVisible(petsclib::PetscLibType, handler::PetscLogHandler, stage::PetscLogStage, isVisible::PetscBool)
    error("PetscLogHandlerStageSetVisible: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerStageSetVisible(petsclib::$UnionPetscLib, handler::PetscLogHandler, stage::PetscLogStage, isVisible::PetscBool )

    @chk ccall(
               (:PetscLogHandlerStageSetVisible, $petsc_library),
               PetscErrorCode,
               (PetscLogHandler, PetscLogStage, PetscBool),
               handler, stage, isVisible,
              )


	return nothing
end 

"""
	PetscLogHandlerStart(petsclib::PetscLibType,h::PetscLogHandler) 
Connect a log handler to PETSc's global logging stream and state.

Logically collective

Input Parameters:
- `h` - a `PetscLogHandler`

Level: developer

-seealso: [](ch_profiling), `PetscLogHandler`, `PetscLogState`, `PetscLogHandlerStop()`, `PetscInitialize()`

# External Links
$(_doc_external("Log/PetscLogHandlerStart"))
"""
function PetscLogHandlerStart(petsclib::PetscLibType, h::PetscLogHandler)
    error("PetscLogHandlerStart: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerStart(petsclib::$UnionPetscLib, h::PetscLogHandler )

    @chk ccall(
               (:PetscLogHandlerStart, $petsc_library),
               PetscErrorCode,
               (PetscLogHandler,),
               h,
              )


	return nothing
end 

"""
	PetscLogHandlerStop(petsclib::PetscLibType,h::PetscLogHandler) 
Disconnect a log handler from PETSc's global logging stream.

Logically collective

Input Parameters:
- `h` - a `PetscLogHandler`

Level: developer

-seealso: [](ch_profiling), `PetscLogHandler`, `PetscLogState`, `PetscLogHandlerStart()`

# External Links
$(_doc_external("Log/PetscLogHandlerStop"))
"""
function PetscLogHandlerStop(petsclib::PetscLibType, h::PetscLogHandler)
    error("PetscLogHandlerStop: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerStop(petsclib::$UnionPetscLib, h::PetscLogHandler )

    @chk ccall(
               (:PetscLogHandlerStop, $petsc_library),
               PetscErrorCode,
               (PetscLogHandler,),
               h,
              )


	return nothing
end 

"""
	PetscLogHandlerView(petsclib::PetscLibType,h::PetscLogHandler, viewer::PetscViewer) 
View the data recorded in a log handler.

Collective

Input Parameters:
- `h`      - the `PetscLogHandler`
- `viewer` - the `PetscViewer`

Level: developer

-seealso: [](ch_profiling), `PetscLogHandler`, `PetscLogView()`

# External Links
$(_doc_external("Log/PetscLogHandlerView"))
"""
function PetscLogHandlerView(petsclib::PetscLibType, h::PetscLogHandler, viewer::PetscViewer)
    error("PetscLogHandlerView: no generated method for these argument types")
end

@for_petsc function PetscLogHandlerView(petsclib::$UnionPetscLib, h::PetscLogHandler, viewer::PetscViewer )

    @chk ccall(
               (:PetscLogHandlerView, $petsc_library),
               PetscErrorCode,
               (PetscLogHandler, PetscViewer),
               h, viewer,
              )


	return nothing
end 

"""
	info::PetscLogClassInfo = PetscLogStateClassGetInfo(petsclib::PetscLibType,state::PetscLogState, clss::PetscLogClass) 
Get the registration information of an class

Not collective

Input Parameters:
- `state` - a `PetscLogState`
- `clss`  - a registered `PetscLogClass`

Output Parameter:
- `info` - the `PetscLogClassInfo` of the class will be copied into info

Level: developer

-seealso: [](ch_profiling), `PetscLogState`, `PetscLogStateClassRegister()`, `PetscLogStateGetClassFromName()`

# External Links
$(_doc_external("Log/PetscLogStateClassGetInfo"))
"""
function PetscLogStateClassGetInfo(petsclib::PetscLibType, state::PetscLogState, clss::PetscLogClass)
    error("PetscLogStateClassGetInfo: no generated method for these argument types")
end

@for_petsc function PetscLogStateClassGetInfo(petsclib::$UnionPetscLib, state::PetscLogState, clss::PetscLogClass )
	info_ = Ref{PetscLogClassInfo}()

    @chk ccall(
               (:PetscLogStateClassGetInfo, $petsc_library),
               PetscErrorCode,
               (PetscLogState, PetscLogClass, Ptr{PetscLogClassInfo}),
               state, clss, info_,
              )

	info = info_[]

	return info
end 

"""
	logclass::PetscLogClass = PetscLogStateClassRegister(petsclib::PetscLibType,state::PetscLogState, name::String, id::PetscClassId) 
Register a class to with a `PetscLogState` used by `PetscLogHandler`s.

Logically collective on `PETSC_COMM_WORLD`

Input Parameters:
- `state` - a `PetscLogState`
- `name`  - the name of a class registered with `PetscClassIdRegister()`
- `id`    - the `PetscClassId` obtained from `PetscClassIdRegister()`

Output Parameter:
- `logclass` - a `PetscLogClass` for this class with this state

Level: developer

-seealso: [](ch_profiling), `PetscLogStateClassGetInfo()`, `PetscLogStateGetClassFromName()`, `PetscLogStateGetClassFromClassId()`

# External Links
$(_doc_external("Log/PetscLogStateClassRegister"))
"""
function PetscLogStateClassRegister(petsclib::PetscLibType, state::PetscLogState, name::String, id::PetscClassId)
    error("PetscLogStateClassRegister: no generated method for these argument types")
end

@for_petsc function PetscLogStateClassRegister(petsclib::$UnionPetscLib, state::PetscLogState, name::String, id::PetscClassId )
	logclass_ = Ref{PetscLogClass}()

    @chk ccall(
               (:PetscLogStateClassRegister, $petsc_library),
               PetscErrorCode,
               (PetscLogState, Ptr{Cchar}, PetscClassId, Ptr{PetscLogClass}),
               state, name, id, logclass_,
              )

	logclass = logclass_[]

	return logclass
end 

"""
	PetscLogStateClassSetActive(petsclib::PetscLibType,state::PetscLogState, stage::PetscLogStage, classid::PetscClassId, isActive::PetscBool) 
Set logging events associated with an event as active or inactive during a logging stage.

Not collective

Input Parameters:
- `state`    - a `PetscLogState`
- `stage`    - a registered `PetscLogStage`, or `PETSC_DEFAULT` for the current stage
- `classid`  - a `PetscClassId`
- `isActive` - if `PETSC_FALSE`, `PetscLogStateEventGetActive()` will return
`PETSC_FALSE` for this stage and all events that were associated
with this class when they were registered (see
`PetscLogStateEventRegister()`).

Level: developer

-seealso: [](ch_profiling), `PetscLogState`, `PetscLogEventGetActive()`, `PetscLogStateEventSetActive()`

# External Links
$(_doc_external("Log/PetscLogStateClassSetActive"))
"""
function PetscLogStateClassSetActive(petsclib::PetscLibType, state::PetscLogState, stage::PetscLogStage, classid::PetscClassId, isActive::PetscBool)
    error("PetscLogStateClassSetActive: no generated method for these argument types")
end

@for_petsc function PetscLogStateClassSetActive(petsclib::$UnionPetscLib, state::PetscLogState, stage::PetscLogStage, classid::PetscClassId, isActive::PetscBool )

    @chk ccall(
               (:PetscLogStateClassSetActive, $petsc_library),
               PetscErrorCode,
               (PetscLogState, PetscLogStage, PetscClassId, PetscBool),
               state, stage, classid, isActive,
              )


	return nothing
end 

"""
	PetscLogStateClassSetActiveAll(petsclib::PetscLibType,state::PetscLogState, classid::PetscClassId, isActive::PetscBool) 
Set logging events associated with an event as active or inactive for all logging stages

Not collective

Input Parameters:
- `state`    - a `PetscLogState`
- `classid`  - a `PetscClassId`
- `isActive` - if `PETSC_FALSE`, `PetscLogStateEventGetActive()` will return
`PETSC_FALSE` for all events that were associated with this class when they
were registered (see `PetscLogStateEventRegister()`).

Level: developer

-seealso: [](ch_profiling), `PetscLogState`, `PetscLogEventGetActive()`, `PetscLogStateClassSetActive()`

# External Links
$(_doc_external("Log/PetscLogStateClassSetActiveAll"))
"""
function PetscLogStateClassSetActiveAll(petsclib::PetscLibType, state::PetscLogState, classid::PetscClassId, isActive::PetscBool)
    error("PetscLogStateClassSetActiveAll: no generated method for these argument types")
end

@for_petsc function PetscLogStateClassSetActiveAll(petsclib::$UnionPetscLib, state::PetscLogState, classid::PetscClassId, isActive::PetscBool )

    @chk ccall(
               (:PetscLogStateClassSetActiveAll, $petsc_library),
               PetscErrorCode,
               (PetscLogState, PetscClassId, PetscBool),
               state, classid, isActive,
              )


	return nothing
end 

"""
	state::PetscLogState = PetscLogStateCreate(petsclib::PetscLibType) 
Create a logging state.

Not collective

Output Parameters:
- `state` - a `PetscLogState`

Level: developer

-seealso: [](ch_profiling), `PetscLogState`, `PetscLogStateDestroy()`

# External Links
$(_doc_external("Log/PetscLogStateCreate"))
"""
function PetscLogStateCreate(petsclib::PetscLibType)
    error("PetscLogStateCreate: no generated method for these argument types")
end

@for_petsc function PetscLogStateCreate(petsclib::$UnionPetscLib)
	state_ = Ref{PetscLogState}()

    @chk ccall(
               (:PetscLogStateCreate, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscLogState},),
               state_,
              )

	state = state_[]

	return state
end 

"""
	PetscLogStateDestroy(petsclib::PetscLibType,state::Union{PetscLogState, Ref{PetscLogState}}) 
Destroy a logging state.

Not collective

Input Parameters:
- `state` - a `PetscLogState`

Level: developer

-seealso: [](ch_profiling), `PetscLogState`, `PetscLogStateCreate()`

# External Links
$(_doc_external("Log/PetscLogStateDestroy"))
"""
function PetscLogStateDestroy(petsclib::PetscLibType, state::Union{PetscLogState, Ref{PetscLogState}})
    error("PetscLogStateDestroy: no generated method for these argument types")
end

@for_petsc function PetscLogStateDestroy(petsclib::$UnionPetscLib, state::Union{PetscLogState, Ref{PetscLogState}} )
	state_ = state isa Base.RefValue ? state : Ref{PetscLogState}(state)

    @chk ccall(
               (:PetscLogStateDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscLogState},),
               state_,
              )


	return nothing
end 

"""
	isActive::PetscBool = PetscLogStateEventGetActive(petsclib::PetscLibType,state::PetscLogState, stage::PetscLogStage, event::PetscLogEvent) 
Check if a logging event is active or inactive during a logging stage.

Not collective

Input Parameters:
- `state` - a `PetscLogState`
- `stage` - a registered `PetscLogStage`, or `PETSC_DEFAULT` for the current stage
- `event` - a registered `PetscLogEvent`

Output Parameter:
- `isActive` - If `PETSC_FALSE`, log handlers should not be notified of the event's beginning or end.

Level: developer

-seealso: [](ch_profiling), `PetscLogState`, `PetscLogEventGetActive()`, `PetscLogStateGetCurrentStage()`, `PetscLogHandler()`

# External Links
$(_doc_external("Log/PetscLogStateEventGetActive"))
"""
function PetscLogStateEventGetActive(petsclib::PetscLibType, state::PetscLogState, stage::PetscLogStage, event::PetscLogEvent)
    error("PetscLogStateEventGetActive: no generated method for these argument types")
end

@for_petsc function PetscLogStateEventGetActive(petsclib::$UnionPetscLib, state::PetscLogState, stage::PetscLogStage, event::PetscLogEvent )
	isActive_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscLogStateEventGetActive, $petsc_library),
               PetscErrorCode,
               (PetscLogState, PetscLogStage, PetscLogEvent, Ptr{PetscBool}),
               state, stage, event, isActive_,
              )

	isActive = isActive_[]

	return isActive
end 

"""
	info::PetscLogEventInfo = PetscLogStateEventGetInfo(petsclib::PetscLibType,state::PetscLogState, event::PetscLogEvent) 
Get the registration information of an event

Not collective

Input Parameters:
- `state` - a `PetscLogState`
- `event` - a registered `PetscLogEvent`

Output Parameter:
- `info` - the `PetscLogEventInfo` of the event will be copied into info

Level: developer

-seealso: [](ch_profiling), `PetscLogState`, `PetscLogStateEventRegister()`, `PetscLogStateGetEventFromName()`

# External Links
$(_doc_external("Log/PetscLogStateEventGetInfo"))
"""
function PetscLogStateEventGetInfo(petsclib::PetscLibType, state::PetscLogState, event::PetscLogEvent)
    error("PetscLogStateEventGetInfo: no generated method for these argument types")
end

@for_petsc function PetscLogStateEventGetInfo(petsclib::$UnionPetscLib, state::PetscLogState, event::PetscLogEvent )
	info_ = Ref{PetscLogEventInfo}()

    @chk ccall(
               (:PetscLogStateEventGetInfo, $petsc_library),
               PetscErrorCode,
               (PetscLogState, PetscLogEvent, Ptr{PetscLogEventInfo}),
               state, event, info_,
              )

	info = info_[]

	return info
end 

"""
	event::PetscLogEvent = PetscLogStateEventRegister(petsclib::PetscLibType,state::PetscLogState, sname::String, id::PetscClassId) 
Register a new event with a logging state

Not collective

Input Parameters:
- `state` - a `PetscLogState`
- `sname` - a unique name
- `id`    - the `PetscClassId` for the type of object most closely associated with this event

Output Parameter:
- `event` - the identifier for the registered event

Level: developer

-seealso: [](ch_profiling), `PetscLogState`, `PetscLogStageRegister()`

# External Links
$(_doc_external("Log/PetscLogStateEventRegister"))
"""
function PetscLogStateEventRegister(petsclib::PetscLibType, state::PetscLogState, sname::String, id::PetscClassId)
    error("PetscLogStateEventRegister: no generated method for these argument types")
end

@for_petsc function PetscLogStateEventRegister(petsclib::$UnionPetscLib, state::PetscLogState, sname::String, id::PetscClassId )
	event_ = Ref{PetscLogEvent}()

    @chk ccall(
               (:PetscLogStateEventRegister, $petsc_library),
               PetscErrorCode,
               (PetscLogState, Ptr{Cchar}, PetscClassId, Ptr{PetscLogEvent}),
               state, sname, id, event_,
              )

	event = event_[]

	return event
end 

"""
	PetscLogStateEventSetActive(petsclib::PetscLibType,state::PetscLogState, stage::PetscLogStage, event::PetscLogEvent, isActive::PetscBool) 
Set a logging event as active or inactive during a logging stage.

Not collective

Input Parameters:
- `state`    - a `PetscLogState`
- `stage`    - a registered `PetscLogStage`, or `PETSC_DEFAULT` for the current stage
- `event`    - a registered `PetscLogEvent`
- `isActive` - if `PETSC_FALSE`, `PetscLogStateEventGetActive()` will return `PETSC_FALSE` for this stage and this event

Level: developer

-seealso: [](ch_profiling), `PetscLogState`, `PetscLogEventGetActive()`, `PetscLogStateGetCurrentStage()`, `PetscLogEventSetActiveAll()`

# External Links
$(_doc_external("Log/PetscLogStateEventSetActive"))
"""
function PetscLogStateEventSetActive(petsclib::PetscLibType, state::PetscLogState, stage::PetscLogStage, event::PetscLogEvent, isActive::PetscBool)
    error("PetscLogStateEventSetActive: no generated method for these argument types")
end

@for_petsc function PetscLogStateEventSetActive(petsclib::$UnionPetscLib, state::PetscLogState, stage::PetscLogStage, event::PetscLogEvent, isActive::PetscBool )

    @chk ccall(
               (:PetscLogStateEventSetActive, $petsc_library),
               PetscErrorCode,
               (PetscLogState, PetscLogStage, PetscLogEvent, PetscBool),
               state, stage, event, isActive,
              )


	return nothing
end 

"""
	PetscLogStateEventSetActiveAll(petsclib::PetscLibType,state::PetscLogState, event::PetscLogEvent, isActive::PetscBool) 
Set logging event as active or inactive for all logging stages

Not collective

Input Parameters:
- `state`    - a `PetscLogState`
- `event`    - a registered `PetscLogEvent`
- `isActive` - if `PETSC_FALSE`, `PetscLogStateEventGetActive()` will return `PETSC_FALSE` for all stages and all events

Level: developer

-seealso: [](ch_profiling), `PetscLogState`, `PetscLogEventGetActive()`

# External Links
$(_doc_external("Log/PetscLogStateEventSetActiveAll"))
"""
function PetscLogStateEventSetActiveAll(petsclib::PetscLibType, state::PetscLogState, event::PetscLogEvent, isActive::PetscBool)
    error("PetscLogStateEventSetActiveAll: no generated method for these argument types")
end

@for_petsc function PetscLogStateEventSetActiveAll(petsclib::$UnionPetscLib, state::PetscLogState, event::PetscLogEvent, isActive::PetscBool )

    @chk ccall(
               (:PetscLogStateEventSetActiveAll, $petsc_library),
               PetscErrorCode,
               (PetscLogState, PetscLogEvent, PetscBool),
               state, event, isActive,
              )


	return nothing
end 

"""
	PetscLogStateEventSetCollective(petsclib::PetscLibType,state::PetscLogState, event::PetscLogEvent, collective::PetscBool) 
Set the collective nature of a logging event

Logically collective

Input Parameters:
- `state`      - a `PetscLogState`
- `event`      - a registered `PetscLogEvent`
- `collective` - if `PETSC_TRUE`, MPI processes synchronize during this event, and `PetscLogHandlerEventSync()` can be used to help measure the delays between when the processes begin the event

Level: developer

-seealso: [](ch_profiling), `PetscLogState`, `PetscLogEventRegister()`

# External Links
$(_doc_external("Log/PetscLogStateEventSetCollective"))
"""
function PetscLogStateEventSetCollective(petsclib::PetscLibType, state::PetscLogState, event::PetscLogEvent, collective::PetscBool)
    error("PetscLogStateEventSetCollective: no generated method for these argument types")
end

@for_petsc function PetscLogStateEventSetCollective(petsclib::$UnionPetscLib, state::PetscLogState, event::PetscLogEvent, collective::PetscBool )

    @chk ccall(
               (:PetscLogStateEventSetCollective, $petsc_library),
               PetscErrorCode,
               (PetscLogState, PetscLogEvent, PetscBool),
               state, event, collective,
              )


	return nothing
end 

"""
	clss::PetscLogClass = PetscLogStateGetClassFromClassId(petsclib::PetscLibType,state::PetscLogState, classid::PetscClassId) 
Get a `PetscLogClass` from the `PetscClassId` it was registered with.

Not collective

Input Parameters:
- `state`   - a `PetscLogState`
- `classid` - a `PetscClassId`

Output Parameter:
- `clss` - the classes's logging id

Level: developer

-seealso: [](ch_profiling), `PetscLogState`, `PetscLogStateClassRegister()`, `PetscLogStateClassGetInfo()`

# External Links
$(_doc_external("Log/PetscLogStateGetClassFromClassId"))
"""
function PetscLogStateGetClassFromClassId(petsclib::PetscLibType, state::PetscLogState, classid::PetscClassId)
    error("PetscLogStateGetClassFromClassId: no generated method for these argument types")
end

@for_petsc function PetscLogStateGetClassFromClassId(petsclib::$UnionPetscLib, state::PetscLogState, classid::PetscClassId )
	clss_ = Ref{PetscLogClass}()

    @chk ccall(
               (:PetscLogStateGetClassFromClassId, $petsc_library),
               PetscErrorCode,
               (PetscLogState, PetscClassId, Ptr{PetscLogClass}),
               state, classid, clss_,
              )

	clss = clss_[]

	return clss
end 

"""
	clss::PetscLogClass = PetscLogStateGetClassFromName(petsclib::PetscLibType,state::PetscLogState, name::String) 
Get a `PetscLogClass` from the name of the class it was registered with.

Not collective

Input Parameters:
- `state` - a `PetscLogState`
- `name`  - the name string of the class

Output Parameter:
- `clss` - the classes's logging id

Level: developer

-seealso: [](ch_profiling), `PetscLogState`, `PetscLogStateClassRegister()`, `PetscLogStateClassGetInfo()`

# External Links
$(_doc_external("Log/PetscLogStateGetClassFromName"))
"""
function PetscLogStateGetClassFromName(petsclib::PetscLibType, state::PetscLogState, name::String)
    error("PetscLogStateGetClassFromName: no generated method for these argument types")
end

@for_petsc function PetscLogStateGetClassFromName(petsclib::$UnionPetscLib, state::PetscLogState, name::String )
	clss_ = Ref{PetscLogClass}()

    @chk ccall(
               (:PetscLogStateGetClassFromName, $petsc_library),
               PetscErrorCode,
               (PetscLogState, Ptr{Cchar}, Ptr{PetscLogClass}),
               state, name, clss_,
              )

	clss = clss_[]

	return clss
end 

"""
	current::PetscLogStage = PetscLogStateGetCurrentStage(petsclib::PetscLibType,state::PetscLogState) 
Get the last stage that was started

Not collective

Input Parameter:
- `state` - a `PetscLogState`

Output Parameter:
- `current` - the last `PetscLogStage` started with `PetscLogStateStagePop()`

Level: developer

-seealso: [](ch_profiling), `PetscLogState`, `PetscLogStateStageRegister()`, `PetscLogStateStagePush()`, `PetscLogStateStagePop()`

# External Links
$(_doc_external("Log/PetscLogStateGetCurrentStage"))
"""
function PetscLogStateGetCurrentStage(petsclib::PetscLibType, state::PetscLogState)
    error("PetscLogStateGetCurrentStage: no generated method for these argument types")
end

@for_petsc function PetscLogStateGetCurrentStage(petsclib::$UnionPetscLib, state::PetscLogState )
	current_ = Ref{PetscLogStage}()

    @chk ccall(
               (:PetscLogStateGetCurrentStage, $petsc_library),
               PetscErrorCode,
               (PetscLogState, Ptr{PetscLogStage}),
               state, current_,
              )

	current = current_[]

	return current
end 

"""
	event::PetscLogEvent = PetscLogStateGetEventFromName(petsclib::PetscLibType,state::PetscLogState, name::String) 
Get a `PetscLogEvent` from the name it was registered with.

Not collective

Input Parameters:
- `state` - a `PetscLogState`
- `name`  - an event's name

Output Parameter:
- `event` - the event's id

Level: developer

-seealso: [](ch_profiling), `PetscLogState`, `PetscLogStateEventRegister()`, `PetscLogStateEventGetInfo()`

# External Links
$(_doc_external("Log/PetscLogStateGetEventFromName"))
"""
function PetscLogStateGetEventFromName(petsclib::PetscLibType, state::PetscLogState, name::String)
    error("PetscLogStateGetEventFromName: no generated method for these argument types")
end

@for_petsc function PetscLogStateGetEventFromName(petsclib::$UnionPetscLib, state::PetscLogState, name::String )
	event_ = Ref{PetscLogEvent}()

    @chk ccall(
               (:PetscLogStateGetEventFromName, $petsc_library),
               PetscErrorCode,
               (PetscLogState, Ptr{Cchar}, Ptr{PetscLogEvent}),
               state, name, event_,
              )

	event = event_[]

	return event
end 

"""
	numClasses::PetscInt = PetscLogStateGetNumClasses(petsclib::PetscLibType,state::PetscLogState) 
Get the number of registered classes in a logging state.

Not collective

Input Parameter:
- `state` - a `PetscLogState`

Output Parameter:
- `numClasses` - the number of registered `PetscLogClass`s

Level: developer

-seealso: [](ch_profiling), `PetscLogState`, `PetscLogStateClassRegister()`

# External Links
$(_doc_external("Log/PetscLogStateGetNumClasses"))
"""
function PetscLogStateGetNumClasses(petsclib::PetscLibType, state::PetscLogState)
    error("PetscLogStateGetNumClasses: no generated method for these argument types")
end

@for_petsc function PetscLogStateGetNumClasses(petsclib::$UnionPetscLib, state::PetscLogState )
	numClasses_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscLogStateGetNumClasses, $petsc_library),
               PetscErrorCode,
               (PetscLogState, Ptr{$PetscInt}),
               state, numClasses_,
              )

	numClasses = numClasses_[]

	return numClasses
end 

"""
	numEvents::PetscInt = PetscLogStateGetNumEvents(petsclib::PetscLibType,state::PetscLogState) 
Get the number of registered events in a logging state.

Not collective

Input Parameter:
- `state` - a `PetscLogState`

Output Parameter:
- `numEvents` - the number of registered `PetscLogEvent`s

Level: developer

-seealso: [](ch_profiling), `PetscLogState`, `PetscLogStateEventRegister()`

# External Links
$(_doc_external("Log/PetscLogStateGetNumEvents"))
"""
function PetscLogStateGetNumEvents(petsclib::PetscLibType, state::PetscLogState)
    error("PetscLogStateGetNumEvents: no generated method for these argument types")
end

@for_petsc function PetscLogStateGetNumEvents(petsclib::$UnionPetscLib, state::PetscLogState )
	numEvents_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscLogStateGetNumEvents, $petsc_library),
               PetscErrorCode,
               (PetscLogState, Ptr{$PetscInt}),
               state, numEvents_,
              )

	numEvents = numEvents_[]

	return numEvents
end 

"""
	numStages::PetscInt = PetscLogStateGetNumStages(petsclib::PetscLibType,state::PetscLogState) 
Get the number of registered stages in a logging state.

Not collective

Input Parameter:
- `state` - a `PetscLogState`

Output Parameter:
- `numStages` - the number of registered `PetscLogStage`s

Level: developer

-seealso: [](ch_profiling), `PetscLogState`, `PetscLogStateStageRegister()`

# External Links
$(_doc_external("Log/PetscLogStateGetNumStages"))
"""
function PetscLogStateGetNumStages(petsclib::PetscLibType, state::PetscLogState)
    error("PetscLogStateGetNumStages: no generated method for these argument types")
end

@for_petsc function PetscLogStateGetNumStages(petsclib::$UnionPetscLib, state::PetscLogState )
	numStages_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscLogStateGetNumStages, $petsc_library),
               PetscErrorCode,
               (PetscLogState, Ptr{$PetscInt}),
               state, numStages_,
              )

	numStages = numStages_[]

	return numStages
end 

"""
	stage::PetscLogStage = PetscLogStateGetStageFromName(petsclib::PetscLibType,state::PetscLogState, name::String) 
Get a `PetscLogStage` from the name it was registered with.

Not collective

Input Parameters:
- `state` - a `PetscLogState`
- `name`  - a stage's name

Output Parameter:
- `stage` - the stage's id

Level: developer

-seealso: [](ch_profiling), `PetscLogState`, `PetscLogStateStageRegister()`, `PetscLogStateStageGetInfo()`

# External Links
$(_doc_external("Log/PetscLogStateGetStageFromName"))
"""
function PetscLogStateGetStageFromName(petsclib::PetscLibType, state::PetscLogState, name::String)
    error("PetscLogStateGetStageFromName: no generated method for these argument types")
end

@for_petsc function PetscLogStateGetStageFromName(petsclib::$UnionPetscLib, state::PetscLogState, name::String )
	stage_ = Ref{PetscLogStage}()

    @chk ccall(
               (:PetscLogStateGetStageFromName, $petsc_library),
               PetscErrorCode,
               (PetscLogState, Ptr{Cchar}, Ptr{PetscLogStage}),
               state, name, stage_,
              )

	stage = stage_[]

	return stage
end 

"""
	isActive::PetscBool = PetscLogStateStageGetActive(petsclib::PetscLibType,state::PetscLogState, stage::PetscLogStage) 
Check if a logging stage is active or inactive.

Not collective

Input Parameters:
- `state` - a `PetscLogState`
- `stage` - a registered `PetscLogStage`

Output Parameter:
- `isActive` - if `PETSC_FALSE`, the state should not send logging events to log handlers during this stage.

Level: developer

-seealso: [](ch_profiling), `PetscLogState`, `PetscLogStageSetActive()`, `PetscLogHandler`, `PetscLogHandlerStart()`, `PetscLogHandlerEventBegin()`, `PetscLogHandlerEventEnd()`

# External Links
$(_doc_external("Log/PetscLogStateStageGetActive"))
"""
function PetscLogStateStageGetActive(petsclib::PetscLibType, state::PetscLogState, stage::PetscLogStage)
    error("PetscLogStateStageGetActive: no generated method for these argument types")
end

@for_petsc function PetscLogStateStageGetActive(petsclib::$UnionPetscLib, state::PetscLogState, stage::PetscLogStage )
	isActive_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscLogStateStageGetActive, $petsc_library),
               PetscErrorCode,
               (PetscLogState, PetscLogStage, Ptr{PetscBool}),
               state, stage, isActive_,
              )

	isActive = isActive_[]

	return isActive
end 

"""
	info::PetscLogStageInfo = PetscLogStateStageGetInfo(petsclib::PetscLibType,state::PetscLogState, stage::PetscLogStage) 
Get the registration information of an stage

Not collective

Input Parameters:
- `state` - a `PetscLogState`
- `stage` - a registered `PetscLogStage`

Output Parameter:
- `info` - the `PetscLogStageInfo` of the stage will be copied into info

Level: developer

-seealso: [](ch_profiling), `PetscLogState`, `PetscLogStateStageRegister()`, `PetscLogStateGetStageFromName()`

# External Links
$(_doc_external("Log/PetscLogStateStageGetInfo"))
"""
function PetscLogStateStageGetInfo(petsclib::PetscLibType, state::PetscLogState, stage::PetscLogStage)
    error("PetscLogStateStageGetInfo: no generated method for these argument types")
end

@for_petsc function PetscLogStateStageGetInfo(petsclib::$UnionPetscLib, state::PetscLogState, stage::PetscLogStage )
	info_ = Ref{PetscLogStageInfo}()

    @chk ccall(
               (:PetscLogStateStageGetInfo, $petsc_library),
               PetscErrorCode,
               (PetscLogState, PetscLogStage, Ptr{PetscLogStageInfo}),
               state, stage, info_,
              )

	info = info_[]

	return info
end 

"""
	PetscLogStateStagePop(petsclib::PetscLibType,state::PetscLogState) 
End a running logging stage.

Not collective

Input Parameter:
- `state` - a `PetscLogState`

Level: developer

-seealso: [](ch_profiling), `PetscLogState`, `PetscLogStateStageRegister()`, `PetscLogStateStagePush()`, `PetscLogStateGetCurrentStage()`

# External Links
$(_doc_external("Log/PetscLogStateStagePop"))
"""
function PetscLogStateStagePop(petsclib::PetscLibType, state::PetscLogState)
    error("PetscLogStateStagePop: no generated method for these argument types")
end

@for_petsc function PetscLogStateStagePop(petsclib::$UnionPetscLib, state::PetscLogState )

    @chk ccall(
               (:PetscLogStateStagePop, $petsc_library),
               PetscErrorCode,
               (PetscLogState,),
               state,
              )


	return nothing
end 

"""
	PetscLogStateStagePush(petsclib::PetscLibType,state::PetscLogState, stage::PetscLogStage) 
Start a new logging stage.

Not collective

Input Parameters:
- `state` - a `PetscLogState`
- `stage` - a registered `PetscLogStage`

Level: developer

-seealso: [](ch_profiling), `PetscLogState`, `PetscLogStateStageRegister()`, `PetscLogStateStagePop()`, `PetscLogStateGetCurrentStage()`

# External Links
$(_doc_external("Log/PetscLogStateStagePush"))
"""
function PetscLogStateStagePush(petsclib::PetscLibType, state::PetscLogState, stage::PetscLogStage)
    error("PetscLogStateStagePush: no generated method for these argument types")
end

@for_petsc function PetscLogStateStagePush(petsclib::$UnionPetscLib, state::PetscLogState, stage::PetscLogStage )

    @chk ccall(
               (:PetscLogStateStagePush, $petsc_library),
               PetscErrorCode,
               (PetscLogState, PetscLogStage),
               state, stage,
              )


	return nothing
end 

"""
	stage::PetscLogStage = PetscLogStateStageRegister(petsclib::PetscLibType,state::PetscLogState, sname::String) 
Register a new stage with a logging state

Not collective

Input Parameters:
- `state` - a `PetscLogState`
- `sname` - a unique name

Output Parameter:
- `stage` - the identifier for the registered stage

Level: developer

-seealso: [](ch_profiling), `PetscLogState`, `PetscLogStateStagePush()`, `PetscLogStateStagePop()`

# External Links
$(_doc_external("Log/PetscLogStateStageRegister"))
"""
function PetscLogStateStageRegister(petsclib::PetscLibType, state::PetscLogState, sname::String)
    error("PetscLogStateStageRegister: no generated method for these argument types")
end

@for_petsc function PetscLogStateStageRegister(petsclib::$UnionPetscLib, state::PetscLogState, sname::String )
	stage_ = Ref{PetscLogStage}()

    @chk ccall(
               (:PetscLogStateStageRegister, $petsc_library),
               PetscErrorCode,
               (PetscLogState, Ptr{Cchar}, Ptr{PetscLogStage}),
               state, sname, stage_,
              )

	stage = stage_[]

	return stage
end 

"""
	PetscLogStateStageSetActive(petsclib::PetscLibType,state::PetscLogState, stage::PetscLogStage, isActive::PetscBool) 
Mark a stage as active or inactive.

Not collective

Input Parameters:
- `state`    - a `PetscLogState`
- `stage`    - a registered `PetscLogStage`
- `isActive` - if `PETSC_FALSE`, `PetscLogStateEventGetActive()` will return `PETSC_FALSE` for all events during this stage

Level: developer

-seealso: [](ch_profiling), `PetscLogState`, `PetscLogStateEventSetActive()`

# External Links
$(_doc_external("Log/PetscLogStateStageSetActive"))
"""
function PetscLogStateStageSetActive(petsclib::PetscLibType, state::PetscLogState, stage::PetscLogStage, isActive::PetscBool)
    error("PetscLogStateStageSetActive: no generated method for these argument types")
end

@for_petsc function PetscLogStateStageSetActive(petsclib::$UnionPetscLib, state::PetscLogState, stage::PetscLogStage, isActive::PetscBool )

    @chk ccall(
               (:PetscLogStateStageSetActive, $petsc_library),
               PetscErrorCode,
               (PetscLogState, PetscLogStage, PetscBool),
               state, stage, isActive,
              )


	return nothing
end 

