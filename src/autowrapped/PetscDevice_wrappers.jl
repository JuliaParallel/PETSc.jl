# autodefined type arguments for class ------
mutable struct _n_PetscDevice end
const PetscDevice = Ptr{_n_PetscDevice}

mutable struct _n_PetscDeviceContext end
const PetscDeviceContext = Ptr{_n_PetscDeviceContext}
# -------------------------------------------------------

"""
	device::PetscDevice = PetscDeviceCreate(petsclib::PetscLibType,type::PetscDeviceType, devid::PetscInt) 

# External Links
$(_doc_external("Sys/PetscDeviceCreate"))
"""
function PetscDeviceCreate(petsclib::PetscLibType, type::PetscDeviceType, devid::PetscInt) end

@for_petsc function PetscDeviceCreate(petsclib::$UnionPetscLib, type::PetscDeviceType, devid::$PetscInt )
	device_ = Ref{PetscDevice}()

    @chk ccall(
               (:PetscDeviceCreate, $petsc_library),
               PetscErrorCode,
               (PetscDeviceType, $PetscInt, Ptr{PetscDevice}),
               type, devid, device_,
              )

	device = device_[]

	return device
end 

"""
	PetscDeviceDestroy(petsclib::PetscLibType,device::PetscDevice) 

# External Links
$(_doc_external("Sys/PetscDeviceDestroy"))
"""
function PetscDeviceDestroy(petsclib::PetscLibType, device::PetscDevice) end

@for_petsc function PetscDeviceDestroy(petsclib::$UnionPetscLib, device::PetscDevice )

    @chk ccall(
               (:PetscDeviceDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscDevice},),
               device,
              )


	return nothing
end 

"""
	PetscDeviceConfigure(petsclib::PetscLibType,device::PetscDevice) 

# External Links
$(_doc_external("Sys/PetscDeviceConfigure"))
"""
function PetscDeviceConfigure(petsclib::PetscLibType, device::PetscDevice) end

@for_petsc function PetscDeviceConfigure(petsclib::$UnionPetscLib, device::PetscDevice )

    @chk ccall(
               (:PetscDeviceConfigure, $petsc_library),
               PetscErrorCode,
               (PetscDevice,),
               device,
              )


	return nothing
end 

"""
	PetscDeviceView(petsclib::PetscLibType,device::PetscDevice, viewer::PetscViewer) 

# External Links
$(_doc_external("Sys/PetscDeviceView"))
"""
function PetscDeviceView(petsclib::PetscLibType, device::PetscDevice, viewer::PetscViewer) end

@for_petsc function PetscDeviceView(petsclib::$UnionPetscLib, device::PetscDevice, viewer::PetscViewer )

    @chk ccall(
               (:PetscDeviceView, $petsc_library),
               PetscErrorCode,
               (PetscDevice, PetscViewer),
               device, viewer,
              )


	return nothing
end 

"""
	type::PetscDeviceType = PetscDeviceGetType(petsclib::PetscLibType,device::PetscDevice) 

# External Links
$(_doc_external("Sys/PetscDeviceGetType"))
"""
function PetscDeviceGetType(petsclib::PetscLibType, device::PetscDevice) end

@for_petsc function PetscDeviceGetType(petsclib::$UnionPetscLib, device::PetscDevice )
	type_ = Ref{PetscDeviceType}()

    @chk ccall(
               (:PetscDeviceGetType, $petsc_library),
               PetscErrorCode,
               (PetscDevice, Ptr{PetscDeviceType}),
               device, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	id::PetscInt = PetscDeviceGetDeviceId(petsclib::PetscLibType,device::PetscDevice) 

# External Links
$(_doc_external("Sys/PetscDeviceGetDeviceId"))
"""
function PetscDeviceGetDeviceId(petsclib::PetscLibType, device::PetscDevice) end

@for_petsc function PetscDeviceGetDeviceId(petsclib::$UnionPetscLib, device::PetscDevice )
	id_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDeviceGetDeviceId, $petsc_library),
               PetscErrorCode,
               (PetscDevice, Ptr{$PetscInt}),
               device, id_,
              )

	id = id_[]

	return id
end 

"""
	PetscDeviceSetDefaultDeviceType(petsclib::PetscLibType,type::PetscDeviceType) 

# External Links
$(_doc_external("Sys/PetscDeviceSetDefaultDeviceType"))
"""
function PetscDeviceSetDefaultDeviceType(petsclib::PetscLibType, type::PetscDeviceType) end

@for_petsc function PetscDeviceSetDefaultDeviceType(petsclib::$UnionPetscLib, type::PetscDeviceType )

    @chk ccall(
               (:PetscDeviceSetDefaultDeviceType, $petsc_library),
               PetscErrorCode,
               (PetscDeviceType,),
               type,
              )


	return nothing
end 

"""
	PetscDeviceInitialize(petsclib::PetscLibType,type::PetscDeviceType) 

# External Links
$(_doc_external("Sys/PetscDeviceInitialize"))
"""
function PetscDeviceInitialize(petsclib::PetscLibType, type::PetscDeviceType) end

@for_petsc function PetscDeviceInitialize(petsclib::$UnionPetscLib, type::PetscDeviceType )

    @chk ccall(
               (:PetscDeviceInitialize, $petsc_library),
               PetscErrorCode,
               (PetscDeviceType,),
               type,
              )


	return nothing
end 

"""
	PetscDeviceMemcpy(petsclib::PetscLibType,dctx::PetscDeviceContext, dest::Cvoid, src::Cvoid, n_std::Csize_t) 

# External Links
$(_doc_external("Sys/PetscDeviceMemcpy"))
"""
function PetscDeviceMemcpy(petsclib::PetscLibType, dctx::PetscDeviceContext, dest::Cvoid, src::Cvoid, n_std::Csize_t) end

@for_petsc function PetscDeviceMemcpy(petsclib::$UnionPetscLib, dctx::PetscDeviceContext, dest::Cvoid, src::Cvoid, n_std::Csize_t )

    @chk ccall(
               (:PetscDeviceMemcpy, $petsc_library),
               PetscErrorCode,
               (PetscDeviceContext, Ptr{Cvoid}, Ptr{Cvoid}, Csize_t),
               dctx, dest, src, n_std,
              )


	return nothing
end 

"""
	PetscDeviceMemset(petsclib::PetscLibType,dctx::PetscDeviceContext, ptr::Cvoid, v::PetscInt, n_std::Csize_t) 

# External Links
$(_doc_external("Sys/PetscDeviceMemset"))
"""
function PetscDeviceMemset(petsclib::PetscLibType, dctx::PetscDeviceContext, ptr::Cvoid, v::PetscInt, n_std::Csize_t) end

@for_petsc function PetscDeviceMemset(petsclib::$UnionPetscLib, dctx::PetscDeviceContext, ptr::Cvoid, v::$PetscInt, n_std::Csize_t )

    @chk ccall(
               (:PetscDeviceMemset, $petsc_library),
               PetscErrorCode,
               (PetscDeviceContext, Ptr{Cvoid}, $PetscInt, Csize_t),
               dctx, ptr, v, n_std,
              )


	return nothing
end 

"""
	PetscDeviceFinalizePackage(petsclib::PetscLibType) 
This function cleans up all components of the `PetscDevice`
package. It is called from `PetscFinalize()`.

-seealso: `PetscFinalize()`, `PetscDeviceInitializePackage()`

# External Links
$(_doc_external("Sys/PetscDeviceFinalizePackage"))
"""
function PetscDeviceFinalizePackage(petsclib::PetscLibType) end

@for_petsc function PetscDeviceFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscDeviceFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscDeviceInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `PetscDevice`
package. It is called on the first call to `PetscDeviceContextCreate()` or
`PetscDeviceCreate()` when using shared or static libraries.

Level: developer

-seealso: `PetscInitialize()`, `PetscDeviceFinalizePackage()`, `PetscDeviceContextCreate()`,
`PetscDeviceCreate()`

# External Links
$(_doc_external("Sys/PetscDeviceInitializePackage"))
"""
function PetscDeviceInitializePackage(petsclib::PetscLibType) end

@for_petsc function PetscDeviceInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscDeviceInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscDeviceContextGetCurrentContext(petsclib::PetscLibType,dctx::PetscDeviceContext) 

# External Links
$(_doc_external("Sys/PetscDeviceContextGetCurrentContext"))
"""
function PetscDeviceContextGetCurrentContext(petsclib::PetscLibType, dctx::PetscDeviceContext) end

@for_petsc function PetscDeviceContextGetCurrentContext(petsclib::$UnionPetscLib, dctx::PetscDeviceContext )

    @chk ccall(
               (:PetscDeviceContextGetCurrentContext, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscDeviceContext},),
               dctx,
              )


	return nothing
end 

"""
	PetscDeviceContextSetCurrentContext(petsclib::PetscLibType,dctx::PetscDeviceContext) 

# External Links
$(_doc_external("Sys/PetscDeviceContextSetCurrentContext"))
"""
function PetscDeviceContextSetCurrentContext(petsclib::PetscLibType, dctx::PetscDeviceContext) end

@for_petsc function PetscDeviceContextSetCurrentContext(petsclib::$UnionPetscLib, dctx::PetscDeviceContext )

    @chk ccall(
               (:PetscDeviceContextSetCurrentContext, $petsc_library),
               PetscErrorCode,
               (PetscDeviceContext,),
               dctx,
              )


	return nothing
end 

"""
	dctx::PetscDeviceContext = PetscDeviceContextCreate(petsclib::PetscLibType) 

# External Links
$(_doc_external("Sys/PetscDeviceContextCreate"))
"""
function PetscDeviceContextCreate(petsclib::PetscLibType) end

@for_petsc function PetscDeviceContextCreate(petsclib::$UnionPetscLib)
	dctx_ = Ref{PetscDeviceContext}()

    @chk ccall(
               (:PetscDeviceContextCreate, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscDeviceContext},),
               dctx_,
              )

	dctx = dctx_[]

	return dctx
end 

"""
	PetscDeviceContextDestroy(petsclib::PetscLibType,dctx::PetscDeviceContext) 

# External Links
$(_doc_external("Sys/PetscDeviceContextDestroy"))
"""
function PetscDeviceContextDestroy(petsclib::PetscLibType, dctx::PetscDeviceContext) end

@for_petsc function PetscDeviceContextDestroy(petsclib::$UnionPetscLib, dctx::PetscDeviceContext )

    @chk ccall(
               (:PetscDeviceContextDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscDeviceContext},),
               dctx,
              )


	return nothing
end 

"""
	PetscDeviceContextSetStreamType(petsclib::PetscLibType,dctx::PetscDeviceContext, type::PetscStreamType) 

# External Links
$(_doc_external("Sys/PetscDeviceContextSetStreamType"))
"""
function PetscDeviceContextSetStreamType(petsclib::PetscLibType, dctx::PetscDeviceContext, type::PetscStreamType) end

@for_petsc function PetscDeviceContextSetStreamType(petsclib::$UnionPetscLib, dctx::PetscDeviceContext, type::PetscStreamType )

    @chk ccall(
               (:PetscDeviceContextSetStreamType, $petsc_library),
               PetscErrorCode,
               (PetscDeviceContext, PetscStreamType),
               dctx, type,
              )


	return nothing
end 

"""
	type::PetscStreamType = PetscDeviceContextGetStreamType(petsclib::PetscLibType,dctx::PetscDeviceContext) 

# External Links
$(_doc_external("Sys/PetscDeviceContextGetStreamType"))
"""
function PetscDeviceContextGetStreamType(petsclib::PetscLibType, dctx::PetscDeviceContext) end

@for_petsc function PetscDeviceContextGetStreamType(petsclib::$UnionPetscLib, dctx::PetscDeviceContext )
	type_ = Ref{PetscStreamType}()

    @chk ccall(
               (:PetscDeviceContextGetStreamType, $petsc_library),
               PetscErrorCode,
               (PetscDeviceContext, Ptr{PetscStreamType}),
               dctx, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	PetscDeviceContextSetDevice(petsclib::PetscLibType,dctx::PetscDeviceContext, device::PetscDevice) 

# External Links
$(_doc_external("Sys/PetscDeviceContextSetDevice"))
"""
function PetscDeviceContextSetDevice(petsclib::PetscLibType, dctx::PetscDeviceContext, device::PetscDevice) end

@for_petsc function PetscDeviceContextSetDevice(petsclib::$UnionPetscLib, dctx::PetscDeviceContext, device::PetscDevice )

    @chk ccall(
               (:PetscDeviceContextSetDevice, $petsc_library),
               PetscErrorCode,
               (PetscDeviceContext, PetscDevice),
               dctx, device,
              )


	return nothing
end 

"""
	PetscDeviceContextGetDevice(petsclib::PetscLibType,dctx::PetscDeviceContext, device::PetscDevice) 

# External Links
$(_doc_external("Sys/PetscDeviceContextGetDevice"))
"""
function PetscDeviceContextGetDevice(petsclib::PetscLibType, dctx::PetscDeviceContext, device::PetscDevice) end

@for_petsc function PetscDeviceContextGetDevice(petsclib::$UnionPetscLib, dctx::PetscDeviceContext, device::PetscDevice )

    @chk ccall(
               (:PetscDeviceContextGetDevice, $petsc_library),
               PetscErrorCode,
               (PetscDeviceContext, Ptr{PetscDevice}),
               dctx, device,
              )


	return nothing
end 

"""
	type::PetscDeviceType = PetscDeviceContextGetDeviceType(petsclib::PetscLibType,dctx::PetscDeviceContext) 

# External Links
$(_doc_external("Sys/PetscDeviceContextGetDeviceType"))
"""
function PetscDeviceContextGetDeviceType(petsclib::PetscLibType, dctx::PetscDeviceContext) end

@for_petsc function PetscDeviceContextGetDeviceType(petsclib::$UnionPetscLib, dctx::PetscDeviceContext )
	type_ = Ref{PetscDeviceType}()

    @chk ccall(
               (:PetscDeviceContextGetDeviceType, $petsc_library),
               PetscErrorCode,
               (PetscDeviceContext, Ptr{PetscDeviceType}),
               dctx, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	PetscDeviceContextSetUp(petsclib::PetscLibType,dctx::PetscDeviceContext) 

# External Links
$(_doc_external("Sys/PetscDeviceContextSetUp"))
"""
function PetscDeviceContextSetUp(petsclib::PetscLibType, dctx::PetscDeviceContext) end

@for_petsc function PetscDeviceContextSetUp(petsclib::$UnionPetscLib, dctx::PetscDeviceContext )

    @chk ccall(
               (:PetscDeviceContextSetUp, $petsc_library),
               PetscErrorCode,
               (PetscDeviceContext,),
               dctx,
              )


	return nothing
end 

"""
	dctxdup::PetscDeviceContext = PetscDeviceContextDuplicate(petsclib::PetscLibType,dctx::PetscDeviceContext) 

# External Links
$(_doc_external("Sys/PetscDeviceContextDuplicate"))
"""
function PetscDeviceContextDuplicate(petsclib::PetscLibType, dctx::PetscDeviceContext) end

@for_petsc function PetscDeviceContextDuplicate(petsclib::$UnionPetscLib, dctx::PetscDeviceContext )
	dctxdup_ = Ref{PetscDeviceContext}()

    @chk ccall(
               (:PetscDeviceContextDuplicate, $petsc_library),
               PetscErrorCode,
               (PetscDeviceContext, Ptr{PetscDeviceContext}),
               dctx, dctxdup_,
              )

	dctxdup = dctxdup_[]

	return dctxdup
end 

"""
	idle::PetscBool = PetscDeviceContextQueryIdle(petsclib::PetscLibType,dctx::PetscDeviceContext) 

# External Links
$(_doc_external("Sys/PetscDeviceContextQueryIdle"))
"""
function PetscDeviceContextQueryIdle(petsclib::PetscLibType, dctx::PetscDeviceContext) end

@for_petsc function PetscDeviceContextQueryIdle(petsclib::$UnionPetscLib, dctx::PetscDeviceContext )
	idle_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscDeviceContextQueryIdle, $petsc_library),
               PetscErrorCode,
               (PetscDeviceContext, Ptr{PetscBool}),
               dctx, idle_,
              )

	idle = idle_[]

	return idle
end 

"""
	PetscDeviceContextWaitForContext(petsclib::PetscLibType,dctxa::PetscDeviceContext, dctxb::PetscDeviceContext) 

# External Links
$(_doc_external("Sys/PetscDeviceContextWaitForContext"))
"""
function PetscDeviceContextWaitForContext(petsclib::PetscLibType, dctxa::PetscDeviceContext, dctxb::PetscDeviceContext) end

@for_petsc function PetscDeviceContextWaitForContext(petsclib::$UnionPetscLib, dctxa::PetscDeviceContext, dctxb::PetscDeviceContext )

    @chk ccall(
               (:PetscDeviceContextWaitForContext, $petsc_library),
               PetscErrorCode,
               (PetscDeviceContext, PetscDeviceContext),
               dctxa, dctxb,
              )


	return nothing
end 

"""
	PetscDeviceContextForkWithStreamType(petsclib::PetscLibType,dctx::PetscDeviceContext, stype::PetscStreamType, n::PetscInt, dsub::PetscDeviceContext) 

# External Links
$(_doc_external("Sys/PetscDeviceContextForkWithStreamType"))
"""
function PetscDeviceContextForkWithStreamType(petsclib::PetscLibType, dctx::PetscDeviceContext, stype::PetscStreamType, n::PetscInt, dsub::PetscDeviceContext) end

@for_petsc function PetscDeviceContextForkWithStreamType(petsclib::$UnionPetscLib, dctx::PetscDeviceContext, stype::PetscStreamType, n::$PetscInt, dsub::PetscDeviceContext )

    @chk ccall(
               (:PetscDeviceContextForkWithStreamType, $petsc_library),
               PetscErrorCode,
               (PetscDeviceContext, PetscStreamType, $PetscInt, PetscDeviceContext),
               dctx, stype, n, dsub,
              )


	return nothing
end 

"""
	PetscDeviceContextFork(petsclib::PetscLibType,dctx::PetscDeviceContext, n::PetscInt, dsub::PetscDeviceContext) 

# External Links
$(_doc_external("Sys/PetscDeviceContextFork"))
"""
function PetscDeviceContextFork(petsclib::PetscLibType, dctx::PetscDeviceContext, n::PetscInt, dsub::PetscDeviceContext) end

@for_petsc function PetscDeviceContextFork(petsclib::$UnionPetscLib, dctx::PetscDeviceContext, n::$PetscInt, dsub::PetscDeviceContext )

    @chk ccall(
               (:PetscDeviceContextFork, $petsc_library),
               PetscErrorCode,
               (PetscDeviceContext, $PetscInt, PetscDeviceContext),
               dctx, n, dsub,
              )


	return nothing
end 

"""
	PetscDeviceContextJoin(petsclib::PetscLibType,dctx::PetscDeviceContext, n::PetscInt, joinMode::PetscDeviceContextJoinMode, dsub::PetscDeviceContext) 

# External Links
$(_doc_external("Sys/PetscDeviceContextJoin"))
"""
function PetscDeviceContextJoin(petsclib::PetscLibType, dctx::PetscDeviceContext, n::PetscInt, joinMode::PetscDeviceContextJoinMode, dsub::PetscDeviceContext) end

@for_petsc function PetscDeviceContextJoin(petsclib::$UnionPetscLib, dctx::PetscDeviceContext, n::$PetscInt, joinMode::PetscDeviceContextJoinMode, dsub::PetscDeviceContext )

    @chk ccall(
               (:PetscDeviceContextJoin, $petsc_library),
               PetscErrorCode,
               (PetscDeviceContext, $PetscInt, PetscDeviceContextJoinMode, PetscDeviceContext),
               dctx, n, joinMode, dsub,
              )


	return nothing
end 

"""
	PetscDeviceContextSynchronize(petsclib::PetscLibType,dctx::PetscDeviceContext) 

# External Links
$(_doc_external("Sys/PetscDeviceContextSynchronize"))
"""
function PetscDeviceContextSynchronize(petsclib::PetscLibType, dctx::PetscDeviceContext) end

@for_petsc function PetscDeviceContextSynchronize(petsclib::$UnionPetscLib, dctx::PetscDeviceContext )

    @chk ccall(
               (:PetscDeviceContextSynchronize, $petsc_library),
               PetscErrorCode,
               (PetscDeviceContext,),
               dctx,
              )


	return nothing
end 

"""
	PetscDeviceContextSetFromOptions(petsclib::PetscLibType,comm::MPI_Comm, dctx::PetscDeviceContext) 

# External Links
$(_doc_external("Sys/PetscDeviceContextSetFromOptions"))
"""
function PetscDeviceContextSetFromOptions(petsclib::PetscLibType, comm::MPI_Comm, dctx::PetscDeviceContext) end

@for_petsc function PetscDeviceContextSetFromOptions(petsclib::$UnionPetscLib, comm::MPI_Comm, dctx::PetscDeviceContext )

    @chk ccall(
               (:PetscDeviceContextSetFromOptions, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, PetscDeviceContext),
               comm, dctx,
              )


	return nothing
end 

"""
	PetscDeviceContextView(petsclib::PetscLibType,dctx::PetscDeviceContext, viewer::PetscViewer) 

# External Links
$(_doc_external("Sys/PetscDeviceContextView"))
"""
function PetscDeviceContextView(petsclib::PetscLibType, dctx::PetscDeviceContext, viewer::PetscViewer) end

@for_petsc function PetscDeviceContextView(petsclib::$UnionPetscLib, dctx::PetscDeviceContext, viewer::PetscViewer )

    @chk ccall(
               (:PetscDeviceContextView, $petsc_library),
               PetscErrorCode,
               (PetscDeviceContext, PetscViewer),
               dctx, viewer,
              )


	return nothing
end 

"""
	PetscDeviceContextViewFromOptions(petsclib::PetscLibType,dctx::PetscDeviceContext, obj::PetscObject, name::String) 

# External Links
$(_doc_external("Sys/PetscDeviceContextViewFromOptions"))
"""
function PetscDeviceContextViewFromOptions(petsclib::PetscLibType, dctx::PetscDeviceContext, obj::PetscObject, name::String) end

@for_petsc function PetscDeviceContextViewFromOptions(petsclib::$UnionPetscLib, dctx::PetscDeviceContext, obj::PetscObject, name::String )

    @chk ccall(
               (:PetscDeviceContextViewFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscDeviceContext, PetscObject, Ptr{Cchar}),
               dctx, obj, name,
              )


	return nothing
end 

"""
	PetscDeviceContextGetStreamHandle(petsclib::PetscLibType,dctx::PetscDeviceContext, handle::Cvoid) 

# External Links
$(_doc_external("Sys/PetscDeviceContextGetStreamHandle"))
"""
function PetscDeviceContextGetStreamHandle(petsclib::PetscLibType, dctx::PetscDeviceContext, handle::Cvoid) end

@for_petsc function PetscDeviceContextGetStreamHandle(petsclib::$UnionPetscLib, dctx::PetscDeviceContext, handle::Cvoid )

    @chk ccall(
               (:PetscDeviceContextGetStreamHandle, $petsc_library),
               PetscErrorCode,
               (PetscDeviceContext, Cvoid),
               dctx, handle,
              )


	return nothing
end 

