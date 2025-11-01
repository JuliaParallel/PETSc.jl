# autodefined type arguments for class ------
mutable struct _n_PetscContainer end
const PetscContainer = Ptr{_n_PetscContainer}

# -------------------------------------------------------
"""
	PetscContainerGetPointer(petsclib::PetscLibType,obj::PetscContainer, ptr::PeCtx) 
Gets the pointer value contained in the container that was provided with `PetscContainerSetPointer()`

Not Collective, No Fortran Support

Input Parameter:
- `obj` - the object created with `PetscContainerCreate()`

Output Parameter:
- `ptr` - the pointer value

Level: advanced

-seealso: `PetscContainerCreate()`, `PetscContainerDestroy()`, `PetscObject`,
`PetscContainerSetPointer()`, `PetscObjectContainerCompose()`, `PetscObjectContainerQuery()`

# External Links
$(_doc_external("Sys/PetscContainerGetPointer"))
"""
function PetscContainerGetPointer(petsclib::PetscLibType, obj::PetscContainer, ptr::PeCtx) end

@for_petsc function PetscContainerGetPointer(petsclib::$UnionPetscLib, obj::PetscContainer, ptr::PeCtx )

    @chk ccall(
               (:PetscContainerGetPointer, $petsc_library),
               PetscErrorCode,
               (PetscContainer, PeCtx),
               obj, ptr,
              )


	return nothing
end 

"""
	PetscContainerSetPointer(petsclib::PetscLibType,obj::PetscContainer, ptr::Cvoid) 
Sets the pointer value contained in the container.

Logically Collective, No Fortran Support

Input Parameters:
- `obj` - the object created with `PetscContainerCreate()`
- `ptr` - the pointer value

Level: advanced

-seealso: `PetscContainerCreate()`, `PetscContainerDestroy()`, `PetscObjectCompose()`, `PetscObjectQuery()`, `PetscObject`,
`PetscContainerGetPointer()`, `PetscObjectContainerCompose()`, `PetscObjectContainerQuery()`

# External Links
$(_doc_external("Sys/PetscContainerSetPointer"))
"""
function PetscContainerSetPointer(petsclib::PetscLibType, obj::PetscContainer, ptr::Cvoid) end

@for_petsc function PetscContainerSetPointer(petsclib::$UnionPetscLib, obj::PetscContainer, ptr::Cvoid )

    @chk ccall(
               (:PetscContainerSetPointer, $petsc_library),
               PetscErrorCode,
               (PetscContainer, Ptr{Cvoid}),
               obj, ptr,
              )


	return nothing
end 

"""
	PetscContainerDestroy(petsclib::PetscLibType,obj::PetscContainer) 
Destroys a PETSc container object.

Collective, No Fortran Support

Input Parameter:
- `obj` - an object that was created with `PetscContainerCreate()`

Level: advanced

-seealso: `PetscContainerCreate()`, `PetscContainerSetCtxDestroy()`, `PetscObject`, `PetscObjectContainerCompose()`, `PetscObjectContainerQuery()`

# External Links
$(_doc_external("Sys/PetscContainerDestroy"))
"""
function PetscContainerDestroy(petsclib::PetscLibType, obj::PetscContainer) end

@for_petsc function PetscContainerDestroy(petsclib::$UnionPetscLib, obj::PetscContainer )

    @chk ccall(
               (:PetscContainerDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscContainer},),
               obj,
              )


	return nothing
end 

"""
	PetscContainerSetCtxDestroy(petsclib::PetscLibType,obj::PetscContainer, des::PetscCtxDestroyFn) 
Sets the destroy function for the data provided to the `PetscContainer` with `PetscContainerSetPointer()`

Logically Collective, No Fortran Support

Input Parameters:
- `obj` - an object that was created with `PetscContainerCreate()`
- `des` - name of the ctx destroy function, see `PetscCtxDestroyFn` for its calling sequence

Level: advanced

-seealso: `PetscContainerDestroy()`, `PetscContainerUserDestroyDefault()`, `PetscMalloc()`, `PetscMalloc1()`, `PetscCalloc()`, `PetscCalloc1()`, `PetscObject`,
`PetscObjectContainerCompose()`, `PetscObjectContainerQuery()`

# External Links
$(_doc_external("Sys/PetscContainerSetCtxDestroy"))
"""
function PetscContainerSetCtxDestroy(petsclib::PetscLibType, obj::PetscContainer, des::PetscCtxDestroyFn) end

@for_petsc function PetscContainerSetCtxDestroy(petsclib::$UnionPetscLib, obj::PetscContainer, des::PetscCtxDestroyFn )

    @chk ccall(
               (:PetscContainerSetCtxDestroy, $petsc_library),
               PetscErrorCode,
               (PetscContainer, Ptr{PetscCtxDestroyFn}),
               obj, des,
              )


	return nothing
end 

"""
	container::PetscContainer = PetscContainerCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Creates a PETSc object that has room to hold a single pointer.

Collective, No Fortran Support

Input Parameter:
- `comm` - MPI communicator that shares the object

Output Parameter:
- `container` - the container created

Level: advanced

-seealso: `PetscContainerDestroy()`, `PetscContainerSetPointer()`, `PetscContainerGetPointer()`, `PetscObjectCompose()`, `PetscObjectQuery()`,
`PetscContainerSetCtxDestroy()`, `PetscObject`, `PetscObjectContainerCompose()`, `PetscObjectContainerQuery()`

# External Links
$(_doc_external("Sys/PetscContainerCreate"))
"""
function PetscContainerCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function PetscContainerCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	container_ = Ref{PetscContainer}()

    @chk ccall(
               (:PetscContainerCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscContainer}),
               comm, container_,
              )

	container = container_[]

	return container
end 

