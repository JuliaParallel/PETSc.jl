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
function PetscContainerCreate(petsclib::PetscLibType, comm::MPI_Comm)
    error("PetscContainerCreate: no generated method for these argument types")
end

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

"""
	PetscContainerDestroy(petsclib::PetscLibType,obj::Union{PetscContainer, Ref{PetscContainer}}) 
Destroys a PETSc container object.

Collective, No Fortran Support

Input Parameter:
- `obj` - an object that was created with `PetscContainerCreate()`

Level: advanced

-seealso: `PetscContainerCreate()`, `PetscContainerSetCtxDestroy()`, `PetscObject`, `PetscObjectContainerCompose()`, `PetscObjectContainerQuery()`

# External Links
$(_doc_external("Sys/PetscContainerDestroy"))
"""
function PetscContainerDestroy(petsclib::PetscLibType, obj::Union{PetscContainer, Ref{PetscContainer}})
    error("PetscContainerDestroy: no generated method for these argument types")
end

@for_petsc function PetscContainerDestroy(petsclib::$UnionPetscLib, obj::Union{PetscContainer, Ref{PetscContainer}} )
	obj_ = obj isa Base.RefValue ? obj : Ref{PetscContainer}(obj)

    @chk ccall(
               (:PetscContainerDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscContainer},),
               obj_,
              )


	return nothing
end 

"""
	ptr::Ptr{Cvoid} = PetscContainerGetPointer(petsclib::PetscLibType,obj::PetscContainer) 
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
function PetscContainerGetPointer(petsclib::PetscLibType, obj::PetscContainer)
    error("PetscContainerGetPointer: no generated method for these argument types")
end

@for_petsc function PetscContainerGetPointer(petsclib::$UnionPetscLib, obj::PetscContainer )
	ptr_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:PetscContainerGetPointer, $petsc_library),
               PetscErrorCode,
               (PetscContainer, Ptr{Cvoid}),
               obj, ptr_,
              )

	ptr = ptr_[]

	return ptr
end 

"""
	PetscContainerSetCtxDestroy(petsclib::PetscLibType,obj::PetscContainer, des::Ptr{Cvoid}) 
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
function PetscContainerSetCtxDestroy(petsclib::PetscLibType, obj::PetscContainer, des::Ptr{Cvoid})
    error("PetscContainerSetCtxDestroy: no generated method for these argument types")
end

@for_petsc function PetscContainerSetCtxDestroy(petsclib::$UnionPetscLib, obj::PetscContainer, des::Ptr{Cvoid} )

    @chk ccall(
               (:PetscContainerSetCtxDestroy, $petsc_library),
               PetscErrorCode,
               (PetscContainer, Ptr{Cvoid}),
               obj, des,
              )


	return nothing
end 

"""
	PetscContainerSetPointer(petsclib::PetscLibType,obj::PetscContainer, ptr::Ptr{Cvoid}) 
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
function PetscContainerSetPointer(petsclib::PetscLibType, obj::PetscContainer, ptr::Ptr{Cvoid})
    error("PetscContainerSetPointer: no generated method for these argument types")
end

@for_petsc function PetscContainerSetPointer(petsclib::$UnionPetscLib, obj::PetscContainer, ptr::Ptr{Cvoid} )

    @chk ccall(
               (:PetscContainerSetPointer, $petsc_library),
               PetscErrorCode,
               (PetscContainer, Ptr{Cvoid}),
               obj, ptr,
              )


	return nothing
end 

