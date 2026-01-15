# autodefined type arguments for class ------
mutable struct _n_PetscOptionItems end
const PetscOptionItems = Ptr{_n_PetscOptionItems}

# -------------------------------------------------------
"""
	PetscObjectSAWsTakeAccess(petsclib::PetscLibType,obj::PetscObject) 
Take access of the data fields that have been published to SAWs
by a `PetscObject` so their values may  be changed in the computation

Collective

Input Parameter:
- `obj` - the `PetscObject` variable. This must be cast with a (`PetscObject`), for example, `PetscObjectSAWSTakeAccess`((`PetscObject`)mat);

Level: advanced

-seealso: `PetscObjectSetName()`, `PetscObjectSAWsViewOff()`, `PetscObjectSAWsGrantAccess()`

# External Links
$(_doc_external("Sys/PetscObjectSAWsTakeAccess"))
"""
function PetscObjectSAWsTakeAccess(petsclib::PetscLibType, obj::PetscObject) end

@for_petsc function PetscObjectSAWsTakeAccess(petsclib::$UnionPetscLib, obj::PetscObject )

    @chk ccall(
               (:PetscObjectSAWsTakeAccess, $petsc_library),
               PetscErrorCode,
               (PetscObject,),
               obj,
              )


	return nothing
end 

"""
	PetscObjectSAWsGrantAccess(petsclib::PetscLibType,obj::PetscObject) 
Grants access of the data fields that have been published to
SAWs called when the changes made during `PetscObjectSAWsTakeAccess()` are complete.

Collective

Input Parameter:
- `obj` - the `PetscObject` variable. This must be cast with a (`PetscObject`), for example, `PetscObjectSAWSRestoreAccess`((`PetscObject`)mat);

Level: advanced

-seealso: `PetscObjectSetName()`, `PetscObjectSAWsViewOff()`, `PetscObjectSAWsTakeAccess()`

# External Links
$(_doc_external("Sys/PetscObjectSAWsGrantAccess"))
"""
function PetscObjectSAWsGrantAccess(petsclib::PetscLibType, obj::PetscObject) end

@for_petsc function PetscObjectSAWsGrantAccess(petsclib::$UnionPetscLib, obj::PetscObject )

    @chk ccall(
               (:PetscObjectSAWsGrantAccess, $petsc_library),
               PetscErrorCode,
               (PetscObject,),
               obj,
              )


	return nothing
end 

"""
	PetscObjectSAWsBlock(petsclib::PetscLibType,obj::PetscObject) 
Blocks the object if `PetscObjectSAWsSetBlock()` has been called

Collective

Input Parameter:
- `obj` - the PETSc variable

Level: advanced

-seealso: `PetscObjectSetName()`, `PetscObjectSAWsViewOff()`, `PetscObjectSAWsSetBlock()`, `PetscSAWsBlock()`

# External Links
$(_doc_external("Sys/PetscObjectSAWsBlock"))
"""
function PetscObjectSAWsBlock(petsclib::PetscLibType, obj::PetscObject) end

@for_petsc function PetscObjectSAWsBlock(petsclib::$UnionPetscLib, obj::PetscObject )

    @chk ccall(
               (:PetscObjectSAWsBlock, $petsc_library),
               PetscErrorCode,
               (PetscObject,),
               obj,
              )


	return nothing
end 

"""
	PetscObjectSAWsSetBlock(petsclib::PetscLibType,obj::PetscObject, flg::PetscBool) 
Sets whether an object will block at `PetscObjectSAWsBlock()`

Collective

Input Parameters:
- `obj` - the PETSc variable
- `flg` - whether it should block

Level: advanced

-seealso: `PetscObjectSetName()`, `PetscObjectSAWsViewOff()`, `PetscObjectSAWsBlock()`, `PetscSAWsBlock()`

# External Links
$(_doc_external("Sys/PetscObjectSAWsSetBlock"))
"""
function PetscObjectSAWsSetBlock(petsclib::PetscLibType, obj::PetscObject, flg::PetscBool) end

@for_petsc function PetscObjectSAWsSetBlock(petsclib::$UnionPetscLib, obj::PetscObject, flg::PetscBool )

    @chk ccall(
               (:PetscObjectSAWsSetBlock, $petsc_library),
               PetscErrorCode,
               (PetscObject, PetscBool),
               obj, flg,
              )


	return nothing
end 

"""
	PetscObjectSAWsViewOff(petsclib::PetscLibType,obj::PetscObject) 

# External Links
$(_doc_external("Sys/PetscObjectSAWsViewOff"))
"""
function PetscObjectSAWsViewOff(petsclib::PetscLibType, obj::PetscObject) end

@for_petsc function PetscObjectSAWsViewOff(petsclib::$UnionPetscLib, obj::PetscObject )

    @chk ccall(
               (:PetscObjectSAWsViewOff, $petsc_library),
               PetscErrorCode,
               (PetscObject,),
               obj,
              )


	return nothing
end 

"""
	PetscObjectViewSAWs(petsclib::PetscLibType,obj::PetscObject, viewer::PetscViewer) 
View the base portion of any object with an SAWs viewer

Collective

Input Parameters:
- `obj`    - the `PetscObject` variable. It must be cast with a (`PetscObject`), for example, `PetscObjectSetName`((`PetscObject`)mat,name);
- `viewer` - the SAWs viewer

Level: advanced

-seealso: [](sec_viewers), `PetscViewer`, `PetscObject`, `PetscObjectSetName()`

# External Links
$(_doc_external("Sys/PetscObjectViewSAWs"))
"""
function PetscObjectViewSAWs(petsclib::PetscLibType, obj::PetscObject, viewer::PetscViewer) end

@for_petsc function PetscObjectViewSAWs(petsclib::$UnionPetscLib, obj::PetscObject, viewer::PetscViewer )

    @chk ccall(
               (:PetscObjectViewSAWs, $petsc_library),
               PetscErrorCode,
               (PetscObject, PetscViewer),
               obj, viewer,
              )


	return nothing
end 

"""
	type::String = PetscObjectGetType(petsclib::PetscLibType,obj::PetscObject) 
Gets the object type of any `PetscObject`.

Not Collective

Input Parameter:
- `obj` - any PETSc object, for example a `Vec`, `Mat` or `KSP`. It must be cast with a (`PetscObject`), for example,
`PetscObjectGetType`((`PetscObject`)mat,&type);

Output Parameter:
- `type` - the object type, for example, `MATSEQAIJ`

Level: advanced

-seealso: `PetscObject`, `PetscClassId`, `PetscObjectGetClassName()`, `PetscObjectGetClassId()`

# External Links
$(_doc_external("Sys/PetscObjectGetType"))
"""
function PetscObjectGetType(petsclib::PetscLibType, obj::PetscObject) end

@for_petsc function PetscObjectGetType(petsclib::$UnionPetscLib, obj::PetscObject )
	type_ = Ref{Ptr{Cchar}}()

    @chk ccall(
               (:PetscObjectGetType, $petsc_library),
               PetscErrorCode,
               (PetscObject, Ptr{Ptr{Cchar}}),
               obj, type_,
              )

	type = unsafe_wrap(Array, type_[], VecGetLocalSize(petsclib, x); own = false)

	return type
end 

"""
	PetscObjectsDump(petsclib::PetscLibType,fd::Libc.FILE, all::PetscBool) 
Prints all the currently existing objects.

Input Parameters:
- `fd`  - file pointer
- `all` - by default only tries to display objects created explicitly by the user, if all is `PETSC_TRUE` then lists all outstanding objects

Options Database Key:
- `-objects_dump <all>` - print information about all the objects that exist at the end of the programs run

Level: advanced

-seealso: `PetscObject`

# External Links
$(_doc_external("Sys/PetscObjectsDump"))
"""
function PetscObjectsDump(petsclib::PetscLibType, fd::Libc.FILE, all::PetscBool) end

@for_petsc function PetscObjectsDump(petsclib::$UnionPetscLib, fd::Libc.FILE, all::PetscBool )

    @chk ccall(
               (:PetscObjectsDump, $petsc_library),
               PetscErrorCode,
               (Ptr{Libc.FILE}, PetscBool),
               fd, all,
              )


	return nothing
end 

"""
	PetscObjectsView(petsclib::PetscLibType,viewer::PetscViewer) 
Prints the currently existing objects.

Logically Collective

Input Parameter:
- `viewer` - must be an `PETSCVIEWERASCII` viewer

Level: advanced

-seealso: `PetscObject`

# External Links
$(_doc_external("Sys/PetscObjectsView"))
"""
function PetscObjectsView(petsclib::PetscLibType, viewer::PetscViewer) end

@for_petsc function PetscObjectsView(petsclib::$UnionPetscLib, viewer::PetscViewer )

    @chk ccall(
               (:PetscObjectsView, $petsc_library),
               PetscErrorCode,
               (PetscViewer,),
               viewer,
              )


	return nothing
end 

"""
	PetscObjectsGetObject(petsclib::PetscLibType,name::String, obj::PetscObject, classname::String) 
Get a pointer to a named object

Not Collective

Input Parameter:
- `name` - the name of an object

Output Parameters:
- `obj`       - the object or `NULL` if there is no object, optional, pass in `NULL` if not needed
- `classname` - the name of the class of the object, optional, pass in `NULL` if not needed

Level: advanced

-seealso: `PetscObject`

# External Links
$(_doc_external("Sys/PetscObjectsGetObject"))
"""
function PetscObjectsGetObject(petsclib::PetscLibType, name::String, obj::PetscObject, classname::String) end

@for_petsc function PetscObjectsGetObject(petsclib::$UnionPetscLib, name::String, obj::PetscObject, classname::String )
	classname_ = Ref(pointer(classname))

    @chk ccall(
               (:PetscObjectsGetObject, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{PetscObject}, Ptr{Ptr{Cchar}}),
               name, obj, classname_,
              )


	return nothing
end 

"""
	PetscObjectSetPrintedOptions(petsclib::PetscLibType,obj::PetscObject) 
indicate to an object that it should behave as if it has already printed the help for its options so it will not display the help message

Input Parameter:
- `obj` - the `PetscObject`

Level: developer

-seealso: `PetscOptionsInsert()`, `PetscObject`

# External Links
$(_doc_external("Sys/PetscObjectSetPrintedOptions"))
"""
function PetscObjectSetPrintedOptions(petsclib::PetscLibType, obj::PetscObject) end

@for_petsc function PetscObjectSetPrintedOptions(petsclib::$UnionPetscLib, obj::PetscObject )

    @chk ccall(
               (:PetscObjectSetPrintedOptions, $petsc_library),
               PetscErrorCode,
               (PetscObject,),
               obj,
              )


	return nothing
end 

"""
	PetscObjectInheritPrintedOptions(petsclib::PetscLibType,pobj::PetscObject, obj::PetscObject) 
If the child object is not on the MPI rank 0 process of the parent object and the child is sequential then the child gets it set.

Input Parameters:
- `pobj` - the parent object
- `obj`  - the `PetscObject`

Level: developer

-seealso: `PetscOptionsInsert()`, `PetscObjectSetPrintedOptions()`, `PetscObject`

# External Links
$(_doc_external("Sys/PetscObjectInheritPrintedOptions"))
"""
function PetscObjectInheritPrintedOptions(petsclib::PetscLibType, pobj::PetscObject, obj::PetscObject) end

@for_petsc function PetscObjectInheritPrintedOptions(petsclib::$UnionPetscLib, pobj::PetscObject, obj::PetscObject )

    @chk ccall(
               (:PetscObjectInheritPrintedOptions, $petsc_library),
               PetscErrorCode,
               (PetscObject, PetscObject),
               pobj, obj,
              )


	return nothing
end 

"""
	PetscObjectAddOptionsHandler(petsclib::PetscLibType,obj::PetscObject, handle::external, destroy::external, ctx::Cvoid) 
Adds an additional function to check for options when `XXXSetFromOptions()` is called.

Not Collective

Input Parameters:
- `obj`     - the PETSc object
- `handle`  - function that checks for options
- `destroy` - function to destroy `ctx` if provided
- `ctx`     - optional context for check function

Calling sequence of `handle`:
- `obj`                - the PETSc object
- `PetscOptionsObject` - the `PetscOptionItems` object
- `ctx`                - optional context for `handle`

Calling sequence of `destroy`:
- `obj` - the PETSc object
- `ctx` - optional context for `handle`

Level: developer

-seealso: `KSPSetFromOptions()`, `PCSetFromOptions()`, `SNESSetFromOptions()`, `PetscObjectProcessOptionsHandlers()`, `PetscObjectDestroyOptionsHandlers()`,
`PetscObject`

# External Links
$(_doc_external("Sys/PetscObjectAddOptionsHandler"))
"""
function PetscObjectAddOptionsHandler(petsclib::PetscLibType, obj::PetscObject, handle::external, destroy::external, ctx::Cvoid) end

@for_petsc function PetscObjectAddOptionsHandler(petsclib::$UnionPetscLib, obj::PetscObject, handle::external, destroy::external, ctx::Cvoid )

    @chk ccall(
               (:PetscObjectAddOptionsHandler, $petsc_library),
               PetscErrorCode,
               (PetscObject, external, external, Ptr{Cvoid}),
               obj, handle, destroy, ctx,
              )


	return nothing
end 

"""
	PetscObjectProcessOptionsHandlers(petsclib::PetscLibType,obj::PetscObject, PetscOptionsObject::PetscOptionItems) 
Calls all the options handlers attached to an object

Not Collective

Input Parameters:
- `obj`                - the PETSc object
- `PetscOptionsObject` - the options context

Level: developer

-seealso: `KSPSetFromOptions()`, `PCSetFromOptions()`, `SNESSetFromOptions()`, `PetscObjectAddOptionsHandler()`, `PetscObjectDestroyOptionsHandlers()`,
`PetscObject`

# External Links
$(_doc_external("Sys/PetscObjectProcessOptionsHandlers"))
"""
function PetscObjectProcessOptionsHandlers(petsclib::PetscLibType, obj::PetscObject, PetscOptionsObject::PetscOptionItems) end

@for_petsc function PetscObjectProcessOptionsHandlers(petsclib::$UnionPetscLib, obj::PetscObject, PetscOptionsObject::PetscOptionItems )

    @chk ccall(
               (:PetscObjectProcessOptionsHandlers, $petsc_library),
               PetscErrorCode,
               (PetscObject, PetscOptionItems),
               obj, PetscOptionsObject,
              )


	return nothing
end 

"""
	PetscObjectDestroyOptionsHandlers(petsclib::PetscLibType,obj::PetscObject) 
Destroys all the option handlers attached to an object

Not Collective

Input Parameter:
- `obj` - the PETSc object

Level: developer

-seealso: `KSPSetFromOptions()`, `PCSetFromOptions()`, `SNESSetFromOptions()`, `PetscObjectAddOptionsHandler()`, `PetscObjectProcessOptionsHandlers()`,
`PetscObject`

# External Links
$(_doc_external("Sys/PetscObjectDestroyOptionsHandlers"))
"""
function PetscObjectDestroyOptionsHandlers(petsclib::PetscLibType, obj::PetscObject) end

@for_petsc function PetscObjectDestroyOptionsHandlers(petsclib::$UnionPetscLib, obj::PetscObject )

    @chk ccall(
               (:PetscObjectDestroyOptionsHandlers, $petsc_library),
               PetscErrorCode,
               (PetscObject,),
               obj,
              )


	return nothing
end 

"""
	PetscObjectReference(petsclib::PetscLibType,obj::PetscObject) 
Indicates to a `PetscObject` that it is being
referenced by another `PetscObject`. This increases the reference
count for that object by one.

Logically Collective

Input Parameter:
- `obj` - the PETSc object. This must be cast with (`PetscObject`), for example, `PetscObjectReference`((`PetscObject`)mat);

Level: advanced

-seealso: `PetscObjectCompose()`, `PetscObjectDereference()`, `PetscObject`

# External Links
$(_doc_external("Sys/PetscObjectReference"))
"""
function PetscObjectReference(petsclib::PetscLibType, obj::PetscObject) end

@for_petsc function PetscObjectReference(petsclib::$UnionPetscLib, obj::PetscObject )

    @chk ccall(
               (:PetscObjectReference, $petsc_library),
               PetscErrorCode,
               (PetscObject,),
               obj,
              )


	return nothing
end 

"""
	cnt::PetscInt = PetscObjectGetReference(petsclib::PetscLibType,obj::PetscObject) 
Gets the current reference count for a PETSc object.

Not Collective

Input Parameter:
- `obj` - the PETSc object; this must be cast with (`PetscObject`), for example,
`PetscObjectGetReference`((`PetscObject`)mat,&cnt); `obj` cannot be `NULL`

Output Parameter:
- `cnt` - the reference count

Level: advanced

-seealso: `PetscObjectCompose()`, `PetscObjectDereference()`, `PetscObjectReference()`, `PetscObject`

# External Links
$(_doc_external("Sys/PetscObjectGetReference"))
"""
function PetscObjectGetReference(petsclib::PetscLibType, obj::PetscObject) end

@for_petsc function PetscObjectGetReference(petsclib::$UnionPetscLib, obj::PetscObject )
	cnt_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscObjectGetReference, $petsc_library),
               PetscErrorCode,
               (PetscObject, Ptr{$PetscInt}),
               obj, cnt_,
              )

	cnt = cnt_[]

	return cnt
end 

"""
	PetscObjectDereference(petsclib::PetscLibType,obj::PetscObject) 
Indicates to any `PetscObject` that it is being
referenced by one less `PetscObject`. This decreases the reference
count for that object by one.

Collective on `obj` if reference reaches 0 otherwise Logically Collective

Input Parameter:
- `obj` - the PETSc object; this must be cast with (`PetscObject`), for example,
`PetscObjectDereference`((`PetscObject`)mat);

Level: advanced

-seealso: `PetscObjectCompose()`, `PetscObjectReference()`, `PetscObjectDestroy()`, `PetscObject`

# External Links
$(_doc_external("Sys/PetscObjectDereference"))
"""
function PetscObjectDereference(petsclib::PetscLibType, obj::PetscObject) end

@for_petsc function PetscObjectDereference(petsclib::$UnionPetscLib, obj::PetscObject )

    @chk ccall(
               (:PetscObjectDereference, $petsc_library),
               PetscErrorCode,
               (PetscObject,),
               obj,
              )


	return nothing
end 

"""
	PetscObjectRemoveReference(petsclib::PetscLibType,obj::PetscObject, name::String) 

# External Links
$(_doc_external("Sys/PetscObjectRemoveReference"))
"""
function PetscObjectRemoveReference(petsclib::PetscLibType, obj::PetscObject, name::String) end

@for_petsc function PetscObjectRemoveReference(petsclib::$UnionPetscLib, obj::PetscObject, name::String )

    @chk ccall(
               (:PetscObjectRemoveReference, $petsc_library),
               PetscErrorCode,
               (PetscObject, Ptr{Cchar}),
               obj, name,
              )


	return nothing
end 

"""
	PetscObjectCompose(petsclib::PetscLibType,obj::PetscObject, name::String, ptr::PetscObject) 
Associates another PETSc object with a given PETSc object.

Not Collective

Input Parameters:
- `obj`  - the PETSc object; this must be cast with (`PetscObject`), for example,
`PetscObjectCompose`((`PetscObject`)mat,...);
- `name` - name associated with the child object
- `ptr`  - the other PETSc object to associate with the PETSc object; this must also be
cast with (`PetscObject`)

Level: advanced

-seealso: `PetscObjectQuery()`, `PetscContainerCreate()`, `PetscObjectComposeFunction()`, `PetscObjectQueryFunction()`, `PetscContainer`,
`PetscContainerSetPointer()`, `PetscObject`, `PetscObjectContainerCompose()`

# External Links
$(_doc_external("Sys/PetscObjectCompose"))
"""
function PetscObjectCompose(petsclib::PetscLibType, obj::PetscObject, name::String, ptr::PetscObject) end

@for_petsc function PetscObjectCompose(petsclib::$UnionPetscLib, obj::PetscObject, name::String, ptr::PetscObject )

    @chk ccall(
               (:PetscObjectCompose, $petsc_library),
               PetscErrorCode,
               (PetscObject, Ptr{Cchar}, PetscObject),
               obj, name, ptr,
              )


	return nothing
end 

"""
	PetscObjectQuery(petsclib::PetscLibType,obj::PetscObject, name::String, ptr::PetscObject) 
Gets a PETSc object associated with a given object that was composed with `PetscObjectCompose()`

Not Collective

Input Parameters:
- `obj`  - the PETSc object. It must be cast with a (`PetscObject`), for example,
`PetscObjectCompose`((`PetscObject`)mat,...);
- `name` - name associated with child object
- `ptr`  - the other PETSc object associated with the PETSc object, this must be
cast with (`PetscObject`*)

Level: advanced

-seealso: `PetscObjectCompose()`, `PetscObjectComposeFunction()`, `PetscObjectQueryFunction()`, `PetscContainer`
`PetscContainerGetPointer()`, `PetscObject`

# External Links
$(_doc_external("Sys/PetscObjectQuery"))
"""
function PetscObjectQuery(petsclib::PetscLibType, obj::PetscObject, name::String, ptr::PetscObject) end

@for_petsc function PetscObjectQuery(petsclib::$UnionPetscLib, obj::PetscObject, name::String, ptr::PetscObject )

    @chk ccall(
               (:PetscObjectQuery, $petsc_library),
               PetscErrorCode,
               (PetscObject, Ptr{Cchar}, Ptr{PetscObject}),
               obj, name, ptr,
              )


	return nothing
end 

"""
	has::PetscBool = PetscObjectHasFunction(petsclib::PetscLibType,obj::PetscObject, name::String) 
Query if a function is associated with a given object.

Logically Collective

Input Parameters:
- `obj`  - the PETSc object
- `name` - name associated with the child function

Output Parameter:
- `has` - the boolean value

Level: advanced

-seealso: `PetscObject`, `PetscObjectComposeFunction()`, `PetscObjectQueryFunction()`

# External Links
$(_doc_external("Sys/PetscObjectHasFunction"))
"""
function PetscObjectHasFunction(petsclib::PetscLibType, obj::PetscObject, name::String) end

@for_petsc function PetscObjectHasFunction(petsclib::$UnionPetscLib, obj::PetscObject, name::String )
	has_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscObjectHasFunction, $petsc_library),
               PetscErrorCode,
               (PetscObject, Ptr{Cchar}, Ptr{PetscBool}),
               obj, name, has_,
              )

	has = has_[]

	return has
end 

"""
	PetscObjectContainerCompose(petsclib::PetscLibType,obj::PetscObject, name::String, pointer::Cvoid, destroy::PetscCtxDestroyFn) 
Creates a `PetscContainer`, provides all of its values and composes it with a `PetscObject`

Collective

Input Parameters:
- `obj`     - the `PetscObject`
- `name`    - the name for the composed container
- `pointer` - the pointer to the data
- `destroy` - the routine to destroy the container's data, see `PetscCtxDestroyFn` for its calling sequence; use `PetscCtxDestroyDefault()` if a `PetscFree()` frees the data

Level: advanced

-seealso: `PetscContainerCreate()`, `PetscContainerDestroy()`, `PetscContainerSetPointer()`, `PetscContainerGetPointer()`, `PetscObjectCompose()`, `PetscObjectQuery()`,
`PetscContainerSetCtxDestroy()`, `PetscObject`, `PetscObjectContainerQuery()`

# External Links
$(_doc_external("Sys/PetscObjectContainerCompose"))
"""
function PetscObjectContainerCompose(petsclib::PetscLibType, obj::PetscObject, name::String, pointer::Cvoid, destroy::PetscCtxDestroyFn) end

@for_petsc function PetscObjectContainerCompose(petsclib::$UnionPetscLib, obj::PetscObject, name::String, pointer::Cvoid, destroy::PetscCtxDestroyFn )

    @chk ccall(
               (:PetscObjectContainerCompose, $petsc_library),
               PetscErrorCode,
               (PetscObject, Ptr{Cchar}, Ptr{Cvoid}, Ptr{PetscCtxDestroyFn}),
               obj, name, pointer, destroy,
              )


	return nothing
end 

"""
	PetscObjectContainerQuery(petsclib::PetscLibType,obj::PetscObject, name::String, pointer::PeCtx) 
Accesses the pointer in a container composed to a `PetscObject` with `PetscObjectContainerCompose()`

Collective

Input Parameters:
- `obj`  - the `PetscObject`
- `name` - the name for the composed container

Output Parameter:
- `pointer` - the pointer to the data

Level: advanced

-seealso: `PetscContainerCreate()`, `PetscContainerDestroy()`, `PetscContainerSetPointer()`, `PetscContainerGetPointer()`, `PetscObjectCompose()`, `PetscObjectQuery()`,
`PetscContainerSetCtxDestroy()`, `PetscObject`, `PetscObjectContainerCompose()`

# External Links
$(_doc_external("Sys/PetscObjectContainerQuery"))
"""
function PetscObjectContainerQuery(petsclib::PetscLibType, obj::PetscObject, name::String, pointer::PeCtx) end

@for_petsc function PetscObjectContainerQuery(petsclib::$UnionPetscLib, obj::PetscObject, name::String, pointer::PeCtx )

    @chk ccall(
               (:PetscObjectContainerQuery, $petsc_library),
               PetscErrorCode,
               (PetscObject, Ptr{Cchar}, PeCtx),
               obj, name, pointer,
              )


	return nothing
end 

"""
	PetscObjectSetFromOptions(petsclib::PetscLibType,obj::PetscObject) 
Sets generic parameters from user options.

Collective

Input Parameter:
- `obj` - the `PetscObject`

Level: beginner

-seealso: `PetscObjectSetOptionsPrefix()`, `PetscObjectGetOptionsPrefix()`, `PetscObject`

# External Links
$(_doc_external("Sys/PetscObjectSetFromOptions"))
"""
function PetscObjectSetFromOptions(petsclib::PetscLibType, obj::PetscObject) end

@for_petsc function PetscObjectSetFromOptions(petsclib::$UnionPetscLib, obj::PetscObject )

    @chk ccall(
               (:PetscObjectSetFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscObject,),
               obj,
              )


	return nothing
end 

"""
	PetscObjectSetUp(petsclib::PetscLibType,obj::PetscObject) 
Sets up the internal data structures for later use of the object

Collective

Input Parameter:
- `obj` - the `PetscObject`

Level: advanced

-seealso: `PetscObjectDestroy()`, `PetscObject`

# External Links
$(_doc_external("Sys/PetscObjectSetUp"))
"""
function PetscObjectSetUp(petsclib::PetscLibType, obj::PetscObject) end

@for_petsc function PetscObjectSetUp(petsclib::$UnionPetscLib, obj::PetscObject )

    @chk ccall(
               (:PetscObjectSetUp, $petsc_library),
               PetscErrorCode,
               (PetscObject,),
               obj,
              )


	return nothing
end 

"""
	PetscObjectGetName(petsclib::PetscLibType,obj::PetscObject, name::String) 
Gets a string name associated with a PETSc object.

Not Collective unless `obj` has not yet been named

Input Parameters:
- `obj`  - the PETSc variable. It must be cast with a (`PetscObject`), for example,
`PetscObjectGetName`((`PetscObject`)mat,&name);
- `name` - the name associated with `obj`, do not free

Level: intermediate

-seealso: `PetscObjectSetName()`, `PetscObjectName()`, `PetscObject`, `PetscObjectGetId()`

# External Links
$(_doc_external("Sys/PetscObjectGetName"))
"""
function PetscObjectGetName(petsclib::PetscLibType, obj::PetscObject, name::String) end

@for_petsc function PetscObjectGetName(petsclib::$UnionPetscLib, obj::PetscObject, name::String )
	name_ = Ref(pointer(name))

    @chk ccall(
               (:PetscObjectGetName, $petsc_library),
               PetscErrorCode,
               (PetscObject, Ptr{Ptr{Cchar}}),
               obj, name_,
              )


	return nothing
end 

"""
	PetscObjectDestroy(petsclib::PetscLibType,obj::PetscObject) 
Destroys a `PetscObject`, regardless of the type.

Collective

Input Parameter:
- `obj` - any PETSc object, for example a `Vec`, `Mat` or `KSP`. It must be cast with a (`PetscObject`), for example,
`PetscObjectDestroy` (`PetscObject` mat);

Level: beginner

-seealso: `PetscObject`

# External Links
$(_doc_external("Sys/PetscObjectDestroy"))
"""
function PetscObjectDestroy(petsclib::PetscLibType, obj::PetscObject) end

@for_petsc function PetscObjectDestroy(petsclib::$UnionPetscLib, obj::PetscObject )

    @chk ccall(
               (:PetscObjectDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscObject},),
               obj,
              )


	return nothing
end 

"""
	PetscObjectView(petsclib::PetscLibType,obj::PetscObject, viewer::PetscViewer) 
Views a `PetscObject` regardless of the type.

Collective

Input Parameters:
- `obj`    - any PETSc object, for example a `Vec`, `Mat` or `KSP`. It must be cast with a (`PetscObject`), for example,
`PetscObjectView`((`PetscObject`)mat,`viewer`);
- `viewer` - any PETSc viewer

Level: intermediate

-seealso: `PetscObject`, `PetscObjectViewFromOptions()`, `PetscViewer`

# External Links
$(_doc_external("Sys/PetscObjectView"))
"""
function PetscObjectView(petsclib::PetscLibType, obj::PetscObject, viewer::PetscViewer) end

@for_petsc function PetscObjectView(petsclib::$UnionPetscLib, obj::PetscObject, viewer::PetscViewer )

    @chk ccall(
               (:PetscObjectView, $petsc_library),
               PetscErrorCode,
               (PetscObject, PetscViewer),
               obj, viewer,
              )


	return nothing
end 

"""
	PetscObjectViewFromOptions(petsclib::PetscLibType,obj::PetscObject, bobj::PetscObject, optionname::String) 
Processes command line options to determine if/how a `PetscObject` is to be viewed.

Collective

Input Parameters:
- `obj`        - the object
- `bobj`       - optional other object that provides prefix (if `NULL` then the prefix in `obj` is used)
- `optionname` - option string that is used to activate viewing

Options Database Key:
- `-optionname_view [viewertype]:...` - option name and values. In actual usage this would be something like `-mat_coarse_view`

Level: developer

-seealso: `PetscObject`, `PetscObjectView()`, `PetscOptionsCreateViewer()`

# External Links
$(_doc_external("Sys/PetscObjectViewFromOptions"))
"""
function PetscObjectViewFromOptions(petsclib::PetscLibType, obj::PetscObject, bobj::PetscObject, optionname::String) end

@for_petsc function PetscObjectViewFromOptions(petsclib::$UnionPetscLib, obj::PetscObject, bobj::PetscObject, optionname::String )

    @chk ccall(
               (:PetscObjectViewFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscObject, PetscObject, Ptr{Cchar}),
               obj, bobj, optionname,
              )


	return nothing
end 

"""
	same::PetscBool = PetscObjectTypeCompare(petsclib::PetscLibType,obj::PetscObject, type_name::String) 
Determines whether a PETSc object is of a particular type.

Not Collective

Input Parameters:
- `obj`       - a PETSc object, for example a `Vec`, `Mat` or `KSP`. It must be cast with a (`PetscObject`), for example,
`PetscObjectTypeCompare`((`PetscObject`)mat);
- `type_name` - string containing a type name

Output Parameter:
- `same` - `PETSC_TRUE` if the type of `obj` and `type_name` are the same or both `NULL`, else `PETSC_FALSE`

Level: intermediate

-seealso: `PetscObject`, `VecGetType()`, `KSPGetType()`, `PCGetType()`, `SNESGetType()`, `PetscObjectBaseTypeCompare()`, `PetscObjectTypeCompareAny()`, `PetscObjectBaseTypeCompareAny()`, `PetscObjectObjectTypeCompare()`

# External Links
$(_doc_external("Sys/PetscObjectTypeCompare"))
"""
function PetscObjectTypeCompare(petsclib::PetscLibType, obj::PetscObject, type_name::String) end

@for_petsc function PetscObjectTypeCompare(petsclib::$UnionPetscLib, obj::PetscObject, type_name::String )
	same_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscObjectTypeCompare, $petsc_library),
               PetscErrorCode,
               (PetscObject, Ptr{Cchar}, Ptr{PetscBool}),
               obj, type_name, same_,
              )

	same = same_[]

	return same
end 

"""
	same::PetscBool = PetscObjectObjectTypeCompare(petsclib::PetscLibType,obj1::PetscObject, obj2::PetscObject) 
Determines whether two PETSc objects are of the same type

Logically Collective

Input Parameters:
- `obj1` - any PETSc object, for example a `Vec`, `Mat` or `KSP`.
- `obj2` - another PETSc object

Output Parameter:
- `same` - `PETSC_TRUE` if they are the same or both unset, else `PETSC_FALSE`

Level: intermediate

-seealso: `PetscObjectTypeCompare()`, `VecGetType()`, `KSPGetType()`, `PCGetType()`, `SNESGetType()`, `PetscObjectBaseTypeCompare()`, `PetscObjectTypeCompareAny()`, `PetscObjectBaseTypeCompareAny()`


# External Links
$(_doc_external("Sys/PetscObjectObjectTypeCompare"))
"""
function PetscObjectObjectTypeCompare(petsclib::PetscLibType, obj1::PetscObject, obj2::PetscObject) end

@for_petsc function PetscObjectObjectTypeCompare(petsclib::$UnionPetscLib, obj1::PetscObject, obj2::PetscObject )
	same_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscObjectObjectTypeCompare, $petsc_library),
               PetscErrorCode,
               (PetscObject, PetscObject, Ptr{PetscBool}),
               obj1, obj2, same_,
              )

	same = same_[]

	return same
end 

"""
	same::PetscBool = PetscObjectBaseTypeCompare(petsclib::PetscLibType,obj::PetscObject, type_name::String) 
Determines whether a `PetscObject` is of a given base type. For example the base type of `MATSEQAIJPERM` is `MATSEQAIJ`

Not Collective

Input Parameters:
- `obj`       - the object
- `type_name` - string containing a type name

Output Parameter:
- `same` - `PETSC_TRUE` if the object is of the same base type identified by `type_name` or both `NULL`, `PETSC_FALSE` otherwise

Level: intermediate

-seealso: `PetscObject`, `PetscObjectTypeCompare()`, `PetscObjectTypeCompareAny()`, `PetscObjectBaseTypeCompareAny()`

# External Links
$(_doc_external("Sys/PetscObjectBaseTypeCompare"))
"""
function PetscObjectBaseTypeCompare(petsclib::PetscLibType, obj::PetscObject, type_name::String) end

@for_petsc function PetscObjectBaseTypeCompare(petsclib::$UnionPetscLib, obj::PetscObject, type_name::String )
	same_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscObjectBaseTypeCompare, $petsc_library),
               PetscErrorCode,
               (PetscObject, Ptr{Cchar}, Ptr{PetscBool}),
               obj, type_name, same_,
              )

	same = same_[]

	return same
end 

"""
	PetscObjectRegisterDestroy(petsclib::PetscLibType,obj::PetscObject) 
Registers a PETSc object to be destroyed when
`PetscFinalize()` is called.

Logically Collective

Input Parameter:
- `obj` - a PETSc object, for example a `Vec`, `Mat` or `KSP`. It must be cast with a (`PetscObject`), for example,
`PetscObjectRegisterDestroy`((`PetscObject`)mat);

Level: developer

-seealso: `PetscObjectRegisterDestroyAll()`

# External Links
$(_doc_external("Sys/PetscObjectRegisterDestroy"))
"""
function PetscObjectRegisterDestroy(petsclib::PetscLibType, obj::PetscObject) end

@for_petsc function PetscObjectRegisterDestroy(petsclib::$UnionPetscLib, obj::PetscObject )

    @chk ccall(
               (:PetscObjectRegisterDestroy, $petsc_library),
               PetscErrorCode,
               (PetscObject,),
               obj,
              )


	return nothing
end 

"""
	PetscObjectRegisterDestroyAll(petsclib::PetscLibType) 
Frees all the PETSc objects that have been registered
with `PetscObjectRegisterDestroy()`. Called by `PetscFinalize()`

Logically Collective on the individual `PetscObject`s that are being processed

Level: developer

-seealso: `PetscObjectRegisterDestroy()`

# External Links
$(_doc_external("Sys/PetscObjectRegisterDestroyAll"))
"""
function PetscObjectRegisterDestroyAll(petsclib::PetscLibType) end

@for_petsc function PetscObjectRegisterDestroyAll(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscObjectRegisterDestroyAll, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscObjectGetId(petsclib::PetscLibType,obj::PetscObject, id::PetscObjectId) 
get a unique object ID for the `PetscObject`

Not Collective

Input Parameter:
- `obj` - object

Output Parameter:
- `id` - integer ID

Level: developer

-seealso: `PetscObjectStateGet()`, `PetscObjectCompareId()`

# External Links
$(_doc_external("Sys/PetscObjectGetId"))
"""
function PetscObjectGetId(petsclib::PetscLibType, obj::PetscObject, id::PetscObjectId) end

@for_petsc function PetscObjectGetId(petsclib::$UnionPetscLib, obj::PetscObject, id::PetscObjectId )

    @chk ccall(
               (:PetscObjectGetId, $petsc_library),
               PetscErrorCode,
               (PetscObject, Ptr{PetscObjectId}),
               obj, id,
              )


	return nothing
end 

"""
	eq::PetscBool = PetscObjectCompareId(petsclib::PetscLibType,obj::PetscObject, id::PetscObjectId) 
compares the objects ID with a given id

Not Collective

Input Parameters:
- `obj` - object
- `id`  - integer ID

Output Parameter:
- `eq` - the ids are equal

Level: developer

-seealso: `PetscObjectStateGet()`, `PetscObjectGetId()`

# External Links
$(_doc_external("Sys/PetscObjectCompareId"))
"""
function PetscObjectCompareId(petsclib::PetscLibType, obj::PetscObject, id::PetscObjectId) end

@for_petsc function PetscObjectCompareId(petsclib::$UnionPetscLib, obj::PetscObject, id::PetscObjectId )
	eq_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscObjectCompareId, $petsc_library),
               PetscErrorCode,
               (PetscObject, PetscObjectId, Ptr{PetscBool}),
               obj, id, eq_,
              )

	eq = eq_[]

	return eq
end 

"""
	PetscObjectGetOptions(petsclib::PetscLibType,obj::PetscObject, options::PetscOptions) 
Gets the options database used by the object that has been set with `PetscObjectSetOptions()`

Collective

Input Parameter:
- `obj` - any PETSc object, for example a `Vec`, `Mat` or `KSP`.

Output Parameter:
- `options` - the options database

Level: advanced

-seealso: `PetscOptionsCreate()`, `PetscOptionsDestroy()`, `PetscObjectSetOptionsPrefix()`, `PetscObjectAppendOptionsPrefix()`, `PetscObjectPrependOptionsPrefix()`,
`PetscObjectGetOptionsPrefix()`, `PetscObjectSetOptions()`

# External Links
$(_doc_external("Sys/PetscObjectGetOptions"))
"""
function PetscObjectGetOptions(petsclib::PetscLibType, obj::PetscObject, options::PetscOptions) end

@for_petsc function PetscObjectGetOptions(petsclib::$UnionPetscLib, obj::PetscObject, options::PetscOptions )
	options_ = Ref(options.ptr)

    @chk ccall(
               (:PetscObjectGetOptions, $petsc_library),
               PetscErrorCode,
               (PetscObject, Ptr{COptions}),
               obj, options_,
              )

	options.ptr = C_NULL

	return nothing
end 

"""
	PetscObjectSetOptions(petsclib::PetscLibType,obj::PetscObject, options::PetscOptions) 
Sets the options database used by the object. Call immediately after creating the object.

Collective

Input Parameters:
- `obj`     - any PETSc object, for example a `Vec`, `Mat` or `KSP`.
- `options` - the options database, use NULL for default

Level: advanced

-seealso: `PetscOptionsCreate()`, `PetscOptionsDestroy()`, `PetscObjectSetOptionsPrefix()`, `PetscObjectAppendOptionsPrefix()`, `PetscObjectPrependOptionsPrefix()`,
`PetscObjectGetOptionsPrefix()`, `PetscObjectGetOptions()`

# External Links
$(_doc_external("Sys/PetscObjectSetOptions"))
"""
function PetscObjectSetOptions(petsclib::PetscLibType, obj::PetscObject, options::PetscOptions) end

@for_petsc function PetscObjectSetOptions(petsclib::$UnionPetscLib, obj::PetscObject, options::PetscOptions )

    @chk ccall(
               (:PetscObjectSetOptions, $petsc_library),
               PetscErrorCode,
               (PetscObject, COptions),
               obj, options,
              )


	return nothing
end 

"""
	PetscObjectSetOptionsPrefix(petsclib::PetscLibType,obj::PetscObject, prefix::String) 
Sets the prefix used for searching for all
options for the given object in the database.

Collective

Input Parameters:
- `obj`    - any PETSc object, for example a `Vec`, `Mat` or `KSP`.
- `prefix` - the prefix string to prepend to option requests of the object.

Level: advanced

-seealso: `PetscOptionsCreate()`, `PetscOptionsDestroy()`, `PetscObjectAppendOptionsPrefix()`, `PetscObjectPrependOptionsPrefix()`,
`PetscObjectGetOptionsPrefix()`, `TSSetOptionsPrefix()`, `SNESSetOptionsPrefix()`, `KSPSetOptionsPrefix()`

# External Links
$(_doc_external("Sys/PetscObjectSetOptionsPrefix"))
"""
function PetscObjectSetOptionsPrefix(petsclib::PetscLibType, obj::PetscObject, prefix::String) end

@for_petsc function PetscObjectSetOptionsPrefix(petsclib::$UnionPetscLib, obj::PetscObject, prefix::String )

    @chk ccall(
               (:PetscObjectSetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (PetscObject, Ptr{Cchar}),
               obj, prefix,
              )


	return nothing
end 

"""
	PetscObjectAppendOptionsPrefix(petsclib::PetscLibType,obj::PetscObject, prefix::String) 
Appends to the prefix used for searching for options for the given object in the database.

Input Parameters:
- `obj`    - any PETSc object, for example a `Vec`, `Mat` or `KSP`.
- `prefix` - the prefix string to prepend to option requests of the object.

Level: advanced

-seealso: `PetscOptionsCreate()`, `PetscOptionsDestroy()`, `PetscObjectSetOptionsPrefix()`, `PetscObjectPrependOptionsPrefix()`,
`PetscObjectGetOptionsPrefix()`, `TSAppendOptionsPrefix()`, `SNESAppendOptionsPrefix()`, `KSPAppendOptionsPrefix()`

# External Links
$(_doc_external("Sys/PetscObjectAppendOptionsPrefix"))
"""
function PetscObjectAppendOptionsPrefix(petsclib::PetscLibType, obj::PetscObject, prefix::String) end

@for_petsc function PetscObjectAppendOptionsPrefix(petsclib::$UnionPetscLib, obj::PetscObject, prefix::String )

    @chk ccall(
               (:PetscObjectAppendOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (PetscObject, Ptr{Cchar}),
               obj, prefix,
              )


	return nothing
end 

"""
	PetscObjectGetOptionsPrefix(petsclib::PetscLibType,obj::PetscObject, prefix::String) 
Gets the prefix of the `PetscObject` used for searching in the options database

Input Parameter:
- `obj` - any PETSc object, for example a `Vec`, `Mat` or `KSP`.

Output Parameter:
- `prefix` - pointer to the prefix string used is returned

Level: advanced

-seealso: `PetscOptionsCreate()`, `PetscOptionsDestroy()`, `PetscObjectSetOptionsPrefix()`, `PetscObjectAppendOptionsPrefix()`, `PetscObjectPrependOptionsPrefix()`,
`TSGetOptionsPrefix()`, `SNESGetOptionsPrefix()`, `KSPGetOptionsPrefix()`

# External Links
$(_doc_external("Sys/PetscObjectGetOptionsPrefix"))
"""
function PetscObjectGetOptionsPrefix(petsclib::PetscLibType, obj::PetscObject, prefix::String) end

@for_petsc function PetscObjectGetOptionsPrefix(petsclib::$UnionPetscLib, obj::PetscObject, prefix::String )
	prefix_ = Ref(pointer(prefix))

    @chk ccall(
               (:PetscObjectGetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (PetscObject, Ptr{Ptr{Cchar}}),
               obj, prefix_,
              )


	return nothing
end 

"""
	PetscObjectPrependOptionsPrefix(petsclib::PetscLibType,obj::PetscObject, prefix::String) 
Sets the prefix used for searching for options of for this object in the database.

Input Parameters:
- `obj`    - any PETSc object, for example a `Vec`, `Mat` or `KSP`.
- `prefix` - the prefix string to prepend to option requests of the object.

Level: advanced

-seealso: `PetscOptionsCreate()`, `PetscOptionsDestroy()`, `PetscObjectSetOptionsPrefix()`, `PetscObjectAppendOptionsPrefix()`,
`PetscObjectGetOptionsPrefix()`

# External Links
$(_doc_external("Sys/PetscObjectPrependOptionsPrefix"))
"""
function PetscObjectPrependOptionsPrefix(petsclib::PetscLibType, obj::PetscObject, prefix::String) end

@for_petsc function PetscObjectPrependOptionsPrefix(petsclib::$UnionPetscLib, obj::PetscObject, prefix::String )

    @chk ccall(
               (:PetscObjectPrependOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (PetscObject, Ptr{Cchar}),
               obj, prefix,
              )


	return nothing
end 

"""
	comm::MPI_Comm = PetscObjectGetComm(petsclib::PetscLibType,obj) 
Gets the MPI communicator for any `PetscObject` regardless of the type.

Not Collective

Input Parameter:
- `obj` - any PETSc object, for example a `Vec`, `Mat` or `KSP`. It must be cast with a (`PetscObject`), for example,
`PetscObjectGetComm`((`PetscObject`)mat,&comm);

Output Parameter:
- `comm` - the MPI communicator

Level: advanced

-seealso: `PetscObject`, `PetscObjectComm()`

# External Links
$(_doc_external("Sys/PetscObjectGetComm"))
"""
function PetscObjectGetComm(petsclib::PetscLibType, obj) end

@for_petsc function PetscObjectGetComm(petsclib::$UnionPetscLib, obj)
    comm = Ref(MPI_Comm())
    @chk ccall(
               (:PetscObjectGetComm, $petsc_library),
               PetscErrorCode,
               (PetscObject, Ptr{MPI_Comm}),
               obj, comm
              )

	return comm[]
end 

"""
	tab::PetscInt = PetscObjectGetTabLevel(petsclib::PetscLibType,obj::PetscObject) 
Gets the number of tabs that `PETSCVIEWERASCII` output for that object uses

Not Collective

Input Parameter:
- `obj` - any PETSc object, for example a `Vec`, `Mat` or `KSP`. It must be cast with a (`PetscObject`), for example,
`PetscObjectGetTabLevel`((`PetscObject`)mat,&tab);

Output Parameter:
- `tab` - the number of tabs

Level: developer

-seealso: `PetscObjectIncrementTabLevel()`, `PetscObjectSetTabLevel()`, `PETSCVIEWERASCII`, `PetscObject`

# External Links
$(_doc_external("Sys/PetscObjectGetTabLevel"))
"""
function PetscObjectGetTabLevel(petsclib::PetscLibType, obj::PetscObject) end

@for_petsc function PetscObjectGetTabLevel(petsclib::$UnionPetscLib, obj::PetscObject )
	tab_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscObjectGetTabLevel, $petsc_library),
               PetscErrorCode,
               (PetscObject, Ptr{$PetscInt}),
               obj, tab_,
              )

	tab = tab_[]

	return tab
end 

"""
	PetscObjectSetTabLevel(petsclib::PetscLibType,obj::PetscObject, tab::PetscInt) 
Sets the number of tabs that `PETSCVIEWERASCII` output for that object uses

Not Collective

Input Parameters:
- `obj` - any PETSc object, for example a `Vec`, `Mat` or `KSP`. It must be cast with a (`PetscObject`), for example,
`PetscObjectSetTabLevel`((`PetscObject`)mat,tab;
- `tab` - the number of tabs

Level: developer

-seealso: `PetscObjectIncrementTabLevel()`, `PetscObjectGetTabLevel()`

# External Links
$(_doc_external("Sys/PetscObjectSetTabLevel"))
"""
function PetscObjectSetTabLevel(petsclib::PetscLibType, obj::PetscObject, tab::PetscInt) end

@for_petsc function PetscObjectSetTabLevel(petsclib::$UnionPetscLib, obj::PetscObject, tab::$PetscInt )

    @chk ccall(
               (:PetscObjectSetTabLevel, $petsc_library),
               PetscErrorCode,
               (PetscObject, $PetscInt),
               obj, tab,
              )


	return nothing
end 

"""
	PetscObjectIncrementTabLevel(petsclib::PetscLibType,obj::PetscObject, oldobj::PetscObject, tab::PetscInt) 
Increments the number of tabs that `PETSCVIEWERASCII` output for that object use based on
the tablevel of another object. This should be called immediately after the object is created.

Not Collective

Input Parameters:
- `obj`    - any PETSc object where we are changing the tab
- `oldobj` - the object providing the tab, optional pass `NULL` to use 0 as the previous tablevel for `obj`
- `tab`    - the increment that is added to the old objects tab

Level: developer

-seealso: `PETSCVIEWERASCII`, `PetscObjectSetTabLevel()`, `PetscObjectGetTabLevel()`

# External Links
$(_doc_external("Sys/PetscObjectIncrementTabLevel"))
"""
function PetscObjectIncrementTabLevel(petsclib::PetscLibType, obj::PetscObject, oldobj::PetscObject, tab::PetscInt) end

@for_petsc function PetscObjectIncrementTabLevel(petsclib::$UnionPetscLib, obj::PetscObject, oldobj::PetscObject, tab::$PetscInt )

    @chk ccall(
               (:PetscObjectIncrementTabLevel, $petsc_library),
               PetscErrorCode,
               (PetscObject, PetscObject, $PetscInt),
               obj, oldobj, tab,
              )


	return nothing
end 

"""
	PetscObjectGetNewTag(petsclib::PetscLibType,obj::PetscObject, tag::PetscMPIInt) 
Gets a unique new tag from a PETSc object. All
processors that share the object MUST call this routine EXACTLY the same
number of times.  This tag should only be used with the current objects
communicator; do NOT use it with any other MPI communicator.

Collective

Input Parameter:
- `obj` - the PETSc object; this must be cast with a (`PetscObject`), for example,
`PetscObjectGetNewTag`((`PetscObject`)mat,&tag);

Output Parameter:
- `tag` - the new tag

Level: developer

-seealso: `PetscCommGetNewTag()`

# External Links
$(_doc_external("Sys/PetscObjectGetNewTag"))
"""
function PetscObjectGetNewTag(petsclib::PetscLibType, obj::PetscObject, tag::PetscMPIInt) end

@for_petsc function PetscObjectGetNewTag(petsclib::$UnionPetscLib, obj::PetscObject, tag::PetscMPIInt )

    @chk ccall(
               (:PetscObjectGetNewTag, $petsc_library),
               PetscErrorCode,
               (PetscObject, Ptr{PetscMPIInt}),
               obj, tag,
              )


	return nothing
end 

"""
	count::PetscInt,numbering::PetscInt = PetscObjectsListGetGlobalNumbering(petsclib::PetscLibType,comm::MPI_Comm, len::PetscInt, objlist::Vector{PetscObject}) 
computes a global numbering
of `PetscObject`s living on subcommunicators of a given communicator.

Collective.

Input Parameters:
- `comm`    - the `MPI_Comm`
- `len`     - local length of `objlist`
- `objlist` - a list of PETSc objects living on subcomms of comm and containing this comm rank
(subcomm ordering is assumed to be deadlock-free)

Output Parameters:
- `count`     - global number of distinct subcommunicators on objlist (may be > `len`)
- `numbering` - global numbers of objlist entries (allocated by user)

Level: developer

-seealso: `PetscCommDuplicate()`, `PetscObjectDestroy()`

# External Links
$(_doc_external("Sys/PetscObjectsListGetGlobalNumbering"))
"""
function PetscObjectsListGetGlobalNumbering(petsclib::PetscLibType, comm::MPI_Comm, len::PetscInt, objlist::Vector{PetscObject}) end

@for_petsc function PetscObjectsListGetGlobalNumbering(petsclib::$UnionPetscLib, comm::MPI_Comm, len::$PetscInt, objlist::Vector{PetscObject} )
	count_ = Ref{$PetscInt}()
	numbering_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscObjectsListGetGlobalNumbering, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, Ptr{PetscObject}, Ptr{$PetscInt}, Ptr{$PetscInt}),
               comm, len, objlist, count_, numbering_,
              )

	count = count_[]
	numbering = numbering_[]

	return count,numbering
end 

"""
	PetscObjectGetClassId(petsclib::PetscLibType,obj::PetscObject, classid::PetscClassId) 
Gets the classid for any `PetscObject`

Not Collective

Input Parameter:
- `obj` - any PETSc object, for example a `Vec`, `Mat` or `KSP`. It must be cast with a (`PetscObject`), for example,
`PetscObjectGetClassId`((`PetscObject`)mat,&classid);

Output Parameter:
- `classid` - the classid

Level: developer

-seealso: `PetscObject`, `PetscClassId`, `PetscObjectGetClassName()`, `PetscObjectGetType()`

# External Links
$(_doc_external("Sys/PetscObjectGetClassId"))
"""
function PetscObjectGetClassId(petsclib::PetscLibType, obj::PetscObject, classid::PetscClassId) end

@for_petsc function PetscObjectGetClassId(petsclib::$UnionPetscLib, obj::PetscObject, classid::PetscClassId )

    @chk ccall(
               (:PetscObjectGetClassId, $petsc_library),
               PetscErrorCode,
               (PetscObject, Ptr{PetscClassId}),
               obj, classid,
              )


	return nothing
end 

"""
	PetscObjectGetClassName(petsclib::PetscLibType,obj::PetscObject, classname::String) 
Gets the class name for any `PetscObject`

Not Collective

Input Parameter:
- `obj` - any PETSc object, for example a `Vec`, `Mat` or `KSP`. It must be cast with a (`PetscObject`), for example,
`PetscObjectGetClassName`((`PetscObject`)mat,&classname);

Output Parameter:
- `classname` - the class name, for example "Vec"

Level: developer

-seealso: `PetscObject`, `PetscClassId`, `PetscObjectGetType()`, `PetscObjectGetClassId()`

# External Links
$(_doc_external("Sys/PetscObjectGetClassName"))
"""
function PetscObjectGetClassName(petsclib::PetscLibType, obj::PetscObject, classname::String) end

@for_petsc function PetscObjectGetClassName(petsclib::$UnionPetscLib, obj::PetscObject, classname::String )
	classname_ = Ref(pointer(classname))

    @chk ccall(
               (:PetscObjectGetClassName, $petsc_library),
               PetscErrorCode,
               (PetscObject, Ptr{Ptr{Cchar}}),
               obj, classname_,
              )


	return nothing
end 

"""
	PetscObjectSetName(petsclib::PetscLibType,obj::PetscObject, name::String) 
Sets a string name for a PETSc object.

Not Collective

Input Parameters:
- `obj`  - the PETSc object. It must be cast with a (`PetscObject`), for example,
`PetscObjectSetName`((`PetscObject`)mat,name);
- `name` - the name to give obj

Level: advanced

-seealso: `PetscObjectGetName()`, `PetscObjectName()`

# External Links
$(_doc_external("Sys/PetscObjectSetName"))
"""
function PetscObjectSetName(petsclib::PetscLibType, obj::PetscObject, name::String) end

@for_petsc function PetscObjectSetName(petsclib::$UnionPetscLib, obj::PetscObject, name::String )

    @chk ccall(
               (:PetscObjectSetName, $petsc_library),
               PetscErrorCode,
               (PetscObject, Ptr{Cchar}),
               obj, name,
              )


	return nothing
end 

"""
	PetscObjectPrintClassNamePrefixType(petsclib::PetscLibType,obj::PetscObject, viewer::PetscViewer) 
used in the `XXXView()` methods to display information about the class, name, prefix and type of an object

Input Parameters:
- `obj`    - the PETSc object
- `viewer` - `PETSCVIEWERASCII` viewer where the information is printed, function does nothing if the viewer is not `PETSCVIEWERASCII` type

Level: developer

-seealso: `PetscObjectSetName()`, `PetscObjectName()`

# External Links
$(_doc_external("Sys/PetscObjectPrintClassNamePrefixType"))
"""
function PetscObjectPrintClassNamePrefixType(petsclib::PetscLibType, obj::PetscObject, viewer::PetscViewer) end

@for_petsc function PetscObjectPrintClassNamePrefixType(petsclib::$UnionPetscLib, obj::PetscObject, viewer::PetscViewer )

    @chk ccall(
               (:PetscObjectPrintClassNamePrefixType, $petsc_library),
               PetscErrorCode,
               (PetscObject, PetscViewer),
               obj, viewer,
              )


	return nothing
end 

"""
	PetscObjectName(petsclib::PetscLibType,obj::PetscObject) 
Gives `obj` a name if it does not have one

Collective

Input Parameter:
- `obj` - the PETSc object. It must be cast with a (`PetscObject`), for example, `PetscObjectName`((`PetscObject`)mat,name);

Level: developer

-seealso: `PetscObjectGetName()`, `PetscObjectSetName()`

# External Links
$(_doc_external("Sys/PetscObjectName"))
"""
function PetscObjectName(petsclib::PetscLibType, obj::PetscObject) end

@for_petsc function PetscObjectName(petsclib::$UnionPetscLib, obj::PetscObject )

    @chk ccall(
               (:PetscObjectName, $petsc_library),
               PetscErrorCode,
               (PetscObject,),
               obj,
              )


	return nothing
end 

"""
	PetscObjectChangeTypeName(petsclib::PetscLibType,obj::PetscObject, type_name::String) 

# External Links
$(_doc_external("Sys/PetscObjectChangeTypeName"))
"""
function PetscObjectChangeTypeName(petsclib::PetscLibType, obj::PetscObject, type_name::String) end

@for_petsc function PetscObjectChangeTypeName(petsclib::$UnionPetscLib, obj::PetscObject, type_name::String )

    @chk ccall(
               (:PetscObjectChangeTypeName, $petsc_library),
               PetscErrorCode,
               (PetscObject, Ptr{Cchar}),
               obj, type_name,
              )


	return nothing
end 

