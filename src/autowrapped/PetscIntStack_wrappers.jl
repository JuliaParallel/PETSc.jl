
"""
	PetscIntStackDestroy(petsclib::PetscLibType,stack::PetscIntStack) 
This function destroys a stack.

Not Collective, No Fortran Support

Input Parameter:
- `stack` - The stack

Level: developer

-seealso: `PetscIntStackCreate()`, `PetscIntStackEmpty()`, `PetscIntStackPush()`, `PetscIntStackPop()`, `PetscIntStackTop()`

# External Links
$(_doc_external("Sys/PetscIntStackDestroy"))
"""
function PetscIntStackDestroy(petsclib::PetscLibType, stack::PetscIntStack) end

@for_petsc function PetscIntStackDestroy(petsclib::$UnionPetscLib, stack::$PetscIntStack )

    @chk ccall(
               (:PetscIntStackDestroy, $petsc_library),
               PetscErrorCode,
               ($PetscIntStack,),
               stack,
              )


	return nothing
end 

"""
	empty::PetscBool = PetscIntStackEmpty(petsclib::PetscLibType,stack::PetscIntStack) 
This function determines whether any items have been pushed.

Not Collective, No Fortran Support

Input Parameter:
- `stack` - The stack

Output Parameter:
- `empty` - `PETSC_TRUE` if the stack is empty

Level: developer

-seealso: `PetscIntStackCreate()`, `PetscIntStackDestroy()`, `PetscIntStackPush()`, `PetscIntStackPop()`, `PetscIntStackTop()`

# External Links
$(_doc_external("Sys/PetscIntStackEmpty"))
"""
function PetscIntStackEmpty(petsclib::PetscLibType, stack::PetscIntStack) end

@for_petsc function PetscIntStackEmpty(petsclib::$UnionPetscLib, stack::$PetscIntStack )
	empty_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscIntStackEmpty, $petsc_library),
               PetscErrorCode,
               ($PetscIntStack, Ptr{PetscBool}),
               stack, empty_,
              )

	empty = empty_[]

	return empty
end 

"""
	PetscIntStackTop(petsclib::PetscLibType,stack::PetscIntStack, top::Cint) 
This function returns the top of the stack.

Not Collective, No Fortran Support

Input Parameter:
- `stack` - The stack

Output Parameter:
- `top` - The integer on top of the stack

Level: developer

-seealso: `PetscIntStackCreate()`, `PetscIntStackDestroy()`, `PetscIntStackEmpty()`, `PetscIntStackPush()`, `PetscIntStackPop()`

# External Links
$(_doc_external("Sys/PetscIntStackTop"))
"""
function PetscIntStackTop(petsclib::PetscLibType, stack::PetscIntStack, top::Cint) end

@for_petsc function PetscIntStackTop(petsclib::$UnionPetscLib, stack::$PetscIntStack, top::Cint )

    @chk ccall(
               (:PetscIntStackTop, $petsc_library),
               PetscErrorCode,
               ($PetscIntStack, Ptr{Cint}),
               stack, top,
              )


	return nothing
end 

"""
	PetscIntStackPush(petsclib::PetscLibType,stack::PetscIntStack, item::Cint) 
This function pushes an integer on the stack.

Not Collective, No Fortran Support

Input Parameters:
- `stack` - The stack
- `item`  - The integer to push

Level: developer

-seealso: `PetscIntStackCreate()`, `PetscIntStackDestroy()`, `PetscIntStackEmpty()`, `PetscIntStackPop()`, `PetscIntStackTop()`

# External Links
$(_doc_external("Sys/PetscIntStackPush"))
"""
function PetscIntStackPush(petsclib::PetscLibType, stack::PetscIntStack, item::Cint) end

@for_petsc function PetscIntStackPush(petsclib::$UnionPetscLib, stack::$PetscIntStack, item::Cint )

    @chk ccall(
               (:PetscIntStackPush, $petsc_library),
               PetscErrorCode,
               ($PetscIntStack, Cint),
               stack, item,
              )


	return nothing
end 

"""
	PetscIntStackPop(petsclib::PetscLibType,stack::PetscIntStack, item::Cint) 
This function pops an integer from the stack.

Not Collective, No Fortran Support

Input Parameter:
- `stack` - The stack

Output Parameter:
- `item` - The integer popped

Level: developer

-seealso: `PetscIntStackCreate()`, `PetscIntStackDestroy()`, `PetscIntStackEmpty()`, `PetscIntStackPush()`, `PetscIntStackTop()`

# External Links
$(_doc_external("Sys/PetscIntStackPop"))
"""
function PetscIntStackPop(petsclib::PetscLibType, stack::PetscIntStack, item::Cint) end

@for_petsc function PetscIntStackPop(petsclib::$UnionPetscLib, stack::$PetscIntStack, item::Cint )

    @chk ccall(
               (:PetscIntStackPop, $petsc_library),
               PetscErrorCode,
               ($PetscIntStack, Ptr{Cint}),
               stack, item,
              )


	return nothing
end 

"""
	stack::PetscIntStack = PetscIntStackCreate(petsclib::PetscLibType) 
This function creates a stack.

Not Collective, No Fortran Support

Output Parameter:
- `stack` - The stack

Level: developer

-seealso: `PetscIntStackDestroy()`, `PetscIntStackEmpty()`, `PetscIntStackPush()`, `PetscIntStackPop()`, `PetscIntStackTop()`

# External Links
$(_doc_external("Sys/PetscIntStackCreate"))
"""
function PetscIntStackCreate(petsclib::PetscLibType) end

@for_petsc function PetscIntStackCreate(petsclib::$UnionPetscLib)
	stack_ = Ref{$PetscIntStack}()

    @chk ccall(
               (:PetscIntStackCreate, $petsc_library),
               PetscErrorCode,
               (Ptr{$PetscIntStack},),
               stack_,
              )

	stack = stack_[]

	return stack
end 

