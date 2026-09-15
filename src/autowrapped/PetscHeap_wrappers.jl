"""
	PetscHeapAdd(petsclib::PetscLibType,h::PetscHeap, id::PetscInt, val::PetscInt) 
Insert an item into a `PetscHeap`.

Not Collective

Input Parameters:
- `h`   - the `PetscHeap`
- `id`  - the item identifier
- `val` - the value used for heap ordering

Level: developer

-seealso: `PetscHeap`, `PetscHeapCreate()`, `PetscHeapPop()`, `PetscHeapPeek()`, `PetscHeapStash()`, `PetscHeapUnstash()`, `PetscHeapDestroy()`

# External Links
$(_doc_external("Mat/PetscHeapAdd"))
"""
function PetscHeapAdd(petsclib::PetscLibType, h::PetscHeap, id::Integer, val::Integer)
    error("PetscHeapAdd: no generated method for these argument types")
end

@for_petsc function PetscHeapAdd(petsclib::$UnionPetscLib, h::PetscHeap, id::$PetscInt, val::$PetscInt )

    @chk ccall(
               (:PetscHeapAdd, $petsc_library),
               PetscErrorCode,
               (PetscHeap, $PetscInt, $PetscInt),
               h, id, val,
              )


	return nothing
end 

"""
	heap::PetscHeap = PetscHeapCreate(petsclib::PetscLibType,maxsize::PetscInt) 
Creates a `PetscHeap` object, a simple min

Not Collective

Input Parameter:
- `maxsize` - the maximum number of items the heap can hold at once

Output Parameter:
- `heap` - the newly created `PetscHeap` object

Level: developer

-seealso: `PetscHeap`, `PetscHeapAdd()`, `PetscHeapPop()`, `PetscHeapPeek()`, `PetscHeapStash()`, `PetscHeapUnstash()`, `PetscHeapView()`, `PetscHeapDestroy()`

# External Links
$(_doc_external("Mat/PetscHeapCreate"))
"""
function PetscHeapCreate(petsclib::PetscLibType, maxsize::Integer)
    error("PetscHeapCreate: no generated method for these argument types")
end

@for_petsc function PetscHeapCreate(petsclib::$UnionPetscLib, maxsize::$PetscInt )
	heap_ = Ref{PetscHeap}()

    @chk ccall(
               (:PetscHeapCreate, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{PetscHeap}),
               maxsize, heap_,
              )

	heap = heap_[]

	return heap
end 

"""
	PetscHeapDestroy(petsclib::PetscLibType,heap::Union{PetscHeap, Ref{PetscHeap}}) 
Destroys a `PetscHeap` created with `PetscHeapCreate()`.

Not Collective

Input Parameter:
- `heap` - the `PetscHeap` to destroy; set to `NULL` on return

Level: developer

-seealso: `PetscHeap`, `PetscHeapCreate()`

# External Links
$(_doc_external("Mat/PetscHeapDestroy"))
"""
function PetscHeapDestroy(petsclib::PetscLibType, heap::Union{PetscHeap, Ref{PetscHeap}})
    error("PetscHeapDestroy: no generated method for these argument types")
end

@for_petsc function PetscHeapDestroy(petsclib::$UnionPetscLib, heap::Union{PetscHeap, Ref{PetscHeap}} )
	heap_ = heap isa Base.RefValue ? heap : Ref{PetscHeap}(heap)

    @chk ccall(
               (:PetscHeapDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscHeap},),
               heap_,
              )


	return nothing
end 

"""
	id::PetscInt,val::PetscInt = PetscHeapPeek(petsclib::PetscLibType,h::PetscHeap) 
Return the minimum item of a `PetscHeap` without removing it.

Not Collective

Input Parameter:
- `h` - the `PetscHeap`

Output Parameters:
- `id`  - identifier of the minimum item, or `-1` if the heap is empty
- `val` - value of the minimum item, or `PETSC_INT_MIN` if the heap is empty

Level: developer

-seealso: `PetscHeap`, `PetscHeapCreate()`, `PetscHeapAdd()`, `PetscHeapPop()`, `PetscHeapDestroy()`

# External Links
$(_doc_external("Mat/PetscHeapPeek"))
"""
function PetscHeapPeek(petsclib::PetscLibType, h::PetscHeap)
    error("PetscHeapPeek: no generated method for these argument types")
end

@for_petsc function PetscHeapPeek(petsclib::$UnionPetscLib, h::PetscHeap )
	id_ = Ref{$PetscInt}()
	val_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscHeapPeek, $petsc_library),
               PetscErrorCode,
               (PetscHeap, Ptr{$PetscInt}, Ptr{$PetscInt}),
               h, id_, val_,
              )

	id = id_[]
	val = val_[]

	return id,val
end 

"""
	id::PetscInt,val::PetscInt = PetscHeapPop(petsclib::PetscLibType,h::PetscHeap) 
Remove and return the minimum item from a `PetscHeap`.

Not Collective

Input Parameter:
- `h` - the `PetscHeap`

Output Parameters:
- `id`  - identifier of the popped item, or `-1` if the heap is empty
- `val` - value of the popped item, or `PETSC_INT_MIN` if the heap is empty

Level: developer

-seealso: `PetscHeap`, `PetscHeapCreate()`, `PetscHeapAdd()`, `PetscHeapPeek()`, `PetscHeapStash()`, `PetscHeapUnstash()`, `PetscHeapDestroy()`

# External Links
$(_doc_external("Mat/PetscHeapPop"))
"""
function PetscHeapPop(petsclib::PetscLibType, h::PetscHeap)
    error("PetscHeapPop: no generated method for these argument types")
end

@for_petsc function PetscHeapPop(petsclib::$UnionPetscLib, h::PetscHeap )
	id_ = Ref{$PetscInt}()
	val_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscHeapPop, $petsc_library),
               PetscErrorCode,
               (PetscHeap, Ptr{$PetscInt}, Ptr{$PetscInt}),
               h, id_, val_,
              )

	id = id_[]
	val = val_[]

	return id,val
end 

"""
	PetscHeapStash(petsclib::PetscLibType,h::PetscHeap, id::PetscInt, val::PetscInt) 
Set aside an item in a `PetscHeap` for later insertion via `PetscHeapUnstash()`.

Not Collective

Input Parameters:
- `h`   - the `PetscHeap`
- `id`  - the item identifier
- `val` - the value used for heap ordering

Level: developer

-seealso: `PetscHeap`, `PetscHeapCreate()`, `PetscHeapAdd()`, `PetscHeapUnstash()`, `PetscHeapDestroy()`

# External Links
$(_doc_external("Mat/PetscHeapStash"))
"""
function PetscHeapStash(petsclib::PetscLibType, h::PetscHeap, id::Integer, val::Integer)
    error("PetscHeapStash: no generated method for these argument types")
end

@for_petsc function PetscHeapStash(petsclib::$UnionPetscLib, h::PetscHeap, id::$PetscInt, val::$PetscInt )

    @chk ccall(
               (:PetscHeapStash, $petsc_library),
               PetscErrorCode,
               (PetscHeap, $PetscInt, $PetscInt),
               h, id, val,
              )


	return nothing
end 

"""
	PetscHeapUnstash(petsclib::PetscLibType,h::PetscHeap) 
Reinsert all items previously stashed with `PetscHeapStash()` into the heap.

Not Collective

Input Parameter:
- `h` - the `PetscHeap`

Level: developer

-seealso: `PetscHeap`, `PetscHeapCreate()`, `PetscHeapAdd()`, `PetscHeapStash()`, `PetscHeapDestroy()`

# External Links
$(_doc_external("Mat/PetscHeapUnstash"))
"""
function PetscHeapUnstash(petsclib::PetscLibType, h::PetscHeap)
    error("PetscHeapUnstash: no generated method for these argument types")
end

@for_petsc function PetscHeapUnstash(petsclib::$UnionPetscLib, h::PetscHeap )

    @chk ccall(
               (:PetscHeapUnstash, $petsc_library),
               PetscErrorCode,
               (PetscHeap,),
               h,
              )


	return nothing
end 

"""
	PetscHeapView(petsclib::PetscLibType,h::PetscHeap, viewer::PetscViewer) 
View the contents of a `PetscHeap`, including any stashed items.

Not Collective

Input Parameters:
- `h`      - the `PetscHeap`
- `viewer` - a `PetscViewer`, or `NULL` to use `PETSC_VIEWER_STDOUT_SELF`

Level: developer

-seealso: `PetscHeap`, `PetscHeapCreate()`, `PetscHeapAdd()`, `PetscHeapPop()`

# External Links
$(_doc_external("Mat/PetscHeapView"))
"""
function PetscHeapView(petsclib::PetscLibType, h::PetscHeap, viewer::PetscViewer)
    error("PetscHeapView: no generated method for these argument types")
end

@for_petsc function PetscHeapView(petsclib::$UnionPetscLib, h::PetscHeap, viewer::PetscViewer )

    @chk ccall(
               (:PetscHeapView, $petsc_library),
               PetscErrorCode,
               (PetscHeap, PetscViewer),
               h, viewer,
              )


	return nothing
end 

