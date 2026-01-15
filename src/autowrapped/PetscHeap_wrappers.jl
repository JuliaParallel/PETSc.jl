# autodefined type arguments for class ------
mutable struct _n_PetscHeap end
const PetscHeap = Ptr{_n_PetscHeap}

# -------------------------------------------------------
"""
	heap::PetscHeap = PetscHeapCreate(petsclib::PetscLibType,maxsize::PetscInt) 

# External Links
$(_doc_external("Mat/PetscHeapCreate"))
"""
function PetscHeapCreate(petsclib::PetscLibType, maxsize::PetscInt) end

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
	PetscHeapAdd(petsclib::PetscLibType,h::PetscHeap, id::PetscInt, val::PetscInt) 

# External Links
$(_doc_external("Mat/PetscHeapAdd"))
"""
function PetscHeapAdd(petsclib::PetscLibType, h::PetscHeap, id::PetscInt, val::PetscInt) end

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
	id::PetscInt,val::PetscInt = PetscHeapPop(petsclib::PetscLibType,h::PetscHeap) 

# External Links
$(_doc_external("Mat/PetscHeapPop"))
"""
function PetscHeapPop(petsclib::PetscLibType, h::PetscHeap) end

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
	id::PetscInt,val::PetscInt = PetscHeapPeek(petsclib::PetscLibType,h::PetscHeap) 

# External Links
$(_doc_external("Mat/PetscHeapPeek"))
"""
function PetscHeapPeek(petsclib::PetscLibType, h::PetscHeap) end

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
	PetscHeapStash(petsclib::PetscLibType,h::PetscHeap, id::PetscInt, val::PetscInt) 

# External Links
$(_doc_external("Mat/PetscHeapStash"))
"""
function PetscHeapStash(petsclib::PetscLibType, h::PetscHeap, id::PetscInt, val::PetscInt) end

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

# External Links
$(_doc_external("Mat/PetscHeapUnstash"))
"""
function PetscHeapUnstash(petsclib::PetscLibType, h::PetscHeap) end

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
	PetscHeapDestroy(petsclib::PetscLibType,heap::PetscHeap) 

# External Links
$(_doc_external("Mat/PetscHeapDestroy"))
"""
function PetscHeapDestroy(petsclib::PetscLibType, heap::PetscHeap) end

@for_petsc function PetscHeapDestroy(petsclib::$UnionPetscLib, heap::PetscHeap )

    @chk ccall(
               (:PetscHeapDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscHeap},),
               heap,
              )


	return nothing
end 

"""
	PetscHeapView(petsclib::PetscLibType,h::PetscHeap, viewer::PetscViewer) 

# External Links
$(_doc_external("Mat/PetscHeapView"))
"""
function PetscHeapView(petsclib::PetscLibType, h::PetscHeap, viewer::PetscViewer) end

@for_petsc function PetscHeapView(petsclib::$UnionPetscLib, h::PetscHeap, viewer::PetscViewer )

    @chk ccall(
               (:PetscHeapView, $petsc_library),
               PetscErrorCode,
               (PetscHeap, PetscViewer),
               h, viewer,
              )


	return nothing
end 

