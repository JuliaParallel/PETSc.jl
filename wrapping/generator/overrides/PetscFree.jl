# override for PetscFree; C macro: PetscFree(a), a call through the function pointer PetscTrFree
"""
	PetscFree(petsclib::PetscLibType, ptr::Ptr)
Frees memory that PETSc allocated with `PetscMalloc()` and handed to the caller, such as the arrays
some `Create` and `Get` functions return with the note "the caller must free it with `PetscFree()`".
Does nothing on `C_NULL`.

`PetscFree()` is a C macro, so this calls the allocator's free routine `PetscTrFree` directly, which
also honours an allocator installed with `PetscMallocSet()`.

Not Collective

Input Parameter:
- `ptr` - memory allocated by PETSc

Level: beginner

See also: `PetscMalloc()`, `PetscMallocSet()`

# External Links
$(_doc_external("Sys/PetscFree"))
"""
function PetscFree(petsclib::PetscLibType, ptr::Ptr) end

@for_petsc function PetscFree(petsclib::$UnionPetscLib, ptr::Ptr)
	ptr == C_NULL && return nothing
	# PetscTrFree is a global function pointer, read on every call because PetscMallocSet can replace it
	free_ = unsafe_load(cglobal((:PetscTrFree, $petsc_library), Ptr{Cvoid}))
	@chk ccall(
		free_,
		PetscErrorCode,
		(Ptr{Cvoid}, Cint, Cstring, Cstring),
		ptr, Cint(0), "PetscFree", "PETSc.jl",
	)
	return nothing
end
