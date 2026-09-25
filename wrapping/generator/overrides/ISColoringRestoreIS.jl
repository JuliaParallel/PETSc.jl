# override for ISColoringRestoreIS; C signature: ISColoringRestoreIS(ISColoring iscoloring, PetscCopyMode mode, IS* is[])
"""
	ISColoringRestoreIS(petsclib::PetscLibType, iscoloring::ISColoring, mode::PetscCopyMode, is::Union{Ptr, AbstractVector{<:AbstractIS}})
Restores the index sets extracted from the coloring context with `ISColoringGetIS()` using `PETSC_USE_POINTER`

Collective

Input Parameters:
- `iscoloring` - the coloring context
- `mode`       - who retains ownership of the is
- `is`         - the index sets `ISColoringGetIS()` returned, or a pointer to a C array of them

Level: advanced

See also: `ISColoring()`, `IS`, `ISColoringGetIS()`, `ISColoringView()`, `PetscCopyMode`

# External Links
$(_doc_external("IS/ISColoringRestoreIS"))
"""
function ISColoringRestoreIS(petsclib::PetscLibType, iscoloring::ISColoring, mode::PetscCopyMode, is::Union{Ptr, AbstractVector{<:AbstractIS}}) end

@for_petsc function ISColoringRestoreIS(petsclib::$UnionPetscLib, iscoloring::ISColoring, mode::PetscCopyMode, is::Union{Ptr, AbstractVector{<:AbstractIS}})
	# C takes an IS* array: build one from the handles when given the Julia wrappers
	handles = is isa Ptr ? CIS[] : CIS[x.ptr for x in is]
	GC.@preserve handles begin
		is_ = Ref{Ptr{CIS}}(is isa Ptr ? is : pointer(handles))
		@chk ccall(
			(:ISColoringRestoreIS, $petsc_library),
			PetscErrorCode,
			(ISColoring, PetscCopyMode, Ptr{Ptr{CIS}}),
			iscoloring, mode, is_,
		)
	end
	return nothing
end
