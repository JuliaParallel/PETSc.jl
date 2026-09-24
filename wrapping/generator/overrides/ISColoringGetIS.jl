# override for ISColoringGetIS; C signature: ISColoringGetIS(ISColoring iscoloring, PetscCopyMode mode, PetscInt* nn, IS* isis[])
"""
	nn::PetscInt,isis::Vector{IS} = ISColoringGetIS(petsclib::PetscLibType, iscoloring::ISColoring, mode::PetscCopyMode)
Extracts index sets from the coloring context. Each is contains the nodes of one color

Collective

Input Parameters:
- `iscoloring` - the coloring context
- `mode`       - `PETSC_OWN_POINTER` hands the index sets to the caller; any other value leaves them with the coloring

Output Parameters:
- `nn`   - number of index sets in the coloring context
- `isis` - array of index sets

With `PETSC_OWN_POINTER` each `IS` in `isis` is owned by the caller and must be released with `ISDestroy()` (or
`PETSc.destroy!`); the C array that held them is freed here. With any other mode the index sets are borrowed: the
coloring destroys them in `ISColoringDestroy()`, and destroying one here does nothing.

Level: advanced

See also: `ISColoring`, `IS`, `ISColoringRestoreIS()`, `ISColoringView()`, `ISColoringGetColoring()`, `ISColoringGetColors()`

# External Links
$(_doc_external("IS/ISColoringGetIS"))
"""
function ISColoringGetIS(petsclib::PetscLibType, iscoloring::ISColoring, mode::PetscCopyMode) end

@for_petsc function ISColoringGetIS(petsclib::$UnionPetscLib, iscoloring::ISColoring, mode::PetscCopyMode)
	nn_ = Ref{$PetscInt}()
	isis_ = Ref{Ptr{CIS}}(C_NULL)

	@chk ccall(
		(:ISColoringGetIS, $petsc_library),
		PetscErrorCode,
		(ISColoring, PetscCopyMode, Ptr{$PetscInt}, Ptr{Ptr{CIS}}),
		iscoloring, mode, nn_, isis_,
	)

	nn = nn_[]
	own = mode == PETSC_OWN_POINTER
	isis = IS{$PetscLib}[]
	if isis_[] != C_NULL
		for i in 1:nn
			push!(isis, IS(unsafe_load(isis_[], i), petsclib; own))
		end
		# the coloring keeps the array unless the caller takes it
		own && PetscFree(petsclib, isis_[])
	end
	return nn, isis
end
