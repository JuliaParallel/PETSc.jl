# override for PetscSFGetGraphLayout; C signature: PetscSFGetGraphLayout(PetscSF sf, PetscLayout* layout, PetscInt* nleaves, const PetscInt* ilocal[], PetscInt* gremote[])
"""
	layout::PetscLayout,nleaves::PetscInt,ilocal::Union{Nothing, Vector{PetscInt}},gremote::Vector{PetscInt} = PetscSFGetGraphLayout(petsclib::PetscLibType,sf::PetscSF)
Get the global indices and `PetscLayout` that describe this star forest

Collective

Input Parameter:
- `sf` - star forest

Output Parameters:
- `layout`  - `PetscLayout` defining the global space for roots, new on every call: destroy it with `PetscLayoutDestroy`
- `nleaves` - number of leaf vertices on the current process, each of these references a root on any process
- `ilocal`  - locations of leaves in leafdata buffers, borrowed from `sf`, or `nothing` for contiguous storage
- `gremote` - root vertices in global numbering corresponding to the leaves, a Julia copy

Level: intermediate

-seealso: `PetscSF`, `PetscSFSetGraphLayout()`, `PetscSFCreate()`, `PetscSFView()`, `PetscSFSetGraph()`, `PetscSFGetGraph()`

# External Links
$(_doc_external("Vec/PetscSFGetGraphLayout"))
"""
function PetscSFGetGraphLayout(petsclib::PetscLibType, sf::PetscSF) end

@for_petsc function PetscSFGetGraphLayout(petsclib::$UnionPetscLib, sf::PetscSF )
	layout_ = Ref{PetscLayout}(C_NULL)
	nleaves_ = Ref{$PetscInt}()
	ilocal_ = Ref{Ptr{$PetscInt}}(C_NULL)
	gremote_ = Ref{Ptr{$PetscInt}}(C_NULL)

    @chk ccall(
               (:PetscSFGetGraphLayout, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{PetscLayout}, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}),
               sf, layout_, nleaves_, ilocal_, gremote_,
              )

	nleaves = nleaves_[]
	# ilocal == NULL means the leaves are contiguous [0, nleaves)
	ilocal = ilocal_[] == C_NULL ? nothing : unsafe_wrap(Array, ilocal_[], max(nleaves, 0); own = false)
	# PETSc allocates gremote for the caller: copy it and free PETSc's array
	gremote = copy(unsafe_wrap(Array, gremote_[], max(nleaves, 0); own = false))
	PetscFree(petsclib, gremote_[])

	return layout_[],nleaves,ilocal,gremote
end
