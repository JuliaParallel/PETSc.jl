# override for PetscSFGetGraphLayout; C signature: PetscSFGetGraphLayout(PetscSF sf, PetscLayout* layout, PetscInt* nleaves, PetscInt* ilocal[], PetscInt* gremote[])
"""
	nleaves::PetscInt,iloc::Vector{PetscInt},gremote::Vector{PetscInt} = PetscSFGetGraphLayout(petsclib::PetscLibType,sf::PetscSF, layout::PetscLayout) 
Get the global indices and `PetscLayout` that describe this star forest

Collective

Input Parameter:
- `sf` - star forest

Output Parameters:
- `layout`  - `PetscLayout` defining the global space for roots
- `nleaves` - number of leaf vertices on the current process, each of these references a root on any process
- `ilocal`  - locations of leaves in leafdata buffers, or `NULL` for contiguous storage
- `gremote` - root vertices in global numbering corresponding to leaves in ilocal

Level: intermediate

-seealso: `PetscSF`, `PetscSFSetGraphLayout()`, `PetscSFCreate()`, `PetscSFView()`, `PetscSFSetGraph()`, `PetscSFGetGraph()`

# External Links
$(_doc_external("Vec/PetscSFGetGraphLayout"))
"""
function PetscSFGetGraphLayout(petsclib::PetscLibType, sf::PetscSF, layout::PetscLayout) end

@for_petsc function PetscSFGetGraphLayout(petsclib::$UnionPetscLib, sf::PetscSF, layout::PetscLayout )
	nleaves_ = Ref{$PetscInt}()
	iloc_ = Ref{Ptr{$PetscInt}}()
	gremote_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PetscSFGetGraphLayout, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{PetscLayout}, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}),
               sf, layout, nleaves_, iloc_, gremote_,
              )

	nleaves = nleaves_[]
	iloc = unsafe_wrap(Array, iloc_[], VecGetLocalSize(petsclib, x); own = false)
	gremote = unsafe_wrap(Array, gremote_[], VecGetLocalSize(petsclib, x); own = false)

	return nleaves,iloc,gremote
end

