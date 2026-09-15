# override for PetscSFGetGraph; C signature: PetscSFGetGraph(PetscSF sf, PetscInt* nroots, PetscInt* nleaves, PetscInt* ilocal[], PetscSFNode* iremote[])
"""
	nroots::PetscInt,nleaves::PetscInt,iloc::Vector{PetscInt} = PetscSFGetGraph(petsclib::PetscLibType,sf::PetscSF, iremote::Vector{PetscSFNode}) 
Get the graph specifying a parallel star forest

Not Collective

Input Parameter:
- `sf` - star forest

Output Parameters:
- `nroots`  - number of root vertices on the current process (these are possible targets for other process to attach leaves)
- `nleaves` - number of leaf vertices on the current process, each of these references a root on any process
- `ilocal`  - locations of leaves in leafdata buffers (if returned value is `NULL`, it means leaves are in contiguous storage)
- `iremote` - remote locations of root vertices for each leaf on the current process

Level: intermediate

-seealso: `PetscSF`, `PetscSFType`, `PetscSFCreate()`, `PetscSFView()`, `PetscSFSetGraph()`

# External Links
$(_doc_external("Vec/PetscSFGetGraph"))
"""
function PetscSFGetGraph(petsclib::PetscLibType, sf::PetscSF) end

@for_petsc function PetscSFGetGraph(petsclib::$UnionPetscLib, sf::PetscSF )
	nroots_ = Ref{$PetscInt}()
	nleaves_ = Ref{$PetscInt}()
	ilocal_ = Ref{Ptr{$PetscInt}}(C_NULL)
	iremote_ = Ref{Ptr{PetscSFNode}}(C_NULL)
    @chk ccall(
               (:PetscSFGetGraph, $petsc_library),
               PetscErrorCode,
               (PetscSF, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{PetscSFNode}}),
               sf, nroots_, nleaves_, ilocal_, iremote_,
              )
	nroots = nroots_[]
	nleaves = nleaves_[]
	# ilocal == NULL means the leaves are contiguous [0, nleaves)
	ilocal = ilocal_[] == C_NULL ? nothing : unsafe_wrap(Array, ilocal_[], max(nleaves, 0); own = false)
	iremote = iremote_[] == C_NULL ? nothing : unsafe_wrap(Array, iremote_[], max(nleaves, 0); own = false)
	return nroots,nleaves,ilocal,iremote
end

