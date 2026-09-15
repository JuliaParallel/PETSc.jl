# override for DMStagGetOwnershipRanges; C signature: DMStagGetOwnershipRanges(DM dm, PetscInt* lx[], PetscInt* ly[], PetscInt* lz[])
"""
	lx::Vector{PetscInt},ly::Vector{PetscInt},lz::Vector{PetscInt} = DMStagGetOwnershipRanges(petsclib::PetscLibType,dm::AbstractPetscDM) 
get elements per rank in each direction

Not Collective

Input Parameter:
- `dm` - the `DMSTAG` object

Output Parameters:
- `lx` - ownership along x direction (optional)
- `ly` - ownership along y direction (optional)
- `lz` - ownership along z direction (optional)

Level: intermediate

Notes:
These correspond to the optional final arguments passed to `DMStagCreate1d()`, `DMStagCreate2d()`, and `DMStagCreate3d()`.

Arguments corresponding to higher dimensions are ignored for 1D and 2D grids. These arguments may be set to `NULL` in this case.

In C you should not free these arrays, nor change the values in them.
They will only have valid values while the `DMSTAG` they came from still exists (has not been destroyed).

See also: 
=== 
`DMSTAG`, `DMStagSetGlobalSizes()`, `DMStagSetOwnershipRanges()`, `DMStagCreate1d()`, `DMStagCreate2d()`, `DMStagCreate3d()`, `DMDAGetOwnershipRanges()`

# External Links
$(_doc_external("DMStag/DMStagGetOwnershipRanges"))
"""
function DMStagGetOwnershipRanges(petsclib::PetscLibType, dm::AbstractPetscDM) end

@for_petsc function DMStagGetOwnershipRanges(petsclib::$UnionPetscLib, dm::AbstractPetscDM )
	lx_ = Ref{Ptr{$PetscInt}}()
	ly_ = Ref{Ptr{$PetscInt}}()
	lz_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:DMStagGetOwnershipRanges, $petsc_library),
               PetscErrorCode,
               (CDM, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}),
               dm, lx_, ly_, lz_,
              )
    # todo: determine the sizes of these arrays to wrap them properly
	#lx = unsafe_wrap(Array, lx_[], VecGetLocalSize(petsclib, x); own = false)
	#ly = unsafe_wrap(Array, ly_[], VecGetLocalSize(petsclib, x); own = false)
	#lz = unsafe_wrap(Array, lz_[], VecGetLocalSize(petsclib, x); own = false)
    lx = lx_[]
    ly = ly_[]
    lz = lz_[]
	return lx,ly,lz
end

