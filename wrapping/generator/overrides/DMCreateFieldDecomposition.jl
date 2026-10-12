# override for DMCreateFieldDecomposition; C signature: DMCreateFieldDecomposition(DM dm, PetscInt* len, char*** namelist, IS* islist[], DM* dmlist[])
"""
	len::PetscInt,namelist::Vector{String},islist::Vector{IS},dmlist::Vector{PetscDM} = DMCreateFieldDecomposition(petsclib::PetscLibType,dm::AbstractPetscDM) 
Returns a list of `IS` objects defining a decomposition of a problem into subproblems
corresponding to different fields.

Not Collective; No Fortran Support

Input Parameter:
- `dm` - the `DM` object

Output Parameters:
- `len`      - The number of fields (or `NULL` if not requested)
- `namelist` - The name for each field (or `NULL` if not requested)
- `islist`   - The global indices for each field (or `NULL` if not requested)
- `dmlist`   - The `DM`s for each field subproblem (or `NULL`, if not requested; if `NULL` is returned, no `DM`s are defined)

Level: intermediate

Note:
The names come back as Julia `String`s and the C arrays are freed here. Every entry of `islist` and `dmlist` is owned
by the caller and should be destroyed with `ISDestroy()` and `DMDestroy()` (or `PETSc.destroy!`).

See also: `DM`, `DMAddField()`, `DMCreateFieldIS()`, `DMCreateSubDM()`, `DMCreateDomainDecomposition()`, `DMDestroy()`, `DMView()`, `DMCreateInterpolation()`, `DMCreateColoring()`, `DMCreateMatrix()`, `DMCreateMassMatrix()`, `DMRefine()`, `DMCoarsen()`

# External Links
$(_doc_external("DM/DMCreateFieldDecomposition"))
"""
function DMCreateFieldDecomposition(petsclib::PetscLibType, dm::AbstractPetscDM) end

@for_petsc function DMCreateFieldDecomposition(petsclib::$UnionPetscLib, dm::AbstractPetscDM )
	len_ = Ref{$PetscInt}()
	namelist_ = Ref{Ptr{Ptr{Cchar}}}(C_NULL)
	islist_ = Ref{Ptr{CIS}}(C_NULL)
	dmlist_ = Ref{Ptr{CDM}}(C_NULL)
    @chk ccall(
               (:DMCreateFieldDecomposition, $petsc_library),
               PetscErrorCode,
               (CDM, Ptr{$PetscInt}, Ptr{Ptr{Ptr{Cchar}}}, Ptr{Ptr{CIS}}, Ptr{Ptr{CDM}}),
               dm, len_, namelist_, islist_, dmlist_,
              )
	len = len_[]
	# the caller owns the names, the three arrays, every IS and every DM: copy the names, free the C memory
	namelist = String[]
	if namelist_[] != C_NULL
		for i in 1:len
			name_ = unsafe_load(namelist_[], i)
			push!(namelist, unsafe_string(name_))
			PetscFree(petsclib, name_)
		end
		PetscFree(petsclib, namelist_[])
	end
	islist = islist_[] == C_NULL ? IS{$PetscLib}[] : [IS(p, petsclib) for p in unsafe_wrap(Array, islist_[], len; own = false)]
	islist_[] == C_NULL || PetscFree(petsclib, islist_[])
	dmlist = dmlist_[] == C_NULL ? PetscDM{$PetscLib}[] : [PetscDM(p, petsclib) for p in unsafe_wrap(Array, dmlist_[], len; own = false)]
	dmlist_[] == C_NULL || PetscFree(petsclib, dmlist_[])
	return len,namelist,islist,dmlist
end
