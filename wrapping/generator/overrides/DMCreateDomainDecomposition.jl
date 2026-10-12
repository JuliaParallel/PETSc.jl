# override for DMCreateDomainDecomposition; C signature: DMCreateDomainDecomposition(DM dm, PetscInt* n, char*** namelist, IS* innerislist[], IS* outerislist[], DM* dmlist[])
"""
	n::PetscInt,namelist::Vector{String},innerislist::Vector{IS},outerislist::Vector{IS},dmlist::Vector{PetscDM} = DMCreateDomainDecomposition(petsclib::PetscLibType,dm::AbstractPetscDM) 
Returns lists of `IS` objects defining a decomposition of a
problem into subproblems corresponding to restrictions to pairs of nested subdomains.

Not Collective

Input Parameter:
- `dm` - the `DM` object

Output Parameters:
- `n`           - The number of subproblems in the domain decomposition (or `NULL` if not requested), also the length of the four arrays below
- `namelist`    - The name for each subdomain (or `NULL` if not requested)
- `innerislist` - The global indices for each inner subdomain (or `NULL`, if not requested)
- `outerislist` - The global indices for each outer subdomain (or `NULL`, if not requested)
- `dmlist`      - The `DM`s for each subdomain subproblem (or `NULL`, if not requested; if `NULL` is returned, no `DM`s are defined)

Level: intermediate

Note:
The names come back as Julia `String`s and the C arrays are freed here. Every entry of `innerislist`, `outerislist`
and `dmlist` is owned by the caller and should be destroyed with `ISDestroy()` and `DMDestroy()` (or `PETSc.destroy!`).

See also: `DM`, `DMCreateFieldDecomposition()`, `DMDestroy()`, `DMCreateDomainDecompositionScatters()`, `DMView()`, `DMCreateInterpolation()`,
`DMSubDomainHookAdd()`, `DMSubDomainHookRemove()`, `DMCreateColoring()`, `DMCreateMatrix()`, `DMCreateMassMatrix()`, `DMRefine()`, `DMCoarsen()`

# External Links
$(_doc_external("DM/DMCreateDomainDecomposition"))
"""
function DMCreateDomainDecomposition(petsclib::PetscLibType, dm::AbstractPetscDM) end

@for_petsc function DMCreateDomainDecomposition(petsclib::$UnionPetscLib, dm::AbstractPetscDM )
	n_ = Ref{$PetscInt}()
	namelist_ = Ref{Ptr{Ptr{Cchar}}}(C_NULL)
	innerislist_ = Ref{Ptr{CIS}}()
	outerislist_ = Ref{Ptr{CIS}}()
	dmlist_ = Ref{Ptr{CDM}}()

    @chk ccall(
               (:DMCreateDomainDecomposition, $petsc_library),
               PetscErrorCode,
               (CDM, Ptr{$PetscInt}, Ptr{Ptr{Ptr{Cchar}}}, Ptr{Ptr{CIS}}, Ptr{Ptr{CIS}}, Ptr{Ptr{CDM}}),
               dm, n_, namelist_, innerislist_, outerislist_, dmlist_,
              )

	n = n_[]
	# the caller owns the names and the name array: copy the names, free the C memory
	namelist = String[]
	if namelist_[] != C_NULL
		for i in 1:n
			name_ = unsafe_load(namelist_[], i)
			push!(namelist, unsafe_string(name_))
			PetscFree(petsclib, name_)
		end
		PetscFree(petsclib, namelist_[])
	end
	innerislist = innerislist_[] == C_NULL ? IS{$PetscLib}[] : [IS(p, petsclib) for p in unsafe_wrap(Array, innerislist_[], n; own = false)]
	innerislist_[] == C_NULL || PetscFree(petsclib, innerislist_[])
	outerislist = outerislist_[] == C_NULL ? IS{$PetscLib}[] : [IS(p, petsclib) for p in unsafe_wrap(Array, outerislist_[], n; own = false)]
	outerislist_[] == C_NULL || PetscFree(petsclib, outerislist_[])
	dmlist = dmlist_[] == C_NULL ? PetscDM{$PetscLib}[] : [PetscDM(p, petsclib) for p in unsafe_wrap(Array, dmlist_[], n; own = false)]
	dmlist_[] == C_NULL || PetscFree(petsclib, dmlist_[])

	return n,namelist,innerislist,outerislist,dmlist
end
