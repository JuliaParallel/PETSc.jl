# override for PCGASMDestroySubdomains; C signature: PCGASMDestroySubdomains(PetscInt n, IS* iis[], IS* ois[])
"""
	PCGASMDestroySubdomains(petsclib::PetscLibType, n::PetscInt, iis::Union{Ptr, AbstractVector{<:AbstractIS}}, ois = nothing)
Destroys the index sets created with
`PCGASMCreateSubdomains()` or `PCGASMCreateSubdomains2D()`. Should be
called after setting subdomains with `PCGASMSetSubdomains()`.

Collective

Input Parameters:
- `n`   - the number of index sets
- `iis` - the vector of inner subdomains, or the raw pointer to a PETSc array of them
- `ois` - the vector of outer subdomains, or the raw pointer to a PETSc array of them, can be `nothing`

With vectors, each index set is destroyed and left with a null pointer: the creators have already freed PETSc's array.
With pointers, PETSc destroys the index sets and frees both arrays.

Level: intermediate

See also: `PCGASM`, `PCGASMCreateSubdomains()`, `PCGASMSetSubdomains()`

# External Links
$(_doc_external("PC/PCGASMDestroySubdomains"))
"""
function PCGASMDestroySubdomains(petsclib::PetscLibType, n::Integer, iis::Union{Ptr, AbstractVector{<:AbstractIS}}, ois = nothing) end

@for_petsc function PCGASMDestroySubdomains(petsclib::$UnionPetscLib, n::$PetscInt, iis::Union{Ptr, AbstractVector{<:AbstractIS}}, ois = nothing)
	if iis isa AbstractVector
		ois isa Union{Nothing, AbstractVector{<:AbstractIS}} ||
			throw(ArgumentError("iis is a vector, so ois must be a vector or nothing"))
		destroy_index_sets(petsclib, n, iis)
		ois === nothing || isempty(ois) || destroy_index_sets(petsclib, n, ois)
		return nothing
	end
	ois isa Union{Nothing, Ptr} || throw(ArgumentError("iis is a pointer, so ois must be a pointer or nothing"))
	iis_ = Ref{Ptr{CIS}}(iis)
	ois_ = Ref{Ptr{CIS}}(ois === nothing ? C_NULL : ois)

	@chk ccall(
		(:PCGASMDestroySubdomains, $petsc_library),
		PetscErrorCode,
		($PetscInt, Ptr{Ptr{CIS}}, Ptr{Ptr{CIS}}),
		n, iis_, ois_,
	)

	return nothing
end
