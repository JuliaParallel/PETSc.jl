# override for PCASMDestroySubdomains; C signature: PCASMDestroySubdomains(PetscInt n, IS* is[], IS* is_local[])
"""
	PCASMDestroySubdomains(petsclib::PetscLibType, n::PetscInt, is::Union{Ptr, AbstractVector{<:AbstractIS}}, is_local = nothing)
Destroys the index sets created with
`PCASMCreateSubdomains()`. Should be called after setting subdomains with `PCASMSetLocalSubdomains()`.

Collective

Input Parameters:
- `n`        - the number of index sets
- `is`       - the vector of index sets, or the raw pointer to a PETSc array of them
- `is_local` - the vector of local index sets, or the raw pointer to a PETSc array of them, can be `nothing`

With vectors, each index set is destroyed and left with a null pointer: the creators have already freed PETSc's array.
With pointers, PETSc destroys the index sets and frees both arrays.

Level: advanced

See also: `PCASM`, `PCASMCreateSubdomains()`, `PCASMSetLocalSubdomains()`

# External Links
$(_doc_external("PC/PCASMDestroySubdomains"))
"""
function PCASMDestroySubdomains(petsclib::PetscLibType, n::Integer, is::Union{Ptr, AbstractVector{<:AbstractIS}}, is_local = nothing) end

@for_petsc function PCASMDestroySubdomains(petsclib::$UnionPetscLib, n::$PetscInt, is::Union{Ptr, AbstractVector{<:AbstractIS}}, is_local = nothing)
	if is isa AbstractVector
		is_local isa Union{Nothing, AbstractVector{<:AbstractIS}} ||
			throw(ArgumentError("is is a vector, so is_local must be a vector or nothing"))
		destroy_index_sets(petsclib, n, is)
		is_local === nothing || isempty(is_local) || destroy_index_sets(petsclib, n, is_local)
		return nothing
	end
	is_local isa Union{Nothing, Ptr} || throw(ArgumentError("is is a pointer, so is_local must be a pointer or nothing"))
	is_ = Ref{Ptr{CIS}}(is)
	is_local_ = Ref{Ptr{CIS}}(is_local === nothing ? C_NULL : is_local)

	@chk ccall(
		(:PCASMDestroySubdomains, $petsc_library),
		PetscErrorCode,
		($PetscInt, Ptr{Ptr{CIS}}, Ptr{Ptr{CIS}}),
		n, is_, is_local_,
	)

	return nothing
end
