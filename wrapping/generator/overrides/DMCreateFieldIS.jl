# override for DMCreateFieldIS; C signature: DMCreateFieldIS(DM dm, PetscInt* numFields, char** fieldNames, IS* fields[])
"""
	numFields::PetscInt,fieldNames::Cchar,fields::Vector{IS} = DMCreateFieldIS(petsclib::PetscLibType,dm::AbstractPetscDM) 
Creates a set of `IS` objects with the global indices of dofs for each field defined with `DMAddField()`

Not Collective; No Fortran Support

Input Parameter:
- `dm` - the `DM` object

Output Parameters:
- `numFields`  - The number of fields (or `NULL` if not requested)
- `fieldNames` - The name of each field (or `NULL` if not requested)
- `fields`     - The global indices for each field (or `NULL` if not requested)

Level: intermediate

Note:
The user is responsible for freeing all requested arrays. In particular, every entry of `fieldNames` should be freed with
`PetscFree()`, every entry of `fields` should be destroyed with `ISDestroy()`, and both arrays should be freed with
`PetscFree()`.

Developer Note:
It is not clear why both this function and `DMCreateFieldDecomposition()` exist. Having two seems redundant and confusing. This function should
likely be removed.

See also: 
=== 
`DM`, `DMAddField()`, `DMGetField()`, `DMDestroy()`, `DMView()`, `DMCreateInterpolation()`, `DMCreateColoring()`, `DMCreateMatrix()`,
`DMCreateFieldDecomposition()`

# External Links
$(_doc_external("DM/DMCreateFieldIS"))
"""
function DMCreateFieldIS(petsclib::PetscLibType, dm::AbstractPetscDM) end

@for_petsc function DMCreateFieldIS(petsclib::$UnionPetscLib, dm::AbstractPetscDM )
	numFields_ = Ref{$PetscInt}()
	fieldNames_ = Ref{Ptr{Ptr{Cchar}}}(C_NULL)
	fields_ = Ref{Ptr{CIS}}(C_NULL)
    @chk ccall(
               (:DMCreateFieldIS, $petsc_library),
               PetscErrorCode,
               (CDM, Ptr{$PetscInt}, Ptr{Ptr{Ptr{Cchar}}}, Ptr{Ptr{CIS}}),
               dm, numFields_, fieldNames_, fields_,
              )
	numFields = numFields_[]
	fieldNames = String[]
	if fieldNames_[] != C_NULL
		for i in 1:numFields
			push!(fieldNames, unsafe_string(unsafe_load(fieldNames_[], i)))
		end
	end
	fields = IS{$PetscLib}[]
	if fields_[] != C_NULL
		for i in 1:numFields
			push!(fields, IS(unsafe_load(fields_[], i), petsclib))
		end
	end
	return numFields,fieldNames,fields
end

