"""
	t::PetscToken = PetscTokenCreate(petsclib::PetscLibType, a::String, b::Cchar) 
Creates a `PetscToken` used to find tokens in a string

Not Collective; No Fortran Support

Input Parameters:
- `a` - the string to look in
- `b` - the separator character

Output Parameter:
- `t` - the token object

Level: intermediate

See also: `PetscToken`, `PetscTokenFind()`, `PetscTokenDestroy()`

# External Links
$(_doc_external("Sys/PetscTokenCreate"))
"""
function PetscTokenCreate(petsclib::PetscLibType, a::String, b::Cchar)
    error("PetscTokenCreate: no generated method for these argument types")
end

@for_petsc function PetscTokenCreate(petsclib::$UnionPetscLib, a::String, b::Cchar )
	t_ = Ref{PetscToken}()

    @chk ccall(
               (:PetscTokenCreate, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Cchar, Ptr{PetscToken}),
               a, b, t_,
              )

	t = t_[]

	return t
end 

"""
	PetscTokenDestroy(petsclib::PetscLibType, a::Union{PetscToken, Ref{PetscToken}}) 
Destroys a `PetscToken`

Not Collective; No Fortran Support

Input Parameter:
- `a` - pointer to token

Level: intermediate

See also: `PetscToken`, `PetscTokenCreate()`, `PetscTokenFind()`

# External Links
$(_doc_external("Sys/PetscTokenDestroy"))
"""
function PetscTokenDestroy(petsclib::PetscLibType, a::Union{PetscToken, Ref{PetscToken}})
    error("PetscTokenDestroy: no generated method for these argument types")
end

@for_petsc function PetscTokenDestroy(petsclib::$UnionPetscLib, a::Union{PetscToken, Ref{PetscToken}} )
	a_ = a isa Base.RefValue ? a : Ref{PetscToken}(a)

    @chk ccall(
               (:PetscTokenDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscToken},),
               a_,
              )


	return nothing
end 

"""
	result::String = PetscTokenFind(petsclib::PetscLibType, a::PetscToken) 
Locates next "token" in a `PetscToken`

Not Collective; No Fortran Support

Input Parameter:
- `a` - pointer to token

Output Parameter:
- `result` - location of occurrence, `NULL` if not found

Level: intermediate

See also: `PetscToken`, `PetscTokenCreate()`, `PetscTokenDestroy()`

# External Links
$(_doc_external("Sys/PetscTokenFind"))
"""
function PetscTokenFind(petsclib::PetscLibType, a::PetscToken)
    error("PetscTokenFind: no generated method for these argument types")
end

@for_petsc function PetscTokenFind(petsclib::$UnionPetscLib, a::PetscToken )
	result_ = Ref{Ptr{Cchar}}()

    @chk ccall(
               (:PetscTokenFind, $petsc_library),
               PetscErrorCode,
               (PetscToken, Ptr{Ptr{Cchar}}),
               a, result_,
              )

	result = result_[] == C_NULL ? "" : unsafe_string(result_[])

	return result
end 

