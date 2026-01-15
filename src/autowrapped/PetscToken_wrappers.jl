# autodefined type arguments for class ------
mutable struct _n_PetscToken end
const PetscToken = Ptr{_n_PetscToken}

# -------------------------------------------------------
"""
	PetscTokenFind(petsclib::PetscLibType,a::PetscToken, result::String) 
Locates next "token" in a `PetscToken`

Not Collective; No Fortran Support

Input Parameter:
- `a` - pointer to token

Output Parameter:
- `result` - location of occurrence, `NULL` if not found

Level: intermediate

-seealso: `PetscToken`, `PetscTokenCreate()`, `PetscTokenDestroy()`

# External Links
$(_doc_external("Sys/PetscTokenFind"))
"""
function PetscTokenFind(petsclib::PetscLibType, a::PetscToken, result::String) end

@for_petsc function PetscTokenFind(petsclib::$UnionPetscLib, a::PetscToken, result::String )
	result_ = Ref(pointer(result))

    @chk ccall(
               (:PetscTokenFind, $petsc_library),
               PetscErrorCode,
               (PetscToken, Ptr{Ptr{Cchar}}),
               a, result_,
              )


	return nothing
end 

"""
	t::PetscToken = PetscTokenCreate(petsclib::PetscLibType,a::String, b::Cchar) 
Creates a `PetscToken` used to find tokens in a string

Not Collective; No Fortran Support

Input Parameters:
- `a` - the string to look in
- `b` - the separator character

Output Parameter:
- `t` - the token object

Level: intermediate

-seealso: `PetscToken`, `PetscTokenFind()`, `PetscTokenDestroy()`

# External Links
$(_doc_external("Sys/PetscTokenCreate"))
"""
function PetscTokenCreate(petsclib::PetscLibType, a::String, b::Cchar) end

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
	PetscTokenDestroy(petsclib::PetscLibType,a::PetscToken) 
Destroys a `PetscToken`

Not Collective; No Fortran Support

Input Parameter:
- `a` - pointer to token

Level: intermediate

-seealso: `PetscToken`, `PetscTokenCreate()`, `PetscTokenFind()`

# External Links
$(_doc_external("Sys/PetscTokenDestroy"))
"""
function PetscTokenDestroy(petsclib::PetscLibType, a::PetscToken) end

@for_petsc function PetscTokenDestroy(petsclib::$UnionPetscLib, a::PetscToken )

    @chk ccall(
               (:PetscTokenDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscToken},),
               a,
              )


	return nothing
end 

