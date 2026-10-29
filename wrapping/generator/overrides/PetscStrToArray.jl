# override for PetscStrToArray; C signature: PetscStrToArray(const char s[], char sp, int* argc, char*** args)
"""
	args::Vector{String} = PetscStrToArray(petsclib::PetscLibType, s::String, sp::Cchar) 
Converts a string to an array of strings by splitting it at every `sp`.

Not Collective, No Fortran Support

Input Parameters:
- `s`  - pointer to string
- `sp` - separator character

Output Parameter:
- `args` - the pieces of `s`

Level: intermediate

Note:
The pieces are copied into Julia `String`s and the C memory is freed here, so `PetscStrToArrayDestroy()` is not needed.

See also: `PetscStrToArrayDestroy()`, `PetscToken`, `PetscTokenCreate()`

# External Links
$(_doc_external("Sys/PetscStrToArray"))
"""
function PetscStrToArray(petsclib::PetscLibType, s::String, sp::Cchar) end

@for_petsc function PetscStrToArray(petsclib::$UnionPetscLib, s::String, sp::Cchar )
	argc_ = Ref{Cint}()
	args_ = Ref{Ptr{Ptr{Cchar}}}(C_NULL)

    @chk ccall(
               (:PetscStrToArray, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Cchar, Ptr{Cint}, Ptr{Ptr{Ptr{Cchar}}}),
               s, sp, argc_, args_,
              )

	# PETSc allocates with C malloc here, not PetscMalloc, so PetscStrToArrayDestroy (C free) releases it
	args = String[unsafe_string(unsafe_load(args_[], i)) for i in 1:(args_[] == C_NULL ? 0 : argc_[])]
	args_[] == C_NULL || @chk ccall(
		(:PetscStrToArrayDestroy, $petsc_library),
		PetscErrorCode,
		(Cint, Ptr{Ptr{Cchar}}),
		argc_[], args_[],
	)

	return args
end
