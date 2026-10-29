# override for PetscGetArguments; C signature: PetscGetArguments(char*** args)
"""
	args::Vector{String} = PetscGetArguments(petsclib::PetscLibType) 
Gets the command line arguments, without the program name.

Not Collective, No Fortran Support

Output Parameter:
- `args` - the command line arguments

Level: developer

Note:
The arguments are copied into Julia `String`s and the C memory is freed here. `PetscFreeArguments()` is not needed (and
in PETSc 3.25 it frees every other string and then runs past the end of the array).

See also: `PetscFinalize()`, `PetscInitializeFortran()`, `PetscGetArgs()`, `PetscInitialize()`

# External Links
$(_doc_external("Sys/PetscGetArguments"))
"""
function PetscGetArguments(petsclib::PetscLibType) end

@for_petsc function PetscGetArguments(petsclib::$UnionPetscLib)
	args_ = Ref{Ptr{Ptr{Cchar}}}(C_NULL)

    @chk ccall(
               (:PetscGetArguments, $petsc_library),
               PetscErrorCode,
               (Ptr{Ptr{Ptr{Cchar}}},),
               args_,
              )

	# NULL-terminated array of malloc'd strings: copy, then free each string and the array
	args = String[]
	if args_[] != C_NULL
		i = 1
		while (arg_ = unsafe_load(args_[], i)) != C_NULL
			push!(args, unsafe_string(arg_))
			PetscFree(petsclib, arg_)
			i += 1
		end
		PetscFree(petsclib, args_[])
	end

	return args
end
