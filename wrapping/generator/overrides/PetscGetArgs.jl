# override for PetscGetArgs; C signature: PetscGetArgs(int* argc, char*** args)
"""
	argc::Cint,args::Vector{String} = PetscGetArgs(petsclib::PetscLibType) 
Allows you to access the raw command line arguments anywhere
after PetscInitialize() is called but before `PetscFinalize()`.

Not Collective, No Fortran Support

Output Parameters:
- `argc` - count of number of command line arguments
- `args` - the command line arguments

Level: intermediate

Note:
The arguments are copied into Julia `String`s. PETSc keeps ownership of its own copy, which is not freed here.

See also: `PetscFinalize()`, `PetscInitializeFortran()`, `PetscGetArguments()`, `PetscInitialize()`

# External Links
$(_doc_external("Sys/PetscGetArgs"))
"""
function PetscGetArgs(petsclib::PetscLibType) end

@for_petsc function PetscGetArgs(petsclib::$UnionPetscLib)
	argc_ = Ref{Cint}()
	args_ = Ref{Ptr{Ptr{Cchar}}}(C_NULL)

    @chk ccall(
               (:PetscGetArgs, $petsc_library),
               PetscErrorCode,
               (Ptr{Cint}, Ptr{Ptr{Ptr{Cchar}}}),
               argc_, args_,
              )

	argc = argc_[]
	# PETSc keeps ownership of its copy of argv
	args = String[unsafe_string(unsafe_load(args_[], i)) for i in 1:(args_[] == C_NULL ? 0 : argc)]

	return argc,args
end
