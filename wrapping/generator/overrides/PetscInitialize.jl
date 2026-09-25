# override for PetscInitialize; C signature: PetscInitialize(int* argc, char*** args, const char file[], const char help[])
"""
	PetscInitialize(petsclib::PetscLibType, args::Vector{String}, file::Union{Nothing, String} = nothing, help::Union{Nothing, String} = nothing)
Initializes the PETSc database and MPI, with `args` as the command line: `args[1]` is the program name and the rest are
options, one token per entry, such as `["prog", "-ksp_type", "cg"]`. Options given here override those in the
`PETSC_OPTIONS` environment variable.

PETSc keeps pointers to `args` until `PetscFinalize()`, so the strings are held here until the next call.

Collective on `MPI_COMM_WORLD` or `PETSC_COMM_WORLD` if it has been set

Input Parameters:
- `args` - the command line, program name first
- `file` - an optional options file, or `nothing`
- `help` - an optional help message printed by `-help`, or `nothing`

Level: beginner

See also: `PetscFinalize()`, `PetscInitializeNoArguments()`, `PetscGetArgs()`

# External Links
$(_doc_external("Sys/PetscInitialize"))
"""
function PetscInitialize(petsclib::PetscLibType, args::Vector{String}, file::Union{Nothing, String} = nothing, help::Union{Nothing, String} = nothing) end

# the command line each library was last initialized with, kept alive for PETSc
const _petsc_initialize_args = IdDict{Any, Any}()

@for_petsc function PetscInitialize(petsclib::$UnionPetscLib, args::Vector{String}, file::Union{Nothing, String} = nothing, help::Union{Nothing, String} = nothing)
	isempty(args) && throw(ArgumentError("args must start with a program name"))
	args = copy(args)
	argv = Ptr{Cchar}[pointer(a) for a in args]
	push!(argv, C_NULL)
	_petsc_initialize_args[petsclib] = (args, argv)
	argc_ = Ref{Cint}(length(args))
	args_ = Ref{Ptr{Ptr{Cchar}}}(pointer(argv))
	@chk ccall(
		(:PetscInitialize, $petsc_library),
		PetscErrorCode,
		(Ptr{Cint}, Ptr{Ptr{Ptr{Cchar}}}, Ptr{Cchar}, Ptr{Cchar}),
		argc_, args_, something(file, C_NULL), something(help, C_NULL),
	)
	return nothing
end
