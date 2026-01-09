"""
   initialized(petsclib)

Check if `petsclib` is initialized

# External Links
$(_doc_external("Sys/PetscInitialized"))
"""
initialized(petsclib) = LibPETSc.PetscInitialized(petsclib)

const _petsc_program_name = "petsc_julia"
const _lib_handles = IdDict{Any, Tuple{Any, Bool}}()

"""
    initialize([petsclib]; log_view = false, options = String[])

Initialize the `petsclib`. If no `petsclib` is given, all `PETSc.petsclibs`
will be initialized.

Additionally:
 - This will initialize MPI if it has not already been initialized.
 - It will disable the PETSc signal handler (via
   $(_petsc_link("Sys/PetscPopSignalHandler"))
 - Add an [`atexit`](https://docs.julialang.org/en/v1/base/base/#Base.atexit)
   hook to call [`PETSc.finalize`](@ref).

# Arguments
- `log_view::Bool = false`: Enable PETSc's `-log_view` performance logging.
  When enabled, PETSc will output performance statistics at finalization.
- `options::Vector{String} = String[]`: Additional PETSc command-line options.
  These are passed via the `PETSC_OPTIONS` environment variable.

# Examples
```julia
# Basic initialization
PETSc.initialize(petsclib)

# Enable performance logging to stdout
PETSc.initialize(petsclib; log_view = true)

# Write log to a file
PETSc.initialize(petsclib; log_view = true, options = [":logfile.txt"])

# Enable memory logging
PETSc.initialize(petsclib; log_view = true, options = [":logfile.txt", "-log_view_memory"])

# Pass custom PETSc options without logging
PETSc.initialize(petsclib; options = ["-malloc_debug", "-on_error_abort"])
```

# External Links
$(_doc_external("Sys/PetscInitializeNoArguments"))
"""
function initialize(; log_view::Bool = false, options = String[])
    map(petsclib -> initialize(petsclib; log_view, options), petsclibs)
    return nothing
end

function initialize(petsclib; log_view::Bool = false, options = String[])
    if !initialized(petsclib)
        if log_view || !isempty(options)
            cli_opts = _build_petsc_options(log_view, options)
            prev_opts = get(ENV, "PETSC_OPTIONS", "")
            ENV["PETSC_OPTIONS"] = isempty(prev_opts) ? cli_opts : "$prev_opts $cli_opts"
            try
                _ensure_mpi_initialized()
                petsclib.age += 1
                LibPETSc.PetscInitializeNoArguments(petsclib)
                _post_initialize(petsclib)
            finally
                if isempty(prev_opts)
                    delete!(ENV, "PETSC_OPTIONS")
                else
                    ENV["PETSC_OPTIONS"] = prev_opts
                end
            end
        else
            _ensure_mpi_initialized()
            petsclib.age += 1
            LibPETSc.PetscInitializeNoArguments(petsclib)
            _post_initialize(petsclib)
        end
    end
    return nothing
end

"""
   finalize(petsclib)

Finalize the `petsclib`, if no `petsclib` is given then all `PETSc.petsclibs`
will be finalized.

# External Links
$(_doc_external("Sys/PetscFinalize"))
"""
function finalize()
    map(finalize, petsclibs)
    return nothing
end

function finalize(petsclib)
    if !finalized(petsclib)
        petsclib.age += 1
        LibPETSc.PetscFinalize(petsclib)
    end
    return nothing
end

"""
   finalized(petsclib)

Check if `petsclib` is finalized

# External Links
$(_doc_external("Sys/PetscFinalized"))
"""
finalized(petsclib) = LibPETSc.PetscFinalized(petsclib)

function _build_petsc_options(log_view::Bool, options)
    opts = String[]
    if log_view
        push!(opts, "-log_view")
    end
    append!(opts, [String(opt) for opt in options])
    return join(opts, " ")
end

function _ensure_mpi_initialized()
    MPI.Initialized() || MPI.Init()
    return nothing
end

function _post_initialize(petsclib)
    # disable signal handler
    LibPETSc.PetscPopSignalHandler(petsclib)
    atexit(() -> finalize(petsclib))
    return nothing
end

function _ensure_library_handle(petsclib)
    return get!(_lib_handles, petsclib) do
        libref = petsclib.petsc_library
        if libref isa AbstractString
            return (Libdl.dlopen(libref), true)
        else
            return (libref, false)
        end
    end
end

function _library_ptr(lib_handle)
    if lib_handle isa Ptr{Cvoid}
        return lib_handle
    end
    try
        return Base.unsafe_convert(Ptr{Cvoid}, lib_handle)
    catch err
        throw(ArgumentError("Unsupported PETSc library handle type $(typeof(lib_handle))"))
    end
end

function _release_library_handle(petsclib)
    entry = pop!(_lib_handles, petsclib, nothing)
    isnothing(entry) && return nothing
    handle, owned = entry
    owned || return nothing
    Libdl.dlclose(_library_ptr(handle))
    return nothing
end


"""
    scalartype(petsclib::PetscLibType)

return the scalar type for the associated `petsclib`
"""
scalartype(::LibPETSc.PetscLibType{ST}) where {ST} = ST
scalartype(
    ::Type{PetscLib},
) where {PetscLib <: PetscLibType{ST}} where {ST} = ST


"""
    inttype(petsclib::PetscLibType)

return the int type for the associated `petsclib`
"""
inttype(::LibPETSc.PetscLibType{ST, IT}) where {ST, IT} = IT
inttype(
    ::Type{PetscLib},
) where {PetscLib <: PetscLibType{ST, IT}} where {ST, IT} = IT