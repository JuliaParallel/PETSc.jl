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

"""
    SetPetscLib(library_path::String; PetscScalar=Float64, PetscInt=Int64)

Create a custom PETSc library instance from a user-specified shared library path.

This function allows you to use a custom-compiled PETSc library instead of the 
pre-built libraries provided by PETSc_jll. The custom library must be:
- Compiled as a shared/dynamic library (not static)
- Built with the correct scalar type (real vs complex) and integer size
- Compatible with the same MPI installation that MPI.jl is using

# Arguments
- `library_path::String`: Path to the custom PETSc shared library (e.g., "/path/to/libpetsc.so")
- `PetscScalar::Type`: Scalar type used by the library. Options: `Float64`, `Float32`, 
  `Complex{Float64}`, or `Complex{Float32}`. Default: `Float64`
- `PetscInt::Type`: Integer type used by the library. Options: `Int32` or `Int64`. 
  Default: `Int64`

# Returns
- A `PetscLibType` instance that can be used with `initialize()`, `finalize()`, and 
  all other PETSc.jl functions

# Examples
```julia
using PETSc

# For a custom double-precision real PETSc library with 64-bit indices
petsclib = PETSc.SetPetscLib("/path/to/custom/libpetsc.so"; 
                             PetscScalar=Float64, PetscInt=Int64)

# For a single-precision complex library with 32-bit indices  
petsclib = PETSc.SetPetscLib("/opt/petsc/lib/libpetsc.so"; 
                             PetscScalar=Complex{Float32}, PetscInt=Int32)

# Initialize and use the custom library
PETSc.initialize(petsclib)
# ... your code ...
PETSc.finalize(petsclib)
```

# See Also
- [`initialize`](@ref): Initialize a PETSc library
- [`finalize`](@ref): Finalize a PETSc library
"""
function SetPetscLib(library_path::String; PetscScalar::Type=Float64, PetscInt::Type=Int64)
    return LibPETSc.PetscLibType{PetscScalar, PetscInt}(library_path)
end


"""
    check_petsc_wrappers_version(petsclib=nothing)

Load the generated `petsc_wrappers_version.jl` (if present) and compare the
declared wrapper version `PETSC_WRAPPERS_VERSION` with the installed PETSc
version obtained from `LibPETSc.PetscGetVersionNumber` for `petsclib`.

Arguments
- `petsclib`: optional `PetscLibType` or path string. If `nothing`, the
  first available `PETSc.petsclibs[1]` is used.

Returns a named tuple: `(:wrappers_version, :installed_version, :match)`.
`match` is `true` when versions are equal, `false` when they differ, and
`nothing` if either side could not be determined.
"""
function check_petsc_wrappers_version(petsclib=nothing)
    verfile = joinpath(@__DIR__, "autowrapped", "petsc_wrappers_version.jl")

    if !isdefined(@__MODULE__, :PETSC_WRAPPERS_VERSION) && isfile(verfile)
        try
            include(verfile)
        catch err
            @warn "Failed to include petsc_wrappers_version.jl" exception=(err,)
        end
    end

    wrappers_version = isdefined(@__MODULE__, :PETSC_WRAPPERS_VERSION) ? PETSC_WRAPPERS_VERSION : nothing

    if petsclib === nothing
        if isdefined(@__MODULE__, :petsclibs) && !isempty(petsclibs)
            petsclib = petsclibs[1]
        else
            error("No PETSc library available to check installed version")
        end
    end

    if isa(petsclib, String)
        petsclib = SetPetscLib(petsclib)
    end

    installed_version = nothing
    try
        major, minor, subminor, _release = LibPETSc.PetscGetVersionNumber(petsclib)
        installed_version = VersionNumber(Int(major), Int(minor), Int(subminor))
    catch err
        @warn "Failed to query installed PETSc version" exception=(err,)
    end

    match = isnothing(wrappers_version) || isnothing(installed_version) ? nothing : (installed_version == wrappers_version)
    if !isnothing(match) && match === false
        @warn "PETSc wrappers version does not match PETSc version of library; this can cause undesired behavior" wrappers_version=wrappers_version installed_version=installed_version
    end

    return (wrappers_version = wrappers_version, installed_version = installed_version, match = match)
end