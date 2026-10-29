function ensure_library_handle(petsclib)
    return get!(lib_handles, petsclib) do
        libref = petsclib.petsc_library
        if libref isa AbstractString
            return (Libdl.dlopen(libref), true)
        else
            return (libref, false)
        end
    end
end

function library_ptr(lib_handle)
    if lib_handle isa Ptr{Cvoid}
        return lib_handle
    end
    try
        return Base.unsafe_convert(Ptr{Cvoid}, lib_handle)
    catch err
        throw(ArgumentError("Unsupported PETSc library handle type $(typeof(lib_handle))"))
    end
end

function release_library_handle(petsclib)
    entry = pop!(lib_handles, petsclib, nothing)
    isnothing(entry) && return nothing
    handle, owned = entry
    owned || return nothing
    Libdl.dlclose(library_ptr(handle))
    return nothing
end


# PETSc and LibPETSc share these functions: there is one of each
"""
    PetscLibType(library_path::String; PetscScalar=Float64, PetscInt=Int64)

Create a custom PETSc library instance from a user-specified shared library path.

Replaces v0.4's `set_petsclib`, which mutated nothing despite its name and so
was never a `set_*!` (docs/src/man/naming.md §7): it builds and returns a
library handle, which is what a constructor does.

This function allows you to use a custom-compiled PETSc library instead of the
pre-built libraries provided by `PETSc_jll`. The custom library must be compiled as a
shared/dynamic library (not static), built with the matching scalar type and integer
size, and linked against the same MPI installation that `MPI.jl` uses.

On HPC systems, set `JULIA_PETSC_SKIP_JLL=1` before starting Julia to prevent
`PETSc_jll` from being precompiled (its MPI stack is typically incompatible with
cluster MPI). Then call this function in your script to load the cluster library.

# Arguments
- `library_path::String`: Path to the PETSc shared library (e.g. `"/path/to/libpetsc.so"`)
- `PetscScalar::Type`: Scalar type the library was built with. One of `Float64`,
  `Float32`, `Complex{Float64}`, `Complex{Float32}`. Default: `Float64`
- `PetscInt::Type`: Integer type the library was built with. `Int32` or `Int64`.
  Default: `Int64`

# Returns
A `PetscLibType` instance for use with `initialize`, `finalize`, and all PETSc.jl functions.

# Environment-variable alternative
Instead of calling this constructor, you can configure everything before Julia starts:
```
JULIA_PETSC_LIBRARY=/path/to/libpetsc.so   # also suppresses PETSc_jll
JULIA_PETSC_SCALAR=Float64                  # Float32 | ComplexFloat64 | ComplexFloat32
JULIA_PETSC_INT=Int64                       # Int32
```
With these set, `PETSc.getlib(; PetscScalar=Float64, PetscInt=Int64)` returns the
custom library directly.

# Examples
```julia
# Double-precision real, 64-bit indices (typical HPC build)
petsclib = PETSc.LibPETSc.PetscLibType("/path/to/libpetsc.so";
                                       PetscScalar=Float64, PetscInt=Int64)
PETSc.initialize(petsclib)
# ... your code ...
PETSc.finalize(petsclib)

# Single-precision complex, 32-bit indices
petsclib = PETSc.LibPETSc.PetscLibType("/opt/petsc/lib/libpetsc.so";
                                       PetscScalar=Complex{Float32}, PetscInt=Int32)
```

# See Also
- [`initialize`](@ref): Initialize a PETSc library
- [`finalize`](@ref): Finalize a PETSc library
"""
function LibPETSc.PetscLibType(library_path::String; PetscScalar::Type=Float64, PetscInt::Type=Int64)
    petsclib = LibPETSc.PetscLibType{PetscScalar, PetscInt}(library_path)
    try
        check_wrappers_version(petsclib)
    catch err
        @warn "Failed to perform PETSc wrappers version check" exception=(err,)
    end
    return petsclib
end

# The resolved path of a loaded library, whether it was configured as a path or
# as a JLL handle.
function library_path_string(lib)
    return lib.petsc_library isa AbstractString ? lib.petsc_library :
           try
               Libdl.dlpath(Libdl.dlopen(lib.petsc_library))
           catch
               string(lib.petsc_library)
           end
end

"""
    library_info()

Report the PETSc library configuration as a `NamedTuple` (§12).

| Field | Meaning |
|---|---|
| `source` | `:preferences` if `LocalPreferences.toml` names a library, `:jll` for the bundled `PETSc_jll` binaries |
| `path` | the library path in use |
| `scalar` | `PetscScalar` of the preferred library |
| `int` | `PetscInt` of the preferred library |
| `real` | `PetscReal` of the preferred library |

The preferred library is `petsclibs[1]`: the configured one when a preference is
set, and the first of the bundled builds otherwise.

It prints as the `NamedTuple` it is. `PETSc.petsclibs` lists every library
loaded in this session.

```julia
info = library_info()
info.scalar          # Float64
```

# See Also
- [`set_library!`](@ref): configure a custom library persistently
- [`unset_library!`](@ref): revert to `PETSc_jll`
"""
function library_info()
    pref_path = @load_preference("library_path", nothing)
    lib = isempty(petsclibs) ? nothing : petsclibs[1]

    source = pref_path === nothing ? :jll : :preferences
    path = pref_path !== nothing ? pref_path :
           lib === nothing ? nothing : library_path_string(lib)

    return (
        source = source,
        path   = path,
        scalar = lib === nothing ? nothing : lib.PetscScalar,
        int    = lib === nothing ? nothing : lib.PetscInt,
        real   = lib === nothing ? nothing : lib.PetscReal,
    )
end

"""
    set_library!(path; PetscScalar=Float64, PetscInt=Int64)

Persistently configure PETSc.jl to use a custom PETSc shared library.

The path and type configuration are stored in `LocalPreferences.toml` (per-project,
git-ignorable) and take effect on the next Julia session. Recompilation is triggered
automatically — no environment variables are needed.

To revert to the default `PETSc_jll` libraries, call [`unset_library!`](@ref).

# Arguments
- `path`: path to the PETSc shared library (e.g. `"/path/to/libpetsc.so"`)
- `PetscScalar`: scalar type the library was built with (`Float64`, `Float32`,
  `Complex{Float64}`, `Complex{Float32}`). Default: `Float64`
- `PetscInt`: integer type the library was built with (`Int64` or `Int32`).
  Default: `Int64`

# Examples
```julia
PETSc.set_library!(
    "/project/petsc/lib/libpetsc.so";
    PetscScalar = Float64,
    PetscInt    = Int64,
)
# Restart Julia — the new library is used automatically from here on.
```

# See Also
- [`unset_library!`](@ref): remove the preference and revert to `PETSc_jll`
- `PetscLibType(path)`: load a custom library for the current session only
"""
function set_library!(path; PetscScalar::Type=Float64, PetscInt::Type=Int64)
    ispath(path) || throw(ArgumentError("PETSc library not found: $path"))
    @set_preferences!(
        "library_path" => realpath(path),
        "PetscScalar"  => string(PetscScalar),
        "PetscInt"     => string(PetscInt),
    )
    @info "PETSc library configured — restart Julia to use the new library." path PetscScalar PetscInt
    return nothing
end

"""
    unset_library!()

Remove the persistent custom-library preference set by [`set_library!`](@ref),
reverting to the default `PETSc_jll` binaries on the next Julia session.
"""
function unset_library!()
    @delete_preferences!("library_path", "PetscScalar", "PetscInt")
    @info "PETSc library preference removed — restart Julia to revert to PETSc_jll."
    return nothing
end

"""
    set_petscint!(::Type{T}) where {T<:Union{Int32,Int64}}

Persistently select which `PetscInt` width of the `PETSc_jll` libraries is loaded
(default `Int64`). Only one width is registered per process: the Int64 and Int32 library
variants link external packages (hypre, SuperLU_DIST) that export identical symbols with
different integer ABIs, so mixing both widths in one process is unsafe on platforms with a
flat dynamic-linker namespace. `PETSc.petsclibs` then holds the four scalar variants of
that width.

Takes effect on the next Julia session (Preferences.jl recompilation). Has no effect when a
custom library is configured via [`set_library!`](@ref).
"""
function set_petscint!(::Type{T}) where {T<:Union{Int32,Int64}}
    @set_preferences!("PetscInt" => string(T))
    @info "PETSc_jll PetscInt width set — restart Julia to use it." PetscInt = T
    return nothing
end


"""
    check_wrappers_version(petsclib=nothing)

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
function check_wrappers_version(petsclib=nothing)
    verfile = joinpath(@__DIR__, "autowrapped", "petsc_wrappers_version.jl")

    if !isdefined(@__MODULE__, :PETSC_WRAPPERS_VERSION) && isfile(verfile)
        try
            Base.invokelatest(() -> include(verfile))
        catch err
            @warn "Failed to include petsc_wrappers_version.jl" exception=(err,)
        end
    end

    wrappers_version = Base.invokelatest(() -> isdefined(@__MODULE__, :PETSC_WRAPPERS_VERSION) ? getproperty(@__MODULE__, :PETSC_WRAPPERS_VERSION) : nothing)

    if petsclib === nothing
        if isdefined(@__MODULE__, :petsclibs) && !isempty(petsclibs)
            petsclib = petsclibs[1]
        else
            throw(
                ArgumentError(
                    "no PETSc library available to check the installed version",
                ),
            )
        end
    end

    if isa(petsclib, String)
        petsclib = LibPETSc.PetscLibType(petsclib)
    end

    installed_version = nothing
    try
        major, minor, subminor, _release = LibPETSc.PetscGetVersionNumber(petsclib)
        installed_version = VersionNumber(Int(major), Int(minor), Int(subminor))
    catch err
        @warn "Failed to query installed PETSc version" exception=(err,)
    end

    match = isnothing(wrappers_version) || isnothing(installed_version) ? nothing : (installed_version.major == wrappers_version.major && installed_version.minor == wrappers_version.minor)
    if !isnothing(match) && match === false
        @warn "PETSc wrappers version does not match PETSc version of library (major.minor); this can cause undesired behavior" wrappers_version=wrappers_version installed_version=installed_version
    end

    return (wrappers_version = wrappers_version, installed_version = installed_version, match = match)
end
