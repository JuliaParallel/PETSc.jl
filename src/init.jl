"""
    PetscNotInitialized(petsclib)

Thrown when a PETSc object is built or used before its library is initialized.

Call [`initialize`](@ref) on the library first. Every PETSc call needs the
library up, so this is reported as its own type rather than as an
`ArgumentError`, which lets callers catch it specifically.
"""
struct PetscNotInitialized <: Exception
    petsclib::Any
end

function Base.showerror(io::IO, e::PetscNotInitialized)
    print(io, "PetscNotInitialized: the PETSc library ")
    # Name the library by its scalar and integer types rather than by printing
    # the whole handle, which carries an artifact path and drowns the message.
    try
        print(io, "(PetscScalar = ", scalartype(e.petsclib))
        print(io, ", PetscInt = ", inttype(e.petsclib), ") ")
    catch
        print(io, e.petsclib, " ")
    end
    return print(io, "is not initialized. Call PETSc.initialize(petsclib) first.")
end

"""
    check_initialized(petsclib)

Throw [`PetscNotInitialized`](@ref) unless `petsclib` is initialized.
"""
function check_initialized(petsclib)
    isinitialized(petsclib) || throw(PetscNotInitialized(petsclib))
    return nothing
end

"""
   isinitialized(petsclib)

Check if `petsclib` is initialized

# External Links
$(doc_external("Sys/PetscInitialized"))
"""
isinitialized(petsclib) = LibPETSc.PetscInitialized(petsclib)

const petsc_program_name = "petsc_julia"
const lib_handles = IdDict{Any, Tuple{Any, Bool}}()

"""
    initialize([petsclib]; log_view = false, options = String[])

Initialize the `petsclib`. If no `petsclib` is given, all `PETSc.petsclibs`
will be initialized.

Additionally:
 - This will initialize MPI if it has not already been initialized.
 - It will disable the PETSc signal handler (via
   $(petsc_link("Sys/PetscPopSignalHandler"))
 - Add an [`atexit`](https://docs.julialang.org/en/v1/base/base/#Base.atexit)
   hook to call [`PETSc.finalize`](@ref).

# Arguments
- `log_view::Bool = false`: Enable PETSc's `-log_view` performance logging.
  When enabled, PETSc will output performance statistics at finalization.
- `options::Vector{String} = String[]`: Additional PETSc command-line options.
  These are passed via the `PETSC_OPTIONS` environment variable.

# BLAS threads
PETSc's BLAS calls run in Julia's OpenBLAS thread pool, the one
`LinearAlgebra.BLAS.set_num_threads` sizes, for every library in `petsclibs`.
`-blas_num_threads n`, given in `options`, in `PETSC_OPTIONS` or on the command
line, sets that pool to `n` threads. Under MPI, several ranks on one node each
own such a pool, and its busy-waiting threads compete with the ranks for cores:
a solve can run many times slower with the right answer. Pass
`-blas_num_threads 1` there, or set `OPENBLAS_NUM_THREADS=1`.

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

# One BLAS thread per rank under MPI
PETSc.initialize(petsclib; options = ["-blas_num_threads", "1"])
```

# External Links
$(doc_external("Sys/PetscInitializeNoArguments"))
"""
function initialize(; log_view::Bool = false, options = String[])
    map(petsclib -> initialize(petsclib; log_view, options), petsclibs)
    return nothing
end

function initialize(petsclib; log_view::Bool = false, options = String[])
    if !isinitialized(petsclib)
        # PETSc's signal handler conflicts with Julia's own handlers when using multithreading
        # Appended to a copy: the caller's vector stays as it was.
        options = [String.(options); "-no_signal_handler"]
        cli_opts = build_petsc_options(log_view, options)
        prev_opts = get(ENV, "PETSC_OPTIONS", "")
        ENV["PETSC_OPTIONS"] = isempty(prev_opts) ? cli_opts : "$prev_opts $cli_opts"
        try
            ensure_mpi_initialized()
            petsclib.age += 1
            LibPETSc.PetscInitializeNoArguments(petsclib)
            post_initialize(petsclib)
        finally
            if isempty(prev_opts)
                delete!(ENV, "PETSC_OPTIONS")
            else
                ENV["PETSC_OPTIONS"] = prev_opts
            end
        end
    end
    return nothing
end

"""
   finalize(petsclib)

Finalize the `petsclib`, if no `petsclib` is given then all `PETSc.petsclibs`
will be finalized.

# External Links
$(doc_external("Sys/PetscFinalize"))
"""
function finalize()
    map(finalize, petsclibs)
    return nothing
end

function finalize(petsclib)
    if !isfinalized(petsclib)
        petsclib.age += 1
        LibPETSc.PetscFinalize(petsclib)
        drop_object_states!(typeof(petsclib))
    end
    return nothing
end

"""
   isfinalized(petsclib)

Check if `petsclib` is finalized

# External Links
$(doc_external("Sys/PetscFinalized"))
"""
isfinalized(petsclib) = LibPETSc.PetscFinalized(petsclib)

"""
    isdestroyable(obj, ::Type{PetscLib})

Whether `obj` still refers to a PETSc object this process is allowed to destroy.

It is not destroyable when the library is finalized, when the pointer is already
null, or when the object predates the current initialize/finalize cycle.
`initialize` and `finalize` both bump `petsclib.age`, and every object records
the age it was created under. `PetscFinalize` frees everything it owns,
including the inner communicator, so calling `xxxDestroy` on an object from an
earlier cycle reaches a communicator that no longer exists and aborts inside
MPI with "Invalid communicator". That happens from a GC finalizer, so it
surfaces at an arbitrary later point rather than where the object was dropped.
"""
function isdestroyable(obj, ::Type{PetscLib}) where {PetscLib}
    isfinalized(PetscLib) && return false
    obj.ptr == C_NULL && return false
    return obj.age == getlib(PetscLib).age
end

"""
    owns(obj)

Whether `destroy!` on `obj` destroys the PETSc object it holds, which is the
wrapper's `own` field.

A constructor returns an owning wrapper. A reader returns a borrowed one: the
object belongs to whatever it was read from, and `destroy!` on the borrowed
wrapper does nothing and leaves it usable. That covers the high-level readers
(`pc(ksp)`, `snes(ts)`, `dm(ksp)`, `solution(ksp)` and the rest), every
`LibPETSc` function with `Get` in its name except the few that hand out a new
reference (`MatGetFactor`, `DMLabelGetStratumIS`, ...: their manual pages say
the caller destroys the result), and the objects a callback receives.
"""
owns(obj) = obj.own

function build_petsc_options(log_view::Bool, options)
    opts = String[]
    if log_view
        push!(opts, "-log_view")
    end
    append!(opts, [String(opt) for opt in options])
    return join(opts, " ")
end

function ensure_mpi_initialized()
    MPI.Initialized() || MPI.Init()
    return nothing
end

function post_initialize(petsclib)
    # disable signal handler
    LibPETSc.PetscPopSignalHandler(petsclib)
    _reset_stale_register_flags(petsclib)
    atexit(() -> finalize(petsclib))
    apply_blas_num_threads(petsclib)
    return nothing
end

# Every PETSc_jll build calls BLAS through libblastrampoline, so PETSc's BLAS runs
# in Julia's own OpenBLAS thread pool. PETSc cannot size that pool: it was built
# against a generic BLAS, so `PetscBLASSetNumThreads` only records the number.
# `-blas_num_threads` is therefore forwarded here. The pool is process-wide and
# is left as it is at `finalize`: it serves all Julia code, not PETSc objects.
function apply_blas_num_threads(petsclib)
    global_options = LibPETSc.PetscOptions{typeof(petsclib)}(C_NULL, petsclib.age; own = false)
    n, set = LibPETSc.PetscOptionsGetInt(petsclib, global_options, "", "-blas_num_threads")
    Bool(set) || return nothing
    n >= 1 || throw(ArgumentError("-blas_num_threads must be at least 1, got $n"))
    LinearAlgebra.BLAS.set_num_threads(Int(n))
    return nothing
end

# PETSc 3.25.x: `TaoFinalizePackage` destroys the `TaoTerm` type list but never resets
# `TaoTermRegisterAllCalled`, so after a finalize/initialize cycle `TaoCreate` fails with
# "Unable to find requested TaoTerm type callbacks". Reset the flag so the list is rebuilt.
# The flag is an internal symbol: reachable on Linux/macOS (ELF/Mach-O export everything),
# not on Windows, where Tao therefore only works in the first initialize/finalize cycle.
const _taoterm_resettable = Ref{Union{Nothing,Bool}}(nothing)

"""
    tao_usable_after_reinitialize()

Whether `Tao` objects can be created after `finalize` followed by `initialize` with the current
PETSc binaries (false on Windows with PETSc 3.25.x, see `_reset_stale_register_flags`).
"""
tao_usable_after_reinitialize() = _taoterm_resettable[] !== false

function _reset_stale_register_flags(petsclib)
    handle, _ = ensure_library_handle(petsclib)
    lib = library_ptr(handle)
    # PETSc 3.25.x: these packages destroy their type lists at PetscFinalize without resetting
    # the RegisterAll flag (TaoFinalizePackage for Tao and TaoTerm; TSTrajectory likewise).
    # TSIRKRegisterAllCalled is `static`, but the public TSIRKRegisterDestroy resets it.
    # Two more are `static` with no reset at all, so they cannot be fixed from here and their
    # types work only in the first initialize/finalize cycle of a process:
    # KSPMatRegisterAllCalled (the LMVM matrix types behind the `lmvm` Tao types) and
    # KSPGuessRegisterAllCalled (`-ksp_guess_type`). All six belong in an upstream fix.
    LibPETSc.TSIRKRegisterDestroy(petsclib)
    ok = true
    for sym in (:TaoRegisterAllCalled, :TaoTermRegisterAllCalled, :TSTrajectoryRegisterAllCalled)
        p = Libdl.dlsym_e(lib, sym)
        if p == C_NULL
            ok = false
        else
            # the flags are PetscBool (one byte): a wider store overwrites the neighbouring globals
            unsafe_store!(Ptr{LibPETSc.PetscBool}(p), LibPETSc.PETSC_FALSE)
        end
    end
    if !ok && _taoterm_resettable[] === nothing
        # see tao_usable_after_reinitialize(); a warning here would repeat in every process
        @debug "PETSc 3.25.x loses its Tao/TaoTerm types at PetscFinalize and the workaround cannot be applied " *
               "with these binaries (internal symbols not exported): Tao objects can only be created before the first finalize"
    end
    _taoterm_resettable[] = ok
    return nothing
end

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
    LibraryInfo

The shape of the `NamedTuple` [`library_info`](@ref) returns, and the type its
`show` method is written for.
"""
const LibraryInfo = NamedTuple{(:source, :path, :scalar, :int, :real)}

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

v0.4 printed a report and returned `nothing`, so the name promised data it never
handed back. The report is unchanged — it is now the `show` method — and the
values are reachable from code.

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

# The report v0.4's `library_info` printed. It is written for the exact field
# names `library_info` returns, so it cannot claim a `NamedTuple` this package
# does not own.
function Base.show(io::IO, ::MIME"text/plain", info::LibraryInfo)
    if info.source === :preferences
        println(io, "Source  : LocalPreferences.toml")
        println(io, "Path    : ", info.path)
        println(io, "Scalar  : ", info.scalar)
        println(io, "Int     : ", info.int)
    else
        println(io, "Source  : PETSc_jll (default precompiled binaries)")
    end

    println(io, "\nLoaded libraries (this session):")
    for lib in petsclibs
        println(
            io,
            "  [$(lib.PetscScalar), $(lib.PetscInt)]: ",
            library_path_string(lib),
        )
    end
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

# ─────────────────────────────────────────────────────────────────────────────
# Type names (docs/src/man/naming.md §3.1)
#
# PETSc registers type names at runtime as strings, so they cannot become enums.
# They are `Symbol` at the Julia API and `String` at the C boundary, and the
# conversion happens here, once, where the call meets C.
#
# The generated `*GetType` wrappers answer `""` (and `MatGetType` answers
# `"(not set)"`) when PETSc has no type for the object yet. A reader that can
# decline to answer is documented rather than papered over with a default
# (§15.1), so both spellings come back as `nothing`.
# ─────────────────────────────────────────────────────────────────────────────

type_name_symbol(s::AbstractString) =
    (isempty(s) || s == "(not set)") ? nothing : Symbol(s)

"""
    set_type!(obj, type::Symbol)

Set the PETSc implementation `obj` uses, for example `set_type!(ksp, :gmres)`.

Defined for `PetscVec`, `PetscMat`, `KSP`, `PC`, `SNES`, `TS` and `AbstractPetscDM`.
The `Symbol` is converted to a `String` at the C boundary (§3.1).

`set_type!(obj, "gmres")` still works in v0.5 and warns; it is a `MethodError`
in v0.6.
"""
function set_type! end
