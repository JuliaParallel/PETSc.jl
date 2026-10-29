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
- `options::Vector{String} = String[]`: Additional PETSc command-line options,
  handed to `PetscInitialize` as its command line. They override the
  `PETSC_OPTIONS` environment variable.

# BLAS threads
PETSc's BLAS calls run in Julia's OpenBLAS thread pool, the one
`LinearAlgebra.BLAS.set_num_threads` sizes, for every library in `petsclibs`.
`-blas_num_threads n`, given in `options` or in `PETSC_OPTIONS`, sets that pool
to `n` threads.

Without that option, when several MPI ranks share a node, `initialize` sets the
pool to one thread. Each rank owns a pool, and its busy-waiting threads would
otherwise compete with the ranks for cores: a solve can run many times slower
with the right answer. A set `OPENBLAS_NUM_THREADS` or `OMP_NUM_THREADS` also
leaves the pool as it is. Serial runs, and runs with one rank per node, are
unaffected. This changes Julia's process-wide BLAS setting, which `finalize`
does not restore. `initialize` sets the pool only once, so a later
`BLAS.set_num_threads(n)` takes effect as usual.

# Examples
```julia
# Basic initialization
PETSc.initialize(petsclib)

# Enable performance logging to stdout
PETSc.initialize(petsclib; log_view = true)

# Write log to a file. PETSc reads `:file:format`, so on Windows give a path
# without a drive letter, such as a relative one
PETSc.initialize(petsclib; log_view = true, options = [":logfile.txt"])

# Enable memory logging
PETSc.initialize(petsclib; log_view = true, options = [":logfile.txt", "-log_view_memory"])

# Pass custom PETSc options without logging
PETSc.initialize(petsclib; options = ["-malloc_debug", "-on_error_abort"])

# Four BLAS threads per rank, overriding the one-thread default for ranks sharing a node
PETSc.initialize(petsclib; options = ["-blas_num_threads", "4"])
```

# External Links
$(doc_external("Sys/PetscInitialize"))
"""
function initialize(; log_view::Bool = false, options = String[])
    map(petsclib -> initialize(petsclib; log_view, options), petsclibs)
    return nothing
end

function initialize(petsclib; log_view::Bool = false, options = String[])
    if !isinitialized(petsclib)
        ensure_mpi_initialized()
        petsclib.age += 1
        LibPETSc.PetscInitialize(petsclib, build_petsc_options(log_view, options))
        post_initialize(petsclib)
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
        empty!(options_stack(typeof(petsclib)))
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

# The command line PETSc is initialized with: the options go in as arguments,
# not through `PETSC_OPTIONS`, which PETSc on Windows cannot see once set from
# Julia. An entry holding several words, such as "-ksp_type cg", is split.
# PETSc's signal handler conflicts with Julia's own, notably under multithreading.
function build_petsc_options(log_view::Bool, options)
    args = [petsc_program_name]
    log_view && push!(args, "-log_view")
    for opt in options
        append!(args, split(String(opt)))
    end
    push!(args, "-no_signal_handler")
    return args
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
# `-blas_num_threads` is therefore forwarded here. Without it, ranks sharing a
# node get one thread each, unless the environment already sizes the pool. The
# pool is process-wide and is left as it is at `finalize`: it serves all Julia
# code, not PETSc objects.
function apply_blas_num_threads(petsclib)
    global_options = LibPETSc.PetscOptions{typeof(petsclib)}(C_NULL, petsclib.age; own = false)
    n, set = LibPETSc.PetscOptionsGetInt(petsclib, global_options, "", "-blas_num_threads")
    if Bool(set)
        n >= 1 || throw(ArgumentError("-blas_num_threads must be at least 1, got $n"))
        LinearAlgebra.BLAS.set_num_threads(Int(n))
    elseif !any(v -> haskey(ENV, v), BLAS_THREAD_VARIABLES) && ranks_on_node() > 1
        LinearAlgebra.BLAS.set_num_threads(1)
        LibPETSc.PetscBLASSetNumThreads(petsclib, petsclib.PetscInt(1))
    end
    return nothing
end

# The variables OpenBLAS and PETSc read for a thread count
const BLAS_THREAD_VARIABLES = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS")

# How many ranks of MPI.COMM_WORLD share this rank's node. Collective.
function ranks_on_node()
    world = MPI.COMM_WORLD
    MPI.Comm_size(world) == 1 && return 1
    node = MPI.Comm_split_type(world, MPI.COMM_TYPE_SHARED, MPI.Comm_rank(world))
    n = MPI.Comm_size(node)
    MPI.free(node)
    return n
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
