# Custom display for REPL
function Base.show(io::IO, v::AbstractPetscOptions{PetscLib}) where {PetscLib}
    if v.ptr == C_NULL
        println(io, "PETSc Options (null pointer)")
        return
    end
    
    # Try to display options, but handle errors gracefully
    try
        println(io, "PETSc Options database:")
        LibPETSc.PetscOptionsView(PetscLib, v, C_NULL)  # NULL viewer: PETSc prints to stdout
    catch
        println(io, "PETSc Options (not yet initialized)")
    end
    return nothing
end


"""
    PetscOptions(petsclib; kwargs...)

Create a PETSc options database for the given `petsclib`.

Keyword arguments are converted to PETSc options:
- Options with value `nothing` or `true` are set without a value (flags)
- Options with value `false` are left out, and `opts[key] = false` clears an
  option set earlier
- Other values are converted to strings

# Examples
```julia-repl
julia> using PETSc

julia> petsclib = PETSc.petsclibs[1];

julia> PETSc.initialize(petsclib)

julia> opt = PETSc.PetscOptions(
                         petsclib,
                         ksp_monitor = nothing,
                         ksp_view = true,
                         pc_type = "mg",
                         pc_mg_levels = 1,
                         false_opt = false,
                     )
#PETSc Option Table entries:
-ksp_monitor
-ksp_view
-pc_mg_levels 1
-pc_type mg
#End of PETSc Option Table entries


julia> opt["ksp_monitor"]
""

julia> opt["pc_type"]
"mg"

julia> opt["pc_type"] = "ilu"
"ilu"

julia> opt["pc_type"]
"ilu"

julia> opt["false_opt"]
ERROR: KeyError: key "false_opt" not found

julia> opt["ksp_view"] = false
false

julia> opt["ksp_view"]
ERROR: KeyError: key "ksp_view" not found

julia> opt["bad_key"]
ERROR: KeyError: key "bad_key" not found
```

# External Links
$(doc_external("Sys/PetscOptionsCreate"))
"""
function LibPETSc.PetscOptions(petsclib::PetscLibType; kwargs...)
    check_initialized(petsclib)
    opts = LibPETSc.PetscOptionsCreate(petsclib)
    destroy_at_gc!(opts)
    for (k, v) in kwargs
        opts[k] = v
    end

    return opts
end

"""
    destroy!(opts::AbstractPetscOptions)

Free the options database `opts` holds, if this process is still allowed to.

Does nothing when the library has been finalized or re-initialized, or when
`opts` was already destroyed: see [`isdestroyable`](@ref).
Does nothing on a borrowed handle either: see [`owns`](@ref).
"""
destroy!(opts::AbstractPetscOptions{PetscLib}) where {PetscLib} =
    destroy_handle!(LibPETSc.PetscOptionsDestroy, opts, PetscLib)

function Base.setindex!(
    opts::AbstractPetscOptions{PetscLib},
    val,
    key,
) where {PetscLib}
    val === true && (val = nothing)
    # `-flag false` would still read as on through PetscOptionsHasName, so clear it
    if val === false
        LibPETSc.PetscOptionsClearValue(PetscLib, opts, string('-', key))
        return opts
    end

    LibPETSc.PetscOptionsSetValue(
        PetscLib,
        opts,
        string('-', key),
        isnothing(val) ? C_NULL : string(val),
    )

    return opts
end

function Base.getindex(opts::AbstractPetscOptions{PetscLib}, key) where {PetscLib}
    val = LibPETSc.PetscOptionsGetString(PetscLib, opts, C_NULL, string('-', key))
    if val == false
        throw(KeyError(key))
    end
    return val
end

"""
    parse_options(args::Vector{String})

Parse the `args` vector into a `NamedTuple` that can be used as the options for
the PETSc solvers.

```sh
julia --project file.jl -ksp_monitor -pc_type mg -ksp_view -da_refine=1
```
"""
function parse_options(args::Vector{String})
    i = 1
    opts = Dict{Symbol, Union{String, Nothing}}()
    while i <= length(args)
        (length(args[i]) > 1 && args[i][1] == '-') || throw(
            ArgumentError(
                "expected an option starting with '-', got $(repr(args[i]))",
            ),
        )
        if i == length(args) || args[i + 1][1] == '-'
            token = split(args[i][2:end], "=")
            if length(token) == 1
                opts[Symbol(token[1])] = nothing
            elseif length(token) == 2
                opts[Symbol(token[1])] = token[2]
            else
                throw(ArgumentError("invalid argument: $(repr(args[i]))"))
            end
            i = i + 1
        else
            opts[Symbol(args[i][2:end])] = args[i + 1]
            i = i + 2
        end
    end
    return NamedTuple(opts)
end

# PETSc's options stack is a linked list through the pushed objects, so pushing
# one that is already on it breaks the chain: the next pops fail and the global
# database is lost. These mirror the stack per library to refuse that.
const options_stacks = IdDict{Any, Vector{Ptr{Cvoid}}}()
options_stack(::Type{PetscLib}) where {PetscLib} =
    get!(() -> Ptr{Cvoid}[], options_stacks, PetscLib)

function Base.push!(opts::AbstractPetscOptions{PetscLib}) where {PetscLib}
    stack = options_stack(PetscLib)
    opts.ptr in stack && throw(ArgumentError(
        "these options are already pushed; PETSc's options stack can hold an " *
        "options database only once"))
    LibPETSc.PetscOptionsPush(PetscLib, opts)
    push!(stack, opts.ptr)
    return nothing
end

"""
    pop!(opts::PetscOptions)

Pop `opts` off PETSc's options stack, where `push!(opts)` put it. An
`opts` that `push!` put below other databases throws an `ArgumentError` and pops
nothing: pop those first. A database pushed through `LibPETSc.PetscOptionsPush`
is popped as PETSc pops it, since only PETSc knows its place on the stack.
Prefer [`with_options`](@ref), which pairs the two.

# External Links
$(doc_external("Sys/PetscOptionsPop"))
"""
function Base.pop!(opts::AbstractPetscOptions{PetscLib}) where {PetscLib}
    stack = options_stack(PetscLib)
    i = findlast(==(opts.ptr), stack)
    isnothing(i) || i == lastindex(stack) || throw(ArgumentError(
        "these options are not on top of PETSc's options stack: pop the " *
        "$(lastindex(stack) - i) pushed after them first"))
    LibPETSc.PetscOptionsPop(PetscLib)
    isnothing(i) || pop!(stack)
    return nothing
end

# The constructor options of a KSP, SNES or TS go with the solver
function destroy_options!(solver)
    isnothing(solver.opts) || (solver.opts = nothing)
    return nothing
end

# Keys this package set in PETSc's options database for a solver's options block,
# with their values, so that a nested block overrides them instead of taking them
# for the user's.
const layered_keys = IdDict{Any, Dict{String, Union{Nothing, String}}}()

# Whether the default database holds `key` with the value an enclosing block set.
# Checked by value, as PETSc has no call that names the database on top of its
# stack: under a database pushed through `LibPETSc` the key is then the user's.
function set_by_enclosing_block(petsclib, db, ours, key)
    haskey(ours, key) || return false
    # shortcut: PetscOptionsGetString reads at most 255 bytes, so a longer value
    # counts as the user's; read the length first if such values turn up.
    current = LibPETSc.PetscOptionsGetString(petsclib, db, "", key)
    return current == something(ours[key], "")
end

# Run `f()` with the `-key => value` pairs `entries` in PETSc's default options
# database, the top of its stack. A key the database already has wins, unless an
# enclosing block set it. Each key is put back as it was afterwards, even if `f`
# throws, so the database holds only what the user set.
function with_layered_options(f, ::Type{PetscLib}, entries) where {PetscLib}
    petsclib = getlib(PetscLib)
    db = LibPETSc.PetscOptions{PetscLib}(C_NULL, petsclib.age; own = false)
    ours = get!(Dict{String, Union{Nothing, String}}, layered_keys, PetscLib)
    undo = Tuple{String, Bool, Union{Missing, Nothing, String}}[]
    try
        for (key, val) in entries
            mine = set_by_enclosing_block(petsclib, db, ours, key)
            !mine && LibPETSc.PetscOptionsHasName(petsclib, db, "", key) && continue
            push!(undo, (key, mine, get(ours, key, missing)))
            LibPETSc.PetscOptionsSetValue(petsclib, db, key, isnothing(val) ? C_NULL : val)
            ours[key] = val
        end
        return f()
    finally
        for (key, mine, prev) in Iterators.reverse(undo)
            if mine
                LibPETSc.PetscOptionsSetValue(petsclib, db, key, isnothing(prev) ? C_NULL : prev)
            else
                LibPETSc.PetscOptionsClearValue(petsclib, db, key)
            end
            ismissing(prev) ? delete!(ours, key) : (ours[key] = prev)
        end
    end
end

# Constructor options as `-key => value` pairs, read as `setindex!` on a
# PetscOptions reads them: `true` and `nothing` are flags, `false` leaves it out.
option_entries(options) = Pair{String, Union{Nothing, String}}[
    string('-', k) => (v === true || isnothing(v) ? nothing : string(v))
    for (k, v) in pairs(options) if v !== false
]

"""
    with_options(f, opts::PetscOptions)

Run `f()` with `opts` pushed onto PETSc's options stack, so that PETSc reads its
options from `opts` alone, and pop it again even if `f` throws. Returns the result
of `f`. PETSc's stack can hold a database only once, so nesting `with_options` on
the same `opts` throws an `ArgumentError`.

# External Links
$(doc_external("Sys/PetscOptionsPush"))
"""
function with_options(f, opts::AbstractPetscOptions{PetscLib}) where {PetscLib}
    push!(opts)
    try
        return f()
    finally
        pop!(opts)
    end
end

"""
    parse_option(opt::NamedTuple, key::Symbol, default::T)

Parse `opt` similar to `Base.get` but ensures that the returned value is the
same type as the default value. When `T <: NTuple` keys that result in a single
value will be filled into an `NTuple` of the same length as `T`; in the case of
strings it is parsed using `Base.split` with comma delimiter

# Examples
```julia-repl
julia> opt = (tup = (1, 2, 3), string_tup = "1,2,3", string_int = "4", int = 4)
(tup = (1, 2, 3), string_tup = "1,2,3", string_int = "4", int = 4)

julia> parse_option(opt, :int, 7)
4

julia> parse_option(opt, :bad_key, 7)
7

julia> parse_option(opt, :tup, (1, 1, 1))
(1, 2, 3)

julia> parse_option(opt, :string_tup, (1, 1, 1))
tokens = SubString{String}["1", "2", "3"]
(1, 2, 3)

julia> parse_option(opt, :string_int, (1, 1, 1))
tokens = SubString{String}["4"]
(4, 4, 4)

julia> parse_option(opt, :int, (1, 1, 1))
(4, 4, 4)

julia> parse_option(opt, :int, (1., 1., 1.))
(4.0, 4.0, 4.0)
```
"""
function parse_option(opt::NamedTuple, key::Symbol, default::T) where {T}
    v = get(opt, key, default)
    if !(v isa T)
        if T <: String
            return string(v)
        elseif v isa String
            if T <: NTuple
                ET = T.types[1]
                tokens = split(v, ",")
                if length(tokens) == 1
                    return ntuple(_ -> parse(ET, tokens[1]), length(T.types))
                else
                    return ntuple(j -> parse(ET, tokens[j]), length(T.types))
                end
            else
                return parse(T, v)
            end
        else
            if T <: NTuple && !(v isa Tuple)
                ET = T.types[1]
                return ntuple(j -> convert(ET, v), length(T.types))
            else
                return convert(T, v)
            end
        end
    end
    return v
end


# ============================================================================
#   Options on solvers and preconditioners
# ============================================================================

"""
    set_from_options!(obj)

Apply PETSc's options to the `KSP`, `SNES` or `TS` `obj`, read under its prefix
(see [`set_options_prefix!`](@ref)): PETSc's options database
(`PETSc.initialize(; options)`, the `PETSC_OPTIONS` environment variable), and the
options given to its constructor for keys the database does not set. Returns `obj`.

The constructor's `options_at` keyword decides when the solver does this itself;
call it to configure `obj` at another point, for example a nested solver before
its first use. With `options_at = :first_solve`, a call here counts as that first
application.

# External Links
$(doc_external("KSP/KSPSetFromOptions"))
$(doc_external("SNES/SNESSetFromOptions"))
$(doc_external("TS/TSSetFromOptions"))
"""
function set_from_options!(obj::Union{LibPETSc.AbstractKSP, LibPETSc.AbstractSNES, LibPETSc.AbstractTS})
    with_options(() -> apply_options!(obj), obj)
    return obj
end

for (T, setfromoptions) in (
    (:AbstractKSP, :KSPSetFromOptions),
    (:AbstractSNES, :SNESSetFromOptions),
    (:AbstractTS, :TSSetFromOptions),
)
    @eval function apply_options!(obj::LibPETSc.$T{PetscLib}) where {PetscLib}
        LibPETSc.$setfromoptions(getlib(PetscLib), obj)
        obj.options_applied = true
        return obj
    end
end

"""
    with_options(f, obj)

Run `f()` with the options given to the constructor of `obj`, a `KSP`, `SNES` or
`TS`, in PETSc's options database, and return the result of `f`. A key the
database already sets keeps its value there, so `PETSC_OPTIONS` and
`PETSc.initialize(; options)` override the constructor. Afterwards the database
is as it was.

PETSc reads some options only when it creates a sub-object, such as the solvers
of a fieldsplit preconditioner during setup. A call that triggers that outside
`solve!`, for example `KSPSetUp` through `LibPETSc`, sees the constructor options
only inside this block.
"""
function with_options(f, obj::Union{LibPETSc.AbstractKSP{PetscLib}, LibPETSc.AbstractSNES{PetscLib}, LibPETSc.AbstractTS{PetscLib}}) where {PetscLib}
    opts = obj.opts
    isnothing(opts) && return f()
    return with_layered_options(f, PetscLib, opts)
end

# When the constructor options of a KSP, SNES or TS are applied
# ----------------------------------------------------------------------------

const OPTIONS_AT = (:construction, :each_solve, :first_solve)

function check_options_at(options_at::Symbol)
    options_at in OPTIONS_AT || throw(ArgumentError(
        "options_at must be :construction, :each_solve or :first_solve, got :$options_at"))
    return options_at
end

# Whether `XSetFromOptions` runs now. `solving` is true in `solve!` and false in
# `setup!` and `step!`, which apply the options only under `:first_solve`. A
# solver wrapped from a handle has `options_at = :never`: its owner applies them.
function options_due(obj, solving::Bool)
    obj.options_at === :each_solve  && return solving
    obj.options_at === :first_solve && return !obj.options_applied
    return false
end

# Run `f()` with the constructor options in effect, applying them first when due.
# Under `:each_solve` the options belong to `solve!` alone: `setup!` and `step!`
# neither set nor apply them.
function with_solver_options(f, obj; solving::Bool)
    !solving && obj.options_at === :each_solve && return f()
    return with_options(obj) do
        options_due(obj, solving) && apply_options!(obj)
        f()
    end
end

# Apply the options of a freshly built solver if `options_at = :construction`,
# and record the choice for later solves.
function store_options!(obj, options_at::Symbol, options)
    obj.options_at = check_options_at(options_at)
    isempty(options) || (obj.opts = option_entries(options))
    options_at === :construction && set_from_options!(obj)
    return obj
end

"""
    set_from_options!(obj::Union{AbstractPC, AbstractPetscMat, AbstractPetscVec})

Apply PETSc's global options database to the preconditioner, matrix or vector
`obj`, read under its prefix (see [`set_options_prefix!`](@ref)): `-pc_type`,
`-mat_type`, `-vec_type` and the options of the chosen type. Returns `obj`.

A matrix or vector reads its type here, so call it before the object is
preallocated or filled.

# External Links
$(doc_external("PC/PCSetFromOptions"))
$(doc_external("Mat/MatSetFromOptions"))
$(doc_external("Vec/VecSetFromOptions"))
"""
function set_from_options!(obj::Union{LibPETSc.AbstractPC, LibPETSc.AbstractPetscMat, LibPETSc.AbstractPetscVec})
    apply_global_options!(obj)
    return obj
end

apply_global_options!(p::LibPETSc.AbstractPC{PetscLib}) where {PetscLib} =
    LibPETSc.PCSetFromOptions(getlib(PetscLib), p)
apply_global_options!(A::LibPETSc.AbstractPetscMat{PetscLib}) where {PetscLib} =
    LibPETSc.MatSetFromOptions(getlib(PetscLib), A)
apply_global_options!(v::LibPETSc.AbstractPetscVec{PetscLib}) where {PetscLib} =
    LibPETSc.VecSetFromOptions(getlib(PetscLib), v)

"""
    set_options_prefix!(obj, prefix::AbstractString)

Make `obj`, a `KSP`, `SNES`, `TS`, `PC`, DM, matrix or vector, read its options
under `prefix`: with `"inner_"`, a `KSP` reads `-inner_ksp_type`. This is how two
objects of one kind in one program are configured apart. Returns `obj`.

# External Links
$(doc_external("KSP/KSPSetOptionsPrefix"))
$(doc_external("SNES/SNESSetOptionsPrefix"))
$(doc_external("TS/TSSetOptionsPrefix"))
$(doc_external("PC/PCSetOptionsPrefix"))
$(doc_external("DM/DMSetOptionsPrefix"))
$(doc_external("Mat/MatSetOptionsPrefix"))
$(doc_external("Vec/VecSetOptionsPrefix"))
"""
function set_options_prefix! end

"""
    options_prefix(obj)

The prefix `obj`, a `KSP`, `SNES`, `TS`, `PC`, DM, matrix or vector, reads its
options under, or `""` when it has none. A solver PETSc builds inside another
inherits the outer prefix: the `KSP` of a `SNES` prefixed `"outer_"` reads
`-outer_ksp_type`.

# External Links
$(doc_external("KSP/KSPGetOptionsPrefix"))
$(doc_external("SNES/SNESGetOptionsPrefix"))
$(doc_external("TS/TSGetOptionsPrefix"))
$(doc_external("PC/PCGetOptionsPrefix"))
$(doc_external("DM/DMGetOptionsPrefix"))
$(doc_external("Mat/MatGetOptionsPrefix"))
$(doc_external("Vec/VecGetOptionsPrefix"))
"""
function options_prefix end

for (T, setprefix, getprefix) in (
    (:AbstractKSP, :KSPSetOptionsPrefix, :KSPGetOptionsPrefix),
    (:AbstractSNES, :SNESSetOptionsPrefix, :SNESGetOptionsPrefix),
    (:AbstractTS, :TSSetOptionsPrefix, :TSGetOptionsPrefix),
    (:AbstractPC, :PCSetOptionsPrefix, :PCGetOptionsPrefix),
    (:AbstractPetscDM, :DMSetOptionsPrefix, :DMGetOptionsPrefix),
    (:AbstractPetscMat, :MatSetOptionsPrefix, :MatGetOptionsPrefix),
    (:AbstractPetscVec, :VecSetOptionsPrefix, :VecGetOptionsPrefix),
)
    @eval begin
        function set_options_prefix!(obj::LibPETSc.$T{PetscLib}, prefix::AbstractString) where {PetscLib}
            LibPETSc.$setprefix(getlib(PetscLib), obj, String(prefix))
            return obj
        end
        options_prefix(obj::LibPETSc.$T{PetscLib}) where {PetscLib} =
            LibPETSc.$getprefix(getlib(PetscLib), obj)
    end
end
