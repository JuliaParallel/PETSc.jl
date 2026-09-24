

import .LibPETSc: AbstractPetscOptions, PetscOptions, COptions

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

Replaces v0.4's `Options` factory: construction goes through the type
(docs/src/man/naming.md §5.1).

Keyword arguments are converted to PETSc options:
- Options with value `nothing` or `true` are set without a value (flags)
- Options with value `false` are not set
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
ERROR: KeyError: key "bad_key" not found

julia> opt["bad_key"]
ERROR: KeyError: key "bad_key" not found
```

# External Links
$(doc_external("Sys/PetscOptionsCreate"))
"""
function LibPETSc.PetscOptions(petsclib::PetscLibType; kwargs...)
    check_initialized(petsclib)
    opts = LibPETSc.PetscOptionsCreate(petsclib)
    finalizer(destroy!, opts)
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
function destroy!(opts::AbstractPetscOptions{PetscLib}) where {PetscLib}
    owns(opts) || return nothing
    if isdestroyable(opts, PetscLib)
        LibPETSc.PetscOptionsDestroy(PetscLib, opts)
    end
    opts.ptr = C_NULL
    return nothing
end

function Base.setindex!(
    opts::AbstractPetscOptions{PetscLib},
    val,
    key,
) where {PetscLib}
    val === true && (val = nothing)
    val === false && (return opts)

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

function Base.push!(opts::AbstractPetscOptions{PetscLib}) where {PetscLib}
    LibPETSc.PetscOptionsPush(PetscLib, opts)
    return nothing
end

function Base.pop!(opts::AbstractPetscOptions{PetscLib}) where {PetscLib}
    LibPETSc.PetscOptionsPop(PetscLib)
    return nothing
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

Apply PETSc's options to the `KSP`, `SNES` or `TS` `obj`: those given to its
constructor, then the global database, so a command-line option given for the
object's prefix (see [`set_options_prefix!`](@ref)) is read too. Returns `obj`.

`solve!` does this itself; call it to configure `obj` earlier, for example a
nested solver before its first use.

# External Links
$(doc_external("KSP/KSPSetFromOptions"))
$(doc_external("SNES/SNESSetFromOptions"))
$(doc_external("TS/TSSetFromOptions"))
"""
function set_from_options! end

for (T, setfromoptions) in (
    (:AbstractKSP, :KSPSetFromOptions),
    (:AbstractSNES, :SNESSetFromOptions),
    (:AbstractTS, :TSSetFromOptions),
)
    @eval function set_from_options!(obj::LibPETSc.$T{PetscLib}) where {PetscLib}
        opts = obj.opts
        isnothing(opts) || push!(opts)
        try
            LibPETSc.$setfromoptions(getlib(PetscLib), obj)
        finally
            isnothing(opts) || pop!(opts)
        end
        return obj
    end
end

"""
    set_options_prefix!(obj, prefix::AbstractString)

Make `obj`, a `KSP`, `SNES`, `TS` or `PC`, read its options under `prefix`:
with `"inner_"`, a `KSP` reads `-inner_ksp_type`. This is how two solvers in one
program are configured apart. Returns `obj`.

# External Links
$(doc_external("KSP/KSPSetOptionsPrefix"))
$(doc_external("SNES/SNESSetOptionsPrefix"))
$(doc_external("TS/TSSetOptionsPrefix"))
$(doc_external("PC/PCSetOptionsPrefix"))
"""
function set_options_prefix! end

"""
    options_prefix(obj)

The prefix `obj`, a `KSP`, `SNES`, `TS` or `PC`, reads its options under, or
`""` when it has none. A solver PETSc builds inside another inherits the outer
prefix: the `KSP` of a `SNES` prefixed `"outer_"` reads `-outer_ksp_type`.

# External Links
$(doc_external("KSP/KSPGetOptionsPrefix"))
$(doc_external("SNES/SNESGetOptionsPrefix"))
$(doc_external("TS/TSGetOptionsPrefix"))
$(doc_external("PC/PCGetOptionsPrefix"))
"""
function options_prefix end

for (T, setprefix, getprefix) in (
    (:AbstractKSP, :KSPSetOptionsPrefix, :KSPGetOptionsPrefix),
    (:AbstractSNES, :SNESSetOptionsPrefix, :SNESGetOptionsPrefix),
    (:AbstractTS, :TSSetOptionsPrefix, :TSGetOptionsPrefix),
    (:AbstractPC, :PCSetOptionsPrefix, :PCGetOptionsPrefix),
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
