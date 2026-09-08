# ============================================================================
#   PETSc object leak audit
# ============================================================================
#
# Matches object creations against destroy calls by walking the parsed syntax tree. 
# Matching on source text instead would miss multi-statement lines and
# would fire inside comments and string literals.

# Names that create an object the caller owns
# ----------------------------------------------------------------------------

const AUDIT_TYPE_CREATORS =
    Dict(:KSP => "KSP", :SNES => "SNES", :DMDA => "DM", :DMStag => "DM", :DMPlex => "DM")

const AUDIT_NAMED_CREATORS = Dict(
    :DMGlobalVec => "Vec",
    :DMLocalVec => "Vec",
    :DMGetCoordinateDM => "DM",
    :DMStagCreateCompatibleDMStag => "DM",
    :DMCreateMatrix => "Mat",
    :MatCreateVecs => "Vec",
    # high-level factories, whose names carry no Create/Duplicate marker
    :VecSeq => "Vec",
    :MatAIJ => "Mat",
    :MatShell => "Mat",
    :MatSeqAIJ => "Mat",
    :MatSeqDense => "Mat",
)

"""
    audit_creator(name::Symbol) -> Union{Nothing, String}

The kind of object `name` creates, or `nothing` if it creates none. Covers the
`Vec`/`Mat` families (`VecCreateSeq`, `MatDuplicate`, `MatSeqAIJWithArrays`, …)
by shape rather than by listing every member.
"""
function audit_creator(name::Symbol)
    haskey(AUDIT_TYPE_CREATORS, name) && return AUDIT_TYPE_CREATORS[name]
    haskey(AUDIT_NAMED_CREATORS, name) && return AUDIT_NAMED_CREATORS[name]
    s = String(name)
    for (prefix, kind) in (("Vec", "Vec"), ("Mat", "Mat"))
        startswith(s, prefix) || continue
        rest = s[(length(prefix) + 1):end]
        if occursin("Create", rest) ||
           occursin("Duplicate", rest) ||
           occursin("Load", rest) ||
           occursin("WithArray", rest)
            return kind
        end
    end
    return nothing
end

"""
    audit_destroyer(name::Symbol) -> Bool

Whether `name` releases a PETSc object. Covers the high-level `destroy`/`destroy!`
and the low-level `VecDestroy`, `MatDestroy`, `DMDestroy` and friends, which the
previous text-matching version treated as leaks.

`finalizer` counts too. `finalizer(destroy, v)` hands the release to the garbage
collector rather than performing it, but the object is accounted for and must not
read as a leak. The package uses that idiom for sequential objects.
"""
function audit_destroyer(name::Symbol)
    name in (:destroy, :destroy!, :finalizer) && return true
    s = String(name)
    return endswith(s, "Destroy") && length(s) > length("Destroy")
end

# Syntax tree helpers
# ————————————————————————————————————————————————————————————————————————————

"""
    audit_isbroadcast(ex) -> Bool

Whether `ex` is a dot-call such as `destroy!.(vs)`. Those parse as `Expr(:.)`
with a tuple of arguments, not as `Expr(:call)`, so they need recognising
separately. Plain field access `a.b` carries a `QuoteNode` instead and is not a
call.
"""
audit_isbroadcast(ex) =
    ex isa Expr &&
    ex.head === :. &&
    length(ex.args) == 2 &&
    ex.args[2] isa Expr &&
    ex.args[2].head === :tuple

"""
    audit_callee(ex) -> Union{Nothing, Symbol}

The bare name a call expression invokes, discarding any module qualification, so
`f(x)`, `PETSc.f(x)`, `PETSc.LibPETSc.f(x)` and `PETSc.f.(xs)` all yield `:f`.
"""
function audit_callee(ex)
    ex isa Expr || return nothing
    if ex.head === :call
        f = ex.args[1]
    elseif audit_isbroadcast(ex)
        f = ex.args[1]
    else
        return nothing
    end
    f isa Symbol && return f
    while f isa Expr && f.head === :.
        f = f.args[2]
    end
    f isa QuoteNode && (f = f.value)
    return f isa Symbol ? f : nothing
end

"""
    audit_argnames(ex) -> Vector{Symbol}

The plain names a call passes, for either call form, looking through splats and
one level of array or tuple literal. `destroy!(v)`, `destroy!(v...)` and
`destroy!.([v, w])` all name the objects they release.
"""
function audit_argnames(ex)
    args = if ex.head === :call
        ex.args[2:end]
    elseif audit_isbroadcast(ex)
        ex.args[2].args
    else
        return Symbol[]
    end

    names = Symbol[]
    for a in args
        if a isa Expr && a.head === :...
            a = a.args[1]
        end
        if a isa Symbol
            push!(names, a)
        elseif a isa Expr && a.head in (:vect, :tuple)
            append!(names, (x for x in a.args if x isa Symbol))
        end
    end
    return names
end

"""
    audit_targets(lhs) -> Vector{Symbol}

The variables an assignment binds. Handles `x = …`, `x, y = …` and `(x, y) = …`,
skipping anything that is not a plain name.
"""
function audit_targets(lhs)
    lhs isa Symbol && return [lhs]
    if lhs isa Expr && lhs.head in (:tuple, :block)
        return Symbol[a for a in lhs.args if a isa Symbol]
    end
    return Symbol[]
end

"""
    audit_walk(f, ex, line = 0) -> Int

Apply `f(expr, line)` to every `Expr` in `ex`, tracking the source line from the
`LineNumberNode`s the parser emits. Returns the last line seen.
"""
function audit_walk(f, ex, line::Int = 0)
    if ex isa LineNumberNode
        return ex.line
    elseif ex isa Expr
        # Quoted code is data, not execution. A `destroy!(v)` in a macro body
        # releases nothing where it is written, and counting it would report a
        # genuinely leaking file as clean.
        ex.head === :quote && return line
        f(ex, line)
        for a in ex.args
            line = audit_walk(f, a, line)
        end
    end
    return line
end

"""
    audit_hasparseerror(ex) -> Bool

Whether the parsed tree contains an error or incomplete node, which happens when
the file has a syntax error. `Meta.parseall` reports those in the tree rather
than throwing.
"""
function audit_hasparseerror(ex)
    ex isa Expr || return false
    ex.head in (:error, :incomplete) && return true
    return any(audit_hasparseerror, ex.args)
end

# ————————————————————————————————————————————————————————————————————————————

"""
    audit_petsc_file(path::AbstractString; verbose::Bool = true)

Scan a Julia source file for PETSc objects that are created but never destroyed.

Creations are calls to a type constructor (`KSP`, `SNES`, `DMDA`, `DMStag`,
`DMPlex`), to a `Vec`/`Mat` creation routine (`VecCreateSeq`, `MatDuplicate`,
`MatSeqAIJWithArrays`, …), or to one of the DM allocators (`DMGlobalVec`,
`DMLocalVec`, `DMCreateMatrix`, …), through either `PETSc` or `LibPETSc`.
Releases are `destroy`/`destroy!` and the low-level `VecDestroy`-style routines.

Returns a `NamedTuple`:

- `created`: `(line, var, kind)` for each creation, `var === nothing` when the
  result is not assigned
- `destroyed`: `(line, var)` for each release
- `finalized`: `(line, var)` for each `finalize` call
- `leaked`: variables that are created and never released

Pass `verbose = false` to suppress the printed report.

# Notes

Heuristic: it does not follow control flow, scopes, or aliasing, and creations
whose result is not assigned cannot be matched against a release.

# Examples

```julia
julia> report = audit_petsc_file("examples/ex1.jl");

julia> isempty(report.leaked)
true
```
"""
function audit_petsc_file(path::AbstractString; verbose::Bool = true)
    ast = Meta.parseall(read(path, String); filename = path)

    # A file that does not parse yields no creations, which would otherwise be
    # reported as "no leaks" and read as a clean bill of health.
    if audit_hasparseerror(ast)
        @warn "$(path) does not parse; the audit below covers only what could be read"
    end

    created = Tuple{Int, Union{Nothing, Symbol}, String}[]
    destroyed = Tuple{Int, Symbol}[]
    finalized = Tuple{Int, Symbol}[]
    assigned = Base.IdSet{Any}()

    audit_walk(ast) do ex, line
        if ex.head === :(=) && ex.args[2] isa Expr
            rhs = ex.args[2]
            callee = audit_callee(rhs)
            if callee !== nothing
                kind = audit_creator(callee)
                if kind !== nothing
                    push!(assigned, rhs)
                    targets = audit_targets(ex.args[1])
                    if isempty(targets)
                        push!(created, (line, nothing, kind))
                    else
                        for t in targets
                            push!(created, (line, t, kind))
                        end
                    end
                end
            end
            return
        end

        callee = audit_callee(ex)
        callee === nothing && return

        if audit_creator(callee) !== nothing && !(ex in assigned)
            push!(created, (line, nothing, audit_creator(callee)))
        end

        # The object is the first argument high-level (`destroy!(v)`) but the
        # second low-level (`VecDestroy(petsclib, v)`), so record every plain
        # name the call mentions. `petsclib` is never a tracked object, so the
        # extra entries cannot mask a leak.
        if audit_destroyer(callee) || callee === :finalize
            sink = callee === :finalize ? finalized : destroyed
            for a in audit_argnames(ex)
                push!(sink, (line, a))
            end
        end
    end

    created_vars = Set{Symbol}(v for (_, v, _) in created if v !== nothing)
    released = Set{Symbol}(v for (_, v) in destroyed)
    leaked = sort!(collect(setdiff(created_vars, released)))

    # Only report releases of objects this file created; the rest are arguments
    # that happened to sit in a destroy call, such as `petsclib`.
    filter!(((_, v),) -> v in created_vars, destroyed)

    verbose && audit_report(created, destroyed, finalized, leaked)

    return (
        created = created,
        destroyed = destroyed,
        finalized = finalized,
        leaked = leaked,
    )
end

"""
    audit_report(created, destroyed, finalized, leaked)

Print the human-readable form of an [`audit_petsc_file`](@ref) result.
"""
function audit_report(created, destroyed, finalized, leaked)
    println("CREATION statements:")
    for (line, var, kind) in created
        if var === nothing
            println("  line $(line): $(kind)(...) (no assignment)")
        else
            println("  line $(line): $(var) = $(kind)(...)")
        end
    end

    println("DESTROY calls:")
    for (line, var) in destroyed
        println("  line $(line): destroy($(var))")
    end

    if isempty(leaked)
        println("No obvious leaks: all assigned creations have a destroy call.")
    else
        println("POSSIBLY UNDESTROYED objects:")
        for v in leaked
            println("  $(v)")
        end
    end

    if isempty(finalized)
        println("FINALIZE calls: none detected")
        println("Suggestion: add PETSc.finalize(petsclib) at the end of the routine.")
    else
        println("FINALIZE calls:")
        for (line, var) in finalized
            println("  line $(line): finalize($(var))")
        end
    end

    return nothing
end
