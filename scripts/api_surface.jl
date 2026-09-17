#!/usr/bin/env julia
#
# The high-level API surface (docs/src/man/naming.md §1.1).
#
#   julia --project=. scripts/api_surface.jl --check
#   julia --project=. scripts/api_surface.jl --sweeps
#
# `--check` enumerates the surface and fails if a binding is in none of the
# register's three sets (renamed, internal, unchanged-public). The register is
# `scripts/renames.jl`; a name missing from it is added there rather than
# excused here, because five generated files are derived from the same data.
#
# `--sweeps` derives the four lists §1.1 asks for, so that no count in
# `naming.md` has to be remembered:
#
#   1. functions taking `petsclib` first alongside a dispatchable PETSc object (§8)
#   2. functions returning a `NamedTuple` (§12)
#   3. exported names against §13's list
#   4. readers returning a PETSc object without a `doc_borrowed` entry (§3.3)
#
# The module is loaded rather than parsed, because §13's short accessors (`dm`,
# `comm`, `info`, `ds`, `label`, `solution`) cannot be found by a text search.

module APISurface

using PETSc

include(joinpath(@__DIR__, "renames.jl"))

const SRC = normpath(joinpath(@__DIR__, "..", "src"))

# ---------------------------------------------------------------------------
# The surface
# ---------------------------------------------------------------------------

"""
    own_source(file) -> Bool

Whether a method was defined in PETSc's own high-level source: `src/*.jl`, but
not the generated `src/autowrapped/` layer and not the `LibPETSc` plumbing.
"""
function own_source(file)
    path = normpath(String(file))
    startswith(path, SRC) || return false
    rel = relpath(path, SRC)
    startswith(rel, "autowrapped") && return false
    startswith(rel, "LibPETSc") && return false
    return true
end

has_own_method(x) = any(m -> own_source(m.file), methods(x))

# Bindings that exist in every module and name nothing of this package's.
const MODULE_BUILTINS = Set{Symbol}([:PETSc, :eval, :include])

"""
    surface() -> Vector{Symbol}

`names(PETSc; all = true)` filtered to bindings whose methods are defined in
PETSc's own source files. A name whose binding `PETSc` does not own (`ndims`,
`size`, `show`: the Base extensions of §11) is not part of the surface, and
neither is anything that only has generated or `LibPETSc` methods.
"""
function surface()
    out = Symbol[]
    for name in names(PETSc; all = true)
        str = String(name)
        (startswith(str, "#") || name in MODULE_BUILTINS) && continue
        # A Base extension: the binding belongs to Base (or to another module),
        # and §1.1 exempts it by definition.
        Base.which(PETSc, name) === PETSc || continue
        isdefined(PETSc, name) || continue
        value = getproperty(PETSc, name)
        value isa Module && continue
        # Types and functions carry methods; everything else (abstract types,
        # constants) is on the surface as soon as `PETSc` owns the binding.
        if value isa Function || value isa Type
            has_own_method(value) ||
                (value isa Type && isabstracttype(value)) ||
                continue
        end
        push!(out, name)
    end
    return sort!(out; by = String)
end

# ---------------------------------------------------------------------------
# --check
# ---------------------------------------------------------------------------

"""
    registered() -> Set{Symbol}

Every name the register accounts for: both sides of `RENAMES`, the internal set,
the unchanged-public list, the exports, and the `Base` targets.
"""
function registered()
    known = Set{Symbol}()
    for (old, new) in RENAMES
        push!(known, old)
        push!(known, new)
    end
    union!(known, INTERNAL)
    union!(known, UNCHANGED_PUBLIC)
    union!(known, EXPORTED)
    union!(known, BASE_TARGETS)
    return known
end

"""
    unregistered() -> Vector{Symbol}

The bindings on the surface that the register does not account for. Empty is the
only acceptable value; `--check` and the test suite both assert it.
"""
unregistered() = filter(!in(registered()), surface())

"""
    check(io = stdout) -> Bool

Print the offenders and return whether the surface is fully registered.
"""
function check(io::IO = stdout)
    absent = unregistered()
    if isempty(absent)
        println(
            io,
            "api_surface --check: $(length(surface())) bindings, all registered.",
        )
        return true
    end
    println(
        io,
        "api_surface --check: $(length(absent)) binding(s) in no set of scripts/renames.jl:",
    )
    for name in absent
        println(io, "  ", name)
    end
    println(io, "\nAdd each to RENAMES, INTERNAL or UNCHANGED_PUBLIC in scripts/renames.jl.")
    return false
end

# ---------------------------------------------------------------------------
# --sweeps
# ---------------------------------------------------------------------------

const PETSC_OBJECT_TYPES = Any[
    PETSc.LibPETSc.AbstractPetscVec,
    PETSc.LibPETSc.AbstractPetscMat,
    PETSc.LibPETSc.AbstractPetscDM,
    PETSc.LibPETSc.AbstractKSP,
    PETSc.LibPETSc.AbstractSNES,
    PETSc.LibPETSc.AbstractTS,
    PETSc.LibPETSc.AbstractIS,
    PETSc.LibPETSc.AbstractPetscOptions,
    PETSc.LibPETSc.AbstractTao,
    PETSc.LibPETSc.AbstractAO,
]

is_petsc_object_type(T) =
    T isa Type && T !== Union{} && any(A -> T <: A, PETSC_OBJECT_TYPES)

is_petsclib_type(T) = T isa Type && T <: PETSc.LibPETSc.PetscLibType

# The declared argument types of a method, `self` dropped.
function argtypes(m::Method)
    sig = Base.unwrap_unionall(m.sig)
    return sig isa DataType ? collect(sig.parameters)[2:end] : Any[]
end

"""
    petsclib_first() -> Vector

Sweep 1 (§8): methods whose first argument is a `PetscLibType` while a later
argument is a dispatchable PETSc object, so the library could have been read off
the object instead.
"""
function petsclib_first()
    hits = Tuple{Symbol, Method}[]
    for name in surface()
        f = getproperty(PETSc, name)
        (f isa Function || f isa Type) || continue
        for m in methods(f)
            own_source(m.file) || continue
            ts = argtypes(m)
            length(ts) >= 2 || continue
            is_petsclib_type(ts[1]) || continue
            any(is_petsc_object_type, ts[2:end]) || continue
            push!(hits, (name, m))
        end
    end
    return hits
end

# The return types inference can see for every own-source method of `f`.
function return_types_of(f)
    out = Any[]
    for m in methods(f)
        own_source(m.file) || continue
        try
            sig = Base.unwrap_unionall(m.sig)
            sig isa DataType || continue
            append!(
                out,
                Base.return_types(f, Tuple{collect(sig.parameters)[2:end]...}),
            )
        catch
        end
    end
    return out
end

flatten_union(T) = T isa Union ? vcat(flatten_union(T.a), flatten_union(T.b)) : Any[T]

mentions(T, pred) = any(pred, flatten_union(T))

is_namedtuple_type(T) = T isa Type && T !== Union{} && T <: NamedTuple

"""
    namedtuple_returns() -> Vector{Symbol}

Sweep 2 (§12): the v0.5 functions that return a `NamedTuple`. The deprecated
spellings forward to them and the internal helpers are not API, so neither is
counted.
"""
function namedtuple_returns()
    deprecated = Set(first.(RENAMES))
    hits = Symbol[]
    for name in surface()
        (name in deprecated || name in INTERNAL) && continue
        f = getproperty(PETSc, name)
        (f isa Function || f isa Type) || continue
        any(T -> mentions(T, is_namedtuple_type), return_types_of(f)) &&
            push!(hits, name)
    end
    return hits
end

"""
    export_sweep() -> NamedTuple

Sweep 3 (§13): what `PETSc` exports against the register's `EXPORTED` list.
"""
function export_sweep()
    # The module's own name is reported as exported; it is not an export.
    actual = Set(
        n for n in names(PETSc; all = true) if
        n !== :PETSc && Base.isexported(PETSc, n)
    )
    declared = Set(EXPORTED)
    return (
        exported = sort!(collect(actual); by = String),
        unexpected = sort!(collect(setdiff(actual, declared)); by = String),
        absent = sort!(collect(setdiff(declared, actual)); by = String),
    )
end

"""
    borrowed_documented() -> Set{Symbol}

The names whose docstring interpolates `doc_borrowed()`, found by reading the
source: `Docs.doc` on a binding with several docstrings answers with a summary
rather than the text, and §1.1 asks for a grep for the helper anyway.
"""
function borrowed_documented()
    documented = Set{Symbol}()
    for file in readdir(SRC; join = true)
        endswith(file, ".jl") || continue
        lines = readlines(file)
        for (i, line) in pairs(lines)
            occursin("doc_borrowed()", line) || continue
            # Walk to the end of the docstring, then to the definition it documents.
            j = i
            while j <= length(lines) && strip(lines[j]) != "\"\"\""
                j += 1
            end
            j += 1
            while j <= length(lines) && isempty(strip(lines[j]))
                j += 1
            end
            j <= length(lines) || continue
            m = match(r"^\s*(?:function\s+)?([A-Za-z_][A-Za-z0-9_!]*)", lines[j])
            m === nothing || push!(documented, Symbol(m.captures[1]))
        end
    end
    return documented
end

"""
    borrowed_sweep() -> Vector{Symbol}

Sweep 4 (§3.3): readers (names without a trailing `!`, and not creators) that
hand back a PETSc object without a `doc_borrowed` note in their docstring. The
deprecated spellings and the internal helpers are not part of the documented
API and are skipped.
"""
function borrowed_sweep()
    creators = Set(keys(AUDIT_TYPE_CREATORS))
    union!(creators, keys(AUDIT_NAMED_CREATORS))
    deprecated = Set(first.(RENAMES))
    documented = borrowed_documented()
    hits = Symbol[]
    for name in surface()
        str = String(name)
        (endswith(str, "!") || startswith(str, "@")) && continue
        (name in creators || name in deprecated || name in INTERNAL) && continue
        name in documented && continue
        f = getproperty(PETSc, name)
        f isa Function || continue
        any(T -> mentions(T, is_petsc_object_type), return_types_of(f)) || continue
        push!(hits, name)
    end
    return hits
end

function sweeps(io::IO = stdout)
    println(io, "1. petsclib first alongside a dispatchable PETSc object (§8)")
    hits = petsclib_first()
    isempty(hits) && println(io, "   (none)")
    for (name, m) in hits
        println(io, "   ", name, "  ", basename(String(m.file)), ":", m.line)
    end

    println(io, "\n2. functions returning a NamedTuple (§12)")
    nts = namedtuple_returns()
    println(io, "   ", length(nts), " function(s)")
    for name in nts
        println(io, "   ", name)
    end

    println(io, "\n3. exported names against §13")
    es = export_sweep()
    println(io, "   exported: ", join(String.(es.exported), ", "))
    println(
        io,
        "   not in the register's EXPORTED: ",
        isempty(es.unexpected) ? "(none)" : join(String.(es.unexpected), ", "),
    )
    println(
        io,
        "   declared but not exported: ",
        isempty(es.absent) ? "(none)" : join(String.(es.absent), ", "),
    )

    println(io, "\n4. readers returning a PETSc object with no doc_borrowed entry (§3.3)")
    bs = borrowed_sweep()
    isempty(bs) && println(io, "   (none)")
    for name in bs
        println(io, "   ", name)
    end
    return nothing
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    if "--sweeps" in ARGS
        APISurface.sweeps()
    elseif isempty(ARGS) || "--check" in ARGS
        exit(APISurface.check() ? 0 : 1)
    else
        println("usage: api_surface.jl [--check | --sweeps]")
        exit(2)
    end
end
