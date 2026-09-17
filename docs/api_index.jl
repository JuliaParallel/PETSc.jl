# The C-name-to-Julia-name index page (docs/src/man/naming.md §15).
#
# Renaming a thin wrapper costs discoverability: someone who knows
# `PetscFECopyQuadrature` cannot grep for it once it is `copy_quadrature!`. Every
# high-level docstring therefore carries a `doc_external("Section/CName")` entry,
# and this page is built from those entries at docs build time, so the lookup
# table maintains itself and cannot drift from the docstrings.
#
# The page is generated, and `docs/src/man/api_index.md` is ignored by git.

const HIGH_LEVEL_SOURCES = [
    "init.jl",
    "vec.jl",
    "mat.jl",
    "options.jl",
    "ts.jl",
    "ksp.jl",
    "snes.jl",
    "dm.jl",
    "sys.jl",
    "dmda.jl",
    "dmstag.jl",
    "dmplex.jl",
    "audit.jl",
]

# `doc_external("DMSTAG/DMStagGetCorners")`, and the `_doc_external` alias the
# generated layer uses.
const DOC_EXTERNAL_RE = r"_?doc_external\(\"([^\"]+)\"\)"

# The name a definition line defines, whatever wraps it: `function f(`,
# `LibPETSc.@for_petsc function f(`, `macro f(`, `mutable struct F{`, `f(x) = …`.
const DEFINITION_RE = r"^\s*(?:LibPETSc\.@for_petsc\s+)?(?:@for_petsc\s+)?(?:function|macro|(?:mutable\s+)?struct|const|abstract\s+type)?\s*((?:Base\.|LibPETSc\.)?@?[A-Za-z_][A-Za-z0-9_!]*)"

"""
    collect_doc_external(srcdir) -> Dict{String, Vector{String}}

Map each C function named in a `doc_external` entry to the Julia names whose
docstrings mention it. A docstring's entries are attached to the definition that
follows it.
"""
function collect_doc_external(srcdir::AbstractString)
    index = Dict{String, Vector{String}}()
    for file in HIGH_LEVEL_SOURCES
        path = joinpath(srcdir, file)
        isfile(path) || continue
        lines = readlines(path)
        pending = String[]
        i = 1
        while i <= length(lines)
            # Skip `#= … =#` blocks: they hold code that is not compiled, and a
            # docstring inside one documents nothing.
            if startswith(strip(lines[i]), "#=")
                while i <= length(lines) && !occursin("=#", lines[i])
                    i += 1
                end
                i += 1
                continue
            end
            m = match(DOC_EXTERNAL_RE, lines[i])
            if m === nothing
                i += 1
                continue
            end
            # Gather every entry in this docstring, then find the definition.
            empty!(pending)
            while i <= length(lines) && strip(lines[i]) != "\"\"\""
                for e in eachmatch(DOC_EXTERNAL_RE, lines[i])
                    push!(pending, e.captures[1])
                end
                i += 1
            end
            i += 1  # past the closing """
            while i <= length(lines) && isempty(strip(lines[i]))
                i += 1
            end
            i <= length(lines) || break
            d = match(DEFINITION_RE, lines[i])
            d === nothing && continue
            julia_name = d.captures[1]
            for entry in pending
                names = get!(index, entry, String[])
                julia_name in names || push!(names, julia_name)
            end
        end
    end
    return index
end

# The type constructors `PETSc` exports, written on `LibPETSc`'s types.
const REEXPORTED = Set(["PetscVec", "PetscMat", "PetscOptions"])

function reexported_name(n::AbstractString)
    bare = last(split(n, '.'))
    bare in REEXPORTED && return "PETSc." * bare
    occursin('.', n) && return String(n)
    return "PETSc." * n
end

# "DMSTAG/DMStagGetCorners" -> ("DMStagGetCorners", the petsc.org URL)
function split_entry(entry::AbstractString)
    cname = last(split(entry, '/'))
    url = "https://petsc.org/release/docs/manualpages/$entry.html"
    return String(cname), url
end

"""
    write_api_index(path; srcdir = "../src")

Write the C-to-Julia index page.
"""
function write_api_index(
    path::AbstractString;
    srcdir::AbstractString = normpath(joinpath(@__DIR__, "..", "src")),
)
    index = collect_doc_external(srcdir)
    rows = sort!(collect(index); by = p -> lowercase(last(split(first(p), '/'))))

    open(path, "w") do io
        println(io, "# C to Julia name index")
        println(io)
        println(
            io,
            """
            The high-level interface renames PETSc's C functions to Julia
            conventions (see [the naming conventions](naming.md)), which costs
            discoverability: someone who knows `PetscFECopyQuadrature` cannot
            grep for it once it is `copy_quadrature!`. This page is the lookup
            table in the other direction.

            It is generated at docs build time from the `doc_external` entries in
            the high-level docstrings, so it cannot drift from them: a wrapper
            that links a C page appears here, and one that does not is a
            docstring bug rather than a missing table row.

            $(length(rows)) C functions are covered.
            """,
        )
        println(io, "| C function | Julia |")
        println(io, "|---|---|")
        for (entry, names) in rows
            cname, url = split_entry(entry)
            # A method added to `Base` or to `LibPETSc` keeps its own module,
            # except for the type constructors `PETSc` exports (§13): those are
            # written on `LibPETSc`'s types but called as `PetscVec(…)`.
            qualify(n) = "`$(reexported_name(n))`"
            julia = join((qualify(n) for n in sort(names)), ", ")
            println(io, "| [`$cname`]($url) | $julia |")
        end
    end
    return path
end
