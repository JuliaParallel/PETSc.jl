# Mine the PETSc source for the handle arguments a function checks unconditionally, and
# write rules/args_required.toml.
#
#   julia --project=wrapping/generator wrapping/generator/mine_required.jl PETSC_DIR API_SNAPSHOT
#
# A debug PETSc rejects a null object passed where a function calls
# `PetscValidHeaderSpecific(x, CLASSID, k)` (or `PetscValidHeader`, `...SpecificType`) on its
# k-th argument. An optimized build, such as PETSc_jll, skips the check and crashes instead,
# so the generated wrapper makes the same check. Only calls at the body's top level count
# (two-space indent): one inside `if (x) { ... }` is for an optional argument. An argument the
# body also tests (`if (x)`, `if (!x)`, `x ?`, `x &&`) is skipped, as is a dereference such as
# `*v` in `XDestroy`, which returns early on NULL. A function that passes an argument on to
# one that validates it, `PetscCall(VecScaleAsync_Private(x, alpha, NULL))`, validates it too.
include(joinpath(@__DIR__, "src", "PetscWrapGen.jl"))
using .PetscWrapGen

const VALIDATE = r"^  PetscValidHeader(?:SpecificType|Specific)?\(\s*([^,]+?)\s*,(?:[^,]*,)?\s*(\d+)\s*[,)]"m
const DEFINITION = r"^PetscErrorCode\s+(\w+)\(([^)]*)\)\s*\n\{\n(.*?)\n\}"ms
const FORWARD = r"^  PetscCall\((\w+)\(([^()]*)\)\);"m

# The body tests `x` itself, so NULL is allowed there
tests_null(body, x) = occursin(Regex("if \\(!?" * x * "\\)|\\(!?" * x * "\\) ?\\?|!?\\b" * x * " (\\?|&&|\\|\\|)|(&&|\\|\\|) !?" * x * "\\b(?!->|\\[)"), body)

params(sig) = [String(something(match(r"(\w+)\s*(\[\])?$", strip(p)), (captures = [""],)).captures[1]) for p in split(sig, ",")]

# function name => sorted argument numbers (1-based) PETSc validates unconditionally,
# directly or through a function it forwards the argument to
function mine_source(petsc_dir::AbstractString)
    defs = Dict{String,Tuple{Vector{String},String}}()
    for (root, _, files) in walkdir(joinpath(petsc_dir, "src")), f in files
        endswith(f, ".c") || continue
        (occursin("/tests", root) || occursin("/tutorials", root) || occursin("/ftn-", root)) && continue
        for m in eachmatch(DEFINITION, read(joinpath(root, f), String))
            defs[m.captures[1]] = (params(m.captures[2]), m.captures[3])
        end
    end
    found = Dict{String,Set{Int}}()
    for (name, (_, body)) in defs, v in eachmatch(VALIDATE, body)
        x, k = strip(v.captures[1]), parse(Int, v.captures[2])
        (occursin(r"^\w+$", x) && !tests_null(body, x)) && push!(get!(found, name, Set{Int}()), k)
    end
    changed = true
    while changed
        changed = false
        for (name, (ps, body)) in defs, c in eachmatch(FORWARD, body)
            callee = get(found, c.captures[1], nothing)
            callee === nothing && continue
            cargs = strip.(split(c.captures[2], ","))
            for k in callee
                k <= length(cargs) || continue
                j = findfirst(==(cargs[k]), ps)
                (j === nothing || tests_null(body, ps[j])) && continue
                changed |= !(j in get!(found, name, Set{Int}()))
                push!(found[name], j)
            end
        end
    end
    return Dict(k => sort!(collect(v)) for (k, v) in found)
end

function main(petsc_dir, snapshot)
    api = load_api(snapshot)
    found = mine_source(petsc_dir)
    out = joinpath(@__DIR__, "rules", "args_required.toml")
    n = 0
    open(out, "w") do io
        println(io, "# Handle arguments PETSc checks for NULL in a debug build, MINED from the PETSc source")
        println(io, "# by mine_required.jl: do not edit by hand, use rules/args.toml instead.")
        for fname in sort!(collect(keys(found)))
            fn = get(api.functions, fname, nothing)
            fn === nothing && continue
            for k in found[fname]
                k <= length(fn.args) || continue
                a = fn.args[k]
                (a.stars == 0 && !isempty(a.name)) || continue
                println(io, "\n[$fname.$(a.name)]\nrequired = true")
                n += 1
            end
        end
    end
    println("wrote $n required arguments to $out")
end

main(ARGS...)
