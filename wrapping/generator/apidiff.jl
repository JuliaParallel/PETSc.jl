# Compare two API snapshots and report what a regeneration will change.
#
#   julia --project=wrapping/generator wrapping/generator/apidiff.jl api/petsc-OLD.json api/petsc-NEW.json
#
# Reports new, removed and signature-changed functions (restricted to the ones the generator
# wraps), and rules/overrides that refer to functions or arguments absent from the new snapshot.
include(joinpath(@__DIR__, "src", "PetscWrapGen.jl"))
using .PetscWrapGen
using TOML

old = load_api(ARGS[1]); new = load_api(ARGS[2])
r = load_rules(joinpath(@__DIR__, "rules"))
wrapped(api) = Set(n for n in keys(api.functions) if !occursin("_", n) && !haskey(r.exclude, n))
wo, wn = wrapped(old), wrapped(new)
sig(fn) = [(a.typename, a.stars, a.array, a.isconst, a.name) for a in fn.args]

added = sort!(collect(setdiff(wn, wo))); removed = sort!(collect(setdiff(wo, wn)))
changed = sort!([n for n in intersect(wo, wn) if sig(old.functions[n]) != sig(new.functions[n])])
println("PETSc $(old.version) -> $(new.version): $(length(wo)) -> $(length(wn)) wrapped functions")
println("  new: $(length(added))   removed: $(length(removed))   signature changed: $(length(changed))")
println("\n--- removed:"); foreach(println, removed)
println("\n--- signature changed:")
for n in changed
    o = join(("$(a.typename)$("*"^a.stars) $(a.name)$(a.array ? "[]" : "")" for a in old.functions[n].args), ", ")
    m = join(("$(a.typename)$("*"^a.stars) $(a.name)$(a.array ? "[]" : "")" for a in new.functions[n].args), ", ")
    println(n, "\n    old: ", o, "\n    new: ", m)
end
println("\n--- new:"); foreach(println, added)

println("\n--- rules referring to functions/arguments missing in the new snapshot:")
stale = 0
for (fn, argrules) in r.args
    if !haskey(new.functions, fn)
        println("  args.toml [$fn.*]: function not in snapshot"); global stale += 1; continue
    end
    names = Set(PetscWrapGen.rename_arg(r, a.name) for a in new.functions[fn].args)
    for a in keys(argrules)
        a in names || (println("  args.toml [$fn.$a]: argument not in snapshot"); global stale += 1)
    end
end
for fn in keys(r.exclude)
    haskey(new.functions, fn) || (println("  files.toml [exclude] $fn: function not in snapshot"); global stale += 1)
end
ovdir = joinpath(@__DIR__, "overrides")
for f in readdir(ovdir)
    endswith(f, ".jl") || continue
    n = f[1:end-3]
    haskey(new.functions, n) || continue     # hand-written extras are fine
    first_line = first(eachline(joinpath(ovdir, f)))
    m = match(r"C signature: \w+\((.*)\)\s*$", first_line)
    m === nothing && continue
    cur = join(("$(a.typename)$("*"^a.stars) $(a.name)$(a.array ? "[]" : "")" for a in new.functions[n].args), ", ")
    cur == m.captures[1] || (println("  overrides/$f: C signature changed\n      was: $(m.captures[1])\n      now: $cur"); global stale += 1)
end
println("$stale stale entries")
