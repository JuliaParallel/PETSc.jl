# Copy hand-written wrapper blocks from a golden autowrapped directory into overrides/NAME.jl.
#
#   julia --project=wrapping/generator wrapping/generator/make_overrides.jl GOLDEN_DIR NAME...
#
# Each override file starts with a comment recording the C signature from the API snapshot, so
# that a later generator run can flag overrides whose C signature changed.
include(joinpath(@__DIR__, "src", "PetscWrapGen.jl"))
using .PetscWrapGen
include(joinpath(@__DIR__, "src", "blocks.jl"))

golden = ARGS[1]
names = ARGS[2:end]
blocks = load_blocks(golden)
api = load_api(joinpath(@__DIR__, "api", "petsc-3.24.0.json"))
mkpath(joinpath(@__DIR__, "overrides"))
for n in names
    key = "fn:$n"
    haskey(blocks, key) || (println("no golden block for $n"); continue)
    sig = haskey(api.functions, n) ? join(("$(a.typename)$("*"^a.stars) $(a.name)$(a.array ? "[]" : "")" for a in api.functions[n].args), ", ") : "<not in API snapshot>"
    open(joinpath(@__DIR__, "overrides", "$n.jl"), "w") do io
        println(io, "# override for $n; C signature: $n($sig)")
        print(io, rstrip(blocks[key].text), "\n\n")
    end
end
println("wrote $(length(names)) overrides")
