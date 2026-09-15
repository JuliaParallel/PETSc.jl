module PetscWrapGen

using JSON3, TOML

include("api.jl")
include("docs.jl")
include("types.jl")
include("classify.jl")
include("render.jl")
include("support.jl")
include("driver.jl")

export generate, load_api, load_rules, build_docindex, function_docs

end
