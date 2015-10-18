module PETSc

import MPI
#generated_path = joinpath(Pkg.dir("PETSc"), "src/auto2/generated")
#push!(LOAD_PATH, generated_path)

include(joinpath("generated", "C.jl"))
using .C
include("petsc_com.jl")
include("vec.jl")
include("mat.jl")
include("error.jl")
include("ksp.jl")


#=
function __init__()

    PetscInitialize()
end
=#

end
