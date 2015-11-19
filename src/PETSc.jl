module PETSc

import MPI
include(joinpath("generated", "C.jl"))
using .C
include("petsc_com.jl")
include("options.jl")
include("vec.jl")
include("mat.jl")
include("error.jl")
include("pc.jl")
include("ksp.jl")
include("is.jl")

end
