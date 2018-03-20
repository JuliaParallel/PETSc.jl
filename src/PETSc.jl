module PETSc

import MPI
export PetscInt
include(joinpath("generated", "C.jl"))
using .C
include("petsc_com.jl")
include("options.jl")
include("is.jl")
include("vec.jl")
include("mat.jl")
include("vec_scatter.jl")
include("pc.jl")
include("ksp.jl")
include("mapping.jl")
include("ts.jl")
end
