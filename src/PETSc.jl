module PETSc

import MPI
export PetscInt
include(joinpath("generated", "C.jl"))
using .C
include("petsc_com.jl")
include("options.jl")
include("vec.jl")
include("is.jl")
include("mat.jl")
include("pc.jl")
include("ksp.jl")
include("mapping.jl")
end
