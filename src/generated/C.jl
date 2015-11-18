module C

using Compat
export PetscInt

include("defs.jl")
include("libPETSc_commonRealDouble.jl")
include("PETScRealDouble.jl")
include("PETScRealSingle.jl")
include("PETScComplexDouble.jl")
include("defs2.jl")

end
