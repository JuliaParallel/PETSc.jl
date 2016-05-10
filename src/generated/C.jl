module C

using Compat
export PetscInt

include("defs.jl")
include("libPETSc_commonRealDouble.jl")
if have_petsc[1]
  include("PETScRealDouble.jl")
end
if have_petsc[2]
  include("PETScRealSingle.jl")
end
if have_petsc[3]
  include("PETScComplexDouble.jl")
end
include("error.jl")
include("defs2.jl")
include("c_funcs.jl")
end
