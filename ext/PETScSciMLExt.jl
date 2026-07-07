module PETScSciMLExt

using PETSc
using SciMLBase

# Algorithm types live in PETSc proper (`src/sciml_algorithms.jl`)
# so users can write `PETSc.TSRK("3bs")`.
using PETSc: PETScTSAlgorithm, TSRK, TSRosW, TSImplicit, TSARKIMEX, TSGeneric

include("sciml/options.jl")
include("sciml/integrator.jl")
include("sciml/interface.jl")
include("sciml/retcode.jl")
include("sciml/helpers.jl")
include("sciml/rhs_callback.jl")
include("sciml/ifunction_callback.jl")
include("sciml/solve.jl")

end # module
