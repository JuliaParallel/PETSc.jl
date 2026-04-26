# Algorithm types live in PETSc proper (`src/sciml_algorithms.jl`) so users
# can write `PETSc.TSRK("3bs")` without going through `Base.get_extension`.
# This file just imports them into the extension's namespace.
using PETSc: PETScTSAlgorithm, TSRK, TSRosW, TSImplicit, TSARKIMEX, TSGeneric
