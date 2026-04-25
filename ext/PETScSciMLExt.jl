module PETScSciMLExt

using PETSc
using SciMLBase
using DiffEqBase

# This extension is loaded when both SciMLBase and DiffEqBase are present in
# the user's environment (loading OrdinaryDiffEq satisfies both).
#
# Implementation is staged following PLAN_INTERFACE.md:
#   Step 1 — algorithm types, integrator, options, interface contract
#   Step 2 — TSRK (explicit Runge-Kutta) end-to-end solve
#   Step 3 — TSRosW (Rosenbrock-W)
#   Step 4 — TSImplicit (BEULER / CN / Theta / BDF)
#   Step 5 — TSARKIMEX (SplitODEProblem)
#   Step 6 — save_everystep, saveat
#   Step 7 — discrete callbacks, terminate!
#   Step 8 — lifecycle (finalizers, destroy)
#   Step 9 — polish

end # module
