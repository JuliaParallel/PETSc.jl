# Standalone script: run as a subprocess and check it exits cleanly.
#
# Regression guard for the at-exit finalizer crash: integrators created via
# `init` but never `solve!`-ed / `destroy`-ed are reclaimed by GC finalizers,
# which can fire at process exit *after* the `atexit` `PetscFinalize` /
# `MPI_Finalize` hooks. `_destroy_petsc!` must then be a no-op (guarded by
# `PETSc.finalized`), otherwise `TSDestroy` calls MPI after finalize and
# aborts under MPICH (PETSC_ERR_MPI = 98), failing the whole test run.
#
# The parent test asserts this process exits 0 and prints no
# "error in running finalizer" / "after finalizing" banners.
using PETSc
using SciMLBase

prob = ODEProblem((du, u, p, t) -> (du[1] = -u[1]; nothing), [1.0], (0.0, 1.0))

# Abandon several integrators across both callback paths (explicit RHS and the
# implicit IFunction path) without cleaning them up.
for _ in 1:5
    init(prob, PETSc.TSRK("3bs"); dt = 0.1)
    init(prob, PETSc.TSImplicit("beuler", ["-snes_fd"]); dt = 0.1)
end

# Force a GC pass so finalizers are queued, then let the process exit; the
# atexit PetscFinalize must turn the remaining finalizers into no-ops.
GC.gc()
GC.gc()
println("finalizer_safety: reached clean exit")
