# provide some most commonly used options, leave rest as low level
# common options

export SNES

###########################################################################

# solver context
type SNES{T}
  p::C.SNES{T}
  function SNES(p::C.SNES{T})
    o = new(p)
    finalizer(o, SNESDestroy)
    return o
  end
end

comm{T}(a::SNES{T}) = MPI.Comm(C.PetscObjectComm(T, a.p.pobj))

function SNESDestroy{T}(o::SNES{T})
  PetscFinalized(T) || C.SNESDestroy(Ref(o.p))
end

"""
    SNES(A, PA=A, FormJacobian; kws...)

Create a SNES solver object that can be used to solve non-linear equations
with the matrix `A`, where `PA` (defaults to `A`) is used to construct the
default preconditioner.  FormJacobian builds the Jacobian matrix when it
is called, and FormFunction builds the function to solve when it is called.

The keyword options are zero or more of the following:

These control the solver and preconditioner characteristics:
* `snes_type="a"`: use SNES algorithm `a`

The following keyword options control the stopping criteria for
iterative solvers:
* `snes_rtol=x`: `x` is relative decrease in residual norm
* `snes_atol=x`: `x` is absolute decrease in residual norm
* `snes_stol=x`: `x` is absolute decrease in residual norm between solution steps
* `snes_trtol=x`: `x` is the trust region tolerance
* `snes_max_it=n`: `n` is the max number of iterations
* `snes_converged_use_initial_residual_norm=true`: use initial residual norm for computing relative convergence
* `snes_convergence_test=:default` or `:skip`: use the default convergence test (tolerances and `max_it`) or skip convergence tests and run until `max_it` is reached
* `snes_lag_jacobian=true`: lag the calculation of the Jacobian by one iteration (trades off reduced communication for an additional iteration)
* `snes_lag_preconditioner=true`: lag the calculation of the preconditioner by one iteration (trades off reduced communication for an additional iteration)

The following options control output that monitors the progress of the
solver (default none).
* `snes_monitor=filename`: print the residual norm at each iteration to `filename` (`""` for `STDOUT`)
* `snes_monitor_short=filename`: print preconditioned residual norm with fewer digits
* `snes_monitor_solution=true`: plot solution graphically
* `snes_monitor_lg_residualnorm=true`: plot preconditioned residual norm graphically
* `snes_monitor_lg_range=true`: plot preconditioned residual norm and range of residual values
* `snes_monitor_cancel=true`: remove any hardwired monitor routines

In addition, if default preconditioner is being used,
then any of the preconditioner options (see `PC`) can be specified to control
this preconditioner (e.g. `pc_type`).
"""
function SNES{T}(FormJacobian, A::Mat{T}, PA::Mat{T}=A; kws...)
  snes_c = Ref{C.SNES{T}}()
  chk(C.SNESCreate(comm(A), snes_c))
  snes = snes_c[]
  vptr = Ref{Void}()
  #lmao this is broken
  const jacfpointer = cfunction(FormJacobian,PetscErrorCode,(C.SNES{T},C.Vec{T},C.Mat{T},Ptr{Void}))
  chk(C.SNESSetJacobian(snes, A.p, PA.p, jacfpointer, C_NULL))
  #chk(C.SNESSetFunction(snes, A.p, PA.p, [FormFunction], vptr))
  withoptions(T, kws) do
    chk(C.SNESSetFromOptions(snes))
  end
  return SNES{T}(snes)
end

# Retrieve a reference to the matrix in the KSP object, as a raw C.Mat
# pointer.  Note that we should not wrap this in a Mat object, or
# call MatDestroy on it, without incrementing its reference count first!
function _snes_ksp{T}(snes::SNES{T})
  k = Ref{C.KSP{T}}()
  chk(C.SNESGetKSP(snes.p,k))
  return k[]
end
