# provide some most commonly used options, leave rest as low level
# common options: Orthogonilization type, KSP type, PC type
# use PC context created as part of KSP

export KSP, petscview

###########################################################################

# solver context
type KSP{T}
  p::C.KSP{T}
  function KSP(p::C.KSP{T})
    o = new(p)
    finalizer(o, KSPDestroy)
    return o
  end
end

comm{T}(a::KSP{T}) = MPI.Comm(C.PetscObjectComm(T, a.p.pobj))

function KSPDestroy{T}(o::KSP{T})
  PetscFinalized(T) || C.KSPDestroy(Ref(o.p))
end

function petscview{T}(o::KSP{T})
  viewer = C.PetscViewer{T}(C_NULL)
  chk(C.KSPView(o.p, viewer))
end

"""
    KSP(A::Mat, PA=A; kws...)
    KSP(pc::PC; kws...)

Create a KSP solver object that can be used to solve equations `Ax=b` with
the matrix `A`, where `PA` (defaults to `A`) is used to construct the
default preconditioner.  Alternatively, you can supply a preconditioner
object (`PC`).

The keyword options are zero or more of the following:

These control the solver and preconditioner characteristics:
* `ksp_type="a"`: use KSP algorithm `a`
* `ksp_pc_side=n`: set preconditioner side to `PETSc.C.PC_LEFT`, `PETSc.C.PC_RIGHT`, or `PETSc.C.PC_SYMMETRIC`
* `ksp_reuse_preconditioner=true`: use initial preconditioner and don't ever compute a new one
* `ksp_diagonal_scale=true`: symmetrically diagonally scale `A` before solving (note that this *changes* `A` and the right-hand side in a solve, unless you also set `ksp_diagonal_scale_fix=true`)
* `ksp_diagonal_scale_fix=true`: undo diagonal scaling after solve
* `ksp_knoll=true`: use preconditioner applied to `b` for initial guess
* `ksp_constant_null_space=true`: add constant null space to Krylov solver matrix
* `ksp_initial_guess_nonzero=true`: use the contents of initial `x` instead of zero for initial guess
* `ksp_fischer_guess="model,size"`: use Fischer initial guess generator (`model=1` or `2`) for repeated linear solves with subspace of dimension `size`

The following keyword options control the stopping criteria for
iterative solvers:
* `ksp_rtol=x`: `x` is relative decrease in residual norm
* `ksp_atol=x`: `x` is absolute decrease in residual norm
* `ksp_divtol=x`: `x` is amount residual can increase before method is considered to be diverging
* `ksp_max_it=n`: `n` is the max number of iterations
* `ksp_converged_use_initial_residual_norm=true`: use initial residual norm for computing relative convergence
* `ksp_converged_use_min_initial_residual_norm=true`: use min of initial residual norm and `b` for computing relative convergence
* `ksp_error_if_not_converged=true`: generate error if solver does not converge
* `ksp_convergence_test=:default` or `:skip`: use the default convergence test (tolerances and `max_it`) or skip convergence tests and run until `max_it` is reached
* `ksp_norm_type=n`: in residual tests, use norm type `n`, one of default (`PETSc.C.KSP_NORM_DEFAULT`), none (`PETSc.C.KSP_NORM_NONE`), of the preconditioned residual (`PETSc.C.KSP_NORM_PRECONDITIONED`), the true residual (`PETSc.C.KSP_NORM_UNPRECONDITIONED`), or the "natural" norm (`PETSc.C.KSP_NORM_NATURAL`)
* `ksp_check_norm_iteration=n`: compute residual norm starting on iteration `n`
* `ksp_lag_norm=true`: lag the calculation of the residual norm by one iteration (trades off reduced communication for an additional iteration)

The following options control output that monitors the progress of the
solver (default none).
* `ksp_monitor=filename`: print the residual norm at each iteration to `filename` (`""` for `STDOUT`)
* `ksp_monitor_short=filename`: print preconditioned residual norm with fewer digits
* `ksp_monitor_range=filename`: prints the percentage of residual elements that are more then 10% of the maximum value
* `ksp_monitor_true_residual=filename`: print true residual norm
* `ksp_monitor_singular_value=filename`: print extreme singular values (via Lanczos or Arnoldi process as the linear system is solved)
* `ksp_monitor_solution=true`: plot solution graphically
* `ksp_monitor_lg_residualnorm=true`: plot preconditioned residual norm graphically
* `ksp_monitor_lg_true_residualnorm=true`: plot preconditioned and true residual norm graphically
* `ksp_monitor_lg_range=true`: plot preconditioned residual norm and range of residual values
* `ksp_monitor_cancel=true`: remove any hardwired monitor routines
* `ksp_compute_singularvalues=true`: print extreme singular values (via Lanczos or Arnoldi process as the linear system is solved)

In addition, if default preconditioner is being used,
then any of the preconditioner options (see `PC`) can be specified to control
this preconditioner (e.g. `pc_type`).
"""
function KSP{T}(pc::PC{T}; kws...)
  ksp_c = Ref{C.KSP{T}}()
  chk(C.KSPCreate(comm(pc), ksp_c))
  ksp = ksp_c[]
  chk(C.KSPSetPC(ksp, pc.p))
  withoptions(T, kws) do
    chk(C.KSPSetFromOptions(ksp))
  end
  return KSP{T}(ksp)
end

KSP{T}(A::Mat{T}, PA::Mat{T}=A; kws...) = KSP(PC(A, PA; kws...))

# Retrieve a reference to the matrix in the KSP object, as a raw C.Mat
# pointer.  Note that we should not wrap this in a Mat object, or
# call MatDestroy on it, without incrementing its reference count first!
function _ksp_A{T}(ksp::KSP{T})
  a = Ref{C.Mat{T}}()
  pa = Ref{C.Mat{T}}()
  chk(C.KSPGetOperators(ksp.p, a, pa))
  return a[]
end

# x = A \ b
function Base.A_ldiv_B!{T}(ksp::KSP{T}, b::Vec{T}, x::Vec{T})
  assemble(_ksp_A(ksp))
  assemble(b)
  assemble(x)

  chk(C.KSPSolve(ksp.p, b.p, x.p))

  reason = Ref{Cint}()
  chk(C.KSPGetConvergedReason(ksp.p, reason))
  reason[] < 0 && warn("KSP solve did not converge")

  return x
end

function Base.size(ksp::KSP)
  m = Ref{PetscInt}()
  n = Ref{PetscInt}()
  chk(C.MatGetSize(_ksp_A(ksp), m, n))
  (Int(m[]), Int(n[]))
end
Base.size{T}(ksp::KSP{T}, dim::Integer) = dim > 2 ? 1 : size(ksp)[dim]

import Base: \
(\){T}(ksp::KSP{T}, b::Vec{T}) = A_ldiv_B!(ksp, b, similar(b, size(ksp, 2)))
