# ============================================================================
#   Convergence of KSP, SNES and TS solves
# ============================================================================

"""
    ConvergenceError <: Exception

Thrown by [`solve!`](@ref) and [`step!`](@ref) when the solver stopped without
converging, unless they are called with `check = false`. The solution vector is
left as PETSc left it.

# Fields
- `solver`: the `KSP`, `SNES` or `TS` that was solved
- `reason`: its negative converged reason (`KSPConvergedReason`,
  `SNESConvergedReason` or `TSConvergedReason`), as [`converged_reason`](@ref) reports it
- `iterations`: the iteration count, or the step number for a `TS`
- `function_norm`: the last residual norm for a `SNES`, `nothing` otherwise
"""
struct ConvergenceError{S, R} <: Exception
    solver::S
    reason::R
    iterations::Int
    function_norm::Union{Nothing, Real}
end

function Base.showerror(io::IO, e::ConvergenceError)
    name = nameof(typeof(e.solver))
    count = e.solver isa AbstractTS ? "at step $(e.iterations)" : "after $(e.iterations) iterations"
    print(io, "ConvergenceError: $name did not converge: $(e.reason), ", reason_text(e.reason), ", ", count)
    isnothing(e.function_norm) || print(io, " (function norm ", e.function_norm, ")")
    print(io, ".\n", convergence_hint(e.solver))
end

const KSP_REASON_TEXT = Dict(
    LibPETSc.KSP_DIVERGED_NULL => "the solver hit a null vector",
    LibPETSc.KSP_DIVERGED_ITS => "it reached the maximum number of iterations (-ksp_max_it)",
    LibPETSc.KSP_DIVERGED_DTOL => "the residual grew past the divergence tolerance (-ksp_divtol)",
    LibPETSc.KSP_DIVERGED_BREAKDOWN => "the Krylov method broke down",
    LibPETSc.KSP_DIVERGED_BREAKDOWN_BICG => "the BiCG method broke down",
    LibPETSc.KSP_DIVERGED_NONSYMMETRIC => "the method needs a symmetric operator or preconditioner",
    LibPETSc.KSP_DIVERGED_INDEFINITE_PC => "the method needs a positive definite preconditioner",
    LibPETSc.KSP_DIVERGED_NANORINF => "the residual became NaN or Inf",
    LibPETSc.KSP_DIVERGED_INDEFINITE_MAT => "the method needs a positive definite operator",
    LibPETSc.KSP_DIVERGED_PC_FAILED => "the preconditioner failed, for example a zero pivot in a factorization",
)

const SNES_REASON_TEXT = Dict(
    LibPETSc.SNES_DIVERGED_FUNCTION_DOMAIN => "the residual was evaluated outside its domain",
    LibPETSc.SNES_DIVERGED_FUNCTION_COUNT => "it reached the maximum number of residual evaluations (-snes_max_funcs)",
    LibPETSc.SNES_DIVERGED_LINEAR_SOLVE => "the linear solver failed",
    LibPETSc.SNES_DIVERGED_FUNCTION_NANORINF => "the residual norm became NaN or Inf",
    LibPETSc.SNES_DIVERGED_MAX_IT => "it reached the maximum number of iterations (-snes_max_it)",
    LibPETSc.SNES_DIVERGED_LINE_SEARCH => "the line search failed",
    LibPETSc.SNES_DIVERGED_INNER => "an inner solver failed",
    LibPETSc.SNES_DIVERGED_LOCAL_MIN => "it reached a local minimum of the residual norm",
    LibPETSc.SNES_DIVERGED_DTOL => "the residual grew past the divergence tolerance (-snes_divergence_tolerance)",
    LibPETSc.SNES_DIVERGED_JACOBIAN_DOMAIN => "the Jacobian was evaluated outside its domain",
    LibPETSc.SNES_DIVERGED_TR_DELTA => "the trust region became too small",
    LibPETSc.SNES_DIVERGED_USER => "a convergence test stopped it",
    LibPETSc.SNES_DIVERGED_OBJECTIVE_DOMAIN => "the objective was evaluated outside its domain",
    LibPETSc.SNES_DIVERGED_OBJECTIVE_NANORINF => "the objective became NaN or Inf",
)

const TS_REASON_TEXT = Dict(
    LibPETSc.TS_DIVERGED_NONLINEAR_SOLVE => "the nonlinear solve of a step failed more often than allowed (-ts_max_snes_failures)",
    LibPETSc.TS_DIVERGED_STEP_REJECTED => "the step was rejected more often than allowed (-ts_max_reject)",
)

reason_text(r::LibPETSc.KSPConvergedReason) = get(KSP_REASON_TEXT, r, "see the PETSc manual page for KSPConvergedReason")
reason_text(r::LibPETSc.SNESConvergedReason) = get(SNES_REASON_TEXT, r, "see the PETSc manual page for SNESConvergedReason")
reason_text(r::LibPETSc.TSConvergedReason) = get(TS_REASON_TEXT, r, "see the PETSc manual page for TSConvergedReason")

convergence_hint(::AbstractKSP) =
    "The options -ksp_converged_reason and -ksp_monitor_true_residual show how the residual evolved."
convergence_hint(::AbstractSNES) =
    "The options -snes_converged_reason, -snes_monitor and -ksp_converged_reason show where it stopped."
convergence_hint(::AbstractTS) =
    "The options -ts_monitor and -snes_converged_reason show which step failed and why."

convergence_error(ksp::AbstractKSP, reason) =
    ConvergenceError(ksp, reason, Int(iteration_number(ksp)), nothing)
convergence_error(snes::AbstractSNES, reason) =
    ConvergenceError(snes, reason, Int(iteration_number(snes)), function_norm(snes))
convergence_error(ts::AbstractTS, reason) =
    ConvergenceError(ts, reason, Int(step_number(ts)), nothing)

# Runs the solve `f` on `solver` and then acts on how it ended: `check = true` throws a
# ConvergenceError, `false` does nothing.
function with_convergence_check(f, solver, check::Bool)
    try
        f()
    catch err
        # A solver set to error out when it fails (a TS by default, or
        # -ksp_error_if_not_converged) raises a PetscError; `check = true` turns it
        # into a ConvergenceError, `check = false` lets it stand.
        err isa LibPETSc.PetscError && check || rethrow()
        reason = converged_reason(solver)
        Int(reason) < 0 || rethrow()
        throw(convergence_error(solver, reason))
    end
    check || return nothing
    reason = converged_reason(solver)
    Int(reason) < 0 && throw(convergence_error(solver, reason))
    return nothing
end
