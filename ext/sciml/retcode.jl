@generated function _ts_converged_reason(
    petsclib::PETSc.LibPETSc.PetscLibType{ST, IT, LT},
    ts,
) where {ST, IT, LT}
    libsym = PETSc.LibPETSc.petsclibs[findfirst(
        l -> typeof(l) == PETSc.LibPETSc.PetscLibType{ST, IT, LT},
        PETSc.LibPETSc.petsclibs,
    )].petsc_library
    quote
        reason_ref = Ref(PETSc.LibPETSc.TS_CONVERGED_ITERATING)
        err = ccall(
            (:TSGetConvergedReason, $libsym),
            PETSc.LibPETSc.PetscErrorCode,
            (PETSc.LibPETSc.CTS, Ptr{PETSc.LibPETSc.TSConvergedReason}),
            ts,
            reason_ref,
        )
        iszero(err) || error("TSGetConvergedReason failed with code $err")
        return reason_ref[]
    end
end

# `TSSetTolerances` accepts NULL vectors when only scalar tolerances are
# desired, but the autowrapped Julia signature requires a `PetscVec` for the
# vector arguments. ccall directly so we can pass `C_NULL` for the per-
# component tolerances.
@generated function _ts_set_scalar_tolerances!(
    petsclib::PETSc.LibPETSc.PetscLibType{ST, IT, LT},
    ts,
    abstol,
    reltol,
) where {ST, IT, LT}
    libsym = PETSc.LibPETSc.petsclibs[findfirst(
        l -> typeof(l) == PETSc.LibPETSc.PetscLibType{ST, IT, LT},
        PETSc.LibPETSc.petsclibs,
    )].petsc_library
    PetscReal = ST <: Complex ? real(ST) : ST
    quote
        err = ccall(
            (:TSSetTolerances, $libsym),
            PETSc.LibPETSc.PetscErrorCode,
            (PETSc.LibPETSc.CTS, $PetscReal, PETSc.LibPETSc.CVec,
             $PetscReal, PETSc.LibPETSc.CVec),
            ts, $PetscReal(abstol), C_NULL, $PetscReal(reltol), C_NULL,
        )
        iszero(err) || error("TSSetTolerances failed with code $err")
        return nothing
    end
end

# Translate PETSc's `TSConvergedReason` into the SciML `ReturnCode` the solution
# is tagged with. The mapping is not one-to-one, so the non-obvious cases are
# spelled out below.
function _petsc_retcode(petsclib, ts)
    reason = _ts_converged_reason(petsclib, ts)
    if reason == PETSc.LibPETSc.TS_CONVERGED_TIME ||
       reason == PETSc.LibPETSc.TS_CONVERGED_USER ||
       reason == PETSc.LibPETSc.TS_CONVERGED_EVENT
        # Reached the final time, or a user/event callback asked to stop on
        # purpose — all of these are a successful, intended termination.
        return SciMLBase.ReturnCode.Success
    elseif reason == PETSc.LibPETSc.TS_CONVERGED_ITS
        # Hit PETSc's own max-step limit before reaching `tf`. SciML's closest
        # equivalent is `MaxIters` (the step budget, not the solution, ran out).
        return SciMLBase.ReturnCode.MaxIters
    elseif reason == PETSc.LibPETSc.TS_CONVERGED_ITERATING
        # Still "iterating" means the integrator stopped while mid-integration
        # without reaching `tf` and without any converged/diverged verdict —
        # i.e. it never actually finished, so report a generic failure.
        return SciMLBase.ReturnCode.Failure
    elseif reason == PETSc.LibPETSc.TS_DIVERGED_NONLINEAR_SOLVE
        # The implicit (SNES) solve diverged: this is a solver failure, not a
        # statement about the ODE's dynamics, so map to the generic `Failure`.
        return SciMLBase.ReturnCode.Failure
    elseif reason == PETSc.LibPETSc.TS_DIVERGED_STEP_REJECTED
        # The adaptive controller kept rejecting steps (it could not satisfy the
        # error tolerance even at the minimum step). That is the classic
        # signature of a stiff/blowing-up trajectory, which SciML spells
        # `Unstable` rather than a bare `Failure`.
        return SciMLBase.ReturnCode.Unstable
    else
        # Any other (or future) diverged reason: be conservative and report a
        # generic failure rather than silently claiming success.
        return SciMLBase.ReturnCode.Failure
    end
end
