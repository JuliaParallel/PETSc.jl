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

function _petsc_retcode(petsclib, ts)
    reason = _ts_converged_reason(petsclib, ts)
    if reason == PETSc.LibPETSc.TS_CONVERGED_TIME ||
       reason == PETSc.LibPETSc.TS_CONVERGED_USER ||
       reason == PETSc.LibPETSc.TS_CONVERGED_EVENT
        return SciMLBase.ReturnCode.Success
    elseif reason == PETSc.LibPETSc.TS_CONVERGED_ITS
        return SciMLBase.ReturnCode.MaxIters
    elseif reason == PETSc.LibPETSc.TS_CONVERGED_ITERATING
        return SciMLBase.ReturnCode.Failure
    elseif reason == PETSc.LibPETSc.TS_DIVERGED_NONLINEAR_SOLVE
        return SciMLBase.ReturnCode.Failure
    elseif reason == PETSc.LibPETSc.TS_DIVERGED_STEP_REJECTED
        return SciMLBase.ReturnCode.Unstable
    else
        return SciMLBase.ReturnCode.Failure
    end
end
