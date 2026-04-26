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
