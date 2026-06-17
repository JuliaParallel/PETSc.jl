mutable struct IFunctionCtx{F, P, SZ, Lib}
    f::F
    p::P
    sizeu::SZ
    petsclib::Lib
    nf::Int   # cumulative count of user-RHS evaluations for `sol.stats.nf`
    err::Any  # first exception thrown by the user RHS, surfaced after TSStep
end
IFunctionCtx(f, p, sizeu, petsclib) = IFunctionCtx(f, p, sizeu, petsclib, 0, nothing)

function _petsc_ifunction!(
    ::PETSc.LibPETSc.CTS,
    t::Float64,
    u_ptr::PETSc.LibPETSc.CVec,
    udot_ptr::PETSc.LibPETSc.CVec,
    F_ptr::PETSc.LibPETSc.CVec,
    ctx_ptr::Ptr{Cvoid},
)::PETSc.LibPETSc.PetscErrorCode
    ctx = unsafe_pointer_to_objref(ctx_ptr)::IFunctionCtx
    # See the note in `_petsc_rhs!`: never let a Julia exception unwind across
    # the C boundary. Stash the first error and return a nonzero code so the
    # implicit solve / TSStep aborts and `step!` can rethrow the real error.
    try
        petsclib = ctx.petsclib
        u = PETSc.VecPtr(petsclib, u_ptr, false)
        udot = PETSc.VecPtr(petsclib, udot_ptr, false)
        Fv = PETSc.VecPtr(petsclib, F_ptr, false)
        PETSc.withlocalarray!(
            (u, udot, Fv);
            read = (true, true, false),
            write = (false, false, true),
        ) do u_array, udot_array, F_array
            F_reshaped = reshape(F_array, ctx.sizeu)
            u_reshaped = reshape(u_array, ctx.sizeu)
            udot_reshaped = reshape(udot_array, ctx.sizeu)
            ctx.f(F_reshaped, u_reshaped, ctx.p, t)
            @. F_reshaped = udot_reshaped - F_reshaped
        end
        ctx.nf += 1
        return PETSc.LibPETSc.PetscErrorCode(0)
    catch e
        ctx.err === nothing && (ctx.err = e)
        return _PETSC_CALLBACK_ERRCODE
    end
end

const _PETSC_IFUNCTION_PTR = Ref{Ptr{Cvoid}}(C_NULL)

function _petsc_ifunction_ptr()
    _PETSC_IFUNCTION_PTR[] == C_NULL && (_PETSC_IFUNCTION_PTR[] = @cfunction(
        _petsc_ifunction!,
        PETSc.LibPETSc.PetscErrorCode,
        (
            PETSc.LibPETSc.CTS,
            Float64,
            PETSc.LibPETSc.CVec,
            PETSc.LibPETSc.CVec,
            PETSc.LibPETSc.CVec,
            Ptr{Cvoid},
        ),
    ))
    return _PETSC_IFUNCTION_PTR[]
end
