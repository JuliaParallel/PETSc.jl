mutable struct RHSCtx{F, P, SZ, Lib}
    f::F
    p::P
    sizeu::SZ
    petsclib::Lib
    nf::Int   # cumulative count of user-RHS evaluations for `sol.stats.nf`
    err::Any  # first exception thrown by the user RHS, surfaced after TSStep
end
RHSCtx(f, p, sizeu, petsclib) = RHSCtx(f, p, sizeu, petsclib, 0, nothing)

function _petsc_rhs!(
    ::PETSc.LibPETSc.CTS,
    t::Float64,
    u_ptr::PETSc.LibPETSc.CVec,
    f_ptr::PETSc.LibPETSc.CVec,
    ctx_ptr::Ptr{Cvoid},
)::PETSc.LibPETSc.PetscErrorCode
    ctx = unsafe_pointer_to_objref(ctx_ptr)::RHSCtx
    # A Julia exception must never unwind through this `@cfunction` boundary
    # back into PETSc's C stack — that is undefined behaviour and typically
    # segfaults. Catch everything, stash the first error on the context, and
    # return a nonzero PETSc error code so the surrounding `TSStep` aborts.
    # `step!` rethrows the stored exception so the user sees their original
    # error instead of an opaque PETSc one. See `_take_callback_error!`.
    try
        petsclib = ctx.petsclib
        u = PETSc.VecPtr(petsclib, u_ptr, false)
        fv = PETSc.VecPtr(petsclib, f_ptr, false)
        PETSc.withlocalarray!(
            (u, fv);
            read = (true, false),
            write = (false, true),
        ) do u_array, f_array
            ctx.f(reshape(f_array, ctx.sizeu), reshape(u_array, ctx.sizeu), ctx.p, t)
        end
        ctx.nf += 1
        return PETSc.LibPETSc.PetscErrorCode(0)
    catch e
        ctx.err === nothing && (ctx.err = e)
        return _PETSC_CALLBACK_ERRCODE
    end
end

const _PETSC_RHS_PTR = Ref{Ptr{Cvoid}}(C_NULL)

function _petsc_rhs_ptr()
    _PETSC_RHS_PTR[] == C_NULL && (_PETSC_RHS_PTR[] = @cfunction(
        _petsc_rhs!,
        PETSc.LibPETSc.PetscErrorCode,
        (
            PETSc.LibPETSc.CTS,
            Float64,
            PETSc.LibPETSc.CVec,
            PETSc.LibPETSc.CVec,
            Ptr{Cvoid},
        ),
    ))
    return _PETSC_RHS_PTR[]
end
