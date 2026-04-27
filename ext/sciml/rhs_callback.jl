mutable struct RHSCtx{F, P, SZ, Lib}
    f::F
    p::P
    sizeu::SZ
    petsclib::Lib
    nf::Int   # cumulative count of user-RHS evaluations for `sol.stats.nf`
end
RHSCtx(f, p, sizeu, petsclib) = RHSCtx(f, p, sizeu, petsclib, 0)

function _petsc_rhs!(
    ::PETSc.LibPETSc.CTS,
    t::Float64,
    u_ptr::PETSc.LibPETSc.CVec,
    f_ptr::PETSc.LibPETSc.CVec,
    ctx_ptr::Ptr{Cvoid},
)::PETSc.LibPETSc.PetscErrorCode
    ctx = unsafe_pointer_to_objref(ctx_ptr)::RHSCtx
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
