mutable struct RHSCtx{F, P, SZ, Lib}
    f::F
    p::P
    sizeu::SZ
    petsclib::Lib
    nf::Int   # cumulative count of user-RHS evaluations for `sol.stats.nf`
    err::Any  # first exception thrown by the user RHS, surfaced after TSStep
end
RHSCtx(f, p, sizeu, petsclib) = RHSCtx(f, p, sizeu, petsclib, 0, nothing)

# `t` is typed `Real` (not a concrete `Float64`) so a single body specializes
# for whichever `PetscReal` the registered `@cfunction` is built for — see
# `_petsc_rhs_ptr`.
function _petsc_rhs!(
    ::PETSc.LibPETSc.CTS,
    t::Real,
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

# PETSc passes the time argument across the C ABI as `PetscReal`, which is
# `Float64` or `Float32` depending on how the library was built. `@cfunction`
# requires literal argument types, so we register (and cache) a separate
# function pointer per real type, keyed by `PetscReal`. A mismatch between the
# ABI type and the registered signature would corrupt the call, so the type is
# selected from `lib.PetscReal` at registration time.
const _PETSC_RHS_PTR = IdDict{DataType, Ptr{Cvoid}}()

_petsc_rhs_ptr(::Type{Float64}) = get!(_PETSC_RHS_PTR, Float64) do
    @cfunction(
        _petsc_rhs!,
        PETSc.LibPETSc.PetscErrorCode,
        (
            PETSc.LibPETSc.CTS,
            Float64,
            PETSc.LibPETSc.CVec,
            PETSc.LibPETSc.CVec,
            Ptr{Cvoid},
        ),
    )
end

_petsc_rhs_ptr(::Type{Float32}) = get!(_PETSC_RHS_PTR, Float32) do
    @cfunction(
        _petsc_rhs!,
        PETSc.LibPETSc.PetscErrorCode,
        (
            PETSc.LibPETSc.CTS,
            Float32,
            PETSc.LibPETSc.CVec,
            PETSc.LibPETSc.CVec,
            Ptr{Cvoid},
        ),
    )
end
