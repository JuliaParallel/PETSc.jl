mutable struct IFunctionCtx{F, P, SZ, Lib}
    f::F
    p::P
    sizeu::SZ
    petsclib::Lib
    nf::Int   # cumulative count of user-RHS evaluations for `sol.stats.nf`
    err::Any  # first exception thrown by the user RHS, surfaced after TSStep
end
IFunctionCtx(f, p, sizeu, petsclib) = IFunctionCtx(f, p, sizeu, petsclib, 0, nothing)

# `t` is typed `Real` (not a concrete `Float64`) so a single body specializes
# for whichever `PetscReal` the registered `@cfunction` is built for — see
# `_petsc_ifunction_ptr`.
function _petsc_ifunction!(
    ::PETSc.LibPETSc.CTS,
    t::Real,
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

# See `_petsc_rhs_ptr`: one cached function pointer per `PetscReal` because
# `@cfunction` needs literal argument types and the time argument's ABI type
# must match the library's `PetscReal`.
const _PETSC_IFUNCTION_PTR = IdDict{DataType, Ptr{Cvoid}}()

_petsc_ifunction_ptr(::Type{Float64}) = get!(_PETSC_IFUNCTION_PTR, Float64) do
    @cfunction(
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
    )
end

_petsc_ifunction_ptr(::Type{Float32}) = get!(_PETSC_IFUNCTION_PTR, Float32) do
    @cfunction(
        _petsc_ifunction!,
        PETSc.LibPETSc.PetscErrorCode,
        (
            PETSc.LibPETSc.CTS,
            Float32,
            PETSc.LibPETSc.CVec,
            PETSc.LibPETSc.CVec,
            PETSc.LibPETSc.CVec,
            Ptr{Cvoid},
        ),
    )
end
