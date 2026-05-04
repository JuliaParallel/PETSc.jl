function _check_isinplace(prob)
    SciMLBase.isinplace(prob) || throw(ArgumentError(
        "PETSc.jl time-stepping wrappers only support in-place ODEProblems " *
        "(f!(du, u, p, t)). Wrap your function in an in-place form or use a " *
        "different solver.",
    ))
end

# The PETSc TS C callback signatures we register via `@cfunction` use
# `Float64` for the time argument. PETSc libraries built with
# `PetscReal = Float32` would pass a Float32 across the C ABI, which would
# not match. Detect this up-front so users get a clear error rather than a
# bus error when the callback fires.
function _check_petscreal(lib)
    lib.PetscReal === Float64 || throw(ArgumentError(
        "PETSc.jl SciML extension currently only supports PetscReal = Float64. " *
        "Got PetscReal = $(lib.PetscReal). Pass `petsclib = PETSc.getlib(PetscScalar = Float64)` " *
        "or use a PETSc build with PetscReal = Float64.",
    ))
end

_pick_petsclib(prob, petsclib) = petsclib
function _pick_petsclib(prob, ::Nothing)
    T = eltype(prob.u0)
    T <: Complex && throw(ArgumentError(
        "PETSc.jl SciML extension currently only supports real-valued ODE " *
        "problems (got eltype(u0) = $T). Pass a real `u0`, or — if you have " *
        "a PETSc build with a complex `PetscScalar` and the matching support " *
        "in this extension lands — supply `petsclib` explicitly.",
    ))
    return PETSc.getlib(PetscScalar = T)
end

function _setfromoptions!(petsclib, ts, petsc_options::AbstractVector{<:AbstractString})
    isempty(petsc_options) && return nothing
    opts = PETSc.Options(petsclib; PETSc.parse_options(String.(petsc_options))...)
    push!(opts)
    try
        PETSc.LibPETSc.TSSetFromOptions(petsclib, ts)
    finally
        pop!(opts)
        PETSc.destroy(opts)
    end
    return nothing
end
_setfromoptions!(petsclib, ts, ::Nothing) = nothing

function _sync_petsc_to_julia!(integ::PETScTSIntegrator)
    PETSc.withlocalarray!(integ.u_petsc; read = true, write = false) do arr
        copyto!(integ.u, reshape(arr, integ.sizeu))
    end
    return nothing
end

function _sync_julia_to_petsc!(integ::PETScTSIntegrator)
    PETSc.withlocalarray!(integ.u_petsc; read = false, write = true) do arr
        copyto!(arr, vec(integ.u))
    end
    return nothing
end

# Reject SciML keywords that this extension does not actually honour. Letting
# the open-ended `kwargs...` sink swallow standard knobs like `adaptive` /
# `dtmin` / `progress` would silently break the usual SciML solver contract,
# so any unsupported key fails loudly with a clear, named error.
const _SUPPORTED_SCIML_KWARGS = (
    :save_everystep, :save_on, :save_start, :save_end, :save_discretes,
    :saveat, :tstops, :callback, :initialize_save,
    :reltol, :abstol,
    :dt, :dtmin, :dtmax, :adaptive,
    :maxiters, :petsclib,
    :verbose,
)

function _reject_unsupported_kwargs(kwargs)
    for key in keys(kwargs)
        key in _SUPPORTED_SCIML_KWARGS && continue
        throw(ArgumentError(
            "PETSc.jl SciML extension: keyword argument `$(key)` is not " *
            "supported. Supported keywords are: " *
            join(_SUPPORTED_SCIML_KWARGS, ", ") * ". " *
            "Pass equivalent PETSc CLI flags via the algorithm's " *
            "`petsc_options` argument if a matching SciML knob is missing.",
        ))
    end
    return nothing
end
