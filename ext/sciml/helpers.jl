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
_pick_petsclib(prob, ::Nothing) = PETSc.getlib(PetscScalar = real(eltype(prob.u0)))

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
