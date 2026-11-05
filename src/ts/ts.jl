function LibPETSc.TSSetRHSFunction(
    petsclib::LibPETSc.PetscLibType,
    ts::LibPETSc.AbstractTS,
    ::Nothing,
    fptr::Ptr{Cvoid},
    ctx::Ptr{Cvoid} = C_NULL,
) end

LibPETSc.@for_petsc function LibPETSc.TSSetRHSFunction(
    petsclib::$UnionPetscLib,
    ts::LibPETSc.AbstractTS,
    ::Nothing,
    fptr::Ptr{Cvoid},
    ctx::Ptr{Cvoid} = C_NULL,
)
    LibPETSc.@chk ccall(
        (:TSSetRHSFunction, $petsc_library),
        LibPETSc.PetscErrorCode,
        (LibPETSc.CTS, LibPETSc.CVec, Ptr{Cvoid}, Ptr{Cvoid}),
        ts,
        C_NULL,
        fptr,
        ctx,
    )
    return nothing
end

function LibPETSc.TSSetIFunction(
    petsclib::LibPETSc.PetscLibType,
    ts::LibPETSc.AbstractTS,
    ::Nothing,
    fptr::Ptr{Cvoid},
    ctx::Ptr{Cvoid} = C_NULL,
) end

LibPETSc.@for_petsc function LibPETSc.TSSetIFunction(
    petsclib::$UnionPetscLib,
    ts::LibPETSc.AbstractTS,
    ::Nothing,
    fptr::Ptr{Cvoid},
    ctx::Ptr{Cvoid} = C_NULL,
)
    LibPETSc.@chk ccall(
        (:TSSetIFunction, $petsc_library),
        LibPETSc.PetscErrorCode,
        (LibPETSc.CTS, LibPETSc.CVec, Ptr{Cvoid}, Ptr{Cvoid}),
        ts,
        C_NULL,
        fptr,
        ctx,
    )
    return nothing
end


"""
    TSMonitorSet(petsclib, ts, monitor::Ptr{Cvoid}, ctx = C_NULL, mdestroy = C_NULL)

Defaults for the context and destroy-callback arguments of the generated `TSMonitorSet`.
"""
LibPETSc.TSMonitorSet(petsclib::LibPETSc.UnionPetscLibType, ts::LibPETSc.AbstractTS, monitor::Ptr{Cvoid}) =
    LibPETSc.TSMonitorSet(petsclib, ts, monitor, C_NULL, C_NULL)
LibPETSc.TSMonitorSet(petsclib::LibPETSc.UnionPetscLibType, ts::LibPETSc.AbstractTS, monitor::Ptr{Cvoid}, ctx::Ptr{Cvoid}) =
    LibPETSc.TSMonitorSet(petsclib, ts, monitor, ctx, C_NULL)

"""
    TSARKIMEXRegister(
        petsclib,
        name::String,
        order,
        s,
        At,
        bt,
        ct,
        A,
        b,
        c,
        bembedt,
        bembed,
        pinterp,
        binterpt,
        binterp,
    )

Julia-friendly overload for registering a custom `TSARKIMEX` tableau. Optional
PETSc arrays may be passed as `nothing`, which is translated to `NULL`.

For the stage tables `At` and `A`, PETSc expects flat vectors in row-major
order, matching the layout used by C arrays. If you start from a Julia matrix,
do not pass `vec(A)` directly since Julia stores matrices column-major; flatten
row-by-row instead, for example with `vec(permutedims(A))`.
"""
# (the docstring stub is the generated one in TS_wrappers.jl; these methods add `nothing` for optional arrays)
LibPETSc.@for_petsc function LibPETSc.TSARKIMEXRegister(
    petsclib::$UnionPetscLib,
    name::String,
    order::Integer,
    s::Integer,
    At::AbstractVector,
    bt::Union{Nothing, AbstractVector},
    ct::Union{Nothing, AbstractVector},
    A::AbstractVector,
    b::Union{Nothing, AbstractVector},
    c::Union{Nothing, AbstractVector},
    bembedt::Union{Nothing, AbstractVector},
    bembed::Union{Nothing, AbstractVector},
    pinterp::Integer,
    binterpt::Union{Nothing, AbstractVector},
    binterp::Union{Nothing, AbstractVector},
)
    At_vals = $PetscReal.(At)
    A_vals = $PetscReal.(A)
    bt_vals = bt === nothing ? Ptr{$PetscReal}(C_NULL) : $PetscReal.(bt)
    ct_vals = ct === nothing ? Ptr{$PetscReal}(C_NULL) : $PetscReal.(ct)
    b_vals = b === nothing ? Ptr{$PetscReal}(C_NULL) : $PetscReal.(b)
    c_vals = c === nothing ? Ptr{$PetscReal}(C_NULL) : $PetscReal.(c)
    bembedt_vals =
        bembedt === nothing ? Ptr{$PetscReal}(C_NULL) : $PetscReal.(bembedt)
    bembed_vals =
        bembed === nothing ? Ptr{$PetscReal}(C_NULL) : $PetscReal.(bembed)
    binterpt_vals =
        binterpt === nothing ? Ptr{$PetscReal}(C_NULL) : $PetscReal.(binterpt)
    binterp_vals =
        binterp === nothing ? Ptr{$PetscReal}(C_NULL) : $PetscReal.(binterp)

    LibPETSc.@chk ccall(
        (:TSARKIMEXRegister, $petsc_library),
        LibPETSc.PetscErrorCode,
        (
            LibPETSc.TSARKIMEXType,
            $PetscInt,
            $PetscInt,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            Ptr{$PetscReal},
            $PetscInt,
            Ptr{$PetscReal},
            Ptr{$PetscReal},
        ),
        name,
        $PetscInt(order),
        $PetscInt(s),
        At_vals,
        bt_vals,
        ct_vals,
        A_vals,
        b_vals,
        c_vals,
        bembedt_vals,
        bembed_vals,
        $PetscInt(pinterp),
        binterpt_vals,
        binterp_vals,
    )
    return nothing
end

# ============================================================================
#   High-level TS interface
# ============================================================================

# The Julia side of a TS, kept with the PETSc object
mutable struct TSState <: ObjectState
    rhs_function!::Any
    rhs_jacobian!::Any
    ifunction!::Any
    ijacobian!::Any
    monitor::Any
    pre_step!::Any
    post_step!::Any
    opts::Any
    options_at::Symbol
    options_applied::Bool
    alive::Bool
end
TSState() = TSState(
    nothing,
    nothing,
    nothing,
    nothing,
    nothing,
    nothing,
    nothing,
    nothing,
    :never,
    false,
    true,
)
state_type(::Type{<:TS}) = TSState

# `ts.opts` and the callbacks reach the state
@inline Base.getproperty(ts::TS, name::Symbol) =
    (name === :ptr || name === :age || name === :own) ? getfield(ts, name) :
    forward_getproperty(ts, name)
@inline function Base.setproperty!(ts::TS, name::Symbol, value)
    (name === :ptr || name === :age || name === :own) &&
        return Base.setfield!(ts, name, convert(fieldtype(typeof(ts), name), value))
    return forward_setproperty!(ts, name, value)
end
Base.propertynames(ts::TS, private::Bool = false) = forward_propertynames(typeof(ts))


function Base.show(io::IO, ts::AbstractTS{PetscLib}) where {PetscLib}
    if ts.ptr == C_NULL
        print(io, "PETSc TS (null pointer)")
    else
        print(io, "PETSc TS object")
    end
    return nothing
end

"""
    TS(petsclib, comm::MPI.Comm; prefix = "", options_at = :first_solve, options...)

Create a PETSc time stepper on the communicator `comm`.

`options_at` says when the options are applied, through `TSSetFromOptions`: at the first
`setup!`, `solve!` or `step!` (`:first_solve`, the default), at every `solve!`
(`:each_solve`) or in the constructor (`:construction`). PETSc's options database
wins on a key both set. See [When the options apply](@ref) for the details.

`exact_final_time` defaults to `TS_EXACTFINALTIME_MATCHSTEP`, so the last step
lands on the time set by [`set_max_time!`](@ref). 
PETSc's own default is `TS_EXACTFINALTIME_UNSPECIFIED`, which integrates 
to the wrong time without reporting an error.

If `comm` has size 1 the garbage collector calls [`destroy!`](@ref).
Otherwise destruction is the caller's responsibility.

# External Links
$(doc_external("TS/TSCreate"))
$(doc_external("TS/TSSetExactFinalTime"))
"""
function TS(
    petsclib::PetscLib,
    comm::MPI.Comm;
    prefix::AbstractString = "",
    exact_final_time::LibPETSc.TSExactFinalTimeOption = LibPETSc.TS_EXACTFINALTIME_MATCHSTEP,
    options_at::Symbol = :first_solve,
    options...,
) where {PetscLib}
    check_initialized(getlib(PetscLib))

    petsclib = getlib(PetscLib)
    return own_solver!(LibPETSc.TSCreate(petsclib, comm), comm) do ts
        isempty(prefix) || LibPETSc.TSSetOptionsPrefix(petsclib, ts, String(prefix))
        LibPETSc.TSSetExactFinalTime(petsclib, ts, exact_final_time)
        store_options!(ts, options_at, options)
    end
end

"""
    set_exact_final_time!(ts::AbstractTS, option)

Choose how the final step meets the time set by [`set_max_time!`](@ref):
`TS_EXACTFINALTIME_MATCHSTEP`, `TS_EXACTFINALTIME_INTERPOLATE` or `TS_EXACTFINALTIME_STEPOVER`.

# External Links
$(doc_external("TS/TSSetExactFinalTime"))
"""
function set_exact_final_time!(
    ts::AbstractTS{PetscLib},
    option::LibPETSc.TSExactFinalTimeOption,
) where {PetscLib}
    LibPETSc.TSSetExactFinalTime(getlib(PetscLib), ts, option)
    return ts
end

"""
    set_adapt_type!(ts::AbstractTS, type::Symbol)

Set the timestep adaptivity controller, for example `:none` to hold the step size fixed, 
or `:basic` for the default error-based controller.

# External Links
$(doc_external("TS/TSAdaptSetType"))
"""
function set_adapt_type!(ts::AbstractTS{PetscLib}, type::Symbol) where {PetscLib}
    petsclib = getlib(PetscLib)
    LibPETSc.TSAdaptSetType(
        petsclib,
        LibPETSc.TSGetAdapt(petsclib, ts),
        String(type),
    )
    return ts
end

"""
    destroy!(ts::AbstractTS)

Destroy `ts` and release the options database attached to it.

The call is a no-op when the library has been finalized or when `ts` predates
the current initialize/finalize cycle, so a stale handle never reaches `TSDestroy`.
Does nothing on a borrowed handle: see [`owns`](@ref).

# External Links
$(doc_external("TS/TSDestroy"))
"""
function destroy!(ts::AbstractTS{PetscLib}) where {PetscLib}
    owns(ts) && destroy_options!(ts)
    return destroy_handle!(LibPETSc.TSDestroy, ts, PetscLib)
end

"""
    type_name(obj)

The name of the PETSc implementation `obj` currently uses.

The reader is `type_name` rather than `type`, because `type` is the conventional name
for the argument of the matching setter and would be shadowed by it, and because what
comes back is the name of a PETSc implementation rather than a Julia type.

Declared here, with methods on `TS`, `KSP`, `SNES`, `PC`, `PetscVec`, `PetscMat` and
`AbstractPetscDM`.
Every method answers with a `Symbol`, or `nothing` while no type has been set.

# External Links
$(doc_external("TS/TSGetType"))
$(doc_external("KSP/KSPGetType"))
$(doc_external("DM/DMGetType"))
"""
function type_name end

# `type_name(ts)` answers with a `Symbol`, or `nothing` when no type has been set yet:
# PETSc reports an unset type as a null string, which the generated `TSGetType` cannot
# convert, and asking a freshly created stepper for its type is reasonable.
LibPETSc.@for_petsc function type_name(ts::AbstractTS{$PetscLib})
    r_type = Ref{Ptr{Cchar}}(C_NULL)
    LibPETSc.@chk ccall(
        (:TSGetType, $petsc_library),
        LibPETSc.PetscErrorCode,
        (LibPETSc.CTS, Ptr{Ptr{Cchar}}),
        ts,
        r_type,
    )
    return r_type[] == C_NULL ? nothing : Symbol(unsafe_string(r_type[]))
end

"""
    set_type!(ts::AbstractTS, type::Symbol)

Set the time-stepping method, for example `:bdf`, `:rk` or `:arkimex`.

# External Links
$(doc_external("TS/TSSetType"))
"""
function set_type!(ts::AbstractTS{PetscLib}, type::Symbol) where {PetscLib}
    LibPETSc.TSSetType(getlib(PetscLib), ts, String(type))
    return ts
end

"""
    set_problem_type!(ts::AbstractTS, type)

Declare the problem as `LibPETSc.TS_LINEAR` or `LibPETSc.TS_NONLINEAR`.

# External Links
$(doc_external("TS/TSSetProblemType"))
"""
function set_problem_type!(
    ts::AbstractTS{PetscLib},
    type::LibPETSc.TSProblemType,
) where {PetscLib}
    LibPETSc.TSSetProblemType(getlib(PetscLib), ts, type)
    return ts
end

"""
    equation_type(ts::AbstractTS)

The form of the equations `ts` integrates, a `LibPETSc.TSEquationType`
(`TS_EQ_UNSPECIFIED` until [`set_equation_type!`](@ref) is called).

# External Links
$(doc_external("TS/TSGetEquationType"))
"""
equation_type(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetEquationType(getlib(PetscLib), ts)

"""
    set_equation_type!(ts::AbstractTS, type)

Declare the form of the equations, for example
`LibPETSc.TS_EQ_DAE_IMPLICIT_INDEX1` for an index-1 DAE written as
``F(t, u, du/dt) = 0``. Some integrators read it to treat the algebraic
components correctly. Returns `ts`.

# External Links
$(doc_external("TS/TSSetEquationType"))
"""
function set_equation_type!(
    ts::AbstractTS{PetscLib},
    type::LibPETSc.TSEquationType,
) where {PetscLib}
    LibPETSc.TSSetEquationType(getlib(PetscLib), ts, type)
    return ts
end

"""
    dm(ts::AbstractTS)

The DM attached to `ts`, [`narrow`](@ref)ed to its flavour.

$(doc_borrowed())

# External Links
$(doc_external("TS/TSGetDM"))
"""
dm(ts::AbstractTS{PetscLib}) where {PetscLib} =
    narrow(LibPETSc.TSGetDM(getlib(PetscLib), ts))

"""
    set_dm!(ts::AbstractTS, dm::AbstractPetscDM)

Attach `dm` to `ts`.

# External Links
$(doc_external("TS/TSSetDM"))
"""
function set_dm!(
    ts::AbstractTS{PetscLib},
    dm::AbstractPetscDM{PetscLib},
) where {PetscLib}
    LibPETSc.TSSetDM(getlib(PetscLib), ts, dm)
    return ts
end

"""
    solution(ts::AbstractTS)

The solution vector held by `ts`, as a `PetscVec`.

$(doc_borrowed())

# External Links
$(doc_external("TS/TSGetSolution"))
"""
solution(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetSolution(getlib(PetscLib), ts)

"""
    set_solution!(ts::AbstractTS, u::AbstractPetscVec)

Set the initial condition of `ts`.

# External Links
$(doc_external("TS/TSSetSolution"))
"""
function set_solution!(
    ts::AbstractTS{PetscLib},
    u::AbstractPetscVec{PetscLib},
) where {PetscLib}
    LibPETSc.TSSetSolution(getlib(PetscLib), ts, u)
    return ts
end

# Time and step controls
# ----------------------------------------------------------------------------

"""
    current_time(ts::AbstractTS)

The time `ts` has reached.

# External Links
$(doc_external("TS/TSGetTime"))
"""
current_time(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetTime(getlib(PetscLib), ts)

"""
    prev_time(ts::AbstractTS)

The time at the start of the last step `ts` took.

# External Links
$(doc_external("TS/TSGetPrevTime"))
"""
prev_time(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetPrevTime(getlib(PetscLib), ts)

"""
    set_time!(ts::AbstractTS, t)

Set the current time of `ts`.

# External Links
$(doc_external("TS/TSSetTime"))
"""
function set_time!(ts::AbstractTS{PetscLib}, t) where {PetscLib}
    LibPETSc.TSSetTime(getlib(PetscLib), ts, PetscLib.PetscReal(t))
    return ts
end

"""
    timestep(ts::AbstractTS)

The current step size.

# External Links
$(doc_external("TS/TSGetTimeStep"))
"""
timestep(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetTimeStep(getlib(PetscLib), ts)

"""
    set_timestep!(ts::AbstractTS, dt)

Set the step size.

# External Links
$(doc_external("TS/TSSetTimeStep"))
"""
function set_timestep!(ts::AbstractTS{PetscLib}, dt) where {PetscLib}
    LibPETSc.TSSetTimeStep(getlib(PetscLib), ts, PetscLib.PetscReal(dt))
    return ts
end

"""
    max_time(ts::AbstractTS)

The time at which integration stops.

# External Links
$(doc_external("TS/TSGetMaxTime"))
"""
max_time(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetMaxTime(getlib(PetscLib), ts)

"""
    set_max_time!(ts::AbstractTS, t)

Set the time at which integration stops.

# External Links
$(doc_external("TS/TSSetMaxTime"))
"""
function set_max_time!(ts::AbstractTS{PetscLib}, t) where {PetscLib}
    LibPETSc.TSSetMaxTime(getlib(PetscLib), ts, PetscLib.PetscReal(t))
    return ts
end

"""
    max_steps(ts::AbstractTS)

The step count at which integration stops.

# External Links
$(doc_external("TS/TSGetMaxSteps"))
"""
max_steps(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetMaxSteps(getlib(PetscLib), ts)

"""
    set_max_steps!(ts::AbstractTS, n)

Set the step count at which integration stops.

# External Links
$(doc_external("TS/TSSetMaxSteps"))
"""
function set_max_steps!(ts::AbstractTS{PetscLib}, n) where {PetscLib}
    LibPETSc.TSSetMaxSteps(getlib(PetscLib), ts, PetscLib.PetscInt(n))
    return ts
end

"""
    set_max_snes_failures!(ts::AbstractTS, n)

Allow `n` failed nonlinear solves over the whole run before [`solve!`](@ref)
stops with `TS_DIVERGED_NONLINEAR_SOLVE` in [`converged_reason`](@ref). A failed
solve makes the adaptor retry the step with a smaller time step. The default
is 1; `LibPETSc.PETSC_UNLIMITED` removes the limit. Returns `ts`.

A solve stopped by [`set_function_domain_error!`](@ref) does not count here:
it rejects the step, which [`set_max_step_rejections!`](@ref) limits.

# External Links
$(doc_external("TS/TSSetMaxSNESFailures"))
"""
function set_max_snes_failures!(ts::AbstractTS{PetscLib}, n::Integer) where {PetscLib}
    LibPETSc.TSSetMaxSNESFailures(getlib(PetscLib), ts, PetscLib.PetscInt(n))
    return ts
end

"""
    set_max_step_rejections!(ts::AbstractTS, n)

Allow the adaptor to reject `n` attempts at a single step before
[`solve!`](@ref) stops with `TS_DIVERGED_STEP_REJECTED`. The default is 10;
`LibPETSc.PETSC_UNLIMITED` removes the limit. Returns `ts`.

# External Links
$(doc_external("TS/TSSetMaxStepRejections"))
"""
function set_max_step_rejections!(ts::AbstractTS{PetscLib}, n::Integer) where {PetscLib}
    LibPETSc.TSSetMaxStepRejections(getlib(PetscLib), ts, PetscLib.PetscInt(n))
    return ts
end

"""
    set_error_if_step_fails!(ts::AbstractTS, flag::Bool)

With `flag = true`, PETSc's default, a step that fails for good makes
[`solve!`](@ref) and [`step!`](@ref) throw a `PetscError`. With `false` they
return normally and [`converged_reason`](@ref) says why the run stopped, so
check it. Returns `ts`.

# External Links
$(doc_external("TS/TSSetErrorIfStepFails"))
"""
function set_error_if_step_fails!(ts::AbstractTS{PetscLib}, flag::Bool) where {PetscLib}
    LibPETSc.TSSetErrorIfStepFails(getlib(PetscLib), ts, LibPETSc.PetscBool(flag))
    return ts
end

"""
    step_number(ts::AbstractTS)

The number of steps taken so far.

# External Links
$(doc_external("TS/TSGetStepNumber"))
"""
step_number(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetStepNumber(getlib(PetscLib), ts)

"""
    set_step_number!(ts::AbstractTS, n)

Set the number of steps taken so far, for example to continue a run from a
checkpoint so that monitors and [`set_max_steps!`](@ref) count on from there.
Set the time with [`set_time!`](@ref) alongside it. Returns `ts`.

# External Links
$(doc_external("TS/TSSetStepNumber"))
"""
function set_step_number!(ts::AbstractTS{PetscLib}, n::Integer) where {PetscLib}
    LibPETSc.TSSetStepNumber(getlib(PetscLib), ts, PetscLib.PetscInt(n))
    return ts
end

"""
    tolerances(ts::AbstractTS)

Local truncation error tolerances, as `(; atol, rtol, vatol, vrtol)`.

`vatol` and `vrtol` are the per-component tolerance vectors, as `PetscVec`s
holding a null pointer when only the scalar tolerances were set.

$(doc_borrowed())

# External Links
$(doc_external("TS/TSGetTolerances"))
"""
function tolerances(ts::AbstractTS{PetscLib}) where {PetscLib}
    petsclib = getlib(PetscLib)
    atol, vatol, rtol, vrtol = LibPETSc.TSGetTolerances(petsclib, ts)
    return (; atol, rtol, vatol, vrtol)
end

"""
    set_tolerances!(ts::AbstractTS; atol, rtol, vatol, vrtol)

Set the local truncation error tolerances. A keyword left at `nothing` keeps
the value `ts` currently has.

Pass `vatol` or `vrtol` to give per-component tolerances.

# External Links
$(doc_external("TS/TSSetTolerances"))
"""
function set_tolerances!(
    ts::AbstractTS{PetscLib};
    atol = nothing,
    rtol = nothing,
    vatol::Union{Nothing, AbstractPetscVec{PetscLib}} = nothing,
    vrtol::Union{Nothing, AbstractPetscVec{PetscLib}} = nothing,
) where {PetscLib}
    petsclib = getlib(PetscLib)
    PetscReal = PetscLib.PetscReal
    cur_atol, _, cur_rtol, _ = LibPETSc.TSGetTolerances(petsclib, ts)
    null_vec = LibPETSc.PetscVec(petsclib)
    LibPETSc.TSSetTolerances(
        petsclib,
        ts,
        isnothing(atol) ? cur_atol : PetscReal(atol),
        isnothing(vatol) ? null_vec : vatol,
        isnothing(rtol) ? cur_rtol : PetscReal(rtol),
        isnothing(vrtol) ? null_vec : vrtol,
    )
    return ts
end

"""
    converged_reason(ts::AbstractTS)

Why the integration stopped, as a `LibPETSc.TSConvergedReason`.

# External Links
$(doc_external("TS/TSGetConvergedReason"))
"""
converged_reason(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetConvergedReason(getlib(PetscLib), ts)

"""
    solve_time(ts::AbstractTS)

The time reached by the last [`solve!`](@ref).

This is the time the integration actually stopped at, which is not the time
asked for by [`set_max_time!`](@ref) unless the final step was made to land on it; 
see [`set_exact_final_time!`](@ref).

# External Links
$(doc_external("TS/TSGetSolveTime"))
"""
solve_time(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetSolveTime(getlib(PetscLib), ts)

"""
    snes_iterations(ts::AbstractTS)

Total number of nonlinear iterations taken so far, summed over the steps.

# External Links
$(doc_external("TS/TSGetSNESIterations"))
"""
snes_iterations(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetSNESIterations(getlib(PetscLib), ts)

"""
    ksp_iterations(ts::AbstractTS)

Total number of linear iterations taken so far, summed over the steps.

# External Links
$(doc_external("TS/TSGetKSPIterations"))
"""
ksp_iterations(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetKSPIterations(getlib(PetscLib), ts)

"""
    step_rejections(ts::AbstractTS)

Number of steps the adaptivity controller has rejected.

# External Links
$(doc_external("TS/TSGetStepRejections"))
"""
step_rejections(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetStepRejections(getlib(PetscLib), ts)

"""
    snes_failures(ts::AbstractTS)

Number of failed nonlinear solves.

# External Links
$(doc_external("TS/TSGetSNESFailures"))
"""
snes_failures(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetSNESFailures(getlib(PetscLib), ts)

"""
    snes(ts::AbstractTS)

The nonlinear solver `ts` steps with.

$(doc_borrowed())

Only the implicit methods build one. Asking an explicit method for its `SNES`
creates an unused solver rather than reporting an error.

# External Links
$(doc_external("TS/TSGetSNES"))
"""
snes(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetSNES(getlib(PetscLib), ts)

"""
    ksp(ts::AbstractTS)

The linear solver `ts` steps with.

$(doc_borrowed())

PETSc only offers this for a problem declared `TS_LINEAR` with [`set_problem_type!`](@ref), 
and raises `PETSC_ERR_ARG_WRONG` otherwise. The linear solver of a nonlinear problem 
belongs to its `SNES`, so reach it through [`snes`](@ref).

# External Links
$(doc_external("TS/TSGetKSP"))
"""
ksp(ts::AbstractTS{PetscLib}) where {PetscLib} =
    LibPETSc.TSGetKSP(getlib(PetscLib), ts)

"""
    setup!(ts::AbstractTS)

Complete the setup of `ts`. [`solve!`](@ref) calls this, so it is only needed
when the setup must happen at a controlled point. Applies the constructor options
first under `options_at = :first_solve` (see [`TS`](@ref)).

# External Links
$(doc_external("TS/TSSetUp"))
"""
function setup!(ts::AbstractTS{PetscLib}) where {PetscLib}
    with_solver_options(ts; solving = false) do
        capture_callback_errors(() -> LibPETSc.TSSetUp(getlib(PetscLib), ts))
    end
    return ts
end

"""
    TSSolve(petsclib, ts, ::Nothing)

Integrate `ts` using the solution vector already set on it.

The generated binding requires a vector. PETSc reads the solution set by
`TSSetSolution` when it is handed `NULL` instead, which is what this overload passes.

# External Links
$(doc_external("TS/TSSolve"))
"""
function LibPETSc.TSSolve(
    petsclib::LibPETSc.PetscLibType,
    ts::LibPETSc.AbstractTS,
    ::Nothing,
) end

LibPETSc.@for_petsc function LibPETSc.TSSolve(
    petsclib::$UnionPetscLib,
    ts::LibPETSc.AbstractTS,
    ::Nothing,
)
    LibPETSc.@chk ccall(
        (:TSSolve, $petsc_library),
        LibPETSc.PetscErrorCode,
        (LibPETSc.CTS, LibPETSc.CVec),
        ts,
        C_NULL,
    )
    return nothing
end

"""
    solve!(u::AbstractPetscVec, ts::AbstractTS; check = true)
    solve!(ts::AbstractTS; check = true)

Integrate `ts`, starting from `u` and returning it. 
The second form uses the solution vector already set on `ts` and returns `ts`.

`u` is registered with `TSSetSolution` before the solve. 
Passing it to `TSSolve` alone leaves the stepper reading an uninitialized solution vector,
which the implicit methods see as an initial condition of zero.

The options passed to the [`TS`](@ref) constructor are applied at the first
`setup!`, `solve!` or `step!` by default; its `options_at` keyword says when.

A failed step stops the run. With `errorifstepfailed` set, PETSc's default (see
[`set_error_if_step_fails!`](@ref)), the failure is thrown as a [`ConvergenceError`](@ref),
or as a `PetscError` with `check = false`. Without it the run returns and
[`converged_reason`](@ref) says why it stopped; `check` decides what happens then:

- `check = true` (the default) throws a [`ConvergenceError`](@ref).
- `check = false` returns normally, and the caller reads `converged_reason`.

# External Links
$(doc_external("TS/TSSolve"))
$(doc_external("TS/TSSetSolution"))
"""
function solve!(
    u::AbstractPetscVec{PetscLib},
    ts::AbstractTS{PetscLib};
    check::Bool = true,
) where {PetscLib}
    LibPETSc.TSSetSolution(PetscLib, ts, u)
    with_solver_options(ts; solving = true) do
        with_convergence_check(ts, check) do
            capture_callback_errors(() -> LibPETSc.TSSolve(PetscLib, ts, u))
        end
    end
    return u
end

function solve!(ts::AbstractTS{PetscLib}; check::Bool = true) where {PetscLib}
    with_solver_options(ts; solving = true) do
        with_convergence_check(ts, check) do
            capture_callback_errors(() -> LibPETSc.TSSolve(PetscLib, ts, nothing))
        end
    end
    return ts
end

"""
    step!(ts::AbstractTS; check = true)

Take a single step. PETSc refuses to step unless a max time or a max number of
steps is set ([`set_max_time!`](@ref), [`set_max_steps!`](@ref)), but `step!` stops
at neither: calling it again steps past them. The step that crosses the max time
can be shortened to land on it, and the time step stays that
short afterwards. The options given to the constructor are applied only under
`options_at = :first_solve`, when it is the first call to do so (see [`TS`](@ref)).

`check` treats a failed step as [`solve!`](@ref) treats a failed run.

# External Links
$(doc_external("TS/TSStep"))
"""
function step!(ts::AbstractTS{PetscLib}; check::Bool = true) where {PetscLib}
    with_solver_options(ts; solving = false) do
        with_convergence_check(ts, check) do
            capture_callback_errors(() -> LibPETSc.TSStep(getlib(PetscLib), ts))
        end
    end
    return ts
end

"""
    reset!(ts::AbstractTS)

Release the work vectors and matrices `ts` allocated, keeping the callbacks
and the options.

The clock is not part of that state: the time and the step count survive, so
set them with [`set_time!`](@ref) and [`set_max_time!`](@ref) before integrating again.

# External Links
$(doc_external("TS/TSReset"))
"""
function reset!(ts::AbstractTS{PetscLib}) where {PetscLib}
    LibPETSc.TSReset(getlib(PetscLib), ts)
    return ts
end

"""
    interpolate!(u::AbstractPetscVec, ts::AbstractTS, t)

Fill `u` with the solution interpolated to time `t`, and return it.

Only the methods that keep a dense output can do this, 
and `t` must lie inside the step just taken.

# External Links
$(doc_external("TS/TSInterpolate"))
"""
function interpolate!(
    u::AbstractPetscVec{PetscLib},
    ts::AbstractTS{PetscLib},
    t,
) where {PetscLib}
    LibPETSc.TSInterpolate(getlib(PetscLib), ts, PetscLib.PetscReal(t), u)
    return u
end

"""
    interpolate(ts::AbstractTS, t)

The solution interpolated to time `t`, in a new vector the caller owns. The
out-of-place [`interpolate!`](@ref), with the same limits.

# External Links
$(doc_external("TS/TSInterpolate"))
"""
interpolate(ts::AbstractTS, t) = interpolate!(similar(solution(ts)), ts, t)

# Callbacks
# ----------------------------------------------------------------------------
#
# Each setter stores the Julia function in the TS's state and hands PETSc a
# `@cfunction` trampoline plus a pointer to that state as the context, so the
# closure lives as long as the PETSc object (src/callbacks.jl). Every
# trampoline runs the user's function through `run_callback`.

"""
    set_rhs_function!(f!, ts::AbstractTS, r = nothing)

Set the right-hand side ``G`` of an explicit problem ``du/dt = G(t, u)``.

`f!` is called as `f!(F, ts, t, u)`, filling the vector `F`.

`r` is an optional template vector for the residual.

$(doc_callback())

# External Links
$(doc_external("TS/TSSetRHSFunction"))

The callback comes first, so `do` block syntax works.
"""
function set_rhs_function! end

mutable struct TSSetRHSFunctionFn{PetscLib, PetscReal} end
function (::TSSetRHSFunctionFn{PetscLib, PetscReal})(
    ts_ptr::CTS,
    t::PetscReal,
    u_ptr::CVec,
    F_ptr::CVec,
    ctx::Ptr{Cvoid},
) where {PetscLib, PetscReal}
    return run_callback("rhs_function!") do
        ts = unsafe_pointer_to_objref(ctx)::TSState
        actual_ts = TS{PetscLib}(ts_ptr, getlib(PetscLib).age; own = false)
        u = PetscVec{PetscLib}(u_ptr, getlib(PetscLib).age; own = false)
        F = PetscVec{PetscLib}(F_ptr, getlib(PetscLib).age; own = false)
        ts.rhs_function!(F, actual_ts, t, u)
    end
end

LibPETSc.@for_petsc function set_rhs_function!(
    f!,
    ts::AbstractTS{$PetscLib},
    r::Union{Nothing, AbstractPetscVec{$PetscLib}} = nothing,
)
    fptr = @cfunction(
        TSSetRHSFunctionFn{$PetscLib, $PetscReal}(),
        LibPETSc.PetscErrorCode,
        (CTS, $PetscReal, CVec, CVec, Ptr{Cvoid})
    )
    ts.rhs_function! = f!
    LibPETSc.TSSetRHSFunction($PetscLib, ts, r, fptr, state_pointer(ts))
    return ts
end

"""
    set_rhs_jacobian!(updateJ!, ts::AbstractTS, A, P = A)

Set the Jacobian of the right-hand side ``G``.

`updateJ!` is called as `updateJ!(A, P, ts, t, u)`, filling the Jacobian `A`
and the preconditioning matrix `P`.

$(doc_callback())

# External Links
$(doc_external("TS/TSSetRHSJacobian"))

The callback comes first, so `do` block syntax works.
"""
function set_rhs_jacobian! end

mutable struct TSSetRHSJacobianFn{PetscLib, PetscReal} end
function (::TSSetRHSJacobianFn{PetscLib, PetscReal})(
    ts_ptr::CTS,
    t::PetscReal,
    u_ptr::CVec,
    A_ptr::CMat,
    P_ptr::CMat,
    ctx::Ptr{Cvoid},
) where {PetscLib, PetscReal}
    return run_callback("rhs_jacobian!") do
        ts = unsafe_pointer_to_objref(ctx)::TSState
        actual_ts = TS{PetscLib}(ts_ptr, getlib(PetscLib).age; own = false)
        u = PetscVec{PetscLib}(u_ptr, getlib(PetscLib).age; own = false)
        A = PetscMat{PetscLib}(A_ptr, getlib(PetscLib).age; own = false)
        P = PetscMat{PetscLib}(P_ptr, getlib(PetscLib).age; own = false)
        ts.rhs_jacobian!(A, P, actual_ts, t, u)
    end
end

LibPETSc.@for_petsc function set_rhs_jacobian!(
    updateJ!,
    ts::AbstractTS{$PetscLib},
    A::AbstractPetscMat{$PetscLib},
    P::AbstractPetscMat{$PetscLib} = A,
)
    fptr = @cfunction(
        TSSetRHSJacobianFn{$PetscLib, $PetscReal}(),
        LibPETSc.PetscErrorCode,
        (CTS, $PetscReal, CVec, CMat, CMat, Ptr{Cvoid})
    )
    ts.rhs_jacobian! = updateJ!
    LibPETSc.TSSetRHSJacobian(
        $PetscLib,
        ts,
        A,
        P,
        fptr,
        state_pointer(ts),
    )
    return ts
end

"""
    set_ifunction!(f!, ts::AbstractTS, r = nothing)

Set the residual ``F`` of an implicit problem ``F(t, u, du/dt) = 0``.

`f!` is called as `f!(F, ts, t, u, u_t)`, filling the vector `F`.

$(doc_callback())

# External Links
$(doc_external("TS/TSSetIFunction"))

The callback comes first, so `do` block syntax works.
"""
function set_ifunction! end

mutable struct TSSetIFunctionFn{PetscLib, PetscReal} end
function (::TSSetIFunctionFn{PetscLib, PetscReal})(
    ts_ptr::CTS,
    t::PetscReal,
    u_ptr::CVec,
    udot_ptr::CVec,
    F_ptr::CVec,
    ctx::Ptr{Cvoid},
) where {PetscLib, PetscReal}
    return run_callback("ifunction!") do
        ts = unsafe_pointer_to_objref(ctx)::TSState
        actual_ts = TS{PetscLib}(ts_ptr, getlib(PetscLib).age; own = false)
        u = PetscVec{PetscLib}(u_ptr, getlib(PetscLib).age; own = false)
        u_t = PetscVec{PetscLib}(udot_ptr, getlib(PetscLib).age; own = false)
        F = PetscVec{PetscLib}(F_ptr, getlib(PetscLib).age; own = false)
        ts.ifunction!(F, actual_ts, t, u, u_t)
    end
end

LibPETSc.@for_petsc function set_ifunction!(
    f!,
    ts::AbstractTS{$PetscLib},
    r::Union{Nothing, AbstractPetscVec{$PetscLib}} = nothing,
)
    fptr = @cfunction(
        TSSetIFunctionFn{$PetscLib, $PetscReal}(),
        LibPETSc.PetscErrorCode,
        (CTS, $PetscReal, CVec, CVec, CVec, Ptr{Cvoid})
    )
    ts.ifunction! = f!
    LibPETSc.TSSetIFunction($PetscLib, ts, r, fptr, state_pointer(ts))
    return ts
end

"""
    set_ijacobian!(updateJ!, ts::AbstractTS, A, P = A)

Set the Jacobian of the implicit residual ``F``.

`updateJ!` is called as `updateJ!(A, P, ts, t, u, u_t, shift)` and should fill
`A` with ``dF/du + shift * dF/du_t``.

$(doc_callback())

# External Links
$(doc_external("TS/TSSetIJacobian"))

The callback comes first, so `do` block syntax works.
"""
function set_ijacobian! end

mutable struct TSSetIJacobianFn{PetscLib, PetscReal} end
function (::TSSetIJacobianFn{PetscLib, PetscReal})(
    ts_ptr::CTS,
    t::PetscReal,
    u_ptr::CVec,
    udot_ptr::CVec,
    shift::PetscReal,
    A_ptr::CMat,
    P_ptr::CMat,
    ctx::Ptr{Cvoid},
) where {PetscLib, PetscReal}
    return run_callback("ijacobian!") do
        ts = unsafe_pointer_to_objref(ctx)::TSState
        actual_ts = TS{PetscLib}(ts_ptr, getlib(PetscLib).age; own = false)
        u = PetscVec{PetscLib}(u_ptr, getlib(PetscLib).age; own = false)
        u_t = PetscVec{PetscLib}(udot_ptr, getlib(PetscLib).age; own = false)
        A = PetscMat{PetscLib}(A_ptr, getlib(PetscLib).age; own = false)
        P = PetscMat{PetscLib}(P_ptr, getlib(PetscLib).age; own = false)
        ts.ijacobian!(A, P, actual_ts, t, u, u_t, shift)
    end
end

LibPETSc.@for_petsc function set_ijacobian!(
    updateJ!,
    ts::AbstractTS{$PetscLib},
    A::AbstractPetscMat{$PetscLib},
    P::AbstractPetscMat{$PetscLib} = A,
)
    fptr = @cfunction(
        TSSetIJacobianFn{$PetscLib, $PetscReal}(),
        LibPETSc.PetscErrorCode,
        (CTS, $PetscReal, CVec, CVec, $PetscReal, CMat, CMat, Ptr{Cvoid})
    )
    ts.ijacobian! = updateJ!
    LibPETSc.TSSetIJacobian(
        $PetscLib,
        ts,
        A,
        P,
        fptr,
        state_pointer(ts),
    )
    return ts
end

"""
    set_monitor!(f, ts::AbstractTS)

Call `f` before the first step and after every accepted step.

`f` is called as `f(ts, step, t, u)`, where `step` counts the steps taken, `t`
is the time reached and `u` holds the solution there.
`u` is owned by `ts` and must not be destroyed.

Only one monitor can be set this way; a second call replaces the first. 
The monitors PETSc installs from the options database, such as `-ts_monitor`, 
are unaffected.

$(doc_callback())

# External Links
$(doc_external("TS/TSMonitorSet"))

The callback comes first, so `do` block syntax works.
"""
function set_monitor! end

mutable struct TSMonitorSetFn{PetscLib, PetscInt, PetscReal} end
function (::TSMonitorSetFn{PetscLib, PetscInt, PetscReal})(
    ts_ptr::CTS,
    step::PetscInt,
    t::PetscReal,
    u_ptr::CVec,
    ctx::Ptr{Cvoid},
) where {PetscLib, PetscInt, PetscReal}
    return run_callback("monitor") do
        ts = unsafe_pointer_to_objref(ctx)::TSState
        actual_ts = TS{PetscLib}(ts_ptr, getlib(PetscLib).age; own = false)
        u = PetscVec{PetscLib}(u_ptr, getlib(PetscLib).age; own = false)
        ts.monitor(actual_ts, step, t, u)
    end
end

LibPETSc.@for_petsc function set_monitor!(f, ts::AbstractTS{$PetscLib})
    fptr = @cfunction(
        TSMonitorSetFn{$PetscLib, $PetscInt, $PetscReal}(),
        LibPETSc.PetscErrorCode,
        (CTS, $PetscInt, $PetscReal, CVec, Ptr{Cvoid})
    )
    ts.monitor = f
    LibPETSc.TSMonitorSet($PetscLib, ts, fptr, state_pointer(ts))
    return ts
end

"""
    set_pre_step!(f!, ts::AbstractTS)

Call `f!(ts)` at the start of every step [`solve!`](@ref) takes, before the
step is attempted. Unlike a monitor, `f!` may change the state of the run,
for example the time step. [`step!`](@ref) does not call it.

A second call replaces the first. Returns `ts`.

$(doc_callback())

# External Links
$(doc_external("TS/TSSetPreStep"))
"""
function set_pre_step! end

"""
    set_post_step!(f!, ts::AbstractTS)

Call `f!(ts)` after every step [`solve!`](@ref) accepts. Unlike a monitor, `f!` may change the state of the run: the solution from
[`solution`](@ref), auxiliary fields, the time step. [`step!`](@ref) does not
call it.

A second call replaces the first. Returns `ts`.

$(doc_callback())

# External Links
$(doc_external("TS/TSSetPostStep"))
"""
function set_post_step! end

# PETSc passes the pre- and post-step hooks no context, only the TS: the state
# is found from the TS pointer
mutable struct TSStepHookFn{PetscLib, Hook} end
function (::TSStepHookFn{PetscLib, Hook})(ts_ptr::CTS) where {PetscLib, Hook}
    # the return value is ignored: these hooks never took an error code. The
    # lookup is inside run_callback too, so that nothing throws through PETSc.
    run_callback(String(Hook)) do
        actual_ts = TS{PetscLib}(ts_ptr, getlib(PetscLib).age; own = false)
        state = object_state(actual_ts)::TSState
        getfield(state, Hook)(actual_ts)
        return nothing
    end
end

LibPETSc.@for_petsc function set_pre_step!(f!, ts::AbstractTS{$PetscLib})
    fptr = @cfunction(TSStepHookFn{$PetscLib, :pre_step!}(), LibPETSc.PetscErrorCode, (CTS,))
    ts.pre_step! = f!
    LibPETSc.TSSetPreStep($PetscLib, ts, fptr)
    return ts
end

LibPETSc.@for_petsc function set_post_step!(f!, ts::AbstractTS{$PetscLib})
    fptr = @cfunction(TSStepHookFn{$PetscLib, :post_step!}(), LibPETSc.PetscErrorCode, (CTS,))
    ts.post_step! = f!
    LibPETSc.TSSetPostStep($PetscLib, ts, fptr)
    return ts
end
