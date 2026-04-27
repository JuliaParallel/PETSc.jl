function _destroy_petsc!(integ::PETScTSIntegrator)
    if integ.ts.ptr != C_NULL
        PETSc.LibPETSc.TSDestroy(integ.petsclib, integ.ts)
    end
    if integ.u_petsc.ptr != C_NULL
        PETSc.destroy(integ.u_petsc)
    end
    return nothing
end

# Explicit, deterministic cleanup. Idempotent — repeated calls are no-ops.
PETSc.destroy(integ::PETScTSIntegrator) = _destroy_petsc!(integ)

# Pick the integrator's internal time type. We always promote integer `tspan`s
# to a floating-point type so PETSc-side `Float64` times do not get truncated
# back into `Int` on assignment.
_pick_tType(tspan) = float(eltype(tspan))

function _check_tspan(t0, tf)
    isfinite(t0) && isfinite(tf) || throw(ArgumentError(
        "PETSc.jl SciML extension: tspan endpoints must be finite, got ($t0, $tf).",
    ))
    tf == t0 && throw(ArgumentError(
        "PETSc.jl SciML extension: zero-length tspan ($t0, $tf) is not supported. " *
        "Build the trivial solution (the initial state) directly in user code.",
    ))
    tf > t0 || throw(ArgumentError(
        "PETSc.jl SciML extension: backward integration (tspan = ($t0, $tf)) is " *
        "not yet supported. PETSc TS adaptivity does not handle negative dt " *
        "reliably; reverse the problem yourself or open an issue if you need this.",
    ))
    return nothing
end

# SciML default tolerances applied when the user passes only one of
# `reltol` / `abstol`. This matches OrdinaryDiffEq's defaults so a partial
# specification (e.g. `reltol = 1e-10`) actually reaches PETSc instead of
# being silently dropped.
const _SCIML_DEFAULT_RELTOL = 1e-3
const _SCIML_DEFAULT_ABSTOL = 1e-6

# Builds the bare PETSc TS skeleton shared by every algorithm: pick the
# library, allocate the solution vector, set time bounds and a maybe-supplied
# initial step. Algorithm-specific TS type, subtype, and callback registration
# happen in `_setup_petsc_algorithm!` afterwards.
function _common_ts_setup(
    prob, dt, maxiters, petsclib, reltol, abstol,
    adaptive, dtmin, dtmax,
)
    _check_isinplace(prob)
    lib = _pick_petsclib(prob, petsclib)
    _check_petscreal(lib)
    PETSc.initialized(lib) || PETSc.initialize(lib)

    u0 = copy(prob.u0)
    tType = _pick_tType(prob.tspan)
    t0 = tType(prob.tspan[1])
    tf = tType(prob.tspan[2])
    _check_tspan(t0, tf)
    tdir = tType(sign(tf - t0))

    ts = PETSc.LibPETSc.TSCreate(lib, PETSc.LibPETSc.PETSC_COMM_SELF)
    u_v = PETSc.VecSeq(lib, length(u0))
    try
        PETSc.withlocalarray!(u_v; read = false, write = true) do arr
            copyto!(arr, vec(u0))
        end
        PETSc.LibPETSc.TSSetSolution(lib, ts, u_v)
        PETSc.LibPETSc.TSSetTime(lib, ts, lib.PetscReal(t0))
        PETSc.LibPETSc.TSSetMaxTime(lib, ts, lib.PetscReal(tf))
        PETSc.LibPETSc.TSSetMaxSteps(lib, ts, lib.PetscInt(maxiters))
        PETSc.LibPETSc.TSSetExactFinalTime(
            lib, ts, PETSc.LibPETSc.TS_EXACTFINALTIME_MATCHSTEP,
        )
        # Validate `dt` against the SciML `dtmin` / `dtmax` bounds *before*
        # installing it: PETSc applies the initial step verbatim and only
        # consults `TSAdaptSetStepLimits` for subsequent step proposals, so
        # silently allowing `dt = 0.5, dtmax = 0.2` would let the very first
        # step violate the user-supplied bound.
        _check_dt_bounds(dt, dtmin, dtmax)
        if dt !== nothing
            PETSc.LibPETSc.TSSetTimeStep(lib, ts, lib.PetscReal(dt))
        end
        _apply_tolerances!(lib, ts, reltol, abstol)
        _apply_adaptivity!(lib, ts, adaptive, dtmin, dtmax)
    catch
        PETSc.LibPETSc.TSDestroy(lib, ts)
        PETSc.destroy(u_v)
        rethrow()
    end

    return (lib, ts, u_v, u0, tType, t0, tdir)
end

# Forward `reltol` / `abstol` to PETSc when at least one is set. Filling the
# missing side from SciML's defaults is the convention upstream wrappers use:
# `solve(prob, alg; reltol = 1e-10)` should reach the adaptive controller
# rather than be silently ignored.
function _apply_tolerances!(lib, ts, reltol, abstol)
    (reltol === nothing && abstol === nothing) && return nothing
    rt = reltol === nothing ? _SCIML_DEFAULT_RELTOL : reltol
    at = abstol === nothing ? _SCIML_DEFAULT_ABSTOL : abstol
    (rt isa Real && at isa Real) || throw(ArgumentError(
        "PETSc.jl SciML extension: only scalar `reltol` / `abstol` are " *
        "supported. Got types $(typeof(reltol)) / $(typeof(abstol)).",
    ))
    _ts_set_scalar_tolerances!(lib, ts, at, rt)
    return nothing
end

# Reject contradictory inputs where the user-supplied initial `dt` would
# violate the user-supplied `dtmin` / `dtmax` bounds. PETSc installs the
# initial step verbatim, so a clean `ArgumentError` is the only way to
# uphold the SciML "all steps respect the bound" contract from step one.
function _check_dt_bounds(dt, dtmin, dtmax)
    dt === nothing && return nothing
    if dtmax !== nothing && dt > dtmax
        throw(ArgumentError(
            "PETSc.jl SciML extension: initial `dt = $dt` exceeds `dtmax = $dtmax`. " *
            "Pass `dt <= dtmax`, or omit `dt` to let PETSc choose the initial step.",
        ))
    end
    if dtmin !== nothing && dt < dtmin
        throw(ArgumentError(
            "PETSc.jl SciML extension: initial `dt = $dt` is below `dtmin = $dtmin`. " *
            "Pass `dt >= dtmin`, or omit `dt` to let PETSc choose the initial step.",
        ))
    end
    return nothing
end

# Map SciML's `adaptive` / `dtmin` / `dtmax` knobs onto the PETSc TSAdapt
# controller. `adaptive = false` selects PETSc's "none" adapter (fixed
# `dt`); `dtmin` / `dtmax` set step limits via `TSAdaptSetStepLimits`.
function _apply_adaptivity!(lib, ts, adaptive, dtmin, dtmax)
    (adaptive === true && dtmin === nothing && dtmax === nothing) && return nothing
    adapt = PETSc.LibPETSc.TSGetAdapt(lib, ts)
    if adaptive === false
        PETSc.LibPETSc.TSAdaptSetType(lib, adapt, "none")
    end
    if dtmin !== nothing || dtmax !== nothing
        hmin = dtmin === nothing ? zero(lib.PetscReal) : lib.PetscReal(dtmin)
        hmax = dtmax === nothing ? lib.PetscReal(Inf) : lib.PetscReal(dtmax)
        PETSc.LibPETSc.TSAdaptSetStepLimits(lib, adapt, hmin, hmax)
    end
    return nothing
end

function _make_integrator(
    alg, u0, tType, t0, tdir, dt, prob,
    opts, sol, lib, ts, u_v, cb_ctx,
)
    integ = PETScTSIntegrator(
        alg,
        u0,
        copy(u0),
        t0,
        t0,
        tType(something(dt, zero(tType))),
        prob.p,
        opts,
        false,           # u_modified
        false,           # derivative_discontinuity
        tdir,
        size(u0),
        sol,
        nothing,
        lib,
        ts,
        u_v,
        cb_ctx,
        false,
        false,
        SciMLBase.ReturnCode.Default,
    )
    finalizer(_destroy_petsc!, integ)
    return integ
end

function _register_rhs!(lib, ts, prob, u0)
    cb_ctx = RHSCtx(prob.f.f, prob.p, size(u0), lib)
    PETSc.LibPETSc.TSSetRHSFunction(
        lib, ts, nothing, _petsc_rhs_ptr(), pointer_from_objref(cb_ctx),
    )
    return cb_ctx
end

function _register_ifunction!(lib, ts, prob, u0)
    cb_ctx = IFunctionCtx(prob.f.f, prob.p, size(u0), lib)
    PETSc.LibPETSc.TSSetIFunction(
        lib, ts, nothing, _petsc_ifunction_ptr(), pointer_from_objref(cb_ctx),
    )
    return cb_ctx
end

# Per-algorithm hooks. Each returns the callback context object that needs to
# stay live for the lifetime of the integrator.
function _setup_petsc_algorithm!(lib, ts, prob, u0, alg::TSRK)
    PETSc.LibPETSc.TSSetType(lib, ts, "rk")
    PETSc.LibPETSc.TSRKSetType(lib, ts, alg.subtype)
    return _register_rhs!(lib, ts, prob, u0)
end

function _setup_petsc_algorithm!(lib, ts, prob, u0, alg::TSRosW)
    PETSc.LibPETSc.TSSetType(lib, ts, "rosw")
    PETSc.LibPETSc.TSRosWSetType(lib, ts, alg.subtype)
    PETSc.LibPETSc.TSSetProblemType(lib, ts, PETSc.LibPETSc.TS_NONLINEAR)
    return _register_ifunction!(lib, ts, prob, u0)
end

const _IMPLICIT_SUBTYPES = ("beuler", "cn", "theta", "bdf")

function _check_implicit_subtype(s)
    s in _IMPLICIT_SUBTYPES || throw(ArgumentError(
        "TSImplicit subtype must be one of " *
        join(_IMPLICIT_SUBTYPES, ", ") * "; got $(repr(s))",
    ))
end

function _setup_petsc_algorithm!(lib, ts, prob, u0, alg::TSImplicit)
    _check_implicit_subtype(alg.subtype)
    PETSc.LibPETSc.TSSetType(lib, ts, alg.subtype)
    if alg.subtype == "theta"
        PETSc.LibPETSc.TSThetaSetTheta(lib, ts, lib.PetscReal(alg.theta))
    end
    PETSc.LibPETSc.TSSetProblemType(lib, ts, PETSc.LibPETSc.TS_NONLINEAR)
    return _register_ifunction!(lib, ts, prob, u0)
end

# Per-side context registration for IMEX. The integrator field stores both
# contexts as a NamedTuple so they remain GC-rooted via the integrator.
function _register_rhs_with_f!(lib, ts, f, prob, u0)
    cb_ctx = RHSCtx(f, prob.p, size(u0), lib)
    PETSc.LibPETSc.TSSetRHSFunction(
        lib, ts, nothing, _petsc_rhs_ptr(), pointer_from_objref(cb_ctx),
    )
    return cb_ctx
end

function _register_ifunction_with_f!(lib, ts, f, prob, u0)
    cb_ctx = IFunctionCtx(f, prob.p, size(u0), lib)
    PETSc.LibPETSc.TSSetIFunction(
        lib, ts, nothing, _petsc_ifunction_ptr(), pointer_from_objref(cb_ctx),
    )
    return cb_ctx
end

function _setup_petsc_algorithm!(lib, ts, prob, u0, alg::TSGeneric)
    _check_tsgeneric_type(alg)
    PETSc.LibPETSc.TSSetType(lib, ts, alg.ts_type)
    if alg.explicit
        # Explicit PETSc TS types (e.g. `"euler"`, `"ssp"`) reject an
        # IFunction at TSStep time — register the RHS instead.
        return _register_rhs!(lib, ts, prob, u0)
    end
    PETSc.LibPETSc.TSSetProblemType(lib, ts, PETSc.LibPETSc.TS_NONLINEAR)
    return _register_ifunction!(lib, ts, prob, u0)
end

# PETSc TS types known to require the explicit (RHS) calling convention.
# Using one of these without `explicit = true` would otherwise dump a full
# PETSc error banner at `TSSetUp` time before the wrapper rethrows; catching
# it on the Julia side gives the user a clear, actionable message instead.
const _EXPLICIT_ONLY_TS_TYPES = ("euler", "ssp")

function _check_tsgeneric_type(alg::TSGeneric)
    if !alg.explicit && alg.ts_type in _EXPLICIT_ONLY_TS_TYPES
        throw(ArgumentError(
            "PETSc.jl SciML extension: TSGeneric ts_type $(repr(alg.ts_type)) " *
            "is an explicit-only PETSc TS family and must be constructed " *
            "with `explicit = true`, e.g. " *
            "`PETSc.TSGeneric($(repr(alg.ts_type)); explicit = true)`. " *
            "The default `explicit = false` registers an IFunction, which " *
            "PETSc rejects for explicit-only types at `TSSetUp` time.",
        ))
    end
    return nothing
end

function _setup_petsc_algorithm!(lib, ts, prob, u0, alg::TSARKIMEX)
    PETSc.LibPETSc.TSSetType(lib, ts, "arkimex")
    PETSc.LibPETSc.TSARKIMEXSetType(lib, ts, alg.subtype)
    PETSc.LibPETSc.TSSetProblemType(lib, ts, PETSc.LibPETSc.TS_NONLINEAR)

    if prob.problem_type isa SciMLBase.SplitODEProblem
        f1 = prob.f.f1     # implicit / stiff
        f2 = prob.f.f2     # explicit / non-stiff
        rhs_ctx = _register_rhs_with_f!(lib, ts, f2, prob, u0)
        ifunc_ctx = _register_ifunction_with_f!(lib, ts, f1, prob, u0)
        return (rhs = rhs_ctx, ifunc = ifunc_ctx)
    else
        return _register_ifunction!(lib, ts, prob, u0)
    end
end

# Mirror the OrdinaryDiffEq initialization contract:
#
# 1. Pessimistically mark `u` as modified before calling `initialize!`. A
#    callback initializer that does not mutate `u` is expected to call
#    `DiffEqBase.u_modified!(integ, false)`; otherwise we conservatively assume
#    it did and resync the PETSc Vec.
# 2. Force an initialize-time save when any discrete callback requests
#    `save_positions[2]` (the post-event side, which corresponds to "after the
#    callback ran"). Duplicate suppression in `step!` keeps `t0` from being
#    recorded twice when `save_start = true`.
function initialize_callbacks!(integ::PETScTSIntegrator, cb_set)
    integ.u_modified = true
    DiffEqBase.initialize!(cb_set, integ.u, integ.t, integ)
    if integ.u_modified
        _sync_julia_to_petsc!(integ)
        PETSc.LibPETSc.TSSetSolution(integ.petsclib, integ.ts, integ.u_petsc)
        integ.u_modified = false
    end

    if integ.opts.save_on && _any_initialize_save(cb_set) &&
       (isempty(integ.sol.t) || last(integ.sol.t) != integ.t)
        push!(integ.sol.t, integ.t)
        push!(integ.sol.u, copy(integ.u))
    end
    return nothing
end

_any_initialize_save(cb_set) =
    any(cb -> cb.save_positions[2], cb_set.discrete_callbacks)

function SciMLBase.__init(
    prob::SciMLBase.AbstractODEProblem,
    alg::PETScTSAlgorithm;
    save_everystep::Bool = false,
    save_on::Bool = true,
    save_start::Bool = true,
    save_end::Bool = true,
    saveat = (),
    tstops = (),
    callback = nothing,
    reltol = nothing,
    abstol = nothing,
    dt = nothing,
    dtmin = nothing,
    dtmax = nothing,
    adaptive::Bool = true,
    maxiters::Integer = Int(1e5),
    petsclib = nothing,
    kwargs...,
)
    _reject_unsupported_kwargs(kwargs)
    cb_set = DiffEqBase.CallbackSet(callback)
    if !isempty(cb_set.continuous_callbacks)
        throw(ArgumentError(
            "PETSc.jl SciML extension: ContinuousCallbacks are not yet " *
            "supported. Their `initialize` and `finalize` hooks would still " *
            "run while event detection is silently dropped, so we reject " *
            "them up front. Use DiscreteCallback or wrap the event detection " *
            "in PETSc's TSSetEventHandler manually.",
        ))
    end

    (lib, ts, u_v, u0, tType, t0, tdir) = _common_ts_setup(
        prob, dt, maxiters, petsclib, reltol, abstol,
        adaptive, dtmin, dtmax,
    )

    if !isempty(_as_time_iter(tstops, Float64))
        @warn "PETSc.jl SciML extension: `tstops` is not yet honoured. " *
              "PETSc adapts step sizes internally; pass `dt` and " *
              "`-ts_adapt_type none` via `petsc_options` to force fixed steps."
    end

    opts = _build_opts(
        tType, saveat, tstops, tdir, prob.tspan;
        save_everystep, save_on, save_start, save_end,
        callback = cb_set,
        reltol, abstol, maxiters,
    )

    sol = SciMLBase.build_solution(
        prob, alg, tType[], typeof(u0)[];
        retcode = SciMLBase.ReturnCode.Default,
        stats = SciMLBase.DEStats(0),
    )

    try
        cb_ctx = _setup_petsc_algorithm!(lib, ts, prob, u0, alg)
        _setfromoptions!(lib, ts, alg.petsc_options)

        integ = _make_integrator(
            alg, u0, tType, t0, tdir, dt, prob,
            opts, sol, lib, ts, u_v, cb_ctx,
        )

        initialize_callbacks!(integ, cb_set)

        return integ
    catch
        # Setup failed before _make_integrator wired up the finalizer; clean up
        # the half-constructed PETSc objects so they do not leak.
        if ts.ptr != C_NULL
            PETSc.LibPETSc.TSDestroy(lib, ts)
        end
        if u_v.ptr != C_NULL
            PETSc.destroy(u_v)
        end
        rethrow()
    end
end

function SciMLBase.__solve(
    prob::SciMLBase.AbstractODEProblem,
    alg::PETScTSAlgorithm,
    args...;
    kwargs...,
)
    integ = SciMLBase.__init(prob, alg; kwargs...)
    return SciMLBase.solve!(integ)
end

function SciMLBase.step!(integ::PETScTSIntegrator)
    integ.done && return nothing

    if !integ.initialized
        PETSc.LibPETSc.TSSetUp(integ.petsclib, integ.ts)
        integ.initialized = true
        if integ.opts.save_start &&
           (isempty(integ.sol.t) || last(integ.sol.t) != integ.t)
            push!(integ.sol.t, integ.t)
            push!(integ.sol.u, copy(integ.u))
        end
    end

    integ.uprev .= integ.u
    integ.tprev = integ.t

    GC.@preserve integ begin
        PETSc.LibPETSc.TSStep(integ.petsclib, integ.ts)
    end

    _sync_petsc_to_julia!(integ)
    integ.t = typeof(integ.t)(PETSc.LibPETSc.TSGetTime(integ.petsclib, integ.ts))
    integ.dt = typeof(integ.dt)(PETSc.LibPETSc.TSGetTimeStep(integ.petsclib, integ.ts))

    handle_callbacks!(integ)

    if integ.u_modified
        integ.u_modified = false
        _sync_julia_to_petsc!(integ)
        PETSc.LibPETSc.TSSetSolution(integ.petsclib, integ.ts, integ.u_petsc)
    end
    integ.derivative_discontinuity = false

    if integ.done
        # `terminate!` was triggered by a callback. Tell PETSc not to keep
        # stepping past the current time and exit.
        PETSc.LibPETSc.TSSetMaxTime(
            integ.petsclib, integ.ts, integ.petsclib.PetscReal(integ.t),
        )
        return nothing
    end

    reason = _ts_converged_reason(integ.petsclib, integ.ts)
    if reason != PETSc.LibPETSc.TS_CONVERGED_ITERATING
        integ.retcode = _petsc_retcode(integ.petsclib, integ.ts)
        integ.done = true
    elseif integ.tdir * (integ.t - integ.sol.prob.tspan[2]) >= 0
        integ.retcode = SciMLBase.ReturnCode.Success
        integ.done = true
    end

    return nothing
end

function SciMLBase.solve!(integ::PETScTSIntegrator)
    while !integ.done
        SciMLBase.step!(integ)
    end
    DiffEqBase.finalize!(integ.opts.callback, integ.u, integ.t, integ)
    if integ.opts.save_end &&
       (isempty(integ.sol.t) || last(integ.sol.t) != integ.t)
        push!(integ.sol.t, integ.t)
        push!(integ.sol.u, copy(integ.u))
    end
    integ.sol = SciMLBase.solution_new_retcode(integ.sol, integ.retcode)
    _destroy_petsc!(integ)
    return integ.sol
end
