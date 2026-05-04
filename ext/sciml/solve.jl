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
        # Screen all step-control input on the Julia side before touching
        # PETSc: bad `dt` / `dtmin` / `dtmax` would otherwise either be
        # silently accepted (for some pathological values) or fall through
        # to a raw PETSc error banner.
        _validate_step_control(dt, dtmin, dtmax)
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
    _validate_tolerance(:reltol, reltol, rt)
    _validate_tolerance(:abstol, abstol, at)
    _ts_set_scalar_tolerances!(lib, ts, at, rt)
    return nothing
end

# Reject scalar tolerances that PETSc's `TSSetTolerances` would refuse:
# non-finite or strictly negative values. `0` is allowed (it disables that
# side of the `atol + rtol * |u|` test) so the policy is "non-negative,
# finite, real". The `user_value` argument is `nothing` when the SciML
# default kicked in, in which case there is nothing to validate.
function _validate_tolerance(name::Symbol, user_value, applied_value)
    user_value === nothing && return nothing
    (isfinite(applied_value) && applied_value >= 0) || throw(ArgumentError(
        "PETSc.jl SciML extension: `$(name) = $(user_value)` is not a valid " *
        "tolerance. Tolerances must be finite and non-negative.",
    ))
    return nothing
end

# Screen `dt`, `dtmin`, and `dtmax` against PETSc's expectations before any
# of them reach the underlying TS / TSAdapt API. Catches the common
# misuse cases (negative / non-finite / contradictory bounds) at the Julia
# boundary instead of letting them fall through to raw PETSc error banners.
function _validate_step_control(dt, dtmin, dtmax)
    _validate_step_size(:dt, dt; allow_zero = false)
    _validate_step_size(:dtmin, dtmin; allow_zero = true)
    _validate_step_size(:dtmax, dtmax; allow_zero = false)
    if dtmin !== nothing && dtmax !== nothing && dtmin > dtmax
        throw(ArgumentError(
            "PETSc.jl SciML extension: `dtmin = $dtmin` exceeds `dtmax = $dtmax`. " *
            "Adaptive step limits must satisfy `dtmin <= dtmax`.",
        ))
    end
    if dt !== nothing
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
    end
    return nothing
end

# Reject negative SciML `maxiters`. With the manual `TSStep` loop, a
# negative cap would otherwise immediately exhaust and turn into a
# zero-step `MaxIters` solve, which mismatches PETSc's own
# `-ts_max_steps -1 = unlimited` interpretation. Keep `maxiters = 0`
# valid (it is the documented "zero-step solve" case from Review-10).
function _validate_maxiters(maxiters)
    maxiters < 0 && throw(ArgumentError(
        "PETSc.jl SciML extension: `maxiters = $(maxiters)` must be " *
        "non-negative. PETSc's `-ts_max_steps -1` means \"unlimited\", but " *
        "the SciML wrapper enforces the cap in its own loop, so a negative " *
        "value would be a zero-step solve. Use `maxiters = 0` if that is " *
        "what you want, or a large positive integer for an effective " *
        "no-cap solve.",
    ))
    return nothing
end

function _validate_step_size(name::Symbol, value; allow_zero::Bool)
    value === nothing && return nothing
    value isa Real || throw(ArgumentError(
        "PETSc.jl SciML extension: `$(name) = $(value)` must be a real scalar, " *
        "got type $(typeof(value)).",
    ))
    # `dtmax = Inf` is the natural SciML spelling of "do not cap the step";
    # the wrapper already maps an omitted `dtmax` to `Inf` internally, so
    # accepting the explicit form keeps the public API consistent with the
    # internal semantics. `dt = Inf` and `dtmin = Inf` remain rejected.
    if name === :dtmax && value == Inf
        return nothing
    end
    isfinite(value) || throw(ArgumentError(
        "PETSc.jl SciML extension: `$(name) = $(value)` must be finite.",
    ))
    if allow_zero
        value < 0 && throw(ArgumentError(
            "PETSc.jl SciML extension: `$(name) = $(value)` must be non-negative.",
        ))
    else
        value > 0 || throw(ArgumentError(
            "PETSc.jl SciML extension: `$(name) = $(value)` must be strictly positive. " *
            "Backward integration is rejected upstream by `_check_tspan`; pass a " *
            "positive step size.",
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
#    `SciMLBase.u_modified!(integ, false)`; otherwise we conservatively assume
#    it did and resync the PETSc Vec.
# 2. Force an initialize-time save when any discrete callback requests
#    `save_positions[2]` (the post-event side, which corresponds to "after the
#    callback ran"). Duplicate suppression in `step!` keeps `t0` from being
#    recorded twice when `save_start = true`.
function initialize_callbacks!(
    integ::PETScTSIntegrator, cb_set, initialize_save::Bool = true,
)
    integ.u_modified = true
    _initialize_callbacks!(cb_set, integ.u, integ.t, integ)
    if integ.u_modified
        _sync_julia_to_petsc!(integ)
        PETSc.LibPETSc.TSSetSolution(integ.petsclib, integ.ts, integ.u_petsc)
        integ.u_modified = false
    end

    # `initialize_save = false` is the upstream SciML knob for "run callback
    # `initialize` hooks but do not append a post-init save record". Honour
    # both that gate and `save_on` here, mirroring OrdinaryDiffEq's
    # `initialize_callbacks!(integrator, initialize_save)` flow.
    if initialize_save && integ.opts.save_on && _any_initialize_save(cb_set) &&
       (isempty(integ.sol.t) || last(integ.sol.t) != integ.t)
        push!(integ.sol.t, integ.t)
        push!(integ.sol.u, copy(integ.u))
    end
    # Forward the SciML discrete-save lifecycle hook so callbacks that
    # populate observable timeseries (MTK-style `saved_clock_partitions`,
    # `initialize_save_discretes`) get their `t0` snapshot. `skip_duplicates`
    # mirrors what OrdinaryDiffEq passes here. The hook is only available
    # in newer SciMLBase releases, so guard with `isdefined`. Skip it
    # entirely when the user opted out via `initialize_save = false`.
    if initialize_save && isdefined(SciMLBase, :save_discretes_if_enabled!)
        SciMLBase.save_discretes_if_enabled!(integ, cb_set; skip_duplicates = true)
    end
    return nothing
end

_any_initialize_save(cb_set) =
    any(cb -> cb.save_positions[2], cb_set.discrete_callbacks)

function SciMLBase.__init(
    prob::SciMLBase.AbstractODEProblem,
    alg::PETScTSAlgorithm;
    save_everystep = nothing,   # default: true when saveat is empty, false otherwise
    save_on::Bool = true,
    save_start::Bool = true,
    save_end::Bool = true,
    save_discretes::Bool = true,
    saveat = (),
    tstops = (),
    callback = nothing,
    initialize_save::Bool = true,
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
    # Pure Julia validation that does not need PETSc state runs first, so
    # bad input fails before the wrapper allocates any PETSc TS / Vec
    # handles. Otherwise the cleanup `try` further down would never see
    # those exceptions and PETSc objects would leak (the autowrapped `TS`
    # has no finalizer of its own).
    _reject_unsupported_kwargs(kwargs)
    _validate_maxiters(maxiters)
    cb_set = SciMLBase.CallbackSet(callback)
    if !isempty(cb_set.continuous_callbacks)
        throw(ArgumentError(
            "PETSc.jl SciML extension: ContinuousCallbacks are not yet " *
            "supported. Their `initialize` and `finalize` hooks would still " *
            "run while event detection is silently dropped, so we reject " *
            "them up front. Use DiscreteCallback or wrap the event detection " *
            "in PETSc's TSSetEventHandler manually.",
        ))
    end

    # `tstops` is documented as a SciML contract for "the integrator must
    # land exactly on these times so step-end callback logic can see them",
    # which means silently ignoring it can skip exact-time discrete
    # callbacks. Reject it up front until the manual `TSStep` loop honours
    # it, so users get a loud error rather than a silently wrong solve.
    # Materialize the iterable first so a one-shot iterator is never
    # consumed by the non-empty check.
    tstops_materialized = _materialize_times(tstops)
    if !_is_empty_times(tstops_materialized)
        throw(ArgumentError(
            "PETSc.jl SciML extension: `tstops` is not yet honoured. " *
            "Accepting it silently would risk skipping exact-time discrete " *
            "callbacks. Either remove the `tstops` keyword, or use " *
            "PETSc's `TSSetEventHandler` directly. To force fixed-step " *
            "integration, pass `adaptive = false` (or " *
            "`-ts_adapt_type none` via the algorithm's `petsc_options`).",
        ))
    end

    # Pre-validate `saveat` here too. `_build_opts` would otherwise throw
    # on bad `saveat` only after `_common_ts_setup` had already allocated
    # PETSc objects that have no finalizer of their own. Materialize the
    # iterable exactly once so a stateful iterator survives both this
    # validation pass *and* `_build_opts`'s own consumption.
    saveat_materialized = _materialize_times(saveat)
    _validate_saveat(saveat_materialized)
    # Follow the SciML convention: save_everystep defaults to true when saveat
    # is empty (record every step), and false when saveat times are given (only
    # the requested times matter, not every intermediate step).
    _save_everystep::Bool = save_everystep === nothing ?
        _is_empty_times(saveat_materialized) : Bool(save_everystep)

    (lib, ts, u_v, u0, tType, t0, tdir) = _common_ts_setup(
        prob, dt, maxiters, petsclib, reltol, abstol,
        adaptive, dtmin, dtmax,
    )

    try
        opts = _build_opts(
            tType, saveat_materialized, tstops_materialized, tdir, prob.tspan;
            save_everystep = _save_everystep, save_on, save_start, save_end, save_discretes,
            callback = cb_set,
            reltol, abstol, maxiters,
        )

        # `DEStats()` defaults every counter to `-1`, which SciML reads as
        # "unknown / not reported". `_populate_stats!` overwrites the
        # fields we can map accurately (`naccept`, `nreject`,
        # `nnonliniter`, `nf`) before `solve!` returns; everything else
        # stays at the sentinel so users can distinguish "no work
        # happened" from "we don't track this".
        sol = SciMLBase.build_solution(
            prob, alg, tType[], typeof(u0)[];
            retcode = SciMLBase.ReturnCode.Default,
            stats = SciMLBase.DEStats(),
        )

        cb_ctx = _setup_petsc_algorithm!(lib, ts, prob, u0, alg)
        _setfromoptions!(lib, ts, alg.petsc_options)

        integ = _make_integrator(
            alg, u0, tType, t0, tdir, dt, prob,
            opts, sol, lib, ts, u_v, cb_ctx,
        )

        initialize_callbacks!(integ, cb_set, initialize_save)

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

# `PETScTSAlgorithm` does not inherit from `SciMLBase.AbstractODEAlgorithm` so
# that the algorithm structs can be defined in `src/` without a hard SciMLBase
# dependency. That means SciMLBase's generic CommonSolve.solve / init dispatch
# (which requires AbstractSciMLAlgorithm) does not match. Provide thin glue
# methods here (inside the extension, where SciMLBase is already available) so
# that the standard `solve(prob, alg; kwargs...)` / `init(prob, alg; kwargs...)`
# entry points reach our `__solve` / `__init` implementations.
function SciMLBase.solve(
    prob::SciMLBase.AbstractODEProblem,
    alg::PETScTSAlgorithm,
    args...;
    kwargs...,
)
    SciMLBase.__solve(prob, alg, args...; kwargs...)
end

function SciMLBase.init(
    prob::SciMLBase.AbstractODEProblem,
    alg::PETScTSAlgorithm,
    args...;
    kwargs...,
)
    SciMLBase.__init(prob, alg, args...; kwargs...)
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

    # Enforce the step cap *before* `TSStep` so `maxiters = 0` (or
    # `-ts_max_steps 0`) yields a true zero-step solve instead of
    # advancing once and reporting the cap afterwards.
    if _ts_step_count(integ) >= _effective_maxiters(integ)
        integ.retcode = SciMLBase.ReturnCode.MaxIters
        integ.done = true
        return nothing
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
    elseif _ts_step_count(integ) >= _effective_maxiters(integ)
        # PETSc's own `TSSetMaxSteps` (and the equivalent `-ts_max_steps`
        # option) is honoured inside `TSSolve`, but the extension drives
        # `TSStep` directly so the cap has to be enforced in this loop.
        # `_effective_maxiters` takes the min of the SciML-side `maxiters`
        # and PETSc's `TSGetMaxSteps`, so an algorithm-side
        # `-ts_max_steps` flag is honoured too.
        integ.retcode = SciMLBase.ReturnCode.MaxIters
        integ.done = true
    end

    return nothing
end

# Read PETSc's accepted-step count via the autowrapped getter. Wrapped in a
# helper so the `step!` and `solve!` paths both see the same value.
_ts_step_count(integ::PETScTSIntegrator) =
    Int(PETSc.LibPETSc.TSGetStepNumber(integ.petsclib, integ.ts))

# Combine the SciML-side `maxiters` and whatever PETSc currently has stored
# as its `TSSetMaxSteps` into a single effective cap. This is what makes
# algorithm-side `-ts_max_steps` reach the manual `TSStep` loop too.
function _effective_maxiters(integ::PETScTSIntegrator)
    petsc_max = Int(PETSc.LibPETSc.TSGetMaxSteps(integ.petsclib, integ.ts))
    return min(integ.opts.maxiters, petsc_max)
end

# Pull whatever PETSc bookkeeping is meaningful for the current TS family
# back into the SciML `DEStats` object. Only counters we can map accurately
# are written; everything else stays at its `DEStats` initialiser value
# (`-1`, the SciML "unknown" sentinel — see
# `SciMLBase.DEStats(x = -1) = DEStats(x, x, ..., 0.0)`).
#
# - `naccept` ← `TSGetStepNumber`
# - `nreject` ← `TSGetStepRejections`
# - `nnonliniter` ← `TSGetSNESIterations` (zero on explicit families)
# - `nf` / `nf2` ← user-callback hit counts maintained on the
#         RHS / IFunction contexts. PETSc itself does not expose a uniform
#         "RHS calls" counter, so we tally evaluations in the C-callback.
#         For split IMEX problems the implicit (`ifunc`) stream lands in
#         `nf` and the explicit (`rhs`) stream in `nf2`, matching SciML
#         convention.
#
# `nsolve` is intentionally left at the `DEStats` sentinel: PETSc's
# `TSGetKSPIterations` returns linear *iteration* counts, while SciML's
# `nsolve` is the number of linear *solves*, so they are not equivalent.
function _populate_stats!(integ::PETScTSIntegrator)
    stats = integ.sol.stats
    stats === nothing && return nothing
    stats.naccept = _ts_step_count(integ)
    stats.nreject = Int(PETSc.LibPETSc.TSGetStepRejections(integ.petsclib, integ.ts))
    stats.nnonliniter = Int(PETSc.LibPETSc.TSGetSNESIterations(integ.petsclib, integ.ts))
    _populate_nf!(stats, integ.cb_ctx)
    return nothing
end

# Single-stream callback contexts: the user RHS is the only function being
# called, so all evaluations roll up into `nf`.
_populate_nf!(stats, ctx::Union{RHSCtx, IFunctionCtx}) = (stats.nf = ctx.nf; nothing)

# Split IMEX context: keep the implicit / explicit streams separate so
# users can distinguish stiff vs. non-stiff function work. This matches
# upstream OrdinaryDiffEq, which uses `nf` for `f1` and `nf2` for `f2` on
# `SplitODEProblem`s.
function _populate_nf!(stats, ctx::NamedTuple)
    haskey(ctx, :ifunc) && (stats.nf = ctx.ifunc.nf)
    haskey(ctx, :rhs) && (stats.nf2 = ctx.rhs.nf)
    return nothing
end

_populate_nf!(_, _) = nothing

function SciMLBase.solve!(integ::PETScTSIntegrator)
    while !integ.done
        SciMLBase.step!(integ)
    end
    _finalize_callbacks!(integ.opts.callback, integ.u, integ.t, integ)
    if integ.opts.save_end &&
       (isempty(integ.sol.t) || last(integ.sol.t) != integ.t)
        push!(integ.sol.t, integ.t)
        push!(integ.sol.u, copy(integ.u))
    end
    # Forward the SciML discrete-save lifecycle hook so callbacks that
    # populate observable timeseries get their final snapshot. Mirrors
    # OrdinaryDiffEq's end-of-solve sequence: `finalize!` first, then the
    # final discrete save. The hook is only available in newer SciMLBase
    # releases, so guard with `isdefined`.
    if isdefined(SciMLBase, :save_final_discretes!)
        SciMLBase.save_final_discretes!(integ, integ.opts.callback)
    end
    _populate_stats!(integ)
    integ.sol = SciMLBase.solution_new_retcode(integ.sol, integ.retcode)
    _destroy_petsc!(integ)
    return integ.sol
end
