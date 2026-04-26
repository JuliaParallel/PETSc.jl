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

# Builds the bare PETSc TS skeleton shared by every algorithm: pick the
# library, allocate the solution vector, set time bounds and a maybe-supplied
# initial step. Algorithm-specific TS type, subtype, and callback registration
# happen in `_setup_petsc_algorithm!` afterwards.
function _common_ts_setup(prob, alg, dt, maxiters, petsclib, reltol, abstol)
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
        if dt !== nothing
            PETSc.LibPETSc.TSSetTimeStep(lib, ts, lib.PetscReal(dt))
        end
        if reltol !== nothing && abstol !== nothing
            (reltol isa Real && abstol isa Real) || throw(ArgumentError(
                "PETSc.jl SciML extension: only scalar `reltol` / `abstol` are " *
                "supported. Got types $(typeof(reltol)) / $(typeof(abstol)).",
            ))
            _ts_set_scalar_tolerances!(lib, ts, abstol, reltol)
        end
    catch
        PETSc.LibPETSc.TSDestroy(lib, ts)
        PETSc.destroy(u_v)
        rethrow()
    end

    return (lib, ts, u_v, u0, tType, t0, tdir)
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
    PETSc.LibPETSc.TSSetType(lib, ts, alg.ts_type)
    PETSc.LibPETSc.TSSetProblemType(lib, ts, PETSc.LibPETSc.TS_NONLINEAR)
    return _register_ifunction!(lib, ts, prob, u0)
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

# Run callback initialization and, if a callback initializer mutated `u`, push
# the modified state into the PETSc Vec so the first `TSStep` starts from the
# right initial condition.
#
# We deliberately do not push into `sol.t` / `sol.u` here. The first call to
# `step!` saves `t0` exactly once (when `save_start = true`) using the
# possibly-modified `integ.u`, which keeps the trajectory free of duplicate
# timestamps and routes start-time saving through a single code path.
function initialize_callbacks!(integ::PETScTSIntegrator, cb_set)
    DiffEqBase.initialize!(cb_set, integ.u, integ.t, integ)
    if integ.u_modified
        _sync_julia_to_petsc!(integ)
        PETSc.LibPETSc.TSSetSolution(integ.petsclib, integ.ts, integ.u_petsc)
        integ.u_modified = false
    end
    return nothing
end

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
    maxiters::Integer = Int(1e5),
    petsclib = nothing,
    verbose::Bool = false,
    kwargs...,
)
    (lib, ts, u_v, u0, tType, t0, tdir) =
        _common_ts_setup(prob, alg, dt, maxiters, petsclib, reltol, abstol)

    cb_set = DiffEqBase.CallbackSet(callback)
    if !isempty(cb_set.continuous_callbacks)
        @warn "PETSc.jl SciML extension: ContinuousCallbacks are not yet " *
              "supported and will be ignored. Use DiscreteCallback or wrap " *
              "the event detection in PETSc's TSSetEventHandler manually."
    end
    if !isempty(_as_time_iter(tstops, Float64))
        @warn "PETSc.jl SciML extension: `tstops` is not yet honoured. " *
              "PETSc adapts step sizes internally; pass `dt` and " *
              "`-ts_adapt_type none` via `petsc_options` to force fixed steps."
    end

    opts = _build_opts(
        tType, saveat, tstops, tdir, prob.tspan;
        save_everystep, save_on, save_start, save_end,
        callback = cb_set,
        reltol, abstol, maxiters, verbose,
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
        if integ.opts.save_start
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
