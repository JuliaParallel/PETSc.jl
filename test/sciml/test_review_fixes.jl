using Test
using PETSc
using SciMLBase

ext = Base.get_extension(PETSc, :PETScSciMLExt)
@assert ext !== nothing
const TSRK = ext.TSRK
const TSImplicit = ext.TSImplicit
const TSARKIMEX = ext.TSARKIMEX

# Counters for hook-forwarding tests. The shim methods below are more specific
# than the SciMLBase generics and intercept calls for `PETScTSIntegrator`, so
# tests can count actual invocations of the SciML hook API rather than relying
# on callback lifecycle counters that would pass even if the hooks were removed.
# Each shim increments its counter and then `invoke`s the real SciMLBase generic
# so actual hook behavior (e.g. observable timeseries writes) still runs.
const _save_discretes_hook_count = Ref(0)
const _save_final_discretes_hook_count = Ref(0)
if isdefined(SciMLBase, :save_discretes_if_enabled!) &&
   isdefined(SciMLBase, :save_final_discretes!)
    function SciMLBase.save_discretes_if_enabled!(
        integ::ext.PETScTSIntegrator, cb::SciMLBase.CallbackSet; kw...
    )
        _save_discretes_hook_count[] += 1
        invoke(SciMLBase.save_discretes_if_enabled!,
               Tuple{SciMLBase.DEIntegrator, SciMLBase.CallbackSet},
               integ, cb; kw...)
    end
    function SciMLBase.save_final_discretes!(
        integ::ext.PETScTSIntegrator, cb::SciMLBase.CallbackSet; kw...
    )
        _save_final_discretes_hook_count[] += 1
        invoke(SciMLBase.save_final_discretes!,
               Tuple{SciMLBase.DEIntegrator, SciMLBase.CallbackSet},
               integ, cb; kw...)
    end
end

function decay!(du, u, p, t)
    du[1] = -u[1]
    return nothing
end

@testset "Review-driven fixes" begin
    u0 = [1.0]
    tspan = (0.0, 1.0)
    prob = ODEProblem(decay!, u0, tspan)

    # ── Review-1 #6 / Review-2 #11 ───────────────────────────────────────────
    @testset "Backward / zero-length tspan are rejected with clear errors" begin
        prob_bw = ODEProblem(decay!, [exp(-1.0)], (1.0, 0.0))
        prob_zr = ODEProblem(decay!, u0, (0.0, 0.0))
        @test_throws ArgumentError solve(prob_bw, TSRK("3bs"); dt = 0.1)
        @test_throws ArgumentError solve(prob_zr, TSRK("3bs"); dt = 0.1)
    end

    # ── Review-2 #3 ─────────────────────────────────────────────────────────
    @testset "Integrator exposes derivative_discontinuity field" begin
        integ = init(prob, TSRK("3bs"); dt = 0.1)
        @test hasfield(typeof(integ), :derivative_discontinuity)
        @test integ.derivative_discontinuity == false
        # u_modified! and the discontinuity field are independent.
        SciMLBase.u_modified!(integ, true)
        @test integ.u_modified == true
        @test integ.derivative_discontinuity == false
        # Direct write should also work for any SciMLBase code path that
        # touches the field.
        integ.derivative_discontinuity = true
        @test integ.derivative_discontinuity == true
        PETSc.destroy(integ)
    end

    @testset "DiscreteCallback that never fires still completes the solve" begin
        # Review-2 #3 explicitly asks for this regression: SciMLBase reads
        # derivative_discontinuity in both the fires-and-doesn't-fire paths.
        cb = DiscreteCallback(
            (u, t, integ) -> false,        # never fires
            integ -> nothing,
        )
        sol = solve(prob, TSRK("3bs"); dt = 0.1, callback = cb)
        @test sol.retcode == ReturnCode.Success
        @test sol.t[end] ≈ 1.0
    end

    # ── Review-1 #2 / Review-2 #4 ───────────────────────────────────────────
    @testset "Callback initialize that mutates u is propagated to PETSc" begin
        # If the callback initializer rewrites u0, the first PETSc step must
        # start from the rewritten value — not from the original u0.
        function init_cb!(cb, u, t, integ)
            u[1] = 5.0
            SciMLBase.u_modified!(integ, true)
            return nothing
        end
        cb = DiscreteCallback(
            (u, t, integ) -> false,
            integ -> nothing;
            initialize = init_cb!,
        )
        sol = solve(prob, TSRK("3bs"); dt = 0.1, callback = cb)
        @test sol.retcode == ReturnCode.Success
        # u(t) = 5 * exp(-t) at t = 1 is 5/e ≈ 1.84.
        @test sol.u[end][1] ≈ 5 * exp(-1.0) atol = 1e-2
    end

    # ── Review-1 #5 / Review-2 #6 ───────────────────────────────────────────
    @testset "SciMLBase-equivalent finalize! is called at end of solve" begin
        finalized = Ref(false)
        function finalize_cb!(cb, u, t, integ)
            finalized[] = true
            return nothing
        end
        cb = DiscreteCallback(
            (u, t, integ) -> false,
            integ -> nothing;
            finalize = finalize_cb!,
        )
        sol = solve(prob, TSRK("3bs"); dt = 0.1, callback = cb)
        @test sol.retcode == ReturnCode.Success
        @test finalized[]
    end

    # ── Review-1 #3 / Review-2 #5 ───────────────────────────────────────────
    @testset "reltol / abstol reach PETSc adaptive controller" begin
        # Same problem, two different tolerance settings: at coarse tolerance
        # PETSc takes fewer adaptive steps than at fine tolerance.
        prob_long = ODEProblem(decay!, u0, (0.0, 5.0))
        sol_loose = solve(prob_long, TSRK("5dp"); dt = 0.1, reltol = 1e-2, abstol = 1e-2,
                          save_everystep = true)
        sol_tight = solve(prob_long, TSRK("5dp"); dt = 0.1, reltol = 1e-10, abstol = 1e-10,
                          save_everystep = true)
        @test sol_loose.retcode == ReturnCode.Success
        @test sol_tight.retcode == ReturnCode.Success
        @test length(sol_tight.t) > length(sol_loose.t)
    end

    @testset "Vector tolerances raise ArgumentError" begin
        @test_throws ArgumentError solve(
            prob, TSRK("5dp");
            dt = 0.1, reltol = [1e-6], abstol = [1e-6],
        )
    end

    # ── Review-2 #10 ────────────────────────────────────────────────────────
    @testset "Failed __init does not leak PETSc objects (subsequent solves work)" begin
        @test_throws ArgumentError solve(prob, TSImplicit("does-not-exist"); dt = 0.1)
        # The failed solve allocated a TS and a Vec via _common_ts_setup, which
        # the catch-block in __init must have destroyed before rethrowing. The
        # next valid solve in the same Julia session should succeed.
        sol = solve(prob, TSRK("3bs"); dt = 0.1)
        @test sol.retcode == ReturnCode.Success
    end

    # ── Review-3 #2 ─────────────────────────────────────────────────────────
    @testset "Complex-valued ODEProblem is rejected with a clear error" begin
        prob_c = ODEProblem(decay!, ComplexF64[1 + 0im], tspan)
        err = try
            solve(prob_c, TSRK("3bs"); dt = 0.1)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("real-valued", err.msg)
        @test occursin("ComplexF64", err.msg) || occursin("Complex", err.msg)
    end

    # ── Review-3 #3 ─────────────────────────────────────────────────────────
    @testset "Initialize callback that mutates u does not duplicate t0" begin
        function init_cb!(cb, u, t, integ)
            u[1] = 5.0
            SciMLBase.u_modified!(integ, true)
            return nothing
        end
        cb = DiscreteCallback(
            (u, t, integ) -> false,
            integ -> nothing;
            initialize = init_cb!,
        )
        sol = solve(
            prob, TSRK("3bs");
            dt = 0.1, callback = cb,
            save_everystep = true,           # forces a full trajectory record
        )
        @test sol.retcode == ReturnCode.Success
        # The trajectory has no duplicate timestamps and starts from the
        # mutated initial state.
        @test length(unique(sol.t)) == length(sol.t)
        @test sol.t[1] ≈ 0.0
        @test sol.u[1][1] ≈ 5.0
    end

    # ── Review-4 #1 ─────────────────────────────────────────────────────────
    @testset "Initialize that mutates u without u_modified! still propagates" begin
        # The pessimistic-modified contract: a callback that mutates `u` but
        # forgets to call `SciMLBase.u_modified!(integ, true)` must still
        # affect the first PETSc step. This is what OrdinaryDiffEq does.
        function init_cb!(cb, u, t, integ)
            u[1] = 5.0
            return nothing
        end
        cb = DiscreteCallback(
            (u, t, integ) -> false,
            integ -> nothing;
            initialize = init_cb!,
        )
        sol = solve(prob, TSRK("3bs"); dt = 0.1, callback = cb)
        @test sol.retcode == ReturnCode.Success
        # u(t) = 5 * exp(-t) at t = 1 is 5/e ≈ 1.84.
        @test sol.u[end][1] ≈ 5 * exp(-1.0) atol = 1e-2
    end

    @testset "Initialize-time save records t0 when save_start=false but cb wants it" begin
        # A DiscreteCallback with default save_positions = (true, true) and
        # save_start = false, save_end = false should still record exactly one
        # t0 entry, matching upstream OrdinaryDiffEq behavior.
        cb = DiscreteCallback(
            (u, t, integ) -> false,
            integ -> nothing,
        )
        sol = solve(
            prob, TSRK("3bs");
            dt = 0.1, callback = cb,
            save_start = false, save_end = false, save_everystep = false,
        )
        @test sol.retcode == ReturnCode.Success
        @test length(sol.t) == 1
        @test sol.t[1] ≈ 0.0
        @test sol.u[1][1] ≈ 1.0

        # And with save_positions = (false, false) the initialize-time save
        # must NOT happen — the user explicitly opted out.
        cb_nosave = DiscreteCallback(
            (u, t, integ) -> false,
            integ -> nothing;
            save_positions = (false, false),
        )
        sol2 = solve(
            prob, TSRK("3bs");
            dt = 0.1, callback = cb_nosave,
            save_start = false, save_end = false, save_everystep = false,
        )
        @test sol2.retcode == ReturnCode.Success
        @test isempty(sol2.t)
    end

    @testset "save_start=true with initialize-saving cb still records t0 only once" begin
        cb = DiscreteCallback(
            (u, t, integ) -> false,
            integ -> nothing,
        )
        sol = solve(
            prob, TSRK("3bs");
            dt = 0.1, callback = cb,
            save_start = true, save_end = false, save_everystep = true,
        )
        @test sol.retcode == ReturnCode.Success
        @test count(==(0.0), sol.t) == 1
    end

    @testset "save_on = false suppresses all trajectory output even with init mutation" begin
        function init_cb!(cb, u, t, integ)
            u[1] = 5.0
            SciMLBase.u_modified!(integ, true)
            return nothing
        end
        cb = DiscreteCallback(
            (u, t, integ) -> false,
            integ -> nothing;
            initialize = init_cb!,
        )
        sol = solve(
            prob, TSRK("3bs");
            dt = 0.1, callback = cb,
            save_on = false, save_start = false, save_end = false,
        )
        @test sol.retcode == ReturnCode.Success
        @test isempty(sol.t) && isempty(sol.u)
    end

    # ── Review-2 #12 ────────────────────────────────────────────────────────
    @testset "AbstractVector petsc_options constructors work for every alg" begin
        # Tuples, SubStrings, generic AbstractVectors should all coerce.
        @test PETSc.TSRK("3bs", ("-ts_max_steps", "100")).petsc_options == ["-ts_max_steps", "100"]
        @test PETSc.TSRosW("ra34pw2", ["-snes_fd"]).petsc_options == ["-snes_fd"]
        @test PETSc.TSImplicit("beuler", ("-snes_fd",)).petsc_options == ["-snes_fd"]
        @test PETSc.TSImplicit("theta", 0.5, ("-snes_fd",)).petsc_options == ["-snes_fd"]
        @test PETSc.TSARKIMEX("2e", ("-snes_fd",)).petsc_options == ["-snes_fd"]
        @test PETSc.TSGeneric("alpha", ("-snes_fd",)).petsc_options == ["-snes_fd"]
        # SubString round-trip via split.
        opts_split = split("-snes_fd -ts_max_steps 100")
        @test PETSc.TSRK("3bs", opts_split).petsc_options ==
              ["-snes_fd", "-ts_max_steps", "100"]
    end

    # ── Review-5 #1 ─────────────────────────────────────────────────────────
    @testset "Single-sided reltol still reaches PETSc" begin
        # `solve(...; reltol = 1e-10)` (no `abstol`) used to be silently
        # ignored — the wrapper required *both* sides to be set. Now the
        # missing side is filled with SciML's default and PETSc's adaptive
        # controller actually responds.
        prob_long = ODEProblem(decay!, u0, (0.0, 5.0))
        sol_default = solve(prob_long, TSRK("5dp"); dt = 0.1, save_everystep = true)
        sol_rel = solve(
            prob_long, TSRK("5dp");
            dt = 0.1, reltol = 1e-10, save_everystep = true,
        )
        @test length(sol_rel.t) > length(sol_default.t)
    end

    @testset "Single-sided abstol still reaches PETSc" begin
        # The decay problem decays exponentially toward zero, so the
        # `atol + rtol * |u|` threshold becomes dominated by `atol` once `u`
        # is small. Tightening `abstol` alone (with `reltol` defaulted)
        # therefore changes the adaptive step count — provided abstol is
        # actually forwarded to PETSc, which used to require both sides set.
        prob_decay = ODEProblem(decay!, u0, (0.0, 30.0))
        sol_default = solve(
            prob_decay, TSRK("5dp"); dt = 0.1, save_everystep = true,
        )
        sol_abs = solve(
            prob_decay, TSRK("5dp");
            dt = 0.1, abstol = 1e-14, save_everystep = true,
        )
        @test length(sol_abs.t) > length(sol_default.t)
    end

    # ── Review-5 #2 ─────────────────────────────────────────────────────────
    @testset "Unsupported solve keywords are rejected with ArgumentError" begin
        # Anything not on the explicit allowlist should fail loudly. Pick a
        # set of common SciML knobs that this extension does NOT honour.
        for bad in (:progress, :progress_steps, :alias_u0, :internalnorm,
                    :force_dtmin, :unstable_check)
            err = try
                solve(prob, TSRK("3bs"); dt = 0.1, (; bad => true)...)
                nothing
            catch e
                e
            end
            @test err isa ArgumentError
            @test occursin(string(bad), err.msg)
        end
    end

    @testset "adaptive = false disables PETSc's adaptive controller" begin
        # Fixed-step PETSc takes evenly-spaced dt steps; the default adaptive
        # solver completes in many fewer steps for this trivial problem.
        prob_long = ODEProblem(decay!, u0, (0.0, 1.0))
        sol_adapt = solve(
            prob_long, TSRK("5dp");
            dt = 0.1, save_everystep = true,
        )
        sol_fixed = solve(
            prob_long, TSRK("5dp");
            dt = 0.1, adaptive = false, save_everystep = true,
        )
        @test sol_adapt.retcode == ReturnCode.Success
        @test sol_fixed.retcode == ReturnCode.Success
        # Fixed-step at dt = 0.1 over [0, 1] is 10 internal steps + start/end.
        @test length(sol_fixed.t) >= 10
        # And the adjacent intervals should be (almost) the same fixed dt.
        diffs = diff(sol_fixed.t)
        @test all(d -> isapprox(d, 0.1; atol = 1e-12), diffs[1:(end - 1)])
    end

    @testset "dtmax caps the step size of the adaptive controller" begin
        prob_long = ODEProblem(decay!, u0, (0.0, 5.0))
        sol_capped = solve(
            prob_long, TSRK("5dp");
            dt = 0.1, dtmax = 0.2, save_everystep = true,
        )
        @test sol_capped.retcode == ReturnCode.Success
        # Every internal interval must respect the cap (modulo the final
        # match-step trim, which can be smaller).
        for d in diff(sol_capped.t)
            @test d <= 0.2 + 1e-12
        end
    end

    @testset "TSGeneric explicit = true accepts euler / ssp" begin
        # Previously TSGeneric always registered an IFunction, so explicit
        # PETSc TS types failed at TSStep with a raw PETSc error. With
        # `explicit = true` the RHS path is selected and the solve succeeds.
        prob_short = ODEProblem(decay!, u0, (0.0, 1.0))
        sol_euler = solve(
            prob_short, PETSc.TSGeneric("euler"; explicit = true); dt = 0.1,
        )
        @test sol_euler.retcode == ReturnCode.Success
        @test sol_euler.t[end] ≈ 1.0
        @test sol_euler.u[end][1] ≈ exp(-1.0) atol = 1e-1

        sol_ssp = solve(
            prob_short, PETSc.TSGeneric("ssp"; explicit = true); dt = 0.1,
        )
        @test sol_ssp.retcode == ReturnCode.Success
        @test sol_ssp.t[end] ≈ 1.0
        @test sol_ssp.u[end][1] ≈ exp(-1.0) atol = 1e-1

        # Default `explicit = false` for an explicit-only TS type fails up
        # front with a clear Julia-side `ArgumentError` instead of a noisy
        # PETSc banner from `TSSetUp`.
        for ts_type in ("euler", "ssp")
            err = try
                solve(prob_short, PETSc.TSGeneric(ts_type); dt = 0.1)
                nothing
            catch e
                e
            end
            @test err isa ArgumentError
            @test occursin("explicit = true", err.msg)
            @test occursin(ts_type, err.msg)
        end
    end

    @testset "TSGeneric positional petsc_options work with explicit kwarg" begin
        alg = PETSc.TSGeneric("euler", ["-ts_max_steps", "100"]; explicit = true)
        @test alg.ts_type == "euler"
        @test alg.explicit == true
        @test alg.petsc_options == ["-ts_max_steps", "100"]
    end

    # ── Review-6 #1 ─────────────────────────────────────────────────────────
    @testset "Initial dt outside [dtmin, dtmax] is rejected" begin
        # PETSc installs the initial `dt` verbatim and only consults
        # `TSAdaptSetStepLimits` for subsequent step proposals, so the wrapper
        # has to validate `dt` against the user's bounds itself.
        err = try
            solve(prob, TSRK("5dp"); dt = 0.5, dtmax = 0.2)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("dtmax", err.msg)

        err = try
            solve(prob, TSRK("5dp"); dt = 0.1, dtmin = 0.2)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("dtmin", err.msg)

        # The valid case — `dt` inside the allowed range — still succeeds.
        sol = solve(
            prob, TSRK("5dp");
            dt = 0.1, dtmin = 0.05, dtmax = 0.2, save_everystep = true,
        )
        @test sol.retcode == ReturnCode.Success
        for d in diff(sol.t)
            @test d <= 0.2 + 1e-12
        end
    end

    # ── Review-6 #2 ─────────────────────────────────────────────────────────
    @testset "verbose is no longer in the supported-keyword set" begin
        # Previously `verbose = true` was silently accepted but unused; the
        # extension now treats it like any other unsupported keyword.
        err = try
            solve(prob, TSRK("3bs"); dt = 0.1, verbose = true)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("verbose", err.msg)
    end

    # ── Review-6 #3 ─────────────────────────────────────────────────────────
    @testset "TSGeneric without explicit = true rejects euler / ssp upfront" begin
        # The Julia-side validator must fire before any PETSc setup runs, so
        # the regular test output stays free of raw PETSc error banners.
        for ts_type in ("euler", "ssp")
            err = try
                solve(prob, PETSc.TSGeneric(ts_type); dt = 0.1)
                nothing
            catch e
                e
            end
            @test err isa ArgumentError
            @test occursin("explicit = true", err.msg)
        end
    end

    # ── Review-7 #1 ─────────────────────────────────────────────────────────
    @testset "Invalid dt values are rejected with ArgumentError" begin
        # Negative `dt` would otherwise let a forward solve effectively step
        # backward and still report Success.
        @test_throws ArgumentError solve(prob, TSRK("3bs"); dt = -0.1)
        # Zero `dt` produces a degenerate trajectory under PETSc's adaptive
        # controller — also reject up front.
        @test_throws ArgumentError solve(prob, TSRK("3bs"); dt = 0.0)
        # Non-finite `dt` previously fell through to a raw PETSc banner.
        @test_throws ArgumentError solve(prob, TSRK("3bs"); dt = Inf)
        @test_throws ArgumentError solve(prob, TSRK("3bs"); dt = NaN)
    end

    # ── Review-7 #2 ─────────────────────────────────────────────────────────
    @testset "Invalid dtmin / dtmax values are rejected with ArgumentError" begin
        # Inverted bounds.
        err = try
            solve(prob, TSRK("5dp"); dtmin = 0.2, dtmax = 0.1)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("dtmin", err.msg)
        @test occursin("dtmax", err.msg)

        # Negative bounds.
        @test_throws ArgumentError solve(prob, TSRK("5dp"); dtmin = -0.1)
        @test_throws ArgumentError solve(prob, TSRK("5dp"); dtmax = -0.1)

        # Non-finite bounds.
        @test_throws ArgumentError solve(prob, TSRK("5dp"); dtmin = NaN)
        @test_throws ArgumentError solve(prob, TSRK("5dp"); dtmax = NaN)
        @test_throws ArgumentError solve(prob, TSRK("5dp"); dtmin = Inf)

        # Zero `dtmax` must also be rejected (no positive step would satisfy
        # it). Zero `dtmin` is allowed because it is the default lower bound.
        @test_throws ArgumentError solve(prob, TSRK("5dp"); dtmax = 0.0)
        sol = solve(prob, TSRK("5dp"); dt = 0.1, dtmin = 0.0)
        @test sol.retcode == ReturnCode.Success
    end

    # ── Review-7 #3 ─────────────────────────────────────────────────────────
    @testset "Invalid reltol / abstol values are rejected with ArgumentError" begin
        @test_throws ArgumentError solve(
            prob, TSRK("5dp"); dt = 0.1, reltol = -1e-3, abstol = 1e-6,
        )
        @test_throws ArgumentError solve(
            prob, TSRK("5dp"); dt = 0.1, reltol = 1e-3, abstol = -1e-6,
        )
        @test_throws ArgumentError solve(
            prob, TSRK("5dp"); dt = 0.1, reltol = -1e-3, abstol = -1e-6,
        )
        @test_throws ArgumentError solve(
            prob, TSRK("5dp"); dt = 0.1, reltol = NaN,
        )
        @test_throws ArgumentError solve(
            prob, TSRK("5dp"); dt = 0.1, abstol = Inf,
        )
        # Zero tolerances are allowed (matching PETSc's policy: zero on one
        # side simply disables that part of the `atol + rtol * |u|` test).
        sol = solve(prob, TSRK("5dp"); dt = 0.1, reltol = 0.0, abstol = 1e-6)
        @test sol.retcode == ReturnCode.Success
    end

    # ── Review-8 #1 ─────────────────────────────────────────────────────────
    @testset "Invalid scalar saveat values are rejected" begin
        # Previously these silently produced an empty save schedule with a
        # `Success` retcode. Now they fail loudly at the Julia boundary.
        @test_throws ArgumentError solve(
            prob, TSRK("3bs"); dt = 0.1, saveat = 0.0,
        )
        @test_throws ArgumentError solve(
            prob, TSRK("3bs"); dt = 0.1, saveat = -0.1,
        )
        @test_throws ArgumentError solve(
            prob, TSRK("3bs"); dt = 0.1, saveat = Inf,
        )
        @test_throws ArgumentError solve(
            prob, TSRK("3bs"); dt = 0.1, saveat = NaN,
        )
    end

    @testset "Iterable saveat with non-finite entries is rejected" begin
        @test_throws ArgumentError solve(
            prob, TSRK("3bs"); dt = 0.1, saveat = [0.25, NaN, 0.75],
        )
        @test_throws ArgumentError solve(
            prob, TSRK("3bs"); dt = 0.1, saveat = (0.25, Inf),
        )
    end

    # ── Review-9 #1 ─────────────────────────────────────────────────────────
    @testset "maxiters caps the manual TSStep loop and reports MaxIters" begin
        # Fixed-step integration over [0, 1] with `dt = 0.1` would otherwise
        # run for 10 steps; `maxiters = 1` must stop it after one accepted
        # step and surface a `ReturnCode.MaxIters`.
        sol = solve(
            prob, TSRK("3bs");
            dt = 0.1, adaptive = false, maxiters = 1, save_everystep = true,
        )
        @test sol.retcode == ReturnCode.MaxIters
        @test sol.t[end] < 1.0
        @test sol.stats.naccept >= 1
    end

    @testset "PETSc-side -ts_max_steps is honoured by the step-count cap" begin
        # The cap is enforced as `min(opts.maxiters, TSGetMaxSteps(ts))`,
        # so passing `-ts_max_steps 1` through `petsc_options` produces the
        # same stop-early behaviour as the SciML `maxiters` knob.
        alg = PETSc.TSRK("3bs", ["-ts_max_steps", "1"])
        sol = solve(prob, alg; dt = 0.1, adaptive = false, save_everystep = true)
        @test sol.retcode == ReturnCode.MaxIters
        @test sol.t[end] ≈ 0.1
        @test sol.stats.naccept == 1
    end

    # ── Review-9 #2 ─────────────────────────────────────────────────────────
    @testset "sol.stats reflects the actual number of steps taken" begin
        sol = solve(
            prob, TSRK("3bs");
            dt = 0.1, adaptive = false, save_everystep = true,
        )
        @test sol.retcode == ReturnCode.Success
        # Fixed-step at dt = 0.1 over [0, 1] takes 10 accepted steps.
        @test sol.stats.naccept == 10
        @test sol.stats.nreject == 0
        # Explicit RK has no SNES work.
        @test sol.stats.nnonliniter == 0
        # `nf` is incremented in the C-callback; an explicit RK must call
        # the user RHS multiple times per step (3bs has 3 stages plus FSAL),
        # so the counter must be strictly positive.
        @test sol.stats.nf > sol.stats.naccept
        # `nsolve` is left at the SciML "unknown" sentinel because
        # `TSGetKSPIterations` reports linear iterations, not solves.
        @test sol.stats.nsolve == -1
    end

    @testset "Implicit solve populates SNES iteration count in stats" begin
        sol = solve(
            prob, TSImplicit("beuler", ["-snes_fd"]);
            dt = 0.1, adaptive = false,
        )
        @test sol.retcode == ReturnCode.Success
        @test sol.stats.naccept == 10
        # An implicit method must do at least one nonlinear solve per step.
        @test sol.stats.nnonliniter > 0
        @test sol.stats.nf > 0
    end

    # ── Review-10 #1 ────────────────────────────────────────────────────────
    @testset "maxiters = 0 yields a zero-step solve, not one step" begin
        # The cap is now checked *before* `TSStep`, so `maxiters = 0` must
        # not advance time, regardless of what `dt` would otherwise produce.
        sol = solve(
            prob, TSRK("3bs");
            dt = 0.1, adaptive = false, maxiters = 0, save_everystep = true,
        )
        @test sol.retcode == ReturnCode.MaxIters
        @test sol.t[end] == 0.0
        @test sol.stats.naccept == 0

        # Even when `dt` is large enough that one step would land on `tf`,
        # the cap must still pre-empt the step and report MaxIters.
        sol_one_shot = solve(
            prob, TSRK("3bs");
            dt = 1.0, adaptive = false, maxiters = 0, save_everystep = true,
        )
        @test sol_one_shot.retcode == ReturnCode.MaxIters
        @test sol_one_shot.t[end] == 0.0
        @test sol_one_shot.stats.naccept == 0
    end

    @testset "petsc_options -ts_max_steps 0 yields a zero-step solve" begin
        alg = PETSc.TSRK("3bs", ["-ts_max_steps", "0"])
        sol = solve(prob, alg; dt = 0.1, adaptive = false, save_everystep = true)
        @test sol.retcode == ReturnCode.MaxIters
        @test sol.t[end] == 0.0
        @test sol.stats.naccept == 0
    end

    # ── Review-14 #1 ────────────────────────────────────────────────────────
    @testset "Front-end validation runs before any PETSc TS allocation" begin
        # Both `tstops` and bad `saveat` should now fail before
        # `_common_ts_setup` allocates a `TS` / `Vec` pair, so repeated
        # invalid solves cannot accumulate live PETSc state.
        for _ in 1:5
            @test_throws ArgumentError solve(
                prob, TSRK("3bs"); dt = 0.1, tstops = [0.4],
            )
            @test_throws ArgumentError solve(
                prob, TSRK("3bs"); dt = 0.1, saveat = NaN,
            )
        end
        # Subsequent valid solve still succeeds — proves nothing leaked
        # leaves the PETSc state in a bad shape.
        sol = solve(prob, TSRK("3bs"); dt = 0.1)
        @test sol.retcode == ReturnCode.Success
    end

    # ── Review-14 #2 ────────────────────────────────────────────────────────
    @testset "Callback initialize / finalize hooks fire around the solve" begin
        # The lifecycle calls (`SciMLBase-equivalent initialize!` /
        # `SciMLBase-equivalent finalize!`) are part of the basic discrete-callback
        # contract and must run regardless of SciMLBase version.
        init_runs = Ref(0)
        finalize_runs = Ref(0)
        cb = DiscreteCallback(
            (u, t, integ) -> false,
            integ -> nothing;
            initialize = (cb, u, t, integ) -> (init_runs[] += 1; nothing),
            finalize = (cb, u, t, integ) -> (finalize_runs[] += 1; nothing),
        )
        sol = solve(prob, TSRK("3bs"); dt = 0.1, callback = cb)
        @test sol.retcode == ReturnCode.Success
        @test init_runs[] == 1
        @test finalize_runs[] == 1
    end

    # The new discrete-save lifecycle hooks
    # (`SciMLBase.save_discretes_if_enabled!` / `save_final_discretes!`)
    # only exist in newer SciMLBase releases. The wrapper guards them with
    # `isdefined`; mirror that gate here so the test file describes what
    # is actually being verified on the loaded compat set.
    if isdefined(SciMLBase, :save_discretes_if_enabled!) &&
       isdefined(SciMLBase, :save_final_discretes!)
        @testset "SciMLBase discrete-save lifecycle hooks are forwarded" begin
            # Reset the shim counters defined at file scope. The shims are more
            # specific than the SciMLBase generics, so they intercept the exact
            # call sites in solve.jl and let us pin the forwarding path directly.
            _save_discretes_hook_count[] = 0
            _save_final_discretes_hook_count[] = 0
            sol = solve(prob, TSRK("3bs"); dt = 0.1)
            @test sol.retcode == ReturnCode.Success
            # `save_discretes_if_enabled!` is called once during initialization;
            # `save_final_discretes!` is called once at the end of `solve!`.
            @test _save_discretes_hook_count[] == 1
            @test _save_final_discretes_hook_count[] == 1
        end
    else
        @info "Skipping `SciMLBase discrete-save lifecycle hooks are forwarded`: " *
              "loaded SciMLBase does not expose `save_discretes_if_enabled!` / " *
              "`save_final_discretes!`."
    end

    # ── Review-15 #1 ────────────────────────────────────────────────────────
    @testset "initialize_save = false suppresses the post-init save record" begin
        # `initialize_save` is the upstream SciML knob that controls
        # whether the integrator appends a save record immediately after
        # callback `initialize!` runs. With it `false`, the `initialize`
        # hook still fires but no `t0` row is added on top of the regular
        # `save_start` row. With it `true` (the default), the row is added
        # exactly once even on top of `save_start = true`.
        init_ran = Ref(false)
        # A discrete callback with default `save_positions = (true, true)`
        # asks for a post-init save row. Silence everything else so we can
        # observe whether the post-init save fired.
        cb = DiscreteCallback(
            (u, t, integ) -> false,
            integ -> nothing;
            initialize = (cb, u, t, integ) -> (init_ran[] = true; nothing),
        )
        sol_off = solve(
            prob, TSRK("3bs");
            dt = 0.1, callback = cb,
            save_start = false, save_end = false, save_everystep = false,
            initialize_save = false,
        )
        @test sol_off.retcode == ReturnCode.Success
        @test init_ran[]              # the `initialize` hook still fires
        @test isempty(sol_off.t)      # but the post-init save was suppressed

        # Default `initialize_save = true` keeps the post-init save row.
        init_ran[] = false
        sol_on = solve(
            prob, TSRK("3bs");
            dt = 0.1, callback = cb,
            save_start = false, save_end = false, save_everystep = false,
        )
        @test sol_on.retcode == ReturnCode.Success
        @test init_ran[]
        @test length(sol_on.t) == 1
        @test sol_on.t[1] ≈ 0.0
    end

    # ── Review-16 #1 ────────────────────────────────────────────────────────
    # `save_discretes` was added to SciMLBase's `allowedkeywords` whitelist in
    # 2.120.0. On older stacks the kwarg is rejected by SciMLBase before it
    # reaches the PETSc extension, so we gate the test on that floor.
    if pkgversion(SciMLBase) >= v"2.120.0"
        @testset "save_discretes is accepted and stored in integrator opts" begin
            # `save_discretes` controls whether `_apply_discrete_callback!`
            # records discrete observable state after a callback fires. The wrapper
            # must store it in `DEOptions` so the integrator can observe it through the
            # standard `integrator.opts.save_discretes` path.
            integ_on = init(prob, TSRK("3bs"); dt = 0.1, save_discretes = true)
            @test integ_on.opts.save_discretes == true
            PETSc.destroy(integ_on)

            integ_off = init(prob, TSRK("3bs"); dt = 0.1, save_discretes = false)
            @test integ_off.opts.save_discretes == false
            PETSc.destroy(integ_off)
        end
    else
        @info "Skipping `save_discretes` test: loaded SciMLBase $(pkgversion(SciMLBase)) " *
              "is older than 2.120.0 (the first release to include :save_discretes in " *
              "allowedkeywords)."
    end

    # ── Review-13 #1 ────────────────────────────────────────────────────────
    @testset "tstops kwarg is rejected with ArgumentError" begin
        # `tstops` carries a strict SciML contract — the integrator must
        # land on those times so step-end callback logic can see them.
        # Silently ignoring it would skip exact-time discrete callbacks,
        # so the wrapper rejects it up front until the manual `TSStep`
        # loop honours it natively.
        err = try
            solve(prob, TSRK("3bs"); dt = 0.1, tstops = [0.4, 0.6])
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("tstops", err.msg)
    end

    # ── Review-13 #2 ────────────────────────────────────────────────────────
    @testset "Stateful saveat iterators survive validation" begin
        # `Iterators.Stateful` is a one-shot iterator: a previous version
        # of the wrapper validated `saveat` by iterating it once, which
        # consumed the iterator before `_expand_saveat` could read it
        # again. We now materialize iterables to a `Vector` exactly once
        # so this case round-trips correctly.
        saveit = Iterators.Stateful([0.25, 0.5, 0.75])
        sol = solve(
            prob, TSRK("3bs");
            dt = 0.1, saveat = saveit,
            save_start = false, save_end = false,
        )
        @test sol.retcode == ReturnCode.Success
        @test all(t -> any(s -> isapprox(s, t; atol = 1e-12), sol.t),
                  (0.25, 0.5, 0.75))
    end

    # ── Review-12 #1 ────────────────────────────────────────────────────────
    @testset "TSARKIMEX with SplitODEProblem populates both nf and nf2" begin
        # The implicit (`f1`) stream should land in `stats.nf`; the
        # explicit (`f2`) stream should land in `stats.nf2`. We pin the
        # mapping by counting `f1` and `f2` calls in user closures and
        # asserting `stats.nf == f1_calls[]` and `stats.nf2 == f2_calls[]`,
        # so a future regression that swaps the two streams is caught.
        f1_calls = Ref(0)
        f2_calls = Ref(0)
        f1!(du, u, p, t) = (f1_calls[] += 1; du[1] = -u[1]; nothing)
        f2!(du, u, p, t) = (f2_calls[] += 1; du[1] = cos(t); nothing)
        prob_split = SplitODEProblem(f1!, f2!, [1.0], (0.0, 1.0))
        sol = solve(prob_split, TSARKIMEX("2e", ["-snes_fd"]); dt = 0.1, adaptive = false)
        @test sol.retcode == ReturnCode.Success
        @test sol.stats.naccept == 10
        @test f1_calls[] > 0
        @test f2_calls[] > 0
        # The promised mapping: `nf` <- f1 (implicit), `nf2` <- f2 (explicit).
        @test sol.stats.nf == f1_calls[]
        @test sol.stats.nf2 == f2_calls[]
        # And the two streams must produce different totals here, so a
        # broken implementation that wrote the same value into both fields
        # would be caught.
        @test f1_calls[] != f2_calls[]
    end

    # ── Review-11 #1 ────────────────────────────────────────────────────────
    @testset "Negative maxiters is rejected with ArgumentError" begin
        # Without explicit validation, a negative SciML `maxiters` would
        # turn the manual `TSStep` loop into a zero-step `MaxIters` solve —
        # but PETSc's `-ts_max_steps -1` means "unlimited", so the two
        # spellings would disagree. Reject the negative SciML form so the
        # asymmetry can never bite.
        err = try
            solve(prob, TSRK("3bs"); dt = 0.1, maxiters = -1)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("maxiters", err.msg)

        # `maxiters = 0` remains a valid zero-step solve (Review-10).
        sol = solve(prob, TSRK("3bs"); dt = 0.1, adaptive = false, maxiters = 0)
        @test sol.retcode == ReturnCode.MaxIters
        @test sol.stats.naccept == 0
    end

    # ── Review-10 #2 ────────────────────────────────────────────────────────
    @testset "Unsupported DEStats fields stay at the SciML \"unknown\" sentinel" begin
        sol = solve(prob, TSRK("3bs"); dt = 0.1, adaptive = false)
        # We populate `naccept`, `nreject`, `nnonliniter`, and `nf`. Every
        # other counter should remain at `-1` so users can distinguish
        # "we don't track this" from "no work happened".
        @test sol.stats.naccept == 10
        @test sol.stats.nreject == 0
        @test sol.stats.nnonliniter == 0
        @test sol.stats.nf > 0
        @test sol.stats.nsolve == -1
        @test sol.stats.nf2 == -1
        @test sol.stats.nw == -1
        @test sol.stats.njacs == -1
    end

    # ── Review-8 #2 ─────────────────────────────────────────────────────────
    @testset "dtmax = Inf is accepted and behaves like an omitted dtmax" begin
        # `Inf` is the natural SciML spelling of "no upper cap", and the
        # wrapper already uses `Inf` as the default internally. Reject the
        # spelling-asymmetry the previous validator introduced.
        sol_inf = solve(
            prob, TSRK("5dp"); dt = 0.1, dtmax = Inf, save_everystep = true,
        )
        sol_nothing = solve(
            prob, TSRK("5dp"); dt = 0.1, save_everystep = true,
        )
        @test sol_inf.retcode == ReturnCode.Success
        @test sol_nothing.retcode == ReturnCode.Success
        @test length(sol_inf.t) == length(sol_nothing.t)

        # `dt = Inf` and `dtmin = Inf` remain rejected — only `dtmax = Inf`
        # has the "no cap" meaning, so only it is special-cased.
        @test_throws ArgumentError solve(prob, TSRK("5dp"); dt = Inf)
        @test_throws ArgumentError solve(prob, TSRK("5dp"); dtmin = Inf)
    end
end
