using Test
using PETSc
using SciMLBase
using DiffEqBase

ext = Base.get_extension(PETSc, :PETScSciMLExt)
@assert ext !== nothing
const TSRK = ext.TSRK
const TSImplicit = ext.TSImplicit

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
        DiffEqBase.u_modified!(integ, true)
        @test integ.u_modified == true
        @test integ.derivative_discontinuity == false
        # Direct write should also work for any SciMLBase code path that
        # touches the field.
        integ.derivative_discontinuity = true
        @test integ.derivative_discontinuity == true
        PETSc.destroy(integ)
    end

    @testset "DiscreteCallback that never fires still completes the solve" begin
        # Review-2 #3 explicitly asks for this regression: DiffEqBase reads
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
            DiffEqBase.u_modified!(integ, true)
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
    @testset "DiffEqBase.finalize! is called at end of solve" begin
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
            DiffEqBase.u_modified!(integ, true)
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
        # forgets to call `DiffEqBase.u_modified!(integ, true)` must still
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
            save_start = false, save_end = false,
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
            save_start = false, save_end = false,
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
            DiffEqBase.u_modified!(integ, true)
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

        # Default `explicit = false` for an explicit-only TS type still
        # fails — but with a PETSc-side error, not silently. Catch any
        # exception (PETSc emits a plain `ErrorException` from `@chk`).
        @test_throws Exception solve(
            prob_short, PETSc.TSGeneric("euler"); dt = 0.1,
        )
    end

    @testset "TSGeneric positional petsc_options work with explicit kwarg" begin
        alg = PETSc.TSGeneric("euler", ["-ts_max_steps", "100"]; explicit = true)
        @test alg.ts_type == "euler"
        @test alg.explicit == true
        @test alg.petsc_options == ["-ts_max_steps", "100"]
    end
end
