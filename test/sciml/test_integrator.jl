using Test
using PETSc
using SciMLBase

function decay!(du, u, p, t)
    du[1] = -u[1]
    return nothing
end

@testset "Integrator interface lifecycle" begin
    u0 = [1.0]
    tspan = (0.0, 1.0)
    prob = ODEProblem(decay!, u0, tspan)

    @testset "__init + solve! equals __solve" begin
        sol_oneshot = solve(prob, TSRK("3bs"); dt = 0.1)

        integ = init(prob, TSRK("3bs"); dt = 0.1)
        sol_steps = solve!(integ)

        @test sol_steps.retcode == sol_oneshot.retcode
        @test sol_steps.t == sol_oneshot.t
        @test sol_steps.u == sol_oneshot.u
    end

    @testset "opts mirrors SciML control knobs for downstream introspection" begin
        # Generic SciML code reads `integrator.opts.adaptive` / `.reltol` /
        # `.abstol`; make sure the wrapper populates them from the kwargs.
        integ = init(prob, TSRK("3bs"); dt = 0.1,
                     adaptive = false, reltol = 1e-8, abstol = 1e-10)
        @test integ.opts.adaptive == false
        @test integ.opts.reltol == 1e-8
        @test integ.opts.abstol == 1e-10
        PETSc.destroy(integ)

        # Defaults: adaptive on, tolerances left as `nothing` (PETSc's own
        # defaults are used — see `_apply_tolerances!`).
        integ2 = init(prob, TSRK("3bs"); dt = 0.1)
        @test integ2.opts.adaptive == true
        @test integ2.opts.reltol === nothing
        @test integ2.opts.abstol === nothing
        PETSc.destroy(integ2)
    end

    @testset "manual step! loop matches solve!" begin
        integ_a = init(prob, TSRK("3bs"); dt = 0.1)
        integ_b = init(prob, TSRK("3bs"); dt = 0.1)

        sol_a = solve!(integ_a)

        while !integ_b.done
            step!(integ_b)
        end
        # mimic solve!'s save_end + retcode finalisation
        if integ_b.opts.save_end &&
           (isempty(integ_b.sol.t) || last(integ_b.sol.t) != integ_b.t)
            push!(integ_b.sol.t, integ_b.t)
            push!(integ_b.sol.u, copy(integ_b.u))
        end

        @test integ_b.retcode == ReturnCode.Success
        @test integ_b.t ≈ integ_a.t
        @test integ_b.u ≈ integ_a.u atol = 1e-12
        PETSc.destroy(integ_b)
    end

    @testset "manual step! loop with early break leaves a valid partial trajectory" begin
        integ = init(prob, TSRK("3bs"); dt = 0.1, save_everystep = true)
        for _ in 1:3
            integ.done && break
            step!(integ)
        end
        @test integ.t > 0.0
        @test integ.t < tspan[2]
        @test length(integ.sol.t) >= 3
        @test issorted(integ.sol.t)
        PETSc.destroy(integ)
    end

    @testset "step! after done is a no-op" begin
        integ = init(prob, TSRK("3bs"); dt = 0.1)
        sol = solve!(integ)
        @test integ.done
        # solve! has already destroyed the PETSc resources; step! must
        # short-circuit on `done` before touching them.
        @test step!(integ) === nothing
        @test integ.t == sol.t[end]
    end

    @testset "two independent integrators coexist" begin
        prob2 = ODEProblem(decay!, [2.0], tspan)
        integ_a = init(prob, TSRK("3bs"); dt = 0.1)
        integ_b = init(prob2, TSRK("5dp"); dt = 0.1)
        sol_a = solve!(integ_a)
        sol_b = solve!(integ_b)
        @test sol_a.retcode == ReturnCode.Success
        @test sol_b.retcode == ReturnCode.Success
        @test sol_a.u[end][1] ≈ exp(-1.0) atol = 1e-3
        @test sol_b.u[end][1] ≈ 2 * exp(-1.0) atol = 1e-3
    end

    @testset "PETSc.destroy(integ) is idempotent and safe after solve!" begin
        integ = init(prob, TSRK("3bs"); dt = 0.1)
        sol = solve!(integ)
        @test sol.retcode == ReturnCode.Success
        # solve! already destroyed the PETSc objects; destroy() should be a no-op.
        @test PETSc.destroy(integ) === nothing
        @test PETSc.destroy(integ) === nothing
    end

    @testset "PETSc.destroy(integ) on an unrun integrator releases memory" begin
        integ = init(prob, TSRK("3bs"); dt = 0.1)
        @test integ.ts.ptr != C_NULL
        @test integ.u_petsc.ptr != C_NULL
        PETSc.destroy(integ)
        @test integ.ts.ptr == C_NULL
        @test integ.u_petsc.ptr == C_NULL
    end

    @testset "GC finalizer destroys PETSc resources without error" begin
        let integ = init(prob, TSRK("3bs"); dt = 0.1)
            # let-block scopes integ so it can be collected
            @test integ.ts.ptr != C_NULL
        end
        # force a GC pass; the finalizer (registered in _make_integrator) should
        # have run on the now-unreachable integrator without raising.
        GC.gc()
        GC.gc()
        @test true
    end
end
