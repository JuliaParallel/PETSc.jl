using Test
using PETSc
using SciMLBase

# ── Reference problem (PETSc TS tutorial ex51) ───────────────────────────────
#   u1' = cos(t),   u2' = sin(u2)
# Analytical solution:
#   u1(t) = sin(t),         u2(t) = 2 atan(exp(t) tan(0.5))
function ex51_rhs!(du, u, p, t)
    du[1] = cos(t)
    du[2] = sin(u[2])
    return nothing
end

ex51_exact(t) = [sin(t), 2 * atan(exp(t) * tan(0.5))]

@testset "TSRK explicit Runge-Kutta" begin
    @testset "ex51 with TSRK(\"3bs\") matches analytical solution" begin
        u0 = [0.0, 1.0]
        tspan = (0.0, 1.0)
        prob = ODEProblem(ex51_rhs!, u0, tspan)
        sol = solve(prob, TSRK("3bs"); dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        @test length(sol.t) >= 2
        @test sol.t[end] ≈ tspan[2]
        @test sol.u[end] ≈ ex51_exact(tspan[2]) atol = 1e-3
    end

    @testset "ex51 with TSRK(\"5dp\") matches with tighter tolerance" begin
        u0 = [0.0, 1.0]
        tspan = (0.0, 1.0)
        prob = ODEProblem(ex51_rhs!, u0, tspan)
        sol = solve(prob, TSRK("5dp"); dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        @test sol.u[end] ≈ ex51_exact(tspan[2]) atol = 1e-3
    end

    @testset "petsc_options on the algorithm overrides the subtype" begin
        u0 = [0.0, 1.0]
        tspan = (0.0, 1.0)
        prob = ODEProblem(ex51_rhs!, u0, tspan)
        sol_direct = solve(prob, TSRK("5dp"); dt = 0.1)
        sol_via_opts = solve(
            prob,
            TSRK("3bs", ["-ts_rk_type", "5dp"]);
            dt = 0.1,
        )
        @test sol_via_opts.retcode == ReturnCode.Success
        @test sol_via_opts.u[end] ≈ sol_direct.u[end] atol = 1e-10
    end

    @testset "scalar exponential decay u' = -u" begin
        u0 = [1.0]
        tspan = (0.0, 1.0)
        prob = ODEProblem((du, u, p, t) -> (du[1] = -u[1]), u0, tspan)
        sol = solve(prob, TSRK("3bs"); dt = 0.05)
        @test sol.retcode == ReturnCode.Success
        @test sol.u[end][1] ≈ exp(-1.0) atol = 1e-3
    end

    @testset "out-of-place ODEProblem is rejected with informative error" begin
        u0 = [1.0]
        tspan = (0.0, 1.0)
        prob = ODEProblem((u, p, t) -> -u, u0, tspan)
        @test_throws ArgumentError solve(prob, TSRK("3bs"); dt = 0.1)
    end

    @testset "init / step! / solve! produce the same result as solve" begin
        u0 = [0.0, 1.0]
        tspan = (0.0, 1.0)
        prob = ODEProblem(ex51_rhs!, u0, tspan)
        sol_oneshot = solve(prob, TSRK("5dp"); dt = 0.1)

        integ = init(prob, TSRK("5dp"); dt = 0.1)
        sol_steps = solve!(integ)
        @test sol_steps.retcode == ReturnCode.Success
        @test sol_steps.u[end] ≈ sol_oneshot.u[end] atol = 1e-12
    end

    @testset "manual step! loop reaches tspan[2]" begin
        u0 = [0.0, 1.0]
        tspan = (0.0, 1.0)
        prob = ODEProblem(ex51_rhs!, u0, tspan)
        integ = init(prob, TSRK("5dp"); dt = 0.1)
        steps = 0
        while !integ.done && steps < 1000
            step!(integ)
            steps += 1
        end
        @test integ.retcode == ReturnCode.Success
        @test integ.t ≈ tspan[2]
        @test integ.u ≈ ex51_exact(tspan[2]) atol = 1e-3
    end
end
