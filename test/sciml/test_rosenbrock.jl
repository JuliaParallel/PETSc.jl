using Test
using PETSc
using SciMLBase

# ── Reference problem (PETSc TS tutorial ex51) ───────────────────────────────
function ex51_rhs!(du, u, p, t)
    du[1] = cos(t)
    du[2] = sin(u[2])
    return nothing
end

ex51_exact(t) = [sin(t), 2 * atan(exp(t) * tan(0.5))]

# ── Van der Pol oscillator (stiff for large mu) ───────────────────────────────
function vdp!(du, u, p, t)
    mu = p[1]
    du[1] = u[2]
    du[2] = mu * ((1 - u[1]^2) * u[2] - u[1])
    return nothing
end

@testset "TSRosW Rosenbrock-W" begin
    @testset "ex51 with TSRosW(\"ra34pw2\") matches analytical solution" begin
        u0 = [0.0, 1.0]
        tspan = (0.0, 1.0)
        prob = ODEProblem(ex51_rhs!, u0, tspan)
        sol = solve(prob, TSRosW("ra34pw2", ["-snes_fd"]); dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        @test sol.t[end] ≈ tspan[2]
        @test sol.u[end] ≈ ex51_exact(tspan[2]) atol = 1e-3
    end

    @testset "ex51 with TSRosW(\"rodas3\") converges similarly" begin
        u0 = [0.0, 1.0]
        tspan = (0.0, 1.0)
        prob = ODEProblem(ex51_rhs!, u0, tspan)
        sol = solve(prob, TSRosW("rodas3", ["-snes_fd"]); dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        @test sol.u[end] ≈ ex51_exact(tspan[2]) atol = 1e-2
    end

    @testset "Van der Pol stiff (mu = 1000) reaches tspan[2] with -snes_fd" begin
        u0 = [2.0, 0.0]
        tspan = (0.0, 1.0)
        prob = ODEProblem(vdp!, u0, tspan, [1000.0])
        sol = solve(prob, TSRosW("ra34pw2", ["-snes_fd"]); dt = 1e-3)
        @test sol.retcode == ReturnCode.Success
        @test sol.t[end] ≈ tspan[2]
        # final state should be bounded for the limit cycle
        @test all(isfinite, sol.u[end])
        @test abs(sol.u[end][1]) < 5
    end

    @testset "out-of-place ODEProblem is rejected with informative error" begin
        u0 = [1.0]
        tspan = (0.0, 1.0)
        prob = ODEProblem((u, p, t) -> -u, u0, tspan)
        @test_throws ArgumentError solve(prob, TSRosW("ra34pw2", ["-snes_fd"]); dt = 0.1)
    end

    @testset "init / step! / solve! parity with solve" begin
        u0 = [0.0, 1.0]
        tspan = (0.0, 1.0)
        prob = ODEProblem(ex51_rhs!, u0, tspan)
        sol_oneshot = solve(prob, TSRosW("ra34pw2", ["-snes_fd"]); dt = 0.1)

        integ = init(prob, TSRosW("ra34pw2", ["-snes_fd"]); dt = 0.1)
        sol_steps = solve!(integ)
        @test sol_steps.retcode == ReturnCode.Success
        @test sol_steps.u[end] ≈ sol_oneshot.u[end] atol = 1e-12
    end
end
