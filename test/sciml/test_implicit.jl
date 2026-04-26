using Test
using PETSc
using SciMLBase
using DiffEqBase
using DataStructures

ext = Base.get_extension(PETSc, :PETScSciMLExt)
@assert ext !== nothing
const TSImplicit = ext.TSImplicit

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

@testset "Step 4 — TSImplicit fully implicit methods" begin
    @testset "TSImplicit(\"beuler\") on ex51" begin
        u0 = [0.0, 1.0]
        tspan = (0.0, 1.0)
        prob = ODEProblem(ex51_rhs!, u0, tspan)
        sol = solve(prob, TSImplicit("beuler", ["-snes_fd"]); dt = 0.01)
        @test sol.retcode == ReturnCode.Success
        @test sol.t[end] ≈ tspan[2]
        @test sol.u[end] ≈ ex51_exact(tspan[2]) atol = 5e-2
    end

    @testset "TSImplicit(\"cn\") on ex51 (2nd-order)" begin
        u0 = [0.0, 1.0]
        tspan = (0.0, 1.0)
        prob = ODEProblem(ex51_rhs!, u0, tspan)
        sol = solve(prob, TSImplicit("cn", ["-snes_fd"]); dt = 0.05)
        @test sol.retcode == ReturnCode.Success
        @test sol.u[end] ≈ ex51_exact(tspan[2]) atol = 1e-2
    end

    @testset "TSImplicit(\"theta\", 0.5) reaches second-order accuracy" begin
        u0 = [0.0, 1.0]
        tspan = (0.0, 1.0)
        prob = ODEProblem(ex51_rhs!, u0, tspan)
        sol = solve(prob, TSImplicit("theta", 0.5, ["-snes_fd"]); dt = 0.05)
        @test sol.retcode == ReturnCode.Success
        @test sol.u[end] ≈ ex51_exact(tspan[2]) atol = 1e-2
    end

    @testset "TSImplicit(\"theta\", 1.0) matches TSImplicit(\"beuler\")" begin
        u0 = [0.0, 1.0]
        tspan = (0.0, 1.0)
        prob = ODEProblem(ex51_rhs!, u0, tspan)
        sol_theta = solve(prob, TSImplicit("theta", 1.0, ["-snes_fd"]); dt = 0.01)
        sol_beuler = solve(prob, TSImplicit("beuler", ["-snes_fd"]); dt = 0.01)
        @test sol_theta.retcode == ReturnCode.Success
        @test sol_theta.u[end] ≈ sol_beuler.u[end] atol = 1e-8
    end

    @testset "TSImplicit(\"bdf\") on stiff Van der Pol (mu = 1000)" begin
        u0 = [2.0, 0.0]
        tspan = (0.0, 1.0)
        prob = ODEProblem(vdp!, u0, tspan, [1000.0])
        sol = solve(prob, TSImplicit("bdf", ["-snes_fd"]); dt = 1e-3)
        @test sol.retcode == ReturnCode.Success
        @test sol.t[end] ≈ tspan[2]
        @test all(isfinite, sol.u[end])
    end

    @testset "Unknown TSImplicit subtype raises ArgumentError" begin
        u0 = [0.0, 1.0]
        tspan = (0.0, 1.0)
        prob = ODEProblem(ex51_rhs!, u0, tspan)
        @test_throws ArgumentError solve(
            prob, TSImplicit("does-not-exist", ["-snes_fd"]); dt = 0.1,
        )
    end

    @testset "init / step! / solve! parity for TSImplicit" begin
        u0 = [0.0, 1.0]
        tspan = (0.0, 1.0)
        prob = ODEProblem(ex51_rhs!, u0, tspan)
        sol_oneshot = solve(prob, TSImplicit("cn", ["-snes_fd"]); dt = 0.05)
        integ = init(prob, TSImplicit("cn", ["-snes_fd"]); dt = 0.05)
        sol_steps = solve!(integ)
        @test sol_steps.retcode == ReturnCode.Success
        @test sol_steps.u[end] ≈ sol_oneshot.u[end] atol = 1e-12
    end
end
