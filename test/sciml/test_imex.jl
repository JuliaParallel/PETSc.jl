using Test
using PETSc
using SciMLBase
using DiffEqBase
using DataStructures

ext = Base.get_extension(PETSc, :PETScSciMLExt)
@assert ext !== nothing
const TSARKIMEX = ext.TSARKIMEX
const TSImplicit = ext.TSImplicit

# ── Linear stiff IMEX problem ────────────────────────────────────────────────
# u' = -u (implicit) + cos(t) (explicit)
# analytical: u(t) = (1/2) * (cos(t) + sin(t)) + (u0 - 1/2) * exp(-t)
function linear_implicit!(du, u, p, t)
    du[1] = -u[1]
    return nothing
end

function linear_explicit!(du, u, p, t)
    du[1] = cos(t)
    return nothing
end

linear_imex_exact(t, u0) = 0.5 * (cos(t) + sin(t)) + (u0 - 0.5) * exp(-t)

# ── Van der Pol (mu = 100) split ─────────────────────────────────────────────
# f1 (implicit, stiff): mu * ((1 - u1^2) * u2 - u1) on the second component
# f2 (explicit, non-stiff): u2 on the first component
function vdp_implicit!(du, u, p, t)
    mu = p[1]
    du[1] = 0.0
    du[2] = mu * ((1 - u[1]^2) * u[2] - u[1])
    return nothing
end

function vdp_explicit!(du, u, p, t)
    du[1] = u[2]
    du[2] = 0.0
    return nothing
end

@testset "Step 5 — TSARKIMEX with SplitODEProblem" begin
    @testset "Linear IMEX with TSARKIMEX(\"2e\") matches analytical" begin
        u0 = [1.0]
        tspan = (0.0, 1.0)
        prob = SplitODEProblem(linear_implicit!, linear_explicit!, u0, tspan)
        sol = solve(prob, TSARKIMEX("2e", ["-snes_fd"]); dt = 0.05)
        @test sol.retcode == ReturnCode.Success
        @test sol.t[end] ≈ tspan[2]
        @test sol.u[end][1] ≈ linear_imex_exact(tspan[2], u0[1]) atol = 1e-3
    end

    @testset "Linear IMEX with TSARKIMEX(\"3\") matches analytical" begin
        u0 = [1.0]
        tspan = (0.0, 1.0)
        prob = SplitODEProblem(linear_implicit!, linear_explicit!, u0, tspan)
        sol = solve(prob, TSARKIMEX("3", ["-snes_fd"]); dt = 0.05)
        @test sol.retcode == ReturnCode.Success
        @test sol.u[end][1] ≈ linear_imex_exact(tspan[2], u0[1]) atol = 1e-3
    end

    @testset "Van der Pol mu=100 split with TSARKIMEX reaches tspan[2]" begin
        u0 = [2.0, 0.0]
        tspan = (0.0, 1.0)
        prob = SplitODEProblem(vdp_implicit!, vdp_explicit!, u0, tspan, [100.0])
        sol = solve(prob, TSARKIMEX("3", ["-snes_fd"]); dt = 1e-2)
        @test sol.retcode == ReturnCode.Success
        @test sol.t[end] ≈ tspan[2]
        @test all(isfinite, sol.u[end])
    end

    @testset "Degenerate split (f2 ≡ 0) reduces to implicit-only" begin
        u0 = [0.0, 1.0]
        tspan = (0.0, 1.0)
        f1!(du, u, p, t) = (du[1] = cos(t); du[2] = sin(u[2]); nothing)
        f2!(du, u, p, t) = (du .= 0.0; nothing)
        prob_split = SplitODEProblem(f1!, f2!, u0, tspan)
        sol_split = solve(prob_split, TSARKIMEX("2e", ["-snes_fd"]); dt = 0.05)
        @test sol_split.retcode == ReturnCode.Success
        # The split with zero explicit part should match a non-split solve of
        # the same total RHS treated as implicit.
        prob_full = ODEProblem(f1!, u0, tspan)
        sol_full = solve(prob_full, TSARKIMEX("2e", ["-snes_fd"]); dt = 0.05)
        @test sol_full.retcode == ReturnCode.Success
        @test sol_split.u[end] ≈ sol_full.u[end] atol = 1e-8
    end

    @testset "init / step! / solve! parity for TSARKIMEX" begin
        u0 = [1.0]
        tspan = (0.0, 1.0)
        prob = SplitODEProblem(linear_implicit!, linear_explicit!, u0, tspan)
        sol_oneshot = solve(prob, TSARKIMEX("2e", ["-snes_fd"]); dt = 0.05)
        integ = init(prob, TSARKIMEX("2e", ["-snes_fd"]); dt = 0.05)
        sol_steps = solve!(integ)
        @test sol_steps.retcode == ReturnCode.Success
        @test sol_steps.u[end] ≈ sol_oneshot.u[end] atol = 1e-12
    end
end
