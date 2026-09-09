using Test
using PETSc
using SciMLBase

# ── Reference problem (PETSc TS tutorial ex51), Float32 variant ──────────────
#   u1' = cos(t),   u2' = sin(u2)
function ex51_rhs_f32!(du, u, p, t)
    du[1] = cos(t)
    du[2] = sin(u[2])
    return nothing
end

ex51_exact_f32(t) = Float32[sin(t), 2 * atan(exp(t) * tan(0.5))]

@testset "Float32 in-place ODE problems" begin
    @testset "auto-selected petsclib is the Float32 library" begin
        u0 = Float32[0.0, 1.0]
        tspan = (0.0f0, 1.0f0)
        prob = ODEProblem(ex51_rhs_f32!, u0, tspan)
        integ = init(prob, TSRK("3bs"); dt = 0.1f0)
        @test integ.petsclib.PetscScalar === Float32
        @test integ.petsclib.PetscReal === Float32
        PETSc.destroy(integ)
    end

    @testset "explicit RK on ex51 keeps Float32 state and matches solution" begin
        u0 = Float32[0.0, 1.0]
        tspan = (0.0f0, 1.0f0)
        prob = ODEProblem(ex51_rhs_f32!, u0, tspan)
        sol = solve(prob, TSRK("5dp"); dt = 0.1f0)
        @test sol.retcode == ReturnCode.Success
        @test eltype(sol.u[end]) === Float32
        @test eltype(sol.t) === Float32
        @test sol.t[end] ≈ tspan[2]
        @test sol.u[end] ≈ ex51_exact_f32(1.0f0) atol = 1.0f-3
    end

    @testset "the time argument reaches the user RHS as Float32" begin
        seen_T = Ref{Any}(nothing)
        f!(du, u, p, t) = (seen_T[] = typeof(t); du[1] = -u[1]; nothing)
        prob = ODEProblem(f!, Float32[1.0], (0.0f0, 1.0f0))
        sol = solve(prob, TSRK("3bs"); dt = 0.05f0)
        @test sol.retcode == ReturnCode.Success
        @test seen_T[] === Float32
    end

    @testset "scalar exponential decay u' = -u" begin
        u0 = Float32[1.0]
        prob = ODEProblem((du, u, p, t) -> (du[1] = -u[1]), u0, (0.0f0, 1.0f0))
        sol = solve(prob, TSRK("3bs"); dt = 0.05f0)
        @test sol.retcode == ReturnCode.Success
        @test eltype(sol.u[end]) === Float32
        @test sol.u[end][1] ≈ exp(-1.0f0) atol = 1.0f-3
    end

    @testset "implicit (IFunction) path with -snes_fd" begin
        u0 = Float32[0.0, 1.0]
        prob = ODEProblem(ex51_rhs_f32!, u0, (0.0f0, 1.0f0))
        sol = solve(prob, TSImplicit("beuler", ["-snes_fd"]); dt = 0.01f0)
        @test sol.retcode == ReturnCode.Success
        @test eltype(sol.u[end]) === Float32
        @test sol.u[end] ≈ ex51_exact_f32(1.0f0) atol = 5.0f-2
    end

    @testset "IMEX split problem stays Float32" begin
        f1!(du, u, p, t) = (du[1] = -u[1]; nothing)   # stiff / implicit
        f2!(du, u, p, t) = (du[1] = cos(t); nothing)  # non-stiff / explicit
        prob = SplitODEProblem(f1!, f2!, Float32[1.0], (0.0f0, 1.0f0))
        sol = solve(prob, TSARKIMEX("2e", ["-snes_fd"]); dt = 0.05f0)
        @test sol.retcode == ReturnCode.Success
        @test eltype(sol.u[end]) === Float32
    end
end
