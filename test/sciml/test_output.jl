using Test
using PETSc
using SciMLBase
using DiffEqBase
using DataStructures

ext = Base.get_extension(PETSc, :PETScSciMLExt)
@assert ext !== nothing
const TSRK = ext.TSRK

# u' = -u with analytical solution exp(-t) starting from u0 = 1.
function decay!(du, u, p, t)
    du[1] = -u[1]
    return nothing
end

@testset "Step 6 — save_everystep / saveat / save_end" begin
    u0 = [1.0]
    tspan = (0.0, 1.0)
    prob = ODEProblem(decay!, u0, tspan)

    @testset "default: only start/end states saved" begin
        sol = solve(prob, TSRK("3bs"); dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        @test length(sol.t) == 2
        @test sol.t[1] ≈ 0.0
        @test sol.t[end] ≈ 1.0
        @test sol.u[end][1] ≈ exp(-1) atol = 1e-3
    end

    @testset "save_everystep = true populates the trajectory" begin
        sol = solve(prob, TSRK("3bs"); dt = 0.1, save_everystep = true)
        @test sol.retcode == ReturnCode.Success
        @test length(sol.t) > 2
        @test issorted(sol.t)
        # spot-check a couple of intermediate values against the analytical solution
        for k in 1:length(sol.t)
            @test sol.u[k][1] ≈ exp(-sol.t[k]) atol = 5e-3
        end
    end

    @testset "saveat (Vector) hits exactly the requested times" begin
        saveat = [0.25, 0.5, 0.75]
        sol = solve(prob, TSRK("3bs"); dt = 0.1, saveat)
        @test sol.retcode == ReturnCode.Success
        # save_start (t0=0.0) + saveat (3 entries) + save_end (tf=1.0)
        @test sol.t[1] ≈ 0.0
        @test sol.t[end] ≈ 1.0
        # the requested saveat times should appear, sorted
        for ts in saveat
            @test any(t -> isapprox(t, ts; atol = 1e-12), sol.t)
        end
        for k in 1:length(sol.t)
            @test sol.u[k][1] ≈ exp(-sol.t[k]) atol = 5e-3
        end
    end

    @testset "saveat without save_start / save_end gives just those times" begin
        saveat = [0.25, 0.5, 0.75]
        sol = solve(
            prob, TSRK("3bs");
            dt = 0.1, saveat,
            save_start = false, save_end = false,
        )
        @test sol.retcode == ReturnCode.Success
        @test length(sol.t) == length(saveat)
        @test all(any(t -> isapprox(t, ts; atol = 1e-12), sol.t) for ts in saveat)
    end

    @testset "backward saveat heap orders times in integration direction" begin
        # PETSc TS does not support negative-dt forward integration on RK
        # (TSAdaptChoose rejects it). Verify directly that the saveat heap is
        # built with `tdir * t` so backward solves would dispense times in
        # decreasing order — the actual mechanism under `savevalues!`.
        opts = ext._build_opts(
            Float64, [0.25, 0.5, 0.75], (), -1.0, (1.0, 0.0);
            save_everystep = false, save_on = true,
            save_start = true, save_end = true, callback = nothing,
            reltol = 1e-3, abstol = 1e-6, maxiters = 1000, verbose = false,
        )
        ordered = Float64[]
        while !isempty(opts.saveat)
            push!(ordered, pop!(opts.saveat) / -1.0)
        end
        @test ordered == [0.75, 0.5, 0.25]
    end

    @testset "scalar saveat" begin
        sol = solve(prob, TSRK("3bs"); dt = 0.1, saveat = 0.5)
        @test sol.retcode == ReturnCode.Success
        @test any(t -> isapprox(t, 0.5; atol = 1e-12), sol.t)
    end
end
