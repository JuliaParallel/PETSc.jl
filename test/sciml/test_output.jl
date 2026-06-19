using Test
using PETSc
using SciMLBase

# `ext` is needed below for the internal `ext._build_opts` helper; the
# algorithm types themselves come from `using PETSc` (they are exported).
ext = Base.get_extension(PETSc, :PETScSciMLExt)
@assert ext !== nothing

# u' = -u with analytical solution exp(-t) starting from u0 = 1.
function decay!(du, u, p, t)
    du[1] = -u[1]
    return nothing
end

@testset "Step 6 — save_everystep / saveat / save_end" begin
    u0 = [1.0]
    tspan = (0.0, 1.0)
    prob = ODEProblem(decay!, u0, tspan)

    @testset "default (no saveat): save_everystep = true, all steps saved" begin
        sol = solve(prob, TSRK("3bs"); dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        @test length(sol.t) > 2   # every step is recorded by default
        @test issorted(sol.t)
        @test sol.t[1] ≈ 0.0
        @test sol.t[end] ≈ 1.0
        @test sol.u[end][1] ≈ exp(-1) atol = 1e-3
        # spot-check intermediate values against the analytical solution
        for k in 1:length(sol.t)
            @test sol.u[k][1] ≈ exp(-sol.t[k]) atol = 5e-3
        end
    end

    @testset "save_everystep = false saves only start and end" begin
        sol = solve(prob, TSRK("3bs"); dt = 0.1, save_everystep = false)
        @test sol.retcode == ReturnCode.Success
        @test length(sol.t) == 2
        @test sol.t[1] ≈ 0.0
        @test sol.t[end] ≈ 1.0
        @test sol.u[end][1] ≈ exp(-1) atol = 1e-3
    end

    @testset "default with saveat: save_everystep = false (only requested times)" begin
        saveat_times = [0.25, 0.5, 0.75]
        sol = solve(prob, TSRK("3bs"); dt = 0.1, saveat = saveat_times)
        @test sol.retcode == ReturnCode.Success
        # With saveat provided, save_everystep defaults to false:
        # only start (0.0) + saveat times + end (1.0) are saved, no intermediates.
        @test length(sol.t) == 1 + length(saveat_times) + 1
        @test sol.t[1] ≈ 0.0
        @test sol.t[end] ≈ 1.0
        for ts in saveat_times
            @test any(t -> isapprox(t, ts; atol = 1e-12), sol.t)
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
            save_start = true, save_end = true, save_discretes = true,
            callback = nothing,
            adaptive = true, reltol = 1e-3, abstol = 1e-6, maxiters = 1000,
        )
        ordered = Float64[]
        while !isempty(opts.saveat)
            push!(ordered, pop!(opts.saveat) / -1.0)
        end
        @test ordered == [0.75, 0.5, 0.25]
    end

    @testset "scalar saveat is a spacing, not a single time (SciML semantics)" begin
        # tspan = (0, 1), saveat = 0.25 should save at 0.25, 0.5, 0.75, 1.0
        sol = solve(prob, TSRK("3bs"); dt = 0.1, saveat = 0.25)
        @test sol.retcode == ReturnCode.Success
        for ts in (0.25, 0.5, 0.75, 1.0)
            @test any(t -> isapprox(t, ts; atol = 1e-12), sol.t)
        end
    end

    @testset "save_everystep + saveat keeps trajectory sorted" begin
        # Forces a saveat time inside a step (saveat = 0.25 with dt = 0.4
        # means PETSc's first step ends at 0.4; without the drain-first fix
        # the trajectory would record 0.4 before 0.25).
        sol = solve(
            prob, TSRK("3bs");
            dt = 0.4, save_everystep = true, saveat = [0.25],
        )
        @test sol.retcode == ReturnCode.Success
        @test issorted(sol.t)
        # 0.25 (interpolated) appears between save_start (0.0) and the next
        # step endpoint.
        @test any(t -> isapprox(t, 0.25; atol = 1e-12), sol.t)
        # and the trajectory has no duplicate timestamps.
        @test length(unique(sol.t)) == length(sol.t)
    end

    @testset "integer tspan is promoted to a float type internally" begin
        # The integrator must not truncate PETSc's Float64 step times back
        # into Int when prob.tspan has an integer eltype.
        prob_i = ODEProblem(decay!, [1.0], (0, 1))
        sol = solve(prob_i, TSRK("3bs"); dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        @test sol.t[end] ≈ 1.0
        @test sol.u[end][1] ≈ exp(-1) atol = 1e-3
    end
end
