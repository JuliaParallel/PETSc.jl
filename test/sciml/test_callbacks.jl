using Test
using PETSc
using SciMLBase

# u' = -u, exact: exp(-t)
function decay!(du, u, p, t)
    du[1] = -u[1]
    return nothing
end

@testset "Discrete callbacks and terminate!" begin
    u0 = [1.0]
    tspan = (0.0, 1.0)
    prob = ODEProblem(decay!, u0, tspan)

    @testset "DiscreteCallback + terminate! stops the integration early" begin
        condition = (u, t, integ) -> t >= 0.5
        affect! = integ -> terminate!(integ)
        cb = DiscreteCallback(condition, affect!)
        sol = solve(prob, TSRK("3bs"); dt = 0.1, callback = cb)
        @test sol.retcode == ReturnCode.Terminated
        @test sol.t[end] < tspan[2]
        @test sol.t[end] >= 0.5
    end

    @testset "DiscreteCallback that modifies state changes the final solution" begin
        # Without any callback, u(1) = exp(-1).
        sol_ref = solve(prob, TSRK("3bs"); dt = 0.1)

        # Callback fires once at t >= 0.5 and resets u to 0. After firing it
        # raises a flag in the closure so it does not fire again.
        fired = Ref(false)
        condition = (u, t, integ) -> !fired[] && t >= 0.5
        function affect!(integ)
            integ.u[1] = 0.0
            fired[] = true
            SciMLBase.u_modified!(integ, true)
        end
        cb = DiscreteCallback(condition, affect!)
        sol = solve(prob, TSRK("3bs"); dt = 0.1, callback = cb)
        @test sol.retcode == ReturnCode.Success
        @test sol.t[end] ≈ tspan[2]
        # After the reset, u continues from ~0; final value should be near 0.
        @test sol.u[end][1] < sol_ref.u[end][1] / 10
    end

    @testset "CallbackSet with two discrete callbacks" begin
        seen_a = Ref(false)
        seen_b = Ref(false)
        cb_a = DiscreteCallback(
            (u, t, integ) -> t >= 0.3,
            integ -> (seen_a[] = true; nothing),
        )
        cb_b = DiscreteCallback(
            (u, t, integ) -> t >= 0.7,
            integ -> (seen_b[] = true; nothing),
        )
        cbs = CallbackSet(cb_a, cb_b)
        sol = solve(prob, TSRK("3bs"); dt = 0.1, callback = cbs)
        @test sol.retcode == ReturnCode.Success
        @test seen_a[]
        @test seen_b[]
    end

    @testset "ContinuousCallback is rejected with ArgumentError" begin
        cc = ContinuousCallback(
            (u, t, integ) -> u[1] - 0.5,
            integ -> nothing,
        )
        err = try
            solve(prob, TSRK("3bs"); dt = 0.1, callback = cc)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("ContinuousCallback", err.msg)

        # `initialize` / `finalize` hooks must not run when the callback is
        # rejected — otherwise users get hidden side effects from an API the
        # extension claims is unsupported.
        init_ran = Ref(false)
        finalize_ran = Ref(false)
        cc_hooks = ContinuousCallback(
            (u, t, integ) -> u[1] - 0.5,
            integ -> nothing;
            initialize = (cb, u, t, integ) -> (init_ran[] = true; nothing),
            finalize = (cb, u, t, integ) -> (finalize_ran[] = true; nothing),
        )
        @test_throws ArgumentError solve(
            prob, TSRK("3bs"); dt = 0.1, callback = cc_hooks,
        )
        @test !init_ran[]
        @test !finalize_ran[]

        # A CallbackSet that mixes discrete and continuous callbacks should
        # also be rejected, since the continuous half cannot be honoured.
        cb_d = DiscreteCallback((u, t, integ) -> false, integ -> nothing)
        cbs = CallbackSet(cb_d, cc)
        @test_throws ArgumentError solve(
            prob, TSRK("3bs"); dt = 0.1, callback = cbs,
        )
    end

    @testset "tstops kwarg is rejected with ArgumentError" begin
        # Silently ignoring `tstops` could skip exact-time discrete
        # callbacks, so the wrapper now refuses it until the manual
        # `TSStep` loop can honour it.
        err = try
            solve(prob, TSRK("3bs"); dt = 0.1, tstops = [0.4, 0.6])
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("tstops", err.msg)
    end
end
