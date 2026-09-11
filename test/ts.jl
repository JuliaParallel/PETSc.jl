using Test
using PETSc
using MPI
using LinearAlgebra: I, inv
using Logging: Logging, NullLogger, with_logger

if !Sys.iswindows()
    MPI.Initialized() || MPI.Init()
end

# High-level TS interface.
#
# The reference values are the exact solutions of the *discrete* scheme, not of
# the ODE, so the tolerances only have to cover round-off and a broken callback
# still shows up immediately.

@testset "TS" begin
    comm = Sys.iswindows() ? LibPETSc.PETSC_COMM_SELF : MPI.COMM_SELF

    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscReal = petsclib.PetscReal
        PetscInt = petsclib.PetscInt
        rtol = PetscReal == Float32 ? 1e-4 : 1e-9

        # du/dt = -u, u(0) = 1. Backward Euler with a fixed step gives
        # u_n = (1 + dt)^-n exactly, which is what the implicit paths below are
        # measured against.
        decay(dt, n) = PetscReal(1 / (1 + dt)^n)

        @testset "construction and introspection" begin
            ts = PETSc.TS(petsclib, comm)
            @test ts isa PETSc.LibPETSc.TS
            @test ts.ptr != C_NULL
            @test occursin("TS", sprint(show, ts))
            # No method has been chosen yet.
            @test PETSc.type(ts) === nothing

            PETSc.set_type!(ts, :beuler)
            @test PETSc.type(ts) === :beuler
            PETSc.set_type!(ts, :bdf)
            @test PETSc.type(ts) === :bdf

            @test PETSc.comm(ts) isa MPI.Comm
            @test PETSc.dm(ts).ptr != C_NULL

            # Nothing is read from the options database until `solve!`.
            ts_opt = PETSc.TS(petsclib, comm; ts_type = "rk")
            @test PETSc.type(ts_opt) === nothing
            PETSc.destroy!(ts_opt)

            PETSc.destroy!(ts)
            @test ts.ptr == C_NULL
        end

        @testset "time and step controls" begin
            ts = PETSc.TS(petsclib, comm)
            PETSc.set_type!(ts, :beuler)

            PETSc.set_time!(ts, 0.25)
            @test PETSc.current_time(ts) ≈ 0.25
            PETSc.set_timestep!(ts, 0.01)
            @test PETSc.timestep(ts) ≈ 0.01
            PETSc.set_max_time!(ts, 2.0)
            @test PETSc.max_time(ts) ≈ 2.0
            PETSc.set_max_steps!(ts, 17)
            @test PETSc.max_steps(ts) == 17
            @test PETSc.step_number(ts) == 0

            @test PETSc.converged_reason(ts) ==
                  PETSc.LibPETSc.TS_CONVERGED_ITERATING

            PETSc.destroy!(ts)
        end

        @testset "tolerances" begin
            ts = PETSc.TS(petsclib, comm)
            PETSc.set_type!(ts, :bdf)

            PETSc.set_tolerances!(ts; atol = 1e-9, rtol = 1e-7)
            tol = PETSc.tolerances(ts)
            @test tol.atol ≈ 1e-9
            @test tol.rtol ≈ 1e-7
            # Only the scalar tolerances were set, so the vector ones are null.
            @test tol.vatol.ptr == C_NULL
            @test tol.vrtol.ptr == C_NULL

            vatol = PETSc.VecSeq(petsclib, 2)
            PETSc.set_tolerances!(ts; vatol = vatol)
            @test PETSc.tolerances(ts).vatol.ptr != C_NULL

            PETSc.destroy!(ts)
            PETSc.destroy(vatol)
        end

        @testset "explicit right-hand side" begin
            dt, n = 0.01, 100

            # `do` block form: the callback comes first, so this works.
            ts = PETSc.TS(petsclib, comm)
            PETSc.set_type!(ts, :rk)
            PETSc.set_adapt_type!(ts, :none)
            u = PETSc.VecSeq(petsclib, 1)
            u[1] = PetscScalar(1)
            PETSc.assemble!(u)

            calls = Ref(0)
            PETSc.set_rhs_function!(ts) do F, _ts, _t, x
                calls[] += 1
                PETSc.withlocalarray!(
                    (x, F);
                    read = (true, false),
                    write = (false, true),
                ) do xa, Fa
                    Fa[1] = -xa[1]
                end
                return 0
            end

            PETSc.set_time!(ts, 0.0)
            PETSc.set_timestep!(ts, dt)
            PETSc.set_max_time!(ts, 1.0)
            PETSc.solve!(u, ts)

            @test calls[] > 0
            @test PETSc.converged_reason(ts) ==
                  PETSc.LibPETSc.TS_CONVERGED_TIME
            @test PETSc.step_number(ts) == n
            @test PETSc.solve_time(ts) ≈ 1.0
            # A Runge-Kutta method at this step size is well inside 1e-3 of the
            # true solution; a callback that never ran would not be.
            @test u[1] ≈ exp(-one(PetscReal)) rtol = 1e-3

            PETSc.destroy!(ts)
            PETSc.destroy(u)

            # Reversed argument order, for a callback that is already a value.
            ts2 = PETSc.TS(petsclib, comm)
            PETSc.set_type!(ts2, :rk)
            PETSc.set_adapt_type!(ts2, :none)
            u2 = PETSc.VecSeq(petsclib, 1)
            u2[1] = PetscScalar(1)
            PETSc.assemble!(u2)
            PETSc.set_rhs_function!(ts2, (F, _ts, _t, x) -> (F[1] = -x[1]; 0))
            PETSc.set_time!(ts2, 0.0)
            PETSc.set_timestep!(ts2, dt)
            PETSc.set_max_time!(ts2, 1.0)
            PETSc.solve!(u2, ts2)
            @test u2[1] ≈ exp(-one(PetscReal)) rtol = 1e-3

            PETSc.destroy!(ts2)
            PETSc.destroy(u2)
        end

        @testset "implicit residual and Jacobian" begin
            dt, n = 0.1, 10

            ts = PETSc.TS(petsclib, comm)
            PETSc.set_type!(ts, :beuler)
            PETSc.set_adapt_type!(ts, :none)

            u = PETSc.VecSeq(petsclib, 1)
            u[1] = PetscScalar(1)
            PETSc.assemble!(u)

            J = PETSc.MatSeqAIJ(petsclib, 1, 1, PetscInt(1))
            J[1, 1] = PetscScalar(1)
            PETSc.assemble!(J)

            PETSc.set_ifunction!(ts) do F, _ts, _t, x, xdot
                PETSc.withlocalarray!(
                    (x, xdot, F);
                    read = (true, true, false),
                    write = (false, false, true),
                ) do xa, xda, Fa
                    Fa[1] = xda[1] + xa[1]
                end
                return 0
            end

            # F = u_t + u, so dF/du + shift dF/du_t = 1 + shift.
            PETSc.set_ijacobian!(ts, J) do A, _P, _ts, _t, _x, _xdot, shift
                A[1, 1] = PetscScalar(shift + 1)
                PETSc.assemble!(A)
                return 0
            end

            PETSc.set_time!(ts, 0.0)
            PETSc.set_timestep!(ts, dt)
            PETSc.set_max_time!(ts, 1.0)
            PETSc.solve!(u, ts)

            @test PETSc.step_number(ts) == n
            @test u[1] ≈ decay(dt, n) rtol = rtol
            @test PETSc.snes_iterations(ts) >= n
            @test PETSc.step_rejections(ts) == 0
            @test PETSc.snes_failures(ts) == 0

            PETSc.destroy!(ts)
            PETSc.destroy(u)
            PETSc.destroy(J)
        end

        @testset "right-hand side Jacobian" begin
            dt, n = 0.1, 10

            ts = PETSc.TS(petsclib, comm)
            PETSc.set_type!(ts, :beuler)
            PETSc.set_adapt_type!(ts, :none)

            u = PETSc.VecSeq(petsclib, 1)
            u[1] = PetscScalar(1)
            PETSc.assemble!(u)

            J = PETSc.MatSeqAIJ(petsclib, 1, 1, PetscInt(1))
            J[1, 1] = PetscScalar(1)
            PETSc.assemble!(J)

            PETSc.set_rhs_function!(ts, (F, _ts, _t, x) -> (F[1] = -x[1]; 0))
            jac_calls = Ref(0)
            PETSc.set_rhs_jacobian!(ts, J) do A, _P, _ts, _t, _x
                jac_calls[] += 1
                A[1, 1] = PetscScalar(-1)
                PETSc.assemble!(A)
                return 0
            end

            PETSc.set_time!(ts, 0.0)
            PETSc.set_timestep!(ts, dt)
            PETSc.set_max_time!(ts, 1.0)
            PETSc.solve!(u, ts)

            @test jac_calls[] > 0
            @test u[1] ≈ decay(dt, n) rtol = rtol

            PETSc.destroy!(ts)
            PETSc.destroy(u)
            PETSc.destroy(J)
        end

        @testset "coupled system exercises the shift" begin
            # du/dt = A u with A = [0 1; -1 0], so the Jacobian has off-diagonal
            # entries and a wrong shift convention cannot pass unnoticed.
            dt, n = 0.05, 20
            A = [0.0 1.0; -1.0 0.0]
            step = inv(I - dt * A)
            expected = (step^n) * [1.0, 0.0]

            ts = PETSc.TS(petsclib, comm)
            PETSc.set_type!(ts, :beuler)
            PETSc.set_adapt_type!(ts, :none)

            u = PETSc.VecSeq(petsclib, 2)
            u[1] = PetscScalar(1)
            u[2] = PetscScalar(0)
            PETSc.assemble!(u)

            J = PETSc.MatSeqAIJ(petsclib, 2, 2, PetscInt(2))
            for i in 1:2, j in 1:2
                J[i, j] = PetscScalar(0)
            end
            PETSc.assemble!(J)

            PETSc.set_ifunction!(ts) do F, _ts, _t, x, xdot
                PETSc.withlocalarray!(
                    (x, xdot, F);
                    read = (true, true, false),
                    write = (false, false, true),
                ) do xa, xda, Fa
                    Fa[1] = xda[1] - xa[2]
                    Fa[2] = xda[2] + xa[1]
                end
                return 0
            end

            PETSc.set_ijacobian!(ts, J) do M, _P, _ts, _t, _x, _xdot, shift
                M[1, 1] = PetscScalar(shift)
                M[1, 2] = PetscScalar(-1)
                M[2, 1] = PetscScalar(1)
                M[2, 2] = PetscScalar(shift)
                PETSc.assemble!(M)
                return 0
            end

            PETSc.set_time!(ts, 0.0)
            PETSc.set_timestep!(ts, dt)
            PETSc.set_max_time!(ts, n * dt)
            PETSc.solve!(u, ts)

            @test u[1] ≈ expected[1] rtol = rtol
            @test u[2] ≈ expected[2] rtol = rtol

            PETSc.destroy!(ts)
            PETSc.destroy(u)
            PETSc.destroy(J)
        end

        @testset "user context" begin
            dt, n = 0.1, 10
            ts = PETSc.TS(petsclib, comm)
            PETSc.set_type!(ts, :beuler)
            PETSc.set_adapt_type!(ts, :none)

            @test PETSc.user_ctx(ts) === nothing
            PETSc.set_user_ctx!(ts, (rate = PetscReal(2),))
            @test PETSc.user_ctx(ts).rate == 2

            u = PETSc.VecSeq(petsclib, 1)
            u[1] = PetscScalar(1)
            PETSc.assemble!(u)

            seen = Ref(false)
            # The method taking a trailing context is the one that must be
            # picked, since `user_ctx` is set.
            PETSc.set_rhs_function!(ts) do F, _ts, _t, x, ctx
                seen[] = true
                F[1] = -ctx.rate * x[1]
                return 0
            end

            PETSc.set_time!(ts, 0.0)
            PETSc.set_timestep!(ts, dt)
            PETSc.set_max_time!(ts, 1.0)
            PETSc.solve!(u, ts)

            @test seen[]
            @test u[1] ≈ PetscReal(1 / (1 + 2dt)^n) rtol = rtol

            PETSc.destroy!(ts)
            PETSc.destroy(u)
        end

        @testset "monitor" begin
            dt, n = 0.1, 10
            ts = PETSc.TS(petsclib, comm)
            PETSc.set_type!(ts, :beuler)
            PETSc.set_adapt_type!(ts, :none)

            u = PETSc.VecSeq(petsclib, 1)
            u[1] = PetscScalar(1)
            PETSc.assemble!(u)
            PETSc.set_rhs_function!(ts, (F, _ts, _t, x) -> (F[1] = -x[1]; 0))

            steps = Int[]
            times = Float64[]
            PETSc.set_monitor!(ts) do _ts, step, t, x
                push!(steps, Int(step))
                push!(times, Float64(t))
                return 0
            end

            PETSc.set_time!(ts, 0.0)
            PETSc.set_timestep!(ts, dt)
            PETSc.set_max_time!(ts, 1.0)
            PETSc.solve!(u, ts)

            # Once before the first step and once after each one.
            @test steps == collect(0:n)
            @test times[1] ≈ 0.0
            @test times[end] ≈ 1.0

            PETSc.destroy!(ts)
            PETSc.destroy(u)
        end

        @testset "solve! without a vector" begin
            dt, n = 0.1, 10
            ts = PETSc.TS(petsclib, comm)
            PETSc.set_type!(ts, :beuler)
            PETSc.set_adapt_type!(ts, :none)

            u = PETSc.VecSeq(petsclib, 1)
            u[1] = PetscScalar(1)
            PETSc.assemble!(u)
            PETSc.set_solution!(ts, u)
            @test PETSc.solution(ts).ptr == u.ptr

            PETSc.set_rhs_function!(ts, (F, _ts, _t, x) -> (F[1] = -x[1]; 0))
            PETSc.set_time!(ts, 0.0)
            PETSc.set_timestep!(ts, dt)
            PETSc.set_max_time!(ts, 1.0)

            @test PETSc.solve!(ts) === ts
            @test u[1] ≈ decay(dt, n) rtol = rtol

            PETSc.destroy!(ts)
            PETSc.destroy(u)
        end

        @testset "options are applied at solve time" begin
            dt, n = 0.1, 10
            # `ts_type` only takes effect in `solve!`, after the callbacks are
            # attached.
            ts = PETSc.TS(
                petsclib,
                comm;
                ts_type = "beuler",
                ts_adapt_type = "none",
            )
            u = PETSc.VecSeq(petsclib, 1)
            u[1] = PetscScalar(1)
            PETSc.assemble!(u)
            PETSc.set_rhs_function!(ts, (F, _ts, _t, x) -> (F[1] = -x[1]; 0))
            PETSc.set_time!(ts, 0.0)
            PETSc.set_timestep!(ts, dt)
            PETSc.set_max_time!(ts, 1.0)
            PETSc.solve!(u, ts)

            @test PETSc.type(ts) === :beuler
            @test u[1] ≈ decay(dt, n) rtol = rtol

            PETSc.destroy!(ts)
            PETSc.destroy(u)
        end

        @testset "stepping by hand" begin
            dt = 0.1
            ts = PETSc.TS(petsclib, comm)
            PETSc.set_type!(ts, :beuler)
            PETSc.set_adapt_type!(ts, :none)

            u = PETSc.VecSeq(petsclib, 1)
            u[1] = PetscScalar(1)
            PETSc.assemble!(u)
            PETSc.set_solution!(ts, u)
            PETSc.set_rhs_function!(ts, (F, _ts, _t, x) -> (F[1] = -x[1]; 0))
            PETSc.set_time!(ts, 0.0)
            PETSc.set_timestep!(ts, dt)
            PETSc.set_max_time!(ts, 1.0)

            PETSc.setup!(ts)
            PETSc.step!(ts)
            @test PETSc.step_number(ts) == 1
            @test PETSc.current_time(ts) ≈ dt
            @test u[1] ≈ decay(dt, 1) rtol = rtol

            w = PETSc.VecSeq(petsclib, 1)
            PETSc.interpolate!(w, ts, dt / 2)
            @test w.ptr != C_NULL
            @test isfinite(abs(w[1]))

            PETSc.reset!(ts)

            PETSc.destroy!(ts)
            PETSc.destroy(u)
            PETSc.destroy(w)
        end

        @testset "sub-solvers" begin
            ts = PETSc.TS(petsclib, comm)
            PETSc.set_type!(ts, :beuler)

            s = PETSc.snes(ts)
            @test s isa PETSc.LibPETSc.PetscSNES
            @test s.ptr != C_NULL
            # The handle we are given must survive being asked for twice; the
            # generated three-argument form nulls it instead.
            @test PETSc.snes(ts).ptr == s.ptr

            # A KSP is only offered for a problem declared linear.
            PETSc.set_problem_type!(ts, PETSc.LibPETSc.TS_LINEAR)
            k = PETSc.ksp(ts)
            @test k isa PETSc.LibPETSc.PetscKSP
            @test k.ptr != C_NULL

            PETSc.destroy!(ts)
        end

        @testset "destruction" begin
            ts = PETSc.TS(petsclib, comm; ts_type = "beuler")
            @test !isnothing(ts.opts)
            PETSc.destroy!(ts)
            @test ts.ptr == C_NULL
            @test isnothing(ts.opts)
            # Destroying twice is a no-op, not an error.
            PETSc.destroy!(ts)
            @test ts.ptr == C_NULL

            # `destroy` keeps working for code written against the rest of the
            # package.
            ts2 = PETSc.TS(petsclib, comm)
            PETSc.destroy(ts2)
            @test ts2.ptr == C_NULL
        end

        @testset "a callback that throws is reported, not fatal" begin
            ts = PETSc.TS(petsclib, comm)
            PETSc.set_type!(ts, :rk)
            PETSc.set_adapt_type!(ts, :none)
            u = PETSc.VecSeq(petsclib, 1)
            u[1] = PetscScalar(1)
            PETSc.assemble!(u)
            PETSc.set_rhs_function!(ts) do _F, _ts, _t, _x
                error("deliberate failure from a user callback")
            end
            PETSc.set_time!(ts, 0.0)
            PETSc.set_timestep!(ts, 0.1)
            PETSc.set_max_time!(ts, 1.0)

            # The trampoline logs the Julia error and reports failure to PETSc,
            # which unwinds and leaves `@chk` to raise. Both it and PETSc's own
            # traceback are silenced here to keep the test output readable.
            @test_throws PETSc.LibPETSc.PetscError with_logger(NullLogger()) do
                redirect_stderr(devnull) do
                    PETSc.solve!(u, ts)
                end
            end

            PETSc.destroy!(ts)
            PETSc.destroy(u)
        end

        PETSc.finalize(petsclib)
    end
end
