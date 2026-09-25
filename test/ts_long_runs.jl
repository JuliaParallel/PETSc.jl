# test/ts_long_runs.jl
# What a long TS run needs: the equation type, restarting the step count,
# pre- and post-step hooks that may change the run, and the limits on failed
# nonlinear solves and rejected steps.

using Test
using PETSc
using MPI

MPI.Initialized() || MPI.Init()

@testset "TS for long runs" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    comm = MPI.COMM_SELF
    PetscScalar = petsclib.PetscScalar
    LibPETSc = PETSc.LibPETSc

    # du/dt = 0 with forward Euler and a fixed step of 0.1
    function still_ts()
        ts = PETSc.TS(petsclib, comm)
        PETSc.set_type!(ts, :euler)
        PETSc.set_adapt_type!(ts, :none)
        PETSc.set_rhs_function!(ts) do F, _ts, _t, _u
            fill!(F, 0)
        end
        PETSc.set_time!(ts, 0.0)
        PETSc.set_timestep!(ts, 0.1)
        PETSc.set_max_time!(ts, 100.0)
        return ts
    end

    @testset "equation type" begin
        ts = PETSc.TS(petsclib, comm)
        @test PETSc.equation_type(ts) == LibPETSc.TS_EQ_UNSPECIFIED
        @test PETSc.set_equation_type!(ts, LibPETSc.TS_EQ_DAE_IMPLICIT_INDEX1) === ts
        @test PETSc.equation_type(ts) == LibPETSc.TS_EQ_DAE_IMPLICIT_INDEX1
        PETSc.destroy!(ts)
    end

    @testset "set_step_number! continues the count" begin
        ts = still_ts()
        u = PETSc.PetscVec(petsclib, PetscScalar[1])
        @test PETSc.set_step_number!(ts, 10) === ts
        PETSc.set_max_steps!(ts, 13)
        steps = Int[]
        PETSc.set_monitor!((_ts, step, _t, _u) -> push!(steps, step), ts)
        PETSc.solve!(u, ts)
        @test steps == 10:13
        @test PETSc.step_number(ts) == 13
        PETSc.destroy!(ts)
        PETSc.destroy!(u)
    end

    @testset "pre- and post-step hooks" begin
        ts = still_ts()
        u = PETSc.PetscVec(petsclib, PetscScalar[1])
        PETSc.set_max_steps!(ts, 3)
        starts = Float64[]
        @test PETSc.set_pre_step!(t -> push!(starts, PETSc.current_time(t)), ts) === ts
        # the post-step hook changes the state: add 1 to the solution after every step
        @test PETSc.set_post_step!(ts) do t
            x = PETSc.solution(t)
            x[1] = x[1] + 1
            PETSc.assemble!(x)
        end === ts
        PETSc.solve!(u, ts)
        @test starts ≈ [0.0, 0.1, 0.2]
        @test u[1] == 4

        # step! runs neither hook
        PETSc.step!(ts)
        @test length(starts) == 3
        @test u[1] == 4
        PETSc.destroy!(ts)
        PETSc.destroy!(u)
    end

    @testset "hooks see user_ctx; a throwing hook is rethrown" begin
        ts = still_ts()
        u = PETSc.PetscVec(petsclib, PetscScalar[1])
        PETSc.set_max_steps!(ts, 2)
        count = Ref(0)
        PETSc.set_user_ctx!(ts, count)
        PETSc.set_post_step!((_t, ctx) -> ctx[] += 1, ts)
        PETSc.solve!(u, ts)
        @test count[] == 2

        PETSc.set_max_steps!(ts, 4)
        PETSc.set_post_step!(_t -> throw(DomainError(0.0, "deliberate")), ts)
        @test_throws DomainError redirect_stderr(devnull) do
            PETSc.solve!(u, ts)
        end
        PETSc.destroy!(ts)
        PETSc.destroy!(u)
    end

    @testset "limits on failed solves and rejected steps" begin
        # backward Euler on a residual `residual!(F, x)` with Jacobian `J = dF/dx`
        function implicit_ts(residual!, jacobian; options...)
            ts = PETSc.TS(petsclib, comm; options...)
            PETSc.set_type!(ts, :beuler)
            PETSc.set_ifunction!(ts) do F, t, _t, x, _xdot
                residual!(F, x, t)
                PETSc.assemble!(F)
            end
            J = LibPETSc.MatCreateSeqDense(petsclib, comm, petsclib.PetscInt(1), petsclib.PetscInt(1), zeros(PetscScalar, 1))
            PETSc.set_ijacobian!(ts, J) do A, _P, _t, _tt, x, _xdot, _shift
                A[1, 1] = jacobian(x[1])
                PETSc.assemble!(A)
            end
            PETSc.set_time!(ts, 0.0)
            PETSc.set_timestep!(ts, 0.1)
            PETSc.set_max_time!(ts, 1.0)
            return ts, J
        end
        u = PETSc.PetscVec(petsclib, PetscScalar[1])

        # exp(x) = 0 has no root, so every nonlinear solve fails
        no_root!(F, x, _t) = (F[1] = exp(x[1]))
        ts, J = implicit_ts(no_root!, exp; snes_max_it = 3)
        @test PETSc.set_error_if_step_fails!(ts, false) === ts
        @test PETSc.set_max_snes_failures!(ts, 3) === ts
        PETSc.solve!(u, ts)
        @test PETSc.converged_reason(ts) == LibPETSc.TS_DIVERGED_NONLINEAR_SOLVE
        @test PETSc.snes_failures(ts) == 3

        # with PETSc's default, the same failure throws
        PETSc.set_error_if_step_fails!(ts, true)
        PETSc.set_time!(ts, 0.0)
        @test_throws LibPETSc.PetscError redirect_stderr(devnull) do
            PETSc.solve!(u, ts)
        end
        PETSc.destroy!(ts)
        PETSc.destroy!(J)

        # a domain error rejects the step instead of counting as a failed solve
        function domain_error!(F, x, t)
            F[1] = x[1]
            PETSc.set_function_domain_error!(PETSc.snes(t))
        end
        ts, J = implicit_ts(domain_error!, x -> one(x))
        PETSc.set_error_if_step_fails!(ts, false)
        @test PETSc.set_max_step_rejections!(ts, 2) === ts
        PETSc.solve!(u, ts)
        @test PETSc.converged_reason(ts) == LibPETSc.TS_DIVERGED_STEP_REJECTED
        @test PETSc.snes_failures(ts) == 0
        @test PETSc.step_rejections(ts) > 2
        PETSc.destroy!(ts)
        PETSc.destroy!(J)
        PETSc.destroy!(u)
    end

    PETSc.finalize(petsclib)
end
