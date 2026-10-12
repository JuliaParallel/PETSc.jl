# test/solvers/convergence.jl
# What solve! does when a solver fails: ConvergenceError with `check = true`,
# nothing with `check = false`, and without `check` a warning once per call site.

using Test
using PETSc
using MPI
using Logging
using SparseArrays: spdiagm

MPI.Initialized() || MPI.Init()
isdefined(Main, :PETScTestUtils) || include(joinpath(@__DIR__, "..", "testutils.jl"))
using .PETScTestUtils: without_petsc_traceback

# The log records `f` emits, with maxlog applied as it is for a user
function collect_logs(f)
    logger = Test.TestLogger(min_level = Logging.Warn)
    with_logger(f, logger)
    return logger.logs
end

@testset "ConvergenceError and solve!(…; check)" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    comm = MPI.COMM_SELF
    PetscScalar = petsclib.PetscScalar
    PetscInt = petsclib.PetscInt
    LibPETSc = PETSc.LibPETSc

    n = 10
    laplacian = spdiagm(-1 => -ones(PetscScalar, n - 1), 0 => fill(PetscScalar(2), n), 1 => -ones(PetscScalar, n - 1))

    @testset "KSP" begin
        A = PETSc.PetscMat(petsclib, laplacian)
        b = PETSc.PetscVec(petsclib, ones(PetscScalar, n))
        x = PETSc.PetscVec(petsclib, zeros(PetscScalar, n))

        # one CG iteration cannot reach rtol = 1e-12
        failing = PETSc.KSP(A; ksp_type = "cg", pc_type = "none", ksp_max_it = 1, ksp_rtol = 1e-12)
        err = try
            PETSc.solve!(x, failing, b; check = true)
            nothing
        catch e
            e
        end
        @test err isa PETSc.ConvergenceError
        @test err.reason == LibPETSc.KSP_DIVERGED_ITS
        @test err.iterations == 1
        @test err.solver === failing
        @test isnothing(err.function_norm)
        message = sprint(showerror, err)
        @test occursin("KSP did not converge: KSP_DIVERGED_ITS", message)
        @test occursin("maximum number of iterations (-ksp_max_it)", message)
        @test occursin("-ksp_converged_reason", message)

        @test isempty(collect_logs(() -> PETSc.solve!(x, failing, b; check = false)))
        @test PETSc.converged_reason(failing) == LibPETSc.KSP_DIVERGED_ITS

        # without `check`: one warning per call site, pointing at the caller
        solve_twice() = for _ in 1:2
            PETSc.solve!(x, failing, b)
        end
        logs = collect_logs(solve_twice)
        @test length(logs) == 1
        @test logs[1].level == Logging.Warn
        @test occursin("From v0.6 this throws", logs[1].message)
        @test logs[1].file == @__FILE__
        # a second call site warns again
        logs = collect_logs(() -> (PETSc.solve!(x, failing, b); PETSc.solve!(x, failing, b)))
        @test length(logs) == 1
        @test length(collect_logs(() -> failing \ b)) == 1

        # a converged solve is silent whatever `check` says
        working = PETSc.KSP(A; ksp_type = "cg", pc_type = "none", ksp_rtol = 1e-10)
        for check in (nothing, true, false)
            @test isempty(collect_logs(() -> PETSc.solve!(x, working, b; check)))
        end
        @test x[:] ≈ laplacian \ ones(PetscScalar, n)
        foreach(PETSc.destroy!, (working, failing, x, b, A))
    end

    @testset "SNES" begin
        # x₁² + x₁x₂ = 3, x₁x₂ + x₂² = 6 from (2, 3): one Newton step is not enough
        snes = PETSc.SNES(petsclib, comm; pc_type = "lu", snes_max_it = 1, snes_rtol = 1e-12)
        r = PETSc.PetscVec(petsclib, zeros(PetscScalar, 2))
        PETSc.set_function!(snes, r) do fx, _snes, x
            fx[1] = x[1]^2 + x[1] * x[2] - 3
            fx[2] = x[1] * x[2] + x[2]^2 - 6
        end
        J = LibPETSc.MatCreateSeqDense(petsclib, comm, PetscInt(2), PetscInt(2), zeros(PetscScalar, 4))
        PETSc.set_snes_jacobian!(snes, J) do J, _snes, x
            J[1, 1] = 2x[1] + x[2]
            J[1, 2] = x[1]
            J[2, 1] = x[2]
            J[2, 2] = x[1] + 2x[2]
            PETSc.assemble!(J)
        end
        start() = PETSc.PetscVec(petsclib, PetscScalar[2, 3])

        x = start()
        err = try
            PETSc.solve!(x, snes; check = true)
            nothing
        catch e
            e
        end
        @test err isa PETSc.ConvergenceError
        @test err.reason == LibPETSc.SNES_DIVERGED_MAX_IT
        @test err.iterations == 1
        @test err.function_norm > 0
        @test occursin("SNES did not converge: SNES_DIVERGED_MAX_IT", sprint(showerror, err))

        @test isempty(collect_logs(() -> PETSc.solve!(start(), snes; check = false)))
        @test length(collect_logs(() -> PETSc.solve!(start(), snes))) == 1
        PETSc.destroy!(snes)
        PETSc.destroy!(J)
    end

    @testset "TS" begin
        # backward Euler on exp(x) = 0, which has no root, so every step fails
        function failing_ts()
            ts = PETSc.TS(petsclib, comm; snes_max_it = 3)
            PETSc.set_type!(ts, :beuler)
            PETSc.set_ifunction!(ts) do F, _ts, _t, x, _xdot
                F[1] = exp(x[1])
                PETSc.assemble!(F)
            end
            J = LibPETSc.MatCreateSeqDense(petsclib, comm, PetscInt(1), PetscInt(1), zeros(PetscScalar, 1))
            PETSc.set_ijacobian!(ts, J) do A, _P, _ts, _t, x, _xdot, _shift
                A[1, 1] = exp(x[1])
                PETSc.assemble!(A)
            end
            PETSc.set_time!(ts, 0.0)
            PETSc.set_timestep!(ts, 0.1)
            PETSc.set_max_time!(ts, 1.0)
            PETSc.set_max_snes_failures!(ts, 1)
            return ts, J
        end
        u = PETSc.PetscVec(petsclib, PetscScalar[1])

        # the run stops without an error from PETSc: `check` decides
        ts, J = failing_ts()
        PETSc.set_error_if_step_fails!(ts, false)
        err = try
            PETSc.solve!(u, ts; check = true)
            nothing
        catch e
            e
        end
        @test err isa PETSc.ConvergenceError
        @test err.reason == LibPETSc.TS_DIVERGED_NONLINEAR_SOLVE
        @test occursin("TS did not converge", sprint(showerror, err))
        PETSc.set_time!(ts, 0.0)
        @test isempty(collect_logs(() -> PETSc.solve!(u, ts; check = false)))
        PETSc.set_time!(ts, 0.0)
        @test length(collect_logs(() -> PETSc.solve!(u, ts))) == 1
        PETSc.set_time!(ts, 0.0)
        @test_throws PETSc.ConvergenceError PETSc.step!(ts; check = true)
        PETSc.destroy!(ts)
        PETSc.destroy!(J)

        # PETSc's default raises the failure itself: a PetscError, or a
        # ConvergenceError with `check = true`
        ts, J = failing_ts()
        @test_throws LibPETSc.PetscError without_petsc_traceback(() -> PETSc.solve!(u, ts), petsclib)
        PETSc.set_time!(ts, 0.0)
        @test_throws PETSc.ConvergenceError without_petsc_traceback(
            () -> PETSc.solve!(u, ts; check = true), petsclib,
        )
        PETSc.destroy!(ts)
        PETSc.destroy!(J)
    end

    PETSc.finalize(petsclib)
end
