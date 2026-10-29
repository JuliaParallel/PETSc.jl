# test/solvers/convergence.jl
# What solve! does when a solver fails: ConvergenceError by default and with
# `check = true`, nothing with `check = false`.

using Test
using PETSc
using MPI
using Logging
using SparseArrays: spdiagm
using LinearAlgebra: LinearAlgebra

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
            PETSc.solve!(x, failing, b)
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

        @test_throws PETSc.ConvergenceError PETSc.solve!(x, failing, b; check = true)
        @test_throws PETSc.ConvergenceError PETSc.solve!(failing)
        @test_throws PETSc.ConvergenceError failing \ b
        @test_throws PETSc.ConvergenceError LinearAlgebra.ldiv!(x, failing, b)
        # `nothing` meant "warn" before the default became `true`
        @test_throws TypeError PETSc.solve!(x, failing, b; check = nothing)

        # a converged solve returns quietly whatever `check` says
        working = PETSc.KSP(A; ksp_type = "cg", pc_type = "none", ksp_rtol = 1e-10)
        for check in (true, false)
            @test isempty(collect_logs(() -> PETSc.solve!(x, working, b; check)))
        end
        @test isempty(collect_logs(() -> PETSc.solve!(x, working, b)))
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
        J = PETSc.PetscMat(petsclib, zeros(PetscScalar, 2, 2))
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
            PETSc.solve!(x, snes)
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
        @test PETSc.converged_reason(snes) == LibPETSc.SNES_DIVERGED_MAX_IT
        @test_throws PETSc.ConvergenceError PETSc.solve!(start(), snes; check = true)
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
            J = PETSc.PetscMat(petsclib, zeros(PetscScalar, 1, 1))
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
            PETSc.solve!(u, ts)
            nothing
        catch e
            e
        end
        @test err isa PETSc.ConvergenceError
        @test err.reason == LibPETSc.TS_DIVERGED_NONLINEAR_SOLVE
        @test occursin("TS did not converge", sprint(showerror, err))
        PETSc.set_time!(ts, 0.0)
        @test isempty(collect_logs(() -> PETSc.solve!(u, ts; check = false)))
        @test PETSc.converged_reason(ts) == LibPETSc.TS_DIVERGED_NONLINEAR_SOLVE
        PETSc.set_time!(ts, 0.0)
        @test_throws PETSc.ConvergenceError PETSc.step!(ts)
        PETSc.set_time!(ts, 0.0)
        @test_throws PETSc.ConvergenceError PETSc.step!(ts; check = true)
        PETSc.set_time!(ts, 0.0)
        PETSc.step!(ts; check = false)
        @test PETSc.converged_reason(ts) == LibPETSc.TS_DIVERGED_NONLINEAR_SOLVE
        PETSc.destroy!(ts)
        PETSc.destroy!(J)

        # PETSc's default raises the failure itself: a ConvergenceError, or the
        # PetscError with `check = false`
        ts, J = failing_ts()
        @test_throws PETSc.ConvergenceError without_petsc_traceback(() -> PETSc.solve!(u, ts), petsclib)
        PETSc.set_time!(ts, 0.0)
        @test_throws LibPETSc.PetscError without_petsc_traceback(
            () -> PETSc.solve!(u, ts; check = false), petsclib,
        )
        PETSc.destroy!(ts)
        PETSc.destroy!(J)
    end

    PETSc.finalize(petsclib)
end
