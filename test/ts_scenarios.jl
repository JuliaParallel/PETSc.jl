using Test
using PETSc
using MPI

# Two scenarios for the high-level TS interface.
#
# The first solves a problem the package already solves through the low-level
# bindings, and checks the two agree to the last digit. The second is shaped
# the way the low-level layer makes awkward.
#
# `examples/ex16.jl` is loaded into its own module because `test/ts_ex16.jl`
# also includes it, and its `@cfunction` pointers are `const`: a second
# include into the same namespace would be a redefinition error, and the order
# the two test files run in should not matter.
module Ex16Reference
include(joinpath(dirname(dirname(@__FILE__)), "examples", "ex16.jl"))
end

@testset "High-level TS scenarios" begin
    petsclib = PETSc.getlib(PetscScalar = Float64)
    PETSc.initialize(petsclib)
    PetscScalar = petsclib.PetscScalar
    PetscInt = petsclib.PetscInt
    comm = Sys.iswindows() ? LibPETSc.PETSC_COMM_SELF : MPI.COMM_SELF

    @testset "van der Pol: same answer as the low-level ex16" begin
        # `examples/ex16.jl` is PETSc's TS tutorial ex16 written against the
        # low-level bindings: four `@cfunction` blocks, a mutable context
        # struct, `pointer_from_objref`, and a `GC.@preserve` region wrapped
        # around every registration. 
        # Here we is the same problem, the same method and the same step size 
        # through the high-level interface, and the assertion is that the two agree exactly.
        mu = 1000.0
        imex = true
        final_time = 0.5
        dt = 0.01

        ts = PETSc.TS(
            petsclib,
            comm;
            exact_final_time = PETSc.LibPETSc.TS_EXACTFINALTIME_STEPOVER,
        )
        PETSc.set_type!(ts, :beuler)
        PETSc.set_problem_type!(ts, PETSc.LibPETSc.TS_NONLINEAR)
        PETSc.set_user_ctx!(ts, (; mu, imex))

        u = PETSc.VecSeq(petsclib, 2)
        PETSc.withlocalarray!(u; read = false, write = true) do ua
            ua[1] = 2.0
            ua[2] = -2.0 / 3.0 + 10.0 / (81.0 * mu) - 292.0 / (2187.0 * mu * mu)
        end
        PETSc.set_solution!(ts, u)

        jac = PETSc.MatSeqAIJ(petsclib, 2, 2, PetscInt(2))
        jac[1, [1, 2]] = PetscScalar.([1.0, 1.0])
        jac[2, [1, 2]] = PetscScalar.([1.0, 1.0])
        PETSc.assemble!(jac)
        PETSc.LibPETSc.MatZeroEntries(petsclib, jac)

        # G(u, t), the explicit half of the IMEX split.
        PETSc.set_rhs_function!(ts) do F, _ts, _t, x, ctx
            PETSc.withlocalarray!(
                (x, F);
                read = (true, false),
                write = (false, true),
            ) do xa, Fa
                Fa[1] = ctx.imex ? xa[2] : 0.0
                Fa[2] = 0.0
            end
            return 0
        end

        # F(u_t, u, t), the implicit half.
        PETSc.set_ifunction!(ts) do F, _ts, _t, x, xdot, ctx
            PETSc.withlocalarray!(
                (x, xdot, F);
                read = (true, true, false),
                write = (false, false, true),
            ) do xa, xda, Fa
                Fa[1] = xda[1] + (ctx.imex ? 0.0 : xa[2])
                Fa[2] =
                    xda[2] -
                    ctx.mu *
                    ((1.0 - xa[1] * xa[1]) * xa[2] - xa[1])
            end
            return 0
        end

        PETSc.set_ijacobian!(ts, jac) do A, _P, _ts, _t, x, _xdot, shift, ctx
            x1, x2 = PETSc.withlocalarray!(x; read = true, write = false) do xa
                (xa[1], xa[2])
            end
            A[1, 1] = PetscScalar(shift)
            A[1, 2] = PetscScalar(ctx.imex ? 0.0 : 1.0)
            A[2, 1] = PetscScalar(ctx.mu * (2.0 * x1 * x2 + 1.0))
            A[2, 2] = PetscScalar(shift - ctx.mu * (1.0 - x1 * x1))
            PETSc.assemble!(A)
            return 0
        end

        PETSc.set_timestep!(ts, dt)
        PETSc.set_max_time!(ts, final_time)

        PETSc.solve!(u, ts)

        solution = copy(u[:])
        steps = PETSc.step_number(ts)
        solve_time = PETSc.solve_time(ts)

        reference = Ex16Reference.solve_ex16(;
            petsclib,
            mu,
            imex,
            final_time,
            dt,
            options = String[],
            verbose = false,
        )

        # Same arithmetic in the same order, so this is exact, not merely close.
        @test steps == reference.steps
        @test solve_time == reference.final_time
        @test solution == reference.solution

        # And the answer is the one ex16 is known to give.
        @test solve_time ≈ 0.5 atol = 100 * eps(Float64)
        @test length(solution) == 2
        @test all(isfinite, solution)

        @test PETSc.converged_reason(ts) ==
              PETSc.LibPETSc.TS_CONVERGED_TIME
        @test PETSc.snes_iterations(ts) >= steps

        PETSc.destroy!(ts)
        PETSc.destroy(u)
        PETSc.destroy(jac)
    end

    @testset "damping sweep, closures and a recorded trajectory" begin
        # du/dt = [-k 1; -1 -k] u from u(0) = [1, 0], whose solution is
        # exp(-k t) * [cos t, -sin t]: a spiral that decays at rate k.
        #
        # Three things here are what the low-level layer makes awkward. The
        # callbacks close over `k` directly, so the sweep needs no context
        # struct, no `pointer_from_objref` and no `GC.@preserve`. The monitor
        # appends to an ordinary Julia vector, which stays rooted because `ts`
        # holds the closure. And the whole stepper is built and thrown away
        # once per parameter, which with hand-written `@cfunction` pointers
        # means either a context struct mutated between runs or a fresh set of
        # trampolines each time.
        dt = 0.005
        nsteps = 200
        tfinal = dt * nsteps

        for k in (0.0, 0.5, 2.0)
            ts = PETSc.TS(petsclib, comm)
            # Crank-Nicolson: second order, and exactly norm-preserving on the
            # undamped problem, which a wrong shift in the Jacobian would break.
            PETSc.set_type!(ts, :cn)
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

            # F(t, u, u_t) = u_t - A(k) u
            PETSc.set_ifunction!(ts) do F, _ts, _t, x, xdot
                PETSc.withlocalarray!(
                    (x, xdot, F);
                    read = (true, true, false),
                    write = (false, false, true),
                ) do xa, xda, Fa
                    Fa[1] = xda[1] + k * xa[1] - xa[2]
                    Fa[2] = xda[2] + xa[1] + k * xa[2]
                end
                return 0
            end

            # dF/du + shift dF/du_t = shift*I - A(k)
            PETSc.set_ijacobian!(ts, J) do A, _P, _ts, _t, _x, _xdot, shift
                A[1, 1] = PetscScalar(shift + k)
                A[1, 2] = PetscScalar(-1)
                A[2, 1] = PetscScalar(1)
                A[2, 2] = PetscScalar(shift + k)
                PETSc.assemble!(A)
                return 0
            end

            trajectory = NTuple{3, Float64}[]
            PETSc.set_monitor!(ts) do _ts, step, t, x
                push!(trajectory, (Float64(t), Float64(x[1]), Float64(x[2])))
                return 0
            end

            PETSc.set_time!(ts, 0.0)
            PETSc.set_timestep!(ts, dt)
            PETSc.set_max_time!(ts, tfinal)
            PETSc.solve!(u, ts)

            decay = exp(-k * tfinal)
            @test u[1] ≈ decay * cos(tfinal) rtol = 1e-3
            @test u[2] ≈ -decay * sin(tfinal) rtol = 1e-3

            # One monitor call before the first step, one after each of them.
            @test length(trajectory) == nsteps + 1
            @test trajectory[1] == (0.0, 1.0, 0.0)
            @test trajectory[end][1] ≈ tfinal
            @test trajectory[end][2] ≈ u[1]
            @test trajectory[end][3] ≈ u[2]

            radii = [hypot(p[2], p[3]) for p in trajectory]
            if k == 0
                # Crank-Nicolson on a skew-symmetric right-hand side is a
                # Cayley transform, so the radius is conserved to round-off.
                @test all(r -> isapprox(r, 1.0; rtol = 1e-10), radii)
            else
                @test all(<(0), diff(radii))
                @test radii[end] ≈ decay rtol = 1e-3
            end

            PETSc.destroy!(ts)
            PETSc.destroy(u)
            PETSc.destroy(J)
        end
    end

    PETSc.finalize(petsclib)
end
