# test/solvers/solver_api.jl
# Readers and setters on KSP, SNES, TS and PC that a nested or logged solve
# needs: iteration counts, reasons, operators, tolerances, options prefixes,
# set_from_options! and comm, also on the DM, Mat, Vec and IS that have them.

using Test
using PETSc
using MPI

MPI.Initialized() || MPI.Init()

@testset "solver readers and setters" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    comm = MPI.COMM_SELF
    PetscScalar = petsclib.PetscScalar
    PetscInt = petsclib.PetscInt
    LibPETSc = PETSc.LibPETSc

    @testset "KSP: operators, prefix, options, readers" begin
        A = PETSc.PetscMat(petsclib, PetscScalar[4 -1 0; -1 4 -1; 0 -1 4])
        b = PETSc.PetscVec(petsclib, PetscScalar[3, 2, 3])
        x = PETSc.PetscVec(petsclib, PetscScalar[0, 0, 0])

        # a nested solver, built the way an inner solve is: no DM, no constructor options
        ksp = LibPETSc.KSPCreate(petsclib, comm)
        @test PETSc.set_operators!(ksp, A) === ksp
        @test PETSc.options_prefix(ksp) == ""
        @test PETSc.set_options_prefix!(ksp, "inner_") === ksp
        @test PETSc.options_prefix(ksp) == "inner_"

        # the prefixed option is read from the global database
        opts = PETSc.PetscOptions(petsclib; inner_ksp_type = "cg", inner_pc_type = "jacobi")
        push!(opts)
        @test PETSc.set_from_options!(ksp) === ksp
        pop!(opts)
        @test PETSc.type_name(ksp) === :cg
        @test PETSc.type_name(PETSc.pc(ksp)) === :jacobi

        PETSc.solve!(x, ksp, b)
        @test x[:] ≈ PetscScalar[1, 1, 1]
        @test PETSc.iteration_number(ksp) > 0
        @test Integer(PETSc.converged_reason(ksp)) > 0

        p = PETSc.pc(ksp)
        @test PETSc.set_options_prefix!(p, "inner_") === p
        @test PETSc.options_prefix(p) == "inner_"

        foreach(PETSc.destroy!, (ksp, opts, x, b, A))
    end

    @testset "KSP: set_dm_active!" begin
        da = PETSc.DMDA(petsclib, comm, (PETSc.DM_BOUNDARY_NONE,), (8,), 1; stencil_width = 1)
        ksp = PETSc.KSP(da)
        @test PETSc.set_dm_active!(ksp, false) === ksp
        @test PETSc.set_dm_active!(ksp, :rhs, true) === ksp
        @test_throws ArgumentError PETSc.set_dm_active!(ksp, :matrix, true)
        PETSc.destroy!(ksp)
        PETSc.destroy!(da)
    end

    @testset "SNES: readers and domain error" begin
        function residual!(fx, snes, x)
            PETSc.with_local_array!(fx, x; read = (false, true), write = (true, false)) do f, u
                f[1] = u[1]^2 + u[1] * u[2] - 3
                f[2] = u[1] * u[2] + u[2]^2 - 6
            end
        end
        function jacobian!(J, snes, x)
            PETSc.with_local_array!(x; write = false) do u
                J[1, 1] = 2u[1] + u[2]
                J[1, 2] = u[1]
                J[2, 1] = u[2]
                J[2, 2] = u[1] + 2u[2]
            end
            PETSc.assemble!(J)
        end

        snes = PETSc.SNES(petsclib, comm; snes_rtol = 1e-10)
        @test PETSc.set_options_prefix!(snes, "outer_") === snes
        @test PETSc.options_prefix(PETSc.ksp(snes)) == "outer_"
        r = PETSc.PetscVec(petsclib, PetscScalar[0, 0])
        J = PETSc.PetscMat(petsclib, zeros(PetscScalar, 2, 2))
        PETSc.set_function!(residual!, snes, r)
        PETSc.set_snes_jacobian!(jacobian!, snes, J)
        x = PETSc.PetscVec(petsclib, PetscScalar[2, 3])
        PETSc.solve!(x, snes)
        @test x[:] ≈ PetscScalar[1, 2] rtol = 1e-6

        @test PETSc.iteration_number(snes) > 0
        @test PETSc.ksp_iterations(snes) >= PETSc.iteration_number(snes)
        @test PETSc.function_norm(snes) < 1e-6
        @test Integer(PETSc.converged_reason(snes)) > 0
        k = PETSc.ksp(snes)
        @test k isa LibPETSc.KSP && !PETSc.owns(k)
        s = PETSc.solution(snes)
        @test s.ptr == x.ptr && !PETSc.owns(s)
        foreach(PETSc.destroy!, (snes, J, r))

        # a residual that always reports a domain error: solve! returns, the reason says why
        snes = PETSc.SNES(petsclib, comm)
        r = PETSc.PetscVec(petsclib, PetscScalar[0, 0])
        PETSc.set_function!(snes, r) do fx, snes, u
            residual!(fx, snes, u)
            @test PETSc.set_function_domain_error!(snes) === snes
        end
        J = PETSc.PetscMat(petsclib, zeros(PetscScalar, 2, 2))
        PETSc.set_snes_jacobian!(jacobian!, snes, J)
        x2 = PETSc.PetscVec(petsclib, PetscScalar[2, 3])
        PETSc.solve!(x2, snes; check = false)
        @test PETSc.converged_reason(snes) == LibPETSc.SNES_DIVERGED_FUNCTION_DOMAIN
        foreach(PETSc.destroy!, (snes, J, r, x2, x))
    end

    @testset "TS: prev_time, prefix, constructor options" begin
        ts = PETSc.TS(petsclib, comm; ts_type = "rk")
        @test PETSc.type_name(ts) === nothing
        @test PETSc.set_from_options!(ts) === ts
        @test PETSc.type_name(ts) === :rk      # the constructor's options were applied
        @test PETSc.set_options_prefix!(ts, "run_") === ts
        @test PETSc.options_prefix(ts) == "run_"

        PETSc.set_adapt_type!(ts, :none)
        u = PETSc.PetscVec(petsclib, PetscScalar[1])
        PETSc.set_rhs_function!(ts) do F, _ts, _t, x
            F[1] = -x[1]
            PETSc.assemble!(F)
        end
        PETSc.set_time!(ts, 0.0)
        PETSc.set_timestep!(ts, 0.1)
        PETSc.set_max_steps!(ts, 3)
        PETSc.solve!(u, ts)
        @test PETSc.current_time(ts) ≈ 0.3
        @test PETSc.prev_time(ts) ≈ 0.2

        # inside the last step, the out-of-place form matches interpolate!
        mid = PETSc.interpolate(ts, 0.25)
        @test PETSc.owns(mid) && mid.ptr != u.ptr
        @test mid[1] ≈ exp(-0.25) rtol = 1e-2       # the RK error at Δt = 0.1 is about 1e-3
        last = PETSc.interpolate(ts, 0.3)
        @test last[1] ≈ u[1]                         # the end of the step is the solution
        PETSc.destroy!(last)
        @test PETSc.interpolate!(similar(u), ts, 0.25)[1] == mid[1]
        foreach(PETSc.destroy!, (mid, ts, u))
    end

    @testset "KSP and SNES tolerances" begin
        A = PETSc.PetscMat(petsclib, PetscScalar[4 -1 0; -1 4 -1; 0 -1 4])
        b = PETSc.PetscVec(petsclib, PetscScalar[3, 2, 3])
        x = PETSc.PetscVec(petsclib, PetscScalar[0, 0, 0])

        ksp = PETSc.KSP(A; ksp_rtol = 1e-7, ksp_type = "gmres", pc_type = "none")
        t = PETSc.tolerances(ksp)
        @test keys(t) == (:rtol, :atol, :divtol, :max_it)
        @test t.rtol == 1e-7 && t.max_it isa Int
        @test PETSc.set_tolerances!(ksp; atol = 1e-12, max_it = 7) === ksp
        @test PETSc.tolerances(ksp) == (; rtol = 1e-7, atol = 1e-12, divtol = t.divtol, max_it = 7)

        # the limit reaches the solve: one iteration, then it stops
        PETSc.set_tolerances!(ksp; rtol = 1e-30, atol = 0, max_it = 1)
        PETSc.solve!(x, ksp, b; check = false)
        @test PETSc.iteration_number(ksp) == 1
        @test PETSc.converged_reason(ksp) == LibPETSc.KSP_DIVERGED_ITS

        # u² = 2, from u = 1
        snes = PETSc.SNES(petsclib, comm)
        r = PETSc.PetscVec(petsclib, PetscScalar[0])
        J = PETSc.PetscMat(petsclib, zeros(PetscScalar, 1, 1))
        PETSc.set_function!(snes, r) do f, _, u
            f[1] = u[1]^2 - 2
            PETSc.assemble!(f)
        end
        PETSc.set_snes_jacobian!(snes, J) do J, _, u
            J[1, 1] = 2u[1]
            PETSc.assemble!(J)
        end
        t = PETSc.tolerances(snes)
        @test keys(t) == (:atol, :rtol, :stol, :max_it, :max_funcs)
        @test PETSc.set_tolerances!(snes; rtol = 1e-9, max_funcs = 50) === snes
        @test PETSc.tolerances(snes) == (; t..., rtol = 1e-9, max_funcs = 50)

        PETSc.set_tolerances!(snes; atol = 0, rtol = 0, stol = 0, max_it = 2)
        u = PETSc.PetscVec(petsclib, PetscScalar[1])
        PETSc.solve!(u, snes; check = false)
        @test PETSc.iteration_number(snes) == 2
        @test PETSc.converged_reason(snes) == LibPETSc.SNES_DIVERGED_MAX_IT

        foreach(PETSc.destroy!, (snes, u, J, r, ksp, x, b, A))
    end

    @testset "set_from_options!, options prefix and comm on PC, DM, Mat, Vec, IS" begin
        opts = PETSc.PetscOptions(petsclib; mine_pc_type = "jacobi", mine_mat_type = "dense",
                                  mine_vec_type = "mpi")

        p = LibPETSc.PCCreate(petsclib, comm)
        M = LibPETSc.MatCreate(petsclib, comm)
        LibPETSc.MatSetSizes(petsclib, M, PetscInt(2), PetscInt(2), PetscInt(2), PetscInt(2))
        v = LibPETSc.VecCreate(petsclib, comm)
        LibPETSc.VecSetSizes(petsclib, v, PetscInt(3), PetscInt(3))
        for obj in (p, M, v)
            @test PETSc.options_prefix(obj) == ""
            @test PETSc.set_options_prefix!(obj, "mine_") === obj
            @test PETSc.options_prefix(obj) == "mine_"
        end

        push!(opts)
        @test all(obj -> PETSc.set_from_options!(obj) === obj, (p, M, v))
        pop!(opts)
        @test PETSc.type_name(p) === :jacobi
        @test LibPETSc.MatGetType(petsclib, M) == "seqdense"
        @test LibPETSc.VecGetType(petsclib, v) == "mpi"

        da = PETSc.DMDA(petsclib, comm, (PETSc.DM_BOUNDARY_NONE,), (8,), 1; stencil_width = 1)
        @test PETSc.set_options_prefix!(da, "grid_") === da
        @test PETSc.options_prefix(da) == "grid_"

        is = LibPETSc.ISCreateStride(petsclib, comm, PetscInt(3), PetscInt(0), PetscInt(1))
        for obj in (p, is)
            @test MPI.Comm_size(PETSc.comm(obj)) == MPI.Comm_size(comm)
        end

        LibPETSc.ISDestroy(petsclib, is)
        foreach(PETSc.destroy!, (da, v, M, p, opts))
    end

    PETSc.finalize(petsclib)
end
