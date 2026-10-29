# test/solvers/options_timing.jl
# When a KSP, SNES or TS applies the options given to its constructor
# (`options_at`), and `with_options` on solvers and on a PetscOptions.

using Test
using PETSc
using MPI
using SparseArrays: sparse

MPI.Initialized() || MPI.Init()

@testset "options timing" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    comm = MPI.COMM_SELF
    PetscScalar = petsclib.PetscScalar
    LibPETSc = PETSc.LibPETSc

    M = sparse(PetscScalar[4 -1; -1 4])

    @testset "with_options on a PetscOptions" begin
        opts = PETSc.PetscOptions(petsclib; ksp_type = "cg")
        ksp = PETSc.KSP(PETSc.PetscMat(petsclib, M))
        @test PETSc.with_options(() -> PETSc.set_from_options!(ksp), opts) === ksp
        @test PETSc.type_name(ksp) === :cg

        # nested on the same options: no second push, so the stack stays intact
        @test PETSc.with_options(opts) do
            PETSc.with_options(() -> 42, opts)
        end == 42

        # popped even when the block throws
        @test_throws ErrorException PETSc.with_options(() -> error("boom"), opts)
        @test isempty(PETSc.options_stack(typeof(petsclib)))

        # a second push of the same options would break PETSc's stack
        push!(opts)
        @test_throws "already pushed" push!(opts)
        pop!(opts)
        @test isempty(PETSc.options_stack(typeof(petsclib)))
        foreach(PETSc.destroy!, (ksp, opts))
    end

    @testset "with_options on a solver" begin
        ksp = PETSc.KSP(PETSc.PetscMat(petsclib, M); ksp_type = "cg")
        @test PETSc.with_options(() -> :result, ksp) === :result
        # the constructor options are visible inside, not outside
        other = PETSc.KSP(PETSc.PetscMat(petsclib, M))
        LibPETSc.KSPSetFromOptions(petsclib, other)
        @test PETSc.type_name(other) === :gmres
        PETSc.with_options(() -> LibPETSc.KSPSetFromOptions(petsclib, other), ksp)
        @test PETSc.type_name(other) === :cg
        PETSc.destroy!(other)
        # nesting set_from_options! inside does not push twice
        @test PETSc.with_options(() -> PETSc.set_from_options!(ksp), ksp) === ksp
        PETSc.destroy!(ksp)

        bare = PETSc.KSP(PETSc.PetscMat(petsclib, M))
        @test PETSc.with_options(() -> 1, bare) == 1
        PETSc.destroy!(bare)
    end

    @testset "options_at is checked" begin
        A = PETSc.PetscMat(petsclib, M)
        @test_throws "options_at must be :construction, :each_solve or :first_solve" PETSc.KSP(
            A; options_at = :never, ksp_type = "cg")
        @test_throws ArgumentError PETSc.SNES(petsclib, comm; options_at = :sometimes)
        @test_throws ArgumentError PETSc.TS(petsclib, comm; options_at = :later)
        PETSc.destroy!(A)
    end

    # Each solver is built with `pc_type = ilu` in its options, and the PC is set
    # to Jacobi in code before each of two solves. Which one is in effect after
    # each solve shows when the options were applied. `:default` passes no
    # `options_at`, which must behave as `:first_solve`.
    expected = Dict(
        :construction => (:jacobi, :jacobi),
        :each_solve => (:ilu, :ilu),
        :first_solve => (:ilu, :jacobi),
        :default => (:ilu, :jacobi),
    )
    at_keyword(at) = at === :default ? (;) : (; options_at = at)

    function pc_types_after_two_solves(solve, pc)
        types = Symbol[]
        for _ in 1:2
            PETSc.set_type!(pc, :jacobi)
            solve()
            push!(types, PETSc.type_name(pc))
        end
        return Tuple(types)
    end

    @testset "KSP, options_at = :$at" for at in (:default, :construction, :each_solve, :first_solve)
        A = PETSc.PetscMat(petsclib, M)
        b = PETSc.PetscVec(petsclib, PetscScalar[3, 3])
        x = PETSc.PetscVec(petsclib, zeros(PetscScalar, 2))
        ksp = PETSc.KSP(A; at_keyword(at)..., pc_type = "ilu", ksp_rtol = 1e-12)
        @test pc_types_after_two_solves(() -> PETSc.solve!(x, ksp, b; check = true), PETSc.pc(ksp)) ==
            expected[at]
        @test x[:] ≈ PetscScalar[1, 1]
        foreach(PETSc.destroy!, (ksp, x, b, A))
    end

    # F(x) = M x − (3, 3), root (1, 1)
    function linear_snes(at)
        snes = PETSc.SNES(petsclib, comm; at_keyword(at)..., pc_type = "ilu", snes_rtol = 1e-12)
        r = PETSc.PetscVec(petsclib, zeros(PetscScalar, 2))
        PETSc.set_function!(snes, r) do fx, _snes, x
            fx[1] = 4x[1] - x[2] - 3
            fx[2] = -x[1] + 4x[2] - 3
        end
        J = PETSc.PetscMat(petsclib, M)
        PETSc.set_snes_jacobian!((J, _snes, _x) -> nothing, snes, J)
        return snes, r, J
    end

    @testset "SNES, options_at = :$at" for at in (:default, :construction, :each_solve, :first_solve)
        snes, r, J = linear_snes(at)
        x = PETSc.PetscVec(petsclib, zeros(PetscScalar, 2))
        solve() = (fill!(x, 0); PETSc.solve!(x, snes; check = true))
        @test pc_types_after_two_solves(solve, PETSc.pc(PETSc.ksp(snes))) == expected[at]
        @test x[:] ≈ PetscScalar[1, 1]
        foreach(PETSc.destroy!, (snes, x, r, J))
    end

    # du/dt = −M u by backward Euler, steps of 0.1. PETSc shifts the Jacobian in
    # place, so the callback fills it on every call.
    function linear_ts(at)
        ts = PETSc.TS(petsclib, comm; at_keyword(at)..., pc_type = "ilu")
        PETSc.set_type!(ts, :beuler)
        PETSc.set_adapt_type!(ts, :none)
        PETSc.set_rhs_function!(ts) do F, _ts, _t, u
            F[1] = -4u[1] + u[2]
            F[2] = u[1] - 4u[2]
        end
        J = PETSc.PetscMat(petsclib, -M)
        PETSc.set_rhs_jacobian!(ts, J) do A, _P, _ts, _t, _u
            A[1, 1], A[1, 2], A[2, 1], A[2, 2] = -4, 1, 1, -4
            PETSc.assemble!(A)
        end
        PETSc.set_time!(ts, 0.0)
        PETSc.set_timestep!(ts, 0.1)
        return ts, J
    end

    @testset "TS, options_at = :$at" for at in (:default, :construction, :each_solve, :first_solve)
        ts, J = linear_ts(at)
        u = PETSc.PetscVec(petsclib, PetscScalar[1, 1])
        t_end = Ref(0.0)
        function solve()
            t_end[] += 0.2
            PETSc.set_max_time!(ts, t_end[])
            PETSc.solve!(u, ts; check = true)
        end
        @test pc_types_after_two_solves(solve, PETSc.pc(PETSc.ksp(PETSc.snes(ts)))) == expected[at]
        @test PETSc.type_name(ts) === :beuler
        foreach(PETSc.destroy!, (ts, u, J))
    end

    @testset ":first_solve applies at the first setup!" begin
        A = PETSc.PetscMat(petsclib, M)
        ksp = PETSc.KSP(A; options_at = :first_solve, pc_type = "ilu")
        @test PETSc.setup!(ksp) === ksp
        @test PETSc.type_name(PETSc.pc(ksp)) === :ilu
        b = PETSc.PetscVec(petsclib, PetscScalar[3, 3])
        x = PETSc.PetscVec(petsclib, zeros(PetscScalar, 2))
        @test pc_types_after_two_solves(() -> PETSc.solve!(x, ksp, b; check = true), PETSc.pc(ksp)) ==
            (:jacobi, :jacobi)
        foreach(PETSc.destroy!, (ksp, x, b, A))

        snes, r, J = linear_snes(:first_solve)
        @test PETSc.setup!(snes) === snes
        @test PETSc.type_name(PETSc.pc(PETSc.ksp(snes))) === :ilu
        foreach(PETSc.destroy!, (snes, r, J))

        ts, J = linear_ts(:first_solve)
        u = PETSc.PetscVec(petsclib, PetscScalar[1, 1])
        PETSc.set_solution!(ts, u)
        p = PETSc.pc(PETSc.ksp(PETSc.snes(ts)))
        PETSc.set_type!(p, :jacobi)
        @test PETSc.setup!(ts) === ts
        @test PETSc.type_name(p) === :ilu
        foreach(PETSc.destroy!, (ts, u, J))

        # set_from_options! counts as the first application
        ts, J = linear_ts(:first_solve)
        u = PETSc.PetscVec(petsclib, PetscScalar[1, 1])
        PETSc.set_solution!(ts, u)
        PETSc.set_from_options!(ts)
        p = PETSc.pc(PETSc.ksp(PETSc.snes(ts)))
        @test PETSc.type_name(p) === :ilu
        PETSc.set_type!(p, :jacobi)
        PETSc.setup!(ts)
        @test PETSc.type_name(p) === :jacobi
        foreach(PETSc.destroy!, (ts, u, J))
    end

    @testset "step! applies the options only under :first_solve" begin
        for (at, pc_type) in ((:default, :ilu), (:first_solve, :ilu), (:each_solve, :jacobi))
            ts, J = linear_ts(at)
            u = PETSc.PetscVec(petsclib, PetscScalar[1, 1])
            PETSc.set_solution!(ts, u)
            PETSc.set_max_time!(ts, 1.0)
            p = PETSc.pc(PETSc.ksp(PETSc.snes(ts)))
            PETSc.set_type!(p, :jacobi)
            PETSc.step!(ts; check = true)
            @test PETSc.type_name(p) === pc_type
            foreach(PETSc.destroy!, (ts, u, J))
        end
    end

    @test isempty(PETSc.options_stack(typeof(petsclib)))
end
