# test/bang_returns.jl
# A `!` function returns the object it mutates (docs/src/man/naming.md §7), with
# the exceptions §7 lists: releases and package-state setters return `nothing`,
# and a do-block function returns the block's result.

using Test
using PETSc
using MPI
using LinearAlgebra: mul!, ldiv!

MPI.Initialized() || MPI.Init()

@testset "! functions return their object" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    comm = MPI.COMM_SELF
    PetscScalar = petsclib.PetscScalar

    @testset "vectors and matrices" begin
        v = PETSc.PetscVec(petsclib, PetscScalar[1, 2, 3])
        @test PETSc.set_type!(v, :seq) === v
        @test setindex!(v, PetscScalar(4), 1) === v
        @test fill!(v, PetscScalar(1)) === v
        @test PETSc.assemble!(v) === v
        @test PETSc.with_local_array!(a -> sum(a), v) == 3

        A = PETSc.PetscMat(petsclib, PetscScalar[2 0 0; 0 2 0; 0 0 2])
        @test setindex!(A, PetscScalar(3), 1, 1) === A
        @test PETSc.assemble!(A) === A
        y = PETSc.PetscVec(petsclib, PetscScalar[0, 0, 0])
        @test mul!(y, A, v) === y

        ksp = PETSc.KSP(A)
        @test PETSc.set_type!(ksp, :cg) === ksp
        p = PETSc.pc(ksp)
        @test PETSc.set_type!(p, :jacobi) === p
        x = PETSc.PetscVec(petsclib, PetscScalar[0, 0, 0])
        @test PETSc.solve!(x, ksp, y) === x
        @test ldiv!(x, ksp, y) === x

        foreach(PETSc.destroy!, (ksp, x, y, A, v))
    end

    @testset "nonlinear and time-stepping solvers" begin
        snes = PETSc.SNES(petsclib, comm)
        r = PETSc.PetscVec(petsclib, PetscScalar[0, 0])
        @test PETSc.set_type!(snes, :newtonls) === snes
        @test PETSc.set_function!((F, s, x) -> nothing, snes, r) === snes
        PETSc.destroy!(snes)
        PETSc.destroy!(r)

        ts = PETSc.TS(petsclib, comm)
        @test PETSc.set_type!(ts, :beuler) === ts
        @test PETSc.set_time!(ts, 0.0) === ts
        @test PETSc.set_timestep!(ts, 0.1) === ts
        @test PETSc.set_max_steps!(ts, 3) === ts
        @test PETSc.set_max_time!(ts, 1.0) === ts
        PETSc.destroy!(ts)
    end

    @testset "DMs and options" begin
        da = PETSc.DMDA(petsclib, comm, (PETSc.DM_BOUNDARY_NONE,), (8,), 1, 1)
        @test PETSc.set_uniform_coordinates!(da, (0.0,), (1.0,)) === da
        g = PETSc.global_vec(da)
        l = PETSc.local_vec(da)
        @test PETSc.global_to_local!(l, da, g) === l
        @test PETSc.local_to_global!(g, da, l) === g
        foreach(PETSc.destroy!, (l, g, da))

        opts = PETSc.PetscOptions(petsclib)
        @test setindex!(opts, 3, :ksp_max_it) === opts
        @test PETSc.destroy!(opts) === nothing
    end

    PETSc.finalize(petsclib)
end
