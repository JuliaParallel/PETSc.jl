# test/blas_threads.jl
# PETSc's BLAS runs in Julia's OpenBLAS pool, so `-blas_num_threads` has to size
# that pool: PETSc's own setter cannot reach it.

using Test
using PETSc
using MPI
using LinearAlgebra: BLAS

MPI.Initialized() || MPI.Init()

@testset "BLAS threads" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.isinitialized(petsclib) && PETSc.finalize(petsclib)
    n0 = BLAS.get_num_threads()

    @testset "-blas_num_threads sizes Julia's pool" begin
        options = ["-blas_num_threads", "2"]
        PETSc.initialize(petsclib; options)
        @test BLAS.get_num_threads() == 2
        @test options == ["-blas_num_threads", "2"]     # the caller's vector is untouched
        PETSc.finalize(petsclib)
        @test BLAS.get_num_threads() == 2               # finalize leaves the pool alone
    end

    @testset "without the option the pool keeps its size" begin
        BLAS.set_num_threads(3)
        PETSc.initialize(petsclib)
        @test BLAS.get_num_threads() == 3
        PETSc.finalize(petsclib)
    end

    @testset "a count below one is rejected" begin
        @test_throws ArgumentError PETSc.initialize(petsclib; options = ["-blas_num_threads", "0"])
        PETSc.finalize(petsclib)
    end

    BLAS.set_num_threads(n0)
end
