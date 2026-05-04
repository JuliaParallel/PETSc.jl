using Test
using PETSc
using SciMLBase
using DiffEqBase
using LinearAlgebra

ext = Base.get_extension(PETSc, :PETScSciMLExt)
@assert ext !== nothing

include(joinpath(dirname(dirname(@__DIR__)), "examples", "ex51_sciml.jl"))

@testset "ex51_sciml example" begin
    @testset "default algorithm (TSRK 5dp), fixed step, reaches final time and is accurate" begin
        result = solve_ex51(; verbose = false)
        @test result.final_time ≈ 1.0 atol = 100 * eps(Float64)
        @test length(result.solution) == 2
        @test result.error < 5.0e-8
    end

    @testset "TSRK 3bs is less accurate than 5dp (fixed step)" begin
        result_3bs = solve_ex51(; alg = PETSc.TSRK("3bs"), verbose = false)
        result_5dp = solve_ex51(; alg = PETSc.TSRK("5dp"), verbose = false)
        @test result_3bs.final_time ≈ 1.0 atol = 100 * eps(Float64)
        @test result_3bs.error < 5.0e-4
        @test result_5dp.error < result_3bs.error
    end

    @testset "repeated solves give the same result (fixed step)" begin
        result_1 = solve_ex51(; verbose = false)
        result_2 = solve_ex51(; verbose = false)
        @test result_1.error ≈ result_2.error rtol = 1e-12
    end

    @testset "smaller dt improves accuracy (fixed step)" begin
        result_coarse = solve_ex51(; dt = 0.25, verbose = false)
        result_fine   = solve_ex51(; dt = 0.05, verbose = false)
        @test result_fine.error < result_coarse.error
    end

    @testset "adaptive = true is accepted and reaches final time" begin
        result = solve_ex51(; adaptive = true, verbose = false)
        @test result.final_time ≈ 1.0 atol = 100 * eps(Float64)
        @test length(result.solution) == 2
        @test result.error < 5.0e-4
    end
end
