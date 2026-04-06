using Test
using PETSc

include(joinpath(dirname(@__DIR__), "examples", "ex51.jl"))

@testset "TS ex51 example" begin
    petsclib = PETSc.getlib(PetscScalar = Float64)
    try
        result_default_1 = solve_ex51(; petsclib, options = String[], verbose = false)
        result_3bs = solve_ex51(;
            petsclib,
            options = ["-ts_type", "rk", "-ts_rk_type", "3bs"],
            verbose = false,
        )

        result_default_2 = solve_ex51(; petsclib, options = String[], verbose = false)
        result_5dp = solve_ex51(;
            petsclib,
            options = ["-ts_type", "rk", "-ts_rk_type", "5dp"],
            verbose = false,
        )

        @test result_default_1.final_time ≈ 1.0 atol = 100 * eps(Float64)
        @test length(result_default_1.solution) == 2
        @test result_default_1.error < 5.0e-4

        @test result_3bs.final_time ≈ 1.0 atol = 100 * eps(Float64)
        @test length(result_3bs.solution) == 2
        @test result_3bs.error < 5.0e-5

        @test result_default_2.final_time ≈ 1.0 atol = 100 * eps(Float64)
        @test length(result_default_2.solution) == 2
        @test result_default_2.error ≈ result_default_1.error rtol = 1e-12

        @test result_5dp.final_time ≈ 1.0 atol = 100 * eps(Float64)
        @test length(result_5dp.solution) == 2
        @test result_5dp.error < 5.0e-8
    finally
        if PETSc.initialized(petsclib) && !PETSc.finalized(petsclib)
            PETSc.finalize(petsclib)
        end
    end
end
