using Test
using PETSc

include(joinpath(dirname(@__DIR__), "examples", "ex51.jl"))

@testset "TS ex51 example" begin
    petsclib = PETSc.getlib(PetscScalar = Float64)
    try
        result1 = solve_ex51(; petsclib, options = String[], verbose = false)
        result2 = solve_ex51(; petsclib, options = String[], verbose = false)

        @test result1.final_time ≈ 1.0 atol = 100 * eps(Float64)
        @test length(result1.solution) == 2
        @test result1.error < 0.5

        @test result2.final_time ≈ 1.0 atol = 100 * eps(Float64)
        @test length(result2.solution) == 2
        @test result2.error < 0.5
    finally
        if PETSc.initialized(petsclib) && !PETSc.finalized(petsclib)
            PETSc.finalize(petsclib)
        end
    end
end
