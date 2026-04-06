using Test
using PETSc

include(joinpath(dirname(@__DIR__), "examples", "ex51_implicit.jl"))

@testset "TS ex51 implicit example" begin
    petsclib = PETSc.getlib(PetscScalar = Float64)
    PETSc.initialize(petsclib)
    try
        result_stage_1 = solve_ex51_implicit(;
            petsclib,
            options = ["-ts_type", "irk", "-ts_irk_type", "gauss", "-ts_irk_nstages", "1"],
            save_trajectory = false,
            verbose = false,
        )
        result_stage_2 = solve_ex51_implicit(;
            petsclib,
            options = ["-ts_type", "irk", "-ts_irk_type", "gauss", "-ts_irk_nstages", "2"],
            save_trajectory = false,
            verbose = false,
        )
        result_stage_3 = solve_ex51_implicit(;
            petsclib,
            options = ["-ts_type", "irk", "-ts_irk_type", "gauss", "-ts_irk_nstages", "3"],
            save_trajectory = false,
            verbose = false,
        )
        result_stage_4 = solve_ex51_implicit(;
            petsclib,
            options = ["-ts_type", "irk", "-ts_irk_type", "gauss", "-ts_irk_nstages", "4"],
            save_trajectory = false,
            verbose = false,
        )
        result_analytic = solve_ex51_implicit(;
            petsclib,
            options = String[],
            jacobian_mode = :analytic,
            save_trajectory = false,
            verbose = false,
        )

        @test result_stage_1.final_time ≈ 1.0 atol = 100 * eps(Float64)
        @test length(result_stage_1.solution) == 2
        @test result_stage_1.error < 5.0e-3

        @test result_stage_2.final_time ≈ 1.0 atol = 100 * eps(Float64)
        @test length(result_stage_2.solution) == 2
        @test result_stage_2.error < 5.0e-6

        @test result_stage_3.final_time ≈ 1.0 atol = 100 * eps(Float64)
        @test length(result_stage_3.solution) == 2
        @test result_stage_3.error < 5.0e-9

        @test result_stage_4.final_time ≈ 1.0 atol = 100 * eps(Float64)
        @test length(result_stage_4.solution) == 2
        @test result_stage_4.error < 5.0e-10

        @test result_analytic.final_time ≈ 1.0 atol = 100 * eps(Float64)
        @test length(result_analytic.solution) == 2
        @test result_analytic.error < 5.0e-5

        err = try
            solve_ex51_implicit(;
                petsclib,
                options = ["-snes_mf"],
                save_trajectory = false,
                verbose = false,
            )
            nothing
        catch ex
            ex
        end
        @test err isa ArgumentError
        @test occursin("TSIRK/Gauss", sprint(showerror, err))
    finally
        if PETSc.initialized(petsclib) && !PETSc.finalized(petsclib)
            PETSc.finalize(petsclib)
        end
    end
end
