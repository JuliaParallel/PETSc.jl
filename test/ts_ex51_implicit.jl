using Test
using PETSc

include(joinpath(dirname(@__DIR__), "examples", "ex51_implicit.jl"))

@testset "TS ex51 implicit example" begin
    petsclib = PETSc.getlib(PetscScalar = Float64)
    PETSc.initialize(petsclib)
    try
        result_fd_1 = solve_ex51_implicit(;
            petsclib,
            options = String[],
            save_trajectory = false,
            verbose = false,
        )
        result_fd_2 = solve_ex51_implicit(;
            petsclib,
            options = String[],
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

        @test result_fd_1.final_time ≈ 1.0 atol = 100 * eps(Float64)
        @test length(result_fd_1.solution) == 2
        @test result_fd_1.error < 5.0e-5

        @test result_fd_2.final_time ≈ 1.0 atol = 100 * eps(Float64)
        @test result_fd_2.error ≈ result_fd_1.error rtol = 1e-12

        @test result_analytic.final_time ≈ 1.0 atol = 100 * eps(Float64)
        @test length(result_analytic.solution) == 2
        @test result_analytic.error < 5.0e-5
        @test result_analytic.error ≈ result_fd_1.error rtol = 1e-6

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
