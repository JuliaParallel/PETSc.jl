using Test
using PETSc

include(joinpath(dirname(@__DIR__), "examples", "ex16.jl"))

function capture_stdout_to_string(f::Function)
    path, io = mktemp()
    close(io)
    try
        open(path, "w") do out
            redirect_stdout(out) do
                f()
            end
        end
        return read(path, String)
    finally
        isfile(path) && rm(path)
    end
end

@testset "TS ex16 example" begin
    petsclib = PETSc.getlib(PetscScalar = Float64)
    PETSc.initialize(petsclib)
    try
        function solve_arkimex(method::String)
            solve_ex16(;
                petsclib,
                options = [
                    "-ts_type",
                    "arkimex",
                    "-ts_arkimex_type",
                    method,
                    "-ts_adapt_type",
                    "none",
                ],
                verbose = false,
            )
        end

        result_default_1 = solve_ex16(;
            petsclib,
            options = String[],
            verbose = false,
        )
        result_default_2 = solve_ex16(;
            petsclib,
            options = String[],
            verbose = false,
        )
        result_no_imex = solve_ex16(;
            petsclib,
            imex = false,
            options = String[],
            verbose = false,
        )
        result_mu_100 = solve_ex16(;
            petsclib,
            mu = 100.0,
            options = String[],
            verbose = false,
        )
        result_myark2 = solve_arkimex("myark2")
        extra_arkimex_results = Dict(
            "ars122" => solve_arkimex("ars122"),
            "ars443" => solve_arkimex("ars443"),
            "3" => solve_arkimex("3"),
            "4" => solve_arkimex("4"),
            "5" => solve_arkimex("5"),
        )
        monitor_output = capture_stdout_to_string() do
            solve_ex16(;
                petsclib,
                monitor = true,
                options = String[],
                verbose = false,
            )
        end

        @test result_default_1.final_time ≈ 0.5 atol = 100 * eps(Float64)
        @test result_default_1.steps > 0
        @test length(result_default_1.solution) == 2
        @test all(isfinite, result_default_1.solution)

        @test result_default_2.final_time ≈ 0.5 atol = 100 * eps(Float64)
        @test result_default_2.steps == result_default_1.steps
        @test result_default_2.solution ≈ result_default_1.solution rtol = 1e-12

        @test result_no_imex.final_time ≈ 0.5 atol = 100 * eps(Float64)
        @test result_no_imex.steps > 0
        @test !result_no_imex.imex
        @test length(result_no_imex.solution) == 2
        @test all(isfinite, result_no_imex.solution)

        @test result_mu_100.final_time ≈ 0.5 atol = 100 * eps(Float64)
        @test result_mu_100.steps > 0
        @test result_mu_100.mu ≈ 100.0
        @test length(result_mu_100.solution) == 2

        @test result_myark2.final_time ≈ 0.5 atol = 100 * eps(Float64)
        @test result_myark2.steps > 0
        @test length(result_myark2.solution) == 2
        @test all(isfinite, result_myark2.solution)

        for (method, result) in extra_arkimex_results
            @testset "ARKIMEX $method" begin
                @test result.final_time ≈ 0.5 atol = 100 * eps(Float64)
                @test result.steps > 0
                @test length(result.solution) == 2
                @test all(isfinite, result.solution)
            end
        end

        @test occursin("[0.0]", monitor_output)
        @test occursin("[0.5]", monitor_output)
        @test length(split(chomp(monitor_output), "\n")) >= 6
    finally
        if PETSc.initialized(petsclib) && !PETSc.finalized(petsclib)
            PETSc.finalize(petsclib)
        end
    end
end
