using PETSc
using Test
using MPI

using .PETScTestUtils: find_sources

@testset "examples" begin
  examples_dir = joinpath(pkgdir(PETSc), "examples")
  examples = find_sources(examples_dir)
  filter!(file -> readline(file) != "# EXCLUDE FROM TESTING", examples)

  @testset "$(basename(example))" for example in examples
    @show example
    code = """
    $(Base.load_path_setup_code())
    include($(repr(example)))
    """
    cmd = `$(Base.julia_cmd()) --startup-file=no -e $code`
    @debug "Testing $example" Text(code) cmd
    @test success(pipeline(cmd, stderr=stderr))
  end

  # The runs above pass no options; ex19 takes PETSc options on its command line
  @testset "ex19.jl with options" begin
    example = joinpath(examples_dir, "dmda", "ex19.jl")
    code = """
    $(Base.load_path_setup_code())
    include($(repr(example)))
    """
    args = `-da_grid_x 9 -da_grid_y 9 -pc_type mg -pc_mg_levels 2 -snes_converged_reason`
    # `--` ends Julia's own options, so that it does not parse -da_grid_x
    cmd = `$(Base.julia_cmd()) --startup-file=no -e $code -- $args`
    @test success(pipeline(cmd, stderr=stderr))
  end

end
