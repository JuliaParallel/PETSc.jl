using PETSc
using Test
using MPI

using .PETScTestUtils: find_sources

@testset "examples" begin
  examples_dir = joinpath(@__DIR__, "..", "examples")
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

end
