using Test
using MPI

using .PETScTestUtils: find_sources

@testset "mpi examples" begin
  examples_dir = joinpath(@__DIR__, "..", "examples")
  examples = find_sources(examples_dir)
  filter!(file -> readline(file) == "# INCLUDE IN MPI TEST", examples)

  @testset "$(basename(example))" for example in examples
    @info "MPI example $example"
    code = """
    $(Base.load_path_setup_code())
    include($(repr(example)))
    """
    cmd = `$(mpiexec()) -n 4 $(Base.julia_cmd()) --startup-file=no -e $code`
    run(cmd)
    @debug "Testing $example" Text(code) cmd
    @test success(pipeline(cmd, stderr=stderr))
  end
 
end
