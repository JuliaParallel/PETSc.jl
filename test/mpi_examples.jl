using Test
using MPI

function find_sources(path::String, sources=String[])
  if isdir(path)
    for entry in readdir(path)
      find_sources(joinpath(path, entry), sources)
    end
  elseif endswith(path, ".jl")
    push!(sources, path)
  end
  sources
end

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
    success(pipeline(cmd, stderr=stderr))
  end
end
