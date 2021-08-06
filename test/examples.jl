using Test

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

@testset "examples" begin
  examples_dir = joinpath(@__DIR__, "..", "examples")
  stale_examples_dir = joinpath(@__DIR__, "..", "examples", "stale")

  examples = find_sources(examples_dir)
  stale_examples = find_sources(stale_examples_dir)

  filter!(file -> file ∉ stale_examples && readline(file) != "# EXCLUDE FROM TESTING", examples)

  @testset "$(basename(example))" for example in examples
    code = """
    $(Base.load_path_setup_code())
    include($(repr(example)))
    """
    cmd = `$(Base.julia_cmd()) --startup-file=no -e $code`
    @debug "Testing $example" Text(code) cmd
    @test success(pipeline(cmd, stderr=stderr))
  end

end
