# This is modeled on a combination of ideas from
#  - https://github.com/JuliaGPU/KernelAbstractions.jl/blob/545457ad3918b4d5aa687d6423f3e33909f368b6/test/examples.jl
#  - https://github.com/lcw/Bennu.jl/blob/54ddea41ad04d024fbd24687124ab43461609585/test/runtests.jl

using Test
using Pkg

function find_sources(path::String, sources = String[])
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

    # XXX: Not sure why this test fails on windows...
    if Sys.iswindows()
        push!(stale_examples, joinpath(examples_dir, "laplacian.jl"))
    end

    filter!(
        file ->
            file ∉ stale_examples && readline(file) != "# EXCLUDE FROM TESTING",
        examples,
    )

    julia = Base.julia_cmd()
    base_dir = joinpath(@__DIR__, "..")

    old_pwd = pwd()
    try
        mktempdir() do tmp_dir
            cd(tmp_dir)
            examples_project = Pkg.Types.projectfile_path(examples_dir)
            tmp_project = Pkg.Types.projectfile_path(tmp_dir)
            cp(examples_project, tmp_project)

            inst_cmd =
            `$julia --project=$tmp_project -e "import Pkg; Pkg.develop(path=raw\"$base_dir\"); Pkg.instantiate();"`
            @test success(pipeline(inst_cmd, stderr=stderr, stdout=stdout))

            cmd(script) =
            `$julia --project=$tmp_project -e "include(raw\"$script\")"`

            @testset "$(basename(example))" for example in examples
                @test success(pipeline(cmd(example), stderr = stderr))
            end
        end
    finally
        cd(old_pwd)
    end
end
