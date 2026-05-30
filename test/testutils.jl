module PETScTestUtils

export find_sources

function find_sources(path::String, sources = String[])
    if isdir(path)
        for entry in readdir(path)
            find_sources(joinpath(path, entry), sources)
        end
    elseif endswith(path, ".jl") && !startswith(basename(path), "_")
        push!(sources, path)
    end
    return sources
end

end # module PETScTestUtils
