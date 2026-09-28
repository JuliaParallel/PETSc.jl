module PETScTestUtils

using Libdl
using PETSc: LibPETSc

export find_sources, without_petsc_traceback

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

# Runs `f()` with PETSc's return-only error handler, so an error raised on purpose
# fails the call as usual without printing PETSc's traceback to stderr
function without_petsc_traceback(f, petsclib)
    handler = Libdl.dlsym(Libdl.dlopen(petsclib.petsc_library), :PetscReturnErrorHandler)
    LibPETSc.PetscPushErrorHandler(petsclib, handler, C_NULL)
    try
        return f()
    finally
        LibPETSc.PetscPopErrorHandler(petsclib)
    end
end

end # module PETScTestUtils
