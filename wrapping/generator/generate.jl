# Regenerate src/autowrapped from a PETSc source tree.
#
#   julia --project=wrapping/generator wrapping/generator/generate.jl --petsc-dir PATH [--api wrapping/api/petsc-X.Y.Z.json] [--out DIR]
#
# Without --api the snapshot is (re)created with getapi_dump.py from PETSC_DIR.
include(joinpath(@__DIR__, "src", "PetscWrapGen.jl"))
using .PetscWrapGen

function main(args)
    petsc_dir = api = ""
    out = joinpath(dirname(dirname(@__DIR__)), "src", "autowrapped")
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--petsc-dir"; petsc_dir = args[i+1]; i += 2
        elseif a == "--api"; api = args[i+1]; i += 2
        elseif a == "--out"; out = args[i+1]; i += 2
        else error("unknown argument $a") end
    end
    isempty(petsc_dir) && error("--petsc-dir is required")
    if isempty(api)
        v = match(r"PETSC_VERSION_(MAJOR|MINOR|SUBMINOR)\s+(\d+)", read(joinpath(petsc_dir, "include", "petscversion.h"), String))
        api = joinpath(@__DIR__, "api", "petsc-snapshot.json")
        run(`python3 $(joinpath(@__DIR__, "getapi_dump.py")) $petsc_dir $api`)
    end
    generate(; api_json = api, petsc_dir = petsc_dir, outdir = out, wrapping_dir = @__DIR__)
end

main(ARGS)
