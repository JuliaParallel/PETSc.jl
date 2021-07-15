#! /bin/bash julia --project generator.jl
using Pkg
using Pkg.Artifacts
using Clang.Generators
using Clang.Generators.JLLEnvs
using PETSc_jll
using MPICH_jll

cd(@__DIR__)

# headers
PETSc_toml = joinpath(dirname(pathof(PETSc_jll)), "..", "Artifacts.toml")
PETSc_dir = Pkg.Artifacts.ensure_artifact_installed("PETSc", PETSc_toml)

MPICH_toml = joinpath(dirname(pathof(MPICH_jll)), "..", "Artifacts.toml")
MPICH_dir = Pkg.Artifacts.ensure_artifact_installed("MPICH", MPICH_toml)

petsc_include_dir = joinpath(PETSc_dir, "include") |> normpath
petsc_h = joinpath(petsc_include_dir, "petsc.h")
@assert isfile(petsc_h)

mpich_include_dir = joinpath(MPICH_dir, "include") |> normpath
mpi_h = joinpath(mpich_include_dir, "mpi.h")
@assert isfile(mpi_h)

# load common option
options = load_options(joinpath(@__DIR__, "generator.toml"))

# run generator for all platforms
for target in JLLEnvs.JLL_ENV_TRIPLES[1:1]
    @info "processing $target"

    options["general"]["output_file_path"] =
        joinpath(@__DIR__, "..", "lib", "$target.jl")

    args = get_default_args(target)
    push!(args, "-I$petsc_include_dir")
    push!(args, "-isystem$mpich_include_dir")

    header_files = [petsc_h]

    ctx = create_context(header_files, args, options)

    build!(ctx)
end
