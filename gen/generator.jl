#! /bin/bash julia --project generator.jl

# This is modeled on:
# https://github.com/JuliaGPU/oneAPI.jl/blob/c040cba550659f82931ba054e6c321444e53c124/res/wrap.jl
# https://github.com/JuliaLang/SuiteSparse.jl/blob/76856153eef26c008f13520ffa12288e214fe02c/gen/generator.jl

using Pkg
using Pkg.Artifacts
using PETSc_jll
using MPICH_jll

include("wrap.jl")

cd(@__DIR__)

# headers
PETSc_toml = joinpath(dirname(pathof(PETSc_jll)), "..", "Artifacts.toml")
PETSc_dir = Pkg.Artifacts.ensure_artifact_installed("PETSc", PETSc_toml)

MPICH_toml = joinpath(dirname(pathof(MPICH_jll)), "..", "Artifacts.toml")
MPICH_dir = Pkg.Artifacts.ensure_artifact_installed("MPICH", MPICH_toml)

# Copy all petsc headers to a local directory
# NOTE: Doing this by hand rather than leaving them at the original location, so I can later add docstrings to the headers
petsc_jll_include_dir = joinpath(PETSc_dir, "lib", "petsc", "double_real_Int64", "include") |> normpath
if !isdir("headers");  mkdir("headers") end
for file in readdir(petsc_jll_include_dir)
    cp(joinpath(petsc_jll_include_dir, file), joinpath("headers",file), force=true)
end

petsc_include_dir = joinpath(pwd(),"headers") |> normpath
#petsc_h = joinpath(petsc_include_dir, "petsc.h")

# for debugging, only wrap a few files, so we don't generate that many files
petscdmstag_h   = joinpath(petsc_include_dir, "petscdmstag.h")
petscsystypes_h = joinpath(petsc_include_dir, "petscsystypes.h")
#petsc_h_files   = [petscdmstag_h, petscsystypes_h]
petsc_h_files   = [petscdmstag_h]


petsc_h = [joinpath(petsc_include_dir, "petsc.h")]
@assert all(isfile.(petsc_h))

mpich_include_dir = joinpath(MPICH_dir, "include") |> normpath
mpi_h = joinpath(mpich_include_dir, "mpi.h")
@assert isfile(mpi_h)

# run generator for all platforms
output_file = joinpath(@__DIR__, "..", "lib", "petsc_library.jl")
wrap(output_file, petsc_include_dir, mpich_include_dir, petsc_h_files)
