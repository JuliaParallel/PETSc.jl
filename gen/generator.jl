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

petsc_include_dir = joinpath(PETSc_dir, "lib", "petsc", "double_real_Int64", "include") |> normpath
petsc_h = joinpath(petsc_include_dir, "petsc.h")
@assert isfile(petsc_h)

mpich_include_dir = joinpath(MPICH_dir, "include") |> normpath
mpi_h = joinpath(mpich_include_dir, "mpi.h")
@assert isfile(mpi_h)

# run generator for all platforms
output_file = joinpath(@__DIR__, "..", "lib", "petsc_library.jl")
wrap(output_file, petsc_include_dir, mpich_include_dir)
