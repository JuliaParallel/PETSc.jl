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
#=
petsc_jll_include_dir = joinpath(PETSc_dir, "lib", "petsc", "double_real_Int64", "include") |> normpath
if !isdir("headers");  mkdir("headers") end
for file in readdir(petsc_jll_include_dir)
    cp(joinpath(petsc_jll_include_dir, file), joinpath("headers",file), force=true)
end
=#

petsc_include_dir = joinpath(pwd(),"headers") |> normpath
#petsc_h = joinpath(petsc_include_dir, "petsc.h")

# for debugging, only wrap a few files, so we don't generate that many files
petscdmstag_h   = joinpath(petsc_include_dir, "petscdmstag.h")
petscmacros_h   = joinpath(petsc_include_dir, "petscmacros.h")
petscconf_h   = joinpath(petsc_include_dir, "petscconf.h")
petscfix_h   = joinpath(petsc_include_dir, "petscfix.h")
petscdm_h   = joinpath(petsc_include_dir, "petscdm.h")
petscvec_h   = joinpath(petsc_include_dir, "petscvec.h")
petscvec_h   = joinpath(petsc_include_dir, "petscvec.h")


petscconf_poison_h          = joinpath(petsc_include_dir, "petscconf_poison.h")
petscversion_h              = joinpath(petsc_include_dir, "petscversion.h")
petscadvancedmacros_h       = joinpath(petsc_include_dir, "petsc/private/petscadvancedmacros.h")
petscsystypes_h = joinpath(petsc_include_dir, "petscsystypes.h")
petscmath_h                 = joinpath(petsc_include_dir, "petscmath.h")
petscerror_h                = joinpath(petsc_include_dir, "petscerror.h")
petscviewertypes_h      = joinpath(petsc_include_dir, "petscviewertypes.h")
petscoptions_h          = joinpath(petsc_include_dir, "petscoptions.h")
petsclog_h          = joinpath(petsc_include_dir, "petsclog.h")
petsctime_h         = joinpath(petsc_include_dir, "petsctime.h")
petscbt_h               = joinpath(petsc_include_dir, "petscbt.h")
petscstring_h           = joinpath(petsc_include_dir, "petscstring.h")
petsclogtypes_h         = joinpath(petsc_include_dir, "petsclogtypes.h")
petsclogdeprecated_h    = joinpath(petsc_include_dir, "petsclogdeprecated.h")
petscsectiontypes_h     = joinpath(petsc_include_dir, "petscsectiontypes.h")
petscistypes_h          = joinpath(petsc_include_dir, "petscistypes.h")
petscdrawtypes_h        = joinpath(petsc_include_dir, "petscdrawtypes.h")


petscsys_h   = joinpath(petsc_include_dir, "petscsys.h")
petscsftypes_h   = joinpath(petsc_include_dir, "petscsftypes.h")
petscis_h   = joinpath(petsc_include_dir, "petscis.h")
petscdevicetypes_h   = joinpath(petsc_include_dir, "petscdevicetypes.h")
petscviewer_h   = joinpath(petsc_include_dir, "petscviewer.h")
petscmacros_h =     joinpath(petsc_include_dir, "petscmacros.h")

#petscsystypes_h = joinpath(petsc_include_dir, "petscsystypes.h")
#petsc_h_files   = [petscdmstag_h, petscsystypes_h]
petsc_h_files   = [petscdmstag_h, 
                    #petscconf_h,
                    #petscconf_poison_h,   
                    #petscfix_h,
                    ##petscversion_h     ,         
                    ##petscadvancedmacros_h,       
                    #petscsystypes_h,
                    #petscmath_h,
                    #petscerror_h,
                    #petscviewertypes_h,
                    #petscoptions_h,
                    #petsclog_h,
                    #petsctime_h ,
                    #petscbt_h   ,
                    #petscstring_h    ,
                    #petsclogtypes_h     ,
                    #petsclogdeprecated_h  ,
                    #petscsectiontypes_h   ,
                    #petscistypes_h      ,
                    #petscdrawtypes_h,
                    #petscsys_h,
                    #petscis_h,
                    #petscsftypes_h,
                    #petscdevicetypes_h,
                    #petscviewer_h,
                    #petscmacros_h
                    ]


petsc_h = [joinpath(petsc_include_dir, "petsc.h")]
@assert all(isfile.(petsc_h))

mpich_include_dir = joinpath(MPICH_dir, "include") |> normpath
mpi_h = joinpath(mpich_include_dir, "mpi.h")
@assert isfile(mpi_h)

# run generator for all platforms
output_file = joinpath(@__DIR__, "..", "lib", "petsc_library.jl")
wrap(output_file, petsc_include_dir, mpich_include_dir, petsc_h_files)
