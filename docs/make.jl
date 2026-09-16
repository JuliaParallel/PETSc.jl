using Documenter, PETSc

# The maintainer guide for regenerating the LibPETSc wrappers lives next to the generator;
# mirror it into the manual so it is published with the docs (docs/src/man/wrapping.md is ignored by git).
cp(joinpath(@__DIR__, "..", "wrapping", "WRAPPING.md"), joinpath(@__DIR__, "src", "man", "wrapping.md"); force = true)

makedocs(;
    modules=[PETSc],
    sitename="PETSc.jl",
    checkdocs=:exports,  # Only check exported functions, skip LibPETSc internals
    warnonly=true,  # Warn but don't error for any documentation issues
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        size_threshold_warn = nothing,  # Disable size warnings for large low-level API pages
        size_threshold = nothing,  # Disable size errors for large low-level API pages
    ),
    pages=[
        "Home" => "index.md",
        "Installation" => "man/installation.md",
        "Getting Started" => "man/getting_started.md",
        "High-level interface" => Any[
            "Vec" =>  "man/vec.md",
            "Mat" =>  "man/mat.md",
            "DM" =>  "man/dm.md",
            "DMDA" =>  "man/dmda.md",
            "DMStag" =>  "man/dmstag.md",
            "DMPlex" =>  "man/dmplex.md",
            "KSP" =>  "man/ksp.md",
            "SNES" =>  "man/snes.md",
            "TS" =>  "man/ts.md",
        ],
        "Low-level interface (LibPETSc)" => Any[
            "Introduction" =>  "man/lowlevel_intro.md",
            "Vec" =>  "man/vec_lowlevel.md",
            "Mat" =>  "man/mat_lowlevel.md",
            "DM" => Any[
                "DM" =>  "man/dm_lowlevel.md",
                "DMDA" =>  "man/dmda_lowlevel.md",
                "DMPlex" =>  "man/dmplex_lowlevel.md",
                "DMStag" =>  "man/dmstag_lowlevel.md",
                "DMSwarm" =>  "man/dmswarm_lowlevel.md",
                "DMForest" =>  "man/dmforest_lowlevel.md",
                "DMNetwork" =>  "man/dmnetwork_lowlevel.md",
                "DMShell and others" =>  "man/dmshell_lowlevel.md",
            ],
            "KSP" =>  "man/ksp_lowlevel.md",
            "PC (Preconditioners)" =>  "man/pc_lowlevel.md",
            "SNES" =>  "man/snes_lowlevel.md",
            "TS (Time Stepping)" =>  "man/ts_lowlevel.md",
            "Tao (Optimization)" =>  "man/tao_lowlevel.md",
            "IS (Index Sets)" =>  "man/is_lowlevel.md",
            "PetscViewer (I/O)" =>  "man/petscviewer_lowlevel.md",
            "PetscSection (DOF Layout)" =>  "man/petscsection_lowlevel.md",
            "PetscSF (Communication)" =>  "man/petscsf_lowlevel.md",
            "AO (Application Ordering)" =>  "man/ao_lowlevel.md",
            "Discretization (PetscFE, PetscDS, ...)" =>  "man/discretization_lowlevel.md",
            "PetscOptions" =>  "man/petscoptions_lowlevel.md",
            "PetscObject, Logging, Devices" =>  "man/petscobject_lowlevel.md",
            "PetscDraw (Graphics)" =>  "man/petscdraw_lowlevel.md",
            "Sys (Runtime utilities)" =>  "man/sys_lowlevel.md",
            "Utilities (Random, Layout, ...)" =>  "man/utilities_lowlevel.md",
            "PF, Partitioner, Regressor, ..." =>  "man/otherclasses_lowlevel.md",
        ],
        "Utilities" => "man/utilities.md",
        "Running on HPC Systems" => "man/hpc.md",
        "GPU Support (CUDA)" => "man/gpu.md",
        "FAQ"  => "man/FAQ.md",
        "Naming Conventions" => "man/naming.md",
        "Contributing"  => "man/contributing.md",
        "Regenerating the wrappers" => "man/wrapping.md",
        "Funding" => "man/funding.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaParallel/PETSc.jl.git",
    branch = "gh-pages",
    target = "build",
    devbranch = "main",
    devurl = "dev",
    forcepush=true,
    push_preview = true
)
