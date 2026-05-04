# PETSc.jl

[![Build Status](https://github.com/JuliaParallel/PETSc.jl/workflows/CI/badge.svg)](https://github.com/JuliaParallel/PETSc.jl/actions/workflows/ci.yml)
[![doc stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaparallel.github.io/PETSc.jl/stable/)
[![doc dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaparallel.github.io/PETSc.jl/dev/)
[![DOI](https://zenodo.org/badge/38933145.svg)](https://doi.org/10.5281/zenodo.18274809)


PETSc.jl provides an interface to the Portable, Extensible Toolkit for Scientific Computation ([PETSc](https://petsc.org)) library, allowing the combination of Julia features (such as automatic differentiation) with the PETSc's infrastructure, including linear, nonlinear, and optimization solvers, timesteppers, domain management (DM), and more, in a distributed-memory (MPI) environment.

This package comprises three main components:

1. An automatically generated, low-level interface for large parts of the PETSc API (see `PETSc.LibPETSc`).
2. A curated, high-level, more Julian interface for selected functionality.
3. A package extension based on [SciMLBase.jl](https://github.com/SciML/SciMLBase.jl) that allows solving problems such as `ODEProblem`s with PETSc's algorithms.

The low-level interface covers nearly the entire PETSc API, but may be awkward to work with and likely requires previous experience with PETSc to use effectively.
The high-level interface is designed to be more familiar and convenient for Julia users, and allows, for example, to set matrix entries with `A[1,2] = 3.0`, rather than having to call `LibPETSc.MatSetValue`.
However, it exposes only a small portion of the functionality of the underlying library.
The SciML package extension is work in progress and currently supports only (in-place) `ODEProblem`s.

## Installation
This package can be added with the Julia command:
```julia
julia>]add PETSc
```
The installation can be tested with
```julia
julia>]test PETSc
```

## PETSc binaries

By default, the package uses a pre-built binary of PETSc (see [PETSc_jll](https://github.com/JuliaBinaryWrappers/PETSc_jll.jl)) along with a default installation of `MPI.jl`, so you don't have to install it on your machine.

If you want to use the package with a custom PETSc build, use `set_library!` to configure it once — the path is stored persistently in `LocalPreferences.toml` and no environment variables are needed:

```julia
using PETSc
PETSc.set_library!("/path/to/custom/libpetsc.so"; PetscScalar=Float64, PetscInt=Int64)
# Restart Julia — PETSc_jll is not loaded and your library is used automatically.
```

To revert to the bundled binaries: `PETSc.unset_library!()`. To check the current configuration: `PETSc.library_info()`.

To get an overview of available precompiled libraries:
```julia
julia>using PETSc
julia>[PETSc.petsclibs...]
```

## Windows users
The package currently does not work on Windows, mainly because `MicrosoftMPI_jll` does not function when used along with the precompiled version used in `PETSc_jll`.
Windows users are therefore advised to install the [Windows Subsystem for Linux](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux) (WSL) and run PETSc through there.

## Getting started
Have a look at the [documentation](https://juliaparallel.org/PETSc.jl/stable/), at the [examples](./examples/) directory or at the tests in the [test](./test) directory.
We do keep the tests up to date, so that is a good starting point.

Note, that we do not have tests in place for the whole library at this stage.
The best supported parts are `DMDA`,`DMStag`, `KSP`,`SNES`,`Vec` and `Mat` interfaces, while other parts such as `DMPlex` do not have a high-level interface or tests yet.
Users will thus have to rely on the low-level interface.

## SciML integration

`PETSc.jl` ships a package extension (`PETScSciMLExt`) that lets you solve in-place `ODEProblem`s with PETSc's TS time integrators through the standard SciML interface. The extension activates automatically when `SciMLBase` is loaded — `DiffEqBase` and `OrdinaryDiffEq` are not required.

```julia
using PETSc, SciMLBase

f!(du, u, p, t) = (du[1] = -u[1]; nothing)
prob = ODEProblem(f!, [1.0], (0.0, 1.0))

sol = solve(prob, PETSc.TSRK("3bs"); dt = 0.1)                                # explicit RK
sol = solve(prob, PETSc.TSRosW("ra34pw2", ["-snes_fd"]); dt = 1e-3)           # Rosenbrock-W
sol = solve(prob, PETSc.TSImplicit("bdf", ["-snes_fd"]); dt = 1e-3)           # BDF / theta / CN / BEuler
sol = solve(SplitODEProblem(f1!, f2!, u0, tspan),
            PETSc.TSARKIMEX("2e", ["-snes_fd"]); dt = 0.05)                   # IMEX

integrator = init(prob, PETSc.TSRK("5dp"); dt = 0.05)                         # step!/solve! interface
step!(integrator); sol = solve!(integrator)
```

`OrdinaryDiffEq` also works as the trigger since it re-exports `SciMLBase`, but the extension only depends on `SciMLBase` — no `DiffEqBase` is needed.

Per-solver PETSc command-line options are passed on the algorithm itself (a `Vector{String}` of raw tokens, e.g. `["-snes_fd", "-ts_max_steps", "100"]`). Discrete callbacks and `terminate!` are supported; `ContinuousCallback`s and `tstops` are rejected with `ArgumentError` so silently-skipped exact-time callbacks cannot bite users. Standard SciML control knobs `adaptive`, `dtmin`, and `dtmax` are mapped onto PETSc's `TSAdapt` controller; unrecognized solve keywords are rejected with a clear `ArgumentError` rather than silently dropped. The extension currently requires `PetscReal = Float64`, real-valued in-place ODE problems, and forward integration (`tspan[1] < tspan[2]`). See the docstrings of `PETSc.TSRK`, `PETSc.TSRosW`, `PETSc.TSImplicit`, `PETSc.TSARKIMEX`, and `PETSc.TSGeneric` for details.

Known callback-lifecycle gaps (compared to OrdinaryDiffEq):

- The SciML *discrete-save* lifecycle hooks `SciMLBase.save_discretes_if_enabled!` (after `initialize!`) and `SciMLBase.save_final_discretes!` (after `finalize!`) are invoked when the loaded SciMLBase version provides them, so `DiscreteCallback` machinery that records observable state at `t0` and at the end of the solve interoperates with the standard SciML lifecycle. Callbacks that depend on per-step `save_discretes_if_enabled!` calls *during* the solve are still not driven by the manual `TSStep` loop.
- `ContinuousCallback`s are rejected with `ArgumentError`. Wrap event detection through PETSc's `TSSetEventHandler` directly if you need it.
