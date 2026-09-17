# Installation

## Using pre-built libraries

The easiest way to install the package is:
```julia
julia> ]
(@v1.12) pkg> add PETSc
```
which will install a pre-built PETSc library (`PETSc_jll`) as well as `MPI.jl` on your system. This will work both in serial and in parallel on your machine.

!!! note "Windows Users"
    Since `PETSc_jll` 3.25.4 the Windows binaries are built with MPI (`MicrosoftMPI_jll`), and
    PETSc.jl runs natively on Windows in serial and in parallel; the test suite, including the
    MPI tests, runs there in CI. One limitation of PETSc 3.25.x on Windows: `Tao` objects can
    only be created before the first `PETSc.finalize`, see `PETSc.tao_usable_after_reinitialize()`.
    The [Windows Subsystem for Linux](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux)
    remains an alternative if you need a custom PETSc build.

## Using a custom PETSc build

Sometimes, you may be interested in a PETSc installation that comes with additional external packages, or that you compiled yourself. Ensure the library is compiled as a **dynamic** (not static) library.

Use `set_library!` to configure the path once — it is stored in `LocalPreferences.toml` and no environment variables are needed afterwards:

```julia
using PETSc
PETSc.set_library!(
    "/path/to/custom/libpetsc.so";
    PetscScalar = Float64,
    PetscInt    = Int64,
)
# Restart Julia — PETSc_jll is not loaded and your library is used automatically.
```

To revert to the bundled binaries: `PETSc.unset_library!()`.

To check what is currently configured: `PETSc.library_info()`.

For a one-off session without changing persistent settings, use the `PetscLibType` constructor directly:

```julia
petsclib = PETSc.LibPETSc.PetscLibType("/path/to/custom/libpetsc.so";
                              PetscScalar=Float64, PetscInt=Int64)
PETSc.initialize(petsclib, log_view=true)
# ... your code ...
PETSc.finalize(petsclib)
```

## HPC systems

On many high-performance clusters, you will need to use the cluster's MPI installation. There are a number of options to do so - see the [Running on HPC Systems](hpc.md) page for details.
