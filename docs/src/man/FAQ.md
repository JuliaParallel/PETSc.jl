# Frequently Asked Questions


## 1. Can I use my own PETSc library?
Yes. You do need to compile PETSc as a dynamic (shared) library, and `MPI.jl` must be configured to be compatible with the MPI used to compile PETSc. If you run this on an HPC system and don't know the configuration options, compile one of the PETSc examples and run it with `-log_view`; the output lists all configuration options used on that machine.

Please note that the version of PETSc should be compatible with the version used for the wrappers.

The recommended approach is `set_library!`, which stores the path persistently in `LocalPreferences.toml` — no environment variables needed:

```julia
using PETSc
PETSc.set_library!("/path/to/libpetsc.so"; PetscScalar=Float64, PetscInt=Int64)
# Restart Julia — the custom library is used automatically from here on.
```

To revert to the bundled `PETSc_jll` binaries:
```julia
PETSc.unset_library!()
```

For a one-off session without changing the persistent preference, use `set_petsclib` directly:
```julia
petsclib = PETSc.set_petsclib("/path/to/libpetsc.so"; PetscScalar=Float64, PetscInt=Int64)
```

## 2. Help, my code crashes?
That is very possible. If you provide a *short* minimum working example (MWE), feel free to open an issue on the github repo, so we can check it out.

## 3. What about the garbage collector in Julia?
Using the GC in combination with MPI code is a tricky business. The users are therefore responsible to free PETSc objects in the code, as shown in the various examples. We have some help for that:
1. The julia function `PETSc.audit_petsc_file("path/to/your/file.jl")` which scans your julia file and tries to guess whether objects are destroyed.
2. You can initialize the PETSc library with `log_view=true`. At the end of the code, it will give an  

## 4. Is it compatible with GPUs?
The precompiled `PETSc_jll` binaries are currently not compatible with GPUs. You can, however, use your own PETSc compilation and things should be fine. 

## 5. I really like this package!
Thanks a lot, we appreciate if you give us a star on github!

## 6. How do I cite this in publications?
We are planning to submit a JOSS paper once the current release is sufficiently stable. Meanwhile, you can cite the [zenodo release](https://doi.org/10.5281/zenodo.18274809).