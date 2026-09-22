# Upgrading from 0.4

PETSc.jl 0.5 is a breaking release. Most of the high-level API has new names, the low-level `LibPETSc` functions return their outputs, and the package needs Julia 1.12 and PETSc 3.25. This page is the order in which to update code written for 0.4. Every change is listed in the [Release notes](release_notes.md), and every renamed function in the rename table at the end of [Naming Conventions](naming.md).

## 1. Requirements

- **Julia 1.12 or newer.** On an older Julia, the package manager keeps installing 0.4.
- **PETSc 3.25.** `PETSc_jll` 3.25 is installed automatically. A library configured with `PETSc.set_library!` must be a PETSc 3.25 build; `PETSc.check_wrappers_version()` warns when it is not.
- **One `PetscInt` width per process.** `PETSc.petsclibs` holds four libraries (`Float64`, `Float32`, `ComplexF64`, `ComplexF32`) with `Int64` indices. Code that indexes `petsclibs[5:8]`, or asks `PETSc.getlib` for an `Int32` library, fails. To use `Int32` indices, call `PETSc.set_petscint!(Int32)` once and restart Julia (see [Installation](installation.md)).

## 2. Run your code and fix the warnings

About a hundred functions were renamed. The old names still work until 0.6 and print a warning the first time each call site runs, for example

```
┌ Warning: getcorners is deprecated, use corners
```

These warnings appear without `--depwarn=yes`. Run your code or its test suite once, and replace every name a warning reports.

## 3. Qualify names that are no longer exported

`using PETSc` now brings in only the types you construct:

```julia
LibPETSc
DMDA, DMStag, DMPlex
PetscVec, PetscMat, PetscOptions
KSP, SNES, TS
petsclibs
```

Everything else, including functions 0.4 exported such as `set_library!`, `library_info` and `audit_petsc_file`, fails with an `UndefVarError`. Write `PETSc.set_library!(…)`, or import the names you use: `using PETSc: set_library!, initialize, finalize`.

## 4. Changes that print no warning

These changes have no deprecation warning, because the name stayed the same or dispatch cannot tell the old call from the new one. Check your code for each of them.

**Callback setters take the callback first.** 0.4 accepted both orders; 0.5 accepts only the callback-first one, which is the one `do` syntax needs. This applies to `set_function!`, `set_snes_jacobian!`, `set_convergence_test!`, `set_compute_rhs!`, `set_compute_operators!`, `set_rhs_function!`, `set_rhs_jacobian!`, `set_ifunction!`, `set_ijacobian!`, `set_monitor!` and `add_coarsen_hook!`.

```julia
# 0.4
PETSc.setfunction!(snes, f!, r)

# 0.5
PETSc.set_function!(f!, snes, r)
PETSc.set_function!(snes, r) do fx, snes, x
    # ...
end
```

**Type names are `Symbol`s.** `type_name` (0.4: `gettype` or `type`) returns a `Symbol`, or `nothing` when PETSc has not set a type yet. A comparison with a `String` is now `false` without any error.

```julia
# 0.4
PETSc.gettype(snes) == "newtonls"

# 0.5
PETSc.type_name(snes) === :newtonls
```

**DM queries return the DM's own dimension.** 0.4 padded every answer to three dimensions. In 0.5, `corners`, `ghost_corners`, `local_indices`, `global_indices`, `info` and `size` return values of the DM's dimension: for a 2D DMDA, `corners(da).lower` is a `CartesianIndex{2}` and `corners(da).size` an `NTuple{2,Int}`, so `corners(da).size[3]` throws a `BoundsError`. `info` renames `dof` to `ndofs` and `mpi_proc_size` to `procs`, and drops `s`, which duplicated `stencil_width`.

**DMs have their own types.** `DMDA`, `DMStag` and `DMPlex` are concrete types below `LibPETSc.AbstractPetscDM`. A method annotated `::PetscDM` no longer matches them; annotate with `LibPETSc.AbstractPetscDM` or with the concrete type. A DM that comes back from PETSc without a known type (for example `PETSc.dm(ksp)`) can be turned into the typed one with `PETSc.narrow(dm)`.

**Errors have specific types.** Invalid arguments raise `ArgumentError` or `DimensionMismatch` where 0.4 failed an `@assert`, and creating a PETSc object with `PetscVec`, `PetscMat`, `KSP`, `SNES`, `TS` or a DM constructor on a library that is not initialized raises `PETSc.PetscNotInitialized`. Code that caught `AssertionError` must catch these instead.

**Handles you do not own.** Readers such as `PETSc.dm(ksp)`, `PETSc.solution(snes)` and `PETSc.local_coordinates(dm)` return a handle that belongs to the object they were called on. `destroy!` on such a handle does nothing (`PETSc.owns` tells you which is which). Destroying it was never correct, and code that did so should stop.

**Reordered arguments.** Functions that write into a vector take that vector first, and the DM second:

```julia
# 0.4
PETSc.dm_global_to_local!(gvec, lvec, dm)
PETSc.dm_local_to_global!(lvec, gvec, dm)

# 0.5
PETSc.global_to_local!(lvec, dm, gvec)
PETSc.local_to_global!(gvec, dm, lvec)
```

The 0.4 names still work and warn, but only when called directly: a 0.4-order call passed through `invoke`, or through a variable holding the function, is not translated. The same applies to `project_function!` and `project_field!`.

## 5. Low-level `LibPETSc` code

The `LibPETSc` wrappers are regenerated for PETSc 3.25 and follow a new calling convention:

- **Outputs are returned.** A C output argument is no longer passed in as a `Ref`:

  ```julia
  # 0.4
  s = Ref{LibPETSc.PetscSection}()
  LibPETSc.PetscSectionCreate(petsclib, comm, s)
  section = s[]

  # 0.5
  section = LibPETSc.PetscSectionCreate(petsclib, comm)
  ```

  Several outputs come back as a tuple, for example `nroots, nleaves, ilocal, iremote = LibPETSc.PetscSFGetGraph(petsclib, sf)`.
- **Type names are `String`s.** `LibPETSc.PCSetType(petsclib, pc, "ilu")` works directly, and the registered names are constants such as `LibPETSc.PCILU` and `LibPETSc.KSPGMRES`. `XGetType` functions return a `String`.
- **Renamed handle types.** `PetscKSP` and `PetscSNES` are now `KSP` and `SNES` (`AbstractPetscKSP` and `AbstractPetscSNES` are `AbstractKSP` and `AbstractSNES`). The old names still work until 0.6.
- **Functions PETSc defines only in its headers are gone.** About 95 macros and inline functions with no symbol in the library (`VecSetValue`, `MatSetValue`, `PetscStrcmp`, `PetscTime`, …) are no longer wrapped. Use the array versions (`VecSetValues`, `MatSetValues`) or plain Julia.
- **`PetscBool` is one byte**, as it is in PETSc since 3.24. Code that uses `LibPETSc.PetscBool` needs no change. Code that writes a PETSc boolean as `Cint` or `Int32` in its own `ccall` or `unsafe_store!` must use `LibPETSc.PetscBool` instead.

The [Release notes](release_notes.md) list the remaining changes, and [`wrapping/DEVIATIONS.md`](https://github.com/JuliaParallel/PETSc.jl/blob/main/wrapping/DEVIATIONS.md) lists every category of changed wrapper.

## 6. Looking up a name

- **An old high-level name:** the rename table at the end of [Naming Conventions](naming.md).
- **A PETSc C function:** the [C to Julia name index](api_index.md).
- **Why a name is what it is:** the rules in [Naming Conventions](naming.md).
