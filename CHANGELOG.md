# Changelog

## Unreleased

- `TSIRK` types (`-ts_irk_type gauss` and the rest) are found again after `finalize` followed by `initialize`. PETSc 3.25 leaves their registration flag set at finalize; `initialize` now resets it.
- `PetscOptions(petsclib)` throws `PetscNotInitialized` on a library that is not initialized, like every other high-level constructor. It used to succeed and return an options database that `destroy!` skipped once `initialize` ran.
- `LibPETSc.PC` is a Julia type, `PC{PetscLib}`, like `KSP` and `IS`, instead of a raw pointer alias. `KSPGetPC`, `PCCreate` and the other functions that hand out a PC return it, and every `PC*` function takes an `AbstractPC`. Code that passes a PC from one call to the next is unaffected; code that annotates a variable as `LibPETSc.PC` expecting a `Ptr`, or passes raw pointers, needs `pc.ptr`.
- `PETSc.pc(ksp)` returns the preconditioner, borrowed from `ksp`: `destroy!` on it does nothing. `set_type!`, `type_name` and `set_fieldsplit_is!` work on it.
- `set_shell_apply!(apply!, p)` and `set_shell_setup!(setup!, p)` write a `:shell` preconditioner in Julia. `apply!(y, p, x)` writes into `y`.
- Callback closures, the user context and the options applied at `solve!` are kept with the PETSc object rather than on the wrapper (docs/src/man/naming.md §18.3). A callback set through a borrowed handle, such as `set_function!(f!, snes(ts), r)`, lives as long as the solver. `snes.f!`, `snes.user_ctx`, `ts.opts` and the other former fields stay readable and writable. `user_ctx(obj)` and `set_user_ctx!(obj, ctx)` read and set the context of a `SNES` or `TS`.
- An exception thrown in a callback (SNES residual, Jacobian or convergence test, KSP right-hand side or operators, TS right-hand side, implicit function, Jacobians or monitor, `MatShell` multiply, shell preconditioner) comes out of the `solve!`, `step!` or `setup!` that ran it as the original exception. SNES, KSP and `MatShell` callbacks used to unwind through PETSc's C code, which is undefined behaviour. When the work was started through `LibPETSc` directly, the exception is logged and the `LibPETSc` call throws `PetscError`.
- **Behaviour change:** an exception in a `TS` callback comes out of `solve!` or `step!` as itself. It used to be logged and surface as a `PetscError`.
- A callback's return value is ignored. A nonzero `Integer` still fails the call, as the PETSc error code 0.5.0 read it as, and warns once; from v0.6 it is ignored. `return 0` keeps working.
- `KSP` and `MatShell` attach a finalizer on a one-process communicator, like `SNES` and `TS`.
- Every `!` function returns the object it mutates (docs/src/man/naming.md §7): `set_type!(ksp, :cg)` returns `ksp`, `solve!(x, ksp, b)` returns `x`, `assemble!(A)` returns `A`, `set_function!(f!, snes, r)` returns `snes`, and so on for about 60 functions that returned `nothing`. `set_function!` returned `0`. `setindex!`, `fill!`, `mul!` and `copyto!` on PETSc objects return the object, as their Base contracts say. `destroy!`, `restore_local_arrays!`, `set_library!` and `push!`/`pop!` on `PetscOptions` still return `nothing`; `with_local_array!` returns what its block returns.
- **Behaviour change:** `add_boundary!` and `add_natural_boundary!` return `dm` instead of PETSc's boundary number; `LibPETSc.DMAddBoundary` still returns it.
- `save_vtk!` and `vtk_merge_tensor!` are renamed `save_vtk` and `vtk_merge_tensor`: they write a file and mutate no argument. The old names warn until v0.6.
- `MPIPreferences` is a test-only dependency. `src/` never loaded it; install it yourself to select an MPI binary, as the HPC guide describes.
- Every handle knows whether it owns its PETSc object. `PetscVec`, `PetscMat`, `PetscDM`, `KSP`, `PC`, `SNES`, `TS`, `PetscOptions`, `IS`, `PF`, `Tao` and `AO` carry `ptr`, `age` and `own`, and `PETSc.owns(obj)` reads the field. Constructors take `own` as a keyword, `true` by default.
- **Behaviour change:** a handle returned by a `LibPETSc` function with `Get` in its name is borrowed, following PETSc's convention, and `destroy!` on it does nothing. That covers `snes(ts)`, `ksp(ts)`, `pc(ksp)` and the objects a callback receives, which could previously be destroyed from under the solver still using them. The `Get` functions that hand the caller a new reference (`MatGetFactor`, `DMLabelGetStratumIS`, `DMGetStratumIS`, `MatGetOrdering`, `MatGetOwnershipIS` and 31 more, listed in `wrapping/generator/rules/ownership.toml`) return an owner as before. `LibPETSc.XDestroy` still frees any handle it is given.
- `destroy!(pc)` destroys a `PC` the caller created with `LibPETSc.PCCreate`; on the borrowed one from `pc(ksp)` it still does nothing.
- `set_convergence_test!` no longer overwrites `snes.user_ctx`. The closure is kept on a field of its own, so residual and Jacobian callbacks that take `user_ctx` keep receiving it.

## v0.5.0

PETSc.jl 0.5 wraps PETSc 3.25.4 (PETSc_jll 3.25) and requires Julia 1.12. The low-level `LibPETSc` bindings are regenerated by a new rules-based generator (`wrapping/`, see `wrapping/WRAPPING.md`) that can be rerun for every PETSc release; the hand edits that used to live in `src/autowrapped/` are expressed as rules and overrides. This changes the calling convention of many `LibPETSc` functions.

### Breaking changes (low-level `LibPETSc`)

- **Outputs are returned, not written into caller-supplied handles or `Ref`s.**
  `PCCreate(petsclib, comm)` returns the `PC`; `KSPGetPC(petsclib, ksp)` returns the `PC`; `DMPlexDistribute(petsclib, dm, overlap)` returns `(sf, dmParallel)`; `ISSorted(petsclib, is)` returns a `PetscBool`; `PetscSFGetGraph(petsclib, sf)` returns `(nroots, nleaves, ilocal, iremote)`. About 480 functions changed this way (`wrapping/DEVIATIONS.md` lists every category).
- **`XDestroy` takes the handle or a `Ref` to it**:
  `PetscViewerDestroy(petsclib, viewer)`.
- **String enums are Julia `String`s.**
  `PCSetType(petsclib, pc, "ilu")` works directly; `Base.unsafe_convert(Ptr{Int8}, "ilu")` is no longer accepted, and the hand-written `String` overloads (`src/string_wrappers*.jl`) are gone. The registered names are constants: `LibPETSc.PCMG == "mg"`, `LibPETSc.KSPGMRES`, `LibPETSc.MATSEQAIJ`, ... `XGetType` returns a `String` (`""` when unset).
- **`PetscKSP` -> `KSP`, `PetscSNES` -> `SNES`**, `AbstractPetscKSP` -> `AbstractKSP`, `AbstractPetscSNES` -> `AbstractSNES`; `AbstractPETScMemBackend` -> `AbstractPetscMemBackend` (PR #260 naming conventions). The old names remain as deprecated aliases until v0.6.
- Input arguments accept abstract types (`AbstractPetscVec`, `AbstractVector{<:Number}`, `Integer`), returned handles are concrete (PR #263). A call with unsupported argument types throws instead of silently returning `nothing`.
- Callbacks and contexts are `Ptr{Cvoid}`; `const T *x` inputs are `Vector{T}`; `char[]` inputs are `String`; `const char *x[]` outputs return a `String`.
- PETSc-owned output arrays with a documented length come back as `Vector`s (arrays of handles as `Vector{IS}` etc.); those without one return the raw pointer.
- `PetscSplitOwnership*`, `PetscSortRemoveDups*` and the `nmax` of `PetscOptionsGet*Array` take and return their in/out scalar.
- 95 header-inline functions and macros without a symbol in `libpetsc` (`PetscStrcmp`, `PetscTime`, `VecSetValue`, `MatSetValue`, `PetscOptionsBegin`, ...) are no longer wrapped.
- Deprecated enum values are not emitted. `KSPSetDMActive` takes the `KSPDMActive` flag (PETSc 3.25).

### High-level API renamed (naming conventions)

The high-level interface follows `docs/src/man/naming.md` from v0.5 on. Roughly a hundred names change, and the full rename table is at the end of that page — what follows is what the changes are and how to move.

**The register and the shims.** `scripts/renames.jl` is the register: a plain list of `old => new` pairs plus the internal and unchanged-public sets. `src/deprecations.jl`, `src/public_names.jl`, `src/audit_names.jl` and `test/test_deprecations.jl` are generated from it by `scripts/generate_renames.jl`, so they cannot drift apart, and `scripts/api_surface.jl --check` (run by the test suite) fails if a binding is in none of its sets.
Every renamed name keeps a forwarding shim that **warns once per call site** — through `@warn`, not `Base.depwarn`, so it prints whatever `--depwarn` is set to — and the shims are **removed in v0.6**.

**Exports.** `PETSc` now exports only types and construction entry points:

```julia
export LibPETSc
export DMDA, DMStag, DMPlex
export PetscVec, PetscMat, PetscOptions
export KSP, SNES, TS
export petsclibs
```

v0.4 exported twelve functions and no types at all. `audit_petsc_file`, `set_petsclib`, `set_library!`, `unset_library!`, `library_info`, `AbstractPetscMemBackend`, `AbstractPETScMemBackend`, `determine_memtype`, `get_petsc_arrays`, `restore_petsc_arrays` and `dmda_star_fd_coloring` lose their export (`HostBackend` was exported but never defined).
No shim can help here: the replacements are not exported either, so `using PETSc` code qualifies the call (`PETSc.set_library!`) or imports the name. The rest of the API is marked `public`, so `names(PETSc)` reports it without exporting it.

**Construction goes through the type.** `PetscVec(petsclib, …)` replaces `VecSeq` and `as_petsc_vec`; `PetscMat(petsclib, …)` replaces `MatSeqAIJ`, `MatSeqDense`, `MatCreateSeqAIJ`, `MatSeqAIJWithArrays` and `MatAIJ`; `PetscOptions` replaces `Options`; `PetscLibType(path; …)` replaces `set_petsclib`. `KSP`, `SNES` and `TS` are types rather than factory functions, so `ksp isa PETSc.KSP` holds.

**A typed DM hierarchy.** `DMDA{L,N}`, `DMStag{L,N}` and `DMPlex{L}` are concrete types under `LibPETSc.AbstractPetscDM`, not one type with a runtime string flavour. Code annotated `::PetscDM` no longer matches — use `AbstractPetscDM`. `PETSc.narrow(dm)` turns a low-level handle into the typed one; it queries PETSc, so its return type is a wide `Union` and hot code should narrow once behind a function barrier. The nine DM-flavour `@assert`s are gone: dispatch enforces what they checked.

**Borrowed handles.** A reader that hands back a PETSc object owned by another object — `dm(ksp)`, `solution(snes)`, `local_coordinates(dm)`, `tolerances(ts)`'s vectors — returns a borrowed handle: no finalizer, and `destroy!` on it is a no-op (`PETSc.owns` tells the two apart). Nothing changed at runtime; destroying such a handle was corrupting the owner's already. `destroy` is now `destroy!`, per the mutation convention.

**Type names are `Symbol`.** `type_name` on a Vec, Mat, KSP, SNES, TS or DM returns a `Symbol`, or `nothing` when PETSc has no type for the object yet, so `type_name(ksp) == "gmres"` is now false; compare against `:gmres`. The new `set_type!(obj, :gmres)` covers Vec, Mat, KSP, SNES and DM as well as TS, and `set_type!(obj, "gmres")` warns until v0.6. This break has no shim on the reader side.

**Argument order, and `petsclib` dropped.** The written vector leads, the DM follows it, and the library is recovered from the object: `project_function!(X, dm, time, funcs, ctxs, mode)`, `project_field!(X, dm, time, U, funcs, mode)`, `global_to_local!(lvec, dm, gvec, mode)`, `local_to_global!(gvec, dm, lvec, mode)`, `l2diff(dm, time, funcs, ctxs, X)`, `add_boundary!(dm, …)`, `add_natural_boundary!(dm, …)`, `set_snes_local_fem!(dm)`, `save_vtk!(vec, filename)`, `star_fd_coloring(da)`. The v0.4 spellings forward and warn, but a call passed through `invoke` or a function reference is not caught.

**Callback setters take the callback first, only.** `set_function!`, `set_snes_jacobian!`, `set_convergence_test!`, `set_compute_rhs!`, `set_compute_operators!`, `set_rhs_function!`, `set_rhs_jacobian!`, `set_ifunction!`, `set_ijacobian!`, `set_monitor!` and `add_coarsen_hook!` no longer accept the subject-first order v0.4 also offered, so `do` syntax is always available. No shim: dispatch cannot tell the two orders apart.

**Dimension-correct returns.** `corners`, `ghost_corners`, `local_indices`, `global_indices`, `info` and `size` answer with the DM's own dimension: `lower`/`upper` are `CartesianIndex{N}`, `size`/`nextra` are `NTuple{N,Int}`, and `center`/`vertex` are keyed `x`, `y` (, `z`) by dimension. v0.4 padded everything to three, so `corners(dm2d).size[3]` returned `1` and now throws a `BoundsError`. This is also a performance fix: the tuples are built with `ntuple(…, Val(N))` and infer concretely instead of allocating. `info` drops the duplicate `s` field, renames `dof` to `ndofs` and `mpi_proc_size` to `procs`, and `ghost_corners(::DMStag)` has no `nextra` field: `DMStagGetGhostCorners` never reported one, though the v0.4 docstring promised it.

**`ownership_range(A)` is 1-based.** That was already the default; the positional `ownership_range(A, false)` still returns PETSc's numbering, warns, and is a `MethodError` in v0.6. `set_values!` spells its index parameters `rows_0b`/`cols_0b` and `star_fd_coloring` returns `row_coo_local_0b`, `col_coo_local_0b`, `perturb_cols_1b`, `coo_idxs_1b` and `local_rows_1b`, so every bulk index vector says which base it uses.

**Smaller changes.** `library_info()` returns `(; source, path, scalar, int, real)` and prints its old report from a `show` method. Argument problems raise `ArgumentError`, `DimensionMismatch` or `PetscNotInitialized` rather than `AssertionError` or a bare `error`. `PETSc.KSP(…)` / `PETSc.SNES(…)` construct the `LibPETSc.KSP` / `LibPETSc.SNES` types.

The rename table, and the reasoning behind each rule, are in [`docs/src/man/naming.md`](docs/src/man/naming.md); the C-function-to-Julia-name lookup is the generated "C to Julia name index" page in the manual.

### Changed

- `ForwardDiff`, `UnicodePlots`, `Statistics` and `Pkg` are no longer dependencies of the package: nothing in `src/` used them, and the tests and examples that do now pull them through `[extras]` and `examples/Project.toml`. Installing PETSc.jl resolves 86 packages instead of 121, which cuts cold precompilation by roughly a factor of five.

- Only one `PetscInt` width of `PETSc_jll` is loaded per process (#241): `PETSc.petsclibs` holds the four scalar variants of `Int64` by default, `PETSc.set_petscint!(Int32)` switches to the `Int32` libraries on the next session. Loading both widths cross-binds the identically named HYPRE/SuperLU_DIST symbols of the two integer ABIs. Code indexing `petsclibs[5:8]` breaks.

### Added

- PetscSF communication: `PetscSFBcastBegin/End`, `PetscSFReduceBegin/End`, `PetscSFFetchAndOpBegin/End` (hand-written; `getAPI.py` rejects `MPI_Datatype` functions).
- `PetscDraw` and `TSMonitorLGCtx` handles work (they were empty placeholder structs).
- Windows CI with MPI (PETSc_jll 3.25.4 ships MPI-enabled Windows binaries).
- Documentation pages for every wrapper file, and the maintainer guide "Regenerating the wrappers".
- Test coverage: `test/wrapper_quality.jl` (every generated method infers a concrete return type; representatives are allocation free), `test/wrapper_leaks.jl`, an ambiguity budget, `test/low_level_petscsf.jl`.

### Fixed

- `PetscBool` is one byte, matching PETSc's `typedef bool PetscBool` (since 3.24). The previous 32-bit type read three bytes PETSc never wrote, so a `PETSC_FALSE` output could come back as true (#268), `Vector{PetscBool}` arguments had the wrong stride and the `MatFactorInfo`, `PetscFEGeom` and `PetscEventPerfInfo` struct layouts were off.

- Re-initialising PETSc after `finalize` works with Tao, TaoTerm and TSTrajectory (PETSc 3.25.x does not reset their `RegisterAllCalled` flags; on Windows the internal symbols are not exported, see `PETSc.tao_usable_after_reinitialize()`).
- Re-initialising PETSc no longer overwrites three bytes of PETSc's memory: the `TaoRegisterAllCalled`, `TaoTermRegisterAllCalled` and `TSTrajectoryRegisterAllCalled` flags were reset with a 4-byte store, while `PetscBool` is 1 byte.
- `unsafe_local_array` no longer touches a destroyed vector or a finalized library.
- Double free of DMDA coordinate vectors; `PCMGSetLevels` reading uninitialised communicators.
- Block assignment `A[rows, cols] = block` on a `PetscMat` works (#248, #271). It threw a `MethodError`, and the values were laid out column by column where `MatSetValues` reads them row by row. A block whose size does not match the index ranges raises a `DimensionMismatch`.
- `PetscOptions` records the initialize/finalize cycle it was created in, like `PetscVec`, `PetscMat`, `KSP`, `SNES`, `TS` and the DMs, so `destroy!` and its finalizer no longer call `PetscOptionsDestroy` on an object from an earlier cycle (#270).

### Removed

- The old PythonCall-based generator (`wrapping/generatejuliabindings.jl`), `src/deprecated/`, `src/startup.jl`, `src/string_wrappers*.jl`.
- `LibPETSc.petsc_version` and `PETSc.SUPPORTED_PETSC_VERSIONS`, added in 0.4.14 for running on PETSc 3.22 and 3.25. 0.5 binds PETSc 3.25 only; `PETSc.check_wrappers_version(petsclib).installed_version` returns the library version.
