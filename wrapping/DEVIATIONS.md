# Deliberate differences between the regenerated wrappers and the hand-edited baseline

The baseline is `src/autowrapped/` as of commit `7dbff00` (PETSc 3.24.0, after PR #263).
`golden_diff.jl --categorize` classifies every wrapper that differs from it. The categories below
are differences by design; anything the tool reports as `other` is either whitespace/hand
formatting or still to be reviewed. Counts are from the run of 2026-09-15 (6086 wrappers).

| category | count | what changed | why |
|---|---|---|---|
| `direction-new-output` | 481 | a `T *x` argument the old heuristics treated as an input (and passed the caller's handle, nulling it afterwards) is now an output and is **returned** | the manual page lists it as an output; the old code threw the result away or corrupted the caller's handle (PR #258 fixed six of these by hand) |
| `enum-direction-fix` | 302 | enum/opaque outputs are returned instead of taken as by-value inputs | same: `Type *out` was refused as an output when `Type` was not a scalar |
| `voidptr-fix` | 415 | `void *ctx` is `ctx::Ptr{Cvoid}` instead of `ctx::Cvoid` | `::Cvoid` could never be satisfied by a caller |
| `pointer-input-as-array` | 185 | `const T *x` inputs are `x::Vector{T}` instead of `x::T` | a pointer to an array was passed a scalar |
| `string-arg-fix` | 128 | `const char s[]` inputs are `s::String` instead of `s::Vector{Cchar}` | consistent with the majority of the baseline and usable |
| `placeholder-removed` | 118 | PETSc-owned output arrays without a size rule return the raw pointer (`Ptr{T}`) | the baseline called `VecGetLocalSize(petsclib, x)` with an undefined `x` (UndefVarError). Since 2026-09-16 `rules/args.toml` carries `size` rules for the ~145 arrays whose length the manual page documents (they come back as `Vector`s, handles as `Vector{IS}` etc., `char**` as `Vector{String}`); about 190 outputs with no documented length (workspaces, tabulations, CAD data, cones of whole meshes, file pointers) still return the raw pointer |
| `type-null-guard` | 59 | `*GetType` returns `""` when PETSc returns NULL | baseline segfaulted on `unsafe_string(C_NULL)` (`MatGetType` had the guard by hand) |
| `enum-out-fix` | 47 | enum outputs whose name contains `Type` are returned as the enum, not `unsafe_string`d | `unsafe_string` on an enum is wrong |
| `nullable-form` | 40 | `Union{Ptr, X}` instead of `Union{Ptr,X}` / `Union{Ptr{X}, Ptr{Nothing}}` | one spelling |
| `fnptr-ptrcvoid-fix` | 39 | callback arguments are `Ptr{Cvoid}` instead of the opaque `XFn` struct | an `XFn` struct could never hold an `@cfunction` pointer |
| `handle-ccall-fix` | 32 | `Ptr{CIS}`/`Ptr{Ptr{CIS}}` in ccall tuples instead of `Ptr{IS}` | `IS` is the Julia struct, `CIS` the C handle |
| `writeback-fix` | 31 | remaining `x.ptr = C_NULL` after a Get became `x.ptr = x_[]` | result was discarded |
| `handle-return-fix` | 13 | `Create*` returning `IS`/`SF` construct the Julia handle instead of returning a raw `Ref{IS}()[]` | consistent handle construction |
| `restore-input-fix` | 15 | `*Restore*` functions take the array back as an input | Restore never outputs |
| `scalar-byref-input` | 15 | `PetscInt *n` inputs of Restore functions are `n::PetscInt` passed as `Ref(n)` | previously passed by value |
| `byref-new` | 6 | `XDestroy(X*)` for opaque handles uses the `Union{X, Ref{X}}` + `Base.RefValue` idiom | the baseline had a hand variant (`Ref(x)`) |
| `byref-missing` | 8 | `DMGetLocalToGlobalMapping` and friends return the object instead of taking `Union{X, Ref{X}}` | return style (see `direction-new-output`) |
| `handle-out-byvalue-fix` | 6 | handle outputs no longer passed by value into a `Ptr{CX}` slot | segfault |
| `ni-alloc-removed` | 1 | caller-allocated output arrays without a known length are inputs | `Vector{T}(undef, ni)` with undefined `ni` |
| `direction-golden-output` | 27 | baseline treated a by-value input as an output (e.g. `DMDASetElementType`, `MatDenseReplaceArray`) | baseline bug |
| `docs-only` | 163 | docstring text differs | the baseline was generated from a different checkout; also the fixed manual-page block termination |
| doc links | (all) | `_doc_external("TS/...")` instead of `"Ts/..."` | manual section from the makefile, so the links resolve |

Other differences from the baseline that are not per-function:

- the docstring stub `function X(petsclib::PetscLibType, ...)` has loosened scalar types and throws
  instead of silently returning `nothing` (this surfaced two test bugs in `test/test_dmstag.jl`);
- `PetscObject` arguments are untyped, so any handle (including `VecPtr`, `MatPtr`) is accepted;
- `MPI_Comm` outputs are read through the C handle and wrapped in `MPI.Comm`;
- `*Restore*` functions accept the array a `Get` returned, or the raw pointer if the size was unknown;
- arrays of handles (`const Vec vecs[]`) are `Vector{<:AbstractPetscVec}` and converted element-wise;
- `XCreate(..., X *x)` returns `x` even when the manual page lists it as an input (`PetscSectionCreate`);
- non-const `T *x` scalar pointers are outputs even when listed as inputs (`TSIRKGetNumStages`);
- `void *ctx` documented as an output returns the pointer (`MatShellGetContext`);
- deprecated enum aliases are skipped rather than truncating the enum (`SNESConvergedReason`);
- string-enum arguments (`PCType`, `MatType`, ...) are `String`s (the baseline needed
  `Base.unsafe_convert(Ptr{Int8}, "ilu")` or the hand-written `src/string_wrappers.jl`, now removed)
  and `senums_wrappers.jl` defines the registered names as constants (`LibPETSc.PCMG == "mg"`);
- `const char *x[]` outputs (`PetscObjectGetType`, `KSPGetOptionsPrefix`, ...) return a `String`
  instead of a raw pointer; `PetscObjectGetName` returns the name (the manual page lists it as an input);
- `direction = "inout"` rules make `PetscSplitOwnership*`, `PetscSortRemoveDups*` and the `nmax` of
  `PetscOptionsGet*Array` take and return the scalar (the baseline passed an uninitialised `Ref`);
- `PCMGSetLevels` takes the optional `comms` (`C_NULL`) instead of returning garbage;
- `PetscSFBcastBegin/End`, `PetscSFReduceBegin/End`, `PetscSFFetchAndOpBegin/End` exist as
  hand-written extras (they are absent from the API snapshot);
- `PetscDraw` and `TSMonitorLGCtx` are opaque pointer handles; the baseline's
  `mutable struct PetscDraw end` placeholder made every `PetscDraw*` call fail (`Ref{PetscDraw}()`
  is an undefined reference and the handle was passed as a Julia object pointer);
- 95 header-inline functions and macros that have no symbol in `libpetsc` (`PetscStrcmp`,
  `PetscTime`, `VecSetValue`, `MatSetValue`, `PetscOptionsBegin`, ...) are excluded; the baseline
  had wrappers that failed with "could not load symbol";

- type names are mapped with the old substring replacement by default (`fix_substring_replacements = false` in `types.toml`) so that `PetscPoCintFn`-style names are reproduced; flip the flag to get correct names (they are opaque placeholders, so nothing else changes);
- the opaque type declarations that used to sit at the top of each file are collected in `opaque_types.jl`;
- functions are sorted by name within each file;
- the `XFn` placeholder structs are no longer declared (callbacks are `Ptr{Cvoid}`);
- `petsc_wrappers_version.jl` no longer records a machine path;
- `PC_wrappers.jl` is generated and included (the baseline excluded it).

Tests adapted: `test/dmplex.jl`, `test/mat.jl`, `test/snes.jl`, `test/test_dmstag.jl`; example
`examples/ex62b.jl`. Also `src/mat.jl` (`MatShellGetContext`), `src/ksp.jl`/`src/ts.jl`/`src/dm.jl`
(solution and coordinate accessors return borrowed `VecPtr` handles, `destroy!` is a no-op on them),
`src/options.jl` (NULL viewer for `PetscOptionsView`).

High-level code adapted for the return convention: `dm.jl` (`DMGetCoordinatesLocal`),
`dmplex.jl` (`DMPlexDistribute`, `DMClone`, `DMCoarsenHookAdd`, `DMGetStratumIS`, `DMGetCoarseDM`),
`ts.jl` (hand overloads of `TSSetRHSFunction`, `TSSetIFunction`, `TSSetIJacobian`,
`TSSetRHSJacobian`, `TSGetAdapt`, `TSIRKGetNumStages`, `TSGetSolution`, `TSGetSNES`, `TSGetKSP`,
`TSMonitorSet` removed; the generated ones now have the same behaviour).

## Type renames (PR #260 naming conventions)

`PetscKSP` -> `KSP`, `PetscSNES` -> `SNES`, `AbstractPetscKSP` -> `AbstractKSP`,
`AbstractPetscSNES` -> `AbstractSNES` (prologue and `rules/types.toml`), and
`AbstractPETScMemBackend` -> `AbstractPetscMemBackend` (high-level). The old names remain as
aliases until v0.6. The high-level `KSP(...)`/`SNES(...)` factories are now constructor methods of
the LibPETSc types, as §5.1 of the naming document asks.
