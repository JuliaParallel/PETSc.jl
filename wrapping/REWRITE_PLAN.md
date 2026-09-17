# Plan: a reproducible generator for `LibPETSc` (v0.5 rewrite)

Status (2026-09-15, end of day): M0-M3 essentially done on the local `v0.5` branch. The generator
lives in `wrapping/generator/` (see `WRAPPING.md`), `src/autowrapped/` is regenerated from PETSc
3.24.0, the high-level layer and tests are adapted to the return convention (see `DEVIATIONS.md`),
`test/wrapper_quality.jl` and `.github/workflows/wrappers.yml` exist. All serial test files pass
locally with the regenerated wrappers (examples and MPI runs pending at the time of writing). M5 done: `src/autowrapped` is
generated from PETSc 3.25.4 and `PETSc_jll = "3.25"` (full suite run against it in progress at the
time of writing). Open: M4 hygiene items (delete the old generator files in `wrapping/`,
`src/deprecated/`, `src/startup.jl`; revisit the remaining `other` category in `DEVIATIONS.md`).

## 1. What the analysis found

### 1.1 The generator and its output have drifted apart

- `src/autowrapped/` holds 57 `*_wrappers.jl` files with **6102 `@for_petsc` methods** for PETSc **3.24.0**.
- The checked-in `wrapping/generatejuliabindings.jl` is **not** the version that produced them. Its
  array-output branch emits `Ref{Ptr{T}}` in the ccall tuple and a `PetscArray(...)` return, which appears
  in only ~20 functions (`VecGetArray{1,2,3,4}d*`). The other ~6080 functions were produced by an earlier
  version (commit `0bbf53a`) whose template uses `Ptr{...}` for every star in the ccall tuple and
  `unsafe_wrap(Array, p_[], VecGetLocalSize(petsclib, x); own = false)` as a size placeholder.
  The new generator must target the `0bbf53a` template, not the checked-in one.
- The generator has a type-renaming bug: `r"int" => "Cint"` and `"char" => "Cchar"` are applied inside
  identifiers, producing `PetscPoCintFn`, `DMPoCintLocationType`, `PCRiCchardsonConvergedReason`, `hCsize_t`
  in five wrapper files. All are opaque placeholders, so fixing the names is behaviour-neutral.
- Docstring lookup (`find_doc_strings.jl`) walks all ~2800 C files once **per function**. That is
  O(functions x files) and is the main reason a full regeneration takes so long.

### 1.2 Roughly 600 wrappers were edited by hand, in nine recurring patterns

About 88-90 % of the wrappers are byte-identical to the `0bbf53a` template. The rest fall into these
categories (counts approximate, from the analysis of `src/autowrapped/` and PRs #237, #254-#261):

| Cat. | Pattern | Count | Example |
|---|---|---|---|
| A | `XDestroy(X*)` for non-custom opaque handles: `x::Union{X, Ref{X}}`, `x_ = x isa Base.RefValue ? x : Ref{X}(x)` | 76 (+~90 `Get*` with same idiom) | `VecTaggerDestroy`, `DMGetLocalToGlobalMapping` |
| B | `Get*` writing a handle back into a passed object: `x.ptr = x_[]` instead of `x.ptr = C_NULL` | ~100 (+13 rewritten to return a new object) | `SNESGetSolution`, `DMGetCoordinates`, `MatGetLocalToGlobalMapping` |
| C | PETSc-owned output arrays with a computed size (sibling call, other output, `n+1`, dense `(m,n)`) and `Ref{Ptr{T}}(C_NULL)` + NULL guard | ~70 | `DMPlexGetCone`, `ISGetIndices`, `MatDenseGetArray`, `MatGetRowIJ`, `DMCreateFieldIS` |
| D | Caller-allocated output arrays sized by the real length argument instead of `ni` | 6 (+`PetscOptionsGetString` rewrite) | `ISLocalToGlobalMappingApply` |
| E | Nullable inputs: `Union{Ptr, Vector{T}}`, `Union{PetscVec, Ptr}`, `Union{Ptr, String}`, `Union{XFn, Ptr}`, ctx `Ptr{Cvoid}` | ~100 signatures | `MatMPIAIJSetPreallocation`, `KSPSolve`, `SNESSetFunction` |
| F | Abstract handle types in input positions so `VecPtr`/`MatPtr`/`MatShell`/typed DMs pass | 23 by hand; **all** inputs in PR #263 | `VecDestroy`, `MatMult`, `KSPSetOperators` |
| G | `PetscArray`/`OffsetArray` based Get/Restore pairs with "already restored" guard | ~20 | `VecGetArray2d`, `DMStagVecGetArray` (fully hand-written) |
| H | Strings: `VecType/MatType = Cstring`, `MatGetType` NULL guard, `char**` outputs as `Ref{Ptr{Cchar}}`, `PetscMemType` not stringified | ~15 | `PetscObjectGetType` |
| I | Struct out-args by reference | 1 | `MatGetInfo(..., info::Ref{MatInfo})` |

Support files: `petsc_library.jl` = `wrapping/prologue.jl` plus hand edits (extra fields on `PetscKSP`,
`PetscSNES`, `TS`; `Base.unsafe_convert(::Type{Ptr{CIS}}, ::AbstractIS)`; `PC_wrappers.jl` commented out).
`typedefs_wrappers.jl` has a hand-written 31-line `primitive type PetscBool 32 end` block. `petscarray.jl`
is fully hand-written. `enums_wrappers.jl`, `senums_wrappers.jl`, `struct_wrappers.jl` are generated with 1-2
hand lines each.

### 1.3 Known-wrong wrappers that nobody has touched yet (~560)

These are template output that is incorrect but has not been fixed. The new generator will reproduce
them first (to prove fidelity) and fix them afterwards through rules, so the fix is a reviewable diff:

| Problem | Count |
|---|---|
| `Vector{T}(undef, ni);  # CHECK SIZE!!` with a wrong or nonexistent `ni` | 96 |
| `unsafe_wrap(..., VecGetLocalSize(petsclib, x))` size placeholder | 189 |
| Non-Destroy functions that null a passed handle and throw the result away (`MatDuplicate`, `MatConvert`, `SNESGetKSP`, ...) | 236 |
| `char **` outputs typed as bare `Cchar` | 39 |
| `PC_wrappers.jl` excluded entirely because `PC` is a raw `Ptr{_n_PC}` | 1 file |

### 1.4 What the rest of the package depends on (must not change)

- The prologue contract: `ptr` is the first field of every handle struct, `age::Int`, constructors
  `PetscX{L}(ptr, age=0)` and `PetscX(ptr, lib, age=lib.age)`, `CVec`/`CMat`/... aliases, the abstract types,
  `PetscViewer = Ptr{Cvoid}`, `KSPConvergedReason`, `PetscMemType`, and the export list in `src/LibPETSc.jl`.
- `@for_petsc` semantics (`$PetscLib`, `$UnionPetscLib`, `$PetscInt`, `$petsc_library`, eval in caller module);
  ~35 high-level definitions in `snes.jl`, `ksp.jl`, `mat.jl`, `dmplex.jl`, `ts.jl` use it.
- 272 `LibPETSc.*` names used by `src/`, plus ~220 more used only by tests/examples; the test suite is the
  executable form of that list.
- Return conventions: creators return a wrapped handle with `petsclib.age`; multiple outputs are positional
  tuples in C argument order (`DMDAGetInfo` 13-tuple, `DMStagGetCorners` 9-tuple); `VecGetArray*` return a
  plain `Vector`; `*AndMemType` return `(array, mtype)`; `PetscOptionsGetString` returns `false` or `String`.
- Hand-written overloads in `src/ts.jl` and `src/string_wrappers*.jl` extend `LibPETSc.X` functions; a
  regenerated method with the same signature would silently overwrite them.
- Docs: `docs/src/man/*_lowlevel.md` pull docstrings with `@autodocs Pages=["autowrapped/Vec_wrappers.jl", ...]`,
  so file names and the `_doc_external` docstring footer matter.
- Dead code that can go: `src/startup.jl`, `src/deprecated/` (22k lines, never included),
  `PetscViennaCLIndices_wrappers.jl` (never included).

### 1.5 PR #263

One commit: every input argument of a custom handle type becomes its abstract supertype
(`PetscVec` -> `AbstractPetscVec`, `Vector{PetscVec}` -> `Vector{<:AbstractPetscVec}`), return positions keep
the concrete type, the old generator gets an `abstract_arg_type` helper, and `test/wrapper_signatures.jl`
fails if any LibPETSc method still takes a concrete wrapper type. This supersedes category F and becomes a
plain generator rule. **Merge #263 into `v0.5` before freezing the golden baseline** (section 3).

### 1.6 Environment

- No PETSc source tree with `getAPI.py` exists on this machine (local trees are <= 3.22.5; `getAPI.py`
  first appears in 3.23). Network works; a shallow clone of `v3.24.0` is ~150 MB, the tarball 16 MB.
- `getAPI.py` in 3.24.x: `config/utils/getAPI.py`, `getAPI()` with cwd = PETSC_DIR, returns a 9-tuple.
  In **3.25.x it moved to `lib/petsc/bin/getAPI.py`**, takes the directory as argument, and returns a
  10-tuple (new `functiontypedefs`). The generator must handle both layouts.
- Julia 1.10/1.11/1.12/1.13 installed; `python3` 3.10 (stdlib only, which is all `getAPI.py` needs).
  PythonCall 0.9.28 pinned in `wrapping/Manifest.toml` is not in the depot.
- PETSc_jll 3.25.4 artifacts are in the depot; `Project.toml` still pins `PETSc_jll = "3.22"` while the
  wrappers say 3.24.0.

### 1.7 `REGENERATING.md` from the PR #263 discussion

filoferra attached `REGENERATING.md` to PR #263 (kept in this repository until the old generator was removed; see git history). It
matches categories A-F and H above and adds these points, which become explicit rules or checks:

- Output-ness must come from the man-page `Input Parameters:` / `Output Parameters:` sections, not from
  the star count. The star heuristic caused PR #254 (`TSGetConvergedReason`, `TSGetTolerances` took outputs
  as inputs) and the hand bindings for `TSGetSNES`/`TSGetKSP` in `src/ts.jl`. Once the generator gets this
  right, those hand bindings in `src/ts.jl` can be retired (they would otherwise collide).
- Never emit `x isa Ref ? ...`; `Ptr{T} <: Ref{T}`, so it must be `isa Base.RefValue` (82 places today).
  Regression check: `grep -rn "isa Ref ?" src/autowrapped/` must be empty.
- The `VecGetLocalSize(petsclib, x)` placeholder is a hard `UndefVarError` in 126 functions where `x` is not
  an argument; the rule is "size from the paired count output or the documented size query, otherwise emit
  nothing and hand-write", never a call that cannot run.
- Type-name aliases: 2 are `Cstring` (`MatType`, `VecType`), 58 are `Ptr{Cchar}`. The `AbstractString`
  overloads in `src/string_wrappers*.jl` stay hand-written unless the generator emits one per alias; the
  plan keeps them hand-written in M2 and revisits in M4.
- Struct fields in `struct_wrappers.jl` (`JacActionCtx`, `DMDALocalInfo`, `TSMonitorDMDARayCtx`) stay
  concrete; only function inputs are widened.
- Hand-written multi-line wrappers also exist in `DMaddons_wrappers.jl` (`DMProjectFunction`,
  `DMComputeL2Diff`), and `Vecs_wrappers.jl` documents `VecSetValues` with a bare signature instead of a
  stub. Both go to `overrides/`.
- Prologue drift, five items: `computerhs!`/`computeops!`/`opts` on `PetscKSP`; `user_ctx`/`opts` on
  `PetscSNES`; `Base.unsafe_convert(::Type{Ptr{CIS}}, ::AbstractIS)`; `const DMLabel = Ptr{Cvoid}`; the
  commented-out `include("PC_wrappers.jl")`. These move back into `prologue.jl` in M0. `local_types.jl` is
  unused and is deleted in M3.
- Ambiguity budget: `length(detect_ambiguities(PETSc; recursive = true))` is 161 after #263; the
  regenerate-and-diff CI job also asserts this number does not grow.

### 1.8 Target version: PETSc 3.25.4

The wrappers to ship are for **PETSc 3.25.4**, matching the PETSc_jll binaries being built for it
(2026-09-15). Consequences:

- `getapi_dump.py` treats the 3.25 layout (`lib/petsc/bin/getAPI.py`, `getAPI(dir)`, 10-tuple with
  `functiontypedefs`) as the primary path; the 3.24 layout is only needed for the fidelity step.
- The fidelity step (section 3) still runs on the 3.24.0 source, because that is the only PETSc version the
  current hand-edited output corresponds to; it is the sole way to prove the rules reproduce the hand fixes.
  It costs one 16 MB tarball. Immediately afterwards the release workflow (section 5) is run for 3.25.4, and
  `apidiff.jl petsc-3.24.0.json petsc-3.25.4.json` explains every function that differs.
- Milestone M5 therefore moves directly after M4 (regenerated output committed) and before the hygiene fixes
  in M4-bis, and `Project.toml` gets `PETSc_jll = "3.25"` in the same PR. The full suite runs against the
  new binaries, which also validates the wrappers on Int32/Float32/complex variants of 3.25.4.
- The GitLab tag `v3.25.4` exists, and PETSc_jll 3.25.4 artifacts (with `include/petscversion.h`) are
  already in the local depot, so the target library is available for testing today.

### 1.9 `docs/src/man/naming.md` (PR #260)

PR #260 proposes the v0.5 high-level naming conventions. It declares `LibPETSc` out of scope (verbatim C
names), but four of its rules reach into the prologue and the generator:

- **Type renames in the prologue** (§5.2, §5.5): `PetscKSP` -> `KSP`, `PetscSNES` -> `SNES`,
  `AbstractPetscKSP` -> `AbstractKSP`, `AbstractPetscSNES` -> `AbstractSNES`; `PetscVec`, `PetscMat`,
  `PetscDM`, `PetscOptions` and the bare `TS`, `IS`, `AO`, `PF`, `Tao` stay. Every KSP/SNES wrapper
  signature changes with it, so the custom handle table (`types.toml`: C name -> Julia struct, C alias,
  abstract type) is the single place holding these names. The rename is then one table edit plus a
  regeneration, done in its own PR after fidelity, together with the prologue change and the high-level
  shims. It must not be mixed into the fidelity step.
- **Typed DM hierarchy** (§5.3, §5.4, closing paragraph): `DMDA{L,N}`, `DMStag{L,N}`, `DMPlex{L}` subtype
  `AbstractPetscDM`, and `PetscDM{L}` stays as the low-level handle returned by creators. The generator
  therefore emits `AbstractPetscDM` for inputs and `PetscDM` for returns, which PR #263 already does.
- **Borrowed handles** (§3.3): readers return `own = false`. The low-level wrappers construct handles
  without finalizers already; the generator can additionally emit a `_doc_borrowed` marker in the docstring
  of Get-functions that hand back a PETSc-owned object (the "writeback" and non-Create "handle out" kinds),
  so the high-level `_doc_borrowed` CI check has something to build on. Optional, decide in M4.
- **Type names are Strings at the C boundary** (§3.1): the `AbstractString` overloads in
  `src/string_wrappers*.jl` stay hand-written; `Symbol` handling is purely high-level.

Also consistent with §17: `src/deprecated/` is an older unrelated API and can be deleted; the new
`src/deprecations.jl` is a different file.

## 2. Design of the new generator

Goal: `julia wrapping/generate.jl --petsc-dir <src>` reproduces `src/autowrapped/` exactly, and rerunning it
on the next PETSc release changes only the functions whose C API or docstring changed.

### 2.1 Pipeline

```
PETSc source --python3 getapi_dump.py--> api/petsc-<ver>.json   (API snapshot: functions, args, enums, ...)
PETSc source --docindex.jl------------> docstrings index        (one pass over src/**/*.c, in memory)
api.json + docs + rules/*.toml + overrides/*.jl --generate.jl--> src/autowrapped/*.jl
api/petsc-<old>.json vs api/petsc-<new>.json --apidiff.jl-----> report: new / removed / changed functions,
                                                                 rules that no longer match anything
```

- **Drop PythonCall.** A ~30-line `wrapping/getapi_dump.py` imports `getAPI` from either location, runs it,
  and writes plain JSON. The Julia side needs stdlib `TOML` for rules plus one JSON package (JSON3) in
  `wrapping/Project.toml`. This removes CondaPkg, makes the run reproducible, and the JSON snapshot doubles
  as the API manifest for release-to-release diffs.
- **Docstring index in one pass**: scan all `/*@ ... @*/` blocks once into `Dict{String,Vector{String}}`.
  Reuse the existing cleaning logic from `find_doc_strings.jl` (input/output parameter extraction,
  note removal, backticks), but as pure functions on the indexed block.

### 2.2 Emission model

Each function goes through: classify arguments -> apply rules -> render with the `0bbf53a` template.
Argument classification produces one of a small set of **kinds**, each with a fixed rendering:

| Kind | Julia signature | ccall type | init | extract |
|---|---|---|---|---|
| scalar in | `x::T` | `T` | | |
| scalar out | (return) | `Ptr{T}` | `x_ = Ref{T}()` | `x = x_[]` |
| type-name out | (return) | `Ptr{XType}` | `Ref{XType}()` | `unsafe_string(x_[])` |
| handle in (custom) | `x::AbstractPetscVec` | `CVec` | | |
| handle out create | (return) | `Ptr{CVec}` | `Ref{CVec}()` | `x = PetscVec(x_[], petsclib)` |
| handle destroy (custom) | `x::AbstractPetscVec` | `Ptr{CVec}` | `Ref(x.ptr)` | `x.ptr = C_NULL` |
| handle writeback (custom) | `x::AbstractPetscVec` | `Ptr{CVec}` | `Ref(x.ptr)` | `x.ptr = x_[]` |
| opaque by-ref | `x::Union{X, Ref{X}}` | `Ptr{X}` | `x isa Base.RefValue ? x : Ref{X}(x)` | |
| caller-alloc array out | (return) | `Ptr{T}` | `Vector{T}(undef, <len>)` | |
| petsc-owned array out | (return) | `Ptr{Ptr{T}}` | `Ref{Ptr{T}}(C_NULL)` | `unsafe_wrap(Array, x_[], <size>; own=false)` + guard |
| array in | `x::Vector{T}` | `Ptr{T}` | | |
| array in restore | `x::Vector{T}` | `Ptr{Ptr{T}}` | `Ref(pointer(x))` | |
| string in | `x::String` | `Ptr{Cchar}` | | |
| char buffer + len | `x::Vector{Cchar}`, `len::Csize_t` | `Ptr{Cchar}` | | |
| `char**` out | (return) | `Ptr{Ptr{Cchar}}` | `Ref{Ptr{Cchar}}()` | `unsafe_string` |
| `char***` out | (return) | `Ptr{Ptr{Ptr{Cchar}}}` | `Ref{...}(C_NULL)` | loop `unsafe_string(unsafe_load(p,i))` |
| function pointer | `f::XFn` (or `Union{XFn,Ptr}`) | `Ptr{XFn}` | | |
| void ctx | `ctx::Cvoid` / `Ptr{Cvoid}` | `Ptr{Cvoid}` | | |
| struct by ref | `s::Ref{S}` | `Ptr{S}` | | |
| verbatim | whole body from `overrides/<Name>.jl` | | | |

Nullable (`Union{..., Ptr}`) is an orthogonal flag on any input kind. Every kind is rendered by exactly one
function, so the template is in one place and formatting is deterministic (tabs in body, 4-space `@chk`,
trailing space after `end`, same docstring header as today).

### 2.3 Rules live in `wrapping/rules/`, never in the output

All hand knowledge becomes declarative TOML, keyed by function name (and argument name), so it survives
regeneration and can be checked against the API snapshot:

- `files.toml`: class -> output file map (keep current file names for the docs), the "Sys" catch-all,
  excluded functions (with the reason), the include order for `petsc_library.jl`.
- `types.toml`: C -> Julia type map (word-boundary regexes, fixing the `PoCint` bug), custom handle types
  (`Vec -> PetscVec/CVec/AbstractPetscVec`, ...; add `PC` here later), string-enum overrides
  (`VecType = Cstring`), extra hand-declared types that must not be auto-declared.
- `args.toml`: per-function argument overrides: `direction = "in"|"out"|"writeback"|"destroy"`,
  `nullable = true`, `size = "DMPlexGetConeSize(petsclib, dm, p)"`, `size = "(Int(m), Int(n))"`,
  `len = "N"`, `byref = true`, `nullinit = true`, `guard = "done == PETSC_TRUE"`. This is where
  categories A-E, H, I are encoded. Generic rules handle the common cases (any `XDestroy(X*)`,
  any `Get*` with a `Type *out`, `XType*` outputs, `char**`) so the per-function table stays small.
- `overrides/*.jl`: verbatim bodies for the ~30 functions whose shape cannot be derived from the API
  (`DMStagVecGetArray/RestoreArray`, `VecGetArray{1..4}d*`, `MatGetRowIJ/RestoreRowIJ`, `MatSeqAIJGetArray`,
  `MatDenseGetArray*`, `PetscOptionsGetString`, `DMCreateFieldIS`, `ISColoringGetColors`, `TSGetTolerances`,
  `MatGetType`, `MatMPIAIJGetSeqAIJ`, `MatGetLocalToGlobalMapping`). The generator still emits the docstring
  from the C source and checks that the override's recorded C signature matches the snapshot, so a changed
  function is flagged.
- `prologue.jl` and `petscarray.jl`: hand-written, copied verbatim (the hand edits currently in
  `petsc_library.jl` move back into `prologue.jl` so there is one source of truth). The `PetscBool` block
  moves out of `typedefs_wrappers.jl` into the prologue.

Every rule entry that references a function or argument absent from the snapshot is reported as stale;
this keeps the rule set honest across releases.

### 2.4 Determinism and change detection

- Functions sorted by name within each file; header type declarations sorted; no timestamps or machine
  paths anywhere in the output (the current version file records `/Users/kausb/Downloads/petsc`; drop it).
- `petsc_wrappers_version.jl` stays, generated from `include/petscversion.h`.
- `apidiff.jl` compares two API snapshots and prints new, removed and signature-changed functions, and
  which rules/overrides are affected. This is the "what changed in 3.25" report to read before regenerating.
- Idempotency: running the generator twice yields a zero diff; a CI job (16 MB tarball of the pinned PETSc
  tag) regenerates and fails on any diff against the committed files.

### 2.5 Quality checks: type stability, allocations, leaks

Every argument kind has one renderer, so one representative wrapper per kind (times all `petsclibs`)
covers the template; functions with rules or overrides get their own entry. These run against the scratch
output in steps 2-3 and join the test suite as `test/wrapper_quality.jl` from step 4.

- **Inference**: `Test.@inferred` on the representatives; return types must be concrete. JET `@report_opt`
  on the same set to catch runtime dispatch in bodies (`Union{X, Ref{X}}` handling, `PetscArray` fields
  `data::AbstractArray` and untyped `ptr`, which should become concrete type parameters).
- **Ambiguities**: `detect_ambiguities(PETSc; recursive = true)` held at 161.
- **Allocations**: `@allocated` after warm-up is 0 for scalar-returning wrappers, one array for
  array-returning ones (PR #263 baseline: 0 per call).
- **Leaks, PETSc side**: run the low-level tests with `-malloc_dump -objects_dump` (or `PetscMallocDump`,
  `PetscObjectsDump` via the wrappers) in a separate CI job; fail if objects survive `PetscFinalize`.
  For write-back and by-reference kinds, `PetscObjectGetReference` must return to its starting value after
  the paired Restore/Destroy.
- **Leaks, Julia side**: `test_destroy.jl` ownership/`age` checks for every emitted create/destroy pair,
  and `audit.jl` over tests and examples.

## 3. Reproduction strategy (proving fidelity)

1. Merge PR #263 into `v0.5`. Freeze the resulting `src/autowrapped/` as the **golden baseline**
   (a copy under the scratch dir, not committed).
2. Get PETSc `v3.24.0` source (shallow clone or tarball into `~/Software/PETSc/petsc-3.24.0`) and dump
   `api/petsc-3.24.0.json`.
3. Implement the core (2.1-2.2) with an empty rule set, generate into a scratch directory, diff against
   the golden baseline. Expect ~600 differing functions plus the `PoCint` names.
4. Add rules and overrides (2.3) until the diff is empty, except for a short **documented deviation list**
   (renamed `PoCint` placeholders, dropped machine path, normalised formatting of hand-edited bodies).
   Each deviation is listed in `wrapping/DEVIATIONS.md`.
5. Run the full test suite (`] test`) and the docs build on the regenerated output, including
   `test/wrapper_signatures.jl`.
6. Commit the generator, rules, snapshot and regenerated output to `v0.5`.

Acceptance for this phase: zero diff (modulo the deviation list), all tests green, docs build green,
generator runtime a few minutes rather than hours.

## 4. After fidelity: the fixes that become cheap

Each of these is a rule change plus a reviewable diff, no hand editing:

- The 236 non-Destroy functions that null a passed handle: switch the default for `Type *out` on non-Destroy
  functions to "writeback" (or "create" when the docstring marks it as output and the name says
  Create/Duplicate/Convert). Review the diff class by class.
- The 189 `VecGetLocalSize` placeholders and 96 `ni` sizes: fill `args.toml` sizes from the C docstring
  ("length n", "of size m") where the length argument is obvious; emit `error("size unknown")` rather than
  a silently wrong size for the rest.
- 39 `char**` outputs: generic rule.
- Promote `PC` to a custom handle type and re-include `PC_wrappers.jl`; the docs already reference PC
  functions.
- Remove `src/startup.jl`, `src/deprecated/`, orphaned wrapper files; bump `PETSc_jll` compat to match the
  wrapper version and make `check_wrappers_version` an error rather than a warning on a major mismatch.

## 5. Release workflow (the stated long-term goal)

For PETSc `3.X`:

```
python3 wrapping/getapi_dump.py <petsc-3.X-src> wrapping/api/petsc-3.X.json
julia --project=wrapping wrapping/apidiff.jl wrapping/api/petsc-3.24.0.json wrapping/api/petsc-3.X.json
   # read: new/removed/changed functions, stale rules, overrides whose C signature changed
julia --project=wrapping wrapping/generate.jl --petsc-dir <petsc-3.X-src>
git diff --stat src/autowrapped      # only changed/new/removed functions appear
] test
```

The wrapper file names, the prologue contract and the rule schema stay stable, so a release update is a
rule review plus a diff review. Supporting both `getAPI.py` layouts (3.24 vs 3.25+) is part of
`getapi_dump.py` from day one, since 3.25.x is the next target (JLL 3.25.4 is already in the depot).

## 6. Milestones

| # | Deliverable | Acceptance |
|---|---|---|
| M0 | PETSc 3.24.0 source, `getapi_dump.py`, `api/petsc-3.24.0.json`, PR #263 merged, golden copy | JSON loads in Julia; snapshot lists 6102+ functions |
| M1 | `generate.jl` core: type maps, kinds, template, docstring index, file map, header blocks | Generates all files; diff vs golden ~600 functions; runtime < 5 min |
| M2 | Rules + overrides encoding categories A-I and section 1.7 | Zero diff modulo `DEVIATIONS.md`; tests + docs green; idempotent; no `isa Ref ?`; ambiguities == 161 |
| M3 | `apidiff.jl`, stale-rule check, CI regenerate-and-diff job, `wrapping/README.md`, delete old generator, `local_types.jl`, PythonCall manifest | CI green on v0.5 |
| M4 | Hygiene fixes from section 4 (each its own PR) | Tests green, reviewed diffs |
| M5 | Update to PETSc 3.25.4 using the workflow in section 5 (runs right after M3, before M4) | Only changed functions in the diff; `PETSc_jll = "3.25"`; suite green on the new binaries |

## 7. Decisions recommended (say so if you disagree)

1. **Replace PythonCall with a `python3` subprocess and a JSON snapshot.** Simpler environment, and the
   snapshot is the manifest needed for release diffs.
2. **Target the `0bbf53a` template** for all functions; treat the checked-in generator's `PetscArray`
   array branch as an override for the `VecGetArray{1..4}d*` family only.
3. **Fix the `PoCint`/`RiCchardson` names** rather than reproduce them, and list them in `DEVIATIONS.md`.
4. **Normalise formatting** of the hand-edited bodies to the template (tabs, blank lines) as part of M2, so
   later diffs are purely semantic. The first commit's diff is larger but every later one smaller.
5. **Merge #263 first**, then build the golden baseline from it.
6. **Commit the API snapshot** (`wrapping/api/petsc-<ver>.json`) if it is under ~10 MB, otherwise gzip it.
