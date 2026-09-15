# How the LibPETSc wrappers are generated

Read this before touching anything under `src/autowrapped/` or `wrapping/generator/`.
It is written for the next person (or agent) who has to regenerate the low-level wrappers
for a new PETSc release, or fix a wrapper that is wrong.

## The one rule

**Never edit files in `src/autowrapped/` by hand.** They are generated, and the next
regeneration overwrites them. Every fix goes into `wrapping/generator/` (a rule, an override,
or the generator code) and is then regenerated. The history of this package shows what happens
otherwise: four PRs of hand fixes were silently lost the moment someone reran the old script.

## Layout

```
wrapping/
  WRAPPING.md            this file
  REWRITE_PLAN.md        analysis and plan of the 2026 rewrite (background, milestones)
  REGENERATING.md        notes from PR #263 on the bugs of the old generator (background)
  generatejuliabindings.jl, find_doc_strings.jl, local_types.jl, Project/Manifest.toml
                         the OLD PythonCall-based generator; obsolete, kept until deleted
  generator/             the generator (a Julia project, deps: JSON3, TOML)
    generate.jl            entry point
    getapi_dump.py         runs PETSc's getAPI.py and writes an API snapshot (JSON)
    api/petsc-X.Y.Z.json   API snapshots (one per PETSc version wrapped)
    rules/                 declarative rules (TOML), the place for fixes
      files.toml             class -> output file, exclusions, include order
      types.toml             type maps, handle types, keyword renames, string-enum overrides
      args.toml              per-function, per-argument overrides (hand-maintained)
      args_mined.toml        the same, mined once from the old hand-edited wrappers (do not edit)
    overrides/NAME.jl      verbatim replacement for one wrapper (last resort)
    prologue.jl            hand-written head of petsc_library.jl (handle structs, MPI, ...)
    petscarray.jl          hand-written PetscArray type (copied verbatim)
    structs.jl             hand-maintained struct_wrappers.jl (copied verbatim)
    petscbool.jl           hand-written PetscBool primitive type (appended to typedefs)
    src/                   generator code (see below)
    golden_diff.jl         compare two autowrapped directories function by function
    bootstrap_rules.jl     mine rules from a hand-edited autowrapped directory (one-off)
    make_overrides.jl      copy functions from an autowrapped directory into overrides/
```

## Pipeline

```
PETSc source tree ──getapi_dump.py──▶ api/petsc-X.Y.Z.json       (functions, args, enums, structs)
PETSc source tree ──docs.jl─────────▶ manual-page index          (one pass over src/**/*.c)
snapshot + index + rules + overrides ──generate.jl──▶ src/autowrapped/*.jl
```

Run it:

```sh
python3 wrapping/generator/getapi_dump.py /path/to/petsc-X.Y.Z wrapping/generator/api/petsc-X.Y.Z.json
julia --project=wrapping/generator wrapping/generator/generate.jl \
      --petsc-dir /path/to/petsc-X.Y.Z --api wrapping/generator/api/petsc-X.Y.Z.json
```

Only the *source* tree is needed (no `configure`, no build): the 16 MB release tarball from
`https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-X.Y.Z.tar.gz`.
A full run takes a few seconds. `--out DIR` writes somewhere else than `src/autowrapped`.

`getAPI.py` lives in `config/utils/` up to PETSc 3.24 and in `lib/petsc/bin/` from 3.25; it
needs the current directory to be the PETSc tree in both cases. `getapi_dump.py` handles both.
The snapshot is deterministic (sorted keys), so two dumps of one tree are identical and two
releases can be diffed.

## Generator code (`generator/src/`)

| file | does |
|---|---|
| `api.jl` | loads the JSON snapshot into `API` (functions, classes, enums, senums, typedefs, structs) |
| `docs.jl` | indexes the `/*@ ... @*/` manual pages, extracts Input/Output parameter lists, cleans the text into a Julia docstring, finds the manual section (`DMPlex`, `PC`, ...) from the `makefile` next to the C file |
| `types.jl` | `Rules` (from the TOML files), C -> Julia type mapping, `abstract_arg_type` |
| `classify.jl` | decides for every argument what kind it is and emits its code fragments |
| `render.jl` | writes docstring + untyped stub + `@for_petsc` method for one function |
| `support.jl` | enums, string enums, typedefs, opaque type declarations, version file, `petsc_library.jl` |
| `driver.jl` | file planning, overrides, and the `generate()` entry point |
| `blocks.jl` | splits autowrapped files into named blocks and categorises differences (used by the tools) |

### What decides input vs output

1. The manual page. If the argument is listed under `Input Parameters` it is an input; if
   under `Output Parameters` (and it is a pointer) it is an output. This is the primary rule.
2. Otherwise a single pointer (`T *x`) is guessed to be an output, unless the function is not a
   `Create`/`Duplicate`/`*Type*` function and `T` is not a scalar/enum/string type, or the
   function name contains `Restore` or `Copy` (the old generator's heuristics).
3. `[FunctionName.arg] direction = "in" | "out"` in `args.toml` overrides both.

### Argument kinds (how each C shape is wrapped)

| C shape | Julia signature | ccall | body |
|---|---|---|---|
| `T x` scalar | `x::T` | `T` | |
| `T *x` output, scalar/enum | return `x::T` | `Ptr{T}` | `x_ = Ref{T}()` ... `x = x_[]` |
| `XType *x` output (string enum) | return `x::String` | `Ptr{XType}` | `unsafe_string`, `""` when NULL |
| `Vec x` (handle) | `x::AbstractPetscVec` | `CVec` | via `unsafe_convert` |
| `Vec *x` output | return `x::PetscVec` | `Ptr{CVec}` | `x_ = Ref{CVec}()` ... `x = PetscVec(x_[], petsclib)` |
| `Vec *x` in `XDestroy` | `x::AbstractPetscVec` | `Ptr{CVec}` | `x_ = Ref(x.ptr)` ... `x.ptr = C_NULL` |
| `Vec *x` not documented as output, not Destroy | `x::AbstractPetscVec` | `Ptr{CVec}` | `x_ = Ref(x.ptr)` ... `x.ptr = x_[]` (write-back) |
| `Opaque *x` in `XDestroy` | `x::Union{X, Ref{X}}` | `Ptr{X}` | `x_ = x isa Base.RefValue ? x : Ref{X}(x)` |
| `const T *x` / `T x[]` input | `x::Vector{T}` | `Ptr{T}` | |
| `const char *s` / `char s[]` input | `s::String` | `Ptr{Cchar}` | |
| `char buf[], size_t len` | `buf::Vector{Cchar}, len::Csize_t` | `Ptr{Cchar}, Csize_t` | caller allocates |
| `char **s` output | return `s::String` | `Ptr{Ptr{Cchar}}` | `unsafe_string` |
| `T *x[]` output (PETSc-owned array) | return `x::Vector{T}` if a `size` rule exists, else `x::Ptr{T}` | `Ptr{Ptr{T}}` | `unsafe_wrap(Array, x_[], size; own = false)` or the raw pointer |
| `T x[]` output (caller-allocated) | return `x::Vector{T}` if a `len` rule exists or an argument `ni` exists, else input `x::Vector{T}` | `Ptr{T}` | `x = Vector{T}(undef, len)` |
| `T *x[]` input (Restore*) | `x::Vector{T}` | `Ptr{Ptr{T}}` | `x_ = Ref(pointer(x))` |
| `T *n` input in Restore* | `n::T` | `Ptr{T}` | `n_ = Ref{T}(n)` |
| `T **x` output, not an array | return `x::Ptr{T}` | `Ptr{Ptr{T}}` | raw pointer |
| function pointer (`XFn *f`, `void (*f)(...)`) | `f::Ptr{Cvoid}` | `Ptr{Cvoid}` | pass an `@cfunction` pointer |
| `XFn **f` output | return `f::Ptr{Cvoid}` | `Ptr{Ptr{Cvoid}}` | |
| `void *ctx` | `ctx::Ptr{Cvoid}` | `Ptr{Cvoid}` | |
| `void **ctx` output | return `ctx::Ptr{Cvoid}` | `Ptr{Ptr{Cvoid}}` | |

Input handle arguments always take the abstract type (`AbstractPetscVec`), so `VecPtr`, `MatShell`
and the typed DM hierarchy pass; return positions use the concrete type (`PetscVec`).
`test/wrapper_signatures.jl` enforces this.

### Rules (`generator/rules/`)

`files.toml`: one `[[file]]` per output file with `classes = [...]` (or `standalone = true` for
functions not attached to a class), optional `include = false` to generate but not include.
`[exclude]` lists functions never wrapped, with a reason. Functions containing `_` are skipped.

`types.toml`: the C -> Julia replacements (`fix_substring_replacements` switches from the old
substring replacement, which turned `PetscPointFn` into `PetscPoCintFn`, to word-boundary
replacement), `[[handles]]` (C name, Julia struct, abstract type, C alias), `[rename_args]`
(Julia keywords used as C argument names), `[senum_overrides]` (`VecType = "Cstring"`), and
`[predeclared]` names the generator must not declare as opaque types.

`args.toml` (hand-maintained) and `args_mined.toml` (written by `bootstrap_rules.jl`, never edited;
`args.toml` wins on conflicts): per argument, keyed `[FunctionName.argname]`:

| key | meaning |
|---|---|
| `direction = "in"/"out"` | force the classification |
| `nullable = true` | input also accepts a `Ptr` (so `C_NULL` can be passed) |
| `byref = true` | opaque handle passed as `Union{X, Ref{X}}` |
| `size = "expr"` | length or dims to `unsafe_wrap` a PETSc-owned output array; may use other arguments/outputs |
| `prelude = "code"` | line(s) emitted before the wrap, e.g. `m, n = MatGetLocalSize(petsclib, A)` |
| `nullinit = true` | initialise the output `Ref` with `C_NULL` |
| `len = "expr"` | length of a caller-allocated output `Vector` |

Rules that name a function or argument absent from the snapshot are reported when generating
(TODO: `apidiff.jl`, see REWRITE_PLAN.md).

### Overrides (`generator/overrides/NAME.jl`)

A file `NAME.jl` replaces the generated block for `NAME` verbatim (docstring, stub and
`@for_petsc` method). Its first line records the C signature it was written against, so a
changed signature can be flagged. Files whose name is not a function in the snapshot are
collected into `extra_wrappers.jl` (macros such as `PETSC_VIEWER_STDOUT_WORLD`, hand-written
helpers such as `DMProjectFunction`). Use an override only when no rule can express the wrapper
(multi-dimensional `PetscArray` views, `MatGetRowIJ`, `PetscOptionsGetString`, ...).

## Checking a regeneration

`golden_diff.jl` compares two autowrapped directories block by block, independent of the order of
functions in a file:

```sh
julia --project=wrapping/generator wrapping/generator/golden_diff.jl OLD_DIR NEW_DIR --categorize
julia --project=wrapping/generator wrapping/generator/golden_diff.jl OLD_DIR NEW_DIR --show VecGetArray
```

`--categorize` buckets every differing wrapper (docs-only, voidptr-fix, writeback-fix, ...) and
writes the lists to `cat_*.txt` in `$GOLDEN_DIFF_OUT` (default: temp dir). The categories that
represent deliberate differences from the old hand-edited files are listed in `DEVIATIONS.md`.

After regenerating, always run

```sh
grep -rn "isa Ref ?" src/autowrapped/                    # must be empty (Ptr{T} <: Ref{T})
grep -rn "VecGetLocalSize(petsclib, x)" src/autowrapped/ # placeholder sizes: must be empty
julia --project=. -e 'using Pkg; Pkg.test()'             # includes test/wrapper_signatures.jl
```

## Things that bite

- The manual page of a PETSc function sometimes ends with `*/` instead of `@*/`; the indexer
  stops a block at `*/` or at the next `/*@` so the following page is not swallowed.
- The first docstring line is the text after the first `-` of `Name - description`; a
  description containing `-` is truncated (`Normalizes a vector by its 2` for `2-norm`). This
  reproduces the old output; fix it in `docs.jl` (`_finish_block`) if you want.
- `getAPI.py` gives `void (*f0)(...)` arguments the type name `void(f0` and an empty name; the
  classifier turns these into `f0::Ptr{Cvoid}`.
- C++ scoped names (`moab::Range`, `std::size_t`) become `moab_Range`, `Csize_t`.
- Argument names that are Julia keywords (`begin`, `end`, `function`, `global`, `local`, ...)
  get a trailing underscore or a short alias (see `[rename_args]`).
- `@for_petsc` bodies get `$PetscScalar`, `$PetscReal`, `$PetscInt`, `$PetscComplex` by plain
  substring replacement (`dispatch_types`), exactly like the old generator.
- `struct_wrappers.jl` is hand-maintained (`generator/structs.jl`): field order must match the
  C struct. Check `api/petsc-X.Y.Z.json` (`structs`) when moving to a new release.
- The prologue (`generator/prologue.jl`) is the only copy of the handle structs. The old
  `wrapping/prologue.jl` and `local_types.jl` are gone; do not resurrect them.
