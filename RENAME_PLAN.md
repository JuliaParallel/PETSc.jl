# Plan: high-level API rename (docs/src/man/naming.md) for v0.5

Status: planning, 2026-09-16. Inventory of the code against `naming.md` done (three read-only
sweeps over `src/*.jl`, tests, examples and docs). Nothing implemented yet beyond §5.2/§5.5 (the
`KSP`/`SNES` type names and constructors) and the `AbstractPetscDM` widening of the generated layer.

## What the inventory found

- About 130 high-level bindings across `dm.jl`, `dmda.jl`, `dmstag.jl`, `dmplex.jl`, `vec.jl`,
  `mat.jl`, `ksp.jl`, `snes.jl`, `ts.jl`, `init.jl`, `options.jl`, `sys.jl`, `audit.jl`;
  roughly 100 of them change name, signature or both.
- Heaviest call-site moves: `destroy` -> `destroy!` (351 sites in tests and examples),
  `setfield!` -> `set_field!` (52), `DMGlobalVec`/`DMLocalVec` (79), `DMDA`/`DMStag`/`DMPlex`
  constructors (236, unchanged names but new types), `withlocalarray!` (63), the `@petsc_*_fn`
  macros (114), `dm_project_*`/`add_boundary!`/`vtk_save!` (drop `petsclib`, reorder; ~80).
- `ts.jl` already follows the document except: `type` (should be `type_name`), a hand
  `destroy` alias, and subject-first duplicates of all five callback setters.
- Typed exceptions (§14) are done except the nine DM-flavour `@assert`s, which the typed DM
  hierarchy deletes. `PetscNotInitialized` and `PetscError` exist.
- None of the migration tooling exists: no `scripts/`, `src/deprecations.jl`, `@renamed`,
  `test/test_deprecations.jl`, `_doc_borrowed`, `api_surface.jl`; CI does not run
  `--depwarn=yes`.
- `PetscDM` is one concrete type with a runtime string flavour; `DMDA`/`DMStag`/`DMPlex` are
  factory functions. Thirteen sites return 3-padded tuples; five examples and `test/dmda.jl`
  index the padding.
- Latent bugs to fix on the way: `DMStag(dm, dof; options...)` slurps `options` positionally;
  `getghostcorners_dmstag` documents a `nextra` field it does not return; `getlocalcoordinatearray`
  computes unused ghost corners; `ownershiprange` indexes plain integers with `[]`;
  `Base.size(dm)` has a wrong error message; `HostBackend` is exported but never defined.

## Decisions the document leaves open (settled by Boris, 2026-09-16)

Decided: 1) export the nine types plus `petsclibs`, keep `export LibPETSc`; 2) `public`
unconditionally; 3) the names missing from the table are kept (renamed as proposed below);
4) `set_type!(obj, ::Symbol)` is added for Vec, Mat, KSP, SNES; 5) `MatAT` is deleted.
The original questions follow for reference.


1. Export list: §13 includes `petsclibs`, §13.1 does not. Whether `LibPETSc` stays exported is
   not stated. Proposal: export the nine types plus `petsclibs` and keep `export LibPETSc`.
2. §13.1's `@static if VERSION >= v"1.11"` gate is dead: `Project.toml` requires Julia 1.12.
   Proposal: `public` unconditionally.
3. Names missing from the rename table (would fail the register check): `owns`,
   `tao_usable_after_reinitialize`, `_taoterm_resettable`, `_reset_stale_register_flags`,
   `setconvergencetest!` (+ `Fn_SNESSetConvergenceTest`, `SNESConvergenceTestBox`), `setup!`,
   `set_dm!`/`comm`/`set_from_options!` on TS, `HostBackend`, `_MATSEQAIJ_WITHARRAYS_STORAGE`,
   `_PETSC_ERR_LIB`, `KSP(petsclib, comm, S::SparseMatrixCSC)`. Proposal: `setconvergencetest!`
   -> `set_convergence_test!` (callback-first only), `Fn_SNESSetConvergenceTest` ->
   `SNESSetConvergenceTestFn`; the rest marked internal or unchanged.
4. `set_type!(ksp, :gmres)` / `set_type!(snes, ...)` / `set_type!(v|A, ...)` do not exist today
   (only the TS one); they are additions, not renames. `type_name` must replace two spellings
   (`type` on Vec/Mat/KSP/TS, `gettype` on SNES/DM).
5. `info(dm)` field `mpi_proc_size` -> `procs` is implied by §12 but not in the table.
6. §16 says nine exports are dropped; the real list is twelve (also the three MemBackend types).
7. `MatAT` has no uses; rename to `MatOrTranspose` or delete.

## Steps (CI green after each; one commit per step on `v0.5`)

### Step 1 — register and shims, pure renames

- `scripts/renames.jl`: the register (`old => new` pairs, internal set, public set, creator and
  destroyer sets), generating `src/deprecations.jl` (`@renamed old new`, `@warn ... maxlog = 1`),
  `test/test_deprecations.jl`, the `public` declaration, and the `audit.jl` name sets.
- `scripts/api_surface.jl --check | --sweeps` per §1.1.
- Rename every name in the table whose meaning is unchanged (`destroy!`, `corners`, `info`,
  `set_field!`, `with_local_array!`, `ghost_update!`, `set_values!`, `isinitialized`, `comm`,
  `parse_option`, `audit_file`, `check_wrappers_version`, macros, `Fn_*` -> `*Fn`,
  `MatOrTranspose`, internals losing `_`). `finalizer(destroy!, ...)`.
- Convert tests, examples and docs to the new names. CI runs with `--depwarn=yes`.
- Merge duplicates: `DMGlobalVec`/`dm_create_global_vec` -> `global_vec`, local likewise,
  `mat_null_space_destroy!` -> `destroy!`, `vtk_save!`/`vtk_save_fields!` -> `save_vtk!`.

### Step 2 — constructors and types

- `PetscVec(...)`, `PetscMat(...)`, `PetscOptions(...)`, `PetscLibType(path; ...)` constructors
  replacing `VecSeq`, `MatSeqAIJ`, `MatSeqDense`, `MatCreateSeqAIJ`, `MatSeqAIJWithArrays`,
  `as_petsc_vec`, `Options`, `set_petsclib` (shims forward).
- Typed DM hierarchy: `DMDA{L,N}`, `DMStag{L,N}`, `DMPlex{L} <: AbstractPetscDM{L}` with `own`
  field; `narrow`; `dm(ksp)`/`dm(snes)`/`dm(ts)` call it; the nine flavour asserts and the three
  string dispatches go. Generated wrappers already take `AbstractPetscDM`.
- `own` flag on every wrapper, `destroy!` consults it, readers return `own = false`,
  `_doc_borrowed` helper and its CI grep.

### Step 3 — semantic breaks without shims (§16)

- Dimension-correct `corners`, `ghost_corners`, `local_indices`, `global_indices`, `info`,
  `Base.size` built with `ntuple(..., Val(N))`; consumers in `dmda.jl`, five examples and
  `test/dmda.jl` adapted. `info` drops `s`, `dof` -> `ndofs`, `mpi_proc_size` -> `procs`.
- `type_name` returns `Symbol`; new `set_type!(obj, ::Symbol)` for Vec/Mat/KSP/SNES; the
  `Base.show` string comparison in `mat.jl` and the string-comparing tests adapted.
- Argument order: `project_function!(X, dm, ...)`, `project_field!(X, dm, ...)`,
  `global_to_local!(lvec, dm, gvec)`, `local_to_global!(gvec, dm, lvec)`; drop `petsclib` from
  the eight §8 functions; forwarding shims for the four reorders.
- Callback-first only: delete the subject-first methods in `snes.jl`, `ksp.jl`, `ts.jl`
  (12 call sites in tests and examples, 4 in docs move).
- `ownership_range(A)` 1-based only (positional `false` warns); `set_values!` parameter names
  `rows_0b`/`cols_0b`; `star_fd_coloring` fields get `_0b`/`_1b`.

### Step 4 — exports, docs, tooling

- Export list per decision 1; `public` for the register; remove the twelve old exports.
- `library_info` returns `(; source, path, scalar, int, real)` with a `show` method;
  `error(...)` on user input in `init.jl`/`options.jl` -> `ArgumentError`.
- Docs: `dmplex.md`, `dmstag.md`, `mat.md`, `vec.md`, `snes.md`, `getting_started.md`,
  `installation.md`, `FAQ.md`, `hpc.md`, README converted; `naming.md` corrected where the
  inventory disagrees (decisions 2, 3, 5, 6, §14 assert count, `src/deprecated/` reference);
  the C-name -> Julia-name index page from `_doc_external` entries (§15).
- `CHANGELOG.md` section for the rename; `audit_file` name sets generated from the register.

## How the work is run

Steps are implemented by subagents (Opus) in an isolated worktree, one step at a time, each
followed by the full test suite locally and CI after the push. Step 1 is the largest and
mostly mechanical; steps 2 and 3 carry the design risk and get a review pass before merging.
