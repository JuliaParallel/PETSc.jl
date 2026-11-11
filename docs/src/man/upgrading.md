# Upgrading from 0.4 and 0.5

PETSc.jl 0.5 is a breaking release. Most of the high-level API has new names, the low-level `LibPETSc` functions return their outputs, and the package needs Julia 1.12 and PETSc 3.25. This page is the order in which to update code written for 0.4. 0.6 removes the 0.5 warnings that name each replacement, so upgrade to 0.5.5 first, follow this page, and then move to 0.6 as described in [Upgrading from 0.5 to 0.6](@ref). Every change is listed in the [Release notes](release_notes.md), and every renamed function in the rename [appendix](#Appendix:-the-0.4-to-0.5-renames) below.

## 1. Requirements

- **Julia 1.12 or newer.** On an older Julia, the package manager keeps installing 0.4.
- **PETSc 3.25.** `PETSc_jll` 3.25 is installed automatically. A library configured with `PETSc.set_library!` must be a PETSc 3.25 build; `PETSc.check_wrappers_version()` warns when it is not.
- **One `PetscInt` width per process.** `PETSc.petsclibs` holds four libraries (`Float64`, `Float32`, `ComplexF64`, `ComplexF32`) with `Int64` indices. Code that indexes `petsclibs[5:8]`, or asks `PETSc.getlib` for an `Int32` library, fails. To use `Int32` indices, call `PETSc.set_petscint!(Int32)` once and restart Julia (see [Installation](installation.md)).

## 2. Run your code on 0.5.5 and fix the warnings

About a hundred functions were renamed. On 0.5.5 the old names still work and print a warning the first time each call site runs, for example

```
┌ Warning: getcorners is deprecated, use corners
```

These warnings appear without `--depwarn=yes`. Pin `PETSc = "0.5.5"` in your `[compat]`, run your code or its test suite once, and replace every name a warning reports. On 0.6 an old name is an `UndefVarError` and an old argument order a `MethodError`, with no hint at the replacement.

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

**Callback setters take the callback first.** 0.4 accepted both orders; 0.5 accepts only the callback-first one, which is the one `do` syntax needs. This applies to `set_function!`, `set_snes_jacobian!`, `set_convergence_test!`, `set_compute_rhs!`, `set_compute_operators!`, `set_rhs_function!`, `set_rhs_jacobian!`, `set_ifunction!`, `set_ijacobian!` and `set_monitor!`. `add_coarsen_hook!` takes C function pointers, which `do` can't supply, so it takes the DM first: `add_coarsen_hook!(dm, hook)`.

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

On 0.5.5 the 0.4 names still work and warn, but only when called directly: a 0.4-order call passed through `invoke`, or through a variable holding the function, is not translated. The same applies to `project_function!` and `project_field!`.

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
- **Renamed handle types.** `PetscKSP` and `PetscSNES` are now `KSP` and `SNES` (`AbstractPetscKSP` and `AbstractPetscSNES` are `AbstractKSP` and `AbstractSNES`). The old names work up to 0.5.5 and are gone in 0.6.
- **Functions PETSc defines only in its headers are gone.** About 95 macros and inline functions with no symbol in the library (`VecSetValue`, `MatSetValue`, `PetscStrcmp`, `PetscTime`, …) are no longer wrapped. Use the array versions (`VecSetValues`, `MatSetValues`) or plain Julia.
- **`PetscBool` is one byte**, as it is in PETSc since 3.24. Code that uses `LibPETSc.PetscBool` needs no change. Code that writes a PETSc boolean as `Cint` or `Int32` in its own `ccall` or `unsafe_store!` must use `LibPETSc.PetscBool` instead.

- **`PetscFE`, `DMLabel` and `MatNullSpace` are handle types** (0.6), as `PC` became one in 0.5: `PetscFE{PetscLib}` and the others carry `ptr`, `age` and `own`, and the wrappers take them where they took a `Ptr`. Code that passes a raw `Ptr` to these wrappers, or uses what `LibPETSc.DMGetLabel` and the `Create` functions return as a `Ptr`, gets a `MethodError`. To wrap a pointer, write `LibPETSc.PetscFE(ptr, petsclib; own)`.

The [Release notes](release_notes.md) list the remaining changes, and [`wrapping/DEVIATIONS.md`](https://github.com/JuliaParallel/PETSc.jl/blob/main/wrapping/DEVIATIONS.md) lists every category of changed wrapper.

## 6. Within 0.5: the 0.5.4 spellings

0.5.4 moves options to keywords and puts arguments in the same order across constructors. Every old form works up to 0.5.5 and warns, so step 2 finds them, and is gone in 0.6; this is the list to update by hand.

| Up to 0.5.3 | From 0.5.4 |
|---|---|
| `set_values!(A, rows, cols, vals, ADD_VALUES)` | `set_values!(A, rows, cols, vals; mode = ADD_VALUES)` |
| `ghost_update!(v, ADD_VALUES, SCATTER_REVERSE)` | `ghost_update!(v; mode = ADD_VALUES, scatter = SCATTER_REVERSE)` |
| `global_to_local!(l, dm, g, ADD_VALUES)`, also `local_to_global!`, `local_to_local!`, `project_function!`, `project_field!` | `…; mode = ADD_VALUES` |
| `zero_rows!(A, rows, 1.0)` | `zero_rows!(A, rows; diag = 1.0)` |
| `DMStag(petsclib, comm, bnd, dims, dof, 1, :box)`, the same for `DMDA` | `DMStag(petsclib, comm, bnd, dims, dof; stencil_width = 1, stencil_type = :box)`; both default to these values |
| `setfromoptions`, `dmsetup`, `dmsetfromoptions` on a DM constructor | `set_from_options`, `setup` |
| `add_coarsen_hook!(hook, dm)` | `add_coarsen_hook!(dm, hook)` |
| `distribute!(dm)` | `distribute(dm)`: it returns a new DM and leaves `dm` as it is |
| `VecPtr(petsclib, ptr, false)`, also `MatPtr` and `PetscVec(petsclib, ptr, own)` | `PetscVec(petsclib, ptr; own = false)`, `PetscMat(petsclib, ptr; own = false)`; `own` defaults to `false` |
| `PetscMat(petsclib, rowptr, colval, nzval; comm)` | `PetscMat(petsclib, comm, rowptr, colval, nzval)` |
| `MatShell(petsclib, f, comm, m, n, M, N)` | `MatShell(f, petsclib, comm, m, n; global_rows = M, global_cols = N)`, or `MatShell(petsclib, comm, m, n) do y, x … end` |

0.5.4 also adds `PetscVec(petsclib, comm, n; local_size)` for a vector distributed over `comm`.

## 7. Within 0.5: the 0.5.5 deprecations

0.5.5 is the last release before 0.6. Everything below works on 0.5.5 and warns once, so step 2 finds it, and is removed in 0.6; this is the list to update by hand.

| Up to 0.5.4 | From 0.5.5 |
|---|---|
| `VecPtr(petsclib, ptr; own)`, `MatPtr(petsclib, ptr; own)` | `PetscVec(petsclib, ptr; own)`, `PetscMat(petsclib, ptr; own)` (new) |
| A SNES or TS callback that takes a trailing context argument; `user_ctx(snes)`, `user_ctx(ts)`, `set_user_ctx!` | A closure carries the context: `set_rhs_function!(ts) do F, ts, t, u … end` reads it from the enclosing scope |
| `local_arrays(…)` with `restore_local_arrays!` | `with_local_arrays!(f!, g_fx, l_x)` |
| `local_coordinate_array(dm)` | `with_local_coordinate_array(f, dm)`; the array is handed back when `f` returns or throws |
| `@residual_fn f0 Nf`, `@jacobian_fn`, `@bd_fn`, `@simple_fn` | `@residual_fn petsclib f0 Nf`: the library comes first, as a top-level `const` |
| `vtk_merge_tensor`, `create_split_boundary_labels!` | Copies in `examples/dmplex/vtk_tensor.jl` and `examples/dmplex/ex69.jl` |
| `isapprox(v, w)` on two PetscVecs | Compares entry by entry up to 0.5.5 and norm-wise from 0.6, ‖v − w‖ ≤ max(atol, rtol·max(‖v‖, ‖w‖)); say which you mean with `all(isapprox.(v[:], w[:]))` or `norm(v - w) ≤ …` |

0.5.5 also adds `x + y`, `x - y`, `α * x`, `dot`, `axpy!` and `copy` on PetscVecs and PetscMats, `MatShell` operations (`mult_transpose`, `diagonal`, `mult_add`) and AMD GPU support.

### What changes in 0.6 without a warning

These can't warn without firing on correct code, so check them now: see [Upgrading from 0.5 to 0.6](@ref) below. On 0.5.5, `options_at = :first_solve` and `check = true` on `solve!` give the 0.6 defaults early.

## 8. Looking up a name

- **An old high-level name:** the [Appendix: the 0.4 to 0.5 renames](@ref) below.
- **A PETSc C function:** the [C to Julia name index](api_index.md).
- **Why a name is what it is:** the [API conventions](conventions.md).

## Upgrading from 0.5 to 0.6

0.6 removes every 0.5 deprecation and changes some behavior without a warning. The changes below give no error, or an error that does not name the fix. The [Release notes](release_notes.md) list all of them.

- **Removed deprecations.** Every name, argument order and keyword that warned on 0.5.5 is an `UndefVarError` or `MethodError` on 0.6. Run your code on 0.5.5 first and fix every warning (step 2).
- **Pointwise functions get `PtrVector`s.** The array arguments of a function made with `@simple_fn`, `@residual_fn`, `@jacobian_fn` or `@bd_fn` are `PETSc.PtrVector`s, not `Vector`s. An argument annotated `::Vector`, `::DenseVector` or `::StridedVector`, or a `reshape` annotated `::Matrix`, throws a `MethodError`. Annotate `::AbstractVector` (or `::AbstractMatrix`), or leave the annotation out.
- **A callback's return value is ignored**, a nonzero `Integer` too. Throw an exception to report a failure.
- **`user_ctx` is removed**, with `set_user_ctx!` and the SNES and TS callbacks that take a trailing context. A callback with only the trailing-context method throws a `MethodError` when called. Let a closure carry the context.
- **`solve!` and `step!` throw a `ConvergenceError`** when the solver stops without converging. Pass `check = false` where the code reads `converged_reason` after a failed solve.
- **Options are applied at the first `setup!` or `solve!`** of a KSP, SNES or TS, not in the constructor or at every solve. Until then a KSP has no type: code that needs it earlier, such as `set_fieldsplit_is!` on a PC chosen by `pc_type`, calls `set_type!` itself or passes `options_at = :construction`. Pass `options_at = :each_solve` to re-read them at every solve.
- **The options database wins over constructor keywords.** KSP, SNES, TS and the DM constructors read `PETSC_OPTIONS` and the options given to `PETSc.initialize` together with their keywords, and the database wins on a key both set. Remove the key from the database where the keyword should apply.
- **`w = 2 .* x` returns a `PetscVec`**, not a `Vector`. On more than one rank, destroy it with `destroy!`. Code that used the `Vector` (`2 .* x == y_julia`, `findall(x .> 0)`) indexes with `x[:]` first.
- **`size(v)` is the tuple `(n,)`.** `size(v) == n` is `false` without an error; write `length(v) == n`.
- **`v == w` compares the entries** of two PetscVecs, collectively, and `hash` throws on a PetscVec or PetscMat. Compare handles with `===`, and key a `Dict` of handles with an `IdDict`.
- **`A == B` compares the entries** of two PetscMats, so a stored zero no longer makes them unequal. `PETSc.equal(A, B)` is the old `MatEqual` test.
- **`copy_disc!(dst, src)` takes the destination first.** A call in the old order runs and copies the wrong way: swap the arguments.
- **A matrix of values in `set_values!(A, rows_0b, cols_0b, vals)` is the block**, as in `A[rows, cols] = vals`. On 0.5 a non-symmetric block landed transposed. Code that passed `permutedims(B)` now passes `B`.
- **`LibPETSc.MatStencil` fields are `k, j, i, c`.** Positional construction is unchanged, but `s.i` now reads the third field. Write `s.k` where the code read `s.i`, and the reverse.
- **`A * B` of two `PetscMat`s is a `PetscMat`**, not a `Matrix`. Call `collect(A * B)` for a Julia matrix, and `destroy!` the product on more than one rank.
- **`unsafe_local_array` and `wrap_local_array` are removed.** Use `with_local_array!(f!, v)`, or `acquire_local_array` with `release_local_array`.
- **`set_convergence_test!` passes `PetscReal` norms**, not `Float64`. This matters only for a Float32 library with a callback whose arguments are annotated `Float64`.
- **`distribute(dm)` on one rank returns a clone of `dm`** with its fields, not a DM with a null pointer, so `d = distribute(dm); destroy!(dm); dm = d` needs no special case. A `d.ptr == C_NULL` check is now always false.
- **FE, label and null-space functions lose `petsclib`.** The handles carry the library, so `set_name!(obj, name)`, `copy_quadrature!(dst_fe, src_fe)` (destination first, it was `(petsclib, src, dst)`), `compose_constant_nullspace!(fe)` and `destroy!(nullspace)` take no `petsclib` or communicator. `label(dm, name)` returns a `DMLabel`, or `nothing` where it returned `C_NULL`: `c_label != C_NULL` is now always true without an error, so test `isnothing(c_label)`. `set_field!`, `add_boundary!` and `add_natural_boundary!` take a `PetscFE` and a `DMLabel` where they took `Ptr`s; `set_field!(dm, i, fe; label = C_NULL)` is `label = nothing` (the default).
- **`isapprox(v, w)` compares norms**, ‖v − w‖ ≤ max(atol, rtol·max(‖v‖, ‖w‖)), as for Julia arrays, instead of entry by entry. For the old test write `all(isapprox.(v[:], w[:]))`.
- **`PetscVec` is an `AbstractVector` and `PetscMat` an `AbstractMatrix`.** Your methods for `AbstractVector` or `AbstractMatrix` now receive them. If such a method reads entry by entry, each read is a C call, and an entry owned by another rank is an error: work on the local part with `with_local_array!`, or add a method for `PetscVec`.
- **`PetscInt`, `PetscScalar` and `PetscReal` are properties of a library, not of its type.** `typeof(petsclib).PetscInt` throws an error; use `PETSc.inttype(T)`, `PETSc.scalartype(T)` and `PETSc.realtype(T)`, which take a library or its type.
- **The positional constructors that took callbacks are removed**, such as `SNES{PetscLib}(ptr, age, f!, …)` and `KSP{PetscLib}(ptr, age, computerhs!, …)`. Wrap the pointer with `SNES{PetscLib}(ptr, age)` and set the callbacks with `set_function!`, `set_compute_rhs!` and the other setters.
- **`audit_file` is removed** from the package. The leak audit is the script `scripts/audit.jl` in the repository: `julia scripts/audit.jl file.jl`.

## Appendix: the 0.4 to 0.5 renames

What each 0.4 high-level name became in 0.5, by the file that defined it. On 0.5.5 every old name still works and warns with its replacement; on 0.6 it is gone.

### `dm.jl`

| v0.4 | v0.5 |
|---|---|
| `destroy` | `destroy!` |
| `getinfo` | `info` (drops the duplicate `s` field, `dof` becomes `ndofs` and `mpi_proc_size` becomes `procs`) |
| `getcorners`, `getcorners_dmda` | `corners` |
| `getghostcorners`, `getghostcorners_dmda` | `ghost_corners` |
| `dm_local_to_global`, `dm_local_to_global!` | `local_to_global`, `local_to_global!` |
| `dm_global_to_local`, `dm_global_to_local!` | `global_to_local`, `global_to_local!` |
| `setuniformcoordinates_dmda!` | `set_uniform_coordinates!` |
| `coordinatesDMLocalVec` | `local_coordinates` |
| `getlocalcoordinatearray` | `local_coordinate_array`; removed in 0.6, use `with_local_coordinate_array(f, dm)` |
| `MatAIJ` | `PetscMat` constructor |
| `DMGlobalVec` | `global_vec` (merged, see below) |
| `DMLocalVec` | `local_vec` (merged, see below) |
| `getdimension` | `Base.ndims` |
| `setfromoptions!` | `set_from_options!` |

`DMGlobalVec` in `dm.jl` and `dm_create_global_vec` in `dmplex.jl` both call `DMCreateGlobalVector`; the local pair duplicates likewise.
Each pair collapses to one name, so two shims point at each replacement.

### `dmstag.jl`

| v0.4 | v0.5 |
|---|---|
| `getcorners_dmstag` | `corners` (method on `DMStag`) |
| `getghostcorners_dmstag` | `ghost_corners` (method on `DMStag`; it has no `nextra` field, which `DMStagGetGhostCorners` never reported although the v0.4 docstring promised it) |
| `local_indices_dmstag` | `local_indices` |
| `global_indices_dmstag` | `global_indices` |
| `setuniformcoordinates_stag!` | `set_uniform_coordinates!` |
| `DMStagDOF_Slot` | `dof_slot` |

### `dmda.jl`

| v0.4 | v0.5 |
|---|---|
| `reshapelocalarray` | `reshape_local_array` |
| `localinteriorlinearindex` | `local_interior_linear_index` |
| `dmda_star_fd_coloring` | `star_fd_coloring(da::DMDA{L,2})` (drops `petsclib`, loses its export, and every index field gains a `_0b` or `_1b` suffix) |
| `ndofs` | unchanged (closed list) |

### `dmplex.jl`

| v0.4 | v0.5 |
|---|---|
| `isplexsimplex` | `issimplex` |
| `plexdistribute!` | `distribute` |
| `petsc_setname!` | `set_name!` |
| `getds` | `ds` |
| `createds!` | `create_ds!` |
| `getlabel` | `label` |
| `dm_project_function!` | `project_function!` |
| `dm_project_field!` | `project_field!` |
| `dm_compute_l2diff` | `l2diff` |
| `dm_create_global_vec` | `global_vec` |
| `dm_create_local_vec` | `local_vec` |
| `dm_set_auxiliary_vec!` | `set_auxiliary_vec!` |
| `dm_coarsen_hook_add!` | `add_coarsen_hook!` |
| `dm_copy_disc!` | `copy_disc!` |
| `dm_get_coarse` | `coarse_dm` |
| `fe_copy_quadrature!` | `copy_quadrature!` |
| `mat_null_space_create` | `mat_nullspace_create` (free function, prefix kept) |
| `mat_set_null_space!` | `set_nullspace!` |
| `mat_null_space_destroy!` | `destroy!` |
| `fe_create_default`, `fe_create_lagrange` | unchanged (free functions) |
| `vtk_save!`, `vtk_save_fields!` | `save_vtk!` (drops `petsclib`); `save_vtk` from 0.5.1 |
| `vtk_merge_tensor!` | unchanged (free function, no PETSc object); `vtk_merge_tensor` from 0.5.1; deprecated in 0.5.5, removed in 0.6 (a copy is in examples/dmplex/vtk_tensor.jl) |
| `setfield!` | `set_field!` (resolves the `Base.setfield!` shadow) |
| `dmclone` | `clone` |
| `plex_set_snes_local_fem!` | `set_snes_local_fem!` |
| `snes_set_jacobian_null_space!` | `set_jacobian_nullspace!` |
| `fe_compose_constant_null_space!` | `compose_constant_nullspace!` |
| `add_boundary!`, `add_natural_boundary!` | unchanged |
| `create_split_boundary_labels!` | unchanged; deprecated in 0.5.5, removed in 0.6 (a copy is in examples/dmplex/ex69.jl) |
| `set_constants!`, `set_exact_solution!` | unchanged |
| `set_residual!`, `set_jacobian!`, `set_jacobian_preconditioner!` | unchanged |
| `@petsc_residual_fn`, `@petsc_jacobian_fn` | `@residual_fn`, `@jacobian_fn` |
| `@petsc_bd_fn`, `@petsc_simple_fn` | `@bd_fn`, `@simple_fn` |

### `ksp.jl`, `snes.jl`

| v0.4 | v0.5 |
|---|---|
| `PetscKSP`, `PetscSNES` (types) | `KSP`, `SNES` |
| `KSP(…)`, `SNES(…)` (factories) | type constructors |
| `getDM` | `dm` |
| `setDM!` | `set_dm!` |
| `get_solution` | `solution` |
| `gettype` | `type_name` |
| `destroy` | `destroy!` |
| `setcomputeoperators!` | `set_compute_operators!` |
| `setcomputerhs!` | `set_compute_rhs!` |
| `setfunction!` | `set_function!` |
| `setjacobian!` | `set_snes_jacobian!` (see note) |
| `setconvergencetest!` | `set_convergence_test!` (callback first) |
| `Fn_SNESSetConvergenceTest` | `SNESSetConvergenceTestFn`; internal and without the old spelling from 0.6 |
| — | `narrow`, new |
| — | `set_type!(obj, ::Symbol)` for Vec, Mat, KSP, SNES and DM, new: v0.4 had it on TS only |

`setjacobian!` and `dmplex.jl`'s `set_jacobian!` share an English word and nothing else. They stay separate names under the prefix rule, because the callback-first rule puts the callback in argument 1 and so leaves `snes` unable to carry the name.

`setfunction!` and `setjacobian!` each ship two argument orders in v0.4 (`src/snes.jl`), subject-first and callback-first. Only callback-first survives, which is a break with no shim: see [Changes that print no warning](#4.-Changes-that-print-no-warning).

### `ts.jl`

No v0.4 column: PETSc.jl had no high-level TS, so every name here is new and written to these rules from the start.
They are listed because the rename table is the register of what is public, not only of what changed.

| v0.5 | Kind |
|---|---|
| `TS` | type, exported |
| `set_rhs_function!`, `set_rhs_jacobian!` | callback setters, callback first |
| `set_ifunction!`, `set_ijacobian!` | callback setters, likewise |
| `set_monitor!` | callback setter, likewise |
| `set_user_ctx!`, `user_ctx` | writer and reader; deprecated in 0.5.5, removed in 0.6 |
| `set_solution!` | writer |
| `set_time!`, `set_timestep!`, `set_max_time!`, `set_max_steps!` | writers |
| `current_time`, `timestep`, `max_time`, `max_steps` | noun readers |
| `set_tolerances!`, `tolerances` | writer and reader, on `TS`, `KSP` and `SNES`. On a `TS` `tolerances` returns `(; atol, rtol, vatol, vrtol)`, PETSc parameter names, and its two vectors are **borrowed**; on a `KSP` `(; rtol, atol, divtol, max_it)` and on a `SNES` `(; atol, rtol, stol, max_it, max_funcs)` |
| `set_type!`, `type_name` | shared with `ksp.jl` and `snes.jl`, `Symbol`. `type_name` replaces `ts.jl`'s v0.4 spelling `type`, the one name in this file the document did change |
| `set_adapt_type!` | writer, `Symbol`. Prefixed: the target is the `TSAdapt` fetched inside, not the `TS` |
| `set_problem_type!`, `set_exact_final_time!` | writers taking PETSc **enums** (`TSProblemType`, `TSExactFinalTimeOption`), per the enum rule, not `Symbol` |
| `solve!`, `step!`, `interpolate!`, `reset!` | mutation |
| `converged_reason`, `solve_time`, `step_number` | noun readers |
| `snes_iterations`, `ksp_iterations`, `snes_failures` | noun readers, prefixed: the counts belong to the sub-solver |
| `step_rejections` | noun reader, count, belongs to the `TS` itself |
| `snes`, `ksp`, `dm`, `solution` | readers returning **borrowed** handles |
| `destroy!` | mutation, shared |

### `vec.jl`, `mat.jl`

| v0.4 | v0.5 |
|---|---|
| `VecSeq` | `PetscVec` constructor |
| `MatCreateSeqAIJ`, `MatSeqAIJ`, `MatSeqDense`, `MatSeqAIJWithArrays`, `MatAIJ` | `PetscMat` constructor |
| `VecPtr`, `MatPtr` | removed in 0.6: `PetscVec(petsclib, ptr; own)` and `PetscMat(petsclib, ptr; own)` wrap a raw handle, borrowed by default |
| `unsafe_localarray`, `unsafe_local_array` | removed in 0.6: `with_local_array!`, or `acquire_local_array` with `release_local_array` |
| `wrap_localarray` | removed in 0.6, as above |
| `acquire_petsc_local_array` | `acquire_local_array` |
| `release_petsc_local_array` | `release_local_array` |
| `get_petsc_arrays` | `local_arrays`; removed in 0.6, use `with_local_arrays!(f!, g_fx, l_x)` |
| `restore_petsc_arrays` | `restore_local_arrays!`; removed in 0.6 with `local_arrays` |
| `withlocalarray!` | `with_local_array!` |
| `ghostupdate!` | `ghost_update!` |
| `ghostupdatebegin!`, `ghostupdateend!` | `ghost_update_begin!`, `ghost_update_end!` |
| `ownershiprange` | `ownership_range` (drops the positional `base_one`) |
| `setvalues!` | `set_values!` (stays 0-based; `row0idxs`/`col0idxs` become `rows_0b`/`cols_0b`) |
| `addindex!` | `add_index!` |
| `determine_memtype` | `memtype` |
| `as_petsc_vec` | `PetscVec` constructor |
| `array_type`, `memtype_backend` | unchanged |
| `assemble!`, `setup!` | unchanged |
| — | `owns`, new: it answers whether the wrapper owns its handle, which is what `destroy!` consults |
| `destroy` | `destroy!` |

### `init.jl`

| v0.4 | v0.5 |
|---|---|
| `initialized` | `isinitialized` |
| `finalized` | `isfinalized` |
| `check_petsc_wrappers_version` | `check_wrappers_version` (returns a `NamedTuple`) |
| `scalartype`, `inttype` | unchanged (closed list) |
| `initialize`, `finalize` | unchanged, not exported |
| `tao_usable_after_reinitialize` | unchanged (a platform predicate, not a PETSc accessor) |
| — | `equal(A, B)`, new in 0.6.0: PETSc's `MatEqual` and `VecEqual`, since `==` compares values as for Julia arrays |
| — | `with_destroy(f)`, `track!`, `untrack!`, new in 0.5.3: destroy a group of objects together, in reverse order, when a block ends |

### Types

| v0.4 | v0.5 |
|---|---|
| `AbstractPetscKSP`, `AbstractPetscSNES` | `AbstractKSP`, `AbstractSNES` |
| `AbstractPETScMemBackend` | `AbstractPetscMemBackend` |
| `Fn_KSPComputeRHS`, `Fn_KSPComputeOperators` | `KSPComputeRHSFn`, `KSPComputeOperatorsFn`; internal and without the old spellings from 0.6 |
| `Fn_SNESSetFunction`, `Fn_SNESSetJacobian` | `SNESSetFunctionFn`, `SNESSetJacobianFn`; internal and without the old spellings from 0.6 |
| `TSSetRHSFunctionFn`, `TSSetRHSJacobianFn`, `TSSetIFunctionFn`, `TSSetIJacobianFn`, `TSMonitorSetFn` | internal from 0.6 |
| `MatAT` | deleted (no uses; it would have been `MatOrTranspose`) |
| `MatShell`, `MatOp` | unchanged |
| `AbstractPetscDS`, `PetscDS` | unchanged |
| `DMStagGetIndices` | removed (already deprecated in v0.4) |

### `options.jl`, `sys.jl`

| v0.4 | v0.5 |
|---|---|
| `Options` (factory) | `PetscOptions` constructor |
| `getcomm` | `comm` |
| `typedget` | `parse_option` |
| `parse_options` | unchanged |
| — | `with_options(f, obj)` and `with_options(f, opts)`, new in 0.5.3: run `f` with options pushed, callback first, returning the block's result |
| — | `options_at = :construction \| :each_solve \| :first_solve` on the `KSP`, `SNES` and `TS` constructors, new in 0.5.3: a `Symbol` keyword |

`typedget` operates on a plain `NamedTuple`, not on `PetscOptions`, and coerces the looked-up value to the type of the supplied default.
`parse_option` names that, and pairs with `parse_options(args)` in the same file.

### `library.jl`, `audit.jl`

| v0.4 | v0.5 |
|---|---|
| `set_petsclib` | `PetscLibType` constructor |
| `library_info` | returns a `NamedTuple` |
| `audit_petsc_file` | `audit_file` |
| `set_library!`, `unset_library!` | unchanged |

`library_info` printed a report and returned `nothing`, so its name promised data it never handed back.
It now returns the values, so code and tests can read them:

```julia
info = library_info()
info.scalar          # Float64
```
