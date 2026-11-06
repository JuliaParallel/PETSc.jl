# Upgrading from 0.4 and 0.5

PETSc.jl 0.5 is a breaking release. Most of the high-level API has new names, the low-level `LibPETSc` functions return their outputs, and the package needs Julia 1.12 and PETSc 3.25. This page is the order in which to update code written for 0.4. 0.6 removes the 0.5 warnings that name each replacement, so upgrade to 0.5.5 first, follow this page, and then move to 0.6 as described in [Upgrading from 0.5 to 0.6](@ref). Every change is listed in the [Release notes](release_notes.md), and every renamed function in the rename table at the end of [Naming Conventions](naming.md).

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

- **An old high-level name:** the rename table at the end of [Naming Conventions](naming.md).
- **A PETSc C function:** the [C to Julia name index](api_index.md).
- **Why a name is what it is:** the rules in [Naming Conventions](naming.md).

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
