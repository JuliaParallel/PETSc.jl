# API conventions

These rules govern the **high-level interface**: everything reachable as `PETSc.foo`. The low-level layer (`LibPETSc.*`) mirrors the PETSc C API name for name and is out of scope; it keeps its C names so that PETSc's own documentation stays usable. Every addition to the high-level layer follows them, including thin wrappers around a single C call.

## Scope and the tiebreak

| Layer | Namespace | Naming |
|---|---|---|
| Low level | `PETSc.LibPETSc` | Verbatim PETSc C names (`DMStagGetCorners`) |
| High level | `PETSc` | These conventions |

Two goals pull against each other: being idiomatic Julia, and staying recognizable to someone reading the PETSc manual. When they conflict:

> **Julia idiom decides structure. PETSc decides vocabulary.**

Dispatch, `!`, argument order and extending `Base` follow Julia. Which *word* to use follows PETSc. So `DMDAGetDof` becomes `ndofs(dm)`: Julian in shape, PETSc's term for the thing.

### The register

The high-level surface is `names(PETSc; all = true)`, filtered to bindings whose methods are defined in PETSc.jl's own source. Since only types are exported (see [Exports](#Exports)), "unexported" cannot separate the API from internals, so a register does: `scripts/renames.jl` lists every public name and the names that are internal. `scripts/api_surface.jl --check` fails when a binding is in neither.

The register is data rather than prose because three things read it: the `public` declaration in `src/public_names.jl` is generated from it, `scripts/audit.jl` takes its creator and destroyer name sets from it, and `--check` compares against it. Kept by hand, they drift.

`scripts/api_surface.jl --sweeps` lists the functions that break the rules below: those taking `petsclib` first beside a dispatchable object, those returning a `NamedTuple`, the exports, readers returning a PETSc object without a `doc_borrowed` entry, and `!` functions that can return `nothing`. Read a count from the sweep, not from prose.

**Internal helpers** may be renamed freely, with no shim. A leading underscore marks one thing only: an inner worker that shares its name with the wrapper around it, as `Base` does with `_growend!`. `_mul!` is the case here, and `_local_arrays` and `_restore_local_arrays!`, the per-backend workers the GPU extensions extend. It is not a visibility marker.

**Macros** follow the rules for functions: snake_case, and no prefix that repeats the module name (`@residual_fn`, not `@petsc_residual_fn`).

**Base extensions** keep their Base names by definition; [Extending Base and stdlib](#Extending-Base-and-stdlib) governs which ones may exist.

## Case

Default to `snake_case` for any multi-word name. Run words together only for names on the closed list below.

```julia
set_type!(ksp, :gmres)
ghost_corners(dm)
ownership_range(A)
with_local_coordinate_array(f, dm)
```

**Closed run-together list**:

- Base-style predicates: `issimplex`, `issymmetric`, `ishermitian`, `isassembled`, `isinitialized`, `isfinalized`
- PETSc single tokens: `ndofs`, `nextra`, `l2diff`, `nullspace`, `seqaij`, `seqdense`
- Type accessors mirroring `Base.eltype`: `scalartype`, `inttype`, `realtype`

The default requires no judgement and the exceptions are countable. That matters more than matching `Base`, which is itself split between `setindex!` and `set_zero_subnormals`.

Never use `camelCase`. `PascalCase` means "this is a type", without exception (see [Types](#Types)).

## Accessors and setters

Readers are nouns, without PETSc's `Get`: `corners(dm)`, `info(dm)`, `dm(ksp)`, `solution(ksp)`, `comm(v)`. Writers are `set_*!`: `set_dm!(snes, dm)`, `set_type!(ksp, :gmres)`, `set_values!(A, …)`. This targets PETSc's C `Get`/`Set` prefixes; Base's `get(collection, key, default)` idiom is unaffected.

The reader for `XXXGetType` is `type_name`, not `type`, while the setter stays `set_type!`. `type` is the conventional name for a setter's argument, so a `type` reader would be shadowed by its own writers (`set_type!`, `set_adapt_type!` and `set_problem_type!` all take a `type`). And what comes back is the *name* of a PETSc implementation, not a Julia type: `type_name(ksp) === :gmres` says that, where `type(ksp)` invites the reader to expect `typeof`.

### Enums and type names

PETSc has two kinds of named value, and they are represented differently.

**Enumerations** are fixed sets, and the wrapped layer exposes them as Julia enums. Use the enum: `DMStag(petsclib, comm, (DM_BOUNDARY_PERIODIC, DM_BOUNDARY_NONE), …)`.

**Type names** (`KSPType`, `MatType`, `DMType`, `PCType`, …) cannot be enums: PETSc registers them at runtime as strings, so plugins add ones the package has never heard of. They are `Symbol` at the Julia API and `String` at the C boundary:

```julia
set_type!(ksp, :gmres)
type_name(ksp)                             # :gmres, not "gmres"
PetscMat(petsclib, m, n; type = :seqaij)
DMStag(…; stencil_type = :box)
```

A `Symbol` costs nothing to compare, needs no fixed set, and matches `parse_options`, which returns a `NamedTuple` with `Symbol` keys. The conversion to `String` happens once, where the call meets C. A reader that can decline to answer returns `nothing`: `type_name(ts)` before the options have been applied.

### Readers whose names are also variable names

`dm`, `snes` and `ksp` are accessors and also the usual variable names for what they return. Inside the package `dm = dm(ksp)` shadows the function for the rest of the scope, so bind the result to something else: `d = dm(ksp)`, `s = snes(ts)`. Users are unaffected, because none of these is exported and a local binding cannot shadow `PETSc.dm`.

### What an accessor hands back

A reader returning a PETSc object returns a **borrowed handle**: the object belongs to whatever it was asked of, no finalizer is attached, and `destroy!` on it does nothing. `solution(ts)`, `snes(ts)`, `dm(ksp)` and the vectors of `tolerances(ts)` are all borrowed. Constructors return owned objects (see [Ownership](#Ownership)).

The docstring says so through `doc_borrowed`, a helper beside `doc_external`, so that the sweep can check it: a reader returning a PETSc object with no `doc_borrowed` entry is reported unless it is listed as owning its result (the constructors, and `narrow`).

## Object prefixes

Drop `dm_`, `mat_`, `vec_`, `plex_` when the first argument already carries the type; dispatch does the work: `project_function!(vec, dm, f)`, `set_nullspace!(A, ns)`, `distribute(dm)`.

Keep a prefix only for free functions, where no argument identifies the domain. These take `petsclib` and `comm` rather than a PETSc object: `fe_create_default(petsclib, comm, …)`, `mat_nullspace_create(petsclib, comm; has_const = true)`. Check this against the argument list, not the name: `save_vtk(vec, filename)` takes a `PetscVec`, so it is not a free function.

### When the subject is not the first argument

> **Keep a prefix when the object that identifies the operation is not the first argument**, either because a callback occupies that slot (see [Callbacks come first](#Callbacks-come-first)) or because the target is a sub-object reached inside the call.

```julia
set_snes_jacobian!(updateJ!, snes, J, PJ)   # callback first, so snes is not argument 1
set_adapt_type!(ts, :basic)                 # the target is the TSAdapt inside, not the TS
snes_iterations(ts)                         # the count belongs to the sub-solver
```

The alternative for a sub-object is a reader with no prefix. A sub-object that has a type of its own takes that form: `set_type!(pc(ksp), :jacobi)`, not `set_pc_type!(ksp, :jacobi)`. `set_adapt_type!` keeps its prefix because `TSAdapt` has no type.

`set_snes_jacobian!` and `PetscDS`'s `set_jacobian!` would not collide by arity, but merging a callback registration with seven function pointers behind one generic function hides two unrelated operations under one name.

## Types

> **PascalCase always names a type. Construction always goes through the type.**

### Prefixes

Use PETSc's bare class name. Add a `Petsc` prefix only when the bare name is too generic to export safely:

```julia
# prefixed
PetscVec  PetscMat  PetscOptions  PetscDS
PetscDM   # low-level handle only

# bare
DMDA  DMStag  DMPlex
KSP  SNES  TS  Tao  IS  AO  PF
```

`Vec` and `Mat` stay prefixed: they are exported, and bare `Vec` and `Mat` would clash across the ecosystem.

### The DM type is a Julia type

The DM type and, where methods depend on it, the dimension are type parameters, so ordinary dispatch applies and each method has one concrete return type:

```julia
abstract type AbstractPetscDM{PetscLib} end

DMDA{PetscLib,N}   <: AbstractPetscDM{PetscLib}
DMStag{PetscLib,N} <: AbstractPetscDM{PetscLib}
DMPlex{PetscLib}   <: AbstractPetscDM{PetscLib}

corners(dm::DMDA{L,N})   where {L,N} = …   # (lower, upper, size)
corners(dm::DMStag{L,N}) where {L,N} = …   # (lower, upper, size, nextra)
```

So there are **no type suffixes** on function names: one `corners`, with a method per DM type, never `corners_dmstag`.

A property becomes a type parameter only where a method dispatches on it or a return type is shaped by it. The DM type qualifies, because `corners`, `local_indices` and the coordinate handling differ by it. The dimension qualifies for DMDA and DMStag, whose corners and creation are written per dimension, but not for DMPlex: nothing dispatches on it, and `DMPlex(petsclib, comm)` leaves it unset until setup, so a plex reports it through `ndims(dm)`. The Vec and Mat type only affects construction, so it stays a keyword (see [Keyword versus positional](#Keyword-versus-positional)).

### DMs of unknown provenance

`PetscDM{PetscLib}` is the low-level handle, since `LibPETSc.DMCreate` and friends must return something before the DM type is known. It is not part of the high-level API.

A DM obtained from another object (`dm(ksp)`, a `DMPlex` read from file) has its type and dimension known only at runtime. `narrow(dm)` queries `DMGetType` and `DMGetDimension` and returns the typed handle; accessors returning a DM PETSc owns call it themselves. Its return type is a `Union` of more than nine members, past `MAX_UNION_SPLITTING`, so the call is a dynamic dispatch. One dispatch is cheap; an abstractly typed DM carried through a hot loop is not, so narrow once behind a function barrier:

```julia
function work(ksp)
    d = dm(ksp)     # dynamic dispatch here, once
    _work(d)        # specialized on the concrete type
end
```

`narrow` returns a second handle onto the same PETSc object, borrowed. A DM type with no Julia type of its own comes back unchanged.

### Abstract, callback and wrapper types

**Abstract types** are `Abstract` followed by the concrete name [Prefixes](#Prefixes) gives, so the prefix is decided once: `AbstractPetscVec`, `AbstractPetscMat`, `AbstractPetscDM`, `AbstractKSP`, `AbstractPC`, `AbstractSNES`, `AbstractTS`.

**Callback types** take an `Fn` suffix, following the generated layer's `DMDATSRHSFunctionLocalFn`: `KSPComputeRHSFn`, `SNESSetFunctionFn`. They are the trampolines the callback setters hand to `@cfunction`, and internal (see [Public surface](#Public-surface)).

**Wrapper and alias types** keep bare PETSc names. The prefix rule asks whether *this name* is too generic to export, not whether its class is: `Mat` needs the prefix, `MatShell` does not.

## Constructors

Construction goes through the type, dispatching on the argument types:

```julia
A = PetscMat(petsclib, S::SparseMatrixCSC)
A = PetscMat(petsclib, m, n, nnz)
A = PetscMat(petsclib, m, n; type = :dense)
A = PetscMat(petsclib, rowptr, colval, nzval)

v = PetscVec(petsclib, 10)
v = PetscVec(petsclib, x::Vector)
```

`petsclib` appears in constructors because no object exists yet to carry it. Elsewhere it appears only in free functions (see [Argument order](#Argument-order)).

## Mutation

A function gets `!` if and only if it mutates one of its arguments, **or** package-global state. Configuring an opaque PETSc object counts as mutation:

```julia
corners(dm)              # reads
ndofs(dm)                # reads

destroy!(dm)
assemble!(A)
set_type!(ksp, :gmres)
solve!(x, ksp, b)

set_library!(path)       # mutates package state, no argument mutated
```

The global-state clause keeps the bang on `set_library!` and `unset_library!`, following `Random.seed!`.

A `!` function returns the object it mutates, as `push!`, `fill!`, `mul!` and `copyto!` do: `solve!(x, ksp, b)` returns `x`, `set_type!(ksp, :gmres)` returns `ksp`, and `set_function!(f!, snes, r)` returns `snes`, the object being configured, even though a callback precedes it. Adding to an object returns the object too: `add_boundary!(dm, …)` returns `dm`. The exceptions:

- A function that releases what it is given returns `nothing`, like `close`: `destroy!`.
- A function that mutates package state only returns `nothing`: `set_library!`, and `push!` and `pop!` on a `PetscOptions`, which change PETSc's global options stack.
- A function whose first argument is a `do` block returns the block's result, like `open(f, file)`: `with_local_array!`, `with_destroy`.
- A function that files an object away for later returns that object, like `get!`, so it can wrap a constructor: `track!(list, obj)` returns `obj`.

A function that writes a file and mutates no argument takes no `!`, like `write` and FileIO's `save`: `save_vtk`.

## Argument order

Subject first, and mutated arguments before read-only ones, matching `mul!(C, A, B)` and `copyto!(dest, src)`:

```julia
solve!(x, ksp, b)                  # x is written
project_function!(vec, dm, f)      # vec is written
global_to_local!(lvec, dm, gvec)   # lvec is written
corners(dm)                        # readers: subject first
```

`petsclib` never leads a high-level call: the object carries `PetscLib` as a type parameter, so passing it again is redundant. The exceptions are constructors and the free functions of [Object prefixes](#Object-prefixes).

### Callbacks come first

A function-valued argument goes first, ahead of the subject, so `do` syntax works. This outranks subject-first, and follows `map(f, c)` and `open(f, path)`:

```julia
set_function!(snes, x) do f, x
    …
end
```

The rule is for Julia functions. A C function pointer from `@cfunction` cannot be a `do` block, so a setter that takes only pointers keeps the subject first: `set_residual!(ds, field, f0, f1)`, `add_coarsen_hook!(dm, hook)`. The order of a callback's own arguments is in [Callback signatures](#Callback-signatures).

## Keyword versus positional

Positional arguments are the ones the object cannot exist without. Keyword arguments are the ones with a sensible default:

```julia
DMStag(petsclib, comm, boundary, dims, dof; stencil_width = 1, stencil_type = :box)
PetscMat(petsclib, m, n; type = :dense)                       # defaults to :seqaij
set_values!(A, rows_0b, cols_0b, vals; mode = ADD_VALUES)     # InsertMode and ScatterMode
ghost_update!(v; mode = ADD_VALUES, scatter = SCATTER_REVERSE)
```

Never take a keyword argument that changes the return type. `type = :dense` is fine because every variant returns a `PetscMat`; a keyword selecting the DM type would not be, which is why that is a type parameter.

## Abbreviations

Follow PETSc's own abbreviations, so a user reading the PETSc manual guesses the Julia name: `ndofs` (`DMDAGetDof`), `l2diff` (`DMComputeL2Diff`), `nextra` (`DMStagGetCorners`), and `seqaij`, `dense`, `ds`, `fe`, `is`, `pc`, `ksp`, `snes`, `ts`. Do not invent abbreviations PETSc does not use, and do not expand ones it does.

## Extending Base and stdlib

Add methods to `Base` and `LinearAlgebra` wherever the PETSc operation satisfies the existing contract, so that generic Julia code works on PETSc objects:

```julia
Base.size, Base.length, Base.ndims, Base.eltype, Base.axes
Base.similar, Base.copy, Base.fill!, Base.iterate, Base.copyto!
Base.getindex, Base.setindex!, Base.show
Base.:(==), Base.isequal, Base.sum, Base.maximum, Base.minimum,
Base.any, Base.all, Base.count, Base.mapreduce, Base.collect, Base.isapprox
Base.:+, Base.:-, Base.:*, Base.:/
LinearAlgebra.norm, LinearAlgebra.mul!, LinearAlgebra.dot,
LinearAlgebra.axpy!, LinearAlgebra.axpby!, LinearAlgebra.rmul!, LinearAlgebra.lmul!,
LinearAlgebra.issymmetric, LinearAlgebra.ishermitian
```

Do not shadow a Base name with something unrelated: the TS reader is `current_time`, because `Base.time` exists.

### Blocklist

Do **not** extend these. The word matches but the contract does not, and a silent mismatch is worse than an unfamiliar name:

| Name | Why not | Use instead |
|---|---|---|
| `LinearAlgebra.nullspace` | Returns a basis matrix; PETSc's `MatNullSpace` is a solver hint object | `set_nullspace!`, `mat_nullspace_create` |
| `Base.values` | PETSc "values" are a setter concept (`MatSetValues`), not a view of contents | `v[:]`, `A[:, :]`, `collect` |
| `Base.view` | `PetscViewer` writes an object to a stream; it is not an array view | `Base.show`, or `LibPETSc.VecView` and `MatView` |
| `Base.setfield!` | Core builtin that writes a struct field; PETSc's attaches a finite element to a DM | `set_field!` |
| `Base.empty` | Returns an empty collection of the same type; a DM is not a collection | the type constructor |

An addition needs the same three columns: the name, the contract mismatch, and the replacement.

## Return values

Functions returning several related values return a `NamedTuple`. **Index-shaped returns** (corners, ghost corners, local and global index ranges) take their field names from a closed vocabulary:

| Field | Type | Meaning |
|---|---|---|
| `lower`, `upper` | `CartesianIndex{N}` | 1-based, inclusive |
| `size` | `NTuple{N,Int}` | Local extent |
| `center`, `vertex` | `NamedTuple` of `UnitRange` | Staggered index ranges, keyed `x`, `y`, `z` |
| `x`, `y`, `z` | `UnitRange` | Per-axis range inside `center`/`vertex` |
| `nextra` | `NTuple{N,Int}` | DMStag partial elements |

Extend the vocabulary rather than inventing a synonym: a local extent is `size`, never `dims` or `extent`.

Every other `NamedTuple` takes its field names from the PETSc parameters it reports, in `snake_case`, one field per value: `tolerances(ts)` is `(; atol, rtol, vatol, vrtol)` after `TSGetTolerances`. Where PETSc spells one parameter differently across types, the field takes the name of its option without the prefix, so a keyword means the same everywhere: `tolerances(ksp)` is `(; rtol, atol, divtol, max_it)` after `-ksp_divtol` and `-ksp_max_it`, although `KSPGetTolerances` says `dtol` and `maxits`.

Results are **dimension-correct**: a 2D DM returns `CartesianIndex{2}`, 2-tuples, and `center` and `vertex` keyed `(x = …, y = …)` with no `z`. Build them with `ntuple(…, Val(N))` from the type parameter and index the axis names by it, never by splatting a collection: splatting hides the length from the compiler, so the fields would infer as `Any` and every call would allocate. Written this way `local_indices(dm2d).center` infers as `@NamedTuple{x::UnitRange{Int}, y::UnitRange{Int}}`.

Two returns stay abstract by design: `corners` has three fields on a `DMDA` and four on a `DMStag`, so code generic over `AbstractPetscDM` sees a small union, and `dm(ksp)` returns a wide one (see [DMs of unknown provenance](#DMs-of-unknown-provenance)).

### Index base

> An index the high-level layer accepts or returns **into Julia data** is 1-based. `LibPETSc` is 0-based.

The rule covers positions in an array and index ranges over one: `corners`, `ghost_corners`, `local_indices`, `global_indices`, `ownership_range`, `dof_slot`. It does not cover identifiers PETSc issues and expects back unchanged, such as the field number `set_field!` gave a field, or a label value passed to `add_boundary!`. Those keep PETSc's numbering, so they agree with every PETSc example.

**Bulk index arrays keep PETSc's base.** The scalar and slice interface is 1-based, because `Base.getindex` and `Base.setindex!` have no choice. A vector of indices handed to C, or handed back from C, keeps the base C uses and says so in its name:

| Function | Base | Naming |
|---|---|---|
| `set_values!` rows and columns (`MatSetValues`, `MatSetValuesStencil`) | 0-based | `rows_0b`, `cols_0b` |
| `star_fd_coloring` returned index vectors | mixed | `row_coo_local_0b`, `perturb_cols_1b`, … |
| `A[i, j] = v`, `A[I::CartesianIndex, J] = v` | 1-based | Base's contract, converts internally |

`set_values!`'s stencil method takes `MatStencil` structs whose layout must match PETSc's C struct, so a 1-based form would rebuild the caller's vector on a hot path. `A[i, j] = v` is the 1-based route, and `set_values!` the bulk one for code that already holds PETSc indices.

## Exports

Only types and construction entry points are exported. Verbs and accessors stay qualified:

```julia
export LibPETSc
export DMDA, DMStag, DMPlex
export PetscVec, PetscMat, PetscOptions
export KSP, SNES, TS
export petsclibs
```

`petsclibs` is exported because every entry point takes one, and `LibPETSc` because the low-level layer is reached through it. `initialize` and `finalize` are not exported: `Base.finalize` exists, and shadowing it would be worse than typing `PETSc.finalize`.

Two reasons to stay narrow. Exporting `solve!` would make `using PETSc, LinearSolve` ambiguous, a common combination. And adding an export later is not breaking, while removing one is.

### The `public` keyword

`public` marks a name as API without exporting it, which is the distinction a type-only export list needs. The declaration is `src/public_names.jl`, generated from `scripts/renames.jl` minus the exported names, since a name cannot be both. `names(PETSc)` then reports the API directly, which matters for the short accessors (`dm`, `comm`, `info`, `ds`, `label`, `solution`) that a text search cannot verify.

## Errors

Argument problems raise standard Julia exceptions. `@assert` is reserved for invariants that cannot fail unless the package itself is wrong, and never validates user input: it carries no useful message and is not guaranteed to run.

| Condition | Exception |
|---|---|
| Wrong argument value or type | `ArgumentError` |
| Mismatched sizes or shapes | `DimensionMismatch` |
| Library not initialized | `PetscNotInitialized` |
| Error returned by PETSc itself | `PetscError` |
| Solver stopped without converging, unless `check = false` | `ConvergenceError` |

Every high-level constructor that creates a PETSc object checks that its library is initialized before calling PETSc, so a missing `initialize` is reported as `PetscNotInitialized` and not as the `PetscError` PETSc would raise; `test/core/errors.jl` checks each of them. `PetscOptions` is included although `PetscOptionsCreate` works before `PetscInitialize`: the object records the library's `age`, which `initialize` advances, so an options database created first would already be stale and `destroy!` would skip it.

```julia
length(A) == prod(sz) ||
    throw(DimensionMismatch("array has $(length(A)) entries, expected $(prod(sz))"))
isinitialized(lib) || throw(PetscNotInitialized(lib))
```

## Docstrings

A renamed wrapper costs discoverability: a user who knows `PetscFECopyQuadrature` cannot grep for `copy_quadrature!`. So every high-level docstring links the C function it wraps through `doc_external`:

```julia
"""
    corners(dm::DMStag)

…

# External Links
$(doc_external("DMSTAG/DMStagGetCorners"))
"""
```

The [C to Julia name index](api_index.md) is generated from those entries.

### Say what goes wrong

Where a failure mode is known, the docstring states it. A link to the C page says what a function does, rarely what it does when misused, and PETSc is a library where misuse often succeeds quietly. `ksp(ts)` says it fails unless the problem was declared `TS_LINEAR`, `snes(ts)` that an explicit method gets an unused solver, and `type_name(ts)` that it answers `nothing` until the options are applied. This is a documentation rule, not a CI rule: whether a failure mode is known cannot be checked mechanically.

## Changing the API

A rename keeps a forwarding shim until the next breaking release, which removes it. Old and new names are recorded in `scripts/renames.jl`.

A shim calls `@warn … maxlog = 1`, not `Base.@deprecate`: `@deprecate` prints nothing unless Julia runs with `--depwarn=yes`, and the default is `no`, so a user would see silence and then a `MethodError`.

```julia
function getcorners(args...; kwargs...)
    @warn "getcorners is deprecated, use corners" maxlog = 1
    return corners(args...; kwargs...)
end
```

A shim forwards only when the old call and the new one mean the same thing:

| Kind | Shim |
|---|---|
| Name changed, meaning unchanged | forwards, warns once |
| Arguments reordered, types differ | forwards, warns once |
| Argument dropped | keeps the argument, warns when it is passed |
| Argument types unchanged, **meaning changed** | throws, naming the replacement |

In the last row dispatch cannot tell the old call from the new one, and a forwarding shim corrupts data without a sound. Reconsider such a change first; if it survives, it throws rather than forwards.

Removing a shim must not turn a call into a silent no-op. The DM constructors pass unknown keywords on as PETSc options, so a removed keyword throws an `ArgumentError` instead of becoming an option nothing reads.

The test suite uses only the current names, every shim gets a test that it warns and forwards, and `test/core/api_surface.jl` checks that no old name in the register is still defined.

## Handles, callbacks and their state

The sections above decide what things are called. This one decides how a wrapper relates to the PETSc object behind it, and what happens when PETSc calls back into Julia. The tiebreak mirrors the first one:

> **Behavior follows what a Julia user expects of a Julia package. PETSc decides what is possible.**

An object is freed by the garbage collector, an exception thrown in a callback comes out of the call that ran it, and every object of a kind behaves like every other. Where PETSc's model rules that out, PETSc wins and the docstring says why: destruction on more than one process is collective, so there it stays explicit.

### A wrapper is a view of the PETSc object

Several wrappers can point at one PETSc object: `pc(ksp)` builds a new `PC` on every read, and `snes(ts)` a new `SNES`. Anything stored on a wrapper is lost with that wrapper, while the object lives on. So a wrapper holds only what identifies the object and who may destroy it:

| Field | Meaning |
|---|---|
| `ptr` | the PETSc object |
| `age` | the initialize/finalize cycle it was created in (`isdestroyable`) |
| `own` | whether `destroy!` on this wrapper destroys the object |

Everything else Julia keeps about an object lives with the PETSc object (see [Callback state lives with the PETSc object](#Callback-state-lives-with-the-PETSc-object)).

### Ownership

Every constructor returns `own = true` and, on a one-process communicator, attaches a finalizer, so the object is freed when it is collected. On more than one process `XDestroy` is collective and a finalizer runs at a different moment on each rank, so there the caller calls `destroy!`. Readers return `own = false` and attach nothing, and `destroy!` on what they return does nothing. `owns(obj)` answers the field. `with_destroy` keeps that bookkeeping for a group of objects: it destroys every object `track!` registered, in reverse order, when its block ends.

`LibPETSc` follows PETSc's own convention, which the function name carries. A handle returned by a function with `Get` in its name is a reference its owner keeps, so the wrapper has `own = false`: `KSPGetPC`, `TSGetSNES`, `DMGetCoordinates`. Every other function that returns a handle (`Create`, `Duplicate`, `Clone`, `Convert`, `DMPlexDistribute`, …) hands out a new reference, and the wrapper has `own = true`. The few `Get` functions that also hand out a new reference, such as `MatGetFactor`, `DMLabelGetStratumIS` and `MatGetOrdering`, are listed in `wrapping/generator/rules/ownership.toml`. A function whose ownership is unclear stays off that list: a borrowed handle the caller should have destroyed only leaks, while an owned one that was not the caller's would be destroyed twice.

`own` governs `destroy!` only. The generated `XDestroy` frees whatever it is given, as in C, and sets the wrapper's `ptr` to `C_NULL`, so freeing a high-level object through `LibPETSc` and letting its finalizer run afterwards is safe.

Objects passed to a callback are borrowed, and valid only for the duration of the call.

### Callback state lives with the PETSc object

Each PETSc object that Julia calls back from has one state box, held by a `PetscContainer` composed onto the object. PETSc destroys the container together with the object, and that ends the box. It follows that:

- A callback set through any wrapper lives as long as the PETSc object, whichever wrapper was used and whether it is still reachable: `set_function!(f!, snes(ts), r)` is safe.
- Setting a callback again replaces the previous one, and the old closure can be collected.
- `finalize` drops the boxes of its library, since none of its objects can call back afterwards.
- The context pointer PETSc passes to a trampoline is the box, never a wrapper.

The box is internal, and there is no separate user context: a callback's closure captures what it needs.

### Exceptions in callbacks

A Julia exception must not unwind through PETSc's C frames. Every trampoline runs the user's function under a guard that catches the exception, stores it, and returns an error code, so that PETSc unwinds its own stack. The high-level call that started the work (`solve!`, `step!`, `setup!`) then rethrows the original exception rather than a `PetscError`: a `DomainError` thrown in a residual comes out of `solve!` as that `DomainError`. The backtrace shown starts at the rethrow; the frames inside the callback are logged at debug level (`JULIA_DEBUG=PETSc`).

When the work was started through `LibPETSc` directly, no high-level call is there to rethrow. The exception is logged with its backtrace, and the `LibPETSc` call throws `PetscError`.

The DMPlex pointwise functions (`@residual_fn` and the other macros) are exempt. They are plain C-callable functions, called once per quadrature point, with no Julia state, and they must not throw.

### Callback signatures

- The callback is the setter's first argument (see [Callbacks come first](#Callbacks-come-first)).
- A callback that fills an argument takes that argument first, then the object, then what it reads: `f!(F, snes, x)`, `updateJ!(J, P, snes, x)`, `apply!(y, pc, x)`, `f!(F, ts, t, u)`. A callback that only observes takes the object first: `monitor(ts, step, t, u)`, `test!(snes, it, xnorm, gnorm, fnorm)`.
- A callback succeeds by returning and fails by throwing. Its return value is ignored, unless it carries a meaning its setter's docstring states, such as the convergence reason `test!` returns.
- Each callback has one calling form, with no trailing user context: a closure captures what it needs.

### Public surface

A name is public when users write it. Trampoline types (the `…Fn` callables handed to `@cfunction`), state boxes and their helpers are internal, and registered as such in `scripts/renames.jl`. They stay reachable as `PETSc.SNESSetFunctionFn` and the like, without a stability promise.
